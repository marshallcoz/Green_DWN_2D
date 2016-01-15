program main
! Calculo de la funcion de Green de desplazamientos y esfuerzos en 2D
! dada una fuente puntual en un semiespacio con estratos planos.
! El calculo se hace con el metodo de la matriz global en cada
! numero de onda (DWN).

! main.F90
! variables.f90    :  definir modelo y receptores
! Green.F90        :  calculo de funcion de Green en frecuencia

! Modulos con definiciones del tipo de varible
use modelVars, only : nfrec,nmax,NPTSTIME,come ! Definicion del modelo
use sourceVars, only : currentiFte,nFuentes,skipdir
use matrizGlobal
implicit none
integer :: J,dir
logical :: saveASmatlab,saveAstxt,sheetOfReceivers

print*,"";print*,"";print*,"";print*,"";print*,"";print*,""
print*,"  Funcion de Green de Desplazamientos y Tracciones  "
print*," para receptores en semiespacio con estratos planos "

call setupmodel !leer del archivo input.txt
call checarWisdom(2*nfrec,2*nmax,NPTSTIME) ! prepar FFTw
call predimensionar !variables de salida y de trabajo
call preparePointerTable !reconocer receptores a misma profundidad
DO J = NFREC+1,1,-1
    call preProcessAtThisfrec(J)

!···· campo en medio estratificado ····································
    do currentiFte = 1,nFuentes !para cada una de las incidencias
       call makeGANU0 !gamma y nu con el número de onda de ésta incidencia plana

       do dir= 1,3 !x,y,z direction of force application
        if(dir .eq. 2) then ! SH
           if(skipdir(dir)) cycle
        else ! P_SV
           if(skipdir(1) .and. skipdir(3)) cycle
        end if! dir
        if(.not. skipdir(dir)) then
              call diffField_at_iz(0,dir,J,cOME)
        end if
       end do !dir
     end do ! currentiFte
end do ! loop de frecuencia

saveASmatlab = .false.
saveAstxt = .false.
sheetOfReceivers = .true.
call guardarResultados(saveASmatlab,saveAstxt,sheetOfReceivers)

call vaciarWisdom
print*,"done"
end program main
