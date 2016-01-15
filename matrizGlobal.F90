module matrizGlobal
public :: deltaij,porLoMenosUnoEsEstr,makeGANU0,diffField_at_iz
contains

subroutine diffField_at_iz(i_zF,dir_j,J,cOME_in)
!#define ver 1
      ! esta función es llamada con cada profundidad donde hay
      ! por lo menos una fuente.
      use gloVars, only: z0,UI,UR,PWfrecReal!,plotFKS
      use resultVars, only : pota,Punto,nZs,MecaElem,FFres
      use refSolMatrixVars, only : B,Ak
      use modelVars, only : NMAX,k_vec,dk,vecNK,SpliK,OME
      use wavelets
      !use dislin
      use sourceVars, only: Po,iFte=>currentiFte
      use modelVars, only:N,Z,alfa0,beta0,alfa,beta
      use, intrinsic :: iso_c_binding
      use debugStuff
      implicit none
      integer, intent(in) :: i_zF,dir_j,J
      complex*16, intent(in),target  :: cOME_in
      logical,pointer :: intf
      integer, pointer :: ef
      real*8,target :: k
      real*8, pointer :: zf,xf,pt_k
      complex*16, dimension(:,:), allocatable, target :: auxK,savedAuxK
      complex*16, target  :: cOME,alf,bet
      complex*16, pointer :: pt_cOME_i
      integer, dimension(:),allocatable,target :: ipivA
      integer, dimension(:),pointer :: pt_ipivA
      complex*16, dimension(:),allocatable,target :: workA
      complex*16, dimension(:),pointer :: pt_workA
      complex*16,dimension(:,:),pointer :: pointA
      type(Punto), pointer :: pXi,p_X
      logical :: auxLogic ,isPW
      integer :: ik,tam,itabla_z,itabla_x,iMec,mecS,mecE,&
                 nXis,n_Xs,iXi,dj,pos,ne!porLoMenosUnoEsEstr
#ifdef ver
      character(LEN=32) :: arg
      real :: result,lastresult
      real, dimension(2) :: tarray
#endif
      allocate(auxK(2*nmax,5)); allocate(savedAuxK(2*nmax,5))
      isPW = .false.
      if (i_zF .eq. 0) then
      if (Po(iFte)%tipofuente .eq. 1) then
      isPW = .true.
      end if
      end if

      cOME = cOME_in
      dj = dir_j; if(dj .eq. 3) dj = 2
      if (i_zF .eq. 0) then
         itabla_x = 3 ! En la tabla (pota) de índices: la fuente real -> (0,3)
!        i_FuenteFinal = 1 ! Cantidad de fuentes reales
         if (isPW) then ! onda plana incidente·p
            ! con incidencia de onda plana no usamos atenuación
           if (PWfrecReal) then
             cOME = OME * UR
             alf = alfa0(N+1)
             bet = beta0(N+1)
           else
             cOME = cOME_in
             alf = alfa(N+1)
             bet = beta(N+1)
           end if
                                      k = real(cOME/bet)*sin(Po(iFte)%gamma) !SV,SH
           if(Po(iFte)%PW_pol .eq. 2) k = real(cOME/alf)*sin(Po(iFte)%gamma) ! P
         end if! ·································································n
      else; itabla_x = 2 + pota(i_zF,1) + 1 !    la primer fuente virtual
      end if       !        a esa profundidad
      ! ------- para cada una de las fuentes en la tabla -------------------------
       call asociar(pXi,iFte,i_zF,itabla_x) ! asociar apuntador a fuente [pXi]
       xf=>pXi%center%x;zf=>pXi%center%z;ef=>pXi%layer;intf=>pXi%isOnInterface
#ifdef ver
       call ETIME(tarray, result)
       print*,"pXi: x",pXi%center%x,"z",pXi%center%z,&
       "isOnInterface",intf, result
       lastresult = result
#endif
      ! Si es la fuente real y es una onda plana no se usa el DWN. Se calcula para
      ! el número de onda horizontal asociado al ángulo de incidencia ············
      if (isPW) then ! onda plana incidente   ·
          if (dir_j .eq. 2) then! SH                                             ·o
            tam = 2*N+1; if (Z(0) .lt. 0.0) tam = tam + 1!                       ·n
            pointA => Ak(1:tam,1:tam,0) !indice 0 reservado para onda plana      ·d
            pt_k => k; pt_come_i => cOME!                                        ·a
            allocate(ipivA(tam)); allocate(workA((tam)*(tam)))!                  ·
            pt_ipivA => ipivA; pt_workA => workA!                                ·
            call globalmatrix_SH(pointA,pt_k,0)!                                 ·p
            call inverseA(pointA,pt_ipivA,pt_workA,tam)!                         ·l
            call SHvectorB_ondaplana(B(:,0),pxi%gamma)!                          ·a
            B(:,0) = matmul(Ak(1:tam,1:tam,0),B(:,0))!                           ·n
          else!  P-SV                                                            ·a
            tam = 4*N+2; if (Z(0) .lt. 0.0) tam = tam + 2!                       ·
            pointA => Ak(1:tam,1:tam,0) !indice 0 reservado para onda plana      ·
            pt_k => k; pt_come_i => cOME!                                        ·
            allocate(ipivA(tam)); allocate(workA((tam)*(tam)))!                  ·
            pt_ipivA => ipivA; pt_workA => workA!                                ·
            call gloMat_PSV(pointA,pt_k,0)!                                      ·
            call inverseA(pointA,pt_ipivA,pt_workA,tam)!                         ·
            call PSVvectorB_ondaplana(B(:,0),pxi%gamma)!                         ·
            B(:,0) = matmul(Ak(1:tam,1:tam,0),B(:,0))!                           ·
          end if!                                                                ·
          pos = 0; ne = 2*nmax+1
      else ! · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·
      ! La fuente es cilíndrica, puede tratarse de la fuente real o una virtual ··
        pos = min(int(vecNK(J)*SpliK),nmax); ne = 2*nmax-(pos-2)                    !
        B(:,pos:ne) = 0                                                           !
        if (dir_j .eq. 2) then                                                   !
         tam = 2*N+1; if (Z(0) .lt. 0.0) tam = tam + 1                           !o
!        if ((i_zF .eq. 0) .and. (pXi%region .eq. 2)) return

         ! si la fuente está sobre una interfaz ...............
!        if (pXi%isOnInterface) then
         do ik = 1,pos+1                                                            !n
           call SHvectorB_force(i_zF,B(:,ik),tam,pXi,cOME,k_vec(ik),ik)             !d
           B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !a
         end do                                                                  !
         do ik = ne,2*NMAX                                                       !c
           call SHvectorB_force(i_zF,B(:,ik),tam,pXi,cOME,k_vec(ik),ik)             !i
           B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !l
         end do                                                                  !i
!        else ! la fuente está entre interfaces ..............
!        stop "fuerza entre interfaces no implementado para SH"
!          do ik = 1,pos+1
!            stop 1990
!            call eGAeNU(i_zF,ik,pXI,dj)
!          end do!
!          do ik = ne,2*NMAX
!            call eGAeNU(i_zF,ik,pXI,dj)
!          end do
!        end if
        else!  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .     !n
         tam = 4*N+2; if (Z(0) .lt. 0.0) tam = tam + 2                           !d
!         if ((i_zF .eq. 0) .and. (pXi%region .eq. 2)) return

         ! si la fuente está sobre una interfaz ...............
         if (pXi%isOnInterface) then
           do ik = 1,pos+1                                                         !r
             call PSVvectorB_force(i_zF,B(:,ik),tam,pXi,dir_j,cOME,k_vec(ik),ik)   !i
!            stop 2046
             B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !c
           end do                                                                  !a
           do ik = ne,2*NMAX                                                       !
             call PSVvectorB_force(i_zF,B(:,ik),tam,pXi,dir_j,cOME,k_vec(ik),ik)   !
             B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !
           end do
         else ! la fuente está entre interfaces ..............
           do ik = 1,pos+1
             call eGAeNU(i_zF,ik,pXI,dj)
           end do!
           do ik = ne,2*NMAX
             call eGAeNU(i_zF,ik,pXI,dj)
           end do
         end if
        end if                                                                   !
      end if! · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·
#ifdef ver
       call ETIME(tarray, result)
       print*,"vector B,  x=inv(A)B (todas k) ",result-lastresult
       lastresult = result
       print*,"resultados para cada Z donde hay receptores:"
#endif
      ! resultados para cada Z donde hay receptores en el medio estratificado ...
      do itabla_z = 1,nZs
      !(en la tabla todos los puntos son receptores)
       auxLogic = porLoMenosUnoEsEstr(itabla_z)
#ifdef ver
       print*,"  itabla_z",itabla_z,"de",nZs," -----------------------------"
       print*,"   pota(itabla_z,:)=",pota(itabla_z,:),auxLogic
#endif
       !
       if (auxLogic .eqv. .false.) cycle
        ! (el primer receptor de la tabla --------| a esa profundidad)
            call asociar(p_X, 0, itabla_z, 3)
#ifdef ver
       call ETIME(tarray, result)
       print*,"  p_x%center%z",p_x%center%z,"[m]  e=",p_X%layer
       lastresult = result
#endif
        !
            if (p_X%layer .eq. N+2) cycle
            if (dir_j .eq. 2) then; mecS = 1; mecE = 3 !V,s32,s12
            else;                   mecS = 1; mecE = 5 !W,U,s33,s31,s11
            end if

      ! ... elementos mecánicos a la profundidad p_X%center%z ........
      ! .... usando los coeficientes de las ondas en el estrato ......
!     savedauxk = z0
!     savedauxk(po+1:ne-1,:) = 0
      if (isPW) then ! onda plana·············
          if (dir_j .eq. 2) then!                                               ·
                 savedauxK(1,1:3) = SHdiffByStrata(B(:,0), &
                              p_X%center%z, p_X%layer,cOME,k,0)
          else !                                                                ·
                 savedauxk(1,1:5) = PSVdiffByStrata(B(:,0), &!                  ·a
                              p_X%center%z, p_X%layer,cOME,k,0)!                ·a
          end if!                                                               ·
      else ! onda plana incidente / onda cilíndrica circular ····················
        do ik = 1,pos+1                                                         !
          if (dir_j .eq. 2) then !SH
              savedauxK(ik,1:3) = SHdiffByStrata(B(:,ik), &
                              p_X%center%z, p_X%layer,cOME,k_vec(ik),ik)
          else !PSV                                                             !a
              savedauxk(ik,1:5) = PSVdiffByStrata(B(:,ik), &                    !c
                              p_X%center%z, p_X%layer,cOME,k_vec(ik),ik)        !i
          end if                                                                !l
        end do ! ik                                                             !i
        do ik = ne,2*Nmax                                                       !n
          if (dir_j .eq. 2) then !SH                                            !d
             savedauxK(ik,1:3) = SHdiffByStrata(B(:,ik), &
                                 p_X%center%z, p_X%layer,cOME,k_vec(ik),ik)
          else !PSV                                                             !
             savedauxk(ik,1:5) = PSVdiffByStrata(B(:,ik), &                     !
                                 p_X%center%z, p_X%layer,cOME,k_vec(ik),ik)     !
          end if                                                                !
        end do ! ik                                                             !
      end if! onda cilíndrica circular ··········································
      ! fase horizontal (fuente..receptor)
      ! para cada fuente a la profundidad iz ...................
        if (i_zF .eq. 0) then
           nXis = 1 ! la fuente real es sólo una
        else
           nXis = pota(i_zF,2) ! fuente virtual (pueden ser varias
        end if                 !                 a la misa profundidad)
#ifdef ver
       call ETIME(tarray, result)
       print*,"  obtenidos el mecanicos (todos los k)",result-lastresult
       lastresult = result
       print*,"   ahora la fase horizontal:"
       print*,"   pota(i_zF,:)=",pota(i_zF,:)
#endif
        do iXi = 1,nXis
          if (i_zF .ne. 0) then ! fuentes virtuales con igual Z
            itabla_x = 2 + pota(i_zF,1) + iXi
            call asociar(pXi, 0, i_zF, itabla_x)
            xf => pXi%center%x
          end if! i_zF ne 0
#ifdef ver
      print*,"     faseHorzDeFuente iXi",ixi,"de",nXis
      print*,"     fuente virtual pXi%center%x  %z", &
                   xf,pXi%center%z
      print*,"     cada horz de receptores"
#endif
      ! para cada X receptor a la profundidad itabla_z ......................
        n_Xs = pota(itabla_z,1) + pota(itabla_z,2)
        do itabla_x =1,n_Xs
          call asociar(p_x, 0, itabla_z, 2+ itabla_x)

!         if (i_zF .ne. 0) then
!         if (p_X%pointIndex .eq. 3 .or. p_X%pointIndex .eq. 106) then
!          print*,"hit",p_X%pointIndex
!         end if;end if!

!          if (p_x%isboundary .eqv. .false.) then
!              if (p_x%region .ne. 1) cycle; end if!'estr'
#ifdef ver
       print*,"        faseHorzDeReceptor itabla_x",itabla_x,"de",n_Xs
       print*,"        p_x%center%x",p_x%center%x
#endif
      ! fase horizontal
            ! reponer auxK original sin fase horizontal
          auxK(1:pos+1,mecS:mecE)      = savedAuxK(1:pos+1,mecS:mecE)
          auxK(ne:2*nmax,mecS:mecE) = savedAuxK(ne:2*nmax,mecS:mecE)
            ! agregar información fase horizontal de fuente y receptor
          do imec = mecS,mecE !.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
!          if (p_x%guardarMovieSiblings .eqv. .false.) then
            if (isPW) then
!               if (p_x%guardarMovieSiblings .eqv. .false.) then
                auxk(1,imec) = auxk(1,imec) * &           ! onda plana      ·
                exp(-UI*k*(p_x%center%x - xf))            !                 ·
                CYCLE ! imec                              !                 ·
!               end if
            end if ! ························································
!           print*,k_vec(:)
            do ik = 1,pos+1                                                 !
                auxk(ik,imec) = auxk(ik,imec) * &                           !
                exp(cmplx(0.0_8, (-1.0_8)*k_vec(ik)*(p_x%center%x - xf), 8))!c
            end do !  ik                                                    !i
            do ik = ne,2*Nmax                                               !l
                auxk(ik,imec) = auxk(ik,imec) * &                           !i
                exp(cmplx(0.0_8, (-1.0_8)*k_vec(ik)*(p_x%center%x - xf), 8))!n
            end do !  ik                                                    !d
#ifdef ver
       call ETIME(tarray, result)                                           !i
       print*,"        agregado fase horizontal",result-lastresult          !c
       lastresult = result                                                  !a
#endif
                ! guardar el FK por fuente real ---------------------       !
!                if(plotFKS) then
                 if (i_zF .eq. 0) then
                 if (p_x%guardarFK .eqv. .true.) then                       !
                 if (dir_j .eq. 2) then
                    if ( imec .eq. 1) then
                      p_x%FK(J,1:nmax,3) = &                                !
                      p_x%FK(J,1:nmax,3) + auxK(1:nmax,imec)
                    end if
                 else
                    if( imec .le. 2) then
                      p_x%FK(J,1:pos,imec) = auxK(1:pos,imec)
                      p_x%FK(J,pos+1:nmax,imec) = 0
                    end if
                end if;end if;end if!;end if!-------------------------       !
#ifdef ver
       call ETIME(tarray, result)
       print*,"        guardado el integrando",result-lastresult            !
       lastresult = result
#endif
      ! K -> X  .........................................................   !

#ifdef ver
        write(arg,'(a,I0,I0,a)') "k_at",p_x%pointIndex,imec,".m"
        OPEN(739,FILE=trim(arg),FORM="FORMATTED",ACTION='WRITE')
        write(arg,'(a,I0,a,I0)') "kat",p_x%pointIndex,"_",imec
        CALL scripToMatlabMNmatrixZ(2*nmax,1,auxK(:,imec),arg,739)
        close(739)
#endif
             auxK(1,iMec) = sum(auxK(1:pos+1,iMec))+sum(auxK(ne:2*nmax,iMec))
             auxK(1,iMec) = auxK(1,iMec)*dk

#ifdef ver
       call ETIME(tarray, result)                                           !
       print*,"        fork non movie", result-lastresult                   !
       lastresult = result                                                  !
#endif
!         end if  ! is a movie                                               !
      end do !imec !.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
              !
              !                    |        RECEPTOR
              !                    | allpoint  |   boupoint
              !   ----------------------------------------------
              !   fuente real      |     W,U   | G,T: term. indep.
              !   fuente virtual   |      G    | G,T: ibem.mat.
          ! Se ha calculado el campo difractado por estratificación,
          ! falta el campo directo (analítico)
!              if(p_x%isBoundary) then
!                if(i_zF .eq. 0) then ! pXi es la fuente real
!                  ! term. indep.
!                  call fill_termindep(auxK(1,mecS:mecE),&
!                                      come,mecS,mecE,p_x,pXi,dir_j)
!                else ! pXi es una fuente virtual
!                  ! func. Green tracciones para matriz IBEM
!                  call fill_ibemMat(i_zF,auxK(1,mecS:mecE),&
!                                    come,mecS,mecE,p_x,pXi,dir_j)
!                end if
!              else ! p_x no es un punto de colocación, es un receptor
                call fill_diffbyStrata(i_zF,J,auxK(1:2*nmax,mecS:mecE),&
                                       come,mecS,mecE,p_x,pXi,dir_j)
!              end if ! isBoundary
#ifdef ver
      print*,"La función de Green resultado se distribuye";print*,""
#endif
          end do !itabla_x
        end do !iXi
      end do !itabla_z
#ifdef ver
       stop "2921 end of diffField_at_iz"
#endif
end subroutine diffField_at_iz

subroutine asociar(PX, i_Fuente ,itabla_z, itabla_x)
      use resultVars, only : pota,Punto,allpoints
      use sourcevars, only : Po
      implicit none
      type(Punto), pointer :: PX
      integer, intent(in) :: itabla_x, itabla_z, i_Fuente
        nullify(PX)
        if (itabla_z .ne. 0) then
      !    if (pota(itabla_z, itabla_x) .gt. 0) then
             PX => allpoints(pota(itabla_z, itabla_x))
      !    else
      !       PX => boupoints(abs(pota(itabla_z, itabla_x)))
      !    end if
        else
          if (i_Fuente .eq. 0) stop "i_fuente = 0 en asociar"
          PX => Po(i_Fuente)
        end if
end subroutine asociar

subroutine makeGANU (J)
      use modelVars, only : vecNK,SpliK,nmax,cOME,k_vec, &
         gamma=>gamma_E,nu=>nu_E,eta=>eta_E,ALFA,BETA,N
!     use dislin
      implicit none
      integer :: J,po,ne,e,ii,ik,ikI(2),ikF(2)
      real*8  :: k
      complex*16 :: omeAlf,omeBet

      ! gamma y nu en esta frecuencia
       po = min(int(vecNK(J)*SpliK),nmax); ne = 2*nmax-(po-2)
       ikI(1) = 1  ! positivos
       ikF(1) = po+1 !
       ikI(2) = ne     ! negativos
       ikF(2) = 2*nmax !
       do e = 1,N+1
          omeAlf = cOME**2.0/ALFA(e)**2.0
          omeBet = cOME**2.0/BETA(e)**2.0

       ! de ondas cilíndricas
       do ii = 1,2 ! los num de onda positivos, negativos
       do ik = ikI(ii),ikF(ii) !  todos los num de onda
          k = k_vec(ik)
          ! algunas valores constantes para todo el estrato
          gamma(ik,e) = sqrt(omeAlf - k**2.0)
          nu(ik,e) = sqrt(omeBet - k**2.0)
          ! Se debe cumplir que la parte imaginaria del número de onda
          ! vertical debe ser menor que cero. La parte imaginaria contri-
          ! buye a tener ondas planas inhomogéneas con decaimiento expo-
          ! nencial a medida que z crece.
          if(aimag(gamma(ik,e)).gt.0.0_8)gamma(ik,e) = conjg(gamma(ik,e))
          if(aimag(nu(ik,e)).gt.0.0_8)nu(ik,e)=conjg(nu(ik,e))
          eta(ik,e) = 2.0*gamma(ik,e)**2.0 - cOME**2.0 / BETA(e)**2.0
       end do ! ik
       end do ! ii
!      call qplot(real(k_vec,4),real(aimag(gamma(:,e)),4), 2*nmax)
!      stop
       end do ! e
end subroutine makeGANU

subroutine makeGANU0
      use modelVars, only : come,ome, &
         gamma=>gamma_E,nu=>nu_E,eta=>eta_E, &
        N,alfa0,beta0,alfa,beta
      use sourceVars, only : PoFte=>Po,iFte=>currentiFte
      use glovars, only: Ur,PWfrecReal
      implicit none
      integer :: e,ik
      real*8  :: k
       k = 0
       do e = 1,N+1
       ! de la onda plana
       ik = 0
       if (PoFte(iFte)%PW_pol .eq. 1) then
        if (PWfrecReal) then
        k = real(OME/beta0(e)*sin(PoFte(iFte)%gamma))
        else
        k = real(COME/beta(e)*sin(PoFte(iFte)%gamma))
        end if
      elseif (PoFte(iFte)%PW_pol .eq. 2) then
        if (PWfrecReal) then
        k = real(OME/alfa0(e)*sin(PoFte(iFte)%gamma))
        else
        k = real(cOME/alfa(e)*sin(PoFte(iFte)%gamma))
        end if
      elseif (PoFte(iFte)%PW_pol .eq. 3) then
        if (PWfrecReal) then
        k = real(OME/beta0(e)*sin(PoFte(iFte)%gamma))
        else
        k = real(COME/beta(e)*sin(PoFte(iFte)%gamma))
        end if
      end if!
       if (PWfrecReal) then
       gamma(ik,e) = sqrt(UR*OME**2.0/ALFA0(e)**2.0 - k**2.0)
!      print*,gamma(ik,e),"gamma(ik,e)"
       nu(ik,e) = sqrt(UR*OME**2.0/BETA0(e)**2.0 - k**2.0)
       if(aimag(gamma(ik,e)).gt.0.0_8)gamma(ik,e) = conjg(gamma(ik,e))
!      print*,gamma(ik,e),"gamma(ik,e)"
       if(aimag(nu(ik,e)).gt.0.0_8)nu(ik,e)=conjg(nu(ik,e))
       eta(ik,e) = 2.0*gamma(ik,e)**2.0 - UR*OME**2.0 / BETA0(e)**2.0
       else
       gamma(ik,e) = sqrt(cOME**2.0/ALFA(e)**2.0 - k**2.0)
       nu(ik,e) = sqrt(cOME**2.0/BETA(e)**2.0 - k**2.0)
       if(aimag(gamma(ik,e)).gt.0.0_8)gamma(ik,e) = conjg(gamma(ik,e))
       if(aimag(nu(ik,e)).gt.0.0_8)nu(ik,e)=conjg(nu(ik,e))
       eta(ik,e) = 2.0*gamma(ik,e)**2.0 - cOME**2.0 / BETA(e)**2.0
       end if
       end do
end subroutine makeGANU0

subroutine gloMat_PSV(this_A,k,ik)
      ! Calcular para +k. Una vez invertida la matrix, hay paridades para -k.
      use modelVars, only : N,Z,AMU!,BETA,ALFA,LAMBDA
      use gloVars, only : UI,UR,Z0
!      use debugStuff
      use modelVars, only : gamma_E,nu_E,eta_E,nmax
      use refSolMatrixVars, only : subMatD0,subMatS0
      implicit none
      complex*16,    intent(inout), dimension(:,:),pointer :: this_A
      real*8,     intent(in),pointer     :: k
      integer :: ik,neik
      real*8     :: k2
      complex*16, dimension(2,4) :: subMatD,subMatS
      complex*16, dimension(4) :: diagMat
      complex*16 :: gamma,nu,xi,eta,ega,enu,ck
      integer    :: i,iR,iC,e,bord
      this_A = 0
      ck = -k*UR
      iR= 0;iC= 0;i=1
      if (Z(0) .lt. 0.0) then !Half-Space por arriba de z=0
        i = 0;iR= -2; iC= -2
      end if

      DO e = i,N+1
!         if (ik .ne. 0) then
          gamma = gamma_E(ik,e)
          nu = nu_E(ik,e)
!         else
!         gamma = sqrt(cOME_i**2.0_8/ALFA(e)**2.0_8 - k**2.0_8)
!         nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
!         if(aimag(gamma).gt.0.0)gamma = conjg(gamma)
!         if(aimag(nu).gt.0.0)nu= conjg(nu)
!         end if
          xi = k**2.0 - nu**2.0
!         eta = 2.0*gamma**2.0 - cOME_i**2.0 / BETA(e)**2.0
          eta = eta_E(ik,e)
!         print*,ik,e,k,gamma,nu,eta," L2184"
          ! en fortran los elementos se indican por columnas:
          subMatD0(1:2,1:4,e,ik) = RESHAPE((/ -gamma, ck, ck,nu,&
                                               gamma, ck, ck,-nu /), &
                           (/ 2,4 /))
          subMatD0(1:2,1:4,e,ik) = subMatD0(1:2,1:4,e,ik) * UI

          k2 = 2.0*k
!         subMatS0 = RESHAPE((/ xi,-k2*gamma,-k2*nu,-xi,&
!                              xi,k2*gamma,k2*nu,-xi /),&
!                          (/2,4/))
          subMatS0(1:3,1:4,e,ik) = RESHAPE(&
                        (/ xi,      -k2*gamma,     eta,     &
                          -k2*nu,     -xi,        k2*nu,   &
                           xi,       k2*gamma,     eta,       &
                           k2*nu,     -xi,       -k2*nu /),&
                           (/3,4/))
          subMatS0(1:3,1:4,e,ik) = amu(e) * subMatS0(1:3,1:4,e,ik)
          ! y para k negativo aprovechando la simetria de gamma y nu
          neik = 2*nmax - (ik-2)
          if (neik .le. 2*nmax) then
          subMatD0(1:2,1:4,e,neik) = RESHAPE((/ -gamma, -ck, -ck,nu,&
                                               gamma, -ck, -ck,-nu /), &
                           (/ 2,4 /))
          subMatD0(1:2,1:4,e,neik) = subMatD0(1:2,1:4,e,neik) * UI
          ! y para el k negativo aprovechando simetria de gamma,nu,xi,eta,
          subMatS0(1:3,1:4,e,neik) = RESHAPE(&
                        (/ xi,      k2*gamma,     eta,     &
                           k2*nu,     -xi,       -k2*nu,   &
                           xi,      -k2*gamma,     eta,       &
                          -k2*nu,     -xi,        k2*nu /),&
                           (/3,4/))
          subMatS0(1:3,1:4,e,neik) = amu(e) * subMatS0(1:3,1:4,e,neik)
          end if



          ! la profundidad z de la frontera superior del estrato
!         z_i = Z(e)   ! e=1  ->  z = z0 = 0
!         z_f = Z(e+1) ! e=1  ->  z = z1
        do bord = 0,1
          if ((e .eq. 0) .and. (bord .eq. 0)) cycle
          if (e+bord > N+1) exit
          ! si 1+0;1+1;2+0;[2+1] > 2
        if (bord .eq. 0) then !---->---->---->---->---->---->
          if (e /= N+1) then !(radiation condition lower HS)
            ega = exp(-UI * gamma * (Z(e+1)-Z(e)))
            enu = exp(-UI * nu * (Z(e+1)-Z(e)))
          else
            ega = Z0
            enu = Z0
          end if
!         diagMat = RESHAPE((/ UR, Z0, Z0, Z0, &
!                              Z0, UR, Z0, Z0, &
!                              Z0, Z0, ega, Z0, &
!                              Z0, Z0, Z0, enu /), &
!                          (/ 4,4 /))
          diagMat = (/ UR, UR, ega, enu /)
        else !bord .eq. 1 -----------------------------------
          if (e /= 0) then !(radiation condition upper HS)
            ega = exp(-UI * gamma * (Z(e+1)-Z(e)))
            enu = exp(-UI * nu * (Z(e+1)-Z(e)))
          else
            ega = Z0
            enu = Z0
          end if !<----<----<----<----<----<----<----<----<--
!         diagMat = RESHAPE((/ ega, Z0, Z0, Z0, &
!                              Z0, enu, Z0, Z0, &
!                              Z0, Z0, UR, Z0, &
!                              Z0, Z0, Z0, UR /), &
!                          (/ 4,4 /))
          diagMat = (/ ega, enu, UR, UR /)
        end if
          ! desplazamientos
!         subMatD = matmul(subMatD0(1:2,1:4,e,ik),diagMat)
          subMatD(1:2,1) = subMatD0(1:2,1,e,ik)*diagMat(1)
          subMatD(1:2,2) = subMatD0(1:2,2,e,ik)*diagMat(2)
          subMatD(1:2,3) = subMatD0(1:2,3,e,ik)*diagMat(3)
          subMatD(1:2,4) = subMatD0(1:2,4,e,ik)*diagMat(4)
          ! esfuerzos
!         subMatS = matmul(subMatS0(1:2,1:4,e,ik),diagMat)
          subMatS(1:2,1) = subMatS0(1:2,1,e,ik)*diagMat(1)
          subMatS(1:2,2) = subMatS0(1:2,2,e,ik)*diagMat(2)
          subMatS(1:2,3) = subMatS0(1:2,3,e,ik)*diagMat(3)
          subMatS(1:2,4) = subMatS0(1:2,4,e,ik)*diagMat(4)
        !ensamble de la macro columna i
          !evaluadas en el borde SUPERIOR del layer i
          if (bord == 0 .AND. e /= 0 ) then
           if (e /= N+1) then !(radiation condition)
            if (e/=i) then !(only stress bound.cond. in the surface
             this_A( iR-1 : iR   , iC+1 : iC+4 ) = -subMatD
            end if
             this_A( iR+1 : iR+2 , iC+1 : iC+4 ) = -subMatS
           else
           if (e/=i) then !(only stress bound.cond. in the surface
             this_A( iR-1 : iR   , iC+1 : iC+2 ) = -subMatD(:,1:2)
           end if
             this_A( iR+1 : iR+2 , iC+1 : iC+2 ) = -subMatS(:,1:2)
            ! exit
           end if
          end if

          !evaluadas en el borde INFERIOR del layer i
          if (bord == 1 .AND. e /= N+1 ) then ! cond de radiación en HS
           if (e /= 0) then !(radiation condition upward)
            ! ambas ondas
            this_A( iR+3 : iR+4 , iC+1 : iC+4 ) = subMatD
            this_A( iR+5 : iR+6 , iC+1 : iC+4 ) = subMatS
           else
            ! solo onda hacia arriba
            this_A( iR+3 : iR+4 , iC+3 : iC+4 ) = subMatD(:,3:4)
            this_A( iR+5 : iR+6 , iC+3 : iC+4 ) = subMatS(:,3:4)
           end if
          end if
        end do !bord loop del borde i superior o nferior
          iR= iR+4
          iC= iC+4
      END DO !{e} loop de las macro columnas para cada estrato
!            call showMNmatrixZ(size(this_A,1),&
!            size(this_A,2),this_A,"A    ",6)
!            stop "gloMat_PSV"
end subroutine gloMat_PSV

subroutine globalmatrix_SH(this_A,k,ik)
      use modelVars, only : N,Z,AMU!,LAMBDA!,BETA,ALFA
      use gloVars, only : UI,UR,Z0
!      use debugStuff
      use modelVars, only : nu_E
      implicit none

      complex*16,    intent(inout), dimension(:,:),pointer :: this_A
      real*8,     intent(in),pointer     :: k
      integer :: ik
!     complex*16, intent(in),pointer     :: cOME_i

      real*8     :: z_i
      complex*16, dimension(1,2) :: subMatD, subMatS, subMatD0, subMatS0
      complex*16, dimension(2,2) :: diagMat
      complex*16 :: nu,enuN,enuP
      integer    :: i,iR,iC,e,bord
      this_A = 0
      z_i = k ! (nada mas para que no chiste)
      iR= 0;iC= 0;i=1
      if (Z(0) .lt. 0.0) then !Half-Space por arriba de z=0
        i = 0;iR= -1; iC= -1
      end if

      DO e = i,N+1!cada estrato
          nu = nu_E(ik,e)!; print*,e,N,nu
          ! algunas valores constantes para todo el estrato
!         nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
          ! Se debe cumplir que la parte imaginaria del número de onda
          ! vertical debe ser menor que cero. La parte imaginaria contri-
!         ! buye a tener ondas planas inhomogéneas con decaimiento expo-
          ! nencial a medida que z es mayor que cero.
!         if(aimag(nu).gt.0.0)nu= conjg(nu)

          ! en fortran los elementos se indican por columnas:
          subMatD0 = RESHAPE((/ UR,UR /), (/ 1,2 /))
          subMatS0 = RESHAPE((/ -UI*nu,UI*nu /),(/1,2/))
          subMatS0 = amu(e) * subMatS0

          ! la profundidad z de la frontera superior del estrato
!         z_i = Z(e)   ! e=1  ->  z = z0 = 0
!         z_f = Z(e+1) ! e=1  ->  z = z1
        do bord = 0,1
          if ((e .eq. 0) .and. (bord .eq. 0)) cycle
          if (e+bord > N+1) then ! si 1+0;1+1;2+0;[2+1] > 2
            exit
          end if
          ! la profundidad z de la frontera superior del estrato
                z_i = Z(e+bord)   ! e=1 , bord=0  ->  z = z0 = 0
                                  ! e=1 , bord=1  ->  z = Z1 = h1
          !downward waves
          if (e /= 0) then !(radiation condition upper HS)
            enuN = exp(-UI * nu * (z_i-Z(e)))
          else
            enuN = Z0
          end if
          !upward waves
          if (e /= N+1) then !(radiation condition lower HS)
            enuP = exp(UI * nu * (z_i-Z(e+1)))
          else
            enuP = Z0
          end if
          diagMat = RESHAPE((/ enuN, Z0,   &
                               Z0,   enuP /), &
                              (/ 2,2 /))
          subMatD = matmul(subMatD0,diagMat)
          subMatS = matmul(subMatS0,diagMat)
!         call showMNmatrixZ(2,2,diagmat,"diag ",6)
!         call showMNmatrixZ(1,2,submatd,"  D  ",6)
!         call showMNmatrixZ(1,2,submats,"  S  ",6)
        !ensamble de la macro columna i
          !evaluadas en el borde SUPERIOR del layer i
          if (bord == 0 .AND. e /= 0 ) then
           if (e /= N+1) then !(radiation condition downward)
            if (e/=i) then !(only stress bound.cond. in the surface
             this_A( iR   : iR   , iC+1 : iC+2 ) = -subMatD!; print*,"b1"
            end if !(e==2,3,...)
             this_A( iR+1 : iR+1 , iC+1 : iC+2 ) = -subMatS!; print*,"b2"
           else
           if (e/=i) then !(only stress bound.cond. in the surface
             this_A( iR   : iR   , iC+1 : iC+1 ) = -subMatD(:,1:1)!; print*,"b3"
           end if
             this_A( iR+1 : iR+1 , iC+1 : iC+1 ) = -subMatS(:,1:1)!; print*,"b4"
!            exit
           end if
          end if
          !evaluadas en el borde INFERIOR del layer i
          if (bord == 1 .AND. e /= N+1 ) then ! cond de radiación en HS
           if (e /= 0) then !(radiation condition upward)
            ! ambas ondas
            this_A( iR+2 : iR+2 , iC+1 : iC+2 ) = subMatD!; print*,"b5",iR+2, iC+1
            this_A( iR+3 : iR+3 , iC+1 : iC+2 ) = subMatS
           else
            ! solo onda hacia arriba
            this_A( iR+2 : iR+2 , iC+2 : iC+2 ) = subMatD(:,2:2)!; print*,"b5",iR+2, iC+1
            this_A( iR+3 : iR+3 , iC+2 : iC+2 ) = subMatS(:,2:2)
           end if
!           this_A( iR+2 : iR+2 , iC+1 : iC+2 ) = subMatD!; print*,"b5",iR+2, iC+1
!           this_A( iR+3 : iR+3 , iC+1 : iC+2 ) = subMatS
          end if
!         print*,"bord=",bord,"e=",e,"IR=",IR,"IC=",IC
!         call showMNmatrixZ(2*N+2,2*N+2,this_A,"Ash  ",6)
        end do !bord loop del borde i superior o nferior
          iR= iR+2
          iC= iC+2
      END DO !{e} loop de las macro columnas para cada estrato
!         call showMNmatrixZ(size(this_A,1),size(this_A,2),this_A,"Ash  ",6)
!         stop 5248
!     print*,"glomat",ik,"line2789"
end subroutine globalmatrix_SH

subroutine inverseA(A,ipiv,work,n)
!     use debugstuff
      integer, intent(in) :: n
      complex*16, dimension(:,:), intent(inout),pointer :: A
      integer, dimension(:), intent(inout),pointer :: ipiv
      complex*16, dimension(:),intent(inout),pointer :: work
      integer :: info
      integer :: lwork
      lwork = n*n
      call zgetrf(n,n,A,n,ipiv,info)
      if(info .ne. 0) stop "inverseA :Problem at LU factorization of matrix Try more Q"
      ! Para más de 3 estratos !--------
      ! 1) convertir a almacenamiento bandeado
!       2) call ZGBTRF( n, n, 5, 5, A, n, IPIV, INFO )
!     if(info .ne. 0) stop "Problem at LU factorization of matrix "
      ! luego con cada B
!       3) call ZGBTRS( 'N', N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
      call zgetri(n,a,n,ipiv,work,lwork,info)
      if(info .ne. 0) stop "inverseA :Problem at inverse of matrix "
end subroutine inverseA

subroutine intrplr_gloMat(k0,n,pt_cOME_i,pt_ipivA,pt_workA)
      use refSolMatrixVars, only : Ak !Ak(#,#,ik:2*nmax)
      use modelVars, only : k_vec, NMAX
      use fitting
      use modelVars, only: alfa,beta,Nestr => N
      use modelVars, only : gamma=>gamma_E,nu=>nu_E
      use sourceVars, only  : PSV
      implicit none
!      interface   ! no se requiere la interfaz si estan en el mismo modulo
!      include 'interfaz.f'
!      end interface
      integer :: i,k0,k1,k2,n,tam,r,c,ik,e
!     integer :: d != 4
      real*8 :: ratio,step
      real*8, dimension(n)     :: xdat
      integer,dimension(n)     :: idat
      complex*16, dimension(n) :: ydat
      complex*16, pointer :: pt_cOME_i
      integer, dimension(:), pointer :: pt_ipivA
      complex*16, dimension(:),pointer :: pt_workA
      complex*16, dimension(:,:), pointer :: pointAp
      real*8, pointer :: pt_k
      complex*16 :: m
      ! armar los vectores de datos
      ratio = 1.0 ! mayor o igual a 1
      tam = size(Ak,1)
!     step = ((nmax+1) - (k0))/(n-1) !iguales
      if (int(n/2) .lt. real(n)/2.) then
      step = (int((n-1)/2)) * n  !impar
      else
      step = (int((n-1)/2) + 1) * n !par
      end if
      step = ((nmax+1) - (k0))/(ratio*step) !crecim geometrico
      idat(1) = k0
      do i=2,n-1
         idat(i) =  idat(i-1) + int((i-1)*step)
         xdat(i) = k_vec(idat(i))
      end do
      idat(n) = nmax+1
      xdat(n) = k_vec(idat(n))
!     do i=1,n;print*,i,idat(i);end do;stop
      do i=2,n
         ik = idat(i)
         pointAp => Ak(1:tam,1:tam,ik)
         pt_k => k_vec(ik)
            do e=1,Nestr+1
              gamma(ik,e) = sqrt(pt_cOME_i**2.0/ALFA(e)**2.0 - pt_k**2.0)
              nu(ik,e) = sqrt(pt_cOME_i**2.0/BETA(e)**2.0 - pt_k**2.0)
              if(aimag(gamma(ik,e)).gt.0.0_8)gamma(ik,e) = conjg(gamma(ik,e))
              if(aimag(nu(ik,e)).gt.0.0_8)nu(ik,e)=conjg(nu(ik,e))
            end do ! e
         if (PSV) then
         call gloMat_PSV(pointAp,pt_k,ik)
         else
         call globalmatrix_SH(pointAp,pt_k,ik)
         end if
         call inverseA(pointAp,pt_ipivA,pt_workA,tam)
      end do

      ! interpolar cada elemento (rectas)
      do r=1,tam
      do c=1,tam
        ydat(1) = Ak(r,c,idat(1))
!       print*,ydat(1)
        k1 = idat(1)+1
        do i=2,n
           k2 = idat(i) !k0 + int((i-1)*step)
           ydat(i) = Ak(r,c,k2)
           m = (ydat(i)-ydat(i-1))/(xdat(i)-xdat(i-1))
           do ik=k1,k2
             Ak(r,c,ik) = ydat(i-1) + m * (k_vec(ik) - xdat(i-1))
           end do
           k1 = k2+1
        end do
!       print*,""
!       do i=k0-2,nmax+1
!         print*,i,k_vec(i),Ak(r,c,i)
!       end do
!       stop
      end do !c
      end do !r
end subroutine intrplr_gloMat

!subroutine intrplr_zpoly_gloMat

      ! de los elementos a_ij
      ! pares en k:  i={1,3,5...} con j={1,3,4...} ; i={2,4,6...} con j={2,4,6...}
      ! impares  k:  i={1,3,5...} con j={2,4,6...} ; i={2,4,6...} con j={1,3,5...}
subroutine parImpar_gloMat
      use refSolMatrixVars, only : Ak !Ak(#,#,ik:2*nmax)
      use modelVars, only : NMAX
!      use debugStuff
      implicit none
      integer :: r,c,ik,ikp,tam
      integer,dimension(:,:),allocatable :: Comal
      logical :: s !empieza en priemr columna
      tam = size(Ak,1)
      allocate(Comal(tam,tam))
      Comal =1
      ! los odd
      s = .false.
      do r = 1,tam
      if (s) then
      Comal(r,1:tam:2) = -1
      else
      Comal(r,2:tam:2) = -1
      end if
      s = (.not. s)
      end do
!     call showMNmatrixI(tam,tam, Comal," Que ",6)
      ! negativos
      do ik=nmax+2,2*NMAX
        ikp = (2*NMAX +2-ik)
!     print*,ik,ikp!;cycle
!     call scripToMatlabMNmatrixZ(tam,tam,Ak(1:tam,1:tam,ikp),"lado positivo",6)
        do r=1,tam
        do c=1,tam
        Ak(r,c,ik) = Ak(r,c,ikp) * Comal(r,c)
        end do
        end do
!     call scripToMatlabMNmatrixZ(tam,tam,Ak(1:tam,1:tam,ik),"lado negativo",6)
!     stop
      end do
!     stop "parImpar_gloMat"
end subroutine parImpar_gloMat

! G_stra - termIndPSV
subroutine PSVvectorB_ondaplana(this_B,gamma)
      use modelVars, only : n,lambda0,amu0,lambda,amu,alfa0,beta0,Z,beta,alfa
      use glovars, only:UI,z0,PWfrecReal
      use sourceVars, only: Po,iFte=>currentiFte
      use modelVars, only : cOME,ome
      !use debugStuff
      implicit none
      complex*16, intent(inout), dimension(1:4*N+2) :: this_B
!     complex*16, intent(in)    :: come ! no trae amortiguamiento
      real*8, intent(in) :: gamma
      integer :: i,e
      real*8,dimension(1:2) :: theta
      complex*16 :: kx,kz,U,W,c,la,am
      real*8 :: z_loc
      !     Colocamos la onda incindente en la interfaz
      !     con el semiespacio de abajo.
      theta = 0; c = z0
      e = Po(iFte)%layer! N+1
      z_loc = Z(N+1)

      if (Po(iFte)%PW_pol .eq. 1) then
        if (PWfrecReal) then
        c = beta0(e) !SV
        else
        c = beta(e) !SV
        end if
        theta(1) = cos(gamma)
        theta(2) = sin(gamma)
      elseif (Po(iFte)%PW_pol .eq. 2) then
        if (PWfrecReal) then
        c = alfa0(e) !SV
        else
        c = alfa(e) !SV
        end if
        theta(1) = sin(gamma)
        theta(2) = -cos(gamma)
      end if
      !
      if (PWfrecReal) then
      kx = ome/c * sin(gamma)
      kz = sqrt((ome/c)**2 - kx**2)
      if(aimag(kz).gt.0.0_8) kz=conjg(kz)
      la = LAMBDA0(e)
      am = AMU0(e)
      else
      kx = come/c * sin(gamma)
      kz = sqrt((come/c)**2 - kx**2)
      if(aimag(kz).gt.0.0_8) kz=conjg(kz)
      la = LAMBDA(e)
      am = AMU(e)
      end if
      U = (theta(1))* exp(UI * kz * (z_loc))
      W = (theta(2))* exp(UI * kz * (z_loc))

      i=0
      this_B(1:4*N+2) = Z0
      if (Z(0) .lt. 0.0) then ! Semiespacio en z<0 ···········
        i = 2                                                !
        this_B(1+4*(e-1)-2 + i) = W !  w                     !
        this_B(1+4*(e-1)-1 + i) = U !  u                     !
      end if                                                 !
      ! ······················································
      if (e .ne. 1) then                                     !
        this_B(1+4*(e-1)-2 + i) = W !  w                     !
        this_B(1+4*(e-1)-1 + i) = U !  u                     !
      end if                                                 !
      !.......................................................

      ! Tracciones en la frontera de la región de la fuente.........
      this_B(1+4*(e-1)   + i) = UI * ( &                           !
                            ( W * kz * (la + 2.0 * am)) & !
                          - ( U * kx * la)) ! szz           !
      this_B(1+4*(e-1)+1 + i) = UI * am &                      !
                          * ( kz * U - kx * W ) ! szx              !
      !                   sxx = UI * ( &                           !
      !                   - ( U * kx * (LAMBDA(e) + 2.0*AMU(e))) & !
      !                   + ( W * kz * LAMBDA(e)))                 !
      !............................................................!
!     call showMNmatrixZ(4*N+2,1, this_B ,"  B  ",6)
end subroutine PSVvectorB_ondaplana


subroutine PSVvectorB_force(i_zF,this_B,tam,pXi,direction,cOME,k,ik)
      use modelVars, only: N,Z,AMU,LAMBDA,RHO!,NPAR,BETA,ALFA
      use gloVars, only : UR,UI,PI,Z0
!     use debugStuff
      use resultvars, only : Punto
      use sourceVars, only: Po,iFte=>currentiFte!use sourceVars, only: tipofuente
      use modelVars, only : gamma_E,nu_E
      implicit none

      integer, intent(in) :: i_zF,tam,ik
      complex*16, intent(inout), dimension(tam) :: this_B
      integer,    intent(in)    :: direction
      real*8,     intent(in)    :: k
      complex*16, intent(in)    :: cOME
      type(Punto),intent(in),target    :: pXi

      integer, pointer :: e_f
      real*8, pointer :: z_f
      logical, pointer :: fisInterf

      integer :: iIf,nInterf
      real    :: SGNz
      complex*16 :: gamma,nu,DEN,L2M!, argum!,sincmod
      real*8     :: errT = 0.0001_8
!     complex*16 :: omeAlf,omeBet

      ! una para cada interfaz (la de arriba [1] y la de abajo [2])
      real*8,     dimension(2) :: z_loc
      complex*16, dimension(2) :: egamz,enuz,gamz,nuz
      complex*16, dimension(2) :: sincGamma, sincNu

                                  !  1 para Greeni1 (horizontal),
                                  !  3 para Greeni3 (vertical)
      complex*16, dimension(2) :: G11,G31,G33
      complex*16, dimension(2) :: s331,s131,s333,s313

      !real*8 :: a
      !real*8, pointer :: cose,seno
      integer :: el_tipo_de_fuente

      this_B = Z0
      sincGamma = Z0
      sincNu = Z0
      e_f => pXi%layer
      z_f => pXi%center%z
      fisInterf => pXi%isOnInterface
      nInterf = 2
      if (fisInterf) then
         nInterf = 1
         go to 349
      end if
      z_loc(1:2) = 0.0_8
      if (e_f .ne. 0) z_loc(1) = Z(e_f) - z_f !downward (-)
      if (e_f .ne. N+1)  z_loc(2) = Z(e_f+1) - z_f !upward (+)

      el_tipo_de_fuente = 2 ! fuente segmento (para ibem)
      if  ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "PSVvectorB_force iFte=0"
      if (i_zF .eq. 0) el_tipo_de_fuente = Po(iFte)%tipofuente

      DEN = 4.0*PI*RHO(e_f)*cOME**2.0
      L2M = LAMBDA(e_f) + 2.0*AMU(e_f)
      gamma = gamma_E(ik,e_f)
      nu = nu_E(ik,e_f)

      do iIf = 1,2
          egamz(iIf) = exp(-UI*gamma*ABS(z_loc(iIf)))
          enuz(iIf) = exp(-UI*nu*ABS(z_loc(iIf)))

          gamz(iIf) = gamma*ABS(z_loc(iIf))
          nuz(iIf) = nu*ABS(z_loc(iIf))
      end do
      if (el_tipo_de_fuente .le. 1) then ! fuente puntual 0 u onda plana 1
        sincGamma(1) = UR * egamz(1)
        sincGamma(2) = UR * egamz(2)
        sincNu(1) = UR * enuz(1)
        sincNu(2) = UR * enuz(2)
!      elseif (el_tipo_de_fuente .eq. 2) then ! la fuente es un segmento
!      ! ahora estamos involucrando información del receptor en el vector de fuente
!      ! en particular, si el receptor (la interfaz) está arriba o abajo de la fuente.
!
!        A = pxi%length * 0.5_8  ! 2a=lenght
!        cose => pxi%cosT
!        seno => pxi%sinT
!
!      ! en cada interfaz
!        ! caso (receptor abajo de la fuente) : interfaz de abajo (2)
!        argum = gamma*A*seno + UR * (k*A*cose)
!        sincGamma(2) = sincmod(argum,gamz(2))*(2*A)
!
!        argum = nu*A*seno + UR * (k*A*cose)
!        sincNu(2) = sincmod(argum,nuz(2))*(2*A)
!
!        ! caso (receptor arriba de la fuente) : interfaz de arriba (1)
!        argum = -gamma*A*seno + UR * (k*A*cose)
!        sincGamma(1) = sincmod(argum,gamz(1))*(2*A)
!
!        argum = -nu*A*seno + UR * (k*A*cose)
!        sincNu(1) = sincmod(argum,nuz(1))*(2*A)

      end if !fuente tipo 0 o 1

      do iIf = 1,2
          if (abs(z_loc(iIf)) > errT ) then
            SGNz = real(z_loc(iIf) / ABS(z_loc(iIf)),4)
          else
            SGNz = 0.0
          end if

      G31(iIf) = -UI/DEN * SGNz*k*(sincGamma(iIf) &
                           - sincNu(iIf))

      if (direction .eq. 1) then !  G31, G11, S331, S131

      G11(iIf) = -UI/DEN * (k**2.0/gamma*sincGamma(iIf) &
                           +nu*sincNu(iIf))

      s331(iIf) = -UR/DEN * ( &
                    (k*gamma*L2M + lambda(e_f)*k**3.0/gamma)* sincGamma(iIf)&
                  + (-2.0*amu(e_f)*k*nu)* sincNu(iIf) &
                    ) !

      s131(iIf) = -UR/DEN * amu(e_f)*SGNz * ( &
                    (2.0*k**2.0)* sincGamma(iIf) &
                    + (nu**2.0-k**2.0)* sincNu(iIf) &
                    ) !
      end if
      !
      if (direction .eq. 3) then ! G33,G31,S333,S313

      G33(iIf) = -UI/DEN * (gamma*sincGamma(iIf) &
                           +k**2.0/nu*sincNu(iIf))

      s333(iIf) = -UR/DEN * SGNz * ( &
                    (gamma**2.0*L2M + k**2.0*lambda(e_f))* sincGamma(iIf) &
                  + (2.0*amu(e_f)*k**2.0)* sincNu(iIf) &
                    ) !

      s313(iIf) = -UR/DEN * amu(e_f) * ( &
                    (2.0*k*gamma)* sincGamma(iIf) &
                  - (k/nu*(nu**2.0-k**2.0))* sincNu(iIf) &
                    ) !
      end if
      end do !iIf interface

 349  if (fisInterf) then
      if(direction .eq. 1) then
        s331(1) = 0
        S131(1) = + cmplx(1.0 / (2.0 * PI ),0.0,8)
      elseif (direction .eq. 3) then
        s333(1) = + cmplx(1.0 / (2.0 * PI ),0.0,8)
        S313(1) = 0
      end if
        G33(1) = 0
        G31(1) = 0
        G11(1) = 0
      end if

      ! El vector de términos independientes genera el campo difractado
      if (z(0) .gt. 0.0) then !--------------------------
      if (direction .eq. 1) then ! fuerza HORIZONTAL
      !                     =      (1) interfaz de arriba
       if (e_f .ne. 1) then
        this_B(1+4*(e_f-1)-2) = G31(1)!  w
        this_B(1+4*(e_f-1)-1) = G11(1)!  u
       end if
        this_B(1+4*(e_f-1)  ) = S331(1)! szz
        this_B(1+4*(e_f-1)+1) = S131(1)! szx   ! delta

      if (.not. fisInterf) then ! la fuerza no en la interfaz
      !                     =      (2) interfaz de abajo
       if (e_f .ne. N+1) then
        this_B(1+4*(e_f-1)+2) = - G31(2)!  w
        this_B(1+4*(e_f-1)+3) = - G11(2)!  u
        this_B(1+4*(e_f-1)+4) = - S331(2)! szz
        this_B(1+4*(e_f-1)+5) = - S131(2)! szx
       end if
      end if

      elseif (direction .eq. 3) then ! fuerza VERTICAL
      !                     =     (1) interfaz de arriba
       if (e_f .ne. 1) then
        this_B(1+4*(e_f-1)-2) = G33(1)!  w
        this_B(1+4*(e_f-1)-1) = G31(1)!  u
       end if
        this_B(1+4*(e_f-1)  ) = S333(1)! szz   ! delta
        this_B(1+4*(e_f-1)+1) = S313(1)! szx

      if (.not. fisInterf) then
      !                     =    (2) interfaz de abajo (con Fza en el HS no entra aqui)
       if (e_f .ne. N+1) then
        this_B(1+4*(e_f-1)+2) = - G33(2)!  w
        this_B(1+4*(e_f-1)+3) = - G31(2)!  u
        this_B(1+4*(e_f-1)+4) = - S333(2)! szz
        this_B(1+4*(e_f-1)+5) = - S313(2)! szx
       end if
      end if
      end if ! direction
      else !Z(0) < 0 un semiespacio arriba --------------------------
      if (direction .eq. 1) then ! fuerza HORIZONTAL
      !                     =      (1) interfaz de arriba
       if (e_f .ne. 0) then
        this_B(1+4*(e_f-1)  ) = G31(1)!  w
        this_B(1+4*(e_f-1)+1) = G11(1)!  u
        this_B(1+4*(e_f-1)+2) = S331(1)! szz
        this_B(1+4*(e_f-1)+3) = S131(1)! szx   ! delta
       end if
      !
      if (.not. fisInterf) then ! la fuerza no en la interfaz
      !                     =      (2) interfaz de abajo
       if (e_f .ne. N+1) then
        this_B(1+4*(e_f-1)+4) = - G31(2)!  w
        this_B(1+4*(e_f-1)+5) = - G11(2)!  u
        this_B(1+4*(e_f-1)+6) = - S331(2)! szz
        this_B(1+4*(e_f-1)+7) = - S131(2)! szx
       end if
      end if

      elseif (direction .eq. 3) then ! fuerza VERTICAL
      !                     =     (1) interfaz de arriba
       if (e_f .ne. 0) then
        this_B(1+4*(e_f-1)  ) = G33(1)!  w
        this_B(1+4*(e_f-1)+1) = G31(1)!  u
        this_B(1+4*(e_f-1)+2) = S333(1)! szz   ! delta
        this_B(1+4*(e_f-1)+3) = S313(1)! szx
       end if
      !
      if (.not. fisInterf) then
      !                     =    (2) interfaz de abajo (con Fza en el HS no entra aqui)
       if (e_f .ne. N+1) then
        this_B(1+4*(e_f-1)+4) = - G33(2)!  w
        this_B(1+4*(e_f-1)+5) = - G31(2)!  u
        this_B(1+4*(e_f-1)+6) = - S333(2)! szz
        this_B(1+4*(e_f-1)+7) = - S313(2)! szx
       end if
      end if
      end if ! direction
      end if
      ! cuand la fuente es una onda plana, sólo un número de onda
      ! es distinto de cero, k = omega/c y no se hace fft en k
!     print*,"e_f=",e_f
!     print*,this_B
end subroutine PSVvectorB_force


subroutine PSVpaGaNU(J)
      !use modelVars, only : cOME,vecNK,k_vec,nmax,SpliK,gamma_E,nu_E
      use modelVars
      use gloVars, only : UR,UI,PI
      use refSolMatrixVars, only : BparaGa,BparaNu
      implicit none
      integer, intent(in) :: J
      complex*16, dimension(:,:,:,:),pointer :: this_B
      complex*16 :: gamma,nu,DEN,L2M!,omeAlf,omeBet
      integer :: ii,ik,e,Ga_o_Nu,dir,pos,neg,ikI(2),ikF(2)
      complex*16, dimension(2) :: G11,G31,G33
      complex*16, dimension(2) :: s331,s131,s333,s313
      real*8 :: k
       pos = min(int(vecNK(J)*SpliK),nmax); neg = 2*nmax-(pos-2)
       ikI(1) = 1
       ikF(1) = pos+1
       ikI(2) = neg
       ikF(2) = 2*nmax
       do ii=1,2 ! los num de onda positivos, negativos
       do ik = ikI(ii),ikF(ii) !  todos los num de onda
       k = k_vec(ik)
       do e=1,N+1 ! estrato que contiene la fuerza

         DEN = 4.0*PI*RHO(e)*cOME**2.0
!        omeAlf = cOME**2.0/ALFA(e)**2.0
!        omeBet = cOME**2.0/BETA(e)**2.0
         L2M = LAMBDA(e) + 2.0*AMU(e)

!         ! algunas valores constantes para todo el estrato
!         gamma = sqrt(omeAlf - k**2.0)
!         nu = sqrt(omeBet - k**2.0)
!         ! Se debe cumplir que la parte imaginaria del número de onda
!         ! vertical debe ser menor que cero. La parte imaginaria contri-
!         ! buye a tener ondas planas inhomogéneas con decaimiento expo-
!         ! nencial a medida que z crece.
!         if(aimag(gamma).gt.0.0_8)gamma = conjg(gamma)
!         if(aimag(nu).gt.0.0_8)nu=conjg(nu)
          gamma = gamma_E(ik,e)
          nu = nu_E(ik,e)

         ! the green function indexes go for the part
         ! that is multiplied by exp(gamma) : 1
         ! and                by exp(nu) : 2
         G31(1) = -UI/DEN * k
         G31(2) =  UI/DEN * k

         G11(1) = -UI/DEN * k**2.0/gamma
         G11(2) = -UI/DEN * nu

         S331(1) = -UR/DEN * (k*gamma*L2M + lambda(e)*k**3.0/gamma)
         S331(2) = -UR/DEN * (-2.0*amu(e)*k*nu)

         S131(1) = -UR/DEN * amu(e)* (2.0*k**2.0)
         S131(2) = -UR/DEN * amu(e)* (nu**2.0-k**2.0)

         G33(1) = -UI/DEN * gamma
         G33(2) = -UI/DEN * (k**2.0/nu)

         s333(1) = -UR/DEN * (gamma**2.0*L2M + k**2.0*lambda(e))
         s333(2) = -UR/DEN * (2.0*amu(e)*k**2.0)

         S313(1) = -UR/DEN * amu(e) * (2.0*k*gamma)
         S313(2) =  UR/DEN * amu(e) * (k/nu*(nu**2.0-k**2.0))

       do Ga_o_Nu = 1,2
         if (Ga_o_Nu .eq. 1) this_B=>BparaGa
         if (Ga_o_Nu .eq. 2) this_B=>BparaNu
         ! para dir 1
         dir = 1
       if (e .ne. 1) then
        this_B(1+4*(e-1)-2,e,ik,dir) = G31(Ga_o_Nu) * (-1)!  w
        this_B(1+4*(e-1)-1,e,ik,dir) = G11(Ga_o_Nu)!  u
       end if
        this_B(1+4*(e-1)  ,e,ik,dir) = S331(Ga_o_Nu)! szz
        this_B(1+4*(e-1)+1,e,ik,dir) = S131(Ga_o_Nu) * (-1)! szx   ! delta
      !                     =      (2) interfaz de abajo
       if (e .ne. N+1) then
        this_B(1+4*(e-1)+2,e,ik,dir) = - G31(Ga_o_Nu)* (1)!  w
        this_B(1+4*(e-1)+3,e,ik,dir) = - G11(Ga_o_Nu)!  u
        this_B(1+4*(e-1)+4,e,ik,dir) = - S331(Ga_o_Nu)! szz
        this_B(1+4*(e-1)+5,e,ik,dir) = - S131(Ga_o_Nu) * (1)! szx
       end if

         ! para dir 3
         dir = 2
       if (e .ne. 1) then
        this_B(1+4*(e-1)-2,e,ik,dir) = G33(Ga_o_Nu)!  w
        this_B(1+4*(e-1)-1,e,ik,dir) = G31(Ga_o_Nu)* (-1)!  u
       end if
        this_B(1+4*(e-1)  ,e,ik,dir) = S333(Ga_o_Nu)* (-1)! szz   ! delta
        this_B(1+4*(e-1)+1,e,ik,dir) = S313(Ga_o_Nu)! szx
      !                     =    (2) interfaz de abajo
       if (e .ne. N+1) then
        this_B(1+4*(e-1)+2,e,ik,dir) = - G33(Ga_o_Nu)!  w
        this_B(1+4*(e-1)+3,e,ik,dir) = - G31(Ga_o_Nu)* (1)!  u
        this_B(1+4*(e-1)+4,e,ik,dir) = - S333(Ga_o_Nu)* (1)! szz
        this_B(1+4*(e-1)+5,e,ik,dir) = - S313(Ga_o_Nu)! szx
       end if

       end do !Ga_o_Nu
       end do ! e
       end do ! ik
       end do !ii
end subroutine PSVpaGANU

subroutine PSVMatAporGaNU(J) ! multiplicar marices A y (casi) B
      !use modelVars, only : vecNK, nmax,SpliK
      use modelVars!, only:N
      use refSolMatrixVars, only : BparaGa,BparaNu,CoefparGa,CoefparNu,Ak
      implicit none
      integer, intent(in) :: J
      complex*16, dimension(:,:,:,:),pointer :: this_B
      complex*16, dimension(:,:,:,:,:),pointer ::this_coef
      integer :: ii,ik,e,Ga_o_Nu,dir,pos,neg,ikI(2),ikF(2),tam,&
                coIu,coFu,coId,coFd
       coId=0;coFd=0
       tam = 4*N+2
       pos = min(int(vecNK(J)*SpliK),nmax); neg = 2*nmax-(pos-2)
       ikI(1) = 1
       ikF(1) = pos+1
       ikI(2) = neg
       ikF(2) = 2*nmax
       do ii=1,2 ! los num de onda positivos, negativos
       do ik = ikI(ii),ikF(ii) !  todos los num de onda
       do Ga_o_Nu = 1,2
         if (Ga_o_Nu .eq. 1) then
           this_B=>BparaGa ;  this_coef=>CoefparGa ; end if!
         if (Ga_o_Nu .eq. 2) then
           this_B=>BparaNu ;  this_coef=>CoefparNu ; end if
         do dir = 1,2 ! para cada dirección de la fuerza
           do e=1,N+1 ! estrato que contiene la fuerza
               coIu = 1+4*(e-1)   ! sz
             if (e .ne. 1) then
               coIu = 1+4*(e-1)-2 ! w
               coFu = 1+4*(e-1)-1 ! u
             end if ! e!= 1
               coFu = 1+4*(e-1)+1 ! sx
             if (e .ne. N+1) then
               coId = 1+4*(e-1)+2 ! w
!              1+4*(e-1)+3        ! u
!              1+4*(e-1)+4        ! sz
               coFd = 1+4*(e-1)+5 ! sx
             end if ! e!= HS
!            print*,ii,ik,Ga_o_nu,dir,e;print*,"   ",coIu,coFu,coId,coFd
               this_coef(1:tam,e,ik,dir,1) = matmul(Ak(1:tam,coIu:coFu,ik),&
                  this_B(coIu:coFu,e,ik,dir))
             if (e .ne. N+1) then
               this_coef(1:tam,e,ik,dir,2) = matmul(Ak(1:tam,coId:coFd,ik),&
                  this_B(coId:coFd,e,ik,dir))
             end if
           end do ! e
         end do ! dir
       end do ! Ga_o_Nu
       end do ! ik
       end do ! ii
!      print*,CoefparGa(:,1,10,1,1);stop "PSVMatAporGaNU"
end subroutine PSVMatAporGaNU

      ! ****  la fuente cilíndrica está entre interfaces de los estratos ****
subroutine  eGAeNU(i_zF,ik,pXI,dj)
      !hacer sincGamma y sincNu para cada interfaz y devolver B final
        use modelVars, only : k_vec,gamma_E,nu_E!cOME,
        use resultVars, only : Punto
        use modelVars
        use sourceVars, only: Po,iFte=>currentiFte
        use glovars, only : z0,UR,UI
        use refSolMatrixVars, only : B, CoefparGa, CoefparNu
        implicit none
        integer :: ik,i_zF,dj
        type(Punto), pointer :: pXi
        integer, pointer :: e_f
        real*8, pointer :: z_f
        logical, pointer :: fisInterf
        integer :: iIf
        complex*16 :: gamma,nu!,argum,sincmod

      ! una para cada interfaz (la de arriba [1] y la de abajo [2])
        real*8,     dimension(2) :: z_loc
        complex*16, dimension(2) :: egamz,enuz,gamz,nuz
        complex*16, dimension(2) :: sincGamma, sincNu
        real*8 :: k!,a
        !real*8, pointer :: cose,seno
        integer :: el_tipo_de_fuente,tam

      sincGamma=z0; sincNu=z0;z_loc=0
      tam = 4*N+2
      k = k_vec(ik)
      e_f => pXi%layer
      z_f => pXi%center%z
      fisInterf => pXi%isOnInterface
      if (fisInterf) stop "eGAeNU: force on interface"
      if (e_f .ne. 0) z_loc(1) = Z(e_f) - z_f !downward (-)    (interfaz de arriba)
      if (e_f .ne. N+1)  z_loc(2) = Z(e_f+1) - z_f !upward (+) (interfaz de abajo )

      el_tipo_de_fuente = 2 ! fuente segmento (para ibem)
      if ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "eGAeNU iFte=0"
      if (i_zF .eq. 0) el_tipo_de_fuente = Po(iFte)%tipofuente
      gamma = gamma_E(ik,e_f)
      nu = nu_E(ik,e_f)

      do iIf = 1,2
          egamz(iIf) = exp(-UI*gamma*ABS(z_loc(iIf)))
          enuz(iIf) = exp(-UI*nu*ABS(z_loc(iIf)))

          gamz(iIf) = gamma*ABS(z_loc(iIf))
          nuz(iIf) = nu*ABS(z_loc(iIf))
      end do
      if (el_tipo_de_fuente .le. 1) then ! fuente puntual 0 u onda plana 1
        sincGamma(1) = UR * egamz(1)
        sincGamma(2) = UR * egamz(2)
        sincNu(1) = UR * enuz(1)
        sincNu(2) = UR * enuz(2)
!      elseif (el_tipo_de_fuente .eq. 2) then ! la fuente es un segmento
!      ! ahora estamos involucrando información del receptor en el vector de fuente
!      ! en particular, si el receptor (la interfaz) está arriba o abajo de la fuente.
!        A = pxi%length * 0.5_8  ! 2a=lenght
!        cose => pxi%cosT
!        seno => pxi%sinT
!
!      ! en cada interfaz
!        ! caso (receptor abajo de la fuente) : interfaz de abajo (2)
!        argum = gamma*A*seno + UR * (k*A*cose)
!        sincGamma(2) = sincmod(argum,gamz(2))*(2*A)
!
!        argum = nu*A*seno + UR * (k*A*cose)
!        sincNu(2) = sincmod(argum,nuz(2))*(2*A)
!
!        ! caso (receptor arriba de la fuente) : interfaz de arriba (1)
!        argum = -gamma*A*seno + UR * (k*A*cose)
!        sincGamma(1) = sincmod(argum,gamz(1))*(2*A)
!
!        argum = -nu*A*seno + UR * (k*A*cose)
!        sincNu(1) = sincmod(argum,nuz(1))*(2*A)

      end if !fuente tipo 0 o 1

      !#< r obtener vector de Coeficientes de ondas final !#>
      ! multiplicar por los coeficientes de cada interfaz para Ga y Nu
      B(1:tam,ik) = CoefparGa(1:tam,e_f,ik,dj,1)* sincGamma(1) + &
                    CoefparNu(1:tam,e_f,ik,dj,1)* sincNu(1)
      if (e_f .ne. N+1) then
      B(1:tam,ik) = B(1:tam,ik) + &
                    CoefparGa(1:tam,e_f,ik,dj,2)* sincGamma(2) + &
                    CoefparNu(1:tam,e_f,ik,dj,2)* sincNu(2)
      end if
!     print*,"eGAeNU",ik,sum(B(1:tam,ik))!;stop
end subroutine  eGAeNU

! G_stra - termIndSH
subroutine SHvectorB_ondaplana(this_B,gamma)
      use modelVars, only : n,amu0,amu,beta0,Z,beta
      use glovars, only:UR,UI,z0,PWfrecReal
      use modelVars, only : cOME,ome
      !use debugStuff
      use sourceVars, only: Po,iFte=>currentiFte
      implicit none
      complex*16, intent(inout), dimension(1:4*N+2) :: this_B
      real*8, intent(in) :: gamma
      integer :: e
      complex*16 :: kx,kz,V,am,c
      real*8 :: z_loc

      e = Po(iFte)%layer

      if (PWfrecReal) then
        c = beta0(e) !SV
        kx = ome/c * sin(gamma)
!       kz = ome/c * cos(gamma)
        kz = sqrt((ome/c)**2 - kx**2)
        if(aimag(kz).gt.0.0_8) kz=conjg(kz)
        am = AMU0(e)
      else
        c = beta(e) !SV
        kx = come/c * sin(gamma)
!       kz = come/c * cos(gamma)
        kz = sqrt((come/c)**2 - kx**2)
        if(aimag(kz).gt.0.0_8) kz=conjg(kz)
        am = AMU(e)
      end if

      this_B(1:2*N+1) = Z0
      if (Z(0) .lt. 0.0) then ! Semiespacio en z<0 ···········
!       i = 1                                                !
!       this_B(1+2*(e-1)-1 + i) = V !  v                     !
        stop "SHvectorB_ondaplana Semiespacio arriba no implementado"
      end if                                                 !
      ! ······················································

       z_loc = Z(N+1) ! La incidencia en la última interfaz

       ! Tracciones libres en z=0
       V = UR * exp(UI * kz * (z_loc))
       this_B(1+2*(e-1)  ) = UI * am * V * (- kx * 0 + kz * 1)

      ! Continuidad en z = z1
      if (N .ge. 1) then
       this_B(1+2*(e-1)-1) = V
      end if
end subroutine SHvectorB_ondaplana

subroutine SHvectorB_force(i_zF,this_B,tam,pXi,cOME,k,ik)
      use modelVars, only : N,Z,AMU!,LAMBDA,RHO!,NPAR,BETA,ALFA
      use gloVars, only : UR,UI,PI,Z0
!     use debugStuff
      use resultvars, only : Punto
      use sourceVars, only: Po,iFte=>currentiFte! use sourceVars, only : tipofuente!,PW_theta
      use modelVars, only : nu_E
      implicit none

      integer, intent(in) :: i_zF,tam,ik
      complex*16,    intent(inout), dimension(tam) :: this_B
      real*8,     intent(in)    :: k
      complex*16, intent(in)    :: cOME
      type(Punto),intent(in),target    :: pXi

      integer, pointer :: e_f
      real*8, pointer :: z_f!,seno!,cose
      logical, pointer :: fisInterf
      integer :: iIf,nInterf, el_tipo_de_fuente
      real    :: SGNz
      complex*16 :: nu,DEN, argum!,sincmod
      real*8     :: errT = 0.0001_8
      ! una para cada interfaz (la de arriba [1] y la de abajo [2])
      real*8,     dimension(2) :: z_loc
      complex*16, dimension(2) :: G22,S22,enuz,nuz,sincNu
      real*8 :: a
      a = k ! no current use of dummy variable k
      argum = cOME ! (ana mas para que no chiste)
      this_B = Z0
      sincNu = Z0
      e_f => pXi%layer
      z_f => pXi%center%z
      fisInterf => pXi%isOnInterface
      nInterf = 2
      if (fisInterf) then
         nInterf = 1
         go to 459
      end if
      z_loc(1:2) = 0.0_8
      if (e_f .ne. 0) z_loc(1) = Z(e_f) - z_f !downward (-)
      if (e_f .ne. N+1) z_loc(2) = Z(e_f+1) - z_f !upward (+)

      el_tipo_de_fuente = 2 ! fuente segmento (para ibem)
      if ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "SHvectorForce iFte=0"
      if (i_zF .eq. 0) el_tipo_de_fuente = Po(iFte)%tipofuente !(puntual u onda plana)
      !  0  puntual
      !  2  segmento

      G22=Z0;S22=z0
      DEN = (4.0*AMU(e_f)*UI*PI)
!     omeBet = cOME**2.0/BETA(e_f)**2.0
!     nu = sqrt(omeBet - k**2.0)
!     if(aimag(nu).gt.0.0_8)nu= conjg(nu)
      nu = nu_E(ik,e_f)

      do iIf = 1,2
          enuz(iIf) = exp(-UI*nu*ABS(z_loc(iIf)))
          nuz(iIf) = nu*ABS(z_loc(iIf))
      end do

      if (el_tipo_de_fuente .le. 1) then ! fuente puntual 0 u onda plana 1
        sincNu(1) = UR * enuz(1)
        sincNu(2) = UR * enuz(2)
!      elseif (el_tipo_de_fuente .eq. 2) then ! la fuente es un segmento
!      ! ahora estamos involucrando información del receptor en el vector de fuente
!      ! en particular, si el receptor (la interfaz) está arriba o abajo de la fuente.
!
!        A = pxi%length * 0.5_8 ! 2a=lenght
!        cose => pxi%cosT
!        seno => pxi%sinT
!
!      ! en cada interfaz
!        ! caso (receptor abajo de la fuente) : interfaz de abajo (2)
!        argum = nu*a*seno + UR * (k*a*cose)
!        sincNu(2) = sincmod(argum,nuz(2))*A*2.0_8
!!       !sincNu(2) = enuz(2) * (sin(k*A)/(k*A)*(2*A)) !no se inclina
!
!        ! caso (receptor arriba de la fuente) : interfaz de arriba (1)
!        argum = -nu*a*seno + UR * (k*a*cose)
!        sincNu(1) = sincmod(argum,nuz(1))*A*2.0_8
!!       !sincNu(1) = enuz(1) * (sin(k*A)/(k*A)*(2*A)) !no se inclina
      end if !fuente tipo 2

      ! en cada interfaz (1 arriba) y (2 abajo)
      do iIf = 1,2! nInterf
          if (abs(z_loc(iIf)) > errT ) then
            SGNz = real(z_loc(iIf) / ABS(z_loc(iIf)),4)
          else
            SGNz = 0.0
          end if

          G22(iIf) = UR/DEN * sincNu(iIf) / nu
          S22(iIf) = -UR / (4.0_8*pi) * sincNu(iIf) * SGNz

      end do !iIf interf
  459 if (fisInterf) then
          S22(1) = + cmplx(1.0 / (2.0 * PI ),0.0,8)
          G22(1) = 0
      end if

      ! El vector de términos independientes genera el campo difractado
      if (z(0) .gt. 0.0) then !--------------------------
      !                     =      (1) interfaz de arriba
       if (e_f .ne. 1) then
        this_B(1+2*(e_f-1)-1) = G22(1)
       end if
        this_B(1+2*(e_f-1)  ) = S22(1)

       if (.not. fisInterf) then ! la fuerza no calló en la interfaz
      !                     =      (2) interfaz de abajo
        if (e_f .ne. N+1) then
         this_B(1+2*(e_f-1)+1) = - G22(2)
         this_B(1+2*(e_f-1)+2) = - S22(2)
        end if
       end if
      else !Z(0) < 0 un semiespacio arriba
      !                     =      (1) interfaz de arriba
       if (e_f .ne. 0) then
        this_B(1+2*(e_f-1)  ) = G22(1)
        this_B(1+2*(e_f-1)+1) = S22(1)
       end if !
      !                     =      (2) interfaz de abajo
       if (.not. fisInterf) then ! la fuerza no calló en la interfaz
        if (e_f .ne. N+1) then
         this_B(1+2*(e_f-1)+2) = - G22(2)
         this_B(1+2*(e_f-1)+3) = - S22(2)
        end if
       end if
      end if
!     print*,"e_f=",e_f
!     print*,this_B;stop 3759
end subroutine SHvectorB_force

FUNCTION SINCMOD(ARG,ARGZ)
      use glovars, only : UI
      COMPLEX*16 :: SINCMOD,ARG, ARGZ
      SINCMOD=EXP(-UI*ARGZ)
      IF(ABS(ARG).LE.0.0001_8)RETURN
      SINCMOD=(EXP(UI*(ARG-ARGZ))-EXP(-UI*(ARG+ARGZ)) )/2.0_8/UI/ARG
END function SINCMOD
! G_stra - results
      ! coeficientes de las ondas planas emitidias en cada interfaz
      ! para representar el campo difractado por estratigrafía
function PSVdiffByStrata(coefOndas_PSV,z_i,e,cOME_i,k,ik)
      use modelVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : UI,UR,Z0,PWfrecReal!,PrintNum,verbose
      use resultVars, only : MecaElem
      use modelVars, only : gamma_E,nu_E,ome
      use refSolMatrixVars, only : subMatD0,subMatS0
      implicit none

      complex*16, dimension(1:5) :: PSVdiffByStrata
      real*8, intent(in)           ::  z_i,k
      complex*16, intent(in)       :: cOME_i
      integer, intent(in)          :: e,ik
      complex*16, dimension(1:4*N+2),intent(in) :: coefOndas_PSV
      complex*16 :: gamma,nu,xi,eta
      complex*16 :: egammaN,enuN,egammaP,enuP
      complex*16, dimension(2,4) :: subMatD
      complex*16, dimension(3,4) :: subMatS
      complex*16, dimension(1:4) :: coeffsPSV
      integer :: i

      PSVdiffByStrata = 0
      if (ik .ne. 0) then
       gamma = gamma_E(ik,e)
       nu = nu_E(ik,e)
       xi = k**2.0 - nu**2.0
       eta = 2.0*gamma**2.0 - cOME_i**2.0 / BETA(e)**2.0
!     print*,ik,gamma,nu
      else !ik = 0
       if (PWfrecReal) then
          gamma = sqrt(UR*OME**2.0_8/ALFA0(e)**2.0_8 - k**2.0_8)
          nu = sqrt(UR*OME**2.0_8/BETA0(e)**2.0_8 - k**2.0_8)
          if(aimag(gamma).gt.0.0)gamma = conjg(gamma)
          if(aimag(nu).gt.0.0)nu= conjg(nu)
          xi = k**2.0 - nu**2.0
          eta = 2.0*gamma**2.0 - OME**2.0 / BETA0(e)**2.0
       else
          gamma = sqrt(cOME_i**2.0_8/ALFA(e)**2.0_8 - k**2.0_8)
          nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
          if(aimag(gamma).gt.0.0)gamma = conjg(gamma)
          if(aimag(nu).gt.0.0)nu= conjg(nu)
          xi = k**2.0 - nu**2.0
          eta = 2.0*gamma**2.0 - cOME_i**2.0 / BETA(e)**2.0
       end if
      end if
      !downward waves
        if (e /= 0) then !(radiation condition upper HS)
          egammaN = exp(-UI * gamma * (z_i-Z(e)))
          enuN = exp(-UI * nu * (z_i-Z(e)))
        else
          egammaN = Z0
          enuN = Z0
        end if
      !upward waves
        if (e /= N+1) then !(radiation condition)
          egammaP = exp(UI * gamma * (z_i-Z(e+1)))
          enuP = exp(UI * nu * (z_i-Z(e+1)))
        else
          egammaP = Z0
          enuP = Z0
        end if

      !la matrix diagonal
      !coeficientes de las ondas en el estrato
      i = 0
      if (z(0) .lt. 0.0) i = 2
        if (e .eq. 0) then! semiespacio de arriba
          coeffsPSV(1:2) = (/z0,z0/)
          coeffsPSV(3:4) = coefOndas_PSV(1 : 2)
        elseif (e .eq. N+1) then ! semiespacio de abajo
          coeffsPSV(1:2) = coefOndas_PSV(4*(e-1)+1+i : 4*(e-1)+2+i)
          coeffsPSV(3:4) = (/z0,z0/)
        else! estrato
          coeffsPSV(1:4) = coefOndas_PSV(4*(e-1)+1+i : 4*(e-1)+4+i)
        end if

      ! desplazamientos
        ! {W}
        ! {U}
        subMatD = subMatD0(1:2,1:4,e,ik)
!       print*,ik,e,sum(subMatD0(1:2,1:4,e,ik))
        subMatD(:,1) = subMatD(:,1)*egammaN*coeffsPSV(1)
        subMatD(:,2) = subMatD(:,2)*enuN*coeffsPSV(2)
        subMatD(:,3) = subMatD(:,3)*egammaP*coeffsPSV(3)
        subMatD(:,4) = subMatD(:,4)*enuP*coeffsPSV(4)

        PSVdiffByStrata(1) = sum(subMatD(1,:)) !W
        PSVdiffByStrata(2) = sum(subMatD(2,:)) !U
!       print*,ik,"e",egammaN,enuN
!       print*,ik,"e",egammaP,enuP
!       print*,ik,coeffsPSV(1)
!       print*,ik,coeffsPSV(2)
!       print*,ik,coeffsPSV(3)
!       print*,ik,coeffsPSV(4)
!       print*,ik,"(1)_",PSVdiffByStrata(1)
!       print*,ik,"(2)_",PSVdiffByStrata(2)


      ! esfuerzos
        subMatS = subMatS0(1:3,1:4,e,ik)
        subMatS(:,1) = subMatS(:,1)*egammaN*coeffsPSV(1)
        subMatS(:,2) = subMatS(:,2)*enuN*coeffsPSV(2)
        subMatS(:,3) = subMatS(:,3)*egammaP*coeffsPSV(3)
        subMatS(:,4) = subMatS(:,4)*enuP*coeffsPSV(4)

        PSVdiffByStrata(3) = sum(subMatS(1,:)) !s33
        PSVdiffByStrata(4) = sum(subMatS(2,:)) !s31
        PSVdiffByStrata(5) = sum(subMatS(3,:)) !s11

!       if (sum(abs(PSVdiffByStrata)) .gt. 10000) then
!       print*,ik,e
!       print*,"coeffsPSV ",coeffsPSV(1:4) !ok
!       print*,"subMatD0(1:2,1:4,e,ik)",subMatD0(1:2,1:4,e,ik)
!       print*,"subMatS0(1:3,1:4,e,ik)",subMatS0(1:3,1:4,e,ik)
!       print*,"egammaN",egammaN
!       print*,"enuN",enuN
!       print*,"egammaP",egammaP
!       print*,"enuP",enuP
!       print*,"   W,U",PSVdiffByStrata(1:2)
!       print*,"   str",PSVdiffByStrata(3:5)
!
!       call system("echo -e ""\033[31m problem \033[0m here""")
!       end if
end function PSVdiffByStrata

function SHdiffByStrata(coefOndas_SH,z_i,e,cOME_i,k,ik)
      use modelVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : UI,UR,Z0,PWfrecReal
      use modelVars, only : nu_E,ome
!     use resultVars, only : MecaElem

      implicit none
      complex*16, dimension(1:3)   :: SHdiffByStrata
      real*8, intent(in)           :: z_i,k
      complex*16, intent(in)       :: cOME_i
      integer, intent(in)          :: e,ik
      complex*16, dimension(2*N+1),intent(in) :: coefOndas_SH
      complex*16 :: nu
      complex*16 :: enuN,enuP
      complex*16, dimension(1,2) :: subMatDsh
      complex*16, dimension(2,2) :: subMatSsh
      complex*16, dimension(2,2) :: diagMatSH
      complex*16, dimension(2,1) :: coeffsSH
      complex*16, dimension(1,1) :: resDsh
      complex*16, dimension(2,1) :: resSsh
      integer :: i
      SHdiffByStrata = z0
      ! algunas valores constantes para todo el estrato
      !nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
      !if(aimag(nu).gt.0.0)nu= conjg(nu)
      if (ik .ne. 0) then
       nu = come_i ! (nada mas para que no chiste)
       nu = nu_E(ik,e)
      else ! ik = 0
       if (PWfrecReal) then
          nu = sqrt(UR*OME**2.0_8/BETA0(e)**2.0_8 - k**2.0_8)
          if(aimag(nu).gt.0.0)nu= conjg(nu)
       else
          nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
          if(aimag(nu).gt.0.0)nu= conjg(nu)
       end if
!      nu = sqrt(UR*OME**2.0_8/BETA0(e)**2.0_8 - k**2.0_8)
!      if(aimag(nu).gt.0.0)nu= conjg(nu)
      end if
          !downward waves
          if (e /= 0) then !(radiation condition upper HS)
             enuN = exp(-UI * nu * (z_i-Z(e)))! enuN = exp(-UI * nu * abs(z_i-Z(e)))
          else
            enuN = Z0
          end if
          !upward waves
          if (e /= N+1) then !(radiation condition)
            enuP = exp(UI * nu * (z_i-Z(e+1)))! enuP = exp(-UI * nu * abs(z_i-Z(e+1)))
          else
            enuP = Z0
          end if

          !la matrix diagonal
          diagMatSH = RESHAPE((/ enuN, Z0,   &
                                 Z0,   enuP /), &
                              (/ 2,2 /))
      !coeficientes de las ondas en el estrato
      i = 0
      if (z(0) .lt. 0.0) i = 1
        if (e .eq. 0) then! semiespacio de arriba
          coeffsSH(1:1,1) = z0
          coeffsSH(2:2,1) = coefOndas_SH(1:1)
        elseif (e .eq. N+1) then! semiespacio de abajo
          coeffsSH(1:1,1) = coefOndas_SH(2*(e-1)+1+i : 2*(e-1)+1+i)
          coeffsSH(2:2,1) = z0
        else! estrato
          coeffsSH(1:2,1) = coefOndas_SH(2*(e-1)+1+i : 2*(e-1)+2+i)
        end if

      ! desplazamientos
!     if (mecStart .eq. 1)then
        subMatDsh = RESHAPE((/ UR,UR /), (/ 1,2 /))
        subMatDsh = matmul(subMatDsh,diagMatSH)
        resDsh = matmul(subMatDsh, coeffsSH)
        SHdiffByStrata(1) = resDsh(1,1) !V
!     end if ! desplazamientos

      ! esfuerzos
!     if (mecEnd .eq. 3) then
        subMatSsh = RESHAPE((/ -UI*nu,-UI*k,UI*nu,-UI*k/),(/2,2/))
        subMatSsh = amu(e) * subMatSsh
        subMatSsh = matmul(subMatSsh,diagMatSH)
        resSsh = matmul(subMatSsh, coeffsSH)
        SHdiffByStrata(2) = resSsh(1,1) !s32
        SHdiffByStrata(3) = resSsh(2,1) !s12
!     end if ! esfuerzos

end function SHdiffByStrata

! G_full PSV
Subroutine FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,mecS,mecE)
!#define ver 1
      ! Sanchez-Sesma y Campillo, 1991.  Mismo resultado que:
      ! Kaussel, Fundamental solutions in elastodynamics... pag 38
      use modelVars ,only : alfa0,beta0,Lambda0,AMU0,alfa,beta,amu,lambda,rho!,N!,minWL
      use gloVars, only:UI,UR,one,z0,PWfrecReal
      !use hank !     use specfun
      use resultvars, only : Punto,FFres
      use sourceVars, only: Po,iFte=>currentiFte
      use modelVars, only : OME
      !use Gquadrature, only : Gquad_n,Gquad_f
      implicit none
!      interface
!        subroutine greenexPSV(greenex,dir_j,p_X,pXi,estrato,cOME)
!         use resultvars, only : Punto,FFres
!         type (FFres), intent(out) :: greenex
!         integer,intent(in) :: dir_j, estrato
!         type(Punto),intent(in), pointer :: p_X,pXi
!         complex*16,intent(in) :: cOME
!        end subroutine greenexPSV
!      end interface
      type (FFres), intent(out) :: FF
      integer,    intent(in)  :: mecS,mecE, i_zF
      type(Punto),intent(in), pointer :: p_X,pXi
      complex*16, intent(in) :: cOME
      integer, intent(in) :: dir_j
      real*8 :: r,gamma(2)
      complex*16 :: A,B,C,Dqr,Dkr,kx,kz
      complex*16 :: omeP,omeS
      complex*16 :: H0s,H1s,H2s,H0p,H1p,H2p !Hankel
      complex*16 :: szz,szx,sxx
      complex*16 :: la,am
      integer :: i,j,iGqDistan
      integer, pointer :: e
      real*8 :: nX(2)
      integer :: iGq,nGq
      real*8, pointer :: xf,zf,GqC
      !real*8 :: deltaij !se omite la declaracion porque esta en elmismo modulo
      real*8,dimension(2) :: theta
      integer :: el_tipo_de_fuente
      integer, target :: estrato
      logical :: shouldI,XinoEstaEnInterfaz,usarGreenex,estratosIguales
#ifdef ver
      print*,""
      print*,"-------------------------------------------"
      print*,"i_zF=",i_zF
      print*,"px",p_x%center,"e=",p_x%layer
      print*,"pxi",pXi%center,"e=",pXi%layer
      print*,"dir_j=",dir_j
#endif
      FF%W=z0; FF%U=z0;FF%Tz=z0;FF%Tx=z0;e=>p_X%layer
      FF%sxx = z0;FF%szx = z0;FF%szz = z0
      estratosIguales = .false.
      XinoEstaEnInterfaz = .false.
      usarGreenex = .false.
      shouldI = .false.
      c = z0; theta(1) = 0; theta(2) = 0
      if (p_x%layer .eq. pXi%layer) estratosIguales = .true.
      if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
      j = dir_j ! <***********  dir_j = 3  (vertical)
      if(j .eq. 3) j = 2 ! para coincidir con los indicies
      el_tipo_de_fuente = 2 !(fuente segmento)
      nx(1) = p_X%normal%x; nx(2) = p_X%normal%z
      xf => one; zf => one

      if (i_zF .eq. 0) then ! es la fuente real
         if ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "FFPSV iFte=0" !failsafe
         el_tipo_de_fuente = Po(iFte)%tipofuente !(puntual:0, onda plana:1, segmento:2)
!         if (pXi%region .eq. 2) then !en la inclusión R
!            shouldI = .true.
!            XinoEstaEnInterfaz = .true.
!            estrato = N+2
!            e => estrato
!         else !pXi%region = 1 E
            if (p_x%layer .eq. pXi%layer) shouldI = .true.
            if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
            estrato = p_x%layer
            e => p_x%layer
!         end if
!      else ! i_zf .ne. 0  una fuente virtual de ibem
!         el_tipo_de_fuente = 2 !(fuente segmento)
!
!         !#< r XinoEstaEnInterfaz = .true. !(se ha encargado la geometría)
!         if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
!         !#>
!         if (i_zF .eq. -1) then !(campo refractado en inclusion R)
!            shouldI = .true.
!            estrato = N+2 !(para tomar las propiedades de la inclusión)
!            e => estrato
!         else !(normal: campo refractado en medio estratificado E)
!            if (p_x%layer .eq. pXi%layer) then
!               shouldI = .true. ! si están en el mismo estrato
!               estrato = p_x%layer
!               e => p_x%layer
!            end if
!         end if
      end if
#ifdef ver
      print*,"shouldI=",ShouldI," estrato=",estrato
      print*,"el_tipo_de_fuente=",el_tipo_de_fuente
#endif
      if (shouldI) then
      if ((i_zF .eq. 0) .and. (el_tipo_de_fuente .eq. 1)) then ! onda plana !
#ifdef ver
      print*,"onda plana Pol=",Po(iFte)%PW_pol
#endif
!        if (p_x%isOnInterface .eqv. .true.) return  ! creo
         if (Po(iFte)%PW_pol .eq. 1) then !SV
            if (PWfrecReal) then
               c = UR*beta0(e)
            else
               c = beta(e)
            end if
            theta(1) = cos(Po(iFte)%gamma)
            theta(2) = sin(Po(iFte)%gamma)
         elseif (Po(iFte)%PW_pol .eq. 2) then
            if (PWfrecReal) then !P
                c = UR*alfa0(e)
            else
                c = alfa(e) !SV
            end if
            theta(1) = sin(Po(iFte)%gamma)
            theta(2) = -cos(Po(iFte)%gamma)
         end if!
         if (PWfrecReal) then
            kx = UR*real(ome/c * sin(Po(iFte)%gamma))
!           kz = UR*real(ome/c * cos(Po(iFte)%gamma))
            kz = sqrt((ome/c)**2 - kx**2)
            if(aimag(kz).gt.0.0_8) kz=conjg(kz)
            la = UR*real(LAMBDA0(e))
            am = UR*real(AMU0(e))
         else
            kx = UR*real(come/c * sin(Po(iFte)%gamma))
!           kz = UR*real(come/c * cos(Po(iFte)%gamma))
            kz = sqrt((come/c)**2 - kx**2)
            if(aimag(kz).gt.0.0_8) kz=conjg(kz)
            la = UR*real(LAMBDA(e))
            am = UR*real(AMU(e))
         end if
!       print*,"kxkzlaam",kx,kz,la,am
        FF%U = (theta(1))* exp(UI * kz * (p_x%center%z - 0)) &         !
              * exp(-UI * kx * (p_x%center%x - Po(iFte)%center%x))          !
        FF%W = (theta(2))* exp(UI * kz * (p_x%center%z - 0)) &         !
              * exp(-UI * kx * (p_x%center%x - Po(iFte)%center%x))          !
        szz = UI * ( &                                                      !
                    ( FF%W * kz * (la + 2.0* am)) &                         !
                  - ( FF%U * kx * la))                                      !
        szx = UI * am * ( kz * FF%U - kx * FF%W )                           !
        sxx = UI * ( &                                                      !
                  - ( FF%U * kx * (la + 2.0*am)) &                          !
                  + ( FF%W * kz * la))                                      !
!        FF%Tx = sxx * nx(1) + szx * nx(2)                                   !
!        FF%Tz = szx * nx(1) + szz * nx(2)                                   !
        FF%sxx = sxx                                                        !
        FF%szx = szx                                                        !
        FF%szz = szz                                                        !
#ifdef ver
      print*,ome,"ome"
      print*,come,"come"
      print*,c,"c"
      print*,kz,"kz"
      print*,kx,"kx"
      print*,la,"la"
      print*,am,"am"
      print*,theta,"theta"
      print*,FF%W,"FF%W"
      print*,FF%U,"FF%U"
      print*,FF%Tx,"FF%Tx"
      print*,FF%Tz,"FF%Tz"
      print*,FF%sxx,"FF%sxx"
      print*,FF%szx,"FF%szx"
      print*,FF%szz,"FF%szz"
      stop
#endif
        return                                                              !
      end if! fin onda plana -´-´-´-´-´-´-´-´-´-´-´-´-´-´--´-´ fin onda plana
      if (XinoEstaEnInterfaz .eqv. .true.) then
         xf => pXi%center%x
         zf => pXi%center%z
         r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.)

!        ! para funciones de Greeen en frontera con cont. de desplazamiento
!        if ((p_x%isboundary .eqv. .true.) .and. &
!            (pXi%isboundary .eqv. .true.) .and. &
!            (p_x%boundaryIndex .eq. pXi%boundaryIndex)) then
!            usarGreenex = .true.; end if
!        ! para campo cercano por fuente segmento
!        if ((i_zF .eq. 0) .and. &
!            (el_tipo_de_fuente .eq. 2) .and. &      !fuente real
!            (abs(r) .lt. pXi%length/2)) then
!            usarGreenex = .true.; end if
!        if ((pXi%isboundary .eqv. .true.) .and. &   !fuente virtual
!            (abs(r) .lt. pXi%length/2)) then
!            usarGreenex = .true.; end if
!       print*,usarGreenex
!        ! si fuente real y cilindrica entonces no
!        if (el_tipo_de_fuente .eq. 0) usarGreenex = .false.
!       print*,usarGreenex
!        ! si es un puntno para sobredeterminar el sistema entonces no
!        if (p_X%isOD) usarGreenex = .false.
!       print*,usarGreenex
!       !#< r
!       usarGreenex = .false.
!       !#>
#ifdef ver
      print*,"onda cilíndrica r=",r, "usarGreenex=",usarGreenex
#endif
      !Con integración de Gauss
      if (usarGreenex .eqv. .false.) then
         if (el_tipo_de_fuente .eq. 0) then
            nGq = 1
            GqC => ONE
            iGqDistan = 0
!         else !.eq. 2 (fuente segmento)
!            if (abs(r) .le. minWL(e)) then
!            ! si está muy cerca:
!            nGq = Gquad_n
!            iGqDistan = 1
!            else if (abs(r) .le. 2*minWL(e)) then
!            ! si está mas o menos cerca
!            nGq = Gquad_f
!            iGqDistan = 2
!            else
!            ! está lejos
!            nGq = 1
!            GqC => ONE
!            iGqDistan = 0
!            end if!
         end if
         do iGq = 1,nGq
!            if (nGq .gt. 1) then
!                xf => pXi%Gq_xXx_coords(iGq,1,iGqDistan)
!                zf => pXi%Gq_xXx_coords(iGq,2,iGqDistan)
!                GqC => pXi%Gq_xXx_C(iGq,iGqDistan)
!            end if
            r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.)
#ifdef ver
      print*,"iGq=",iGq,"xf,zf,r,Gqc=",xf,zf,r,Gqc
#endif
      gamma(1) = (p_X%center%x - xf) / r ! gamma x
      gamma(2) = (p_X%center%z - zf) / r ! gamma z
      omeP = cOME * r / alfa(e)
      omeS = cOME * r / beta(e)
      ! funcs de Hankel de segunda especie
      call hankels(omeP,H0p,H1p)
      H2p = -H0p + 2./omeP * H1p
      call hankels(omeS,H0s,H1s)
      H2s = -H0s + 2./omeS * H1s

      A = H0p/alfa(e)**2. + H0s/beta(e)**2.
      B = H2p/alfa(e)**2. - H2s/beta(e)**2.
#ifdef ver
      print*,"gamma,omeP,omeS,A,B",gamma(1),gamma(2)
      print*,omeP,omeS
      print*,A,B
#endif
      ! desplazamientos
      if (mecS .eq. 1) then
      ! W
      i = 2
      FF%W = FF%W + (-UI/8.0/rho(e)*&
      (A*deltaij(i,j)-(2.*gamma(i)*gamma(j)-deltaij(i,j))*B)) * GqC
      ! U
      i = 1
      FF%U = FF%U + (-UI/8.0/rho(e)*&
      (A*deltaij(i,j)-(2.*gamma(i)*gamma(j)-deltaij(i,j))*B)) * GqC
      end if!mecs1

      if (mecE .eq. 5) then
      nx(1) = p_X%normal%x; nx(2) = p_X%normal%z
      ! tracciones
      Dqr = omeP*H1p
      Dkr = omeS*H1s
      C = Dqr/alfa(e)**2. - Dkr/beta(e)**2.

      ! TZ
      i = 2
      FF%Tz = FF%Tz + (&
      amu(e)*UI /(2.*rho(e)*r)*(&
      (B+(lambda(e)*Dqr)/(2.*amu(e)*alfa(e)**2.))*gamma(j)*nX(i) &
      + (B + Dkr/(2.*beta(e)**2.))* &
      (gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4.*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))&
      ) * GqC

      ! TX
      i = 1
      FF%Tx = FF%Tx + ( &
      amu(e)*UI /(2.*rho(e)*r)*(&
      (B+(lambda(e)*Dqr)/(2.*amu(e)*alfa(e)**2.))*gamma(j)*nX(i) &
      + (B + Dkr/(2.*beta(e)**2.))* &
      (gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4.*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))&
      ) * GqC

      ! sxx
      i = 1
      nX(1) = 1; nX(2) = 0
      FF%sxx = FF%sxx + ( &
      amu(e)*UI /(2.*rho(e)*r)*(&
      (B+(lambda(e)*Dqr)/(2.*amu(e)*alfa(e)**2.))*gamma(j)*nX(i) &
      + (B + Dkr/(2.*beta(e)**2.))* &
      (gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4.*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))&
      ) * GqC
      ! szx
      i = 1
      nX(1) = 0; nX(2) = 1
      FF%szx = FF%szx + ( &
      amu(e)*UI /(2.*rho(e)*r)*(&
      (B+(lambda(e)*Dqr)/(2.*amu(e)*alfa(e)**2.))*gamma(j)*nX(i) &
      + (B + Dkr/(2.*beta(e)**2.))* &
      (gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4.*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))&
      ) * GqC
      ! szz
      i = 2
      nX(1) = 0; nX(2) = 1
      FF%szz = FF%szz + (&
      amu(e)*UI /(2.*rho(e)*r)*(&
      (B+(lambda(e)*Dqr)/(2.*amu(e)*alfa(e)**2.))*gamma(j)*nX(i) &
      + (B + Dkr/(2.*beta(e)**2.))* &
      (gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4.*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))&
      ) * GqC

      end if !mece5
      end do ! iGq
      else !. usarGreenex ..................................................
#ifdef ver
      print*,"Greenex"
      stop
#endif
 !       call greenexPSV(FF,dir_j,p_X,pXi,estrato,cOME) ! W,U
!     print*,FF%W,"FF%W"
!     print*,FF%U,"FF%U"
!     print*,FF%Tx,"FF%Tx"
!     print*,FF%Tz,"FF%Tz"
!     print*,"Greenex";stop
!        if (p_X%isOD) then
!          print*,p_x%center,"  ->  ",pXi%center
!          print*,"tx=",FF%Tx
!          print*,"tz=",FF%Tz
!          stop "OD in FFpsv"
!          if (dir_j .eq. 1) then; FF%Tx=0.5*UR; FF%Tz=z0;end if!
!          if (dir_j .eq. 3) then; FF%Tx=z0;     FF%Tz=0.5*UR; end if
!        end if
      end if ! greenex o gaussiana
      end if ! XinoEstaEnInterfaz: on the interface
#ifdef ver
      print*,FF%W,"FF%W"
      print*,FF%U,"FF%U"
      print*,FF%Tx,"FF%Tx"
      print*,FF%Tz,"FF%Tz"
      print*,""
!     stop "FFpsv"
#endif
      end if !should I?
end subroutine FFpsv

!  subroutine greenexPSV

! G_full SH
subroutine FFsh(i_zF,FF,p_X,pXi,cOME,mecS,mecE)
!#define ver 1
      use modelVars ,only : amu,amu0,beta,beta0!,N!,minWL
      use gloVars, only:UR,UI,z0,ONE,PWfrecReal
      !use hank !use specfun
      use resultvars, only : Punto,FFres
      use sourceVars, only: Po,iFte=>currentiFte!use sourceVars, only: tipofuente
      use modelVars, only : OME
      !use Gquadrature, only : Gquad_n,Gquad_f
      implicit none
!      interface
!        subroutine greenexSH(FF,p_X,pXi,e,cOME)
!         use resultvars, only : Punto,FFres
!         type (FFres), intent(out) :: FF
!         integer,intent(in) :: e
!         type(Punto),intent(in), pointer :: p_X,pXi
!         complex*16,intent(in) :: cOME
!        end subroutine greenexSH
!      end interface
      type (FFres), intent(out) :: FF !V,Ty
      integer,    intent(in)     :: mecS,mecE, i_zF
      type(Punto),intent(in), pointer :: p_X,pXi
      complex*16, intent(in)     :: cOME

      real*8 :: r,gamma(2)!,longonda
      complex*16 :: omeS,c,am
      complex*16 :: H0s,H1s,kx,kz,szy,sxy
      integer, pointer :: e
      integer :: iGq,nGq,iGqDistan
      real*8 ::  nX(2)
      real*8, pointer :: xf,zf,GqC
      integer :: el_tipo_de_fuente
      integer, target :: estrato
      logical :: shouldI,XinoEstaEnInterfaz,usarGreenex
#ifdef ver
      print*,""
      print*,"-------------------------------------------"
      print*,"i_zF=",i_zF
      print*,"px",p_x%center,"e=",p_x%layer
      print*,"pxi",pXi%center,"e=",pXi%layer
#endif
      FF%V=z0;FF%Ty=z0;C=z0;e=>p_X%layer
!     estratosIguales = .false.
      XinoEstaEnInterfaz = .false.
      usarGreenex = .false.
      shouldI = .false.
!     if (p_x%layer .eq. pXi%layer) estratosIguales = .true.
!     if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
!     el_tipo_de_fuente = 2 ! fuente segmento para IBEM
      nx(1) = p_X%normal%x;nx(2) = p_X%normal%z !para Tracciones
      xf => one;zf => one

      if (i_zF .eq. 0) then
         if ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "FFSH iFte=0"
         if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
         el_tipo_de_fuente = Po(iFte)%tipofuente !(puntual u onda plana)
!            if (pXi%region .eq. 2) then !en la inclusión R
!               shouldI = .true.
!               XinoEstaEnInterfaz = .true.
!               estrato = N+2
!               e => estrato
!            else ! en el exterior E
               if (p_x%layer .eq. pXi%layer) shouldI = .true.
               if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
               estrato = p_x%layer
               e => p_x%layer
!            end if
!      else !i_zf .ne. 0  una fuente virtual de ibem
!         el_tipo_de_fuente = 2 !(fuente segmento)
!         !#< r XinoEstaEnInterfaz = .true. !(se ha encargado la geometría)
!         if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
!         !#>
!         if (i_zF .eq. -1) then !(Fza distr en segmento en región R)
!!           if (p_X%region .eq. 2) then
!               !nx(1) =- p_X%normal%x;nx(2) =- p_X%normal%z
!               shouldI = .true.
!               estrato = N+2 !(para tomar las propiedades de la inclusión)
!               e => estrato
!!              print*,beta(1:N+2)
!!              print*,amu(1:N+2)
!!           end
!         else !(Fza en medio estratificado E)
!            if (p_x%layer .eq. pXi%layer) then
!               shouldI = .true. ! si están en el mismo estrato
!               estrato = p_x%layer
!               e => p_x%layer
!            end if
!         end if
      end if!
#ifdef ver
      print*,"shouldI=",ShouldI," estrato=",estrato
      print*,"el_tipo_de_fuente=",el_tipo_de_fuente
#endif
      if (shouldI) then
      if ((i_zF .eq. 0) .and. (el_tipo_de_fuente .eq. 1)) then ! onda plana !
#ifdef ver
      print*,"onda plana"
#endif
            if (PWfrecReal) then
               c = UR*beta0(e)
               kx = UR*real(ome/c * sin(Po(iFte)%gamma))
!              kz = UR*real(ome/c * cos(Po(iFte)%gamma))
               kz = sqrt((ome/c)**2 - kx**2)
               if(aimag(kz).gt.0.0_8) kz=conjg(kz)
               am = UR*real(AMU0(e))
            else
               c = beta(e)
               kx = UR*real(come/c * sin(Po(iFte)%gamma))
!              kz = UR*real(come/c * cos(Po(iFte)%gamma))
               kz = sqrt((come/c)**2 - kx**2)
               if(aimag(kz).gt.0.0_8) kz=conjg(kz)
               am = UR*real(AMU(e))
            end if
        ! las expresiones se encuentran en todos lados, por ejemplo:
        ! Gil-Zepeda et al 2003
        ! SV:
        FF%V = (UR)* exp(-UI*(kx * (p_x%center%x - Po(iFte)%center%x) - &
                              kz * (p_x%center%z - 0)))
!       szy, sxy
        sxy = am * FF%V * (-UI * kx)
        szy = am * FF%V * (UI * kz)
        FF%sxy = sxy
        FF%szy = szy
!        FF%Ty = sxy * nx(1) + szy * nx(2)
!       !#< r
!       FF%V = 0
!       FF%Ty = 0
!       print*,"lin 4713"
!       !#>
#ifdef ver
      print*,FF%V,"FF%V"
!      print*,FF%Ty,"FF%Ty"
#endif
        return
      end if! onda plana
      if (XinoEstaEnInterfaz .eqv. .true.) then !..............................
        xf => pXi%center%x
        zf => pXi%center%z
        r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.)

!        ! para funciones de Greeen en frontera con cont. de desplazamiento
!        if ((p_x%isboundary .eqv. .true.) .and. &
!            (pXi%isboundary .eqv. .true.) .and. &
!            (p_x%boundaryIndex .eq. pXi%boundaryIndex)) then
!            usarGreenex = .true.
!            end if
!!       print*,usarGreenex
!        ! para campo cercano por fuente segmento
!        if ((i_zF .eq. 0) .and. &
!            (el_tipo_de_fuente .eq. 2) .and. &      !fuente real
!            (abs(r) .lt. pXi%length/2)) then
!            usarGreenex = .true.
!            end if
!!       print*,usarGreenex
!!       if ((p_X%isboundary .eqv. .false.) .and. &
!!           (pXi%isboundary .eqv. .true.) .and. &   !fuente virtual
!!           (abs(r) .lt. pXi%length/2)) usarGreenex = .true.
!        if ((pXi%isboundary .eqv. .true.) .and. &   !fuente virtual
!            (abs(r) .lt. pXi%length/2))then
!             usarGreenex = .true.
!             end if
!
!!       print*,usarGreenex
!        ! si fuente real y cilindrica entonces no
!        if (el_tipo_de_fuente .eq. 0) usarGreenex = .false.
!!       print*,usarGreenex
!        ! si es un puntno para sobredeterminar el sistema entonces no
!        if (p_X%isOD) usarGreenex = .false.
!       print*,usarGreenex
!       !#< r
!       usarGreenex = .false.
!       !#>
#ifdef ver
      print*,"onda cilíndrica r=",r, "usarGreenex=",usarGreenex
#endif
      !Con integración de Gauss
      if (usarGreenex .eqv. .false.) then
         if (el_tipo_de_fuente .eq. 0) then                             !
           nGq = 1                                                      !
           GqC => ONE
           iGqDistan = -1000
!         else !.eq. 2 (fuente segmento)
!            if (abs(r) .le. minWL(e)) then
!            ! si está muy cerca:
!            nGq = Gquad_n
!            iGqDistan = 1
!            else if (abs(r) .le. 2*minWL(e)) then
!            ! si está mas o menos cerca
!            nGq = Gquad_f
!            iGqDistan = 2
!            else
!            ! está lejos
!            nGq = 1
!            GqC => ONE
!            iGqDistan = -1000
!            end if!                                                 !
!!        else ! eq. 2; print*,"fuente segmento "                        !
!!          if (pXi%isBoundary) nGq = Gquad_n ! IBEM                     !
!!          if (pXi%isSourceSegmentForce) nGq = Gquad_n                  !
!!          if (nGq .ne. Gquad_n) stop "chin 6600"                       !
         end if                                                         !
         do iGq = 1,nGq !..........................................     !
!            if (nGq .gt. 1) then                                  !     !
!               xf => pXi%Gq_xXx_coords(iGq,1,iGqDistan)           !     !
!               zf => pXi%Gq_xXx_coords(iGq,2,iGqDistan)           !     !
!               GqC => pXi%Gq_xXx_C(iGq,iGqDistan)                 !     !
!            end if                                                !     !
      r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.)     !     !
#ifdef ver
      print*,"iGq=",iGq,"xf,zf,r,Gqc=",xf,zf,r,Gqc
#endif
      gamma(1) = (p_X%center%x - xf) / r ! gamma x                !     !
      gamma(2) = (p_X%center%z - zf) / r ! gamma z                !     !
      omeS = cOME * r / beta(e)                                   !     !
      call hankels(omeS,H0s,H1s) ! Hankel de segunda especie      !     !
#ifdef ver
      print*,"omeS=",omeS
      print*,"H0s,H1s=",H0s,H1s
#endif
      if(mecS .eq. 1) FF%V = FF%V + (-UI/(4.*amu(e))*H0s) * GqC   !     !
      if(mecE .eq. 3) then                                        !     !
     FF%szy = FF%szy + (UI*cOME*(p_x%center%z-zf)) &             !     !
                        /(4.0_8*beta(e)*r)*H1s * C               !     !
     FF%sxy = FF%sxy + (UI*cOME*(p_x%center%x-xf))/ &            !     !
                         (4.0_8*beta(e)*r)*H1s * C               !     !
!      FF%Ty = FF%Ty + (UI*cOME/beta(e)/4.0*H1s* &                 !     !
!               (gamma(1)*nx(1) + gamma(2)*nx(2))) * GqC           !     !
      end if !mec3                                                !     !
      end do ! iGq ................................................     !
      else !. usarGreenex ...............................
#ifdef ver
      print*,"Greenex"
      stop
#endif
!        call greenexSH(FF,p_X,pXi,e,cOME)
      end if ! greenex o gaussiana .....................................!
      end if ! XinoEstaEnInterfaz: on the interface
#ifdef ver
      print*,FF%V,"FF%V"
      print*,FF%Ty,"FF%Ty"
      print*,""
!     stop "FFsh"
#endif
      end if ! should I? ....................................................
end subroutine FFsh

! subroutine greenexSH

subroutine fill_diffbyStrata(i_zf,J,auxK,come,mecS,mecE,p_x,pXi,dir_j)
      use resultvars, only:Punto,FFres!, nIpts
      use modelVars, only : NMAX,SpliK!,ome!k_vec,vecNK,,dk
!      use meshvars, only: npixX,MeshMaxX,MeshMinX
!      use wavelets !fork
!      use modelVars, only:N!beta0,beta,alfa0,alfa,
      use sourceVars, only: PoFte =>  Po, iFte => currentiFte!tipofuente, PW_pol
      use glovars, only : UR,UI,PWfrecReal
!      use peli, only : fotogramas_Region
      implicit none
!      interface
!       include 'interfazFF.f'
!      end interface
      type(Punto), pointer :: p_X,pXi!,p_Xmov
!      type(Punto),target :: p_xaux
      integer :: i_zf,mecS,mecE,J,mecaElemEnd!,po!,ne!i,ik,iMec,
      integer, intent(in) :: dir_j
      complex*16, dimension(2*nmax,mecS:mecE), target :: auxK
!      complex*16, dimension(2*nmax,mecS:mecE) :: auxKmo
!      real*8 :: mov_x
!      complex*16 :: c!,kx!,omega,bet,alf
      complex*16, intent(in)  :: cOME
      type(FFres),target :: FF
      type(FFres) :: SolRef
      real*8 :: nf(3)
!      logical :: PW
      nf(1) = pXi%normal%x
      nf(2) = 1!pXi%normal%y
      nf(3) = pXi%normal%z
!     print*,i_zf,dir_j
      mecaElemEnd = 2 !PSV
      if (dir_j .eq. 2) mecaElemEnd = 1 !SH
      if (i_zF .eq. 0) then
        if ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "fill_diffbyStra iFte=0"
        if (PoFte(iFte)%tipofuente .eq. 1) nf(1:3) = 1.0 !fuente real y onda plana
        if (abs(nf(dir_j)) .lt. 0.0001) return !fuente real fuerza; componente nulo
      end if
      ! diffracted field due to the stratification U,V,W or G
      if (dir_j .eq. 2) then ! SH
       call FFsh(i_zf,FF,p_X,pXi,cOME,1,3) !incidencia directa
       !SolRef%Ty = ((auxk(1,3)* p_x%normal%x + auxk(1,2)* p_x%normal%z) + FF%Ty)
       SolRef%sxy = FF%sxy + auxk(1,3)
       SolRef%szy = FF%szy + auxk(1,2)
       SolRef%V = (auxk(1,1) + FF%V)
!      !#< r
!      print*,p_X%center
!      print*,auxk(1,1) ,abs(auxk(1,1)),FF%V,abs(FF%V),abs(auxk(1,1)+FF%V)
!      FF%V = 0
!      auxk(1,1) = 0
!      !#>

!       if(i_zf .eq.0) then
!         print*,p_x%resp(J,iFte)%V,auxk(1,1),FF%V,nf(dir_j);stop 4898
!         if (pXi%region .ne. p_x%region) return
          p_x%resp(J,iFte)%V = &
          p_x%resp(J,iFte)%V + SolRef%V !(auxk(1,1) + FF%V) !* nf(dir_j) ! V
          p_x%resp(J,iFte)%sxy = &
          p_x%resp(J,iFte)%sxy + SolRef%sxy
          p_x%resp(J,iFte)%szy = &
          p_x%resp(J,iFte)%szy + SolRef%szy
!          p_x%resp(J,iFte)%Ty = &
!          p_x%resp(J,iFte)%Ty + SolRef%Ty !&
!         ((auxk(1,3)* p_x%normal%x + auxk(1,2)* p_x%normal%z) + FF%Ty) !* nf(dir_j)
          !                     s12                       s32
!       else
!          pXi%G(p_x%pointIndex,3,dir_j) = &
!          pXi%G(p_x%pointIndex,3,dir_j) + SolRef%V !auxK(1,1) + FF%V ! V
!          pXi%G(p_x%pointIndex,6,dir_j) = &
!          pXi%G(p_x%pointIndex,6,dir_j) + SolRef%Ty !&
!!         (auxk(1,3)* p_x%normal%x + auxk(1,2)* p_x%normal%z) + FF%Ty ! Ty
!          !                     s12                       s32
!       end if
      else !PSV
!     if (pXi%region .ne. p_x%region) return
       call FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,1,5)  !incidencia directa
!      print*,"x= ",p_x%center
!      print*,p_x%resp(J,iFte)%W
!      print*,auxK(1,1)
!      print*,FF%W
!      print*," "
!       if(i_zf .eq. 0) then
!        if (pXi%region .ne. p_x%region) return
         p_x%resp(J,iFte)%W = p_x%resp(J,iFte)%W + (auxK(1,1) + FF%W) * nf(dir_j) !W
         p_x%resp(J,iFte)%U = p_x%resp(J,iFte)%U + (auxK(1,2) + FF%U) * nf(dir_j) !U
!         p_x%resp(J,iFte)%Tz = p_x%resp(J,iFte)%Tz + &
!         ((auxk(1,4)* p_x%normal%x + auxk(1,3)* p_x%normal%z) + FF%Tz) * nf(dir_j)
!         !      s31                       s33
!         p_x%resp(J,iFte)%Tx = p_x%resp(J,iFte)%Tx + &
!         ((auxk(1,5)* p_x%normal%x + auxk(1,4)* p_x%normal%z) + FF%Tx) * nf(dir_j)
!         !      s11                       s31
         p_x%resp(J,iFte)%sxx = p_x%resp(J,iFte)%sxx + &
         (auxk(1,5) + FF%sxx) * nf(dir_j) !sxx
         p_x%resp(J,iFte)%szx = p_x%resp(J,iFte)%szx + &
         (auxk(1,4) + FF%szx) * nf(dir_j) !szx
         p_x%resp(J,iFte)%szz = p_x%resp(J,iFte)%szz + &
         (auxk(1,3) + FF%szz) * nf(dir_j) !szz
!       else
!         pXi%G(p_X%pointIndex,1,dir_j) = &
!         pXi%G(p_X%pointIndex,1,dir_j) + (auxK(1,1) + FF%W) !W
!         pXi%G(p_X%pointIndex,2,dir_j) = &
!         pXi%G(p_X%pointIndex,2,dir_j) + (auxK(1,2) + FF%U) !U
!         pXi%G(p_X%pointIndex,4,dir_j) = &
!         pXi%G(p_X%pointIndex,4,dir_j) + &
!         ((auxk(1,4)* p_x%normal%x + auxk(1,3)* p_x%normal%z) + FF%Tz) !Tz
!         pXi%G(p_X%pointIndex,5,dir_j) = &
!         pXi%G(p_X%pointIndex,5,dir_j) + &
!         ((auxk(1,5)* p_x%normal%x + auxk(1,4)* p_x%normal%z) + FF%Tx) !Tx
!
!         pXi%G(p_X%pointIndex,7,dir_j) = pXi%G(p_X%pointIndex,7,dir_j) + &
!         auxk(1,5) + FF%sxx !sxx
!         pXi%G(p_X%pointIndex,8,dir_j) = pXi%G(p_X%pointIndex,8,dir_j) + &
!         auxk(1,4) + FF%szx !szx
!         pXi%G(p_X%pointIndex,9,dir_j) = pXi%G(p_X%pointIndex,9,dir_j) + &
!         auxk(1,3) + FF%szz !szz
!       end if
      end if !dir_j
!      end if ! guardarMovieSiblings
end subroutine fill_diffbyStrata


function deltaij(i,j)
      integer :: i,j
      real*8 :: deltaij
      ! i : dirección del receptor {1;2} x,z
      ! j : dirección de la fuente {1;2} x,z
      deltaij = real(0.,8)
      if (i .eq. j) then
        deltaij = real(1.,8)
      end if
      return
end function deltaij

function porLoMenosUnoEsEstr(itabla_z)
      use resultVars, only : pota,Punto,allpoints
      implicit none
      type(Punto), pointer :: PX
      logical :: porLoMenosUnoEsEstr
      integer :: itabla_z ,itabla_x

      porLoMenosUnoEsEstr = .true.
      if (pota(itabla_z,2) .gt. 0) return !si hay una frontera

      porLoMenosUnoEsEstr = .false.
      do itabla_x = 3,2+ pota(itabla_z,1) ! si hay un receptor en medio estratificado
       nullify(PX)
       PX => allpoints(pota(itabla_z,itabla_x))
       !if (PX%region .eq. 1) &
       porLoMenosUnoEsEstr = .true. !'estr'
      end do
end function porLoMenosUnoEsEstr

SUBROUTINE HANKELS(Z,H0,H1)
!#define comparar
!#ifdef comparar
!     use specfun
!     use glovars,only:UI
!     integer :: NM,K
!     complex*16, dimension(0:10) :: CBJ,CDJ,CBY,CDY
!#endif
!     Z = COMPLEX ARGUMENT
!
!     COMPUTE SECOND KIND HANKEL FUNCTIONS H0 AND H1
!
      COMPLEX*16 :: Z,H0,H1,C,A,E,E2,ZH,P
      real*8 :: X,Y,R,PHI,J,AR
      X=REAL(Z)
      Y=AIMAG(Z)
      R=SQRT(X*X+Y*Y)
      PHI=ATAN2(Y,X)
      IF(R.LE.10.0)GO TO 20
      J=2.0*R
      C=(0.0,0.1250)/Z
      K=2
      P=C*C
      A=4.5*P
      P=7.5*P
      H0=1.0+C+A
      H1=1.0-3.0*C-P
10    I=4*K
      K=K+1
      DI=I
      DK=K
      A=A*C*(DI+1.0/DK)
      P=P*C*(DI-3.0/DK)
      H0=H0+A
      H1=H1-P
      AR=ABS(REAL(P))+ABS(AIMAG(P))
      IF(AR.GT.1.E-16.AND.K.LT.J)GO TO 10
      AR=0.785398163397448-X-PHI/2.0
      E=0.0
      IF(Y.GT.-160.0) E=0.7978845608028650/SQRT(R)*EXP(Y)*CMPLX(COS(AR),SIN(AR),8)
!     IF(X.EQ.0.0)E=CMPLX(0.0,AIMAG(E))
      IF(abs(X) .lt. 0.00001)E=CMPLX(0.0,AIMAG(E),8)
      H0=H0*E
      H1=H1*E*(0.0,1.0)
      GO TO 23
20    ZH=Z/2.0
      C=-ZH*ZH
      E=CMPLX(0.0,0.3183098861837910,8)
      E2=E*2.0
      A=1.0-E2*(0.5772156649015330+LOG(R/2.0))+PHI*0.636619772367582
      P=1.0
      K=1
      H0=A
      H1=A+E*(1.0-1.0/C)
25    A=A+E2/K
      P=P*C
      H0=H0+A*P
      K=K+1
      P=P/(K*K)
      H1=H1+(A*K+E)*P
      IF(ABS(REAL(P))+ABS(AIMAG(P)).GT.1.E-16) GO TO 25
      H1=H1*ZH
!     IF(X.NE.0.0)GO TO 23
      IF(abs(X) .gt. 0.00001) GO TO 23
      H0=CMPLX(0.0,AIMAG(H0),8)
      H1=CMPLX(REAL(H1),0.0,8)

23    K=K

!#ifdef comparar
!     print*,""
!     print*,"Z=",Z
!     print*,"H_0^(2)", H0
!     print*,"H_1^(2)", H1
!     call CJYNB(2,Z,NM,CBJ,CDJ,CBY,CDY)
!       WRITE(*,*)
!       WRITE(*,*)'   n     Re[Jn(z)]       Im[Jn(z)]',&
!                 '       Re[Yn(z)]      Im[Yn(z)]'
!       WRITE(*,*)' -------------------------------------',&
!                 '-------------------------------'
!       DO K=0,NM,1
!          WRITE(*,'(1X,I4,4D16.8)')K,CBJ(K),CBY(K)
!       end do
!       WRITE(*,*)
!       WRITE(*,*)'   n     Re[H_n^(2)(z)]       Im[H_n^(2)(z)]'
!       WRITE(*,*)' -------------------------------------',&
!                 '-------------------------------'
!       DO K=0,NM,1
!          WRITE(*,'(1X,I4,2D26.16)')K,CBJ(K)-UI*CBY(K)
!       end do
!     stop "HANK"
!#endif
      RETURN
      END SUBROUTINE HANKELS

end module matrizGlobal
