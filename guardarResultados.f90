subroutine guardarResultados(saveASmatlab,saveAstxt,sl)
use debugstuff
use resultVars, only : allpoints,nPts
use sourcevars, only: PSV,SH,currentiFte,nfuentes,Po
use modelVars, only: NFREC,NPTSTIME,dfrec,t_vec,Uo
use glovars, only: z0,pi
implicit none
integer :: ip,i,comp
character(LEN=100) :: yax
logical, intent(in) :: saveASmatlab,saveAstxt,sl
real :: t0
complex*16, dimension(NPTSTIME) :: S
complex*16, dimension(:,:,:), allocatable :: ShR
character(LEN=32)  :: name
character(LEN=100) :: titleN
character(LEN=3),dimension(8) :: cname

!    W_to_t  guarda la respuesta en frecuencia y tiempo en los formatos indicados
   allocate(t_vec(NPTSTIME))
   allocate(Uo(NPTSTIME))
   call system("mkdir out")
   call chdir("out")
   CALL getcwd(yax)
      write(6,'(a,/,a)') 'Saving files At:',TRIM(yax)
   call system('rm *.m *.txt')

   if (sl) allocate(ShR(nPts,8,NPTSTIME))

   do currentiFte=1,nFuentes
      call sourceAmplitudeFunction
      if (sl) ShR = 0; cname(1:8) = ''
      ! para la correción del tiempo inicial de la señal
      t0 = Po(currentiFte)%t0
      t_vec(1:NPTSTIME) = z0
      t_vec(1) = exp(cmplx(0.0,-0.01* dfrec * t0*(2*pi),8))
      do i = 2,nfrec+1
        t_vec(i) = exp(cmplx(0.0,-(i-1)* dfrec * t0*(2*pi),8))
      end do!
      do i = NPTSTIME-NFREC+2,NPTSTIME
        t_vec(i) = exp(cmplx(0.0,-(i-NPTSTIME-1)* dfrec * t0*(2*pi),8))
      end do

      if (PSV) then
        do iP = 1,nPts
        !#< r ____________   W   _________________________________________ !#>
          write(yAx,'(a)') '$u_3$ [m]'
!         print*,ip,"---------------------",allpoints(iP)%resp(:,currentiFte)%W; print*," "
          call W_to_t(S,allpoints(iP)%resp(:,currentiFte)%W,'w--',yax,iP,&
                    allpoints(iP)%center%x,&
                    allpoints(iP)%center%z,currentiFte,saveASmatlab,saveAstxt)
          if (sl) shR(iP,3,:) = S; cname(3) = '_w_'
        !#< r ____________   U   _________________________________________ !#>
          write(yAx,'(a)') '$u_1$ [m]'
          call W_to_t(S,allpoints(iP)%resp(:,currentiFte)%U,'u--',yax,iP,&
                    allpoints(iP)%center%x,&
                    allpoints(iP)%center%z,currentiFte,saveASmatlab,saveAstxt)
          if (sl) shR(iP,1,:) = S; cname(1) = '_u_'
        !#< r ____________   sxx   _________________________________________ !#>
          write(yAx,'(a)') '$\sigma_{xx}$ [m]'
          call W_to_t(S,allpoints(iP)%resp(:,currentiFte)%sxx,'sxx',yax,iP,&
                    allpoints(iP)%center%x,&
                    allpoints(iP)%center%z,currentiFte,saveASmatlab,saveAstxt)
          if (sl) shR(iP,4,:) = S; cname(4) = 'sxx'
        !#< r ____________   szx   _________________________________________ !#>
          write(yAx,'(a)') '$\sigma_{zx}$ [m]'
          call W_to_t(S,allpoints(iP)%resp(:,currentiFte)%szx,'szx',yax,iP,&
                    allpoints(iP)%center%x,&
                    allpoints(iP)%center%z,currentiFte,saveASmatlab,saveAstxt)
          if (sl) shR(iP,5,:) = S; cname(5) = 'szx'
        !#< r ____________   szz   _________________________________________ !#>
          write(yAx,'(a)') '$\sigma_{zz}$ [m]'
          call W_to_t(S,allpoints(iP)%resp(:,currentiFte)%szz,'szz',yax,iP,&
                    allpoints(iP)%center%x,&
                    allpoints(iP)%center%z,currentiFte,saveASmatlab,saveAstxt)
          if (sl) shR(iP,6,:) = S; cname(6) = 'szz'
        end do
      end if !psv
      if (SH) then
        do iP = 1,nPts
          !#< r ____________   V   _________________________________________ !#>
          write(yAx,'(a)') '$u_2$ [m]'
          call W_to_t(S,allpoints(iP)%resp(:,currentiFte)%V,'v--',yax,iP,&
                    allpoints(iP)%center%x,&
                    allpoints(iP)%center%z,currentiFte,saveASmatlab,saveAstxt)
          if (sl) shR(iP,2,:) = S; cname(2) = '_v_'
          write(yAx,'(a)') '$\sigma_{xy}$ [m]'
          call W_to_t(S,allpoints(iP)%resp(:,currentiFte)%sxy,'sxy',yax,iP,&
                    allpoints(iP)%center%x,&
                    allpoints(iP)%center%z,currentiFte,saveASmatlab,saveAstxt)
          if (sl) shR(iP,7,:) = S; cname(7) = 'sxy'
          write(yAx,'(a)') '$\sigma_{zy}$ [m]'
          call W_to_t(S,allpoints(iP)%resp(:,currentiFte)%szy,'szy',yax,iP,&
                    allpoints(iP)%center%x,&
                    allpoints(iP)%center%z,currentiFte,saveASmatlab,saveAstxt)
          if (sl) shR(iP,8,:) = S; cname(8) = 'szy'
!          call W_to_t(allpoints(iP)%resp(:,currentiFte)%Ty,'Ty-',yax,iP,&
!                    allpoints(iP)%center%x,&
!                    allpoints(iP)%center%z,saveASmatlab,saveAstxt)
        end do !ip
      end if !sh

      if (sl) then
      do comp=1,8
       if (cname(comp) .eq. '') cycle
       write(titleN,'(a,I0,a,a,a)') 'In_',currentiFte,'_',cname(comp),'.m'
       OPEN(3555,FILE=trim(titleN),FORM="FORMATTED",ACTION='WRITE')
       !write(name,'(a,I0,a,I0)') 'inci_',currentiFte,'_',comp
        write(name,'(a)') 'var'
       call scripToMatlabMNmatrixZ(nPts,NPTSTIME,shR(:,comp,:),name,3555)
       close(3555)
      end do
      end if
    end do !insidencia
end subroutine guardarResultados

subroutine W_to_t(Sout,W,nombre,yax,iP,x_i,z_i,currentiFte,saveASmatlab,saveAstxt)
      use modelVars, only : NFREC, NPTSTIME, OMEI, &
                            t_vec, dt,Uo
      use glovars, only : z0
      use resultvars, only : allpoints
      use sourceVars, only : Po
      use wavelets
      use debugStuff
      implicit none
      complex*16, dimension(Nfrec+1), intent(in) :: W !espectro
      character(LEN=3)   :: nombre ! del componente
      character(LEN=100) :: titleN,yAx, CTIT ! texto para grafica
      integer ,intent(in) :: iP ! indice del receptor
      real*8, intent(in) :: x_i,z_i ! coordenadas del receptor
      integer, intent(in) :: currentiFte ! indice de la incidencia actual
      logical, intent(in) :: saveASmatlab,saveAstxt

      character(LEN=32)  :: name
      complex*16, dimension(NPTSTIME) :: S,Sout
      integer :: i,n_maxtime

      write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x_i,' , ',z_i,')'

      S = z0
      S(1:nfrec+1)= W(1:nfrec+1:+1)
      S(NPTSTIME-NFREC+2:NPTSTIME) = conjg(W(nfrec:2:-1))

      ! (0) espectro sin corregir
      if (saveAstxt .eqv. .true.) then
        write(titleN,'(I0,a,I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') &
               iP,'_in_',currentiFte,'_rawfrec_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', &
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].txt'
        write(name,'(a,I0,a,a,I0)') 'rawf_IP',iP,trim(nombre),'_in_',currentiFte
        OPEN(3201,FILE=trim(titleN),FORM="FORMATTED",ACTION='WRITE')
        write(3201,'(a,a,/,a,EN26.9,a,EN26.9,a,/,a)') '%',yAx, &
         '%(', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')','%% '
        call showMNmatrixZ(int(nfrec +1),1,S(1:int(nfrec +1)),name,3201)
        close (3201)
       end if
!        call plotXYcomp(S(1:int(nfrec)),real(DFREC,4), &
!         int(nfrec),titleN, 'frec[hz] ',yAx, CTIT ,1200,800,0.0)
       if (saveASmatlab .eqv. .true.) then
        write(titleN,'(I0,a,I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') &
               iP,'_in_',currentiFte,'_rawfrec_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', &
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].m'
        write(name,'(a,I0,a,a,I0)') 'rawf_IP',iP,trim(nombre),'_in_',currentiFte
        OPEN(3202,FILE=trim(titleN),FORM="FORMATTED",ACTION='WRITE')
        write(3202,'(a,a,/,a,EN26.9,a,EN26.9,a,/,a)') '%',yAx, &
         '%(', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')','%% '
        call scripToMatlabMNmatrixZ(int(nfrec +1),1,S(1:int(nfrec +1)),name,3202)
        close (3202)
       end if

      !(1) señal en tiempo


!     if(ip.eq.3) then
!     print*,"*******"
!     DO i=1,NPTSTIME
!     print*,i,S(i),t_vec(i)
!     end do
!     end if


      S = S * t_vec ! t0 tiempo inicial
      S = S * Uo ! conv con fucion de amplitud
      S = FFTW(NPTSTIME,S,+1,1/(NPTSTIME*dt)) !backward (pasar al tiempo)

      ! (2) remover efecto de la frecuencia imaginaria (método DWN)
      !  La frecuencia imaginaria se utiliza como un amortiguamiento artificial
      !  que tiene el efecto de atenuar la amplitud de los arrivos de las
      !  fuentes periódicas en el espacio
      !  (surgen de truncar y discretizar una integral infinita)
      S = S * exp(- OMEI * Dt*((/(i,i=0, NPTSTIME-1)/))) !

      ! tiempo maximo para graficar
         n_maxtime = int(Po(currentiFte)%maxtime/dt)
         if(Po(currentiFte)%maxtime .lt. dt) n_maxtime = 2*nfrec
         if(Po(currentiFte)%maxtime .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
         S(n_maxtime+1: NPTSTIME) = z0;
         Sout = S

      !  (3) GRAFICAR en el tiempo
      if (saveAstxt .eqv. .true.) then
        write(titleN,'(I0,a,I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') &
               iP,'_in_',currentiFte,'_sism_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', &
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].txt'
        write(name,'(a,I0,a,a,I0)') 'Sism_IP',iP,trim(nombre),'_in_',currentiFte
        OPEN(3203,FILE=trim(titleN),FORM="FORMATTED",ACTION='WRITE')
        write(3203,'(a,a,/,a,EN26.9,a,EN26.9,a,/,a)') '%',yAx, &
         '%(', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')','%% '
        call showMNmatrixZ(NPTSTIME,1,S(1:int(NPTSTIME)),name,3203)
        close (3203)
       end if

       if (saveASmatlab .eqv. .true.) then
        write(titleN,'(I0,a,I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') &
               iP,'_in_',currentiFte,'_sism_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', &
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].m'
        write(name,'(a,I0,a,a,I0)') 'Sism_IP',iP,trim(nombre),'_in_',currentiFte
        OPEN(3204,FILE=trim(titleN),FORM="FORMATTED",ACTION='WRITE')
        write(3204,'(a,a,/,a,EN26.9,a,EN26.9,a,/,a)') '%',yAx, &
         '%(', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')','%% '
        call scripToMatlabMNmatrixZ(NPTSTIME,1,S(1:int(NPTSTIME)),name,3204)
        close (3204)
       end if


      !  (4) espectro correcto
!       ! tapper ?

        ! fft
        S = FFTW(NPTSTIME,S(1:int(NPTSTIME)),-1,dt) !forward

        if (saveAstxt .eqv. .true.) then
        write(titleN,'(I0,a,I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') &
               iP,'_in_',currentiFte,'_f_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', &
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].txt'
        write(name,'(a,I0,a,a,I0)') 'f_IP',iP,trim(nombre),'_in_',currentiFte
        OPEN(3205,FILE=trim(titleN),FORM="FORMATTED",ACTION='WRITE')
        write(3205,'(a,a,/,a,EN26.9,a,EN26.9,a,/,a)') '%',yAx, &
         '%(', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')','%% '
        call showMNmatrixZ(int(NPTSTIME/2+1),1,S(1:int(NPTSTIME/2+1)),name,3205)
        close (3205)
       end if

       if (saveASmatlab .eqv. .true.) then
        write(titleN,'(I0,a,I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') &
               iP,'_in_',currentiFte,'_f_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', &
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].m'
        write(name,'(a,I0,a,a,I0)') 'f_IP',iP,trim(nombre),'_in_',currentiFte
        OPEN(3206,FILE=trim(titleN),FORM="FORMATTED",ACTION='WRITE')
        write(3206,'(a,a,/,a,EN26.9,a,EN26.9,a,/,a)') '%',yAx, &
         '%(', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')','%% '
        call scripToMatlabMNmatrixZ(int(NPTSTIME/2+1),1,S(1:int(NPTSTIME/2+1)),name,3206)
        close (3206)
       end if

!     ! grafica logaritmica
!           write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') &
!              'fL_',nombre,iP,'[', &
!              int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', &
!              int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
!           logflag = 'logx     '
!           write(xAx,'(a)') "frec[Hz]"
!           call plotSpectrum(S(1: NPTSTIME/2),real(DFREC,4), NPTSTIME/2, NPTSTIME/2, &
!               titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))

end subroutine W_to_t

subroutine sourceAmplitudeFunction
      use wavelets !las funciones: ricke
      use modelVars, only : Dt,Uo,dt!,maxtime
      use modelVars, only : nfrec, NPTSTIME,OMEI!DFREC
      use gloVars, only : ve => verbose,Ur,Ui,z0,Printnum
      use sourceVars, only:Po,iFte=>currentiFte
      !use ploteo10pesos
      implicit none
      integer  :: i,nval
!      character(LEN=9)          :: logflag
!      character(LEN=100)        :: titleN,xAx,yAx,CTIT

!      CHARACTER(len=32) :: arg
      logical :: inquire, lexist!, argumA
      real*8 :: val1,val2
      complex*16, dimension(NPTSTIME) :: S
      !complex*16 :: FFTW!(NPTSTIME)
      integer :: n_maxtime
!      argumA = .false.
!      CALL get_command_argument(1, arg)
!      IF (LEN_TRIM(arg) .ne. 0) then
!        if (trim(arg) .eq. '-a') then
!        iFte = 1
!        ve = 6
!        argumA = .true.
!        end if
!      end if

       n_maxtime = int(Po(iFte)%maxtime/dt)
       if(Po(iFte)%maxtime .lt. dt) n_maxtime = 2*nfrec
       if(Po(iFte)%maxtime .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME

      ! grafica del factor de correción
      S = exp(- OMEI * Dt*((/(i,i=0, NPTSTIME-1)/)))
!      write(titleN,'(a)') 'x-LevantonDWN.pdf'
!      write(CTIT,'(a)') 'exp(-omei t)'
!      write(xAx,'(a)') 't[sec]'
!      write(yAx,'(a)') ''
!        call plotXYcomp(S(1:n_maxtime),real(Dt,4), &
!         int(n_maxtime),titleN,xAx,yAx, CTIT ,1200,800,0.0)


!     real :: factor
      ! Amplitude function of incident wave
      ! prepare the signal we will use as incident wave amplitude.
      ! el tamaño del ricker es 2*NFREC porque incluirá el conjugado

! 153
      Uo(:)=z0

      if(Po(iFte)%ampfunction .eq. 1) then
        call ricker(Uo(:),Po(iFte)%Ts,Po(iFte)%Tp) ! Ricker wavelet saved on Uo
!          if (ve .ge. 1) then
!            write(Printnum,'(a)')'   Incident wave amplitude function: Ricker'
!            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-t.pdf'
!            write(CTIT,'(a)') 'WaveAmplitude of Ricker wavelet'
!            xAx = 'time[sec]'
!            write(yAx,'(a)') 'amplitude'
!            call plotXYcomp(Uo(1:n_maxtime),real(Dt,4), &
!                 n_maxtime,titleN,xAx,yAx,CTIT,1200,800,0.0)
!          end if
        ! forward
        Uo(:) = FFTW(NPTSTIME,Uo(:),-1,Dt)
        Uo(1) = 0
!          if (ve .ge. 1) then
!            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-f.pdf'
!            xAx = 'frec[Hz] '
!            write(yAx,'(a)') 'amplitude'
!!           logflag = 'logx     '
!!           logflag = 'none     '
!            call plotXYcomp(Uo(:),real(DFREC,4), &
!                 size(Uo(:)),titleN,xAx,yAx,CTIT,1200,800,0.0)
!!           call plotSpectrum(Uo(:,iFte),real(DFREC,4), size(Uo(:,iFte)),int(size(Uo(:,iFte))/2), &
!!           titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
!          end if
     !********************************************************************************
      elseif(Po(iFte)%ampfunction .eq. 2) then ! Gaussian
        call gaussian(Uo(:),NPTSTIME,Po(iFte)%sigGaus)
!          if (ve .ge. 1) then
!           write(Printnum,'(a)')'   Incident wave amplitude function: Gaussian'
!            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-f.pdf'
!            write(CTIT,'(a)') 'WaveAmplitude of Gaussian wavelet'
!            write(xAx,'(a)') 'frec[Hz] '
!            write(yAx,'(a)') 'amplitude'
!!           call plotXYcomp(Uo(:,iFte),real(DFREC,4), &
!!                n_maxtime,titleN,xAx,yAx,CTIT,1200,800,0.0)
!!           print*,Uo(:,iFte)
!            call plotXYcomp(Uo(:),real(DFREC,4), &
!                 size(Uo(:)),titleN,xAx,yAx,CTIT,1200,800,0.0)
!!           call plotSpectrum(Uo(:,iFte),real(DFREC,4), size(Uo(:,iFte)),int(size(Uo(:,iFte))/2), &
!!           titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
!          end if
      elseif(Po(iFte)%ampfunction .gt. 20) then ! GaussianCúbico
        call gaussian(Uo(:),NPTSTIME,Po(iFte)%sigGaus) !shift =  int((int(sigGaus) - sigGaus) * 1000)
        if (int((int(Po(iFte)%sigGaus) - Po(iFte)%sigGaus) * 1000) .gt. 0) &
        write(Printnum,'(a,I0)')'   One sided Gaussian filter with shift in points:',&
         int((int(Po(iFte)%sigGaus) - Po(iFte)%sigGaus) * 1000)
        Uo(:) = Uo(:) ** (Po(iFte)%ampfunction-20)
!          if (ve .ge. 1) then
!           write(Printnum,'(a,I0)')'   Incident wave amplitude function: Gaussian ',int((Po(iFte)%ampfunction-20))
!            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-f.pdf'
!            write(CTIT,'(a,I0,a)') 'WaveAmplitude of Gaussian',int((Po(iFte)%ampfunction-20)),' wavelet'
!            write(xAx,'(a)') 'frec[Hz] '
!            write(yAx,'(a)') 'amplitude'
!!           call plotXYcomp(Uo(:,iFte),real(DFREC,4), &
!!                n_maxtime,titleN,xAx,yAx,CTIT,1200,800,0.0)
!!           print*,Uo(:,iFte)
!            call plotXYcomp(Uo(:),real(DFREC,4), &
!                 size(Uo(:)),titleN,xAx,yAx,CTIT,1200,800,0.0)
!!           call plotSpectrum(Uo(:,iFte),real(DFREC,4), size(Uo(:,iFte)),int(size(Uo(:,iFte))/2), &
!!           titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
!          end if
      elseif(Po(iFte)%ampfunction .eq. 3) then ! inDispl.txt
     !********************************************************************************
!        CALL chdir("..")
!        CALL chdir("ins")
        inquire(file="inAmplitude.txt",exist=lexist)
        if (lexist) then
        OPEN(77,FILE="inAmplitude.txt",FORM="FORMATTED")
        else
        write(6,'(a)') 'There is a missing input file. '
        stop 'Check "inAmplitude.txt" on this directory'
        end if
        READ(77,*)
        READ(77,*) inquire
        READ(77,*) nval
        if (nval .gt. NPTSTIME) then
        print*, "**** warning inAmplitude trimmed to ", NPTSTIME
        nval = NPTSTIME
        end if
        Uo(:) = 0
        if ( inquire) then
        do i = 1,nval
        READ(77,*) val1
        Uo(i) = UR * val1
        end do
        close(77)
!          if (ve .ge. 1) then
!            CALL chdir("..")
!            call chdir(trim(adjustl(rutaOut)))
!            write(Printnum,'(a)')'   Incident wave amplitude function from file'
!            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-t.pdf'
!            write(CTIT,'(a)') 'WaveAmplitude'
!            xAx = 'time[sec]'
!            write(yAx,'(a)') 'amplitude'
!            call plotXYcomp(Uo(1:NPTSTIME),real(Dt,4), &
!                 n_maxtime,titleN,xAx,yAx,CTIT,1200,800,0.0)
!          end if

        ! forward
        Uo(:) = FFTW(NPTSTIME,Uo(:),-1,Dt)

!          if (ve .ge. 1) then
!            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-f.pdf'
!            xAx = 'frec[Hz] '
!            write(yAx,'(a)') 'amplitude'
!            logflag = 'logx     '
!!           logflag = 'none     '
!            call plotSpectrum(Uo(:),real(DFREC,4), size(Uo(:)),int(size(Uo(:))/2), &
!            titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
!            CALL chdir("..")
!            CALL chdir("ins")
!          end if
        else
        do i = 1,nval
        READ(77,*) val1, val2
        Uo(i) = UR * val1 + UI * val2
        end do
!        if (ve .ge. 1) then
!            CALL chdir("..")
!            call chdir(trim(adjustl(rutaOut)))
!            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-f.pdf'
!            xAx = 'frec[Hz] '
!            write(yAx,'(a)') 'amplitude'
!            logflag = 'logx     '
!!           logflag = 'none     '
!            call plotSpectrum(Uo(:),real(DFREC,4), size(Uo(:)),int(size(Uo(:))/2), &
!            titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
!            CALL chdir("..")
!            CALL chdir("ins")
!          end if
        close(77)
        end if

!        CALL chdir("..")
!        call chdir(trim(adjustl(rutaOut)))

      else  ! DIRAC -----------------------------------------
       if (ve .ge. 1) write(Printnum,'(a)')'   Incident wave amplitude function: Dirac delta'
        Uo(:)=UR
      end if !tipoPulso
      Uo(:) = Uo(:)*Po(iFte)%Escala !escala de la señal
!      write(Printnum,'(a)') &
!       "---------------------------------------------------------------------------------"
!      if (argumA .eqv. .true.) then
!      if (iFte .eq. nfuentes) then
!      stop "argumento -a"
!      else
!      iFte = iFte + 1
!      go to 153
!      end if
!      end if
end subroutine sourceAmplitudeFunction
