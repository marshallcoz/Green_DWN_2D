subroutine setupmodel
! Definir el modelo a partir de input.txt
use gloVars
use modelVars
use sourceVars
use fitting
use resultVars
implicit none
! variables auxiliares:
logical :: lexist
character(21) :: time,auxtx
CHARACTER(len=400) :: path
real :: H, ALF, BET, RO
integer :: J,i
real*8 :: BEALF,frac,minBeta,k1_3,kmax
integer :: thelayeris,ths_layer
logical :: tellisoninterface, ths_isoninterface
real    :: xfsource,zfsource,l
real*8  :: nxfsource,nzfsource, PW_theta
real :: Escala
real :: mastimere,Ts,Tp,T0
integer :: ampfunction
real :: sigGaus
integer :: PW_pol
integer :: tipofuente,nnsabanas,thisnsab,iIndex,nIpts,nSabanapts
real :: xini,deltax,zini,deltaz,dx,dz
real*8 :: escalax,escalay,offsetx,offsety



CALL getcwd(path)
write(PrintNum,'(a,/,a)') 'At:',TRIM(path)
call date_and_time(TIME=time); write(PrintNum,'(a,a)') "hhmmss.sss = ",time
      inquire(file="input.txt",exist=lexist)
      if (lexist) then
        OPEN(35,FILE="input.txt",FORM="FORMATTED")
      else
        write(PrintNum,'(a)') 'There is a missing input file. '
        stop 'Check "input.txt" on Working directory'
      end if
READ(35,*) !########## del modelo #################################################################
READ(35,*) DFREC
READ(35,*) NFREC
READ(35,*) NPTSTIME
         NPTSTIME = int(log10(real(NPTSTIME)) / log10(2.)+0.99 ) ! porqueFFT
         NPTSTIME = 2** NPTSTIME     !adjuste to a 2**N value
         if (NPTSTIME .lt. 2*NFREC) then
           NPTSTIME=2*NFREC
           print*,"WARNING, modified NPTSTIME to" , NPTSTIME
         end if
         Dt = (1.0) / (real(NPTSTIME) * DFREC)
READ(35,*) Qq
READ(35,'(L1)') useAzimi
READ(35,*) NK
READ(35,*) cKbeta
READ(35,*) TW
READ(35,*) periodicdamper
READ(35,*) !########### del medio estratificado ###################################################
READ(35,'(I2)')N   !NUMBER OF LAYERS
      ALLOCATE (Z(0:N+1)); Z(1)=real(0,8)
      ALLOCATE (RHO(0:N+2));     ALLOCATE(ANU(0:N+2))
      ALLOCATE (BETA0  (0:N+2)); ALLOCATE (ALFA0(0:N+2))
      ALLOCATE (BETA   (0:N+2)); ALLOCATE (ALFA (0:N+2))
      ALLOCATE (LAMBDA0(0:N+2)); ALLOCATE (AMU0 (0:N+2))
      ALLOCATE (LAMBDA (0:N+2)); ALLOCATE (AMU  (0:N+2))
write(PrintNum,'(a,I0)')'Number of layers= ',N
write(PrintNum,'(a)')'    zini     zfin      alpha    beta    mu       rho      lambda      nu '
READ(35,*) ! H(m)    Alfa(m/s) Beta(m/s) Dens(T/m3)
DO J=1,N
         READ(35,*) H, ALF, BET, RO
         Z(J+1)=Z(J)+real(H)
         AMU0(J)=RO*BET**2.0
         BETA0(J)=BET
         ALFA0(J)=ALF
         RHO(J) = RO
         LAMBDA0(J)=RHO(J)*ALF**2.0 - real(2.)*real(AMU0(J))
!        BEALF=SQRT((0.5-ANU)/(1.0-ANU)) !IF POISSON RATIO IS GIVEN
         BEALF = beta0(J)/alfa0(J)
         anu(J) = (bealf**2 - 0.5)/(-1 + bealf**2)
!        ALFA(J)=BET/BEALF
          write(PrintNum,&
          '(F7.1,A,F7.1,2x, F7.1,2x, F7.1,2x, E8.2,2x, E8.2,2x, E8.2,2x, E8.2)') &
          Z(J),' - ',Z(J+1),ALFA0(J),BETA0(J),real(AMU0(J)),&
          RHO(J), real(LAMBDA0(J)),real(anu(J))
      END DO
      READ(35,*) H, ALF, BET, RO
      AMU0(N+1)=RO*BET**2
      BETA0(N+1)=BET
      ALFA0(N+1)=ALF
      RHO(N+1) = RO
      LAMBDA0(N+1)=RHO(n+1)*ALF**2 - real(2.)*real(AMU0(n+1))
      BEALF = beta0(N+1)/alfa0(N+1)
      anu(N+1) = (bealf**2 - 0.5)/(-1 + bealf**2)
         write(PrintNum,&
         '(F7.1,A,           F7.1,2x, F7.1,2x, E8.2,2x, E8.2,2x, E8.2,2x, E8.2)') &
         Z(1),' -     inf. ',ALFA0(N+1),BETA0(N+1),real(AMU0(N+1)),&
         RHO(N+1), real(LAMBDA0(N+1)),real(anu(N+1))

!#< r ----  dk para ver estar CKbeta veces por arriba de polo Rayleigh------- !#>
      frac = 1./3.
      minBeta = minval(beta0(1:N+1));
      kmax = 0.9* (2*pi*DFREC*NFREC) / minBeta * cKbeta !1.5
      DK = kmax / nk
      k1_3 = (2*pi*DFREC*NFREC * frac) / minBeta * cKbeta

      allocate(vecNk(NFREC+1))
      ! cantidad de numeros de onda {donde se invierte la matriz} en cada frecuencia
      ! el resto hasta nmax+1(o hasta SpliK) se interpola
      vecNk = (/(i,i=1, NFREC+1)/)
      call splineIn(vecNk, & ! vector x [int]
                    vecNk, & ! vector interpolado [int]
                    int(nfrec * frac),& ! n puntos [x0 ... x1-1]
                    NFREC+1-int(nfrec * frac), & ! [x1 ...  x2 ]
                    1,                  int(0.55*nK), &  !x0,y0
                    int(nfrec * frac)+1,int(0.6*nk), &  !x1,y1
                    NFREC+1,            nk) !x2,y2
      NMAX = vecNK(NFREC)
      NMAX = int(log10(real(NMAX)) / log10(2.)+0.99 )
      NMAX = 2**NMAX     !adjuste to a 2**N value

      ! complex frecuency to "smooth" the poles and be able to integrate
        ! Bouchon (2003) OMEI entre -pi/T y -2pi/T ; T= 2pi/DFREC tiempo total
        ! tengo TW la ventana de tiempo de inter�s.
      OMEI = - periodicdamper*PI/TW ! se cumple que exp(omei TW) << 1

       write(PrintNum,'(/,a)') &
       "--- frecuency --------------------------------------------------------------------"
       write(PrintNum,'(a,F9.7,/,a,F15.5)') "   dt = (1.0) / (real(NPTSTIME) * DFREC) = ",Dt,"   Tmax=",Dt* NPTSTIME
       write(PrintNum,'(a,F8.1)') '   Atenuation Q = ',Qq
       write(PrintNum,'(a,I0,a,F8.4,a,F12.4,a,/)') &
           '   N. frequencies: ',NFREC,'  @',DFREC,'Hertz :. Fmax = ', &
           NFREC*DFREC,'Hertz'

       write(PrintNum,'(a)') &
       "--- DWN -------------------------------------------------------------------------"
       write(PrintNum,'(a,F15.7)') "   kmax = 0.9* (2*pi*DFREC*NFREC) / minBeta * cKbeta = ",kmax
       write(PrintNum,'(a,F9.7)') "   DK = kmax / nk = ",DK
       write(PrintNum,'(a,I0,a,I0,a,I0,a,I0,a,I0,a,I0,a)') "   nk a 3pt spline with (",&
          1,",",int(0.55*nK),"), (",&
          int(nfrec * frac)+1,",",int(0.6*nk),"), (",&
          NFREC+1,",",nk,")"
       write(PrintNum,'(a,I0)') '   nk average = ',int(sum(vecNK)/(NFREC+1))
       write(PrintNum,'(a,I0)') '   nmax (each sign): ',NMAX
       write(PrintNum,'(a,F12.7)') "   delta X = ", real(pi / (NMAX * DK),4)
       write(PrintNum,'(a)') "   Periodicidad de la fuente (metodo dwn):"
       write(PrintNum,'(a,EN14.2E2,a)') "   L = 2*pi/DK = ",2*pi/DK, " m"
       write(PrintNum,'(a,EN19.5E2,a)') &
       "   L/maxalfa = ",(2*pi/DK)/maxval(abs(ALFA0(1:N+1))), " seconds"
       write(PrintNum,'(a,EN19.5E2,a)') &
       "   L/maxbeta = ",(2*pi/DK)/maxval(abs(BETA0(1:N+1))), " seconds"
       write(PrintNum,'(a,EN19.5E2,a)') &
       "   L/minalfa = ",(2*pi/DK)/minval(abs(ALFA0(1:N+1))), " seconds"
       write(PrintNum,'(a,EN19.5E2,a)') &
       "   L/minbeta = ",(2*pi/DK)/minval(abs(BETA0(1:N+1))), " seconds"
       write(PrintNum,'(a,E10.2,a)') '   Frec. Imaginary part: - periodicdamper*PI/TW = ',OMEI,' rad/s'
       write(PrintNum,'(/,a)') &
       "--- SOURCE ----------------------------------------------------------------------"
READ(35,*) !########### de la incidencia ##########################################################
      READ(35,'(I1)') i
      SH = .false.; PSV = .true.; skipdir=.true.
       if (i .eq. 2) then ! SH
            SH = .true.; PSV = .false.
            skipdir(1) = .true.
            skipdir(2) = .false.
            skipdir(3) = .true.
       end if

      READ(35,'(I1)') nFuentes
      READ(35,*) !  X    Z  | nx   nz  |  th    l    Esca  Afn  T0   Tmax   Ts    Tp    Gsg   PWp  Knd
      allocate(Po(nFuentes))
do i=1,nFuentes
       READ(35,*) xfsource, zfsource, &
                  nxfsource, nzfsource, &
                  PW_theta, l,&
                  Escala, ampfunction,&
                  t0,mastimere, Ts, Tp, sigGaus,&
                  PW_pol
       Po(i)%center%x = xfsource
       Po(i)%center%z = zfsource
       Po(i)%normal%x = nxfsource
       Po(i)%normal%z = nzfsource
       Po(i)%cosT = cos((360-PW_theta)*pi/180.0)
       Po(i)%sinT = sin((360-PW_theta)*pi/180.0)
       Po(i)%gamma = PW_theta*pi/180.0 !clockwise desde el eje z (hacia abajo)

       Po(i)%Escala=Escala
       Po(i)%ampfunction=ampfunction
       Po(i)%maxtime = mastimere
       Po(i)%Ts=Ts
       Po(i)%Tp=Tp
       Po(i)%t0=T0;
       Po(i)%sigGaus=sigGaus
       Po(i)%PW_pol=PW_pol
     if (PSV .eqv. .true.) then
       if (PW_pol .eq. 0)then
          !fuerza puntual P_SV,
          if (abs(nxfsource) > 0.01) skipdir(1) = .false.
          if (abs(nzfsource) > 0.01) skipdir(3) = .false.
          if (abs(nxfsource)+abs(nzfsource) < 0.01) then
          stop "normales nulas"
          end if
          tipoFuente = 0
       elseif (PW_pol .eq. 1) then
          !SV
          skipdir(1) = .false.
          tipoFuente = 1
       elseif (PW_pol .eq. 2) then
          !P
          skipdir(1) = .false.
          tipoFuente = 1
       else
          stop "PW_pol incorrecto"
       end if
     else !SH
       if (PW_pol .eq. 3) then
          !SH
          tipoFuente = 1
        elseif (PW_pol .eq. 0) then
          !fuerza puntual SH
          tipoFuente = 0
        else
           stop "PW_pol incorrecto"
        end if
     end if
       Po(i)%tipoFuente=tipoFuente
       if (tipoFuente .eq. 1) then ! onda plana
            ! se coloca en la interfaz con el semiespacio
            Po(i)%center%z = Z(N+1)
            ths_layer = N+1
            ths_isoninterface = .true.
       else ! fuente puntual
           ths_layer = thelayeris(real(zfsource,8))
           ths_isoninterface = tellisoninterface(real(zfsource,8))
       end if
       Po(i)%layer = ths_layer
       Po(i)%isOnInterface = ths_isoninterface
       if (tipoFuente .eq. 1) then ! onda plana
        write(auxtx,'(a)') " plane wave"
      else
        write(auxtx,'(a)') " point force"
      end if
      write(PrintNum,'(a,F8.2,a,F8.2,a,2x,a,F9.2,a,F9.2,a,I0,a,L1,a)') &
      "   (",Po(i)%center%x,",",Po(i)%center%z,")", &
      "n=[",Po(i)%normal%x,",",Po(i)%normal%z,&
      "] e=",Po(i)%layer," intf=",Po(i)%isOnInterface,auxtx
end do
       print*,"   skipdir=",skipdir
       print*,"   PSV= ",PSV
       print*,"   SH=  ",SH

      write(PrintNum,'(/,a)') &
       "--- RECEIVERS -------------------------------------------------------------------"
      READ(35,*) !########### de los receptores #########################################################
      READ(35,*) nIpts
      READ(35,*) nnsabanas, nSabanapts
      READ(35,*) !_____________________SINGLE RECEIVERS__________________________
      READ(35,*) escalax,escalay
      READ(35,*) offsetx,offsety
      READ(35,*) ! X        Z          nx       nz   guardarFK
      allocate(allPoints(nIpts + nSabanapts))
      nPts = nIpts + nSabanapts
      do i=1, nIpts
         READ(35,*) allPoints(i)%center%x, allPoints(i)%center%z, allPoints(i)%guardarFK
           allPoints(i)%center%x = allPoints(i)%center%x * escalax + offsetx
           allPoints(i)%center%z = allPoints(i)%center%z * escalay + offsety
      !encontrar el layer en el que estan o 0 si est� sobre la interfaz
           allPoints(i)%layer = thelayeris(real(allPoints(i)%center%z,8))
           allPoints(i)%isOnInterface = &
                         tellisoninterface(real(allPoints(i)%center%z,8))
      end do
      iIndex = nIpts
      READ(35,*) !___________________RECEIVER LINES_____________________________
      if (nSabanapts .gt. 0) then
      read(35,*) escalax,escalay
      READ(35,*) offsetx,offsety
      read(35,*) !npuntos xini   deltax   zini   delta z  guardarFK
      do j=1,nnsabanas
      read(35,*) thisnsab,xini,deltax,zini,deltaz
      dx = 0.0
      dz = 0.0
      do i=1,thisnsab
        iIndex = iIndex + 1
        allPoints(iIndex)%isSabana = .true.
        allPoints(iIndex)%center%x = (xini + dx)*escalax + offsetx
        allPoints(iIndex)%center%z = (zini + dz)*escalay + offsety
        ths_layer = thelayeris(allPoints(iIndex)%center%z)
        ths_isoninterface = tellisoninterface(allPoints(iIndex)%center%z)
        allPoints(iIndex)%layer = ths_layer
        allPoints(iIndex)%isOnInterface = ths_isoninterface
        dx = dx + deltax
        dz = dz + deltaz
        allPoints(iIndex)%guardarFK = .false.
      end do ! i
      end do ! j
      else
      READ(35,*) !escala
      READ(35,*) !offsetof
      READ(35,*) !encabezado
      end if
      CLOSE(35)
      do i=1,nPts
      write(PrintNum,'(I0,2x,a,f7.2,a,f7.2,a,a,i0,a,l)'),&
        i,"[",allpoints(i)%center%x,",",allpoints(i)%center%z,"]",&
        " l",allpoints(i)%layer,&
        " OnIntfce=",allpoints(i)%isOnInterface
      end do
end subroutine setupmodel

      function thelayeris(zi)
      use modelVars, only : Z,N
      implicit none
      integer :: i,e,thelayeris
      real*8 :: zi
      i = 1
      if (Z(0) .lt. 0.0) i = 0
      do e=i,N
         if(real(zi) .lt. Z(e+1)) then
            exit
         end if
      end do
      thelayeris = e
      end function thelayeris


      function tellisoninterface(zi)
      use modelVars, only : Z,N
      implicit none
      real*8 ::  errT = 0.0001 !0.01_8
      integer :: e
      real*8 :: zi
      logical :: tellisoninterface
      tellisoninterface = .false.
      do e=1,N+1
        if((Z(e)-errT .lt. real(zi)) &
             .and. (real(zi) .lt. Z(e)+errT)) then
           tellisoninterface = .true.
        end if
      end do
      end function tellisoninterface
