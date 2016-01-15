module gloVars
      use, intrinsic :: iso_c_binding
      complex*16, parameter :: UI = cmplx(0.0d0,1.0d0,8), &
                               UR = cmplx(1.0d0,0.0d0,8), &
                               Z0 = cmplx(0.0d0,0.0d0,8)
      real*8, target :: ONE = 1.0_8
      real*8, parameter :: PI = real(4.0d0*ATAN(1.0d0),8)
      CHARACTER(len=400) :: rutaOut
      type(C_PTR),save  :: planNmaxF,planNmaxB,&
                           planNfrecF,planNfrecB,&
                           planNtimeF,planNtimeB
      integer,parameter :: verbose = 3
      integer,parameter :: PrintNum = 6
      logical,parameter :: PWfrecReal=.true.
end module gloVars

module resultVars
       type FFres
         complex*16 :: U,V,W,Tx,Ty,Tz,szz,szx,sxx,sxy,szy
       end type FFres

       type Punto2d
        real*8 :: x,z
       end type Punto2d

       type MecaElem
        type (Punto2d) :: center
        complex*16, dimension(5) :: Rw !u1,u3,s33,s31,s11
        complex*16, dimension(3) :: Rw_SH !u2,s32,s12
       end type MecaElem

       type Punto
        ! para el receptor
        type(Punto2d) :: center
        type(Punto2d) :: normal
        integer :: layer
        logical :: isOnInterface
        logical :: guardarFK
        logical :: isSabana
        integer :: pointIndex

        ! para la fuente
        real*8  :: gamma,cosT,sinT
        real :: Escala,Ts,Tp,t0,maxtime,sigGaus
        integer :: ampfunction ! 0 dirac; 1 ricker; 2 gaussian; 3 inAmplitude.txt
        integer :: PW_pol,tipofuente
      !                       ,--- f: 1...nfrec+1
      !                       | ,--- k: 1...NMAX+1 / 2*NMAX
      !                       | | ,--- iMec: 1:3 (desps) W,U,V
      !                       | | |
        complex*16, dimension(:,:,:), allocatable :: FK
      !
        complex*16, dimension(:,:,:)   , allocatable :: G

        ! espectro campo total inquirePoints :
      !                        ,--- f: 1...nfrec+1
      !                        | ,--- iMec: 1:2 y 3
        type(FFres),dimension (:,:), allocatable :: resp
!        complex*16, dimension (:,:), allocatable :: S !(traza,componente)
       end type Punto

      ! bondary elements:     ,--- POINT index / x (receptor)
      type (Punto), dimension(:), allocatable, save, target :: allpoints
        integer, save :: nPts

        integer, allocatable, dimension(:,:), save :: pota
        integer :: nZs ! depths at pota
end module resultVars

module sourceVars
        use resultVars, only : Punto
        type (Punto),dimension(:),allocatable, save, target :: Po
        logical, save :: SH,PSV
        integer :: nFuentes
        integer :: currentiFte
        logical :: skipdir(3)
end module sourceVars

module modelVars
      !frequency loop vars:
      integer,save      :: NFREC,NPTSTIME
      integer,dimension (:),allocatable :: vecNK
      complex*16,target :: cOME
      real*8   ,save    :: FREC,DFREC,OME,OMEI,TW
      real*8, parameter :: SpliK=1.25
      real*8 :: Qq
      logical    :: useAzimi
      !Discrete Wave-number:
      real*8   ,save    :: DK,cKbeta,periodicdamper    ! delta k
      real*8, dimension(:), allocatable,target :: k_vec
      integer,save      :: NK,NMAX
      !Soil
      integer ::  N !number of layers from z>= 0. HALF-SPACE at N+1
      real*8,     dimension(:), allocatable :: ALFA0,BETA0,AMU0,LAMBDA0,ANU,Z,RHO
      complex*16 ,dimension(:), allocatable :: ALFA ,BETA ,AMU, LAMBDA
      complex*16, dimension(:,:), allocatable :: gamma_E,nu_E,eta_E
      !Output
      real*8, save :: Dt  !segundos
      complex*16, dimension(:), allocatable :: t_vec
      complex*16, dimension(:), allocatable :: Uo
end module modelVars


      module refSolMatrixVars
        complex*16, save, allocatable :: B(:,:)
        complex*16, dimension(:,:,:), allocatable,target :: Ak
        integer, dimension(:), allocatable :: IPIV

        complex*16, dimension(:,:,:,:), allocatable,target :: BparaGa,BparaNu
        complex*16, dimension(:,:,:,:,:), allocatable,target :: CoefparGa, CoefparNu
        !                     | | | | '-- interfaz de arriba o abajo [2]
        !                     | | | '-- direcci�n de la fuerza [2]
        !                     | | `-- k [2nmax]
        !                     | '-- estrato [N+1]
        !                     '-- onda desde la interfaz [4N+2]

        complex*16, dimension(:,:,:,:),allocatable :: subMatD0,subMatS0
        !                         | '-- ik
        !                         '- estrato
        !                      ' '-- renglon,columna
        complex*16, pointer :: pt_cOME_i
        real*8, pointer :: pt_k
        integer, dimension(:),allocatable,target :: ipivA
        integer, dimension(:),pointer :: pt_ipivA
        complex*16, dimension(:),allocatable,target :: workA
        complex*16, dimension(:),pointer :: pt_workA
        complex*16,dimension(:,:),pointer :: pointA
 !     type(Punto), pointer :: pXi,p_X
      end module refSolMatrixVars

module wavelets
      contains

      !  The Ricker wavelet on the time domain saved on   Uo
      subroutine ricker(Uo,Ts,Tp)
      use gloVars, only : PI,UR
      use modelVars, only : Dt,NPTSTIME
      implicit none
      integer :: i
      complex*16, dimension(NPTSTIME) :: Uo
      real :: Ts,Tp
      real*8 :: A,B
      Uo = 0;
      A = nearest(pi*(-Ts) / Tp,1.0)
      A = nearest(real((A*A-0.5)* exp(- A * A), 8),1.0)
      if (abs(A) .lt. 0.0001) A = 0
      Uo(1) = A * UR

      do i = 2, NPTSTIME/2+1
      ! NEAREST(X, S) returns the processor-representable number
      ! nearest to X in the direction indicated by the sign of S
        A = nearest(pi*(Dt*(i*1.0_8-1.0_8)-Ts) / Tp,1.0)
        A = nearest(A * A,1.0)
        B = nearest(exp(-A),1.0)
        A = nearest((A - 0.5_8) * B ,1.0)
        if (abs(A) .lt. 0.001) A = 0
        Uo(i) = A * UR
      end do
      end subroutine ricker

      subroutine gaussian(Uo,NPTS,sigGausIn)
!     use waveVars, only : Uo,sigGaus
!     use waveNumVars, only : NPTSTIME!,nfrec
      implicit none
      integer :: i,NPTS
      real*8 :: f,s
      complex*16, dimension(NPTS) :: Uo
      real :: sigGaus,sigGausIn
      integer :: shift
      f = 0.0
      sigGaus = int(sigGausIn)
      if (sigGausIn  - sigGaus .gt. 0.001) then
        shift = int((sigGausIn  - sigGaus) * 1000)
      else
        shift = 0
      end if
      !positivos
      s = real(sigGaus/100.0 * NPTS/2,8)
      do i=1, NPTS/2+1
        f = real(i-1,8) ! Hz
        Uo(i) = cmplx(exp(-0.5*(f/s)**2.),0.,8)
      end do
      !negativos
      if (shift .eq. 0) then
      Uo(NPTS/2+2:NPTS) = conjg(Uo(NPTS/2:2:-1))
      else
      if (shift + NPTS/2+1 .gt. NPTS) shift = NPTS/2
      Uo(shift:shift + NPTS/2+1) = Uo(1:NPTS/2+1)
      Uo(1:shift) = 1.0
      Uo(NPTS/2+2:NPTS) = conjg(Uo(NPTS/2:2:-1))
      end if
      end subroutine gaussian

      function FFTW(n,Uin,direccion,escala)
      use gloVars,only:planNmaxF,planNmaxB,&
                           planNfrecF,planNfrecB,&
                           planNtimeF,planNtimeB
      use modelVars, only: nfrec,nmax,NPTSTIME
      use, intrinsic :: iso_c_binding
      include 'fftw3.f03'
      integer, intent(in) :: n,direccion
      complex*16 :: Uin(n)
      complex*16 :: FFTW(n)
      real*8,intent(in) :: escala
      type(C_PTR) :: plan
      ! El plano creado en checarWisdom a partir del wisdom se reusa
      if (direccion .lt. 1) then ! forward -1
      if (n .eq. 2*NFREC) plan = planNfrecF
      if (n .eq. 2*NMAX) plan = planNmaxF
      if (n .eq. NPTSTIME) plan = planNtimeF
      else !backward +1
      if (n .eq. 2*NFREC) plan = planNfrecB
      if (n .eq. 2*NMAX) plan = planNmaxB
      if (n .eq. NPTSTIME) plan = planNtimeB
      end if
      ! ejecutar con plan de reuso
      if (C_ASSOCIATED (plan) .eqv. .false.) stop "fail FFTW_WISDOM_ONLY"
      call fftw_execute_dft(plan,Uin,FFTW)
      ! hacia adelante
      !Uof = Uof * Dt
      ! hacia atr�s
      !Uot = Uot / (NPTSTIME*dt) ! *DFREC
      FFTW = FFTW * escala
      end function FFTW
end module wavelets


module debugStuff
      contains

      subroutine showMNmatrixZ(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n,outpf
      complex*16, dimension(m,n), intent(in) :: MAT
      integer :: i,j
      character(LEN=5), intent(in) :: name

      write(outpf,'(A)') trim(name)
      if (outpf .eq. 6) then
      do i = 1,m
        do j = 1,n
          write(outpf,'(A,E9.2,A)',advance='no') "(",REAL(MAT(i,j)),","
          write(outpf,'(E9.2,A)',advance='no') AIMAG(MAT(i,j)),"i) "
        end do
        write(outpf,'(A)',advance='yes')''
      end do
      else !a un archivo
      do i = 1,m
        do j = 1,n
          write(outpf,'(EN26.9,1X)',advance='no') REAL(MAT(i,j))
          write(outpf,'(EN26.9,3X)',advance='no') AIMAG(MAT(i,j))
        end do
        write(outpf,'(A)',advance='yes')''
      end do
      end if
      end subroutine
      !
      subroutine scripToMatlabMNmatrixZ(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n, outpf
      complex*16, dimension(m,n), intent(in) :: MAT
      integer :: i,j
      character(LEN=32), intent(in) :: name

      write(outpf,'(A,A)') trim(name), " = ["
      do i = 1,m
        write(outpf,'(a)',advance='no') "["
        do j = 1,n
          write(outpf,'(a,EN26.9,a)',advance='no') "(",REAL(MAT(i,j)),")+("
          write(outpf,'(EN26.9,a,2X)',advance='no') AIMAG(MAT(i,j)),")*1i"
        end do
        write(outpf,'(A)',advance='yes')'];'
      end do
      write(outpf,'(a)') "];"
      end subroutine scripToMatlabMNmatrixZ

      subroutine showMNmatrixI(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n, outpf
      integer, dimension(m,n), intent(in) :: MAT
      integer :: i,j
      character(LEN=5), intent(in) :: name

      write(outpf,'(A)') ""
      write(outpf,'(A)') trim(name)
      do i = 1,m
        do j = 1,n
          write(outpf,'(I0,3x)',advance='no') MAT(i,j)
        end do
        write(outpf,'(A)',advance='yes')''
!       write(outpf,'(A)',advance='yes')''
      end do
      write(outpf,'(A)') ""
      end subroutine showMNmatrixI

end module debugStuff

