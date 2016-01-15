      module fitting
      contains

      function splitatY(surf_poly,degree,Y,aIN,bIN) !SINGLE PRES, REAL
      ! there is only one intersection.
      implicit none
!     integer, parameter           :: dp = selected_real_kind(15, 307)
      real, dimension(:),intent(in) :: surf_poly
      real, intent(in) :: Y,aIN,bIN
      real, allocatable, dimension(:) :: surf_poly0
      integer, intent(in) :: degree
      integer :: i
      real :: a,b,splitatY,af,bf,temp,tmp2,tmp, &
                 cx,d,cf,errorTol,s,sf
      logical :: mflag = .true.
      errorTol = real(0.000001,8)

      allocate(surf_poly0(degree+1))
      a = aIN
      b = bIN

      surf_poly0 = surf_poly
      surf_poly0(1) = surf_poly0(1) - Y

      !encontramos el cero de surf_poly0 entre a y b
      af = polyVal(surf_poly0,degree,a)
      bf = polyVal(surf_poly0,degree,b)
      ! usamos el mŽtodo de Brent.
      if (af*bf >= real(0,8)) then
        if (af < bf ) then
          splitatY = af
        else
          splitatY = bf
        end if
        return
      else

        if (abs(af) < abs(bf)) then
          temp = b
          b = a
          a = temp
          temp = bf
          bf = af
          af = temp
        end if
        cx = a
        cf = af
        mflag = .true.
        i = 0
        d = 0.
        do while( (abs(bf) > errortol) .and. (abs(a-b) > errorTol ) )
!       do while((.not.(bf == real(0,8) )) .and. (abs(a-b) > errorTol ))
!          print*,"go_",i
          if ((abs(af-cf) > errorTol) .and. (abs(bf-cf)>errorTol)) then
!          if ((af /= cf) .and. (bf /= cf)) then
          ! interpolaci—n cuadr‡tica inversa
            s = a * bf * cf / (af-bf) / (af-cf) + &
            b*af*cf/(bf-af)/(bf-cf)+cx*af*bf/(cf-af)/(cf-bf)
          else
          ! regla de la secante
            s = b - bf * (b-a)/(bf-af)
          end if
          tmp2 = (3.0*a + b)/4.0
          if ( (.not.(((s > tmp2) .and. (s < b)) .or. &
          ((s < tmp2) .and. (s > b))) ) .or. &
          (mflag .and. ((abs(s-b)) .ge. (abs(b-cx)/2.0 ))) .or. &
          ((.not. (mflag)) .and. ((abs(s-b)) .ge. (abs(cx-d)/2.0 ))) ) then
            s = (a+b) / 2.0
            mflag = .true.
          else
            if ((mflag .and. (abs(b-cx)< errorTol)) .or. &
            ((.not. (mflag)) .and. (abs(cx-d) < errorTol))) then
              s = (a+b) / 2.0
              mflag = .true.
            else
              mflag = .false.
            end if
          end if
           sf = polyVal(surf_poly0,degree,s)
           d = cx
           cx = b
           cf = bf
!          if (af * sf < real(0,8)) then
           if (af * sf < errorTol) then
             b = s
             bf = sf
           else
             a = s
             af = sf
           end if
           if (abs(af) < abs(bf)) then
             tmp = a
             a = b
             b = tmp
             tmp = af
             af = bf
             bf = tmp
           end if
           i = i + 1
           if (i> 1000) then
             !error
             b = real(0.123456789)
             exit
           end if
        end do
      splitatY = b
      end if
      end function splitatY


      ! function evaluation
      function polyVal(surf_poly,degree,X) ! SINGLE, REAL
      implicit none
!     integer, parameter           :: dp = selected_real_kind(15, 307)
      real, dimension(:),intent(in) :: surf_poly
      integer, intent(in) :: degree
      real, intent(in) :: X
      real :: polyVal
      integer :: i

      !surfo_poly are the polynomial coefficients: A0 A1 A2...An
      polyVal = surf_poly(1) !A0
      DO i = 1,degree
      polyVal = polyVal + X**i * surf_poly(i+1)
      end do

      end function

      !fit a polynomial Anx^n + ... + A2x^2 + A1x + A0 = 0
      !of degree 'd' to data: (vx,vy)
      ! coefficients are orthered from the lower power: A0 A1 A2 .. AN
      ! http://rosettacode.org/wiki/Polynomial_regression#Fortran
      function polyfit(vx,vy,Ln,d,verbose,outpf) ! DOUBLE, REAL
      implicit none
      integer, intent(in)               :: verbose,outpf
      integer, intent(in)               :: Ln, d
      real*8, dimension(d+1)              :: polyfit
      real*8, dimension(:), intent(in)    :: vx, vy
      real*8, dimension(:,:), allocatable :: X
      real*8, dimension(:,:), allocatable :: XT
      real*8, dimension(:,:), allocatable :: XTX
      integer :: i, j
      integer :: n, lda, lwork
      integer :: info
      integer, dimension(:), allocatable :: ipiv
      real*8, dimension(:), allocatable :: work

      n = d+1
      lda = n
      lwork = n

      allocate(ipiv(n))
      allocate(work(lwork))
      allocate(XT(n, Ln))
      allocate(X(Ln, n))
      allocate(XTX(n, n))

      if (verbose >= 1) then
        write(outpf,*)"fitting curve.."
      end if !
      if (verbose >= 4) then
       write(outpf,'(a)')"begin curve fit with:"
       write(outpf,'(a,I4)')"vx: of size:",size(vx)
       write (outpf, '(F9.4)') vx
       write(outpf,'(a,I4)')"vy: of size:",size(vy)
       write (outpf, '(F9.4)') vy
      end if

      ! prepare the matrix
      do i = 0, d
       do j = 1, size(vx)
          X(j, i+1) = vx(j)**i
       end do
      end do
      if (verbose >= 4) then
       write(outpf,'(a,I4,a,I4,a)') &
                              "X: size, (",size(X,1),',',size(X,2),')'
       write(outpf,*) X
      end if

      XT  = transpose(X)
      XTX = matmul(XT, X)
      if (verbose >= 4) then
       write(outpf,'(a,I4,a,I4,a)') &
                         "XTX: size, (",size(XTX,1),',',size(XTX,2),')'
       write(outpf,*)XTX
      end if

      ! calls to LAPACK subs DGETRF and DGETRI
      ! factorizacion LU
      call DGETRF(n, n, XTX, &
                  lda, ipiv, info)
      if (verbose >= 4) then
       write(outpf,'(a)')"DGETRF (XTX):"
       write (outpf, '(E12.3)') XTX
      end if
      !
      if ( info /= 0 ) then
       write(outpf,*) "problem DGETRF =",info
       stop 1
      end if
      ! inversa de la matriz
      call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
      if (verbose >= 4) then
       write(outpf,'(a)')"DGETRI (XTX):"
       write (outpf, '(F9.4)') XTX
      end if

      if ( info /= 0 ) then
       write(outpf,'(a)') "problem DGETRI =",info
       stop 1
      end if

      polyfit = matmul( matmul(XTX, XT), vy)
      deallocate(ipiv)
      deallocate(work)
      deallocate(X)
      deallocate(XT)
      deallocate(XTX)
      if (verbose >= 4) then
       write(outpf,'(a)')'polyfit='
       write (outpf, '(E12.3)') polyfit
      end if!
      if (verbose .ge. 4) then
      !fit a polynomial Anx^n + ... + A2x^2 + A1x + A0 = 0
      !of degree 'd' to data: (vx,vy)
      ! coefficients are orthered from the lower power: A0 A1 A2 .. AN
       write(outpf,'(a)', ADVANCE = "NO") "0 ="
       do i=1,d+1
          write(outpf,'(ES12.3,a,I0)', ADVANCE = "NO") polyfit(i),"x^",i-1
       end do
      end if
      end function polyfit

      function Zpolyfit(vx,vy,Ln,d) ! complex*16, dimension(d+1)
      use glovars, only : verbose, outpf => PrintNum
      implicit none
      real*8,     dimension(:), intent(in)    :: vx
      complex*16, dimension(:), intent(in)    :: vy
      integer, intent(in)                     :: d, Ln
      complex*16, dimension(d+1)              :: zpolyfit

      complex*16, dimension(:,:), allocatable :: X
      complex*16, dimension(:,:), allocatable :: XT
      complex*16, dimension(:,:), allocatable :: XTX
      integer :: i, j
      integer :: n, lda, lwork
      integer :: info
      integer, dimension(:), allocatable :: ipiv
      complex*16, dimension(:), allocatable :: work

      n = d+1
      lda = n
      lwork = n

      allocate(ipiv(n))
      allocate(work(lwork))
      allocate(XT(n, Ln))
      allocate(X(Ln, n))
      allocate(XTX(n, n))

      if (verbose >= 4) then
        write(outpf,*)"fitting curve.."
      end if !

      ! prepare the matrix
      do i = 0, d
       do j = 1, Ln
          X(j, i+1) = vx(j)**i
       end do
      end do
      XT  = transpose(X)
      XTX = matmul(XT, X)

      ! calls to LAPACK subs DGETRF and DGETRI
      ! factorizacion LU
      call ZGETRF(n, n, XTX, lda, ipiv, info)
      if ( info /= 0 ) then

       print*,Ln
       print*,vx
       write(6,*) "Zpolyfit :problem ZGETRF =",info
       stop
      end if
      ! inversa de la matriz
      call ZGETRI(n, XTX, lda, ipiv, work, lwork, info)
      if ( info /= 0 ) then
       write(6,'(a)') "Zpolyfit :problem ZGETRI =",info
       stop
      end if

      zpolyfit = matmul( matmul(XTX, XT), vy(1:Ln))
      deallocate(ipiv)
      deallocate(work)
      deallocate(X)
      deallocate(XT)
      deallocate(XTX)
      if (verbose .ge. 4) then
      !fit a polynomial Anx^n + ... + A2x^2 + A1x + A0 = 0
      !of degree 'd' to data: (vx,vy)
      ! coefficients are orthered from the lower power: A0 A1 A2 .. AN
       write(outpf,'(a)', ADVANCE = "NO") "0 ="
       do i=1,d+1
          write(outpf,'(a,ES12.3,a,ES12.3,a,I0)', ADVANCE = "NO") &
          'polyfit= ',real(zpolyfit(i))," i",aimag(zpolyfit(i)),"x^",i-1
       end do
      end if
      end function Zpolyfit

      function Zevapol(po,d,x)
      implicit none
      integer, intent(in)                     :: d
      complex*16, dimension(d+1)              :: po
      real*8 :: x
      complex*16 :: Zevapol
      integer :: i
      Zevapol = po(1)
      do i =2,d+1
        Zevapol = Zevapol + po(i)* x**(i-1)
      end do
      end function Zevapol

      !normal vectors to surf_poly at the nXI points (x,y)
      function normalto(x,y,nXI,surf_poly,degree,verbose,outpf)
      implicit none
      integer,intent(in)      :: verbose,outpf
      integer, intent(in)     :: degree
!     integer, parameter      :: dp = selected_real_kind(15, 307)
      real, intent(in), dimension(degree+1) :: surf_poly !surface
      integer :: nXI !number of points
      real, intent(in), dimension(nXI) :: x,y !points coordinates

      integer :: i,j
      real*8, allocatable :: normalto(:,:) !function result value
      real*8, parameter      :: a = real(1,8) ! normal poing up parameter
!     real*8, parameter      :: tolerance = real(0.001,8) !to tell if its vertical
      real*8, dimension(degree) :: fprime !surf_poly derivative coeficients
      real*8 :: fprimeX, x0, mag

      ALLOCATE (normalto(nXI,2))
      if(verbose >= 1)then
       write(outpf,'(a)',advance='no') "  ...getting the normal vectors"
      end if
      ! the normal line to curve f(x) = An x^n + ... + A1 x + A0
      ! is y = y1 - 1 / (f'(x1)) (x - x1)
      ! if we want the normal vectors pointing up (bigger 'y')
      ! we force   y = y1 + a   and we find x in terms of x1

      ! f'(x) coefficients
      do i = degree,1,-1
!      write(outpf,*)'i=',i
       fprime(i)= real(i,8) * surf_poly(i+1)
      end do
      if (verbose >= 4 ) then
       write(outpf,'(a)')'surf_poly=A0,A1,...,AN '
       write(outpf,'(F9.4)') surf_poly
       write(outpf,'(a)')'f_prime_(x) A0 A1 ... An '
       write(outpf,'(ES9.4/)') fprime
      end if

      do i = 1,nXI !for every point
      !the derivative at (xi)
       fprimeX = real(0,8)
       do j = size(fprime),2,-1
         fprimeX = fprimeX + fprime(j) * x(i)**(j-1)
       end do
       fprimeX = fprimeX + fprime(1)
       x0 = x(i) - a * fprimeX

       !normalizar y poner como vector
       mag = sqrt((x(i)-x0)**2 + a**2)

       normalto(i,1)= (x0-x(i))/mag
       normalto(i,2)= a/mag

       if (verbose >= 4) then
         write(outpf,'(A,f7.2,A,f7.2,A,/A,f6.1,A,f6.1,/A,f6.2,A,f6.2)') &
         "(",x(i),",",y(i),")","x0= ",x0,"  y0= ",y(i)+a, &
         "nx=",normalto(i,1),"  ny=",normalto(i,2)
       end if



      end do

      if(verbose >= 1)then
       write(outpf,'(a)',advance='yes') " ... done"
      end if

      ! if by incrementing  y = y1 + a the value becomes infinite, then
      ! it is a vertical surface.                                     TODO


      end function normalto


!     subroutine makeTaperFuncs_cutF_cutK!(dir)
!     use resultvars, only: Hf!,Hk
!     use waveNumVars, only : NFREC,DFREC!,NMAX,DK
!     use dislin
!     implicit none
!     integer, intent(in) :: dir
!     integer, parameter :: NpolF = 19
!     integer, parameter :: NpolK = 20
!     real, parameter :: cutoff_fracFmax = 0.92
!     real, parameter :: cutoff_fracKmax = 0.8
!     integer :: i
!
!     if(.not. allocated(Hf)) allocate(Hf(NFREC+1))
!     if(.not. allocated(Hk)) allocate(Hk(NMAX*2))
!
!        Hf = (/ (i,i=0,nfrec) /) * DFREC
!        Hf = 1.0/sqrt(1.0+(Hf/(cutoff_fracFmax*NFREC*DFREC))**(2*NpolF))


!!        Hf = 1.0_8
!
!        Hk(1:nmax+1) = (/ (cmplx(i,0,8),i=0,nmax) /) * DK
!        Hk(nmax+2:2*nmax) = (/ (cmplx(i,0,8),i=nmax-1,1,-1) /) * DK
!        HK = 1.0/sqrt(1.0+(Hk/(cutoff_fracKmax*NMAX*DK))**(2*NpolK))
!
!        if (dir .eq. 2) then
!          Hk = 1.0_8/1.135_8  !SH
!        else
!!          Hk = 1.0_8/3.385_8 ! dk = 0.00343750
! !         Hk = 1.0_8/2.3_8    !PSV ! 256  Kmax = 14.07 nK = 2048 dk = 0.006875
!  !        Hk = 1.0_8/1.155_8  !PSV ! 512  Kmax = 28.15 nK = 2048 dk = 0.01375
!   !       Hk = 1.0_8/0.595_8  !PSV ! 1024 Kmax = 56.29 nK = 2048 dk = 0.0275
!          HK = HK * 1.0_8/(6100.275482093669_8 * dk * dk &
!                      - 0292.363636363636 * dk &
!                      + 0004.021666666667)
!        end if

!!        CALL SETFIL('tapperK.pdf')
! !       call qplot((/(1.0*i,i=1,nmax*2)/),real(Hk,4),2*nmax)
!  !      stop
!         CALL SETFIL('tapperF.pdf')
!         call qplot((/(1.0*i,i=1,NFREC+1)/),real(Hf,4),NFREC+1)
!     !   stop
!     end subroutine makeTaperFuncs_cutF_cutK

      subroutine spline3p(x,y,n01,n12,x0,y0,x1,y1,x2,y2)
      ! interpolaci—n con spline de 3 puntos
      !    0     1         2    :  nodos
      !    1 2 3 4 5 6 7 8 9    :  vector interpolado
      !    <---> <--------->
      !     n01      n12
      implicit none
      integer, intent(in) :: n01,n12 !puntos entre nodos,incluidos los 3 nodos
      real*8, intent(inout), dimension(n01+n12) :: x
      complex*16, intent(inout), dimension(n01+n12) :: y
      complex*16, intent(in) :: y0,y1,y2
      real*8, intent(in) :: x0,x1,x2
      integer :: i
      complex*16, dimension(3,3) :: A
      complex*16, dimension(3) :: B
      complex*16 :: a1,b1,a2,b2,t
      integer, dimension(3) :: ipiv
      complex*16, dimension(3) :: work
      integer :: info
      integer :: lwork
      A = 0
      A(1,1) = 2 / (x1-x0)
      A(1,2) = 1 / (x1-x0)
      A(2,1) = 1 / (x1-x0)
      A(2,2) = 2*(1/(x1-x0) + 1/(x2-x1))
      A(2,3) = 1 / (x2-x1)
      A(3,2) = 1 / (x2-x1)
      A(3,3) = 2 / (x2-x1)
      B(1) = 3*(y1-y0) / (x1-x0)**2
      B(2) = 3*((y1-y0)/(x1-x0)**2 + (y2-y1)/(x2-x1)**2)
      B(3) = 3*(y2-y1)/ (x2-x1)**2
      ! resolver y obtener ko, k1, k2
      lwork = 3*3
      call zgetrf(3,3,A,3,ipiv,info)
      call zgetri(3,A,3,ipiv,work,lwork,info)
      if(info .ne. 0) stop "Problem at inverse of SPLINE3P matrix "
      B = matmul(A,B)!     print*,B
      a1 = B(1)*(x1-x0) - (y1-y0)!; print*,a1
      b1 = -B(2)*(x1-x0) + (y1-y0)!; print*,b1
      a2 = B(2)*(x2-x1) - (y2-y1)!; print*,a2
      b2 = -B(3)*(x2-x1) + (y2-y1)!; print*,b2
      y(1) = y0
      do i = 2,n01
      t = (x(i)-x0) / (x1-x0)
      ! en puntos intermedios a x0 y x1
      y(i) = (1-t)*y0 + t * y1 + t * (1- t)* (a1*(1-t)+b1*t)
      end do!
      y(n01+1) = y1
      do i = n01+2,n01+n12-1
      t = (x(i)-x1) / (x2-x1)
      ! en puntos intermedios a x1 y x2
      y(i) = (1-t)*y1 + t * y2 + t * (1- t)* (a2*(1-t)+b2*t)
      end do
      y(n01+n12) = y2
      end subroutine spline3p

      subroutine splineIn(x,y,n01,n12,x0,y0,x1,y1,x2,y2)
!     use debugstuff
      ! interpolaci—n con spline de 3 puntos
      !    0     1         2    :  nodos
      !    1 2 3 4 5 6 7 8 9    :  vector interpolado
      !    <---> <--------->
      !     n01      n12
      implicit none
      integer, intent(in) :: n01,n12 !puntos entre nodos,incluidos los 3 nodos
      integer, intent(inout), dimension(n01+n12) :: x
      integer, intent(inout), dimension(n01+n12) :: y
      integer, intent(in) :: y0,y1,y2
      integer, intent(in) :: x0,x1,x2
      integer :: i
      real, dimension(3,3) :: A
      real, dimension(3) :: B
      real :: a1,b1,a2,b2,t
      integer, dimension(3) :: ipiv
      real, dimension(3) :: work
      integer :: info
      integer :: lwork
!     print*,"spline3pIn"
!     print*,n01,n12
!     print*,x0,y0
!     print*,x1,y1
!     print*,x2,y2
      A(1,1) = 2. / (x1-x0)
      A(1,2) = 1. / (x1-x0)
      A(1,3) = 0.
      A(2,1) = 1. / (x1-x0)
      A(2,2) = 2.*(1./(x1-x0) + 1./(x2-x1))
      A(2,3) = 1. / (x2-x1)
      A(3,1) = 0.
      A(3,2) = 1. / (x2-x1)
      A(3,3) = 2. / (x2-x1)
      B(1) = 3.*(y1-y0) / (x1-x0)**2.
      B(2) = 3.*(1.*(y1-y0)/(x1-x0)**2. + 1.*(y2-y1)/(x2-x1)**2.)
      B(3) = 3.*(y2-y1)/ (x2-x1)**2.

      ! resolver y obtener ko, k1, k2
      lwork = 3*3
!     call showMNmatrixR(3,3,A,"  A  ",6)
!     call showMNmatrixR(3,1,B,"  B  ",6)
      call sgetrf(3,3,A,3,ipiv,info)
!     if(info .ne. 0) stop "Problem at LU of SPLINE3P matrix "
      call sgetri(3,A,3,ipiv,work,lwork,info)
!     if(info .ne. 0) stop "Problem at inverse of SPLINE3P matrix "
      B = matmul(A,B)!     print*,B
      a1 = B(1)*(x1-x0) - 1.*(y1-y0)!; print*,a1
      b1 = -B(2)*(x1-x0) + 1.*(y1-y0)!; print*,b1
      a2 = B(2)*(x2-x1) - 1.*(y2-y1)!; print*,a2
      b2 = -B(3)*(x2-x1) + 1.*(y2-y1)!; print*,b2
      y(1) = y0
      do i = 2,n01
      t = 1.*(x(i)-x0) / (x1-x0)
      ! en puntos intermedios a x0 y x1
      y(i) = int((1.-t)*y0 + t * 1.*y1 + t * (1.- t)* (a1*(1.-t)+1.*b1*t))
      end do!
      y(n01+1) = y1
      do i = n01+2,n01+n12-1
      t = 1.*(x(i)-x1) / (x2-x1)
      ! en puntos intermedios a x1 y x2
      y(i) = int((1.-t)*y1 + t * 1.*y2 + t * (1.- t)* (a2*(1.-t)+1.*b2*t))
      end do
      y(n01+n12) = y2
      end subroutine splineIn

      end module fitting
