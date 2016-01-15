subroutine preparePointerTable ! pota
      use resultVars, only : Punto,nPts,allpoints,nZs,pota
      use debugstuff
      use glovars, only : verbose,PrintNum

      ! tabla con las coordenadas Z (sin repetir).
      implicit none
      real*8 :: ratioWL_bola
      integer :: i,j,iP
      logical :: nearby,esnueva
!      type(Punto), pointer :: PX
      integer, allocatable, dimension(:,:) :: auxpota

!     !#< r
!     bola = min(0.01_8*smallestWL,0.01_8)
!     bola = 0.1*smallestWL !#>
      ratioWL_bola = 0.01
      ! si es la primera vez que corre sólo agregamos los allpoints
!      if (firstTime) then
       if (verbose .ge. 4) write(PrintNum,'(a)') " ...generating pointer table"
       allocate(auxpota(npts,npts+2))
       auxpota = 0
       ! siempre agregamos el primer punto [ allpoints(1)% ]
       nZs = 1
       auxpota(1,1) = 1 !nXs allpoints
       auxpota(1,3) = 1 !4,5,6 ... allpoints -7,-8,-9 ... boupoints
       !         ^--- el (3) siempre está
       ! agregamos los demás sin repetir
       if (nPts>=2) then
       do iP = 2,nPts
         esnueva = .true.
         do i = 1,nZs
           if (nearby(allpoints(iP)%center%z, &
               allpoints(auxpota(i,3))%center%z,0.01_8)) then
             !agregar a coordenada Z existente
             auxpota(i,1) = 1 + auxpota(i,1)
             auxpota(i,auxpota(i,1) + 2) = iP
             esnueva = .false.
             exit
           end if
         end do
         if (esnueva) then
           !nueva coordenada Z
           nZs = nZs + 1
           auxpota(nZs,1) = 1
           auxpota(nZs,3) = iP
         end if
       end do
       end if
       j = maxval(auxpota(:,1))

       allocate(PoTa(nZs,2+j))
       Pota = auxpota(1:nZs,1:2+j)
       deallocate(auxpota)

      if (verbose .ge. 4) call showMNmatrixI(nZs,2+j,pota,"po_ta",6)
      if (verbose .ge. 4) print*, "out of preparePointerTable"
      end subroutine preparePointerTable

      function nearby(a,b,bola)
      implicit none
      real*8, intent(in) :: a,b,bola
      logical :: nearby
      nearby = .false.
      if ((b-bola .le. a) .and. (a .le. b+bola)) then
        nearby = .true.
      end if
!     print*, b-bola, "<=", a, "<=",b+bola ," :: ", nearby
      end function nearby


