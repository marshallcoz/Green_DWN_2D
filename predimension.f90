subroutine predimensionar
use glovars, only: z0
use sourceVars, only : SH,PSV,nFuentes
use resultVars
use refSolMatrixVars
use modelVars
implicit none
integer :: i,ik,ip

      if (SH) then
        i = 2*N+1
      end if !
      if (PSV) then
        i = 4*N+2
      end if
      ik = 0 !indice para onda plana
        allocate (Ak(i,i,0:2*nmax)); allocate (B(i,0:2*nmax))
        allocate(ipivA(i)); allocate(workA((i)*(i)))
   !   to store the strata diffracted displacement and stress: W,U,V,...
      do iP=1,nPts
        allocate(allpoints(ip)%resp(NFREC+1,nFuentes))
        if (PSV) then
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%U = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%W = z0
!        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%Tx = z0
!        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%Tz = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%szz = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%szx = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%sxx = z0
        end if
        if (SH) then
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%V = z0
!        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%Ty = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%sxy = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%szy = z0
        end if
         ! Plot integrand FK of Green(l) function
        if(allpoints(iP)%guardarFK)then
         allocate(allpoints(iP)%FK(NFREC+1,NMAX,3)) !W,U,V
         allpoints(iP)%FK = z0
        end if
      end do
      do i=1,Npts
        allpoints(i)%pointIndex = i
      end do ! i
      allocate( gamma_E(0:2*nmax,N+1))
      allocate( nu_E(0:2*nmax,N+1))
      allocate( eta_E(0:2*nmax,N+1))
      allocate( subMatD0(2,4,N+1,0:2*nmax))
      allocate( subMatS0(3,4,N+1,0:2*nmax))
      allocate( k_vec(2*nmax))
      k_vec(1) = real(dk * 0.01,8)
      do ik = 2,NMAX+1
        k_vec(ik) = real(ik-1,8) * dk
      end do!
      do ik = nmax+2,2*NMAX
        k_vec(ik) = (ik - 2*NMAX - 1) * dk
      end do
end subroutine predimensionar
