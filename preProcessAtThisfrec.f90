subroutine preProcessAtThisfrec(J)
use glovars
use modelVars
use refSolMatrixVars
use matrizGlobal
use sourceVars, only: PSV,SH
implicit none
integer, intent(in) :: J
integer :: i,l,ik,tam,k0
complex*16, dimension(:,:), pointer :: pointAp

        FREC=DFREC*real(J-1); if (J .eq. 1)  FREC = 0.5_8 * DFREC  ! Hz
        OME=2.0*PI*FREC !rad/s
        COME = CMPLX(OME, OMEI,8)!periodic sources damping
        COME = COME * cmplx(1.0, -1.0/2.0/Qq,8) !histeretic damping
        pt_come_i => COME
        ik = N+1; i = 1
        if (useAzimi) then
        ! Azimi attenuation: see Kennett, The seismic wavefield, vol1 pag 148-150
          do l=i,ik !N+1 estratos y una inclusión
          ! "Normally we expect that loss in dilatation is very small compared 
          ! with that in shear so that Qp^-1 << Qs^-1 , and then "
!         Qs(l) = Qq !* 4./3. * (beta0(l)/alfa0(l))**2. !; print*,"Qs ",Qs(l)
!         Qp(l) = Qq !; print*,"Qp ",Qp(l)
          ! alfa0 y beta0 son velocidades de referencia a T=1sec
           beta(l) = cmplx((beta0(l))*(1.+1./pi/Qq*log(ome/2./pi)),-1./2./Qq,8)
           alfa(l) = cmplx((alfa0(l))*(1.+1./pi/Qq*log(ome/2./pi)),-1./2./Qq,8)
           aMU(l) = RHO(l) * beta(l)**2.
           Lambda(l) = RHO(l)*alfa(l)**2. - real(2.)*aMU(l)
          end do
        else
          do l=i,ik !N+1 estratos y una inclusión
           beta(l) = beta0(l); alfa(l) = alfa0(l)
           aMU(l) = RHO(l) * beta(l)**2.
           Lambda(l) = RHO(l)*alfa(l)**2. - real(2.)*aMU(l)
          end do
        end if

         ! apuntadores a la matriz global
         pt_ipivA => ipivA
         pt_workA => workA
         call makeGANU (J) !los numeros de onda horizontales para P y S

         if (PSV) then
         tam = size(Ak,1)
       do ik = 1,vecNK(J) ! k positivo (Aorig -> nmax+1)
         pointAp => Ak(1:tam,1:tam,ik)
         pt_k => k_vec(ik)
         call gloMat_PSV(pointAp,pt_k,ik)
         call inverseA(pointAp,pt_ipivA,pt_workA,tam)
       end do ! ik
       k0 = vecNK(J);
       call intrplr_gloMat(k0,15,pt_cOME_i,pt_ipivA,pt_workA)
!      do ik = 1,min(int(vecNK(J)*SpliK),nmax)+1
!        print*,ik,sum(Ak(1:tam,1:tam,ik))
!      end do
!      print*,"sum=",sum(Ak(1:tam,1:tam,1:min(int(vecNK(J)*SpliK),nmax)+1))
!      stop
       call parImpar_gloMat ! k negativo

       ! Funciones de Green para propagar a la frontera de cada estrato
       ! (sin fase vertical, la fase vertical se agrega para cada fuente)
       if(.not. allocated(BparaGa)) allocate(BparaGa(tam,N+1,2*nmax,2));BparaGa=0
       if(.not. allocated(BparaNu)) allocate(BparaNu(tam,N+1,2*nmax,2));BparaNu=0
       call PSVpaGaNU(J)
       ! Multiplicar partes de la m
       if(.not. allocated(CoefparGa)) allocate(CoefparGa(tam,N+1,2*nmax,2,2))
       if(.not. allocated(CoefparNu)) allocate(CoefparNu(tam,N+1,2*nmax,2,2))
       call PSVMatAporGaNU(J)

      end if!psv ............................................
      if (SH) then
         tam = size(Ak,1)
      Do ik = 1,vecNK(J)!k0!nmax+1
      ! k positivo
         pointAp => Ak(1:tam,1:tam,ik)
         pt_k => k_vec(ik)!        print*,ik, 2*nmax - (ik-2)
         call globalmatrix_SH(pointAp,pt_k,ik)
         call inverseA(pointAp,pt_ipivA,pt_workA,tam)
      end do ! ik
       k0 = vecNK(J);
      call intrplr_gloMat(k0,15,pt_cOME_i,pt_ipivA,pt_workA)
      ! k negativo (es simétrico)
         Ak(1:tam,1:tam,Nmax+2:2*nmax) = &
         Ak(1:tam,1:tam,Nmax:2:-1)
      end if!sh


      ! barra de avance
        write(6,'(A)', ADVANCE = "NO") repeat(char(8),60)
        write(6,'(A)', ADVANCE = "NO") repeat(char(8),17) !eta
        write(6,'(A)', ADVANCE = "NO") "["
        write(6,'(A,A)', ADVANCE = "NO") &
        repeat("X",int((58.0/NFREC)*(NFREC+1-J))),&
        repeat("_",58-int((58.0/NFREC)*(NFREC+1-J)))
        write(6,'(A)', ADVANCE = "NO") "]"
end subroutine preProcessAtThisfrec
