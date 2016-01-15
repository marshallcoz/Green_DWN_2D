subroutine checarWisdom(a,b,c)
      use glovars, only : planNmaxF,planNmaxB,&
                             planNfrecF,planNfrecB,&
                             planNtimeF,planNtimeB
      use, intrinsic :: iso_c_binding
      use glovars, only : PrintNum
      include 'fftw3.f03'
      type(C_PTR) :: plan
      integer(C_INT) :: j,a,b,c,flag
      complex*16, dimension(a) :: Ua,Va !frec
      complex*16, dimension(b) :: Ub,Vb !nmax
      complex*16, dimension(c) :: Uc,Vc !npstime
      character(len=100) :: tx
      logical :: lexist
      write(tx,'(a,I0,a,I0,a,I0,a)') "wis",a,"-",b,"-",c,".dat"
      inquire(file=trim(tx),exist=lexist)
      if (lexist .eqv. .false.) then
      flag = FFTW_PATIENT ! FFTW_MEASURE
      plan = fftw_plan_dft_1d(a, Ua,Va, FFTW_FORWARD,flag)
      plan = fftw_plan_dft_1d(a, Ua,Va, FFTW_BACKWARD,flag)
      plan = fftw_plan_dft_1d(b, Ub,Vb, FFTW_FORWARD, flag)
      plan = fftw_plan_dft_1d(b, Ub,Vb, FFTW_BACKWARD, flag)
      plan = fftw_plan_dft_1d(c, Uc,Vc, FFTW_FORWARD, flag)
      plan = fftw_plan_dft_1d(c, Uc,Vc, FFTW_BACKWARD, flag)
      j = fftw_export_wisdom_to_filename(trim(tx) // C_NULL_CHAR)
      call fftw_destroy_plan(plan)
      if (j .eq. 0) then
      print*,"i=",j," (non-zero on success)"
      stop "Error al exportar wisdom en :checarWisdom"
      else
      print*, "Plans created, at"
      print*,trim(tx)
      stop "Re run program"
      end if
      print*,"Se export— wisdom ",trim(tx)
      else
      !importar
      j = fftw_import_wisdom_from_filename(trim(tx) // C_NULL_CHAR)
      if (j .eq. 0) then
      print*,"i=",j," (non-zero on success)"
      stop "Error al importar wisdom en :FFTW"
      else
      write(PrintNum,*)"  "
      write(PrintNum,*)"  le’do wisdom ",trim(tx)
      ! crear planos
      flag = FFTW_WISDOM_ONLY
      planNfrecF = fftw_plan_dft_1d(a, Ua,Va, FFTW_FORWARD,flag)
      planNfrecB = fftw_plan_dft_1d(a, Ua,Va, FFTW_BACKWARD,flag)
      planNmaxF = fftw_plan_dft_1d(b, Ub,Vb, FFTW_FORWARD, flag)
      planNmaxB = fftw_plan_dft_1d(b, Ub,Vb, FFTW_BACKWARD, flag)
      planNtimeF = fftw_plan_dft_1d(c, Uc,Vc, FFTW_FORWARD, flag)
      planNtimeB = fftw_plan_dft_1d(c, Uc,Vc, FFTW_BACKWARD, flag)
      write(PrintNum,*)"  planos creados"
      end if
      end if
      end subroutine checarWisdom


      subroutine vaciarWisdom
      use glovars, only : planNmaxF,planNmaxB,&
                             planNfrecF,planNfrecB,&
                             planNtimeF,planNtimeB
      use, intrinsic :: iso_c_binding
      include 'fftw3.f03'
      call fftw_destroy_plan(planNfrecF)
      call fftw_destroy_plan(planNfrecB)
      call fftw_destroy_plan(planNmaxF)
      call fftw_destroy_plan(planNmaxB)
      call fftw_destroy_plan(planNtimeF)
      call fftw_destroy_plan(planNtimeB)
      !call fftwf_forget_wisdom()
      !call fftwf_cleanup()
      end subroutine vaciarWisdom
