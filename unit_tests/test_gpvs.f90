! This is a test program for UPP.
!
! This program tests the functions and subroutines in GPVS.f.
!
! The function FPVS0() is intentionally skipped because it's
! no longer being used.
!
! FPVS() has mostly been replaced, but is still used in
! CALMICT.f, so it is still being tested here.
!
! Alyson Stahl, 1/2026
program test_gpvs
    use svptbl_mod, only: nx, tbpvs, tbpvs0
    implicit none

    real, parameter :: tol = 1.0e-8
    real, parameter :: XMIN = 180.0, XMAX = 330.0 ! From GPVS.f
    integer, parameter :: npts = 16 ! Number of points to check in svp tables
    integer, parameter :: ntests = 8 ! Number of tests for FPVS
    real, parameter :: TINC = (XMAX - XMIN) / (ntests - 1)
    integer :: i, j, res
    real :: T, SVP
    real, dimension(npts) :: EXP_TBPVS = (/&
        5.1082324717E-06, 3.1216961361E-05, 1.5894723765E-04, 6.9217011333E-04, &
        2.6336372830E-03, 8.9114867151E-03, 2.7213586494E-02, 7.5932152569E-02, &
        1.9561828673E-01, 4.6946316957E-01, 9.8994475603E-01, 1.9149112701E+00, &
        3.5241372585E+00, 6.2015328407E+00, 1.0480564117E+01, 1.7075515747E+01 /)
    real, dimension(npts) :: EXP_TBPVS0 = (/&
        1.2800095647E-05, 7.1109141572E-05, 3.2837456092E-04, 1.2950585224E-03, &
        4.4593093917E-03, 1.3652097434E-02, 3.7726629525E-02, 9.5303639770E-02, &
        2.2244605422E-01, 4.8410034180E-01, 9.8994475603E-01, 1.9149112701E+00, &
        3.5241372585E+00, 6.2015328407E+00, 1.0480564117E+01, 1.7075515747E+01 /)
    real, dimension(ntests) :: EXP_SVP = (/&
        5.1082324717E-06, 1.9790380611E-04, 3.7731970660E-03, 4.2686607689E-02, &
        3.2524418831E-01, 1.5943791866E+00, 5.7352380753E+00, 1.7075515747E+01 /)
        
    interface 
        subroutine GPVS()
        end subroutine GPVS
        real function FPVS(T)
            real, intent(in) :: T
        end function FPVS
    end interface

    ! No input or output. Just initializes tables tbpvs and tbpvs0
    call GPVS()

    ! Check results of TBPVS and TBPVS0 against expected values
    res = 0 
    do i = 1, npts
        j =  1 + (i-1) * 500
        if (abs(tbpvs(j) - EXP_TBPVS(i)) > tol) then
            print *, 'tbpvs Failed for test', i, ': ', &
                        'Expected ', EXP_TBPVS(i), &
                        ' but got ', tbpvs(j)
            res = 1
        end if
        if (abs(tbpvs0(j) - EXP_TBPVS0(i)) > tol) then
            print *, 'tbpvs0 Failed for test', i, ': ', &
                        'Expected ', EXP_TBPVS0(i), &
                        ' but got ', tbpvs0(j)
            res = 1
        end if
        if (res .ne. 0) stop 10
    end do

    do i = 1, ntests
        T = XMIN + TINC * (i - 1)
        SVP = FPVS(T)
        if (abs(SVP - EXP_SVP(i)) > tol) then
            print *, 'FPVS Failed for test', i, ': ', &
                        'Expected ', EXP_SVP(i), &
                        ' but got ', SVP
            stop 20
        end if
    end do

    print *, 'SUCCESS!'
end program test_gpvs