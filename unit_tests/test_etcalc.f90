! This is a test program for UPP.
!
! This program tests the ETCALC() subroutine.
!
! Alyson Stahl, 1/2026
program test_etcalc
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: ntests = 24
    integer :: i, res
    integer, dimension(1:ntests) :: ISOIL
    real, dimension(1:ntests) :: ETA, ETP, ESD, VEGFAC, SMC
    real, dimension(1:ntests) :: CMC, EC, EDIR, ETRANS, ESNOW, SMCDRY, SMCMAX
    real, dimension(1:ntests) :: EXP_CMC, EXP_EC, EXP_EDIR, EXP_ETRANS, EXP_ESNOW, &
                                    EXP_SMCDRY, EXP_SMCMAX

    interface 
        subroutine ETCALC(ETA, ETP, ESD, VEGFAC, ISOIL, SMC, CMC,        &
                          EC, EDIR, ETRANS, ESNOW, SMCDRY, SMCMAX)
        integer, intent(in) :: ISOIL
        real, intent(in) :: ETA, ETP, ESD, VEGFAC, SMC
        real, intent(inout) :: CMC
        real, intent(out) ::  EC, EDIR, ETRANS, ESNOW, SMCDRY, SMCMAX
        end subroutine ETCALC
    end interface

    ! Initialize arrays with default test case values
    ETA = 150.0
    ETP = 200.0
    ESD = 0.0
    VEGFAC = 0.5
    SMC = 0.15
    CMC = 0.00025
    ISOIL = 1
    ! We want to test every soil type. ISOIL = 1 will be used for additional test cases.
    do i = 1, 19
        ISOIL(i) = i
    end do

    ESD(20) = 0.5 ! Test Case where ESD > 0, expect ESNOW = ETA
    SMC(21) = 0.4 ! Test Case where SMC > SMCMAX for ISOIL = 1
    SMC(22) = 0.023 ! Test Case where SMC = SMCDRY for ISOIL = 1
    CMC(23) = 0.006 ! Test Case where CMC > CMCMAX
    ETA(24) = 0.0 ! Test Case where ETRANS = ETA - EDIR - EC < 0 (clipped to ETRANS = 0)

    EXP_CMC = 0.00025
    EXP_CMC(23) = 0.0005

    EXP_EC = 7.0710678101E+01
    EXP_EC(20) = 0.0
    EXP_EC(23) = 100.0

    EXP_ESNOW = 0.0
    EXP_ESNOW(20) = 150.0

    EXP_SMCDRY = (/ 0.023, 0.028, 0.047, 0.084, 0.084, 0.066, &
                0.069, 0.120, 0.103, 0.100, 0.126, 0.135, &
                0.069, 0.000, 0.012, 0.028, 0.135, 0.012, &
                0.023, 0.023, 0.023, 0.023, 0.023, 0.023 /)

    EXP_SMCMAX = (/ 0.395, 0.421, 0.434, 0.476, 0.476, 0.439, &
                0.404, 0.464, 0.465, 0.406, 0.468, 0.457, &
                0.464, 0.000, 0.200, 0.421, 0.457, 0.200, &
                0.395, 0.395, 0.395, 0.395, 0.395, 0.395 /)

    EXP_EDIR = (/ 1.1655249596E+01, 9.6368389130E+00, 7.0835776329E+00, &
                2.8347566128E+00, 2.8347566128E+00, 5.0715522766E+00, &
                5.8462915421E+00, 7.6054668427E-01, 1.6856940985E+00, &
                2.6699137688E+00, 4.9245935678E-01, 2.1700558066E-01, &
                4.2050962448E+00, 1.0000000000E+02, 5.3881855011E+01, &
                9.6368389130E+00, 2.1700558066E-01, 5.3881855011E+01, &
                1.1655249596E+01, 0.0000000000E+00, 1.0000000000E+02, &
                0.0000000000E+00, 1.1655249596E+01, 1.1655249596E+01 /)

    EXP_ETRANS = (/ 6.7634078979E+01, 6.9652481079E+01, 7.2205749512E+01, &
                7.6454559326E+01, 7.6454559326E+01, 7.4217773438E+01, &
                7.3443023682E+01, 7.8528778076E+01, 7.7603622437E+01, &
                7.6619415283E+01, 7.8796859741E+01, 7.9072311401E+01, &
                7.5084228516E+01, 0.0000000000E+00, 2.5407470703E+01, &
                6.9652481079E+01, 7.9072311401E+01, 2.5407470703E+01, &
                6.7634078979E+01, 0.0000000000E+00, 0.0000000000E+00, &
                7.9289321899E+01, 3.8344757080E+01, 0.0000000000E+00 /)
    
    res = 0
    do i = 1, ntests
        call ETCALC(ETA(i),ETP(i),ESD(i),VEGFAC(i),ISOIL(i),SMC(i),CMC(i),        &
     &                  EC(i),EDIR(i),ETRANS(i),ESNOW(i),SMCDRY(i),SMCMAX(i))

        if (abs(CMC(i) - EXP_CMC(i)) > tol) then
            print *, "Test ", i, " failed: CMC = ", CMC(i), " expected ", EXP_CMC(i)
            res = 1
        end if
        if (abs(EC(i) - EXP_EC(i)) > tol) then
            print *, "Test ", i, " failed: EC = ", EC(i), " expected ", EXP_EC(i)
            res = 1
        end if
        if (abs(EDIR(i) - EXP_EDIR(i)) > tol) then
            print *, "Test ", i, " failed: EDIR = ", EDIR(i), " expected ", EXP_EDIR(i)
            res = 1
        end if
        if (abs(ETRANS(i) - EXP_ETRANS(i)) > tol) then
            print *, "Test ", i, " failed: ETRANS = ", ETRANS(i), " expected ", EXP_ETRANS(i)
            res = 1
        end if
        if (abs(ESNOW(i) - EXP_ESNOW(i)) > tol) then
            print *, "Test ", i, " failed: ESNOW = ", ESNOW(i), " expected ", EXP_ESNOW(i)
            res = 1
        end if
        if (abs(SMCDRY(i) - EXP_SMCDRY(i)) > tol) then
            print *, "Test ", i, " failed: SMCDRY = ", SMCDRY(i), " expected ", EXP_SMCDRY(i)
            res = 1
        end if
        if (abs(SMCMAX(i) - EXP_SMCMAX(i)) > tol) then
            print *, "Test ", i, " failed: SMCMAX = ", SMCMAX(i), " expected ", EXP_SMCMAX(i)
            res = 1
        end if

        if (res .ne. 0) stop 10
    end do


    print *, "SUCCESS!"
end program test_etcalc