! This is a test program for UPP.
!
! This program tests the ZENSUN() subroutine.
!
! Note: There is one branch of the ZENSUN() logic (di == 74) that
! appears to be unreachable.
!
! Alyson Stahl, 1/2026
program test_zensun
    use kinds, only: r_kind,i_kind
    implicit none
    
    real, parameter :: tol = 1.0e-8
    integer, parameter :: ntests = 10
    integer :: i, res
    ! Input
    real(r_kind), parameter :: PI = 3.14159265358979323846
    integer(i_kind), dimension(ntests) :: DAY
    real(r_kind), dimension(ntests) :: TIME, LAT, LON
    ! Output
    real(r_kind), dimension(ntests) :: SUN_ZENITH, SUN_AZIMUTH
    real(r_kind), dimension(ntests) :: EXP_SUN_ZENITH, EXP_SUN_AZIMUTH

    interface
        subroutine ZENSUN(DAY, TIME, LAT, LON, PI, SUN_ZENITH, SUN_AZIMUTH)
            use kinds, only: r_kind,i_kind
            integer(i_kind), intent(in) :: DAY
            real(r_kind), intent(in) :: TIME, LAT, LON, PI
            real(r_kind), intent(out) :: SUN_ZENITH, SUN_AZIMUTH
        end subroutine ZENSUN
    end interface

    ! 1) di == 1 lower boundary: tt = 1.0 (nday(1))
    DAY(1)  = 1
    TIME(1) = 0.0
    LAT(1)  = 0.0
    LON(1)  = 0.0

    ! 2) di == 1 upper boundary: tt = 6.0 (nday(2)), first matching interval is di=1
    DAY(2)  = 6
    TIME(2) = 0.0
    LAT(2)  = 45.0
    LON(2)  = 0.0

    ! 3) di == 2 lower boundary just above 6: tt ≈ 6.0417
    DAY(3)  = 6
    TIME(3) = 1.0
    LAT(3)  = 0.0
    LON(3)  = 30.0

    ! 4) di == 2 typical interior: tt in [6,11]
    DAY(4)  = 8
    TIME(4) = 12.0
    LAT(4)  = -30.0
    LON(4)  = -60.0

    ! 5) di in [3,72] typical interior: mid-year day
    DAY(5)  = 100
    TIME(5) = 6.0
    LAT(5)  = 50.0
    LON(5)  = 10.0

    ! 6) di in [3,72] exact LOWTRAN point: tt = 171.0 (summer solstice)
    DAY(6)  = 171
    TIME(6) = 0.0
    LAT(6)  = 23.5
    LON(6)  = 0.0

    ! 7) di in [3,72] equinox region: tt = 266.0 (fall equinox)
    DAY(7)  = 266
    TIME(7) = 12.0
    LAT(7)  = 0.0
    LON(7)  = 0.0

    ! 8) di == 73 lower boundary: tt = 361.0
    DAY(8)  = 361
    TIME(8) = 0.0
    LAT(8)  = 0.0
    LON(8)  = 0.0

    ! 9) di == 73 interior: tt in [361,366]
    DAY(9)  = 363
    TIME(9) = 18.0
    LAT(9)  = -60.0
    LON(9)  = 120.0

    ! 10) di == 73 upper boundary: tt = 366.0 (nday(74))
    DAY(10)  = 366
    TIME(10) = 0.0
    LAT(10)  = 10.0
    LON(10)  = -90.0
    
    EXP_SUN_ZENITH = (/ &
         1.5711235046E+02, 1.5745872498E+02, 1.3195394897E+02, 5.5247062683E+01, &
         7.8017562866E+01, 1.3316987610E+02, 1.8452537060E+00, 1.5675169373E+02, &
         1.0507037354E+02, 1.8356952667E+01 /)
         
    EXP_SUN_AZIMUTH = (/ &
        -1.7825204468E+02, -3.2312366962E+00, 1.2082527924E+02, 9.7385551453E+01, &
         9.2414581299E+01, -4.0916731954E-01, -8.7924919128E+01, -1.7962820435E+02, &
        -9.5679412842E+01, 8.6554504395E+01 /)

    res = 0
    do i = 1, ntests
        call ZENSUN(DAY(i), TIME(i), LAT(i), LON(i), PI, SUN_ZENITH(i), SUN_AZIMUTH(i))
        if ( abs(SUN_ZENITH(i) - EXP_SUN_ZENITH(i)) > tol ) then
            print *, "Test failed for SUN_ZENITH(", i, "): ", &
                     "Expected = ", EXP_SUN_ZENITH(i), &
                     ", Computed = ", SUN_ZENITH(i)
            res = 1
        end if
        if ( abs(SUN_AZIMUTH(i) - EXP_SUN_AZIMUTH(i)) > tol ) then
            print *, "Test failed for SUN_AZIMUTH(", i, "): ", &
                     "Expected = ", EXP_SUN_AZIMUTH(i), &
                     ", Computed = ", SUN_AZIMUTH(i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10

    print *, "SUCCESS!"
end program test_zensun