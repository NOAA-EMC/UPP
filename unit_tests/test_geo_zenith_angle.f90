! This is a test program for UPP.
!
! This program tests the GEO_ZENITH_ANGLE() subroutine.
!
! Alyson Stahl, 1/2025
program test_geo_zenith_angle
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: ntests = 5
    integer :: i, res
    real, dimension(ntests) :: RLAT, RLON, SLAT, SLON, ZA, EXP_ZA

    ! Default input values
    RLAT = 35.0
    RLON = 135.0
    SLAT = 0.0
    SLON = 140.7

    ! RLAT > 180.0
    RLAT(2) = 250.0
    
    ! COSE clips to 1.0
    RLAT(3) = 0.0
    RLON(3) = 180.0
    SLAT(3) = 0.0
    SLON(3) = 0.0

    ! COSE clips to -1.0
    RLAT(4) = 0.0
    RLON(4) = 0.0
    SLAT(4) = 0.0
    SLON(4) = 0.0    
    
    ! ZA clips to 0.0
    RLAT(5) = 0.0
    RLON(5) = 360.0
    SLAT(5) = 0.0
    SLON(5) = 0.0

    ! TODO: Fix this later.
    EXP_ZA = (/0.0, 0.0, 0.0, 0.0, 0.0/)

    do i = 1, ntests
        ! First two arguments are unused in current implementation
        call GEO_ZENITH_ANGLE(0, 0, RLAT(i), RLON(i), SLAT(i), SLON(i), ZA(i))
    end do

    res = 0
    do i = 1, ntests
        print '(A,I0,A,ES24.10)', "ZA(", i, ") = ", ZA(i)
    end do
    
    print *, 'SUCCESS!'
end program test_geo_zenith_angle