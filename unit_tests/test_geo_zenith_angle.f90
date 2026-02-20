! This is a test program for UPP.
!
! This program tests the GEO_ZENITH_ANGLE() subroutine.
!
! Alyson Stahl, 1/2026
program test_geo_zenith_angle
    implicit none

    real, parameter :: tol_abs = 1.0e-6, tol_rel = 1.0e-6
    integer, parameter :: ntests = 5
    integer :: i, res
    real, dimension(ntests) :: RLAT, RLON, SLAT, SLON, ZA, EXP_ZA

    interface 
        subroutine GEO_ZENITH_ANGLE(i, j, RLAT, RLON, SLAT, SLON, ZA)
            integer, intent(in) :: i, j
            real, intent(in) :: RLAT, RLON, SLAT, SLON
            real, intent(out) :: ZA
        end subroutine GEO_ZENITH_ANGLE
    end interface

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

    EXP_ZA = (/ 41.077945709,  117.58034515, 180.0, 0.0, 0.0 /)

    do i = 1, ntests
        ! First two arguments are unused in current implementation
        call GEO_ZENITH_ANGLE(0, 0, RLAT(i), RLON(i), SLAT(i), SLON(i), ZA(i))
    end do

    res = 0
    do i = 1, ntests
        if (abs(ZA(i) - EXP_ZA(i)) > max(tol_abs, tol_rel*abs(EXP_ZA(i)))) then
            print *, 'Test ', i, ' failed: computed ZA = ', ZA(i), &
                    ', expected ZA = ', EXP_ZA(i)
            res = 1
        end if
    end do
    
    if (res .ne. 0) stop 10
    
    print *, 'SUCCESS!'
end program test_geo_zenith_angle