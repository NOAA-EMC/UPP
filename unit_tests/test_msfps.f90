! This is a test program for UPP.
!
! This program tests the MSFPS() subroutine.
!
! Alyson Stahl, 3/2026
program test_msfps
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: ntests = 10
    integer :: i, res
    real :: LAT(ntests), TRUELAT1(ntests), MSF(ntests), EXP_MSF(ntests)

    interface
        subroutine MSFPS(LAT, TRUELAT1, MSF)
            real, intent(in) :: LAT, TRUELAT1
            real, intent(out) :: MSF
        end subroutine MSFPS
    end interface

    ! Test Case 1: 0 < TRUELAT1 < 90, LAT = TRUELAT1
    TRUELAT1(1) = 60.0
    LAT(1) = TRUELAT1(1)
    EXP_MSF(1) = 1.0
    
    ! Test Case 2: -90 < TRUELAT1 < 0, LAT = TRUELAT1
    TRUELAT1(2) = -60.0
    LAT(2) = TRUELAT1(2)
    EXP_MSF(2) = 1.0

    ! Test Case 3: TRUELAT1 = 0, LAT = 0
    TRUELAT1(3) = 0.0
    LAT(3) = 0.0
    EXP_MSF(3) = 1.0

    ! Test Case 4: 0 < TRUELAT1 < 90, LAT = 90
    TRUELAT1(4) = 30.0
    LAT(4) = 90.0
    EXP_MSF(4) = 0.75

    ! Test Case 5: -90 < TRUELAT1 < 0, LAT = -90
    TRUELAT1(5) = -30.0
    LAT(5) = -90.0
    EXP_MSF(5) = 0.75

    ! Test Case 6: 0 <= LAT < TRUELAT1 < 90
    TRUELAT1(6) = 45.0
    LAT(6) = 30.0
    EXP_MSF(6) = 1.1380711794

    ! Test Case 7: -90 < TRUELAT1 < LAT <= 0
    TRUELAT1(7) = -45.0
    LAT(7) = -30.0
    EXP_MSF(7) = 1.1380711794

    ! Test Case 8: 0 < TRUELAT1 < LAT <= 90
    TRUELAT1(8) = 60.0
    LAT(8) = 80.0
    EXP_MSF(8) = 9.4015425444E-01

    ! Test Case 9: -90 <= LAT < TRUELAT1 < 0
    TRUELAT1(9) = -60.0
    LAT(9) = -80.0
    EXP_MSF(9) = 9.4015425444E-01

    ! Test Case 10: 0 < TRUELAT1 < 90, LAT ~ -90
    TRUELAT1(10) = 60.0
    LAT(10) = -89.0
    EXP_MSF(10) = 1.2253116211E+04

    res = 0
    do i = 1, ntests
        call MSFPS(LAT(i), TRUELAT1(i), MSF(i))
        if (abs(MSF(i) - EXP_MSF(i)) > tol) then
            print *, 'Test case ', i, ' failed: MSF = ', MSF(i), ' expected ', EXP_MSF(i)
            res = 1
        end if
    end do
    
    if (res .ne. 0) stop 10
    
    print *, 'SUCCESS!'
end program test_msfps