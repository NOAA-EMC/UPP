! This is a test program for UPP.
!
! This program tests the XMYTW_POST() function.
!
! Alyson Stahl, 7/2026
program test_xmytw_post
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: ntests = 11
    !
    integer :: i, res
    real*4 :: T(ntests), TD(ntests), P(ntests)
    real*4 :: XMYTW_POST_OUT(ntests), EXP_XMYTW_POST(ntests)

    interface
        function XMYTW_POST(T, TD, P)
            real*4, intent(in) :: T, TD, P
        end function XMYTW_POST
    end interface

    ! Use same value of P for most test cases
    P = 0.0 

    ! Test Case 1: td > t (returns average of t and td)
    T(1) = 20.0
    TD(1) = 30.0
    EXP_XMYTW_POST(1) = 25.0

    ! Test Case 2: td = t (boundary case, returns average of t and td)
    T(2) = 20.0
    TD(2) = 20.0
    EXP_XMYTW_POST(2) = 20.0


    ! Test Case 3: td < t < 100.0 (standard case that returns value in Celsius)
    T(3) = 90.0
    TD(3) = 80.0
    EXP_XMYTW_POST(3) = 80.0

    ! Test Case 4: t == 100.0 (boundary case, should return in Kelvin)
    T(4) = 100.0
    TD(4) = 90.0
    EXP_XMYTW_POST(4) = 95.0

    ! Test Case 5: t > 100.0 & td < t (standard case that returns value in Kelvin)
    T(5) = 300.0
    TD(5) = 290.0
    EXP_XMYTW_POST(5) = 290.0

    ! Note: The input values below were calculated and validated in Python. They were specifically chosen
    ! for branch coverage and to create every possible return condition. If future changes to the code
    ! intentionally break these conditions, test cases should be updated accordingly to maintain 
    ! coverage.

    ! Test Case 6: Returns average in Kelvin with ed < -14.0
    T(6) = 200.0
    TD(6) =  150.0
    EXP_XMYTW_POST(6) = 175.0

    ! Test Case 7: Returns average in Kelvin with ed > 7.0
    T(7) = 600.0
    TD(7) = 500.0
    EXP_XMYTW_POST(7) = 550.0

    ! Test Case 8: Returns average in Kelvin with ew < -14.0
    T(8) = 4300.0
    TD(8) = 300.0
    EXP_XMYTW_POST(8) = 2300.0

    ! Test Case 9: Returns average in Kelvin with ew > 7.0
    T(9) = 1700.0
    TD(9) = 300.0
    EXP_XMYTW_POST(9) = 1000.0

    ! Test Case 10: Get ew < -14.0 within iteration - returns average in Kelvin
    T(10) = 4000.0
    TD(10) = 320.0
    P(10) = 500.0
    EXP_XMYTW_POST(10) = 2160.0

    ! Test Case 11: Get ew > 7.0 within iteration - returns average in Kelvin
    T(11) = 4000.0
    TD(11) = 320.0
    P(11) = 10.0
    EXP_XMYTW_POST(11) = 2160.0

    res = 0
    do i = 1, ntests
        XMYTW_POST_OUT(i) = XMYTW_POST(T(i), TD(i), P(i))
        if (abs(XMYTW_POST_OUT(i) - EXP_XMYTW_POST(i)) > tol) then
            print *, "Test case ", i, " failed: Returned ", XMYTW_POST_OUT(i), &
                ", but expected ", EXP_XMYTW_POST(i)
            res = 1
        end if
    end do
    
    if (res .ne. 0) stop 10 

    print *, "SUCCESS!"
end program test_xmytw_post