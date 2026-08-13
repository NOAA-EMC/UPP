! This is a test program for UPP.
!
! This program tests the ESAT() function.
!
! Alyson Stahl, 7/2026
program test_esat
    implicit none

    real, parameter :: tol = 1.0e-5
    integer, parameter :: ntests = 10
    ! Value from the ESAT function to avoid floating underflow
    real*4, parameter :: DEFAULT_ESAT = 3.777647E-05
    !
    integer :: i, res
    real*4 :: T(ntests), ESAT_OUT(ntests), EXP_ESAT(ntests)
    real*4 :: FLAG, FLG, k
    
    interface 
        real*4 function ESAT(T,FLAG,FLG)
            real*4, intent(in) :: T, FLAG, FLG
        end function ESAT
    end interface

    ! NOTE: FLG is not actually used by the function
    FLAG = 0.
    FLG = 0.

    ! Test Case 1: T < -273.15 => T + 273.15 < 0. (returns flag value)
    T(1) = -274.0
    EXP_ESAT(1) = FLAG

    ! Test Case 2: -273.15 <= T < -100 => T + 273.15 < 173.15 (set to default value to 
    ! avoid floating underflow)
    T(2) = -200.0
    EXP_ESAT(2) = DEFAULT_ESAT

    ! Test Case 3: -100 <= T < 100 => T + 273.15 > 173.15
    T(3) = 50.0
    k = T(3) + 273.15
    EXP_ESAT(3) = exp(26.660820-0.0091379024*k-6106.3960/k)

    ! Test Case 4: 100 <= T < 173.15 (T not adjusted, returns default value to 
    ! avoid floating underflow)
    T(4) = 150.0
    EXP_ESAT(4) = DEFAULT_ESAT

    ! Test Case 5: 173.15 <= T <= 373.15 (T not adjusted)
    T(5) = 200.0
    k = T(5)
    EXP_ESAT(5) = exp(26.660820-0.0091379024*k-6106.3960/k)

    ! Test Case 6: T > 373.15 (T not adjusted, returns flag value)
    T(6) = 400.0
    EXP_ESAT(6) = FLAG

    ! Test Case 7: T = -273.15 => T + 273.15 = 0 (Boundary Case)
    T(7) = -273.15
    EXP_ESAT(7) = DEFAULT_ESAT

    ! Test Case 8: T = -100 => T + 273.15 = 173.15 (Boundary Case)
    T(8) = -100.0
    k = T(8) + 273.15
    EXP_ESAT(8) = exp(26.660820-0.0091379024*k-6106.3960/k)

    ! Test Case 9: T = 100. (Boundary Case)
    T(9) = 100.0
    EXP_ESAT(9) = DEFAULT_ESAT

    ! Test Case 10: T = 373.15 (Boundary Case)
    T(10) = 373.15
    k = T(10)
    EXP_ESAT(10) = exp(26.660820-0.0091379024*k-6106.3960/k)

    res = 0
    do i = 1, ntests
        ESAT_OUT(i) = ESAT(T(i), FLAG, FLG)
        if (abs(ESAT_OUT(i) - EXP_ESAT(i)) > tol) then
            print *, "Test case ", i, " failed: Returned ", ESAT_OUT(i), &
                ", but expected ", EXP_ESAT(i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10

    print *, "SUCCESS!"
end program test_esat