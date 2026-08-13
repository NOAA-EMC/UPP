! This is a test program for UPP.
!
! This program tests the TDOFESAT() function.
!
! Alyson Stahl, 7/2026
program test_tdofesat
    implicit none

    real, parameter :: tol = 1.0e-5
    integer, parameter :: ntests = 7
    ! Values from the TDOFESAT function to avoid underflow/overflow
    real, parameter :: lim1 = 3.777647E-05, lim2 = 980.5386
    real*4, parameter :: DEFAULT_TDOFESAT = 173.15
    !
    integer :: i, res
    real*4 :: ES(ntests), TDOFESAT_OUT(ntests), EXP_TDOFESAT(ntests)
    real*4 :: FLAG, FLG, b
    
    interface 
        real*4 function TDOFESAT(ES,FLAG,FLG)
            real*4, intent(in) :: ES, FLAG, FLG
        end function TDOFESAT
    end interface

    ! NOTE: FLG is not actually used by the function
    FLAG = 0.
    FLG = 0.

    ! Test Case 1: ES < 0. (returns flag value)
    ES(1) = -1.0
    EXP_TDOFESAT(1) = FLAG

    ! Test Case 2: 0 <= ES < lim1 (returns default value)
    ES(2) = 1.0e-6
    EXP_TDOFESAT(2) = DEFAULT_TDOFESAT

    ! Test Case 3: lim1 <= ES <= lim2 (returns computed value)
    ES(3) = 200.0
    b = 26.66082 - alog(ES(3))
    EXP_TDOFESAT(3) = (b-sqrt(b*b-223.1986)) / 0.0182758048

    ! Test Case 4: ES > lim2 (returns flag value)
    ES(4) = 1000.0
    EXP_TDOFESAT(4) = FLAG

    ! Test Case 5: ES = 0. (Boundary Case)
    ES(5) = 0.0
    EXP_TDOFESAT(5) = DEFAULT_TDOFESAT

    ! Test Case 6: ES = lim1 (Boundary Case)
    ES(6) = lim1
    b = 26.66082 - alog(ES(6))
    EXP_TDOFESAT(6) = (b-sqrt(b*b-223.1986)) / 0.0182758048

    ! Test Case 7: ES = lim2 (Boundary Case)
    ES(7) = lim2
    b = 26.66082 - alog(ES(7))
    EXP_TDOFESAT(7) = (b-sqrt(b*b-223.1986)) / 0.0182758048

    res = 0
    do i = 1, ntests
        TDOFESAT_OUT(i) = TDOFESAT(ES(i), FLAG, FLG)
        if (abs(TDOFESAT_OUT(i) - EXP_TDOFESAT(i)) > tol) then
            print *, "Test case ", i, " failed: Returned ", TDOFESAT_OUT(i), &
                ", but expected ", EXP_TDOFESAT(i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10
    
    print *, "SUCCESS!"
end program test_tdofesat