! This is a test program for UPP.
!
! This program tests some of the smaller functions and subroutines in 
! the upp_physics module.
!
! Tested functions/subroutines:
!   - TVIRTUAL()
!   - FPVSNEW()
!
! Alyson Stahl, 5/2026
program test_upp_physics
    use upp_physics, only: FPVSNEW, TVIRTUAL
    implicit none

    real, parameter :: tol = 1.0e-8
    ! Number of test points for TVIRTUAL(). 
    integer, parameter :: N_TV = 5
    ! Number of test points for FPVSNEW().
    integer, parameter :: N_FPVS = 5
    !
    real :: T(N_TV), Q(N_TV), TV(N_TV), EXP_TV(N_TV)
    real :: T_FPVS(N_FPVS), EXP_SVP(N_FPVS), SVP(N_FPVS)
    !
    integer :: i, res

    print *, "Testing TVIRTUAL()..."

    ! No branches or edge cases in TVIRTUAL(), so just covering a range a reasonable inputs.
    T(1) = 0.0
    Q(1) = 0.0
    EXP_TV(1) = 0.0

    T(2) = 220.0
    Q(2) = 0.0
    EXP_TV(2) = 220.0

    T(3) = 280.0
    Q(3) = 0.005
    EXP_TV(3) = 280.0 * (1 + 0.608 * 0.005)

    T(4) = 295.0
    Q(4) = 0.015
    EXP_TV(4) = 295.0 * (1 + 0.608 * 0.015)

    T(5) = 305.0
    Q(5) = 0.025
    EXP_TV(5) = 305.0 * (1 + 0.608 * 0.025)

    TV = TVIRTUAL(T, Q)

    res = 0
    do i = 1, N_TV
        if (abs(TV(i) - EXP_TV(i)) > tol) then
            print '(A,I0,A,A,ES24.10,A,ES24.10)', "TVIRTUAL() Failed for test ", i, ": ", &
                        "Expected ", EXP_TV(i), &
                        " but got ", TV(i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10

    print *, "Testing FPVSNEW()..."

    ! Test Case 1: x, xp1 < tice
    T_FPVS(1) = 200.0
    EXP_SVP(1) = 1.5914107859E-01

    ! Test Case 2: x == tice
    T_FPVS(2) = 253.16
    EXP_SVP(2) = 1.0326692963E+02

    ! Test Case 3: tice < x, xp1 < tliq
    T_FPVS(3) = 270.0
    EXP_SVP(3) = 4.8179052734E+02

    ! Test Case 4: x == tliq
    T_FPVS(4) = 273.16
    EXP_SVP(4) = 6.1078002930E+02

    ! Test Case 5: x, xp1 > tliq
    T_FPVS(5) = 290.0
    EXP_SVP(5) = 1.9149111328E+03

    SVP = FPVSNEW(T_FPVS)

    res = 0
    do i = 1, N_FPVS
        if (abs(SVP(i) - EXP_SVP(i)) > tol) then
            print *, "FPVSNEW() Failed for test ", i, ": ", &
                        "Expected ", EXP_SVP(i), &
                        " but got ", SVP(i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 20

    print *, "SUCCESS!"
end program test_upp_physics