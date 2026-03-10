! This is a test program for UPP.
!
! This program tests subroutines PARA_RANGE() and PARA_RANGE2() in PARA_RANGE.f.
!
! Alyson Stahl, 2/2026
program test_para_range
    implicit none

    ! Number of tests for PARA_RANGE() and PARA_RANGE2(), respectively.
    integer, parameter :: n_tests1 = 2, n_tests2 = 8
    integer :: i, res
    integer, dimension(n_tests1) :: N1, N2, NPROCS, IRANK
    integer, dimension(n_tests2) :: IM, JM, NX, NY, NRANK
    integer, dimension(n_tests1) :: ISTA1, IEND1, EXP_ISTA1, EXP_IEND1
    integer, dimension(n_tests2) :: ISTA2, IEND2, JSTA2, JEND2, EXP_ISTA2, & 
        EXP_IEND2, EXP_JSTA2, EXP_JEND2

    interface 
        subroutine PARA_RANGE(N1, N2, NPROCS, IRANK, ISTA, IEND)
            integer, intent(in) :: N1, N2, NPROCS, IRANK
            integer, intent(out) :: ISTA, IEND
        end subroutine PARA_RANGE
        subroutine PARA_RANGE2(IM, JM, NX, NY, NRANK, ISTA, IEND, JSTA, JEND)
            integer, intent(in) :: IM, JM, NX, NY, NRANK
            integer, intent(out) :: ISTA, IEND, JSTA, JEND
        end subroutine PARA_RANGE2
    end interface

    ! Testing PARA_RANGE()
    ! Test Case 1: iwork2 <= irank
    N1(1) = 1
    N2(1) = 10
    NPROCS(1) = 4
    IRANK(1) = 2
    EXP_ISTA1(1) = 7
    EXP_IEND1(1) = 8

    ! Test Case: iwork2 > irank
    N1(2) = 1
    N2(2) = 10
    NPROCS(2) = 4
    IRANK(2) = 1
    EXP_ISTA1(2) = 4
    EXP_IEND1(2) = 6

    res = 0
    do i = 1, n_tests1
        call PARA_RANGE(N1(i), N2(i), NPROCS(i), IRANK(i), ISTA1(i), IEND1(i))
        if (ISTA1(i) .ne. EXP_ISTA1(i)) then
            print *, "Test ", i, " Failed: Expected ISTA=", EXP_ISTA1(i), " but got ISTA=", ISTA1(i)
            res = 1
        end if
        if (IEND1(i) .ne. EXP_IEND1(i)) then
            print *, "Test ", i, " Failed: Expected IEND=", EXP_IEND1(i), " but got IEND=", IEND1(i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10

    ! Testing PARA_RANGE2()
    IM = 8
    JM = 8
    NX = 2
    NY = 2

    ! Test Cases 1 to 4: Check all "corners" of the domain with same (IM, JM, NX, NY).
    ! Top-left corner
    NRANK(1) = 0
    EXP_ISTA2(1) = 1
    EXP_IEND2(1) = 4
    EXP_JSTA2(1) = 1
    EXP_JEND2(1) = 4

    ! Top-right corner
    NRANK(2) = 1 ! NX-1
    EXP_ISTA2(2) = 5
    EXP_IEND2(2) = 8
    EXP_JSTA2(2) = 1
    EXP_JEND2(2) = 4

    ! Bottom-left corner
    NRANK(3) = 2 ! NX*(NY-1)
    EXP_ISTA2(3) = 1
    EXP_IEND2(3) = 4
    EXP_JSTA2(3) = 5
    EXP_JEND2(3) = 8

    ! Bottom-right corner
    NRANK(4) = 3 ! NX*NY-1
    EXP_ISTA2(4) = 5
    EXP_IEND2(4) = 8
    EXP_JSTA2(4) = 5
    EXP_JEND2(4) = 8

    ! Test Cases 5 to 7: 1D-degenerate decompositions
    NRANK(5:7) = 0

    ! NX = 1, NY > 1
    NX(5) = 1
    EXP_ISTA2(5) = 1
    EXP_IEND2(5) = 8
    EXP_JSTA2(5) = 1
    EXP_JEND2(5) = 4

    ! NX > 1, NY = 1
    NY(6) = 1
    EXP_ISTA2(6) = 1
    EXP_IEND2(6) = 4
    EXP_JSTA2(6) = 1
    EXP_JEND2(6) = 8
    
    ! NX = 1, NY = 1
    NX(7) = 1
    NY(7) = 1
    EXP_ISTA2(7) = 1
    EXP_IEND2(7) = 8
    EXP_JSTA2(7) = 1
    EXP_JEND2(7) = 8

    ! Test Case 8: IM and JM not divisible by NX and NY, respectively
    IM(8) = 10
    JM(8) = 7
    NRANK(8) = 5
    NX(8) = 4
    NY(8) = 3
    EXP_ISTA2(8) = 4
    EXP_IEND2(8) = 6
    EXP_JSTA2(8) = 4
    EXP_JEND2(8) = 5

    res = 0
    do i = 1, n_tests2
        call PARA_RANGE2(IM(i), JM(i), NX(i), NY(i), NRANK(i), ISTA2(i), IEND2(i), JSTA2(i), JEND2(i))
        if (ISTA2(i) .ne. EXP_ISTA2(i)) then
            print *, "Test ", i, " Failed: Expected ISTA=", EXP_ISTA2(i), " but got ISTA=", ISTA2(i)
            res = 1
        end if
        if (IEND2(i) .ne. EXP_IEND2(i)) then
            print *, "Test ", i, " Failed: Expected IEND=", EXP_IEND2(i), " but got IEND=", IEND2(i)
            res = 1
        end if
        if (JSTA2(i) .ne. EXP_JSTA2(i)) then
            print *, "Test ", i, " Failed: Expected JSTA=", EXP_JSTA2(i), " but got JSTA=", JSTA2(i)
            res = 1
        end if
        if (JEND2(i) .ne. EXP_JEND2(i)) then
            print *, "Test ", i, " Failed: Expected JEND=", EXP_JEND2(i), " but got JEND=", JEND2(i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 20

    print *, "SUCCESS!"
end program test_para_range