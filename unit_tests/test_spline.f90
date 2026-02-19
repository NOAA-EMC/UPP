! This is a test program for UPP.
!
! This program tests the SPLINE() subroutine.
!
! Alyson Stahl, 2/2026
program test_spline
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: JTB = 5
    integer :: i, res
    integer :: NOLD, NNEW
    real :: XOLD(JTB), YOLD(JTB), XNEW(JTB), P(JTB), Q(JTB), Y2(JTB)
    real :: YNEW(JTB), EXP_YNEW_1(JTB), EXP_YNEW_2(JTB)

    interface
        subroutine SPLINE(JTB,NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)
            integer,intent(in) :: JTB,NOLD,NNEW
            real,dimension(JTB),intent(in) ::  XOLD,YOLD,XNEW 
            real,dimension(JTB),intent(inout) :: P,Q,Y2
            real,dimension(JTB),intent(out) ::  YNEW
        end subroutine SPLINE
    end interface

    ! Test Case: Standard case where NOLD > 3
    EXP_YNEW_1 = 0.0
    EXP_YNEW_1(1) =         0.34558820724
    EXP_YNEW_1(2) =         0.61488968134
    EXP_YNEW_1(3) =         2.2107839584
    EXP_YNEW_1(4) =         25.0

    ! Choose NOLD>3 to exercise forward sweep; NNEW>=2 to test reuse/recompute
    NOLD = 5
    NNEW = 4

    ! Ascending, nonuniform XOLD to exercise variable spacing
    XOLD = (/0.0, 1.0, 2.5, 4.0, 5.0/)

    ! Nonlinear YOLD values
    do i = 1, NOLD
        YOLD(i) = XOLD(i)**2
    end do

    ! Natural spline boundary conditions: Y2(1)=0 and Y2(NOLD)=0
    Y2 = 0.0
    P = 0.0
    Q = 0.0

    ! XNEW sequence to cover interior, reuse, recompute, and right endpoint
    XNEW(1) = 0.5          ! interior in first interval (XOLD(1), XOLD(2))
    XNEW(2) = 0.75         ! same interval as previous (reuse coefficients)
    XNEW(3) = 1.5          ! moves to next interval (recompute coefficients)
    XNEW(4) = XOLD(NOLD)   ! right endpoint path (assign YNEW=YOLD(NOLD))

    ! Initialize YNEW so that unused indices maintain consistent values.
    YNEW = 0.0

    call SPLINE(JTB,NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)

    res = 0
    do i = 1, JTB
        if (abs(YNEW(i) - EXP_YNEW_1(i)) > tol) then
            print *, 'YNEW Failed for test', i, ': ', &
                        'Expected ', EXP_YNEW_1(i), &
                        ' but got ', YNEW(i)
            res = 1
        end if
    end do
    
    if (res .ne. 0) stop 10

    ! Test Case: Edge case where NOLD=3
    NOLD = 3

    EXP_YNEW_2 = 0.0
    EXP_YNEW_2(1) = 0.31249994040
    EXP_YNEW_2(2) = 0.58593744040
    EXP_YNEW_2(3) = 2.3333330154
    EXP_YNEW_2(4) = 6.25
    EXP_YNEW_2(5) = 0.0

    ! Reinitialize inout and output variables
    Y2 = 0.0
    P = 0.0
    Q = 0.0
    YNEW = 0.0

    call SPLINE(JTB,NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)

    do i = 1, JTB
        if (abs(YNEW(i) - EXP_YNEW_2(i)) > tol) then
            print *, 'YNEW Failed for test', i, ': ', &
                        'Expected ', EXP_YNEW_2(i), &
                        ' but got ', YNEW(i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 20
    
    print *, 'SUCCESS!'
end program test_spline