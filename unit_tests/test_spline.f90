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
    real :: YNEW(JTB), EXP_YNEW(JTB)

    interface
        subroutine SPLINE(JTB,NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)
            integer,intent(in) :: JTB,NOLD,NNEW
            real,dimension(JTB),intent(in) ::  XOLD,YOLD,XNEW 
            real,dimension(JTB),intent(inout) :: P,Q,Y2
            real,dimension(JTB),intent(out) ::  YNEW
        end subroutine SPLINE
    end interface

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

    call SPLINE(JTB,NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)

    do i = 1, JTB
        print '(A,I0,A,ES24.10)', 'YNEW(', i, ') = ', YNEW(i)
    end do
    print *, 'SUCCESS!'
end program test_spline