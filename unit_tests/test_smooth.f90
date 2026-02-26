! This is a test program for UPP.
!
! This program tests subroutines SMOOTH() and SMOOTHC() in SMOOTH.f.
!
! Alyson Stahl, 2/2026
program test_smooth
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: nx = 4, ny = 4
    integer :: i, j, res
    integer :: IX, IY
    real :: SMTH, HOLD(nx, 2), FIELD(nx, ny), FIELDC(nx, ny)
    real :: EXP_FIELD(nx, ny), EXP_FIELDC(nx, ny)

    interface 
        subroutine SMOOTH(FIELD, HOLD, IX, IY, SMTH)
            integer, intent(in) :: IX, IY
            real, intent(in) :: SMTH
            real, dimension(IX, 2), intent(in) :: HOLD
            real, dimension(IX, IY), intent(inout) :: FIELD
        end subroutine SMOOTH
        subroutine SMOOTHC(FIELD, HOLD, IX, IY, SMTH)
            integer, intent(in) :: IX, IY
            real, intent(in) :: SMTH
            real, dimension(IX, 2), intent(in) :: HOLD
            real, dimension(IX, IY), intent(inout) :: FIELD
        end subroutine SMOOTHC
    end interface

    IX = nx
    IY = ny
    SMTH = 0.5
    HOLD = 0.0

    do j = 1, ny
        do i = 1, nx
            FIELD(i, j) = 280.0 + 2.0 * real(j - 1) + real(i - 1)
            FIELDC(i, j) = 280.0 + 2.0 * real(j - 1) + real(i - 1)
        end do
    end do
    
    ! Test for SMOOTH()

    call SMOOTH(FIELD, HOLD, IX, IY, SMTH)

    print *, "SMOOTH() Results:"
    do i = 1, nx
        do j = 1, ny
            print '("(",I1,",",I1,"): ",ES24.10)', i, j, FIELD(i, j)
        end do
    end do  

    ! Test for SMOOTHC()

    call SMOOTHC(FIELDC, HOLDC, IX, IY, SMTHC)

    print *, "SMOOTHC() Results:"
    do i = 1, nx
        do j = 1, ny
            print '("(",I1,",",I1,"): ",ES24.10)', i, j, FIELDC(i, j)
        end do
    end do  

    print *, "SUCCESS!"
end program test_smooth