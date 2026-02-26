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

    EXP_FIELD = reshape([280.0, 281.0, 282.0, 283.0, &
                    282.0, 283.0, 284.0, 285.0, &
                    284.0, 285.0, 286.0, 287.0, &
                    286.0, 287.0, 288.0, 289.0], [nx, ny]) 
    EXP_FIELDC = reshape([281.0, 281.25, 282.0625, 282.265625, &
                    283.0, 283.0, 284.0, 284.0, &
                    285.0, 285.0, 286.0, 286.0, &
                    287.0, 287.25, 288.0625, 288.265625], [nx, ny])

    do i = 1, nx
        do j = 1, ny
            FIELD(i, j) = 280.0 + 2.0 * real(j - 1) + real(i - 1)
            FIELDC(i, j) = 280.0 + 2.0 * real(j - 1) + real(i - 1)
        end do
    end do
    
    ! Test for SMOOTH()
    
    call SMOOTH(FIELD, HOLD, IX, IY, SMTH)

    res = 0
    do i = 1, nx
        do j = 1, ny
            if (abs(FIELD(i,j) - EXP_FIELD(i,j)) > tol) then
                print *, 'FIELD Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_FIELD(i,j), &
                         ' but got ', FIELD(i,j)
                res = 1
            end if
        end do
    end do  

    if (res .ne. 0) stop 10
    
    ! Test for SMOOTHC()

    HOLD = 0.0 ! Reset HOLD for SMOOTHC test, just in case.

    call SMOOTHC(FIELDC, HOLD, IX, IY, SMTH)

    res = 0
    do i = 1, nx
        do j = 1, ny
            if (abs(FIELDC(i,j) - EXP_FIELDC(i,j)) > tol) then
                print *, 'FIELDC Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_FIELDC(i,j), &
                         ' but got ', FIELDC(i,j)
                res = 1
            end if
        end do
    end do  

    if (res .ne. 0) stop 20
    
    print *, "SUCCESS!"
end program test_smooth