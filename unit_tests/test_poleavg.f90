! This is a test program for UPP.
!
! This program tests the POLEAVG() subroutine.
!
! Alyson Stahl, 2/2026
program test_poleavg
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: nx = 3, ny = 3
    integer :: i, j, res
    integer :: IM, JM, JSTA, JEND
    real :: SMALL, SPVAL
    real :: COSL(nx, ny), VAR(nx, ny), EXP_VAR(nx, ny)

    interface 
        subroutine POLEAVG(IM, JM, JSTA, JEND, SMALL, COSL, SPVAL, VAR)
            integer, intent(in) :: IM, JM, JSTA, JEND
            real, intent(in) :: SMALL, SPVAL
            real, dimension(IM,JSTA:JEND), intent(in) :: COSL
            real, dimension(IM,JSTA:JEND), intent(inout) :: VAR
        end subroutine POLEAVG
    end interface

    SMALL = 1.0e-6
    SPVAL = 9.9e10
    IM = nx
    JM = ny
    COSL = SMALL * 10.0 
    VAR = 1.0

    ! Test Case: jsta > 1 and jend < jm. VAR should be unchanged.
    JSTA = 2
    JEND = 2

    call POLEAVG(IM, JM, JSTA, JEND, SMALL, COSL, SPVAL, VAR)

    res = 0
    do i = 1, nx
        do j = 1, ny
            if (abs(VAR(i,j) - 1.0) > tol) then
                print *, 'Test failed at (', i, ',', j, '): Expected VAR=1.0 but got VAR=', VAR(i,j)
                res = 1
            end if
        end do
    end do
    if (res .ne. 0) stop 10

    ! Test Case: JJ in bounds, but COSL > SMALL at both poles. VAR should be unchanged.
    JSTA = 1
    JEND = ny
    call POLEAVG(IM, JM, JSTA, JEND, SMALL, COSL, SPVAL, VAR)
    
    res = 0
    do i = 1, nx
        do j = 1, ny
            if (abs(VAR(i,j) - 1.0) > tol) then
                print *, 'Test failed at (', i, ',', j, '): Expected VAR=1.0 but got VAR=', VAR(i,j)
                res = 1
            end if
        end do
    end do
    if (res .ne. 0) stop 20

    print *, "SUCCESS!"
end program test_poleavg