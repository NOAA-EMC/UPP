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
                print *, 'Test failed at (', i, ',', j, '): ', &
                    'Expected VAR=1.0 but got VAR=', VAR(i,j)
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
                print *, 'Test failed at (', i, ',', j, '): ', &
                    'Expected VAR=1.0 but got VAR=', VAR(i,j)
                res = 1
            end if
        end do
    end do
    if (res .ne. 0) stop 20

    ! Test Case: COSL < SMALL at both poles, but VAR is SPVAL everywhere. VAR should be unchanged.
    COSL = SMALL / 10.0
    VAR = SPVAL
    call POLEAVG(IM, JM, JSTA, JEND, SMALL, COSL, SPVAL, VAR)
    
    res = 0
    do i = 1, nx
        do j = 1, ny
            if (abs(VAR(i,j) - SPVAL) > tol) then
                print *, 'Test failed at (', i, ',', j, '): ', &
                    'Expected VAR=', SPVAL, ' but got VAR=', VAR(i,j)
                res = 1
            end if
        end do
    end do
    if (res .ne. 0) stop 30

    ! Test Case: COSL < SMALL at both poles and VAR does not contain SPVAL. 
    ! VAR should be set to the average of all values in the column.
    do i = 1, nx
        do j = 1, ny
            VAR(i,j) = real((j-1)*nx + i)
        end do
    end do

    ! Average of both poles
    EXP_VAR(:,1) = 2.0
    EXP_VAR(1,2) = VAR(1,2)
    EXP_VAR(2,2) = VAR(2,2)
    EXP_VAR(3,2) = VAR(3,2) 
    EXP_VAR(:,3) = 8.0

    call POLEAVG(IM, JM, JSTA, JEND, SMALL, COSL, SPVAL, VAR)
    res = 0
    do i = 1, nx
        do j = 1, ny
            if (abs(VAR(i,j) - EXP_VAR(i,j)) > tol) then
                print *, 'Test failed at (', i, ',', j, '): ', &
                    'Expected VAR=', EXP_VAR(i,j), ' but got VAR=', VAR(i,j)
                res = 1
            end if
        end do
    end do
    if (res .ne. 0) stop 40
    
    ! Test Case: Mix of SPVAL and non-SPVAL values in VAR. VAR should be set to the 
    ! average of the non-SPVAL values in the column.
    VAR(1,1) = 1.0
    VAR(2,1) = SPVAL
    VAR(3,1) = 3.0
    VAR(1,3) = SPVAL
    VAR(2,3) = 8.0
    VAR(3,3) = 9.0

    EXP_VAR(:,1) = 2.0
    EXP_VAR(:,3) = 8.5

    call POLEAVG(IM, JM, JSTA, JEND, SMALL, COSL, SPVAL, VAR)
    res = 0
    do i = 1, nx
        do j = 1, ny
            if (abs(VAR(i,j) - EXP_VAR(i,j)) > tol) then
                print *, 'Test failed at (', i, ',', j, '): ', &
                    'Expected VAR=', EXP_VAR(i,j), ' but got VAR=', VAR(i,j)
                res = 1
            end if
        end do
    end do
    if (res .ne. 0) stop 50
    
    print *, "SUCCESS!"
end program test_poleavg