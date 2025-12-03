! This is a test program for UPP.
!
! This program tests the BOUND() subroutine.
!
! Alyson Stahl, 11/2025
program test_bound
    use ctlblk_mod, only: jsta, jend, spval, ista, iend, im, jm
    implicit none
    
    integer, parameter :: npts = 3
    real, parameter :: FMIN = 10.0, FMAX = 20.0
    integer :: i, j, res
    real :: FLD(1:npts, 1:npts), EXPECTED_FLD(1:npts, 1:npts)

    im = npts
    jm = npts
    ista = 1
    iend = npts
    jsta = 1
    jend = npts
    spval = -9999.0

    ! FLD contains values below FMIN, above FMAX, within range, and spval
    FLD = reshape([spval, 5.0, 15.0, 25.0, spval, 10.0, 20.0, 30.0, spval], [im, jm])
    EXPECTED_FLD = reshape([spval, FMIN, 15.0, FMAX, spval, FMIN, FMAX, FMAX, spval], [im, jm])

    print *, "Testing BOUND subroutine..."
    call BOUND(FLD, FMIN, FMAX)

    res = 0
    do j = 1, jm
        do i = 1, im
            if (FLD(i,j) /= EXPECTED_FLD(i,j)) then
                print *, 'FAIL: FLD(',i,',',j,') = ', FLD(i,j), ' expected ', EXPECTED_FLD(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    print *, "SUCCESS!"
end program test_bound
