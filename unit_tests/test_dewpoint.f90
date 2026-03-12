! This is a test program for UPP.
!
! This program tests the DEWPOINT() subroutine.
!
! Alyson Stahl, 11/2025
program test_dewpoint
    use ctlblk_mod, only: jsta, jend, spval, ista, iend, im, jm
    implicit none
    
    integer, parameter :: npts = 2

    ! From DEWPOINT documentation, plausible physical range for dewpoint (K)
    real, parameter :: min_expected = 233.0, max_expected = 315.0
    integer :: i, j, res
    real :: VP(1:npts, 1:npts), TD(1:npts, 1:npts)

    interface
        subroutine DEWPOINT(VP, TD)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real,intent(in) ::  VP(ista:iend,jsta:jend)
            real,intent(out) :: TD(ista:iend,jsta:jend)
        end subroutine DEWPOINT
    end interface

    im = npts
    jm = npts
    ista = 1
    iend = npts
    jsta = 1
    jend = npts
    spval = -9999.0
    
    ! Set up test vapor pressure values (centibars)
    VP = reshape([0.02, 0.5, 2.0, 8.0], [npts, npts])

    print *, 'Testing DEWPOINT subroutine with spval less than min vapor pressure. Expect TD = spval.'

    call DEWPOINT(VP, TD)
    
    res = 0
    do i = 1, npts
        do j = 1, npts
            if (TD(i,j) /= spval) then
                print *, 'FAIL: TD(',i,',',j,') = ', TD(i,j), ' expected spval = ', spval
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    print *, "Testing DEWPOINT subroutine with spval greater than max vapor pressure."

    spval = 9.9e10
    call DEWPOINT(VP, TD)

    res = 0
    do i = 1, npts
        do j = 1, npts
            if (TD(i,j) < min_expected .or. TD(i,j) > max_expected) then
                print *, 'FAIL: TD(',i,',',j,') = ', TD(i,j), ' out of expected range.'
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 11

    print *, "SUCCESS!"
end program test_dewpoint