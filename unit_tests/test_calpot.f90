! This is a test program for UPP.
!
! This program tests the CALPOT() subroutine.
!
! Alyson Stahl, 12/2025
program test_calpot
    use ctlblk_mod, only: jsta, jend, spval, im,  ista, iend
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: npts = 2
    integer :: i, j, res
    real :: P1D(1:npts,1:npts), T1D(1:npts,1:npts)
    real :: THETA(1:npts,1:npts), EXP_THETA(1:npts,1:npts)
    
    interface
        subroutine CALPOT(P1D, T1D, THETA)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, dimension(ista:iend,jsta:jend), intent(in) :: P1D, T1D
            real, dimension(ista:iend,jsta:jend), intent(inout) :: THETA
        end subroutine CALPOT
    end interface

    ! Grid parameters
    jsta = 1
    jend = npts
    ista = 1
    iend = npts
    im = npts
    spval = 9999.0

    P1D = reshape([100000.0, 1.0, 99500.0, 99500.0], [npts, npts])
    T1D = reshape([275.0, 275.0, 300.0, spval], [npts, npts])
    EXP_THETA = reshape([275.0, 0.0, 300.4302368, spval], [npts, npts])

    call CALPOT(P1D, T1D, THETA)

    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(THETA(i,j) - EXP_THETA(i,j)) > tol) then
                print *, 'Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_THETA(i,j), &
                         ' but got ', THETA(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    print *, 'SUCCESS!'
end program test_calpot