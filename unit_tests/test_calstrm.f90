! This is a test program for UPP.
!
! This program tests the CALSTRM() subroutine.
!
! Alyson Stahl, 12/2025
program test_calstrm
    use ctlblk_mod, only: jsta, jend, im, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: npts = 2
    integer :: i, j, res
    real :: Z1D(1:npts,1:npts), STRM(1:npts,1:npts)
    real :: EXP_STRM(1:npts,1:npts)

    interface
        subroutine CALSTRM(Z1D, STRM)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, dimension(ista:iend,jsta:jend), intent(in) :: Z1D
            real, dimension(ista:iend,jsta:jend), intent(inout) :: STRM
        end subroutine CALSTRM
    end interface
    ! Grid parameters
    ista = 1
    iend = npts
    jsta = 1
    jend = npts
    im = 4

    Z1D(1,1) = 500.0    ! Normal height test
    Z1D(1,2) = 1000.0   ! Another level
    Z1D(2,1) = -500.0   ! Negative height test
    Z1D(2,2) = 0.0      ! Zero height test
    STRM = 0.0
    EXP_STRM = reshape([5.24657E7, -5.24657E7, 1.049314E8, 0.0], [npts, npts])

    call CALSTRM(Z1D, STRM)

    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(STRM(i,j) - EXP_STRM(i,j)) > tol) then
                print *, 'Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_STRM(i,j), &
                         ' but got ', STRM(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10
    
    print *, 'SUCCESS!'
end program test_calstrm