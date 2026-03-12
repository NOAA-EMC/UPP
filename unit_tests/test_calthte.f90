! This is a test program for UPP.
!
! This program tests the CALTHTE() subroutine.
!
! Alyson Stahl, 12/2025
program test_calthte
    use ctlblk_mod, only: jsta, jend, im, spval, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: npts = 2
    integer :: i, j, res
    real, dimension(1:npts,1:npts) :: P1D, T1D, Q1D
    real, dimension(1:npts,1:npts) :: THTE, EXP_THTE

    interface
        subroutine CALTHTE(P1D, T1D, Q1D, THTE)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, dimension(ista:iend,jsta:jend), intent(in) :: P1D, T1D, Q1D
            real, dimension(ista:iend,jsta:jend), intent(out) :: THTE
        end subroutine CALTHTE
    end interface
    ! Grid parameters
    jsta = 1
    jend = npts
    im = npts
    ista = 1
    iend = npts

    spval = 9.9e10

    P1D = reshape([100000.0, spval, 99500.0, 99000.0], [npts, npts])
    T1D = reshape([292.0, 293.0, 293.0, 294.0], [npts, npts])
    Q1D = reshape([0.01, 0.011, 0.011, 0.012], [npts, npts])
    EXP_THTE = reshape([320.884246826171875, 0.0, 325.338592529296875, 329.84857177734375], [npts, npts])

    call CALTHTE(P1D, T1D, Q1D, THTE)

    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(THTE(i,j) - EXP_THTE(i,j)) > tol) then
                print *, 'THTE Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_THTE(i,j), &
                         ' but got ', THTE(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    print *, 'SUCCESS!'
end program test_calthte