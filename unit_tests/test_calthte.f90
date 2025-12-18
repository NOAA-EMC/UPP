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
    real :: P1D(1:npts,1:npts), T1D(1:npts,1:npts), Q1D(1:npts,1:npts)
    real :: THTE(1:npts,1:npts), EXP_THTE(1:npts,1:npts)

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

    call CALTHTE(P1D, T1D, Q1D, THTE)

    do j = jsta, jend
        do i = ista, iend
            print *, 'THTE(', i, ',', j, ') = ', THTE(i,j)
        end do
    end do


end program test_calthte