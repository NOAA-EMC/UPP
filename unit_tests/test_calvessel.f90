! This is a test program for UPP.
!
! This program tests the CALVESSEL() subroutine.
!
! Alyson Stahl, 12/2025
program test_calvessel
    use vrbls2d, only: sst, u10h, v10h, tshltr, pshltr
    use masks, only: sm, sice
    use ctlblk_mod, only: jsta, jend, im, spval, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: npts = 3
    integer :: i, j, res
    real :: ICEG(1:npts,1:npts), EXP_ICEG(1:npts,1:npts)

    ! Grid parameters
    jsta = 1
    jend = npts
    im = npts
    ista = 1
    iend = npts

    spval = 9.9e10

    allocate(sm(1:npts,1:npts))
    allocate(sice(1:npts,1:npts))
    allocate(sst(1:npts,1:npts))
    allocate(u10h(1:npts,1:npts))
    allocate(v10h(1:npts,1:npts))
    allocate(tshltr(1:npts,1:npts))
    allocate(pshltr(1:npts,1:npts))

    sice = 0.0
    sm = 1.0
    ! Test Case: Fail mask check
    sm(1,1) = 0.5 

    sst = 280.0
    ! Test Case: sst > 285.15 K
    sst(1,2) = 286.0 

    u10h = 20.0
    v10h = 30.0

    ! Test Case: SQRT(U10H(I,J)**2+V10H(I,J)**2) > 50.0
    u10h(1,3) = 50.0 
    v10h(1,3) = 60.0

    pshltr = 1.1e5
    tshltr = 264
    ! Test Case: Shelter level temp > 273.15 K
    tshltr(2,1) = 273.15

    ! Test Case: ICEG goes below 0
    tshltr(2,2) = 265.0 

    call CALVESSEL(ICEG)

    do j = jsta, jend
        do i = ista, iend
            write(*,'(A,I0,A,I0,A,E30.20)') 'ICEG(', i, ',', j, ') = ', ICEG(i,j)
        end do
    end do


end program test_calvessel