! This is a test program for UPP.
!
! This program tests the CALLCL() subroutine.
!
! Alyson Stahl, 12/2025
program test_callcl
    use vrbls3d, only: alpint, zint
    use vrbls2d, only: fis
    use masks, only: lmh
    use params_mod, only: eps, oneps, d01, h1m12, gi, d00
    use ctlblk_mod, only: jsta, jend, spval, jsta_m, jend_m, im, &
                        ista, iend, ista_m, iend_m
    implicit none

    real, parameter :: tol = 1.0e-6, p0 = 101325.0, height = 8000.0, ztop = 15000.0
    integer, parameter :: npts = 2, nlevs = 60
    integer :: i, j, k, res
    real :: zsfc(1:npts,1:npts), levels(1:nlevs), pint(1:npts,1:npts,1:nlevs)
    real :: P1D(1:npts,1:npts), T1D(1:npts,1:npts), Q1D(1:npts,1:npts)
    real :: PLCL(1:npts,1:npts), ZLCL(1:npts,1:npts)
    real :: EXP_PLCL(1:npts,1:npts), EXP_ZLCL(1:npts,1:npts)

    ! Grid parameters
    jsta = 1
    jend = npts
    ista = 1
    iend = npts
    im = npts
    jsta_m = 1
    jend_m = npts
    ista_m = 1
    iend_m = npts
    
    spval = 9.9e10

    ! Allocate arrays
    allocate(alpint(1:npts,1:npts,1:nlevs))
    allocate(zint(1:npts,1:npts,1:nlevs))
    allocate(fis(1:npts,1:npts))
    allocate(lmh(1:npts,1:npts))

    ! Expected results
    EXP_ZLCL = reshape([486.8708, 9.9e+10, 428.96264648, 385.38275146], [npts, npts])
    EXP_PLCL = reshape([92988.375, 9.9e+10, 93080.39, 93005.65], [npts, npts])

    ! Initialize input arrays
    P1D = reshape([100000.0, spval, 99500.0, 99000.0], [npts, npts])
    T1D = reshape([292.0, 293.0, 293.0, 294.0], [npts, npts])
    Q1D = reshape([0.01, 0.011, 0.011, 0.012], [npts, npts])

    lmh = nlevs
    fis = reshape([1962.0, 2452.5, 2452.5, 2943.0], [npts, npts])
    zsfc = reshape([200.0, 250.0, 250.0, 300.0], [npts, npts])

    ! Calculate some reasonable values for alpint and zint
    do k = 1, nlevs
        levels(k) = real(k - 1) / real(nlevs - 1)
    end do

    do j = jsta, jend
        do i = ista, iend
            do k = 1, nlevs
                zint(i,j,k) = (1.0 - levels(k)) * ztop + levels(k) * zsfc(i,j)
                pint(i,j,k) = p0 * exp(-zint(i,j,k) / height)
                alpint(i,j,k) = log(pint(i,j,k))
            end do
        end do
    end do

    call callcl(P1D, T1D, Q1D, PLCL, ZLCL)

    ! Check results
    res = 0
    do j = jsta, jend
        do i = ista, iend
            if (abs(PLCL(i,j) - EXP_PLCL(i,j)) > tol) then
                print *, 'PLCL Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_PLCL(i,j), &
                         ' but got ', PLCL(i,j)
                res = 1
            end if
        end do
    end do
    do j = jsta, jend
        do i = ista, iend
            if (abs(ZLCL(i,j) - EXP_ZLCL(i,j)) > tol) then
                print *, 'ZLCL Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_ZLCL(i,j), &
                         ' but got ', ZLCL(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10
    
    print *, 'SUCCESS!'
end program test_callcl