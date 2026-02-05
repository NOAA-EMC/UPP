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

    real, parameter :: tol = 1.0e-6
    ! From CALLCL.f
    real, parameter :: D35=3.5, D4805=4.805,  H2840=2840., H55=55., D2845=0.2845, D28=0.28
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

    ! CALLCL() uses lmh to determine the number of levels
    ! Set to nlevs - 1 since loop accesses lmh(i,j) + 1
    lmh = nlevs - 1

    fis = 980.0      ! Geopotential at surface (~100 m AGL if GI≈1/g)
    P1D = 95000.0    ! Surface parcel pressure (Pa), ~950 hPa
    T1D = 300.0      ! Parcel temperature (K), ~27°C
    Q1D = 0.012      ! Specific humidity (kg/kg), ~12 g/kg

    

    call callcl(P1D, T1D, Q1D, PLCL, ZLCL)

    
    print *, 'SUCCESS!'
end program test_callcl