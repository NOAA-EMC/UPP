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
    ! Used to initialize inputs
    real, parameter :: zsfc = 10.0, dz = 200.0, Hscale = 8000.0
    integer, parameter :: npts = 2, nlevs = 60
    ! Used to calculate expected results
    real :: evp, rmx, rkapa, tlcl, dlplcl, dalp
    integer :: i, j, k, res
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

    fis = zsfc/gi        ! Set surface geopotential consistent with ~10 m AGL
    P1D = 100000.0       ! 1000 hPa surface parcel pressure (Pa)
    T1D = 300.0          ! 300 K surface parcel temperature
    Q1D = 0.012          ! 12 g/kg specific humidity (kg/kg)

    do k = 1, nlevs
        zint(:,:,k) = zsfc + (k-1)*dz
        alpint(:,:,k) = LOG(100000.0) - zint(:,:,k)/Hscale
    end do

    evp = 100000.0 * 0.012 / (eps + 0.012 * oneps)  
    rmx = eps * evp / (100000.0 - evp)
    rkapa = 1.0 / (D2845 * (1.0 - D28 * rmx))
    tlcl = H55 + H2840 / (D35*LOG(300.0)-LOG(evp * D01)-D4805)
    EXP_PLCL(1,1) = 100000.0 * (tlcl/300.0)**rkapa
    dlplcl = LOG(EXP_PLCL(1,1)) - alpint(1,1,nlevs)
    dalp = alpint(1,1,nlevs-1) - alpint(1,1,nlevs)
    EXP_ZLCL(1,1) = zint(1,1,nlevs) + dz*dlplcl/dalp - zsfc

    call callcl(P1D, T1D, Q1D, PLCL, ZLCL)

    if (abs(PLCL(1,1) - EXP_PLCL(1,1)) > tol) then
        print *, 'PLCL Test failed: Expected ', EXP_PLCL(1,1), &
                 ' but got ', PLCL(1,1)
        res = 1
    else
        print *, 'PLCL Test passed: ', PLCL(1,1)
    end if

    if (abs(ZLCL(1,1) - EXP_ZLCL(1,1)) > tol) then
        print *, 'ZLCL Test failed: Expected ', EXP_ZLCL(1,1), &
                 ' but got ', ZLCL(1,1)
        res = 1
    else
        print *, 'ZLCL Test passed: ', ZLCL(1,1)
    end if
    print *, 'SUCCESS!'
end program test_callcl