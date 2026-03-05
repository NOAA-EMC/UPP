! This is a test program for UPP.
!
! This program tests the CALVIS_GSD() subroutine.
!
! Alyson Stahl, 12/2025
program test_calvis_gsd
    use vrbls2d, only: sno, si, ustar, z0
    use vrbls3d, only: qqw, qqi, qqs, qqr, qqg, t, pmid, q, u, v, aextc55
    use ctlblk_mod, only: jm, im, jsta_2l, jend_2u, lm, modelname, spval, method_blsn, &
                          ista_2l, iend_2u
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 4, nlevs = 5
    integer :: i, j, k, res
    real, dimension(1:npts,1:npts) :: CZEN, VIS, EXP_VIS

    interface
        subroutine CALVIS_GSD(CZEN, VIS)
            use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u, &
                                  ista, iend, ista_2l, iend_2u
            real, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in) :: CZEN
            real, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(out) :: VIS
        end subroutine CALVIS_GSD
    end interface

    ! Grid parameters
    jm = npts
    im = npts
    jsta_2l = 1
    jend_2u = npts
    ista_2l = 1
    iend_2u = npts
    lm = nlevs

    spval = 9.9e10
    modelname = 'GFS'
    method_blsn = .true.

    allocate(sno(1:npts,1:npts))
    allocate(si(1:npts,1:npts))
    allocate(ustar(1:npts,1:npts))
    allocate(z0(1:npts,1:npts))
    allocate(qqw(1:npts,1:npts,1:nlevs))
    allocate(qqi(1:npts,1:npts,1:nlevs))
    allocate(qqs(1:npts,1:npts,1:nlevs))
    allocate(qqr(1:npts,1:npts,1:nlevs))
    allocate(qqg(1:npts,1:npts,1:nlevs))
    allocate(t(1:npts,1:npts,1:nlevs))
    allocate(pmid(1:npts,1:npts,1:nlevs))
    allocate(q(1:npts,1:npts,1:nlevs))
    allocate(u(1:npts,1:npts,1:nlevs))
    allocate(v(1:npts,1:npts,1:nlevs))
    allocate(aextc55(1:npts,1:npts,1:nlevs))

    ! Initialize arrays with realistic base values
    sno = 0.05         ! m SWE
    si  = 50.0         ! mm snow depth
    ustar = 0.30       ! m s^-1
    z0 = 0.10          ! m, roughness length

    do k = 1, nlevs
        qqw(:,:,k) = 2.0e-5   ! kg/kg cloud water
        qqi(:,:,k) = 2.0e-5   ! kg/kg ice
        qqs(:,:,k) = 2.0e-5   ! kg/kg snow
        qqr(:,:,k) = 1.0e-5   ! kg/kg rain
        qqg(:,:,k) = 1.0e-5   ! kg/kg graupel

        t(:,:,k)    = 270.0          ! K
        pmid(:,:,k) = 100000.0 - 2000.0*real(k-1)  ! Pa
        q(:,:,k)    = 0.004          ! kg/kg specific humidity

        u(:,:,k) = 2.0 + 1.0*real(k-1)  ! m/s, increasing with height
        v(:,:,k) = 1.0 + 0.5*real(k-1)  ! m/s, increasing with height

        aextc55(:,:,k) = 1.0e-4       ! m^-1 aerosol extinction
    end do

    CZEN = 0.5  ! daytime default

    ! (1,1): Clear-air case, no hydrometeors or aerosols; low humidity for large vis
    do k = 1, 3
        qqw(1,1,k) = 0.0; qqi(1,1,k) = 0.0; qqs(1,1,k) = 0.0
        qqr(1,1,k) = 0.0; qqg(1,1,k) = 0.0
        aextc55(1,1,k) = 0.0
    end do
    q(1,1,1) = 0.001
    CZEN(1,1) = 0.5
    z0(1,1) = 0.10

    ! (1,2): Heavy dry snow + graupel at night; strong ustar to trigger BLSN
    do k = 1, 3
        qqs(1,2,k) = 3.0e-4
        qqg(1,2,k) = 2.0e-4
        qqw(1,2,k) = 1.0e-5
        qqi(1,2,k) = 1.0e-5
        qqr(1,2,k) = 0.0
        t(1,2,k)   = 260.0
    end do
    ustar(1,2) = 0.60
    sno(1,2)   = 0.20
    CZEN(1,2)  = 0.05   ! night/low sun

    ! (1,3): Wet snow near freezing; forest (z0>0.7) disables BLSN even if enabled
    do k = 1, 3
        qqs(1,3,k) = 4.0e-4
        qqr(1,3,k) = 5.0e-5
        t(1,3,k)   = 273.5
    end do
    z0(1,3)    = 1.0
    ustar(1,3) = 0.70
    CZEN(1,3)  = 0.30

    ! (2,1): Aerosol-only attenuation (hydrometeors zero, elevated aerosols)
    do k = 1, 3
        qqw(2,1,k) = 0.0; qqi(2,1,k) = 0.0; qqs(2,1,k) = 0.0
        qqr(2,1,k) = 0.0; qqg(2,1,k) = 0.0
        aextc55(2,1,k) = 5.0e-4
    end do
    q(2,1,1) = 0.0035
    CZEN(2,1) = 0.60

    ! (2,2): Strong low-level shear between levels 1 and 4
    u(2,2,1) = 0.0; v(2,2,1) = 0.0
    u(2,2,4) = 12.0; v(2,2,4) = 8.0
    CZEN(2,2) = 0.20

    ! (2,3): High RH clear-air (warm and moist), hydrometeors zero, no aerosols
    do k = 1, 2
        qqw(2,3,k) = 0.0; qqi(2,3,k) = 0.0; qqs(2,3,k) = 0.0
        qqr(2,3,k) = 0.0; qqg(2,3,k) = 0.0
        t(2,3,k)   = 285.0
        aextc55(2,3,k) = 0.0
        q(2,3,k)   = 0.018
    end do
    CZEN(2,3) = 0.50

    ! (3,1): Very low sun angle (night-like)
    CZEN(3,1) = 0.001

    ! Set some array values to spval to test the handling of spval
    t(4,1,lm) = spval
    u(4,2,lm)    = spval
    v(4,3,lm)    = spval
    pmid(4,4,lm) = spval

    EXP_VIS = reshape([ &
        204.72335815, 204.72335815, 459.18145752, spval, &
        147.99403381, 204.72335815, 204.72335815, spval, &
        163.81585693, 204.72335815, 204.72335815, spval, &
        204.72335815, 204.72335815, 204.72335815, spval  &
    ], [npts, npts])

    call CALVIS_GSD(CZEN, VIS)

    ! Deallocate all allocated arrays
    deallocate(sno)
    deallocate(si)
    deallocate(ustar)
    deallocate(z0)
    deallocate(qqw)
    deallocate(qqi)
    deallocate(qqs)
    deallocate(qqr)
    deallocate(qqg)
    deallocate(t)
    deallocate(pmid)
    deallocate(q)
    deallocate(u)
    deallocate(v)
    deallocate(aextc55)
    
    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(VIS(i,j) - EXP_VIS(i,j)) > tol) then
                print *, 'VIS Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_VIS(i,j), &
                         ' but got ', VIS(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    print *, 'SUCCESS!'
end program test_calvis_gsd