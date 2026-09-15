! This is a test program for UPP.
!
! This program tests routines from the IcingSeverity module.
!
! Alyson Stahl, 9/2026
program test_icing_severity
    use IcingSeverity, only: icing_sev
    use CloudLayers, only: clouds_t
    implicit none

    real, parameter :: tol = 1.0e-6
    integer :: res
    
    res = 0

    ! Test Case 1: Convection scenario.
    call test_convection_scenario(res)
    if (res .ne. 0) stop 10

    ! Test Case 2: Precipitation below warmnose scenario.
    call test_precip_below_warmnose(res)
    if (res .ne. 0) stop 20
    
    ! Test Case 3: Precipitation above warmnose scenario.
    call test_precip_above_warmnose(res)
    if (res .ne. 0) stop 30
    print *, "SUCCESS!"

contains
    subroutine test_convection_scenario(res)
        integer, intent(out) :: res
        integer, parameter :: nz = 4
        integer :: imp_physics, prcpType, i
        real :: hgt(nz), rh(nz), t(nz), pres(nz), vv(nz)
        real :: liqCond(nz), iceCond(nz), twp(nz), ice_pot(nz)
        real :: iseverity(nz), expected(nz)
        real :: hcprcp, cape, lx, kx, tott, pc
        type(clouds_t) :: clouds

        res = 0

        imp_physics = 98
        hgt = (/ 3000.0, 2200.0, 1400.0, 800.0 /)
        rh = (/ 92.0, 88.0, 85.0, 82.0 /)
        t = (/ 243.15, 254.15, 269.15, 271.15 /)
        pres = (/ 70000.0, 76000.0, 82000.0, 88000.0 /)
        vv = (/ -0.20, -0.30, -0.10, -0.05 /)
        liqCond = (/ 0.08, 0.10, 0.12, 0.09 /)
        iceCond = (/ 0.04, 0.03, 0.02, 0.01 /)
        twp = (/ 150.0, 220.0, 300.0, 180.0 /)
        ice_pot = (/ 0.40, 0.55, 0.70, 0.0 /)

        hcprcp = 2.0
        cape = 1750.0
        lx = -5.0
        kx = 30.0
        tott = 40.0
        pc = 0.10
        prcpType = 4

        allocate(clouds%layerQ(nz))
        clouds%nLayers = 1
        clouds%wmnIdx = -1
        clouds%avv = -0.25
        clouds%layerQ = (/ 0.20, 0.35, 0.50, 0.10 /)
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        clouds%topIdx(1) = 1
        clouds%baseIdx(1) = 4
        clouds%ctt(1) = 258.15

        expected(1) = 0.41156417
        expected(2) = 0.54625165
        expected(3) = 0.60204083
        expected(4) = 0.0

        call icing_sev(imp_physics, hgt, rh, t, pres, vv, liqCond, iceCond, twp, &
             ice_pot, nz, hcprcp, cape, lx, kx, tott, pc, prcpType, clouds, iseverity)

        do i = 1, nz
            if (abs(iseverity(i) - expected(i)) > tol) then
                print *, "icing_sev() failed for test case ", i, ": expected ", &
                      expected(i), " but got ", iseverity(i)
                res = 1
            end if
        end do

        deallocate(clouds%layerQ)
    end subroutine test_convection_scenario

    subroutine test_precip_below_warmnose(res)
        integer, intent(out) :: res
        integer, parameter :: nz = 4
        integer :: imp_physics, prcpType, i
        real :: hgt(nz), rh(nz), t(nz), pres(nz), vv(nz)
        real :: liqCond(nz), iceCond(nz), twp(nz), ice_pot(nz)
        real :: iseverity(nz), expected(nz)
        real :: hcprcp, cape, lx, kx, tott, pc
        type(clouds_t) :: clouds

        res = 0

        imp_physics = 98
        hgt = (/ 9000.0, 5000.0, 1000.0, 0.0 /)
        rh = (/ 100.0, 100.0, 100.0, 100.0 /)
        t = (/ 270.0, 270.0, 270.0, 270.0 /)
        pres = (/ 70000.0, 75000.0, 80000.0, 85000.0 /)
        vv = (/ -0.10, -0.50, -0.10, -0.10 /)
        liqCond = (/ 0.15, 0.15, 0.15, 0.15 /)
        iceCond = (/ 0.15, 0.15, 0.15, 0.15 /)
        twp = (/ 200.0, 1000.0, 200.0, 200.0 /)
        ice_pot = (/ 0.0, 0.5, 0.0, 0.0 /)

        hcprcp = 0.5
        cape = 500.0
        lx = 2.0
        kx = 15.0
        tott = 15.0
        pc = 0.20
        prcpType = 1

        allocate(clouds%layerQ(nz))
        clouds%nLayers = 1
        clouds%wmnIdx = 2
        clouds%avv = -0.25
        clouds%layerQ = (/ 0.10, 0.10, 0.10, 0.10 /)
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        clouds%topIdx(1) = 1
        clouds%baseIdx(1) = 4
        clouds%ctt(1) = 223.15

        expected = (/ 0.0, 0.9375, 0.0, 0.0 /)

        call icing_sev(imp_physics, hgt, rh, t, pres, vv, liqCond, iceCond, twp, &
             ice_pot, nz, hcprcp, cape, lx, kx, tott, pc, prcpType, clouds, iseverity)

        do i = 1, nz
            if (abs(iseverity(i) - expected(i)) > tol) then
                print *, "icing_sev() failed for below warmnose test case ", i, &
                      ": expected ", expected(i), " but got ", iseverity(i)
                res = 1
            end if
        end do

        deallocate(clouds%layerQ)
    end subroutine test_precip_below_warmnose

    subroutine test_precip_above_warmnose(res)
        integer, intent(out) :: res
        integer, parameter :: nz = 4
        integer :: imp_physics, prcpType, i
        real :: hgt(nz), rh(nz), t(nz), pres(nz), vv(nz)
        real :: liqCond(nz), iceCond(nz), twp(nz), ice_pot(nz)
        real :: iseverity(nz), expected(nz)
        real :: hcprcp, cape, lx, kx, tott, pc
        type(clouds_t) :: clouds

        res = 0

        imp_physics = 98
        hgt = (/ 10000.0, 7000.0, 3000.0, 0.0 /)
        rh = (/ 100.0, 100.0, 100.0, 100.0 /)
        t = (/ 270.0, 270.0, 270.0, 270.0 /)
        pres = (/ 70000.0, 75000.0, 80000.0, 85000.0 /)
        vv = (/ -0.10, -0.10, -0.50, -0.10 /)
        liqCond = (/ 0.15, 0.15, 0.15, 0.15 /)
        iceCond = (/ 0.15, 0.15, 0.15, 0.15 /)
        twp = (/ 200.0, 200.0, 500.0, 200.0 /)
        ice_pot = (/ 0.0, 0.0, 0.5, 0.0 /)

        hcprcp = 0.5
        cape = 500.0
        lx = 2.0
        kx = 15.0
        tott = 15.0
        pc = 0.20
        prcpType = 1

        allocate(clouds%layerQ(nz))
        clouds%nLayers = 1
        clouds%wmnIdx = 1
        clouds%avv = -0.25
        clouds%layerQ = (/ 0.10, 0.10, 0.10, 0.10 /)
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        clouds%topIdx(1) = 2
        clouds%baseIdx(1) = 4
        clouds%ctt(1) = 223.15

        expected = (/ 0.0, 0.0, 0.8, 0.0 /)

        call icing_sev(imp_physics, hgt, rh, t, pres, vv, liqCond, iceCond, twp, &
             ice_pot, nz, hcprcp, cape, lx, kx, tott, pc, prcpType, clouds, iseverity)

        do i = 1, nz
            if (abs(iseverity(i) - expected(i)) > tol) then
                print *, "icing_sev() failed for above warmnose test case ", i, &
                      ": expected ", expected(i), " but got ", iseverity(i)
                res = 1
            end if
        end do

        deallocate(clouds%layerQ)
    end subroutine test_precip_above_warmnose


end program test_icing_severity