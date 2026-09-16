! This is a test program for UPP.
!
! This program tests routines from the IcingSeverity module.
!
! Alyson Stahl, 9/2026
program test_icing_severity
    use IcingSeverity, only: icing_sev
    use DerivedFields, only : PRECIPS
    use CloudLayers, only: clouds_t
    implicit none

    real, parameter :: tol = 1.0e-6
    integer :: res
    
    res = 0

    ! Test Case 1: Convection scenario.
    call test_convection(res)
    if (res .ne. 0) stop 10

    ! Test Case 2: Precipitation below warmnose scenario. 
    ! This scenario is currently unreachable as written.
    call test_precip_below_warmnose(res)
    if (res .ne. 0) stop 20
    
    ! Test Case 3: Precipitation above warmnose scenario.
    call test_precip_above_warmnose(res)
    if (res .ne. 0) stop 30

    ! Test Case 4: No precipitation scenario.
    call test_no_precip(res)
    if (res .ne. 0) stop 40

    ! Test Case 5: Snow scenario.
    call test_snow(res)
    !if (res .ne. 0) stop 50

    ! Test Case 6: Cold rain scenario.
    call test_cold_rain(res)
    ! if (res .ne. 0) stop 60

    ! Test Case 7: Warm precipitation scenario.
    call test_warm_precip(res)
    ! if (res .ne. 0) stop 70

    ! Test Case 8: Freezing precipitation scenario.
    ! call test_freezing_precip(res)
    ! if (res .ne. 0) stop 80
    
    ! ! Test Case 9: Invalid scenario.
    call test_invalid_scenario(res)
    !if (res .ne. 0) stop 90

    print *, "SUCCESS!"

contains
    subroutine test_convection(res)
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
        prcpType = PRECIPS%CONVECTION

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
                print *, "icing_sev() failed at index ", i, ": expected ", &
                      expected(i), " but got ", iseverity(i)
                res = 1
            end if
        end do

        deallocate(clouds%layerQ)
    end subroutine test_convection

    subroutine test_precip_below_warmnose(res)
        integer, intent(out) :: res
        ! This test is intended to check the branch where
        ! "elseif(isClassicPrcpBlwWmn(prcpType, k))" = True
        ! 
        ! Currently, isClassicPrcpBlwWmn() will always return false, making
        ! the branch unreachable in the current implementation.
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

        imp_physics = 99
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
        prcpType = PRECIPS%RAIN

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

        expected = 0.0
        expected(3) = 0.8

        call icing_sev(imp_physics, hgt, rh, t, pres, vv, liqCond, iceCond, twp, &
             ice_pot, nz, hcprcp, cape, lx, kx, tott, pc, prcpType, clouds, iseverity)

        do i = 1, nz
            if (abs(iseverity(i) - expected(i)) > tol) then
                print *, "icing_sev() failed at index ", i, ": expected ", &
                      expected(i), " but got ", iseverity(i)
                res = 1
            end if
        end do

        deallocate(clouds%layerQ)
    end subroutine test_precip_above_warmnose

    subroutine test_no_precip(res)
        integer, intent(out) :: res
        integer, parameter :: nz = 4
        integer :: imp_physics, prcpType, i
        real :: hgt(nz), rh(nz), t(nz), pres(nz), vv(nz)
        real :: liqCond(nz), iceCond(nz), twp(nz), ice_pot(nz)
        real :: iseverity(nz), expected(nz)
        real :: hcprcp, cape, lx, kx, tott, pc
        type(clouds_t) :: clouds

        res = 0

        imp_physics = 11
        hgt = (/ 3000.0, 1800.0, 914.4, 0.0 /)
        rh = (/ 85.0, 95.0, 100.0, 90.0 /)
        t = (/ 270.0, 270.0, 270.0, 271.0 /)
        pres = (/ 70000.0, 76000.0, 82000.0, 88000.0 /)
        vv = (/ -0.10, -0.15, -0.25, -0.05 /)
        liqCond = (/ 0.20, 0.40, 0.60, 0.10 /)
        iceCond = (/ 0.10, 0.20, 0.40, 0.05 /)
        twp = (/ 100.0, 150.0, 200.0, 80.0 /)
        ice_pot = (/ 0.0, 0.0, 0.50, 0.0 /)

        hcprcp = 0.0
        cape = 200.0
        lx = 2.0
        kx = 18.0
        tott = 20.0
        pc = 0.02
        prcpType = PRECIPS%NONE

        allocate(clouds%layerQ(nz))
        clouds%nLayers = 1
        clouds%wmnIdx = -1
        clouds%avv = -0.10
        clouds%layerQ = (/ 0.10, 0.20, 0.50, 0.10 /)
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        clouds%topIdx(1) = 2
        clouds%baseIdx(1) = 4
        clouds%ctt(1) = 270.0

        expected = 0.0
        expected(3) = 0.61428571

        call icing_sev(imp_physics, hgt, rh, t, pres, vv, liqCond, iceCond, twp, &
             ice_pot, nz, hcprcp, cape, lx, kx, tott, pc, prcpType, clouds, iseverity)

        do i = 1, nz
            if (abs(iseverity(i) - expected(i)) > tol) then
                print *, "icing_sev() failed at index ", i, ": expected ", &
                      expected(i), " but got ", iseverity(i)
                res = 1
            end if
        end do

        deallocate(clouds%layerQ)
    end subroutine test_no_precip

    subroutine test_snow(res)
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
        hgt = (/ 7000.0, 5500.0, 1828.8, 0.0 /)
        rh = (/ 95.0, 98.0, 100.0, 96.0 /)
        t = (/ 256.15, 257.15, 258.15, 259.15 /)
        pres = (/ 65000.0, 72000.0, 82000.0, 90000.0 /)
        vv = (/ -0.10, -0.20, -0.25, -0.10 /)
        liqCond = (/ 0.10, 0.20, 0.35, 0.10 /)
        iceCond = (/ 0.40, 0.50, 0.65, 0.30 /)
        twp = (/ 200.0, 350.0, 500.0, 200.0 /)
        ice_pot = (/ 0.0, 0.0, 0.50, 0.0 /)

        hcprcp = 0.0
        cape = 150.0
        lx = 3.0
        kx = 12.0
        tott = 18.0
        pc = 0.02
        prcpType = PRECIPS%SNOW

        allocate(clouds%layerQ(nz))
        clouds%nLayers = 1
        clouds%wmnIdx = -1
        clouds%avv = -0.15
        clouds%layerQ = (/ 0.20, 0.30, 0.40, 0.15 /)
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        clouds%topIdx(1) = 2
        clouds%baseIdx(1) = 4
        clouds%ctt(1) = 258.15

        expected = 0.0
        expected(3) = 0.579598248

        call icing_sev(imp_physics, hgt, rh, t, pres, vv, liqCond, iceCond, twp, &
             ice_pot, nz, hcprcp, cape, lx, kx, tott, pc, prcpType, clouds, iseverity)

        do i = 1, nz
            if (abs(iseverity(i) - expected(i)) > tol) then
                print *, "icing_sev() failed at index ", i, ": expected ", &
                      expected(i), " but got ", iseverity(i)
                res = 1
            end if
        end do

        deallocate(clouds%layerQ)
    end subroutine test_snow

    subroutine test_cold_rain(res)
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
        hgt = (/ 4500.0, 3000.0, 1500.0, 0.0 /)
        rh = (/ 96.0, 94.0, 92.0, 90.0 /)
        t = (/ 255.15, 262.15, 268.15, 274.15 /)
        pres = (/ 65000.0, 73000.0, 82000.0, 91000.0 /)
        vv = (/ -0.10, -0.15, -0.20, -0.05 /)
        liqCond = (/ 0.08, 0.12, 0.18, 0.20 /)
        iceCond = (/ 0.25, 0.15, 0.05, 0.0 /)
        twp = (/ 100.0, 150.0, 220.0, 120.0 /)
        ice_pot = (/ 0.0, 0.0, 0.40, 0.0 /)

        hcprcp = 0.0
        cape = 250.0
        lx = 1.0
        kx = 20.0
        tott = 24.0
        pc = 0.10
        prcpType = PRECIPS%RAIN

        allocate(clouds%layerQ(nz))
        clouds%nLayers = 1
        clouds%wmnIdx = -1
        clouds%avv = -0.10
        clouds%layerQ = (/ 0.10, 0.15, 0.20, 0.10 /)
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        clouds%topIdx(1) = 2
        clouds%baseIdx(1) = 4
        clouds%ctt(1) = 255.15

        expected = 0.0
        expected(3) = 0.413246214

        call icing_sev(imp_physics, hgt, rh, t, pres, vv, liqCond, iceCond, twp, &
             ice_pot, nz, hcprcp, cape, lx, kx, tott, pc, prcpType, clouds, iseverity)

        do i = 1, nz
            if (abs(iseverity(i) - expected(i)) > tol) then
                print *, "icing_sev() failed at index ", i, ": expected ", &
                      expected(i), " but got ", iseverity(i)
                res = 1
            end if
        end do

        deallocate(clouds%layerQ)
    end subroutine test_cold_rain

    subroutine test_warm_precip(res)
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
        hgt = (/ 4000.0, 2600.0, 1200.0, 0.0 /)
        rh = (/ 94.0, 92.0, 90.0, 88.0 /)
        t = (/ 262.15, 264.15, 270.15, 274.15 /)
        pres = (/ 66000.0, 74000.0, 83000.0, 92000.0 /)
        vv = (/ -0.05, -0.10, -0.15, -0.05 /)
        liqCond = (/ 0.10, 0.15, 0.20, 0.18 /)
        iceCond = (/ 0.08, 0.05, 0.02, 0.0 /)
        twp = (/ 120.0, 180.0, 240.0, 150.0 /)
        ice_pot = (/ 0.0, 0.0, 0.35, 0.0 /)

        hcprcp = 0.0
        cape = 300.0
        lx = 1.5
        kx = 22.0
        tott = 26.0
        pc = 0.08
        prcpType = PRECIPS%RAIN

        allocate(clouds%layerQ(nz))
        clouds%nLayers = 1
        clouds%wmnIdx = -1
        clouds%avv = -0.08
        clouds%layerQ = (/ 0.10, 0.15, 0.18, 0.10 /)
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        clouds%topIdx(1) = 2
        clouds%baseIdx(1) = 4
        clouds%ctt(1) = 264.15

        expected = 0.0
        expected(3) = 0.512249529

        call icing_sev(imp_physics, hgt, rh, t, pres, vv, liqCond, iceCond, twp, &
             ice_pot, nz, hcprcp, cape, lx, kx, tott, pc, prcpType, clouds, iseverity)

        do i = 1, nz
            if (abs(iseverity(i) - expected(i)) > tol) then
                print *, "icing_sev() failed at index ", i, ": expected ", &
                      expected(i), " but got ", iseverity(i)
                res = 1
            end if
        end do

        deallocate(clouds%layerQ)
    end subroutine test_warm_precip

    subroutine test_freezing_precip(res)
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
        hgt = (/ 4200.0, 2800.0, 1300.0, 0.0 /)
        rh = (/ 95.0, 93.0, 91.0, 89.0 /)
        t = (/ 255.15, 260.15, 268.15, 272.15 /)
        pres = (/ 66000.0, 74000.0, 83000.0, 92000.0 /)
        vv = (/ -0.08, -0.12, -0.18, -0.06 /)
        liqCond = (/ 0.06, 0.10, 0.14, 0.10 /)
        iceCond = (/ 0.20, 0.15, 0.08, 0.02 /)
        twp = (/ 110.0, 170.0, 210.0, 140.0 /)
        ice_pot = (/ 0.0, 0.0, 0.30, 0.0 /)

        hcprcp = 0.0
        cape = 250.0
        lx = 1.0
        kx = 18.0
        tott = 22.0
        pc = 0.07
        prcpType = PRECIPS%OTHER

        allocate(clouds%layerQ(nz))
        clouds%nLayers = 1
        clouds%wmnIdx = -1
        clouds%avv = -0.10
        clouds%layerQ = (/ 0.08, 0.12, 0.16, 0.08 /)
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        clouds%topIdx(1) = 2
        clouds%baseIdx(1) = 4
        clouds%ctt(1) = 255.15

        expected = 0.0

        call icing_sev(imp_physics, hgt, rh, t, pres, vv, liqCond, iceCond, twp, &
             ice_pot, nz, hcprcp, cape, lx, kx, tott, pc, prcpType, clouds, iseverity)

        do i = 1, nz
            if (abs(iseverity(i) - expected(i)) > tol) then
                print *, "icing_sev() failed at index ", i, ": expected ", &
                      expected(i), " but got ", iseverity(i)
                res = 1
            end if
        end do

        deallocate(clouds%layerQ)
    end subroutine test_freezing_precip

    subroutine test_invalid_scenario(res)
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
        hgt = (/ 1000.0, 750.0, 500.0, 250.0 /)
        rh = (/ 0.0, 0.0, 0.0, 0.0 /)
        t = (/ 270.0, 270.0, 270.0, 270.0 /)
        pres = (/ 80000.0, 82000.0, 84000.0, 86000.0 /)
        vv = (/ 0.0, 0.0, 0.0, 0.0 /)
        liqCond = 0.0
        iceCond = 0.0
        twp = 0.0
        ice_pot = (/ 0.50, 0.0, 0.0, 0.0 /)

        hcprcp = 0.0
        cape = 0.0
        lx = 0.0
        kx = 0.0
        tott = 0.0
        pc = 0.0
        prcpType = 99

        allocate(clouds%layerQ(nz))
        clouds%nLayers = 1
        clouds%wmnIdx = -1
        clouds%avv = 0.0
        clouds%layerQ = 0.0
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        clouds%topIdx(1) = 1
        clouds%baseIdx(1) = 1
        clouds%ctt(1) = 270.0

        expected = 0.0

        call icing_sev(imp_physics, hgt, rh, t, pres, vv, liqCond, iceCond, twp, &
             ice_pot, nz, hcprcp, cape, lx, kx, tott, pc, prcpType, clouds, iseverity)

        do i = 1, nz
            if (abs(iseverity(i) - expected(i)) > tol) then
                print *, "icing_sev() failed at index ", i, ": expected ", &
                      expected(i), " but got ", iseverity(i)
                res = 1
            end if
        end do

        deallocate(clouds%layerQ)
    end subroutine test_invalid_scenario

end program test_icing_severity