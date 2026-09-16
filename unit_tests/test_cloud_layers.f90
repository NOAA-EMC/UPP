! This is a test program for UPP.
!
! This program tests routines from the CloudLayers module.
!
! Alyson Stahl, 9/2026
program test_cloud_layers
    use CloudLayers, only : calc_CloudLayers, clouds_t
    implicit none

    real, parameter :: tol = 1.0e-6
    integer :: res

    res = 0

    ! Test Case 1: Valid & representative multi-layer cloud profile with warm nose. 
    call test_multilayer_and_warmnose(res)
    if (res .ne. 0) stop 10

    ! Test Case 2: Valid no-cloud profile that exercises the num_lyr <= 0 branch.
    call test_no_cloud_layers(res)
    if (res .ne. 0) stop 20

    ! Test Case 3: Valid midlatitude cloud profile with a top cloud whose base is 
    ! explicitly identified.
    call test_top_cloud_base_detection(res)
    if (res .ne. 0) stop 30

    print *, "Success!"

contains

    subroutine test_multilayer_and_warmnose(res)
        integer, intent(inout) :: res
        integer, parameter :: nz = 12
        integer :: i
        integer :: topoK
        real :: xlat, xlon
        real :: rh(nz), t(nz), pres(nz), ept(nz), vv(nz)
        !
        integer :: region, exp_region
        type(clouds_t) :: clouds, exp_clouds

        exp_region = 2
        
        exp_clouds%nLayers = 2
        exp_clouds%wmnIdx = 8
        exp_clouds%avv = -0.75

        exp_clouds%topIdx = 0
        exp_clouds%topIdx(1) =1
        exp_clouds%topIdx(2) = 8

        exp_clouds%baseIdx = 0
        exp_clouds%baseIdx(1) = 3
        exp_clouds%baseIdx(2) = 12

        exp_clouds%ctt = 0.0
        exp_clouds%ctt(1) = 235.0
        exp_clouds%ctt(2) = 274.0
 
        allocate(exp_clouds%layerQ(nz))
        exp_clouds%layerQ = 0.0
        exp_clouds%layerQ(8) = 1.81183428E-01
        exp_clouds%layerQ(9) = 3.33446525E-02
        exp_clouds%layerQ(10) = 4.4320303E-01

        topoK = nz
        xlat = 70.0
        xlon = -150.0

        rh = (/ 95.0, 92.0, 68.0, 65.0, 60.0, 66.0, & 
                68.0, 70.0, 70.0, 100.0, 102.0, 98.0 /)
        t = (/ 235.0, 240.0, 245.0, 250.0, 254.0, 258.0, &
                262.0, 274.0, 256.0, 259.0, 265.0, 268.0 /)
        pres = (/ 30000.0, 35000.0, 40000.0, 45000.0, 50000.0, &
                    60000.0, 70000.0, 76000.0, 82000.0, 88000.0, 94000.0, 100000.0 /)
        ept = (/ 320.0, 315.0, 310.0, 307.0, 305.0, 303.0, &
                301.0, 310.0, 312.0, 309.0, 307.0, 300.0 /)
        vv = (/ -0.05, -0.08, -0.10, -0.12, -0.15, -0.20, &
                -0.25, -0.30, -0.40, -0.55, -0.10, -0.05 /)

        region = -1
        clouds%nLayers = 0
        clouds%wmnIdx = -1
        clouds%avv = 0.0
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        allocate(clouds%layerQ(nz))
        clouds%layerQ = 0.0

        call calc_CloudLayers(rh, t, pres, ept, vv, nz, topoK, xlat, xlon, region, clouds)

        if (clouds%nLayers .ne. exp_clouds%nLayers) then
            print *, "Expected nLayers: ", exp_clouds%nLayers, &
                     " but got: ", clouds%nLayers
            res = 1
        end if

        if (clouds%wmnIdx .ne. exp_clouds%wmnIdx) then
            print *, "Expected wmnIdx: ", exp_clouds%wmnIdx, &
                     " but got: ", clouds%wmnIdx
            res = 1
        end if

        if (abs(clouds%avv - exp_clouds%avv) > tol) then
            print *, "Expected avv: ", exp_clouds%avv, &
                     " but got: ", clouds%avv
            res = 1
        end if

        do i = 1, nz
            if (clouds%topIdx(i) .ne. exp_clouds%topIdx(i)) then
                print *, "Expected topIdx(", i, "): ", exp_clouds%topIdx(i), &
                         " but got: ", clouds%topIdx(i)
                res = 1
            end if
            if (clouds%baseIdx(i) .ne. exp_clouds%baseIdx(i)) then
                print *, "Expected baseIdx(", i, "): ", exp_clouds%baseIdx(i), &
                         " but got: ", clouds%baseIdx(i)
                res = 1
            end if
            if (abs(clouds%ctt(i) - exp_clouds%ctt(i)) > tol) then
                print *, "Expected ctt(", i, "): ", exp_clouds%ctt(i), &
                         " but got: ", clouds%ctt(i)
                res = 1
            end if
            if (abs(clouds%layerQ(i) - exp_clouds%layerQ(i)) > tol) then
                print *, "Expected layerQ(", i, "): ", exp_clouds%layerQ(i), &
                         " but got: ", clouds%layerQ(i)
                res = 1
            end if
        end do

    end subroutine test_multilayer_and_warmnose

    subroutine test_no_cloud_layers(res)
        integer, intent(inout) :: res
        integer, parameter :: nz = 8
        integer :: i
        integer :: topoK
        real :: xlat, xlon
        real :: rh(nz), t(nz), pres(nz), ept(nz), vv(nz)
        !
        integer :: region, exp_region
        type(clouds_t) :: clouds, exp_clouds

        exp_region = 1

        exp_clouds%nLayers = 0
        exp_clouds%wmnIdx = -1
        exp_clouds%avv = 0.0

        exp_clouds%topIdx = 0
        exp_clouds%baseIdx = 0
        exp_clouds%ctt = 0.0

        allocate(exp_clouds%layerQ(nz))
        exp_clouds%layerQ = 0.0

        topoK = nz
        xlat = 15.0
        xlon = -75.0

        rh = (/ 55.0, 58.0, 60.0, 62.0, 64.0, 66.0, 68.0, 70.0 /)
        t = (/ 238.0, 243.0, 248.0, 253.0, 258.0, 263.0, 268.0, 272.0 /)
        pres = (/ 30000.0, 40000.0, 50000.0, 60000.0, 70000.0, 80000.0, 90000.0, 100000.0 /)
        ept = (/ 300.0, 302.0, 304.0, 306.0, 308.0, 310.0, 312.0, 314.0 /)
        vv = (/ -0.05, -0.06, -0.08, -0.10, -0.12, -0.14, -0.16, -0.18 /)
        region = -1
        clouds%nLayers = 0
        clouds%wmnIdx = -1
        clouds%avv = 0.0
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        allocate(clouds%layerQ(nz))
        clouds%layerQ = 0.0

        call calc_CloudLayers(rh, t, pres, ept, vv, nz, topoK, xlat, xlon, region, clouds)

        if (clouds%nLayers .ne. exp_clouds%nLayers) then
            print *, "Expected nLayers: ", exp_clouds%nLayers, &
                     " but got: ", clouds%nLayers
            res = 1
        end if

        if (clouds%wmnIdx .ne. exp_clouds%wmnIdx) then
            print *, "Expected wmnIdx: ", exp_clouds%wmnIdx, &
                     " but got: ", clouds%wmnIdx
            res = 1
        end if

        if (abs(clouds%avv - exp_clouds%avv) > tol) then
            print *, "Expected avv: ", exp_clouds%avv, &
                     " but got: ", clouds%avv
            res = 1
        end if

        do i = 1, nz
            if (clouds%topIdx(i) .ne. exp_clouds%topIdx(i)) then
                print *, "Expected topIdx(", i, "): ", exp_clouds%topIdx(i), &
                         " but got: ", clouds%topIdx(i)
                res = 1
            end if
            if (clouds%baseIdx(i) .ne. exp_clouds%baseIdx(i)) then
                print *, "Expected baseIdx(", i, "): ", exp_clouds%baseIdx(i), &
                         " but got: ", clouds%baseIdx(i)
                res = 1
            end if
            if (abs(clouds%ctt(i) - exp_clouds%ctt(i)) > tol) then
                print *, "Expected ctt(", i, "): ", exp_clouds%ctt(i), &
                         " but got: ", clouds%ctt(i)
                res = 1
            end if
            if (abs(clouds%layerQ(i) - exp_clouds%layerQ(i)) > tol) then
                print *, "Expected layerQ(", i, "): ", exp_clouds%layerQ(i), &
                         " but got: ", clouds%layerQ(i)
                res = 1
            end if
        end do

    end subroutine test_no_cloud_layers

    subroutine test_top_cloud_base_detection(res)
        integer, intent(inout) :: res
        integer, parameter :: nz = 6
        integer :: i
        integer :: topoK
        real :: xlat, xlon
        real :: rh(nz), t(nz), pres(nz), ept(nz), vv(nz)
        !
        integer :: region, exp_region
        type(clouds_t) :: clouds, exp_clouds

        exp_region = 2

        exp_clouds%nLayers = 1
        exp_clouds%wmnIdx = -1
        exp_clouds%avv = 0.0

        exp_clouds%topIdx = 0
        exp_clouds%topIdx(1) = 3

        exp_clouds%baseIdx = 0
        exp_clouds%baseIdx(1) = 5

        exp_clouds%ctt = 0.0
        exp_clouds%ctt(1) = 262.0

        allocate(exp_clouds%layerQ(nz))
        exp_clouds%layerQ = 0.0

        topoK = nz
        xlat = 40.0
        xlon = -100.0

        rh = (/ 60.0, 70.0, 78.0, 78.0, 70.0, 65.0 /)
        t = (/ 250.0, 256.0, 262.0, 262.0, 262.0, 270.0 /)
        pres = (/ 50000.0, 60000.0, 70000.0, 80000.0, 90000.0, 100000.0 /)
        ept = (/ 300.0, 302.0, 305.0, 305.0, 305.0, 306.0 /)
        vv = (/ -0.05, -0.08, -0.10, -0.12, -0.14, -0.16 /)

        region = -1
        clouds%nLayers = 0
        clouds%wmnIdx = -1
        clouds%avv = 0.0
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0
        allocate(clouds%layerQ(nz))
        clouds%layerQ = 0.0

        call calc_CloudLayers(rh, t, pres, ept, vv, nz, topoK, xlat, xlon, region, clouds)

        if (clouds%nLayers .ne. exp_clouds%nLayers) then
            print *, "Expected nLayers: ", exp_clouds%nLayers, &
                     " but got: ", clouds%nLayers
            res = 1
        end if

        if (clouds%wmnIdx .ne. exp_clouds%wmnIdx) then
            print *, "Expected wmnIdx: ", exp_clouds%wmnIdx, &
                     " but got: ", clouds%wmnIdx
            res = 1
        end if

        if (abs(clouds%avv - exp_clouds%avv) > tol) then
            print *, "Expected avv: ", exp_clouds%avv, &
                     " but got: ", clouds%avv
            res = 1
        end if

        do i = 1, nz
            if (clouds%topIdx(i) .ne. exp_clouds%topIdx(i)) then
                print *, "Expected topIdx(", i, "): ", exp_clouds%topIdx(i), &
                         " but got: ", clouds%topIdx(i)
                res = 1
            end if
            if (clouds%baseIdx(i) .ne. exp_clouds%baseIdx(i)) then
                print *, "Expected baseIdx(", i, "): ", exp_clouds%baseIdx(i), &
                         " but got: ", clouds%baseIdx(i)
                res = 1
            end if
            if (abs(clouds%ctt(i) - exp_clouds%ctt(i)) > tol) then
                print *, "Expected ctt(", i, "): ", exp_clouds%ctt(i), &
                         " but got: ", clouds%ctt(i)
                res = 1
            end if
            if (abs(clouds%layerQ(i) - exp_clouds%layerQ(i)) > tol) then
                print *, "Expected layerQ(", i, "): ", exp_clouds%layerQ(i), &
                         " but got: ", clouds%layerQ(i)
                res = 1
            end if
        end do

    end subroutine test_top_cloud_base_detection

end program test_cloud_layers