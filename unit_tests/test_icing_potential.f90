! This is a test program for UPP.
!
! This program tests routines from the IcingPotential module.
!
! Alyson Stahl, 9/2026
program test_icing_potential
    use IcingPotential, only : icing_pot
    use CloudLayers, only : clouds_t
    implicit none

    real, parameter :: tol = 1.0e-6
    integer :: res

    res = 0

    ! Test Case 1: Representative multilayer cloud column spanning cold, mixed-phase, 
    ! and warm-edge icing conditions.
    call test_mixed_phase_cloud_column(res)
    if (res .ne. 0) stop 10

    print *, "Success!"

contains
    subroutine test_mixed_phase_cloud_column(res)
        integer, intent(inout) :: res
        integer, parameter :: nz = 9
        integer :: i
        real :: hgt(nz), rh(nz), t(nz), liqCond(nz), vv(nz)
        type(clouds_t) :: clouds
        !
        real :: ice_pot(nz), exp_ice_pot(nz)

        exp_ice_pot = 0.0
        exp_ice_pot(1) = 0.2564970255
        exp_ice_pot(2) = 0.0
        exp_ice_pot(3) = 0.0
        exp_ice_pot(4) = 0.0
        exp_ice_pot(5) = 0.6958929896
        exp_ice_pot(6) = 1.0
        exp_ice_pot(7) = 0.0573721007
        exp_ice_pot(8) = 0.0
        exp_ice_pot(9) = 0.0


        clouds%nLayers = 3
        clouds%wmnIdx = -1
        clouds%avv = 0.0
        clouds%topIdx = 0
        clouds%baseIdx = 0
        clouds%ctt = 0.0

        clouds%topIdx(1) = 1
        clouds%baseIdx(1) = 2
        clouds%ctt(1) = 220.0

        clouds%topIdx(2) = 3
        clouds%baseIdx(2) = 4
        clouds%ctt(2) = 250.0

        clouds%topIdx(3) = 5
        clouds%baseIdx(3) = 9
        clouds%ctt(3) = 265.0

        hgt = (/ 9000.0, 8000.0, 7000.0, 6000.0, 5000.0, 4000.0, 3000.0, 2000.0, 1000.0 /)
        rh = (/ 100.0, 45.0, 80.0, 96.0, 90.0, 97.0, 60.0, 40.0, 100.0 /)
        t = (/ 250.0, 274.0, 245.0, 274.0, 260.0, 265.0, 270.0, 274.0, 240.0 /)
        liqCond = (/ 0.30, 0.20, 0.20, 0.30, 0.10, 0.30, 0.0005, 0.25, 0.30 /)
        vv = (/ -0.50, 0.20, -0.10, 0.10, -0.40, -0.60, -0.10, 0.20, 0.00 /)

        call icing_pot(hgt, rh, t, liqCond, vv, nz, clouds, ice_pot)

        do i = 1, nz
            if (abs(ice_pot(i) - exp_ice_pot(i)) > tol) then
                print *, "Expected ice_pot(", i, "): ", exp_ice_pot(i), & 
                        " but got: ", ice_pot(i)
                res = 1
            end if
        end do
    end subroutine test_mixed_phase_cloud_column

end program test_icing_potential