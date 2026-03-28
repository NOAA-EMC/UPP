! This is a test program for UPP.
!
! This program tests the SET_LVLSXML() subroutine.
!
! Alyson Stahl, 3/2026
program test_set_lvlsxml
    use xml_perl_data, only: param_t
    use ctlblk_mod, only: lsm, spl, nsoil, isf_surface_physics, nfd, htfd, &
                        petabnd, nbnd, ifi_nflight, ifi_flight_levels, komax
    use soil, only: SLDPTH, SLLEVEL
    use rqstfld_mod, only : mxlvl, LVLS, LVLSXML
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: nlvls = 25, ntests = 24
    integer, parameter :: KPV = 5, KTH = 5
    ! 
    integer :: i, j, res
    integer :: IFLD
    type(param_t) :: PARAM(1:ntests)
    integer :: IREC(1:ntests), EXP_IREC(1:ntests)
    real :: PV(1:KPV), TH(1:KTH)
    real :: EXP_LVLS(1:nlvls, 1:ntests), EXP_LVLSXML(1:nlvls, 1:ntests)
    ! In some test cases, param%level and/or param%level2 will be updated by SET_LVLSXML().
    ! We want to verify that the updated values are correct and that no unexpected 
    ! updates are occurring.
    real :: EXP_LEVEL(1:nlvls, 1:ntests), EXP_LEVEL2(1:nlvls, 1:ntests)

    interface
        subroutine SET_LVLSXML(PARAM, IFLD, IREC, KPV, PV, KTH, TH)
            use xml_perl_data, only: param_t
            type(param_t), intent(inout)    :: PARAM
            integer, intent(inout)          :: IREC
            integer, intent(in)             :: IFLD, KPV, KTH
            real, intent(in)                :: PV(1:KPV), TH(1:KTH)
        end subroutine SET_LVLSXML
    end interface

    ! Initialize input variables
    ! Note that some of these variables may not be used most of the test cases.
    IREC = 0

    ! Land surface model
    ! Used in Test Cases 1 and 2
    lsm = nlvls

    ! Number of model soil levels
    ! Used in Test Cases 4 and 5
    nsoil = nlvls

    ! Surface physics scheme option in model run
    ! Used in Test Cases 4 and 5
    isf_surface_physics = 1

    ! Number of flight levels
    ! Used in Test Cases 8 and 9
    ifi_nflight = nlvls

    allocate(ifi_flight_levels(1:ifi_nflight))
    allocate(SLDPTH(1:nsoil))
    allocate(SLLEVEL(1:nsoil))
    allocate(LVLSXML(1:mxlvl, 1:ntests))
    do i = 1, ntests
        if (i < ntests-1) then
            allocate(PARAM(i)%level(1))
            allocate(PARAM(i)%level2(1))
        else
            allocate(PARAM(i)%level(1:nlvls))
            allocate(PARAM(i)%level2(1:nlvls))
        end if
        PARAM(i)%level = 0.0
        PARAM(i)%level2 = 0.0
    end do

    ! Specified pressure levels
    ! Used in Test Case 1
    do i = 1, nlvls
        spl(i) = 100000.0 - real(i - 1) * 4000.0
    end do

    ! Thickness of each soil layer
    ! Used in Test Case 4
    SLDPTH = 0.0
    SLDPTH(1) = 0.10
    SLDPTH(2) = 0.20
    SLDPTH(3) = 0.30

    ! Soil level
    ! Used in Test Case 5
    SLLEVEL = 0.0
    SLLEVEL(1) = 0.5
    SLLEVEL(2) = 1.0
    SLLEVEL(3) = 1.5
    
    ! Potential vorticity levels
    ! Used in Test Case 6
    PV = 0.0

    ! Isentropic levels
    ! Used in Test Case 7
    TH = 0.0
    TH(1) = 290.0
    TH(2) = 300.0
    TH(3) = 310.0
    TH(4) = 0.0
    TH(5) = 320.0

    ! Flight levels in feet (usually provided by libIFI)
    ! Used in Test Case 8
    do i = 1, ifi_nflight
        ifi_flight_levels(i) = 10000.0 + 1000.0 * real(i)
    end do

    ! Initialize output arrays
    LVLS = 0.0
    LVLSXML = 0.0
    EXP_LVLS = 0.0
    EXP_LVLSXML = 0.0
    EXP_LEVEL = 0.0
    EXP_LEVEL2 = 0.0
    EXP_IREC = 0

    ! Test Case 1:
    ! Fixed surface 1 type: isobaric_sfc
    ! Short name does not contain "ON_ICAO_STD_SFC"
    PARAM(1)%fixed_sfc1_type = 'isobaric_sfc'

    do i = 1, nlvls
        EXP_LVLS(i, 1) = 1
        if (i == 1) then
            PARAM(1)%level(nlvls) = spl(1)
            EXP_LVLSXML(i, 1) = nlvls
        else
            PARAM(1)%level(i-1) = spl(i)
            EXP_LVLSXML(i, 1) = i - 1
        end if
        EXP_LEVEL(i, 1) = PARAM(1)%level(i)
    end do

    ! Test Case 2:
    ! Fixed surface 1 type: isobaric_sfc
    ! Short name contains "ON_ICAO_STD_SFC"
    PARAM(2)%fixed_sfc1_type = 'isobaric_sfc'
    PARAM(2)%shortname = 'ON_ICAO_STD_SFC'

    do i = 1, nlvls
        EXP_LVLS(i, 2) = 1
        EXP_LVLSXML(i, 2) = i
    end do

    ! Test Case 3:
    ! Fixed surface 1 type: hybrid_lvl
    PARAM(3)%fixed_sfc1_type = 'hybrid_lvl'
    do i = 1, nlvls
        PARAM(3)%level(i) = real(i) * 2.0
        EXP_LEVEL(i, 3) = real(i) * 2.0

        if (mod(i, 2) == 0) then
            EXP_LVLS(i, 3) = 1
            EXP_LVLSXML(i, 3) = i / 2
        end if
    end do

    ! Test Case 4:
    ! Fixed surface 1 type: depth_bel_land_sfc
    ! Fixed surface 2 type: depth_bel_land_sfc
    ! ISF_SURFACE_PHYSICS != 3
    PARAM(4)%fixed_sfc1_type = 'depth_bel_land_sfc'
    PARAM(4)%fixed_sfc2_type = 'depth_bel_land_sfc'

    PARAM(4)%level2(1) = 30.0
    PARAM(4)%level2(2) = 60.0
    PARAM(4)%level2(3) = 10.0

    EXP_LEVEL2(1, 4) = 30.0
    EXP_LEVEL2(2, 4) = 60.0
    EXP_LEVEL2(3, 4) = 10.0

    EXP_LVLS(1, 4)   = 1
    EXP_LVLSXML(1, 4) = 3
    EXP_LVLS(2, 4)   = 1
    EXP_LVLSXML(2, 4) = 1
    EXP_LVLS(3, 4)   = 1
    EXP_LVLSXML(3, 4) = 2

    ! Test Case 5:
    ! Fixed surface 1 type: depth_bel_land_sfc
    ! Fixed surface 2 type: depth_bel_land_sfc
    ! ISF_SURFACE_PHYSICS == 3
    PARAM(5)%fixed_sfc1_type = 'depth_bel_land_sfc'
    PARAM(5)%fixed_sfc2_type = 'depth_bel_land_sfc'
    isf_surface_physics = 3

    PARAM(5)%level(1) = 100.0
    PARAM(5)%level(2) = 150.0
    PARAM(5)%level(3) = 50.0

    EXP_LEVEL(1, 5) = 100.0
    EXP_LEVEL(2, 5) = 150.0
    EXP_LEVEL(3, 5) = 50.0

    EXP_LVLS(1, 5)   = 1
    EXP_LVLSXML(1, 5) = 3
    EXP_LVLS(2, 5)   = 1
    EXP_LVLSXML(2, 5) = 1
    EXP_LVLS(3, 5)   = 1
    EXP_LVLSXML(3, 5) = 2
    
    ! Test Case 6:
    ! Fixed surface 1 type: pot_vort_sfc
    PARAM(6)%fixed_sfc1_type = 'pot_vort_sfc'
    allocate(PARAM(6)%scale_fact_fixed_sfc1(nlvls))

    do i = 1, nlvls
        PARAM(6)%level(i) = real(i) * 1.0e-6
        EXP_LEVEL(i, 6) = real(i) * 1.0e-6
        if (i <= 3) then
            PARAM(6)%scale_fact_fixed_sfc1(i) = 5
        else
            PARAM(6)%scale_fact_fixed_sfc1(i) = 7
        end if
    end do

    PV(1) = PARAM(6)%level(1) * 10.0**(-1.0 * real(PARAM(6)%scale_fact_fixed_sfc1(1) - 6))
    PV(2) = PARAM(6)%level(3) * 10.0**(-1.0 * real(PARAM(6)%scale_fact_fixed_sfc1(3) - 6)) + 1.0e-4
    PV(3) = PARAM(6)%level(4) * 10.0**(-1.0 * real(PARAM(6)%scale_fact_fixed_sfc1(4) - 6))
    PV(4) = 0.0
    PV(5) = PARAM(6)%level(5) * 10.0**(-1.0 * real(PARAM(6)%scale_fact_fixed_sfc1(5) - 6)) + 1.0e-4

    EXP_LVLS(1, 6) = 1
    EXP_LVLSXML(1, 6) = 1
    EXP_LVLS(3, 6) = 1
    EXP_LVLSXML(3, 6) = 4

    ! Test Case 7:
    ! Fixed surface 1 type: isentropic_lvl
    PARAM(7)%fixed_sfc1_type = 'isentropic_lvl'

    PARAM(7)%level(1) = 290.0
    PARAM(7)%level(2) = 305.0
    PARAM(7)%level(3) = 310.0
    PARAM(7)%level(4) = 295.0
    PARAM(7)%level(5) = 325.0

    do i = 1, 5
        EXP_LEVEL(i, 7) = PARAM(7)%level(i)
    end do

    EXP_LVLS(1, 7)   = 1
    EXP_LVLSXML(1, 7) = 1
    EXP_LVLS(3, 7)   = 1
    EXP_LVLSXML(3, 7) = 3

    ! Test Case 8:
    ! Fixed surface 1 type: spec_alt_above_mean_sea_lvl
    ! Short name contains "SPECIFIC_IFI_FLIGHT_LEVEL"
    PARAM(8)%fixed_sfc1_type = 'spec_alt_above_mean_sea_lvl'
    PARAM(8)%shortname = 'SPECIFIC_IFI_FLIGHT_LEVEL'

    do j = 1, nlvls
        if (j == nlvls) then
            PARAM(8)%level(j) = ifi_flight_levels(j) + 50.0
        else
            i = mod(j, 24) + 1
            PARAM(8)%level(j) = ifi_flight_levels(i)
            EXP_LVLS(i, 8) = 1
            EXP_LVLSXML(i, 8) = j
        end if
        EXP_LEVEL(j, 8) = PARAM(8)%level(j)
    end do

    ! Test Case 9:
    ! Fixed surface 1 type: spec_alt_above_mean_sea_lvl
    ! Short name contains "IFI_FLIGHT_LEVEL"
    PARAM(9)%fixed_sfc1_type = 'spec_alt_above_mean_sea_lvl'
    PARAM(9)%shortname = 'IFI_FLIGHT_LEVEL'

    do i = 1, ifi_nflight
        EXP_LVLS(i, 9)   = 1
        EXP_LVLSXML(i, 9) = i
    end do

    ! Test Case 10:
    ! Fixed surface 1 type: spec_alt_above_mean_sea_lvl
    ! Short name contains "GTG_ON_SPEC_ALT_ABOVE_MEAN_SEA_LVL"
    PARAM(10)%fixed_sfc1_type = 'spec_alt_above_mean_sea_lvl'
    PARAM(10)%shortname = 'GTG_ON_SPEC_ALT_ABOVE_MEAN_SEA_LVL'

    do i = 1, nlvls
        EXP_LVLS(i, 10)   = 1
        EXP_LVLSXML(i, 10) = i
    end do

    ! Test Case 11:
    ! Fixed surface 1 type: spec_alt_above_mean_sea_lvl
    ! Short name is not set
    PARAM(11)%fixed_sfc1_type = 'spec_alt_above_mean_sea_lvl'

    PARAM(11)%level(1) = HTFD(1)
    PARAM(11)%level(2) = HTFD(8)

    EXP_LEVEL(1, 11) = HTFD(1)
    EXP_LEVEL(2, 11) = HTFD(8)

    EXP_LVLS(1, 11)   = 2
    EXP_LVLSXML(1, 11) = 1
    EXP_LVLS(8, 11)   = 1
    EXP_LVLSXML(8, 11) = 2

    ! Test Case 12:
    ! Fixed surface 1 type: spec_prec_above_grnd
    ! Short name == "MIXED_LAYER_CAPE_ON_SPEC_PRES_ABOVE_GRND"
    PARAM(12)%fixed_sfc1_type = 'spec_prec_above_grnd'
    PARAM(12)%shortname = 'MIXED_LAYER_CAPE_ON_SPEC_PRES_ABOVE_GRND'

    EXP_LVLSXML(1, 12) = 1
    EXP_LEVEL(1, 12) = nint(PETABND(3)+15.)*100
    EXP_LEVEL2(1, 12) = nint(PETABND(1)-15.)*100
    EXP_IREC(12) = 1

    ! Test Case 13:
    ! Fixed surface 1 type: spec_prec_above_grnd
    ! Short name == "MIXED_LAYER_CIN_ON_SPEC_PRES_ABOVE_GRND"
    PARAM(13)%fixed_sfc1_type = 'spec_prec_above_grnd'
    PARAM(13)%shortname = "MIXED_LAYER_CIN_ON_SPEC_PRES_ABOVE_GRND"

    EXP_LVLSXML(1, 13) = 1
    EXP_LEVEL(1, 13) = nint(PETABND(3)+15.)*100
    EXP_LEVEL2(1, 13) = nint(PETABND(1)-15.)*100
    EXP_IREC(13) = 1

    ! Test Case 14:
    ! Fixed surface 1 type: spec_prec_above_grnd
    ! Short name == "UNSTABLE_CAPE_ON_SPEC_PRES_ABOVE_GRND"
    PARAM(14)%fixed_sfc1_type = 'spec_prec_above_grnd'
    PARAM(14)%shortname = "UNSTABLE_CAPE_ON_SPEC_PRES_ABOVE_GRND"

    EXP_LVLSXML(1, 14) = 1
    EXP_LEVEL(1, 14) = 25500
    EXP_LEVEL2(1, 14) = 0
    EXP_IREC(14) = 1

    ! Test Case 15:
    ! Fixed surface 1 type: spec_prec_above_grnd
    ! Short name == "UNSTABLE_CIN_ON_SPEC_PRES_ABOVE_GRND"
    PARAM(15)%fixed_sfc1_type = 'spec_prec_above_grnd'
    PARAM(15)%shortname = "UNSTABLE_CIN_ON_SPEC_PRES_ABOVE_GRND"

    EXP_LVLSXML(1, 15) = 1
    EXP_LEVEL(1, 15) = 25500
    EXP_LEVEL2(1, 15) = 0
    EXP_IREC(15) = 1
    
    ! Test Case 16:
    ! Fixed surface 1 type: spec_prec_above_grnd
    ! Short name == "BEST_CAPE_ON_SPEC_PRES_ABOVE_GRND"
    PARAM(16)%fixed_sfc1_type = 'spec_prec_above_grnd'
    PARAM(16)%shortname = "BEST_CAPE_ON_SPEC_PRES_ABOVE_GRND"

    EXP_LVLSXML(1, 16) = 1
    EXP_LEVEL(1, 16) = nint(PETABND(NBND)+15.)*100
    EXP_LEVEL2(1, 16) = nint(PETABND(1)-15.)*100
    EXP_IREC(16) = 1

    ! Test Case 17:
    ! Fixed surface 1 type: spec_prec_above_grnd
    ! Short name == "BEST_CIN_ON_SPEC_PRES_ABOVE_GRND"
    PARAM(17)%fixed_sfc1_type = 'spec_prec_above_grnd'
    PARAM(17)%shortname = "BEST_CIN_ON_SPEC_PRES_ABOVE_GRND"

    EXP_LVLSXML(1, 17) = 1
    EXP_LEVEL(1, 17) = nint(PETABND(NBND)+15.)*100
    EXP_LEVEL2(1, 17) = nint(PETABND(1)-15.)*100
    EXP_IREC(17) = 1

    ! Test Case 18:
    ! Fixed surface 1 type: spec_prec_above_grnd
    ! Short name is not listed above
    PARAM(18)%fixed_sfc1_type = 'spec_prec_above_grnd'

    PARAM(18)%level = 0.0
    PARAM(18)%level(1) = 25500.0
    PARAM(18)%level(2) = 5000.0
    PARAM(18)%level(3) = (PETABND(2) + 15.0) * 100.0
    PARAM(18)%level(4) = (PETABND(4) + 15.0) * 100.0
    PARAM(18)%level(5) = (PETABND(6) + 15.0) * 100.0

    do i = 1, 5
        EXP_LEVEL(i, 18) = PARAM(18)%level(i)
    end do

    EXP_LVLS(1, 18)    = 1
    EXP_LVLSXML(1, 18) = 1
    EXP_LVLS(2, 18)    = 1
    EXP_LVLSXML(2, 18) = 3
    EXP_LVLS(4, 18)    = 1
    EXP_LVLSXML(4, 18) = 4
    EXP_LVLS(6, 18)    = 1
    EXP_LVLSXML(6, 18) = 5

    ! Test Case 19:
    ! Fixed surface 1 type: 'spec_hgt_lvl_above_grnd'
    ! Short name contains "SPEC_HGT_LVL_ABOVE_GRND_FDHGT"
    PARAM(19)%fixed_sfc1_type = 'spec_hgt_lvl_above_grnd'
    PARAM(19)%shortname = 'SPEC_HGT_LVL_ABOVE_GRND_FDHGT'

    PARAM(19)%level = 0.0
    PARAM(19)%level(1) = HTFD(2)
    PARAM(19)%level(2) = HTFD(5)
    PARAM(19)%level(3) = HTFD(10)
    PARAM(19)%level(4) = HTFD(1)
    PARAM(19)%level(5) = HTFD(8)

    do i = 1, 5
        EXP_LEVEL(i, 19) = PARAM(19)%level(i)
    end do

    EXP_LVLS(1, 19)    = 1
    EXP_LVLSXML(1, 19) = 4
    EXP_LVLS(2, 19)    = 1
    EXP_LVLSXML(2, 19) = 1
    EXP_LVLS(5, 19)    = 1
    EXP_LVLSXML(5, 19) = 2
    EXP_LVLS(8, 19)    = 1
    EXP_LVLSXML(8, 19) = 5
    EXP_LVLS(10, 19)   = 1
    EXP_LVLSXML(10, 19)= 3

    ! Test Case 20:
    ! Fixed surface 1 type: 'spec_hgt_lvl_above_grnd'
    ! Short name does not contain "SPEC_HGT_LVL_ABOVE_GRND_FDHGT"
    PARAM(20)%fixed_sfc1_type = 'spec_hgt_lvl_above_grnd'

    do i = 1, nlvls
        EXP_LVLS(i, 20)    = 1
        EXP_LVLSXML(i, 20) = i
    end do

    ! Test Case 21:
    ! Short name  == 'TMP_ON_SIGMA_LVL_HPC'
    PARAM(21)%shortname = 'TMP_ON_SIGMA_LVL_HPC'

    PARAM(21)%level = 0.0
    PARAM(21)%level(1) = 8000.0
    PARAM(21)%level(2) = 7000.0
    PARAM(21)%level(3) = 9000.0
    PARAM(21)%level(4) = 7500.0
    PARAM(21)%level(5) = 8500.0

    do i = 1, 5
        EXP_LEVEL(i, 21) = PARAM(21)%level(i)
    end do

    EXP_LVLS(1, 21)    = 1
    EXP_LVLSXML(1, 21) = 2
    EXP_LVLS(2, 21)    = 1
    EXP_LVLSXML(2, 21) = 4
    EXP_LVLS(3, 21)    = 1
    EXP_LVLSXML(3, 21) = 1
    EXP_LVLS(4, 21)    = 1
    EXP_LVLSXML(4, 21) = 5
    EXP_LVLS(5, 21)    = 1
    EXP_LVLSXML(5, 21) = 3

    ! Test Case 22:
    ! Short name containing 'SIGMA_LVLS'
    PARAM(22)%shortname = 'SIGMA_LVLS'

    PARAM(22)%level = 0.0
    PARAM(22)%level(1) = 4550.0
    PARAM(22)%level(2) = 530.0
    PARAM(22)%level(3) = 7585.0
    PARAM(22)%level(4) = 9835.0
    PARAM(22)%level(5) = 2605.0

    do i = 1, 5
        EXP_LEVEL(i, 22) = PARAM(22)%level(i)
    end do

    EXP_LVLS(1, 22)    = 1
    EXP_LVLSXML(1, 22) = 2
    EXP_LVLS(3, 22)    = 1
    EXP_LVLSXML(3, 22) = 5
    EXP_LVLS(5, 22)    = 1
    EXP_LVLSXML(5, 22) = 1
    EXP_LVLS(10, 22)   = 1
    EXP_LVLSXML(10, 22)= 3
    EXP_LVLS(20, 22)   = 1
    EXP_LVLSXML(20, 22)= 4

    ! Test Case 23:
    ! Fixed surface 1 type: spec_prec_above_grnd
    ! nlevels == 1
    PARAM(23)%fixed_sfc1_type = 'spec_prec_above_grnd'

    EXP_LVLS(1, 23)    = 1
    EXP_LVLSXML(1, 23) = 1

    ! Test Case 24:
    ! Unrecognized fixed surface type and short name
    EXP_LVLS(1,24) = 1
    EXP_LVLSXML(1,24) = 1
    EXP_IREC(24) = 1

    res = 0
    do IFLD = 1, ntests
        call SET_LVLSXML(PARAM(IFLD), IFLD, IREC(IFLD), KPV, PV, KTH, TH)

        print * , 'Test Case ', IFLD, ': IREC = ', IREC(IFLD)

        do i = 1, nlvls
            if (abs(LVLS(i, IFLD) - EXP_LVLS(i, IFLD)) > tol) then
                print *, 'LVLS Test failed at (', i, ',', IFLD, '): ', &
                         'Expected ', EXP_LVLS(i, IFLD), &
                         ' but got ', LVLS(i, IFLD)
                res = 1
            end if
            if (abs(LVLSXML(i, IFLD) - EXP_LVLSXML(i, IFLD)) > tol) then
                print *, 'LVLSXML Test failed at (', i, ',', IFLD, '): ', &
                         'Expected ', EXP_LVLSXML(i, IFLD), &
                         ' but got ', LVLSXML(i, IFLD)
                res = 1
            end if
            if (abs(PARAM(IFLD)%level(i) - EXP_LEVEL(i, IFLD)) > tol) then
                print *, 'PARAM%LEVEL Test failed at (', i, ',', IFLD, '): ', &
                         'Expected ', EXP_LEVEL(i, IFLD), &
                         ' but got ', PARAM(IFLD)%level(i)
                res = 1
            end if
            if (abs(PARAM(IFLD)%level2(i) - EXP_LEVEL2(i, IFLD)) > tol) then
                print *, 'PARAM%LEVEL2 Test failed at (', i, ',', IFLD, '): ', &
                         'Expected ', EXP_LEVEL2(i, IFLD), &
                         ' but got ', PARAM(IFLD)%level2(i)
                res = 1
            end if
        end do
    end do
    print *, 'SUCCESS!'
end program test_set_lvlsxml