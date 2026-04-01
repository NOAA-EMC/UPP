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
    type(param_t), pointer :: PARAM(:)
    integer :: IFLD, IREC(1:ntests), EXP_IREC(1:ntests)
    real :: PV(1:KPV), TH(1:KTH)
    integer :: EXP_LVLS(1:nlvls, 1:ntests), EXP_LVLSXML(1:nlvls, 1:ntests)
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

    ! Initialize input variables with default values. These will be modified as needed for each test case.
    ! Note that some of these variables may not be used in most of the test cases.
    IREC = 0
    PV = 0.0
    TH = 0.0

    isf_surface_physics = 1 

    lsm = nlvls
    nsoil = nlvls
    ifi_nflight = nlvls
    
    allocate(ifi_flight_levels(1:ifi_nflight))
    allocate(SLDPTH(1:nsoil))
    allocate(SLLEVEL(1:nsoil))
    allocate(LVLSXML(1:mxlvl, ntests))

    nullify(PARAM)
    allocate(PARAM(1:ntests))
    do i = 1, ntests
        allocate(PARAM(i)%level(nlvls))
        allocate(PARAM(i)%level2(nlvls))
        allocate(PARAM(i)%scale_fact_fixed_sfc1(nlvls))

        PARAM(i)%level  = 0.0
        PARAM(i)%level2 = 0.0
        PARAM(i)%scale_fact_fixed_sfc1 = 0.0
    end do

    spl = 0.0 
    ifi_flight_levels = 0.0
    SLDPTH = 0.0 
    SLLEVEL = 0.0 

    ! Initialize output arrays
    LVLS = 0
    LVLSXML = 0
    EXP_LVLS = 0
    EXP_LVLSXML = 0
    EXP_LEVEL = 0.0
    EXP_LEVEL2 = 0.0
    EXP_IREC = 0

    ! Test Case 1:
    ! Fixed surface 1 type: isobaric_sfc
    ! Short name does not contain "ON_ICAO_STD_SFC"
    IFLD = 1
    PARAM(IFLD)%fixed_sfc1_type = 'isobaric_sfc'

    ! Specified pressure levels
    do i = 1, nlvls
        spl(i) = 100000.0 - real(i - 1) * 4000.0
        if (i == 1) then
            PARAM(IFLD)%level(nlvls) = spl(1)
            EXP_LEVEL(nlvls, IFLD) = spl(1)
            EXP_LVLSXML(i, IFLD) = nlvls
        else
            PARAM(IFLD)%level(i-1) = spl(i)
            EXP_LEVEL(i-1, IFLD) = spl(i)
            EXP_LVLSXML(i, IFLD) = i - 1
        end if
    end do

    EXP_LVLS(:, IFLD) = 1
    EXP_IREC(IFLD) = nlvls

    call SET_LVLSXML(PARAM(IFLD), IFLD, IREC(IFLD), KPV, PV, KTH, TH)

    ! Test Case 2:
    ! Fixed surface 1 type: isobaric_sfc
    ! Short name contains "ON_ICAO_STD_SFC"
    IFLD = 2
    PARAM(IFLD)%fixed_sfc1_type = 'isobaric_sfc'
    PARAM(IFLD)%shortname = 'ON_ICAO_STD_SFC'

    do i = 1, nlvls
        EXP_LVLS(i, IFLD) = 1
        EXP_LVLSXML(i, IFLD) = i
    end do

    EXP_IREC(IFLD) = nlvls

    call SET_LVLSXML(PARAM(IFLD), IFLD, IREC(IFLD), KPV, PV, KTH, TH)

    ! Test Case 3:
    ! Fixed surface 1 type: hybrid_lvl
    IFLD = 3
    PARAM(IFLD)%fixed_sfc1_type = 'hybrid_lvl'

    EXP_IREC(IFLD) = nlvls

    do i = 1, nlvls
        PARAM(IFLD)%level(i) = real(i) * 2.0
        EXP_LEVEL(i, IFLD) = real(i) * 2.0

        if (mod(i, 2) == 0) then
            EXP_LVLS(i, IFLD) = 1
            EXP_LVLSXML(i, IFLD) = i / 2
        end if
    end do

    call SET_LVLSXML(PARAM(IFLD), IFLD, IREC(IFLD), KPV, PV, KTH, TH)

    ! Test Case 4:
    ! Fixed surface 1 type: depth_bel_land_sfc
    ! Fixed surface 2 type: depth_bel_land_sfc
    ! ISF_SURFACE_PHYSICS != 3
    IFLD = 4
    PARAM(IFLD)%fixed_sfc1_type = 'depth_bel_land_sfc'
    PARAM(IFLD)%fixed_sfc2_type = 'depth_bel_land_sfc'

    PARAM(IFLD)%level2(1) = 30.0
    PARAM(IFLD)%level2(2) = 60.0
    PARAM(IFLD)%level2(3) = 10.0

    SLDPTH = 0.0
    SLDPTH(1) = 0.10
    SLDPTH(2) = 0.20
    SLDPTH(3) = 0.30

    EXP_IREC(IFLD) = 3

    EXP_LEVEL2(1, IFLD) = 30.0
    EXP_LEVEL2(2, IFLD) = 60.0
    EXP_LEVEL2(3, IFLD) = 10.0

    EXP_LVLS(1, IFLD)   = 1
    EXP_LVLSXML(1, IFLD) = 3
    EXP_LVLS(2, IFLD)   = 1
    EXP_LVLSXML(2, IFLD) = 1
    EXP_LVLS(3, IFLD)   = 1
    EXP_LVLSXML(3, IFLD) = 2

    call SET_LVLSXML(PARAM(IFLD), IFLD, IREC(IFLD), KPV, PV, KTH, TH)

    ! Test Case 5:
    ! Fixed surface 1 type: depth_bel_land_sfc
    ! Fixed surface 2 type: depth_bel_land_sfc
    ! ISF_SURFACE_PHYSICS == 3
    IFLD = 5
    PARAM(IFLD)%fixed_sfc1_type = 'depth_bel_land_sfc'
    PARAM(IFLD)%fixed_sfc2_type = 'depth_bel_land_sfc'
    isf_surface_physics = 3

    PARAM(IFLD)%level(1) = 100.0
    PARAM(IFLD)%level(2) = 150.0
    PARAM(IFLD)%level(3) = 50.0

    SLLEVEL = 0.0
    SLLEVEL(1) = 0.5
    SLLEVEL(2) = 1.0
    SLLEVEL(3) = 1.5
    
    EXP_IREC(IFLD) = 25

    EXP_LEVEL(1, IFLD) = 100.0
    EXP_LEVEL(2, IFLD) = 150.0
    EXP_LEVEL(3, IFLD) = 50.0

    EXP_LVLS(1, IFLD)   = 1
    EXP_LVLSXML(1, IFLD) = 3
    EXP_LVLS(2, IFLD)   = 1
    EXP_LVLSXML(2, IFLD) = 1
    EXP_LVLS(3, IFLD)   = 1
    EXP_LVLSXML(3, IFLD) = 2
    EXP_LVLS(4, IFLD)   = 1
    EXP_LVLSXML(4, IFLD) = nlvls

    call SET_LVLSXML(PARAM(IFLD), IFLD, IREC(IFLD), KPV, PV, KTH, TH)

    ! Test Case 6:
    ! Fixed surface 1 type: pot_vort_sfc
    IFLD = 6
    PARAM(IFLD)%fixed_sfc1_type = 'pot_vort_sfc'

    do i = 1, nlvls
        PARAM(IFLD)%level(i) = real(i) * 1.0e-6
        EXP_LEVEL(i, IFLD) = real(i) * 1.0e-6
        if (i <= 3) then
            PARAM(IFLD)%scale_fact_fixed_sfc1(i) = 5
        else
            PARAM(IFLD)%scale_fact_fixed_sfc1(i) = 7
        end if
    end do

    PV(1) = PARAM(IFLD)%level(1) * 10.0**(-1.0 * real(PARAM(IFLD)%scale_fact_fixed_sfc1(1) - 6))
    PV(2) = PARAM(IFLD)%level(3) * 10.0**(-1.0 * real(PARAM(IFLD)%scale_fact_fixed_sfc1(3) - 6)) + 1.0e-4
    PV(3) = PARAM(IFLD)%level(4) * 10.0**(-1.0 * real(PARAM(IFLD)%scale_fact_fixed_sfc1(4) - 6))
    PV(4) = 0.0
    PV(5) = PARAM(IFLD)%level(5) * 10.0**(-1.0 * real(PARAM(IFLD)%scale_fact_fixed_sfc1(5) - 6)) + 1.0e-4

    EXP_IREC(IFLD) = 2
    EXP_LVLS(1, IFLD) = 1
    EXP_LVLSXML(1, IFLD) = 1
    EXP_LVLS(3, IFLD) = 1
    EXP_LVLSXML(3, IFLD) = 4

    call SET_LVLSXML(PARAM(IFLD), IFLD, IREC(IFLD), KPV, PV, KTH, TH)

    res = 0
    do j = 1, ntests
        print *, 'Checking Test Case ', j
        if (IREC(j) .ne. EXP_IREC(j)) then
            print *, 'Test Case ', j, ': IREC = ', IREC(j), ' but expected ', EXP_IREC(j)
            res = 1
        end if
        do i = 1, nlvls
            if (LVLS(i, j) .ne. EXP_LVLS(i, j)) then
                print *, 'Test Case ', j, ': LVLS(', i, ') = ', LVLS(i, j), ' but expected ', EXP_LVLS(i, j)
                res = 1
            end if
            if (LVLSXML(i, j) .ne. EXP_LVLSXML(i, j)) then
                print *, 'Test Case ', j, ': LVLSXML(', i, ') = ', LVLSXML(i, j), ' but expected ', EXP_LVLSXML(i, j)
                res = 1
            end if
            if (PARAM(j)%level(i) .ne. EXP_LEVEL(i, j)) then
                print *, 'Test Case ', j, ': PARAM(', j, ')%level(', i, ') = ', PARAM(j)%level(i), &
                         ' but expected ', EXP_LEVEL(i, j)
                res = 1
            end if
            if (PARAM(j)%level2(i) .ne. EXP_LEVEL2(i, j)) then
                print *, 'Test Case ', j, ': PARAM(', j, ')%level2(', i, ') = ', PARAM(j)%level2(i), &
                         ' but expected ', EXP_LEVEL2(i, j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10
    print *, 'SUCCESS!'

contains


end program test_set_lvlsxml

