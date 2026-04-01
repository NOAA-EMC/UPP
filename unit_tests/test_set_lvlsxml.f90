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
    integer, parameter :: nlvls = 25
    integer, parameter :: KPV = 5, KTH = 5, IFLD = 1
    ! 
    integer :: i, j, res
    type(param_t) :: PARAM, PARAM_ONE_LEVEL
    integer :: IREC, EXP_IREC
    real :: PV(1:KPV), TH(1:KTH)
    integer :: EXP_LVLS(1:nlvls), EXP_LVLSXML(1:nlvls)
    ! In some test cases, param%level and/or param%level2 will be updated by SET_LVLSXML().
    ! We want to verify that the updated values are correct and that no unexpected 
    ! updates are occurring.
    real :: EXP_LEVEL(1:nlvls), EXP_LEVEL2(1:nlvls)

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
    allocate(LVLSXML(1:mxlvl, 1))

    allocate(PARAM%level(1:nlvls))
    allocate(PARAM%level2(1:nlvls))
    allocate(PARAM%scale_fact_fixed_sfc1(1:nlvls))

    allocate(PARAM_ONE_LEVEL%level(1))
    allocate(PARAM_ONE_LEVEL%level2(1))

    spl = 0.0 
    ifi_flight_levels = 0.0
    SLDPTH = 0.0 
    SLLEVEL = 0.0 

    PARAM%level = 0.0
    PARAM%level2 = 0.0
    PARAM%scale_fact_fixed_sfc1 = 0.0

    PARAM_ONE_LEVEL%level = 0.0
    PARAM_ONE_LEVEL%level2 = 0.0

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
    PARAM%fixed_sfc1_type = 'isobaric_sfc'
    PARAM%shortname = ''
    PARAM%level = 0.0
    PARAM%level2 = 0.0
    PARAM%scale_fact_fixed_sfc1 = 0.0

    ! Specified pressure levels
    do i = 1, nlvls
        spl(i) = 100000.0 - real(i - 1) * 4000.0
        if (i == 1) then
            PARAM%level(nlvls) = spl(1)
            EXP_LEVEL(nlvls) = spl(1)
            EXP_LVLSXML(i) = nlvls
        else
            PARAM%level(i-1) = spl(i)
            EXP_LEVEL(i-1) = spl(i)
            EXP_LVLSXML(i) = i - 1
        end if
    end do

    EXP_LEVEL2 = 0.0
    EXP_LVLS = 1
    EXP_IREC = nlvls

    res = 0
    call SET_LVLSXML(PARAM, IFLD, IREC, KPV, PV, KTH, TH)

    call check_test_case(1, IREC, EXP_IREC, PARAM, EXP_LVLS, EXP_LVLSXML, &
                    EXP_LEVEL, EXP_LEVEL2, res)

    if (res .ne. 0) stop 10

    ! Test Case 2:
    ! Fixed surface 1 type: isobaric_sfc
    ! Short name contains "ON_ICAO_STD_SFC"
    PARAM%fixed_sfc1_type = 'isobaric_sfc'
    PARAM%shortname = 'ON_ICAO_STD_SFC'
    PARAM%level = 0.0
    PARAM%level2 = 0.0
    PARAM%scale_fact_fixed_sfc1 = 0.0

    IREC = 0
    EXP_LEVEL = 0.0
    EXP_LEVEL2 = 0.0
    EXP_LVLS = 1
    EXP_IREC = nlvls

    do i = 1, nlvls
        EXP_LVLSXML(i) = i
    end do

    res = 0
    call SET_LVLSXML(PARAM, IFLD, IREC, KPV, PV, KTH, TH)

    call check_test_case(2, IREC, EXP_IREC, PARAM, EXP_LVLS, EXP_LVLSXML, &
                    EXP_LEVEL, EXP_LEVEL2, res)

    if (res .ne. 0) stop 20

    print *, 'SUCCESS!'

contains

    subroutine check_test_case(num, irec, exp_irec, param, &
                               exp_lvls, exp_lvlsxml,      &
                               exp_level, exp_level2, res)
        implicit none

        integer,       intent(in)    :: num
        integer,       intent(in)    :: irec, exp_irec
        type(param_t), intent(in)    :: param
        integer,       intent(in)    :: exp_lvls(nlvls),  exp_lvlsxml(nlvls)
        real,          intent(in)    :: exp_level(nlvls), exp_level2(nlvls)
        integer,       intent(inout) :: res

        integer :: i

        if (irec .ne. exp_irec) then
            print *, 'Test Case ', num, ' Failed: IREC = ', irec, &
                     ' Expected: ', exp_irec
            res = 1
        end if

        do i = 1, nlvls
            if (LVLS(i, IFLD) .ne. exp_lvls(i)) then
                print *, 'Test Case ', num, ' Failed: LVLS(', i, ') = ', &
                         LVLS(i, IFLD), ' Expected: ', exp_lvls(i)
                res = 1
            end if
            if (LVLSXML(i, IFLD) .ne. exp_lvlsxml(i)) then
                print *, 'Test Case ', num, ' Failed: LVLSXML(', i, ') = ', &
                         LVLSXML(i, IFLD), ' Expected: ', exp_lvlsxml(i)
                res = 1
            end if
            if (abs(param%level(i) - exp_level(i)) > tol) then
                print *, 'Test Case ', num, ' Failed: PARAM%level(', i, ') = ', &
                         param%level(i), ' Expected: ', exp_level(i)
                res = 1
            end if
            if (abs(param%level2(i) - exp_level2(i)) > tol) then
                print *, 'Test Case ', num, ' Failed: PARAM%level2(', i, ') = ', &
                         param%level2(i), ' Expected: ', exp_level2(i)
                res = 1
            end if
        end do

    end subroutine check_test_case

end program test_set_lvlsxml

