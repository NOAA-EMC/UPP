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
    integer, parameter :: nlvls = 25, ntests = 2
    integer, parameter :: KPV = 5, KTH = 5
    ! 
    integer :: i, j, res
    type(param_t) :: PARAM(1:ntests)
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
    allocate(LVLSXML(1:mxlvl, 1))

    do i = 1, ntests
        allocate(PARAM(i)%level(1:nlvls))
        allocate(PARAM(i)%level2(1:nlvls))
        allocate(PARAM(i)%scale_fact_fixed_sfc1(1:nlvls))
        PARAM(i)%level = 0.0
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

    EXP_LEVEL2(:, IFLD) = 0.0
    EXP_LVLS(:, IFLD) = 1
    EXP_IREC(IFLD) = nlvls

    res = 0
    call SET_LVLSXML(PARAM(IFLD), IFLD, IREC(IFLD), KPV, PV, KTH, TH)

    if (IREC(IFLD) .ne. EXP_IREC(IFLD)) then
        print *, 'Test Case 1 Failed: IREC = ', IREC(IFLD), ' Expected: ', EXP_IREC(IFLD)
        res = 1
    end if

    do i = 1, nlvls
        if (LVLS(i, IFLD) .ne. EXP_LVLS(i, IFLD)) then
            print *, 'Test Case 1 Failed: LVLS(', i, ') = ', LVLS(i, IFLD), &
                     ' Expected: ', EXP_LVLS(i, IFLD)
            res = 1
        end if
        if (LVLSXML(i, IFLD) .ne. EXP_LVLSXML(i, IFLD)) then
            print *, 'Test Case 1 Failed: LVLSXML(', i, ') = ', LVLSXML(i, IFLD), &
                     ' Expected: ', EXP_LVLSXML(i, IFLD)
            res = 1
        end if
        if (abs(PARAM(IFLD)%level(i) - EXP_LEVEL(i, IFLD)) > tol) then
            print *, 'Test Case 1 Failed: PARAM%level(', i, ') = ', PARAM(IFLD)%level(i), &
                     ' Expected: ', EXP_LEVEL(i, IFLD)
            res = 1
        end if
        if (abs(PARAM(IFLD)%level2(i) - EXP_LEVEL2(i, IFLD)) > tol) then
            print *, 'Test Case 1 Failed: PARAM%level2(', i, ') = ', PARAM(IFLD)%level2(i), &
                     ' Expected: ', EXP_LEVEL2(i, IFLD)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10

    ! Test Case 2:
    ! Fixed surface 1 type: isobaric_sfc
    ! Short name contains "ON_ICAO_STD_SFC"
    IFLD = 2
    PARAM(IFLD)%fixed_sfc1_type = 'isobaric_sfc'
    PARAM(IFLD)%shortname = 'ON_ICAO_STD_SFC'

    IREC(IFLD) = 0
    EXP_LEVEL(:, IFLD) = 0.0
    EXP_LEVEL2(:, IFLD) = 0.0
    EXP_LVLS(:, IFLD) = 1
    EXP_IREC(IFLD) = nlvls

    do i = 1, nlvls
        EXP_LVLSXML(i, IFLD) = i
    end do

    res = 0
    call SET_LVLSXML(PARAM(IFLD), IFLD, IREC(IFLD), KPV, PV, KTH, TH)

    if (IREC(IFLD) .ne. EXP_IREC(IFLD)) then
        print *, 'Test Case 2 Failed: IREC = ', IREC(IFLD), ' Expected: ', EXP_IREC(IFLD)
        res = 1
    end if

    do i = 1, nlvls
        if (LVLS(i, IFLD) .ne. EXP_LVLS(i, IFLD)) then
            print *, 'Test Case 2 Failed: LVLS(', i, ') = ', LVLS(i, IFLD), &
                     ' Expected: ', EXP_LVLS(i, IFLD)
            res = 1
        end if
        if (LVLSXML(i, IFLD) .ne. EXP_LVLSXML(i, IFLD)) then
            print *, 'Test Case 2 Failed: LVLSXML(', i, ') = ', LVLSXML(i, IFLD), &
                     ' Expected: ', EXP_LVLSXML(i, IFLD)
            res = 1
        end if
        if (abs(PARAM(IFLD)%level(i) - EXP_LEVEL(i, IFLD)) > tol) then
            print *, 'Test Case 2 Failed: PARAM%level(', i, ') = ', PARAM(IFLD)%level(i), &
                     ' Expected: ', EXP_LEVEL(i, IFLD)
            res = 1
        end if
        if (abs(PARAM(IFLD)%level2(i) - EXP_LEVEL2(i, IFLD)) > tol) then
            print *, 'Test Case 2 Failed: PARAM%level2(', i, ') = ', PARAM(IFLD)%level2(i), &
                     ' Expected: ', EXP_LEVEL2(i, IFLD)
            res = 1
        end if
    end do
    if (res .ne. 0) stop 20

    print *, 'SUCCESS!'
end program test_set_lvlsxml