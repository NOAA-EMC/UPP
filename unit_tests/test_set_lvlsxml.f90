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

    IFLD = 1
    call SET_LVLSXML(PARAM(IFLD), IFLD, IREC(IFLD), KPV, PV, KTH, TH)


    res = 0
    do i = 1, 1
        print * , 'Test Case ', i, ': IREC = ', IREC(i)

        do j = 1, nlvls
            if (abs(LVLS(j, i) - EXP_LVLS(j, i)) > tol) then
                print *, 'LVLS Test failed at (', j, ',', i, '): ', &
                         'Expected ', EXP_LVLS(j, i), &
                         ' but got ', LVLS(j, i)
                res = 1
            end if
            if (abs(LVLSXML(j, i) - EXP_LVLSXML(j, i)) > tol) then
                print *, 'LVLSXML Test failed at (', j, ',', i, '): ', &
                         'Expected ', EXP_LVLSXML(j, i), &
                         ' but got ', LVLSXML(j, i)
                res = 1
            end if
            if (abs(PARAM(i)%level(j) - EXP_LEVEL(j, i)) > tol) then
                print *, 'PARAM%LEVEL Test failed at (', j, ',', i, '): ', &
                         'Expected ', EXP_LEVEL(j, i), &
                         ' but got ', PARAM(i)%level(j)
                res = 1
            end if
            if (abs(PARAM(i)%level2(j) - EXP_LEVEL2(j, i)) > tol) then
                print *, 'PARAM%LEVEL2 Test failed at (', j, ',', i, '): ', &
                         'Expected ', EXP_LEVEL2(j, i), &
                         ' but got ', PARAM(i)%level2(j)
                res = 1
            end if
        end do
    end do
    print *, 'SUCCESS!'
end program test_set_lvlsxml