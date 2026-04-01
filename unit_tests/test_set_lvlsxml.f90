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
    real :: EXP_LVLS(1:nlvls, 1), EXP_LVLSXML(1:nlvls, 1)
    ! In some test cases, param%level and/or param%level2 will be updated by SET_LVLSXML().
    ! We want to verify that the updated values are correct and that no unexpected 
    ! updates are occurring.
    real :: EXP_LEVEL(1:nlvls, 1), EXP_LEVEL2(1:nlvls, 1)

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
    PARAM%fixed_sfc1_type = 'isobaric_sfc'
    PARAM%level = 0.0
    PARAM%level2 = 0.0
    PARAM%scale_fact_fixed_sfc1 = 0.0

    ! Specified pressure levels
    do i = 1, nlvls
        spl(i) = 100000.0 - real(i - 1) * 4000.0
        if (i == 1) then
            PARAM%level(nlvls) = spl(1)
        else
            PARAM%level(i-1) = spl(i)
        end if
    end do

    call SET_LVLSXML(PARAM, IFLD, IREC, KPV, PV, KTH, TH)

    print *, "IREC = ", IREC
    do i = 1, nlvls
        print *, "LVLS(", i, ") = ", LVLS(i, IFLD)
    end do
    do i = 1, nlvls
        print *, "LVLSXML(", i, ") = ", LVLSXML(i, IFLD)
    end do
    do i = 1, nlvls
        print *, "PARAM%level(", i, ") = ", PARAM%level(i)
    end do
    do i = 1, nlvls
        print *, "PARAM%level2(", i, ") = ", PARAM%level2(i)
    end do
    
    print *, 'SUCCESS!'
end program test_set_lvlsxml