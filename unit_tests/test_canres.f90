! This is a test program for UPP.
!
! This program tests the CANRES() subroutine.
!
! Alyson Stahl, 1/2026
program test_canres
    use ctlblk_mod, only: novegtype, nsoil, ivegsrc
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 4
    integer :: res
    ! Input variables
    integer :: IVEG, ISOIL
    real :: SOLAR, SFCTMP, Q2, SFCPRS
    real, dimension(1:npts) :: SMC, SLDPTH
    ! Output variables
    integer :: NROOTS, EXP_NROOTS_1 = 3, EXP_NROOTS_2 = 4
    real :: RCT, RCS, RCQ, RCSOIL, GC, RC, SMCWLT, SMCREF, RSMIN
    real :: EXP_RCT_1 = 9.8559999466E-01, EXP_RCS_1 = 6.3962256908E-01, EXP_RCQ_1 = 7.4187535048E-01, & 
            EXP_RCSOIL_1 = 4.5774650574E-01, EXP_GC_1 = 3.8059044164E-03, EXP_RC_1 = 2.6274963379E+02, &
            EXP_SMCWLT_1 = 2.3000000045E-02, EXP_SMCREF_1 = 2.3600000143E-01, EXP_RSMIN_1 = 225.0
    real :: EXP_RCT_2 = 9.8559999466E-01, EXP_RCS_2 = 8.9405411482E-01, EXP_RCQ_2 = 9.9999997474E-05, & 
            EXP_RCSOIL_2 = 5.5774652958E-01, EXP_GC_2 = 7.9999997979E-04, EXP_RC_2 = 1250.0, &
            EXP_SMCWLT_2 = 2.3000000045E-02, EXP_SMCREF_2 = 2.3600000143E-01, EXP_RSMIN_2 = 100.0

    interface
        subroutine CANRES(SOLAR, SFCTMP, Q2, SFCPRS, SMC, GC, RC, IVEG, ISOIL, RSMIN, &
                          NROOTS, SMCWLT, SMCREF, RCS, RCQ, RCT, RCSOIL, SLDPTH)
            use ctlblk_mod, only: novegtype, nsoil, ivegsrc
            integer , intent(in) :: IVEG, ISOIL
            real, intent(in) :: SOLAR, SFCTMP, Q2, SFCPRS
            real, dimension(nsoil), intent(in) :: SMC, SLDPTH
            integer, intent(out) :: NROOTS
            real, intent(out) :: RCT, RCS, RCQ, RCSOIL, GC, RC, SMCWLT, SMCREF, RSMIN
        end subroutine CANRES
    end interface
    
    res = 0 ! Initialize to no errors

    ivegsrc = 1 ! Test cases where veg type is IGBP
    nsoil = npts
    novegtype = 20
    SOLAR = 600.0
    SFCTMP = 295.0
    Q2 = 0.008
    SFCPRS = 100000.0
    IVEG = 6      ! IGBP: IROOT(6) = 3 -> NROOTS > 1
    ISOIL = 1     ! Soil type: SAND (SMCREF=0.236, SMCWLT=0.023)

    SMC(1) = 0.50 ! GX > 1.0, will be clipped to 1.0
    SMC(2) = 0.00 ! GX < 0.0, will be clipped to 0.0
    SMC(3) = 0.15 ! 0.0 < GX < 1.0
    SMC(4) = 0.15 

    SLDPTH(1) = 0.10
    SLDPTH(2) = 0.30
    SLDPTH(3) = 0.60
    SLDPTH(4) = 0.0
    
    call CANRES(SOLAR, SFCTMP, Q2, SFCPRS, SMC, GC, RC, IVEG, ISOIL, RSMIN, &
                NROOTS, SMCWLT, SMCREF, RCS, RCQ, RCT, RCSOIL, SLDPTH)

    if (NROOTS /= EXP_NROOTS_1) then
        print *, "ERROR: NROOTS expected ", EXP_NROOTS_1, " but got ", NROOTS
        res = 1
    end if
    if (abs(RCT - EXP_RCT_1) > tol) then
        print *, "ERROR: RCT expected ", EXP_RCT_1, " but got ", RCT
        res = 1
    end if
    if (abs(RCS - EXP_RCS_1) > tol) then
        print *, "ERROR: RCS expected ", EXP_RCS_1, " but got ", RCS
        res = 1
    end if
    if (abs(RCQ - EXP_RCQ_1) > tol) then
        print *, "ERROR: RCQ expected ", EXP_RCQ_1, " but got ", RCQ
        res = 1
    end if
    if (abs(RCSOIL - EXP_RCSOIL_1) > tol) then
        print *, "ERROR: RCSOIL expected ", EXP_RCSOIL_1, " but got ", RCSOIL
        res = 1
    end if
    if (abs(GC - EXP_GC_1) > tol) then
        print *, "ERROR: GC expected ", EXP_GC_1, " but got ", GC
        res = 1
    end if
    if (abs(RC - EXP_RC_1) > tol) then
        print *, "ERROR: RC expected ", EXP_RC_1, " but got ", RC
        res = 1
    end if
    if (abs(SMCWLT - EXP_SMCWLT_1) > tol) then
        print *, "ERROR: SMCWLT expected ", EXP_SMCWLT_1, " but got ", SMCWLT
        res = 1
    end if
    if (abs(SMCREF - EXP_SMCREF_1) > tol) then
        print *, "ERROR: SMCREF expected ", EXP_SMCREF_1, " but got ", SMCREF
        res = 1
    end if
    if (abs(RSMIN - EXP_RSMIN_1) > tol) then
        print *, "ERROR: RSMIN expected ", EXP_RSMIN_1, " but got ", RSMIN
        res = 1
    end if

    if (res .ne. 0) stop 10

    ivegsrc = 0 ! Test case where veg type is USGS
    novegtype = 24
    SOLAR = 900.0
    SFCTMP = 301.0
    Q2 = 0.050   
    SFCPRS = 101300.0
    IVEG = 11       ! USGS: IROOT(11) = 4 -> NROOTS > 1
    
    SMC(1) = 0.00 ! GX < 0.0
    SMC(2) = 0.50 ! GX > 1.0
    SMC(3) = 0.15 ! 0.0 < GX < 1.0
    SMC(4) = 0.15 

    SLDPTH(1) = 0.20 
    SLDPTH(2) = 0.20
    SLDPTH(3) = 0.60
    SLDPTH(4) = 0.0

    call CANRES(SOLAR, SFCTMP, Q2, SFCPRS, SMC, GC, RC, IVEG, ISOIL, RSMIN, &
                NROOTS, SMCWLT, SMCREF, RCS, RCQ, RCT, RCSOIL, SLDPTH)

    if (NROOTS /= EXP_NROOTS_2) then
        print *, "ERROR: NROOTS expected ", EXP_NROOTS_2, " but got ", NROOTS
        res = 1
    end if
    if (abs(RCT - EXP_RCT_2) > tol) then
        print *, "ERROR: RCT expected ", EXP_RCT_2, " but got ", RCT
        res = 1
    end if
    if (abs(RCS - EXP_RCS_2) > tol) then
        print *, "ERROR: RCS expected ", EXP_RCS_2, " but got ", RCS
        res = 1
    end if
    if (abs(RCQ - EXP_RCQ_2) > tol) then
        print *, "ERROR: RCQ expected ", EXP_RCQ_2, " but got ", RCQ
        res = 1
    end if
    if (abs(RCSOIL - EXP_RCSOIL_2) > tol) then
        print *, "ERROR: RCSOIL expected ", EXP_RCSOIL_2, " but got ", RCSOIL
        res = 1
    end if
    if (abs(GC - EXP_GC_2) > tol) then
        print *, "ERROR: GC expected ", EXP_GC_2, " but got ", GC
        res = 1
    end if
    if (abs(RC - EXP_RC_2) > tol) then
        print *, "ERROR: RC expected ", EXP_RC_2, " but got ", RC
        res = 1
    end if
    if (abs(SMCWLT - EXP_SMCWLT_2) > tol) then
        print *, "ERROR: SMCWLT expected ", EXP_SMCWLT_2, " but got ", SMCWLT
        res = 1
    end if
    if (abs(SMCREF - EXP_SMCREF_2) > tol) then
        print *, "ERROR: SMCREF expected ", EXP_SMCREF_2, " but got ", SMCREF
        res = 1
    end if
    if (abs(RSMIN - EXP_RSMIN_2) > tol) then
        print *, "ERROR: RSMIN expected ", EXP_RSMIN_2, " but got ", RSMIN
        res = 1
    end if

    if (res .ne. 0) stop 20

    print *, "SUCCESS!"
end program test_canres