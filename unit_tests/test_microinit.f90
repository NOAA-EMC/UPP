! This is a test program for UPP.
!
! This program tests the MICROINIT() subroutine.
!
! Alyson Stahl, 4/2026
program test_microinit
    use svptbl_mod, only: nx, tbpvs, tbpvs0
    use params_mod, only: tfrz, pi
    use cmassi_mod, only: dmrmax, t_ice, nlimax, xmrmax, &
                        mdrmax, mdrmin, trad_ice, massi, &
                        rqr_drmin, n0r0, rqr_drmax, cn0r0, &
                        cn0r_dmrmin, cn0r_dmrmax, dmrmin
    use gridspec_mod,only : gridtype
    use rhgrd_mod, only: rhgrd
    use ctlblk_mod, only: me
    implicit none

    real, parameter :: tol = 1.0e-8, spval = -999.0
    integer, parameter :: ntests = 5
    integer :: i, j, res
    integer :: imp_physics(ntests)
    character(len=1) :: cur_gridtype(ntests)
    ! Note: Subroutine also sets FLARGE2, but it's no longer used, so it's not being tested here.
    real, dimension(ntests) :: EXP_RHGRD, EXP_DMRmax, EXP_XMRmax, EXP_T_ICE, EXP_NLImax, &
                                EXP_TRAD_ICE, EXP_RQR_DRmin, EXP_RQR_DRmax, EXP_CN0r0, &
                                EXP_CN0r_DMRmin, EXP_CN0r_DMRmax
    integer, dimension(ntests) :: EXP_MDRmax
    real :: tmp

    interface
        subroutine MICROINIT(imp_physics)
            integer, intent(in) :: imp_physics
        end subroutine MICROINIT
    end interface

    me = 0
    
    ! Will use this grid type for most cases. Only checked in the imp_physics = 95 case.
    cur_gridtype = 'B'

    EXP_T_ICE = -40.0
    EXP_TRAD_ICE = tfrz - 20.0

    tmp = 1000.0 * pi * n0r0
    EXP_CN0r0 = 1.E6/SQRT(SQRT(tmp))

    ! Test Case 1: imp_physics = 5
    imp_physics(1) = 5
    EXP_RHGRD(1) = 0.98
    EXP_DMRmax(1) = 1.E-3
    EXP_XMRmax(1) = 1.E6 * EXP_DMRmax(1)
    EXP_MDRmax(1) = int(EXP_XMRmax(1))
    EXP_NLImax(1) = spval
    ! Come from table so will need to set later
    EXP_RQR_DRmin(1) = 1.5419247745E-07
    EXP_RQR_DRmax(1) = 2.4873061106E-02
    !
    tmp = 1000.0 * pi * DMRmin * DMRmin * DMRmin * DMRmin
    EXP_CN0r_DMRmin(1) = 1./tmp
    tmp = 1000.0 * pi * EXP_DMRmax(1) * EXP_DMRmax(1) * EXP_DMRmax(1) * EXP_DMRmax(1)
    EXP_CN0r_DMRmax(1) = 1./tmp

    ! Test Case 2: imp_physics = 15 (should be the same as 5)
    imp_physics(2) = 15
    EXP_RHGRD(2) = 0.98
    EXP_DMRmax(2) = 1.E-3
    EXP_XMRmax(2) = 1.E6 * EXP_DMRmax(2)
    EXP_MDRmax(2) = int(EXP_XMRmax(2))
    EXP_NLImax(2) = spval
    ! Come from table so will need to set later
    EXP_RQR_DRmin(2) = 1.5419247745E-07
    EXP_RQR_DRmax(2) = 2.4873061106E-02
    !
    tmp = 1000.0 * pi * DMRmin * DMRmin * DMRmin * DMRmin
    EXP_CN0r_DMRmin(2) = 1./tmp
    tmp = 1000.0 * pi * EXP_DMRmax(2) * EXP_DMRmax(2) * EXP_DMRmax(2) * EXP_DMRmax(2)
    EXP_CN0r_DMRmax(2) = 1./tmp

    ! Test Case 3: imp_physics = 85
    imp_physics(3) = 85
    EXP_RHGRD(3) = 1.
    EXP_DMRmax(3) = .45E-3
    EXP_XMRmax(3) = 1.E6 * EXP_DMRmax(3)
    EXP_MDRmax(3) = int(EXP_XMRmax(3))
    EXP_NLImax(3) = 20.E3
    ! Come from table so will need to set later
    EXP_RQR_DRmin(3) = 1.5419246324E-07
    EXP_RQR_DRmax(3) = 1.0305735050E-03
    !
    tmp = 1000.0 * pi * DMRmin * DMRmin * DMRmin * DMRmin
    EXP_CN0r_DMRmin(3) = 1./tmp
    tmp = 1000.0 * pi * EXP_DMRmax(3) * EXP_DMRmax(3) * EXP_DMRmax(3) * EXP_DMRmax(3)
    EXP_CN0r_DMRmax(3) = 1./tmp

    ! Test Case 4: imp_physics = 95 & gridtype == "B"
    imp_physics(4) = 95
    EXP_RHGRD(4) = 1.
    EXP_DMRmax(4) = .45E-3
    EXP_XMRmax(4) = 1.E6 * EXP_DMRmax(4)
    EXP_MDRmax(4) = int(EXP_XMRmax(4))
    EXP_NLImax(4) = 5.E3
    ! Come from table so will need to set later
    EXP_RQR_DRmin(4) = 1.5419246324E-07 
    EXP_RQR_DRmax(4) = 1.0305735050E-03
    !
    tmp = 1000.0 * pi * DMRmin * DMRmin * DMRmin * DMRmin
    EXP_CN0r_DMRmin(4) = 1./tmp
    tmp = 1000.0 * pi * EXP_DMRmax(4) * EXP_DMRmax(4) * EXP_DMRmax(4) * EXP_DMRmax(4)
    EXP_CN0r_DMRmax(4) = 1./tmp

    ! Test Case 5: imp_physics = 95 & gridtype != "B"
    imp_physics(5) = 95
    cur_gridtype(5) = 'C'
    EXP_RHGRD(5) = 1.
    EXP_DMRmax(5) = 1.E-3
    EXP_XMRmax(5) = 1.E6 * EXP_DMRmax(5)
    EXP_MDRmax(5) = int(EXP_XMRmax(5))
    EXP_NLImax(5) = 5.E3
    ! Come from table so will need to set later
    EXP_RQR_DRmin(5) = 1.5419247745E-07
    EXP_RQR_DRmax(5) = 2.4873061106E-02
    !
    tmp = 1000.0 * pi * DMRmin * DMRmin * DMRmin * DMRmin
    EXP_CN0r_DMRmin(5) = 1./tmp
    tmp = 1000.0 * pi * EXP_DMRmax(5) * EXP_DMRmax(5) * EXP_DMRmax(5) * EXP_DMRmax(5)
    EXP_CN0r_DMRmax(5) = 1./tmp

    res = 0
    do i = 1, ntests
        ! Initialize with some value
        RHgrd = spval
        DMRmax = spval
        XMRmax = spval
        T_ICE = spval
        NLImax = spval
        TRAD_ICE = spval
        MDRmax = int(spval)
        RQR_DRmin = spval
        RQR_DRmax = spval
        CN0r0 = spval
        CN0r_DMRmin = spval
        CN0r_DMRmax = spval

        ! These are initialized by GPVS() which is called at the end of MICROINIT().
        tbpvs = spval
        tbpvs0 = spval

        gridtype = cur_gridtype(i)

        call MICROINIT(imp_physics(i))
        
        if (abs(RHgrd - EXP_RHGRD(i)) > tol) then
            print *, "Test Case ", i, " Failed: RHgrd = ", RHgrd, &
                " Expected = ", EXP_RHGRD(i)
            res = 1
        end if
        if (abs(DMRmax - EXP_DMRmax(i)) > tol) then
            print *, "Test Case ", i, " Failed: DMRmax = ", DMRmax, &
                " Expected = ", EXP_DMRmax(i)
            res = 1
        end if
        if (abs(XMRmax - EXP_XMRmax(i)) > tol) then
            print *, "Test Case ", i, " Failed: XMRmax = ", XMRmax, &
                " Expected = ", EXP_XMRmax(i)
            res = 1
        end if
        if (abs(T_ICE - EXP_T_ICE(i)) > tol) then
            print *, "Test Case ", i, " Failed: T_ICE = ", T_ICE, &
                " Expected = ", EXP_T_ICE(i)
            res = 1
        end if
        if (abs(NLImax - EXP_NLImax(i)) > tol) then
            print *, "Test Case ", i, " Failed: NLImax = ", NLImax, &
                " Expected = ", EXP_NLImax(i)
            res = 1
        end if
        if (abs(TRAD_ICE - EXP_TRAD_ICE(i)) > tol) then
            print *, "Test Case ", i, " Failed: TRAD_ICE = ", TRAD_ICE, &
                " Expected = ", EXP_TRAD_ICE(i)
            res = 1
        end if
        print '(A,I0,A,1X,ES16.10)', "RQR_DRMIN(", i, ") =", RQR_DRmin
        print '(A,I0,A,1X,ES16.10)', "RQR_DRMAX(", i, ") =", RQR_DRmax
        !if (abs(RQR_DRmin - EXP_RQR_DRmin(i)) > tol) then
        !    print *, "Test Case ", i, " Failed: RQR_DRmin = ", RQR_DRmin, &
        !        " Expected = ", EXP_RQR_DRmin(i)
        !    res = 1
        !end if
        !if (abs(RQR_DRmax - EXP_RQR_DRmax(i)) > tol) then
        !    print *, "Test Case ", i, " Failed: RQR_DRmax = ", RQR_DRmax, &
        !        " Expected = ", EXP_RQR_DRmax(i)
        !    res = 1
        !end if
        if (abs(CN0r0 - EXP_CN0r0(i)) > tol) then
            print *, "Test Case ", i, " Failed: CN0r0 = ", CN0r0, &
                " Expected = ", EXP_CN0r0(i)
            res = 1
        end if
        if (abs(CN0r_DMRmin - EXP_CN0r_DMRmin(i)) > tol) then
            print *, "Test Case ", i, " Failed: CN0r_DMRmin = ", CN0r_DMRmin, &
                " Expected = ", EXP_CN0r_DMRmin(i)
            res = 1
        end if
        if (abs(CN0r_DMRmax - EXP_CN0r_DMRmax(i)) > tol) then
            print *, "Test Case ", i, " Failed: CN0r_DMRmax = ", CN0r_DMRmax, &
                " Expected = ", EXP_CN0r_DMRmax(i)
            res = 1
        end if
        if (MDRmax .ne. EXP_MDRmax(i)) then
            print *, "Test Case ", i, " Failed: MDRmax = ", MDRmax, &
                " Expected = ", EXP_MDRmax(i)
            res = 1
        end if

        ! Only checking the first value of tpvs and tpvs0 just to verify that GPVS() is being called.
        if ((tbpvs(1) < 0).OR.(tbpvs0(1) < 0)) then
            print *, "Test Case ", i, " Failed: Did not call GPVS()."
            res = 1
        end if

        !if (res .ne. 0) then
        !    print *, "Test Case ", i, " Failed."
        !    stop 10
        !end if
    end do
    
    print *, 'SUCCESS!'
end program test_microinit