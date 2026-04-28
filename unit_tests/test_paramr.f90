! This is a test program for UPP.
!
! This program tests the PARAMR() subroutine.
!
! The PARAMR() subroutine sets all of the microphysics constants in the PMICRPH_mod module (100+ in total).
! The PMICRPH_mod module is only referenced by the MDLFLD() subroutine in UPP source code and only a fraction of 
! all of the parameters are used. This test only checks those parameters used by MDLFLD(). However, parameters
! can always be added to this test as needed.
!
! Alyson Stahl, 4/2026
program test_paramr
    use params_mod, only: pi
    use pmicrph_mod, only: r1, const1r, qr0, delqr0, const2r, ron, topr, son,&
            tops, dsnow, drain,const_ng1, const_ng2, gon, topg, dgraupel

    implicit none

    real, parameter :: tol = 1.0e-8
    integer :: res
    real :: EXP_r1 = 1.E-15
    real :: EXP_const1r = 4.9E8
    real :: EXP_qr0 = 0.0002
    real :: EXP_delqr0 = 0.0001
    real :: EXP_const2r = 5.1E8
    real :: EXP_ron = 8.E6
    real :: EXP_topr = pi * 8.E9
    real :: EXP_son = 2.E7
    real :: EXP_tops = pi * 2.e9
    real :: EXP_dsnow = 100.
    real :: EXP_drain = 1000.
    real :: EXP_const_ng1 = (1.57**(1./0.52))*((400.*pi)**(12./13.))
    real :: EXP_const_ng2 = -(12./13.)
    real :: EXP_gon = 4.E6
    real :: EXP_topg = pi * 1.6E9
    real :: EXP_dgraupel = 400.
    interface
        subroutine PARAMR()
        end subroutine PARAMR
    end interface

    call PARAMR()

    res = 0
    if (abs(r1 - EXP_r1) > tol) then
        print *, 'Error: r1 = ', r1, ' does not match expected value of ', EXP_r1
        res = 1
    end if
    if (abs(const1r - EXP_const1r) > tol) then
        print *, 'Error: const1r = ', const1r, ' does not match expected value of ', EXP_const1r
        res = 1
    end if
    if (abs(qr0 - EXP_qr0) > tol) then
        print *, 'Error: qr0 = ', qr0, ' does not match expected value of ', EXP_qr0
        res = 1
    end if
    if (abs(delqr0 - EXP_delqr0) > tol) then
        print *, 'Error: delqr0 = ', delqr0, ' does not match expected value of ', EXP_delqr0
        res = 1
    end if
    if (abs(const2r - EXP_const2r) > tol) then
        print *, 'Error: const2r = ', const2r, ' does not match expected value of ', EXP_const2r
        res = 1
    end if
    if (abs(ron - EXP_ron) > tol) then
        print *, 'Error: ron = ', ron, ' does not match expected value of ', EXP_ron
        res = 1
    end if
    if (abs(topr - EXP_topr) > tol) then
        print *, 'Error: topr = ', topr, ' does not match expected value of ', EXP_topr
        res = 1
    end if
    if (abs(son - EXP_son) > tol) then
        print *, 'Error: son = ', son, ' does not match expected value of ', EXP_son
        res = 1
    end if
    if (abs(tops - EXP_tops) > tol) then
        print *, 'Error: tops = ', tops, ' does not match expected value of ', EXP_tops
        res = 1
    end if
    if (abs(dsnow - EXP_dsnow) > tol) then
        print *, 'Error: dsnow = ', dsnow, ' does not match expected value of ', EXP_dsnow
        res = 1
    end if
    if (abs(drain - EXP_drain) > tol) then
        print *, 'Error: drain = ', drain, ' does not match expected value of ', EXP_drain
        res = 1
    end if
    if (abs(const_ng1 - EXP_const_ng1) > tol) then
        print *, 'Error: const_ng1 = ', const_ng1, ' does not match expected value of ', EXP_const_ng1
        res = 1
    end if
    if (abs(const_ng2 - EXP_const_ng2) > tol) then
        print *, 'Error: const_ng2 = ', const_ng2, ' does not match expected value of ', EXP_const_ng2
        res = 1
    end if
    if (abs(gon - EXP_gon) > tol) then
        print *, 'Error: gon = ', gon, ' does not match expected value of ', EXP_gon
        res = 1
    end if
    if (abs(topg - EXP_topg) > tol) then
        print *, 'Error: topg = ', topg, ' does not match expected value of ', EXP_topg
        res = 1
    end if
    if (abs(dgraupel - EXP_dgraupel) > tol) then
        print *, 'Error: dgraupel = ', dgraupel, ' does not match expected value of ', EXP_dgraupel
        res = 1
    end if

    if (res .ne. 0) stop 10

    print *, 'SUCCESS!'
end program test_paramr