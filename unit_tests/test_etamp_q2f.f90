! This is a test program for UPP.
!
! This program tests the ETAMP_Q2F() subroutine.
!
! Alyson Stahl, 1/2026
program test_etamp_q2f
    use ctlblk_mod, only:   lm,jsta,jend,jsta_2l,jend_2u,&
                            ista,iend,ista_2l,iend_2u
    implicit none

    real, parameter :: tol = 1.0e-8
    ! From etamp_q2f.f
    real, parameter :: t_ice=-40., t0c=273.15, t_icek=233.15, epsq=1.e-12
    integer, parameter :: npts = 2, nlevs = 1
    integer :: i, j, res
    real, dimension(1:npts,1:npts,1:nlevs) :: QRIMEF, QQW, QQR, QQI, T
    real, dimension(1:npts,1:npts,1:nlevs) :: F_RAIN, F_ICE, F_RIMEF, CWM
    real, dimension(1:npts,1:npts,1:nlevs) :: EXP_F_RAIN, EXP_F_ICE, EXP_F_RIMEF, EXP_CWM

    interface 
        subroutine ETAMP_Q2F(QRIMEF, QQI, QQR, QQW, CWM, F_RAIN, F_ICE, F_RIMEF, T)
            use ctlblk_mod, only: lm, jsta_2l, jend_2u, ista_2l, iend_2u
            real, intent(in), dimension(ista_2l:iend_2u,jsta_2l:jend_2u,lm) :: &
                QRIMEF, QQW, QQR, QQI, T
            real, intent(out), dimension(ista_2l:iend_2u,jsta_2l:jend_2u,lm) :: &
                F_RAIN, F_ICE, F_RIMEF, CWM
        end subroutine ETAMP_Q2F
    end interface

    ! Grid parameters
    jsta = 1
    jend = npts
    jsta_2l = 1
    jend_2u = npts
    ista = 1
    iend = npts
    ista_2l = 1
    iend_2u = npts
    lm = nlevs
    
    QRIMEF = 0.015
    QQI = 0.004
    QQR = 0.011
    QQW = 0.003
    T = 260.0

    ! Test Case: F_ICE = 0 & F_FRIMEF = 1 with QQI <= EPSQ
    QQI(1,1,1) = epsq / 10.0

    ! Test Case: F_ICE = 1 & F_FRIMEF = 1 with QQI <= EPSQ
    QQI(1,2,1) = epsq / 10.0
    T(1,2,1) = t_icek - 5.0

    ! Test Case: F_RIMEF = 100.0 with QRIMEF > 100 * QQI &
    ! F_ICE = 1 and F_RAIN = 0 with QQR = 0 and QQW = 0
    QRIMEF(2,1,1) = 0.5
    QQR(2,1,1) = 0.0
    QQW(2,1,1) = 0.0

    ! Test Case: F_RAIN = 1.0 with QQW < 0
    QQW(2,2,1) = -0.0001

    EXP_F_RAIN = reshape([7.8571426868E-01, 0.0, 7.8571426868E-01, 1.0], &
                            [npts, npts, nlevs])
    EXP_F_ICE = reshape([0.0, 1.0, 1.0, 2.6845636964E-01], &
                            [npts, npts, nlevs])
    EXP_F_RIMEF = reshape([1.0, 100.0, 1.0, 3.7499997616], & 
                            [npts, npts, nlevs])
    EXP_CWM = reshape([1.4000000432E-02, 4.0000001900E-03, &
                        1.4000000432E-02, 1.4900000766E-02], &
                        [npts, npts, nlevs])
    
    call ETAMP_Q2F(QRIMEF, QQI, QQR, QQW, CWM, F_RAIN, F_ICE, F_RIMEF, T)

    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(F_RAIN(i,j,1) - EXP_F_RAIN(i,j,1)) > tol) then
                print *, "F_RAIN test failed at (", i, ",", j, "): ", &
                         "Expected ", EXP_F_RAIN(i,j,1), " but got ", F_RAIN(i,j,1)
                res = 1
            end if
            if (abs(F_ICE(i,j,1) - EXP_F_ICE(i,j,1)) > tol) then
                print *, "F_ICE test failed at (", i, ",", j, "): ", &
                         "Expected ", EXP_F_ICE(i,j,1), " but got ", F_ICE(i,j,1)
                res = 1
            end if
            if (abs(F_RIMEF(i,j,1) - EXP_F_RIMEF(i,j,1)) > tol) then
                print *, "F_RIMEF test failed at (", i, ",", j, "): ", &
                         "Expected ", EXP_F_RIMEF(i,j,1), " but got ", F_RIMEF(i,j,1)
                res = 1
            end if
            if (abs(CWM(i,j,1) - EXP_CWM(i,j,1)) > tol) then
                print *, "CWM test failed at (", i, ",", j, "): ", &
                         "Expected ", EXP_CWM(i,j,1), " but got ", CWM(i,j,1)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10
    print *, "SUCCESS!"
end program test_etamp_q2f
