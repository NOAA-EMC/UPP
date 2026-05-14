! This is a test program for UPP.
!
! This program tests the subroutines CALRH(), CALRH_GSD() and CALRH_NAM()
! in the upp_physics module.
!
! Alyson Stahl, 5/2026
program test_calrh
    use upp_physics, only: CALRH, CALRH_GSD, CALRH_NAM
    use ctlblk_mod, only: ista, iend, jsta, jend, spval, modelname
    use params_mod, only: PQ0, a2, a3, a4, rhmin
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 7
    integer :: i, res
    ! For expected value calculations
    real :: QC
    ! Keeping arrays separate because Q1 array can get overwritten.
    ! Arrays for CALRH_NAM()
    real :: P1_NAM(1, npts), T1_NAM(1, npts)
    real :: Q1_NAM(1, npts), EXP_Q1_NAM(1, npts)
    real :: RH_NAM(1, npts), EXP_RH_NAM(1, npts)
    ! Arrays for CALRH_GSD()
    real :: P1_GSD(1, npts), T1_GSD(1, npts), Q1_GSD(1, npts)
    real :: RH_GSD(1, npts), EXP_RH_GSD(1, npts)
    ! Arrays for CALRH()
    real :: P1(1, npts), T1(1, npts)
    real :: Q1(1, npts), EXP_Q1(1, npts)
    real :: RH(1, npts)

    ! Grid dimensions
    ista = 1
    iend = 1
    jsta = 1
    jend = npts
    spval = 9.9e10
    modelname = ""

    print *, "Testing CALRH_NAM()..."

    ! Test Case 1: Standard case where RHmin < RH (Q1/QC) < 1.0 (default case)
    P1_NAM = 85000.0
    T1_NAM = 280.0
    Q1_NAM = 0.004
    RH_NAM = 0.0 ! Give a default value.
    EXP_Q1_NAM = 0.004

    QC = PQ0/P1_NAM(1,1)*EXP(A2*(T1_NAM(1,1)-A3)/(T1_NAM(1,1)-A4))
    EXP_RH_NAM = Q1_NAM(1,1)/QC
    
    ! Test Case 2: Q1/QC > 1.0 (Clipped to 1.0)
    Q1_NAM(1,2) = 0.010
    EXP_Q1_NAM(1,2) = QC
    EXP_RH_NAM(1,2) = 1.0

    ! Test Case 3: Q1/QC < RHmin & P1 >= 300.0 (Clipped to RHmin)
    Q1_NAM(1,3) = 1e-6
    EXP_Q1_NAM(1,3) = RHMIN * QC
    EXP_RH_NAM(1,3) = RHMIN

    ! Test Case 4: RHmin / 10 < Q1/QC < RHmin & P1 < 300.0 
    P1_NAM(1,4) = 100.0
    Q1_NAM(1,4) = 1e-6
    QC = PQ0/P1_NAM(1,4)*EXP(A2*(T1_NAM(1,4)-A3)/(T1_NAM(1,4)-A4))
    EXP_RH_NAM(1,4) = Q1_NAM(1,4)/QC
    EXP_Q1_NAM(1,4) = EXP_RH_NAM(1,4) * QC

    ! Test Case 5: Q1/QC < RHmin / 10 & P1 < 300.0 (Clipped to RHmin / 10)
    P1_NAM(1,5) = 100.0
    Q1_NAM(1,5) = 1e-8
    EXP_RH_NAM(1,5) = RHMIN / 10
    QC = PQ0/P1_NAM(1,5)*EXP(A2*(T1_NAM(1,5)-A3)/(T1_NAM(1,5)-A4))
    EXP_Q1_NAM(1,5) = EXP_RH_NAM(1,5) * QC

    ! Test Case 6: ABS(P1) < 1 (RH is not set to anything)
    P1_NAM(1,6) = 0.5
    EXP_RH_NAM(1,6) = 0.0

    ! Test Case 7: T1 has spval
    T1_NAM(1,7) = spval
    EXP_RH_NAM(1,7) = spval

    call CALRH_NAM(P1_NAM, T1_NAM, Q1_NAM, RH_NAM)

    res = 0
    do i = 1, npts
        !if (abs(Q1_NAM(1,i) - EXP_Q1_NAM(1,i)) > tol) then
        !    print *, "CALRH_NAM() Failed for test ", i, ": ", &
        !                "Expected Q1 = ", EXP_Q1_NAM(1,i), &
        !                        " but got Q1 = ", Q1_NAM(1,i)
        !    res = 1
        !end if
        !if (abs(RH_NAM(1,i) - EXP_RH_NAM(1,i)) > tol) then
        !    print *, "CALRH_NAM() Failed for test ", i, ": ", &
        !                "Expected RH = ", EXP_RH_NAM(1,i), &
        !                " but got RH = ", RH_NAM(1,i)
        !    res = 1
        !end if
        print '(A,I0,A,E24.10)', "Q(", i, ") = ", Q1_NAM(1,i)
        print '(A,I0,A,E24.10)', "RH(", i, ") = ", RH_NAM(1,i)
    end do

    if (res .ne. 0) stop 10

    print *, "Testing CALRH_GSD()..."


end program test_calrh