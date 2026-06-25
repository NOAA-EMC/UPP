! This is a test program for UPP.
!
! This program tests the subroutines CALRH(), CALRH_GSD() and CALRH()
! in the upp_physics module.
!
! Alyson Stahl, 5/2026
program test_calrh
    use upp_physics, only: CALRH_NAM, CALRH_GSD, CALRH
    use ctlblk_mod, only: ista, iend, jsta, jend, spval, modelname
    use params_mod, only: PQ0, a2, a3, a4, rhmin
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: npts = 9
    integer :: i, res
    ! For expected value calculations
    real :: QC
    real :: P1(1, npts), T1(1, npts), RH(1, npts)
    ! Keeping Q1 arrays separate because it can get overwritten.
    real :: Q1_NAM(1, npts), Q1_GSD(1, npts), Q1(1, npts)
    real :: EXP_Q1_NAM(1, npts), EXP_Q1_GSD(1, npts)
    real :: EXP_RH_NAM(1, npts), EXP_RH_GSD(1, npts)

    ! Grid dimensions
    ista = 1
    iend = 1
    jsta = 1
    jend = npts
    spval = 9.9e10
    modelname = "RAPR"

    ! Test Cases completely cover CALRH_NAM() and CALRH_GSD(). CALRH() just calls
    ! CALRH_NAM() or CALRH_GSD() depending on the model, so the tests just check it 
    ! calls the correct one. No need to construct separate inputs.

    print *, "Testing CALRH_NAM()..."

    ! Test Case 1: Standard case where RHmin < RH (Q1/QC) < 1.0 (default case)
    P1 = 85000.0
    T1 = 280.0
    Q1_NAM = 0.004
    RH = 0.0 ! Give a default value.

    EXP_Q1_NAM = 0.004
    QC = PQ0/P1(1,1)*EXP(A2*(T1(1,1)-A3)/(T1(1,1)-A4))
    EXP_RH_NAM = Q1_NAM(1,1)/QC
    
    ! Test Case 2: Q1/QC > 1.0 (Clipped to 1.0)
    Q1_NAM(1,2) = 0.010
    EXP_Q1_NAM(1,2) = QC
    EXP_RH_NAM(1,2) = 1.0

    ! Test Case 3: Q1/QC < RHmin & P1 >= 300.0 (Clipped to RHmin)
    Q1_NAM(1,3) = 1e-10
    EXP_Q1_NAM(1,3) = RHMIN * QC
    EXP_RH_NAM(1,3) = RHMIN

    ! Test Case 4: RHmin / 10 < Q1/QC < RHmin & P1 < 300.0 
    P1(1,4) = 100.0
    Q1_NAM(1,4) = 1e-6
    QC = PQ0/P1(1,4)*EXP(A2*(T1(1,4)-A3)/(T1(1,4)-A4))
    EXP_RH_NAM(1,4) = Q1_NAM(1,4)/QC
    EXP_Q1_NAM(1,4) = EXP_RH_NAM(1,4) * QC

    ! Test Case 5: Q1/QC < RHmin / 10 & P1 < 300.0 (Clipped to RHmin / 10)
    P1(1,5) = 100.0
    Q1_NAM(1,5) = 1e-8
    EXP_RH_NAM(1,5) = RHMIN / 10
    QC = PQ0/P1(1,5)*EXP(A2*(T1(1,5)-A3)/(T1(1,5)-A4))
    EXP_Q1_NAM(1,5) = EXP_RH_NAM(1,5) * QC

    ! Test Case 6: ABS(P1) < 1 (RH is not set to anything)
    P1(1,6) = 0.5
    EXP_RH_NAM(1,6) = 0.0

    ! Test Case 7: T1 has spval
    T1(1,7) = spval
    EXP_RH_NAM(1,7) = spval

    ! Test Case 8: P1 has spval (not explicitly handled by CALRH_NAM)
    P1(1,8) = spval
    QC = PQ0/P1(1,8)*EXP(A2*(T1(1,8)-A3)/(T1(1,8)-A4))
    EXP_RH_NAM(1,8) = 1.0
    EXP_Q1_NAM(1,8) = QC

    ! Test Case 9: Q1 has spval (not explicitly handled by CALRH_NAM)
    Q1_NAM(1,9) = spval
    QC = PQ0/P1(1,9)*EXP(A2*(T1(1,9)-A3)/(T1(1,9)-A4))
    EXP_RH_NAM(1,9) = 1.0
    EXP_Q1_NAM(1,9) = QC

    ! Copy Q1_NAM to other Q1 arrays before it gets overwritten.
    do i = 1, npts
        Q1(1,i) = Q1_NAM(1,i)
        Q1_GSD(1,i) = Q1_NAM(1,i)

        ! Q1 not overwritten in CALRH_GSD(), but copy it to EXP_Q1_GSD for consistency.
        EXP_Q1_GSD(1,i) = Q1_GSD(1,i)
    end do

    call CALRH_NAM(P1, T1, Q1_NAM, RH)

    res = 0
    do i = 1, npts
        if (abs(Q1_NAM(1,i) - EXP_Q1_NAM(1,i)) > tol) then
            print *, "CALRH_NAM() Failed for test ", i, ": ", &
                        "Expected Q1 = ", EXP_Q1_NAM(1,i), &
                                " but got Q1 = ", Q1_NAM(1,i)
            res = 1
        end if
        if (abs(RH(1,i) - EXP_RH_NAM(1,i)) > tol) then
            print *, "CALRH_NAM() Failed for test ", i, ": ", &
                        "Expected RH = ", EXP_RH_NAM(1,i), &
                                " but got RH = ", RH(1,i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10

    print *, "Testing CALRH_GSD()..."

    RH = 0.0 ! Reset

    EXP_RH_GSD(1,1) = 5.5023592710E-01
    EXP_RH_GSD(1,2) = 1.0000000000E+00
    EXP_RH_GSD(1,3) = 1.3789341224E-08
    EXP_RH_GSD(1,4) = 1.6222745103E-07
    EXP_RH_GSD(1,5) = 1.6222753141E-09
    EXP_RH_GSD(1,6) = 3.2366820051E-06
    EXP_RH_GSD(1,7) = spval
    EXP_RH_GSD(1,8) = spval
    EXP_RH_GSD(1,9) = spval

    call CALRH_GSD(P1, T1, Q1_GSD, RH)

    res = 0
    do i = 1, npts
        if (abs(Q1_GSD(1,i) - EXP_Q1_GSD(1,i)) > tol) then
            print *, "CALRH_GSD() Failed for test ", i, ": ", &
                        "Expected Q1 = ", EXP_Q1_GSD(1,i), &
                                " but got Q1 = ", Q1_GSD(1,i)
            res = 1
        end if
        if (abs(RH(1,i) - EXP_RH_GSD(1,i)) > tol) then
            print *, "CALRH_GSD() Failed for test ", i, ": ", &
                        "Expected RH = ", EXP_RH_GSD(1,i), &
                                " but got RH = ", RH(1,i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 20

    print *, "Testing CALRH()..."

    ! Test Case 1: modelname == 'RAPR' (Calls CALRH_GSD())
    RH = 0.0 ! Reset
    call CALRH(P1, T1, Q1, RH)

    res = 0
    do i = 1, npts
        if (abs(Q1(1,i) - EXP_Q1_GSD(1,i)) > tol) then
            print *, "CALRH() Failed for test ", i, ": ", &
                        "Expected Q1 = ", EXP_Q1_GSD(1,i), &
                                " but got Q1 = ", Q1(1,i)
            res = 1
        end if
        if (abs(RH(1,i) - EXP_RH_GSD(1,i)) > tol) then
            print *, "CALRH() Failed for test ", i, ": ", &
                        "Expected RH = ", EXP_RH_GSD(1,i), &
                                " but got RH = ", RH(1,i)
            res = 1
        end if
    end do

    if (res .ne. 0) then
        print *, "CALRH() failed for modelname == 'RAPR'."
        stop 30
    end if

    ! Test Case 2: modelname != 'RAPR' (Calls CALRH_NAM())
    modelname = "Not_RAPR"
    RH = 0.0 ! Reset
    call CALRH(P1, T1, Q1, RH)

    res = 0
    do i = 1, npts
        if (abs(Q1(1,i) - EXP_Q1_NAM(1,i)) > tol) then
            print *, "CALRH() Failed for test ", i, ": ", &
                        "Expected Q1 = ", EXP_Q1_NAM(1,i), &
                                " but got Q1 = ", Q1(1,i)
            res = 1
        end if
        if (abs(RH(1,i) - EXP_RH_NAM(1,i)) > tol) then
            print *, "CALRH() Failed for test ", i, ": ", &
                        "Expected RH = ", EXP_RH_NAM(1,i), &
                                " but got RH = ", RH(1,i)
            res = 1
        end if
    end do

    if (res .ne. 0) then
        print *, "CALRH() failed for modelname != 'RAPR'."
        stop 40
    end if

    print *, "SUCCESS!"
end program test_calrh