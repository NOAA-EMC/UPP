! This is a test program for UPP.
!
! This program tests the TABLE() subroutine.
!
! Alyson Stahl, 6/2026
program test_table
    use table_upp_mod, only: TABLE
    implicit none

    real, parameter :: tol = 1.0e-6, rel_tol = 1.0e-6
    integer, parameter :: ITB=076, JTB=134, ntests = 3
    integer :: i, j, res
    !
    character(len=*), parameter :: data_file_prefix = 'data/ref_table_case'
    character(len=*), parameter :: data_file_suffix = '.txt'
    !
    real :: PT, THL
    real :: PTBL(ITB,JTB), TTBL(JTB,ITB)
    real :: QS0(JTB), SQS(JTB), STHE(ITB), THE0(ITB)
    real :: RDQ, RDTH, RDP, RDTHE, PL
    !
    real :: EXP_PTBL(ITB,JTB,ntests), EXP_TTBL(JTB,ITB,ntests)
    real :: EXP_QS0(JTB,ntests), EXP_SQS(JTB,ntests), EXP_STHE(ITB,ntests), EXP_THE0(ITB,ntests)
    real :: EXP_RDQ(ntests), EXP_RDTH(ntests), EXP_RDP(ntests), EXP_RDTHE(ntests), EXP_PL(ntests)

    ! Load expected data.
    do i = 1, ntests
        call load_reference_data(i, EXP_PTBL(:,:,i), EXP_TTBL(:,:,i), EXP_QS0(:,i), &
                                  EXP_SQS(:,i), EXP_STHE(:,i), EXP_THE0(:,i))
    end do

    ! PL = PT
    ! RDQ = KPM - 1 = ITB - 1 = 75
    ! RDTH = 1./DTH -> DTH = (THH - THL) / REAL(KTHM - 1) = (THH - THL) / REAL(JTB - 1) = (365 - THL) / 133
    ! RDP = 1./DP -> DP = (PH - PL) / REAL(KPM - 1) = (PH - PL) / REAL(ITB - 1) = (105000.- PL) / 75
    ! RDTHE = 1./DTHE -> DTHE = 1 / REAL(KTHM - 1) = 1 / (JTB - 1) = 1 / 133

    ! Test Case 1: Standard case. Same input values used in INITPOST_GFS_NEMS_MPIIO().
    PT = 10000.0
    THL = 210.0

    EXP_RDQ(1) = 75.0
    EXP_RDTH(1) = 133.0 / 155.0
    EXP_RDP(1) = 75.0 / 95000.0
    EXP_RDTHE(1) = 1.0 / (1.0 / 133.0) ! Prevents precision related errors
    EXP_PL(1) = PT

    print *, 'Running Test Case 1: PT = 10000.0, THL = 210.0'
    call TABLE(PTBL, TTBL, PT, RDQ, RDTH, RDP, RDTHE, PL, THL, QS0, SQS, STHE, THE0)

    res = 0

    if (abs(RDQ - EXP_RDQ(1)) > tol) then
        print *, 'Test Case 1 Failed: RDQ = ', RDQ, ' Expected: ', EXP_RDQ(1)
        res = 1
    end if
    if (abs(RDTH - EXP_RDTH(1)) > tol) then
        print *, 'Test Case 1 Failed: RDTH = ', RDTH, ' Expected: ', EXP_RDTH(1)
        res = 1
    end if
    if (abs(RDP - EXP_RDP(1)) > tol) then
        print *, 'Test Case 1 Failed: RDP = ', RDP, ' Expected: ', EXP_RDP(1)
        res = 1
    end if
    if (abs(RDTHE - EXP_RDTHE(1)) > tol) then
        print *, 'Test Case 1 Failed: RDTHE = ', RDTHE, ' Expected: ', EXP_RDTHE(1)
        res = 1
    end if
    if (abs(PL - EXP_PL(1)) > tol) then
        print *, 'Test Case 1 Failed: PL = ', PL, ' Expected: ', EXP_PL(1)
        res = 1
    end if

    do i = 1, ITB
        do j = 1, JTB
            if (abs(PTBL(i,j) - EXP_PTBL(i,j,1)) / EXP_PTBL(i,j,1) > rel_tol) then
                print '(A,I0,A,I0,A,ES24.16,A,ES24.16)', 'Test Case 1 Failed: PTBL(', i, ',', j, ') = ', PTBL(i,j), &
                         ' Expected: ', EXP_PTBL(i,j,1)
                res = 1
            end if
            if (abs(TTBL(j,i) - EXP_TTBL(j,i,1)) / EXP_TTBL(j,i,1) > rel_tol) then
                print '(A,I0,A,I0,A,ES24.16,A,ES24.16)', 'Test Case 1 Failed: TTBL(', j, ',', i, ') = ', TTBL(j,i), &
                         ' Expected: ', EXP_TTBL(j,i,1)
                res = 1
            end if
        end do
        if (abs(STHE(i) - EXP_STHE(i,1)) / EXP_STHE(i,1) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 1 Failed: STHE(', i, ') = ', STHE(i), ' Expected: ', EXP_STHE(i,1)
            res = 1
        end if
        if (abs(THE0(i) - EXP_THE0(i,1)) / EXP_THE0(i,1) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 1 Failed: THE0(', i, ') = ', THE0(i), ' Expected: ', EXP_THE0(i,1)
            res = 1
        end if
    end do

    do j = 1, JTB
        if (abs(QS0(j) - EXP_QS0(j,1)) / EXP_QS0(j,1) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 1 Failed: QS0(', j, ') = ', QS0(j), ' Expected: ', EXP_QS0(j,1)
            res = 1
        end if
        if (abs(SQS(j) - EXP_SQS(j,1)) / EXP_SQS(j,1) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 1 Failed: SQS(', j, ') = ', SQS(j), ' Expected: ', EXP_SQS(j,1)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10

    ! Test Case 2: PT = 0.0 (reaches the p <= 0.0 branch)
    PT = 0.0
    THL = 210.0

    EXP_RDQ(2) = 75.0
    EXP_RDTH(2) = 133.0 / 155.0
    EXP_RDP(2) = 75.0 / 105000.0 
    EXP_RDTHE(2) = 1.0 / (1.0 / 133.0) ! Prevents precision related errors
    EXP_PL(2) = PT

    print *, 'Running Test Case 2: PT = 0.0, THL = 210.0'
    call TABLE(PTBL, TTBL, PT, RDQ, RDTH, RDP, RDTHE, PL, THL, QS0, SQS, STHE, THE0)

    res = 0

    if (abs(RDQ - EXP_RDQ(2)) > tol) then
        print *, 'Test Case 2 Failed: RDQ = ', RDQ, ' Expected: ', EXP_RDQ(2)
        res = 1
    end if
    if (abs(RDTH - EXP_RDTH(2)) > tol) then
        print *, 'Test Case 2 Failed: RDTH = ', RDTH, ' Expected: ', EXP_RDTH(2)
        res = 1
    end if
    if (abs(RDP - EXP_RDP(2)) > tol) then
        print *, 'Test Case 2 Failed: RDP = ', RDP, ' Expected: ', EXP_RDP(2)
        res = 1
    end if
    if (abs(RDTHE - EXP_RDTHE(2)) > tol) then
        print *, 'Test Case 2 Failed: RDTHE = ', RDTHE, ' Expected: ', EXP_RDTHE(2)
        res = 1
    end if
    if (abs(PL - EXP_PL(2)) > tol) then
        print *, 'Test Case 2 Failed: PL = ', PL, ' Expected: ', EXP_PL(2)
        res = 1
    end if

    do i = 1, ITB
        do j = 1, JTB
            if (abs(PTBL(i,j) - EXP_PTBL(i,j,2)) / EXP_PTBL(i,j,2) > rel_tol) then
                print '(A,I0,A,I0,A,ES24.16,A,ES24.16)', 'Test Case 2 Failed: PTBL(', i, ',', j, ') = ', PTBL(i,j), &
                         ' Expected: ', EXP_PTBL(i,j,2)
                res = 1
            end if
            if (abs(TTBL(j,i) - EXP_TTBL(j,i,2)) / EXP_TTBL(j,i,2) > rel_tol) then
                print '(A,I0,A,I0,A,ES24.16,A,ES24.16)', 'Test Case 2 Failed: TTBL(', j, ',', i, ') = ', TTBL(j,i), &
                         ' Expected: ', EXP_TTBL(j,i,2)
                res = 1
            end if
        end do
        if (abs(STHE(i) - EXP_STHE(i,2)) / EXP_STHE(i,2) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 2 Failed: STHE(', i, ') = ', STHE(i), ' Expected: ', EXP_STHE(i,2)
            res = 1
        end if
        if (abs(THE0(i) - EXP_THE0(i,2)) / EXP_THE0(i,2) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 2 Failed: THE0(', i, ') = ', THE0(i), ' Expected: ', EXP_THE0(i,2)
            res = 1
        end if
    end do

    do j = 1, JTB
        if (abs(QS0(j) - EXP_QS0(j,2)) / EXP_QS0(j,2) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 2 Failed: QS0(', j, ') = ', QS0(j), ' Expected: ', EXP_QS0(j,2)
            res = 1
        end if
        if (abs(SQS(j) - EXP_SQS(j,2)) / EXP_SQS(j,2) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 2 Failed: SQS(', j, ') = ', SQS(j), ' Expected: ', EXP_SQS(j,2)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 20

    ! Test Case 3: Low Pressure (reaches DENOM <= EPS branch)
    PT = 100.0
    THL = 210.0

    EXP_RDQ(3) = 75.0
    EXP_RDTH(3) = 133.0 / 155.0
    EXP_RDP(3) = 75.0 / 104900.0
    EXP_RDTHE(3) = 1.0 / (1.0 / 133.0) ! Prevents precision related errors
    EXP_PL(3) = PT

    print *, 'Running Test Case 3: PT = 100.0, THL = 210.0'
    call TABLE(PTBL, TTBL, PT, RDQ, RDTH, RDP, RDTHE, PL, THL, QS0, SQS, STHE, THE0)

    res = 0

    if (abs(RDQ - EXP_RDQ(3)) > tol) then
        print *, 'Test Case 3 Failed: RDQ = ', RDQ, ' Expected: ', EXP_RDQ(3)
        res = 1
    end if
    if (abs(RDTH - EXP_RDTH(3)) > tol) then
        print *, 'Test Case 3 Failed: RDTH = ', RDTH, ' Expected: ', EXP_RDTH(3)
        res = 1
    end if
    if (abs(RDP - EXP_RDP(3)) > tol) then
        print *, 'Test Case 3 Failed: RDP = ', RDP, ' Expected: ', EXP_RDP(3)
        res = 1
    end if
    if (abs(RDTHE - EXP_RDTHE(3)) > tol) then
        print *, 'Test Case 3 Failed: RDTHE = ', RDTHE, ' Expected: ', EXP_RDTHE(3)
        res = 1
    end if
    if (abs(PL - EXP_PL(3)) > tol) then
        print *, 'Test Case 3 Failed: PL = ', PL, ' Expected: ', EXP_PL(3)
        res = 1
    end if

    do i = 1, ITB
        do j = 1, JTB
            if (abs(PTBL(i,j) - EXP_PTBL(i,j,3)) / EXP_PTBL(i,j,3) > rel_tol) then
                print '(A,I0,A,I0,A,ES24.16,A,ES24.16)', 'Test Case 3 Failed: PTBL(', i, ',', j, ') = ', PTBL(i,j), &
                         ' Expected: ', EXP_PTBL(i,j,3)
                res = 1
            end if
            if (abs(TTBL(j,i) - EXP_TTBL(j,i,3)) / EXP_TTBL(j,i,3) > rel_tol) then
                print '(A,I0,A,I0,A,ES24.16,A,ES24.16)', 'Test Case 3 Failed: TTBL(', j, ',', i, ') = ', TTBL(j,i), &
                         ' Expected: ', EXP_TTBL(j,i,3)
                res = 1
            end if
        end do
        if (abs(STHE(i) - EXP_STHE(i,3)) / EXP_STHE(i,3) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 3 Failed: STHE(', i, ') = ', STHE(i), ' Expected: ', EXP_STHE(i,3)
            res = 1
        end if
        if (abs(THE0(i) - EXP_THE0(i,3)) / EXP_THE0(i,3) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 3 Failed: THE0(', i, ') = ', THE0(i), ' Expected: ', EXP_THE0(i,3)
            res = 1
        end if
    end do

    do j = 1, JTB
        if (abs(QS0(j) - EXP_QS0(j,3)) / EXP_QS0(j,3) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 3 Failed: QS0(', j, ') = ', QS0(j), ' Expected: ', EXP_QS0(j,3)
            res = 1
        end if
        if (abs(SQS(j) - EXP_SQS(j,3)) / EXP_SQS(j,3) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Case 3 Failed: SQS(', j, ') = ', SQS(j), ' Expected: ', EXP_SQS(j,3)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 30

    print *, 'SUCCESS!'

contains

    subroutine load_reference_data(case_num, ptbl_out, ttbl_out, qs0_out, sqs_out, sthe_out, the0_out)
        integer, intent(in) :: case_num
        real, intent(out) :: ptbl_out(ITB,JTB), ttbl_out(JTB,ITB)
        real, intent(out) :: qs0_out(JTB), sqs_out(JTB), sthe_out(ITB), the0_out(ITB)
        real :: temp_2d(ITB*JTB)
        character(len=100) :: filename, header_line
        integer :: unit_num, j
        
        write(filename, '(a,i1,a)') data_file_prefix, case_num, data_file_suffix
        open(newunit=unit_num, file=filename, status='old', action='read')
        
        ! Skip 3 header lines
        read(unit_num, '(a)') header_line
        read(unit_num, '(a)') header_line
        read(unit_num, '(a)') header_line
        
        do j = 1, ITB*JTB
            read(unit_num, *) temp_2d(j)
        end do
        ptbl_out = reshape(temp_2d, [ITB, JTB])
        
        do j = 1, JTB*ITB
            read(unit_num, *) temp_2d(j)
        end do
        ttbl_out = reshape(temp_2d, [JTB, ITB])
        
        do j = 1, JTB
            read(unit_num, *) qs0_out(j)
        end do
        
        do j = 1, JTB
            read(unit_num, *) sqs_out(j)
        end do
        
        do j = 1, ITB
            read(unit_num, *) sthe_out(j)
        end do
        
        do j = 1, ITB
            read(unit_num, *) the0_out(j)
        end do
        
        close(unit_num)
    end subroutine load_reference_data

end program test_table