! This is a test program for UPP.
!
! This program tests the TABLEQ() subroutine.
!
! Alyson Stahl, 6/2026
program test_tableq
    use tableq_upp_mod, only: TABLEQ
    implicit none

    real, parameter :: tol = 1.0e-6, rel_tol = 1.0e-6
    integer, parameter :: ITB=152, JTB=440
    integer :: i, j, res
    !
    character(len=*), parameter :: data_file_name = 'data/ref_tableq.txt'
    !
    real :: PL, THL
    real :: TTBLQ(JTB,ITB)
    real :: STHE(ITB), THE0(ITB)
    real :: RDP, RDTHE
    !
    real :: EXP_TTBLQ(JTB,ITB), EXP_STHE(ITB), EXP_THE0(ITB)
    real :: EXP_RDP, EXP_RDTHE

    ! RDP = 1./DP -> DP = (PH - PL) / REAL(KPM - 1) = (PH - PL) / REAL(ITB - 1)
    ! RDTHE = 1./DTHE -> DTHE = 1 / REAL(KTHM - 1) = 1 / (JTB - 1)

    ! Same values used in INITPOST_GFS_NEMS_MPIIO()
    PL = 70000.0
    THL = 210.0

    EXP_RDP = 151.0 / 35000.0
    EXP_RDTHE = 1.0 / (1.0 / 439.0) ! Prevents precision related errors

    ! Load expected data
    call load_tableq_reference_data(data_file_name, EXP_TTBLQ, EXP_STHE, EXP_THE0)

    call TABLEQ(TTBLQ, RDP, RDTHE, PL, THL, STHE, THE0)

    res = 0

    if (abs(RDP - EXP_RDP) > tol) then
        print *, 'Test Failed: RDP = ', RDP, ' Expected: ', EXP_RDP
        res = 1
    end if
    if (abs(RDTHE - EXP_RDTHE) > tol) then
        print *, 'Test Failed: RDTHE = ', RDTHE, ' Expected: ', EXP_RDTHE
        res = 1
    end if
    do i = 1, JTB
        do j = 1, ITB
            if (abs(TTBLQ(i,j) - EXP_TTBLQ(i,j)) / EXP_TTBLQ(i,j) > rel_tol) then
                print '(A,I0,A,I0,A,ES24.16,A,ES24.16)', 'Test Failed: TTBLQ(', i, ',', j, ') = ', TTBLQ(i,j), &
                         ' Expected: ', EXP_TTBLQ(i,j)
                res = 1
            end if
        end do
    end do
    do i = 1, ITB
        if (abs(STHE(i) - EXP_STHE(i)) / EXP_STHE(i) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Failed: STHE(', i, ') = ', STHE(i), ' Expected: ', EXP_STHE(i)
            res = 1
        end if
        if (abs(THE0(i) - EXP_THE0(i)) / EXP_THE0(i) > rel_tol) then
            print '(A,I0,A,ES24.16,A,ES24.16)', 'Test Failed: THE0(', i, ') = ', THE0(i), ' Expected: ', EXP_THE0(i)
            res = 1
        end if
    end do
    if (res .ne. 0) stop 10

    print *, 'SUCCESS!'

contains

    subroutine load_tableq_reference_data(filename, ttblq_out, sthe_out, the0_out)
        character(len=*), intent(in) :: filename
        real, intent(out) :: ttblq_out(JTB,ITB), sthe_out(ITB), the0_out(ITB)
        real :: temp_2d(JTB*ITB)
        character(len=100) :: header_line
        integer :: unit_num, i, j
        
        open(newunit=unit_num, file=filename, status='old', action='read')
        
        ! Skip 3 header lines
        read(unit_num, '(a)') header_line
        read(unit_num, '(a)') header_line
        read(unit_num, '(a)') header_line
        
        do i = 1, JTB*ITB
            read(unit_num, *) temp_2d(i)
        end do
        ttblq_out = reshape(temp_2d, [JTB, ITB])
        
        do i = 1, ITB
            read(unit_num, *) sthe_out(i)
        end do
        
        do i = 1, ITB
            read(unit_num, *) the0_out(i)
        end do
        
        close(unit_num)
    end subroutine load_tableq_reference_data

    subroutine write_tableq_reference_data(filename, pl_in, thl_in, ttblq_in, sthe_in, the0_in)
        character(len=*), intent(in) :: filename
        real, intent(in) :: pl_in, thl_in
        real, intent(in) :: ttblq_in(JTB,ITB), sthe_in(ITB), the0_in(ITB)
        integer :: unit_num, i, j
        
        open(newunit=unit_num, file=filename, status='replace', action='write')
        
        write(unit_num, '(a)') '# Reference data for TABLEQ() subroutine test'
        write(unit_num, '(a,f10.1,a,f10.1)') '# Test case with PL=', pl_in, ', THL=', thl_in
        write(unit_num, '(a,i0,a,i0,a,i0,a)') &
            '# Format: TTBLQ(', JTB*ITB, ' values), STHE(', ITB, '), THE0(', ITB, ')'
        
        do i = 1, ITB
            do j = 1, JTB
                write(unit_num, '(es24.16)') ttblq_in(j,i)
            end do
        end do
        
        do i = 1, ITB
            write(unit_num, '(es24.16)') sthe_in(i)
        end do
        
        do i = 1, ITB
            write(unit_num, '(es24.16)') the0_in(i)
        end do
        
        close(unit_num)
    end subroutine write_tableq_reference_data
    
end program test_tableq