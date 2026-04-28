! This is a test program for UPP.
!
! This program tests the CALFLTCND() subroutine.
!
! Alyson Stahl, 4/2026
program test_calfltcnd
    use vrbls2d, only: vis
    use ctlblk_mod, only: jsta, jend, im, spval, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 16
    integer :: i, res
    real :: CEILING(1, npts)
    real :: FLTCND(1, npts), EXP_FLTCND(1, npts)
    real :: ft_to_m, mi_to_m
    
    ft_to_m = 0.3048 ! Conversion factor from feet to meters
    mi_to_m = 1609.34 ! Conversion factor from miles to meters

    interface
        subroutine CALFLTCND(CEILING,FLTCND)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, dimension(ista:iend,jsta:jend), intent(in) :: CEILING
            real, dimension(ista:iend,jsta:jend), intent(inout) :: FLTCND
        end subroutine CALFLTCND
    end interface

    ! Grid parameters
    ista = 1
    iend = 1
    jsta = 1
    jend = npts
    spval = 9.9e10

    allocate(vis(1,1:npts))

    ! NOTE: This subroutine converts CEILING from m to ft and VIS from m to miles.

    ! Fight Condition 1 Cases: CEILING < 500.0 ft (152.4 m) OR VIS < 1.0 mi (1609.34 m)
    ! Test Case 1: CEILING < 152.4 m & VIS < 1609.34 m
    CEILING(1, 1) = 100.0
    vis(1, 1) = 1000.0
    EXP_FLTCND(1, 1) = 1.0

    ! Test Case 2: CEILING < 152.4 m & VIS > 1609.34 m
    CEILING(1, 2) = 100.0
    vis(1, 2) = 2000.0
    EXP_FLTCND(1, 2) = 1.0

    ! Test Case 3: CEILING > 152.4 m & VIS < 1609.34 m
    CEILING(1, 3) = 200.0
    vis(1, 3) = 1000.0
    EXP_FLTCND(1, 3) = 1.0

    ! Flight Condition 2 Cases: 500 ft (152.4 m) <= CEILING < 1000 ft (304.8 m) 
    ! OR 1.0 mi (1609.34 m) <= VIS < 3.0 mi (4828.02 m)
    ! Test Case 4: CEILING = 152.4 m & VIS = 1609.34 m
    CEILING(1, 4) = 152.4
    vis(1, 4) = 1609.34
    EXP_FLTCND(1, 4) = 2.0

    ! Test Case 5: 152.4 m < CEILING < 304.8 m & VIS > 4828.02 m
    CEILING(1, 5) = 200.0
    vis(1, 5) = 5000.0
    EXP_FLTCND(1, 5) = 2.0

    ! Test Case 6: CEILING > 304.8 m & 1609.34 m < VIS < 4828.02 m
    CEILING(1, 6) = 400.0
    vis(1, 6) = 2000.0
    EXP_FLTCND(1, 6) = 2.0

    ! Flight Condition 3 Cases: 1000 ft (304.8 m) <= CEILING < 3,000 ft (914.4 m) 
    ! OR 3.0 mi (4828.02 m) <= VIS < 5.0 mi (8046.72 m)
    ! Test Case 7: CEILING = 304.8 m & VIS = 4828.02 m
    CEILING(1, 7) = 304.8
    vis(1, 7) = 4828.02
    EXP_FLTCND(1, 7) = 3.0

    ! Test Case 8: 304.8 m < CEILING < 914.4 m & VIS > 8046.72 m
    CEILING(1, 8) = 400.0
    vis(1, 8) = 10000.0
    EXP_FLTCND(1, 8) = 3.0

    ! Test Case 9: CEILING > 914.4 m & 4828.02 m < VIS < 8046.72 m
    CEILING(1, 9) = 1000.0
    vis(1, 9) = 5000.0
    EXP_FLTCND(1, 9) = 3.0

    ! Test Case 10: CEILING = 914.4 m & VIS = 8046.72 m
    CEILING(1, 10) = 914.4
    vis(1, 10) = 8046.72
    EXP_FLTCND(1, 10) = 3.0

    ! Flight Condition 4 Cases: CEILING > 3,000 ft (914.4 m) OR VIS > 5.0 mi (8046.72 m)
    ! Test Case 11: CEILING > 914.4 m & VIS > 8046.72 m
    CEILING(1, 11) = 1000.0
    vis(1, 11) = 10000.0
    EXP_FLTCND(1, 11) = 4.0

    ! Test Case 12: CEILING < 914.4 m & VIS > 8046.72 m
    CEILING(1, 12) = 900.0
    vis(1, 12) = 10000.0
    EXP_FLTCND(1, 12) = 4.0

    ! Test Case 13: CEILING > 914.4 m & VIS < 8046.72 m
    CEILING(1, 13) = 1000.0
    vis(1, 13) = 5000.0
    EXP_FLTCND(1, 13) = 4.0

    ! Test Case 14: CEILING = spval & VIS = spval
    CEILING(1, 14) = spval
    vis(1, 14) = spval
    EXP_FLTCND(1, 14) = spval

    ! Test Case 15: CEILING = spval & VIS != spval
    CEILING(1, 15) = spval
    vis(1, 15) = 1000.0
    EXP_FLTCND(1, 15) = spval

    ! Test Case 16: CEILING != spval & VIS = spval
    CEILING(1, 16) = 100.0
    vis(1, 16) = spval
    EXP_FLTCND(1, 16) = spval

    call CALFLTCND(CEILING, FLTCND)

    res = 0
    do i = 1, npts
        if (abs(FLTCND(1,i) - EXP_FLTCND(1,i)) > tol) then
            print *, "Test Case ", i, " failed: Expected FLTCND = ", EXP_FLTCND(1,i), &
                     " but got ", FLTCND(1,i)
            res = 1
        end if
    end do

    deallocate(vis)
    
    if (res .ne. 0) stop 10

    print *, "SUCCESS!"
end program test_calfltcnd
