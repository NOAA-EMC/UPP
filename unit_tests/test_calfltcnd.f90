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
    integer, parameter :: npts = 14
    integer :: i, res
    real :: CEILING(1, npts)
    real :: FLTCND(1, npts), EXP_FLTCND(1, npts)
    real :: ft_to_m, mi_to_m

    interface
        subroutine CALFLTCND(CEILING,FLTCND)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, dimension(ista:iend,jsta:jend), intent(in) :: CEILING
            real, dimension(ista:iend,jsta:jend), intent(inout) :: FLTCND
        end subroutine CALFLTCND
    end interface

    ! See converstion factors used in CALFLTCND() subroutine
    ft_to_m = 1/3.2808
    mi_to_m = 1609.0

    ! Grid parameters
    ista = 1
    iend = 1
    jsta = 1
    jend = npts
    spval = 9.9e10

    allocate(vis(1,1:npts))

    ! NOTE: This subroutine converts CEILING from m to ft and VIS from m to miles.

    ! Fight Condition 1 Cases: CEILING < 500.0 ft OR VIS < 1.0 mi
    ! Test Case 1: CEILING < 500.0 ft & VIS < 1.0 mi
    CEILING(1, 1) = 400.0 * ft_to_m
    vis(1, 1) = 0.9 * mi_to_m
    EXP_FLTCND(1, 1) = 1.0

    ! Test Case 2: CEILING < 500.0 ft & VIS > 1.0 mi
    CEILING(1, 2) = 400.0 * ft_to_m
    vis(1, 2) = 1.1 * mi_to_m
    EXP_FLTCND(1, 2) = 1.0

    ! Test Case 3: CEILING > 500.0 ft & VIS < 1.0 mi
    CEILING(1, 3) = 600.0 * ft_to_m
    vis(1, 3) = 0.9 * mi_to_m
    EXP_FLTCND(1, 3) = 1.0

    ! Flight Condition 2 Cases: 500 ft <= CEILING < 1000 ft 
    ! OR 1.0 mi <= VIS < 3.0 mi
    ! Test Case 4: CEILING = 500.0 ft & VIS = 1.0 mi
    CEILING(1, 4) = 500.0 * ft_to_m
    vis(1, 4) = 1.0 * mi_to_m
    EXP_FLTCND(1, 4) = 2.0

    ! Test Case 5: 500 ft < CEILING < 1000 ft & VIS > 3.0 mi
    CEILING(1, 5) = 700.0 * ft_to_m
    vis(1, 5) = 3.1 * mi_to_m
    EXP_FLTCND(1, 5) = 2.0

    ! Test Case 6: CEILING > 1000 ft & 1.0 mi < VIS < 3.0 mi
    CEILING(1, 6) = 1200.0 * ft_to_m
    vis(1, 6) = 2.0 * mi_to_m
    EXP_FLTCND(1, 6) = 2.0

    ! Flight Condition 3 Cases: 1000 ft <= CEILING < 3,000 ft
    ! OR 3.0 mi <= VIS < 5.0 mi 
    ! Test Case 7: CEILING = 1000 ft & VIS = 3.0 mi
    CEILING(1, 7) = 1000.0 * ft_to_m
    vis(1, 7) = 3.0 * mi_to_m
    EXP_FLTCND(1, 7) = 3.0

    ! Test Case 8: 1000 ft < CEILING < 3000 ft & VIS > 5.0 mi
    CEILING(1, 8) = 2000.0 * ft_to_m
    vis(1, 8) = 6.0 * mi_to_m
    EXP_FLTCND(1, 8) = 3.0

    ! Test Case 9: CEILING > 3000 ft & 3.0 mi < VIS < 5.0 mi
    CEILING(1, 9) = 3500.0 * ft_to_m
    vis(1, 9) = 3.0 * mi_to_m
    EXP_FLTCND(1, 9) = 3.0

    ! Test Case 10: CEILING = 3000 ft & VIS = 5.0 mi
    CEILING(1, 10) = 3000.0 * ft_to_m
    vis(1, 10) = 5.0 * mi_to_m
    EXP_FLTCND(1, 10) = 3.0

    ! Flight Condition 4 Cases: CEILING > 3,000 ft OR VIS > 5.0 mi 
    ! Test Case 11: CEILING > 3000 ft & VIS > 5.0 mi
    CEILING(1, 11) = 3500.0 * ft_to_m
    vis(1, 11) = 6.0 * mi_to_m
    EXP_FLTCND(1, 11) = 4.0

    ! Test Case 12: CEILING = spval & VIS = spval
    CEILING(1, 12) = spval
    vis(1, 12) = spval
    EXP_FLTCND(1, 12) = spval

    ! Test Case 13: CEILING = spval & VIS != spval
    CEILING(1, 13) = spval
    vis(1, 13) = 1.0 * mi_to_m
    EXP_FLTCND(1, 13) = spval

    ! Test Case 14: CEILING != spval & VIS = spval
    CEILING(1, 14) = 100.0 * ft_to_m
    vis(1, 14) = spval
    EXP_FLTCND(1, 14) = spval

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
