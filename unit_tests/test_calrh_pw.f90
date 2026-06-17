! This is a test program for UPP.
!
! This program tests the CALRH_PW() subroutine in the upp_physics module.
!
! Alyson Stahl, 5/2026
program test_calrh_pw
    use vrbls3d, only: q, pmid, t
    use params_mod, only: g
    use ctlblk_mod, only: lm, ista, iend, jsta, jend, spval
    use upp_physics, only: CALRH_PW
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: npts = 5, nlevs = 30
    real :: RHPW(1, npts), EXP_RHPW(1, npts)
    integer :: i, res

    ! Grid Dimensions
    ista = 1
    iend = 1
    jsta = 1
    jend = npts
    lm = nlevs
    spval = 9.9e10

    allocate(q(ista:iend, jsta:jend, 1:lm))
    allocate(pmid(ista:iend, jsta:jend, 1:lm))
    allocate(t(ista:iend, jsta:jend, 1:lm))

    ! Test Case 1: Default case where none of the input values are spval. 
    do i = 1, lm
        pmid(1, :, i) = 100000.0 - real(i-1) / real(lm-1) * 95000.0
        t(1, :, i)    = 288.0 - 6.5e-3 * real(i-1) / real(lm-1) * 11000.0
        q(1, :, i)    = 0.005
    end do

    EXP_RHPW = 9.110954285E+01

    ! Test Case 2: Result of max(sh, Qs) evaluates to sh in the calculation of pw_sat.
    t(1, 2, :) = 220.0
    EXP_RHPW(1, 2) = 100.0

    ! Test Case 3: t == spval & q != spval
    t(1, 3, :) = spval
    EXP_RHPW(1, 3) = spval

    ! Test Case 4: t != spval & q == spval
    q(1, 4, :) = spval
    EXP_RHPW(1, 4) = spval

    ! Test Case 5: t == spval & q == spval
    t(1, 5, :) = spval
    q(1, 5, :) = spval
    EXP_RHPW(1, 5) = spval

    res = 0
    call CALRH_PW(RHPW)

    do i = 1, npts
        if (abs(RHPW(1, i) - EXP_RHPW(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected RHPW = ", EXP_RHPW(1, i), &
                     " but got RHPW = ", RHPW(1, i)
            res = 1
        end if
    end do

    deallocate(q, pmid, t)

    if (res .ne. 0) stop 10

    print *, "SUCCESS!"
end program test_calrh_pw