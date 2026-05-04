! This is a test program for UPP.
!
! This program tests the CALGUSTCONV() subroutine.
!
! Alyson Stahl, 4/2026
program test_calgustconv
    use vrbls2d , only: u10, v10, ustar
    use ctlblk_mod, only: ista, iend, jsta, jend, ista_2l, iend_2u, jsta_2l, jend_2u, spval
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 5
    integer :: j, res
    real :: SPEED850(1, 1:npts), SPEED950(1, 1:npts)
    real :: GUSTCONV(1, 1:npts), EXP_GUSTCONV(1, 1:npts)

    interface
        subroutine CALGUSTCONV(SPEED850,SPEED950,GUSTCONV)
            use ctlblk_mod, only: ista_2l, iend_2u, jsta_2l, jend_2u
            real, intent(in) :: SPEED850(ista_2l:iend_2u,jsta_2l:jend_2u)
            real, intent(in) :: SPEED950(ista_2l:iend_2u,jsta_2l:jend_2u)
            real, intent(inout) :: GUSTCONV(ista_2l:iend_2u,jsta_2l:jend_2u)
        end subroutine CALGUSTCONV
    end interface

    ! Grid dimensions
    ista = 1
    iend = 1
    jsta = 1
    jend = npts
    ista_2l = 1
    iend_2u = 1
    jsta_2l = 1
    jend_2u = npts
    spval = 9.9e10

    allocate(u10(ista_2l:iend_2u,jsta_2l:jend_2u))
    allocate(v10(ista_2l:iend_2u,jsta_2l:jend_2u))
    allocate(ustar(ista_2l:iend_2u,jsta_2l:jend_2u))

    u10 = 12.0
    v10 = 16.0
    ustar = 0.7

    SPEED950 = 26.0
    SPEED850 = 34.0

    ! Test Case 1: Typical case
    EXP_GUSTCONV(1,1) = 24.5

    ! Test Case 2: SPEED950 > SPEED850 (Negative wind speed difference clipped to 0)
    SPEED950(1,2) = 35.0
    EXP_GUSTCONV(1,2) = 22.1

    ! Test Case 3: SPEED850 = SPEED950 (Wind speed difference is 0)
    SPEED950(1,3) = 34.0
    EXP_GUSTCONV(1,3) = 22.1

    ! Test Case 4:  u10 and v10 have spval
    u10(1,4) = spval
    v10(1,4) = spval
    EXP_GUSTCONV(1,4) = spval

    ! Test Case 5: SPEED850 and SPEED950 have spval
    SPEED850(1,5) = spval
    SPEED950(1,5) = spval
    EXP_GUSTCONV(1,5) = spval

    call CALGUSTCONV(SPEED850,SPEED950,GUSTCONV)

    deallocate(u10, v10, ustar)

    res = 0
    do j = 1, npts
        if (abs(GUSTCONV(1,j) - EXP_GUSTCONV(1,j)) > tol) then
            print *, "Test Case ", j, " Failed: GUSTCONV = ", GUSTCONV(1,j), &
                " Expected = ", EXP_GUSTCONV(1,j)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10
    
    print *, 'SUCCESS!'
end program test_calgustconv