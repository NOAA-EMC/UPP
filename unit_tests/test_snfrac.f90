! This is a test program for UPP.
!
! This program tests the SNFRAC() subroutine.
!
! Alyson Stahl, 2/2026
program test_snfrac
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: ntests = 10
    real, parameter :: SALP = 4.0 ! From SNFRAC.f
    integer :: i, res
    real :: RSNOW ! For calculating expected values
    integer :: IVEGx(ntests)
    real, dimension(ntests) :: SNEQV, SNCOVR, EXP_SNCOVR

    interface
        subroutine SNFRAC(SNEQV, IVEGx, SNCOVR)
            integer,intent(in) :: IVEGx
            real,intent(in) ::  SNEQV
            real,intent(out) ::  SNCOVR
        end subroutine SNFRAC
    end interface

    ! Test Case: Standard case, where 0 < SNEQV < 1, for different vegetation types. 
    ! Covers all distinct values of SNUP (threshold depth)
    ! SNUP(1) = 0.080
    IVEGx(1) = 1
    SNEQV(1) = 0.06
    RSNOW = SNEQV(1) / 0.08
    EXP_SNCOVR(1) = 1. - (EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))

    ! SNUP(6) = 0.020
    IVEGx(2) = 6
    SNEQV(2) = 0.018
    RSNOW = SNEQV(2) / 0.02
    EXP_SNCOVR(2) = 1. - (EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))

    ! SNUP(8) = 0.060
    IVEGx(3) = 8
    SNEQV(3) = 0.048
    RSNOW = SNEQV(3) / 0.06
    EXP_SNCOVR(3) = 1. - (EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))

    ! SNUP(9) = 0.040
    IVEGx(4) = 9
    SNEQV(4) = 0.025
    RSNOW = SNEQV(4) / 0.04
    EXP_SNCOVR(4) = 1. - (EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))

    ! SNUP(11) = 0.01
    IVEGx(5) = 11
    SNEQV(5) = 0.005
    RSNOW = SNEQV(5) / 0.01
    EXP_SNCOVR(5) = 1. - (EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))

    ! SNUP(15) = 0.013
    IVEGx(6) = 15
    SNEQV(6) = 0.0052
    RSNOW = SNEQV(6) / 0.013
    EXP_SNCOVR(6) = 1. - (EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))

    ! Test Case: Invalid inputs for IVEGx, defaulting to IVEGx = 1
    IVEGx(7) = 0
    SNEQV(7) = 0.04
    RSNOW = SNEQV(7) / 0.08
    EXP_SNCOVR(7) = 1. - (EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))

    IVEGx(8) = 21
    SNEQV(8) = 0.05
    RSNOW = SNEQV(8) / 0.08
    EXP_SNCOVR(8) = 1. - (EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))

    ! Test Case: SNEQV exceeds threshold depth, clipping SNCOVR to 1.0
    IVEGx(9) = 1
    SNEQV(9) = 0.1
    EXP_SNCOVR(9) = 1.0

    ! Test Case: Invalid input for SNEQV, clipping SNCOVR to 0.0
    IVEGx(10) = 1
    SNEQV(10) = -0.1
    EXP_SNCOVR(10) = 0.0

    res = 0
    do i = 1, ntests
        call SNFRAC(SNEQV(i), IVEGx(i), SNCOVR(i))
        if (abs(SNCOVR(i) - EXP_SNCOVR(i)) > tol) then
            print *, 'SNCOVR Failed for test', i, ': ', &
                        'Expected ', EXP_SNCOVR(i), &
                        ' but got ', SNCOVR(i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10
    
    print *, 'SUCCESS!'
end program test_snfrac