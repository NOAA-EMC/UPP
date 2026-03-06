! This is a test program for UPP.
!
! This program tests the SNFRAC_GFS() subroutine.
!
! Alyson Stahl, 2/2026
program test_snfrac_gfs
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: ntests = 5
    real, parameter :: SALP = 2.6 ! From SNFRAC_GFS.f
    integer :: i, res
    real :: RSNOW ! For calculating expected values
    integer :: IVEG(ntests)
    real, dimension(ntests) :: SNEQV, SNCOVR, EXP_SNCOVR

    interface
        subroutine SNFRAC_GFS(SNEQV, IVEG, SNCOVR)
            integer,intent(in) :: IVEG
            real,intent(in) ::  SNEQV
            real,intent(out) ::  SNCOVR
        end subroutine SNFRAC_GFS
    end interface

    ! Test Case: Standard case, where 0 < SNEQV < 1, for different vegetation types. 
    ! Covers all distinct values of SNUP (threshold depth)
    ! SNUP(1) = 0.080
    IVEG(1) = 1     
    SNEQV(1) = 0.06
    RSNOW = SNEQV(1) / 0.08
    EXP_SNCOVR(1) = 1. - (EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))

    ! SNUP(7) = 0.040
    IVEG(2) = 7
    SNEQV(2) = 0.02
    RSNOW = SNEQV(2) / 0.04
    EXP_SNCOVR(2) = 1. - (EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))

    ! SNUP(13) = 0.025
    IVEG(3) = 13
    SNEQV(3) = 0.02
    RSNOW = SNEQV(3) / 0.025
    EXP_SNCOVR(3) = 1. - (EXP(-SALP*RSNOW) - RSNOW*EXP(-SALP))

    ! Test Case: SNEQV exceeds threshold depth, clipping SNCOVR to 1.0
    IVEG(4) = 1
    SNEQV(4) = 0.1
    EXP_SNCOVR(4) = 1.0

    ! Test Case: Invalid input for SNEQV, clipping SNCOVR to 0.0
    IVEG(5) = 1
    SNEQV(5) = -0.1
    EXP_SNCOVR(5) = 0.0

    res = 0
    do i = 1, ntests
        call SNFRAC_GFS(SNEQV(i), IVEG(i), SNCOVR(i))
        if (abs(SNCOVR(i) - EXP_SNCOVR(i)) > tol) then
            print *, 'SNCOVR Failed for test', i, ': ', &
                        'Expected ', EXP_SNCOVR(i), &
                        ' but got ', SNCOVR(i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10
    
    print *, 'SUCCESS!'
end program test_snfrac_gfs