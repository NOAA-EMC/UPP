! This is a test program for UPP.
!
! This program tests the CALICING() subroutine.
!
! Alyson Stahl, 4/2026
program test_calicing
    use ctlblk_mod, only: jsta, jend, spval, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 11
    integer :: i, res
    real :: T1(1, npts), RH(1, npts), OMGA(1, npts)
    real :: ICING(1, npts), EXP_ICING(1, npts)

    interface
        subroutine CALICING(T1,RH,OMGA, ICING)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, dimension(ista:iend,jsta:jend), intent(in) :: T1,RH,OMGA
            real, dimension(ista:iend,jsta:jend), intent(inout) :: ICING 
        end subroutine CALICING
    end interface

    ! Grid parameters
    ista = 1
    iend = 1
    jsta = 1
    jend = npts
    spval = 9.9e10

    ! Test Case 1: OMGA < 0 & 251 < T1 < 273 & RH > 70 (expect ICING = 1)
    T1 = 260.0
    RH = 80.0
    OMGA = -0.1
    EXP_ICING(1,1) = 1.0

    ! Test Case 2: OMGA > 0 (expect ICING = 0)
    OMGA(1,2) = 0.1
    EXP_ICING(1,2) = 0.0

    ! Test Case 3: T1 < 251 (expect ICING = 0)
    T1(1,3) = 250.0
    EXP_ICING(1,3) = 0.0

    ! Test Case 4: T1 > 273 (expect ICING = 0)
    T1(1,4) = 274.0
    EXP_ICING(1,4) = 0.0

    ! Test Case 5: RH < 70 (expect ICING = 0)
    RH(1,5) = 60.0
    EXP_ICING(1,5) = 0.0

    ! Test Case 6: OMGA < 0 & T1 == 251 & RH > 70 (expect ICING = 1)
    T1(1,6) = 251.0
    EXP_ICING(1,6) = 1.0

    ! Test Case 7: OMGA < 0 & T1 == 273 & RH > 70 (expect ICING = 1)
    T1(1,7) = 273.0
    EXP_ICING(1,7) = 1.0

    ! Test Case 8: OMGA < 0 & 251 < T1 < 273 & RH == 70 (expect ICING = 1)
    RH(1,8) = 70.0
    EXP_ICING(1,8) = 1.0

    ! Test Case 9: OMGA > spval (expect ICING = spval)
    OMGA(1,9) = spval
    EXP_ICING(1,9) = spval

    ! Test Case 10: T1 > spval (expect ICING = spval)
    T1(1,10) = spval
    EXP_ICING(1,10) = spval

    ! Test Case 11: RH > spval (expect ICING = spval)
    RH(1,11) = spval
    EXP_ICING(1,11) = spval

    call CALICING(T1, RH, OMGA, ICING)

    res = 0
    do i = 1, npts
        if (abs(ICING(1,i) - EXP_ICING(1,i)) > tol) then
            print *, "Test Case ", i, " failed: ICING = ", ICING(1,i), &
                     " but expected ", EXP_ICING(1,i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10
    
    print *, "SUCCESS!"
end program test_calicing
