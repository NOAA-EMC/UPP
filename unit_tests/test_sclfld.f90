! This is a test program for UPP.
!
! This program tests the SCLFLD() subroutine.
!
! Alyson Stahl, 2/2026
program test_sclfld
    use ctlblk_mod, only: jsta, jend, spval, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 2
    integer :: i, j, res
    integer :: IMO, JMO
    real :: SCALE
    real, dimension(1:npts,1:npts) :: FLD, EXP_FLD

    interface
        subroutine SCLFLD(FLD, SCALE, IMO,JMO)
            use ctlblk_mod, only: jsta, jend, ista, iend
            integer,intent(in) :: IMO, JMO
            REAL,intent(in) ::  SCALE
            REAL,dimension(ista:iend,jsta:jend), intent(inout) :: FLD
        end subroutine SCLFLD
    end interface

    ! Grid parameters
    jsta = 1
    jend = npts
    ista = 1
    iend = npts
    spval = 9.9e10

    SCALE = 100.0

    ! Set for sake of completeness. These values aren't actually used by SCLFLD.
    IMO = npts
    JMO = npts

    ! Test Case: Standard cases where FLD is scaled.
    FLD(1,1) = 1.0
    FLD(1,2) = 2.0
    FLD(2,1) = 3.0

    ! Test Case: FLD is filled with spval, should remain unchanged.
    FLD(2,2) = spval

    EXP_FLD(1,1) = 100.0
    EXP_FLD(1,2) = 200.0
    EXP_FLD(2,1) = 300.0
    EXP_FLD(2,2) = spval

    call SCLFLD(FLD, SCALE, IMO, JMO)

    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(FLD(i,j) - EXP_FLD(i,j)) > tol) then
                print *, 'FLD Failed for test', i, ': ', &
                            'Expected ', EXP_FLD(i,j), &
                            ' but got ', FLD(i,j)
                res = 1
            end if
        end do
    end do
    
    if (res .ne. 0) stop 10
    
    print *, 'SUCCESS!'
end program test_sclfld