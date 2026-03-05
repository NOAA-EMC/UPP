! This is a test program for UPP.
!
! This program tests the ICAOHEIGHT() subroutine.
!
! Alyson Stahl, 1/2026
program test_icaoheight
    use ctlblk_mod, only: jsta, jend, ista, iend, spval
    implicit none

    real, parameter :: tol = 1.0e-8
    ! From ICAOHEIGHT.f
    real, parameter :: Press_Bot = 101325., Press_Mid = 22632., Press_Top = 5474.87
    integer, parameter :: npts = 3
    integer :: i, j, res
    real, dimension(1:npts,1:npts) :: MAXWP, MAXWICAOZ, EXP_MAXWICAOZ

    interface
        subroutine ICAOHEIGHT(MAXWP, MAXWICAOZ)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, intent(in) :: MAXWP(ista:iend,jsta:jend)
            real, intent(out) :: MAXWICAOZ(ista:iend,jsta:jend)
        end subroutine ICAOHEIGHT
    end interface

    ! Grid parameters
    jsta = 1
    jend = npts
    ista = 1
    iend = npts
    spval = 9.9e10

    MAXWP(1,1) = 500.0              ! 0 < MAXWP < 1000
    MAXWP(1,2) = spval              ! MAXWP = spval
    MAXWP(1,3) = Press_Bot + 500.0  ! Press_Bot < MAXWP < spval
    MAXWP(2,1) = Press_Mid + 500.0  ! Press_Mid < MAXWP < Press_Bot
    MAXWP(2,2) = Press_Mid          ! MAXWP = Press_Mid
    MAXWP(2,3) = Press_Top + 500.0  ! Press_Top < MAXWP < Press_Mid
    MAXWP(3,1) = Press_Top          ! MAXWP = Press_Top
    MAXWP(3,2) = 3500.0             ! 1000.0 < MAXWP < Press_Top
    MAXWP(3,3) = 2750.0             ! 1000.0 < MAXWP < Press_Top

    EXP_MAXWICAOZ = reshape([3.1054496094E+04, 1.0861050781E+04,  2.0000000000E+04, &
                            spval, 1.1000000000E+04, 2.2855916016E+04, 0.0, &
                            1.9445695312E+04, 2.4410888672E+04], [npts,npts])

    call ICAOHEIGHT(MAXWP, MAXWICAOZ)

    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(MAXWICAOZ(i,j) - EXP_MAXWICAOZ(i,j)) > tol) then
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    print *, 'SUCCESS!'
end program test_icaoheight