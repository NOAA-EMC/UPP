! This is a test program for UPP.
!
! This program tests the ICAOHEIGHT() subroutine.
!
! Alyson Stahl, 1/2025
program test_icaoheight
    use ctlblk_mod, only: jsta, jend, spval, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-8
    ! From ICAOHEIGHT.f
    real, parameter :: Press_Bot = 101325., Press_Mid = 22632., Press_Top = 5474.87
    integer, parameter :: npts = 2
    integer :: i, j, res
    real, dimension(1:npts,1:npts) :: MAXWP, MAXWICAOZ, EXP_MAXWICAOZ

    ! Grid parameters
    jsta = 1
    jend = npts
    ista = 1
    iend = npts
    spval = 9.9e10

    MAXWP(1,1) = 500.0       
    MAXWP(1,2) = Press_Bot + 500.0 
    MAXWP(2,1) = spval 
    MAXWP(2,2) = Press_Top + 500.0

    call ICAOHEIGHT(MAXWP, MAXWICAOZ)

    res = 0
    do j = jsta, jend
        do i = ista, iend
            print '(A,I0,A,I0,A,ES24.10)', 'MAXWICAOZ(', i, ',', j, ') = ', MAXWICAOZ(i,j)
        end do
    end do

    if (res .ne. 0) stop 10

    print *, 'SUCCESS!'
end program test_icaoheight