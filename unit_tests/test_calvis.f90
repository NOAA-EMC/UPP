! This is a test program for UPP.
!
! This program tests the CALVIS() subroutine.
!
! Alyson Stahl, 12/2025
program test_calvis
    use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u, spval, &
                          ista, iend, ista_2l, iend_2u
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 3
    integer :: i, j, res
    real, dimension(1:npts,1:npts) :: QV,QC,QR,QI,QS,TT,PP
    real, dimension(1:npts,1:npts) :: VIS, EXP_VIS

    ! Grid parameters
    jsta = 1
    jend = npts
    jsta_2l = jsta
    jend_2u = jend
    ista = 1
    iend = npts
    ista_2l = ista
    iend_2u = iend
    spval = 9.9e10

    ! TODO: Set baseline values for the input arrays. These values should be reasonable "real life" values. 
    ! The values should also guarantee that CONST1/BETAV < 24.135. All should be LESS than spval. 
    QV = 0.01       ! water vapor mixing ratio (kg/kg)
    QC = 1.0e-3     ! cloud water mixing ratio (kg/kg)
    QR = 0.0        ! rain water mixing ratio (kg/kg)
    QI = 0.0        ! cloud ice mixing ratio (kg/kg)
    QS = 0.0        ! snow mixing ratio (kg/kg)
    TT = 280.0      ! temperature (K)
    PP = 101325.0   ! pressure (Pa)

    QV(1,1) = spval
    QC(2,1) = spval
    QR(3,1) = spval
    QI(1,2) = spval
    QS(2,2) = spval
    TT(3,2) = spval
    PP(1,3) = spval

    ! TODO: Set some array values at (2,3) that will give VIS(2,3) = 1.E3 * 24.135

    ! The assignments go here.
    QV(2,3) = 0.005
    QC(2,3) = 0.0
    QR(2,3) = 0.0
    QI(2,3) = 0.0
    QS(2,3) = 0.0
    TT(2,3) = 280.0
    PP(2,3) = 101325.0

    call CALVIS(QV, QC, QR, QI, QS, TT, PP, VIS)

    res = 0
    do j = jsta, jend
        do i = ista, iend
            ! TODO: Print out VIS at each (i,j) in ESw.d form 10th decimal place
            write(*,'(ES24.10)') VIS(i,j)
        end do
    end do

    if (res .ne. 0) stop 10

    print *, 'SUCCESS!'
end program test_calvis