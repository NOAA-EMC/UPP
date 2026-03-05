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

    interface
        subroutine CALVIS(QV, QC, QR, QI, QS, TT, PP, VIS)
            use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u, &
                                ista, iend, ista_2l, iend_2u
            real, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in) :: &
                QV, QC, QR, QI, QS, TT, PP
            real, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(inout) :: VIS
        end subroutine CALVIS
    end interface

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

    ! Test case where CONST1/BETAV < 24.135
    QV = 0.01       ! water vapor mixing ratio (kg/kg)
    QC = 1.0e-3     ! cloud water mixing ratio (kg/kg)
    QR = 0.0        ! rain water mixing ratio (kg/kg)
    QI = 0.0        ! cloud ice mixing ratio (kg/kg)
    QS = 0.0        ! snow mixing ratio (kg/kg)
    TT = 280.0      ! temperature (K)
    PP = 101325.0   ! pressure (Pa)

    ! Set some array values to spval to test the handling of spval
    QV(1,1) = spval
    QC(2,1) = spval
    QR(3,1) = spval
    QI(1,2) = spval
    QS(2,2) = spval
    TT(3,2) = spval
    PP(1,3) = spval

    ! Test case where CONST1/BETAV > 24.135
    QV(2,3) = 0.005
    QC(2,3) = 0.0
    QR(2,3) = 0.0
    QI(2,3) = 0.0
    QS(2,3) = 0.0
    TT(2,3) = 280.0
    PP(2,3) = 101325.0

    EXP_VIS = reshape([spval, spval, spval, spval, spval, spval, spval, 24135.0, 22.361997604], [npts, npts])
    
    call CALVIS(QV, QC, QR, QI, QS, TT, PP, VIS)

    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(VIS(i,j) - EXP_VIS(i,j)) > tol) then
                print *, 'VIS Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_VIS(i,j), &
                         ' but got ', VIS(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    print *, 'SUCCESS!'
end program test_calvis