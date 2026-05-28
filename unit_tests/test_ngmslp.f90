! This is a test program for UPP.
!
! This program tests the NGMSLP() subroutine.
!
! Alyson Stahl, 5/2026
program test_ngmslp
    use vrbls3d,    only: zint, pint, t, q, zmid
    use vrbls2d,    only: slp, fis, z1000
    use masks,      only: lmh
    use ctlblk_mod, only: jsta, jend, im, jm, spval, ista, iend, lm
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 4, nlevs = 30
    integer :: i, k, res
    real :: EXP_SLP(1, npts), EXP_Z1000(1, npts)

    interface
        subroutine NGMSLP()
        end subroutine NGMSLP
    end interface

    ! Grid Dimensions
    ista = 1
    iend = 1
    jsta = 1
    jend = npts
    lm = nlevs
    spval = 9.9e10

    allocate(zint(ista:iend, jsta:jend, 1:lm+1))
    allocate(pint(ista:iend, jsta:jend, 1:lm+1))
    allocate(t(ista:iend, jsta:jend, 1:lm))
    allocate(q(ista:iend, jsta:jend, 1:lm))
    allocate(zmid(ista:iend, jsta:jend, 1:lm))
    allocate(fis(ista:iend, jsta:jend))
    allocate(lmh(ista:iend, jsta:jend))

    ! Output arrays
    allocate(slp(ista:iend, jsta:jend))
    allocate(z1000(ista:iend, jsta:jend))

    ! Test Case 1: Default test case where none of the input values are spval.
    lmh = real(lm)
    fis(1, :) = 980.0
    
    do k = 1, lm + 1
        pint(1, :, k) = 100000.0 - real(k - 1) / real(lm) * 95000.0
        zint(1, :, k) = real(k - 1) / real(lm) * 25000.0
    end do
    
    do k = 1, lm
        zmid(1, :, k) = 0.5 * (zint(1, :, k) + zint(1, :, k+1))
        t(1, :, k) = 280.0 - real(k - 1) / real(lm - 1) * 60.0
        q(1, :, k) = 0.01 * exp(-real(k - 1) / real(lm - 1) * 3.0)
    end do

    EXP_SLP = 1.444795625E+05
    EXP_Z1000 = 1.229971484E+04

    ! Test Case 2: ((TAUSL>TAUCR).AND.(TAUSFC<=TAUCR))
    zint(1, 2, lm+1) = 8000.0
    fis(1, 2) = 8000.0 * 9.81
    zmid(1, 2, lm) = 8200.0
    t(1, 2, lm) = 260.0
    q(1, 2, lm) = 0.01

    EXP_SLP(1,2) =  1.342745410E+04
    EXP_Z1000(1,2) = -8.556410156E+04

    ! Test Case 3: ((TAUSL>TAUCR).AND.(TAUSFC>TAUCR))
    zint(1, 3, lm+1) = 100.0
    fis(1, 3) = 100.0 * 9.81
    zmid(1, 3, lm) = 500.0
    t(1, 3, lm) = 310.0
    q(1, 3, lm) = 0.02

    EXP_SLP(1,3) =  5.056930664E+03
    EXP_Z1000(1,3) = -1.667697344E+05

    ! Test Case 4: spval case
    pint(1, 4, lm+1) = spval

    EXP_SLP(1,4) =  spval
    EXP_Z1000(1,4) = spval

    res = 0
    call NGMSLP()

    do i = 1, npts
        if (abs(slp(1, i) - EXP_SLP(1, i)) > tol) then
            print *, "Error: SLP(1, ", i, ") = ", slp(1, i), " != EXP_SLP(1, ", i, ") = ", EXP_SLP(1, i)
            res = 1
        end if
        if (abs(z1000(1, i) - EXP_Z1000(1, i)) > tol) then
            print *, "Error: Z1000(1, ", i, ") = ", z1000(1, i), " != EXP_Z1000(1, ", i, ") = ", EXP_Z1000(1, i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10

    deallocate(zint, pint, t, q, zmid, fis, slp, z1000)

    print *, "SUCCESS!"
end program test_ngmslp