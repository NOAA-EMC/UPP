! This is a test program for UPP.
!
! This program tests the LFMFLD() subroutine.
!
! Alyson Stahl, 5/2026
program test_lfmfld
    use vrbls3d, only: pint, alpint, zint, t, q, cwm
    use masks, only: lmh
    use params_mod, only: d00, d50, pq0, a2, a3, a4, h1, d01, gi
    use ctlblk_mod, only: jsta, jend, modelname, spval, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 5, nlevs = 30
    integer :: i, j, res
    integer :: lm
    real :: RH3310(1, npts), RH6610(1, npts), RH3366(1, npts), PW3310(1, npts)
    real :: EXP_RH3310(1, npts), EXP_RH6610(1, npts), EXP_RH3366(1, npts), EXP_PW3310(1, npts)
    real :: EXP_RH3310_GFS(1, npts), EXP_RH6610_GFS(1, npts), EXP_RH3366_GFS(1, npts), EXP_PW3310_GFS(1, npts)
    interface
        subroutine LFMFLD(RH3310,RH6610,RH3366,PW3310)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, dimension(ista:iend,jsta:jend), intent(inout) :: RH3310, RH6610, RH3366
            real, dimension(ista:iend,jsta:jend), intent(inout) :: PW3310
        end subroutine LFMFLD
    end interface

    ! Grid Dimensions
    ista = 1
    iend = 1
    jsta = 1
    jend = npts
    lm = nlevs
    spval = 9.9e10
    modelname = ""

    allocate(pint(ista:iend, jsta:jend, 1:lm+1))
    allocate(alpint(ista:iend, jsta:jend, 1:lm+1))
    allocate(zint(ista:iend, jsta:jend, 1:lm+1))
    allocate(t(ista:iend, jsta:jend, 1:lm))
    allocate(q(ista:iend, jsta:jend, 1:lm))
    allocate(cwm(ista:iend, jsta:jend, 1:lm))
    allocate(lmh(ista:iend, jsta:jend))

    ! Test Case 1: Default case where none of the RH values are spval.
    lmh = real(lm)
    do i = 1, lm + 1
        pint(1, :, i)   = 5000.0 + real(i-1) / real(lm) * 95000.0
        alpint(1, :, i) = log(pint(1, :, i))
        zint(1, :, i)   = 20000.0 - real(i-1) / real(lm) * 20000.0
    end do
    do i = 1, lm
        t(1, :, i)   = 288.0 - 6.5e-3 * 0.5 * (zint(1, :, i) + zint(1, :, i+1))
        q(1, :, i)   = 0.01 * exp(-0.5 * (zint(1, :, i) + zint(1, :, i+1)) / 2000.0)
        cwm(1, :, i) = 1.0e-4
    end do

    EXP_RH3310 = 7.553577423E-01
    EXP_RH6610 = 7.212698460E-01
    EXP_RH3366 = 7.943152785E-01
    EXP_PW3310 = 7.344522476E+00

    EXP_RH3310_GFS = 8.572704792E-01
    EXP_RH6610_GFS = 7.627038360E-01
    EXP_RH3366_GFS = 9.653466344E-01
    EXP_PW3310_GFS = 7.344522476E+00

    ! Test Case 2: RH6610 != spval, RH3310 != spval, RH3366 == spval
    lmh(1, 2) = 2.0
    EXP_RH3310(1, 2) = 1.0
    EXP_RH6610(1, 2) = 1.0
    EXP_RH3366(1, 2) = spval
    EXP_PW3310(1, 2) = 3.245311603E-02

    EXP_RH3310_GFS(1, 2) = 1.0
    EXP_RH6610_GFS(1, 2) = 1.0
    EXP_RH3366_GFS(1, 2) = spval
    EXP_PW3310_GFS(1, 2) = 3.245311603E-02

    ! Test Case 3: RH6610 == spval, RH3310 != spval, RH3366 != spval
    lmh(1, 3)       = 3.0
    pint(1, 3, 2)   = 40000.0
    pint(1, 3, 3)   = 80000.0
    pint(1, 3, 4)   = 100000.0
    alpint(1, 3, 2) = log(40000.0)
    alpint(1, 3, 3) = log(80000.0)
    alpint(1, 3, 4) = log(100000.0)

    EXP_RH3310(1, 3) = 1.0
    EXP_RH6610(1, 3) = spval
    EXP_RH3366(1, 3) = 1.0
    EXP_PW3310(1, 3) = 4.107992351E-01

    EXP_RH3310_GFS(1, 3) = 1.0
    EXP_RH6610_GFS(1, 3) = spval
    EXP_RH3366_GFS(1, 3) = 1.0
    EXP_PW3310_GFS(1, 3) = 4.107992351E-01

    ! Test Case 4: RH6610 == spval, RH3310 == spval, RH3366 == spval
    lmh(1, 4) = 1.0

    EXP_RH3310(1, 4) = spval
    EXP_RH6610(1, 4) = spval
    EXP_RH3366(1, 4) = spval
    EXP_PW3310(1, 4) = 0.0

    EXP_RH3310_GFS(1, 4) = spval
    EXP_RH6610_GFS(1, 4) = spval
    EXP_RH3366_GFS(1, 4) = spval
    EXP_PW3310_GFS(1, 4) = 0.0

    ! Test Case 5: Executes the ``` if (RH < D01) ``` branch.
    q(1, 5, :) = 0.0

    EXP_RH3310(1, 5) = 9.999999776E-03
    EXP_RH6610(1, 5) = 9.999998845E-03
    EXP_RH3366(1, 5) = 1.000000071E-02
    EXP_PW3310(1, 5) = 4.841997027E-01

    EXP_RH3310_GFS(1, 5) = 9.999999776E-03
    EXP_RH6610_GFS(1, 5) = 9.999998845E-03
    EXP_RH3366_GFS(1, 5) = 1.000000071E-02
    EXP_PW3310_GFS(1, 5) = 4.841997027E-01

    print *, "Testing LFMFLD() with modelname != 'GFS'..."

    res = 0
    call LFMFLD(RH3310, RH6610, RH3366, PW3310)

    do i = 1, npts
        if (abs(RH3310(1, i) - EXP_RH3310(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected RH3310 = ", EXP_RH3310(1, i), &
                     " but got RH3310 = ", RH3310(1, i)
            res = 1
        end if
        if (abs(RH6610(1, i) - EXP_RH6610(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected RH6610 = ", EXP_RH6610(1, i), &
                     " but got RH6610 = ", RH6610(1, i)
            res = 1
        end if
        if (abs(RH3366(1, i) - EXP_RH3366(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected RH3366 = ", EXP_RH3366(1, i), &
                     " but got RH3366 = ", RH3366(1, i)
            res = 1
        end if
        if (abs(PW3310(1, i) - EXP_PW3310(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected PW3310 = ", EXP_PW3310(1, i), &
                     " but got PW3310 = ", PW3310(1, i)
            res = 1
        end if
    end do

    if (res .ne. 0) then
        deallocate(pint, alpint, zint, t, q, cwm, lmh)
        stop 10
    end if

    print *, "Testing LFMFLD() with modelname == 'GFS'..."
    modelname = "GFS"
    res = 0
    call LFMFLD(RH3310, RH6610, RH3366, PW3310)

    do i = 1, npts
        if (abs(RH3310(1, i) - EXP_RH3310_GFS(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected RH3310 = ", EXP_RH3310_GFS(1, i), &
                     " but got RH3310 = ", RH3310(1, i)
            res = 1
        end if
        if (abs(RH6610(1, i) - EXP_RH6610_GFS(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected RH6610 = ", EXP_RH6610_GFS(1, i), &
                     " but got RH6610 = ", RH6610(1, i)
            res = 1
        end if
        if (abs(RH3366(1, i) - EXP_RH3366_GFS(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected RH3366 = ", EXP_RH3366_GFS(1, i), &
                     " but got RH3366 = ", RH3366(1, i)
            res = 1
        end if
        if (abs(PW3310(1, i) - EXP_PW3310_GFS(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected PW3310 = ", EXP_PW3310_GFS(1, i), &
                     " but got PW3310 = ", PW3310(1, i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 20

    deallocate(pint, alpint, zint, t, q, cwm, lmh)
    print *, "SUCCESS!"
end program test_lfmfld