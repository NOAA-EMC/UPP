! This is a test program for UPP.
!
! This program tests the FRZLVL() subroutine.
!
! Alyson Stahl, 5/2026
program test_frzlvl
      use vrbls3d, only: pint, t, zmid, q, pmid
      use vrbls2d, only: fis, tshltr, pshltr, qshltr
      use masks, only: lmh
      use ctlblk_mod, only: jsta, jend, spval, lm, modelname, im, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: npts = 7, nlevs = 30
    integer :: i, j, k, res
    real :: ZFRZ(1, npts), RHFRZ(1, npts), PFRZL(1, npts)
    real :: EXP_ZFRZ(1, npts), EXP_RHFRZ(1, npts), EXP_PFRZL(1, npts)
    real :: EXP_ZFRZ_GFS(1, npts), EXP_RHFRZ_GFS(1, npts), EXP_PFRZL_GFS(1, npts)

    interface 
        subroutine FRZLVL(ZFRZ,RHFRZ,PFRZL)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, dimension(ista:iend,jsta:jend), intent(out) :: ZFRZ, RHFRZ, PFRZL
        end subroutine FRZLVL
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
    allocate(t(ista:iend, jsta:jend, 1:lm))
    allocate(zmid(ista:iend, jsta:jend, 1:lm))
    allocate(q(ista:iend, jsta:jend, 1:lm))
    allocate(pmid(ista:iend, jsta:jend, 1:lm))
    allocate(fis(ista:iend, jsta:jend))
    allocate(tshltr(ista:iend, jsta:jend))
    allocate(pshltr(ista:iend, jsta:jend))
    allocate(qshltr(ista:iend, jsta:jend))
    allocate(lmh(ista:iend, jsta:jend))

    ! Test Case 1: Default case where none of the input values are spval.
    lmh = real(lm)
    fis(1, :) = 0.0
    tshltr(1, :) = 288.0
    pshltr(1, :) = 100000.0
    qshltr(1, :) = 0.010
    
    ! Initialize 3D arrays: Level lm is at bottom, level 1 is at top
    do k = 1, lm + 1
        pint(1, :, k) = 5000.0 + real(lm + 1 - k) / real(lm) * 95000.0
    end do
    
    do k = 1, lm
        pmid(1, :, k) = 0.5 * (pint(1, :, k) + pint(1, :, k+1))
        ! Temperature decreases from 288 K at surface (k=lm) to ~218 K at top (k=1)
        ! Using 6.5 K/km standard lapse rate over ~11 km
        t(1, :, k) = 288.0 - 6.5e-3 * real(lm - k) / real(lm - 1) * 11000.0
        ! Height increases from ~0 m at k=lm to ~20000 m at k=1
        zmid(1, :, k) = real(lm - k) / real(lm - 1) * 20000.0
        ! Specific humidity decreases exponentially with height
        q(1, :, k) = 0.010 * exp(-real(lm - k) / real(lm - 1) * 3.0)
    end do

    EXP_ZFRZ(1, 1) =  4.153851562E+03
    EXP_RHFRZ(1, 1) =  3.624264598E-01
    EXP_PFRZL(1, 1) =  2.565234961E+04

    EXP_ZFRZ_GFS(1, 1) =  4.153851562E+03
    EXP_RHFRZ_GFS(1, 1) =  3.591695130E-01
    EXP_PFRZL_GFS(1, 1) =  2.565234961E+04

    ! Test Case 2: TSHLTR and PSHLTR are spval with other default input values.
    tshltr(1, 2) = spval
    pshltr(1, 2) = spval

    EXP_ZFRZ(1, 2) =  4.153851562E+03
    EXP_RHFRZ(1, 2) =  3.624264598E-01
    EXP_PFRZL(1, 2) =  2.565234961E+04

    EXP_ZFRZ_GFS(1, 2) =  4.153851562E+03
    EXP_RHFRZ_GFS(1, 2) =  3.591695130E-01
    EXP_PFRZL_GFS(1, 2) =  2.565234961E+04

    ! Test Case 3: QSHLTR is spval with other default input values.
    qshltr(1, 3) = spval

    EXP_ZFRZ(1, 3) =  4.153851562E+03
    EXP_RHFRZ(1, 3) =  3.624264598E-01
    EXP_PFRZL(1, 3) =  2.565234961E+04

    EXP_ZFRZ_GFS(1, 3) =  4.153851562E+03
    EXP_RHFRZ_GFS(1, 3) =  3.591695130E-01
    EXP_PFRZL_GFS(1, 3) =  2.565234961E+04

    ! Test Case 4: TSFC < TFRZ
    tshltr(1, 4) = 270.0

    EXP_ZFRZ(1, 4) = -4.826144409E+02
    EXP_RHFRZ(1, 4) =  1.000000000E+00
    EXP_PFRZL(1, 4) =  1.000000000E+05

    EXP_ZFRZ_GFS(1, 4) = -4.826144409E+02
    EXP_RHFRZ_GFS(1, 4) =  1.000000000E+00
    EXP_PFRZL_GFS(1, 4) =  1.000000000E+05

    ! Test Case 5: In the L do loop labeled 10, T(I,J,L) <= TFRZ at L == LLMH.
    tshltr(1, 5) = 275.0
    t(1, 5, lm) = 270.0
    
    EXP_ZFRZ(1, 5) =  1.259997606E+00
    EXP_RHFRZ(1, 5) =  1.458200514E-01
    EXP_PFRZL(1, 5) =  5.535744629E+03

    EXP_ZFRZ_GFS(1, 5) =  1.259997606E+00
    EXP_RHFRZ_GFS(1, 5) =  1.397437900E-01
    EXP_PFRZL_GFS(1, 5) =  5.535744629E+03

    ! Test Case 6: QSHLTR and PSHLTR are spval with the T(I,J,LLMH) <= TFRZ case.
    qshltr(1, 6) = spval
    pshltr(1, 6) = spval
    t(1, 6, lm) = 272.0
    zmid(1, 6, lm) = 1000.0

    EXP_ZFRZ(1, 6) = 8.230778198E+02
    EXP_RHFRZ(1, 6) =  1.651607752E-01
    EXP_PFRZL(1, 6) =  6.269973633E+03

    EXP_ZFRZ_GFS(1, 6) = 8.230778198E+02
    EXP_RHFRZ_GFS(1, 6) =  1.590846628E-01
    EXP_PFRZL_GFS(1, 6) =  6.269973633E+03

    ! Test Case 7: T(I,J,L) > TFRZ at all levels.
    tshltr(1, 7) = 300.0
    do k = 1, lm
        t(1, 7, k) = 300.0 - 2.0 * real(lm - k) / real(lm - 1) * 11000.0 / 1000.0
    end do

    EXP_ZFRZ(1, 7) =  0.0
    EXP_RHFRZ(1, 7) =  0.0
    EXP_PFRZL(1, 7) =  5000.0

    EXP_ZFRZ_GFS(1, 7) =  0.0
    EXP_RHFRZ_GFS(1, 7) =  0.0
    EXP_PFRZL_GFS(1, 7) =  5000.0

    print *, "Testing FRZLVL() with modelname != 'GFS' or 'RAPR'..."
    res = 0
    call FRZLVL(ZFRZ, RHFRZ, PFRZL)

    do i = 1, npts
        if (abs(ZFRZ(1, i) - EXP_ZFRZ(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected ZFRZ = ", EXP_ZFRZ(1, i), &
                     " but got ZFRZ = ", ZFRZ(1, i)
            res = 1
        end if
        if (abs(RHFRZ(1, i) - EXP_RHFRZ(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected RHFRZ = ", EXP_RHFRZ(1, i), &
                     " but got RHFRZ = ", RHFRZ(1, i)
            res = 1
        end if
        if (abs(PFRZL(1, i) - EXP_PFRZL(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected PFRZL = ", EXP_PFRZL(1, i), &
                     " but got PFRZL = ", PFRZL(1, i)
            res = 1
        end if
    end do

    if (res .ne. 0) then
        deallocate(pint, t, zmid, q, pmid, fis, tshltr, pshltr, qshltr, lmh)
        stop 10
    end if

    print *, "Testing FRZLVL() with modelname == 'GFS'..."
    modelname = "GFS"
    res = 0
    call FRZLVL(ZFRZ, RHFRZ, PFRZL)

    do i = 1, npts
        if (abs(ZFRZ(1, i) - EXP_ZFRZ_GFS(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected ZFRZ = ", EXP_ZFRZ_GFS(1, i), &
                     " but got ZFRZ = ", ZFRZ(1, i)
            res = 1
        end if
        if (abs(RHFRZ(1, i) - EXP_RHFRZ_GFS(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected RHFRZ = ", EXP_RHFRZ_GFS(1, i), &
                     " but got RHFRZ = ", RHFRZ(1, i)
            res = 1
        end if
        if (abs(PFRZL(1, i) - EXP_PFRZL_GFS(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected PFRZL = ", EXP_PFRZL_GFS(1, i), &
                     " but got PFRZL = ", PFRZL(1, i)
            res = 1
        end if
    end do

    if (res .ne. 0) then
        deallocate(pint, t, zmid, q, pmid, fis, tshltr, pshltr, qshltr, lmh)
        stop 20
    end if

    print *, "Testing FRZLVL() with modelname == 'RAPR'. Will give same result as 'GFS'..."
    modelname = "RAPR"
    res = 0
    call FRZLVL(ZFRZ, RHFRZ, PFRZL)

    do i = 1, npts
        if (abs(ZFRZ(1, i) - EXP_ZFRZ_GFS(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected ZFRZ = ", EXP_ZFRZ_GFS(1, i), &
                     " but got ZFRZ = ", ZFRZ(1, i)
            res = 1
        end if
        if (abs(RHFRZ(1, i) - EXP_RHFRZ_GFS(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected RHFRZ = ", EXP_RHFRZ_GFS(1, i), &
                     " but got RHFRZ = ", RHFRZ(1, i)
            res = 1
        end if
        if (abs(PFRZL(1, i) - EXP_PFRZL_GFS(1, i)) > tol) then
            print *, "Test Case ", i, " FAILED: Expected PFRZL = ", EXP_PFRZL_GFS(1, i), &
                     " but got PFRZL = ", PFRZL(1, i)
            res = 1
        end if
    end do
    
    if (res .ne. 0) stop 30

    deallocate(pint, t, zmid, q, pmid, fis, tshltr, pshltr, qshltr, lmh)
    print *, "SUCCESS!"
end program test_frzlvl