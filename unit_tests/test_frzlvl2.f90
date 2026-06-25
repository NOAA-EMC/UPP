! This is a test program for UPP.
!
! This program tests the FRZLVL2() subroutine.
!
! Alyson Stahl, 5/2026
program test_frzlvl2
      use vrbls3d, only: pint, t, zmid, pmid, q, zint, alpint
      use vrbls2d, only: fis, tshltr, pshltr, qz0, qs, qshltr
      use masks, only: lmh, sm
      use ctlblk_mod, only: jsta, jend, spval, lm, modelname, im, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: npts = 7, nlevs = 30
    integer :: i, j, k, res
    real :: ISOTHERM_1, ISOTHERM_2
    real :: ZFRZ(1, npts), RHFRZ(1, npts), PFRZL(1, npts)
    real :: EXP_ZFRZ_1(1, npts), EXP_RHFRZ_1(1, npts), EXP_PFRZL_1(1, npts)
    real :: EXP_ZFRZ_2(1, npts), EXP_RHFRZ_2(1, npts), EXP_PFRZL_2(1, npts)
    real :: EXP_ZFRZ_GFS_1(1, npts), EXP_RHFRZ_GFS_1(1, npts), EXP_PFRZL_GFS_1(1, npts)
    real :: EXP_ZFRZ_GFS_2(1, npts), EXP_RHFRZ_GFS_2(1, npts), EXP_PFRZL_GFS_2(1, npts)

    interface 
        subroutine FRZLVL2(ISOTHERM,ZFRZ,RHFRZ,PFRZL)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, intent(in) :: ISOTHERM
            real, dimension(ista:iend,jsta:jend), intent(out) :: ZFRZ, RHFRZ, PFRZL
        end subroutine FRZLVL2
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
    allocate(zint(ista:iend, jsta:jend, 1:lm+1))
    allocate(alpint(ista:iend, jsta:jend, 1:lm+1))
    allocate(fis(ista:iend, jsta:jend)) 
    allocate(tshltr(ista:iend, jsta:jend))
    allocate(pshltr(ista:iend, jsta:jend))
    allocate(qz0(ista:iend, jsta:jend))
    allocate(qs(ista:iend, jsta:jend))
    allocate(qshltr(ista:iend, jsta:jend))
    allocate(lmh(ista:iend, jsta:jend))
    allocate(sm(ista:iend, jsta:jend))


    ! This unit test will call the FRZLVL2() subroutine 6 times:
    ! - modelname == ""
    ! - modelname == "GFS"
    ! - modelname == "RAPR"
    ! For each of the 3 modelname cases, there will be 2 test cases:
    ! 1) ISOTHERM_1 = 230 K (< TSFC)
    ! 2) ISOTHERM_2 = 245 K (> TSFC)
    ISOTHERM_1 = 230.0
    ISOTHERM_2 = 245.0

    ! Test Case 1: Default case where none of the input values are spval and LICE == LLMH.
    lmh = real(lm)
    fis(1, :) = 0.0
    sm(1, :) = 0.5
    tshltr(1, :) = 355.0
    pshltr(1, :) = 24400.0
    qshltr(1, :) = 0.001
    qz0(1, :) = 0.001
    qs(1, :) = 0.001
    
    do k = 1, lm + 1
        pint(1, :, k) = 5000.0 + real(lm + 1 - k) / real(lm) * 20000.0
        zint(1, :, k) = real(lm + 1 - k) / real(lm) * 18600.0
        alpint(1, :, k) = log(pint(1, :, k))
    end do
    
    do k = 1, lm
        pmid(1, :, k) = 0.5 * (pint(1, :, k) + pint(1, :, k+1))
        zmid(1, :, k) = 0.5 * (zint(1, :, k) + zint(1, :, k+1))
        t(1, :, k) = 238.0 - 20.0 * real(lm - k) / real(lm - 1)
        q(1, :, k) = 0.001 * exp(-real(lm - k) / real(lm - 1) * 2.0)
    end do
    
    EXP_ZFRZ_1(1, 1) = 0.000000000E+00
    EXP_RHFRZ_1(1, 1) = 6.934353113E-01
    EXP_PFRZL_1(1, 1) =  5.666664551E+03

    EXP_ZFRZ_GFS_1(1, 1) = 0.000000000E+00
    EXP_RHFRZ_GFS_1(1, 1) = 1.000000000E+00
    EXP_PFRZL_GFS_1(1, 1) = 5.666664551E+03

    EXP_ZFRZ_2(1, 1) = 0.000000000E+00
    EXP_RHFRZ_2(1, 1) = 1.000000000E+00
    EXP_PFRZL_2(1, 1) =  2.440000000E+04

    EXP_ZFRZ_GFS_2(1, 1) = 0.000000000E+00
    EXP_RHFRZ_GFS_2(1, 1) = 1.000000000E+00
    EXP_PFRZL_GFS_2(1, 1) = 2.440000000E+04


    ! Test Case 2: QSHLTR == spval with other default input values.
    qshltr(1, 2) = spval

    EXP_ZFRZ_1(1, 2) = 0.000000000E+00
    EXP_RHFRZ_1(1, 2) = 6.934353113E-01
    EXP_PFRZL_1(1, 2) =  5.666664551E+03

    EXP_ZFRZ_GFS_1(1, 2) = 0.000000000E+00
    EXP_RHFRZ_GFS_1(1, 2) = 1.000000000E+00
    EXP_PFRZL_GFS_1(1, 2) = 5.666664551E+03

    EXP_ZFRZ_2(1, 2) = 0.000000000E+00
    EXP_RHFRZ_2(1, 2) = 3.072937727E-01
    EXP_PFRZL_2(1, 2) =  5.333333008E+03

    EXP_ZFRZ_GFS_2(1, 2) = 0.000000000E+00
    EXP_RHFRZ_GFS_2(1, 2) = 4.265226424E-01
    EXP_PFRZL_GFS_2(1, 2) = 5.333333008E+03

    ! Test Case 3: QSHLTR, TSHLTR, and PSHLTR are spval with other default input values.
    qshltr(1, 3) = spval
    tshltr(1, 3) = spval
    pshltr(1, 3) = spval

    EXP_ZFRZ_1(1, 3) = 1.540770020E+03
    EXP_RHFRZ_1(1, 3) = 6.934353113E-01
    EXP_PFRZL_1(1, 3) =  5.666664551E+03

    EXP_ZFRZ_GFS_1(1, 3) = 1.540770020E+03
    EXP_RHFRZ_GFS_1(1, 3) = 1.000000000E+00
    EXP_PFRZL_GFS_1(1, 3) = 5.666664551E+03

    EXP_ZFRZ_2(1, 3) = 0.000000000E+00
    EXP_RHFRZ_2(1, 3) = 2.320170552E-01
    EXP_PFRZL_2(1, 3) =  5.333333008E+03

    EXP_ZFRZ_GFS_2(1, 3) = 0.000000000E+00
    EXP_RHFRZ_GFS_2(1, 3) = 3.142645955E-01
    EXP_PFRZL_GFS_2(1, 3) = 5.333333008E+03

    ! Test Case 4: TSHLTR and PSHLTR are spval with other default input values.
    tshltr(1, 4) = spval
    pshltr(1, 4) = spval

    EXP_ZFRZ_1(1, 4) = 1.540770020E+03
    EXP_RHFRZ_1(1, 4) = 6.934353113E-01
    EXP_PFRZL_1(1, 4) =  5.666664551E+03

    EXP_ZFRZ_GFS_1(1, 4) = 1.540770020E+03
    EXP_RHFRZ_GFS_1(1, 4) = 1.000000000E+00
    EXP_PFRZL_GFS_1(1, 4) = 5.666664551E+03

    EXP_ZFRZ_2(1, 4) = 0.000000000E+00
    EXP_RHFRZ_2(1, 4) = 1.000000000E+00
    EXP_PFRZL_2(1, 4) =  9.900000051E+10

    EXP_ZFRZ_GFS_2(1, 4) = 0.000000000E+00
    EXP_RHFRZ_GFS_2(1, 4) = 1.000000000E+00
    EXP_PFRZL_GFS_2(1, 4) = 9.900000051E+10

    ! Test Case 5: LICE clipped to L < LLMH due to PMID(1,5,L) >= PUCAP
    pint(1, 5, 15) = 35000.0
    pint(1, 5, 16) = 32000.0
    pmid(1, 5, 15) = 0.5 * (pint(1, 5, 15) + pint(1, 5, 16))
    alpint(1, 5, 15) = log(pint(1, 5, 15))
    alpint(1, 5, 16) = log(pint(1, 5, 16))
    t(1, 5, 15) = 228.0
    t(1, 5, 16) = 232.0

    EXP_ZFRZ_1(1, 5) = 9.300000000E+03
    EXP_RHFRZ_1(1, 5) = 1.000000000E+00
    EXP_PFRZL_1(1, 5) =  2.239790430E+04

    EXP_ZFRZ_GFS_1(1, 5) = 9.300000000E+03
    EXP_RHFRZ_GFS_1(1, 5) = 1.000000000E+00
    EXP_PFRZL_GFS_1(1, 5) = 2.239790430E+04

    EXP_ZFRZ_2(1, 5) = 0.000000000E+00
    EXP_RHFRZ_2(1, 5) = 1.000000000E+00
    EXP_PFRZL_2(1, 5) =  2.440000000E+04

    EXP_ZFRZ_GFS_2(1, 5) = 0.000000000E+00
    EXP_RHFRZ_GFS_2(1, 5) = 1.000000000E+00
    EXP_PFRZL_GFS_2(1, 5) = 2.440000000E+04

    ! Test Case 6: Previous LICE clipped to L < LLMH case with QSHLTR, TSHLTR, and PSHLTR all spval.
    pint(1, 6, 15) = 35000.0
    pint(1, 6, 16) = 32000.0
    pmid(1, 6, 15) = 0.5 * (pint(1, 6, 15) + pint(1, 6, 16))
    alpint(1, 6, 15) = log(pint(1, 6, 15))
    alpint(1, 6, 16) = log(pint(1, 6, 16))
    t(1, 6, 15) = 228.0
    t(1, 6, 16) = 232.0
    qshltr(1, 6) = spval
    tshltr(1, 6) = spval
    pshltr(1, 6) = spval

    EXP_ZFRZ_1(1, 6) = 9.300000000E+03
    EXP_RHFRZ_1(1, 6) = 1.000000000E+00
    EXP_PFRZL_1(1, 6) =  2.239790430E+04

    EXP_ZFRZ_GFS_1(1, 6) = 9.300000000E+03
    EXP_RHFRZ_GFS_1(1, 6) = 1.000000000E+00
    EXP_PFRZL_GFS_1(1, 6) = 2.239790430E+04

    EXP_ZFRZ_2(1, 6) = 0.000000000E+00
    EXP_RHFRZ_2(1, 6) = 2.320170552E-01
    EXP_PFRZL_2(1, 6) =  5.333333008E+03

    EXP_ZFRZ_GFS_2(1, 6) = 0.000000000E+00
    EXP_RHFRZ_GFS_2(1, 6) = 3.142645955E-01
    EXP_PFRZL_GFS_2(1, 6) = 5.333333008E+03

    ! Test Case 7: FIS is spval
    fis(1, 7) = spval

    EXP_ZFRZ_1(1,7) = spval 
    EXP_RHFRZ_1(1,7) = spval 
    EXP_PFRZL_1(1,7) = 0.0

    EXP_ZFRZ_2(1,7) = spval 
    EXP_RHFRZ_2(1,7) = spval 
    EXP_PFRZL_2(1,7) = 0.0
    EXP_ZFRZ_GFS_1(1,7) = spval 
    EXP_RHFRZ_GFS_1(1,7) = spval 
    EXP_PFRZL_GFS_1(1,7) = 0.0
    EXP_ZFRZ_GFS_2(1,7) = spval 
    EXP_RHFRZ_GFS_2(1,7) = spval 
    EXP_PFRZL_GFS_2(1,7) = 0.0

    print *, "Testing FRZLVL2() with modelname != 'GFS' and ISOTHERM < TSFC..."
    res = 0
    ! Set default outputs
    ZFRZ = 0.0
    RHFRZ = 0.0
    PFRZL = 0.0
    call FRZLVL2(ISOTHERM_1, ZFRZ, RHFRZ, PFRZL)

    do i = 1, npts
        if (abs(ZFRZ(1, i) - EXP_ZFRZ_1(1, i)) > tol) then
            print *, "ZFRZ(1, ", i, ") = ", ZFRZ(1, i), " does not match expected value of ", EXP_ZFRZ_1(1, i)
            res = 1
        end if
        if (abs(RHFRZ(1, i) - EXP_RHFRZ_1(1, i)) > tol) then
            print *, "RHFRZ(1, ", i, ") = ", RHFRZ(1, i), " does not match expected value of ", EXP_RHFRZ_1(1, i)
            res = 1
        end if
        if (abs(PFRZL(1, i) - EXP_PFRZL_1(1, i)) > tol) then
            print *, "PFRZL(1, ", i, ") = ", PFRZL(1, i), " does not match expected value of ", EXP_PFRZL_1(1, i)
            res = 1
        end if
    end do

    if (res .ne. 0) then
        deallocate(pint, t, zmid, q, pmid, zint, alpint, fis, tshltr, pshltr, qz0, qs, qshltr, lmh, sm)
        stop 10
    end if

    print *, "Testing FRZLVL2() with modelname != 'GFS' and ISOTHERM > TSFC..."
    ! Set default outputs
    ZFRZ = 0.0
    RHFRZ = 0.0
    PFRZL = 0.0
    call FRZLVL2(ISOTHERM_2, ZFRZ, RHFRZ, PFRZL)

    do i = 1, npts
        if (abs(ZFRZ(1, i) - EXP_ZFRZ_2(1, i)) > tol) then
            print *, "ZFRZ(1, ", i, ") = ", ZFRZ(1, i), " does not match expected value of ", EXP_ZFRZ_2(1, i)
            res = 1
        end if
        if (abs(RHFRZ(1, i) - EXP_RHFRZ_2(1, i)) > tol) then
            print *, "RHFRZ(1, ", i, ") = ", RHFRZ(1, i), " does not match expected value of ", EXP_RHFRZ_2(1, i)
            res = 1
        end if
        if (abs(PFRZL(1, i) - EXP_PFRZL_2(1, i)) > tol) then
            print *, "PFRZL(1, ", i, ") = ", PFRZL(1, i), " does not match expected value of ", EXP_PFRZL_2(1, i)
            res = 1
        end if
    end do

    if (res .ne. 0) then
        deallocate(pint, t, zmid, q, pmid, zint, alpint, fis, tshltr, pshltr, qz0, qs, qshltr, lmh, sm)
        stop 20
    end if

    print *, "Testing FRZLVL2() with modelname == 'GFS' and ISOTHERM < TSFC..."
    modelname = "GFS"
    ! Set default outputs
    ZFRZ = 0.0
    RHFRZ = 0.0
    PFRZL = 0.0
    call FRZLVL2(ISOTHERM_1, ZFRZ, RHFRZ, PFRZL)

    do i = 1, npts
        if (abs(ZFRZ(1, i) - EXP_ZFRZ_GFS_1(1, i)) > tol) then
            print *, "ZFRZ(1, ", i, ") = ", ZFRZ(1, i), " does not match expected value of ", EXP_ZFRZ_GFS_1(1, i)
            res = 1
        end if
        if (abs(RHFRZ(1, i) - EXP_RHFRZ_GFS_1(1, i)) > tol) then
            print *, "RHFRZ(1, ", i, ") = ", RHFRZ(1, i), " does not match expected value of ", EXP_RHFRZ_GFS_1(1, i)
            res = 1
        end if
        if (abs(PFRZL(1, i) - EXP_PFRZL_GFS_1(1, i)) > tol) then
            print *, "PFRZL(1, ", i, ") = ", PFRZL(1, i), " does not match expected value of ", EXP_PFRZL_GFS_1(1, i)
            res = 1
        end if
    end do

    if (res .ne. 0) then
        deallocate(pint, t, zmid, q, pmid, zint, alpint, fis, tshltr, pshltr, qz0, qs, qshltr, lmh, sm)
        stop 30
    end if

    print *, "Testing FRZLVL2() with modelname == 'GFS' and ISOTHERM > TSFC..."
    ! Set default outputs
    ZFRZ = 0.0
    RHFRZ = 0.0
    PFRZL = 0.0
    call FRZLVL2(ISOTHERM_2, ZFRZ, RHFRZ, PFRZL)

    do i = 1, npts
        if (abs(ZFRZ(1, i) - EXP_ZFRZ_GFS_2(1, i)) > tol) then
            print *, "ZFRZ(1, ", i, ") = ", ZFRZ(1, i), " does not match expected value of ", EXP_ZFRZ_GFS_2(1, i)
            res = 1
        end if
        if (abs(RHFRZ(1, i) - EXP_RHFRZ_GFS_2(1, i)) > tol) then
            print *, "RHFRZ(1, ", i, ") = ", RHFRZ(1, i), " does not match expected value of ", EXP_RHFRZ_GFS_2(1, i)
            res = 1
        end if
        if (abs(PFRZL(1, i) - EXP_PFRZL_GFS_2(1, i)) > tol) then
            print *, "PFRZL(1, ", i, ") = ", PFRZL(1, i), " does not match expected value of ", EXP_PFRZL_GFS_2(1, i)
            res = 1
        end if
    end do

    if (res .ne. 0) then
        deallocate(pint, t, zmid, q, pmid, zint, alpint, fis, tshltr, pshltr, qz0, qs, qshltr, lmh, sm)
        stop 40
    end if
    
    print *, "Testing FRZLVL2() with modelname == 'RAPR' and ISOTHERM < TSFC. Will give same result as modelname == 'GFS' case..."
    modelname = "RAPR"
    ! Set default outputs
    ZFRZ = 0.0
    RHFRZ = 0.0
    PFRZL = 0.0
    call FRZLVL2(ISOTHERM_1, ZFRZ, RHFRZ, PFRZL)

    do i = 1, npts
        if (abs(ZFRZ(1, i) - EXP_ZFRZ_GFS_1(1, i)) > tol) then
            print *, "ZFRZ(1, ", i, ") = ", ZFRZ(1, i), " does not match expected value of ", EXP_ZFRZ_GFS_1(1, i)
            res = 1
        end if
        if (abs(RHFRZ(1, i) - EXP_RHFRZ_GFS_1(1, i)) > tol) then
            print *, "RHFRZ(1, ", i, ") = ", RHFRZ(1, i), " does not match expected value of ", EXP_RHFRZ_GFS_1(1, i)
            res = 1
        end if
        if (abs(PFRZL(1, i) - EXP_PFRZL_GFS_1(1, i)) > tol) then
            print *, "PFRZL(1, ", i, ") = ", PFRZL(1, i), " does not match expected value of ", EXP_PFRZL_GFS_1(1, i)
            res = 1
        end if
    end do

    if (res .ne. 0) then
        deallocate(pint, t, zmid, q, pmid, zint, alpint, fis, tshltr, pshltr, qz0, qs, qshltr, lmh, sm)
        stop 50
    end if

    print *, "Testing FRZLVL2() with modelname == 'RAPR' and ISOTHERM > TSFC. Will give same result as modelname == 'GFS' case..."
    ! Set default outputs
    ZFRZ = 0.0
    RHFRZ = 0.0
    PFRZL = 0.0
    call FRZLVL2(ISOTHERM_2, ZFRZ, RHFRZ, PFRZL)

    do i = 1, npts
        if (abs(ZFRZ(1, i) - EXP_ZFRZ_GFS_2(1, i)) > tol) then
            print *, "ZFRZ(1, ", i, ") = ", ZFRZ(1, i), " does not match expected value of ", EXP_ZFRZ_GFS_2(1, i)
            res = 1
        end if
        if (abs(RHFRZ(1, i) - EXP_RHFRZ_GFS_2(1, i)) > tol) then
            print *, "RHFRZ(1, ", i, ") = ", RHFRZ(1, i), " does not match expected value of ", EXP_RHFRZ_GFS_2(1, i)
            res = 1
        end if
        if (abs(PFRZL(1, i) - EXP_PFRZL_GFS_2(1, i)) > tol) then
            print *, "PFRZL(1, ", i, ") = ", PFRZL(1, i), " does not match expected value of ", EXP_PFRZL_GFS_2(1, i)
            res = 1
        end if
    end do
    
    if (res .ne. 0) then
        deallocate(pint, t, zmid, q, pmid, zint, alpint, fis, tshltr, pshltr, qz0, qs, qshltr, lmh, sm)
        stop 60
    end if

    deallocate(pint, t, zmid, q, pmid, zint, alpint, fis, tshltr, pshltr, qz0, qs, qshltr, lmh, sm)

    print *, "SUCCESS!"
end program test_frzlvl2