! This is a test program for UPP.
!
! This program tests the MDL2SIGMA2() subroutine.
!
! Alyson Stahl, 3/2026
program test_mdl2sigma2
    use vrbls3d, only:  pint, pmid, t, zint, q
    use masks, only: lmh
    use params_mod, only: pq0, a2, a3, a4, rgamog
    use ctlblk_mod, only: pt, jsta_2l, jend_2u, spval, lp1, lm, jsta, jend,&
                        grib, cfld, datapd, fld_info, &
                        ista, iend, ista_2l, iend_2u
    use rqstfld_mod, only: iget, lvls, iavblfld, lvlsxml
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: nx = 2, ny = 2, nlevs = 3
    integer :: i, j, k, res
    integer :: EXP_IFLD = 296, EXP_LVL = 1000
    real :: EXP_DATAPD(1:nx, 1:ny)

    interface
        subroutine MDL2SIGMA2()
        end subroutine MDL2SIGMA2
    end interface

    ! Grid Dimensions
    ista     = 1
    iend     = nx
    ista_2l  = 1
    iend_2u  = nx
    jsta     = 1
    jend     = ny
    jsta_2l  = 1
    jend_2u  = ny
    lm       = nlevs
    lp1      = nlevs + 1

    spval = 9.9e10
    grib = 'grib2'
    pt   = 1.0e4
    cfld = 0

    allocate(pint(ista_2l:iend_2u, jsta_2l:jend_2u, 1:lm+1))
    allocate(pmid(ista_2l:iend_2u, jsta_2l:jend_2u, 1:lm))
    allocate(t(ista_2l:iend_2u, jsta_2l:jend_2u, 1:lm))
    allocate(zint(ista_2l:iend_2u, jsta_2l:jend_2u, 1:lm+1))
    allocate(q(ista_2l:iend_2u, jsta_2l:jend_2u, 1:lm))
    allocate(lmh(ista:iend, jsta:jend))
    allocate(datapd(iend-ista+1, jend-jsta+1, 1))
    allocate(fld_info(1))
    allocate(lvlsxml(5, 1))

    cfld = 0

    do i = 1, nx
        do j = 1, ny
            do k = 1, lm+1
                pint(i, j, k) = pt + (real(k-1)/real(lm))*(1.0e5 - pt)
                zint(i, j, k) = real(k-1)*1000.0
            end do
            do k = 1, lm
                pmid(i, j, k) = 0.5*(pint(i, j, k) + pint(i, j, k+1))
                t(i, j, k) = 290.0 - 6.0*real(k-1)
                q(i, j, k) = 0.010/real(k)
            end do
            lmh(i, j) = real(lm)
        end do
    end do

    iget(:) = 0
    lvls(:,:) = 0
    lvls(1:5,1) = 0
    lvls(1,1) = 1000
    iavblfld(1) = 296
    lvlsxml(:,:) = 0
    lvlsxml(1,1) = 1000

    datapd(:,:,:) = 0.0

    ! Set to some default values.
    fld_info(1)%ifld = 9999
    fld_info(1)%lvl = 9999

    ! Testing with IGET(296) = 0. Should skip entire subroutine.
    call MDL2SIGMA2()

    res = 0
    if (cfld .ne. 0) then
        print *, "Expected cfld = 0, got ", cfld
        res = 1
    end if
    if (fld_info(1)%ifld .ne. 9999) then
        print *, "Expected fld_info(1)%ifld = 9999, got ", fld_info(1)%ifld
        res = 1
    end if
    if (fld_info(1)%lvl .ne. 9999) then
        print *, "Expected fld_info(1)%lvl = 9999, got ", fld_info(1)%lvl
        res = 1
    end if

    do i = 1, nx
        do j = 1, ny
            if (datapd(i, j, 1) .ne. 0.0) then
                print *, "Expected datapd(", i, ",", j, ",1) = 0.0, got ", datapd(i, j, 1)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    ! Testing with IGET(296) = 1. Should execute subroutine and populate datapd, fld_info.
    iget(296) = 1
    EXP_DATAPD(1,1) =  280.09765625

    ! Test Case: Sigma pressure deeper than the model bottom interface. 0.01 <= RHL <= 1.0
    pmid(1, 2, :) = 0.5*pt
    EXP_DATAPD(1,2) = 321.78915405

    ! Test Case: Sigma pressure deeper than the model bottom interface. RHL < 0.01
    pmid(2, 1, :) = 0.5*pt
    q(2, 1, 1:2) = 1.0e-10
    EXP_DATAPD(2,1) =  321.78912354

    ! Test Case: Sigma pressure deeper than the model bottom interface. RHL > 1.0
    ! Also meets condition where TMT0 < -20.0.
    pmid(2, 2, :) = 0.5*pt
    q(2, 2, 1:2) = 1.0
    t(2, 2, 1) = 200.0
    t(2, 2, 2) = 200.0
    EXP_DATAPD(2,2) = 224.24328613

    call MDL2SIGMA2()

    res = 0
    if (cfld .ne. 1) then
        print *, "Expected cfld = 1, got ", cfld
        res = 1
    end if
    if (fld_info(1)%ifld .ne. EXP_IFLD) then
        print *, "Expected fld_info(1)%ifld = ", EXP_IFLD, " got ", fld_info(1)%ifld
        res = 1
    end if
    if (fld_info(1)%lvl .ne. EXP_LVL) then
        print *, "Expected fld_info(1)%lvl = ", EXP_LVL, " got ", fld_info(1)%lvl
        res = 1
    end if
    
    do i = 1, nx
        do j = 1, ny
            if (abs(datapd(i, j, 1) - EXP_DATAPD(i, j)) > tol) then
                print *, "Expected datapd(", i, ",", j, ",1) = ", EXP_DATAPD(i, j), " got ", datapd(i, j, 1)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 20

    deallocate(pint)
    deallocate(pmid)
    deallocate(t)
    deallocate(zint)
    deallocate(q)
    deallocate(lmh)
    deallocate(datapd)
    deallocate(fld_info)
    deallocate(lvlsxml)

    print *, 'SUCCESS!'
end program test_mdl2sigma2