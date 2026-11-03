! This is a test program for UPP.
!
! This program tests the CALCAPE() subroutine.
!
! Alyson Stahl, 8/2026
program test_calcape
    use upp_physics, only: CALCAPE
    use table_upp_mod, only: TABLE
    use tableq_upp_mod, only: TABLEQ
    use lookup_mod, only: thl, plq, ptbl, ttbl, rdq, rdth, rdp, rdthe, pl, qs0, sqs, sthe, &
                                                the0, ttblq, rdpq, rdtheq, stheq, the0q
    use vrbls3d,    only: pmid, t, q, zint
    use vrbls2d,    only: teql,ieql,tshltr,pshltr,qshltr
    use masks,      only: lmh
    use ctlblk_mod, only: jsta_2l, jend_2u, lm, jsta, jend, spval, &
                                                ista_2l, iend_2u, ista, iend, capecin_2m, pt
    implicit none

    real, parameter :: tol = 1.0e-6
    integer :: res

    ! Adjust the coarse pressure table so the lowest interpolated saturation
    ! pressure can reach the TPSPK guard in this test.
    pt = 0.0
    thl = 210.0
    plq = 70000.0

    call TABLE(ptbl, ttbl, pt, rdq, rdth, rdp, rdthe, pl, thl, qs0, sqs, sthe, the0)
    call TABLEQ(ttblq, rdpq, rdtheq, plq, thl, stheq, the0q)

    res = 0

    ! Test Case 1: ITYPE = 1
    call test_itype_1(res)
    if (res .ne. 0) stop 10

    ! Test Case 2: ITYPE = 2
    call test_itype_2(res)
    if (res .ne. 0) stop 20

    print *, "SUCCESS!"

contains

    subroutine test_itype_1(res)
        integer, intent(inout) :: res
        integer, parameter :: ny = 5, nz = 4
        integer :: itype
        real :: dpbnd
        integer, dimension(1,ny) :: l1d
        real, dimension(1,ny)  :: p1d, t1d
        real, dimension(1,ny) :: q1d, cape, cins, pparc, zeql, thund
        real, dimension(1,ny) :: exp_q1d, exp_cape, exp_cins, exp_pparc, exp_zeql, exp_thund
        integer :: j, k

        ista = 1
        iend = 1
        jsta = 1
        jend = ny
        ista_2l = ista
        iend_2u = iend
        jsta_2l = jsta
        jend_2u = jend
        lm = nz
        spval = 9.9e10

        allocate(pmid(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(t(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(q(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(zint(ista_2l:iend_2u,jsta_2l:jend_2u,lm+1))
        allocate(teql(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(ieql(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(tshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(pshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(qshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(lmh(ista_2l:iend_2u,jsta_2l:jend_2u))

        itype = 1
        dpbnd = 25000.0
        capecin_2m = .true.

        pmid = spval
        t = spval
        q = spval
        zint = spval
        teql = spval
        ieql = 0
        tshltr = spval
        pshltr = spval
        qshltr = spval
        lmh = real(lm)

        do j = 1, ny
            pmid(1,j,:) = (/ 50000.0, 65000.0, 80000.0, 95000.0 /)
            zint(1,j,:) = (/ 6500.0, 4500.0, 2500.0, 1000.0, 0.0 /) 
            pshltr(1,j) = 96000.0
        end do
        
        ! Test Case 1: Unstable boundary-layer parcel meant to support 
        ! buoyancy ascent and nonzero CAPE
        t(1,1,:) = (/ 250.0, 262.0, 286.0, 301.0 /)
        q(1,1,:) = (/ 3.0e-4, 1.0e-3, 8.0e-3, 1.5e-2 /)
        tshltr(1,1) = 303.0
        qshltr(1,1) = 1.6e-2

        ! Test Case 2: Stable boundary-layer parcel meant to suppress buoyancy.
        t(1,2,:) = (/ 247.0, 255.0, 263.0, 271.0 /)
        q(1,2,:) = (/ 1.0e-4, 4.0e-4, 1.2e-3, 3.0e-3 /)
        tshltr(1,2) = 272.0
        qshltr(1,2) = 3.0e-3

        ! Test Case 3: Stress case meant to push potential temp below lookup table
        ! lower bound.
        t(1,3,:) = (/ 180.0, 190.0, 200.0, 205.0 /)
        q(1,3,:) = (/ 0.0, 0.0, 0.0, 0.0 /)
        tshltr(1,3) = 180.0
        qshltr(1,3) = 0.0

        ! Test Case 4: Stress case meant to push potential temp above lookup table
        ! upper bound.
        t(1,4,:) = (/ 360.0, 375.0, 390.0, 400.0 /)
        q(1,4,:) = (/ 2.0e-3, 3.0e-3, 4.0e-3, 5.0e-3 /)
        tshltr(1,4) = 400.0
        qshltr(1,4) = 5.0e-3

        ! Test Case 5: Stress case meant to exercise guard for very low TPSPK value.
        t(1,5,:) = (/ 180.0, 190.0, 200.0, 205.0 /)
        q(1,5,:) = (/ 0.0, 0.0, 0.0, 0.0 /)
        tshltr(1,5) = 180.0
        qshltr(1,5) = 0.0

        ! Dummy variables for ITYPE = 1
        p1d = spval
        t1d = spval
        q1d = spval

        l1d = lm
        
        cape = spval
        cins = spval
        pparc = spval
        zeql = spval
        thund = spval

        exp_q1d = q1d

        call CALCAPE(itype, dpbnd, p1d, t1d, q1d, l1d, cape, cins, pparc, zeql, thund)

        do j = 1, ny
            print *, "cape(1,", j, ") = ", cape(1,j)
        end do

        do j = 1, ny
            print *, "cins(1,", j, ") = ", cins(1,j)
        end do

        do j = 1, ny
            print *, "pparc(1,", j, ") = ", pparc(1,j)
        end do

        do j = 1, ny
            print *, "zeql(1,", j, ") = ", zeql(1,j)
        end do

        do j = 1, ny
            print *, "thund(1,", j, ") = ", thund(1,j)
        end do

        deallocate(pmid)
        deallocate(t)
        deallocate(q)
        deallocate(zint)
        deallocate(teql)
        deallocate(ieql)
        deallocate(tshltr)
        deallocate(pshltr)
        deallocate(qshltr)
        deallocate(lmh)

    end subroutine test_itype_1

    subroutine test_itype_2(res)
        integer, intent(inout) :: res
        integer, parameter :: ny = 7, nz = 4
        integer :: itype
        real :: dpbnd
        integer, dimension(1,ny) :: l1d
        real, dimension(1,ny)  :: p1d, t1d
        real, dimension(1,ny) :: q1d, cape, cins, pparc, zeql, thund
        real, dimension(1,ny) :: exp_q1d, exp_cape, exp_cins, exp_pparc, exp_zeql, exp_thund
        integer :: j, k

        ista = 1
        iend = 1
        jsta = 1
        jend = ny
        ista_2l = ista
        iend_2u = iend
        jsta_2l = jsta
        jend_2u = jend
        lm = nz
        spval = 9.9e10

        allocate(pmid(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(t(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(q(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(zint(ista_2l:iend_2u,jsta_2l:jend_2u,lm+1))
        allocate(teql(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(ieql(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(tshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(pshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(qshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(lmh(ista_2l:iend_2u,jsta_2l:jend_2u))

        itype = 2
        dpbnd = 25000.0
        capecin_2m = .false.

        pmid = spval
        t = spval
        q = spval
        zint = spval
        teql = spval
        ieql = 0
        tshltr = spval
        pshltr = spval
        qshltr = spval
        lmh = real(lm)

        do j = 1, ny
            pmid(1,j,:) = (/ 50000.0, 65000.0, 80000.0, 95000.0 /)
            zint(1,j,:) = (/ 6500.0, 4500.0, 2500.0, 1000.0, 0.0 /)
            t(1,j,:) = (/ 250.0, 262.0, 286.0, 301.0 /)
            q(1,j,:) = (/ 3.0e-4, 1.0e-3, 8.0e-3, 1.5e-2 /)
            pshltr(1,j) = 96000.0
            tshltr(1,j) = 303.0
            qshltr(1,j) = 1.6e-2
        end do

        p1d = 95000.0

        ! Test Case 1: Typical warm, moist parcel for ITYPE = 2.
        t1d(1,1) = 301.0
        q1d(1,1) = 1.5e-2

        ! Test Case 2: Cooler, drier parcel meant to be less buoyant.
        t1d(1,2) = 271.0
        q1d(1,2) = 3.0e-3
        t(1,2,:) = (/ 247.0, 255.0, 263.0, 271.0 /)
        q(1,2,:) = (/ 1.0e-4, 4.0e-4, 1.2e-3, 3.0e-3 /)
        tshltr(1,2) = 272.0
        qshltr(1,2) = 3.0e-3

        ! Test Case 3: Force the low-end Q1D clamp.
        t1d(1,3) = 295.0
        q1d(1,3) = -1.0

        ! Test Case 4: Force the high-end Q1D clamp.
        t1d(1,4) = 295.0
        q1d(1,4) = 1.0e6

        ! Test Case 5: Stress case meant to push ITTBK below the table range.
        t1d(1,5) = 180.0
        q1d(1,5) = 1.0e-3

        ! Test Case 6: Stress case meant to push ITTBK above the table range.
        t1d(1,6) = 400.0
        q1d(1,6) = 5.0e-3

        ! Test Case 7: Stress case with a very low parcel pressure.
        p1d(1,7) = 500.0
        t1d(1,7) = 220.0
        q1d(1,7) = 1.0e-12

        l1d = lm

        cape = spval
        cins = spval
        pparc = spval
        zeql = spval
        thund = spval

        call CALCAPE(itype, dpbnd, p1d, t1d, q1d, l1d, cape, cins, pparc, zeql, thund)

        do j = 1, ny
            print *, "cape(1,", j, ") = ", cape(1,j)
        end do

        do j = 1, ny
            print *, "cins(1,", j, ") = ", cins(1,j)
        end do

        do j = 1, ny
            print *, "pparc(1,", j, ") = ", pparc(1,j)
        end do

        do j = 1, ny
            print *, "zeql(1,", j, ") = ", zeql(1,j)
        end do

        do j = 1, ny
            print *, "thund(1,", j, ") = ", thund(1,j)
        end do

        deallocate(pmid)
        deallocate(t)
        deallocate(q)
        deallocate(zint)
        deallocate(teql)
        deallocate(ieql)
        deallocate(tshltr)
        deallocate(pshltr)
        deallocate(qshltr)
        deallocate(lmh)

    end subroutine test_itype_2
    
end program test_calcape