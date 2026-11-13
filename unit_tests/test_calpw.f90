! This is a test program for UPP.
!
! This program tests the CALPW() subroutine.
!
! Alyson Stahl, 11/2026
program test_calpw
    use vrbls3d,    only: q, qqw, qqi, qqr, qqs, cwm, qqg, t, rswtt,    &
                        train, tcucn, mcvg, pmid, o3, ext, pint, rlwtt, &
                        taod5503d,sca, asy, zmid
    use vrbls4d,    only: smoke, fv3dust, coarsepm
    use masks,      only: htm
    use ctlblk_mod, only: lm, jsta, jend, spval, ista, iend
    use params_mod, only: tfrz
    implicit none

    real, parameter :: tol = 1.0e-6
    integer :: res

    interface 
        subroutine CALPW(PW,IDECID)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, intent(inout) :: PW(ista:iend,jsta:jend)
            integer, intent(in) :: IDECID
        end subroutine CALPW
    end interface

    res = 0
    call test_idecid_le_1(res)
    ! call test_idecid_2(res)
    ! call test_idecid_3(res)
    ! call test_idecid_4(res)
    ! call test_idecid_5(res)
    ! call test_idecid_6(res)
    ! call test_idecid_7(res)
    ! call test_idecid_8(res)
    ! call test_idecid_9(res)
    ! call test_idecid_10(res)
    ! call test_idecid_11(res)
    ! call test_idecid_12(res)
    ! call test_idecid_13(res)
    ! call test_idecid_14(res)
    ! call test_idecid_15(res)
    ! call test_idecid_17(res)
    ! call test_idecid_18(res)
    ! call test_idecid_19(res)
    ! call test_idecid_20(res)
    ! call test_idecid_21(res)
    ! call test_idecid_22(res)
    ! call test_idecid_23(res)
    ! call test_idecid_24(res)
    ! call test_idecid_25(res)

    print *, 'SUCCESS!'

contains
    subroutine test_idecid_le_1(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(q(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        q = 0.005
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        q(1,3,:) = spval
        EXP_PW(1,3) = spval

        ! Testing with IDECID = 1
        IDECID = 1
        PW = 0.0
        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 1. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        ! Testing with IDECID < 1. Should give same result as IDECID = 1.
        IDECID = 0
        PW = 0.0
        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID < 1. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(q)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_le_1

    subroutine test_idecid_2(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 2
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(qqw(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        qqw = 0.005
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        qqw(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 2. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(qqw)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_2

    subroutine test_idecid_3(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 3
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(qqi(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        qqi = 0.005
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        qqi(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 3. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(qqi)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_3

    subroutine test_idecid_4(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 4
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(qqr(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        qqr = 0.005
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        qqr(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 4. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(qqr)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_4

    subroutine test_idecid_5(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 5
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(qqs(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        qqs = 0.005
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        qqs(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 5. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(qqs)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_5

    subroutine test_idecid_6(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 6
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(cwm(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        cwm = 0.04
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        cwm(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 6. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(cwm)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_6

    subroutine test_idecid_7(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 5
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 7
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(t(ista:iend, jsta:jend, 1:lm))
        allocate(qqw(ista:iend, jsta:jend, 1:lm))
        allocate(qqr(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        t = tfrz - 10.0
        qqw = 0.005
        qqr = 0.003
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case for T < TFRZ. No spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: Representative case for T >= TFRZ. No spval input or output.
        t(1,2,:) = tfrz + 10.0
        EXP_PW(1,2) = 0.0

        ! Test Case 3: PINT >= spval
        pint(1,3,:) = spval
        EXP_PW(1,3) = spval

        ! Test Case 4: QQW >= spval
        qqw(1,4,:) = spval
        EXP_PW(1,4) = spval

        ! Test Case 5: QQR >= spval
        qqr(1,5,:) = spval
        EXP_PW(1,5) = spval
    
        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 7. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(t)
        deallocate(qqw)
        deallocate(qqr)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_7

    subroutine test_idecid_8(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 5
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 8
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(t(ista:iend, jsta:jend, 1:lm))
        allocate(qqi(ista:iend, jsta:jend, 1:lm))
        allocate(qqs(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        t = tfrz + 10.0
        qqi = 0.004
        qqs = 0.002
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case for T > TFRZ. No spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: Representative case for T <= TFRZ. No spval input or output.
        t(1,2,:) = tfrz - 10.0
        EXP_PW(1,2) = 0.0

        ! Test Case 3: PINT >= spval
        pint(1,3,:) = spval
        EXP_PW(1,3) = spval

        ! Test Case 4: QQI >= spval
        qqi(1,4,:) = spval
        EXP_PW(1,4) = spval

        ! Test Case 5: QQS >= spval
        qqs(1,5,:) = spval
        EXP_PW(1,5) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 8. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(t)
        deallocate(qqi)
        deallocate(qqs)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_8

    subroutine test_idecid_9(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 9
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(rswtt(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        rswtt = 0.002
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        rswtt(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 9. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(rswtt)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_9

    subroutine test_idecid_10(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 10
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(rlwtt(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        rlwtt = -0.001
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        rlwtt(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 10. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(rlwtt)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_10

    subroutine test_idecid_11(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 11
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(train(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        train = 0.005
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        train(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 11. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(train)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_11

    subroutine test_idecid_12(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 12
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(tcucn(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        tcucn = -0.005
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        tcucn(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 12. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(tcucn)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_12

    subroutine test_idecid_13(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 13
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(mcvg(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        mcvg = 0.0004
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        mcvg(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 13. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(mcvg)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_13

    subroutine test_idecid_14(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 14
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(q(ista:iend, jsta:jend, 1:lm))
        allocate(t(ista:iend, jsta:jend, 1:lm))
        allocate(pmid(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        q = 0.005
        t = 280.0
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do
        do L = 1, lm
            pmid(:,:,L) = 0.5 * (pint(:,:,L) + pint(:,:,L+1))
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        q(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 14. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(q)
        deallocate(t)
        deallocate(pmid)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_14

    subroutine test_idecid_15(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 15
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(o3(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        o3 = 0.00004
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        o3(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 15. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(o3)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_15

    subroutine test_idecid_16(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 16
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(qqg(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        qqg = 0.005
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        qqg(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 16. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(qqg)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_16

    subroutine test_idecid_17(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 17
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(ext(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        ext = 24.0
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        ext(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 17. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(ext)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_17

    subroutine test_idecid_18(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 2
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 18
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(smoke(ista:iend, jsta:jend, 1:lm, 1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        smoke = 10.0
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 18. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(pint)
        deallocate(smoke)
        deallocate(htm)
    end subroutine test_idecid_18

    subroutine test_idecid_19(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 19
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(taod5503d(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))

        taod5503d = 0.1
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        taod5503d(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 19. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(taod5503d)
        deallocate(pint)
    end subroutine test_idecid_19

    subroutine test_idecid_20(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 20
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(sca(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        sca = 10.0
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        sca(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 20. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(sca)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_20

    subroutine test_idecid_21(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 21
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(asy(ista:iend, jsta:jend, 1:lm))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        asy = 35.0
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        ! Test Case 3: Qdum >= spval
        asy(1,3,:) = spval
        EXP_PW(1,3) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 21. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(asy)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_21

    subroutine test_idecid_22(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 2
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 22
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(fv3dust(ista:iend, jsta:jend, 1:lm, 1))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        fv3dust = 10.0
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 22. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(fv3dust)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_22

    subroutine test_idecid_23(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 2
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 23
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(coarsepm(ista:iend, jsta:jend, 1:lm, 1))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        coarsepm = 10.0
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 23. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(coarsepm)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_23

    subroutine test_idecid_24(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 2
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 24
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(coarsepm(ista:iend, jsta:jend, 1:lm, 1))
        allocate(fv3dust(ista:iend, jsta:jend, 1:lm, 1))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        fv3dust = 10.0
        coarsepm = 10.0
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 24. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(coarsepm)
        deallocate(fv3dust)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_24

    subroutine test_idecid_25(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 2
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 25
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(coarsepm(ista:iend, jsta:jend, 1:lm, 1))
        allocate(fv3dust(ista:iend, jsta:jend, 1:lm, 1))
        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(zmid(ista:iend, jsta:jend, 1:lm))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        fv3dust = 10.0
        coarsepm = 10.0
        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do
        do L = 1, lm
            zmid(:,:,L) = 1000.0 * L
        end do

        ! Test Case 1: Representative case with no spval input or output.
        EXP_PW(1,1) = 0.0

        ! Test Case 2: PINT >= spval
        pint(1,2,:) = spval
        EXP_PW(1,2) = spval

        call CALPW(PW, IDECID)

        do J = 1, ntests
            if (abs(PW(1,J) - EXP_PW(1,J)) > tol) then
                print *,  "Test case ", J, " failed for IDECID == 25. Expected PW = ", PW(1,J), &
                        "but got PW = ", PW
                res = 1
            end if
        end do

        deallocate(coarsepm)
        deallocate(fv3dust)
        deallocate(zmid)
        deallocate(pint)
        deallocate(htm)
    end subroutine test_idecid_25

    subroutine test_idecid_undefined(res)
        integer, intent(inout) :: res
        integer, parameter :: ntests = 3
        integer, parameter :: nlevs = 3
        !
        integer :: IDECID, L, J
        real :: PW(1,1:ntests), EXP_PW(1,1:ntests)

        IDECID = 26
        PW = 0.0

        ista = 1
        iend = 1
        jsta = 1
        jend = ntests
        lm = nlevs
        spval = 9.9e10

        allocate(pint(ista:iend, jsta:jend, 1:lm+1))
        allocate(htm(ista:iend, jsta:jend, 1:lm))

        htm = 1.0
        do L = 1, lm+1
            pint(:,:,L) = 70000.0 + 10000.0 * (L - 1)
        end do

    end subroutine test_idecid_undefined
end program test_calpw