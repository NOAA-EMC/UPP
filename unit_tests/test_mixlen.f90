! This is a test program for UPP.
!
! This program tests the MIXLEN() subroutine.
!
! Alyson Stahl, 3/2026
program test_mixlen
    use vrbls3d, only: zint, pmid, t, q2
    use masks, only: lmh, htm
    use ctlblk_mod, only: jsta, jend, jsta_m, jend_m, jsta_2l, jend_2u, lm, &
                          lm1, spval, ista, iend, ista_m, iend_m, ista_2l, iend_2u
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: nx = 2, ny = 2, nlevs = 3
    integer :: i, j, k, res
    real :: EL0(nx, ny), EL(nx, ny, nlevs), EXP_EL(nx, ny, nlevs)
    
    interface
        subroutine MIXLEN(EL0,EL)
            use ctlblk_mod, only: ista_2l, iend_2u, jsta_2l, jend_2u, lm
            real, intent(in) :: EL0(ista_2l:iend_2u,jsta_2l:jend_2u)
            real, intent(out) :: EL(ista_2l:iend_2u,jsta_2l:jend_2u,lm)
        end subroutine MIXLEN
    end interface

    ! Grid parameters
    ista = 1
    iend = nx
    ista_m = 1
    iend_m = nx
    ista_2l = 1
    iend_2u = nx
    jsta = 1
    jend = ny
    jsta_m = 1
    jend_m = ny
    jsta_2l = 1
    jend_2u = ny
    lm = nlevs
    lm1 = nlevs - 1
    spval = 9.9e10

    allocate(zint(ista_2l:iend_2u, jsta_2l:jend_2u, 1:lm+1))
    allocate(pmid(ista_2l:iend_2u, jsta_2l:jend_2u, 1:lm))
    allocate(t(ista_2l:iend_2u, jsta_2l:jend_2u, 1:lm))
    allocate(q2(ista_2l:iend_2u, jsta_2l:jend_2u, 1:lm))
    allocate(htm(ista:iend, jsta:jend, 1:lm))
    allocate(lmh(ista:iend, jsta:jend))
    
    lmh(:,:) = real(lm)
    htm(:,:,:) = 1.0
    EL0(:,:) = 100.0

    do i = 1, nx
        do j = 1, ny
            do k = 1, lm + 1
                zint(i, j, k) = real(lm + 1 - k) * 100.0
            end do
        end do
    end do

    do i = 1, nx
        do j = 1, ny
            do k = 1, lm
                pmid(i, j, k) = 1.0e5 - real(k - 1) * 5.0e3
                t(i, j, k)    = 290.0 - 5.0 * real(k - 1)
                q2(i, j, k)   = 0.1
            end do
        end do
    end do

    EXP_EL(:,:,1) = 10.445438385
    EXP_EL(:,:,2) = 11.520626068
    EXP_EL(:,:,3) = 0.0

    ! Test Case: HGT(1,2) = ZINT(1,2,lm+1) = spval and T(1,2,1) = spval
    ! Expect EL(1,2,1) = EL(1,2,lm) = spval
    ZINT(1,2,lm+1) = spval
    T(1,2,1) = spval

    EXP_EL(1,2,1) = spval
    EXP_EL(1,2,2) = 540.0
    EXP_EL(1,2,3) = spval

    ! Test Case: ZIAG >= CPBLT * EL0(2,1)
    EL0(2,1)     = 1.0
    ZINT(2,1,1)  = 150.0
    ZINT(2,1,2)  = 40.0
    ZINT(2,1,3)  = 10.0
    HTM(2,1,2)   = 0.0
    HTM(2,1,3)   = 0.0

    EXP_EL(2,1,1) = 16.0
    EXP_EL(2,1,2) = 4.0
    EXP_EL(2,1,3) = 0.0

    call MIXLEN(EL0, EL)

    deallocate(zint)
    deallocate(pmid)
    deallocate(t)
    deallocate(q2)
    deallocate(htm)
    deallocate(lmh)

    res = 0
    do i = 1, nx
        do j = 1, ny
            do k = 1, lm
                if (abs(EL(i, j, k) - EXP_EL(i, j, k)) > tol) then
                    print *, 'Test failed at (', i, j, k, '): EL = ', EL(i, j, k), &
                             ' but expected ', EXP_EL(i, j, k)
                    res = 1
                end if
            end do
        end do
    end do

    if (res .ne. 0) stop 10
    
    print *, 'SUCCESS!'
end program test_mixlen