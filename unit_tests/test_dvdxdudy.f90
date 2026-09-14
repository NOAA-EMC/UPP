! This is a test program for UPP.
!
! This program tests the dvdxdudy() subroutine.
!
! Alyson Stahl, 8/2026
program test_dvdxdudy
    use upp_math, only: dvdxdudy, DDVDX, DDUDY, UUAVG
    use masks,        only: dx, dy
    use ctlblk_mod,   only: jsta_2l, jend_2u, jsta_m, jend_m, spval,&
                                ista_2l, iend_2u, ista_m, iend_m
    use gridspec_mod, only: gridtype
    implicit none

    real, parameter :: tol = 1.0e-6
    ! Minimum grid size must be 3x3. Size chosen for some variety in input values.
    integer, parameter :: NX = 5, NY = 5
    integer :: i, j, res
    real :: UWND(1:NX,1:NY), VWND(1:NX,1:NY)
    real :: EXP_DDVDX(1:NX,1:NY), EXP_DDUDY(1:NX,1:NY), EXP_UUAVG(1:NX,1:NY)

    ! Grid dimensions
    ista_2l = 1
    iend_2u = NX
    ista_m = 2
    iend_m = NX - 1
    jsta_2l = 1
    jend_2u = NY
    jsta_m = 2
    jend_m = NY - 1
    spval = 9.9e10

    allocate(dx(ista_2l:iend_2u,jsta_2l:jend_2u))
    allocate(dy(ista_2l:iend_2u,jsta_2l:jend_2u))
    allocate(DDVDX(ista_2l:iend_2u,jsta_2l:jend_2u))
    allocate(DDUDY(ista_2l:iend_2u,jsta_2l:jend_2u))
    allocate(UUAVG(ista_2l:iend_2u,jsta_2l:jend_2u))

    ! Test Case 1: Typical case for grid type A without SPVAL in UWND or VWND. 
    ! Expect calculated values for DDVDX, DDUDY, & UUAVG except where dx < 1e-5 or
    ! dy < 1e-5, in which case SPVAL is expected.

    gridtype = 'A'
    dx = 1000.00
    dy = 2000.00
    EXP_DDVDX = 0.0
    EXP_DDUDY = 0.0
    EXP_UUAVG = 0.0

    do i = 1, NX
        do j = 1, NY
            UWND(i,j) = 12.00 + 1.50 * real(i - 1) + 6.00 * real(j - 1)
            VWND(i,j) = 8.00 + 4.00 * real(i - 1) + 0.50 * real(j - 1)
        end do
    end do

    do i = 2, NX-1
        do j = 2, NY-1
            EXP_DDVDX(i,j) = 0.004
            EXP_DDUDY(i,j) = 0.003 
            EXP_UUAVG(i,j) = UWND(i,j)
        end do
    end do

    dx(2,2) = 1.0e-6
    EXP_DDVDX(2,2) = spval
    EXP_DDUDY(2,2) = spval
    EXP_UUAVG(2,2) = spval

    dy(2,3) = 1.0e-6
    EXP_DDVDX(2,3) = spval
    EXP_DDUDY(2,3) = spval
    EXP_UUAVG(2,3) = spval

    call dvdxdudy(UWND, VWND)

    res = 0
    do i = 1, NX
        do j = 1, NY
            if (abs(EXP_DDVDX(i,j) - DDVDX(i,j)) > tol) then
                print *, "Expected DDVDX(", i, ",", j, ") = ", EXP_DDVDX(i,j), &
                         ", but got ", DDVDX(i,j)
                res = 1
            end if
            if (abs(EXP_DDUDY(i,j) - DDUDY(i,j)) > tol) then
                print *, "Expected DDUDY(", i, ",", j, ") = ", EXP_DDUDY(i,j), &
                         ", but got ", DDUDY(i,j)
                res = 1
            end if
            if (abs(EXP_UUAVG(i,j) - UUAVG(i,j)) > tol) then
                print *, "Expected UUAVG(", i, ",", j, ") = ", EXP_UUAVG(i,j), &
                         ", but got ", UUAVG(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    ! Test Case 2: Grid type A with SPVAL. Should return all SPVAL.
    UWND = spval
    VWND = spval
    dx = 1000.0
    dy = 2000.0
    do i = 2, NX-1
        do j = 2, NY-1
            EXP_DDVDX(i,j) = spval
            EXP_DDUDY(i,j) = spval
            EXP_UUAVG(i,j) = spval
        end do
    end do

    call dvdxdudy(UWND, VWND)

    res = 0
    do i = 1, NX
        do j = 1, NY
            if (abs(EXP_DDVDX(i,j) - DDVDX(i,j)) > tol) then
                print *, "Expected DDVDX(", i, ",", j, ") = ", EXP_DDVDX(i,j), &
                         ", but got ", DDVDX(i,j)
                res = 1
            end if
            if (abs(EXP_DDUDY(i,j) - DDUDY(i,j)) > tol) then
                print *, "Expected DDUDY(", i, ",", j, ") = ", EXP_DDUDY(i,j), &
                         ", but got ", DDUDY(i,j)
                res = 1
            end if
            if (abs(EXP_UUAVG(i,j) - UUAVG(i,j)) > tol) then
                print *, "Expected UUAVG(", i, ",", j, ") = ", EXP_UUAVG(i,j), &
                         ", but got ", UUAVG(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 20

    ! Test Case 3: Typical case for grid type E without SPVAL in UWND or VWND. 
    ! Expect calculated values for DDVDX, DDUDY, & UUAVG

    gridtype = 'E'
    dx = 1000.00
    dy = 2000.00
    do i = 1, NX
        do j = 1, NY
            UWND(i,j) = 12.00 + 1.50 * real(i - 1) + 6.00 * real(j - 1)
            VWND(i,j) = 8.00 + 4.00 * real(i - 1) + 0.50 * real(j - 1)
        end do
    end do

    do i = 2, NX-1
        do j = 2, NY-1
            EXP_DDVDX(i,j) = 0.002
            EXP_DDUDY(i,j) = 0.003
            if (MOD(j,2) == 0) then
                EXP_UUAVG(i,j) = UWND(i,j) + 0.375
            else
                EXP_UUAVG(i,j) = UWND(i,j) - 0.375
            end if
        end do
    end do

    call dvdxdudy(UWND, VWND)

    res = 0
    do i = 1, NX
        do j = 1, NY
            if (abs(EXP_DDVDX(i,j) - DDVDX(i,j)) > tol) then
                print *, "Expected DDVDX(", i, ",", j, ") = ", EXP_DDVDX(i,j), &
                         ", but got ", DDVDX(i,j)
                res = 1
            end if
            if (abs(EXP_DDUDY(i,j) - DDUDY(i,j)) > tol) then
                print *, "Expected DDUDY(", i, ",", j, ") = ", EXP_DDUDY(i,j), &
                         ", but got ", DDUDY(i,j)
                res = 1
            end if
            if (abs(EXP_UUAVG(i,j) - UUAVG(i,j)) > tol) then
                print *, "Expected UUAVG(", i, ",", j, ") = ", EXP_UUAVG(i,j), &
                         ", but got ", UUAVG(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 30

    ! Test Case 4: Grid type E with SPVAL. Should return all SPVAL.
    UWND = spval
    VWND = spval
    dx = 1000.0
    dy = 2000.0
    do i = 2, NX-1
        do j = 2, NY-1
            EXP_DDVDX(i,j) = spval
            EXP_DDUDY(i,j) = spval
            EXP_UUAVG(i,j) = spval
        end do
    end do

    call dvdxdudy(UWND, VWND)

    res = 0
    do i = 1, NX
        do j = 1, NY
            if (abs(EXP_DDVDX(i,j) - DDVDX(i,j)) > tol) then
                print *, "Expected DDVDX(", i, ",", j, ") = ", EXP_DDVDX(i,j), &
                         ", but got ", DDVDX(i,j)
                res = 1
            end if
            if (abs(EXP_DDUDY(i,j) - DDUDY(i,j)) > tol) then
                print *, "Expected DDUDY(", i, ",", j, ") = ", EXP_DDUDY(i,j), &
                         ", but got ", DDUDY(i,j)
                res = 1
            end if
            if (abs(EXP_UUAVG(i,j) - UUAVG(i,j)) > tol) then
                print *, "Expected UUAVG(", i, ",", j, ") = ", EXP_UUAVG(i,j), &
                         ", but got ", UUAVG(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 40

    ! Test Case 5: Typical case for grid type B without SPVAL in UWND or VWND. 
    ! Expect calculated values for DDVDX, DDUDY, & UUAVG

    gridtype = 'B'
    dx = 1000.00
    dy = 2000.00
    do i = 1, NX
        do j = 1, NY
            UWND(i,j) = 12.00 + 1.50 * real(i - 1) + 6.00 * real(j - 1)
            VWND(i,j) = 8.00 + 4.00 * real(i - 1) + 0.50 * real(j - 1)
        end do
    end do


    do i = 2, NX-1
        do j = 2, NY-1
            EXP_DDVDX(i,j) = 0.004
            EXP_DDUDY(i,j) = 0.003
            EXP_UUAVG(i,j) = UWND(i,j) - 3.75
        end do
    end do

    call dvdxdudy(UWND, VWND)

    res = 0
    do i = 1, NX
        do j = 1, NY
            if (abs(EXP_DDVDX(i,j) - DDVDX(i,j)) > tol) then
                print *, "Expected DDVDX(", i, ",", j, ") = ", EXP_DDVDX(i,j), &
                         ", but got ", DDVDX(i,j)
                res = 1
            end if
            if (abs(EXP_DDUDY(i,j) - DDUDY(i,j)) > tol) then
                print *, "Expected DDUDY(", i, ",", j, ") = ", EXP_DDUDY(i,j), &
                         ", but got ", DDUDY(i,j)
                res = 1
            end if
            if (abs(EXP_UUAVG(i,j) - UUAVG(i,j)) > tol) then
                print *, "Expected UUAVG(", i, ",", j, ") = ", EXP_UUAVG(i,j), &
                         ", but got ", UUAVG(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 50

    ! Test Case 6: Grid type B with SPVAL. Should return all SPVAL.
    UWND = spval
    VWND = spval
    dx = 1000.0
    dy = 2000.0
    do i = 2, NX-1
        do j = 2, NY-1
            EXP_DDVDX(i,j) = spval
            EXP_DDUDY(i,j) = spval
            EXP_UUAVG(i,j) = spval
        end do
    end do

    call dvdxdudy(UWND, VWND)

    res = 0
    do i = 1, NX
        do j = 1, NY
            if (abs(EXP_DDVDX(i,j) - DDVDX(i,j)) > tol) then
                print *, "Expected DDVDX(", i, ",", j, ") = ", EXP_DDVDX(i,j), &
                         ", but got ", DDVDX(i,j)
                res = 1
            end if
            if (abs(EXP_DDUDY(i,j) - DDUDY(i,j)) > tol) then
                print *, "Expected DDUDY(", i, ",", j, ") = ", EXP_DDUDY(i,j), &
                         ", but got ", DDUDY(i,j)
                res = 1
            end if
            if (abs(EXP_UUAVG(i,j) - UUAVG(i,j)) > tol) then
                print *, "Expected UUAVG(", i, ",", j, ") = ", EXP_UUAVG(i,j), &
                         ", but got ", UUAVG(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 60

    print *, "SUCCESS!"
end program test_dvdxdudy