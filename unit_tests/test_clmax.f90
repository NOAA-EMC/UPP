! This is a test program for UPP.
!
! This program tests the CLMAX() subroutine.
!
! Alyson Stahl, 1/2026
program test_clmax
    use vrbls3d, only: zint, q2, pint
    use masks, only: lmh, sm
    use params_mod, only: EPSQ2
    use ctlblk_mod, only: jsta, jend, lm, spval, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 2, nlevs = 3
    integer :: i, j, res
    real, dimension(1:npts,1:npts) :: SQZ, SQ, RQ2L, RQ2H
    real, dimension(1:npts,1:npts) :: EL0, EXP_EL0

    interface
        subroutine CLMAX(EL0, SQZ, SQ, RQ2L, RQ2H)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, dimension(ista:iend,jsta:jend), intent(inout) :: &
                EL0, SQZ, SQ, RQ2L, RQ2H
        end subroutine CLMAX
    end interface
    
    ! Grid parameters
    jsta = 1
    jend = npts
    ista = 1
    iend = npts
    lm = npts
    spval = 9.9e6

    ! Allocate arrays
    allocate(lmh(1:npts,1:npts))
    allocate(q2(1:npts,1:npts,1:nlevs))
    allocate(zint(1:npts,1:npts,1:nlevs+1))
    allocate(pint(1:npts,1:npts,1:nlevs+1))
    allocate(sm(1:npts,1:npts))

    ! Expected results
    EXP_EL0 = reshape([11.0, 300.0, spval, 5.430163193E+01], [npts, npts])

    lmh = 1.0        
    q2  = 0.4        
    sm = 0.0

    ! Interface heights (m), monotonically increasing: 0, 100, 500, 1500
    zint(:,:,1) = 0.0
    zint(:,:,2) = 100.0
    zint(:,:,3) = 500.0
    zint(:,:,4) = 1500.0

    ! Interface pressures (Pa), decreasing with height (scale height ≈ 8 km)
    pint(:,:,1) = 100000.0
    pint(:,:,2) = 98760.0
    pint(:,:,3) = 93940.0
    pint(:,:,4) = 83470.0

    ! Test Case: q2 <= EPSQ2 at all levels. Should expect EL0 = ELMIN
    q2(1,1,:) = EPSQ2 * 0.5

    ! Test Case: Expect EL0 = spval
    zint(1,2,:) = spval

    ! Test Case: Expect EL0 = EL0M
    zint(2,1,1) = 0.0
    zint(2,1,2) = 2000.0
    zint(2,1,3) = 4000.0
    zint(2,1,4) = 6000.0

    call CLMAX(EL0, SQZ, SQ, RQ2L, RQ2H)

    deallocate(lmh)
    deallocate(q2)
    deallocate(zint)
    deallocate(pint)
    deallocate(sm)

    ! Check results
    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(EL0(i,j) - EXP_EL0(i,j)) > tol) then
                print *, 'EL0 Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_EL0(i,j), &
                         ' but got ', EL0(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    print *, "SUCCESS!"
end program test_clmax
