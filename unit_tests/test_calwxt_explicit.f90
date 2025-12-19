! This is a test program for UPP.
!
! This program tests the CALWXT_EXPLICIT_POST() subroutine.
!
! Alyson Stahl, 12/2025
program test_calwxt_explicit
    use ctlblk_mod, only: jsta, jend, pthresh, im, jsta_2l, jend_2u,  &
                          lm, ista, iend, ista_2l, iend_2u
    implicit none

    integer, parameter :: npts = 3 , nlevs = 1
    integer :: i, j, res
    real, dimension(1:npts,1:npts) :: LMH, THS, PREC, SR
    real, dimension(1:npts,1:npts,1:nlevs) :: PMID, F_RIMEF
    integer, dimension(1:npts,1:npts) :: IWX, EXP_IWX

    ! Grid parameters
    jsta = 1
    jsta_2l = 1
    jend = npts
    jend_2u = npts
    im = npts
    ista = 1
    iend = npts
    ista_2l = 1
    iend_2u = npts
    lm = nlevs
    pthresh = 1.0e-6

    ! Initialize arrays
    LMH(:,:) = 1.0            ! use lowest level everywhere
    PMID(:,:,1) = 100000.0    ! ~sea-level pressure in Pa
    THS(:,:) = 290.0          ! potential temperature in K
    SR(:,:) = 0.2             ! mostly rain branch (SR < 0.5)
    F_RIMEF(:,:,1) = 0.0      ! low riming by default
    PREC(:,:) = 1.0e-4        ! precip present (> pthresh)

    ! Case IWX = 0: no precip (PREC <= PTHRESH)
    PREC(1,1) = pthresh

    ! Case IWX = 4: freezing rain (SR < 0.5, Tskin < 273.15 K)
    THS(1,2) = 268.0

    ! Case IWX = 1: snow (SR >= 0.5, F_RIMEF < 10)
    SR(2,2) = 0.8
    F_RIMEF(2,2,1) = 0.0

    ! Case IWX = 2: sleet (SR >= 0.5, F_RIMEF >= 10)
    SR(2,3) = 0.8
    F_RIMEF(2,3,1) = 15.0

    EXP_IWX = reshape([0, 8, 8,  4, 1, 8,  8, 2, 8], [npts, npts])

    call CALWXT_EXPLICIT_POST(LMH, THS, PMID, PREC, SR, F_RIMEF, IWX)

    res = 0
    do j = jsta, jend
        do i = ista, iend
            if (IWX(i,j) /= EXP_IWX(i,j)) then
                print *, 'IWX Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_IWX(i,j), &
                         ' but got ', IWX(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    print *, 'SUCCESS!'
end program test_calwxt_explicit