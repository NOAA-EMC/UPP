! This is a test program for UPP.
!
! This program tests the CALWXT_DOMINANT_POST() subroutine.
!
! Alyson Stahl, 12/2025
program test_calwxt_dominant
    use ctlblk_mod, only: jsta, jend, pthresh, im, jsta_2l, jend_2u, &
                          ista, iend, ista_2l, iend_2u
    implicit none

    integer, parameter :: nalg = 5 ! Defined in CALWXT_DOMINANT_POST() subroutine
    integer, parameter :: npts = 3
    integer :: i, j, res
    real, dimension(1:npts,1:npts) :: PREC, DOMS, DOMR, DOMZR, DOMIP
    real, dimension(1:npts,1:npts,nalg) :: RAIN, SNOW, SLEET, FREEZR
    real, dimension(1:npts,1:npts) :: EXP_DOMS, EXP_DOMR, EXP_DOMZR, EXP_DOMIP

    interface
        subroutine CALWXT_DOMINANT_POST(PREC, RAIN, FREEZR, SLEET, SNOW, DOMR, &
                                        DOMZR, DOMIP, DOMS)
            use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u, &
                                  ista, iend, ista_2l, iend_2u
            integer, parameter :: nalg = 5
            real, intent(in) :: PREC(ista_2l:iend_2u,jsta_2l:jend_2u)
            real, dimension(ista:iend,jsta:jend,nalg), intent(in) :: &
                RAIN, SNOW, SLEET, FREEZR
            real, dimension(ista:iend,jsta:jend), intent(inout) ::  &
                DOMS, DOMR, DOMZR, DOMIP
        end subroutine CALWXT_DOMINANT_POST
    end interface

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
    pthresh = 1.0e-6

    ! Initialize input arrays
    PREC = 1.0e-4 
    RAIN = 0.0
    FREEZR = 0.0
    SLEET = 0.0
    SNOW = 0.0

    ! Test Case: PREC <= PTHRESH
    PREC(1,1) = pthresh 

    ! (2,1): TOTSN > TOTIP and TOTSN > TOTZR and TOTSN >= TOTR -> DOMS
    SNOW(2,1,1) = 1.0
    SNOW(2,1,2) = 1.0
    SNOW(2,1,3) = 1.0
    SLEET(2,1,4) = 1.0
    RAIN(2,1,5) = 1.0

    ! (3,1): TOTSN > TOTIP and TOTSN > TOTZR and TOTSN < TOTR -> DOMR
    SNOW(3,1,1) = 1.0
    RAIN(3,1,2) = 1.0
    RAIN(3,1,3) = 1.0

    ! (1,2): TOTSN > TOTIP and TOTSN <= TOTZR and TOTZR >= TOTR -> DOMZR
    SNOW(1,2,1) = 1.0
    FREEZR(1,2,2) = 1.0
    FREEZR(1,2,3) = 1.0
    RAIN(1,2,4) = 1.0

    ! (2,2): TOTSN > TOTIP and TOTSN <= TOTZR and TOTZR < TOTR -> DOMR
    SNOW(2,2,1) = 1.0
    FREEZR(2,2,2) = 1.0
    RAIN(2,2,3) = 1.0
    RAIN(2,2,4) = 1.0

    ! (3,2): TOTSN <= TOTIP and TOTIP > TOTZR and TOTIP >= TOTR -> DOMIP
    SLEET(3,2,1) = 1.0
    SLEET(3,2,2) = 1.0
    FREEZR(3,2,3) = 1.0
    RAIN(3,2,4) = 1.0

    ! (1,3): TOTSN <= TOTIP and TOTIP > TOTZR and TOTIP < TOTR -> DOMR
    SLEET(1,3,1) = 1.0
    RAIN(1,3,2) = 1.0
    RAIN(1,3,3) = 1.0

    ! (2,3): TOTSN <= TOTIP and TOTIP <= TOTZR and TOTZR >= TOTR -> DOMZR
    SNOW(2,3,1) = 1.0
    SLEET(2,3,2) = 1.0
    FREEZR(2,3,3) = 1.0
    FREEZR(2,3,4) = 1.0
    RAIN(2,3,5) = 1.0

    ! (3,3): TOTSN <= TOTIP and TOTIP <= TOTZR and TOTZR < TOTR -> DOMR
    SNOW(3,3,1) = 1.0
    SLEET(3,3,2) = 1.0
    FREEZR(3,3,3) = 1.0
    RAIN(3,3,4) = 1.0
    RAIN(3,3,5) = 1.0

    ! Default output is 0.0
    EXP_DOMS = 0.0
    EXP_DOMR = 0.0
    EXP_DOMZR = 0.0
    EXP_DOMIP = 0.0

    ! Expected dominant type per case
    EXP_DOMS(2,1) = 1.0
    EXP_DOMR(3,1) = 1.0
    EXP_DOMZR(1,2) = 1.0
    EXP_DOMR(2,2) = 1.0
    EXP_DOMIP(3,2) = 1.0
    EXP_DOMR(1,3) = 1.0
    EXP_DOMZR(2,3) = 1.0
    EXP_DOMR(3,3) = 1.0

    call CALWXT_DOMINANT_POST(PREC, RAIN, FREEZR, SLEET, SNOW, DOMR, DOMZR, DOMIP, DOMS)

    res = 0
    do i = 1, npts
        do j = 1, npts
            if (DOMS(i,j) /= EXP_DOMS(i,j)) then
                print *, 'DOMS Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_DOMS(i,j), &
                         ' but got ', DOMS(i,j)
                res = 1
            end if
            if (DOMR(i,j) /= EXP_DOMR(i,j)) then
                print *, 'DOMR Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_DOMR(i,j), &
                         ' but got ', DOMR(i,j)
                res = 1
            end if
            if (DOMZR(i,j) /= EXP_DOMZR(i,j)) then
                print *, 'DOMZR Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_DOMZR(i,j), &
                         ' but got ', DOMZR(i,j)
                res = 1
            end if
            if (DOMIP(i,j) /= EXP_DOMIP(i,j)) then
                print *, 'DOMIP Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_DOMIP(i,j), &
                         ' but got ', DOMIP(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10
    
    print *, "SUCCESS!"
end program test_calwxt_dominant