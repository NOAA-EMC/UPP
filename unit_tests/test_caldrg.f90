! This is a test program for UPP.
!
! This program tests the CALDRG() subroutine.
!
! Alyson Stahl, 11/2025
program test_caldrg
    use vrbls3d, only: uh, vh
    use vrbls2d, only: uz0, vz0, ustar, u10, v10
    use masks, only: lmh
    use ctlblk_mod, only: jsta, jend, jsta_m, jend_m, modelname, spval, im, jm,  &
                        jsta_2l, jend_2u, ista, iend, ista_m, iend_m, ista_2l, iend_2u
    use gridspec_mod, only: gridtype
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: npts = 5
    integer, parameter :: nlevs = 60 ! arbitrarily larger than lmh for buffer
    integer :: i, j, res
    real :: DRAGCO(1:npts, 1:npts)

    ! Expected results for comparison
    real :: EXP_DRAGCO_A(1:npts, 1:npts), EXP_DRAGCO_E(1:npts, 1:npts), EXP_DRAGCO_B(1:npts, 1:npts),  &
            EXP_DRAGCO_B_NMM(1:npts, 1:npts), EXP_DRAGCO_INVALID(1:npts, 1:npts)

    interface
        subroutine CALDRG(DRAGCO)
            use ctlblk_mod, only: jsta_2l, jend_2u, ista_2l, iend_2u
            real, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(inout) :: DRAGCO
        end subroutine CALDRG
    end interface

    ! Grid dimensions and bounds
    im = npts
    jm = npts
    ista = 1
    iend = npts
    jsta = 1
    jend = npts
    ista_2l = 1
    iend_2u = npts
    jsta_2l = 1
    jend_2u = npts
    ista_m = 2
    iend_m = 4
    jsta_m = 2
    jend_m = 4

    spval = -9999.0

    modelname = 'GFS'

    ! Initialize expected results
    EXP_DRAGCO_A = reshape((/ 0.0000000000000, 0.0002941176471, 0.0002941176471, 0.0002941176471, 0.0002941176471, &
                            0.0002941176471, 0.0011764705882, 0.0073529411765, 0.0188235294118, 0.0002941176471, &
                            0.0002941176471, 0.0026470588235, 0.0105882352941, 0.0238235294118, 0.0002941176471, &
                            0.0002941176471, 0.0047058823529, 0.0144117647059, 0.0294117647059, 0.0002941176471, &
                            0.0002941176471, 0.0002941176471, 0.0002941176471, 0.0002941176471, 0.0002941176471 /), (/npts,npts/))
    EXP_DRAGCO_E = reshape((/ 0.0, 0.0, 0.0, 0.0, 0.0, &
                            0.0, 0.0007920792079, 0.0049504950495, 0.0126732673267, 0.0, &
                            0.0, 0.0017821782178, 0.0071287128713, 0.0160396039604, 0.0, &
                            0.0, 0.0031683168317, 0.0097029702970, 0.0198019801980, 0.0, &
                            0.0, 0.0, 0.0, 0.0, 0.0 /), (/npts,npts/))
    EXP_DRAGCO_B = reshape((/ 0.0, 0.0, 0.0, 0.0, 0.0, &
                            0.0, 0.0007920792079, 0.0049504950495, 0.0126732673267, 0.0, &
                            0.0, 0.0017821782178, 0.0071287128713, 0.0160396039604, 0.0, &
                            0.0, 0.0031683168317, 0.0097029702970, 0.0198019801980, 0.0, &
                            0.0, 0.0, 0.0, 0.0, 0.0 /), (/npts,npts/))
    EXP_DRAGCO_B_NMM = RESHAPE((/ 0.0, 0.0, 0.0, 0.0, 0.0, &
                            0.0, 0.0007920792079, 0.0049504950495, 0.0126732673267, 0.0, &
                            0.0, 0.0017821782178, 0.0071287128713, 0.0160396039604, 0.0, &
                            0.0, 0.0031683168317, 0.0097029702970, 0.0198019801980, 0.0, &
                            0.0, 0.0, 0.0, 0.0, 0.0 /), (/npts,npts/))
    EXP_DRAGCO_INVALID = spval

    ! Allocate arrays with full grid dimensions
    allocate(ustar(ista_2l:iend_2u, jsta_2l:jend_2u))
    allocate(u10(ista_2l:iend_2u, jsta_2l:jend_2u))
    allocate(v10(ista_2l:iend_2u, jsta_2l:jend_2u))
    allocate(lmh(ista_2l:iend_2u, jsta_2l:jend_2u))
    allocate(uh(ista_2l:iend_2u, jsta_2l:jend_2u, 1:nlevs))
    allocate(vh(ista_2l:iend_2u, jsta_2l:jend_2u, 1:nlevs))
    allocate(uz0(ista_2l:iend_2u, jsta_2l:jend_2u))
    allocate(vz0(ista_2l:iend_2u, jsta_2l:jend_2u))

    lmh = 40.0           ! Model levels

    !ustar = reshape([0.3, 0.2, spval, 0.5, 0.4, 0.3, 0.1, 0.6, 0.7], [npts,npts])

    ustar = 0.1
    ustar(2,2) = 0.2
    ustar(2,3) = 0.3
    ustar(2,4) = 0.4
    ustar(3,2) = 0.5
    ustar(3,3) = 0.6
    ustar(3,4) = 0.7
    ustar(4,2) = 0.8
    ustar(4,3) = 0.9
    ustar(4,4) = 1.0
    ustar(1,1) = spval

    uh = 8.0
    vh = 6.0
    u10 = 5.0
    v10 = 3.0
    uz0 = 3.0
    vz0 = 3.0

    print *, "Testing CALDRG subroutine..."
    print *, "Testing with gridtype = 'A'."

    gridtype = 'A'
    call CALDRG(DRAGCO)

    ! Compare results for gridtype 'A'
    res = 0
    do j = jsta, jend
        do i = ista, iend
            if (abs(DRAGCO(i,j) - EXP_DRAGCO_A(i,j)) > tol) then
                print *, "Mismatch at (", i, ",", j, "): ", DRAGCO(i,j), " vs ", EXP_DRAGCO_A(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    print *, "Testing with gridtype = 'E'."
    gridtype = 'E'
    call CALDRG(DRAGCO)

    ! Compare results for gridtype 'E'
    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(DRAGCO(i,j) - EXP_DRAGCO_E(i,j)) > tol) then
                print *, "Mismatch at (", i, ",", j, "): ", DRAGCO(i,j), " vs ", EXP_DRAGCO_E(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 20

    print *, "Testing with gridtype = 'B'."

    gridtype = 'B'
    call CALDRG(DRAGCO)

    ! Compare results for gridtype 'B'
    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(DRAGCO(i,j) - EXP_DRAGCO_B(i,j)) > tol) then
                print *, "Mismatch at (", i, ",", j, "): ", DRAGCO(i,j), " vs ", EXP_DRAGCO_B(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 30

    print *, "Testing with gridtype = 'B' and modelname = 'NMM'."

    modelname = 'NMM'
    call CALDRG(DRAGCO)

    ! Compare results for gridtype 'B' with NMM model
    res = 0
    do i = 1, npts
        do j = 1, npts
            if (abs(DRAGCO(i,j) - EXP_DRAGCO_B_NMM(i,j)) > tol) then
                print *, "Mismatch at (", i, ",", j, "): ", DRAGCO(i,j), " vs ", EXP_DRAGCO_B_NMM(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 40

    print *, "Testing with invalid grid type. Expect DRAGCO = spval."
    gridtype = 'X'
    call CALDRG(DRAGCO)

    ! Compare results for invalid grid type
    res = 0
    do i = 1, npts
        do j = 1, npts
            if (DRAGCO(i,j) /= EXP_DRAGCO_INVALID(i,j)) then
                print *, "Mismatch at (", i, ",", j, "): ", DRAGCO(i,j), " vs ", EXP_DRAGCO_INVALID(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 50
    
    print *, "SUCCESS!"
end program test_caldrg
