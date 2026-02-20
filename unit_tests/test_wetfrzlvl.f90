! This is a test program for UPP.
!
! This program tests the WETFRZLVL() subroutine.
!
! Alyson Stahl, 1/2026
program test_wetfrzlvl
    use vrbls3d, only: pint, zint, t
    use vrbls2d, only:  fis, thz0, ths
    use masks, only: lmh, sm
    use params_mod, only: p1000, capa
    use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u, lm, spval, &
                        ista, iend, ista_2l, iend_2u 
    implicit none
    
    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 3, nlevs = 30
    integer :: i, j, k, res
    real :: z_sfc, z_top, dz, p0, H, z_mid, t_sfc
    real, dimension(1:npts,1:npts,1:nlevs+1) :: TWET
    real, dimension(1:npts,1:npts) :: ZWET, EXP_ZWET

    interface
        subroutine WETFRZLVL(TWET,ZWET)
            use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u, lm, &
                                ista, iend, ista_2l, iend_2u 
            real, intent(in) :: TWET(ista_2l:iend_2u,jsta_2l:jend_2u,lm)
            real, intent(out) :: ZWET(ista:iend,jsta:jend)
        end subroutine WETFRZLVL
    end interface

    ! Grid parameters
    lm =  nlevs+1
    jsta = 1
    jend = npts
    jsta_2l = jsta
    jend_2u = jend
    ista = 1
    iend = npts
    ista_2l = ista
    iend_2u = iend
    spval = 9.9e10

    allocate(sm(ista_2l:iend_2u, jsta_2l:jend_2u))
    allocate(fis(ista_2l:iend_2u, jsta_2l:jend_2u))
    allocate(thz0(ista_2l:iend_2u, jsta_2l:jend_2u))
    allocate(ths(ista_2l:iend_2u, jsta_2l:jend_2u))
    allocate(lmh(ista_2l:iend_2u, jsta_2l:jend_2u))
    allocate(zint(ista_2l:iend_2u, jsta_2l:jend_2u, nlevs+2))
    allocate(pint(ista_2l:iend_2u, jsta_2l:jend_2u, nlevs+2))
    allocate(t(ista_2l:iend_2u, jsta_2l:jend_2u, nlevs+1))

    ! Initialize default inputs
    z_sfc = 200.0
    z_top = 5200.0
    dz    = (z_top - z_sfc)/real(nlevs+2)
    p0    = 100000.0
    H     = 7500.0

    sm   = 0.5
    fis  = 9.81*z_sfc
    thz0 = 290.0
    ths  = 288.0
    lmh  = real(nlevs)

    do i = ista_2l, iend_2u
        do j = jsta_2l, jend_2u
            do k = 1, nlevs+2
                zint(i,j,k) = z_top - real(k-1)*dz
                pint(i,j,k) = p0 * exp( - zint(i,j,k) / H )
            end do
            do k = 1, nlevs+1
                t(i,j,k)    = 240.0 + (260.0 - 240.0) * real(k-1) / real(nlevs-1)
                TWET(i,j,k) = 260.0 + (273.15 - 260.0) * real(k-1) / real(nlevs-1)
            end do
        end do
    end do

    ! Expected output for default test case
    EXP_ZWET = 3.9467242432E+02

    ! Test Case: FIS = spval
    FIS(1,1) = spval
    EXP_ZWET(1,1) = spval

    ! Test Case: Freezing level at ground level
    thz0(1,2) = 270.0
    ths(1,2)  = 270.0
    EXP_ZWET(1,2) = -1.0882454834E+03

    ! Test Case: Freezing level above heighest model level
    thz0(1,3) = 290.0
    ths(1,3)  = 290.0
    do k = 1, nlevs+1
        TWET(1,3,k) = 275.0
    end do
    EXP_ZWET(1,3) = z_sfc

    ! Test Case: DELT = 0 branch
    sm(2,1)   = 0.4
    thz0(2,1) = 247.5
    ths(2,1)  = 300.0
    t_sfc = sm(2,1)*thz0(2,1) + (1.0 - sm(2,1))*ths(2,1) * ( pint(2,1,nlevs+1) / p1000 )**capa
    t(2,1,nlevs) = t_sfc
    EXP_ZWET(2,1) = 5.6424841309E+02

    ! Test Case: ZWET clipped to ZU if ZWET > ZU
    t_sfc = sm(2,2)*thz0(2,2) + (1.0 - sm(2,2))*ths(2,2) * ( pint(2,2,nlevs+1) / p1000 )**capa
    t(2,2,nlevs) = t_sfc - 0.5
    EXP_ZWET(2,2) = 5.9062500000E+02

    ! Test Case: ZWET clipped to ZU if -ZWET > ZU
    t_sfc = sm(2,3)*thz0(2,3) + (1.0 - sm(2,3))*ths(2,3) * ( pint(2,3,nlevs+1) / p1000 )**capa
    t(2,3,nlevs) = t_sfc + 1.2
    EXP_ZWET(2,3) = 5.9062500000E+02

    ! Test Case: TWET <= TFRZ below model top
    do k = 1, nlevs+1
        TWET(3,1,k) = 260.0 + (285.0 - 260.0) * real(k-1) / real(nlevs-1)
    end do
    EXP_ZWET(3,1) = 2.7384375000E+03

    call WETFRZLVL(TWET, ZWET)

    deallocate(sm)
    deallocate(fis)
    deallocate(thz0)
    deallocate(ths)
    deallocate(lmh)
    deallocate(zint)
    deallocate(pint)
    deallocate(t)

    res = 0
    do i = 1, npts
        do j = 1, npts
            if ( abs(ZWET(i,j) - EXP_ZWET(i,j)) > tol ) then
                print *, "Test failed at (i,j)=(", i, ",", j, "): ", &
                         "Expected ZWET = ", EXP_ZWET(i,j), &
                         ", Computed ZWET = ", ZWET(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10

    print *, "SUCCESS!"
end program test_wetfrzlvl