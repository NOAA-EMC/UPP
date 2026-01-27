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
    real, dimension(1:npts,1:npts,1:nlevs) :: TWET
    real, dimension(1:npts,1:npts) :: ZWET, EXP_ZWET

    ! Grid parameters
    lm =  nlevs
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
    allocate(zint(ista_2l:iend_2u, jsta_2l:jend_2u, nlevs+1))
    allocate(pint(ista_2l:iend_2u, jsta_2l:jend_2u, nlevs+1))
    allocate(t(ista_2l:iend_2u, jsta_2l:jend_2u, nlevs))

    ! Initialize default inputs
    sm   = 0.5          ! mixed land/sea fraction
    fis  = 981.0        ! surface geopotential (m^2 s^-2) -> ~100 m MSL
    thz0 = 290.0        ! potential temperature near surface (K)
    ths  = 288.0        ! surface skin potential temperature (K)
    lmh  = real(nlevs)

    do i = ista_2l, iend_2u
        do j = jsta_2l, jend_2u
            ! Set up physically realistic vertical profiles
            ! Define surface and top heights, and compute uniform layer spacing
            z_sfc = 100.0     ! meters MSL
            z_top = 16000.0   ! meters MSL
            dz    = (z_top - z_sfc)/real(nlevs)
            p0    = 100000.0  ! Pa
            H     = 8000.0    ! scale height (m)
            t_sfc = 288.0     ! near-surface air temperature (K)
            do k= 1, nlevs+1
                zint(i,j,k) = z_top - dz*real(k-1)
                pint(i,j,k) = p0*exp(-zint(i,j,k)/H)
            end do
            do k= 1, nlevs
                z_mid       = 0.5*(zint(i,j,k) + zint(i,j,k+1))
                t(i,j,k)    = t_sfc - 6.5*((z_mid - z_sfc)/1000.0)   ! lapse rate 6.5 K/km
                TWET(i,j,k) = t(i,j,k) - 1.0                          ! wet-bulb slightly below T
            end do
        end do
    end do

    EXP_ZWET = 2.2307702637E+03  ! Expected output for the default test case
    
    ! Test Case: FIS = spval, expect ZWET = spval
    FIS(1,1) = spval
    EXP_ZWET(1,1) = spval

    ! Test Case:  tsfc < tfrz (tfrz = 273.15 K)
    thz0(1,2) = 270.0
    ths(1,2)  = 270.0
    EXP_ZWET(1,2) = -5.3279840088E+02

    ! Test Case: TWET = tfrz at top level.
    ! T = TSFC at (i,j) = (2,1)
    do k = 25, nlevs
        TWET(1,3,k) = 273.15
        TWET(2,1,k) = 273.15
    end do
    t_sfc = sm(2,1) * thz0(2,1) + (1.0 - sm(2,1)) * ths(2,1)  &
                    * (pint(2,1,nlevs+1)/p1000)**capa

    t(2,1,nlevs) = t_sfc


    ! TODO: Replace ??? with code to set up a test case at (i,j) = (3,1), with the appropriate 
    ! vertical profile, such that 
    ! TWET(3,1,nlev) = 273.15 K and ZWET(3,1) < ZU and -ZWET(3,1) < ZU where:
    ! ZU = 0.5*(ZINT(3,1,nlevs)+ZINT(3,1,nlevs+1))
    sm(3,1)    = 0.5
    thz0(3,1)  = 270.0
    ths(3,1)   = 270.0
    TWET(3,1,nlevs) = 273.15
    
    call WETFRZLVL(TWET, ZWET)

    do i = ista, iend
        do j = jsta, jend
            print '(A,I0,A,I0,A,ES24.10)', "Point (", i, ",", j, "): WETFRZLVL = ", ZWET(i,j)
        end do
    end do

    print *, "SUCCESS!"
end program test_wetfrzlvl