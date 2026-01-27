! This is a test program for UPP.
!
! This program tests the WETFRZLVL() subroutine.
!
! Alyson Stahl, 1/2026
program test_wetfrzlvl
    use vrbls3d, only: pint, zint, t
    use vrbls2d, only:  fis, thz0, ths
    use masks, only: lmh, sm
    use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u, lm, spval, &
                        ista, iend, ista_2l, iend_2u 
    implicit none
    
    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 3, nlevs = 4
    integer :: i, j, k, res
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

    sm   = 0.0
    fis  = 0.0
    thz0 = 285.0
    ths  = 285.0
    lmh  = real(nlevs)

    do i = ista_2l, iend_2u
        do j = jsta_2l, jend_2u
            do k= 1, nlevs+1
                zint(i,j,k) = (k-1)*1000.0      ! interface heights [m]
                pint(i,j,k) = 100000.0 - 20000.0*(k-1)  ! interface pressures [Pa]
            end do
            do k= 1, nlevs
                t(i,j,k)    = 285.0 - 5.0*(k-1) ! layer temperatures [K]
                TWET(i,j,k) = t(i,j,k) - 2.0     ! wet-bulb approx [K]
            end do
        end do
    end do

    print *, "SUCCESS!"
end program test_wetfrzlvl