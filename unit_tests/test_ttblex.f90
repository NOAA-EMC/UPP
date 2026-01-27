! This is a test program for UPP.
!
! This program tests the TTBLEX() subroutine.
!
! Alyson Stahl, 1/2026
program test_ttblex
    use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u,  &
                        ista, iend, ista_2l, iend_2u
    implicit none
    
    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 3, ni = 50, nj = 50
    integer :: i, j, res
    ! Inputs
    integer :: ITB, JTB, KARR(1:npts,1:npts)
    real, dimension(1:ni,1:nj) :: TTBL
    real, dimension(1:npts,1:npts) :: PMIDL, THESP
    real, dimension(1:ni) :: THE0, STHE
    real :: PL, RDP, RDTHE
    ! Outputs
    real, dimension(1:npts,1:npts) :: TREF, QQ, PP, EXP_TREF, EXP_QQ, EXP_PP
    integer, dimension(1:npts,1:npts) :: IPTB, ITHTB, EXP_IPTB, EXP_ITHTB

    ! Grid parameters
    jsta = 1
    jend = npts
    jsta_2l = jsta
    jend_2u = jend
    ista = 1
    iend = npts
    ista_2l = ista
    iend_2u = iend
    
    ITB = ni
    JTB = nj
    KARR = 1
    
    ! Initialize temperature table (TTBL) in Kelvin
    do j = 1, nj
        do i = 1, ni
            TTBL(i,j) = 200.0 + 0.6*j + 0.3*i
        end do
    end do

    ! Initialize mid-layer pressure (PMIDL) in Pascals
    do j = 1, npts
        do i = 1, npts
            PMIDL(i,j) = 60000.0 + 5000.0*(i-2) + 3000.0*(j-2)
        end do
    end do

    ! Initialize saturation potential temperature (THESP) in Kelvin
    do j = 1, npts
        do i = 1, npts
            THESP(i,j) = 305.0 + 0.5*(i-2) + 0.3*(j-2)
        end do
    end do

    ! Initialize theta table base and scale (K)
    do i = 1, ni
        THE0(i) = 300.0 + 0.2*(i-1)
        STHE(i) = 40.0
    end do

    ! Pressure table base and scaling
    PL = 10000.0      ! Pa
    RDP = 5.0e-4      ! 1/Pa

    ! Theta table reciprocal scaling
    RDTHE = 5.0e-2    ! 1/K

    call TTBLEX(TREF, TTBL, ITB, JTB, KARR, PMIDL, PL, QQ, PP, RDP, THE0, &
                STHE, RDTHE, THESP, IPTB, ITHTB)

    print *, "SUCCESS!"
end program test_ttblex