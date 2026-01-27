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
    PL = 10000.0      ! Pa
    RDP = 5.0e-4      ! 1/Pa
    RDTHE = 5.0e-2    ! 1/K

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

    ! TODO: For one (i,j) point in (1:npts,1:npts), modify the input arrays
    ! so that IPTB and ITHTB will fall below 1. Replace the ??? with appropriate code.
    PMIDL(1,1) = PL - 3000.0      ! force TPK <= -1 -> IPTB < 1 before clamp
    THESP(1,1) = THE0(1) - 5.0    ! ensure THESP < BTHK for negative TTHK
    STHE(1)    = 0.1              ! small scale to make TTHK <= -1 -> ITHTB < 1

    call TTBLEX(TREF, TTBL, ITB, JTB, KARR, PMIDL, PL, QQ, PP, RDP, THE0, &
                STHE, RDTHE, THESP, IPTB, ITHTB)

    do i = 1, npts
        do j = 1, npts
            print '(A,I0,A,I0,A,ES24.10)', "TREF(", i, ",", j, ") = ", TREF(i,j)
            print '(A,I0,A,I0,A,ES24.10)', "QQ(", i, ",", j, ")   = ", QQ(i,j)
            print '(A,I0,A,I0,A,ES24.10)', "PP(", i, ",", j, ")   = ", PP(i,j)
            print '(A,I0,A,I0,A,I0)',      "IPTB(", i, ",", j, ") = ", IPTB(i,j)
            print '(A,I0,A,I0,A,I0)',      "ITHTB(", i, ",", j, ")= ", ITHTB(i,j)
        end do
    end do
    print *, "SUCCESS!"
end program test_ttblex