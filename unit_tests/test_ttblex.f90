! This is a test program for UPP.
!
! This program tests the TTBLEX() subroutine.
!
! Alyson Stahl, 1/2026
program test_ttblex
    use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u,  &
                        ista, iend, ista_2l, iend_2u
    implicit none
    
    ! spval for tests
    real, parameter :: tol = 1.0e-8, spval = 9.9e10
    integer, parameter :: npts = 2, ni = 50, nj = 50
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

    interface
        subroutine TTBLEX(TREF, TTBL, ITB, JTB, KARR, PMIDL,            &
                          PL, QQ, PP, RDP, THE0, STHE, RDTHE, THESP,    &
                          IPTB, ITHTB)
            use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u, &
                                ista, iend, ista_2l, iend_2u
            integer, intent(in) :: ITB,JTB
            integer, intent(in) :: KARR(ista:iend,jsta:jend)
            real, dimension(JTB,ITB), intent(in) :: TTBL
            real, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in) :: PMIDL
            real, dimension(ista:iend,jsta:jend), intent(in) :: THESP
            real, dimension(ITB),  intent(in) :: THE0,STHE
            real, intent(in) :: PL,RDP,RDTHE
            integer, dimension(ista:iend,jsta:jend), intent(out) :: IPTB,ITHTB
            real, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(out) :: TREF
            real, dimension(ista:iend,jsta:jend), intent(out) :: QQ,PP
        end subroutine TTBLEX
    end interface

    ! Grid parameters
    jsta = 1
    jend = npts
    jsta_2l = jsta
    jend_2u = jend
    ista = 1
    iend = npts
    ista_2l = ista
    iend_2u = iend
    
    ! Initialize inputs
    ITB = ni
    JTB = nj
    KARR = 1
    PL = 0.0
    RDP = 5.0e-4
    RDTHE = 1.0

    do i = 1, ni
        do j = 1, nj
            TTBL(i,j) = 220.0 + 0.6*i + 0.4*j
        end do
    end do

    do i = 1, npts
        do j = 1, npts
            PMIDL(i,j) = 55000.0 + 137.0*i + 53.0*j
            THESP(i,j) = 310.0 + 0.7*j + 0.3*i
        end do
    end do
    
    do i = 1, ni
        THE0(i) = 300.0 + 0.8*(i-1)
        STHE(i) = 12.5
    end do

    ! Expected Outputs
    EXP_TREF = reshape([2.2080192566E+02, spval, 2.6900000000E+02, 2.3158871460E+02], [npts, npts])
    EXP_QQ   = reshape([-5.0000000745E-02, spval, 0.0000000000E+00, 6.9000053406E-01], [npts, npts])
    EXP_PP   = reshape([-2.9679930210E-01, spval, 0.0000000000E+00, -8.1216067076E-01], [npts, npts])
    EXP_IPTB = reshape([1, 0, 49, 28], [npts, npts])
    EXP_ITHTB= reshape([1, 0, 49, 1], [npts, npts])

    ! Test Case: IPTB and ITHTB clipped to 1
    PMIDL(1,1) = PL - 100.0
    THESP(1,1) = THE0(1) - 0.3*STHE(1)

    ! Test Case: IPTB and ITHTB clipped to ITB-1 and JTB-1, respectively
    PMIDL(1,2) = PL + (ITB + 10)/RDP
    THESP(1,2) = THE0(ITB-1) + STHE(ITB-1)*(JTB + 0.2)

    ! Test Case: KARR = 0. This won't compute anything and will leave output arrays as is.
    ! Output arrays set to default value before call for predictable output. 
    ! (spval for reals, 0 for integers)
    KARR(2,1) = 0
    TREF = spval
    QQ   = spval
    PP   = spval
    IPTB = 0
    ITHTB= 0

    call TTBLEX(TREF, TTBL, ITB, JTB, KARR, PMIDL, PL, QQ, PP, RDP, THE0, &
                STHE, RDTHE, THESP, IPTB, ITHTB)

    res = 0
    do i = 1, npts
        do j = 1, npts
            if ( abs(TREF(i,j) - EXP_TREF(i,j)) > tol ) then
                print *, "Test failed at TREF (i,j)=(", i, ",", j, "): ", &
                         "Expected TREF = ", EXP_TREF(i,j), &
                         ", Computed TREF = ", TREF(i,j)
                res = 1
            end if
            if ( abs(QQ(i,j) - EXP_QQ(i,j)) > tol ) then
                print *, "Test failed at QQ (i,j)=(", i, ",", j, "): ", &
                         "Expected QQ = ", EXP_QQ(i,j), &
                         ", Computed QQ = ", QQ(i,j)
                res = 1
            end if
            if ( abs(PP(i,j) - EXP_PP(i,j)) > tol ) then
                print *, "Test failed at PP (i,j)=(", i, ",", j, "): ", &
                         "Expected PP = ", EXP_PP(i,j), &
                         ", Computed PP = ", PP(i,j)
                res = 1
            end if
            if ( IPTB(i,j) /= EXP_IPTB(i,j) ) then
                print *, "Test failed at IPTB (i,j)=(", i, ",", j, "): ", &
                         "Expected IPTB = ", EXP_IPTB(i,j), &
                         ", Computed IPTB = ", IPTB(i,j)
                res = 1
            end if
            if ( ITHTB(i,j) /= EXP_ITHTB(i,j) ) then
                print *, "Test failed at ITHTB (i,j)=(", i, ",", j, "): ", &
                         "Expected ITHTB = ", EXP_ITHTB(i,j), &
                         ", Computed ITHTB = ", ITHTB(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 10
    
    print *, "SUCCESS!"
end program test_ttblex