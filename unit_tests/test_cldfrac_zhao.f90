! This is a test program for UPP.
!
! This program tests the progcld1() subroutine in CLDFRAC_ZHAO.f
!
! Alyson Stahl, 1/2026
program test_cldfrac_zhao
    use kinds, only: r_kind
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: nx = 3, ny = 3
    integer :: i, j, res
    integer :: IX, NLAY, IFLIP
    real(r_kind), dimension(1:nx,1:ny) :: PLYR, TLYR, QLYR, QSTL, CLW
    real(r_kind), dimension(1:nx,1:ny) :: CLDTOT, EXP_CLDTOT

    interface
        subroutine PROGCLD1(PLYR, TLYR, QLYR, QSTL, CLW, IX, NLAY, IFLIP, CLDTOT)
            use kinds, only: r_kind
            integer, intent(in) :: IX, NLAY, IFLIP
            real(kind=r_kind), dimension(IX, NLAY), intent(in) :: PLYR, TLYR, &
                QLYR, QSTL, CLW
            real(kind=r_kind), dimension(IX, NLAY), intent(out) :: CLDTOT
        end subroutine PROGCLD1
    end interface
    
    IX = nx
    NLAY = ny

    PLYR(1:nx,1) = 400.0
    PLYR(1:nx,2) = 650.0
    PLYR(1:nx,3) = 900.0

    TLYR(1:nx,1) = 240.0
    TLYR(1:nx,2) = 260.0
    TLYR(1:nx,3) = 280.0

    QSTL = 0.0100
    QLYR = 0.0090
    CLW  = 2.0e-6

    EXP_CLDTOT = 2.166432329E-02

    ! Test Case: will result in CLDTOT = 0.0
    CLW(1,1) = 1.0e-7
    EXP_CLDTOT(1,1) = 0.0

    ! Test Case: values set such that tem1 evaluates to 2000.0. Will result in CLDTOT = 0.0
    QSTL(1,2) = 1.0
    QLYR(1,2) = 0.0  
    EXP_CLDTOT(1,2) = 0.0

    ! Test Case: values set such that value evaluates to 50.0.
    QLYR(2,1) = QSTL(2,1)
    CLW(2,1)  = 3.0e-5
    EXP_CLDTOT(2,1) = 1.0

    ! Test Case: values set such that tem1 evaluates to 2000/0.001
    QLYR(3,2) = QSTL(3,2)  
    EXP_CLDTOT(3,2) = 9.816843271E-01
    
    ! Test Case: IFLIP = 0 (input data from toa to sfc)
    IFLIP = 0
    call PROGCLD1(PLYR, TLYR, QLYR, QSTL, CLW, IX, NLAY, IFLIP, CLDTOT)

    res = 0
    do i = 1, nx
        do j = 1, ny
            if (abs(CLDTOT(i,j) - EXP_CLDTOT(i,j)) > tol) then
                print *, 'CLDTOT Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_CLDTOT(i,j), &
                         ' but got ', CLDTOT(i,j)
                res = 1
            end if
        end do
    end do
    
    if (res .ne. 0) stop 10

    ! Test Case: IFLIP = 1 (input data from sfc to toa)
    IFLIP = 1
    call PROGCLD1(PLYR, TLYR, QLYR, QSTL, CLW, IX, NLAY, IFLIP, CLDTOT)
    
    do i = 1, nx
        do j = 1, ny
            if (abs(CLDTOT(i,j) - EXP_CLDTOT(i,j)) > tol) then
                print *, 'CLDTOT Test failed at (', i, ',', j, '): ', &
                         'Expected ', EXP_CLDTOT(i,j), &
                         ' but got ', CLDTOT(i,j)
                res = 1
            end if
        end do
    end do
    
    if (res .ne. 0) stop 20

    print *, "SUCCESS!"
end program test_cldfrac_zhao