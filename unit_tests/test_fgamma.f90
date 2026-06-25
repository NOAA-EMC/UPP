! This is a test program for UPP.
!
! This program tests the FGAMMA() subroutine.
!
! Alyson Stahl, 1/2026
program test_fgamma
    implicit none

    real, parameter :: tol = 1.0e-6
    real, parameter :: XBIG = 35.040E0, XMININ = 1.18E-38, EPS = 1.19E-7, &
                        XINF = 3.4E38
    integer, parameter :: ntests = 12
    integer :: i, res
    real, dimension(1:ntests) :: X, FX, EXP_FX

    interface
        real function fGAMMA(X)
            real, intent(in) :: X
        end function fGAMMA
    end interface

    ! Initialize input array X with test values
    X(1) = 0.0 ! Y = 0 with RES = 0
    X(2) = -1.0 ! Y < 0 with RES = 0
    X(3) = -1.2 ! Y <= 0 with RES != 0 and PARITY = .TRUE.
    X(4) = -0.5 ! Y <= 0 with RES != 0 and PARITY = .FALSE.
    X(5) = 5.9e-39 ! 0 < Y < EPS and Y < XMININ
    X(6) = 1.18e-38 ! 0 < Y < EPS and Y = XMININ
    X(7) = 0.5 ! In range [EPS, 12) with 0 < Y < 1
    X(8) = 1.5 ! In range [EPS, 12) with 1 < Y < 2
    X(9) = 5.0 ! In range [EPS, 12) with 2 < Y < 12
    X(10) = 12.0 ! 12 <= Y <= XBIG
    X(11) = 35.04 ! Y = XBIG
    X(12) = 35.1 ! Y > XBIG

    EXP_FX = (/ XINF, XINF, 4.8509559631E+00, &
                -3.5449078083E+00, XINF, 8.4745768940E+37, &
                1.7724539042E+00, 8.8622695208E-01, 2.4000000000E+01, &
                3.9916800000E+07, 3.4016304551E+38, XINF/)
                
    do i = 1, ntests
        FX(i) = fGAMMA(X(i))
    end do

    res = 0
    do i = 1, ntests
        if (abs(FX(i) - EXP_FX(i)) > tol) then
            print *, "Test ", i, " failed: FGAMMA(", X(i), ") = ", FX(i), &
                     ", expected ", EXP_FX(i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10 

    print *, "SUCCESS!"
end program test_fgamma