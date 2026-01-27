! This is a test program for UPP.
!
! This program tests the ZENSUN() subroutine.
!
! Alyson Stahl, 1/2026
program test_zensun
    use kinds, only: r_kind,i_kind
    implicit none
    
    real, parameter :: tol = 1.0e-8
    integer, parameter :: ntests = 5
    integer :: i, res
    ! Input
    real(r_kind), parameter :: PI = 3.14159265358979323846
    integer(i_kind), dimension(ntests) :: DAY
    real(r_kind), dimension(ntests) :: TIME, LAT, LON
    ! Output
    real(r_kind), dimension(ntests) :: SUN_ZENITH, SUN_AZIMUTH

    print *, "SUCCESS!"
end program test_zensun