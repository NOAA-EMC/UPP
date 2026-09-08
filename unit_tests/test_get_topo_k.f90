! This is a test program for UPP.
!
! This program tests the getTopoK() function from GFIP3.f
!
! Alyson Stahl, 9/2026
program test_get_topo_k
    implicit none
    integer, parameter :: nz = 10
    !
    integer :: i, res
    integer :: topo_k, exp_topo_k
    real :: alt, hgt(nz)

    interface
        function getTopoK(hgt, alt, nz)
            real, intent(in) :: hgt(nz)
            real, intent(in) :: alt
            integer, intent(in) :: nz
            integer :: getTopoK
        end function getTopoK
    end interface

    ! Test Case 1: Typical case where hgt(nz) < alt
    alt = 800.0
    hgt = (/ 12000.0, 10000.0, 8500.0, 7000.0, 5600.0, 4300.0, 3100.0, 2000.0, 1200.0, 100.0 /)

    exp_topo_k = 9
    res = 0

    topo_k = getTopoK(hgt, alt, nz)

    if (topo_k .ne. exp_topo_k) stop 10

    ! Test Case 2: hgt(nz) >= alt (expect topo_k = nz)
    alt = 50.0
    hgt = (/ 12000.0, 10000.0, 8500.0, 7000.0, 5600.0, 4300.0, 3100.0, 2000.0, 1200.0, 100.0 /)

    exp_topo_k = nz

    topo_k = getTopoK(hgt, alt, nz)

    if (topo_k .ne. exp_topo_k) stop 20

    print *, "Success!"

end program test_get_topo_k