!> @file
!
!! This module is set up to ensure the random_number call from mersenne_twister
!! module of W3EMC returns right values (0.-1.0) whenever mersenne_twister is
!! linked with real 4 or real 8. 
!! PROGRAM HISTORY LOG:
!! - 20-12-19 Wen Meng - Initial Code
!!
  module upp_right_mersenne_twister
    interface right_random_number
      module procedure random_number_4
      module procedure random_number_8
    end interface
    contains
    subroutine random_number_4(rn)
      implicit none
      real(kind=4) :: rn(:)
      call random_number(rn)
    end subroutine random_number_4
    subroutine random_number_8(rn)
      implicit none
      real(kind=8) :: rn(:)
      call random_number(rn)
    end subroutine random_number_8
  end module upp_right_mersenne_twister
