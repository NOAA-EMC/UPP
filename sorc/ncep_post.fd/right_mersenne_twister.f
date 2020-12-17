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
  end module
