!> @file
!>
!> @brief This module, native_endianness, was written by Dusan Jovic and has been adapted to GSI for internal translation
!> of WRF ARW and NMM binary restart files as required to match the machine native 
!> endian storage format. The original code only converted from big-endian to little-endian.
!> There are no restrictions in this version.
!> This is required for these two types of files, because they are read/written to using mpi-io,
!> which has no compiler option for automatic switching to machine native endian format
!> for fortran unformatted read/write.
!>
!> @author Parrish wx22 @date 2012-10-11
!> @note functions included: is_little_endian - no argument--returns true for little-endian machine, false for big-endian machine
!>
!> @note variables included: byte_swap - false if machine and wrf binary file are same endian, true if different
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2012-10-11 | Parrish | Initial. Copy/modify original module native_endianness provided by Dusan Jovic, NCEP/EMC 2012
!> 2012-10-19 | parrish | Additional modifications to improve efficiency.  Remove interface and make to_native_endianness to work only with integer(4) arguments. Put to_native_endianness_i4 outside module.
!>
!> @author Parrish wx22 @date 2012-10-11
      module native_endianness


 use kinds, only: i_byte,i_long
 implicit none

 private

 public byte_swap
 public is_little_endian

 !> false if machine and wrf binary file are same endian, true if different
 logical byte_swap

 contains

 logical function is_little_endian()
!> is_little_endian() tests to see if machine is little-endian. Returns true for little-endian, false for big-endian.
!> @return boolean value (true or false)
!> 
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2012-10-11 | Parrish | Add doc block
!>
!> @author Parrish wx22 @date 2012-10-11
   implicit none

   integer(i_byte) :: i1
   integer(i_long) :: i2

   i1 = 1
   i2 = 0
   i2 = transfer(i1,i2)

   is_little_endian = (i1 == i2)

 end function is_little_endian

 end module native_endianness

!----------------------------------------------------------------------
! convert 4-byte integer scalar from big-endian to native-endian
!----------------------------------------------------------------------

 subroutine to_native_endianness_i4(i4,num)
!> to_native_endianness_i4() is to swap bytes of argument.
!>
!> @param[in] i4 Input 4 byte integer array.
!> @param[in] num Length of array i4.  (NOTE:  type of num must be i_llong (8 byte integer) )
!> @param[out] i4 Output 4 byte integer array with bytes in reverse order.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2012-10-11 | Parrish | Add doc block
!> 2012-10-19 | Parrish | Additional modifications to improve efficiency.  Remove interface and make to_native_endianness to work only with integer(4) arguments. Put to_native_endianness_i4 outside module.
!>
!> @author Parrish wx22 @date 2012-10-11
 use kinds, only: i_byte,i_long,i_llong
 implicit none

 integer(i_llong), intent(in) :: num
 integer(i_long), intent(inout) :: i4(num)

 integer(i_byte), dimension(4) :: byte_arr, byte_arr_tmp
 integer(i_long) :: i,n

 do n=1,num
    byte_arr_tmp = transfer (i4(n), byte_arr)
    byte_arr(1)=byte_arr_tmp(4)
    byte_arr(2)=byte_arr_tmp(3)
    byte_arr(3)=byte_arr_tmp(2)
    byte_arr(4)=byte_arr_tmp(1)
    i4(n) = transfer (byte_arr, i4(n))
 end do

 return

 end subroutine to_native_endianness_i4
