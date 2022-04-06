!> @file
!> @brief retrieve_index gets record number of desired variable.
!>
!> By examining previously generated inventory of wrf binary restart file,
!> find record number that contains the header record for variable
!> identified by input character variable "string".
!>
!> @param[in] string Mnemonic for variable desired.
!> @param[in] varname_all List of all mnemonics obtained from inventory of file.
!> @param[in] nrecs Total number of sequential records counted in wrf binary restart file.
!> @param[out] index Desired record number.
!> @param[out] iret Return status, set to 0 if variable was found, non-zero if not.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2004-11-29 | Parrish | Initial
!>
!> @author Parrish np22 @date 2004-11-29
      subroutine retrieve_index(index,string,varname_all,nrecs,iret)


  implicit none

  integer,intent(out)::iret
  integer,intent(in)::nrecs
  integer,intent(out):: index
  character(*), intent(in):: string
  character(132),intent(in)::varname_all(nrecs)

  integer i

  iret=0

  do i=1,nrecs
   if(trim(string) == trim(varname_all(i))) then
      index=i
      return
   end if
  end do

  write(6,*)' problem reading wrf nmm binary file, rec id "',trim(string),'" not found'

  iret=-1

end subroutine retrieve_index
