!> @file
!> @brief PARA_RANGE sets up decomposition values.
!>
!> This subroutine sets up decomposition values.
!>
!> @param[in] N1 First interate value.
!> @param[in] N2 Last interate value.
!> @param[in] NPROCS Number of MPI tasks.
!> @param[in] IRANK My taks ID.
!> @param[out] ISTA First loop value.
!> @param[out] IEND Last loop value.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2000-01-06 | Jim Tuccillo | Initial
!>
!> @author Jim Tuccillo IBM @date 2000-01-06
      SUBROUTINE PARA_RANGE (N1,N2,NPROCS,IRANK,ISTA,IEND)

      implicit none
      integer,intent(in)  ::  n1,n2,nprocs,irank
      integer,intent(out) ::  ista,iend
      integer iwork1, iwork2

      iwork1 = ( n2 - n1 + 1 ) / nprocs
      iwork2 = mod ( n2 - n1 + 1, nprocs )
      ista   = irank * iwork1 + n1 + min ( irank, iwork2 )
      iend   = ista + iwork1 - 1
      if ( iwork2 > irank ) iend = iend + 1
      return
      end

