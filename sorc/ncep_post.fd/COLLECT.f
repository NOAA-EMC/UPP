!> @file
!> @brief Subroutine that collect gathers from all MPI tasks.
!>
!> @param[in] A Array being gathered.
!> @param[out] B gathered array - only valid on task 0.
!>
!> Gather "A" from all MPI tasks onto task 0.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2000-01-06 | Jim Tuccillo | Initial
!>
!> @author Jim Tuccillo IBM @date 2000-01-06
!-----------------------------------------------------------------------
!> COLLECT() Subroutine that collect gathers from all MPI tasks.
!>
!> @param[in] A Array being gathered.
!> @param[out] B gathered array - only valid on task 0.
!-----------------------------------------------------------------------
      SUBROUTINE COLLECT (A, B) 


      use ctlblk_mod, only: num_procs, jsta, icnt, idsp, mpi_comm_comp, im, jm, me
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      implicit none
!
      include 'mpif.h'
!
      real,intent(in)  :: a ( im, jm ) 
      real,intent(out) :: b ( im, jm ) 
!     integer i, j
      integer ierr
!
      if ( num_procs <= 1 ) then
         b = a
      else
         call mpi_gatherv(a(1,jsta),icnt(me),MPI_REAL,              &
     &                    b,icnt,idsp,MPI_REAL,0,MPI_COMM_COMP,ierr)
      end if
      end               
