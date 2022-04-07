!> @file
!> @brief Subroutine that collect gathers from all MPI tasks.
!>
!> @param[in] A Array being gathered.
!> @param[out] A gathered array - only valid on task 0.
!>
!> Gather "A" from all MPI tasks onto task 0.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2000-01-06 | Jim Tuccillo | Initial
!>
!> @author Jim Tuccillo IBM @date 2000-01-06
      SUBROUTINE COLLECT_LOC ( A, B ) 


      use CTLBLK_mod, only: num_procs, jsta, icnt, idsp, mpi_comm_comp, im,&
              jsta_2l, jend_2u, jm, me
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      include 'mpif.h'
      real, dimension(im,jsta_2l:jend_2u), intent(in) :: a
      real, dimension(im,jm), intent(out) :: b
      integer ierr
!
      if ( num_procs <= 1 ) then
         b = a
      else
         call mpi_gatherv(a(1,jsta),icnt(me),MPI_REAL,   &
     &    b,icnt,idsp,MPI_REAL,0,MPI_COMM_COMP, ierr )
      
      end if

      end               
