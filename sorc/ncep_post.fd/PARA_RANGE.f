!> @file
!> @brief para_range() sets up decomposition values.
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
!!
!! USAGE:    CALL PARA_RANGE2(N1,N2,NX,NY,NRANK,ISTA,IEND,JSTA,JEND)(A)
!!   INPUT ARGUMENT LIST:
!!     N1 - LAAT INTERATE VALUE I dimension
!!     N2 - LAST INTERATE VALUE  J dimension
!!     NX  NUMBER OF subdomains in Z  dimension
!!     NY  NUMBER OF subdomains  in Y dimension
!!        NX * NY should be the total number of MPI procs
!!     NRANK - MY TAKS ID
!!
!!   OUTPUT ARGUMENT LIST:
!!     ISTA - FIRST LOOP VALUE I
!!     IEND - LAST LOOP VALUE I
!!     JSTA - FIRST LOOP VALUE J
!!     JEND - LAST LOOP VALUE J
!!
!!   OUTPUT FILES:
!!     STDOUT  - RUN TIME STANDARD OUT.
!!
!!   SUBPROGRAMS CALLED:
!!     UTILITIES:
!!       NONE
!!     LIBRARY:
!!
!!   ATTRIBUTES:
!!     LANGUAGE: FORTRAN
!!     MACHINE : IBM RS/6000 SP
!!
     subroutine para_range2(im,jm,nx,ny,nrank,ista,iend,jsta,jend)
         jx=nrank/nx
         ix=nrank-(jx*nx)
           call para_range(1,im,nx,ix,ista,iend)
           call para_range(1,jm,ny,jx,jsta,jend)
            print 101,n,ix,jx,ista,iend,jsta,jend
 101   format(16i8)
          return
          end


