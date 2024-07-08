!> @file
!> @brief Subroutines in this file set up decomposition values for 1D and 2D decomposition.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2000-01-06 | Jim Tuccillo | Initial
!>
!> @author Jim Tuccillo IBM @date 2000-01-06
!-----------------------------------------------------------------
!> @brief Sets up decomposition values.
!>
!> @param[in] N1 First interate value.
!> @param[in] N2 Last interate value.
!> @param[in] NPROCS Number of MPI tasks.
!> @param[in] IRANK My taks ID.
!> @param[out] ISTA First loop value.
!> @param[out] IEND Last loop value.
!>
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

! ----------------------------------------------------------------------------------------------------
!> @brief para_range2() sets up 2D decomposition values
!> @param[in] N1 - LAAT INTERATE VALUE I dimension
!> @param[in] N2 - LAST INTERATE VALUE  J dimension
!> @param[in] NX NUMBER OF subdomains in Z dimension (NX * NY should be the total number of MPI procs)
!> @param[in] NY NUMBER OF subdomains in Y dimension (NX * NY should be the total number of MPI procs)
!> @param[in] NRANK - MY TASK ID
!> @param[out] ISTA - FIRST LOOP VALUE I
!> @param[out] IEND - LAST LOOP VALUE I
!> @param[out] JSTA - FIRST LOOP VALUE J
!> @param[out] JEND - LAST LOOP VALUE J
     subroutine para_range2(im,jm,nx,ny,nrank,ista,iend,jsta,jend)

      implicit none
      integer,intent(in)  ::  im,jm,nx,ny,nrank
      integer,intent(out) ::  ista,iend,jsta,jend
      integer             ::  ix,jx

         jx=nrank/nx
         ix=nrank-(jx*nx)
           call para_range(1,im,nx,ix,ista,iend)
           call para_range(1,jm,ny,jx,jsta,jend)
!            print 101,n,ix,jx,ista,iend,jsta,jend
! 101   format(16i8)
          return
          end


