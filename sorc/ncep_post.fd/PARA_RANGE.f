!> @file
!
!> SUBPROGRAM:    PARA_RANGE  SET UP DECOMPOSITION VALUES
!!   PRGRMMR: TUCCILLO        ORG: IBM
!!
!! ABSTRACT:
!!     SETS UP DECOMOSITION VALUES
!!
!! PROGRAM HISTORY LOG:
!!   00-01-06  TUCCILLO - ORIGINAL
!!
!! USAGE:    CALL PARA_RANGE (N1,N2,NPROCS,IRANK,ISTA,IEND)(A)
!!   INPUT ARGUMENT LIST:
!!     N1 - FIRST INTERATE VALUE
!!     N2 - LAST INTERATE VALUE
!!     NPROCS - NUMBER OF MPI TASKS
!!     IRANK - MY TAKS ID
!!
!!   OUTPUT ARGUMENT LIST:
!!     ISTA - FIRST LOOP VALUE
!!     IEND - LAST LOOP VALUE
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
      SUBROUTINE PARA_RANGE2 (N1,N2,i1,i2,NPROCS,IRANK,ISTAJ,IENDJ,isx,iex)

      implicit none
      integer,intent(in)  ::  n1,n2,nprocs,irank,i1,i2
      integer,intent(out) ::  istaj,iendj,isx,iex
      integer iwork1, iwork2

      iwork1 = ( n2 - n1 + 1 ) / nprocs
      iwork2 = mod ( n2 - n1 + 1, nprocs )
      istaj   = irank * iwork1 + n1 + min ( irank, iwork2 )
      iendj   = istaj + iwork1 - 1
      if ( iwork2 > irank ) iendj = iendj + 1
      isx=i1
      iex=i2
      print 101,' GWVX para_range2 irank,iwork1,iwork2,istaj,iendj,i1,i2,isx,iex',irank,iwork1,iwork2,istaj,iendj,i1,i2,isx,iex
 101   format( a70,11i8)
      return
      end

