!!!@PROCESS NOEXTCHK
      SUBROUTINE MPI_FIRST()
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    MPI_FIRST   SET UP MESSGAE PASSING INFO
!   PRGRMMR: TUCCILLO        ORG: IBM
!
! ABSTRACT:
!     SETS UP MESSAGE PASSING INFO
!   .
!
! PROGRAM HISTORY LOG:
!   00-01-06  TUCCILLO - ORIGINAL
!   01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-06-19  MIKE BALDWIN - WRF VERSION
!   11-12-16  SARAH LU - MODIFIED TO INITIALIZE AEROSOL FIELDS
!   12-01-07  SARAH LU - MODIFIED TO INITIALIZE AIR DENSITY/LAYER THICKNESS
!
! USAGE:    CALL MPI_FIRST
!   INPUT ARGUMENT LIST:
!
!   OUTPUT ARGUMENT LIST:
!
!   OUTPUT FILES:
!     STDOUT  - RUN TIME STANDARD OUT.
!
!   SUBPROGRAMS CALLED:
!       PARA_RANGE
!     UTILITIES:
!       NONE
!     LIBRARY:
!       COMMON - CTLBLK.comm
!
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : IBM RS/6000 SP
!$$$
!
      use vrbls4d
      use vrbls3d
      use vrbls2d
      use soil
      use masks
!
      use params_mod
      use ctlblk_mod
!- - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - 
      implicit none
!
      include 'mpif.h'
!
      integer ierr,i,jsx,jex
!
      if ( me .eq. 0 ) then
!        print *, ' NUM_PROCS = ',num_procs
      end if

      if ( num_procs .gt. 1024 ) then
         print *, ' too many MPI tasks, max is 1024, stopping'
         call mpi_abort(MPI_COMM_WORLD,1,ierr)
         stop
      end if
!
!     error check
!
      if ( num_procs .gt. JM/2 ) then
         print *, ' too many MPI tasks, max is ',jm/2,' stopping'
         call mpi_abort(MPI_COMM_WORLD,1,ierr)
         stop
      end if
!
!     global loop ranges
!
      call para_range(1,jm,num_procs,me,  &
        jsta,jend)
      jsta_m  = jsta
      jsta_m2 = jsta
      jend_m  = jend
      jend_m2 = jend
      if ( me .eq. 0 ) then
         jsta_m  = 2
         jsta_m2 = 3
      end if
      if ( me .eq. num_procs - 1 ) then
         jend_m  = jm - 1
         jend_m2 = jm - 2
      end if
!
!     neighbors
!
      iup = me + 1
      idn = me - 1
      if ( me .eq. 0 ) then
         idn = MPI_PROC_NULL
      end if
      if ( me .eq. num_procs - 1 ) then
         iup = MPI_PROC_NULL
      end if
!
!     print *, ' ME, NUM_PROCS = ',me,num_procs
!     print *, ' ME, JSTA, JSTA_M, JSTA_M2 = ',me,jsta,jsta_m,jsta_m2
!     print *, ' ME, JEND, JEND_M, JEND_M2 = ',me,jend,jend_m,jend_m2
!     print *, ' ME, IUP, IDN = ',me,iup,idn
!
!     counts, disps for gatherv and scatterv
!
      do i = 0, num_procs - 1
         call para_range(1,jm,num_procs,i,jsx,jex) 
         icnt(i) = (jex-jsx+1)*im
         idsp(i) = (jsx-1)*im
         if ( me .eq. 0 ) then
           print *, ' i, icnt(i),idsp(i) = ',i,icnt(i),      &
            idsp(i)
         end if
      end do
!
!     extraction limits -- set to two rows    
!
      jsta_2l = max(jsta - 2,  1 )
      jend_2u = min(jend + 2, jm )
! special for c-grid v
      jvend_2u = min(jend + 2, jm+1 )
! special for c-grid v
!     print *, ' me, jvend_2u = ',me,jvend_2u
!
!     allocate arrays
!
!
!     FROM VRBLS3D
!
      print *, ' me, jsta_2l, jend_2u = ',me,jsta_2l, jend_2u,  &
               'jvend_2u=',jvend_2u,'im=',im,'jm=',jm,'lm=',lm, &
               'lp1=',lp1

      end
