!> @file
!
!> SUBPROGRAM:    PROCESS     DRIVER FOR MAJOR POST ROUTINES.
!!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-21
!!
!! ABSTRACT:
!!     THIS ROUTINE CALLS THE MAJOR POST PROCESSOR ROUTINES.
!!     THESE ROUTINES ARE
!!        MDLFLD  - CALCULATE NMC SLP, SET BELOW SURFACE FIELDS,
!!                  AND POSTS DATA ON MODEL SURFACES.
!!        MDL2P   - POSTS DATA ON ISOBARIC SURFACES.
!!        SURFCE  - POSTS SOUNDING DATA, SURFACE BASED FIELDS,
!!                  AND STATIC OR FIXED FIELDS.
!!        CLDRAD  - POST SOUNDING/CLOUD/RADIATION FIELDS.
!!        MISCLN  - POST MISCELLANEOUS (SPECIAL) FIELDS.
!!        FIXED   - POST FIXED FIELDS.
!!
!! PROGRAM HISTORY LOG:
!!   92-12-21  RUSS TREADON
!!   98-06-01  T BLACK - CONVERSION OF POST FROM 1-D TO 2-D
!!   00-01-05  JIM TUCCILLO - MPI VERSION
!!   01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!!   02-06-19  MIKE BALDWIN - WRF VERSION
!!   11-02-04  Jun Wang - add grib2 option
!!
!! USAGE:    CALL PROCESS
!!   INPUT ARGUMENT LIST:
!!     NONE
!!
!!   OUTPUT ARGUMENT LIST:
!!     NONE
!!
!!   OUTPUT FILES:
!!     NONE
!!
!!   SUBPROGRAMS CALLED:
!!     UTILITIES:
!!       MDLFLD   - POST DATA MDL SURFACES.
!!       MDL2P    - POST DATA ON PRESSURE SURFACES.
!!       SURFCE   - POST SURFACE BASED FIELDS.
!!       CLDRAD   - POST SOUNDING/CLOUD/RADIATION FIELDS.
!!       MISCLN   - POST MISCELLANEOUS FIELDS.
!!       FIXED    - POST FIXED FIELDS.
!!     LIBRARY:
!!       COMMON   - OUTGRD
!!
!!   ATTRIBUTES:
!!     LANGUAGE: FORTRAN
!!     MACHINE : CRAY C-90
!!
      SUBROUTINE PROCESS(kth,kpv,th,pv,iostatusD3D)
!
!----------------------------------------------------------------------------
!
      use mpi, only: mpi_wtime

      use CTLBLK_mod, only: cfld, etafld2_tim, eta2p_tim, mdl2sigma_tim, surfce2_tim,&
                            mdl2agl_tim, mdl2std_tim, mdl2thandpv_tim, calrad_wcloud_tim,&
                            cldrad_tim, miscln_tim, fixed_tim, ntlfld, me
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
!------ DECLARE VARIABLES.
!
      integer,intent(in) :: kth
      integer,intent(in) :: kpv
      integer,intent(in) :: iostatusD3D
      real,intent(in)    :: th(kth)
      real,intent(in)    :: pv(kpv)
      real(kind=8)       :: btim
      CHARACTER*6           DATSET,PROJ
      LOGICAL               NORTH
!
!
!****************************************************************************
!     START SUBROUTINE PROCESS.
!
      cfld=0
      if(me==0) write(0,*) "PROCESS starts"
!
!     COMPUTE/POST FIELDS ON MDL SURFACES.
!
      btim = mpi_wtime()
      CALL MDLFLD
      if(me==0) write(0,*) "PROCESS MDLFLD done"
      ETAFLD2_tim = ETAFLD2_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON PRESSURE SURFACES.
      btim = mpi_wtime()
      CALL MDL2P(iostatusD3D)
      if(me==0) write(0,*) "PROCESS MDL2P done"
      ETA2P_tim = ETA2P_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON SIGMA SURFACES.
      btim = mpi_wtime()
      CALL MDL2SIGMA
      if(me==0) write(0,*) "PROCESS MDL2SIGMA done"
      CALL MDL2SIGMA2
      if(me==0) write(0,*) "PROCESS MDL2SIGMA2 done"
      MDL2SIGMA_tim = MDL2SIGMA_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON AGL SURFCES.
      btim = mpi_wtime()
      CALL MDL2AGL
      if(me==0) write(0,*) "PROCESS MDL2AGL done"
      MDL2AGL_tim = MDL2AGL_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST SURFACE RELATED FIELDS.
      btim = mpi_wtime()
      CALL SURFCE
      if(me==0) write(0,*) "PROCESS SURFCE done"
      SURFCE2_tim = SURFCE2_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST SOUNDING AND CLOUD RELATED FIELDS.
      btim = mpi_wtime()
      CALL CLDRAD
      if(me==0) write(0,*) "PROCESS CLDRAD done"
      CLDRAD_tim = CLDRAD_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST TROPOPAUSE DATA, FD LEVEL FIELDS,
!     FREEZING LEVEL HEIGHT AND RH, BOUNDARY LAYER FIELDS,
!     AND LFM-NGM LOOK-ALIKE FIELDS.
      btim = mpi_wtime()
      CALL MISCLN
      if(me==0) write(0,*) "PROCESS MISCLN done"
      MISCLN_tim = MISCLN_tim +(mpi_wtime() - btim)

!     COMPUTE/POST TROPOPAUSE DATA, FD LEVEL FIELDS,
!     FREEZING LEVEL HEIGHT AND RH, BOUNDARY LAYER FIELDS,
!     AND LFM-NGM LOOK-ALIKE FIELDS.
      btim = mpi_wtime()
      CALL MDL2STD_P
      if(me==0) write(0,*) "PROCESS MDL2STD_P done"
      MDL2STD_tim = MDL2STD_tim +(mpi_wtime() - btim)
!
!     POST FIXED FIELDS.
      btim = mpi_wtime()
      CALL FIXED
      if(me==0) write(0,*) "PROCESS FIXED done"
      FIXED_tim =  FIXED_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON SIGMA SURFACES.
      btim = mpi_wtime()
      CALL MDL2THANDPV(kth,kpv,th,pv)
      if(me==0) write(0,*) "PROCESS MDL2THANDPV done"
      MDL2THANDPV_tim = MDL2THANDPV_tim +(mpi_wtime() - btim)
!
!     POST RADIANCE AND BRIGHTNESS FIELDS.
      btim = mpi_wtime()
      CALL CALRAD_WCLOUD
      if(me==0) write(0,*) "PROCESS CALRAD_WCLOUD done"
      CALRAD_WCLOUD_tim = CALRAD_WCLOUD_tim +(mpi_wtime() - btim)
!
!     END OF ROUTINE.
!
      NTLFLD=cfld
      if(me==0)print *,'nTLFLD=',NTLFLD
      if(me==0) write(0,*) "PROCESS done"
!
      RETURN
      END
