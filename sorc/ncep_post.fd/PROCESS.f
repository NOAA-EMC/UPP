!> @file
!> @brief PROCESS is a driver for major post routines.
!>
!> This routine calls the major post processor routines.
!> <pre>
!> These routines are
!> MDLFLD  - Calculate NMC SLP, set below surface fields,
!>           and posts data on model surfaces.
!> MDL2P   - Posts data on isobaric surfaces.
!> SURFCE  - Posts sounding data  surface based fields,
!>           and static or fixed fields.
!> CLDRAD  - Post sounding/cloud/radiation fields.
!> MISCLN  - Post miscellaneous (special) fields.
!> FIXED   - Post fixed fields.
!> </pre>
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-21 | Russ Treadon | Initial
!> 1998-06-01 | T Black      | Conversion from 1-D to 2-D
!> 2000-01-05 | Jim Tuccillo | MPI Version
!> 2001-10-25 | H CHUANG     | Modified to process hybrid model output
!> 2002-06-19 | Mike Baldwin | WRF Version
!> 2011-02-04 | Jun Wang     | Add grib2 option
!>
!> @author Russ Treadon W/NP2 @date 1992-12-21
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
!
!     COMPUTE/POST FIELDS ON MDL SURFACES.
!
      btim = mpi_wtime()
      CALL MDLFLD
      ETAFLD2_tim = ETAFLD2_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON PRESSURE SURFACES.
      btim = mpi_wtime()
      CALL MDL2P(iostatusD3D)
      ETA2P_tim = ETA2P_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON SIGMA SURFACES.
      btim = mpi_wtime()
      CALL MDL2SIGMA
      CALL MDL2SIGMA2
      MDL2SIGMA_tim = MDL2SIGMA_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON AGL SURFCES.
      btim = mpi_wtime()
      CALL MDL2AGL
      MDL2AGL_tim = MDL2AGL_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST SURFACE RELATED FIELDS.
      btim = mpi_wtime()
      CALL SURFCE
      SURFCE2_tim = SURFCE2_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST SOUNDING AND CLOUD RELATED FIELDS.
      btim = mpi_wtime()
      CALL CLDRAD
      CLDRAD_tim = CLDRAD_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST TROPOPAUSE DATA, FD LEVEL FIELDS,
!     FREEZING LEVEL HEIGHT AND RH, BOUNDARY LAYER FIELDS,
!     AND LFM-NGM LOOK-ALIKE FIELDS.
      btim = mpi_wtime()
      CALL MISCLN
      MISCLN_tim = MISCLN_tim +(mpi_wtime() - btim)

!     COMPUTE/POST TROPOPAUSE DATA, FD LEVEL FIELDS,
!     FREEZING LEVEL HEIGHT AND RH, BOUNDARY LAYER FIELDS,
!     AND LFM-NGM LOOK-ALIKE FIELDS.
      btim = mpi_wtime()
      CALL MDL2STD_P
      MDL2STD_tim = MDL2STD_tim +(mpi_wtime() - btim)
!
!     POST FIXED FIELDS.
      btim = mpi_wtime()
      CALL FIXED
      FIXED_tim =  FIXED_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON SIGMA SURFACES.
      btim = mpi_wtime()
      CALL MDL2THANDPV(kth,kpv,th,pv)
      MDL2THANDPV_tim = MDL2THANDPV_tim +(mpi_wtime() - btim)
!
!     POST RADIANCE AND BRIGHTNESS FIELDS.
      btim = mpi_wtime()
      CALL CALRAD_WCLOUD
      CALRAD_WCLOUD_tim = CALRAD_WCLOUD_tim +(mpi_wtime() - btim)
!
!     END OF ROUTINE.
!
      NTLFLD=cfld
      if(me==0)print *,'nTLFLD=',NTLFLD
!
      RETURN
      END
