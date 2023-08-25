!> @file
!> @brief process() is a driver for major post routines.
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
!> 2023-01-24 | Sam Trahan   | run IFI and compute its runtime
!> 2023-08-24 | Yali Mao     | Remove running MDL2STD_P
!>
!> @author Russ Treadon W/NP2 @date 1992-12-21
!----------------------------------------------------------------------------
!> process() is a driver for major post routines.
!>
!> @param[in] kth integer Number of isentropic levels.
!> @param[in] kpv integer Number of potential vorticity levels.
!> @param[in] th real Isentropic levels (K).
!> @param[in] pv real Potential vorticity (PV units).
!> @param[in] iostatusD3D integer No longer used/supported.
!> 
      SUBROUTINE PROCESS(kth,kpv,th,pv,iostatusD3D)
!
!----------------------------------------------------------------------------
!
      use mpi, only: mpi_wtime
      use upp_ifi_mod, only: run_ifi
      use CTLBLK_mod, only: cfld, etafld2_tim, eta2p_tim, mdl2sigma_tim, surfce2_tim,&
                            mdl2agl_tim, mdl2std_tim, mdl2thandpv_tim, calrad_wcloud_tim,&
                            cldrad_tim, miscln_tim, fixed_tim, ntlfld, me, run_ifi_tim
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
      if(me==0) write(*,*) "PROCESS starts"
!
!     COMPUTE/POST FIELDS ON MDL SURFACES.
!
      btim = mpi_wtime()
      CALL MDLFLD
      if(me==0) write(*,*) "PROCESS MDLFLD done"
      ETAFLD2_tim = ETAFLD2_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON PRESSURE SURFACES.
      btim = mpi_wtime()
      CALL MDL2P(iostatusD3D)
      if(me==0) write(*,*) "PROCESS MDL2P done"
      ETA2P_tim = ETA2P_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON SIGMA SURFACES.
      btim = mpi_wtime()
      CALL MDL2SIGMA
      if(me==0) write(*,*) "PROCESS MDL2SIGMA done"
      CALL MDL2SIGMA2
      if(me==0) write(*,*) "PROCESS MDL2SIGMA2 done"
      MDL2SIGMA_tim = MDL2SIGMA_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON AGL SURFCES.
      btim = mpi_wtime()
      CALL MDL2AGL
      if(me==0) write(*,*) "PROCESS MDL2AGL done"
      MDL2AGL_tim = MDL2AGL_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST SURFACE RELATED FIELDS.
      btim = mpi_wtime()
      CALL SURFCE
      if(me==0) write(*,*) "PROCESS SURFCE done"
      SURFCE2_tim = SURFCE2_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST SOUNDING AND CLOUD RELATED FIELDS.
      btim = mpi_wtime()
      CALL CLDRAD
      if(me==0) write(*,*) "PROCESS CLDRAD done"
      CLDRAD_tim = CLDRAD_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST TROPOPAUSE DATA, FD LEVEL FIELDS,
!     FREEZING LEVEL HEIGHT AND RH, BOUNDARY LAYER FIELDS,
!     AND LFM-NGM LOOK-ALIKE FIELDS.
      btim = mpi_wtime()
      CALL MISCLN
      if(me==0) write(*,*) "PROCESS MISCLN done"
      MISCLN_tim = MISCLN_tim +(mpi_wtime() - btim)

!     COMPUTE/POST TROPOPAUSE DATA, FD LEVEL FIELDS,
!     FREEZING LEVEL HEIGHT AND RH, BOUNDARY LAYER FIELDS,
!     AND LFM-NGM LOOK-ALIKE FIELDS.
!      btim = mpi_wtime()
!      CALL MDL2STD_P
!      if(me==0) write(*,*) "PROCESS MDL2STD_P done"
!      MDL2STD_tim = MDL2STD_tim +(mpi_wtime() - btim)

!     POST FIXED FIELDS.
      btim = mpi_wtime()
      CALL FIXED
      if(me==0) write(*,*) "PROCESS FIXED done"
      FIXED_tim =  FIXED_tim +(mpi_wtime() - btim)
!
!     COMPUTE/POST FIELDS ON SIGMA SURFACES.
      btim = mpi_wtime()
      CALL MDL2THANDPV(kth,kpv,th,pv)
      if(me==0) write(*,*) "PROCESS MDL2THANDPV done"
      MDL2THANDPV_tim = MDL2THANDPV_tim +(mpi_wtime() - btim)
!
!     POST RADIANCE AND BRIGHTNESS FIELDS.
      btim = mpi_wtime()
      CALL CALRAD_WCLOUD
      if(me==0) write(*,*) "PROCESS CALRAD_WCLOUD done"
      CALRAD_WCLOUD_tim = CALRAD_WCLOUD_tim +(mpi_wtime() - btim)
!
!     IN-FLIGHT ICING PRODUCTS
      btim = mpi_wtime()
      CALL RUN_IFI
      RUN_IFI_tim = RUN_IFI_tim +(mpi_wtime()-btim)
!
!     END OF ROUTINE.
!
      NTLFLD=cfld
      if(me==0)print *,'nTLFLD=',NTLFLD
      if(me==0) write(*,*) "PROCESS done"
!
      RETURN
      END
