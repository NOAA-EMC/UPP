!> @file
!> @brief module: CTLBLK sets default parameters that are used throughout the UPP code
!>
!> ABSTRACT: 
!> This module is replacing the CTLBLK.comm, all the comm block is removed.
!> 
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!>  2011-02    | Jun Wang | ADD variables for grib2
!>  2011-12-14 | SARAH LU | ADD AER FILENAME
!>  2011-12-23 | SARAH LU | ADD NBIN FOR DU, SS, OC, BC, SU
!>  2021-09-30 | JESSE MENG | 2D DECOMPOSITION
!>  2022-09-22 | Li(Kate) Zhang | Add option for NASA GOCART as "nasa_on", add NBIN for NO3 and NH4
!>  2022-11-08 | Kai Wang | Replace aqfcmaq_on with aqf_on
!>  2023-01-24 | Sam Trahan | IFI flight levels, runtime of IFI, and record of the bucket time
!>  2023-03-21 | Jesse Meng | Add slrutah_on option to use U Utah SLR
!>  2023-04-04 | Li(Kate Zhang) | Add namelist optoin for CCPP-Chem (UFS-Chem) and 2D diag. output (d2d_chem) for GEFS-Aerosols and CCPP-Chem model.
!>  2023-04-17 | Eric James | Adding 160 and 320 m above ground to HTFD for RRFS output.
!>  2023-08-16 | Yali Mao   | Add gtg_on logical option
!>  2023-11-24 | Eric James | Add method_blsn logical option
!-----------------------------------------------------------------------
!> @defgroup CTLBLK_mod Sets default parameters that are used throughout the UPP code
!-----------------------------------------------------------------------
  module CTLBLK_mod
!
  implicit none
!
  type field_info
    integer ifld        !< Field number in post control file
    integer lvl         !< _____.
    integer lvl1        !< _____.
    integer lvl2        !< _____.
    integer ntrange     !< _____.
    integer tinvstat    !< _____.
  end type
!
  integer, parameter :: komax=70          !< _____.
  integer, parameter :: LSMDEF=46         !< Default number of pressure levels.
  integer,PARAMETER  :: NFD=20            !< Default number of flight level heights in geopotential meters. 
  integer,PARAMETER  :: NBND=6            !< Default number of ETA boundary layers.
  REAL,  PARAMETER   :: QMIN = 1.E-15     !< A minimum specific humidity value.
!
  integer :: novegtype     !< Number of vegetation types based on vegetation classification.
!
  character(len=256) :: fileName                   !< Name of input dynamics file; name of full 3-D model output file.
  character(len=256) :: fileNameFlux               !< Name of input physics file; name of 2-D model output file with physics and surface fields.
  character(len=256) :: fileNameD3D                !< _____.
  character(len=256) :: fileNameAER                !< _____.
  character(len=256) :: fileNameFlat               !< Input configuration text file defining the requested fields.
  character(len=19)  :: DateStr                    !< Time stamp being processed (e.g., 2022-08-02_19:00:00).
  character(len=4)   :: MODELNAME                  !< Model name used by UPP internally (e.g., FV3R for LAM, GFS for GFS, NCAR for WRF).
  character(len=4)   :: SUBMODELNAME               !< Name of submodel for output differing from parent domain; used to treat a subset of model output in a special way; typically used only for outputting RTMA-specific fields now; previously used with HWRF and NMM to identify and process moving nests.
  character(len=8)   :: FULLMODELNAME              !< No longer used/supported.
  character(len=20)  :: IOFORM                     !< Input file format.
  character(len=4)   :: VTIMEUNITS                 !< Valid time units.
! 
  character(5) :: grib                          !< Grib type (Note that UPP only supports Grib2 currently).
  type(field_info),allocatable :: fld_info(:)   !< _____.
  integer :: cfld                               !< _____.
  integer :: ntlfld                             !< _____.
  integer :: npset                              !< _____.
  real*8 :: gdsdegr                             !< _____.
  real,allocatable :: datapd(:,:,:)             !< _____.
!
!> Logicals to turn on/off different post-processing packages/output depending on model output.
  logical :: gocart_on     !< Turn on option to process the aerosol/chemical tracers related output from GEFS-Aerosols model (GOCART).
  logical :: gccpp_on      !< Turn on option to process the aerosol/chemical tracers related output from UFS-Chem (CCPP-Chem) model.
  logical :: nasa_on       !< Turn on option to process the aerosol/chemical tracers related output from UFS-Aerosols model (NASA GOCART).
  logical :: d3d_on        !< _____.
  logical :: hyb_sigp      !< _____.
  logical :: rdaod         !< Turn on option to process the AOD from GFS scheme.
  logical :: d2d_chem      !< Turn on option to process the 2D aerosol/chemical tracers.
  logical :: aqf_on        !< Turn on Air Quality Forecasting (CMAQ-based).
  logical :: slrutah_on    !< Calculate snow to liquid ratio (SLR) using method from University of Utah.
  logical :: gtg_on        !< Turn on GTG (Graphical Turbulence Guidance)
  logical :: method_blsn   !< Turn on blowing snow effect on visibility diagnostic
!
  logical :: SIGMA      !< No longer used/supported.
  logical :: RUN        !< No longer used/supported.
  logical :: FIRST      !< No longer used/supported.
  logical :: RESTRT     !< Indicates whether it is a restart run.
  logical :: global     !< _____.
  logical :: SMFLAG     !< Smoothing flag for isobaric output.
!
  integer :: IDAT(5)             !< Array storing input month, day, year, hour, min of file being processed (parsed from DateStr)
  integer :: IHRST               !< Hour of file being processed (parsed from DateStr).
  integer :: IFHR                !< Forecast hour (lead time).
  integer :: IMIN                !< Minute of file being processed (parsed from DateStr).
  integer :: ifmin               !< Forecast minute.
  integer :: imp_physics         !< Microphysics option used in the model run.
  integer :: icu_physics         !< Cumulus physics option in the model run.
  integer :: iSF_SURFACE_PHYSICS !< Surface physics scheme option in model run.
  integer :: DataHandle          !< _____. 
  integer :: NPREC               !< _____. 
  integer :: NPHS                !< _____. 
  integer :: ISEC                !< Seconds of file being processed (not parsed from DateStr, hard-coded set to 0).
  integer :: icount_calmict      !< _____. 
  integer :: ivegsrc             !< Flag for vegetation classification source (0=USGS, 1=IGBP, 2=UMD)

!
!> @ingroup CTLBLK_mod 
!> @{
!> No longer used/supported.
   integer :: NFCST,NBC,LIST,IOUT,NTSTM,                 &
             NRADS,NRADL,NDDAMP,IDTAD,NBOCO,NSHDE,NCP,IMDLTY
!> @}
!
  real :: DT            !< Model time step in seconds.
  real :: SDAT(3)       !< Array of month, day, year of restart run.
  real :: AVRAIN        !< Counter for summing latent heating from grid microphysics.
  real :: AVCNVC        !< Counter for summing latent heating from convection.
  real :: DTQ2          !< Model physics time step in seconds.
  real :: PT            !< Model top requested by CMAQ.
  real :: PDTOP         !< Pressure thickness requested by CMAQ.
  real :: SPL(komax)    !< _____.
  real :: ALSL(komax)   !< _____.
  real :: PREC_ACC_DT   !< _____.
  real :: PT_TBL        !< _____.
  real :: PREC_ACC_DT1  !< _____.
  real :: spval         !< _____.
! real :: SPVAL=9.9e10                                     ! Moorthi
!
  integer :: NUM_PROCS        !< The number of MPI ranks available to the post processor.
  integer :: ME               !< MPI rank.
  integer :: JSTA             !< Start latitude on a task subdomain.
  integer :: JEND             !< End latitude on a task subdomain.
  integer :: ISTA             !< Start longitude latitude on a task subdomain.
  integer :: IEND             !< End longitude on a task subdomain.
  integer :: JSTA_M           !< Beginning latitude loop index in subdomain for halo depth 1.
  integer :: JEND_M           !< Ending latitude loop index in subdomain for halo depth 1.
  integer :: JSTA_M2          !< Second latitude below begin latitude of subdomain for halo depth 2 (in NGMFLD.f).
  integer :: JEND_M2          !< Second latitude above end latitude of subdomain for halo depth 2 (in NGMFLD.f).
  integer ::  ISTA_M          !< Beginning longitude loop index in subdomain for halo depth 1.
  integer ::  IEND_M          !< Ending longitude loop index in subdomain for halo depth 1.
  integer ::  ISTA_M2         !< Second longitude before begin longitude for halo depth 2 (not used as of 6/22).
  integer ::  IEND_M2         !< Second longitude after end longitude for halo depth 2 (not used as of 6/22).
  integer :: IUP              !< MPI rank containing the first latitude after jend.
  integer :: IDN              !< MPI rank containing the last latitude before jsta.

!> Used for gathers and scatters 
  integer :: ICNT(0:1023)     !< The number of data items to scatter to each MPI rank; it is a NUM_PROCS array.
  integer :: IDSP(0:1023)     !< Displacement in the array to be scattered where the portion of the array to be scattered to each MPI rank begins.
  integer :: ICNT2(0:1023)    !< The number of data items to gather from each MPI rank; it is a NUM_PROCS array.
  integer :: IDSP2(0:1023)    !< Displacement in the array to be gathered where the portion of the array to be gathered from each MPI rank begins.
  integer :: JSTA_2L          !< Start latitude -2 of the subdomain.
  integer :: JEND_2U          !< End latitude +2 of the subdomain.
  integer :: JVEND_2U         !< Defines the upper boundary for the subdomain used on each MPI rank. Includes information from neighboring ranks (halos).   
  integer :: ISTA_2L          !< Start longitude -2 of the subdomain.
  integer :: IEND_2U          !< End longitude +2 of the subdomain.
  integer :: IVEND_2U         !< Defines the right most boundary for the subdomain used on each MPI rank. Includes information from neighboring ranks (halos).
  integer :: NUM_SERVERS      !< An optional variable to support asynchronous writes of post-processed fields; one if there is more than one total MPI task - otherwise zero; note that the asynchronous write code is not in active development or used.
  integer :: MPI_COMM_INTER   !< An MPI communicator defining a subgroup of the MPI ranks used for asynchronous I/O; asynchronous writes are not in active development.
  integer :: MPI_COMM_COMP    !< an MPI communicator defining the subgroup of MPI ranks used to compute post-processed product fields; all current post implementations use all of the ranks so this again supports an unexploited development path in the code.
  integer :: IM               !< Full longitude domain.
  integer :: JM               !< Full latitude domain.
  integer :: LM               !< Number of vertical levels.
  integer :: NSOIL            !< Number of model soil levels (dependent on the land surface model used).
  integer :: LP1              !< LM+1.
  integer :: LM1              !< LM-1.
  integer :: IM_JM            !< The product of IM and JM, which defines the number of points in the full post domain.
  integer :: ileft            !< MPI rank containing the last longitude before ista.
  integer :: iright           !< MPI rank containing the first longitude after iend.
  integer :: ileftb           !< MPI rank containing the last longitude before ista but for cyclic boundary conditions where "last" at the beginning is the other end of the domain (apparently unused and replaced with local calculation).
  integer :: irightb          !< MPI rank containing the first longitude after iend but for cyclic boundary conditions where "first" at the beginning is the other end of the domain (apparently unused and replaced with local calculation).
  integer :: ibsize           !< Defines the size of the buffer used in mpi_scatter and mpi_gather. It is necessary because the post-processed variables are not contiguous in the 2D ista_2l:iend_2u,jsta_2l:jsta_2u arrays, so they have to be stored in a contigous buffer and that buffer is what is scattered or gathered.
  integer :: ibsum            !< No longer supported.
  !comm mpi
  integer :: lsm              !< _____.
  integer :: lsmp1            !< LSM+1.
!
!> @ingroup CTLBLK_mod
!> @{
!> Arrays that store the coordinates of their elements; used to validate communications; 
!> when scattered or otherwise dispersed, the receiving ranks check that the values of 
!> the arrays match the I and J indices of the receiver.
  integer, allocatable :: icoords(:,:)
  integer, allocatable :: ibcoords(:,:)
  real, allocatable :: rcoords(:,:)
  real, allocatable :: rbcoords(:,:)
!> @}

  real, allocatable :: bufs(:)            !< Unused/no longer supported; replaced by rbufs. 
  real, allocatable :: buff(:)            !< Used in the many variables' gather; note that scattering has been replaced with subdomain reads when the fields to be post-processed are read in.
  integer, allocatable :: isxa(:)         !< Array of i start bounds for the subdomain loop on each MPI rank.
  integer, allocatable :: iexa(:)         !< Array of i end bounds for the subdomain loop on each MPI rank.
  integer, allocatable :: jsxa(:)         !< Array of j start bounds for the subdomain loop on each MPI rank.
  integer, allocatable :: jexa(:)         !< Array of j end bounds for the subdomain loop on each MPI rank.
  integer numx                            !< The number of i regions in a 2D decomposition; Each i row is distibuted to numx ranks; numx=1 is the special case of a 1D decomposition in Y only.
  integer, allocatable :: ibufs(:)        !< The buffer used for scatters of the integer coordinate array to each MPI rank.
  real, allocatable :: rbufs(:)           !< The buffer used for scatters of the real coordinate array to each MPI rank; analagous to buff in the state variable scatter. 
!
!comm rad
  real :: ARDSW      !< Shortwave flux accumulation array.
  real :: ARDLW      !< Longwave flux accumulation array.
  real :: ASRFC      !< Surface flux array.
  real :: TSRFC      !< Number of hours in surface flux buckets.
  real :: TRDLW      !< Number of hours in long wave buckets.
  real :: TRDSW      !< Number of hours in shortwave buckets.
  real :: TCLOD      !< Number of hours in cloud fraction average.
  real :: THEAT      !< Number of hours in latent heating bucket.
  real :: TPREC      !< Number of hours in precipitation bucket.
  real :: TMAXMIN    !< _____.
  real :: TD3D       !< _____.
!
  real PTHRESH !< Threshold for precipitation (used to check if there is precipitation, mainly in ptype routines).
!
!> @ingroup CTLBLK_mod
!> @{ Time to execute named routine; note that ETAFLD2 and ETA2P refer to MDLFLD and MDL2P routines respectively.
  real(kind=8) :: ETAFLD2_tim=0.,ETA2P_tim=0.,SURFCE2_tim=0.,          &
                  CLDRAD_tim=0.,MISCLN_tim=0.,FIXED_tim=0.,            &
                  MDL2SIGMA_tim=0.,READxml_tim=0.,MDL2AGL_tim=0.,      &
                  MDL2STD_tim=0.,MDL2THANDPV_tim=0.,                   &
                  CALRAD_WCLOUD_tim=0.,RUN_IFI_TIM=0.     !comm tim_info
!> @}
!
!> @ingroup CTLBLK_mod
!> @{
!> Initialized as 0, but never used.
  real(kind=8) :: time_output=0., time_e2out=0.           !comm jjt
!> @}
!
!> SPLDEF   !< The fixed pressure levels available for output (Pa).
  real :: SPLDEF(LSMDEF) =                                             &
      (/200.,500.,700.,1000.,2000.,3000.                               &
      ,5000.,7000.,7500.,10000.,12500.,15000.,17500.,20000.,22500.     &
      ,25000.,27500.,30000.,32500.,35000.,37500.,40000.,42500.,45000.  &
      ,47500.,50000.,52500.,55000.,57500.,60000.,62500.,65000.         &
      ,67500.,70000.,72500.,75000.,77500.,80000.,82500.,85000.         &
      ,87500.,90000.,92500.,95000.,97500.,100000./)
!
  REAL HTFD(NFD)        !< The fixed flight level heights available for output (gpm).
  REAL PETABND(NBND)    !< The fixed ETA levels available for output.
  REAL SIGBND(NBND)     !< The fixed sigma levels available for output.

! Add GOCART aerosol specification
  integer, parameter :: nbin_du = 5   		!< dust
  integer, parameter :: nbin_ss = 5   		!< sea salt
  integer, parameter :: nbin_oc = 2   		!< organic carbon
  integer, parameter :: nbin_bc = 2   		!< black carbon
  integer, parameter :: nbin_su = 1   		!< sulfate
  integer, parameter :: nbin_no3 = 3   	!< nitrate
  integer, parameter :: nbin_nh4 = 1   	!< NH4
  integer, parameter :: nbin_sm = 1       !< smoke
!
!     SET FD LEVEL HEIGHTS IN GEOPOTENTAL METERS.
      DATA HTFD  / 20.E0,30.E0,40.E0,50.E0,80.E0,100.E0,160.E0,305.E0,320.E0,457.E0,610.E0,   &
           914.E0,1524.E0,1829.E0,2134.E0,2743.E0,3658.E0,4572.E0, &
	   6000.E0,7010.E0/
!
!     SET MIDPOINT "SIGMA" VALUES FOR ETA BOUNDARY LAYERS.
      DATA SIGBND / 0.985,0.955,0.925,0.895,0.865,0.835 /
      DATA PETABND / 15.,45.,75.,105.,135.,165./
!
      real :: ITPREC=-1    !< Precipitation bucket time
!  
      integer :: ifi_nflight = 0                !< Number of flight levels
      real, allocatable :: ifi_flight_levels(:) !< Flight levels in feet, provided by libIFI
!     
!-----------------------------------------------------------------------
  end module CTLBLK_mod
