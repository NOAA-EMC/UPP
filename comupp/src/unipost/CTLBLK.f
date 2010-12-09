  module CTLBLK_mod
!-----------------------------------------------------------------------
! module: CTLBLK
!
! ABSTRACT: 
! this module is replacing the CTLBLK.comm, all the comm block is 
!    removed.
!-----------------------------------------------------------------------
!
  implicit none
!
  integer, parameter :: komax=70
  integer, parameter :: LSMDEF=46             ! default number of p levels
!
  character(len=256) :: fileName,fileNameFlux,fileNameD3D
  character(len=19)  :: DateStr
  character(len=4)   :: MODELNAME
  character(len=20)  :: IOFORM
!      
  logical :: SIGMA,RUN,FIRST,RESTRT
  logical :: global
  logical :: SMFLAG
  integer :: IDAT(5),IHRST, NFCST,NBC,LIST,IOUT,IFHR,NTSTM,            &
             NDDAMP,NPREC,IDTAD,NBOCO,NSHDE,NCP,IMDLTY,NPHS,           &
             NRADS,NRADL,IMIN,ifmin,DataHandle,imp_physics,            &
             icu_physics,iSF_SURFACE_PHYSICS
  real :: DT,SDAT(3),AVRAIN,AVCNVC,DTQ2,PT,PDTOP,                &
          SPL(komax),ALSL(komax),PREC_ACC_DT
  real :: SPVAL=9.9e10
!
  integer :: NUM_PROCS,ME,JSTA,JEND,JSTA_M,JEND_M,                     &
             JSTA_M2,JEND_M2,IUP,IDN,ICNT(0:1023),IDSP(0:1023),        &
             JSTA_2L, JEND_2U,JVEND_2u,NUM_SERVERS, MPI_COMM_INTER,    &
             MPI_COMM_COMP, IM,JM,LM,NSOIL,LP1,LM1,IM_JM,              &
             lsm,lsmp1                                    !comm mpi
!
  real :: ARDSW, ARDLW, ASRFC, TSRFC,TRDLW,TRDSW,TCLOD,THEAT,          &
          TPREC,TMAXMIN,TD3D                              !comm rad
!
  real(kind=8) :: ETAFLD2_tim=0.,ETA2P_tim=0.,SURFCE2_tim=0.,          &
                  CLDRAD_tim=0.,MISCLN_tim=0.,FIXED_tim=0.,            &
                  MDL2SIGMA_tim=0.                        !comm tim_info
!
  real(kind=8) :: time_output=0., time_e2out=0.           !comm jjt
!
  real :: SPLDEF(LSMDEF) =                                             &
      (/200.,500.,700.,1000.,2000.,3000.                               &
      ,5000.,7000.,7500.,10000.,12500.,15000.,17500.,20000.,22500.     &
      ,25000.,27500.,30000.,32500.,35000.,37500.,40000.,42500.,45000.  &
      ,47500.,50000.,52500.,55000.,57500.,60000.,62500.,65000.         &
      ,67500.,70000.,72500.,75000.,77500.,80000.,82500.,85000.         &
      ,87500.,90000.,92500.,95000.,97500.,100000./)
!-----------------------------------------------------------------------
  end module CTLBLK_mod
