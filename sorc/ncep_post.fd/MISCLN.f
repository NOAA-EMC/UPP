!> @file
!
!> SUBPROGRAM:    MISCLN      POSTS MISCELLANEOUS FIELDS
!!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-20
!!     
!! ABSTRACT:
!!     THIS ROUTINE HAS BECOME THE CATCH-ALL FOR MISCELLANEOUS
!!     OUTPUT FIELDS POSTED BY THE ETA POST PROCESSOR.  
!!     CURRENTLY THIS ROUTINE POSTS THE FOLLOWING FIELDS:
!!        (1) TROPOPAUSE LEVEL Z,P, T, U, V, AND VERTICAL WIND SHEAR,
!!        (2) MAX WIND LEVEL Z, P, U, AND V,
!!        (3) FD LEVEL T, Q, U, AND V,
!!        (4) FREEZING LEVEL Z AND RH,
!!        (5) CONSTANT MASS (BOUNDARY) FIELDS,
!!        (6) LFM LOOK-ALIKE FIELDS, AND
!!        (7) NGM LOOK-ALIKE FIELDS.
!!
!!     
!! PROGRAM HISTORY LOG:
!!   92-12-20  RUSS TREADON
!!   93-06-19  RUSS TREADON - ADDED TYPE 2 CAPE POSTING.
!!   94-11-07  MIKE BALDWIN - ADDED HELICITY POSTING.
!!   96-03-26  MIKE BALDWIN - CHANGE ETA BOUNDARY LAYER LABELS FOR GRIB
!!   96-11-19  MIKE BALDWIN - BACK OUT PREVIOUS CHANGE 
!!   97-04-25  MIKE BALDWIN - CHANGE ETA BOUNDARY LAYER LABELS FOR GRIB
!!   97-04-29  GEOFF MANIKIN - ADDED TROPOPAUSE HEIGHT AND
!!                             MAX WIND LEVEL FIELDS
!!   98-06-15  T BLACK       - CONVERSION FROM 1-D TO 2-D
!!   98-07-17  MIKE BALDWIN - REMOVED LABL84
!!   00-01-04  JIM TUCCILLO - MPI VERSION
!!   02-04-23  MIKE BALDWIN - WRF VERSION
!!   11-02-06  JUN WANG     - ADD GRIB2 OPTION
!!   11-10-16  SARAH LU     - ADD FD LEVEL DUST/ASH
!!   12-04-03  Jun Wang     - FIXED LVLSXML for fields at FD height (spec_hgt_lvl_above_grnd)
!!   13-05-3   Shrinivas Moorthi - Fix some bugs and make more efficient code
!!   14-02-21  Shrinivas Moorthi - Add more threading
!!   14-02-26  S Moorthi - threading datapd assignment and some cleanup &
!!                         bug fix
!!   15-11-18  S Moorthi - fixed some logical errors in the helicity and
!!   i                     storm motion part of the code
!!   17-06-01  Y Mao - ADD FD levels for GTG(EDPARM CATEDR MWTURB) and allow 
!!                     levels input from control file
!!   19-09-03  J Meng - ADD CAPE related variables for HRRR
!!   20-03-24  J Meng - remove grib1
!!   20-11-10  J Meng - USE UPP_PHYSICS MODULE
!!   21-03-25  E Colon - 3D-RTMA-specific SPC fields added as output
!!   21-04-01  J Meng - computation on defined points only
!!     
!! USAGE:    CALL MISCLN
!!   INPUT ARGUMENT LIST:
!!
!!   OUTPUT ARGUMENT LIST: 
!!     NONE
!!     
!!   SUBPROGRAMS CALLED:
!!     UTILITIES:
!!       TRPAUS  - COMPUTE TROPOPAUSE LEVEL FIELDS.
!!       CALMXW  - COMPUTE MAX WIND LEVEL FIELDS.
!!       SCLFLD  - SCALE ARRAY ELEMENTS BY CONSTANT.
!!       GRIBIT  - OUTPUT FIELD TO GRIB FILE.
!!       CALPOT  - CALCULATE POTENTIAL TEMPERATURE.
!!       FDLVL   - COMPUTE FD LEVEL DATA (AGL OR MSL).
!!       FRZLVL  - COMPUTE FREEZING LEVEL DATA.
!!       BOUND   - BOUND ARRAY ELEMENTS BETWEEN MINIMUM AND MAXIMUM VALUES.
!!       BNDLYR  - COMPUTE BOUNDARY LAYER FIELDS.
!!       CALDWP  - CALCULATE DEWPOINT TEMPERATURE.
!!       OTLFT   - COMPUTE LIFTED INDEX AT 500MB.
!!       CALLCL  - COMPUTE LCL DATA.
!!       LFMFLD  - COMPUTE LFM LOOK-ALIKE FIELDS.
!!       NGMFLD  - COMPUTE NGM LOOK-ALIKE FIELDS.
!!       CALTHTE - COMPUTE THETA-E.
!!       CALHEL  - COMPUTE HELICITY AND STORM MOTION.
!!
!!     LIBRARY:
!!       COMMON - RQSTFLD
!!                CTLBLK
!!     
!!   ATTRIBUTES:
!!     LANGUAGE: FORTRAN
!!     MACHINE : CRAY C-90
!!
      SUBROUTINE MISCLN

!
!
      use vrbls3d,    only: pmid, uh, vh, t, zmid, zint, pint, alpint, q, omga
      use vrbls3d,    only: catedr,mwt,gtg
      use vrbls2d,    only: pblh, cprate, fis, T500, T700, Z500, Z700,&
                            teql,ieql
      use masks,      only: lmh
      use params_mod, only: d00, d50, h99999, h100, h1, h1m12, pq0, a2, a3, a4,    &
                            rhmin, rgamog, tfrz, small, g
      use ctlblk_mod, only: grib, cfld, fld_info, datapd, im, jsta, jend, jm, jsta_m, jend_m, &
                            nbnd, nbin_du, lm, htfd, spval, pthresh, nfd, petabnd, me,&
                            jsta_2l, jend_2u, MODELNAME, SUBMODELNAME
      use rqstfld_mod, only: iget, lvls, id, iavblfld, lvlsxml
      use grib2_module, only: pset
      use upp_physics, only: FPVSNEW,CALRH_PW,CALCAPE,CALCAPE2,TVIRTUAL,        &
                             CALCAPE1D
      use gridspec_mod, only: gridtype
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
!     SET LOCAL PARAMETERS.  MAKE SURE NFD AND NBND AGREE
!     WITH THE VALUES SET IN SUBROUTINES FDLVL AND BNDLYR,
!     RESPECTIVELY.
      real,PARAMETER :: C2K=273.15
      real,parameter :: con_rd      =2.8705e+2 ! gas constant air    (J/kg/K)
      real,parameter :: con_rv      =4.6150e+2 ! gas constant H2O
      real,parameter :: con_eps     =con_rd/con_rv
      real,parameter :: con_epsm1   =con_rd/con_rv-1
      real,parameter :: cpthresh    =0.000004
      real,PARAMETER :: D1000=1000
      real,PARAMETER :: D1500=1500
      real,PARAMETER :: D2000=2000
      real,PARAMETER :: HCONST=42000000. 
      real,PARAMETER :: K2C=273.16
      REAL,PARAMETER :: DM9999=-9999.0

!     
!     DECLARE VARIABLES.
!     
      LOGICAL NORTH, FIELD1,FIELD2
      LOGICAL, dimension(IM,JSTA:JEND) :: DONE, DONE1

      INTEGER, allocatable ::  LVLBND(:,:,:),LB2(:,:)
!     INTEGER LVLBND(IM,JM,NBND),LB2(IM,JM),LPBL(IM,JM)

      real,dimension(im,jm)        :: GRID1, GRID2
      real,dimension(im,jsta:jend) :: P1D, T1D, Q1D, U1D, V1D, SHR1D, Z1D,   &
                                      RH1D, EGRID1, EGRID2, EGRID3, EGRID4,  &
                                      EGRID5, EGRID6, EGRID7, EGRID8, &
                                      MLCAPE,MLCIN,MLLCL,MUCAPE,MUCIN,MUMIXR, &
                                      FREEZELVL,MUQ1D,SLCL,THE,MAXTHE
      integer,dimension(im,jsta:jend) :: MAXTHEPOS
      real, dimension(:,:,:),allocatable :: OMGBND, PWTBND, QCNVBND,   &
                                            PBND,   TBND,   QBND,      &
                                            UBND,   VBND,   RHBND,     &
                                            WBND,   T7D,    Q7D,       &
                                            U7D,    V6D,    P7D,       &
                                            ICINGFD,GTGFD,CATFD,MWTFD
      real, dimension(:,:,:,:),allocatable :: AERFD

      real, dimension(:,:),allocatable ::   QM8510, RH4710, RH8498,    &
                                            RH4796, RH1847, UST, VST,  &
                                            RH3310, RH6610, RH3366,    &
                                            PW3310, RH4410, RH7294,    &
                                            RH4472,                    &
                                            T78483, T89671, P78483, P89671

      REAL, dimension(:,:,:),allocatable :: HELI
      real, dimension(:,:),  allocatable :: USHR1, VSHR1, USHR6, VSHR6, &
                                            MAXWP, MAXWZ, MAXWU, MAXWV, &
                                            MAXWT
      INTEGER,dimension(:,:),allocatable :: LLOW, LUPP
      INTEGER,dimension(:,:),allocatable :: LLOW1, LUPP1
      REAL, dimension(:,:),allocatable   :: CANGLE,ESHR,UVECT,VVECT,&
                                            EFFUST,EFFVST,FSHR,HTSFC,&
                                            ESRH
!
      integer I,J,jj,L,ITYPE,ISVALUE,LBND,ILVL,IFD,ITYPEFDLVL(NFD),    &
              iget1, iget2, iget3, LLMH,imax,jmax,lmax
      real    DPBND,PKL1,PKU1,FAC1,FAC2,PL,TL,QL,QSAT,RHL,TVRL,TVRBLO, &
              ES1,ES2,QS1,QS2,RH1,RH2,ZSF,DEPTH(2),work1,work2,work3, &
              SCINtmp,MUCAPEtmp,MUCINtmp,MLLCLtmp,ESHRtmp,MLCAPEtmp,STP,&
              FSHRtmp,MLCINtmp,SLCLtmp,LAPSE,SHIP

!     Variables introduced to allow FD levels from control file - Y Mao
      integer :: N,NFDCTL
      REAL, allocatable :: HTFDCTL(:)
      integer, allocatable :: ITYPEFDLVLCTL(:)
      integer IE,IW,JN,JS,IVE(JM),IVW(JM),JVN,JVS
      integer ISTART,ISTOP,JSTART,JSTOP,MIDCAL
      real    dummy(IM,jsta:jend)
      integer idummy(IM,jsta:jend)
!     NEW VARIABLES USED FOR EFFECTIVE LAYER
      INTEGER,dimension(:,:),allocatable :: EL_BASE, EL_TOPS
      LOGICAL,dimension(:,:),allocatable :: FOUND_BASE, FOUND_TOPS
      INTEGER,dimension(:,:),allocatable :: L_THETAE_MAX
      INTEGER,dimension(:,:),allocatable :: CAPE9, CINS9
      CHARACTER(LEN=5)   :: IM_CH, JSTA_CH, JEND_CH, ME_CH
      CHARACTER(LEN=60)  :: EFFL_FNAME
      CHARACTER(LEN=60)  :: EFFL_FNAME2
      INTEGER            :: IREC, IUNIT
      INTEGER            :: IREC2, IUNIT2
      LOGICAL            :: debugprint
      INTEGER            :: EL_SCHEME
      INTEGER            :: LLL
      INTEGER            :: LLCL_PAR, LEQL_PAR
      REAL               :: LMASK, PSFC, CAPE_PAR, CINS_PAR, LPAR0
      REAL, DIMENSION(4) :: PARCEL0
      REAL, DIMENSION(:), ALLOCATABLE  :: TPAR_B, TPAR_T   
      REAL, DIMENSION(:), ALLOCATABLE  :: TPAR_TMP 
      REAL, DIMENSION(:), ALLOCATABLE  :: P_AMB, T_AMB, Q_AMB, ZINT_AMB
      REAL, DIMENSION(:,:,:), ALLOCATABLE  :: TPAR_BASE, TPAR_TOPS

!     
!****************************************************************************
!     START MISCLN HERE.
!     
       debugprint = .TRUE.
       EL_SCHEME = 1

         allocate(USHR1(IM,jsta_2l:jend_2u),VSHR1(IM,jsta_2l:jend_2u), &
                  USHR6(IM,jsta_2l:jend_2u),VSHR6(IM,jsta_2l:jend_2u))
         allocate(UST(IM,jsta_2l:jend_2u),VST(IM,jsta_2l:jend_2u),     &
                  HELI(IM,jsta_2l:jend_2u,2),FSHR(IM,jsta_2l:jend_2u))
!
!      HELICITY AND STORM MOTION.
       iget1 = IGET(162)
       iget2 = -1
       iget3 = -1
       if (iget1 > 0) then
         iget2 = LVLS(1,iget1)
         iget3 = LVLS(2,iget1)
       endif
       IF (iget1 > 0 .OR. IGET(163) > 0 .OR. IGET(164) > 0) THEN
         DEPTH(1) = 3000.0
         DEPTH(2) = 1000.0
         CALL CALHEL(DEPTH,UST,VST,HELI,USHR1,VSHR1,USHR6,VSHR6)
         IF (iget2 > 0) then 
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = HELI(I,J,1)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(iget1)
             fld_info(cfld)%lvl=LVLSXML(1,iget1)
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF

         IF (iget3 > 0) then
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = HELI(I,J,2)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(iget1)
             fld_info(cfld)%lvl=LVLSXML(2,iget1)
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF

         IF (IGET(163) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = UST(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(163))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
         IF (IGET(164) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = VST(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(164))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
       ENDIF

!     UPDRAFT HELICITY

       if (IGET(427) > 0) THEN
         CALL CALUPDHEL(GRID1(1,jsta_2l))
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(427))
!$omp parallel do private(i,j,jj)
           do j=1,jend-jsta+1
             jj = jsta+j-1
             do i=1,im
               datapd(i,j,cfld) = GRID1(i,jj)
             enddo
           enddo
         endif

       ENDIF

! CRA  0-1 KM AND 0-6 KM SHEAR

       IF(IGET(430) > 0 .OR. IGET(431) > 0 .OR. IGET(432) > 0      &
                                           .OR. IGET(433) > 0) THEN

         DEPTH = 6000.0
         CALL CALHEL(DEPTH,UST,VST,HELI,USHR1,VSHR1,USHR6,VSHR6)
! 0-6 km shear magnitude
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               FSHR(I,J) = SQRT(USHR6(I,J)**2+VSHR6(I,J)**2)
             ENDDO
           ENDDO
         IF(IGET(430) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = USHR1(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(430))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
         IF(IGET(431) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = VSHR1(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(431))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
         IF(IGET(432) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = USHR6(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(432))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
         IF(IGET(433) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
                GRID1(I,J) = VSHR6(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(433))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
       ENDIF
!
       if (allocated(ushr1)) deallocate(ushr1)
       if (allocated(vshr1)) deallocate(vshr1)
       if (allocated(ushr6)) deallocate(ushr6)
       if (allocated(vshr6)) deallocate(vshr6)
       if (allocated(ust))   deallocate(ust)
       if (allocated(vst))   deallocate(vst)
       if (allocated(heli))  deallocate(heli)
! CRA
!     
!
!     ***BLOCK 1:  TROPOPAUSE P, Z, T, U, V, AND WIND SHEAR.
!    
      IF ((IGET(054)>0).OR.(IGET(055)>0).OR.       &
          (IGET(056)>0).OR.(IGET(057)>0).OR.       &
          (IGET(177)>0).OR.                           &
          (IGET(058)>0).OR.(IGET(108)>0) ) THEN
! Chuang: Use GFS algorithm per Iredell's and DiMego's decision on unification
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM

              if(PMID(I,J,1)<spval) then
! INPUT
              CALL TPAUSE(LM,PMID(I,J,1:LM),UH(I,J,1:LM)             & 
! INPUT
                         ,VH(I,J,1:LM),T(I,J,1:LM),ZMID(I,J,1:LM)    &
! OUTPUT
                         ,P1D(I,J),U1D(I,J),V1D(I,J),T1D(I,J)        &
! OUTPUT
                         ,Z1D(I,J),SHR1D(I,J))                       ! OUTPUT
              else
                P1D(I,J) = spval
                U1D(I,J) = spval
                V1D(I,J) = spval
                T1D(I,J) = spval
                Z1D(I,J) = spval
                SHR1D(I,J) = spval
              endif

            END DO
          END DO
!
!        TROPOPAUSE PRESSURE.
         IF (IGET(054) > 0) THEN
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = P1D(I,J)
               ENDDO
             ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(054))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF

!        ICAO HEIGHT OF TROPOPAUSE
         IF (IGET(399)>0) THEN
           CALL ICAOHEIGHT(P1D, GRID1(1,jsta))
!            print*,'sample TROPOPAUSE ICAO HEIGHTS',GRID1(im/2,(jsta+jend)/2)
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(399))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
          ENDIF

!        TROPOPAUSE HEIGHT.
         IF (IGET(177) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = Z1D(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(177))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
!
!        TROPOPAUSE TEMPERATURE.
         IF (IGET(055) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = T1D(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(055))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
!
!        TROPOPAUSE POTENTIAL TEMPERATURE.
         IF (IGET(108) > 0) THEN
           CALL CALPOT(P1D,T1D,GRID1(1,jsta))
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(108))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
!     
!        TROPOPAUSE U WIND AND/OR V WIND.
         IF ((IGET(056) > 0).OR.(IGET(057) > 0)) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=U1D(I,J)
               GRID2(I,J)=V1D(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            if(IGET(056)>0) then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(056))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
            if(IGET(057)>0) then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(057))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID2(i,jj)
                enddo
              enddo
            endif
           endif
         ENDIF
!
!        TROPOPAUSE WIND SHEAR.
         IF (IGET(058) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = SHR1D(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(058))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
      ENDIF
!
!
!     ***BLOCK 2:  MAX WIND LEVEL  P, Z, U, AND V
!
!        MAX WIND LEVEL CALCULATIONS
      IF ((IGET(173)>0) .OR. (IGET(174)>0) .OR.                &
          (IGET(175)>0) .OR. (IGET(176)>0)) THEN

          allocate(MAXWP(IM,jsta:jend), MAXWZ(IM,jsta:jend),         &
                   MAXWU(IM,jsta:jend), MAXWV(IM,jsta:jend),MAXWT(IM,jsta:jend))
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
            DO I=1,IM
             MAXWP(I,J)=SPVAL
             MAXWZ(I,J)=SPVAL
             MAXWU(I,J)=SPVAL
             MAXWV(I,J)=SPVAL
            ENDDO
           ENDDO

!            CALL CALMXW(MAXWP,MAXWZ,MAXWU,MAXWV,MAXWT)
! Chuang: Use GFS algorithm per Iredell's and DiMego's decision on unification
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
           loopI:DO I=1,IM
            DO L=1,LM
              IF (ABS(PMID(I,J,L)-SPVAL)<=SMALL .OR. &
                  ABS(UH(I,J,L)-SPVAL)<=SMALL .OR. &
                  ABS(UH(I,J,L)-SPVAL)<=SMALL .OR. &
                  ABS(VH(I,J,L)-SPVAL)<=SMALL .OR. &
                  ABS(T(I,J,L)-SPVAL)<=SMALL .OR. &
                  ABS(ZMID(I,J,L)-SPVAL)<=SMALL) cycle loopI
            ENDDO
! INPUT
            CALL MXWIND(LM,PMID(I,J,1:LM),UH(I,J,1:LM)               &
! INPUT
                       ,VH(I,J,1:LM),T(I,J,1:LM),ZMID(I,J,1:LM)      &
! OUTPUT
                       ,MAXWP(I,J),MAXWU(I,J),MAXWV(I,J)             &
! OUTPUT
                       ,MAXWT(I,J),MAXWZ(I,J))
           ENDDO loopI
          END DO 
!        PRESSURE OF MAX WIND LEVEL
         IF (IGET(173) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = MAXWP(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(173))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
!        ICAO HEIGHT OF MAX WIND LEVEL
         IF (IGET(398)>0) THEN
           CALL ICAOHEIGHT(MAXWP, GRID1(1,jsta))
!            print*,'sample MAX WIND ICAO HEIGHTS',GRID1(im/2,(jsta+jend)/2)
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(398))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
!        HEIGHT OF MAX WIND LEVEL
         IF (IGET(174) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = MAXWZ(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(174))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF

!        MAX WIND LEVEL U WIND AND/OR V WIND.
        IF ((IGET(175) > 0).OR.(IGET(176) > 0)) THEN
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              GRID1(I,J) = MAXWU(I,J)
              GRID2(I,J) = MAXWV(I,J)
            ENDDO
          ENDDO
          if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(175))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(176))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID2(i,jj)
              enddo
            enddo
          endif
        ENDIF
!        TEMPERATURE OF MAX WIND LEVEL
        IF (IGET(314) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=MAXWT(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(314))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
        ENDIF
          deallocate(MAXWP,MAXWZ,MAXWU,MAXWV,MAXWT)
      ENDIF

!
!     ***BLOCK 3-1:  FD LEVEL (selected) T, Q, U, AND V.
!     
      IF ( (IGET(059)>0.or.IGET(586)>0).OR.IGET(911)>0.OR.     &
           (IGET(060)>0.or.IGET(576)>0).OR.                     &
           (IGET(061)>0.or.IGET(577)>0).OR.                     &
           (IGET(601)>0.or.IGET(602)>0.or.IGET(603)>0).OR.      &
           (IGET(604)>0.or.IGET(605)>0).OR.                     &
           (IGET(451)>0.or.IGET(578)>0).OR.IGET(580)>0 ) THEN

         ALLOCATE(T7D(IM,JSTA:JEND,NFD), Q7D(IM,JSTA:JEND,NFD),    &
                  U7D(IM,JSTA:JEND,NFD), V6D(IM,JSTA:JEND,NFD),    &
                  P7D(IM,JSTA:JEND,NFD), ICINGFD(IM,JSTA:JEND,NFD) &
                 ,AERFD(IM,JSTA:JEND,NFD,NBIN_DU))

!
!     DETERMINE WHETHER TO DO MSL OR AGL FD LEVELS
!
         ITYPEFDLVL=1
         DO IFD = 1,NFD
           IF (IGET(059)>0) THEN
            IF (LVLS(IFD,IGET(059))>1) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(911)>0) THEN
            IF (LVLS(IFD,IGET(911))>1) ITYPEFDLVL(IFD)=2
           ENDIF
!for grib2, spec hgt only
           IF (IGET(586)>0) THEN
            IF(LVLS(IFD,IGET(586))>0) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(060)>0) THEN
            IF (LVLS(IFD,IGET(060))>1) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(576)>0) THEN
            IF(LVLS(IFD,IGET(576))>0) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(061)>0) THEN
            IF (LVLS(IFD,IGET(061))>1) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(577)>0) then
            if(LVLS(IFD,IGET(577))>0) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(451)>0) THEN
	    IF (LVLS(IFD,IGET(451))>1) ITYPEFDLVL(IFD)=2
	   ENDIF
           IF (IGET(578)>0) then
            if(LVLS(IFD,IGET(578))>0) ITYPEFDLVL(IFD)=2
           ENDIF
	   
	   IF (IGET(580)>0) then
            if(LVLS(IFD,IGET(580))>1) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(587)>0) then
            if(LVLS(IFD,IGET(587))>0) ITYPEFDLVL(IFD)=2
           ENDIF

	   IF (IGET(601)>0) THEN
            IF (LVLS(IFD,IGET(601))>1) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(602)>0) THEN
            IF (LVLS(IFD,IGET(602))>1) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(603)>0) THEN
            IF (LVLS(IFD,IGET(603))>1) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(604)>0) THEN
            IF (LVLS(IFD,IGET(604))>1) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(605)>0) THEN
            IF (LVLS(IFD,IGET(605))>1) ITYPEFDLVL(IFD)=2
           ENDIF

         ENDDO
!         print *,'call FDLVL with ITYPEFDLVL: ', ITYPEFDLVL,'for tmp,lvls=',LVLS(1:15,iget(59)), &
!          'grib2tmp lvs=',LVLS(1:15,iget(586))

         CALL FDLVL(ITYPEFDLVL,T7D,Q7D,U7D,V6D,P7D,ICINGFD,AERFD)
!     
         DO 10 IFD = 1,NFD
!
!           FD LEVEL TEMPERATURE.
            iget1 = IGET(059)
            iget2 = IGET(586)
            if (iget1 > 0) then
              work1 = LVLS(IFD,iget1)
            else
              work1 = 0.0
            endif
            if (iget2 > 0) then
              work2 = LVLS(IFD,iget2)
            else
              work2 = 0.0
            endif
            IF (IGET1 > 0 .or. IGET2 > 0) THEN
              IF (work1 > 0 .or. work2 > 0) THEN
             
!$omp parallel do private(i,j)
                DO J=JSTA,JEND
                  DO I=1,IM
                    GRID1(I,J) = T7D(I,J,IFD)
                  ENDDO
                ENDDO
                IF(work1 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET1)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET1)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                ENDIF
                IF (work2 > 0) THEN
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET2)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET2)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                ENDIF
              ENDIF
            ENDIF

!           FD LEVEL VIRTUAL TEMPERATURE.
            IF (IGET(911)>0) THEN
              IF (LVLS(IFD,IGET(911))>0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 if ( T7D(I,J,IFD) > 600 ) then
                 GRID1(I,J)=SPVAL
                 else
                 GRID1(I,J)=T7D(I,J,IFD)*(1.+0.608*Q7D(I,J,IFD))
                 endif
                 !print *, "grid value ",T7D(I,J,IFD),Q7D(I,J,IFD),T7D(I,J,IFD)*(1.+0.608*Q7D(I,J,IFD)),GRID1(I,J)
               ENDDO
               ENDDO
               IF(LVLS(IFD,IGET(911))>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(911))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(911))
                   datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
                 endif
               ENDIF
              ENDIF
            ENDIF

!
!           FD LEVEL SPEC HUMIDITY.
            iget1 = IGET(451)
            iget2 = iget(578)
            if (iget1 > 0) then
              work1 = LVLS(IFD,iget1)
            else
              work1 = 0.0
            endif
            if (iget2 > 0) then
              work2 = LVLS(IFD,iget2)
            else
              work2 = 0.0
            endif
            IF (IGET1 > 0 .or. iget2 > 0) THEN
              IF (work1 > 0 .or. work2 > 0)THEN
!$omp parallel do private(i,j)
                DO J=JSTA,JEND
                  DO I=1,IM
                    GRID1(I,J) = Q7D(I,J,IFD)
                  ENDDO
                ENDDO
                if(work1 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET1)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET1)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                endif
                if(work2 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET2)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET2)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                endif
              ENDIF
            ENDIF
!
!           FD LEVEL PRESSURE
            iget1 = IGET(482)
            iget2 = iget(579)
            if (iget1 > 0) then
              work1 = LVLS(IFD,iget1)
            else
              work1 = 0.0
            endif
            if (iget2 > 0) then
              work2 = LVLS(IFD,iget2)
            else
              work2 = 0.0
            endif
            IF (IGET1 > 0 .or. IGET2 > 0) THEN
              IF (work1 > 0 .or. work2 > 0) THEN
!$omp parallel do private(i,j)
                DO J=JSTA,JEND
                  DO I=1,IM
                    GRID1(I,J) = P7D(I,J,IFD)
                  ENDDO
                ENDDO
                if(work1 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET1)
                   fld_info(cfld)%lvl   = LVLSXML(IFD,IGET1)
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                     jj = jsta+j-1
                     do i=1,im
                       datapd(i,j,cfld) = GRID1(i,jj)
                     enddo
                   enddo
                  endif
                endif
                if(work2 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET2)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET2)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                endif
              ENDIF
            ENDIF
!
!           FD LEVEL ICING
            iget1 = IGET(580)
            iget2 = iget(587)
            if (iget1 > 0) then
              work1 = LVLS(IFD,iget1)
            else
              work1 = 0.0
            endif
            if (iget2 > 0) then
              work2 = LVLS(IFD,iget2)
            else
              work2 = 0.0
            endif
            IF (IGET1 > 0 .or. IGET2 > 0) THEN
              IF (work1 > 0 .or. work2 > 0) THEN
!$omp parallel do private(i,j)
                DO J=JSTA,JEND
                  DO I=1,IM
                    GRID1(I,J) = ICINGFD(I,J,IFD)
                  ENDDO
                ENDDO
                if(iget1 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET1)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET1)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                endif
                if(work2 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET2)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET2)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                endif

              ENDIF
            ENDIF
!
!  ADD FD LEVEL DUST/ASH (GOCART)
            IF (IGET(601)>0) THEN                      ! DUST 1
	      IF (LVLS(IFD,IGET(601))>0) THEN
!$omp parallel do private(i,j)
	       DO J=JSTA,JEND
	       DO I=1,IM
	          GRID1(I,J)=AERFD(I,J,IFD,1) 
               ENDDO
               ENDDO
               if(iget(601)>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(601))
                  fld_info(cfld)%lvl=LVLSXML(IFD,IGET(601))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               endif
	      ENDIF
	    ENDIF

            IF (IGET(602)>0) THEN			! DUST 2
	      IF (LVLS(IFD,IGET(602))>0) THEN
!$omp parallel do private(i,j)
	       DO J=JSTA,JEND
	       DO I=1,IM
	          GRID1(I,J)=AERFD(I,J,IFD,2) 
               ENDDO
               ENDDO
               if(iget(602)>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(602))
                  fld_info(cfld)%lvl=LVLSXML(IFD,IGET(602))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               endif
	      ENDIF
	    ENDIF

            IF (IGET(603)>0) THEN			! DUST 3
	      IF (LVLS(IFD,IGET(603))>0) THEN
!$omp parallel do private(i,j)
	       DO J=JSTA,JEND
	       DO I=1,IM
	          GRID1(I,J)=AERFD(I,J,IFD,3) 
               ENDDO
               ENDDO
               if(iget(603)>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(603))
                  fld_info(cfld)%lvl=LVLSXML(IFD,IGET(603))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               endif
	      ENDIF
	    ENDIF

            IF (IGET(604)>0) THEN			! DUST 4
	      IF (LVLS(IFD,IGET(604))>0) THEN
!$omp parallel do private(i,j)
	       DO J=JSTA,JEND
	       DO I=1,IM
	          GRID1(I,J)=AERFD(I,J,IFD,4) 
               ENDDO
               ENDDO
               if(iget(604)>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(604))
                  fld_info(cfld)%lvl=LVLSXML(IFD,IGET(604))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               endif
	      ENDIF
	    ENDIF

            IF (IGET(605)>0) THEN			! DUST 5
	      IF (LVLS(IFD,IGET(605))>0) THEN
!$omp parallel do private(i,j)
	       DO J=JSTA,JEND
	       DO I=1,IM
	          GRID1(I,J)=AERFD(I,J,IFD,5) 
               ENDDO
               ENDDO
               if(iget(605)>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(605))
                  fld_info(cfld)%lvl=LVLSXML(IFD,IGET(605))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               endif
	      ENDIF
	    ENDIF

!
!
!           FD LEVEL U WIND AND/OR V WIND.
            IF ((IGET(060)>0).OR.(IGET(061)>0)) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=U7D(I,J,IFD)
                 GRID2(I,J)=V6D(I,J,IFD)
               ENDDO
               ENDDO
               IF (IGET(060)>0) THEN
                 IF (LVLS(IFD,IGET(060))>0) then
                  if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(060))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(060))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                  endif
                 ENDIF
               ENDIF
               IF (IGET(061)>0) THEN
                 IF (LVLS(IFD,IGET(061))>0) THEN
                  if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(061))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(061))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID2(i,jj)
                    enddo
                  enddo
                  endif
                 ENDIF
               ENDIF
            ENDIF
!
!           FD LEVEL U WIND AND/OR V WIND.
            IF ((IGET(576)>0).OR.(IGET(577)>0)) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = U7D(I,J,IFD)
                   GRID2(I,J) = V6D(I,J,IFD)
                 ENDDO
               ENDDO
               IF (IGET(576)>0) THEN
                 IF (LVLS(IFD,IGET(576))>0) then
                  if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(576))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(576))
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                     jj = jsta+j-1
                     do i=1,im
                       datapd(i,j,cfld) = GRID1(i,jj)
                     enddo
                   enddo
                  endif
                 ENDIF
               ENDIF
               IF (IGET(577)>0) THEN
                 IF (LVLS(IFD,IGET(577))>0) THEN
                  if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(577))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(577))
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                     jj = jsta+j-1
                     do i=1,im
                       datapd(i,j,cfld) = GRID2(i,jj)
                     enddo
                   enddo
                  endif
                 ENDIF
               ENDIF
            ENDIF

 10      CONTINUE
         DEALLOCATE(T7D,Q7D,U7D,V6D,P7D,ICINGFD,AERFD)
      ENDIF

!
!     ***BLOCK 3-2:  FD LEVEL (from control file) GTG
!     
      IF(IGET(467)>0.or.IGET(468)>0.or.IGET(469)>0) THEN
         if(IGET(467)>0) THEN          ! GTG
            N=IAVBLFLD(IGET(467))
            NFDCTL=size(pset%param(N)%level)
            if(allocated(ITYPEFDLVLCTL)) deallocate(ITYPEFDLVLCTL)
            allocate(ITYPEFDLVLCTL(NFDCTL))
            DO IFD = 1,NFDCTL
               ITYPEFDLVLCTL(IFD)=LVLS(IFD,IGET(467))
            enddo
            if(allocated(HTFDCTL)) deallocate(HTFDCTL)
            allocate(HTFDCTL(NFDCTL))
            HTFDCTL=pset%param(N)%level
!           print *, "GTG 467 levels=",pset%param(N)%level
            allocate(GTGFD(IM,JSTA:JEND,NFDCTL))
            call FDLVL_MASS(ITYPEFDLVLCTL,NFDCTL,HTFDCTL,GTG,GTGFD)
!           print *, "GTG 467 Done GTGFD=",me,GTGFD(IM/2,jend,1:NFDCTL)
            DO IFD = 1,NFDCTL
              IF (LVLS(IFD,IGET(467))>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                 DO I=1,IM
                    GRID1(I,J)=GTGFD(I,J,IFD)
                 ENDDO
                 ENDDO
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(467))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(467))
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                         datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                   enddo
                 endif
              ENDIF
            ENDDO
         endif

         if(IGET(468)>0) THEN          ! CAT
            N=IAVBLFLD(IGET(468))
            NFDCTL=size(pset%param(N)%level)
            if(allocated(ITYPEFDLVLCTL)) deallocate(ITYPEFDLVLCTL)
            allocate(ITYPEFDLVLCTL(NFDCTL))
            DO IFD = 1,NFDCTL
               ITYPEFDLVLCTL(IFD)=LVLS(IFD,IGET(468))
            enddo
            if(allocated(HTFDCTL)) deallocate(HTFDCTL)
            allocate(HTFDCTL(NFDCTL))
            HTFDCTL=pset%param(N)%level
            allocate(CATFD(IM,JSTA:JEND,NFDCTL))
            call FDLVL_MASS(ITYPEFDLVLCTL,NFDCTL,HTFDCTL,catedr,CATFD)
            DO IFD = 1,NFDCTL
              IF (LVLS(IFD,IGET(468))>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                 DO I=1,IM
                    GRID1(I,J)=CATFD(I,J,IFD)
                 ENDDO
                 ENDDO
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(468))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(468))
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                         datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                   enddo
                 endif
              ENDIF
            ENDDO
         endif

         if(IGET(469)>0) THEN          ! MWT
            N=IAVBLFLD(IGET(469))
            NFDCTL=size(pset%param(N)%level)
            if(allocated(ITYPEFDLVLCTL)) deallocate(ITYPEFDLVLCTL)
            allocate(ITYPEFDLVLCTL(NFDCTL))
            DO IFD = 1,NFDCTL
               ITYPEFDLVLCTL(IFD)=LVLS(IFD,IGET(469))
            enddo
            if(allocated(HTFDCTL)) deallocate(HTFDCTL)
            allocate(HTFDCTL(NFDCTL))
            HTFDCTL=pset%param(N)%level
            allocate(MWTFD(IM,JSTA:JEND,NFDCTL))
            call FDLVL_MASS(ITYPEFDLVLCTL,NFDCTL,HTFDCTL,MWT,MWTFD)
            DO IFD = 1,NFDCTL
              IF (LVLS(IFD,IGET(469))>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                 DO I=1,IM
                    GRID1(I,J)=MWTFD(I,J,IFD)
                 ENDDO
                 ENDDO
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(469))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(469))
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                         datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                   enddo
                 endif
              ENDIF
            ENDDO
         endif

         if(allocated(GTGFD)) deallocate(GTGFD)
         if(allocated(CATFD)) deallocate(CATFD)
         if(allocated(MWTFD)) deallocate(MWTFD)

         if(allocated(ITYPEFDLVLCTL)) deallocate(ITYPEFDLVLCTL)
         if(allocated(HTFDCTL)) deallocate(HTFDCTL)

      ENDIF
!     
!
!     ***BLOCK 4:  FREEZING LEVEL Z, RH AND P.
!     
      IF ( (IGET(062)>0).OR.(IGET(063)>0) ) THEN
         CALL FRZLVL(Z1D,RH1D,P1D)
!
!        FREEZING LEVEL HEIGHT.
         IF (IGET(062)>0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J)=Z1D(I,J)
                IF (SUBMODELNAME == 'RTMA' .OR. MODELNAME == 'FV3R') THEN
                  FREEZELVL(I,J)=GRID1(I,J)
                ENDIF
              ENDDO
            ENDDO
            CALL BOUND (GRID1,D00,H99999)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(062))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
         ENDIF
!
!        FREEZING LEVEL RELATIVE HUMIDITY.
         IF (IGET(063)>0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J) = RH1D(I,J)
               ENDDO
            ENDDO
            CALL SCLFLD(GRID1,H100,IM,JM)
            CALL BOUND(GRID1,H1,H100)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(063))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
         ENDIF
!
!        FREEZING LEVEL PRESSURE
         IF (IGET(753)>0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J) = P1D(I,J)
              ENDDO
            ENDDO
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(753))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif

         ENDIF
      ENDIF

      IF (IGET(165)>0 .OR. IGET(350)>0.OR. IGET(756)>0) THEN
         CALL FRZLVL2(TFRZ,Z1D,RH1D,P1D)
!
!        HIGHEST FREEZING LEVEL HEIGHT.
          IF (IGET(165)>0)THEN  
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=Z1D(I,J)
               ENDDO
               ENDDO
            CALL BOUND (GRID1,D00,H99999)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(165))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST FREEZING LEVEL RELATIVE HUMIDITY
          IF (IGET(350)>0)THEN  
               GRID1=spval 
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 IF(RH1D(I,J) < spval) GRID1(I,J)=RH1D(I,J)*100.
               ENDDO
               ENDDO
            CALL BOUND (GRID1,H1,H100)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(350))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST FREEZING LEVEL PRESSURE
         IF (IGET(756)>0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J) = P1D(I,J)
              ENDDO
            ENDDO
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(756))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
         ENDIF

      ENDIF
!
! HIGHEST -10 C ISOTHERM VALUES
!
      IF (IGET(776)>0 .OR. IGET(777)>0.OR. IGET(778)>0) THEN 
         CALL FRZLVL2(263.15,Z1D,RH1D,P1D)
!
!        HIGHEST -10C ISOTHERM LEVEL HEIGHT.
          IF (IGET(776)>0)THEN  
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=Z1D(I,J)
               ENDDO
               ENDDO
            CALL BOUND (GRID1,D00,H99999)
            if(grib=='grib2') then 
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(776))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST -10C ISOTHERM RELATIVE HUMIDITY
          IF (IGET(777)>0)THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 IF(RH1D(I,J) < spval) GRID1(I,J)=RH1D(I,J)*100.
               ENDDO
               ENDDO
            CALL BOUND (GRID1,H1,H100)
            if(grib=='grib2') then 
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(777))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST -10C ISOTHERM LEVEL PRESSURE
         IF (IGET(778)>0) THEN 
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=P1D(I,J)
               ENDDO
               ENDDO
            if(grib=='grib2') then 
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(778))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
         ENDIF

      ENDIF
!
! HIGHEST -20 C ISOTHERM VALUES
!
      IF (IGET(779)>0 .OR. IGET(780)>0.OR. IGET(781)>0) THEN
         CALL FRZLVL2(253.15,Z1D,RH1D,P1D)
!
!        HIGHEST -20C ISOTHERM LEVEL HEIGHT.
          IF (IGET(779)>0)THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=Z1D(I,J)
               ENDDO
               ENDDO
            CALL BOUND (GRID1,D00,H99999)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(779))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST -20C ISOTHERM RELATIVE HUMIDITY
          IF (IGET(780)>0)THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 IF(RH1D(I,J) < spval) GRID1(I,J)=RH1D(I,J)*100.
               ENDDO
               ENDDO
            CALL BOUND (GRID1,H1,H100)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(780))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST -20C ISOTHERM LEVEL PRESSURE
         IF (IGET(781)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=P1D(I,J)
               ENDDO
               ENDDO
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(781))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
         ENDIF

      ENDIF
!     
      allocate(PBND(IM,jsta:jend,NBND), TBND(IM,jsta:jend,NBND),    &
               QBND(IM,jsta:jend,NBND), UBND(IM,jsta:jend,NBND),    &
               VBND(IM,jsta:jend,NBND), RHBND(IM,jsta:jend,NBND),   &
               WBND(IM,jsta:jend,NBND))

!
!     ***BLOCK 5:  BOUNDARY LAYER FIELDS.
!     
      IF ( (IGET(067)>0).OR.(IGET(068)>0).OR.       &
           (IGET(069)>0).OR.(IGET(070)>0).OR.       &
           (IGET(071)>0).OR.(IGET(072)>0).OR.       &
           (IGET(073)>0).OR.(IGET(074)>0).OR.       &
           (IGET(088)>0).OR.(IGET(089)>0).OR.       &
           (IGET(090)>0).OR.(IGET(075)>0).OR.       &
           (IGET(109)>0).OR.(IGET(110)>0).OR.       &
           (IGET(031)>0).OR.(IGET(032)>0).OR.       &
           (IGET(573)>0).OR.                           &
           (IGET(107)>0).OR.(IGET(091)>0).OR.       &
           (IGET(092)>0).OR.(IGET(093)>0).OR.       &
           (IGET(094)>0).OR.(IGET(095)>0).OR.       &
           (IGET(096)>0).OR.(IGET(097)>0).OR.       &
           (IGET(098)>0).OR.(IGET(221)>0) ) THEN
!
           allocate(OMGBND(IM,jsta:jend,NBND),PWTBND(IM,jsta:jend,NBND),  &
                    QCNVBND(IM,jsta:jend,NBND),LVLBND(IM,jsta:jend,NBND), &
                    LB2(IM,jsta:jend))

!        COMPUTE ETA BOUNDARY LAYER FIELDS.
         CALL BNDLYR(PBND,TBND,QBND,RHBND,UBND,VBND,      &
                     WBND,OMGBND,PWTBND,QCNVBND,LVLBND)

!$omp parallel do private(i,j)
         DO J=JSTA,JEND
           DO I=1,IM
             EGRID2(i,j) = SPVAL     
           ENDDO
         ENDDO

!     
!        LOOP OVER NBND BOUNDARY LAYERS.
         DO 20 LBND = 1,NBND
!     
!           BOUNDARY LAYER PRESSURE.
            IF (IGET(067)>0) THEN
              IF (LVLS(LBND,IGET(067))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = PBND(I,J,LBND)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(067))
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(067))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
                endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER TEMPERATURE.
            IF (IGET(068)>0) THEN
              IF (LVLS(LBND,IGET(068))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=TBND(I,J,LBND)
               ENDDO
               ENDDO
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(068))
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(068))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER POTENTIAL TEMPERATURE.
            IF (IGET(069)>0) THEN
              IF (LVLS(LBND,IGET(069))>0) THEN
               CALL CALPOT(PBND(1,jsta,LBND),TBND(1,jsta,LBND),GRID1(1,jsta))
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(069))
                 fld_info(cfld)%lvl=LVLSXML(IFD,IGET(069))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
                endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER RELATIVE HUMIDITY.
            IF (IGET(072)>0) THEN
              IF (LVLS(LBND,IGET(072))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RHBND(I,J,LBND)
               ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(072))
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(072))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER DEWPOINT TEMPERATURE.
            IF (IGET(070)>0) THEN
              IF (LVLS(LBND,IGET(070))>0) THEN
               CALL CALDWP(PBND(1,jsta,LBND), QBND(1,jsta,LBND),     &
                           GRID1(1,jsta),     TBND(1,jsta,LBND))
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(070))
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(070))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER SPECIFIC HUMIDITY.
            IF (IGET(071)>0) THEN
              IF (LVLS(LBND,IGET(071))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=QBND(I,J,LBND)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,H1M12,H99999)
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(071))
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(071))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER MOISTURE CONVERGENCE.
            IF (IGET(088)>0) THEN
              IF (LVLS(LBND,IGET(088))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = QCNVBND(I,J,LBND)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(088))
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(088))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER U WIND AND/OR V WIND.
!
            FIELD1=.FALSE.
            FIELD2=.FALSE.
!
            IF(IGET(073)>0)THEN
              IF(LVLS(LBND,IGET(073))>0)FIELD1=.TRUE.
            ENDIF
            IF(IGET(074)>0)THEN
              IF(LVLS(LBND,IGET(074))>0)FIELD2=.TRUE.
            ENDIF
!
            IF(FIELD1.OR.FIELD2)THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = UBND(I,J,LBND)
                   GRID2(I,J) = VBND(I,J,LBND)
                 ENDDO
               ENDDO
!
               IF (IGET(073)>0) THEN
                 IF (LVLS(LBND,IGET(073))>0) then
                   if(grib=='grib2') then
                    cfld=cfld+1
                    fld_info(cfld)%ifld=IAVBLFLD(IGET(073))
                    fld_info(cfld)%lvl=LVLSXML(LBND,IGET(073))
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                   endif
                 ENDIF
               ENDIF
               IF (IGET(074)>0) THEN
                 IF (LVLS(LBND,IGET(074))>0) THEN
                   if(grib=='grib2') then
                    cfld=cfld+1
                    fld_info(cfld)%ifld=IAVBLFLD(IGET(074))
                    fld_info(cfld)%lvl=LVLSXML(LBND,IGET(074))
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID2(i,jj)
                      enddo
                    enddo
                   endif
                 ENDIF
               ENDIF
            ENDIF
!     
!           BOUNDARY LAYER OMEGA.
            IF (IGET(090)>0) THEN
              IF (LVLS(LBND,IGET(090))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = OMGBND(I,J,LBND)
                 ENDDO
               ENDDO
              if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(090))
               fld_info(cfld)%lvl=LVLSXML(LBND,IGET(090))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
              endiF
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER PRECIPITBLE WATER.
            IF (IGET(089)>0) THEN
              IF (LVLS(LBND,IGET(089))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = PWTBND(I,J,LBND)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
              if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(089))
               fld_info(cfld)%lvl=LVLSXML(LBND,IGET(089))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
              endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER LIFTED INDEX.
            IF (IGET(075)>0 .OR. IGET(031)>0 .OR. IGET(573)>0) THEN
             CALL OTLFT(PBND(1,jsta,LBND),TBND(1,jsta,LBND),    &
                    QBND(1,jsta,LBND),GRID1(1,jsta))
             IF(IGET(075)>0)THEN
              IF (LVLS(LBND,IGET(075))>0) THEN
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(075))
                fld_info(cfld)%lvl=LVLSXML(LBND,IGET(075))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
              ENDIF
             END IF
             IF(IGET(031)>0 .or. IGET(573)>0)THEN
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                DO I=1,IM
                  EGRID2(I,J) = MIN(EGRID2(I,J),GRID1(I,J))
                END DO
              END DO
             END IF
            ENDIF
!
!        END OF ETA BOUNDARY LAYER LOOP.
 20      CONTINUE
         deallocate(OMGBND,PWTBND,QCNVBND)
!     
!        BEST LIFTED INDEX FROM BOUNDARY LAYER FIELDS.
!     
         IF (IGET(031)>0 .OR. IGET(573)>0 ) THEN
!           DO J=JSTA,JEND
!            DO I=1,IM
!              EGRID1(I,J) = H99999
!              EGRID2(I,J) = H99999
!            ENDDO
!            ENDDO
!
!            DO 50 LBND = 1,NBND
!               CALL OTLFT(PBND(1,1,LBND),TBND(1,1,LBND),      &
!                    QBND(1,1,LBND),EGRID2)
!               DO J=JSTA,JEND
!               DO I=1,IM
!                 EGRID1(I,J)=AMIN1(EGRID1(I,J),EGRID2(I,J))
!               ENDDO
!               ENDDO
! 50         CONTINUE
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J)=EGRID2(I,J)
              ENDDO
            ENDDO
!	    print*,'writting out best lifted index'

            if (IGET(031)>0) then
              if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(031))
               datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
              endif
            endif
   
            if(IGET(573)> 0 ) THEN
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(573))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            endif

	 END IF
!     
!        BEST BOUNDARY LAYER CAPE AND CINS.
!     
         FIELD1=.FALSE.
         FIELD2=.FALSE.
!
         IF(IGET(032)>0)THEN
           IF(LVLS(2,IGET(032))>0)FIELD1=.TRUE.
         ENDIF
         IF(IGET(107)>0)THEN
           IF(LVLS(2,IGET(107))>0)FIELD2=.TRUE.
         ENDIF
!
         IF(IGET(566)>0)THEN
           FIELD1=.TRUE.
         ENDIF
         IF(IGET(567)>0)THEN
           FIELD2=.TRUE.
         ENDIF
!
         !if(grib=="grib2") print *,'in MISCLN.f,iget(566)=',          &
         !  iget(566), 'iget(567)=',iget(567),'LVLSXML(1,IGET(566)=',  &
         !  LVLSXML(1,IGET(566)),'LVLSXML(1,IGET(567)=',               &
         !  LVLSXML(1,IGET(567)),'field1=',field1,'field2=',field2
!
         IF(FIELD1.OR.FIELD2)THEN
           ITYPE = 2
!
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               EGRID1(I,J) = -H99999
               EGRID2(I,J) = -H99999
             ENDDO
           ENDDO
!
           DO 80 LBND = 1,NBND
           CALL CALTHTE(PBND(1,jsta,LBND),TBND(1,jsta,LBND),        &
                        QBND(1,jsta,LBND),EGRID1)
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               IF (EGRID1(I,J) > EGRID2(I,J)) THEN
                 EGRID2(I,J) = EGRID1(I,J)
                 LB2(I,J)    = LVLBND(I,J,LBND)
                 P1D(I,J)    = PBND(I,J,LBND)
                 T1D(I,J)    = TBND(I,J,LBND)
                 Q1D(I,J)    = QBND(I,J,LBND)
               ENDIF
             ENDDO
           ENDDO
 80        CONTINUE
!
           DPBND = 0.
           CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,LB2,EGRID1,   &
                        EGRID2,EGRID3,EGRID4,EGRID5) 
!
           IF (IGET(566)>0) THEN
! dong add missing value for cape
              GRID1=spval
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                DO I=1,IM
                  IF(T1D(I,J) < spval) GRID1(I,J) = EGRID1(I,J)
                ENDDO
              ENDDO
             CALL BOUND(GRID1,D00,H99999)
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(566))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(566))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
           ENDIF
!
           IF (IGET(567) > 0) THEN
! dong add missing value for cape
              GRID1=spval
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(T1D(I,J) < spval) GRID1(I,J) = - EGRID2(I,J)
               ENDDO
             ENDDO
!
             CALL BOUND(GRID1,D00,H99999)
!
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(T1D(I,J) < spval) GRID1(I,J) = - GRID1(I,J)
               ENDDO
             ENDDO
!
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(567))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(567))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
           ENDIF
         ENDIF
!

!    PBL HEIGHT 
         IF(IGET(221) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = PBLH(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(221))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         END IF
!        BOUNDARY LAYER LIFTING CONDENSATION PRESSURE AND HEIGHT.
!        EGRID1 IS LCL PRESSURE.  EGRID2 IS LCL HEIGHT.
!
         IF ( (IGET(109)>0).OR.(IGET(110)>0) ) THEN
            CALL CALLCL(PBND(1,jsta,1),TBND(1,jsta,1),          &
                        QBND(1,jsta,1),EGRID1,EGRID2)
            IF (IGET(109)>0) THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(TBND(I,J,1) < spval) GRID1(I,J) = EGRID2(I,J)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(109))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
            IF (IGET(110)>0) THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(TBND(I,J,1) < spval) GRID1(I,J) = EGRID1(I,J)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(110))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
         ENDIF
!     
!        NGM BOUNDARY LAYER FIELDS.
!     
         IF ( (IGET(091)>0).OR.(IGET(092)>0).OR.      &
              (IGET(093)>0).OR.(IGET(094)>0).OR.      &
              (IGET(095)>0).OR.(IGET(095)>0).OR.      &
              (IGET(096)>0).OR.(IGET(097)>0).OR.      &
              (IGET(098)>0) ) THEN

              allocate(T78483(im,jsta:jend), T89671(im,jsta:jend), &
                       P78483(im,jsta:jend), P89671(im,jsta:jend))
!
!  COMPUTE SIGMA 0.89671 AND 0.78483 TEMPERATURES
!    INTERPOLATE LINEAR IN LOG P
            IF (IGET(097)>0.OR.IGET(098)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   P78483(I,J) = LOG(PINT(I,J,NINT(LMH(I,J)))*0.78483)
                   P89671(I,J) = LOG(PINT(I,J,NINT(LMH(I,J)))*0.89671)
                 ENDDO
               ENDDO
               DONE =.FALSE.
               DONE1=.FALSE.
!!$omp  parallel do private(fac1,fac2,pkl1,pku1,t78483,t89671)
               DO L=2,LM
                 DO J=JSTA,JEND
                   DO I=1,IM                  
                     PKL1=0.5*(ALPINT(I,J,L)+ALPINT(I,J,L+1))
                     PKU1=0.5*(ALPINT(I,J,L)+ALPINT(I,J,L-1))
!                     IF(I==1 .AND. J==1)PRINT*,'L,P89671,PKL1,PKU1= ', &
!                                             L,P89671(I,J), PKL1, PKU1
                     IF(P78483(I,J) < PKL1.AND.P78483(I,J) > PKU1) THEN
                      FAC1 = (PKL1-P78483(I,J))/(PKL1-PKU1)
                      FAC2 = (P78483(I,J)-PKU1)/(PKL1-PKU1)
                      T78483(I,J) = T(I,J,L)*FAC2 + T(I,J,L-1)*FAC1
                      DONE1(I,J)=.TRUE.
                     ENDIF
                     IF(P89671(I,J) < PKL1.AND.P89671(I,J) > PKU1)THEN
                       FAC1 = (PKL1-P89671(I,J))/(PKL1-PKU1)
                       FAC2 = (P89671(I,J)-PKU1)/(PKL1-PKU1)
                       T89671(I,J) = T(I,J,L)*FAC2 + T(I,J,L-1)*FAC1
                       DONE(I,J)   = .TRUE.
!	  	       IF(I==1 .AND. J==1)PRINT*,'done(1,1)= ',done(1,1)
                     ENDIF
                   ENDDO
                 ENDDO
               ENDDO
!	       print*,'done(1,1)= ',done(1,1)
!$omp parallel do private(i,j,pl,tl,ql,qsat,rhl)
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(.NOT. DONE(I,J)) THEN
                   PL    = PINT(I,J,LM-1)
                   TL   = 0.5*(T(I,J,LM-2)+T(I,J,LM-1))
                   QL   = 0.5*(Q(I,J,LM-2)+Q(I,J,LM-1))
                   QSAT = PQ0/PL *EXP(A2*(TL-A3)/(TL-A4))
!
                   RHL = QL/QSAT
!
                   IF(RHL > 1.)THEN
                    RHL = 1.
                    QL  = RHL*QSAT
                   ENDIF
!
                   IF(RHL < RHmin)THEN
                    RHL = RHmin
                    QL  = RHL*QSAT
                   ENDIF
!
!                  TVRL   = TL*(1.+0.608*QL)
!                  TVRBLO = TVRL*(P89671(I,J)/PL)**RGAMOG
!                  T89671(I,J) = TVRBLO/(1.+0.608*QL)

                   T89671(I,J) = TL * (P89671(I,J)/PL)**RGAMOG
!     

!                   PKL1=0.5*(ALPINT(I,J,LM)+ALPINT(I,J,LM+1))
!                   PKU1=0.5*(ALPINT(I,J,LM-1)+ALPINT(I,J,LM))
!                   T89671(I,J)=T(I,J,LM)+(T(I,J,LM)-T(I,J,LM-1))*
!     +               (P89671(I,J)-PKL1)/(PKL1-PKU1)

!                   print*,'Debug T89671= ',i,j
!     +		     ,P89671(I,J), PKL1, PKU1  
!     +               ,T89671(I,J),T(I,J,LM-1),T(I,J,LM)
                 END IF
 
                 IF(.NOT. DONE1(I,J))THEN
                   PL   = PINT(I,J,LM-1)
                   TL   = 0.5*(T(I,J,LM-2)+T(I,J,LM-1))
                   QL   = 0.5*(Q(I,J,LM-2)+Q(I,J,LM-1))
                   QSAT = PQ0/PL *EXP(A2*(TL-A3)/(TL-A4))
!
                   RHL = QL/QSAT
!
                   IF(RHL>1.)THEN
                    RHL = 1.
                    QL  = RHL*QSAT
                   ENDIF
!
                   IF(RHL<RHmin)THEN
                     RHL = RHmin
                     QL  = RHL*QSAT
                   ENDIF
!
!                  TVRL   = TL*(1.+0.608*QL)
!                  TVRBLO = TVRL*(P78483(I,J)/PL)**RGAMOG
!                  T78483(I,J)  =TVRBLO/(1.+0.608*QL)

                   T78483(I,J)  = TL * (P78483(I,J)/PL)**RGAMOG
!     
                 END IF
 
               END DO
             END DO
!     
!           SIGMA 0.89671 TEMPERATURE
             IF (IGET(097) > 0) THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T(I,J,LM) < spval) GRID1(I,J) = T89671(I,J)
!                  IF(T89671(I,J)>350.)PRINT*,'LARGE T89671 ',   &
!                    I,J,T89671(I,J)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(097))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(097))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
             ENDIF
!     
!           SIGMA 0.78483 TEMPERATURE
             IF (IGET(098)>0) THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T(I,J,LM) < spval) GRID1(I,J) = T78483(I,J)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(098))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(098))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
             ENDIF
             deallocate(T78483, T89671, P78483, P89671)
            ENDIF
!     
!           NGM SIGMA LAYER 0.98230 FIELDS.  THESE FIELDS ARE 
!           THE FIRST ETA LAYER BOUNDARY LAYER FIELDS. 
!     
!     
            IF ( (IGET(091)>0).OR.(IGET(092)>0).OR.      &
                 (IGET(093)>0).OR.(IGET(094)>0).OR.      &
                 (IGET(095)>0).OR.(IGET(095)>0).OR.      &
                 (IGET(096)>0) ) THEN
!     
!     
!              PRESSURE.
               IF (IGET(091)>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                   DO I=1,IM
                     GRID1(I,J) = PBND(I,J,1)
                   ENDDO
                 ENDDO
                 if(grib=='grib2') then
                  cfld=cfld+1
                  fld_info(cfld)%ifld=IAVBLFLD(IGET(091))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               ENDIF
!     
!              TEMPERATURE.
               IF (IGET(092)>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                   DO I=1,IM
                     GRID1(I,J) = TBND(I,J,1)
                   ENDDO
                 ENDDO
                 if(grib=='grib2') then
                  cfld=cfld+1
                  fld_info(cfld)%ifld=IAVBLFLD(IGET(092))
                  fld_info(cfld)%lvl=LVLSXML(1,IGET(092))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               ENDIF
!     
!              SPECIFIC HUMIDITY.
               IF (IGET(093)>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                   DO I=1,IM
                     GRID1(I,J) = QBND(I,J,1)
                   ENDDO
                 ENDDO
                  CALL BOUND(GRID1,H1M12,H99999)
                 if(grib=='grib2') then
                  cfld=cfld+1
                  fld_info(cfld)%ifld=IAVBLFLD(IGET(093))
                  fld_info(cfld)%lvl=LVLSXML(1,IGET(093))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               ENDIF
!     
!              RELATIVE HUMIDITY.
               IF (IGET(094)>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                   DO I=1,IM
                     GRID1(I,J) = RHBND(I,J,1)
                   ENDDO
                 ENDDO
                  CALL SCLFLD(GRID1,H100,IM,JM)
                  CALL BOUND(GRID1,H1,H100)
                if(grib=='grib2') then
                  cfld=cfld+1
                  fld_info(cfld)%ifld=IAVBLFLD(IGET(094))
                  fld_info(cfld)%lvl=LVLSXML(1,IGET(094))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
               endif
               ENDIF
!     
!              U AND/OR V WIND.
               IF ((IGET(095)>0).OR.(IGET(096)>0)) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                   DO I=1,IM
                     GRID1(I,J) = UBND(I,J,1)
                     GRID2(I,J) = VBND(I,J,1)
                   ENDDO
                 ENDDO
                  IF (IGET(095)>0)  then
                    if(grib=='grib2') then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(IGET(095))
                     fld_info(cfld)%lvl=LVLSXML(1,IGET(095))
!$omp parallel do private(i,j,jj)
                     do j=1,jend-jsta+1
                       jj = jsta+j-1
                       do i=1,im
                         datapd(i,j,cfld) = GRID1(i,jj)
                       enddo
                     enddo
                    endif
                  ENDIF
                  IF (IGET(096)>0) then
                    if(grib=='grib2') then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(IGET(096))
                     fld_info(cfld)%lvl=LVLSXML(1,IGET(096))
!$omp parallel do private(i,j,jj)
                     do j=1,jend-jsta+1
                       jj = jsta+j-1
                       do i=1,im
                         datapd(i,j,cfld) = GRID2(i,jj)
                       enddo
                     enddo
                    endif
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!     
!     ENDIF FOR BOUNDARY LAYER BLOCK.
!
      ENDIF
!     
!
!     ***BLOCK 6:  MISCELLANEOUS LAYER MEAN LFM AND NGM FIELDS.
!     
      IF ( (IGET(066)>0).OR.(IGET(081)>0).OR.        &
           (IGET(082)>0).OR.(IGET(104)>0).OR.        &
           (IGET(099)>0).OR.(IGET(100)>0).OR.        &
           (IGET(101)>0).OR.(IGET(102)>0).OR.        &
           (IGET(103)>0) ) THEN
!     
!        LFM "MEAN" RELATIVE HUMIDITIES AND PRECIPITABLE WATER.
!     
         IF ( (IGET(066)>0).OR.(IGET(081)>0).OR.     &
              (IGET(082)>0).OR.(IGET(104)>0) ) THEN
            allocate(RH3310(IM,jsta:jend),RH6610(IM,jsta:jend),          &
                     RH3366(IM,jsta:jend),PW3310(IM,jsta:jend))
            CALL LFMFLD(RH3310,RH6610,RH3366,PW3310)
!     
!           SIGMA 0.33-1.00 MEAN RELATIVE HUMIIDITY.
            IF (IGET(066)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH3310(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
              if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(066))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(066))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
!               print *,'in miscln,RH0.33-1.0,cfld=',cfld,'fld=',  &
!                IAVBLFLD(IGET(066)),'lvl=',fld_info(cfld)%lvl
              endif
            ENDIF
!     
!           SIGMA 0.66-1.00 MEAN RELATIVE HUMIIDITY.
            IF (IGET(081)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH6610(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(081))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(081))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.33-0.66 MEAN RELATIVE HUMIIDITY.
            IF (IGET(082)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH3366(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(082))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(082))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.33-1.00 PRECIPITABLE WATER.
            IF (IGET(104)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = PW3310(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(104))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(104))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
            deallocate(RH3310,RH6610,RH3366,PW3310)
         ENDIF
!     
!        VARIOUS LAYER MEAN NGM SIGMA FIELDS.
!     
         IF ( (IGET(099)>0).OR.(IGET(100)>0).OR.    &
              (IGET(101)>0).OR.(IGET(102)>0).OR.    &
              (IGET(103)>0) ) THEN
            allocate(RH4710(IM,jsta_2l:jend_2u),RH4796(IM,jsta_2l:jend_2u), &
                     RH1847(IM,jsta_2l:jend_2u))
            allocate(RH8498(IM,jsta_2l:jend_2u),QM8510(IM,jsta_2l:jend_2u))

            CALL NGMFLD(RH4710,RH4796,RH1847,RH8498,QM8510)
!     
!           SIGMA 0.47191-1.00000 RELATIVE HUMIDITY.
            IF (IGET(099)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH4710(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(099))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(099))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.47191-0.96470 RELATIVE HUMIDITY.
            IF (IGET(100)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH4796(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(100))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(100))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.18019-0.47191 RELATIVE HUMIDITY.
            IF (IGET(101)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH1847(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(101))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(101))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.84368-0.98230 RELATIVE HUMIDITY.
            IF (IGET(102)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH8498(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
              if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(102))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(102))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
              endif
            ENDIF
!     
!           SIGMA 0.85000-1.00000 MOISTURE CONVERGENCE.
            IF (IGET(103)>0) THEN
            GRID1=spval
!           CONVERT TO DIVERGENCE FOR GRIB
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(QM8510(I,J) < spval) GRID1(I,J) = -1.0*QM8510(I,J)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(103))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(103))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
            deallocate(RH4710,RH4796,RH1847)
            deallocate(RH8498,QM8510)
         ENDIF
      ENDIF

      IF ( (IGET(318)>0).OR.(IGET(319)>0).OR.     &
           (IGET(320)>0))THEN
       allocate(RH4410(IM,jsta:jend),RH7294(IM,jsta:jend),   &
                RH4472(IM,jsta:jend),RH3310(IM,jsta:jend))
       CALL LFMFLD_GFS(RH4410,RH7294,RH4472,RH3310)
!     
!           SIGMA 0.44-1.00 MEAN RELATIVE HUMIIDITY.
            IF (IGET(318)>0) THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(RH4410(I,J) < spval) GRID1(I,J) = RH4410(I,J)*100.
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(318))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(318))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.72-0.94 MEAN RELATIVE HUMIIDITY.
            IF (IGET(319)>0) THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(RH7294(I,J) < spval) GRID1(I,J) = RH7294(I,J)*100.
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(319))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(319))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.44-0.72 MEAN RELATIVE HUMIIDITY.
            IF (IGET(320)>0) THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(RH4472(I,J) < spval) GRID1(I,J)=RH4472(I,J)*100.
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(320))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(320))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
            deallocate(RH4410,RH7294,RH4472,RH3310)
      END IF

! GFS computes sigma=0.9950 T, THETA, U, V from lowest two model level fields 
         IF ( (IGET(321)>0).OR.(IGET(322)>0).OR.     &
              (IGET(323)>0).OR.(IGET(324)>0).OR.     &
              (IGET(325)>0).OR.(IGET(326)>0)) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
	     DO I=1,IM
               EGRID2(I,J) = 0.995*PINT(I,J,LM+1)
               EGRID1(I,J) = LOG(PMID(I,J,LM)/EGRID2(I,J))   &
                           / LOG(PMID(I,J,LM)/PMID(I,J,LM-1))

        IF (MODELNAME == 'GFS' .OR. MODELNAME == 'FV3R' &
           .OR. MODELNAME == 'FV3R') THEN
               EGRID1(I,J) = LOG(PMID(I,J,LM)/EGRID2(I,J))   &
                           / max(1.e-6,LOG(PMID(I,J,LM)/PMID(I,J,LM-1)))
               EGRID1(I,J) =max(-10.0,min(EGRID1(I,J), 10.0))
                   IF ( ABS(PMID(I,J,LM)-PMID(I,J,LM-1)) < 0.5 ) THEN
                     EGRID1(I,J) = -1.
                   ENDIF
        ENDIF

	     END DO
	   END DO
! Temperature	   
	   IF (IGET(321)>0) THEN
             GRID1=spval
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                IF(T(I,J,LM)<spval.and.T(I,J,LM-1)<spval.and.EGRID1(I,J)<spval)&
                 GRID1(I,J) = T(I,J,LM)+(T(I,J,LM-1)-T(I,J,LM)) &
                            * EGRID1(I,J)
               ENDDO
             ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(321))
            fld_info(cfld)%lvl=LVLSXML(1,IGET(321))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
!            print *,'in miscln,iget(321,temp,sigmadata=',maxval(GRID1(1:im,jsta:jend)), &
!             minval(GRID1(1:im,jsta:jend)),'grib=',grib
            ENDIF
! Potential Temperature	    
            IF (IGET(322)>0) THEN
             GRID2=spval
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                IF(T(I,J,LM)<spval.and.T(I,J,LM-1)<spval.and.EGRID1(I,J)<spval)&
                 GRID2(I,J) = T(I,J,LM)+(T(I,J,LM-1)-T(I,J,LM))     &
                            * EGRID1(I,J)
               ENDDO
             ENDDO
             CALL CALPOT(EGRID2,GRID2(1,jsta),GRID1(1,jsta))
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(322))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(322))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF
! RH	    
            IF (IGET(323)>0) THEN
                 GRID1=spval
!$omp parallel do private(i,j,es1,qs1,rh1,es2,qs2,rh2)
             DO J=JSTA,JEND
               DO I=1,IM
                IF(PMID(I,J,LM)<spval.and.PMID(I,J,LM-1)<spval.and.&
                   Q(I,J,LM)<spval.and.Q(I,J,LM-1)<spval) THEN
                 ES1 = min(PMID(I,J,LM),FPVSNEW(T(I,J,LM)))
                 QS1 = CON_EPS*ES1/(PMID(I,J,LM)+CON_EPSM1*ES1)
                 RH1 = Q(I,J,LM)/QS1
                 ES2 = min(PMID(I,J,LM-1),FPVSNEW(T(I,J,LM-1)))
                 QS2 = CON_EPS*ES2/(PMID(I,J,LM-1)+CON_EPSM1*ES2)
                 RH2 = Q(I,J,LM-1)/QS2
                 GRID1(I,J) = (RH1+(RH2-RH1)*EGRID1(I,J))*100.
                 ENDIF
               ENDDO
             ENDDO
             CALL BOUND(GRID1,D00,H100)
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(323))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(323))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF    
! U	   
            IF (IGET(324)>0) THEN
             GRID1=spval
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                IF(UH(I,J,LM)<spval.and.UH(I,J,LM-1)<spval.and.EGRID1(I,J)<spval)&
                 GRID1(I,J) = UH(I,J,LM)+(UH(I,J,LM-1)-UH(I,J,LM))    &
                            * EGRID1(I,J)
               ENDDO
             ENDDO
            if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(324))
             fld_info(cfld)%lvl=LVLSXML(1,IGET(324))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
            ENDIF
! V	   
            IF (IGET(325)>0) THEN
             GRID1=spval
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                IF(VH(I,J,LM)<spval.and.VH(I,J,LM-1)<spval.and.EGRID1(I,J)<spval)&
                 GRID1(I,J) = VH(I,J,LM)+(VH(I,J,LM-1)-VH(I,J,LM))    &
                            * EGRID1(I,J)
               ENDDO
             ENDDO
            if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(325))
             fld_info(cfld)%lvl=LVLSXML(1,IGET(325))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
              jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
            endif
           ENDIF
! OMGA	   
           IF (IGET(326)>0) THEN
             GRID1=spval
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                IF(OMGA(I,J,LM)<spval.and.OMGA(I,J,LM-1)<spval.and.EGRID1(I,J)<spval)&
                GRID1(I,J) = OMGA(I,J,LM)+(OMGA(I,J,LM-1)-OMGA(I,J,LM))&
                            * EGRID1(I,J)
               ENDDO
             ENDDO
            if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(326))
             fld_info(cfld)%lvl=LVLSXML(1,IGET(326))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
              jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
            endif
           ENDIF
         END IF

!       MIXED LAYER CAPE AND CINS
!
         FIELD1=.FALSE.
         FIELD2=.FALSE.
!
         IF(IGET(032)>0)THEN
           IF(LVLS(3,IGET(032))>0)FIELD1=.TRUE.
         ENDIF
         IF(IGET(107)>0)THEN
           IF(LVLS(3,IGET(107))>0)FIELD2=.TRUE.
         ENDIF
!
         IF(IGET(582)>0)THEN
           FIELD1=.TRUE.
         ENDIF
         IF(IGET(583)>0)THEN
           FIELD2=.TRUE.
         ENDIF

         IF(FIELD1.OR.FIELD2)THEN
           ITYPE = 2
!
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               EGRID1(I,J) = -H99999
               EGRID2(I,J) = -H99999
               LB2(I,J)  = (LVLBND(I,J,1) + LVLBND(I,J,2) +           &
                            LVLBND(I,J,3))/3
               P1D(I,J)  = (PBND(I,J,1) + PBND(I,J,2) + PBND(I,J,3))/3
               T1D(I,J)  = (TBND(I,J,1) + TBND(I,J,2) + TBND(I,J,3))/3
               Q1D(I,J)  = (QBND(I,J,1) + QBND(I,J,2) + QBND(I,J,3))/3
             ENDDO
           ENDDO
!
           DPBND = 0.
           CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,LB2,EGRID1,           &
                        EGRID2,EGRID3,EGRID4,EGRID5)
 
           IF (IGET(582)>0) THEN
! dong add missing value for cape
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) THEN
                     GRID1(I,J) = EGRID1(I,J)
                     IF (SUBMODELNAME == 'RTMA' .OR. MODELNAME == 'FV3R') THEN
                      MLCAPE(I,J)=GRID1(I,J)
                     ENDIF
                   ENDIF
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(582))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(582))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                 jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
           ENDIF
           IF (IGET(583)>0) THEN
! dong add missing value for cape
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) GRID1(I,J) = - EGRID2(I,J)
                 ENDDO
               ENDDO
!
               CALL BOUND(GRID1,D00,H99999)
!
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) THEN
                      GRID1(I,J) = - GRID1(I,J)
                      IF (SUBMODELNAME == 'RTMA' .OR. MODELNAME == 'FV3R') THEN
                        MLCIN(I,J) = GRID1(I,J)
                      ENDIF
                   ENDIF
                 ENDDO
               ENDDO
!
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(583))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(583))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
           ENDIF
         ENDIF
     
!        MIXED LAYER LIFTING CONDENSATION PRESSURE AND HEIGHT.
!        EGRID1 IS LCL PRESSURE.  EGRID2 IS LCL HEIGHT.
!
         IF ( (IGET(109)>0).OR.(IGET(110)>0) ) THEN
            CALL CALLCL(P1D,T1D,Q1D,EGRID1,EGRID2)
            IF (IGET(109)>0) THEN
               GRID1=spval
	       DO J=JSTA,JEND
               DO I=1,IM
                 IF(T1D(I,J) < spval) GRID1(I,J)=EGRID2(I,J)
                 IF (SUBMODELNAME == 'RTMA' .OR. MODELNAME == 'FV3R') THEN
                     MLLCL(I,J) = GRID1(I,J)
                 ENDIF
               ENDDO
               ENDDO
!           
!               ID(1:25) = 0
!	       
!	       CALL GRIBIT(IGET(109),1,
!     X              GRID1,IM,JM)
            ENDIF
!	    
!            IF (IGET(110)>0) THEN
!	       DO J=JSTA,JEND
!               DO I=1,IM
!                 GRID1(I,J)=EGRID1(I,J)
!               ENDDO
!               ENDDO
!	       
!               ID(1:25) = 0
!	       
!	       CALL GRIBIT(IGET(110),1,
!     X              GRID1,IM,JM)
!            ENDIF
         ENDIF
!
!       MOST UNSTABLE CAPE-LOWEST 300 MB
!
         FIELD1=.FALSE.
         FIELD2=.FALSE.
!
         IF(IGET(032)>0)THEN
           IF(LVLS(4,IGET(032))>0)FIELD1=.TRUE.
         ENDIF
              
         IF(IGET(107)>0)THEN
           IF(LVLS(4,IGET(107))>0)FIELD2=.TRUE.
         ENDIF
!
         IF(IGET(584)>0)THEN
           FIELD1=.TRUE.
         ENDIF
         IF(IGET(585)>0)THEN
           FIELD2=.TRUE.
         ENDIF
!
         IF(FIELD1.OR.FIELD2)THEN
           ITYPE = 1
!
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               EGRID1(I,J) = -H99999
               EGRID2(I,J) = -H99999
             ENDDO
           ENDDO
              
           DPBND = 300.E2
           CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,LB2,EGRID1,     &
                        EGRID2,EGRID3,EGRID4,EGRID5)
!           IF (SUBMODELNAME == 'RTMA' .OR. MODELNAME == 'FV3R') THEN
!               MUMIXR(I,J) = Q1D(I,J) 
!           ENDIF
           IF (IGET(584)>0) THEN
! dong add missing value to cin
               GRID1 = spval
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) THEN
                      GRID1(I,J) = EGRID1(I,J)
                      IF (SUBMODELNAME == 'RTMA' .OR. MODELNAME == 'FV3R') THEN
                        MUCAPE(I,J) = GRID1(I,J)
                      ENDIF
                   ENDIF
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
!               IF (SUBMODELNAME == 'RTMA' .OR. MODELNAME == 'FV3R') THEN
!                    CALL BOUND(MUCAPE,D00,H99999)
!               ENDIF
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(584))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(584))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif

           ENDIF
                
           IF (IGET(585)>0) THEN
! dong add missing value to cin
               GRID1 = spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) GRID1(I,J) = - EGRID2(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) THEN 
                   GRID1(I,J) = - GRID1(I,J)
                       IF (SUBMODELNAME == 'RTMA' .OR. &
                           MODELNAME == 'FV3R') THEN
                              MUCAPE(I,J) = GRID1(I,J)
                              MUQ1D(I,J) = Q1D(I,J)
                       ENDIF
                   ENDIF
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(585))
                 fld_info(cfld)%lvl=LVLSXML(1,IGET(585))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif

            ENDIF
              
!    EQUILLIBRIUM HEIGHT
           IF (IGET(443)>0) THEN
               GRID1 = spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) GRID1(I,J) = EGRID4(I,J)
                 ENDDO
               ENDDO
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(443))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(443))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
           ENDIF
!Equilibrium Temperature
            IF (IGET(982)>0) THEN
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) =  TEQL(I,J)
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(982))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(982))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF


!      PRESSURE OF LEVEL FROM WHICH 300 MB MOST UNSTABLE CAPE
!      PARCEL WAS LIFTED (eq. PRESSURE OF LEVEL OF HIGHEST THETA-E)
           IF (IGET(246)>0) THEN
              GRID1 = spval
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) GRID1(I,J) = EGRID3(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
!               print *,'in miscln,PLPL=',maxval(grid1(1:im,jsta:jend)),  &
!                 minval(grid1(1:im,jsta:jend))
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(246))
                 fld_info(cfld)%lvl=LVLSXML(1,IGET(246))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
           ENDIF   ! 246

!    GENERAL THUNDER PARAMETER  ??? 458 ???
        IF (IGET(444)>0) THEN
             GRID1 = spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                  IF(CPRATE(I,J) < spval) THEN
                   IF (CPRATE(I,J) > PTHRESH) THEN
                    GRID1(I,J) = EGRID5(I,J)
                   ELSE
                    GRID1(I,J) = 0
                   ENDIF
                  ENDIF
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(444))
                 fld_info(cfld)%lvl=LVLSXML(1,IGET(444))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
        ENDIF
      ENDIF

!
! --- Effective (inflow) Layer (EL)
!
      IF ( EL_SCHEME > 0 ) THEN
        ALLOCATE(EL_BASE(IM,JSTA_2L:JEND_2U))
        ALLOCATE(EL_TOPS(IM,JSTA_2L:JEND_2U))
        ALLOCATE(FOUND_BASE(IM,JSTA_2L:JEND_2U))
        ALLOCATE(FOUND_TOPS(IM,JSTA_2L:JEND_2U))
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            EL_BASE(I,J) = LM
            EL_TOPS(I,J) = LM
            FOUND_BASE(I,J) = .FALSE.
            FOUND_TOPS(I,J) = .FALSE.
          ENDDO
        ENDDO
      END IF
!
      IF ( EL_SCHEME == 1 ) THEN
        ITYPE = 2
        DPBND = 0.

        DO L = LM, 1, -1

!         SET AIR PARCELS FOR LEVEL L
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              EGRID1(I,J) = -H99999
              EGRID2(I,J) = -H99999
              IDUMMY(I,J) = 0
              P1D(I,J)    = PMID(I,J,L)
              T1D(I,J)    = T(I,J,L)
              Q1D(I,J)    = Q(I,J,L)
            ENDDO
          ENDDO

!---      CALCULATE CAPE/CIN FOR ALL AIR PARCELS on LEVEL L
          IF (debugprint) WRITE(1000+ME,'(1x,A,I2.2,2x,A,I3)')          &
           'NEW EL_SCHEME:', EL_SCHEME,                                 &
           '   CALCULATING CAPE/CINS ON LEVEL:',L
          CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,IDUMMY,EGRID1,           &
                          EGRID2,EGRID3,EGRID4,EGRID5)

!---      CHECK CAPE/CIN OF EACH AIR PARCELS WITH EL CRITERIA
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              IF ( .NOT. FOUND_BASE(I,J) ) THEN
                IF ( EGRID1(I,J) >= 100. .AND. EGRID2(I,J) >= -250. ) THEN
                  EL_BASE(I,J) = L
                  FOUND_BASE(I,J) = .TRUE.
                ELSE
                  EL_BASE(I,J) = LM
                  FOUND_BASE(I,J) = .FALSE.
                END IF 
              ELSE
                IF ( .NOT. FOUND_TOPS(I,J) ) THEN
                  IF ( EGRID1(I,J) < 100. .OR. EGRID2(I,J) < -250. ) THEN
                    EL_TOPS(I,J) = L + 1
                    FOUND_TOPS(I,J) = .TRUE.
                  ELSE
                    EL_TOPS(I,J) = LM
                    FOUND_TOPS(I,J) = .FALSE.
                  END IF 
                END IF 
              END IF
            ENDDO
          ENDDO

        END DO ! L
!
      ELSE IF ( EL_SCHEME == 2 ) THEN

!---    SEARCH FOR EL ALONG EACH PROFILE/COLUMN    
        ALLOCATE(L_THETAE_MAX(IM,JSTA:JEND))
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            EGRID1(I,J) = -H99999
            EGRID2(I,J) = -H99999
            L_THETAE_MAX(I,J) = -9999
          ENDDO
        ENDDO

!------ SEARCH FOR PARCEL WITH MAX THETA-E OF EACH PROFILE/COLUMN
        DO L = LM,1,-1
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              EGRID1(I,J) = -H99999
              P1D(I,J)    = PMID(I,J,L)
              T1D(I,J)    = T(I,J,L)
              Q1D(I,J)    = Q(I,J,L)
            ENDDO
          ENDDO
!          CALL CALTHTE(PMID(1,jsta,L),T(1,jsta,L),Q(1,jsta,L),EGRID1)
          CALL CALTHTE(P1D,T1D,Q1D,EGRID1)
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              IF (EGRID1(I,J) > EGRID2(I,J)) THEN
                 EGRID2(I,J) = EGRID1(I,J)
                 L_THETAE_MAX(I,J) = L
              END IF
            END DO
          END DO

        END DO  ! L

!---    SET UP THE PARCELS WIH MAX THETA-E
        ALLOCATE(CAPE9(IM,JSTA:JEND))
        ALLOCATE(CINS9(IM,JSTA:JEND))
!$omp parallel do private(i,j,LLL)
        DO J=JSTA,JEND
          DO I=1,IM
            LLL         = L_THETAE_MAX(I,J)
            P1D(I,J)    = PMID(I,J,LLL)
            T1D(I,J)    = T(I,J,LLL)
            Q1D(I,J)    = Q(I,J,LLL)
          ENDDO
        ENDDO
!---    COMPUTE CAPE/CIN OF PARCEL WITH MAX THETA-E
        ITYPE = 2
        DPBND = 0.
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            EGRID1(I,J) = -H99999
            EGRID2(I,J) = -H99999
            CAPE9(I,J) = -H99999
            CINS9(I,J) = -H99999
            IDUMMY(I,J) = 0
          ENDDO
        ENDDO
        CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,IDUMMY,EGRID1,   &
                     EGRID2,EGRID3,EGRID4,EGRID5) 
!---    SANITY CHECK IF CAPE/CIN OF MAX THETA-E PARCEL OF THIS COLUMN
!$omp parallel do private(i,j,LLL)
        DO J=JSTA,JEND
          DO I=1,IM
            CAPE9(I,J) = EGRID1(I,J)
            CINS9(I,J) = EGRID2(I,j)
          ENDDO
        ENDDO
        
        ALLOCATE(TPAR_B(LM), TPAR_T(LM))
        ALLOCATE(TPAR_TMP(LM))
        ALLOCATE(P_AMB(LM), T_AMB(LM), Q_AMB(LM), ZINT_AMB(LM))
        IF (debugprint) THEN
          ALLOCATE(TPAR_BASE(IM,JSTA:JEND,LM))
          ALLOCATE(TPAR_TOPS(IM,JSTA:JEND,LM))
!$omp  parallel do private(i,j,l)
          DO L=1,LM
            DO J=JSTA,JEND
              DO I=1,IM
                TPAR_BASE(I,J,L) = DM9999
                TPAR_TOPS(I,J,L) = DM9999
              ENDDO
            ENDDO
          ENDDO
        END IF

!$OMP  PARALLEL DO PRIVATE(i,j,l,P_AMB,T_AMB,Q_AMB,ZINT_AMB,LMASK,PSFC, &
!$OMP&                     TPAR_B,TPAR_T,LPAR0,CAPE_PAR,CINS_PAR,       &
!$OMP&                     PARCEL0,LLCL_PAR,LEQL_PAR,TPAR_TMP)
        DO J=JSTA,JEND
          DO I=1,IM

!-----------GET THE AMBIENT PROFILE
            DO L=1,LM
              P_AMB(L) = PMID(I,J,L)
              T_AMB(L) = T(I,J,L)
              Q_AMB(L) = Q(I,J,L)
              ZINT_AMB(L) = ZINT(I,J,L)
            END DO
            LMASK = NINT(LMH(I,J))
            PSFC  = PMID(I,J,NINT(LMH(I,J)))
            EL_BASE(I,J)=LM
            EL_TOPS(I,J)=LM
            FOUND_BASE(I,J) = .FALSE.
            FOUND_TOPS(I,J) = .FALSE.
            TPAR_B(1:LM) = DM9999
            TPAR_T(1:LM) = DM9999
            TPAR_TMP(1:LM) = DM9999

!----------      SANITY  CHECK                                                   
!---------- USING CAPE/CINS OF AIR PARCEL WITH MAX THETA-E FOR SANITY CHECK FIRST
!---------- CAPE4/CINS4 ARE CAPE/CINS CALCULATED FROM CODE ABOVE FOR AIR PARCEL
!---------- WITH MAX THETA-E.
            IF (CAPE9(I,J) >= 100. .OR. CINS9(I,J) >= -250.) THEN  !EL FOR THIS COLUMN
!             EL SHOULD EXIST ALONG THIS COLUMN/PROFILE.
!-------------NEED TO SEARCH FROM BOTTOM TO TOP, COMPUTE CAPE/CINS OF AIR PARCEL
!             INITIATING AT EACH LEVEL UNTIL THE CRITERIA MEETS.

              VLOOP: DO L=LM,1,-1
!---------------CALCULATE CAPE/CIN OF AIR PARCEL INITIALIZED AT LEVEL L
                IF (.NOT. FOUND_BASE(I,J)) THEN  !SEARCH FOR BASE FIRST
                  LPAR0 = FLOAT(L)
                  CAPE_PAR = D00
                  CINS_PAR = D00
                  PARCEL0(1) = LPAR0
                  PARCEL0(2:4) = D00
                  LLCL_PAR = 1
                  LEQL_PAR = LM
                  TPAR_B(1:LM) = DM9999
                  TPAR_T(1:LM) = DM9999

                  WRITE(1000+me,*)"CALCAPE2_1D_base: parcel@I J L:",I,J,L

                  CALL CALCAPE1D(P_AMB,T_AMB,Q_AMB,ZINT_AMB,LPAR0,              &
                                 PSFC,LMASK,CAPE_PAR,CINS_PAR,TPAR_B,          &
                                 PARCEL0,LLCL_PAR,LEQL_PAR)
                  IF (CAPE_PAR >= 100. .AND. CINS_PAR >= -250.) THEN  ! BASE OF EL
                    EL_BASE(I,J) = L
                    FOUND_BASE(I,J) = .TRUE.
                  ELSE
                    TPAR_B(1:LM) = DM9999     ! RESET TPAR IF NOT FOUND BASE
                  ENDIF
                ELSE                           !SEARCH FOR TOP AFTER BASE IS FOUND
                  IF (.NOT. FOUND_TOPS(I,J)) THEN  !SEARCH FOR BASE FIRST
                    LPAR0 = FLOAT(L)
                    CAPE_PAR = D00
                    CINS_PAR = D00
                    TPAR_T(1:LM) = DM9999
                    PARCEL0(1) = LPAR0
                    PARCEL0(2:4) = D00
                    LLCL_PAR = 1
                    LEQL_PAR = LM

                    WRITE(1000+me,*)"CALCAPE2_1D_top: parcel@I J L:",I,J,L

                    CALL CALCAPE1D(P_AMB,T_AMB,Q_AMB,ZINT_AMB,LPAR0,            &
                                   PSFC,LMASK,CAPE_PAR,CINS_PAR,TPAR_T,        &
                                   PARCEL0,LLCL_PAR,LEQL_PAR)
                    IF (CAPE_PAR < 100. .OR. CINS_PAR < -250.) THEN  ! TOP OF EL
                      EL_TOPS(I,J) = L+1
                      FOUND_TOPS(I,J) = .TRUE.
                      IF (EL_TOPS(I,J) == EL_BASE(I,J)) THEN
                        TPAR_T(1:LM) = TPAR_B(1:LM)
                      ELSE IF (EL_TOPS(I,J) <  EL_BASE(I,J) ) THEN
                        TPAR_T(1:LM) = TPAR_TMP(1:LM)
                      ELSE
                        WRITE(0,'(1x,A,A)') "TOP OF EFFECTIVE LAYER IS",        &
                             " LOWER THAN BASE. WRONG! ABORT ..."
                        STOP 9
                      END IF
                      EXIT VLOOP
                    ELSE
                      TPAR_TMP(1:LM) = TPAR_T(1:LM)
                      TPAR_T(1:LM) = DM9999   ! RESET TPAR IF NOT FOUND TOP
                    ENDIF
                  ENDIF    ! FOUND_TOPS

                ENDIF    ! FOUND_BASE OR NOT

              ENDDO VLOOP

              IF (debugprint) THEN
                DO L=1,LM
                  TPAR_BASE(I,J,L) = TPAR_B(L)
                  TPAR_TOPS(I,J,L) = TPAR_T(L)
                END DO
              END IF

            ENDIF        ! IF PASSING SANITY CHECK

            IF ( FOUND_BASE(I,J) /= FOUND_TOPS(I,J) ) THEN
              WRITE(0,'(1x,A,A,A,I6,1x,I6)') "BASE & TOP OF ",                 &
                    " EFFECTIVE LAYER ARE NOT FOUND TOGETHER. WRONG! ",        &
                    " ABORT! ABORT! ...  AT GRID POINT: ", I, J
              STOP 10
            END IF

          ENDDO          ! J
        ENDDO            ! I

        IF(ALLOCATED(L_THETAE_MAX)) DEALLOCATE(L_THETAE_MAX)
        IF(ALLOCATED(CAPE9)) DEALLOCATE(CAPE9)
        IF(ALLOCATED(CINS9)) DEALLOCATE(CINS9)
        IF(ALLOCATED(TPAR_B)) DEALLOCATE(TPAR_B)
        IF(ALLOCATED(TPAR_T)) DEALLOCATE(TPAR_T)
        IF(ALLOCATED(TPAR_TMP)) DEALLOCATE(TPAR_TMP)
        IF(ALLOCATED(P_AMB)) DEALLOCATE(P_AMB)
        IF(ALLOCATED(T_AMB)) DEALLOCATE(T_AMB)
        IF(ALLOCATED(Q_AMB)) DEALLOCATE(Q_AMB)
        IF(ALLOCATED(ZINT_AMB)) DEALLOCATE(ZINT_AMB)

      END IF             ! EL_SCHEME :1 OR 2

      IF (ALLOCATED(FOUND_BASE))    DEALLOCATE(FOUND_BASE)
      IF (ALLOCATED(FOUND_TOPS))    DEALLOCATE(FOUND_TOPS)

      IF (debugprint .AND. EL_SCHEME > 0) THEN
        WRITE(IM_CH,'(I5.5)') IM
        WRITE(JSTA_CH,'(I5.5)') JSTA
        WRITE(JEND_CH,'(I5.5)') JEND
        EFFL_FNAME="EFFL_NEW_"//IM_CH//"_"//JSTA_CH//"_"//JEND_CH               &
                   //".dat"
        EFFL_FNAME2="EFFL_NEW_LVLS_"//IM_CH//"_"//JSTA_CH//"_"//JEND_CH              &
                   //".dat"
        IUNIT=10000+JSTA
        IUNIT2=20000+JSTA
        IREC=0
        IREC2=0
        OPEN(IUNIT,FILE=TRIM(ADJUSTL(EFFL_FNAME)),FORM='FORMATTED')
        IF (EL_SCHEME == 2) THEN
          OPEN(IUNIT2,FILE=TRIM(ADJUSTL(EFFL_FNAME2)),FORM='FORMATTED')
        END IF
!        OPEN(IUNIT,FILE=TRIM(ADJUSTL(EFFL_FNAME)),FORM='UNFORMATTED',          &
!             ACCESS='DIRECT',RECL=4*6)
        DO J=JSTA,JEND
          DO I=1,IM
            IREC = IREC + 1
            IREC2 = IREC2 + 1
!            WRITE(IUNIT,'(1x,I6,2x,I6,2x,I6,2x,I6)')I,J,EL_BASE(I,J),EL_TOPS(I,J)
            WRITE(IUNIT,'(1x,I6,2x,I6,2(2x,I6,2x,F12.3))') I, J,                &
                  EL_BASE(I,J),PMID(I,J,EL_BASE(I,J)),                          &
                  EL_TOPS(I,J),PMID(I,J,EL_TOPS(I,J))
!            WRITE(IUNIT,REC=IREC) I, J,                &
!                  EL_BASE(I,J),PMID(I,J,EL_BASE(I,J)),                          &
!                  EL_TOPS(I,J),PMID(I,J,EL_TOPS(I,J))
            IF (EL_SCHEME == 2) THEN
              WRITE(IUNIT2,'(1x,I6,2x,I6,2(2x,I6,2x,F12.3))') I, J,               &
                    EL_BASE(I,J),PMID(I,J,EL_BASE(I,J)),                          &
                    EL_TOPS(I,J),PMID(I,J,EL_TOPS(I,J))
              DO L=LM,1,-1
                IREC2=IREC2+1
                WRITE(IUNIT2,'(1x,I4,5(2x,F12.3))')                              &
                  LM+1-L, PMID(I,J,L), ZINT(I,J,L), T(I,J,L),                   &
                  TPAR_BASE(I,J,L), TPAR_TOPS(I,J,L)
!                WRITE(IUNIT2,REC=IREC2)                              &
!                  LM+1-L, PMID(I,J,L), ZINT(I,J,L), T(I,J,L),                   &
!                  TPAR_BASE(I,J,L), TPAR_TOPS(I,J,L)
              END DO
            END IF
          END DO
        ENDDO
        CLOSE(IUNIT)

        IF (EL_SCHEME == 2) THEN
          CLOSE(IUNIT2)
        END IF

      ENDIF

      IF(ALLOCATED(TPAR_BASE)) DEALLOCATE(TPAR_BASE)
      IF(ALLOCATED(TPAR_TOPS)) DEALLOCATE(TPAR_TOPS)
!
!       EXPAND HRRR CAPE/CIN RELATED VARIABLES
!       
!    CAPE AND CINS 0-3KM, FOLLOW ML PROCEDURE WITH HEIGHT 0-3KM
!
         FIELD1=.FALSE.
         FIELD2=.FALSE.
!
         IF(IGET(032)>0)THEN
           IF(LVLS(3,IGET(032))>0)FIELD1=.TRUE.
         ENDIF
         IF(IGET(107)>0)THEN
           IF(LVLS(3,IGET(107))>0)FIELD2=.TRUE.
         ENDIF
!
         IF(IGET(950)>0)THEN
           FIELD1=.TRUE.
         ENDIF
         IF(IGET(951)>0)THEN
           FIELD2=.TRUE.
         ENDIF
!
!         IF(FIELD1)ITYPE=2
!         IF(FIELD2)ITYPE=2

         IF(FIELD1.OR.FIELD2)THEN
           ITYPE = 2

!
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               EGRID1(I,J) = -H99999
               EGRID2(I,J) = -H99999
               EGRID3(I,J) = -H99999
               EGRID4(I,J) = -H99999
               EGRID5(I,J) = -H99999
               EGRID6(I,J) = -H99999
               EGRID7(I,J) = -H99999
               EGRID8(I,J) = -H99999
!          ENDDO
!          ENDDO
!          DO J=JSTA,JEND
!          DO I=1,IM
               LB2(I,J)  = (LVLBND(I,J,1) + LVLBND(I,J,2) +           &
                            LVLBND(I,J,3))/3
               P1D(I,J)  = (PBND(I,J,1) + PBND(I,J,2) + PBND(I,J,3))/3
               T1D(I,J)  = (TBND(I,J,1) + TBND(I,J,2) + TBND(I,J,3))/3
               Q1D(I,J)  = (QBND(I,J,1) + QBND(I,J,2) + QBND(I,J,3))/3
             ENDDO
           ENDDO

           DPBND = 0.
           CALL CALCAPE2(ITYPE,DPBND,P1D,T1D,Q1D,LB2,            &
                         EGRID1,EGRID2,EGRID3,EGRID4,EGRID5,     &
                         EGRID6,EGRID7,EGRID8)

!                        CAPE1, CINS2, LFC3,  ESRHL4,ESRHH5,
!                        DCAPE6,DGLD7, ESP8)
!
           IF (IGET(950)>0) THEN
! dong add missing value for cape
              GRID1=spval
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                DO I=1,IM
                  IF(T1D(I,J) < spval) GRID1(I,J) = EGRID1(I,J)
                ENDDO
              ENDDO
             CALL BOUND(GRID1,D00,H99999)
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(950))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(950))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
           ENDIF   !950
!
           IF (IGET(951)>0) THEN
! dong add missing value for cape
              GRID1=spval
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(T1D(I,J) < spval) GRID1(I,J) = - EGRID2(I,J)
               ENDDO
             ENDDO
!
             CALL BOUND(GRID1,D00,H99999)
!
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(T1D(I,J) < spval) GRID1(I,J) = - GRID1(I,J)
               ENDDO
             ENDDO
!
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(951))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(951))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF   !951

!    LFC HEIGHT

            IF (IGET(952)>0) THEN
             GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) GRID1(I,J) = EGRID3(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(952))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(952))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
            ENDIF   !952


!    EFFECTIVE STORM RELATIVE HELICITY AND STORM MOTION.

         allocate(UST(IM,jsta_2l:jend_2u),VST(IM,jsta_2l:jend_2u),     &
                  HELI(IM,jsta_2l:jend_2u,2))
         allocate(LLOW(IM,jsta_2l:jend_2u),LUPP(IM,jsta_2l:jend_2u),   &
                  CANGLE(IM,jsta_2l:jend_2u))
!         allocate(LLOW1(IM,jsta_2l:jend_2u),LUPP1(IM,jsta_2l:jend_2u))

       iget1 = IGET(953)
       iget2 = -1
       iget3 = -1
       if (iget1 > 0) then
         iget2 = LVLS(1,iget1)
         iget3 = LVLS(2,iget1)
       endif
      if(me==0) write(0,*) '953 ',iget1,iget2,iget3
       IF (iget1 > 0 .OR. IGET(162) > 0 .OR. IGET(953) > 0) THEN
         DEPTH(1) = 3000.0
         DEPTH(2) = 1000.0
!$omp parallel do private(i,j)
         DO J=JSTA,JEND
           DO I=1,IM
             LLOW(I,J) = INT(EGRID4(I,J))
             LUPP(I,J) = INT(EGRID5(I,J))
           ENDDO
         ENDDO
!---     OUTPUT EL BASE & TOP COMPUTED BY OLD SCHEME
         IF (debugprint) THEN
           WRITE(IM_CH,'(I5.5)') IM
           WRITE(JSTA_CH,'(I5.5)') JSTA
           WRITE(JEND_CH,'(I5.5)') JEND
           EFFL_FNAME="EFFL_OLD_"//IM_CH//"_"//JSTA_CH//"_"//JEND_CH               &
                      //".dat"
           IUNIT=10000+JSTA
           IREC=0
           OPEN(IUNIT,FILE=TRIM(ADJUSTL(EFFL_FNAME)),FORM='FORMATTED')
           DO J=JSTA,JEND
             DO I=1,IM
               IREC = IREC + 1
!               WRITE(IUNIT,'(1x,I6,2x,I6,2x,I6,2x,I6)')I,J,LLOW(I,J),LUPP(I,J)
               WRITE(IUNIT,'(1x,I6,2x,I6,2(2x,I6,2x,F12.3))') I, J,                &
                    LLOW(I,J),PMID(I,J,LLOW(I,J)),                                 &
                    LUPP(I,J),PMID(I,J,LUPP(I,J))
             END DO
           ENDDO
           CLOSE(IUNIT)
         ENDIF

!---  IF USSING EL BASE & TOP COMPUTED BY NEW SCHEME FOR THE RELATED VARIABLES
         IF ( EL_SCHEME > 0 ) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
                LLOW(I,J) = EL_BASE(I,J)
                LUPP(I,J) = EL_TOPS(I,J)
             ENDDO
           ENDDO
         END IF

!         CALL CALHEL(DEPTH,UST,VST,HELI,USHR1,VSHR1,USHR6,VSHR6)
         CALL CALHEL2(LLOW,LUPP,DEPTH,UST,VST,HELI,CANGLE)
!                                                  
         IF (iget2 > 0) then
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = HELI(I,J,1)
             !  GRID1(I,J) = HELI(I,J,2)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(iget1)
             fld_info(cfld)%lvl=LVLSXML(1,iget1)
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF

       ENDIF   !953


        IF (SUBMODELNAME == 'RTMA' .OR. MODELNAME=='FV3R') THEN  !Start RTMA block

!EL field allocation

         allocate(ESHR(IM,jsta_2l:jend_2u),UVECT(IM,jsta_2l:jend_2u),&
                  VVECT(IM,jsta_2l:jend_2u),HTSFC(IM,jsta_2l:jend_2u))
         allocate(EFFUST(IM,jsta_2l:jend_2u),EFFVST(IM,jsta_2l:jend_2u),&
                  ESRH(IM,jsta_2l:jend_2u))

!       
           DO J=JSTA,JEND
             DO I=1,IM
               MAXTHE(I,J)=-H99999
               THE(I,J)=-H99999
               MAXTHEPOS(I,J)=0
               MUQ1D(I,J) = 0.
             ENDDO
           ENDDO
           DO L=LM,1,-1

             DO J=JSTA,JEND
               DO I=1,IM
                 EGRID1(I,J) = -H99999
                 P1D(I,J)=PMID(I,J,L)
                 T1D(I,J)=T(I,J,L)
                 Q1D(I,J)=Q(I,J,L)
               ENDDO
             ENDDO
             CALL CALTHTE(P1D,T1D,Q1D,EGRID1)
             DO J=JSTA,JEND
               DO I=1,IM
                 THE(I,J)=EGRID1(I,J)
                 IF(THE(I,J)>=MAXTHE(I,J))THEN
                    MAXTHE(I,J)=THE(I,J)
                    MAXTHEPOS(I,J)=L
                    MUQ1D(I,J) = Q(I,J,L)  ! save the Q of air parcel with max theta-e (MU Parcel)
                 ENDIF 
               ENDDO
             ENDDO

           ENDDO

!
!get surface height
        IF(gridtype == 'E')THEN
        JVN =  1
        JVS = -1
        do J=JSTA,JEND
          IVE(J) = MOD(J,2)
          IVW(J) = IVE(J)-1
        enddo
        ISTART = 2
        ISTOP  = IM-1
        JSTART = JSTA_M
        JSTOP  = JEND_M
      ELSE IF(gridtype == 'B')THEN
        JVN = 1
        JVS = 0
        do J=JSTA,JEND
          IVE(J)=1
          IVW(J)=0
        enddo
        ISTART = 2
        ISTOP  = IM-1
        JSTART = JSTA_M
        JSTOP  = JEND_M
      ELSE
        JVN = 0
        JVS = 0
        do J=JSTA,JEND
          IVE(J) = 0
          IVW(J) = 0
        enddo
        ISTART = 1
        ISTOP  = IM
        JSTART = JSTA
        JSTOP  = JEND
      END IF

      IF(gridtype /= 'A') CALL EXCH(FIS(1:IM,JSTA:JEND))
        DO J=JSTART,JSTOP
          DO I=ISTART,ISTOP
            IE = I+IVE(J)
            IW = I+IVW(J)
            JN = J+JVN
            JS = J+JVS
            IF (gridtype=='B')THEN
            HTSFC(I,J)=(0.25/g)*(FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(IE,JN))
            ELSE
            HTSFC(I,J)=(0.25/g)*(FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(I,JS))
            ENDIF
          ENDDO
        ENDDO

!Height of effbot
            IF (IGET(979)>0) THEN
             GRID1=spval
             DO J=JSTA,JEND
               DO I=1,IM
                IF(ZINT(I,J,LLOW(I,J))<spval.and.HTSFC(I,J)<spval)&
                 GRID1(I,J) = ZINT(I,J,LLOW(I,J)) - HTSFC(I,J)
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(979))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(979))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF
!Height of effbot
            IF (IGET(980)>0) THEN
             GRID1=spval
             DO J=JSTA,JEND
               DO I=1,IM
                IF(ZINT(I,J,LUPP(I,J))<spval.and.HTSFC(I,J)<spval)&
                 GRID1(I,J) = ZINT(I,J,LUPP(I,J)) - HTSFC(I,J)
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(980))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(980))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF

!U inflow based to 50% EL shear vector

            IF (IGET(983)>0) THEN
             GRID1=spval
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(LLOW(I,J)<spval.and.LUPP(I,J)<spval.and.&
                 UH(I,J,MIDCAL)<spval.and.UH(I,J,LLOW(I,J))<spval)THEN
                       MIDCAL=INT(LLOW(I,J)+D50*(LUPP(I,J)-LLOW(I,J)))       
                                                            !mid-layer 
                                                            !vertical
                                                            !index
                       UVECT(I,J)=UH(I,J,MIDCAL)-UH(I,J,LLOW(I,J))
                       GRID1(I,J)=UVECT(I,J)
                  ENDIF
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(983))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(983))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF

!V inflow based to 50% EL shear vector
            IF (IGET(984)>0) THEN
             GRID1=spval
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(LLOW(I,J)<spval.and.LUPP(I,J)<spval.and.&
                 VH(I,J,MIDCAL)<spval.and.VH(I,J,LLOW(I,J))<spval)THEN
                       MIDCAL=INT(LLOW(I,J)+D50*(IEQL(I,J)-LLOW(I,J)))
                                                            !mid-layer 
                                                            !vertical
                                                            !index
                       VVECT(I,J)=VH(I,J,MIDCAL)-VH(I,J,LLOW(I,J))
                       GRID1(I,J)=VVECT(I,J)
                 ENDIF
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(984))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(984))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF

!Inflow based (ESFC) to (50%) EL shear magnitude
            IF (IGET(985)>0) THEN
             GRID1=spval
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(UVECT(I,J)<spval.and.VVECT(I,J)<spval) THEN
                       ESHR(I,J)=SQRT((UVECT(I,J)**2)+(VVECT(I,J))**2)
                                                               !effshear
                                                               !calc
                       GRID1(I,J)=ESHR(I,J) !Effective
                 ENDIF
             ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(985))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(985))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF

! Effective Helicity

       CALL CALHEL3(LLOW,LUPP,EFFUST,EFFVST,ESRH)

!U Bunkers Effective right motion
!

            IF (IGET(986)>0) THEN
             GRID1=spval
             DO J=JSTA,JEND
               DO I=1,IM
                  IF(LLOW(I,J)<spval.and.LUPP(I,J)<spval)&
                       GRID1(I,J)=EFFUST(I,J)
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(986))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(986))
! $omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF

!V Bunkers Effective right motion
            IF (IGET(987)>0) THEN
             GRID1=spval
             DO J=JSTA,JEND
               DO I=1,IM
                  IF(LLOW(I,J)<spval.and.LUPP(I,J)<spval)&
                       GRID1(I,J)=EFFVST(I,J)
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(987))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(987))
! $omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF

!Effective layer helicity
            IF (IGET(988)>0) THEN
             GRID1=spval
             DO J=JSTA,JEND
               DO I=1,IM
                  IF(LLOW(I,J)<spval.and.LUPP(I,J)<spval)&
                       GRID1(I,J)=ESRH(I,J)
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(988))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(988))
! $omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF

!Effective Layer Tornado Parameter
            IF (IGET(989)>0) THEN
            DO J=JSTA,JEND
               DO I=1,IM
                IF (MLLCL(I,J)>D2000) THEN
                        MLLCLtmp=D00
                ELSEIF (MLLCL(I,J)<D1000) THEN
                        MLLCLtmp=1.0
                ELSE
                        MLLCLtmp=((D2000-MLLCL(I,J))/D1000)
                ENDIF
                IF (ESHR(I,J)<12.5) THEN
                        ESHRtmp=D00
                ELSEIF (ESHR(I,J)>30.0) THEN
                        ESHRtmp=1.5
                ELSE
                        ESHRtmp=(ESHR(I,J)/20.)
                ENDIF
                IF (MLCIN(I,J)>-50.) THEN
                        MLCINtmp=1.0
                ELSEIF (MLCIN(I,J)<-200.) THEN
                        MLCINtmp=D00
                ELSE
                        MLCINtmp=(200.+MLCIN(I,J))/150.
                ENDIF
                STP=(MLCAPE(I,J)/D1500)*MLLCLtmp*(ESRH(I,J)/150.)*&
                        ESHRtmp*MLCINtmp
                GRID1(I,J) = spval
                IF(LLOW(I,J)<spval.and.LUPP(I,J)<spval) THEN
                IF (STP>0) THEN
                   GRID1(I,J)=STP
                ELSE
                   GRID1(I,J)=D00
                ENDIF
                ENDIF
               ENDDO
            ENDDO
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(989))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(989))
! $omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
           ENDIF

!Fixed Layer Tornado Parameter
            IF (IGET(990)>0) THEN
            DO J=JSTA,JEND
             DO I=1,IM
                 LLMH = NINT(LMH(I,J))
                 P1D(I,J) = PMID(I,J,LLMH)
                 T1D(I,J) = T(I,J,LLMH)
                 Q1D(I,J) = Q(I,J,LLMH)
             ENDDO
            ENDDO
           CALL CALLCL(P1D,T1D,Q1D,EGRID1,EGRID2)
            DO J=JSTA,JEND
             DO I=1,IM
                SLCL(I,J)=EGRID2(I,J)
             ENDDO
            ENDDO
            ITYPE  = 1
            DPBND  = 10.E2
            dummy  = 0.
            idummy = 0
            CALL CALCAPE(ITYPE,DPBND,dummy,dummy,dummy,&
                         idummy,EGRID1,EGRID2,&
                         EGRID3,dummy,dummy)

            DO J=JSTA,JEND
               DO I=1,IM
                IF (SLCL(I,J)>D2000) THEN
                        SLCLtmp=D00
                ELSEIF (SLCL(I,J)<=D1000) THEN
                        SLCLtmp=1.0
                ELSE
                        SLCLtmp=((D2000-SLCL(I,J))/D1000)
                ENDIF
                IF (FSHR(I,J)<12.5) THEN
                        FSHRtmp=D00
                ELSEIF (FSHR(I,J)>30.0) THEN
                        FSHRtmp=1.5
                ELSE
                        FSHRtmp=(FSHR(I,J)/20.)
                ENDIF
                IF (EGRID2(I,J)>-50.) THEN
                        SCINtmp=1.0
                ELSEIF (EGRID2(I,J)<-200.) THEN
                        SCINtmp=D00
                ELSE
                        SCINtmp=((200.+EGRID2(I,J)/150.))
                ENDIF
                STP=(EGRID1(I,J)/D1500)*SLCLtmp*(HELI(I,J,2)/150.)*&
                        FSHRtmp*SCINtmp
                GRID1(I,J) = spval
                IF(T1D(I,J) < spval) THEN
                IF (STP>0) THEN
                   GRID1(I,J)=STP
                ELSE
                   GRID1(I,J)=D00
                ENDIF
                ENDIF
               ENDDO
            ENDDO
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(990))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(990))
! $omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
           ENDIF

!Effective Layer Supercell Parameter
            IF (IGET(991)>0) THEN
            DO J=JSTA,JEND
               DO I=1,IM
                IF (ESHR(I,J)<10.) THEN
                   ESHRtmp=D00
                ELSEIF (ESHR(I,J)>20.0) THEN
                   ESHRtmp=1
                ELSE
                   ESHRtmp=(ESHR(I,J)/20.)
                ENDIF
                IF (MUCIN(I,J)>-40.) THEN
                   MUCINtmp=1.0
                ELSE
                   MUCINtmp=(-40./MUCIN(I,J))
                ENDIF
                STP=(MUCAPE(I,J)/D1000)*(ESRH(I,J)/50.)*&
                        ESHRtmp*MUCINtmp
                GRID1(I,J) = spval
                IF(T1D(I,J) < spval) THEN
                IF (STP>0) THEN
                   GRID1(I,J)=STP
                ELSE
                   GRID1(I,J)=D00
                ENDIF
                ENDIF
               ENDDO
            ENDDO
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(991))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(991))
! $omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
           ENDIF

!Mixed Layer (100 mb) Virtual LFC

           IF (IGET(992)>0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               EGRID1(I,J) = -H99999
               EGRID2(I,J) = -H99999
               EGRID3(I,J) = -H99999
               EGRID4(I,J) = -H99999
               EGRID5(I,J) = -H99999
               EGRID6(I,J) = -H99999
               EGRID7(I,J) = -H99999
               EGRID8(I,J) = -H99999
               LB2(I,J)  = (LVLBND(I,J,1) + LVLBND(I,J,2) +           &
                            LVLBND(I,J,3))/3
               P1D(I,J)  = (PBND(I,J,1) + PBND(I,J,2) + PBND(I,J,3))/3
               T1D(I,J)  = (TVIRTUAL(TBND(I,J,1),QBND(I,J,1)) +       &
                            TVIRTUAL(TBND(I,J,2),QBND(I,J,2)) +       &
                            TVIRTUAL(TBND(I,J,3),QBND(I,J,3)))/3
               Q1D(I,J)  = (QBND(I,J,1) + QBND(I,J,2) + QBND(I,J,3))/3
             ENDDO
           ENDDO

             DPBND = 0.
             ITYPE = 2
! EGRID3 is Virtual LFC
             CALL CALCAPE2(ITYPE,DPBND,P1D,T1D,Q1D,LB2,            &
                           EGRID1,EGRID2,EGRID3,EGRID4,EGRID5,     &
                           EGRID6,EGRID7,EGRID8)

             GRID1=spval
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(T1D(I,J) < spval) GRID1(I,J) = EGRID3(I,J)
               ENDDO
             ENDDO
             CALL BOUND(GRID1,D00,H99999)
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(992))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(992))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
           ENDIF   !992


           IF (IGET(763)>0) THEN
!$omp parallel do private(i,j)
! EGRID3 is Virtual LFC
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = Q1D(I,J)
               ENDDO
             ENDDO
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(763))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(763))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
           ENDIF  

!Hail parameter
            IF (IGET(993)>0) THEN
            GRID1=spval
            DO J=JSTA,JEND
               DO I=1,IM
               LAPSE=-((T700(I,J)-T500(I,J))/((Z700(I,J)-Z500(I,J))))
                SHIP=(MUCAPE(I,J)*D1000*MUQ1D(I,J)*LAPSE*(T500(I,J)-K2C)*FSHR(I,J))/HCONST
                IF (MUCAPE(I,J)<1300.)THEN
                   SHIP=SHIP*(MUCAPE(I,J)/1300.)
                ENDIF
                IF (LAPSE < 5.8)THEN
                   SHIP=SHIP*(LAPSE/5.8)
                ENDIF
                IF (FREEZELVL(I,J) < 2400.)THEN
                   SHIP=SHIP*(FREEZELVL(I,J)/2400.)
                ENDIF
                GRID1(I,J)=SHIP
               ENDDO
            ENDDO
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(993))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(993))
! $omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
           ENDIF

        ENDIF   !END RTMA BLOCK


!    Critical Angle

            IF (IGET(957)>0) THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                    IF(T1D(I,J) < spval ) GRID1(I,J) = CANGLE(I,J)   
           !         IF(EGRID1(I,J)<100. .OR. EGRID2(I,J)>-250.) THEN
           !           GRID1(I,J) = 0.
           !         ENDIF
                 ENDDO
               ENDDO
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(957))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(957))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
            ENDIF   !957

!    Dendritic Layer Depth, -17C < T < -12C

            IF (IGET(955)>0) THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval ) GRID1(I,J) = EGRID7(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(955))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(955))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
            ENDIF   !955

!    Enhanced Stretching Potential

            IF (IGET(956)>0) THEN
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval ) GRID1(I,J) = EGRID8(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(956))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(956))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
            ENDIF   !956

!    Downdraft CAPE

!           ITYPE = 1
!           DO J=JSTA,JEND
!           DO I=1,IM
!               LB2(I,J)  = (LVLBND(I,J,1) + LVLBND(I,J,2) +           &
!                            LVLBND(I,J,3))/3
!               P1D(I,J)  = (PBND(I,J,1) + PBND(I,J,2) + PBND(I,J,3))/3
!               T1D(I,J)  = (TBND(I,J,1) + TBND(I,J,2) + TBND(I,J,3))/3
!               Q1D(I,J)  = (QBND(I,J,1) + QBND(I,J,2) + QBND(I,J,3))/3
!             ENDDO
!           ENDDO

!           DPBND = 400.E2
!           CALL CALCAPE2(ITYPE,DPBND,P1D,T1D,Q1D,LB2,            &
!                         EGRID1,EGRID2,EGRID3,EGRID4,EGRID5,     &
!                         EGRID6,EGRID7,EGRID8)

           IF (IGET(954)>0) THEN
               GRID1 = spval
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                 DO I=1,IM
                  IF(T1D(I,J) < spval) GRID1(I,J) = -EGRID6(I,J)
                 ENDDO
              ENDDO
               CALL BOUND(GRID1,D00,H99999)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(954))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(954))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif

           ENDIF   !954

       if (allocated(ushr1)) deallocate(ushr1)
       if (allocated(vshr1)) deallocate(vshr1)
       if (allocated(ushr6)) deallocate(ushr6)
       if (allocated(vshr6)) deallocate(vshr6)
       if (allocated(ust))   deallocate(ust)
       if (allocated(vst))   deallocate(vst)
       if (allocated(heli))  deallocate(heli)
       if (allocated(llow))  deallocate(llow)
       if (allocated(lupp))  deallocate(lupp)
       if (allocated(llow1))  deallocate(llow1)
       if (allocated(lupp1))  deallocate(lupp1)
       if (allocated(cangle))deallocate(cangle)
       if (allocated(effust))deallocate(effust)
       if (allocated(effvst))deallocate(effvst)
       if (allocated(eshr))  deallocate(eshr)
       if (allocated(uvect)) deallocate(uvect)
       if (allocated(vvect)) deallocate(vvect)
       if (allocated(esrh))  deallocate(esrh)
       if (allocated(htsfc)) deallocate(htsfc)
       if (allocated(fshr))  deallocate(fshr)
       ENDIF

      if (allocated(pbnd))   deallocate(pbnd)
      if (allocated(tbnd))   deallocate(tbnd)
      if (allocated(qbnd))   deallocate(qbnd)
      if (allocated(ubnd))   deallocate(ubnd)
      if (allocated(vbnd))   deallocate(vbnd)
      if (allocated(rhbnd))  deallocate(rhbnd)
      if (allocated(wbnd))   deallocate(wbnd)
      if (allocated(lvlbnd)) deallocate(lvlbnd)
      if (allocated(lb2))    deallocate(lb2)
!    
!
! RELATIVE HUMIDITY WITH RESPECT TO PRECIPITABLE WATER
       IF (IGET(749)>0) THEN
          CALL CALRH_PW(GRID1(1,jsta))
          if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(749))
!$omp parallel do private(i,j,jj)
           do j=1,jend-jsta+1
             jj = jsta+j-1
             do i=1,im
               datapd(i,j,cfld) = GRID1(i,jj)
             enddo
           enddo
          endif
       ENDIF       

 
!     END OF ROUTINE.
!     
      RETURN
      END
