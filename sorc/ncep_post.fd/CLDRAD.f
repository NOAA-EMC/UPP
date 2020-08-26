      SUBROUTINE CLDRAD
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CLDRAD       POST SNDING/CLOUD/RADTN FIELDS
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-08-30       
!     
! ABSTRACT:  THIS ROUTINE COMPUTES/POSTS SOUNDING, CLOUD 
!   RELATED, AND RADIATION FIELDS.  UNDER THE HEADING OF 
!   SOUNDING FIELDS FALL THE THREE ETA MODEL LIFTED INDICES,
!   CAPE, CIN, AND TOTAL COLUMN PRECIPITABLE WATER.
!
!   THE THREE ETA MODEL LIFTED INDICES DIFFER ONLY IN THE
!   DEFINITION OF THE PARCEL TO LIFT.  ONE LIFTS PARCELS FROM
!   THE LOWEST ABOVE GROUND ETA LAYER.  ANOTHER LIFTS MEAN 
!   PARCELS FROM ANY OF NBND BOUNDARY LAYERS (SEE SUBROUTINE
!   BNDLYR).  THE FINAL TYPE OF LIFTED INDEX IS A BEST LIFTED
!   INDEX BASED ON THE NBND BOUNDARY LAYER LIFTED INDICES.
!
!   TWO TYPES OF CAPE/CIN ARE AVAILABLE.  ONE IS BASED ON PARCELS
!   IN THE LOWEST ETA LAYER ABOVE GROUND.  THE OTHER IS BASED 
!   ON A LAYER MEAN PARCEL IN THE N-TH BOUNDARY LAYER ABOVE 
!   THE GROUND.  SEE SUBROUTINE CALCAPE FOR DETAILS.
!
!   THE CLOUD FRACTION AND LIQUID CLOUD WATER FIELDS ARE DIRECTLY
!   FROM THE MODEL WITH MINIMAL POST PROCESSING.  THE LIQUID 
!   CLOUD WATER, 3-D CLOUD FRACTION, AND TEMPERATURE TENDENCIES
!   DUE TO PRECIPITATION ARE NOT POSTED IN THIS ROUTINE.  SEE
!   SUBROUTINE ETAFLD FOR THESE FIELDS.  LIFTING CONDENSATION
!   LEVEL HEIGHT AND PRESSURE ARE COMPUTED AND POSTED IN
!   SUBROUTINE MISCLN.  
!
!   THE RADIATION FIELDS POSTED BY THIS ROUTINE ARE THOSE COMPUTED
!   DIRECTLY IN THE MODEL.
!     
! PROGRAM HISTORY LOG:
!   93-08-30  RUSS TREADON
!   94-08-04  MICHAEL BALDWIN - ADDED OUTPUT OF INSTANTANEOUS SFC
!                               FLUXES OF NET SW AND LW DOWN RADIATION
!   97-04-25  MICHAEL BALDWIN - FIX PDS FOR PRECIPITABLE WATER
!   97-04-29  GEOFF MANIKIN - MOVED CLOUD TOP TEMPS CALCULATION
!                               TO THIS SUBROUTINE.  CHANGED METHOD
!                               OF DETERMINING WHERE CLOUD BASE AND
!                               TOP ARE FOUND AND ADDED HEIGHT OPTION
!                               FOR TOP AND BASE.
!   98-04-29  GEOFF MANIKIN - CHANGED VALUE FOR CLOUD BASE/TOP PRESSURES
!                               AND HEIGHTS FROM SPVAL TO -500
!   98-06-15  T BLACK       - CONVERSION FROM 1-D TO 2-D
!   98-07-17  MIKE BALDWIN  - REMOVED LABL84
!   00-01-04  JIM TUCCILLO  - MPI VERSION
!   00-02-22  GEOFF MANIKIN - CHANGED VALUE FOR CLOUD BASE/TOP PRESSURES
!                               AND HEIGHTS FROM SPVAL TO -500 (WAS NOT IN
!                               PREVIOUS IBM VERSION)
!   01-10-22  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-01-15  MIKE BALDWIN - WRF VERSION
!   05-01-06  H CHUANG - ADD VARIOUS CLOUD FIELDS
!   05-07-07  BINBIN ZHOU - ADD RSM MODEL
!   05-08-30  BINBIN ZHOU - ADD CEILING and FLIGHT CONDITION RESTRICTION
!   10-09-09  GEOFF MANIKIN - REVISED CALL TO CALCAPE
!   11-02-06  Jun Wang - ADD GRIB2 OPTION
!   11-12-14  SARAH LU - ADD AEROSOL OPTICAL PROPERTIES
!   11-12-16  SARAH LU - ADD AEROSOL 2D DIAG FIELDS
!   11-12-23  SARAH LU - CONSOLIDATE ALL GOCART FIELDS TO BLOCK 4
!   11-12-23  SARAH LU - ADD AOD AT ADDITIONAL CHANNELS
!   12-04-03  Jun Wang - Add lftx and GFS convective cloud cover for grib2
!   13-05-06  Shrinivas Moorthi - Add cloud condensate to total precip water
!   13-12-23  LU/Wang  - READ AEROSOL OPTICAL PROPERTIES LUTS to compute dust aod,
!                        non-dust aod, and use geos5 gocart LUTS
!   15-??-??  S. Moorthi - threading, optimization, local dimension
!   19-07-24  Li(Kate) Zhang Merge and update ARAH Lu's work from NGAC into FV3-Chem
!   19-10-30  Bo CUI - Remove "GOTO" statement
!   20-03-25  Jesse Meng - remove grib1
!     
! USAGE:    CALL CLDRAD
!   INPUT ARGUMENT LIST:
!
!   OUTPUT ARGUMENT LIST: 
!     NONE
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       NONE
!     LIBRARY:
!       COMMON   - RQSTFLD
!                  CTLBLK
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : IBM SP
!$$$  
!
      use vrbls4d, only: DUST,SUSO, SALT, SOOT, WASO
      use vrbls3d, only: QQW, QQR, T, ZINT, CFR, QQI, QQS, Q, EXT, ZMID,PMID,&
                         PINT, DUEM, DUSD, DUDP, DUWT, DUSV, SSEM, SSSD,SSDP,&
                         SSWT, SSSV, BCEM, BCSD, BCDP, BCWT, BCSV, OCEM,OCSD,&
                         OCDP, OCWT, OCSV, SCA, ASY,CFR_RAW
      use vrbls2d, only: CLDEFI, CFRACL, AVGCFRACL, CFRACM, AVGCFRACM, CFRACH,&
                         AVGCFRACH, AVGTCDC, NCFRST, ACFRST, NCFRCV, ACFRCV,  &
                         HBOT, HBOTD, HBOTS, HTOP, HTOPD, HTOPS,  FIS, PBLH,  &
                         PBOT, PBOTL, PBOTM, PBOTH, CNVCFR, PTOP, PTOPL,      &
                         PTOPM, PTOPH, TTOPL, TTOPM, TTOPH, PBLCFR, CLDWORK,  &
                         ASWIN, AUVBIN, AUVBINC, ASWIN, ASWOUT,ALWOUT, ASWTOA,&
                         RLWTOA, CZMEAN, CZEN, RSWIN, ALWIN, ALWTOA, RLWIN,   &
                         SIGT4, RSWOUT, RADOT, RSWINC, ASWINC, ASWOUTC,       &
                         ASWTOAC, ALWOUTC, ASWTOAC, AVISBEAMSWIN,             &
                         AVISDIFFSWIN, ASWINTOA, ASWINC, ASWTOAC, AIRBEAMSWIN,&
                         AIRDIFFSWIN, DUSMASS, DUSMASS25, DUCMASS, DUCMASS25, &
                         ALWINC, ALWTOAC, SWDDNI, SWDDIF, SWDNBC, SWDDNIC,    &
                         SWDDIFC, SWUPBC, LWDNBC, LWUPBC, SWUPT,              &
                         TAOD5502D, AERSSA2D, AERASY2D, MEAN_FRP, LWP, IWP,   &
                         TAOD5502D, AERSSA2D, AERASY2D, AVGCPRATE,            &
                         DUSTCB,SSCB,BCCB,OCCB,SULFCB,DUSTPM,SSPM
      use masks,    only: LMH, HTM
      use params_mod, only: TFRZ, D00, H99999, QCLDMIN, SMALL, D608, H1, ROG, &
                            GI, RD, QCONV, ABSCOEFI, ABSCOEF, STBOL, PQ0, A2, &
                            A3, A4
      use ctlblk_mod, only: JSTA, JEND, SPVAL, MODELNAME, GRIB, CFLD,DATAPD,  &
                            FLD_INFO, AVRAIN, THEAT, IFHR, IFMIN, AVCNVC,     &
                            TCLOD, ARDSW, TRDSW, ARDLW, NBIN_DU, TRDLW, IM,   &
                            NBIN_SS, NBIN_OC, NBIN_BC, NBIN_SU, DTQ2,         &
                            JM, LM, gocart_on, me
      use rqstfld_mod, only: IGET, ID, LVLS, IAVBLFLD
      use gridspec_mod, only: dyval, gridtype
      use cmassi_mod,  only: TRAD_ice
      use machine_post,     only: kind_phys
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     SET CELSIUS TO KELVIN CONVERSION.
      REAL,PARAMETER :: C2K=273.15, PTOP_LOW=64200., PTOP_MID=35000.,        &
                        PTOP_HIGH=15000.
!     
!     DECLARE VARIABLES.
!     
!     LOGICAL,dimension(im,jm) ::  NEED
      INTEGER :: lcbot,lctop,jc,ic  !bsf
      INTEGER,dimension(im,jsta:jend) :: IBOTT, IBOTCu, IBOTDCu, IBOTSCu, IBOTGr,   &
                                         ITOPT, ITOPCu, ITOPDCu, ITOPSCu, ITOPGr
      REAL,dimension(im,jm)           :: GRID1
      REAL,dimension(im,jsta:jend)    :: GRID2, EGRID1, EGRID2, EGRID3,      &
                                         CLDP, CLDZ, CLDT, CLDZCu
      REAL,dimension(lm)       :: RHB, watericetotal, pabovesfc
      REAL   :: watericemax, wimin, zcldbase, zcldtop, zpbltop,              &
                rhoice, coeffp, exponfp, const1, cloud_def_p,                &
                pcldbase, rhoair, vovermd, concfp, betav,                    &
                vertvis, tx, tv, pol, esx, es, e, zsf, zcld, frac
      integer   nfog, nfogn(7),npblcld,nlifr, k1, k2, ll, ii, ib, n, jj,     &
                NUMR, NUMPTS
      real,dimension(lm)       :: cldfra
      real                     :: ceiling_thresh_cldfra, cldfra_max, zceil
      real,dimension(im,jm)    :: ceil
!     B ZHOU: For aviation:
      REAL, dimension(im,jsta:jend) :: TCLD, CEILING
      real   CU_ir(LM), q_conv   !bsf
!jw
      integer I,J,L,K,IBOT,ITCLOD,LBOT,LTOP,ITRDSW,ITRDLW,        &
              LLMH,ITHEAT,IFINCR,ITYPE,ITOP,NUM_THICK
      real    DPBND,RRNUM,QCLD,RSUM,TLMH,FACTRS,FACTRL,DP,        &
              OPDEPTH, TMP,QSAT,RHUM,TCEXT,DELZ,DELY,DY_m
!
      real    FULL_CLD(IM,JM)   !-- Must be dimensioned for the full domain
!
      real    dummy(IM,jsta:jend)
      integer idummy(IM,jsta:jend)
!
!     --- Revision added for GOCART ---

!!    GOCART aerosol optical data from GSFC Mie code calculations,
!!    mapped to 7 channels: 0.34, 0.44, 0.55, 0.66, 0.86, 1.63, 11.1 micron
!!    data  wvnum1/0.338, 0.430, 0.545, 0.62, 0.841, 1.628, 11.0/
!!    data  wvnum2/0.342, 0.450, 0.565, 0.67, 0.876, 1.652, 11.2/

      integer, parameter :: KRHLEV = 36 ! num of rh levels for rh-dep components
      integer, parameter :: KCM1 = 5    ! num of rh independent aer species
      integer, parameter :: KCM2 = 5    ! num of rh dependent aer species
      integer, parameter :: NBDSW = 7   ! total num of sw bands
      integer, parameter :: NOAER = 20  ! unit for LUTs file
      integer, parameter :: nAero=KCM2  ! num of aer species in LUTs
      CHARACTER          :: AerosolName(KCM2)*4, AerosolName_rd*4, aerosol_file*30
      CHARACTER          :: AerName_rd*4, AerOpt*3

!   - aerosol optical properties: mass extinction efficiency
      REAL, ALLOCATABLE  :: extrhd_DU(:,:,:), extrhd_SS(:,:,:), &
     &                      extrhd_SU(:,:,:), extrhd_BC(:,:,:), &
     &                      extrhd_OC(:,:,:)

!   - aerosol optical properties: mass scattering efficienc
      REAL, ALLOCATABLE  :: scarhd_DU(:,:,:), scarhd_SS(:,:,:), &
     &                      scarhd_SU(:,:,:), scarhd_BC(:,:,:), &
     &                      scarhd_OC(:,:,:)

!   - aerosol optical properties: asymmetry factor
      REAL, ALLOCATABLE  :: asyrhd_DU(:,:,:), asyrhd_SS(:,:,:), &
     &                      asyrhd_SU(:,:,:), asyrhd_BC(:,:,:), &
     &                      asyrhd_OC(:,:,:)

!   - aerosol optical properties: single scatter albedo
      REAL, ALLOCATABLE  :: ssarhd_DU(:,:,:), ssarhd_SS(:,:,:), &
     &                      ssarhd_SU(:,:,:), ssarhd_BC(:,:,:), &
     &                      ssarhd_OC(:,:,:)

!  --- aerosol optical properties mapped onto specified spectral bands
!   - relative humidity independent aerosol optical properties: du
      real (kind=kind_phys)  :: extrhi(KCM1,NBDSW)         ! extinction coefficient

!   - relative humidity dependent aerosol optical properties: oc, bc, su, ss001-005
      real (kind=kind_phys)  :: extrhd(KRHLEV,KCM2,NBDSW)  ! extinction coefficient
!
      REAL,dimension(im,jsta:jend)  :: P1D,T1D,Q1D,EGRID4
!     REAL, allocatable  :: RH3D(:,:,:)                     ! RELATIVE HUMIDITY
      real,    allocatable:: rdrh(:,:,:)
      integer, allocatable :: ihh(:,:,:)
      REAL                 :: rh3d, DRH0, DRH1, EXT01, EXT02,SCA01,ASY01
      INTEGER              :: IH1, IH2
      INTEGER            :: IOS, INDX, ISSAM, ISSCM, ISUSO, IWASO, ISOOT, NBIN
      REAL               :: CCDRY, CCWET, SSAM, SSCM
      REAL,dimension(im,jsta:jend) :: AOD_DU, AOD_SS, AOD_SU, AOD_OC, AOD_BC, AOD
      REAL,dimension(im,jsta:jend) :: SCA_DU, SCA_SS, SCA_SU, SCA_OC,SCA_BC, SCA2D
      REAL,dimension(im,jsta:jend) :: ASY_DU, ASY_SS, ASY_SU, ASY_OC, ASY_BC,ASY2D
      REAL,dimension(im,jsta:jend) :: ANGST, AOD_440, AOD_860      ! FORANGSTROM EXPONENT
      REAL               :: ANG1, ANG2
      INTEGER            :: INDX_EXT(nAero), INDX_SCA(nAero)
      LOGICAL            :: LAEROPT, LEXT, LSCA, LASY
      LOGICAL            :: LAERSMASS
      REAL, allocatable  :: fPM25_DU(:),fPM25_SS(:)
      REAL, allocatable, dimension(:,:) :: RHOsfc, smass_du_cr,smass_du_fn, &
     &                      smass_ss_cr, smass_ss_fn, smass_oc,smass_bc,    &
     &                      smass_su, smass_cr, smass_fn
      real               :: rPM, dmass
      real (kind=kind_phys), dimension(KRHLEV) :: rhlev
      data  rhlev (:)/  .0, .05, .10, .15, .20, .25, .30, .35,               &
     &                 .40, .45, .50, .55, .60, .65, .70, .75,               &
     &                 .80, .81, .82, .83, .84, .85, .86, .87,               &
     &                 .88, .89, .90, .91, .92, .93, .94, .95,               &
     &                 .96, .97, .98, .99/
!
      data AerosolName    /'DUST', 'SALT', 'SUSO', 'SOOT', 'WASO'/
!     INDEX FOR TOTAL AND SPECIATED AEROSOLS (DU, SS, SU, OC, BC)
      data INDX_EXT       / 610, 611, 612, 613, 614  /
      data INDX_SCA       / 651, 652, 653, 654, 655  /
!     
!
!*************************************************************************
!     START CLDRAD HERE.
!     
!***  BLOCK 1.  SOUNDING DERIVED FIELDS.
!     
!     ETA SURFACE TO 500MB LIFTED INDEX.  TO BE CONSISTENT WITH THE
!     LFM AND NGM POSTING WE ADD 273.15 TO THE LIFTED INDEX
! GSM     WILL NOT ADD 273 TO VALUE FOR RAPID REFRESH TO BE
!           CONSISTENT WITH RUC
!
!     THE BEST (SIX LAYER) AND BOUNDARY LAYER LIFTED INDICES ARE
!     COMPUTED AND POSTED IN SUBROUTINE MISCLN.
!
      IF (IGET(030).GT.0.OR.IGET(572)>0) THEN
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            EGRID1(I,J) = SPVAL
          ENDDO
        ENDDO
!
        CALL OTLIFT(EGRID1)
!
        IF(MODELNAME == 'RAPR') THEN
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              IF(EGRID1(I,J) < SPVAL) GRID1(I,J) = EGRID1(I,J)
            ENDDO
          ENDDO
        ELSE
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              IF(EGRID1(I,J) < SPVAL) GRID1(I,J) = EGRID1(I,J) + TFRZ
            ENDDO
          ENDDO
        ENDIF
!
        if(IGET(030) > 0) then
          if(grib == "grib2" )then
            cfld = cfld+1
            fld_info(cfld)%ifld = IAVBLFLD(IGET(030))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
          endif
        endif
!for GFS
        if(IGET(572) > 0) then
          if(grib == "grib2" )then
            cfld = cfld+1
            fld_info(cfld)%ifld = IAVBLFLD(IGET(572))
!           where(GRID1 /= SPVAL) GRID1 = GRID1-TFRZ
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                if (grid1(i,jj) /= spval) grid1(i,jj) = grid1(i,jj) - tfrz
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
          endif
        endif

      ENDIF
!
!     SOUNDING DERIVED AREA INTEGRATED ENERGIES - CAPE AND CIN.
!       THIS IS THE SFC-BASED CAPE/CIN (lowest 70 mb searched)
!
!           CONVECTIVE AVAILABLE POTENTIAL ENERGY.
      IF ((IGET(032) > 0))THEN
! dong add missing value for cape 
        GRID1 = spval
        IF ( (LVLS(1,IGET(032)).GT.0) )THEN
          ITYPE  = 1
          DPBND  = 10.E2
          dummy  = 0.
          idummy = 0
          CALL CALCAPE(ITYPE,DPBND,dummy,dummy,dummy,idummy,EGRID1,EGRID2, &
                       EGRID3,dummy,dummy)
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              IF(FIS(I,J) < SPVAL) GRID1(I,J) = EGRID1(I,J)
            ENDDO
          ENDDO
          CALL BOUND(GRID1,D00,H99999)
          if(grib == "grib2" )then
            cfld = cfld+1
            fld_info(cfld)%ifld = IAVBLFLD(IGET(032))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
          endif
        END IF
      END IF
!
!           CONVECTIVE INHIBITION.     
      IF ((IGET(107) > 0))THEN
! dong add missing value for cin
        GRID1 = spval
        IF ( (LVLS(1,IGET(107)) > 0) )THEN
          IF ((IGET(032) > 0))THEN
            IF ( (LVLS(1,IGET(032)) > 0) )THEN
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                DO I=1,IM
                  IF(FIS(I,J) < SPVAL) GRID1(I,J) = - EGRID2(I,J)
                ENDDO
              ENDDO
            END IF
          ELSE
            ITYPE  = 1
            DPBND  = 10.E2
            dummy  = 0.
            idummy = 0
            CALL CALCAPE(ITYPE,DPBND,dummy,dummy,dummy,idummy,EGRID1,EGRID2, &
                         EGRID3,dummy,dummy)
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                IF(FIS(I,J) < SPVAL) GRID1(I,J) = - EGRID2(I,J)
              ENDDO
            ENDDO
          END IF   
          CALL BOUND(GRID1,D00,H99999)
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              IF(FIS(I,J) < SPVAL) GRID1(I,J) = - GRID1(I,J)
            ENDDO
          ENDDO
          if(grib == "grib2" )then
            cfld = cfld+1
            fld_info(cfld)%ifld = IAVBLFLD(IGET(107))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
          endif
        END IF ! end for lvls(107)
      END IF ! end of iget(107)	 
      
!!!=======================================================================
!
!     TOTAL COLUMN PRECIPITABLE WATER (SPECIFIC HUMIDITY).
      IF (IGET(080) > 0) THEN
! dong 
         GRID1 = spval
         CALL CALPW(GRID1(1,jsta),1)
         ID(1:25) = 0
          DO J=JSTA,JEND
            DO I=1,IM
              IF(FIS(I,J) >= SPVAL) GRID1(I,J)=spval
            END DO
          END DO
        CALL BOUND(GRID1,D00,H99999)
        if(grib == "grib2" )then
          cfld = cfld + 1
          fld_info(cfld)%ifld = IAVBLFLD(IGET(080))
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
!     E. James - 8 Dec 2017
!     TOTAL COLUMN AOD (TAOD553D FROM HRRR-SMOKE)
!
      IF (IGET(735) > 0) THEN
         CALL CALPW(GRID1(1,jsta),19)
         ID(1:25) = 0
         CALL BOUND(GRID1,D00,H99999)
        if(grib == "grib2" )then
          cfld = cfld + 1
          fld_info(cfld)%ifld = IAVBLFLD(IGET(735))
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
!     E. James - 8 Dec 2017
!     TOTAL COLUMN FIRE SMOKE (tracer_1a FROM HRRR-SMOKE)
!
      IF (IGET(736) > 0) THEN
         CALL CALPW(GRID1(1,jsta),18)
         ID(1:25) = 0
         CALL BOUND(GRID1,D00,H99999)
        if(grib == "grib2" )then
          cfld = cfld + 1
          fld_info(cfld)%ifld = IAVBLFLD(IGET(736))
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
!     TOTAL COLUMN CLOUD WATER
      IF (IGET(200) > 0 .or. IGET(575) > 0) THEN 
       IF (MODELNAME == 'RAPR') THEN
          DO J=JSTA,JEND
            DO I=1,IM
              GRID1(I,J) = LWP(I,J)/1000.0 ! use WRF-diagnosed value
            ENDDO
          ENDDO
       ELSE
        CALL CALPW(GRID1(1,jsta),2)
        IF(MODELNAME == 'GFS')then
! GFS combines cloud water and cloud ice, hoping to seperate them next implementation    
          CALL CALPW(GRID2(1,jsta),3)
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              GRID1(I,J) = GRID1(I,J) + GRID2(I,J)
            ENDDO
          ENDDO
        END IF ! GFS
       END IF ! RAPR
  
        ID(1:25) = 0
        ID(02)   = 129      !--- Parameter Table 129, PDS Octet 4 = 129)
        CALL BOUND(GRID1,D00,H99999)
        if(IGET(200) > 0) then
          if(grib == "grib2" )then
            cfld = cfld + 1
            fld_info(cfld)%ifld = IAVBLFLD(IGET(200))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
          endif
        endif
        if(iget(575) > 0) then
          if(grib == "grib2" )then
            cfld = cfld + 1
            fld_info(cfld)%ifld = IAVBLFLD(IGET(575))
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
!
!     TOTAL COLUMN CLOUD ICE
      IF (IGET(201) > 0) THEN
       IF (MODELNAME == 'RAPR') THEN
          DO J=JSTA,JEND
            DO I=1,IM
              GRID1(I,J) = IWP(I,J)/1000.0 ! use WRF-diagnosed value
            ENDDO
          ENDDO
       ELSE
         CALL CALPW(GRID1(1,jsta),3)
       END IF
         ID(1:25) = 0
         ID(02)   = 129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
        if(grib == "grib2" )then
          cfld = cfld + 1
          fld_info(cfld)%ifld = IAVBLFLD(IGET(201))
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
!     TOTAL COLUMN RAIN 
      IF (IGET(202) > 0) THEN
         CALL CALPW(GRID1(1,jsta),4)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(202))
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
!     TOTAL COLUMN SNOW 
      IF (IGET(203) > 0) THEN
         CALL CALPW(GRID1(1,jsta),5)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(203))
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
! SRD
!     TOTAL COLUMN GRAUPEL
      IF (IGET(428) > 0) THEN
         CALL CALPW(GRID1(1,jsta),16)
         ID(1:25)=0
!         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(428))
!$omp parallel do private(i,j,jj)
          do j=1,jend-jsta+1
            jj = jsta+j-1
            do i=1,im
              datapd(i,j,cfld) = GRID1(i,jj)
            enddo
          enddo
        endif
      ENDIF
! SRD

!     TOTAL COLUMN CONDENSATE 
      IF (IGET(204) > 0) THEN
         CALL CALPW(GRID1(1,jsta),6)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(204))
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
!     TOTAL COLUMN SUPERCOOLED (<0C) LIQUID WATER 
      IF (IGET(285) > 0) THEN
         CALL CALPW(GRID1(1,jsta),7)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(285))
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
!     TOTAL COLUMN MELTING (>0C) ICE
      IF (IGET(286) > 0) THEN
         CALL CALPW(GRID1(1,jsta),8)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(286))
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
!     TOTAL COLUMN SHORT WAVE T TENDENCY
      IF (IGET(290) > 0) THEN
         CALL CALPW(GRID1(1,jsta),9)
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(290))
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
!     TOTAL COLUMN LONG WAVE T TENDENCY
      IF (IGET(291) > 0) THEN
         CALL CALPW(GRID1(1,jsta),10)
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(291))
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
!     TOTAL COLUMN GRID SCALE LATENT HEATING (TIME AVE)
      IF (IGET(292) > 0) THEN
         CALL CALPW(GRID1(1,jsta),11)
         IF(AVRAIN > 0.)THEN
           RRNUM = 1./AVRAIN
         ELSE
           RRNUM = 0.
         ENDIF
!$omp  parallel do
         DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = GRID1(I,J)*RRNUM
           ENDDO
         ENDDO
         ID(1:25)=0
         ITHEAT     = NINT(THEAT)
         IF (ITHEAT .NE. 0) THEN
          IFINCR     = MOD(IFHR,ITHEAT)
         ELSE
          IFINCR=0
         END IF
         ID(19) = IFHR
         IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20) = 3
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITHEAT
         ELSE
          ID(18) = IFHR-IFINCR
         ENDIF
         IF(IFMIN .GE. 1)ID(18)=ID(18)*60
         IF (ID(18).LT.0) ID(18) = 0
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(292))
            if(ITHEAT>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
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
!     TOTAL COLUMN CONVECTIVE LATENT HEATING (TIME AVE)
      IF (IGET(293) > 0) THEN
         CALL CALPW(GRID1(1,jsta),12)
         IF(AVRAIN > 0.)THEN
           RRNUM = 1./AVCNVC
         ELSE
           RRNUM = 0.
         ENDIF
!$omp  parallel do
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = GRID1(I,J)*RRNUM
         ENDDO
         ENDDO
         ID(1:25)=0
         ITHEAT     = NINT(THEAT)
         IF (ITHEAT .NE. 0) THEN
          IFINCR     = MOD(IFHR,ITHEAT)
         ELSE
          IFINCR=0
         END IF
         ID(19) = IFHR
         IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20) = 3
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITHEAT
         ELSE
          ID(18) = IFHR-IFINCR
         ENDIF
         IF(IFMIN .GE. 1)ID(18)=ID(18)*60
         IF (ID(18).LT.0) ID(18) = 0
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(293))
            if(ITHEAT>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
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
!     TOTAL COLUMN moisture convergence
      IF (IGET(295).GT.0) THEN
         CALL CALPW(GRID1(1,jsta),13)
         ID(1:25)=0
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(295))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF
!
!     TOTAL COLUMN RH
      IF (IGET(312).GT.0) THEN
         CALL CALPW(GRID1(1,jsta),14)
         ID(1:25)=0
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(312))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF
!
!     TOTAL COLUMN OZONE
      IF (IGET(299) > 0) THEN
         CALL CALPW(GRID1(1,jsta),15)
         ID(1:25)=0
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(299))
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
!     BOTTOM AND/OR TOP OF SUPERCOOLED (<0C) LIQUID WATER LAYER
      IF (IGET(287).GT.0 .OR. IGET(288).GT.0) THEN
         DO J=JSTA,JEND
            DO I=1,IM
               GRID1(I,J)=-5000.
               GRID2(I,J)=-5000.
!-- Search for the base first, then look for the top if supercooled liquid exists
               LBOT=0
               LM=NINT(LMH(I,J))
               DO L=LM,1,-1
                  QCLD=QQW(I,J,L)+QQR(I,J,L)
                  IF (QCLD.GE.QCLDmin .AND. T(I,J,L).LT.TFRZ) THEN
                     LBOT=L
                     EXIT
                  ENDIF
               ENDDO    !--- End L loop
               IF (LBOT .GT. 0) THEN
!-- Supercooled liquid exists, so get top & bottom heights.  In this case,
!   be conservative and select the lower interface height at the bottom of the
!   layer and the top interface height at the top of the layer.
                  GRID1(I,J)=ZINT(I,J,LBOT+1)
                  DO L=1,LM
                     QCLD=QQW(I,J,L)+QQR(I,J,L)
                     IF (QCLD.GE.QCLDmin .AND. T(I,J,L).LT.TFRZ) THEN
                        LTOP=L
                        EXIT
                     ENDIF
                  ENDDO    !--- End L loop
                  LTOP=MIN(LBOT,LTOP)
                  GRID2(I,J)=ZINT(I,J,LTOP)
               ENDIF    !--- End IF (LBOT .GT. 0)
            ENDDO       !--- End I loop
         ENDDO          !--- End J loop
         IF (IGET(287).GT.0) THEN
           if(grib=="grib2" )then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(287))
             datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
           endif
         ENDIF
         IF (IGET(288).GT.0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J)=GRID2(I,J)
              ENDDO
            ENDDO
           if(grib=="grib2" )then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(288))
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
!     Convective cloud efficiency parameter used in convection ranges
!     from 0.2 (EFIMN in cuparm in model) to 1.0   (Ferrier, Feb '02) 
      IF (IGET(197).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = CLDEFI(I,J)
         ENDDO
         ENDDO
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(197))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF
!
      IF ((MODELNAME=='NMM' .AND. GRIDTYPE=='B') .OR. &
             MODELNAME=='FV3R') THEN
!nmmb_clds1
!   
!-- Initialize low, middle, high, and total cloud cover; 
!   also a method for cloud ceiling height
!
         DO J=JSTA,JEND
           DO I=1,IM
             CFRACL(I,J)=0.
             CFRACM(I,J)=0.
             CFRACH(I,J)=0.
             TCLD(I,J)=0.
           ENDDO
         ENDDO
!
!-- Average cloud fractions over a 10 mi (16.09 km) radius (R), 
!   approximated by a box of the same area = pi*R**2. Final
!   distance (d) is 1/2 of box size, d=0.5*sqrt(pi)*R=14259 m.
!
        if(grib == "grib2" )then
          DY_m=DYVAL*0.1112    !- DY_m in m 
        endif   
        DELY=14259./DY_m
        numr=NINT(DELY)
  !     write (0,*) 'numr,dyval,DY_m=',numr,dyval,DY_m
        DO L=LM,1,-1
          DO J=JSTA,JEND
            DO I=1,IM
              FULL_CLD(I,J)=CFR(I,J,L)    !- 3D cloud fraction (from radiation)
            ENDDO
          ENDDO
          CALL AllGETHERV(FULL_CLD)
          DO J=JSTA,JEND
            DO I=1,IM
              NUMPTS=0
              FRAC=0.
              DO JC=max(1,J-numr),min(JM,J+numr)
                DO IC=max(1,I-numr),min(IM,I+numr)
                  NUMPTS=NUMPTS+1
                  FRAC=FRAC+FULL_CLD(IC,JC)
                ENDDO
              ENDDO
              IF (NUMPTS>0) FRAC=FRAC/REAL(NUMPTS)
              PCLDBASE=PMID(I,J,L)    !-- Using PCLDBASE variable for convenience
              IF (PCLDBASE>=PTOP_LOW) THEN
                CFRACL(I,J)=MAX(CFRACL(I,J),FRAC)
              ELSE IF (PCLDBASE>=PTOP_MID) THEN
                CFRACM(I,J)=MAX(CFRACM(I,J),FRAC)
              ELSE
                CFRACH(I,J)=MAX(CFRACH(I,J),FRAC)
              ENDIF
              TCLD(I,J)=MAX(TCLD(I,J),FRAC)
            ENDDO  ! I
          ENDDO    ! J
        ENDDO      ! L
!end nmmb_clds1
      ELSEIF (MODELNAME=='GFS') THEN
!Initialize for GLOBAL FV3 which has cluod fraction in range from
!0.0 to 1.0
!
!-- Initialize low, middle, high, and total cloud cover;
!   also a method for cloud ceiling height
!
        DO J=JSTA,JEND
          DO I=1,IM
            CFRACL(I,J)=0.
            CFRACM(I,J)=0.
            CFRACH(I,J)=0.
            TCLD(I,J)=0.
          ENDDO
        ENDDO
        DO L=LM,1,-1
          DO J=JSTA,JEND
            DO I=1,IM
              FRAC=CFR(I,J,L) !- 3D cloud fraction at model layers
              PCLDBASE=PMID(I,J,L)    !-- Using PCLDBASE variable for convenience
              IF (PCLDBASE>=PTOP_LOW) THEN
                CFRACL(I,J)=MAX(CFRACL(I,J),FRAC)
              ELSE IF (PCLDBASE>=PTOP_MID) THEN
                CFRACM(I,J)=MAX(CFRACM(I,J),FRAC)
              ELSE
                CFRACH(I,J)=MAX(CFRACH(I,J),FRAC)
              ENDIF
              TCLD(I,J)=MAX(TCLD(I,J),FRAC)
            ENDDO  ! I
          ENDDO    ! J
        ENDDO      ! L
      ENDIF  
!
!***  BLOCK 2.  2-D CLOUD FIELDS.
!
!     LOW CLOUD FRACTION.
      IF (IGET(037) > 0) THEN
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            IF(CFRACL(I,J) < SPVAL) then
              GRID1(I,J) = CFRACL(I,J)*100.
            else
              GRID1(I,J) = spval
            endif
          ENDDO
        ENDDO
       if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(037))
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
!     TIME AVERAGED LOW CLOUD FRACTION.
      IF (IGET(300) > 0) THEN
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            IF(AVGCFRACL(I,J) < SPVAL) then
              GRID1(I,J) = AVGCFRACL(I,J)*100.   
            else
              GRID1(I,J) = spval
            endif
          ENDDO
        ENDDO
        ID(1:25)=0
        ITCLOD     = NINT(TCLOD)
        IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
          IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
        ELSE
          IFINCR     = 0
        endif

        ID(19)  = IFHR
        IF(IFMIN >= 1)ID(19)=IFHR*60+IFMIN  !USE MIN FOR OFF-HR FORECAST
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
           IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
       if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(300))
          if(ITCLOD>0) then
            fld_info(cfld)%ntrange=1
          else
            fld_info(cfld)%ntrange=0
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)
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
!     MIDDLE CLOUD FRACTION.
      IF (IGET(038) > 0) THEN
!       GRID1=SPVAL
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            IF(CFRACM(I,J) < SPVAL) then
              GRID1(I,J) = CFRACM(I,J)*100.
            else
              GRID1(I,J) = spval
            endif
          ENDDO
        ENDDO
       if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(038))
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
!     TIME AVERAGED MIDDLE CLOUD FRACTION.
      IF (IGET(301) > 0) THEN
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            IF(ABS(AVGCFRACM(I,J)-SPVAL)>SMALL)THEN
              GRID1(I,J) = AVGCFRACM(I,J)*100.
            ELSE
              GRID1(I,J) = SPVAL
            END IF 
          ENDDO
        ENDDO
        ID(1:25)=0
        ITCLOD     = NINT(TCLOD)
        IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
          IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
        ELSE
          IFINCR     = 0
        endif

        ID(19)  = IFHR
        IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN  !USE MIN FOR OFF-HR FORECAST
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
           IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(301))
          if(ITCLOD>0) then
            fld_info(cfld)%ntrange=1
          else
            fld_info(cfld)%ntrange=0
          endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
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
!     HIGH CLOUD FRACTION.
      IF (IGET(039).GT.0) THEN
!       GRID1=SPVAL
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            IF(CFRACH(I,J) < SPVAL) then
              GRID1(I,J) = CFRACH(I,J)*100.
            else
              GRID1(I,J) = spval
            endif
          ENDDO
        ENDDO
       if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(039))
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
!     TIME AVERAGED HIGH CLOUD FRACTION.
      IF (IGET(302) > 0) THEN
!       GRID1=SPVAL
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            IF(AVGCFRACH(I,J) < SPVAL) then
              GRID1(I,J) = AVGCFRACH(I,J)*100.
            else
              GRID1(I,J) = spval
            endif
          ENDDO
        ENDDO
        ID(1:25)=0
        ITCLOD     = NINT(TCLOD)
        IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
          IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
        ELSE
          IFINCR     = 0
        endif

        ID(19)  = IFHR
        IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN  !USE MIN FOR OFF-HR FORECAST
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
           IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(302))
          if(ITCLOD>0) then
            fld_info(cfld)%ntrange=1
          else
            fld_info(cfld)%ntrange=0
          endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
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
!     TOTAL CLOUD FRACTION (INSTANTANEOUS).
      IF ((IGET(161) > 0) .OR. (IGET(260) > 0)) THEN
!        GRID1=SPVAL
         IF(MODELNAME=='NCAR' .OR. MODELNAME=='RAPR')THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(i,j)  = SPVAL
               egrid1(i,j)=0.
              do l = 1,LM
               egrid1(i,j)=max(egrid1(i,j),cfr(i,j,l))
              end do
             ENDDO
           ENDDO

         ELSE IF (MODELNAME=='NMM'.OR.MODELNAME=='FV3R' &
           .OR. MODELNAME=='GFS')THEN
           DO J=JSTA,JEND
             DO I=1,IM
!               EGRID1(I,J)=AMAX1(CFRACL(I,J),
!     1                 AMAX1(CFRACM(I,J),CFRACH(I,J)))
!            EGRID1(I,J)=1.-(1.-CFRACL(I,J))*(1.-CFRACM(I,J))*      &  
!     &                 (1.-CFRACH(I,J))
            GRID1(i,j)=SPVAL
            EGRID1(I,J)=TCLD(I,J)
          ENDDO
          ENDDO
         END IF
!$omp parallel do private(i,j)
         DO J=JSTA,JEND
           DO I=1,IM
             IF(ABS(EGRID1(I,J)-SPVAL) > SMALL) THEN
               GRID1(I,J) = EGRID1(I,J)*100.
               TCLD(I,J)  = EGRID1(I,J)*100.         !B ZHOU, PASSED to CALCEILING
             END IF 
           ENDDO
         ENDDO
         IF (IGET(161).GT.0) THEN
          if(grib=="grib2" )then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(161))
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
!     TIME AVERAGED TOTAL CLOUD FRACTION.
      IF (IGET(144) > 0) THEN
!        GRID1=SPVAL
        IF(MODELNAME == 'GFS')THEN
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              IF(ABS(AVGTCDC(I,J)-SPVAL) > SMALL) then
                GRID1(I,J) = AVGTCDC(I,J)*100.
              else
                GRID1(I,J) = spval
              endif
            END DO
          END DO 

        ELSE IF(MODELNAME == 'NMM')THEN
          DO J=JSTA,JEND
            DO I=1,IM
!             RSUM = NCFRST(I,J)+NCFRCV(I,J)
!             IF (RSUM.GT.0.0) THEN
!                EGRID1(I,J)=(ACFRST(I,J)+ACFRCV(I,J))/RSUM
!             ELSE
!                EGRID1(I,J) = D00
!             ENDIF
!ADDED BRAD'S MODIFICATION
              RSUM = D00
              IF (NCFRST(I,J) .GT. 0) RSUM=ACFRST(I,J)/NCFRST(I,J)
              IF (NCFRCV(I,J) .GT. 0)                               &
                RSUM=MAX(RSUM, ACFRCV(I,J)/NCFRCV(I,J))
              GRID1(I,J) = RSUM*100.
            ENDDO
          ENDDO
        END IF 
        IF(MODELNAME == 'NMM' .OR. MODELNAME == 'GFS')THEN
          ID(1:25)= 0
          ITCLOD     = NINT(TCLOD)
          IF(ITCLOD .ne. 0) then
            IFINCR     = MOD(IFHR,ITCLOD)
            IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
          ELSE
            IFINCR     = 0
          endif

          ID(19)  = IFHR
          IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN  !USE MIN FOR OFF-HR FORECAST
          ID(20)  = 3
          IF (IFINCR.EQ.0) THEN
              ID(18)  = IFHR-ITCLOD
          ELSE
              ID(18)  = IFHR-IFINCR
              IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
          ENDIF
          IF (ID(18).LT.0) ID(18) = 0
        ENDIF
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(144))
          if(ITCLOD>0) then
             fld_info(cfld)%ntrange=1
          else
             fld_info(cfld)%ntrange=0
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)
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
!     TIME AVERAGED STRATIFORM CLOUD FRACTION.
         IF (IGET(139).GT.0) THEN
           IF(MODELNAME /= 'NMM')THEN
             GRID1=SPVAL
           ELSE 
            DO J=JSTA,JEND
            DO I=1,IM
               IF (NCFRST(I,J).GT.0.0) THEN
                  GRID1(I,J) = ACFRST(I,J)/NCFRST(I,J)*100.
               ELSE
                  GRID1(I,J) = D00
               ENDIF
            ENDDO
            ENDDO
           END IF 
          IF(MODELNAME.EQ.'NMM')THEN
           ID(1:25)=0
           ITCLOD     = NINT(TCLOD)
	   IF(ITCLOD .ne. 0) then
            IFINCR     = MOD(IFHR,ITCLOD)
	    IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	   ELSE
	    IFINCR     = 0
           endif 
           ID(19)  = IFHR
	   IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
           ID(20)  = 3
           IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITCLOD
           ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
           ENDIF
           IF (ID(18).LT.0) ID(18) = 0
          ENDIF
          if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(139))
            if(ITCLOD>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
          endif
         ENDIF
!    
!     TIME AVERAGED CONVECTIVE CLOUD FRACTION.
         IF (IGET(143).GT.0) THEN
           IF(MODELNAME /= 'NMM')THEN
	    GRID1=SPVAL
	   ELSE  
            DO J=JSTA,JEND
            DO I=1,IM
               IF (NCFRCV(I,J).GT.0.0) THEN
                  GRID1(I,J) = ACFRCV(I,J)/NCFRCV(I,J)*100.
               ELSE
                  GRID1(I,J) = D00
               ENDIF
            ENDDO
            ENDDO
	   END IF
           IF(MODELNAME.EQ.'NMM')THEN 
            ID(1:25)=0
            ITCLOD     = NINT(TCLOD)
	    IF(ITCLOD .ne. 0) then
             IFINCR     = MOD(IFHR,ITCLOD)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	    ELSE
	     IFINCR     = 0
            endif 
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITCLOD
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
          ENDIF
          if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(143))
            if(ITCLOD>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
          endif
         ENDIF
!    
!     CLOUD BASE AND TOP FIELDS 
      IF((IGET(148).GT.0) .OR. (IGET(149).GT.0) .OR.              &
          (IGET(168).GT.0) .OR. (IGET(178).GT.0) .OR.             &
          (IGET(179).GT.0) .OR. (IGET(194).GT.0) .OR.             &
          (IGET(408).GT.0) .OR. (IGET(798).GT.0) .OR.             & 
          (IGET(409).GT.0) .OR. (IGET(406).GT.0) .OR.             &
          (IGET(195).GT.0) .OR. (IGET(260).GT.0) .OR.             &
          (IGET(275).GT.0))  THEN
!
!--- Calculate grid-scale cloud base & top arrays (Ferrier, Feb '02)
!
!--- Rain is not part of cloud, only cloud water + cloud ice + snow
!
        DO J=JSTA,JEND
          DO I=1,IM
!
!--- Various convective cloud base & cloud top levels
!
!     write(0,*)' hbot=',hbot(i,j),' hbotd=',hbotd(i,j),'
!     hbots=',hbots(i,j)&
!  ,' htop=',htop(i,j),' htopd=',htopd(i,j),' htops=',htops(i,j),i,j
! Initilize
            IBOTCu(I,J) = 0
            ITOPCu(I,J) = 100
            IBOTDCu(I,J) = 0
            ITOPDCu(I,J) = 100
            IBOTSCu(I,J) = 0
            ITOPSCu(I,J) = 100
            if (hbot(i,j) .ne. spval) then
              IBOTCu(I,J) = NINT(HBOT(I,J))
            endif
            if (hbotd(i,j) .ne. spval) then
              IBOTDCu(I,J) = NINT(HBOTD(I,J))
            endif
            if (hbots(i,j) .ne. spval) then
              IBOTSCu(I,J) = NINT(HBOTS(I,J))
            endif
            if (htop(i,j) .ne. spval) then
              ITOPCu(I,J) = NINT(HTOP(I,J))
            endif
            if (htopd(i,j) .ne. spval) then
              ITOPDCu(I,J) = NINT(HTOPD(I,J))
            endif
            if (htops(i,j) .ne. spval) then
              ITOPSCu(I,J) = NINT(HTOPS(I,J))
            endif
            IF (IBOTCu(I,J)-ITOPCu(I,J) .LE. 1) THEN
              IBOTCu(I,J) = 0
              ITOPCu(I,J) = 100
            ENDIF
            IF (IBOTDCu(I,J)-ITOPDCu(I,J) .LE. 1) THEN
              IBOTDCu(I,J) = 0
              ITOPDCu(I,J) = 100
            ENDIF
            IF (IBOTSCu(I,J)-ITOPSCu(I,J) .LE. 1) THEN
              IBOTSCu(I,J) = 0
              ITOPSCu(I,J) = 100
            ENDIF
! Convective cloud top height
           ITOP = ITOPCu(I,J)
           IF (ITOP > 0 .AND. ITOP < 100) THEN
!             print *, 'aha ', ITOP
           ENDIF
           IF (ITOP > 0 .AND. ITOP <= NINT(LMH(I,J))) THEN
             CLDZCu(I,J) = ZMID(I,J,ITOP)
           else
             CLDZCu(I,J) = -5000.
           endif

!   !
    !--- Grid-scale cloud base & cloud top levels 
    !
    !--- Grid-scale cloud occurs when the mixing ratio exceeds QCLDmin
    !    or in the presence of snow when RH>=95% or at/above the PBL top.
    !
        if(MODELNAME == 'RAPR') then
            IBOTGr(I,J)=0
            DO L=NINT(LMH(I,J)),1,-1
              QCLD=QQW(I,J,L)+QQI(I,J,L)+QQS(I,J,L)
              IF (QCLD .GE. QCLDmin) THEN
                IBOTGr(I,J)=L
                EXIT
              ENDIF
            ENDDO    !--- End L loop
            ITOPGr(I,J)=100
            DO L=1,NINT(LMH(I,J))
              QCLD=QQW(I,J,L)+QQI(I,J,L)+QQS(I,J,L)
              IF (QCLD .GE. QCLDmin) THEN
                ITOPGr(I,J)=L
                EXIT
              ENDIF
            ENDDO    !--- End L loop
        else
            IBOTGr(I,J) = 0
            ZPBLtop     = PBLH(I,J)+ZINT(I,J,NINT(LMH(I,J))+1)
            DO L=NINT(LMH(I,J)),1,-1
              QCLD = QQW(I,J,L)+QQI(I,J,L)           !- no snow +QQS(I,J,L)
              IF (QCLD >= QCLDmin) THEN
                IBOTGr(I,J) = L
                EXIT
              ENDIF
snow_check:   IF (QQS(I,J,L)>=QCLDmin) THEN
                TMP=T(I,J,L)
                IF (TMP>=C2K) THEN
                  QSAT=PQ0/PMID(I,J,L)*EXP(A2*(TMP-A3)/(TMP-A4))
                ELSE
!-- Use Teten's formula for ice from Murray (1967).  More info at
!   http://faculty.eas.ualberta.ca/jdwilson/EAS372_13/Vomel_CIRES_satvpformulae.html
                  QSAT=PQ0/PMID(I,J,L)*EXP(21.8745584*(TMP-A3)/(TMP-7.66))
                ENDIF
                RHUM=Q(I,J,L)/QSAT
                IF (RHUM>=0.98 .AND. ZMID(I,J,L)>=ZPBLtop) THEN
                IBOTGr(I,J)=L
                EXIT
              ENDIF
              ENDIF  snow_check
            ENDDO    !--- End L loop
            ITOPGr(I,J) = 100
            DO L=1,NINT(LMH(I,J))
              QCLD=QQW(I,J,L)+QQI(I,J,L)+QQS(I,J,L)
              IF (QCLD >= QCLDmin) THEN
                ITOPGr(I,J)=L
                EXIT
              ENDIF
            ENDDO    !--- End L loop
        endif
    !
    !--- Combined (convective & grid-scale) cloud base & cloud top levels 
            IF(MODELNAME .EQ. 'NCAR' .OR. MODELNAME == 'RAPR')THEN
              IBOTT(I,J) = IBOTGr(I,J)
              ITOPT(I,J) = ITOPGr(I,J)
	    ELSE
              IBOTT(I,J) = MAX(IBOTGr(I,J), IBOTCu(I,J))
!	      if(i==200 .and. j==139)print*,'Debug cloud base 1: ',&
!             IBOTGr(I,J),IBOTCu(I,J),ibott(i,j)
              ITOPT(I,J) = MIN(ITOPGr(I,J), ITOPCu(I,J))
	    END IF 
          ENDDO      !--- End I loop
        ENDDO        !--- End J loop
      ENDIF          !--- End IF tests 
!
! CONVECTIVE CLOUD TOP HEIGHT
      IF (IGET(758).GT.0) THEN

          DO J=JSTA,JEND
          DO I=1,IM
              GRID1(I,J) = CLDZCu(I,J)
          ENDDO
          ENDDO
          ID(1:25)=0
          if(grib=="grib2" )then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(758))
               datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
          endif
      ENDIF
!
!-------------------------------------------------
!-----------  VARIOUS CLOUD BASE FIELDS ----------
!-------------------------------------------------
!
!--- "TOTAL" CLOUD BASE FIELDS (convective + grid-scale;  Ferrier, Feb '02)
!
      IF ((IGET(148).GT.0) .OR. (IGET(178).GT.0) .OR.(IGET(260).GT.0) ) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            IBOT=IBOTT(I,J)     !-- Cloud base ("bottoms")
            IF(MODELNAME == 'RAPR') then
               IF (IBOT .LE. 0) THEN
                 CLDP(I,J) = SPVAL
                 CLDZ(I,J) = SPVAL
               ELSE IF (IBOT .LE. NINT(LMH(I,J))) THEN
                 CLDP(I,J) = PMID(I,J,IBOT)
                 IF (IBOT .EQ. LM) THEN
                   CLDZ(I,J) = ZINT(I,J,LM)
                 ELSE
                   CLDZ(I,J) = HTM(I,J,IBOT+1)*T(I,J,IBOT+1)   &
                           *(Q(I,J,IBOT+1)*D608+H1)*ROG*    &
                           (LOG(PINT(I,J,IBOT+1))-LOG(CLDP(I,J)))&
                           +ZINT(I,J,IBOT+1)
                 ENDIF     !--- End IF (IBOT .EQ. LM) ...
               ENDIF       !--- End IF (IBOT .LE. 0) ...
            ELSE
               IF (IBOT>0 .AND. IBOT<=NINT(LMH(I,J))) THEN
                 CLDP(I,J) = PMID(I,J,IBOT)
                 CLDZ(I,J) = ZMID(I,J,IBOT)
               ELSE
                 CLDP(I,J) = -50000.
                 CLDZ(I,J) = -5000.
               ENDIF       !--- End IF (IBOT .LE. 0) ...
            ENDIF
          ENDDO         !--- End DO I loop
        ENDDO           !--- End DO J loop
!   CLOUD BOTTOM PRESSURE
         IF (IGET(148).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = CLDP(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
             if(grib=="grib2" )then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(148))
               datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
             endif
         ENDIF 
!    CLOUD BOTTOM HEIGHT
         IF (IGET(178).GT.0) THEN
!--- Parameter was set to 148 in operational code  (Ferrier, Feb '02)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = CLDZ(I,J)
               ENDDO
               ENDDO
             if(grib=="grib2" )then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(178))
               datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
             endif
         ENDIF
      ENDIF

!    GSD CLOUD CEILING ALGORITHM
!    J. Kenyon, 3 Feb 2017:  formerly described here as
!    "GSD CLOUD BOTTOM HEIGHT".  An alternative (experimental)
!    GSD cloud ceiling algorithm is offered further below.

      IF (IGET(408).GT.0 .OR. IGET(798).GT.0) THEN
!- imported from RUC post
!  -- constants for effect of snow on ceiling
!      Also found in calvis.f
        rhoice = 970.
        coeffp = 10.36
! - new value from Roy Rasmussen - Dec 2003
!        exponfp = 0.7776 
! change consistent with CALVIS_GSD.f
        exponfp = 1.
        const1 = 3.912

        nfog = 0
        do k=1,7
         nfogn(k) = 0
        end do
        npblcld = 0

        Cloud_def_p = 0.0000001

        DO J=JSTA,JEND
  
          DO I=1,IM
!
!- imported from RUC post
            IF(MODELNAME == 'RAPR') then
              CLDZ(I,J) = SPVAL 
              pcldbase = SPVAL
              zcldbase = SPVAL 
            ELSE
              CLDZ(I,J) = -5000.
              pcldbase = -50000.
              zcldbase = -5000.
            ENDIF
          watericemax = -99999.
          do k=1,lm
            LL=LM-k+1
            watericetotal(k) = QQW(i,j,ll) + QQI(i,j,ll)
            watericemax = max(watericemax,watericetotal(k))
          end do
          loop3701:do
!         if (watericemax.lt.cloud_def_p) go to 3701
          if (watericemax.lt.cloud_def_p) exit loop3701

!  Cloud base
!====================

! --- Check out no. of points with thin cloud layers near surface
         do k=2,3
           pabovesfc(k) = pint(i,j,lm) - pint(i,j,lm-k+1)
           if (watericetotal(k).lt.cloud_def_p) then
! --- wimin is watericemin in lowest few levels
             wimin = 100.
             do k1=k-1,1,-1
              wimin = min(wimin,watericetotal(k1))
             end do
             if (wimin.gt.cloud_def_p) then
               nfogn(k)= nfogn(k)+1
             end if
           end if
         end do

         loop3789:do

!        Eliminate fog layers near surface in watericetotal array
         do 1778 k=2,3
! --- Do this only when at least 10 mb (1000 Pa) above surface
!          if (pabovesfc(k).gt.1000.) then
           if (watericetotal(k).lt.cloud_def_p) then
             loop3441:do
             if (watericetotal(1).gt.cloud_def_p) then
               nfog = nfog+1
!               go to 3441
                exit loop3441
             end if
!            go to 3789
             exit loop3789
             exit loop3441
             enddo loop3441
3441         continue
             do k1=1,k-1
               if (watericetotal(k1).ge.cloud_def_p) then
!                print *,'Zero fog',i,j,k1,watericetotal(k1),
!    1               g3(i,j,k1,p_p)/100.
                 watericetotal(k1)=0.
               end if
             end do
           end if
!          go to 3789
           exit loop3789
!          end if
1778     continue

         exit loop3789
         enddo loop3789

3789     continue

!!       At surface?
!commented out 16aug11
!          if (watericetotal(1).gt.cloud_def_p) then
!            zcldbase = zmid(i,j,lm)
!            go to 3788
!          end if
!!       Aloft?

          loop372: do
          do 371 k=2,lm
            k1 = k
!           if (watericetotal(k).gt.cloud_def_p) go to 372
            if (watericetotal(k).gt.cloud_def_p) exit loop372
 371      continue
!         go to 3701
          exit loop3701
          exit loop372
          enddo loop372
 372      continue
        if (k1.le.4) then
! -- If within 4 levels of surface, just use lowest cloud level
!     as ceiling WITHOUT vertical interpolation.
           zcldbase = zmid(i,j,lm-k1+1)
           pcldbase = pmid(i,j,lm-k1+1)
        else   
! -- Use vertical interpolation to obtain cloud level
        zcldbase = zmid(i,j,lm-k1+1) + (cloud_def_p-watericetotal(k1))    &
                 * (zmid(i,j,lm-k1+2)-zmid(i,j,lm-k1+1))                  &
                 / (watericetotal(k1-1) - watericetotal(k1))
        pcldbase = pmid(i,j,lm-k1+1) + (cloud_def_p-watericetotal(k1))    &
                 * (pmid(i,j,lm-k1+2)-pmid(i,j,lm-k1+1))                  &
                 / (watericetotal(k1-1) - watericetotal(k1))
        end if
        zcldbase  = max(zcldbase,FIS(I,J)*GI+5.)

 3788   continue

        loop3743:do

! -- consider lowering of ceiling due to falling snow
!      -- extracted from calvis.f (visibility diagnostic)
          if (QQS(i,j,LM).gt.0.) then
            TV=T(I,J,lm)*(H1+D608*Q(I,J,lm))
            RHOAIR=PMID(I,J,lm)/(RD*TV)
            vovermd = (1.+Q(i,j,LM))/rhoair + QQS(i,j,LM)/rhoice
            concfp = QQS(i,j,LM)/vovermd*1000.
            betav = coeffp*concfp**exponfp + 1.e-10
            vertvis = 1000.*min(90., const1/betav)
            if (vertvis .lt. zcldbase-FIS(I,J)*GI ) then
              zcldbase = FIS(I,J)*GI + vertvis

              loop3742:do
              do 3741 k=2,LM
              k1 = k
!               if (ZMID(i,j,lm-k+1) .gt. zcldbase) go to 3742
                if (ZMID(i,j,lm-k+1) .gt. zcldbase) exit loop3742
 3741         continue
!             go to 3743
              exit loop3743
              exit loop3742
              enddo loop3742
 3742         continue
           pcldbase = pmid(i,j,lm-k1+2) + (zcldbase-ZMID(i,j,lm-k1+2))   &
               *(pmid(i,j,lm-k1+1)-pmid(i,j,lm-k1+2) )                   &
               /(zmid(i,j,lm-k1+1)-zmid(i,j,lm-k1+2) )
            end if
          end if

          exit loop3743
          enddo loop3743
 3743     continue

          exit loop3701
          enddo loop3701

 3701  continue

!new 15 aug 2011
              CLDZ(I,J) = zcldbase
              CLDP(I,J) = pcldbase

! --- Now, do a PBL cloud check.
! --- First, get a PBL-top cloud ceiling, if it exists.
!     This value is the first level under the cloud top if
!       the RH is greater than 95%.   This should help to identify
!       ceilings that the RUC model doesn't quite catch due to
!       vertical resolution.

! - compute relative humidity
         do k=1,LM
        LL=LM-K+1
        Tx=T(I,J,LL)-273.15
        POL = 0.99999683       + TX*(-0.90826951E-02 +                  &
           TX*(0.78736169E-04   + TX*(-0.61117958E-06 +                 &
           TX*(0.43884187E-08   + TX*(-0.29883885E-10 +                 &
           TX*(0.21874425E-12   + TX*(-0.17892321E-14 +                 &
           TX*(0.11112018E-16   + TX*(-0.30994571E-19)))))))))
        esx = 6.1078/POL**8

          ES = esx
          E = PMID(I,J,LL)/100.*Q(I,J,LL)/(0.62197+Q(I,J,LL)*0.37803)
          RHB(k) = 100.*MIN(1.,E/ES)
!
!     COMPUTE VIRTUAL POTENTIAL TEMPERATURE.
!
         enddo

! PBL height is computed in INITPOST.f
! zpbltop is relative to sea level
            ZSF=ZINT(I,J,NINT(LMH(I,J))+1)
            zpbltop = PBLH(I,J)+ZSF

!            PBLH(I,J)= zpbltop - FIS(I,J)*GI
!         print *,'I,J,k1,zmid(i,j,lm-k1+1),zmid(i,j,lm-k1),PBLH(I,J)',
!     1   I,J,k1,zmid(i,j,lm-k1+1),zmid(i,j,lm-k1),PBLH(I,J),RHB(k1)

         loop745:do

         loop744:do
         do k2=3,20
!          if (zpbltop.lt.ZMID(i,j,LM-k2+1)) go to 744
           if (zpbltop.lt.ZMID(i,j,LM-k2+1)) exit loop744
         end do
!        go to 745     ! No extra considerations for PBL-top cloud
         exit loop745  ! No extra considerations for PBL-top cloud

         exit loop744
         enddo loop744
  744    continue

!       print*,'check RH at PBL top, RH,i,j,k2',RHB(k2-1),i,j,k2-1
         if (rhb(k2-1).gt.95. ) then
           zcldbase = ZMID(i,j,LM-k2+2)
!       print*,' PBL cloud ceiling',zcldbase,i,j
           if (CLDZ(i,j).lt.-100.) then
!       print*,'add PBL cloud ceiling',zcldbase,i,j,k2
!     1         ,RHB(k2-1)
             npblcld = npblcld+1
             CLDZ(i,j) = zcldbase
             CLDP(I,J) = PMID(i,j,LM-k2+2)
!            go to 745
             exit loop745
           end if
           if ( zcldbase.lt.CLDZ(I,J)) then
!       print*,' change to PBL cloud ceiling',zcldbase,CLDZ(I,J),i,j,k2
!     1         ,RHB(k2-1)
!cc             npblcld = npblcld+1
             CLDZ(I,J) = zcldbase
           end if
         end if
         exit loop745
         enddo loop745
  745    continue

!- include convective clouds
           IBOT=IBOTCu(I,J)
       loop746:do
       if(IBOT.gt.0) then
!        print *,'IBOTCu(i,j)',i,j,IBOTCu(i,j)
         if(CLDZ(I,J).lt.-100.) then
!        print *,'add convective cloud, IBOT,CLDZ(I,J),ZMID(I,J,IBOT)'
!     1        ,IBOT,CLDZ(I,J),ZMID(I,J,IBOT),i,j
            CLDZ(I,J)=ZMID(I,J,IBOT)
!           GOTO 746
            exit loop746
         else if(ZMID(I,J,IBOT).lt.CLDZ(I,J)) then
!        print *,'change ceiling for convective cloud, CLDZ(I,J),
!     1              ZMID(I,J,IBOT),IBOT,i,j'
!     1        ,IBOT,CLDZ(I,J),ZMID(I,J,IBOT),IBOT,i,j
            CLDZ(I,J)=ZMID(I,J,IBOT)
         endif
       endif

       exit loop746
       enddo loop746
 746     continue

          ENDDO      !--- End I loop
        ENDDO        !--- End J loop

      write(6,*)'No. pts with PBL-cloud  =',npblcld
      write(6,*)'No. pts to eliminate fog =',nfog
      do k=2,7
       write(6,*)'No. pts with fog below lev',k,' =',nfogn(k)
      end do

      nlifr = 0
      DO J=JSTA,JEND
      DO I=1,IM
        zcld = CLDZ(i,j) - FIS(I,J)*GI
        if (CLDZ(i,j).ge.0..and.zcld.lt.160.) nlifr = nlifr+1
      end do
      end do
      write(6,*)'No. pts w/ LIFR ceiling =',nlifr

! GSD CLOUD BOTTOM HEIGHTS
          IF (IGET(408).GT.0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J) = CLDZ(I,J)
              ENDDO
            ENDDO
               if(grib=="grib2" )then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(408))
                 datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
               endif
          ENDIF
!   GSD CLOUD BOTTOM PRESSURE
          IF (IGET(798).GT.0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J) = CLDP(I,J)
              ENDDO
            ENDDO
               if(grib=="grib2" )then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(798))
                 datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
               endif
          ENDIF
      ENDIF   !End of GSD algorithm

! BEGIN EXPERIMENTAL GSD CEILING DIAGNOSTIC...
! J. Kenyon, 4 Feb 2017:  this approach uses model-state cloud fractions
      IF (IGET(487).GT.0) THEN
!       set some constants for ceiling adjustment in snow (retained from legacy algorithm, also in calvis.f)
        rhoice = 970.
        coeffp = 10.36
        exponfp = 1.
        const1 = 3.912
!       set minimum cloud fraction to represent a ceiling
        ceiling_thresh_cldfra = 0.5

        DO J=JSTA,JEND
          DO I=1,IM
              ceil(I,J) = SPVAL
              zceil     = SPVAL
              cldfra_max = 0.
              do k=1,lm
                LL=LM-k+1
                cldfra(k) = cfr_raw(i,j,ll)
                cldfra_max = max(cldfra_max,cldfra(k))              ! determine the column-maximum cloud fraction
              end do
            
              loop4701: do
!             if (cldfra_max .lt. ceiling_thresh_cldfra) go to 4701 ! threshold cloud fraction not found in column, so skip to end
              if (cldfra_max .lt. ceiling_thresh_cldfra) exit loop4701 ! threshold cloud fraction not found in column, so skip to end

!             threshold cloud fraction (possible ceiling) found somewhere in column, so proceed...
!             first, search for and eliminate fog layers near surface (retained from legacy diagnostic)

              loop4789: do
              do 2778 k=2,3
                if (cldfra(k) .lt. ceiling_thresh_cldfra) then   ! these two lines:
                  loop4441:do
                  if (cldfra(1) .gt. ceiling_thresh_cldfra) then ! ...look for surface-based fog beneath less-cloudy layers 
!                   go to 4441   ! found surface-based fog beneath level k
                    exit loop4441   ! found surface-based fog beneath level k
                  end if
!                 go to 4789     ! level k=2,3 has no ceiling, and no fog at surface, so skip ahead
                  exit loop4789     ! level k=2,3 has no ceiling, and no fog at surface, so skip ahead
                  exit loop4441
                  enddo loop4441
4441              continue       
                  do k1=1,k-1    ! now perform the clearing for k=1 up to k-1
                    if (cldfra(k1) .ge. ceiling_thresh_cldfra) then
                      cldfra(k1)=0.
                    end if
                  end do

                end if
!               go to 4789 
                exit loop4789 
2778          continue

              exit loop4789
              enddo loop4789
4789          continue

!             now search aloft...
              loop472:do
              do 471 k=2,lm
                k1 = k
!               if (cldfra(k) .ge. ceiling_thresh_cldfra) go to 472 ! found ceiling
                if (cldfra(k) .ge. ceiling_thresh_cldfra) exit loop472 ! found ceiling
471           continue
!             go to 4701                                            ! no ceiling found
              exit loop4701                                            ! no ceiling found
              exit loop472
              enddo loop472
472           continue
              if (k1 .le. 4) then ! within 4 levels of surface, no interpolation
                 zceil = zmid(i,j,lm-k1+1)
              else                ! use linear interpolation
                 zceil = zmid(i,j,lm-k1+1) + (ceiling_thresh_cldfra-cldfra(k1)) &
                         * (zmid(i,j,lm-k1+2)-zmid(i,j,lm-k1+1))                &
                         / (cldfra(k1-1) - cldfra(k1))
              end if
              zceil = max(zceil,FIS(I,J)*GI+5.)

!         consider lowering of ceiling due to falling snow (retained from legacy diagnostic)
!         ...this is extracted from calvis.f (visibility diagnostic)
              loop4743: do
              if (QQS(i,j,LM).gt.0.) then
                TV=T(I,J,lm)*(H1+D608*Q(I,J,lm))
                RHOAIR=PMID(I,J,lm)/(RD*TV)
                vovermd = (1.+Q(i,j,LM))/rhoair + QQS(i,j,LM)/rhoice
                concfp = QQS(i,j,LM)/vovermd*1000.
                betav = coeffp*concfp**exponfp + 1.e-10
                vertvis = 1000.*min(90., const1/betav)
                if (vertvis .lt. zceil-FIS(I,J)*GI ) then
                  zceil = FIS(I,J)*GI + vertvis
                  loop4742: do
                  do 4741 k=2,LM
                  k1 = k
!                   if (ZMID(i,j,lm-k+1) .gt. zceil) go to 4742
                    if (ZMID(i,j,lm-k+1) .gt. zceil) exit loop4742
4741              continue
!                 go to 4743
                  exit loop4743
                  exit loop4742
                  enddo loop4742
4742              continue
                end if
              end if
              exit loop4743
              enddo loop4743
4743          continue

              exit loop4701
              enddo loop4701
4701          continue
              ceil(I,J) = zceil
          ENDDO      ! i loop
        ENDDO        ! j loop

! proceed to gridding
        DO J=JSTA,JEND
        DO I=1,IM
          GRID1(I,J) = ceil(I,J)
        ENDDO
        ENDDO
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(487))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF ! end of parameter-487 conditional code
! END OF EXPERIMENTAL GSD CEILING DIAGNOSTIC
 
!    B. ZHOU: CEILING
        IF (IGET(260).GT.0) THEN
            CALL CALCEILING(CLDZ,TCLD,CEILING)
            DO J=JSTA,JEND
              DO I=1,IM
               GRID1(I,J) = CEILING(I,J)
              ENDDO
            ENDDO
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(260))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF
!    B. ZHOU: FLIGHT CONDITION RESTRICTION
        IF (IGET(261) > 0) THEN
          CALL CALFLTCND(CEILING,GRID1(1,jsta))
!         DO J=JSTA,JEND
!          DO I=1,IM
!            GRID1(I,J) = FLTCND(I,J)
!          ENDDO
!         ENDDO
          if(grib=="grib2" )then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(261))
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
!---  Convective cloud base pressures (deep & shallow; Ferrier, Feb '02)
!
      IF (IGET(188) > 0) THEN
        IF(MODELNAME == 'GFS')THEN
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              GRID1(I,J) = PBOT(I,J)
            ENDDO
          ENDDO
        ELSE
          DO J=JSTA,JEND
            DO I=1,IM
              IBOT=IBOTCu(I,J)
              IF (IBOT.GT.0 .AND. IBOT.LE.NINT(LMH(I,J))) THEN
                GRID1(I,J) = PMID(I,J,IBOT)
              ELSE
                GRID1(I,J) = -50000.
              ENDIF
            ENDDO
          ENDDO
        END IF
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(188))
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
!---  Deep convective cloud base pressures  (Ferrier, Feb '02)
!
      IF (IGET(192) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            IBOT=IBOTDCu(I,J)
            IF (IBOT.GT.0 .AND. IBOT.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,IBOT)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
      if(grib=="grib2" )then
        cfld=cfld+1
        fld_info(cfld)%ifld=IAVBLFLD(IGET(192))
        datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
       endif
       ENDIF 
!---  Shallow convective cloud base pressures   (Ferrier, Feb '02)
!
      IF (IGET(190) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            IBOT=IBOTSCu(I,J)  
            IF (IBOT.GT.0 .AND. IBOT.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,IBOT)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
      if(grib=="grib2" )then
        cfld=cfld+1
        fld_info(cfld)%ifld=IAVBLFLD(IGET(190))
        datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
       endif
       ENDIF
!---  Base of grid-scale cloudiness   (Ferrier, Feb '02)
!
      IF (IGET(194) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            IBOT=IBOTGr(I,J)
            IF (IBOT.GT.0 .AND. IBOT.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,IBOT)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
      if(grib=="grib2" )then
        cfld=cfld+1
        fld_info(cfld)%ifld=IAVBLFLD(IGET(194))
        datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
       endif
       ENDIF
       
  !---  Base of low cloud 
  !
      IF (IGET(303) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
!             IF(PBOTL(I,J) > SMALL)THEN
	      GRID1(I,J) = PBOTL(I,J)
!	     ELSE
!	      GRID1(I,J) = SPVAL
!	     END IF  
          ENDDO
        ENDDO
        ID(1:25)=0
	ITCLOD     = NINT(TCLOD)
	IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
	  IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	ELSE
	  IFINCR     = 0
        ENDIF 
        ID(19)  = IFHR
	IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
      if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(303))
            if(ITCLOD==0) then
              fld_info(cfld)%ntrange=0
            else
              fld_info(cfld)%ntrange=1
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)

            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
       endif
       ENDIF
  !---  Base of middle cloud  
  !
      IF (IGET(306) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
	     IF(PBOTM(I,J) > SMALL)THEN
	      GRID1(I,J) = PBOTM(I,J)
	     ELSE
	      GRID1(I,J) = SPVAL
	     END IF
          ENDDO
        ENDDO
        ID(1:25)=0
	ITCLOD     = NINT(TCLOD)
	IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
	  IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	ELSE
	  IFINCR     = 0
        ENDIF 
        ID(19)  = IFHR
	IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
      if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(306))
          if(ITCLOD==0) then
              fld_info(cfld)%ntrange=0
          else
              fld_info(cfld)%ntrange=1
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)

            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
      endif
       ENDIF
  !---  Base of high cloud   
  !
      IF (IGET(309) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
	     IF(PBOTH(I,J) > SMALL)THEN
	      GRID1(I,J) = PBOTH(I,J)
	     ELSE
	      GRID1(I,J) = SPVAL
	     END IF
          ENDDO
        ENDDO
        ID(1:25)=0
	ITCLOD     = NINT(TCLOD)
	IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
	  IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	ELSE
	  IFINCR     = 0
        ENDIF 
        ID(19)  = IFHR
	IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
       if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(309))
          if(ITCLOD==0) then
              fld_info(cfld)%ntrange=0
          else
              fld_info(cfld)%ntrange=1
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)
         
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
       endif
       ENDIF
!
!------------------------------------------------
!-----------  VARIOUS CLOUD TOP FIELDS ----------
!------------------------------------------------
!
!--- "TOTAL" CLOUD TOP FIELDS (convective + grid-scale;  Ferrier, Feb '02)
!
      IF ((IGET(149).GT.0) .OR. (IGET(179).GT.0) .OR.                    &
          (IGET(168).GT.0) .OR. (IGET(275).GT.0)) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            ITOP=ITOPT(I,J)
            IF (ITOP>0 .AND. ITOP<=NINT(LMH(I,J))) THEN
              CLDP(I,J) = PMID(I,J,ITOP)
              CLDZ(I,J) = ZMID(I,J,ITOP)
              CLDT(I,J) = T(I,J,ITOP)
            ELSE
              IF(MODELNAME == 'RAPR') then
                CLDP(I,J) = SPVAL
                CLDZ(I,J) = SPVAL
              ELSE
                CLDP(I,J) = -50000.
                CLDZ(I,J) = -5000.
              ENDIF
              CLDT(I,J) = -500.
            ENDIF      !--- End IF (ITOP.GT.0 .AND. ITOP.LE.LMH(I,J)) ...
          ENDDO        !--- End DO I loop
        ENDDO          !--- End DO J loop
!
!   CLOUD TOP PRESSURE
!
         IF (IGET(149).GT.0) THEN
              DO J=JSTA,JEND
              DO I=1,IM
                 GRID1(I,J) = CLDP(I,J)
               ENDDO
               ENDDO
              if(grib=="grib2" )then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(149))
                datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
              endif
         ENDIF
!   CLOUD TOP HEIGHT
!
          IF (IGET(179).GT.0) THEN
              DO J=JSTA,JEND
              DO I=1,IM
                 GRID1(I,J) = CLDZ(I,J)
               ENDDO
               ENDDO
              if(grib=="grib2" )then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(179))
                datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
              endif
         ENDIF
      ENDIF

! GSD COULD TOP HEIGHTS AND PRESSURE
      IF ((IGET(409).GT.0) .OR. (IGET(406).GT.0)) THEN

        Cloud_def_p = 0.0000001

        DO J=JSTA,JEND
          DO I=1,IM
! imported from RUC post
!  Cloud top
          zcldtop = -5000. 
          IF(MODELNAME == 'RAPR') zcldtop = SPVAL
          do k=1,lm
            LL=LM-k+1
            watericetotal(k) = QQW(i,j,ll) + QQI(i,j,ll)
          enddo

          loop3799: do
          if (watericetotal(LM).gt.cloud_def_p) then
            zcldtop = zmid(i,j,1)
!           go to 3799
            exit loop3799
          end if
! in RUC          do 373 k=LM,2,-1
          loop374: do
          do 373 k=LM-1,2,-1
!           if (watericetotal(k).gt.cloud_def_p) go to 374
            if (watericetotal(k).gt.cloud_def_p) exit loop374
 373      continue
!         go to 3799
          exit loop3799
          exit loop374   
          enddo loop374
 374    zcldtop = zmid(i,j,lm-k+1) + (cloud_def_p-watericetotal(k))   &
                 * (zmid(i,j,lm-k)-zmid(i,j,lm-k+1))                &
                 / (watericetotal(k+1) - watericetotal(k))
          exit loop3799
          enddo loop3799
 3799     continue

            ITOP=ITOPT(I,J)
            IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
              CLDP(I,J) = PMID(I,J,ITOP)
              CLDT(I,J) = T(I,J,ITOP)
            ELSE
              CLDP(I,J) = -50000.
              IF(MODELNAME == 'RAPR') CLDP(I,J) = SPVAL
!              CLDZ(I,J) = -5000.
              CLDT(I,J) = -500.
            ENDIF      !--- End IF (ITOP.GT.0 .AND. ITOP.LE.LMH(I,J)) ...

!- include convective clouds
           ITOP=ITOPCu(I,J)
       if(ITOP.lt.lm+1) then
!        print *,'ITOPCu(i,j)',i,j,ITOPCu(i,j)
         if(zcldtop .lt.-100.) then
!        print *,'add convective cloud, ITOP,CLDZ(I,J),ZMID(I,J,ITOP)'
!     1        ,ITOP,zcldtop,ZMID(I,J,ITOP),i,j
            zcldtop=ZMID(I,J,ITOP)
         else if(ZMID(I,J,ITOP).gt.zcldtop) then
!        print *,'change cloud top for convective cloud, zcldtop,
!     1              ZMID(I,J,ITOP),ITOP,i,j'
!     1        ,zcldtop,ZMID(I,J,ITOP),ITOP,i,j
            zcldtop=ZMID(I,J,ITOP)
         endif
       endif

! check consistency of cloud base and cloud top
            if(CLDZ(I,J).gt.-100. .and. zcldtop.lt.-100.) then
              zcldtop = CLDZ(I,J) + 200.
            endif

              CLDZ(I,J) = zcldtop   !  Now CLDZ is cloud top height

          ENDDO        !--- End DO I loop
        ENDDO          !--- End DO J loop
!
!   GSD CLOUD TOP PRESSURE
!
         IF (IGET(406).GT.0) THEN
              DO J=JSTA,JEND
              DO I=1,IM
                 GRID1(I,J) = CLDP(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
              if(grib=="grib2" )then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(406))
                datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
              endif
         ENDIF
!   GSD CLOUD TOP HEIGHT
!
          IF (IGET(409).GT.0) THEN
              DO J=JSTA,JEND
              DO I=1,IM
                 GRID1(I,J) = CLDZ(I,J)
               ENDDO
               ENDDO
              if(grib=="grib2" )then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(409))
                datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
              endif
         ENDIF
       ENDIF   ! end of GSD algorithm
!
!   CLOUD TOP TEMPS
!
          IF (IGET(168).GT.0) THEN 
              DO J=JSTA,JEND
              DO I=1,IM
                 GRID1(I,J) = CLDT(I,J)
               ENDDO
               ENDDO
              if(grib=="grib2" )then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(168))
                datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
              endif
         ENDIF 
!
!huang  CLOUD TOP BRIGHTNESS TEMPERATURE
          IF (IGET(275).GT.0) THEN
             num_thick=0   ! for debug
             GRID1=spval
             DO J=JSTA,JEND
             DO I=1,IM
               opdepth=0.
               llmh=nint(lmh(i,j))
!bsf - start
!-- Add subgrid-scale convective clouds for WRF runs
               do k=1,llmh
                 CU_ir(k)=0.
               enddo
! Chuang: GFS specified non-convective grid points as missing
              if((hbot(i,j)-spval)>small .and. (htop(i,j)-spval)>small)then
               lcbot=nint(hbot(i,j))
               lctop=nint(htop(i,j))
               if (lcbot-lctop > 1) then
                 q_conv=cnvcfr(i,j)*Qconv
                 do k=lctop,lcbot
                   if (t(i,j,k) < TRAD_ice) then
                     CU_ir(k)=abscoefi*q_conv
                   else
                     CU_ir(k)=abscoef*q_conv
                   endif
                 end do   !-- do k = lctop,lcbot
               endif      !-- if (lcbot-lctop > 1) then
              end if ! end of check for meaningful hbot and htop
               do k=1,llmh
!	         if(imp_physics==99 .and. t(i,j,k)<(tfrz-15.))then
!		  qqi(i,j,k)=qqw(i,j,k) ! because GFS only uses cloud water
!		  qqw(i,j,k)=0.
!		 end if 
                 dp=pint(i,j,k+1)-pint(i,j,k)
                 opdepth=opdepth+( CU_ir(k) + abscoef*qqw(i,j,k)+            &
!bsf - end
     &                   abscoefi*( qqi(i,j,k)+qqs(i,j,k) ) )*dp
                 if (opdepth > 1.) exit
               enddo
               if (opdepth > 1.) num_thick=num_thick+1   ! for debug
               k=min(k,llmh)
	       GRID1(I,J)=T(i,j,k)
             ENDDO
             ENDDO
!      print *,'num_points, num_thick = ',(jend-jsta+1)*im,num_thick
!!              k=0
!! 20           opdepthu=opdepthd
!!              k=k+1
!!!              if(k.eq.1) then
!!!               dp=pint(i,j,itop+k)-pmid(i,j,itop)
!!!               opdepthd=opdepthu+(abscoef*(0.75*qqw(i,j,itop)+
!!!     &                  0.25*qqw(i,j,itop+1))+abscoefi*
!!!     &                  (0.75*qqi(i,j,itop)+0.25*qqi(i,j,itop+1)))
!!!     &                        *dp/g
!!!              else
!!               dp=pint(i,j,k+1)-pint(i,j,k)
!!               opdepthd=opdepthu+(abscoef*qqw(i,j,k)+
!!     &                        abscoefi*qqi(i,j,k))*dp
!!!              end if
!!	      
!!              lmhh=nint(lmh(i,j))
!!              if (opdepthd.lt.1..and. k.lt.lmhh) then
!!               goto 20
!!              elseif (opdepthd.lt.1..and. k.eq.lmhh) then
!!	       GRID1(I,J)=T(i,j,lmhh )
!!!               prsctt=pmid(i,j,lmhh)
!!              else
!!!	       GRID1(I,J)=T(i,j,k) 
!!               if(k.eq.1)then
!!	         GRID1(I,J)=T(i,j,k)
!!	       else if(k.eq.lmhh)then
!!	         GRID1(I,J)=T(i,j,k)
!!	       else 	 	 
!!                 fac=(1.-opdepthu)/(opdepthd-opdepthu)
!!	         GRID1(I,J)=(T(i,j,k)+T(i,j,k-1))/2.0+
!!     &             (T(i,j,k+1)-T(i,j,k-1))/2.0*fac 
!!               end if    	       
!!!               prsctt=pf(i,j,k-1)+fac*(pf(i,j,k)-pf(i,j,k-1))
!!!               prsctt=min(prs(i,j,mkzh),max(prs(i,j,1),prsctt))
!!              endif
!!!              do 30 k=2,mkzh
!!!              if (prsctt.ge.prs(i,j,k-1).and.prsctt.le.prs(i,j,k)) then
!!!               fac=(prsctt-prs(i,j,k-1))/(prs(i,j,k)-prs(i,j,k-1))
!!!               ctt(i,j)=tmk(i,j,k-1)+
!!!     &            fac*(tmk(i,j,k)-tmk(i,j,k-1))-celkel
!!!               goto 40
!!!              endif
!!!   30       continue
!!!   40       continue 
!!             END DO
!!	     END DO 
           if(grib=="grib2" )then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(275))
             datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
           endif
         ENDIF

!
!---  Convective cloud top pressures (deep & shallow; Ferrier, Feb '02)
!
      IF (IGET(189) > 0) THEN
        IF(MODELNAME == 'GFS')THEN
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              GRID1(I,J) = PTOP(I,J)
            ENDDO
          ENDDO
        ELSE
          DO J=JSTA,JEND
            DO I=1,IM
              ITOP=ITOPCu(I,J) 
              IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
                GRID1(I,J) = PMID(I,J,ITOP)
              ELSE
                GRID1(I,J) = -50000.
              ENDIF
            ENDDO
          ENDDO
        END IF
        if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(189))
!$omp parallel do private(i,j,jj)
          do j=1,jend-jsta+1
            jj = jsta+j-1
            do i=1,im
              datapd(i,j,cfld) = GRID1(i,jj)
            enddo
          enddo
        endif
      END IF
!
!---  Deep convective cloud top pressures   (Ferrier, Feb '02)
!
      IF (IGET(193) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            ITOP=ITOPDCu(I,J)
            IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,ITOP)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
       if(grib=="grib2" )then
        cfld=cfld+1
        fld_info(cfld)%ifld=IAVBLFLD(IGET(193))
        datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
      endif
      END IF
!---  Shallow convective cloud top pressures  (Ferrier, Feb '02)
!
      IF (IGET(191) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            ITOP=ITOPSCu(I,J)
            IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,ITOP)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
       if(grib=="grib2" )then
        cfld=cfld+1
        fld_info(cfld)%ifld=IAVBLFLD(IGET(191))
        datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
      endif
      END IF
!
!---  Top of grid-scale cloudiness  (Ferrier, Feb '02)
!
      IF (IGET(195) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            ITOP=ITOPGr(I,J)
            IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,ITOP)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
       if(grib=="grib2" )then
        cfld=cfld+1
        fld_info(cfld)%ifld=IAVBLFLD(IGET(195))
        datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
      endif
      END IF
      
  !---  top of low cloud 
  !
      IF (IGET(304) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
	     IF(PTOPL(I,J) > SMALL)THEN
	      GRID1(I,J) = PTOPL(I,J)
	     ELSE
	      GRID1(I,J) = SPVAL
	     END IF
          ENDDO
        ENDDO
        ID(1:25)=0
	ITCLOD     = NINT(TCLOD)
	IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
	  IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	ELSE
	  IFINCR     = 0
        ENDIF 
        ID(19)  = IFHR
	IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
       if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(304))
          if(ITCLOD==0) then
              fld_info(cfld)%ntrange=0
          else
              fld_info(cfld)%ntrange=1
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)
     
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
       endif
       ENDIF
  !---  top of middle cloud  
  !
      IF (IGET(307) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
             GRID1(I,J) = PTOPM(I,J)
          ENDDO
        ENDDO
        ID(1:25)=0
	ITCLOD     = NINT(TCLOD)
	IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
	  IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	ELSE
	  IFINCR     = 0
        ENDIF 
        ID(19)  = IFHR
	IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
       if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(307))
          if(ITCLOD==0) then
              fld_info(cfld)%ntrange=0
          else
              fld_info(cfld)%ntrange=1
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)

            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
       endif
       ENDIF
  !---  top of high cloud   
  !
      IF (IGET(310) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
             GRID1(I,J) = PTOPH(I,J)
          ENDDO
        ENDDO
        ID(1:25)=0
	ITCLOD     = NINT(TCLOD)
	IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
	  IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	ELSE
	  IFINCR     = 0
        ENDIF 
        ID(19)  = IFHR
	IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
       if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(310))
          if(ITCLOD==0) then
              fld_info(cfld)%ntrange=0
          else
              fld_info(cfld)%ntrange=1
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)

            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
       endif
       ENDIF

  !---  T of low cloud top
  !
      IF (IGET(305) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
             GRID1(I,J) = TTOPL(I,J)
          ENDDO
        ENDDO
        ID(1:25)=0
	ITCLOD     = NINT(TCLOD)
	IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
	  IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	ELSE
	  IFINCR     = 0
        ENDIF 
        ID(19)  = IFHR
	IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
       if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(305))
          if(ITCLOD==0) then
              fld_info(cfld)%ntrange=0
          else
              fld_info(cfld)%ntrange=1
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)       

            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
       endif
       ENDIF
  !---  Base of middle cloud  
  !
      IF (IGET(308) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
             GRID1(I,J) = TTOPM(I,J)
          ENDDO
        ENDDO
        ID(1:25)=0
	ITCLOD     = NINT(TCLOD)
	IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
	  IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	ELSE
	  IFINCR     = 0
        ENDIF 
        ID(19)  = IFHR
	IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
       if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(308))
          if(ITCLOD==0) then
              fld_info(cfld)%ntrange=0
          else
              fld_info(cfld)%ntrange=1
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)

            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
       endif
       ENDIF
  !---  Base of high cloud   
  !
      IF (IGET(311) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
             GRID1(I,J) = TTOPH(I,J)
          ENDDO
        ENDDO
        ID(1:25)=0
	ITCLOD     = NINT(TCLOD)
	IF(ITCLOD .ne. 0) then
          IFINCR     = MOD(IFHR,ITCLOD)
	  IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	ELSE
	  IFINCR     = 0
        ENDIF 
        ID(19)  = IFHR
	IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
        ID(20)  = 3
        IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITCLOD
        ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
        ENDIF
        IF (ID(18).LT.0) ID(18) = 0
       if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(311))
          if(ITCLOD==0) then
              fld_info(cfld)%ntrange=0
          else
              fld_info(cfld)%ntrange=1
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
       endif
       ENDIF
!
!--- Convective cloud fractions from modified Slingo (1987)
!
      IF (IGET(196) .GT. 0.or.IGET(570)>0) THEN
          GRID1=SPVAL
          DO J=JSTA,JEND
          DO I=1,IM
            if(CNVCFR(I,J)/=SPVAL)GRID1(I,J)=100.*CNVCFR(I,J)   !-- convert to percent
          ENDDO
          ENDDO
          if(IGET(196)>0) then
            if(grib=="grib2" )then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(196))
             datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
            endif
          elseif(IGET(570)>0) then
            if(grib=="grib2" )then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(570))
             datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
            endif
          endif
      END IF
!
!--- Boundary layer cloud fractions 
!
      IF (IGET(342) .GT. 0) THEN
          GRID1=SPVAL
          DO J=JSTA,JEND
          DO I=1,IM
            if(PBLCFR(I,J)/=SPVAL)GRID1(I,J)=100.*PBLCFR(I,J)   !-- convert to percent
          ENDDO
          ENDDO
          ID(1:25)=0
	  ITCLOD     = NINT(TCLOD)
	  IF(ITCLOD .ne. 0) then
            IFINCR     = MOD(IFHR,ITCLOD)
	    IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	  ELSE
	    IFINCR     = 0
          ENDIF 
          ID(19)  = IFHR
	  IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
          ID(20)  = 3
          IF (IFINCR.EQ.0) THEN
            ID(18)  = IFHR-ITCLOD
          ELSE
            ID(18)  = IFHR-IFINCR
	    IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
          ENDIF
          IF (ID(18).LT.0) ID(18) = 0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(342))
          if(ITCLOD==0) then
              fld_info(cfld)%ntrange=0
          else
              fld_info(cfld)%ntrange=1
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)

            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      END IF
!
!--- Cloud work function 
!
      IF (IGET(313) .GT. 0) THEN
          DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J)=cldwork(I,J)  
          ENDDO
          ENDDO	  
          ID(1:25)=0
	  ITCLOD     = NINT(TCLOD)
	  IF(ITCLOD .ne. 0) then
            IFINCR     = MOD(IFHR,ITCLOD)
	    IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	  ELSE
	    IFINCR     = 0
          ENDIF 
          ID(19)  = IFHR
	  IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
          ID(20)  = 3
          IF (IFINCR.EQ.0) THEN
            ID(18)  = IFHR-ITCLOD
          ELSE
            ID(18)  = IFHR-IFINCR
	    IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
          ENDIF
          IF (ID(18).LT.0) ID(18) = 0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(313))
          if(ITCLOD==0) then
              fld_info(cfld)%ntrange=0
          else
              fld_info(cfld)%ntrange=1
          endif
          fld_info(cfld)%tinvstat=IFHR-ID(18)

            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      END IF      
!
!***  BLOCK 3.  RADIATION FIELDS.
!     
!
!     TIME AVERAGED SURFACE SHORT WAVE INCOMING RADIATION.
         IF (IGET(126).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME == 'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE  
!          print*,'ARDSW in CLDRAD=',ARDSW 
           IF(ARDSW.GT.0.)THEN
             RRNUM=1./ARDSW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
              IF(ASWIN(I,J)/=SPVAL)THEN
	       GRID1(I,J) = ASWIN(I,J)*RRNUM
	      ELSE
	       GRID1(I,J)=ASWIN(I,J)
	      END IF 
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDSW     = NINT(TRDSW)
	    IF(ITRDSW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDSW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	    ELSE
	     IFINCR     = 0
            endif 	    
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDSW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF 
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(126))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
         ENDIF
!
!     TIME AVERAGED SURFACE UV-B INCOMING RADIATION.
         IF (IGET(298).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME == 'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE  
!          print*,'ARDSW in CLDRAD=',ARDSW 
           IF(ARDSW.GT.0.)THEN
             RRNUM=1./ARDSW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
	     IF(AUVBIN(I,J)/=SPVAL)THEN
              GRID1(I,J) = AUVBIN(I,J)*RRNUM
	     ELSE
	      GRID1(I,J) = AUVBIN(I,J)
	     END IF  
           ENDDO
           ENDDO
            ID(1:25)=0
	    ID(02)=129
            ITRDSW     = NINT(TRDSW)
	    IF(ITRDSW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDSW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	    ELSE
	     IFINCR     = 0
            endif 	    
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDSW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF 
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(298))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
         ENDIF
!
!     TIME AVERAGED SURFACE UV-B CLEAR SKY INCOMING RADIATION.
         IF (IGET(297).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME == 'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE  
!          print*,'ARDSW in CLDRAD=',ARDSW 
           IF(ARDSW.GT.0.)THEN
             RRNUM=1./ARDSW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
	     IF(AUVBINC(I,J)/=SPVAL)THEN
              GRID1(I,J) = AUVBINC(I,J)*RRNUM
	     ELSE
	      GRID1(I,J) = AUVBINC(I,J)
	     END IF  
           ENDDO
           ENDDO
            ID(1:25)=0
	    ID(02)=129
            ITRDSW     = NINT(TRDSW)
	    IF(ITRDSW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDSW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDSW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF 
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(297))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
         ENDIF
!
!     TIME AVERAGED SURFACE LONG WAVE INCOMING RADIATION.
         IF (IGET(127).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME == 'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
           IF(ARDLW.GT.0.)THEN
             RRNUM=1./ARDLW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
	    IF(ALWIN(I,J)/=SPVAL)THEN 
             GRID1(I,J) = ALWIN(I,J)*RRNUM
	    ELSE
	     GRID1(I,J)=ALWIN(I,J)
	    END IF  
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDLW     = NINT(TRDLW)
	    IF(ITRDLW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDLW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDLW*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDLW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(127))
            if(ITRDLW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
         ENDIF
!
!     TIME AVERAGED SURFACE SHORT WAVE OUTGOING RADIATION.
         IF (IGET(128).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME == 'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
           IF(ARDSW.GT.0.)THEN
             RRNUM=1./ARDSW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
	    IF(ASWOUT(I,J)/=SPVAL)THEN
             GRID1(I,J) = -1.0*ASWOUT(I,J)*RRNUM
	    ELSE
	     GRID1(I,J)=ASWOUT(I,J)
	    END IF 
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDSW     = NINT(TRDSW)
	    IF(ITRDSW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDSW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDSW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(128))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
         ENDIF
!
!     TIME AVERAGED SURFACE LONG WAVE OUTGOING RADIATION.
         IF (IGET(129).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME == 'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
           IF(ARDLW.GT.0.)THEN
             RRNUM=1./ARDLW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
	    IF(ALWOUT(I,J)/=SPVAL)THEN
             GRID1(I,J) = -1.0*ALWOUT(I,J)*RRNUM
	    ELSE
	     GRID1(I,J)=ALWOUT(I,J)
	    END IF  
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDLW     = NINT(TRDLW)
	    IF(ITRDLW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDLW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDLW*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDLW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(129))
            if(ITRDLW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
         ENDIF
!
!     TIME AVERAGED TOP OF THE ATMOSPHERE SHORT WAVE RADIATION.
         IF (IGET(130).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME == 'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
           IF(ARDSW.GT.0.)THEN
             RRNUM=1./ARDSW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
	    IF(ASWTOA(I,J)/=SPVAL)THEN
             GRID1(I,J) = ASWTOA(I,J)*RRNUM
	    ELSE
	     GRID1(I,J)=ASWTOA(I,J)
	    END IF  
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDSW     = NINT(TRDSW)
	    IF(ITRDSW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDSW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDSW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(130))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
         ENDIF
!
!     TIME AVERAGED TOP OF THE ATMOSPHERE LONG WAVE RADIATION.
         IF (IGET(131).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME == 'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
           IF(ARDLW.GT.0.)THEN
             RRNUM=1./ARDLW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
	    IF(ALWTOA(I,J)/=SPVAL)THEN
             GRID1(I,J) = ALWTOA(I,J)*RRNUM
	    ELSE
	     GRID1(I,J)=ALWTOA(I,J)
	    END IF  
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDLW     = NINT(TRDLW)
            IF(ITRDLW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDLW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDLW*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDLW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(131))
            if(ITRDLW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
         ENDIF
!
!     CURRENT TOP OF THE ATMOSPHERE LONG WAVE RADIATION.
         IF (IGET(274).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	   GRID1=SPVAL
	   ID(1:25)=0
	  ELSE
           DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RLWTOA(I,J)
           ENDDO
           ENDDO
           ID(1:25)=0
	  END IF  
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(274))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
         ENDIF
!
!     CLOUD TOP BRIGHTNESS TEMPERATURE FROM TOA OUTGOING LW.
         IF (IGET(265).GT.0) THEN
	  GRID1=SPVAL
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME == 'RAPR')THEN
	   GRID1=SPVAL
	  ELSE
           DO J=JSTA,JEND
           DO I=1,IM
             IF(RLWTOA(I,J) .LT. SPVAL)                      &
     &         GRID1(I,J) = (RLWTOA(I,J)*STBOL)**0.25
           ENDDO
           ENDDO
	  END IF  
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(265))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
         ENDIF
!     
!     CURRENT INCOMING SW RADIATION AT THE SURFACE.
      IF (IGET(156).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           IF(CZMEAN(I,J).GT.1.E-6) THEN
             FACTRS=CZEN(I,J)/CZMEAN(I,J)
           ELSE
             FACTRS=0.0
           ENDIF
           GRID1(I,J)=RSWIN(I,J)*FACTRS
         ENDDO
         ENDDO
!
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(156))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!     
!     CURRENT INCOMING LW RADIATION AT THE SURFACE.
      IF (IGET(157).GT.0) THEN
! dong add missing value to DLWRF
         GRID1 = spval
         DO J=JSTA,JEND
         DO I=1,IM
          IF(MODELNAME.eq.'RSM' .OR. MODELNAME == 'RAPR') THEN      !add by Binbin: RSM has direct RLWIN output
           GRID1(I,J)=RLWIN(I,J)
          ELSE
           IF(SIGT4(I,J).GT.0.0) THEN
             LLMH=NINT(LMH(I,J))
             TLMH=T(I,J,LLMH)
             FACTRL=5.67E-8*TLMH*TLMH*TLMH*TLMH/SIGT4(I,J)
           ELSE
             FACTRL=0.0
           ENDIF
           IF(RLWIN(I,J) < spval) GRID1(I,J)=RLWIN(I,J)*FACTRL
          ENDIF
         ENDDO
         ENDDO
!
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(157))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!     
!     CURRENT OUTGOING SW RADIATION AT THE SURFACE.
      IF (IGET(141).GT.0) THEN
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
             IF(CZMEAN(I,J).GT.1.E-6) THEN
               FACTRS=CZEN(I,J)/CZMEAN(I,J)
             ELSE
               FACTRS=0.0
             ENDIF
             GRID1(I,J)=RSWOUT(I,J)*FACTRS
           ENDDO
         ENDDO
!
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(141))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

! Instantaneous clear-sky upwelling SW at the surface
      IF (IGET(743).GT.0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J) = SWUPBC(I,J)
          ENDDO
        ENDDO
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(743))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF

!     CURRENT OUTGOING LW RADIATION AT THE SURFACE.
      IF (IGET(142).GT.0) THEN
!$omp parallel do private(i,j)
         DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RADOT(I,J)
           ENDDO
         ENDDO
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(142))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

! Instantaneous clear-sky downwelling LW at the surface
      IF (IGET(744).GT.0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J) = LWDNBC(I,J)
          ENDDO
        ENDDO
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(744))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF

! Instantaneous clear-sky upwelling LW at the surface
      IF (IGET(745).GT.0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J) = LWUPBC(I,J)
          ENDDO
        ENDDO
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(745))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF

! Instantaneous MEAN_FRP
      IF (IGET(740).GT.0) THEN
        print *,"GETTING INTO MEAN_FRP PART"
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J) = MEAN_FRP(I,J)
          ENDDO
        ENDDO
        if(grib=='grib2') then
          print *,"GETTING INTO MEAN_FRP GRIB2 PART"
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(740))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF

!     CURRENT (instantaneous) INCOMING CLEARSKY SW RADIATION AT THE SURFACE.
      IF (IGET(262).GT.0) THEN
!$omp parallel do private(i,j)
         DO J=JSTA,JEND
           DO I=1,IM
             IF(CZMEAN(I,J).GT.1.E-6) THEN
               FACTRS=CZEN(I,J)/CZMEAN(I,J)
             ELSE
               FACTRS=0.0
             ENDIF
             GRID1(I,J) = RSWINC(I,J)*FACTRS
           ENDDO
         ENDDO
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(262))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

! Instantaneous clear-sky downwelling SW at surface (GSD version)
      IF (IGET(742).GT.0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J) = SWDNBC(I,J)
          ENDDO
        ENDDO
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(742))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF

! Instantaneous SWDDNI
      IF (IGET(772).GT.0)THEN
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J) = SWDDNI(I,J)
          ENDDO
        ENDDO
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(772))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF

! Instantaneous clear-sky SWDDNI
      IF (IGET(796).GT.0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J) = SWDDNIC(I,J)
          ENDDO
        ENDDO
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(796))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF

! Instantaneous SWDDIF
      IF (IGET(773).GT.0) THEN
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J) = SWDDIF(I,J)
          ENDDO
        ENDDO
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(773))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF

! Instantaneous clear-sky SWDDIF
      IF (IGET(797).GT.0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J) = SWDDIFC(I,J)
          ENDDO
        ENDDO
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(797))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF

!     TIME AVERAGED INCOMING CLEARSKY SW RADIATION AT THE SURFACE.
      IF (IGET(383).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = ASWINC(I,J)
         ENDDO
	 ENDDO
	 ID(1:25)=0
         ITRDSW     = NINT(TRDSW)
	 IF(ITRDSW .ne. 0) then
           IFINCR     = MOD(IFHR,ITRDSW)
	   IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	 ELSE
	   IFINCR     = 0
         endif
         ID(19)  = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)  = 3
         IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITRDSW
         ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(383))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!     
!     TIME AVERAGED OUTGOING CLEARSKY SW RADIATION AT THE SURFACE.
      IF (IGET(386).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = ASWOUTC(I,J)
         ENDDO
	 ENDDO
	 ID(1:25)=0
         ITRDSW     = NINT(TRDSW)
	 IF(ITRDSW .ne. 0) then
           IFINCR     = MOD(IFHR,ITRDSW)
	   IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	 ELSE
	   IFINCR     = 0
         endif
         ID(19)  = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)  = 3
         IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITRDSW
         ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(386))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

! Instantaneous all-sky outgoing SW flux at the model top
      IF (IGET(719).GT.0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J) = SWUPT(I,J)
          ENDDO
        ENDDO
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(719))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
      ENDIF

!     TIME AVERAGED OUTGOING CLEARSKY SW RADIATION AT THE MODEL TOP
      IF (IGET(387).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = ASWTOAC(I,J)
         ENDDO
	 ENDDO
	 ID(1:25)=0
         ITRDSW     = NINT(TRDSW)
	 IF(ITRDSW .ne. 0) then
           IFINCR     = MOD(IFHR,ITRDSW)
	   IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	 ELSE
	   IFINCR     = 0
         endif
         ID(19)  = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)  = 3
         IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITRDSW
         ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(387))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!     
!     TIME AVERAGED INCOMING SW RADIATION AT THE MODEL TOP
      IF (IGET(388).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = ASWINTOA(I,J)
         ENDDO
	 ENDDO
	 ID(1:25)=0
         ITRDSW     = NINT(TRDSW)
	 IF(ITRDSW .ne. 0) then
           IFINCR     = MOD(IFHR,ITRDSW)
	   IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	 ELSE
	   IFINCR     = 0
         endif
         ID(19)  = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)  = 3
         IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITRDSW
         ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(388))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!     
!     TIME AVERAGED INCOMING CLEARSKY LW RADIATION AT THE SURFACE
      IF (IGET(382).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = ALWINC(I,J)
         ENDDO
	 ENDDO
	 ID(1:25)=0
         ITRDLW     = NINT(TRDLW)
	 IF(ITRDLW .ne. 0) then
           IFINCR     = MOD(IFHR,ITRDLW)
	   IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDLW*60)
	 ELSE
	   IFINCR     = 0
         endif
         ID(19)  = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)  = 3
         IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITRDLW
         ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(382))
            if(ITRDLW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!     
!     TIME AVERAGED OUTGOING CLEARSKY LW RADIATION AT THE SURFACE
      IF (IGET(384).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = ALWOUTC(I,J)
         ENDDO
	 ENDDO
	 ID(1:25)=0
         ITRDLW     = NINT(TRDLW)
	 IF(ITRDLW .ne. 0) then
           IFINCR     = MOD(IFHR,ITRDLW)
	   IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDLW*60)
	 ELSE
	   IFINCR     = 0
         endif
         ID(19)  = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)  = 3
         IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITRDLW
         ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(384))
            if(ITRDLW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!     
!     TIME AVERAGED OUTGOING CLEARSKY LW RADIATION AT THE MODEL TOP
      IF (IGET(385).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = ALWTOAC(I,J)
         ENDDO
	 ENDDO
	 ID(1:25)=0
         ITRDLW     = NINT(TRDLW)
	 IF(ITRDLW .ne. 0) then
           IFINCR     = MOD(IFHR,ITRDLW)
	   IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDLW*60)
	 ELSE
	   IFINCR     = 0
         endif
         ID(19)  = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)  = 3
         IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITRDLW
         ELSE
           ID(18)  = IFHR-IFINCR
	   IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(385))
            if(ITRDLW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!
!     TIME AVERAGED SURFACE VISIBLE BEAM DOWNWARD SOLAR FLUX
      IF (IGET(401).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = AVISBEAMSWIN(I,J)
         ENDDO
         ENDDO
         ID(1:25)=0
         ITRDSW     = NINT(TRDSW)
         IF(ITRDSW .ne. 0) then
           IFINCR     = MOD(IFHR,ITRDSW)
           IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
         ELSE
           IFINCR     = 0
         endif
         ID(19)  = IFHR
         IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)  = 3
         IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITRDSW
         ELSE
           ID(18)  = IFHR-IFINCR
           IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
! CFS labels time ave fields as inst in long range forecast
         IF(ITRDSW < 0)ID(1:25)=0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(401))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!
!     TIME AVERAGED SURFACE VISIBLE DIFFUSE DOWNWARD SOLAR FLUX
      IF (IGET(402).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = AVISDIFFSWIN(I,J)
         ENDDO
         ENDDO
         ID(1:25)=0
         ITRDSW     = NINT(TRDSW)
         IF(ITRDSW .ne. 0) then
           IFINCR     = MOD(IFHR,ITRDSW)
           IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
         ELSE
           IFINCR     = 0
         endif
         ID(19)  = IFHR
         IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)  = 3
         IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITRDSW
         ELSE
           ID(18)  = IFHR-IFINCR
           IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         IF(ITRDSW < 0)ID(1:25)=0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(402))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!
!     TIME AVERAGED SURFACE VISIBLE BEAM DOWNWARD SOLAR FLUX
      IF (IGET(403).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = AIRBEAMSWIN(I,J)
         ENDDO
         ENDDO
         ID(1:25)=0
         ITRDSW     = NINT(TRDSW)
         IF(ITRDSW .ne. 0) then
           IFINCR     = MOD(IFHR,ITRDSW)
           IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
         ELSE
           IFINCR     = 0
         endif
         ID(19)  = IFHR
         IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)  = 3
         IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITRDSW
         ELSE
           ID(18)  = IFHR-IFINCR
           IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         IF(ITRDSW < 0)ID(1:25)=0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(403))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!
!     TIME AVERAGED SURFACE VISIBLE DIFFUSE DOWNWARD SOLAR FLUX
      IF (IGET(404).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = AIRDIFFSWIN(I,J)
         ENDDO
         ENDDO
         ID(1:25)=0
         ITRDSW     = NINT(TRDSW)
         IF(ITRDSW .ne. 0) then
           IFINCR     = MOD(IFHR,ITRDSW)
           IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
         ELSE
           IFINCR     = 0
         endif
         ID(19)  = IFHR
         IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)  = 3
         IF (IFINCR.EQ.0) THEN
           ID(18)  = IFHR-ITRDSW
         ELSE
           ID(18)  = IFHR-IFINCR
           IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         IF(ITRDSW < 0)ID(1:25)=0
         if(grib=="grib2" )then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(404))
            if(ITRDSW>0) then
               fld_info(cfld)%ntrange=1
            else
               fld_info(cfld)%ntrange=0
            endif
            fld_info(cfld)%tinvstat=IFHR-ID(18)
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

      !2D AEROSOL OPTICAL DEPTH AT 550 NM
      IF (IGET(715).GT.0) THEN
         DO J=JSTA,JEND
           DO I=1,IM
             grid1(i,j)=taod5502d(i,j)
           ENDDO
         ENDDO
         if(grib=="grib2" )then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(715))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
   
      !AEROSOL ASYMMETRY FACTOR
      IF (IGET(716).GT.0) THEN
         DO J=JSTA,JEND
           DO I=1,IM
             grid1(i,j)=aerasy2d(i,j)
           ENDDO
         ENDDO
         if(grib=="grib2" )then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(716))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
   
      !AEROSOL SINGLE-SCATTERING ALBEDO
      IF (IGET(717).GT.0) THEN
         DO J=JSTA,JEND
           DO I=1,IM
             grid1(i,j)=aerssa2d(i,j)
           ENDDO
         ENDDO
         if(grib=="grib2" )then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(717))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!
      if (gocart_on) then
!
!***  BLOCK 4. GOCART AEROSOL FIELDS
!
!
!!! ADD AOD AT 550 NM AND OTHER CHANNELS

!! Aerosol Optical Depth (AOD)
!! CALPW(dust mixing ratio in kg/kg) => column density (kg/m2)
!! CALPW(dust mixing ratio in kg/kg * Qext [aerosol extinction efficiency
!! in m2/g] * 1000. [convert m2/g to m2/kg]) =>  AOD (no unit)
!!

!! DETERMINE WHETHER TO COMPUTE AEROSOL OPTICAL PROPERTIES
        LAEROPT = .FALSE.
        DO I = 609, 614   ! TOTAL AND SPECIATED AOD AT 550NM
          IF  ( IGET(I).GT.0 ) LAEROPT = .TRUE.
        ENDDO
        DO I = 623, 628   ! AOD AT MULTI-CHANNELS
          IF  ( IGET(I).GT.0 ) LAEROPT = .TRUE.
        ENDDO
        DO I = 648, 656   ! (SSA, ASY AT 340),(SCA AT 550), ANGSTROM
          IF  ( IGET(I).GT.0 ) LAEROPT = .TRUE.
        ENDDO

!! DETERMINE WHETHER TO COMPUTE INSTANT SURFACE MASS CONC
        LAERSMASS = .FALSE.
        DO I = 690, 698   ! TOTAL AND SPECIATED AEROSOL
          IF  ( IGET(I).GT.0 ) LAERSMASS = .TRUE.
        ENDDO

        IF ( LAEROPT ) THEN
         PRINT *, 'COMPUTE AEROSOL OPTICAL PROPERTIES'

!!! ALLOCATE AEROSOL OPTICAL PROPERTIES
         ALLOCATE ( extrhd_DU(KRHLEV,nbin_du,NBDSW))
         ALLOCATE ( extrhd_SS(KRHLEV,nbin_ss,NBDSW))
         ALLOCATE ( extrhd_SU(KRHLEV,nbin_su,NBDSW))
         ALLOCATE ( extrhd_BC(KRHLEV,nbin_bc,NBDSW))
         ALLOCATE ( extrhd_OC(KRHLEV,nbin_oc,NBDSW))

         ALLOCATE ( scarhd_DU(KRHLEV,nbin_du,NBDSW))
         ALLOCATE ( scarhd_SS(KRHLEV,nbin_ss,NBDSW))
         ALLOCATE ( scarhd_SU(KRHLEV,nbin_su,NBDSW))
         ALLOCATE ( scarhd_BC(KRHLEV,nbin_bc,NBDSW))
         ALLOCATE ( scarhd_OC(KRHLEV,nbin_oc,NBDSW))

         ALLOCATE ( asyrhd_DU(KRHLEV,nbin_du,NBDSW))
         ALLOCATE ( asyrhd_SS(KRHLEV,nbin_ss,NBDSW))
         ALLOCATE ( asyrhd_SU(KRHLEV,nbin_su,NBDSW))
         ALLOCATE ( asyrhd_BC(KRHLEV,nbin_bc,NBDSW))
         ALLOCATE ( asyrhd_OC(KRHLEV,nbin_oc,NBDSW))

         ALLOCATE ( ssarhd_DU(KRHLEV,nbin_du,NBDSW))
         ALLOCATE ( ssarhd_SS(KRHLEV,nbin_ss,NBDSW))
         ALLOCATE ( ssarhd_SU(KRHLEV,nbin_su,NBDSW))
         ALLOCATE ( ssarhd_BC(KRHLEV,nbin_bc,NBDSW))
         ALLOCATE ( ssarhd_OC(KRHLEV,nbin_oc,NBDSW))
         PRINT *, 'aft  AEROSOL allocate, nbin_du=',nbin_du,  &
          'nbin_ss=',nbin_ss,'nbin_su=',nbin_su,'nbin_bc=',     &
          'nbin_oc=',nbin_oc,'nAero=',nAero

!!! READ AEROSOL LUTS
         DO i = 1, nAero
          CLOSE(UNIT=NOAER)
          aerosol_file='optics_luts_'//AerosolName(i)//'.dat'
          open(unit=NOAER, file=aerosol_file, status='OLD', iostat=ios)
          IF (IOS .GT. 0) THEN
            print *,' ERROR! Non-zero iostat for rd_LUTS ', aerosol_file
            stop
          ENDIF
          print *,'i=',i,'read aerosol_file=',trim(aerosol_file),'ios=',ios
!
          IF (AerosolName(i) .EQ. 'DUST') nbin = nbin_du
          IF (AerosolName(i) .EQ. 'SALT') nbin = nbin_ss
          IF (AerosolName(i) .EQ. 'SUSO') nbin = nbin_su
          IF (AerosolName(i) .EQ. 'SOOT') nbin = nbin_bc
          IF (AerosolName(i) .EQ. 'WASO') nbin = nbin_oc
          DO J = 1, NBIN
           read(NOAER,'(2x,a4,1x,i1,1x,a3)')AerName_rd,ib, AerOpt
           IF (AerName_rd .ne. AerosolName(i)) STOP
           IF (j .ne.  ib                        ) STOP
           IF (AerOpt .ne. 'ext'                 ) STOP

           IF (AerosolName(i) .EQ. 'DUST') THEN
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (extrhd_du(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (scarhd_du(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (asyrhd_du(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (ssarhd_du(ii,j,ib), ii=1,KRHLEV)
            enddo

           ELSEIF (AerosolName(i) .EQ. 'SALT') THEN
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (extrhd_ss(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (scarhd_ss(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (asyrhd_ss(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (ssarhd_ss(ii,j,ib), ii=1,KRHLEV)
            enddo

           ELSEIF (AerosolName(i) .EQ. 'SUSO') THEN
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (extrhd_su(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (scarhd_su(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (asyrhd_su(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (ssarhd_su(ii,j,ib), ii=1,KRHLEV)
            enddo

           ELSEIF (AerosolName(i) .EQ. 'SOOT') THEN
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (extrhd_bc(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (scarhd_bc(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (asyrhd_bc(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (ssarhd_bc(ii,j,ib), ii=1,KRHLEV)
            enddo

           ELSEIF (AerosolName(i) .EQ. 'WASO') THEN
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (extrhd_oc(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (scarhd_oc(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (asyrhd_oc(ii,j,ib), ii=1,KRHLEV)
            enddo
            read(NOAER,'(2x,a4)')  AerName_rd
            do ib = 1, NBDSW
             read(NOAER,'(8f10.5)') (ssarhd_oc(ii,j,ib), ii=1,KRHLEV)
            enddo

           ENDIF

         ENDDO        ! j-loop for nbin
        ENDDO        ! i-loop for nAero
        print *,'finish reading coef'

        CLOSE(UNIT=NOAER)

!!! COMPUTES RELATIVE HUMIDITY AND RDRH
!       allocate (RH3D(im,jsta:jend,lm))
        allocate (rdrh(im,jsta:jend,lm))
        allocate (ihh(im,jsta:jend,lm))
        DO L=1,LM                    ! L FROM TOA TO SFC
          LL=LM-L+1                  ! LL FROM SFC TO TOA
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              P1D(I,J) = PMID(I,J,LL)
              T1D(I,J) = T(I,J,LL)
              Q1D(I,J) = Q(I,J,LL)
            ENDDO
          ENDDO
          CALL CALRH_GFS(P1D,T1D,Q1D,EGRID4)
          DO J=JSTA,JEND
            DO I=1,IM
!             RH3D(I,J,LL) = EGRID4(I,J)
              RH3D         = EGRID4(I,J)
!   DETERMINE RDRH (wgt for IH2) and IHH (index for IH2)
!             IF ( RH3D(I,J,LL) > RHLEV(KRHLEV) ) THEN
              IF ( RH3D > RHLEV(KRHLEV) ) THEN
                IH2 = KRHLEV
                IH1 = IH2 - 1
                RDRH(I,J,LL) = 1.
!             ELSEIF ( RH3D(I,J,LL) < RHLEV(1)) THEN
              ELSEIF ( RH3D < RHLEV(1)) THEN
                IH2 = 1
                IH1 = 1
                RDRH(I,J,LL) = 0.
              ELSE
                IH2 = 1
!               DO WHILE ( RH3D(I,J,LL) > RHLEV(IH2))
                DO WHILE ( RH3D > RHLEV(IH2))
                  IH2 = IH2 + 1
                  IF ( IH2 > KRHLEV ) EXIT
                ENDDO
                IH2 = MIN( KRHLEV, IH2 )
                IH1 = MAX( 1, IH2-1 )
                DRH0 = RHLEV(IH2) - RHLEV(IH1)
!               DRH1 = RH3D(I,J,LL) - RHLEV(IH1)
                DRH1 = RH3D - RHLEV(IH1)
                RDRH(I,J,LL) = DRH1 / DRH0
              ENDIF
              IHH(I,J,LL) = IH1
!
            ENDDO
          ENDDO
        ENDDO

!!!
!!! COMPUTE AOD FOR SPECIFIED WAVELENGTHS
        DO IB = 1, NBDSW

! AOD AT 340 NM
        IF (IB .EQ. 1 )  INDX = 623
! AOD AT 440 NM
        IF (IB .EQ. 2 )  INDX = 624
! AOD AT 550 NM
        IF (IB .EQ. 3 )  INDX = 609
! AOD AT 660 NM
        IF (IB .EQ. 4 )  INDX = 625
! AOD AT 860 NM
        IF (IB .EQ. 5 )  INDX = 626
! AOD AT 1630 NM
        IF (IB .EQ. 6 )  INDX = 627
! AOD AT 11100 NM
        IF (IB .EQ. 7 )  INDX = 628

! DETERMINE LEXT AND LSCA (DEFAULT TO F)
        LEXT = .FALSE.
        LSCA = .FALSE.
        LASY = .FALSE.
! -- CHECK WHETHER TOTAL EXT AOD IS REQUESTED
        IF (IGET(INDX).GT.0 ) LEXT =.TRUE.
! -- CHECK WHETHER SPECIATED AOD AT 550 NM IS REQUESTED
        IF ( IB .EQ. 3 ) THEN
          IF (IGET(650).GT.0  ) LSCA =.TRUE.    !TOTAL SCA AOD
          DO I = 1, nAero
            IF (IGET(INDX_EXT(I)).GT.0 ) LEXT = .TRUE.
            IF (IGET(INDX_SCA(I)).GT.0 ) LSCA = .TRUE.
          ENDDO
        ENDIF
! -- CHECK WHETHER ASY AND SSA AT 340NM IS REQUESTED
        IF ( IB .EQ. 1 ) THEN
          IF (IGET(648).GT.0  ) LSCA =.TRUE.
          IF (IGET(649).GT.0  ) LASY =.TRUE.
        ENDIF
! -- CHECK WHETHER ANGSTROM EXPONENT IS REQUESTED
        IF (IGET(656).GT.0 ) THEN
          IF ( IB .EQ. 2 ) LEXT = .TRUE.
          IF ( IB .EQ. 5 ) LEXT = .TRUE.
        ENDIF
        print *,'LEXT=',LEXT,'LSCA=',LSCA,'LASY=',LASY
! SKIP IF POST PRODUCT IS NOT REQUESTED
        IF ( LEXT .OR. LSCA .OR. LASY ) THEN
! COMPUTE DUST AOD
         AOD_DU=SPVAL
         SCA_DU=SPVAL
         ASY_DU=SPVAL
         EXT=0.0
         SCA=0.0
         ASY=0.0
         DO  J=JSTA,JEND
           DO  I=1,IM
             DO  L=1,LM
               DO N=1, NBIN_DU
               EXT01 = EXTRHD_DU(1,N,IB)
               SCA01 = SCARHD_DU(1,N,IB)
               ASY01 = ASYRHD_DU(1,N,IB)
               EXT(I,J,L) = EXT(I,J,L)+1e-9*DUST(I,J,L,N) * EXT01
               SCA(I,J,L) = SCA(I,J,L)+1e-9*DUST(I,J,L,N) * SCA01
               ASY(I,J,L) = ASY(I,J,L)+1e-9*DUST(I,J,L,N) * SCA01*ASY01
               ENDDO
               EXT(I,J,L) = EXT(I,J,L) * 1000.
               SCA(I,J,L) = SCA(I,J,L) * 1000.
               ASY(I,J,L) = ASY(I,J,L) * 1000.
             ENDDO  ! L-loop
           ENDDO    ! I-loop
         ENDDO      ! J-loop
         CALL CALPW(AOD_DU,17)
         CALL CALPW(SCA_DU,20)
         CALL CALPW(ASY_DU,21)
! COMPUTE SULFATE AOD
         AOD_SU=SPVAL
         SCA_SU=SPVAL
         ASY_SU=SPVAL
         EXT=0.0
         SCA=0.0
         ASY=0.0
         DO  J=JSTA,JEND
           DO  I=1,IM
             DO  L=1,LM
               ih1 = ihh(I,J,L)
               ih2 = ih1 + 1
               DO N = 1, NBIN_SU
               EXT01 = EXTRHD_SU(IH1,N,IB)                                &
     &          + RDRH(I,J,L)*(EXTRHD_SU(IH2,N,IB)-EXTRHD_SU(IH1,N,IB))
               SCA01 = SCARHD_SU(IH1,N,IB)                                &
     &          + RDRH(I,J,L)*(SCARHD_SU(IH2,N,IB)-SCARHD_SU(IH1,N,IB))
               ASY01 = ASYRHD_SU(IH1,N,IB)                                &
     &          + RDRH(I,J,L)*(ASYRHD_SU(IH2,N,IB)-ASYRHD_SU(IH1,N,IB))
          EXT(I,J,L) = EXT(I,J,L)+1e-9*SUSO(I,J,L,N) * EXT01
          SCA(I,J,L) = SCA(I,J,L)+1e-9*SUSO(I,J,L,N)*SCA01
          ASY(I,J,L) = ASY(I,J,L)+1e-9*SUSO(I,J,L,N)*SCA01*ASY01

               ENDDO   ! N-loop
               EXT(I,J,L) = EXT(I,J,L) * 1000.
               SCA(I,J,L) = SCA(I,J,L) * 1000.
               ASY(I,J,L) = ASY(I,J,L) * 1000.
             ENDDO  ! L-loop
           ENDDO    ! I-loop
         ENDDO      ! J-loop
         CALL CALPW(AOD_SU,17)
         CALL CALPW(SCA_SU,20)
         CALL CALPW(ASY_SU,21)

! COMPUTE SEA SALT AOD
         AOD_SS=SPVAL
         SCA_SS=SPVAL
         ASY_SS=SPVAL
         EXT=0.0
         SCA=0.0
         ASY=0.0
         DO  J=JSTA,JEND
           DO  I=1,IM
             DO  L=1,LM
               ih1 = ihh (I,J,L)
               ih2 = ih1 + 1
             DO N = 1, NBIN_SS
               EXT01 = EXTRHD_SS(IH1,N,IB) &
     &          + RDRH(I,J,L)*(EXTRHD_SS(IH2,N,IB)-EXTRHD_SS(IH1,N,IB))
               SCA01 = SCARHD_SS(IH1,N,IB) &
     &          + RDRH(I,J,L)*(SCARHD_SS(IH2,N,IB)-SCARHD_SS(IH1,N,IB))
               ASY01 = ASYRHD_SS(IH1,N,IB) &
     &          + RDRH(I,J,L)*(ASYRHD_SS(IH2,N,IB)-ASYRHD_SS(IH1,N,IB))
               EXT(I,J,L) = EXT(I,J,L)+1e-9*SALT(I,J,L,N)*EXT01
               SCA(I,J,L) = SCA(I,J,L)+1e-9*SALT(I,J,L,N)*SCA01
               ASY(I,J,L) = ASY(I,J,L)+1e-9*SALT(I,J,L,N)*SCA01*ASY01
             ENDDO   ! N-loop
               EXT(I,J,L) = EXT(I,J,L) * 1000.
               SCA(I,J,L) = SCA(I,J,L) * 1000.
               ASY(I,J,L) = ASY(I,J,L) * 1000.
             ENDDO  ! L-loop
           ENDDO    ! I-loop
         ENDDO      ! J-loop
         CALL CALPW(AOD_SS,17)
         CALL CALPW(SCA_SS,20)
         CALL CALPW(ASY_SS,21)

! COMPUTE BLACK CARBON AOD
         AOD_BC=SPVAL
         SCA_BC=SPVAL
         ASY_BC=SPVAL
         EXT=0.0
         SCA=0.0
         ASY=0.0
         DO  J=JSTA,JEND
           DO  I=1,IM
             DO  L=1,LM
             ih1 = ihh (I,J,L)
             ih2 = ih1 + 1
             DO N = 1, NBIN_BC
               EXT01 = EXTRHD_BC(IH1,N,IB)  &
     &          + RDRH(I,J,L)*(EXTRHD_BC(IH2,N,IB)-EXTRHD_BC(IH1,N,IB))
               SCA01 = SCARHD_BC(IH1,N,IB)  &
     &          + RDRH(I,J,L)*(SCARHD_BC(IH2,N,IB)-SCARHD_BC(IH1,N,IB))
               ASY01 = ASYRHD_BC(IH1,N,IB)  &
     &          + RDRH(I,J,L)*(ASYRHD_BC(IH2,N,IB)-ASYRHD_BC(IH1,N,IB))
               EXT(I,J,L) = EXT(I,J,L)+1e-9*SOOT(I,J,L,N)*EXT01
               SCA(I,J,L) = SCA(I,J,L)+1e-9*SOOT(I,J,L,N)*SCA01
               ASY(I,J,L) = ASY(I,J,L)+1e-9*SOOT(I,J,L,N)*SCA01*ASY01
             ENDDO   ! N-loop
               EXT(I,J,L) = EXT(I,J,L) * 1000.
               SCA(I,J,L) = SCA(I,J,L) * 1000.
               ASY(I,J,L) = ASY(I,J,L) * 1000.
             ENDDO  ! L-loop
           ENDDO    ! I-loop
         ENDDO      ! J-loop
         CALL CALPW(AOD_BC,17)
         CALL CALPW(SCA_BC,20)
         CALL CALPW(ASY_BC,21)
! COMPUTE ORGANIC CARBON AOD
         AOD_OC=SPVAL
         SCA_OC=SPVAL
         ASY_OC=SPVAL
         EXT=0.0
         SCA=0.0
         ASY=0.0
         DO  J=JSTA,JEND
           DO  I=1,IM
             DO  L=1,LM
             ih1 = ihh (I,J,L)
             ih2 = ih1 + 1
             DO N = 1, NBIN_OC
               EXT01 = EXTRHD_OC(IH1,N,IB) &
     &          + RDRH(I,J,L)*(EXTRHD_OC(IH2,N,IB)-EXTRHD_OC(IH1,N,IB))
               SCA01 = SCARHD_OC(IH1,N,IB) &
     &          + RDRH(I,J,L)*(SCARHD_OC(IH2,N,IB)-SCARHD_OC(IH1,N,IB))
               ASY01 = ASYRHD_OC(IH1,N,IB) &
     &          + RDRH(I,J,L)*(ASYRHD_OC(IH2,N,IB)-ASYRHD_OC(IH1,N,IB))
               EXT(I,J,L) = EXT(I,J,L)+1e-9*WASO(I,J,L,N)*EXT01
               SCA(I,J,L) = SCA(I,J,L)+1e-9*WASO(I,J,L,N)*SCA01
               ASY(I,J,L) = ASY(I,J,L)+1e-9*WASO(I,J,L,N)*SCA01*ASY01
             ENDDO   ! N-loop
               EXT(I,J,L) = EXT(I,J,L) * 1000.
               SCA(I,J,L) = SCA(I,J,L) * 1000.
               ASY(I,J,L) = ASY(I,J,L) * 1000.
             ENDDO  ! L-loop
           ENDDO    ! I-loop
         ENDDO      ! J-loop
         CALL CALPW(AOD_OC,17)
         CALL CALPW(SCA_OC,20)
         CALL CALPW(ASY_OC,21)

! COMPUTE TOTAL AOD
         AOD=SPVAL
         SCA=SPVAL
         ASY=SPVAL
         DO  J=JSTA,JEND
           DO  I=1,IM
             AOD_DU(I,J) = MAX (AOD_DU(I,J), 0.0)
             AOD_BC(I,J) = MAX (AOD_BC(I,J), 0.0)
             AOD_OC(I,J) = MAX (AOD_OC(I,J), 0.0)
             AOD_SU(I,J) = MAX (AOD_SU(I,J), 0.0)
             AOD_SS(I,J) = MAX (AOD_SS(I,J), 0.0)

            SCA_DU(I,J) = MAX (SCA_DU(I,J), 0.0)
            SCA_BC(I,J) = MAX (SCA_BC(I,J), 0.0)
            SCA_OC(I,J) = MAX (SCA_OC(I,J), 0.0)
            SCA_SU(I,J) = MAX (SCA_SU(I,J), 0.0)
            SCA_SS(I,J) = MAX (SCA_SS(I,J), 0.0)

            ASY_DU(I,J) = MAX (ASY_DU(I,J), 0.0)
            ASY_BC(I,J) = MAX (ASY_BC(I,J), 0.0)
            ASY_OC(I,J) = MAX (ASY_OC(I,J), 0.0)
            ASY_SU(I,J) = MAX (ASY_SU(I,J), 0.0)
            ASY_SS(I,J) = MAX (ASY_SS(I,J), 0.0)

            AOD(I,J)    = AOD_DU(I,J) + AOD_BC(I,J) + AOD_OC(I,J) +   &
     &                     AOD_SU(I,J) + AOD_SS(I,J)
            SCA2D(I,J) = SCA_DU(I,J) + SCA_BC(I,J) + SCA_OC(I,J) +    &
     &                 SCA_SU(I,J) + SCA_SS(I,J)
            ASY2D(I,J) = ASY_DU(I,J) + ASY_BC(I,J) + ASY_OC(I,J) +    &
     &                 ASY_SU(I,J) + ASY_SS(I,J)
           ENDDO   ! I-loop
         ENDDO   ! J-loop
! FILL UP AOD_440 AND AOD_860, IF ANGSTROM EXP IS REQUESTED
         IF ( IGET(656) .GT. 0 )  THEN
          IF (IB .EQ. 2 ) THEN              !! AOD AT 440 NM
!$omp parallel do private(i,j)
           DO  J=JSTA,JEND
           DO  I=1,IM
             AOD_440(I,J) = AOD(I,J)
           ENDDO   ! I-loop
           ENDDO   ! J-loop
          ENDIF

          IF (IB .EQ. 5 ) THEN              !! AOD AT 860 NM
!$omp parallel do private(i,j)
           DO  J=JSTA,JEND
           DO  I=1,IM
             AOD_860(I,J) = AOD(I,J)
           ENDDO   ! I-loop
           ENDDO   ! J-loop
          ENDIF
        ENDIF

! WRITE OUT TOTAL AOD
        IF ( IGET(INDX) .GT. 0)  THEN
        ID(1:25)=0
        ID(02)=141
!$omp parallel do private(i,j)
        do j=jsta,jend
          do i=1,im
            GRID1(i,j) = AOD(i,j)
          enddo
        enddo
        CALL BOUND(GRID1,D00,H99999)
        if(grib=="grib2" )then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(INDX))
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
        endif
        ENDIF
!
! WRITE OUT ASY AND SSA AT 340NM
        IF ( IB .EQ. 1 ) THEN                 !!! FOR 340NM ONLY

! AER ASYM FACTOR AT 340 NM
          IF ( IGET(649) .GT. 0 )  THEN
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
          DO I=1,IM
            IF ( SCA2D(I,J) > 0.0 ) THEN
             ASY2D(I,J) = ASY2D(I,J) / SCA2D(I,J)
            ELSE
             ASY2D(I,J) = 0.
            ENDIF
          GRID1(I,J)=ASY2D(I,J)
          ENDDO
          ENDDO
          ID(1:25)=0
          ID(02)=141
          CALL BOUND(GRID1,D00,H99999)
          if(grib=="grib2" )then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(649))
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
          endif
          ENDIF       ! IGET(649)

! AER SINGLE SCATTER ALB AT 340 NM
          IF ( IGET(648) .GT. 0 )  THEN
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
          DO I=1,IM
            IF ( AOD(I,J) > 0.0 ) THEN
             SCA2D(I,J) = SCA2D(I,J) / AOD(I,J)
            ELSE
             SCA2D(I,J) = 1.0
            ENDIF
             GRID1(I,J)=SCA2D(I,J)
          ENDDO
          ENDDO
          ID(1:25)=0
          ID(02)=141
          CALL BOUND(GRID1,D00,H99999)
          if(grib=="grib2" )then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(648))
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
          endif
          ENDIF        ! IGET(648)
        print *,'aft compute sca340'

        ENDIF       ! IB IF-BLOCK (340NM)


! WRITE OUT AOD FOR DU, SU, SS, OC, BC for all wavelengths
! WRITE OUT SPECIATED AEROSOL OPTICAL PROPERTIES
        IF ( IB .EQ. 3 ) THEN                 !!! FOR 550NM ONLY

! WRITE OUT TOTAL SCATTERING AOD
          IF ( IGET(650) .GT. 0 )  THEN
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
          DO I=1,IM
          GRID1(I,J)=SCA2D(I,J)
          ENDDO
          ENDDO
          ID(1:25)=0
          ID(02)=141
          CALL BOUND(GRID1,D00,H99999)
          if(grib=="grib2" )then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(650))
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
          endif
          ENDIF
! LOOP THROUGH EACH SPECIES
          DO II = 1, nAero

! WRITE OUT EXT AOD
            JJ = INDX_EXT(II)
            IF ( IGET(JJ) .GT. 0)  THEN          ! EXT AOD
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
          DO I=1,IM
             IF ( II .EQ. 1 ) GRID1(I,J) = AOD_DU(I,J)
             IF ( II .EQ. 2 ) GRID1(I,J) = AOD_SS(I,J)
             IF ( II .EQ. 3 ) GRID1(I,J) = AOD_SU(I,J)
             IF ( II .EQ. 4 ) GRID1(I,J) = AOD_OC(I,J)
             IF ( II .EQ. 5 ) GRID1(I,J) = AOD_BC(I,J)
          ENDDO
          ENDDO
             ID(1:25)=0
             ID(02)=141
             CALL BOUND(GRID1,D00,H99999)
             if(grib=="grib2" )then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(JJ))
               datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
             endif
            ENDIF

! WRITE OUT SCA AOD
            JJ = INDX_SCA(II)
            IF ( IGET(JJ) .GT. 0)  THEN          ! SCA AOD
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
          DO I=1,IM
             IF ( II .EQ. 1 ) GRID1(I,J) = SCA_DU(I,J)
             IF ( II .EQ. 2 ) GRID1(I,J) = SCA_SS(I,J)
             IF ( II .EQ. 3 ) GRID1(I,J) = SCA_SU(I,J)
             IF ( II .EQ. 4 ) GRID1(I,J) = SCA_OC(I,J)
             IF ( II .EQ. 5 ) GRID1(I,J) = SCA_BC(I,J)
          ENDDO
          ENDDO
             ID(1:25)=0
             ID(02)=141
             CALL BOUND(GRID1,D00,H99999)
             if(grib=="grib2" )then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(JJ))
               datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
             endif
            ENDIF

          ENDDO     ! II DO-LOOP

        ENDIF       ! IB IF-BLOCK (550NM)

        ENDIF       ! LEXT IF-BLOCK
        ENDDO       ! LOOP THROUGH NBDSW CHANNELS
! COMPUTE AND WRITE OUT ANGSTROM EXPONENT
        IF ( IGET(656) .GT. 0 )  THEN
          ANGST=SPVAL
!          ANG2 = LOG ( 0.860 / 0.440 )
          ANG2 = LOG ( 860. / 440. )
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
          DO I=1,IM
            IF (AOD_860(I,J) .GT. 0.) THEN
             ANG1 = LOG( AOD_440(I,J)/AOD_860(I,J) )
             ANGST(I,J) = ANG1 / ANG2
            ENDIF
          GRID1(I,J)=ANGST(I,J)
          ENDDO
          ENDDO
          print *,'output angstrom exp,angst=',maxval(angst(1:im,jsta:jend)), &
            minval(angst(1:im,jsta:jend))
          ID(1:25)=0
          ID(02)=141
          CALL BOUND(GRID1,D00,H99999)
          if(grib=="grib2" )then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(656))
            datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
          endif
        ENDIF       ! ANGSTROM EXPONENT

      ENDIF         ! END OF LAEROPT IF-BLOCK

!! Multiply by 1.E-6 to revert these fields back
      IF (IGET(659).GT.0) THEN
         GRID1=SPVAL
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               GRID1(I,J) = DUEM(I,J,1)*1.E-6
               DO K=2,NBIN_DU
                GRID1(I,J) = GRID1(I,J) + DUEM(I,J,K)*1.E-6
               END DO
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=141
         if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(659))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

      IF (IGET(660).GT.0) THEN
         GRID1=SPVAL
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               GRID1(I,J) = DUSD(I,J,1)*1.E-6
               DO K=2,NBIN_DU
                GRID1(I,J) = GRID1(I,J)+ DUSD(I,J,K)*1.E-6
               END DO
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=141
         if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(660))
          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!! ADD DUST DRY DEPOSITION FLUXES (kg/m2/sec)
!
!      IF (IGET(661).GT.0) THEN
!         DO J = JSTA,JEND
!            DO I = 1,IM
!               GRID1(I,J) = DUDP(I,J,1)*1.E-6
!               DO K=2,NBIN_DU
!                GRID1(I,J) = GRID1(I,J)+ DUDP(I,J,K)*1.E-6
!               END DO
!            END DO
!         END DO
!         ID(1:25) = 0
!         ID(02)=141
!         if(grib=='grib1') then
!          CALL GRIBIT(IGET(661),LVLS(1,IGET(661)),GRID1,IM,JM)
!         elseif(grib=='grib2') then
!          cfld=cfld+1
!          fld_info(cfld)%ifld=IAVBLFLD(IGET(661))
!          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
!         endif
!      ENDIF

!! ADD AEROSOL SURFACE PM25 DUST MASS CONCENTRATION (ug/m3)
      IF (IGET(686).GT.0 ) THEN
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               !GRID1(I,J) = DUSMASS(I,J) * 1.E-6
               GRID1(I,J) = DUSTPM(I,J)   !ug/m3 
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=129
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(686))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

!! ADD DUST WET DEPOSITION FLUXES (kg/m2/sec)
!      IF (IGET(662).GT.0) THEN
!         DO J = JSTA,JEND
!            DO I = 1,IM
!               GRID1(I,J) = DUWT(I,J,1)*1.E-6
!               DO K=2,NBIN_DU
!                GRID1(I,J) = GRID1(I,J)+ DUWT(I,J,K)*1.E-6
!               END DO
!            END DO
!         END DO
!         ID(1:25) = 0
!         ID(02)=141
!         if(grib=='grib1') then
!          CALL GRIBIT(IGET(662),LVLS(1,IGET(662)),GRID1,IM,JM)
!         elseif(grib=='grib2') then
!          cfld=cfld+1
!          fld_info(cfld)%ifld=IAVBLFLD(IGET(662))
!          datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
!         endif
!      ENDIF

!! ADD AEROSOL SURFACE PM25 SEA SALT MASS CONCENTRATION (ug/m3)
      IF (IGET(684).GT.0 ) THEN
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               !GRID1(I,J) = DUSMASS(I,J) * 1.E-6
               GRID1(I,J) = SSPM(I,J)   !ug/m3 
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=129
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(684))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!! ADD AEROSOL SURFACE PM10 MASS CONCENTRATION (ug/m3)
      IF (IGET(619).GT.0 ) THEN
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               !GRID1(I,J) = DUSMASS(I,J) * 1.E-6
               GRID1(I,J) = DUSMASS(I,J)   !ug/m3 
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=129
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(619))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

!! ADD AEROSOL SURFACE PM2.5 MASS CONCENTRATION (ug/m3)
      IF (IGET(620).GT.0 ) THEN
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               !GRID1(I,J) = DUSMASS25(I,J) * 1.E-6
               GRID1(I,J) = DUSMASS25(I,J) ! ug/m3 
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=129
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(620))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!! ADD TOTAL AEROSOL PM10 COLUMN DENSITY (kg/m2) !
      IF (IGET(621).GT.0 ) THEN
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               !GRID1(I,J) = DUCMASS(I,J) * 1.E-6
               GRID1(I,J) = DUCMASS(I,J) * 1.E-9
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=141
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(621))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

!! ADD TOTAL AEROSOL PM2.5 COLUMN DENSITY (kg/m2)  
      IF (IGET(622).GT.0 ) THEN
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               !GRID1(I,J) = DUCMASS25(I,J) * 1.E-6
               GRID1(I,J) = DUCMASS25(I,J) * 1.E-9
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=141
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(622))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

!! ADD DUST PM2.5 COLUMN DENSITY (kg/m2)  
      IF (IGET(646).GT.0 ) THEN
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               GRID1(I,J) = DUSTCB(I,J) * 1.E-9
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=141
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(646))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

!! ADD SEA SALT PM2.5 COLUMN DENSITY (kg/m2)  
      IF (IGET(647).GT.0 ) THEN
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               GRID1(I,J) = SSCB(I,J) * 1.E-9
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=141
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(647))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!! ADD BC COLUMN DENSITY (kg/m2)  
      IF (IGET(616).GT.0 ) THEN
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               GRID1(I,J) = BCCB(I,J) * 1.E-9
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=141
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(616))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

!! ADD OC COLUMN DENSITY (kg/m2)  !
      IF (IGET(617).GT.0 ) THEN
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               GRID1(I,J) = OCCB(I,J) * 1.E-9
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=141
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(617))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF

!! ADD SULF COLUMN DENSITY (kg/m2)  !
      IF (IGET(618).GT.0 ) THEN
!$omp parallel do private(i,j)
         DO J = JSTA,JEND
            DO I = 1,IM
               GRID1(I,J) = SULFCB(I,J) * 1.E-9
            END DO
         END DO
         ID(1:25) = 0
         ID(02)=141
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(618))
           datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
         endif
      ENDIF
!! ADD EMISSION FLUXES,dry depostion, wet/convective depostion (kg/m2/sec)
!! The AER file uses 1.E6 to scale all 2d diagnosis fields
!! Multiply by 1.E-6 to revert these fields back
      IF (IGET(659).GT.0) call wrt_aero_diag(659,nbin_du,duem)
      print *,'aft wrt disg duem'
      IF (IGET(660).GT.0) call wrt_aero_diag(660,nbin_du,dusd)
      IF (IGET(661).GT.0) call wrt_aero_diag(661,nbin_du,dudp)
      IF (IGET(662).GT.0) call wrt_aero_diag(662,nbin_du,duwt)
      IF (IGET(679).GT.0) call wrt_aero_diag(679,nbin_du,dusv)
      print *,'aft wrt disg duwt'

!! wrt SS diag field
      IF (IGET(663).GT.0) call wrt_aero_diag(663,nbin_ss,ssem)
      IF (IGET(664).GT.0) call wrt_aero_diag(664,nbin_ss,sssd)
      IF (IGET(665).GT.0) call wrt_aero_diag(665,nbin_ss,ssdp)
      IF (IGET(666).GT.0) call wrt_aero_diag(666,nbin_ss,sswt)
      IF (IGET(680).GT.0) call wrt_aero_diag(680,nbin_ss,sssv)
      print *,'aft wrt disg sswt'

!! wrt BC diag field
      IF (IGET(667).GT.0) call wrt_aero_diag(667,nbin_bc,bcem)
      IF (IGET(668).GT.0) call wrt_aero_diag(668,nbin_bc,bcsd)
      IF (IGET(669).GT.0) call wrt_aero_diag(669,nbin_bc,bcdp)
      IF (IGET(670).GT.0) call wrt_aero_diag(670,nbin_bc,bcwt)
      IF (IGET(681).GT.0) call wrt_aero_diag(681,nbin_bc,bcsv)
      print *,'aft wrt disg bcwt'

!! wrt OC diag field
      IF (IGET(671).GT.0) call wrt_aero_diag(671,nbin_oc,ocem)
      IF (IGET(672).GT.0) call wrt_aero_diag(672,nbin_oc,ocsd)
      IF (IGET(673).GT.0) call wrt_aero_diag(673,nbin_oc,ocdp)
      IF (IGET(674).GT.0) call wrt_aero_diag(674,nbin_oc,ocwt)
      IF (IGET(682).GT.0) call wrt_aero_diag(682,nbin_oc,ocsv)
      print *,'aft wrt disg ocwt'

!! wrt SU diag field
!      IF (IGET(675).GT.0) call wrt_aero_diag(675,nbin_su,suem)
!      IF (IGET(676).GT.0) call wrt_aero_diag(676,nbin_su,susd)
!      IF (IGET(677).GT.0) call wrt_aero_diag(677,nbin_su,sudp)
!      IF (IGET(678).GT.0) call wrt_aero_diag(678,nbin_su,suwt)
!      print *,'aft wrt disg suwt'
      endif          ! if gocart_on

! CB for WAFS
      if(IGET(473)>0 .or. IGET(474)>0 .or. IGET(475)>0) then
         ! CB cover is derived from CPRAT (same as #272 in SURFCE.f) 
         EGRID1 = SPVAL
         DO J=JSTA,JEND
            DO I=1,IM
               if(AVGCPRATE(I,J) /= SPVAL) then
                  EGRID1(I,J) = AVGCPRATE(I,J)*(1000./DTQ2)
               end if
            END DO
         END DO
         call cb_cover(EGRID1)

         ! CB base(height):derived from convective cloud base (pressure, same as #188 in CLDRAD.f)
         ! CB top(height): derived from convective cloud top (pressure, same as #189 in CLDRAD.f)
         EGRID2 = SPVAL
         EGRID3 = SPVAL
         IF(MODELNAME == 'GFS' .OR. MODELNAME == 'FV3R') then
            DO J=JSTA,JEND
               DO I=1,IM
                  EGRID2(I,J) = PBOT(I,J)
                  EGRID3(I,J) = PTOP(I,J)
               END DO
            END DO
         END IF

         ! Derive CB base and top, relationship among CB fields
         DO J=JSTA,JEND
            DO I=1,IM
               if(EGRID1(I,J)<= 0. .or. EGRID2(I,J)<= 0. .or. EGRID3(I,J) <= 0.) then
                  EGRID1(I,J) = SPVAL
                  EGRID2(I,J) = SPVAL
                  EGRID3(I,J) = SPVAL
               end if
            END DO
         END DO
         DO J=JSTA,JEND
            DO I=1,IM
               IF(EGRID2(I,J) == SPVAL .or. EGRID3(I,J) == SPVAL) cycle
               if(EGRID3(I,J) < 400.*100. .and. &
                  (EGRID2(I,J)-EGRID3(I,J)) > 300.*100) then
                  ! Convert PBOT to height
                  if(EGRID2(I,J) > PMID(I,J,LM)) then
                     EGRID2(I,J) = 0.
                  else
                     do L = LM-1, 1, -1
                        if(EGRID2(I,J) >= PMID(I,J,L)) then
                           if(EGRID2(I,J)-PMID(I,J,L)<0.5) then
                              EGRID2(I,J) = ZMID(I,J,L)
                           else
                              dp = (log(EGRID2(I,J)) - log(PMID(I,J,L)))/ &
                                max(1.e-6,(LOG(PMID(I,J,L+1))-LOG(PMID(I,J,L))))
                              EGRID2(I,J) = ZMID(I,J,L)+(ZMID(I,J,L+1)-ZMID(I,J,L))*dp
                           end if
                           exit
                        end if
                     end do
                  end if
                  ! Convert PTOP to height
                  if(EGRID3(I,J) < PMID(I,J,1)) then
                     EGRID3(I,J) = ZMID(I,J,1)
                  else
                     do L = 2, LM
                        if(EGRID3(I,J) <= PMID(I,J,L)) then
                           if(PMID(I,J,L)-EGRID3(I,J)<0.5) then
                              EGRID3(I,J) = ZMID(I,J,L)
                           else
                              dp = (log(EGRID3(I,J)) - log(PMID(I,J,L)))/ &
                                max(1.e-6,(LOG(PMID(I,J,L))-LOG(PMID(I,J,L-1))))
                              EGRID3(I,J) = ZMID(I,J,L)+(ZMID(I,J,L)-ZMID(I,J,L-1))*dp
                           end if
                           exit
                        end if
                     end do
                  end if
               else
                  EGRID1(I,J) = SPVAL
                  EGRID2(I,J) = SPVAL
                  EGRID3(I,J) = SPVAL
               end if
            END DO
         END DO

         IF(IGET(473) > 0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
               DO I=1,IM
                  GRID1(I,J) = EGRID1(I,J)
               ENDDO
            ENDDO
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(473))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
               enddo
            enddo
         END IF

         IF(IGET(474) > 0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
               DO I=1,IM
                  GRID1(I,J) = EGRID2(I,J)
               ENDDO
            ENDDO
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(474))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
               enddo
            enddo
         END IF

         IF(IGET(475) > 0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
               DO I=1,IM
                  GRID1(I,J) = EGRID3(I,J)
               ENDDO
            ENDDO
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(475))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
               enddo
            enddo
         END IF
      end if 

!
!     END OF ROUTINE.
!
      RETURN
      END

      subroutine cb_cover(cbcov)
!     Calculate CB coverage by using fuzzy logic
!     Evaluate membership of val in a fuzzy set fuzzy.
!     Assume f is in x-log scale
      use ctlblk_mod, only: SPVAL,JSTA,JEND,IM
      implicit none
      real, intent(inout) :: cbcov(IM,JSTA:JEND)

      ! x - convective precipitation [1.0e6*kg/(m2s)]
      ! y - cloud cover fraction, between 0 and 1
      ! These are original values from Slingo (Table 1):
      ! c = -.006 + 0.125*log(p)
      ! x = 1.6 3.6 8.1 18.5 39.0 89.0 197.0 440.0 984.0 10000.0
      ! y = 0.0 0.1 0.2  0.3  0.4  0.5   0.6   0.7   0.8     0.8
      integer, parameter :: NP=10
      real :: x(NP), y(NP)

      integer :: i, j, k
      real :: val, delta

      x = (/ 1.6,3.6,8.1,18.5,39.0,89.0,197.0,440.0,984.0,10000.0 /)   
      y = (/ 0.0,0.1,0.2, 0.3, 0.4, 0.5,  0.6,  0.7,  0.8,    0.8 /)

      x = log(x)

      do j = jsta, jend
      do i = 1, IM
         if(cbcov(i,j) == SPVAL) cycle
         if(cbcov(i,j) <= 0.) then
            cbcov(i,j) = 0.
         else
            val=log(1.0e6*cbcov(i,j))
            if (val <= x(1)) then
               cbcov(i,j) = 0.0
            else if (val >= x(NP)) then
               cbcov(i,j) = 0.0
            else
               do k = 2, NP
                  if (val < x(k)) then
                     delta = x(k) - x(k-1)
                     if (delta <= 0.0) then
                        cbcov(i,j) = y(k-1)
                     else
                        cbcov(i,j) = (y(k) * (val-x(k-1)) + &
                            y(k-1) * (x(k)-val)) / delta
                     end if
                     exit
                  end if
               end do
            end if      
         end if
      end do
      end do
      end subroutine cb_cover

      subroutine wrt_aero_diag(igetfld,nbin,data)
      use ctlblk_mod, only: jsta, jend, SPVAL, im, jm, grib,     &
                  cfld, datapd, fld_info, jsta_2l, jend_2u
      use rqstfld_mod, only: IGET, ID, LVLS, IAVBLFLD
      implicit none
!
      integer igetfld,nbin
      real, dimension(1:im,jsta_2l:jend_2u,nbin) :: data
!
      integer i,j,k
      REAL,dimension(im,jm)    :: GRID1
!
      GRID1=SPVAL
!$omp parallel do private(i,j)
      DO J = JSTA,JEND
        DO I = 1,IM
          grid1(I,J) = data(I,J,1)
          DO K=2,NBIN
            GRID1(I,J) = GRID1(I,J)+ data(I,J,K)
          END DO
        END DO
      END DO
      ID(1:25) = 0
      ID(02)=141
      if(grib=='grib2') then
        cfld=cfld+1
        fld_info(cfld)%ifld=IAVBLFLD(iget(igetfld))
        datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
      endif

      end subroutine wrt_aero_diag
