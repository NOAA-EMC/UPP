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
!
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
      use vrbls3d
      use vrbls2d
      use masks
      use params_mod
      use ctlblk_mod
      use rqstfld_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     SET CELSIUS TO KELVIN CONVERSION.
      real,PARAMETER :: C2K=273.15
!     
!     DECLARE VARIABLES.
!     
      LOGICAL NEED(IM,JM)
      INTEGER L1D(IM,JM) ,lcbot,lctop  !bsf
      INTEGER IBOTT(IM,JM),IBOTCu(IM,JM),IBOTDCu(IM,JM)        &
       ,IBOTSCu(IM,JM),IBOTGr(IM,JM), ITOPT(IM,JM)             &
      ,ITOPCu(IM,JM),ITOPDCu(IM,JM),ITOPSCu(IM,JM)             &
      ,ITOPGr(IM,JM)
      REAL EGRID1(IM,JM),EGRID2(IM,JM),EGRID3(IM,JM)
      REAL GRID1(IM,JM),GRID2(IM,JM),CLDP(IM,JM),              &
              CLDZ(IM,JM),CLDT(IM,JM)                          &
            , RHB(IM,JM,LM)                                    &
            ,watericetotal(LM),watericemax,wimin               &
            ,zcldbase,zcldtop,zpbltop
        real rhoice, coeffp, exponfp, const1, cloud_def_p,     &
             pcldbase, rhoair, vovermd, concfp, betav,         &
             vertvis, tx, tv, pol, esx, es, e, zsf, zcld
        real frac
        real pabovesfc(LM)
        integer nfog, nfogn(7),npblcld                         &
           ,nlifr, k1, k2, ll
      
!     B ZHOU: For aviation:
      REAL  TCLD(IM,JM), CEILING(IM,JM), FLTCND(IM,JM)         &
     &, CU_ir(LM), q_conv   !bsf
!jw
      integer I,J,L,K,IBOT,ITCLOD,LBOT,LTOP,ITRDSW,ITRDLW,     &
              LLMH,ITHEAT,IFINCR,ITYPE,ITOP,NUM_THICK
      real DPBND,RRNUM,QCLD,RSUM,TLMH,FACTRS,FACTRL,DP,        &
           OPDEPTH
      real dummy(IM,JM)
      integer idummy(IM,JM)
!     
!
!*************************************************************************
!     START CLDRAD HERE.
!     
!***  BLOCK 1.  SOUNDING DERIVED FIELDS.
!     
!     ETA SURFACE TO 500MB LIFTED INDEX.  TO BE CONSISTENT WITH THE
!     LFM AND NGM POSTING WE ADD 273.15 TO THE LIFTED INDEX
!
!     THE BEST (SIX LAYER) AND BOUNDARY LAYER LIFTED INDICES ARE
!     COMPUTED AND POSTED IN SUBROUTINE MISCLN.
!
      IF (IGET(030).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           EGRID1(I,J) = SPVAL
         ENDDO
         ENDDO
!
         CALL OTLIFT(EGRID1)
!
         DO J=JSTA,JEND
         DO I=1,IM
           IF(EGRID1(I,J).LT.SPVAL) GRID1(I,J)=EGRID1(I,J) +TFRZ 
         ENDDO
         ENDDO
!
         ID(1:25)=0
         ID(10)  =50
         ID(11)  =100
         CALL GRIBIT(IGET(030),LVLS(1,IGET(030)),GRID1,IM,JM)
      ENDIF
!
!     SOUNDING DERIVED AREA INTEGRATED ENERGIES - CAPE AND CIN.
!       THIS IS THE SFC-BASED CAPE/CIN (lowest 70 mb searched)
      IF ((IGET(032).GT.0).OR.(IGET(107).GT.0)) THEN
         IF ( (LVLS(1,IGET(032)).GT.0) .OR.                         &
              (LVLS(1,IGET(107)).GT.0) ) THEN
            ITYPE = 1
	    DPBND=10.E2
            dummy=0.
            idummy=0
            CALL CALCAPE(ITYPE,DPBND,dummy,dummy,dummy,idummy,EGRID1,EGRID2, &
                 EGRID3,dummy,dummy)
!
!           CONVECTIVE AVAILABLE POTENTIAL ENERGY.
            IF (IGET(032).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = EGRID1(I,J)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               ID(1:25)=0
               CALL GRIBIT(IGET(032),LVLS(1,IGET(032)),GRID1,IM,JM)
            ENDIF
!
!           CONVECTIVE INHIBITION.
            IF (IGET(107).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = -1.*EGRID2(I,J)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = -1.*GRID1(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(107),LVLS(1,IGET(107)),GRID1,IM,JM)
            ENDIF
         ENDIF
      ENDIF
!
!     TOTAL COLUMN PRECIPITABLE WATER (SPECIFIC HUMIDITY).
      IF (IGET(080).GT.0) THEN
         CALL CALPW(GRID1,1)
         ID(1:25)=0
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(080),LVLS(1,IGET(080)),GRID1,IM,JM)
      ENDIF
!     
!     TOTAL COLUMN CLOUD WATER
      IF (IGET(200).GT.0) THEN
         CALL CALPW(GRID1,2)
	 IF(MODELNAME == 'GFS')then
! GFS combines cloud water and cloud ice, hoping to seperate them next implementation	 
	   CALL CALPW(GRID2,3)
	   DO J=JSTA,JEND
           DO I=1,IM
     	     GRID1(I,J)=GRID1(I,J)+GRID2(I,J)
           ENDDO
           ENDDO
	 END IF  
	    
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(200),LVLS(1,IGET(200)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN CLOUD ICE
      IF (IGET(201).GT.0) THEN
         CALL CALPW(GRID1,3)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(201),LVLS(1,IGET(201)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN RAIN 
      IF (IGET(202).GT.0) THEN
         CALL CALPW(GRID1,4)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(202),LVLS(1,IGET(202)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN SNOW 
      IF (IGET(203).GT.0) THEN
         CALL CALPW(GRID1,5)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(203),LVLS(1,IGET(203)),GRID1,IM,JM)
      ENDIF
!
! SRD
!     TOTAL COLUMN GRAUPEL
      IF (IGET(428).GT.0) THEN
         CALL CALPW(GRID1,16)
         ID(1:25)=0
!         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(428),LVLS(1,IGET(428)),GRID1,IM,JM)
      ENDIF
! SRD

!     TOTAL COLUMN CONDENSATE 
      IF (IGET(204).GT.0) THEN
         CALL CALPW(GRID1,6)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(204),LVLS(1,IGET(204)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN SUPERCOOLED (<0C) LIQUID WATER 
      IF (IGET(285).GT.0) THEN
         CALL CALPW(GRID1,7)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(285),LVLS(1,IGET(285)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN MELTING (>0C) ICE
      IF (IGET(286).GT.0) THEN
         CALL CALPW(GRID1,8)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(286),LVLS(1,IGET(286)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN SHORT WAVE T TENDENCY
      IF (IGET(290).GT.0) THEN
         CALL CALPW(GRID1,9)
         ID(1:25)=0
         CALL GRIBIT(IGET(290),LVLS(1,IGET(290)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN LONG WAVE T TENDENCY
      IF (IGET(291).GT.0) THEN
         CALL CALPW(GRID1,10)
         ID(1:25)=0
         CALL GRIBIT(IGET(291),LVLS(1,IGET(291)),GRID1,IM,JM)
      ENDIF            
!
!     TOTAL COLUMN GRID SCALE LATENT HEATING (TIME AVE)
      IF (IGET(292).GT.0) THEN
         CALL CALPW(GRID1,11)
	 IF(AVRAIN.GT.0.)THEN
           RRNUM=1./AVRAIN
         ELSE
           RRNUM=0.
         ENDIF
!$omp  parallel do
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J)=GRID1(I,J)*RRNUM
         ENDDO
         ENDDO
         ID(1:25)=0
	 ITHEAT     = INT(THEAT)
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
         CALL GRIBIT(IGET(292),LVLS(1,IGET(292)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN CONVECTIVE LATENT HEATING (TIME AVE)
      IF (IGET(293).GT.0) THEN
         CALL CALPW(GRID1,12)
	 IF(AVRAIN.GT.0.)THEN
           RRNUM=1./AVCNVC
         ELSE
           RRNUM=0.
         ENDIF
!$omp  parallel do
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J)=GRID1(I,J)*RRNUM
         ENDDO
         ENDDO
         ID(1:25)=0
	 ITHEAT     = INT(THEAT)
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
         CALL GRIBIT(IGET(293),LVLS(1,IGET(293)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN moisture convergence
      IF (IGET(295).GT.0) THEN
         CALL CALPW(GRID1,13)
         ID(1:25)=0
         CALL GRIBIT(IGET(295),LVLS(1,IGET(295)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN RH
      IF (IGET(312).GT.0) THEN
         CALL CALPW(GRID1,14)
         ID(1:25)=0
         CALL GRIBIT(IGET(312),LVLS(1,IGET(312)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN OZONE
      IF (IGET(299).GT.0) THEN
         CALL CALPW(GRID1,15)
         ID(1:25)=0
         CALL GRIBIT(IGET(299),LVLS(1,IGET(299)),GRID1,IM,JM)
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
            ID(1:25)=0
            CALL GRIBIT(IGET(287),LVLS(1,IGET(287)),GRID1,IM,JM)
         ENDIF
         IF (IGET(288).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
               GRID1(I,J)=GRID2(I,J)
            ENDDO
            ENDDO
            ID(1:25)=0
            CALL GRIBIT(IGET(288),LVLS(1,IGET(288)),GRID1,IM,JM)
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
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL GRIBIT(IGET(197),LVLS(1,IGET(197)),GRID1,IM,JM)
      ENDIF
!   
!
!
!***  BLOCK 2.  2-D CLOUD FIELDS.
!
!     LOW CLOUD FRACTION.
      IF (IGET(037).GT.0) THEN
        GRID1=SPVAL	  
        DO J=JSTA,JEND
        DO I=1,IM
	  IF(CFRACL(I,J) < SPVAL) GRID1(I,J) = CFRACL(I,J)*100.
        ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(037),LVLS(1,IGET(037)),GRID1,IM,JM)
      ENDIF
!
!     TIME AVERAGED LOW CLOUD FRACTION.
      IF (IGET(300).GT.0) THEN	
        GRID1=spval  
        DO J=JSTA,JEND
        DO I=1,IM
	  IF(AVGCFRACL(I,J) < SPVAL)                               &  
     &        GRID1(I,J) = AVGCFRACL(I,J)*100.   
        ENDDO
        ENDDO
        ID(1:25)=0
        ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(300),LVLS(1,IGET(300)),GRID1,IM,JM)
      ENDIF      
!     
!     MIDDLE CLOUD FRACTION.
      IF (IGET(038).GT.0) THEN
        GRID1=SPVAL
        DO J=JSTA,JEND
        DO I=1,IM
	   IF(CFRACM(I,J) < SPVAL)GRID1(I,J) = CFRACM(I,J)*100.
        ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(038),LVLS(1,IGET(038)),GRID1,IM,JM)
      ENDIF
!
!     TIME AVERAGED MIDDLE CLOUD FRACTION.
      IF (IGET(301).GT.0) THEN	  
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
        ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(301),LVLS(1,IGET(301)),GRID1,IM,JM)
      ENDIF   
!     
!     HIGH CLOUD FRACTION.
      IF (IGET(039).GT.0) THEN
        GRID1=SPVAL
        DO J=JSTA,JEND
        DO I=1,IM
	   IF(CFRACH(I,J) < SPVAL)GRID1(I,J) = CFRACH(I,J)*100.
        ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(039),LVLS(1,IGET(039)),GRID1,IM,JM)
      ENDIF
!
!     TIME AVERAGED HIGH CLOUD FRACTION.
      IF (IGET(302).GT.0) THEN	  
        DO J=JSTA,JEND
        DO I=1,IM
	  IF(AVGCFRACH(I,J) < SPVAL)GRID1(I,J) = AVGCFRACH(I,J)*100.
        ENDDO
        ENDDO
        ID(1:25)=0
        ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(302),LVLS(1,IGET(302)),GRID1,IM,JM)
      ENDIF   
!     
!     TOTAL CLOUD FRACTION (INSTANTANEOUS).
      IF ((IGET(161).GT.0) .OR. (IGET(260).GT.0)) THEN
         IF(MODELNAME .EQ. 'GFS')THEN
          EGRID1=SPVAL
         ELSE IF(MODELNAME .EQ. 'NCAR' .OR. MODELNAME == 'RAPR')THEN
          DO J=JSTA,JEND
          DO I=1,IM
           frac = 1.
            do L=1,LM
             LL=LM-L+1
              frac=frac*(1.-CFR(I,J,LL))
            enddo
!      if(frac.le. 1.) print *,'frac,i,j',frac,i,j

           EGRID1(I,J) = 1. - frac
          ENDDO
          ENDDO

         ELSE IF (MODELNAME.EQ.'NMM'.OR.MODELNAME.EQ.'RSM')THEN
          DO J=JSTA,JEND
          DO I=1,IM
!           EGRID1(I,J)=AMAX1(CFRACL(I,J),
!     1                 AMAX1(CFRACM(I,J),CFRACH(I,J)))
            EGRID1(I,J)=1.-(1.-CFRACL(I,J))*(1.-CFRACM(I,J))*      &  
     &                 (1.-CFRACH(I,J))
          ENDDO
          ENDDO
         END IF
         DO J=JSTA,JEND
         DO I=1,IM
	    IF(ABS(EGRID1(I,J)-SPVAL).GT.SMALL)THEN
             GRID1(I,J) = EGRID1(I,J)*100.
	     TCLD(I,J)  = EGRID1(I,J)*100.         !B ZHOU, PASSED to CALCEILING
	    END IF 
         ENDDO
         ENDDO
         IF (IGET(161).GT.0) THEN
            ID(1:25)=0
            CALL GRIBIT(IGET(161),LVLS(1,IGET(161)),GRID1,IM,JM)
         ENDIF
      ENDIF
!
!     TIME AVERAGED TOTAL CLOUD FRACTION.
         IF (IGET(144).GT.0) THEN
           IF(MODELNAME == 'GFS')THEN
	    DO J=JSTA,JEND
            DO I=1,IM
	     IF(ABS(AVGTCDC(I,J)-SPVAL).GT.SMALL)                     &
                    GRID1(I,J) = AVGTCDC(I,J)*100.
            END DO
	    END DO 
	    
	   ELSE IF(MODELNAME == 'NMM')THEN
            DO J=JSTA,JEND
            DO I=1,IM
!               RSUM = NCFRST(I,J)+NCFRCV(I,J)
!               IF (RSUM.GT.0.0) THEN
!                  EGRID1(I,J)=(ACFRST(I,J)+ACFRCV(I,J))/RSUM
!               ELSE
!                  EGRID1(I,J) = D00
!               ENDIF
!ADDED BRAD'S MODIFICATION
               RSUM = D00
               IF (NCFRST(I,J) .GT. 0) RSUM=ACFRST(I,J)/NCFRST(I,J)
               IF (NCFRCV(I,J) .GT. 0)                               &
     &            RSUM=MAX(RSUM, ACFRCV(I,J)/NCFRCV(I,J))
               EGRID1(I,J) = RSUM
            ENDDO
            ENDDO
!
            DO J=JSTA,JEND
            DO I=1,IM
              IF(ABS(EGRID1(I,J)-SPVAL).GT.SMALL)                    &
     &              GRID1(I,J) = EGRID1(I,J)*100.
            ENDDO
            ENDDO
	   ELSE
	    GRID1=SPVAL 
	   END IF 
          IF(MODELNAME.EQ.'NMM' .OR. MODELNAME.EQ.'GFS')THEN
           ID(1:25)= 0
           ITCLOD     = INT(TCLOD)
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
           CALL GRIBIT(IGET(144),LVLS(1,IGET(144)),GRID1,IM,JM)
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
                  EGRID1(I,J) = ACFRST(I,J)/NCFRST(I,J)
               ELSE
                  EGRID1(I,J) = D00
               ENDIF
            ENDDO
            ENDDO
!
            DO J=JSTA,JEND
            DO I=1,IM
              IF(ABS(EGRID1(I,J)-SPVAL).GT.SMALL)                        &
     &             GRID1(I,J) = EGRID1(I,J)*100.
            ENDDO
            ENDDO
	   END IF 
          IF(MODELNAME.EQ.'NMM')THEN
           ID(1:25)=0
           ITCLOD     = INT(TCLOD)
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
           CALL GRIBIT(IGET(139),LVLS(1,IGET(139)),GRID1,IM,JM)
         ENDIF
!    
!     TIME AVERAGED CONVECTIVE CLOUD FRACTION.
         IF (IGET(143).GT.0) THEN
           IF(MODELNAME /= 'NMM')THEN
	    EGRID1=SPVAL
	   ELSE  
            DO J=JSTA,JEND
            DO I=1,IM
               IF (NCFRCV(I,J).GT.0.0) THEN
                  EGRID1(I,J) = ACFRCV(I,J)/NCFRCV(I,J)
               ELSE
                  EGRID1(I,J) = D00
               ENDIF
            ENDDO
            ENDDO
!
            DO J=JSTA,JEND
            DO I=1,IM
               IF(ABS(EGRID1(I,J)-SPVAL).GT.SMALL)                &  
     &               GRID1(I,J) = EGRID1(I,J)*100.
            ENDDO
            ENDDO
	   END IF
           IF(MODELNAME.EQ.'NMM')THEN 
            ID(1:25)=0
            ITCLOD     = INT(TCLOD)
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
           CALL GRIBIT(IGET(143),LVLS(1,IGET(143)),GRID1,IM,JM)
         ENDIF
!    
!     CLOUD BASE AND TOP FIELDS 
      IF((IGET(148).GT.0) .OR. (IGET(149).GT.0) .OR.              &
          (IGET(168).GT.0) .OR. (IGET(178).GT.0) .OR.             &
          (IGET(179).GT.0) .OR. (IGET(194).GT.0) .OR.             &
          (IGET(408).GT.0) .OR. (IGET(444).GT.0) .OR.             & 
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
            IBOTCu(I,J)=NINT(HBOT(I,J))
            IBOTDCu(I,J)=NINT(HBOTD(I,J))
            IBOTSCu(I,J)=NINT(HBOTS(I,J))
            ITOPCu(I,J)=NINT(HTOP(I,J))
            ITOPDCu(I,J)=NINT(HTOPD(I,J))
            ITOPSCu(I,J)=NINT(HTOPS(I,J))
            IF (IBOTCu(I,J)-ITOPCu(I,J) .LE. 1) THEN
              IBOTCu(I,J)=0
              ITOPCu(I,J)=100
            ENDIF
            IF (IBOTDCu(I,J)-ITOPDCu(I,J) .LE. 1) THEN
              IBOTDCu(I,J)=0
              ITOPDCu(I,J)=100
            ENDIF
            IF (IBOTSCu(I,J)-ITOPSCu(I,J) .LE. 1) THEN
              IBOTSCu(I,J)=0
              ITOPSCu(I,J)=100
            ENDIF
    !
    !--- Grid-scale cloud base & cloud top levels 
    !
    !--- Grid-scale cloud occurs when the mixing ratio exceeds QCLDmin
    !
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
    !
    !--- Combined (convective & grid-scale) cloud base & cloud top levels 
            IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME == 'RAPR')THEN
	     IBOTT(I,J)=IBOTGr(I,J)
	     ITOPT(I,J)=ITOPGr(I,J)
	    ELSE
             IBOTT(I,J)=MAX(IBOTGr(I,J), IBOTCu(I,J))
!	     if(i==200 .and. j==139)print*,'Debug cloud base 1: ',&
!             IBOTGr(I,J),IBOTCu(I,J),ibott(i,j)
             ITOPT(I,J)=MIN(ITOPGr(I,J), ITOPCu(I,J))
	    END IF 
          ENDDO      !--- End I loop
        ENDDO        !--- End J loop
      ENDIF          !--- End IF tests 
!
!-------------------------------------------------
!-----------  VARIOUS CLOUD BASE FIELDS ----------
!-------------------------------------------------
!
!--- "TOTAL" CLOUD BASE FIELDS (convective + grid-scale;  Ferrier, Feb '02)
!
      IF ((IGET(148).GT.0) .OR. (IGET(178).GT.0)                         &
           .OR.(IGET(260).GT.0) ) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            IBOT=IBOTT(I,J)
            IF (IBOT .LE. 0) THEN
              CLDP(I,J) = -50000.
              CLDZ(I,J) = -5000.
            ELSE IF (IBOT .LE. NINT(LMH(I,J))) THEN
              CLDP(I,J) = PMID(I,J,IBOT)
!	      if(i==200 .and. j==139)print*,'Debug cloud base 2: ',&
!              ibot,PMID(I,J,IBOT)
              IF (IBOT .EQ. LM) THEN
                CLDZ(I,J) = ZINT(I,J,LM)
              ELSE
                CLDZ(I,J) = HTM(I,J,IBOT+1)*T(I,J,IBOT+1)                &  
                           *(Q(I,J,IBOT+1)*D608+H1)*ROG*                 &
                           (LOG(PINT(I,J,IBOT+1))-LOG(CLDP(I,J)))        &
                           +ZINT(I,J,IBOT+1)
              ENDIF     !--- End IF (IBOT .EQ. LM) ...
            ENDIF       !--- End IF (IBOT .LE. 0) ...
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
               CALL GRIBIT(IGET(148),LVLS(1,IGET(148)),GRID1,IM,JM)
         ENDIF 
!    CLOUD BOTTOM HEIGHT
         IF (IGET(178).GT.0) THEN
  !--- Parameter was set to 148 in operational code  (Ferrier, Feb '02)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = CLDZ(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(178),LVLS(1,IGET(178)),GRID1,IM,JM)
         ENDIF
      ENDIF

!    GSD CLOUD BOTTOM HEIGHT
         IF (IGET(408).GT.0 .OR. IGET(444).GT.0) THEN
!- imported from RUC post
!  -- constants for effect of snow on ceiling
!      Also found in calvis.f
        rhoice = 970.
        coeffp = 10.36
        exponfp = 0.7776
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
          CLDZ(I,J) = -5000.

          pcldbase = -50000.
          zcldbase = -5000.
          watericemax = -99999.
          do k=1,lm
            LL=LM-k+1
            watericetotal(k) = QQW(i,j,ll) + QQI(i,j,ll)
            watericemax = max(watericemax,watericetotal(k))
          end do
          if (watericemax.lt.cloud_def_p) go to 3701

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

!        Eliminate fog layers near surface in watericetotal array
         do 1778 k=2,3
! --- Do this only when at least 10 mb (1000 Pa) above surface
!          if (pabovesfc(k).gt.1000.) then
           if (watericetotal(k).lt.cloud_def_p) then
             if (watericetotal(1).gt.cloud_def_p) then
               nfog = nfog+1
               go to 3441
             end if
             go to 3789
3441         continue
             do k1=1,k-1
               if (watericetotal(k1).ge.cloud_def_p) then
!                print *,'Zero fog',i,j,k1,watericetotal(k1),
!    1               g3(i,j,k1,p_p)/100.
                 watericetotal(k1)=0.
               end if
             end do
           end if
           go to 3789
!          end if
1778     continue

3789     continue

!       At surface?
          if (watericetotal(1).gt.cloud_def_p) then
            zcldbase = zmid(i,j,lm)
            go to 3788
          end if
!       Aloft?
          do 371 k=2,lm
            k1 = k
            if (watericetotal(k).gt.cloud_def_p) go to 372
 371      continue
          go to 3701
 372      continue
! -- Use vertical interpolation to obtain cloud level
        zcldbase = zmid(i,j,lm-k1+1) + (cloud_def_p-watericetotal(k1))    &
                 * (zmid(i,j,lm-k1+2)-zmid(i,j,lm-k1+1))                  &
                 / (watericetotal(k1-1) - watericetotal(k1))
        pcldbase = pmid(i,j,lm-k1+1) + (cloud_def_p-watericetotal(k1))    &
                 * (pmid(i,j,lm-k1+2)-pmid(i,j,lm-k1+1))                  &
                 / (watericetotal(k1-1) - watericetotal(k1))

! -- If within 4 levels of surface, just use lowest cloud level
!     as ceiling WITHOUT vertical interpolation.

          if (k1.le.4) then
           zcldbase = zmid(i,j,lm-k1+1)
           pcldbase = pmid(i,j,lm-k1+1)
          end if

 3788   continue

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
              do 3741 k=2,LM
              k1 = k
                if (ZMID(i,j,lm-k+1) .gt. zcldbase) go to 3742
 3741         continue
              go to 3743
 3742         continue
           pcldbase = pmid(i,j,lm-k1+2) + (zcldbase-ZMID(i,j,lm-k1+2))   &
               *(pmid(i,j,lm-k1+1)-pmid(i,j,lm-k1+2) )                   &
               /(zmid(i,j,lm-k1+1)-zmid(i,j,lm-k1+2) )
            end if
          end if
 3743     continue

 3701  continue

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
          RHB(I,J,k) = 100.*AMIN1(1.,E/ES)
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
!     1   I,J,k1,zmid(i,j,lm-k1+1),zmid(i,j,lm-k1),PBLH(I,J),RHB(i,j,k1)

         do k2=3,20
           if (zpbltop.lt.ZMID(i,j,LM-k2+1)) go to 744
         end do
         go to 745     ! No extra considerations for PBL-top cloud

  744    continue
!       print*,'check RH at PBL top, RH,i,j,k2',RHB(i,j,k2-1),i,j,k2-1
         if (rhb(i,j,k2-1).gt.95. ) then
           zcldbase = ZMID(i,j,LM-k2+2)
!       print*,' PBL cloud ceiling',zcldbase,i,j
           if (CLDZ(i,j).lt.-100.) then
!       print*,'add PBL cloud ceiling',zcldbase,i,j,k2
!     1         ,RHB(i,j,k2-1)
             npblcld = npblcld+1
             CLDZ(i,j) = zcldbase 
             CLDP(I,J) = PMID(i,j,LM-k2+2)
             go to 745 
           end if
           if ( zcldbase.lt.CLDZ(I,J)) then
!       print*,' change to PBL cloud ceiling',zcldbase,CLDZ(I,J),i,j,k2
!     1         ,RHB(i,j,k2-1)
!cc             npblcld = npblcld+1
             CLDZ(I,J) = zcldbase
           end if
         end if
  745    continue

!- include convective clouds
           IBOT=IBOTCu(I,J)
       if(IBOT.gt.0) then
!        print *,'IBOTCu(i,j)',i,j,IBOTCu(i,j)
         if(CLDZ(I,J).lt.-100.) then
!        print *,'add convective cloud, IBOT,CLDZ(I,J),ZMID(I,J,IBOT)'
!     1        ,IBOT,CLDZ(I,J),ZMID(I,J,IBOT),i,j
            CLDZ(I,J)=ZMID(I,J,IBOT)
            GOTO 746
         else if(ZMID(I,J,IBOT).lt.CLDZ(I,J)) then
!        print *,'change ceiling for convective cloud, CLDZ(I,J),
!     1              ZMID(I,J,IBOT),IBOT,i,j'
!     1        ,IBOT,CLDZ(I,J),ZMID(I,J,IBOT),IBOT,i,j
            CLDZ(I,J)=ZMID(I,J,IBOT)
         endif 
       endif

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
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = CLDZ(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(408),LVLS(1,IGET(408)),GRID1,IM,JM)
          ENDIF
!   GSD CLOUD BOTTOM PRESSURE
          IF (IGET(444).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = CLDP(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(444),LVLS(1,IGET(444)),GRID1,IM,JM) 
          ENDIF
      ENDIF   !End of GSD algorithm

!    B. ZHOU: CEILING
        IF (IGET(260).GT.0) THEN
            CALL CALCEILING(CLDZ,TCLD,CEILING)
            DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = CEILING(I,J)
             ENDDO
            ENDDO
            ID(1:25)=0
            CALL GRIBIT(IGET(260),LVLS(1,IGET(260)),GRID1,IM,JM)
         ENDIF

!    B. ZHOU: FLIGHT CONDITION RESTRICTION
        IF (IGET(261).GT.0) THEN
            CALL CALFLTCND(CEILING,FLTCND)
            DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = FLTCND(I,J)
             ENDDO
            ENDDO
            ID(1:25)=0
            CALL GRIBIT(IGET(261),LVLS(1,IGET(261)),GRID1,IM,JM)
         ENDIF
!
!---  Convective cloud base pressures (deep & shallow; Ferrier, Feb '02)
!
      IF (IGET(188) .GT. 0) THEN
       IF(MODELNAME .EQ. 'GFS')THEN
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
       ID(1:25)=0
       CALL GRIBIT(IGET(188),LVLS(1,IGET(188)),GRID1,IM,JM)
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
        ID(1:25)=0
        CALL GRIBIT(IGET(192),LVLS(1,IGET(192)),GRID1,IM,JM)
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
        ID(1:25)=0
        CALL GRIBIT(IGET(190),LVLS(1,IGET(190)),GRID1,IM,JM)
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
        ID(1:25)=0
        CALL GRIBIT(IGET(194),LVLS(1,IGET(194)),GRID1,IM,JM)
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
	ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(303),LVLS(1,IGET(303)),GRID1,IM,JM)
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
	ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(306),LVLS(1,IGET(306)),GRID1,IM,JM)
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
	ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(309),LVLS(1,IGET(309)),GRID1,IM,JM)
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
            IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
              CLDP(I,J) = PMID(I,J,ITOP)
              CLDT(I,J) = T(I,J,ITOP)
              IF (ITOP .EQ. LM) THEN
                CLDZ(I,J) = ZINT(I,J,LM)
              ELSE
                CLDZ(I,J) = HTM(I,J,ITOP+1)*T(I,J,ITOP+1)               &
                          *(Q(I,J,ITOP+1)*D608+H1)*ROG*                 &
                           (LOG(PINT(I,J,ITOP+1))-LOG(CLDP(I,J)))       &
                          +ZINT(I,J,ITOP+1)
              ENDIF    !--- End IF (ITOP .EQ. LM) ...
            ELSE
              CLDP(I,J) = -50000.
              CLDZ(I,J) = -5000.
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
               ID(1:25)=0
               CALL GRIBIT(IGET(149),LVLS(1,IGET(149)),GRID1,IM,JM)
         ENDIF
!   CLOUD TOP HEIGHT
!
          IF (IGET(179).GT.0) THEN
              DO J=JSTA,JEND
              DO I=1,IM
                 GRID1(I,J) = CLDZ(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(179),LVLS(1,IGET(179)),GRID1,IM,JM)
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
          do k=1,lm
            LL=LM-k+1
            watericetotal(k) = QQW(i,j,ll) + QQI(i,j,ll)
          enddo

          if (watericetotal(LM).gt.cloud_def_p) then
            zcldtop = zmid(i,j,1)
            go to 3799
          end if
! in RUC          do 373 k=LM,2,-1
          do 373 k=LM-1,2,-1
            if (watericetotal(k).gt.cloud_def_p) go to 374
 373      continue
          go to 3799
 374      zcldtop = zmid(i,j,lm-k+1) + (cloud_def_p-watericetotal(k))   &
                 * (zmid(i,j,lm-k+2)-zmid(i,j,lm-k+1))                &
                 / (watericetotal(k+1) - watericetotal(k))
 3799     continue

            ITOP=ITOPT(I,J)
            IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
              CLDP(I,J) = PMID(I,J,ITOP)
              CLDT(I,J) = T(I,J,ITOP)
            ELSE
              CLDP(I,J) = -50000.
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

! check consistency of cloud bas and cloud top
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
               CALL GRIBIT(IGET(406),LVLS(1,IGET(406)),GRID1,IM,JM)
         ENDIF
!   GSD CLOUD TOP HEIGHT
!
          IF (IGET(409).GT.0) THEN
              DO J=JSTA,JEND
              DO I=1,IM
                 GRID1(I,J) = CLDZ(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(409),LVLS(1,IGET(409)),GRID1,IM,JM)
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
               ID(1:25)=0
               CALL GRIBIT(IGET(168),LVLS(1,IGET(168)),GRID1,IM,JM)
         ENDIF
!
!huang  CLOUD TOP BRIGHTNESS TEMPERATURE
          IF (IGET(275).GT.0) THEN
             num_thick=0   ! for debug
             DO J=JSTA,JEND
             DO I=1,IM
               opdepth=0.
               llmh=nint(lmh(i,j))
!bsf - start
!-- Add subgrid-scale convective clouds for WRF runs
               do k=1,llmh
                 CU_ir(k)=0.
               enddo
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
      print *,'num_points, num_thick = ',(jend-jsta+1)*im,num_thick
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
            ID(1:25)=0
!	    ID(02)=129    ! Parameter Table 129
            CALL GRIBIT(IGET(275),LVLS(1,IGET(275)),GRID1,IM,JM)
         ENDIF

!
!---  Convective cloud top pressures (deep & shallow; Ferrier, Feb '02)
!
      IF (IGET(189) .GT. 0) THEN
       IF(MODELNAME .EQ. 'GFS')THEN
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
        ID(1:25)=0
        CALL GRIBIT(IGET(189),LVLS(1,IGET(189)),GRID1,IM,JM)
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
        ID(1:25)=0
        CALL GRIBIT(IGET(193),LVLS(1,IGET(193)),GRID1,IM,JM) 
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
        ID(1:25)=0
        CALL GRIBIT(IGET(191),LVLS(1,IGET(191)),GRID1,IM,JM)
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
        ID(1:25)=0
        CALL GRIBIT(IGET(195),LVLS(1,IGET(195)),GRID1,IM,JM)
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
	ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(304),LVLS(1,IGET(304)),GRID1,IM,JM)
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
	ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(307),LVLS(1,IGET(307)),GRID1,IM,JM)
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
	ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(310),LVLS(1,IGET(310)),GRID1,IM,JM)
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
	ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(305),LVLS(1,IGET(305)),GRID1,IM,JM)
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
	ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(308),LVLS(1,IGET(308)),GRID1,IM,JM)
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
	ITCLOD     = INT(TCLOD)
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
        CALL GRIBIT(IGET(311),LVLS(1,IGET(311)),GRID1,IM,JM)
       ENDIF
!
!--- Convective cloud fractions from modified Slingo (1987)
!
      IF (IGET(196) .GT. 0) THEN
          GRID1=SPVAL
          DO J=JSTA,JEND
          DO I=1,IM
            if(CNVCFR(I,J)/=SPVAL)GRID1(I,J)=100.*CNVCFR(I,J)   !-- convert to percent
          ENDDO
          ENDDO
          ID(1:25)=0
           CALL GRIBIT(IGET(196),LVLS(1,IGET(196)),GRID1,IM,JM)
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
	  ITCLOD     = INT(TCLOD)
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
          CALL GRIBIT(IGET(342),LVLS(1,IGET(342)),GRID1,IM,JM)
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
	  ITCLOD     = INT(TCLOD)
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
          CALL GRIBIT(IGET(313),LVLS(1,IGET(313)),GRID1,IM,JM)
      END IF      
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
            ITRDSW     = INT(TRDSW)
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
          CALL GRIBIT(IGET(126),LVLS(1,IGET(126)),GRID1,IM,JM)
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
            ITRDSW     = INT(TRDSW)
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
          CALL GRIBIT(IGET(298),LVLS(1,IGET(298)),GRID1,IM,JM)
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
            ITRDSW     = INT(TRDSW)
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
          CALL GRIBIT(IGET(297),LVLS(1,IGET(297)),GRID1,IM,JM)
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
            ITRDLW     = INT(TRDLW)
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
          CALL GRIBIT(IGET(127),LVLS(1,IGET(127)),GRID1,IM,JM)
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
            ITRDSW     = INT(TRDSW)
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
          CALL GRIBIT(IGET(128),LVLS(1,IGET(128)),GRID1,IM,JM)
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
            ITRDLW     = INT(TRDLW)
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
          CALL GRIBIT(IGET(129),LVLS(1,IGET(129)),GRID1,IM,JM)
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
            ITRDSW     = INT(TRDSW)
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
          CALL GRIBIT(IGET(130),LVLS(1,IGET(130)),GRID1,IM,JM)
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
            ITRDLW     = INT(TRDLW)
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
          CALL GRIBIT(IGET(131),LVLS(1,IGET(131)),GRID1,IM,JM)
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
          CALL GRIBIT(IGET(274),LVLS(1,IGET(274)),GRID1,IM,JM)
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
	  ID(1:25)=0
	  ID(02)=129    ! Parameter Table 129
          CALL GRIBIT(IGET(265),LVLS(1,IGET(265)),GRID1,IM,JM)
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
         ID(1:25)=0
         CALL GRIBIT(IGET(156),LVLS(1,IGET(156)),GRID1,IM,JM)
      ENDIF
!     
!     CURRENT INCOMING LW RADIATION AT THE SURFACE.
      IF (IGET(157).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
          IF(MODELNAME.eq.'RSM') THEN      !add by Binbin: RSM has direct RLWIN output
           GRID1(I,J)=RLWIN(I,J)
          ELSE
           IF(SIGT4(I,J).GT.0.0) THEN
             LLMH=NINT(LMH(I,J))
             TLMH=T(I,J,LLMH)
             FACTRL=5.67E-8*TLMH*TLMH*TLMH*TLMH/SIGT4(I,J)
           ELSE
             FACTRL=0.0
           ENDIF
           GRID1(I,J)=RLWIN(I,J)*FACTRL
          ENDIF
         ENDDO
         ENDDO
!
         ID(1:25)=0
         CALL GRIBIT(IGET(157),LVLS(1,IGET(157)),GRID1,IM,JM)
      ENDIF
!     
!     CURRENT OUTGOING SW RADIATION AT THE SURFACE.
      IF (IGET(141).GT.0) THEN
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
         ID(1:25)=0
         CALL GRIBIT(IGET(141),LVLS(1,IGET(141)),GRID1,IM,JM)
      ENDIF
!     
!     CURRENT OUTGOING LW RADIATION AT THE SURFACE.
      IF (IGET(142).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = RADOT(I,J)
               ENDDO
               ENDDO
         ID(1:25)=0
         CALL GRIBIT(IGET(142),LVLS(1,IGET(142)),GRID1,IM,JM)
      ENDIF
!     
!     CURRENT (instantaneous) INCOMING CLEARSKY SW RADIATION AT THE SURFACE.
      IF (IGET(262).GT.0) THEN
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
         ID(1:25)=0
         CALL GRIBIT(IGET(262),LVLS(1,IGET(262)),GRID1,IM,JM)
      ENDIF
!     
!     TIME AVERAGED INCOMING CLEARSKY SW RADIATION AT THE SURFACE.
      IF (IGET(383).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = ASWINC(I,J)
         ENDDO
	 ENDDO
	 ID(1:25)=0
         ITRDSW     = INT(TRDSW)
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
         CALL GRIBIT(IGET(383),LVLS(1,IGET(383)),GRID1,IM,JM)
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
         ITRDSW     = INT(TRDSW)
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
         CALL GRIBIT(IGET(386),LVLS(1,IGET(386)),GRID1,IM,JM)
      ENDIF
!     
!     TIME AVERAGED OUTGOING CLEARSKY SW RADIATION AT THE MODEL TOP
      IF (IGET(387).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = ASWTOAC(I,J)
         ENDDO
	 ENDDO
	 ID(1:25)=0
         ITRDSW     = INT(TRDSW)
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
         CALL GRIBIT(IGET(387),LVLS(1,IGET(387)),GRID1,IM,JM)
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
         ITRDSW     = INT(TRDSW)
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
         CALL GRIBIT(IGET(388),LVLS(1,IGET(388)),GRID1,IM,JM)
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
         ITRDLW     = INT(TRDLW)
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
         CALL GRIBIT(IGET(382),LVLS(1,IGET(382)),GRID1,IM,JM)
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
         ITRDLW     = INT(TRDLW)
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
         CALL GRIBIT(IGET(384),LVLS(1,IGET(384)),GRID1,IM,JM)
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
         ITRDLW     = INT(TRDLW)
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
         CALL GRIBIT(IGET(385),LVLS(1,IGET(385)),GRID1,IM,JM)
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
         ITRDSW     = INT(TRDSW)
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
         CALL GRIBIT(IGET(401),LVLS(1,IGET(401)),GRID1,IM,JM)
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
         ITRDSW     = INT(TRDSW)
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
         CALL GRIBIT(IGET(402),LVLS(1,IGET(402)),GRID1,IM,JM)
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
         ITRDSW     = INT(TRDSW)
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
         CALL GRIBIT(IGET(403),LVLS(1,IGET(403)),GRID1,IM,JM)
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
         ITRDSW     = INT(TRDSW)
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
         CALL GRIBIT(IGET(404),LVLS(1,IGET(404)),GRID1,IM,JM)
      ENDIF
!
!     END OF ROUTINE.
!
      RETURN
      END
