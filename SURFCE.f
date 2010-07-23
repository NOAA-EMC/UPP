      SUBROUTINE SURFCE
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    SURFCE      POST SURFACE BASED FIELDS
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-21       
!     
! ABSTRACT:
!     THIS ROUTINE POSTS SURFACE BASED FIELDS.
!   .     
!     
! PROGRAM HISTORY LOG:
!   92-12-21  RUSS TREADON
!   94-08-04  MICHAEL BALDWIN - ADDED OUTPUT OF SFC FLUXES OF
!                               SENS AND LATENT HEAT AND THETA AT Z0
!   94-11-04  MICHAEL BALDWIN - ADDED INSTANTANEOUS PRECIP TYPE
!   96-03-19  MICHAEL BALDWIN - CHANGE SOIL PARAMETERS
!   96-09-25  MICHAEL BALDWIN - ADDED SNOW RATIO FROM EXPLICIT SCHEME
!   96-10-17  MICHAEL BALDWIN - CHANGED SFCEVP,POTEVP TO ACCUM.  TOOK
!                               OUT -PTRACE FOR ACSNOW,SSROFF,BGROFF.
!   97-04-23  MICHAEL BALDWIN - TOOK OUT -PTRACE FOR ALL PRECIP FIELDS
!   98-06-12  T BLACK         - CONVERSION FROM 1-D TO 2-D
!   98-07-17  MIKE BALDWIN - REMOVED LABL84
!   98-08-18  MIKE BALDWIN - COMPUTE RH OVER ICE
!   98-12-22  MIKE BALDWIN - BACK OUT RH OVER ICE
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   01-10-22  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-06-11  MIKE BALDWIN - WRF VERSION ASSUMING ALL ACCUM VARS
!                            HAVE BUCKETS THAT FILL FROM T=00H ON
!   02-08-28  H CHUANG - COMPUTE FIELDS AT SHELTER LEVELS FOR WRF
!   04-12-09  H CHUANG - ADD ADDITIONAL LSM FIELDS
!   05-07-07  BINBIN ZHOU - ADD RSM MODEL
!   05-08-24  GEOFF MANIKIN - ADDED DOMINANT PRECIP TYPE
!     
! USAGE:    CALL SURFCE
!   INPUT ARGUMENT LIST:
!
!   OUTPUT ARGUMENT LIST: 
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       BOUND    - ENFORCE LOWER AND UPPER LIMITS ON ARRAY ELEMENTS.
!       DEWPOINT - COMPUTE DEWPOINT TEMPERATURE.
!       CALDRG   - COMPUTE SURFACE LAYER DRAG COEFFICENT
!       CALTAU   - COMPUTE SURFACE LAYER U AND V WIND STRESSES.
!
!     LIBRARY:
!       COMMON   - CTLBLK
!                  RQSTFLD
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
!
!     
!     INCLUDE GRID DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
!
      use vrbls3d   
      use vrbls2d   
      use soil
      use masks
      use params_mod
      use ctlblk_mod
      use rqstfld_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      INCLUDE "mpif.h"
!     
!     IN NGM SUBROUTINE OUTPUT WE FIND THE FOLLOWING COMMENT.
!     "IF THE FOLLOWING THRESHOLD VALUES ARE CHANGED, CONTACT
!     TDL/SYNOPTIC-SCALE TECHNIQUES BRANCH (PAUL DALLAVALLE
!     AND JOHN JENSENIUS).  THEY MAY BE USING IT IN ONE OF 
!     THEIR PACKING CODES."  THE THRESHOLD VALUE IS 0.01 INCH
!     OR 2.54E-4 METER.  PRECIPITATION VALUES LESS THAN THIS
!     THRESHOLD ARE SET TO MINUS ONE TIMES THIS THRESHOLD.
      real,PARAMETER :: PTRACE = 0.000254E0
!     
!     SET CELCIUS TO KELVIN AND SECOND TO HOUR CONVERSION.
      real,PARAMETER :: C2K    = 273.15
      integer,PARAMETER :: NALG    = 5
      real,PARAMETER :: SEC2HR = 1./3600.
      integer,parameter :: nosoiltype=9
!     
!     DECLARE VARIABLES.
!     
      INTEGER IWX1(IM,JM),NROOTS(IM,JM),IWX4(IM,JM),IWX5(IM,JM)
      REAL IWX2(IM,JM),IWX3(IM,JM)
      REAL PSFC(IM,JM),TSFC(IM,JM),QSFC(IM,JM),RHSFC(IM,JM)
      REAL ZSFC(IM,JM),THSFC(IM,JM),DWPSFC(IM,JM),EVP(IM,JM)
!      REAL ANCPRC(IM,JM),P1D(IM,JM),T1D(IM,JM),Q1D(IM,JM)
      REAL P1D(IM,JM),T1D(IM,JM),Q1D(IM,JM),ZWET(IM,JM)
      REAL EGRID1(IM,JM),EGRID2(IM,JM),UA(IM,JM),VA(IM,JM)
      REAL GRID1(IM,JM),GRID2(IM,JM),IW(IM,JM),IWM1
!      REAL SLEET(IM,JM),RAIN(IM,JM),FREEZR(IM,JM),SNOW(IM,JM)
      REAL SLEET(IM,JM,NALG),RAIN(IM,JM,NALG),           &
           FREEZR(IM,JM,NALG),SNOW(IM,JM,NALG)
      REAL SLEET1(IM,JM),RAIN1(IM,JM),FREEZR1(IM,JM)     &
          ,SNOW1(IM,JM)
      REAL SLEET2(IM,JM),RAIN2(IM,JM),FREEZR2(IM,JM)     &
          ,SNOW2(IM,JM)
      REAL SLEET3(IM,JM),RAIN3(IM,JM),FREEZR3(IM,JM)     &
          ,SNOW3(IM,JM)
      REAL SLEET4(IM,JM),RAIN4(IM,JM),FREEZR4(IM,JM)     &
          ,SNOW4(IM,JM)
      REAL SLEET5(IM,JM),RAIN5(IM,JM),FREEZR5(IM,JM)     &
          ,SNOW5(IM,JM)
      REAL DOMS(IM,JM),DOMR(IM,JM),DOMIP(IM,JM),DOMZR(IM,JM)

!GSD
      REAL totprcp, snowratio,t2,rainl

      REAL ECAN(IM,JM),EDIR(IM,JM),ETRANS(IM,JM),ESNOW(IM,JM) &
           ,SMCDRY(IM,JM),SMCMAX(IM,JM)
      REAL RSMIN(IM,JM),SMCREF(IM,JM)     &
       ,RCS(IM,JM),RCQ(IM,JM),RCT(IM,JM),RCSOIL(IM,JM)  &
       ,GC(IM,JM)     
      REAL REFSMC(nosoiltype),WLTSMC(nosoiltype)
      DATA REFSMC /0.2484580,0.3678367,0.3981426  &
                  ,0.2820649,0.3204950,0.3606471  &
                  ,0.2924667,0.3012990,0.2484580/
      DATA WLTSMC /2.8506856E-02,0.1196190,0.1385488  &
                  ,4.6867151E-02,9.9965721E-02,0.1030795  &
                  , 6.9077924E-02,6.5861143E-02  &
                  ,2.8506856E-02/ 
!
      integer I,J,IWX,ISNO,ITMAXMIN,IFINCR,ISVALUE,II,JJ,                &
              ITPREC,IIP,IZR,IRAIN,ITSRFC,L,LS,IVEG,LLMH,                &
              IVG,IRTN,ISEED
      real RDTPHS,TLOW,TSFCK,QSAT,DTOP,DBOT,SNEQV,RRNUM,SFCPRS,SFCQ,     &
           RC,SFCTMP,SNCOVR,FACTRS,SOLAR
!     
!****************************************************************************
!
!     START SURFCE.
!
!     
!***  BLOCK 1.  SURFACE BASED FIELDS.
!
!     IF ANY OF THE FOLLOWING "SURFACE" FIELDS ARE REQUESTED,
!     WE NEED TO COMPUTE THE FIELDS FIRST.
!     
      IF ( (IGET(024).GT.0).OR.(IGET(025).GT.0).OR.     &
           (IGET(026).GT.0).OR.(IGET(027).GT.0).OR.     &
           (IGET(028).GT.0).OR.(IGET(029).GT.0).OR.     &
           (IGET(154).GT.0).OR.                         &
           (IGET(034).GT.0).OR.(IGET(076).GT.0) ) THEN
!     
         DO J=JSTA,JEND
         DO I=1,IM
!
!           SCALE ARRAY FIS BY GI TO GET SURFACE HEIGHT.
!            ZSFC(I,J)=FIS(I,J)*GI
            ZSFC(I,J)=ZINT(I,J,LM+1)
!
!           SURFACE PRESSURE.
            PSFC(I,J)=PINT(I,J,NINT(LMH(I,J))+1)
!     
!           SURFACE (SKIN) POTENTIAL TEMPERATURE AND TEMPERATURE.
            THSFC(I,J)=THS(I,J)
            TSFC(I,J) =spval
            IF(THSFC(i,j) /= spval)  &
            TSFC(I,J) =THSFC(I,J)*(PSFC(I,J)/P1000)**CAPA 
!     
!           SURFACE SPECIFIC HUMIDITY, RELATIVE HUMIDITY,
!           AND DEWPOINT.  ADJUST SPECIFIC HUMIDITY IF
!           RELATIVE HUMIDITY EXCEEDS 0.1 OR 1.0.
!
            QSFC(I,J)=QS(I,J)
            QSFC(I,J)=AMAX1(H1M12,QSFC(I,J))
            TSFCK    =TSFC(I,J)
!     
            QSAT=PQ0/PSFC(I,J)*EXP(A2*(TSFCK-A3)/(TSFCK-A4))
            RHSFC(I,J)=QSFC(I,J)/QSAT

            IF (RHSFC(I,J).GT.H1 ) RHSFC(I,J) = H1
            IF (RHSFC(I,J).LT.D00) RHSFC(I,J) = D01
            QSFC(I,J)  = RHSFC(I,J)*QSAT
            EVP(I,J)   = PSFC(I,J)*QSFC(I,J)/(EPS+ONEPS*QSFC(I,J))
            EVP(I,J)   = EVP(I,J)*D001
!     
!mp           ACCUMULATED NON-CONVECTIVE PRECIP.
!mp            IF(IGET(034).GT.0)THEN
!mp              IF(LVLS(1,IGET(034)).GT.0)THEN

!           ACCUMULATED PRECIP (convective + non-convective)
            IF(IGET(087).GT.0)THEN
              IF(LVLS(1,IGET(087)).GT.0)THEN
!	write(6,*) 'acprec, ancprc, cuprec: ', ANCPRC(I,J)+CUPREC(I,J),
!     +		ANCPRC(I,J),CUPREC(I,J)
!                 ACPREC(I,J)=ANCPRC(I,J)+CUPREC(I,J)
              ENDIF
            ENDIF

         ENDDO
         ENDDO
!     
!        INTERPOLATE/OUTPUT REQUESTED SURFACE FIELDS.
!     
!        SURFACE PRESSURE.
         IF (IGET(024).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=PSFC(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(024),LVLS(1,IGET(024)), GRID1,IM,JM)
         ENDIF
!     
!        SURFACE HEIGHT.
         IF (IGET(025).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=ZSFC(I,J)
            ENDDO
            ENDDO
!            CALL BOUND(GRID1,D00,H99999)
            ID(1:25) = 0
            CALL GRIBIT(IGET(025),LVLS(1,IGET(025)), GRID1,IM,JM)
         ENDIF
!     
!        SURFACE (SKIN) TEMPERATURE.
         IF (IGET(026).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=TSFC(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(026),LVLS(1,IGET(026)), GRID1,IM,JM)
         ENDIF
!     
!        SURFACE (SKIN) POTENTIAL TEMPERATURE.
         IF (IGET(027).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=THSFC(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(027),LVLS(1,IGET(027)), GRID1,IM,JM)
         ENDIF
!     
!        SURFACE SPECIFIC HUMIDITY.
         IF (IGET(028).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=QSFC(I,J)
            ENDDO
            ENDDO
            CALL BOUND(GRID1,H1M12,H99999)
            ID(1:25) = 0
            CALL GRIBIT(IGET(028),LVLS(1,IGET(028)), GRID1,IM,JM)
         ENDIF
!     
!        SURFACE DEWPOINT TEMPERATURE.
         IF (IGET(029).GT.0) THEN
            CALL DEWPOINT(EVP,DWPSFC)
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=DWPSFC(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(029),LVLS(1,IGET(029)), GRID1,IM,JM)
         ENDIF
!     
!        SURFACE RELATIVE HUMIDITY.
         IF (IGET(076).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=RHSFC(I,J)*100.
            ENDDO
            ENDDO
            CALL BOUND(GRID1,H1,H100)
            ID(1:25) = 0
            CALL GRIBIT(IGET(076),LVLS(1,IGET(076)), GRID1,IM,JM)
         ENDIF
!     
      ENDIF
!
!     ADDITIONAL SURFACE-SOIL LEVEL FIELDS.
!

      DO L=1,NSOIL
!     SOIL TEMPERATURE.
      IF (IGET(116).GT.0) THEN
        IF (LVLS(L,IGET(116)).GT.0) THEN
          IF(iSF_SURFACE_PHYSICS==3)THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=STC(I,J,L)
            ENDDO
            ENDDO
            ID(1:25)=0
            ID(9)=111
            IF(L==1)ID(9)=1
            print*,'SLLEVEL in SURFCE= ',SLLEVEL
            ID(11)=NINT(SLLEVEL(L)*100.)
            CALL GRIBIT(IGET(116),L,GRID1,IM,JM)
          ELSE
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=STC(I,J,L)
!            IF(SM(I,J)>=1.0 .OR. SICE(I,J)>SMALL)
!     +         GRID1(I,J)=SPVAL
            ENDDO
            ENDDO
            ID(1:25) = 0
            DTOP=0.
            DO LS=1,L-1
              DTOP=DTOP+SLDPTH(LS)
            ENDDO
            DBOT=DTOP+SLDPTH(L)
            ID(10) = NINT(DTOP*100.)
            ID(11) = NINT(DBOT*100.)
            CALL GRIBIT(IGET(116),L,GRID1,IM,JM)
          END IF
        ENDIF
      ENDIF
!
!     SOIL MOISTURE.
      IF (IGET(117).GT.0) THEN
        IF (LVLS(L,IGET(117)).GT.0) THEN
          IF(iSF_SURFACE_PHYSICS==3)THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SMC(I,J,L)
            ENDDO
            ENDDO
            ID(1:25)=0
            ID(9)=111
            IF(L==1)ID(9)=1
            ID(11)=NINT(SLLEVEL(L)*100.)
            CALL GRIBIT(IGET(117),L,GRID1,IM,JM)
          ELSE
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SMC(I,J,L)
!            IF(SM(I,J)>=1.0 .OR. SICE(I,J)>SMALL)
!     +         GRID1(I,J)=SPVAL
             IF(GRID1(I,J)==spval .and.  &
             (sm(i,j)<1.0.and.sice(i,j)<1.0e-10))  &
             print*,'bad smc land point ',i,j
            ENDDO
            ENDDO
            ID(1:25) = 0
            DTOP=0.
            DO LS=1,L-1
              DTOP=DTOP+SLDPTH(LS)
            ENDDO
            DBOT=DTOP+SLDPTH(L)
            ID(10) = NINT(DTOP*100.)
            ID(11) = NINT(DBOT*100.)
            CALL GRIBIT(IGET(117),L,GRID1,IM,JM)
          END IF
        ENDIF
      ENDIF
!      ENDDO
!     ADD LIQUID SOIL MOISTURE
      IF (IGET(225).GT.0) THEN
        IF (LVLS(L,IGET(225)).GT.0) THEN
          IF(iSF_SURFACE_PHYSICS==3)THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SH2O(I,J,L)
            ENDDO
            ENDDO
            ID(1:25)=0
            ID(02) = 130
            ID(9)=111
            IF(L==1)ID(9)=1
            ID(11)=NINT(SLLEVEL(L)*100.)
            CALL GRIBIT(IGET(225),L,GRID1,IM,JM)
          ELSE
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SH2O(I,J,L)
!            IF(SM(I,J)>=1.0 .OR. SICE(I,J)>SMALL)
!     +         GRID1(I,J)=SPVAL
            ENDDO
            ENDDO
            ID(1:25) = 0
            DTOP=0.
            DO LS=1,L-1
              DTOP=DTOP+SLDPTH(LS)
            ENDDO
            DBOT=DTOP+SLDPTH(L)
            ID(10) = NINT(DTOP*100.)
            ID(11) = NINT(DBOT*100.)
            ID(02) = 130
            CALL GRIBIT(IGET(225),L,GRID1,IM,JM)
          END IF
        ENDIF
      ENDIF
! END OF NSOIL LOOP
      ENDDO
!
!     BOTTOM SOIL TEMPERATURE.
      IF (IGET(115).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=TG(I,J) ! change to tg per discussion with LSM group
            ENDDO
            ENDDO
         ID(1:25) = 0
         ISVALUE     = 300
         ID(11) = ISVALUE
         CALL GRIBIT(IGET(115),LVLS(1,IGET(115)),GRID1,IM,JM)
      ENDIF
!
!     SOIL MOISTURE AVAILABILITY
      IF (IGET(171).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             IF(SMSTAV(I,J)/=SPVAL)GRID1(I,J)=SMSTAV(I,J)*100.
            ENDDO
            ENDDO
         ID(1:25) = 0
         ID(10) = 0
         ID(11) = 100
         CALL GRIBIT(IGET(171),LVLS(1,IGET(171)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL SOIL MOISTURE
      IF (IGET(036).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
	     IF(SMSTOT(I,J)/=SPVAL .AND. SM(I,J).GT.SMALL     &
      	       .AND. SICE(I,J).LT.SMALL)THEN
	      GRID1(I,J)=1000.0  ! TEMPORY FIX TO MAKE SURE SMSTOT=1 FOR WATER
	     ELSE  
              GRID1(I,J)=SMSTOT(I,J)
	     END IF 
            ENDDO
            ENDDO
         ID(1:25) = 0
         ID(10) = 0
         ID(11) = 200
         CALL GRIBIT(IGET(036),LVLS(1,IGET(036)),GRID1,IM,JM)
      ENDIF
!
!     PLANT CANOPY SURFACE WATER.
      IF ( IGET(118).GT.0 ) THEN
            GRID1=SPVAL
            DO J=JSTA,JEND
            DO I=1,IM
             IF(CMC(I,J)/=SPVAL)GRID1(I,J)=CMC(I,J)*1000.
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(118),LVLS(1,IGET(118)),GRID1,IM,JM)
      ENDIF
!
!     SNOW WATER EQUIVALENT.
      IF ( IGET(119).GT.0 ) THEN
            GRID1=SPVAL
            DO J=JSTA,JEND
            DO I=1,IM
!             GRID1(I,J)=SNO(I,J)*1000.
             GRID1(I,J)=SNO(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(119),LVLS(1,IGET(119)),GRID1,IM,JM)
      ENDIF
!
!     PERCENT SNOW COVER.
      IF ( IGET(120).GT.0 ) THEN
         GRID1=SPVAL
         DO J=JSTA,JEND
         DO I=1,IM
!             GRID1(I,J)=PCTSNO(I,J)
           SNEQV=SNO(I,J)
           IVEG=IVGTYP(I,J)
           IF(IVEG.EQ.0)IVEG=7
           CALL SNFRAC (SNEQV,IVEG,SNCOVR)
           GRID1(I,J)=SNCOVR*100.
         ENDDO
         ENDDO
         CALL BOUND(GRID1,D00,H100)
         ID(1:25) = 0
         CALL GRIBIT(IGET(120),LVLS(1,IGET(120)),GRID1,IM,JM)
      ENDIF
! ADD SNOW DEPTH
      IF ( IGET(224).GT.0 ) THEN
            ii=im/2
            jj=(jsta+jend)/2
	    GRID1=SPVAL
            DO J=JSTA,JEND
            DO I=1,IM
             IF(SI(I,J)/=SPVAL)GRID1(I,J)=SI(I,J)/1000.  ! SI comes out of WRF in mm
             IF(i==ii.and.j== jj)print*,'sample snow depth in GRIBIT= '  &
                ,si(i,j)             	     
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(224),LVLS(1,IGET(224)),GRID1,IM,JM)
      ENDIF      
! ADD POTENTIAL EVAPORATION
      IF ( IGET(242).GT.0 ) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=POTEVP(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(242),LVLS(1,IGET(242)),GRID1,IM,JM)
      ENDIF
! ADD ICE THICKNESS
      IF ( IGET(349).GT.0 ) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=DZICE(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(349),LVLS(1,IGET(349)),GRID1,IM,JM)
      ENDIF      

! ADD EC,EDIR,ETRANS,ESNOW,SMCDRY,SMCMAX
! ONLY OUTPUT NEW LSM FIELDS FOR NMM AND ARW BECAUSE RSM USES OLD SOIL TYPES
      IF (MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'NMM' .OR. MODELNAME.EQ.'RAPR')THEN
      IF ( IGET(228).GT.0 .OR. IGET(229).GT.0      &
       .OR.IGET(230).GT.0 .OR. IGET(231).GT.0      &
       .OR.IGET(232).GT.0 .OR. IGET(233).GT.0) THEN
        DO J=JSTA,JEND
         DO I=1,IM
! ----------------------------------------------------------------------
!          IF(QWBS(I,J).gt.0.001)print*,'NONZERO QWBS',i,j,QWBS(I,J)
!          IF(abs(SM(I,J)-0.).lt.1.0E-5)THEN
          IF( (abs(SM(I,J)-0.)   .lt. 1.0E-5) .AND.              &
     &        (abs(SICE(I,J)-0.) .lt. 1.0E-5) ) THEN
           CALL ETCALC(QWBS(I,J),POTEVP(I,J),SNO(I,J),VEGFRC(I,J) &
     &  ,  ISLTYP(I,J),SH2O(I,J,1:1),CMC(I,J)                     &
     &  ,  ECAN(I,J),EDIR(I,J),ETRANS(I,J),ESNOW(I,J),SMCDRY(I,J) &
     &  ,  SMCMAX(I,J) )
          ELSE
           ECAN(I,J)=0.
           EDIR(I,J)=0.
           ETRANS(I,J)=0.
           ESNOW(I,J)=0.
           SMCDRY(I,J)=0.
           SMCMAX(I,J)=0.
          END IF
         END DO
        END DO

        IF ( IGET(228).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=ECAN(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(228),LVLS(1,IGET(228)),GRID1,IM,JM)
        ENDIF	

        IF ( IGET(229).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=EDIR(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(229),LVLS(1,IGET(229)),GRID1,IM,JM)
        ENDIF

        IF ( IGET(230).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=ETRANS(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(230),LVLS(1,IGET(230)),GRID1,IM,JM)
        ENDIF
	

        IF ( IGET(231).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=ESNOW(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
	 ID(02)= 130
         CALL GRIBIT(IGET(231),LVLS(1,IGET(231)),GRID1,IM,JM)
        ENDIF	

        IF ( IGET(232).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SMCDRY(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
	 ID(02)= 130
         CALL GRIBIT(IGET(232),LVLS(1,IGET(232)),GRID1,IM,JM)
        ENDIF

        IF ( IGET(233).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SMCMAX(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
	 ID(02)= 130
         CALL GRIBIT(IGET(233),LVLS(1,IGET(233)),GRID1,IM,JM)
        ENDIF
	
      ENDIF
      END IF  ! endif for ncar and nmm options
!
!     
!
!***  BLOCK 2.  SHELTER (2M) LEVEL FIELDS.
!     
!     COMPUTE/POST SHELTER LEVEL FIELDS.
!     
      IF ( (IGET(106).GT.0).OR.(IGET(112).GT.0).OR.     &
           (IGET(113).GT.0).OR.(IGET(114).GT.0).OR.     &
           (IGET(138).GT.0) ) THEN
!
!HC  COMPUTE SHELTER PRESSURE BECAUSE IT WAS NOT OUTPUT FROM WRF       
        IF(MODELNAME .EQ. 'NCAR' .OR. MODELNAME.EQ.'RSM'.OR. MODELNAME.EQ.'RAPR')THEN
         DO J=JSTA,JEND
         DO I=1,IM
          TLOW=T(I,J,NINT(LMH(I,J)))
          PSHLTR(I,J)=PSFC(I,J)*EXP(-0.068283/TLOW)
         END DO
         END DO 
	END IF 
!
!        SHELTER LEVEL TEMPERATURE
         IF (IGET(106).GT.0) THEN
	    GRID1=spval
            DO J=JSTA,JEND
            DO I=1,IM
!             GRID1(I,J)=TSHLTR(I,J)
!HC CONVERT FROM THETA TO T 
             if(tshltr(i,j)/=spval)GRID1(I,J)=TSHLTR(I,J)*(PSHLTR(I,J)*1.E-5)**CAPA
             IF(GRID1(I,J).LT.200)PRINT*,'ABNORMAL 2MT ',i,j,  &
             TSHLTR(I,J),PSHLTR(I,J)
!             TSHLTR(I,J)=GRID1(I,J) 
            ENDDO
            ENDDO
            ID(1:25) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(106),LVLS(1,IGET(106)),GRID1,IM,JM)
         ENDIF
!
!        SHELTER LEVEL SPECIFIC HUMIDITY.
         IF (IGET(112).GT.0) THEN       
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=QSHLTR(I,J)
            ENDDO
            ENDDO
            CALL BOUND (GRID1,H1M12,H99999)
            ID(1:25) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(112),LVLS(1,IGET(112)),GRID1,IM,JM)
         ENDIF
!     
!        SHELTER LEVEL DEWPOINT.
         IF (IGET(113).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM

!tgs The next 4 lines are GSD algorithm for Dew Point computation
!tgs Results is very close to dew point computed in DEWPOINT subroutine
!          qv = max(1.E-5,(QSHLTR(I,J)/(1.-QSHLTR(I,J))))
!          e=PSHLTR(I,J)/100.*qv/(0.62197+qv)
!          DWPT = (243.5*ALOG(E)-440.8)/(19.48-ALOG(E))+273.15
!          if(i.eq.ii.and.j.eq.jj)print*,'Debug: RUC-type DEWPT,i,j'
!     &,   DWPT,i,j,qv,pshltr(i,j),qshltr(i,j)
!          EGRID1(I,J)=DWPT

              EVP(I,J)=PSHLTR(I,J)*QSHLTR(I,J)/(EPS+ONEPS*QSHLTR(I,J))
	      EVP(I,J)=EVP(I,J)*D001
            ENDDO
            ENDDO
            CALL DEWPOINT(EVP,EGRID1)
	    GRID1=spval
            DO J=JSTA,JEND
            DO I=1,IM
             if(qshltr(i,j)/=spval)GRID1(I,J)=EGRID1(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(113),LVLS(1,IGET(113)),GRID1,IM,JM)
         ENDIF
!     
!        SHELTER LEVEL RELATIVE HUMIDITY.
         IF (IGET(114).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             P1D(I,J)=PSHLTR(I,J)
             T1D(I,J)=TSHLTR(I,J)*(PSHLTR(I,J)*1.E-5)**CAPA
             Q1D(I,J)=QSHLTR(I,J)
            ENDDO
            ENDDO
!            CALL CALRH(PSHLTR,TSHLTR,QSHLTR,EGRID1)
            IF(MODELNAME == 'GFS')THEN
              CALL CALRH_GFS(P1D,T1D,Q1D,EGRID1)
            ELSEIF(MODELNAME.EQ.'RAPR')THEN
              CALL CALRH_GSD(P1D,T1D,Q1D,EGRID1)
            ELSE
              CALL CALRH(P1D,T1D,Q1D,EGRID1)
            END IF
            DO J=JSTA,JEND
            DO I=1,IM
             if(qshltr(i,j)/=spval)then
	      GRID1(I,J)=EGRID1(I,J)*100.
	     else	     
	      grid1(i,j)=spval 
	     end if 
            ENDDO
            ENDDO
            CALL BOUND(GRID1,H1,H100)
            ID(1:25) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(114),LVLS(1,IGET(114)),GRID1,IM,JM)
         ENDIF
!     
!        SHELTER LEVEL PRESSURE.
         IF (IGET(138).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=PSHLTR(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(138),LVLS(1,IGET(138)),GRID1,IM,JM)
         ENDIF
!
      ENDIF
!
!        SHELTER LEVEL MAX TEMPERATURE.
         IF (IGET(345).GT.0) THEN       
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=MAXTSHLTR(I,J)
            ENDDO
            ENDDO
	    ID(1:25) = 0
	    ITMAXMIN     = INT(TMAXMIN)
            IF(ITMAXMIN .ne. 0) then
             IFINCR     = MOD(IFHR,ITMAXMIN)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITMAXMIN*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 2
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITMAXMIN
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(345),LVLS(1,IGET(345)),GRID1,IM,JM)
         ENDIF
!
!        SHELTER LEVEL MIN TEMPERATURE.
         IF (IGET(346).GT.0) THEN       
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=MINTSHLTR(I,J)
            ENDDO
            ENDDO
	    ID(1:25) = 0
	    ITMAXMIN     = INT(TMAXMIN)
            IF(ITMAXMIN .ne. 0) then
             IFINCR     = MOD(IFHR,ITMAXMIN)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITMAXMIN*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 2
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITMAXMIN
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(346),LVLS(1,IGET(346)),GRID1,IM,JM)
         ENDIF
!
!        SHELTER LEVEL MAX RH.
         IF (IGET(347).GT.0) THEN       
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=MAXRHSHLTR(I,J)*100.
            ENDDO
            ENDDO
	    ID(1:25) = 0
	    ID(02)=129
	    ITMAXMIN     = INT(TMAXMIN)
            IF(ITMAXMIN .ne. 0) then
             IFINCR     = MOD(IFHR,ITMAXMIN)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITMAXMIN*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 2
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITMAXMIN
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(347),LVLS(1,IGET(347)),GRID1,IM,JM)
         ENDIF
!
!        SHELTER LEVEL MIN RH.
         IF (IGET(348).GT.0) THEN       
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=MINRHSHLTR(I,J)*100.
            ENDDO
            ENDDO
	    ID(1:25) = 0
	    ID(02)=129
	    ITMAXMIN     = INT(TMAXMIN)
            IF(ITMAXMIN .ne. 0) then
             IFINCR     = MOD(IFHR,ITMAXMIN)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITMAXMIN*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 2
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITMAXMIN
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(348),LVLS(1,IGET(348)),GRID1,IM,JM)
         ENDIF
!
!
!     BLOCK 3.  ANEMOMETER LEVEL (10M) WINDS, THETA, AND Q.
!
      IF ( (IGET(064).GT.0).OR.(IGET(065).GT.0) ) THEN
         ID(1:25) = 0
         ISVALUE = 10
         ID(10) = MOD(ISVALUE/256,256)
         ID(11) = MOD(ISVALUE,256)
!
!        ANEMOMETER LEVEL U WIND AND/OR V WIND.
         IF ((IGET(064).GT.0).OR.(IGET(065).GT.0)) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=U10(I,J)
             GRID2(I,J)=V10(I,J)
            ENDDO
            ENDDO
            IF (IGET(064).GT.0) CALL GRIBIT(IGET(064),     &
                 LVLS(1,IGET(064)),GRID1,IM,JM)
            IF (IGET(065).GT.0) CALL GRIBIT(IGET(065),     &
                 LVLS(1,IGET(065)),GRID2,IM,JM)
         ENDIF
      ENDIF
!
!        ANEMOMETER LEVEL (10 M) POTENTIAL TEMPERATURE.
!   NOT A OUTPUT FROM WRF
      IF (IGET(158).GT.0) THEN
         ID(1:25) = 0
         ISVALUE = 10
         ID(10) = MOD(ISVALUE/256,256)
         ID(11) = MOD(ISVALUE,256)
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=TH10(I,J)
            ENDDO
            ENDDO
         CALL GRIBIT(IGET(158),LVLS(1,IGET(158)),GRID1,IM,JM)
       ENDIF
!
!        ANEMOMETER LEVEL (10 M) SPECIFIC HUMIDITY.
!
      IF (IGET(159).GT.0) THEN
         ID(1:25) = 0
         ISVALUE = 10
         ID(10) = MOD(ISVALUE/256,256)
         ID(11) = MOD(ISVALUE,256)
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=Q10(I,J)
            ENDDO
            ENDDO
         CALL GRIBIT(IGET(159),LVLS(1,IGET(159)),GRID1,IM,JM)
       ENDIF
! SRD
!
!        ANEMOMETER LEVEL (10 M) MAX WIND SPEED.
!
      IF (IGET(422).GT.0) THEN
         print *,' SRD ***** outputting WSPD10MAX '
         ID(1:25) = 0
         ISVALUE = 10
         ID(10) = MOD(ISVALUE/256,256)
         ID(11) = MOD(ISVALUE,256)
         ID(20) = 2
         ID(19) = IFHR
         IF (IFHR.EQ.0) THEN
           ID(18) = 0
         ELSE
           ID(18) = IFHR - 1
         ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=WSPD10MAX(I,J)
            ENDDO
            ENDDO
         CALL GRIBIT(IGET(422),LVLS(1,IGET(422)),GRID1,IM,JM)
       ENDIF
!
! SRD

!
!
!
!***  BLOCK 4.  PRECIPITATION RELATED FIELDS.
!MEB 6/17/02  ASSUMING THAT ALL ACCUMULATED FIELDS NEVER EMPTY
!             THEIR BUCKETS.  THIS IS THE EASIEST WAY TO DEAL WITH
!             ACCUMULATED FIELDS.  SHORTER TIME ACCUMULATIONS CAN
!             BE COMPUTED AFTER THE FACT IN A SEPARATE CODE ONCE
!             THE POST HAS FINISHED.  I HAVE LEFT IN THE OLD
!             ETAPOST CODE FOR COMPUTING THE BEGINNING TIME OF
!             THE ACCUMULATION PERIOD IF THIS IS CHANGED BACK
!             TO A 12H OR 3H BUCKET.  I AM NOT SURE WHAT
!             TO DO WITH THE TIME AVERAGED FIELDS, SO
!             LEAVING THAT UNCHANGED.
!     
!     SNOW FRACTION FROM EXPLICIT CLOUD SCHEME.  LABELLED AS
!      'PROB OF FROZEN PRECIP' IN GRIB, 
!      DIDN'T KNOW WHAT ELSE TO CALL IT
      IF (IGET(172).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
              IF (PREC(I,J) .LE. PTHRESH) THEN
                GRID1(I,J)=-50.
              ELSE
                GRID1(I,J)=SR(I,J)*100.
              ENDIF
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(172),LVLS(1,IGET(172)),GRID1,IM,JM)
      ENDIF
!     INSTANTANEOUS PRECIPITATION RATE.
      IF (IGET(167).GT.0) THEN
!MEB need to get physics DT
         RDTPHS=1./(DT * NPHS) 
!MEB need to get physics DT
            DO J=JSTA,JEND
            DO I=1,IM
             IF(MODELNAME .EQ. 'NCAR' .OR. MODELNAME.EQ.'RAPR')THEN
              GRID1(I,J)=PREC(I,J)/DT*1000.
             ELSE IF (MODELNAME .EQ. 'NMM' )THEN
              GRID1(I,J)=PREC(I,J)*RDTPHS*1000.
             ELSE IF (MODELNAME .EQ. 'RSM') THEN    !Add by Binbin 
              GRID1(I,J)=PREC(I,J)
             END IF
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(167),LVLS(1,IGET(167)),GRID1,IM,JM)
      ENDIF
!
!     INSTANTANEOUS CONVECTIVE PRECIPITATION RATE.
!     SUBSTITUTE WITH CUPPT IN WRF FOR NOW
      IF (IGET(249).GT.0) THEN
         RDTPHS=1000./DTQ2     !--- 1000 kg/m**3, density of liquid water
!         RDTPHS=1000./(TRDLW*3600.)
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J)=CPRATE(I,J)*RDTPHS
!           GRID1(I,J)=CUPPT(I,J)*RDTPHS
         ENDDO
         ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(249),LVLS(1,IGET(249)),GRID1,IM,JM)
      ENDIF
!
!     INSTANTANEOUS PRECIPITATION RATE.
      IF (IGET(167).GT.0) THEN
!MEB need to get physics DT
         RDTPHS=1./(DT * NPHS) 
!MEB need to get physics DT
            DO J=JSTA,JEND
            DO I=1,IM
             IF(MODELNAME .EQ. 'NCAR' .OR. MODELNAME.EQ.'RAPR')THEN
              GRID1(I,J)=PREC(I,J)/DT*1000.
             ELSE IF (MODELNAME .EQ. 'NMM' .OR.        &
             MODELNAME .EQ. 'GFS')THEN                  
              GRID1(I,J)=PREC(I,J)*RDTPHS*1000.
             ELSE IF (MODELNAME .EQ. 'RSM') THEN    !Add by Binbin 
              GRID1(I,J)=PREC(I,J)
             END IF
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(167),LVLS(1,IGET(167)),GRID1,IM,JM)
      ENDIF
!
!     TIME-AVERAGED CONVECTIVE PRECIPITATION RATE.
      IF (IGET(272).GT.0) THEN
         RDTPHS=1000./DTQ2     !--- 1000 kg/m**3, density of liquid water
         ID(1:25) = 0
         ITPREC     = NINT(TPREC)
!mp
	 if (ITPREC .ne. 0) then
          IFINCR     = MOD(IFHR,ITPREC)
	  IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	 else
	  IFINCR     = 0
	 endif
!mp
         ID(18)     = 0
         ID(19)     = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)     = 3
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITPREC
         ELSE
          ID(18) = IFHR-IFINCR
	  IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
	 grid1=spval
         DO J=JSTA,JEND
         DO I=1,IM
           if(AVGCPRATE(I,J)/=spval)GRID1(I,J)=AVGCPRATE(I,J)*RDTPHS
         ENDDO
         ENDDO
        
         CALL GRIBIT(IGET(272),LVLS(1,IGET(272)),GRID1,IM,JM)
      ENDIF
!      
!     TIME-AVERAGED PRECIPITATION RATE.
      IF (IGET(271).GT.0) THEN
         RDTPHS=1000./DTQ2     !--- 1000 kg/m**3, density of liquid water
!         RDTPHS=1000./(TRDLW*3600.)
         ID(1:25) = 0
         ITPREC     = NINT(TPREC)
!mp
	 if (ITPREC .ne. 0) then
          IFINCR     = MOD(IFHR,ITPREC)
	  IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	 else
	  IFINCR     = 0
	 endif
!mp
         ID(18)     = 0
         ID(19)     = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)     = 3
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITPREC
         ELSE
          ID(18) = IFHR-IFINCR
	  IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
	 grid1=spval
         DO J=JSTA,JEND
         DO I=1,IM
           if(avgprec(i,j)/=spval)GRID1(I,J)=AVGPREC(I,J)*RDTPHS
         ENDDO
         ENDDO
        
         CALL GRIBIT(IGET(271),LVLS(1,IGET(271)),GRID1,IM,JM)
      ENDIF
!     
!     ACCUMULATED TOTAL PRECIPITATION.
      IF (IGET(087).GT.0) THEN
         ID(1:25) = 0
         ITPREC     = NINT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
	 IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
         ID(18)     = 0
         ID(19)     = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)     = 4
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITPREC
         ELSE
          ID(18) = IFHR-IFINCR
	  IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
	 IF(MODELNAME == 'GFS') THEN
	   DO J=JSTA,JEND
            DO I=1,IM
             IF(AVGPREC(I,J) < SPVAL)THEN
	      GRID1(I,J)=AVGPREC(I,J)*                   &
      	        FLOAT(ID(19)-ID(18))*3600.*1000./DTQ2
             ELSE
	      GRID1(I,J)=SPVAL
	     END IF 
            ENDDO
          ENDDO
	 ELSE   
	  DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J)=ACPREC(I,J)*1000.
            ENDDO
          ENDDO
	 END IF 
!	 IF(IFMIN .GE. 1 .AND. ID(19) .GT. 256)THEN
!	  IF(ITPREC.EQ.3)ID(17)=10
!	  IF(ITPREC.EQ.6)ID(17)=11
!	  IF(ITPREC.EQ.12)ID(17)=12
!	 END IF 
         IF (ID(18).LT.0) ID(18) = 0
!	write(6,*) 'call gribit...total precip'
         CALL GRIBIT(IGET(087),LVLS(1,IGET(087)),GRID1,IM,JM)
      ENDIF
!     
!     ACCUMULATED CONVECTIVE PRECIPITATION.
      IF (IGET(033).GT.0) THEN
         ID(1:25) = 0
         ITPREC     = NINT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
         ID(18)     = 0
         ID(19)     = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)     = 4
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITPREC
         ELSE
          ID(18) = IFHR-IFINCR
          IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
	 IF(MODELNAME == 'GFS') THEN
	   DO J=JSTA,JEND
            DO I=1,IM
             IF(AVGCPRATE(I,J) < SPVAL)THEN
	      GRID1(I,J)=AVGCPRATE(I,J)*                      &
      	        FLOAT(ID(19)-ID(18))*3600.*1000./DTQ2
             ELSE
	      GRID1(I,J)=SPVAL
	     END IF 
            ENDDO
          ENDDO
	 ELSE
	   DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=CUPREC(I,J)*1000.
            ENDDO
           ENDDO
	 END IF 
!	write(6,*) 'call gribit...convective precip'
         CALL GRIBIT(IGET(033),LVLS(1,IGET(033)),GRID1,IM,JM)
      ENDIF
!     
!     ACCUMULATED GRID-SCALE PRECIPITATION.
      IF (IGET(034).GT.0) THEN
            
         ID(1:25) = 0
         ITPREC     = NINT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
         ID(18)     = 0
         ID(19)     = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)     = 4
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITPREC
         ELSE
          ID(18) = IFHR-IFINCR
          IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
	 IF(MODELNAME == 'GFS') THEN
	   DO J=JSTA,JEND
            DO I=1,IM
             IF(AVGCPRATE(I,J) < SPVAL .AND. AVGPREC(I,J) < SPVAL)   &
               then 	     
     	      GRID1(I,J)=( AVGPREC(I,J) - AVGCPRATE(I,J) )*          &
      	        FLOAT(ID(19)-ID(18))*3600.*1000./DTQ2
             ELSE
	      GRID1(I,J)=SPVAL
	     END IF 
            ENDDO
          ENDDO
	 ELSE
	   DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=ANCPRC(I,J)*1000.
            ENDDO
           ENDDO
	 END IF  
!	write(6,*) 'call gribit...grid-scale precip'
         CALL GRIBIT(IGET(034),LVLS(1,IGET(034)), GRID1,IM,JM)
      ENDIF
!     
!     ACCUMULATED LAND SURFACE PRECIPITATION.
      IF (IGET(256).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=LSPA(I,J)*1000.
            ENDDO
            ENDDO
         ID(1:25) = 0
         ITPREC     = NINT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
         ID(18)     = 0
         ID(19)     = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)     = 4
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITPREC
         ELSE
          ID(18) = IFHR-IFINCR
          IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         ID(02)= 130
         CALL GRIBIT(IGET(256),LVLS(1,IGET(256)), GRID1,IM,JM)
      ENDIF
!     
!     ACCUMULATED SNOWFALL.
         IF (IGET(035).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
!             GRID1(I,J)=ACSNOW(I,J)*1000.
             GRID1(I,J)=ACSNOW(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
      print*,'maxval ACCUMULATED SNOWFALL: ', maxval(GRID1)
            CALL GRIBIT(IGET(035),LVLS(1,IGET(035)),GRID1,IM,JM)
         ENDIF
!     
!     ACCUMULATED SNOW MELT.
         IF (IGET(121).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
!             GRID1(I,J)=ACSNOM(I,J)*1000.
             GRID1(I,J)=ACSNOM(I,J)	     
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(121),LVLS(1,IGET(121)),GRID1,IM,JM)
         ENDIF
!     
!     ACCUMULATED SNOWFALL RATE
         IF (IGET(405).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SNOWFALL(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	    IF(ITPREC < 0)ID(1:25)=0
            CALL GRIBIT(IGET(405),LVLS(1,IGET(405)),GRID1,IM,JM)
         ENDIF
!     
!     ACCUMULATED STORM SURFACE RUNOFF.
         IF (IGET(122).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
!             GRID1(I,J)=SSROFF(I,J)*1000.
             GRID1(I,J)=SSROFF(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(122),LVLS(1,IGET(122)), GRID1,IM,JM)
         ENDIF
!     
!     ACCUMULATED BASEFLOW-GROUNDWATER RUNOFF.
         IF (IGET(123).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
!             GRID1(I,J)=BGROFF(I,J)*1000.
             GRID1(I,J)=BGROFF(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(123),LVLS(1,IGET(123)),GRID1,IM,JM)
         ENDIF
!     
!     ACCUMULATED WATER RUNOFF.
         IF (IGET(343).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=RUNOFF(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
	    if (ITPREC .ne. 0) then
             IFINCR     = MOD(IFHR,ITPREC)
             IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	    else
	     IFINCR     = 0
	    endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(343),LVLS(1,IGET(343)),GRID1,IM,JM)
         ENDIF

!     PRECIPITATION BUCKETS - accumulated between output times
!     'BUCKET TOTAL PRECIP '
         IF (IGET(434).GT.0.) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=PCP_BUCKET(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
            if (ITPREC .ne. 0) then
             IFINCR     = MOD(IFHR,ITPREC)
             IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
            else
             IFINCR     = 0
            endif
!mp
           if(MODELNAME.EQ.'NCAR' .OR. MODELNAME.EQ.'RAPR') IFINCR = NINT(PREC_ACC_DT)/60
            ID(18)     = 0
            ID(19)     = IFHR
            IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(434),LVLS(1,IGET(434)),GRID1,IM,JM)
         ENDIF

!     PRECIPITATION BUCKETS - accumulated between output times
!     'BUCKET CONV PRECIP  '
         IF (IGET(435).GT.0.) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=RAINC_BUCKET(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
            if (ITPREC .ne. 0) then
             IFINCR     = MOD(IFHR,ITPREC)
             IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
            else
             IFINCR     = 0
            endif

           if(MODELNAME.EQ.'NCAR' .OR. MODELNAME.EQ.'RAPR') IFINCR = NINT(PREC_ACC_DT)/60
!mp
            ID(18)     = 0
            ID(19)     = IFHR
            IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0

!  print *,'IFMIN,IFHR,ITPREC',IFMIN,IFHR,ITPREC
  print *,'PREC_ACC_DT,ID(18),ID(19)',PREC_ACC_DT,ID(18),ID(19)

            CALL GRIBIT(IGET(435),LVLS(1,IGET(435)),GRID1,IM,JM)
         ENDIF
!     PRECIPITATION BUCKETS - accumulated between output times
!     'BUCKET GRDSCALE PRCP'
         IF (IGET(436).GT.0.) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=RAINNC_BUCKET(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
            if (ITPREC .ne. 0) then
             IFINCR     = MOD(IFHR,ITPREC)
             IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
            else
             IFINCR     = 0
            endif
!mp
           if(MODELNAME.EQ.'NCAR' .OR. MODELNAME.EQ.'RAPR') IFINCR = NINT(PREC_ACC_DT)/60
            ID(18)     = 0
            ID(19)     = IFHR
            IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(436),LVLS(1,IGET(436)),GRID1,IM,JM)
         ENDIF
!     PRECIPITATION BUCKETS - accumulated between output times
!     'BUCKET SNOW  PRECIP '
         IF (IGET(437).GT.0.) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SNOW_BUCKET(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
            if (ITPREC .ne. 0) then
             IFINCR     = MOD(IFHR,ITPREC)
             IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
            else
             IFINCR     = 0
            endif
!mp
           if(MODELNAME.EQ.'NCAR' .OR. MODELNAME.EQ.'RAPR') IFINCR = NINT(PREC_ACC_DT)/60
            ID(18)     = 0
            ID(19)     = IFHR
            IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
      print*,'maxval BUCKET SNOWFALL: ', maxval(GRID1)
            CALL GRIBIT(IGET(437),LVLS(1,IGET(437)),GRID1,IM,JM)
         ENDIF




!     
!     INSTANTANEOUS PRECIPITATION TYPE.
         IF (IGET(160).GT.0 .OR.(IGET(247).GT.0)) THEN

          CALL CALWXT(T,Q,PMID,PINT,HTM,LMH,PREC,ZINT,IWX1    &
             ,ZWET)
          IF (IGET(160).GT.0) THEN 
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=IWX1(I,J)
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW1(I,J)   = ISNO*1.0
              SLEET1(I,J)  = IIP*1.0
              FREEZR1(I,J) = IZR*1.0
              RAIN1(I,J)   = IRAIN*1.0
	      SNOW(I,J,1)    = SNOW1(I,J)
              SLEET(I,J,1)   = SLEET1(I,J)
              FREEZR(I,J,1) = FREEZR1(I,J)
              RAIN(I,J,1)    = RAIN1(I,J)
            ENDDO
            ENDDO
          ENDIF
!     
!     LOWEST WET BULB ZERO HEIGHT
           IF (IGET(247).GT.0) THEN
             DO J=JSTA,JEND
             DO I=1,IM
              GRID1(I,J)=ZWET(I,J)
             ENDDO
             ENDDO
             ID(1:25) = 0
             CALL GRIBIT(IGET(247),LVLS(1,IGET(247)),GRID1,IM,JM)
           ENDIF

!     DOMINANT PRECIPITATION TYPE
!GSM  IF DOMINANT PRECIP TYPE IS REQUESTED, 4 MORE ALGORITHMS
!GSM    WILL BE CALLED.  THE TALLIES ARE THEN SUMMED IN
!GSM    CALWXT_DOMINANT

           IF (IGET(160).GT.0) THEN   
!  RAMER ALGORITHM
            CALL CALWXT_RAMER(T,Q,PMID,PINT,LMH,PREC,IWX2)
          print *,'in SURFCE,me=',me,'IWX2=',IWX2(1:30,JSTA)
               
!     DECOMPOSE IWX2 ARRAY
!
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=NINT(IWX2(I,J))
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW2(I,J)   = ISNO*1.0
              SLEET2(I,J)  = IIP*1.0
              FREEZR2(I,J) = IZR*1.0
              RAIN2(I,J)   = IRAIN*1.0
              SNOW(I,J,2)    = SNOW2(I,J)
              SLEET(I,J,2)   = SLEET2(I,J)
              FREEZR(I,J,2) = FREEZR2(I,J)
              RAIN(I,J,2)    = RAIN2(I,J)
            ENDDO
            ENDDO

! BOURGOUIN ALGORITHM
            ISEED=44641*(INT(SDAT(1)-1)*24*31+INT(SDAT(2))*24+IHRST)+   &
     &            MOD(IFHR*60+IFMIN,44641)+4357
            CALL CALWXT_BOURG(IM,JM,JSTA_2L,JEND_2U,JSTA,JEND,LM,LP1,   &
     &                        ISEED,G,PTHRESH,                          &
     &                        T,Q,PMID,PINT,LMH,PREC,ZINT,IWX3)
          print *,'in SURFCE,me=',me,'IWX3=',IWX3(1:30,JSTA),'PTHRESH=',PTHRESH

!     DECOMPOSE IWX3 ARRAY
!
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=NINT(IWX3(I,J))
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW3(I,J)   = ISNO*1.0
              SLEET3(I,J)  = IIP*1.0
              FREEZR3(I,J) = IZR*1.0
              RAIN3(I,J)   = IRAIN*1.0
              SNOW(I,J,3)    = SNOW3(I,J)
              SLEET(I,J,3)   = SLEET3(I,J)
              FREEZR(I,J,3) = FREEZR3(I,J)
              RAIN(I,J,3)    = RAIN3(I,J)
            ENDDO
            ENDDO

! REVISED NCEP ALGORITHM
            CALL CALWXT_REVISED(T,Q,PMID,PINT,HTM,LMH,PREC,ZINT,   &
                IWX4)
          print *,'in SURFCE,me=',me,'IWX4=',IWX4(1:30,JSTA)
!     DECOMPOSE IWX2 ARRAY
!
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=IWX4(I,J)
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW4(I,J)   = ISNO*1.0
              SLEET4(I,J)  = IIP*1.0
              FREEZR4(I,J) = IZR*1.0
              RAIN4(I,J)   = IRAIN*1.0
              SNOW(I,J,4)    = SNOW4(I,J)
              SLEET(I,J,4)   = SLEET4(I,J)
              FREEZR(I,J,4) = FREEZR4(I,J)
              RAIN(I,J,4)    = RAIN4(I,J)
            ENDDO
            ENDDO
              
! EXPLICIT ALGORITHM (UNDER 18 NOT ADMITTED WITHOUT PARENT 
!     OR GUARDIAN)
 
            IF(imp_physics==5)then
             CALL CALWXT_EXPLICIT(LMH,THS,PMID,PREC,SR,F_RimeF,IWX5)
            else
             IWX5=0
            end if
          print *,'in SURFCE,me=',me,'IWX5=',IWX5(1:30,JSTA)
!     DECOMPOSE IWX2 ARRAY
!
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=IWX5(I,J)
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW5(I,J)   = ISNO*1.0
              SLEET5(I,J)  = IIP*1.0
              FREEZR5(I,J) = IZR*1.0
              RAIN5(I,J)   = IRAIN*1.0
              SNOW(I,J,5)    = SNOW5(I,J)
              SLEET(I,J,5)   = SLEET5(I,J)
              FREEZR(I,J,5) = FREEZR5(I,J)
              RAIN(I,J,5)    = RAIN5(I,J)
            ENDDO
            ENDDO
               
           CALL CALWXT_DOMINANT(PREC,RAIN,FREEZR,SLEET,SNOW,      &
               DOMR,DOMZR,DOMIP,DOMS)
           if ( me.eq.0) print *,'after CALWXT_DOMINANT, no avrg'
           ID(1:25) = 0
!     SNOW.
            ID(8) = 143 
	    grid1=spval
            DO J=JSTA,JEND
            DO I=1,IM
             if(prec(i,j)/=spval)GRID1(I,J)=DOMS(I,J)
            ENDDO
            ENDDO
            CALL GRIBIT(IGET(160),LVLS(1,IGET(160)),GRID1,IM,JM)
!     ICE PELLETS.
            ID(8) = 142 
	    grid1=spval
            DO J=JSTA,JEND
            DO I=1,IM
             if(prec(i,j)/=spval)GRID1(I,J)=DOMIP(I,J)
            ENDDO
            ENDDO
            CALL GRIBIT(IGET(160),LVLS(1,IGET(160)),GRID1,IM,JM)
!     FREEZING RAIN.
            ID(8) = 141
	    grid1=spval 
            DO J=JSTA,JEND
            DO I=1,IM
!             if (DOMZR(I,J) .EQ. 1) THEN
!               PSFC(I,J)=PINT(I,J,NINT(LMH(I,J))+1)
!               print *, 'aha ', I, J, PSFC(I,J)
!               print *, FREEZR(I,J,1), FREEZR(I,J,2),
!     *  FREEZR(I,J,3), FREEZR(I,J,4), FREEZR(I,J,5)
!             endif
             if(prec(i,j)/=spval)GRID1(I,J)=DOMZR(I,J)
            ENDDO
            ENDDO
            CALL GRIBIT(IGET(160),LVLS(1,IGET(160)),GRID1,IM,JM)
!     RAIN.
            ID(8) = 140
	    grid1=spval 
            DO J=JSTA,JEND
            DO I=1,IM
             if(prec(i,j)/=spval)GRID1(I,J)=DOMR(I,J)
            ENDDO
            ENDDO
	    CALL GRIBIT(IGET(160),LVLS(1,IGET(160)),GRID1,IM,JM)
        ENDIF
      ENDIF
!     
!     TIME AVERAGED PRECIPITATION TYPE.
         IF (IGET(317).GT.0) THEN

          CALL CALWXT(T,Q,PMID,PINT,HTM,LMH,AVGPREC,ZINT,IWX1   &
             ,ZWET)
          print *,'in SURFCE,me=',me,'IWX1=',IWX1(1:30,JSTA)
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=IWX1(I,J)
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW1(I,J)   = ISNO*1.0
              SLEET1(I,J)  = IIP*1.0
              FREEZR1(I,J) = IZR*1.0
              RAIN1(I,J)   = IRAIN*1.0
	      SNOW(I,J,1)    = SNOW1(I,J)
              SLEET(I,J,1)   = SLEET1(I,J)
              FREEZR(I,J,1) = FREEZR1(I,J)
              RAIN(I,J,1)    = RAIN1(I,J)
            ENDDO
            ENDDO

!     DOMINANT PRECIPITATION TYPE
!GSM  IF DOMINANT PRECIP TYPE IS REQUESTED, 4 MORE ALGORITHMS
!GSM    WILL BE CALLED.  THE TALLIES ARE THEN SUMMED IN
!GSM    CALWXT_DOMINANT

!  RAMER ALGORITHM
            CALL CALWXT_RAMER(T,Q,PMID,PINT,LMH,AVGPREC,IWX2)
          print *,'in SURFCE,me=',me,'IWX2=',IWX2(1:30,JSTA)
               
!     DECOMPOSE IWX2 ARRAY
!
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=NINT(IWX2(I,J))
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW2(I,J)   = ISNO*1.0
              SLEET2(I,J)  = IIP*1.0
              FREEZR2(I,J) = IZR*1.0
              RAIN2(I,J)   = IRAIN*1.0
              SNOW(I,J,2)    = SNOW2(I,J)
              SLEET(I,J,2)   = SLEET2(I,J)
              FREEZR(I,J,2) = FREEZR2(I,J)
              RAIN(I,J,2)    = RAIN2(I,J)
            ENDDO
            ENDDO

! BOURGOUIN ALGORITHM
            ISEED=44641*(INT(SDAT(1)-1)*24*31+INT(SDAT(2))*24+IHRST)+   &
     &            MOD(IFHR*60+IFMIN,44641)+4357
            CALL CALWXT_BOURG(IM,JM,JSTA_2L,JEND_2U,JSTA,JEND,LM,LP1,   &
     &                        ISEED,G,PTHRESH,                          &
     &                        T,Q,PMID,PINT,LMH,AVGPREC,ZINT,IWX3)
          print *,'in SURFCE,me=',me,'IWX3=',IWX3(1:30,JSTA)

!     DECOMPOSE IWX3 ARRAY
!
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=NINT(IWX3(I,J))
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW3(I,J)   = ISNO*1.0
              SLEET3(I,J)  = IIP*1.0
              FREEZR3(I,J) = IZR*1.0
              RAIN3(I,J)   = IRAIN*1.0
              SNOW(I,J,3)    = SNOW3(I,J)
              SLEET(I,J,3)   = SLEET3(I,J)
              FREEZR(I,J,3) = FREEZR3(I,J)
              RAIN(I,J,3)    = RAIN3(I,J)
            ENDDO
            ENDDO

! REVISED NCEP ALGORITHM
            CALL CALWXT_REVISED(T,Q,PMID,PINT,HTM,LMH,AVGPREC,ZINT,  &
                IWX4)
          print *,'in SURFCE,me=',me,'IWX4=',IWX4(1:30,JSTA)
!     DECOMPOSE IWX2 ARRAY
!
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=IWX4(I,J)
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW4(I,J)   = ISNO*1.0
              SLEET4(I,J)  = IIP*1.0
              FREEZR4(I,J) = IZR*1.0
              RAIN4(I,J)   = IRAIN*1.0
              SNOW(I,J,4)    = SNOW4(I,J)
              SLEET(I,J,4)   = SLEET4(I,J)
              FREEZR(I,J,4) = FREEZR4(I,J)
              RAIN(I,J,4)    = RAIN4(I,J)
            ENDDO
            ENDDO
              
! EXPLICIT ALGORITHM (UNDER 18 NOT ADMITTED WITHOUT PARENT 
!     OR GUARDIAN)
 
            IF(imp_physics==5)then
             CALL CALWXT_EXPLICIT(LMH,THS,PMID,AVGPREC,SR,F_RimeF,IWX5)
            else
             IWX5=0
            end if
          print *,'in SURFCE,me=',me,'IWX5=',IWX5(1:30,JSTA)
!     DECOMPOSE IWX2 ARRAY
!
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=IWX5(I,J)
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW5(I,J)   = ISNO*1.0
              SLEET5(I,J)  = IIP*1.0
              FREEZR5(I,J) = IZR*1.0
              RAIN5(I,J)   = IRAIN*1.0
              SNOW(I,J,5)    = SNOW5(I,J)
              SLEET(I,J,5)   = SLEET5(I,J)
              FREEZR(I,J,5) = FREEZR5(I,J)
              RAIN(I,J,5)    = RAIN5(I,J)
            ENDDO
            ENDDO
               
            print *,'me=',me,'before SNOW=',snow(1:10,JSTA,1:5)
            print *,'me=',me,'before RAIN=',RAIN(1:10,JSTA,1:5)
            print *,'me=',me,'before FREEZR=',FREEZR(1:10,JSTA,1:5)
            print *,'me=',me,'before SLEET=',SLEET(1:10,JSTA,1:5)
           CALL CALWXT_DOMINANT(AVGPREC,RAIN,FREEZR,SLEET,SNOW,    &
               DOMR,DOMZR,DOMIP,DOMS)
     
           ID(1:25) = 0
           ITPREC     = NINT(TPREC)
!mp
	   if (ITPREC .ne. 0) then
            IFINCR     = MOD(IFHR,ITPREC)
	    IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	   else
	    IFINCR     = 0
	   endif
!mp
           ID(18)     = 0
           ID(19)     = IFHR
	   IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
           ID(20)     = 3
           IF (IFINCR.EQ.0) THEN
            ID(18) = IFHR-ITPREC
           ELSE
            ID(18) = IFHR-IFINCR
	    IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
           ENDIF

           if ( me.eq.0) print *,'after CALWXT_DOMINANT, avrg,TPREC=', &
             TPREC,'IFHR=',IFHR,'IFMIN=',IFMIN,'IFINCR=',IFINCR,'ID=',ID
!     SNOW.
            
            ID(8) = 143 
	    grid1=spval
            DO J=JSTA,JEND
            DO I=1,IM
             if(avgprec(i,j)/=spval)GRID1(I,J)=DOMS(I,J)
            ENDDO
            ENDDO
!            print *,'me=',me,'SNOW=',GRID1(1:10,JSTA)
            CALL GRIBIT(IGET(317),LVLS(1,IGET(317)),GRID1,IM,JM)
!     ICE PELLETS.
            ID(8) = 142 
	    ITPREC     = NINT(TPREC)
!mp
	    if (ITPREC .ne. 0) then
             IFINCR     = MOD(IFHR,ITPREC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	    else
	     IFINCR     = 0
	    endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
	     IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
	    grid1=spval
            DO J=JSTA,JEND
            DO I=1,IM
             if(avgprec(i,j)/=spval)GRID1(I,J)=DOMIP(I,J)
            ENDDO
            ENDDO
            CALL GRIBIT(IGET(317),LVLS(1,IGET(317)),GRID1,IM,JM)
!     FREEZING RAIN.
            ID(8) = 141
	    
	    ITPREC     = NINT(TPREC)
!mp
	    if (ITPREC .ne. 0) then
             IFINCR     = MOD(IFHR,ITPREC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	    else
	     IFINCR     = 0
	    endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
	     IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
	    grid1=spval
            DO J=JSTA,JEND
            DO I=1,IM
!             if (DOMZR(I,J) .EQ. 1) THEN
!               PSFC(I,J)=PINT(I,J,NINT(LMH(I,J))+1)
!               print *, 'aha ', I, J, PSFC(I,J)
!               print *, FREEZR(I,J,1), FREEZR(I,J,2),
!     *  FREEZR(I,J,3), FREEZR(I,J,4), FREEZR(I,J,5)
!             endif
             if(avgprec(i,j)/=spval)GRID1(I,J)=DOMZR(I,J)
            ENDDO
            ENDDO
            CALL GRIBIT(IGET(317),LVLS(1,IGET(317)),GRID1,IM,JM)
!     RAIN.
            ID(8) = 140
	    
	    ITPREC     = NINT(TPREC)
!mp
	    if (ITPREC .ne. 0) then
             IFINCR     = MOD(IFHR,ITPREC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	    else
	     IFINCR     = 0
	    endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
	     IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
	    grid1=spval 
            DO J=JSTA,JEND
            DO I=1,IM
             if(avgprec(i,j)/=spval)GRID1(I,J)=DOMR(I,J)
            ENDDO
            ENDDO
	    CALL GRIBIT(IGET(317),LVLS(1,IGET(317)),GRID1,IM,JM)

      ENDIF


! GSD PRECIPITATION TYPE
         IF (IGET(407).GT.0) THEN

            DO J=JSTA,JEND
            DO I=1,IM
!-- snow
               DOMS(I,J) =0.
!-- rain
               DOMR(I,J) =0.
!-- freezing rain
               DOMZR(I,J)=0.
!-- ice pellets
               DOMIP(I,J)=0.
            ENDDO
            ENDDO

            DO J=JSTA,JEND
            DO I=1,IM
!-- TOTPRCP is total 1-hour accumulated precipitation in  [m] 
              totprcp = (RAINC_BUCKET(I,J) + RAINNC_BUCKET(I,J))*1.e-3
              snowratio = 0.0
!               if (totprcp.gt.0.01) 
               if (totprcp.gt.0.001)                               &
              snowratio = snow_bucket(i,j)*1.e-3/totprcp

!              snowratio = SR(i,j)
!-- 2-m temperature
               t2=TSHLTR(I,J)*(PSHLTR(I,J)*1.E-5)**CAPA
!--snow
!--   SNOW is time step non-convective snow [m]
             if(SNOWNC(i,j)/DT .gt. 0.2e-9 .or.                    &
               (totprcp.gt.0.1.and.snowratio.ge.0.25)) DOMS(i,j)=1.

!-- rain/freezing rain
!--   compute RAIN [m/s] from total convective and non-convective precipitation
               rainl = (1. - SR(i,j))*prec(i,j)/DT
!-- in RUC RAIN is in cm/h and the limit is 1.e-3, 
!-- converted to m/s will be 2.8e-9
             if((rainl .gt. 2.8e-9) .or.                           &
               (totprcp.gt.0.1 .and.snowratio.lt.0.75)) then

               if (t2.ge.273.15) then
!--rain
                  DOMR(I,J) = 1.
               else if (tmax(i,j).gt.273.15) then
!  Only allow frz rain if some level above
!               ground is above freezing.
!-- freezing rain
                  DOMZR(I,J) = 1.
               endif
             endif
!-- graupel/ice pellets/snow
!-- GRAUPEL is time step non-convective graupel in [m]
             if(GRAUPELNC(i,j)/DT .gt. 1.e-9) then
               if (t2.le.273.15) then
!              check for max rain mixing ratio
!              if it's > 0.05 g/kg, => ice pellets
               if (qrmax(i,j).gt.0.00005) then
                 if(GRAUPELNC(i,j) .gt. SNOWNC(i,j)) then
!-- ice pellets
                 DOMIP(I,J) = 1.

! -- If graupel is greater than rain,
!        report graupel only
! in RUC --> if (3.6E5*gex2(i,j,8).gt.   gex2(i,j,6)) then
                  if ((GRAUPELNC(i,j)/DT) .gt. rainl) then
                    DOMIP(I,J) = 1.
                    DOMZR(I,J) = 0.
                    DOMR(I,J)  = 0.
! -- If rain is greater than 4x graupel,
!        report rain only
! in RUC -->  else if (gex2(i,j,6).gt.4.*3.6E5*gex2(i,j,8)) then
                  else if (rainl .gt. (4.*GRAUPELNC(i,j)/DT)) then
                    DOMIP(I,J) = 0.
                  end if

                else
!              snow
                  DOMS(i,j)=1.
                end if
              else
!              snow
                DOMS(i,j)=1.
              end if
            else
              if (t2.ge.276.15) then
!              rain
                DOMR(I,J) = 1.
              else
!              snow
                DOMS(i,j)=1.
              end if
            end if
          end if
            ENDDO
            ENDDO

           ID(1:25) = 0
!     SNOW.
            ID(8) = 143
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=DOMS(I,J)
            ENDDO
            ENDDO
      print*,'maxval SNOW: ', maxval(GRID1)
            CALL GRIBIT(IGET(407),LVLS(1,IGET(407)),       &
                 GRID1,IM,JM)
!     ICE PELLETS.
            ID(8) = 142
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=DOMIP(I,J)
!             if (DOMIP(I,J) .EQ. 1) THEN
!               print *, 'ICE PELLETS at I,J ', I, J
!             endif
            ENDDO
            ENDDO
      print*,'maxval ICE PELLETS: ', maxval(GRID1)
            CALL GRIBIT(IGET(407),LVLS(1,IGET(407)),       &
                 GRID1,IM,JM)
!     FREEZING RAIN.
            ID(8) = 141
            DO J=JSTA,JEND
            DO I=1,IM
!             if (DOMZR(I,J) .EQ. 1) THEN
!               PSFC(I,J)=PINT(I,J,NINT(LMH(I,J))+1)
!               print *, 'FREEZING RAIN AT I,J ', I, J, PSFC(I,J)
!             endif
             GRID1(I,J)=DOMZR(I,J)
            ENDDO
            ENDDO
      print*,'maxval FREEZING RAIN: ', maxval(GRID1)
            CALL GRIBIT(IGET(407),LVLS(1,IGET(407)),       &
                 GRID1,IM,JM)
!     RAIN.
            ID(8) = 140
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=DOMR(I,J)
            ENDDO
            ENDDO
      print*,'maxval RAIN: ', maxval(GRID1)
            CALL GRIBIT(IGET(407),LVLS(1,IGET(407)),       &
                 GRID1,IM,JM)

        ENDIF

!     
!
!
!***  BLOCK 5.  SURFACE EXCHANGE FIELDS.
!     
!     TIME AVERAGED SURFACE LATENT HEAT FLUX.
         IF (IGET(042).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. &
             MODELNAME.EQ.'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE  
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
	     IF(SFCLHX(I,J)/=SPVAL)THEN
              GRID1(I,J)=-1.*SFCLHX(I,J)*RRNUM !change the sign to conform with Grib
	     ELSE
	      GRID1(I,J)=SFCLHX(I,J)
	     END IF 
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
	    IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(042),LVLS(1,IGET(042)),GRID1,IM,JM)
          END IF 
         ENDIF
!
!     TIME AVERAGED SURFACE SENSIBLE HEAT FLUX.
         IF (IGET(043).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. &
             MODELNAME.EQ.'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
	     IF(SFCSHX(I,J)/=SPVAL)THEN
              GRID1(I,J) = -1.* SFCSHX(I,J)*RRNUM !change the sign to conform with Grib
	     ELSE
	      GRID1(I,J)=SFCSHX(I,J)
	     END IF  
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
	    IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(043),LVLS(1,IGET(043)),GRID1,IM,JM)
         ENDIF
!     
!     TIME AVERAGED SUB-SURFACE SENSIBLE HEAT FLUX.
         IF (IGET(135).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. &
             MODELNAME.EQ.'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J) = SUBSHX(I,J)*RRNUM
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
            IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(135),LVLS(1,IGET(135)),GRID1,IM,JM)
         ENDIF
!     
!     TIME AVERAGED SNOW PHASE CHANGE HEAT FLUX.
         IF (IGET(136).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME.EQ.'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J) = SNOPCX(I,J)*RRNUM
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
            IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(136),LVLS(1,IGET(136)),GRID1,IM,JM)
         ENDIF
!     
!     TIME AVERAGED SURFACE MOMENTUM FLUX.
         IF (IGET(046).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR.  &
             MODELNAME.EQ.'RAPR')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
	     IF(SFCUVX(I,J)/=SPVAL)THEN
              GRID1(I,J) = SFCUVX(I,J)*RRNUM
	     ELSE
	      GRID1(I,J) = SFCUVX(I,J)
	     END IF  
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
            IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(046),LVLS(1,IGET(046)),GRID1,IM,JM)
         ENDIF
!
!     TIME AVERAGED SURFACE ZONAL MOMENTUM FLUX.
         IF (IGET(269).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. &
             MODELNAME.EQ.'RAPR')THEN
            GRID1=SPVAL
            ID(1:25)=0
          ELSE
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J) = SFCUX(I,J)*RRNUM
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
            IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
             IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
            ELSE
             IFINCR     = 0
            endif
            ID(19)     = IFHR
            IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
               IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
          END IF
          CALL GRIBIT(IGET(269),LVLS(1,IGET(269)),GRID1,IM,JM)
         ENDIF
!
!     TIME AVERAGED SURFACE MOMENTUM FLUX.
         IF (IGET(270).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME.EQ.'RAPR')THEN
            GRID1=SPVAL
            ID(1:25)=0
          ELSE
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J) = SFCVX(I,J)*RRNUM
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
            IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
             IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
            ELSE
             IFINCR     = 0
            endif
            ID(19)     = IFHR
            IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
               IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
          END IF
          CALL GRIBIT(IGET(270),LVLS(1,IGET(270)),GRID1,IM,JM)
         ENDIF
!     
!     ACCUMULATED SURFACE EVAPORATION
         IF (IGET(047).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SFCEVP(I,J)*1000.
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
	 IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
	     IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(047),LVLS(1,IGET(047)),GRID1,IM,JM)
         ENDIF
!     
!     ACCUMULATED POTENTIAL EVAPORATION
         IF (IGET(137).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=POTEVP(I,J)*1000.
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = NINT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
	 IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
	     IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(137),LVLS(1,IGET(137)),GRID1,IM,JM)
         ENDIF
!     
!     ROUGHNESS LENGTH.
      IF (IGET(044).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=Z0(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(044),LVLS(1,IGET(044)),GRID1,IM,JM)
      ENDIF
!     
!     FRICTION VELOCITY.
      IF (IGET(045).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=USTAR(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(045),LVLS(1,IGET(045)),GRID1,IM,JM)
      ENDIF
!     
!     SURFACE DRAG COEFFICIENT.
      IF (IGET(132).GT.0) THEN
         CALL CALDRG(EGRID1)
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=EGRID1(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(132),LVLS(1,IGET(132)),GRID1,IM,JM)
      ENDIF
!     
!     SURFACE U AND/OR V COMPONENT WIND STRESS
      IF ( (IGET(133).GT.0) .OR. (IGET(134).GT.0) ) THEN
         CALL CALTAU(EGRID1,EGRID2)
!     
!        SURFACE U COMPONENT WIND STRESS.
         IF (IGET(133).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=EGRID1(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(133),LVLS(1,IGET(133)),GRID1,IM,JM)
         ENDIF
!     
!        SURFACE V COMPONENT WIND STRESS
         IF (IGET(134).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=EGRID2(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(134),LVLS(1,IGET(134)),GRID1,IM,JM)
         ENDIF
      ENDIF
!     
!     GRAVITY U AND/OR V COMPONENT STRESS
      IF ( (IGET(315).GT.0) .OR. (IGET(316).GT.0) ) THEN
!     
!        GRAVITY U COMPONENT WIND STRESS.
         IF (IGET(315).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=GTAUX(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
	    ITSRFC     = INT(TSRFC)
            IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(315),LVLS(1,IGET(315)),GRID1,IM,JM)
         ENDIF
!     
!        SURFACE V COMPONENT WIND STRESS
         IF (IGET(316).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=GTAUY(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
	    ITSRFC     = INT(TSRFC)
            IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(316),LVLS(1,IGET(316)),GRID1,IM,JM)
         ENDIF
      ENDIF      
!     
!     INSTANTANEOUS SENSIBLE HEAT FLUX
      IF (IGET(154).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
	     IF(MODELNAME.EQ.'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. &
                MODELNAME.EQ.'RAPR')THEN
               GRID1(I,J)=TWBS(I,J)
             ELSE
               GRID1(I,J)=-1.*TWBS(I,J)
             END IF
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(154),LVLS(1,IGET(154)),GRID1,IM,JM)
      ENDIF
!     
!     INSTANTANEOUS LATENT HEAT FLUX
      IF (IGET(155).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
	     IF(MODELNAME.EQ.'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. &
                MODELNAME.EQ.'RAPR')THEN
               GRID1(I,J)=QWBS(I,J)
             ELSE
               GRID1(I,J)=-1.*QWBS(I,J)
             END IF
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(155),LVLS(1,IGET(155)),GRID1,IM,JM)
      ENDIF
!     
!     SURFACE EXCHANGE COEFF
      IF (IGET(169).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SFCEXC(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(169),LVLS(1,IGET(169)),GRID1,IM,JM)
      ENDIF
!     
!     GREEN VEG FRACTION
      IF (IGET(170).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=VEGFRC(I,J)*100.
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(170),LVLS(1,IGET(170)),GRID1,IM,JM)
      ENDIF
!     
!     INSTANTANEOUS GROUND HEAT FLUX
      IF (IGET(152).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=GRNFLX(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(152),LVLS(1,IGET(152)),GRID1,IM,JM)
      ENDIF
!    VEGETATION TYPE
      IF (IGET(218).GT.0) THEN
         DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = FLOAT(IVGTYP(I,J))
           ENDDO
         ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(218),LVLS(1,IGET(218)),GRID1,IM,JM)                                                          
      ENDIF
!
!    SOIL TYPE
      IF (IGET(219).GT.0) THEN
         DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = FLOAT(ISLTYP(I,J))
           ENDDO
         ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(219),LVLS(1,IGET(219)),GRID1,IM,JM)                                                          
      ENDIF
!    SLOPE TYPE
      IF (IGET(223).GT.0) THEN
         DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = FLOAT(ISLOPE(I,J))                                  
           ENDDO
         ENDDO
         ID(1:25) = 0
         ID(02)= 130
         CALL GRIBIT(IGET(223),LVLS(1,IGET(223)),GRID1,IM,JM)
                                                                                
      ENDIF
!      print*,'starting computing canopy conductance'
!
! CANOPY CONDUCTANCE
! ONLY OUTPUT NEW LSM FIELDS FOR NMM AND ARW BECAUSE RSM USES OLD SOIL TYPES
      IF (MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'NMM' .OR. MODELNAME.EQ.'RAPR')THEN
      IF (IGET(220).GT.0 .OR. IGET(234).GT.0               &
     & .OR. IGET(235).GT.0 .OR. IGET(236).GT.0             &
     & .OR. IGET(237).GT.0 .OR. IGET(238).GT.0             &
     & .OR. IGET(239).GT.0 .OR. IGET(240).GT.0             &
     & .OR. IGET(241).GT.0 .OR. IGET(254).GT.0 ) THEN
      print*,'starting computing canopy conductance'     
         DO J=JSTA,JEND
           DO I=1,IM
!             IF(abs(SM(I,J)-0.).lt.1.0E-5)THEN
             IF( (abs(SM(I,J)-0.)   .lt. 1.0E-5) .AND.     &
     &           (abs(SICE(I,J)-0.) .lt. 1.0E-5) ) THEN
              IF(CZMEAN(I,J).GT.1.E-6) THEN
               FACTRS=CZEN(I,J)/CZMEAN(I,J)
              ELSE
               FACTRS=0.0
              ENDIF
!              SOLAR=HBM2(I,J)*RSWIN(I,J)*FACTRS
              LLMH=NINT(LMH(I,J))
	      SOLAR=RSWIN(I,J)*FACTRS
              SFCTMP=T(I,J,LLMH)
              SFCQ=Q(I,J,LLMH)
              SFCPRS=PINT(I,J,LLMH+1)
!              IF(IVGTYP(I,J).EQ.0)PRINT*,'IVGTYP ZERO AT ',I,J
!     &        ,SM(I,J)
              IVG=IVGTYP(I,J)
!              IF(IVGTYP(I,J).EQ.0)IVG=7
!              CALL CANRES(SOLAR,SFCTMP,SFCQ,SFCPRS
!     &        ,SMC(I,J,1:NSOIL),GC(I,J),RC,IVG,ISLTYP(I,J))
!
              CALL CANRES(SOLAR,SFCTMP,SFCQ,SFCPRS            &
     &        ,SH2O(I,J,1:NSOIL),GC(I,J),RC,IVG,ISLTYP(I,J)   &
     &        ,RSMIN(I,J),NROOTS(I,J),SMCWLT(I,J),SMCREF(I,J) &
     &        ,RCS(I,J),RCQ(I,J),RCT(I,J),RCSOIL(I,J),SLDPTH)  
               IF(abs(SMCWLT(I,J)-0.5).lt.1.e-5)print*,       &
     &       'LARGE SMCWLT',i,j,SM(I,J),ISLTYP(I,J),SMCWLT(I,J)
             ELSE
              GC(I,J)=0.
              RSMIN(I,J)=0.
              NROOTS(I,J)=0
              SMCWLT(I,J)=0.
              SMCREF(I,J)=0.
              RCS(I,J)=0.
              RCQ(I,J)=0.
              RCT(I,J)=0.
              RCSOIL(I,J)=0.
             END IF
           ENDDO
         ENDDO
	 
         IF (IGET(220).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = GC(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(220),LVLS(1,IGET(220)),            &
              GRID1,IM,JM)                                                          
         ENDIF	 	     

         IF (IGET(234).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RSMIN(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(234),LVLS(1,IGET(234)),GRID1,IM,JM)                                                          
         ENDIF	
	 
         IF (IGET(235).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = FLOAT(NROOTS(I,J))
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(235),LVLS(1,IGET(235)),GRID1,IM,JM)                                                          
         ENDIF	

         IF (IGET(236).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = SMCWLT(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(236),LVLS(1,IGET(236)),GRID1,IM,JM)                                                          
         ENDIF	

         IF (IGET(237).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = SMCREF(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(237),LVLS(1,IGET(237)),GRID1,IM,JM)                                                          
         ENDIF	

         IF (IGET(238).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RCS(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(238),LVLS(1,IGET(238)),GRID1,IM,JM)                                                          
         ENDIF	

         IF (IGET(239).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RCT(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(239),LVLS(1,IGET(239)),GRID1,IM,JM)                                                          
         ENDIF	

         IF (IGET(240).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RCQ(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(240),LVLS(1,IGET(240)),GRID1,IM,JM)                                                          
         ENDIF	
         
         IF (IGET(241).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RCSOIL(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(241),LVLS(1,IGET(241)),GRID1,IM,JM)                                                          
         ENDIF	
	 
	 print*,'outputting leaf area index= ',XLAI
         IF (IGET(254).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = XLAI
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(254),LVLS(1,IGET(254)),GRID1,IM,JM)                                                          
         ENDIF

      ENDIF
      END IF
      IF(MODELNAME .EQ. 'GFS')THEN
! Outputting wilting point and field capacity for TIGGE
       IF(IGET(236).GT.0)THEN
        DO J=JSTA,JEND
          DO I=1,IM
!            GRID1(I,J) = smcwlt(i,j)
            IF(isltyp(i,j)/=0)THEN
              GRID1(I,J) = WLTSMC(isltyp(i,j))
            ELSE
              GRID1(I,J) = spval
            END IF
          ENDDO
        ENDDO
        ID(1:25) = 0
	ID(02)= 130
        CALL GRIBIT(IGET(236),LVLS(1,IGET(236)),  &
             GRID1,IM,JM)                                                          
       ENDIF
       
       IF(IGET(397).GT.0)THEN
        DO J=JSTA,JEND
          DO I=1,IM
!            GRID1(I,J) = fieldcapa(i,j)
            IF(isltyp(i,j)/=0)THEN
              GRID1(I,J) = REFSMC(isltyp(i,j))
            ELSE
              GRID1(I,J) = spval
            END IF
          ENDDO
        ENDDO
        ID(1:25) = 0
	ID(02)= 130
        CALL GRIBIT(IGET(397),LVLS(1,IGET(397)),   &
             GRID1,IM,JM)                                                          
       ENDIF
      END IF 
      IF(IGET(396).GT.0)THEN
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J) = suntime(i,j)
          ENDDO
        ENDDO
        ID(1:25) = 0
	ID(02)= 133
        CALL GRIBIT(IGET(396),LVLS(1,IGET(396)),   &
             GRID1,IM,JM)                                                          
       ENDIF    
!     
!     END OF ROUTINE
!     
!     
!       MODEL TOP REQUESTED BY CMAQ
      IF (IGET(282).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=PT
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(282),LVLS(1,IGET(282)),GRID1,IM,JM)
      ENDIF
!     
!       PRESSURE THICKNESS REQUESTED BY CMAQ
      IF (IGET(283).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=PDTOP
            ENDDO
            ENDDO
            ID(1:25) = 0
	    IF(ME == 0)THEN 
	     DO L=1,LM
	      IF(PMID(1,1,L).GE.(PDTOP+PT))EXIT
	     END DO
	     PRINT*,'hybrid boundary ',L
            END IF 
            CALL MPI_BCAST(L,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
	    ID(10)=1
	    ID(11)=L  
            CALL GRIBIT(IGET(283),LVLS(1,IGET(283)),GRID1,IM,JM)
      ENDIF
!      
!       SIGMA PRESSURE THICKNESS REQUESTED BY CMAQ
      IF (IGET(273).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=PD(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
	    ID(10)=L+1
	    ID(11)=LM
            CALL GRIBIT(IGET(273),LVLS(1,IGET(273)),GRID1,IM,JM)
      ENDIF
              
      RETURN
      END
