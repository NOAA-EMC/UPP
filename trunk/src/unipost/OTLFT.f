      SUBROUTINE OTLFT(PBND,TBND,QBND,SLINDX)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    OTLFT       COMPUTES LIFTED INDEX
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-03-10       
!     
! ABSTRACT:
!     THIS ROUTINE COMPUTES LIFTS A PARCEL SPECIFIED BY THE
!     PASSED PRESSURE, TEMPERATURE, AND SPECIFIC HUMIDITY TO
!     500MB AND THEN COMPUTES A LIFTED INDEX.  THIS LIFTED 
!     LIFTED INDEX IS THE DIFFERENCE BETWEEN THE LIFTED 
!     PARCEL'S TEMPERATURE AT 500MB AND THE AMBIENT 500MB
!     TEMPERATURE.
!   .     
!     
! PROGRAM HISTORY LOG:
!   93-03-10  RUSS TREADON - MODIFIED OTLIFT2 TO LIFT PARCELS
!                            SPECIFIED BY PASSED P, T, AND Q.
!   98-06-15  T BLACK      - CONVERSION FROM 1-D TO 2-D
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   02-06-17  MIKE BALDWIN - WRF VERSION
!   11-04-12  GEOFF MANIKIN - USE VIRTUAL TEMPERATURE
!     
! USAGE:    CALL OTLFT(PBND,TBND,QBND,SLINDX)
!   INPUT ARGUMENT LIST:
!     PBND     - PARCEL PRESSURE.
!     TBND     - PARCEL TEMPERATURE.
!     QBND     - PARCEL SPECIFIC HUMIDITY.
!
!   OUTPUT ARGUMENT LIST: 
!     SLINDX   - LIFTED INDEX.
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       NONE
!     LIBRARY:
!       COMMON   - CTLBLK
!                  LOOPS
!                  MASKS
!                  PHYS
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
!     
!     
      use vrbls2d
      use vrbls3d
      use lookup_mod
      use ctlblk_mod
      use params_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
!     SET LOCAL PARAMETERS.
       real,PARAMETER :: D8202=.820231E0 , H5E4=5.E4 , P500=50000.
       REAL TVP, ESATP, QSATP
       real,external::FPVSNEW

!     
!     DECLARE VARIABLES.
     real,dimension(IM,JM),intent(in) :: PBND,TBND,QBND
     real,dimension(IM,JM),intent(out) :: SLINDX

!
      integer              &
     & ITTB  (IM,JM),IQTB  (IM,JM),IPTB  (IM,JM),ITHTB (IM,JM)
      real TBT   (IM,JM),QBT   (IM,JM)                    &
     &,APEBT (IM,JM),TTHBT (IM,JM),TTH   (IM,JM),PP    (IM,JM)         &
     &,BQS00 (IM,JM),SQS00 (IM,JM),BQS10 (IM,JM),SQS10 (IM,JM)         &
     &,BQ    (IM,JM),SQ    (IM,JM),TQ    (IM,JM),QQ    (IM,JM)         &
     &,P00   (IM,JM),P10   (IM,JM),P01   (IM,JM),P11   (IM,JM)         &
     &,TPSP  (IM,JM),APESP (IM,JM),TTHES (IM,JM)                       &
     &,PSP   (IM,JM),THBT  (IM,JM),THESP (IM,JM)                       &
     &,P     (IM,JM),TP    (IM,JM),BTH   (IM,JM),STH   (IM,JM)         &
     &,BTHE00(IM,JM),STHE00(IM,JM),BTHE10(IM,JM),STHE10(IM,JM)         &
     &,T00   (IM,JM),T10   (IM,JM),T01   (IM,JM),T11   (IM,JM)         &
     &,PARTMP(IM,JM)
!     
       integer I,J,ITTBK,IQ,IT,IPTBK,ITH,IP
!     
!********************************************************************
!     START OTLFT HERE.
!     
!     ZERO LIFTED INDEX ARRAY.
!
      DO J=JSTA,JEND
      DO I=1,IM
        SLINDX(I,J)=D00
      ENDDO
      ENDDO
!
!--------------FIND EXNER IN BOUNDARY LAYER-----------------------------
!
      DO 130 J=JSTA,JEND
      DO 130 I=1,IM
      TBT(I,J)   = TBND(I,J)
      QBT(I,J)   = QBND(I,J)
      APEBT(I,J) = PBND(I,J)
      APEBT(I,J) = (H10E5/APEBT(I,J))**CAPA
  130 CONTINUE
!
!--------------SCALING POTENTIAL TEMPERATURE & TABLE INDEX--------------
!
      DO 140 J=JSTA,JEND
      DO 140 I=1,IM
      TTHBT(I,J)=TBT(I,J)*APEBT(I,J)
      TTH(I,J)=(TTHBT(I,J)-THL)*RDTH
      QQ(I,J)=TTH(I,J)-AINT(TTH(I,J))
      ITTB(I,J)=INT(TTH(I,J))+1
  140 CONTINUE
!
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
!
      DO 145 J=JSTA,JEND
      DO 145 I=1,IM
      IF(ITTB(I,J).LT.1)THEN
        ITTB(I,J)=1
        QQ(I,J)=D00
      ENDIF
      IF(ITTB(I,J).GE.JTB)THEN
        ITTB(I,J)=JTB-1
        QQ(I,J)=D00
      ENDIF
  145 CONTINUE
!
!--------------BASE AND SCALING FACTOR FOR SPEC. HUMIDITY---------------
!
      DO 150 J=JSTA,JEND
      DO 150 I=1,IM
      ITTBK=ITTB(I,J)
      BQS00(I,J)=QS0(ITTBK)
      SQS00(I,J)=SQS(ITTBK)
      BQS10(I,J)=QS0(ITTBK+1)
      SQS10(I,J)=SQS(ITTBK+1)
  150 CONTINUE
!
!--------------SCALING SPEC. HUMIDITY & TABLE INDEX---------------------
!
      DO 160 J=JSTA,JEND
      DO 160 I=1,IM
      BQ(I,J)=(BQS10(I,J)-BQS00(I,J))*QQ(I,J)+BQS00(I,J)
      SQ(I,J)=(SQS10(I,J)-SQS00(I,J))*QQ(I,J)+SQS00(I,J)
      TQ(I,J)=(QBT(I,J)-BQ(I,J))/SQ(I,J)*RDQ
      PP(I,J)=TQ(I,J)-AINT(TQ(I,J))
      IQTB(I,J)=INT(TQ(I,J))+1
  160 CONTINUE
!
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
!
      DO 165 J=JSTA,JEND
      DO 165 I=1,IM
        IF(IQTB(I,J).LT.1)THEN
        IQTB(I,J)=1
        PP(I,J)=D00
      ENDIF
      IF(IQTB(I,J).GE.ITB)THEN
        IQTB(I,J)=ITB-1
        PP(I,J)=D00
      ENDIF
 165  CONTINUE
!
!--------------SATURATION PRESSURE AT FOUR SURROUNDING TABLE PTS.-------
!
      DO 170 J=JSTA,JEND
      DO 170 I=1,IM
      IQ=IQTB(I,J)
      IT=ITTB(I,J)
      P00(I,J)=PTBL(IQ,IT)
      P10(I,J)=PTBL(IQ+1,IT)
      P01(I,J)=PTBL(IQ,IT+1)
      P11(I,J)=PTBL(IQ+1,IT+1)
  170 CONTINUE
!
!--------------SATURATION POINT VARIABLES AT THE BOTTOM-----------------
!
      DO 180 J=JSTA,JEND
      DO 180 I=1,IM
      TPSP(I,J)=P00(I,J)+(P10(I,J)-P00(I,J))*PP(I,J)    &
                        +(P01(I,J)-P00(I,J))*QQ(I,J)    &
            +(P00(I,J)-P10(I,J)-P01(I,J)+P11(I,J))*PP(I,J)*QQ(I,J)
      IF(TPSP(I,J).LE.D00)TPSP(I,J)=H10E5
      APESP(I,J)=(H10E5/TPSP(I,J))**CAPA
      TTHES(I,J)=TTHBT(I,J)*EXP(ELOCP*QBT(I,J)*APESP(I,J)/TTHBT(I,J))
  180 CONTINUE
!
!-----------------------------------------------------------------------
!
      DO J=JSTA,JEND
      DO I=1,IM
        PSP(I,J)=TPSP(I,J)
        THBT(I,J)=TTHBT(I,J)
        THESP(I,J)=TTHES(I,J)
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
 190  CONTINUE
!
!--------------SCALING PRESSURE & TT TABLE INDEX------------------------
!
      DO 210 J=JSTA,JEND
      DO 210 I=1,IM
      P (I,J)=H5E4
      TP(I,J)=(P(I,J)-PL)*RDP
      QQ(I,J)=TP(I,J)-AINT(TP(I,J))
      IPTB(I,J)=INT(TP(I,J))+1
  210 CONTINUE
!
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
!
      DO 215 J=JSTA,JEND
      DO 215 I=1,IM
      IF(IPTB(I,J).LT.1)THEN
        IPTB(I,J)=1
        QQ(I,J)=D00
      ENDIF
      IF(IPTB(I,J).GE.ITB)THEN
        IPTB(I,J)=ITB-1
        QQ(I,J)=D00
      ENDIF
 215  CONTINUE
!
!--------------BASE AND SCALING FACTOR FOR THE--------------------------
!
      DO 220 J=JSTA,JEND
      DO 220 I=1,IM
      IPTBK=IPTB(I,J)
      BTHE00(I,J)=THE0(IPTBK)
      STHE00(I,J)=STHE(IPTBK)
      BTHE10(I,J)=THE0(IPTBK+1)
      STHE10(I,J)=STHE(IPTBK+1)
  220 CONTINUE
!
!--------------SCALING THE & TT TABLE INDEX-----------------------------
!
      DO 230 J=JSTA,JEND
      DO 230 I=1,IM
      BTH(I,J)=(BTHE10(I,J)-BTHE00(I,J))*QQ(I,J)+BTHE00(I,J)
      STH(I,J)=(STHE10(I,J)-STHE00(I,J))*QQ(I,J)+STHE00(I,J)
      TTH(I,J)=(THESP(I,J)-BTH(I,J))/STH(I,J)*RDTHE
      PP(I,J)=TTH(I,J)-AINT(TTH(I,J))
      ITHTB(I,J)=INT(TTH(I,J))+1
  230 CONTINUE
!
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
!
      DO 235 J=JSTA,JEND
      DO 235 I=1,IM
      IF(ITHTB(I,J).LT.1)THEN
        ITHTB(I,J)=1
        PP(I,J)=D00
      ENDIF
      IF(ITHTB(I,J).GE.JTB)THEN
        ITHTB(I,J)=JTB-1
        PP(I,J)=D00
      ENDIF
 235  CONTINUE
!
!--------------TEMPERATURE AT FOUR SURROUNDING TT TABLE PTS.------------
!
      DO 240 J=JSTA,JEND
      DO 240 I=1,IM
      ITH=ITHTB(I,J)
      IP=IPTB(I,J)
      T00(I,J)=TTBL(ITH,IP)
      T10(I,J)=TTBL(ITH+1,IP)
      T01(I,J)=TTBL(ITH,IP+1)
      T11(I,J)=TTBL(ITH+1,IP+1)
  240 CONTINUE
!
!--------------PARCEL TEMPERATURE AT 500MB----------------------------
!
      DO 300 J=JSTA,JEND
      DO 300 I=1,IM
      IF(TPSP(I,J).GE.H5E4)THEN
        PARTMP(I,J)=(T00(I,J)+(T10(I,J)-T00(I,J))*PP(I,J)     &
                             +(T01(I,J)-T00(I,J))*QQ(I,J)     &
               +(T00(I,J)-T10(I,J)-T01(I,J)+T11(I,J))*PP(I,J)*QQ(I,J))
      ELSE
        PARTMP(I,J)=TBT(I,J)*APEBT(I,J)*D8202
      ENDIF
!
!--------------LIFTED INDEX---------------------------------------------
!
! GSM  THE PARCEL TEMPERATURE AT 500 MB HAS BEEN COMPUTED, AND WE
!       FIND THE MIXING RATIO AT THAT LEVEL WHICH WILL BE THE SATURATION
!       VALUE SINCE WE'RE FOLLOWING A MOIST ADIABAT.    NOTE THAT THE
!       AMBIENT 500 MB SHOULD PROBABLY BE VIRTUALIZED, BUT THE IMPACT
!       OF MOISTURE AT THAT LEVEL IS QUITE SMALL
       ESATP=FPVSNEW(PARTMP(I,J))
       QSATP=EPS*ESATP/(P500-ESATP*ONEPS)
       TVP=PARTMP(I,J)*(1+0.608*QSATP)
       SLINDX(I,J)=T500(I,J)-TVP
  300 CONTINUE
!     
!     END OF ROUTINE.
      RETURN
      END
