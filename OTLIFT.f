      SUBROUTINE OTLIFT(SLINDX)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    OTLIFT      COMPUTES SFC TO 500MB LIFTED INDEX
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-03-10       
!     
! ABSTRACT:
!     THIS ROUTINE COMPUTES A SURFACE TO 500MB LIFTED INDEX.
!     THE LIFTED PARCEL IS FROM THE FIRST ATMOSPHERIC ETA 
!     LAYER (IE, THE ETA LAYER CLOSEST TO THE MODEL GROUND).
!     THE LIFTED INDEX IS THE DIFFERENCE BETWEEN THIS PARCEL'S
!     TEMPERATURE AT 500MB AND THE AMBIENT 500MB TEMPERATURE.
!   .     
!     
! PROGRAM HISTORY LOG:
!   ??-??-??  ??? - SUBROUTINE OTLIFT IN ETA MODEL.
!   93-03-10  RUSS TREADON - ADAPTED OTLIFT FOR USE WITH NEW POST. 
!   98-06-18  T BLACK      - CONVERSION FROM 1-D TO 2-D
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-06-11  MIKE BALDWIN - WRF VERSION
!     
! USAGE:    CALL OTLIFT(SLINDX)
!   INPUT ARGUMENT LIST:
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
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
!     
      use vrbls3d
      use vrbls2d
      use masks
      use lookup_mod
      use ctlblk_mod
!
      implicit none
!
!     SET LOCAL PARAMETERS.
      real,PARAMETER :: D00=0.E0,H10E5=100000.E0,H5E4=5.E4,CAPA=0.28589641    &
     &, D8202=.820231E0
      real,PARAMETER :: ELIVW=2.72E6,CP=1004.6E0,ELOCP=ELIVW/CP
!     
!     DECLARE VARIABLES.
      real,intent(out) :: SLINDX(IM,JM)
      integer & 
         ITTB  (IM,JM),IQTB  (IM,JM),IPTB  (IM,JM),ITHTB (IM,JM)
      real  TBT   (IM,JM),QBT   (IM,JM)            &
     &,APEBT (IM,JM),TTHBT (IM,JM),TTH   (IM,JM),PP    (IM,JM)    &
     &,BQS00 (IM,JM),SQS00 (IM,JM),BQS10 (IM,JM),SQS10 (IM,JM)    &
     &,BQ    (IM,JM),SQ    (IM,JM),TQ    (IM,JM),QQ    (IM,JM)    &
     &,P00   (IM,JM),P10   (IM,JM),P01   (IM,JM),P11   (IM,JM)    &
     &,TPSP  (IM,JM),APESP (IM,JM),TTHES (IM,JM)                  &
     &,PSP   (IM,JM),THBT  (IM,JM),THESP (IM,JM)                  &
     &,P     (IM,JM),TP    (IM,JM),BTH   (IM,JM),STH   (IM,JM)    &
     &,BTHE00(IM,JM),STHE00(IM,JM),BTHE10(IM,JM),STHE10(IM,JM)    &
     &,T00   (IM,JM),T10   (IM,JM),T01   (IM,JM),T11   (IM,JM)    &
     &,PARTMP(IM,JM)
!
       integer I,J,LBTM,ITTBK,IQ,IT,IPTBK,ITH,IP
!     
!***********************************************************************
!     START OTLIFT HERE
!     
!     INTIALIZE LIFTED INDEX ARRAY TO ZERO.
      DO J=JSTA,JEND
      DO I=1,IM
        SLINDX(I,J)=D00
      ENDDO
      ENDDO
!--------------FIND EXNER AT LOWEST LEVEL-------------------------------
      DO J=JSTA,JEND
      DO I=1,IM
        LBTM=NINT(LMH(I,J))
        TBT(I,J)=T(I,J,LBTM)
        QBT(I,J)=Q(I,J,LBTM)
        APEBT(I,J)=PMID(I,J,LBTM)
        APEBT(I,J)=(H10E5/APEBT(I,J))**CAPA
      ENDDO
      ENDDO
!--------------SCALING POTENTIAL TEMPERATURE & TABLE INDEX--------------
      DO J=JSTA,JEND
      DO I=1,IM
        TTHBT(I,J)=TBT(I,J)*APEBT(I,J)
        TTH(I,J)=(TTHBT(I,J)-THL)*RDTH
        QQ(I,J)=TTH(I,J)-AINT(TTH(I,J))
        ITTB(I,J)=INT(TTH(I,J))+1
      ENDDO
      ENDDO
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
      DO J=JSTA,JEND
      DO I=1,IM
        IF(ITTB(I,J).LT.1)THEN
          ITTB(I,J)=1
          QQ(I,J)=D00
        ENDIF
          IF(ITTB(I,J).GE.JTB)THEN
          ITTB(I,J)=JTB-1
          QQ(I,J)=D00
        ENDIF
      ENDDO
      ENDDO
!--------------BASE AND SCALING FACTOR FOR SPEC. HUMIDITY---------------
      DO J=JSTA,JEND
      DO I=1,IM
        ITTBK=ITTB(I,J)
        BQS00(I,J)=QS0(ITTBK)
        SQS00(I,J)=SQS(ITTBK)
        BQS10(I,J)=QS0(ITTBK+1)
        SQS10(I,J)=SQS(ITTBK+1)
      ENDDO
      ENDDO
!--------------SCALING SPEC. HUMIDITY & TABLE INDEX---------------------
      DO J=JSTA,JEND
      DO I=1,IM
        BQ(I,J)=(BQS10(I,J)-BQS00(I,J))*QQ(I,J)+BQS00(I,J)
        SQ(I,J)=(SQS10(I,J)-SQS00(I,J))*QQ(I,J)+SQS00(I,J)
        TQ(I,J)=(QBT(I,J)-BQ(I,J))/SQ(I,J)*RDQ
        PP(I,J)=TQ(I,J)-AINT(TQ(I,J))
        IQTB(I,J)=INT(TQ(I,J))+1
      ENDDO
      ENDDO
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
      DO J=JSTA,JEND
      DO I=1,IM
        IF(IQTB(I,J).LT.1)THEN
          IQTB(I,J)=1
          PP(I,J)=D00
        ENDIF
        IF(IQTB(I,J).GE.ITB)THEN
          IQTB(I,J)=ITB-1
          PP(I,J)=D00
        ENDIF
      ENDDO
      ENDDO
!--------------SATURATION PRESSURE AT FOUR SURROUNDING TABLE PTS.-------
      DO J=JSTA,JEND
      DO I=1,IM
        IQ=IQTB(I,J)
        IT=ITTB(I,J)
        P00(I,J)=PTBL(IQ,IT)
        P10(I,J)=PTBL(IQ+1,IT)
        P01(I,J)=PTBL(IQ,IT+1)
        P11(I,J)=PTBL(IQ+1,IT+1)
      ENDDO
      ENDDO
!--------------SATURATION POINT VARIABLES AT THE BOTTOM-----------------
      DO J=JSTA,JEND
      DO I=1,IM
        TPSP(I,J)=P00(I,J)+(P10(I,J)-P00(I,J))*PP(I,J)       &
                          +(P01(I,J)-P00(I,J))*QQ(I,J)       &
              +(P00(I,J)-P10(I,J)-P01(I,J)+P11(I,J))*PP(I,J)*QQ(I,J)
        IF(TPSP(I,J).LE.D00)TPSP(I,J)=H10E5
        APESP(I,J)=(H10E5/TPSP(I,J))**CAPA
        TTHES(I,J)=TTHBT(I,J)*EXP(ELOCP*QBT(I,J)*APESP(I,J)/TTHBT(I,J))
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
      DO J=JSTA,JEND
      DO I=1,IM
        PSP(I,J)=TPSP(I,J)
        THBT(I,J)=TTHBT(I,J)
        THESP(I,J)=TTHES(I,J)
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
 190  CONTINUE
!--------------SCALING PRESSURE & TT TABLE INDEX------------------------
      DO J=JSTA,JEND
      DO I=1,IM
        P (I,J)=H5E4
        TP(I,J)=(P(I,J)-PL)*RDP
        QQ(I,J)=TP(I,J)-AINT(TP(I,J))
        IPTB(I,J)=INT(TP(I,J))+1
      ENDDO
      ENDDO
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
      DO J=JSTA,JEND
      DO I=1,IM
        IF(IPTB(I,J).LT.1)THEN
          IPTB(I,J)=1
          QQ(I,J)=D00
        ENDIF
        IF(IPTB(I,J).GE.ITB)THEN
          IPTB(I,J)=ITB-1
          QQ(I,J)=D00
        ENDIF
      ENDDO
      ENDDO
!--------------BASE AND SCALING FACTOR FOR THE--------------------------
      DO J=JSTA,JEND
      DO I=1,IM
        IPTBK=IPTB(I,J)
        BTHE00(I,J)=THE0(IPTBK)
        STHE00(I,J)=STHE(IPTBK)
        BTHE10(I,J)=THE0(IPTBK+1)
        STHE10(I,J)=STHE(IPTBK+1)
      ENDDO
      ENDDO
!--------------SCALING THE & TT TABLE INDEX-----------------------------
      DO J=JSTA,JEND
      DO I=1,IM
        BTH(I,J)=(BTHE10(I,J)-BTHE00(I,J))*QQ(I,J)+BTHE00(I,J)
        STH(I,J)=(STHE10(I,J)-STHE00(I,J))*QQ(I,J)+STHE00(I,J)
        TTH(I,J)=(THESP(I,J)-BTH(I,J))/STH(I,J)*RDTHE
        PP(I,J)=TTH(I,J)-AINT(TTH(I,J))
        ITHTB(I,J)=INT(TTH(I,J))+1
      ENDDO
      ENDDO
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
      DO J=JSTA,JEND
      DO I=1,IM
        IF(ITHTB(I,J).LT.1)THEN
          ITHTB(I,J)=1
          PP(I,J)=D00
        ENDIF
        IF(ITHTB(I,J).GE.JTB)THEN
          ITHTB(I,J)=JTB-1
          PP(I,J)=D00
        ENDIF
      ENDDO
      ENDDO
!--------------TEMPERATURE AT FOUR SURROUNDING TT TABLE PTS.------------
      DO J=JSTA,JEND
      DO I=1,IM
        ITH=ITHTB(I,J)
        IP=IPTB(I,J)
        T00(I,J)=TTBL(ITH,IP)
        T10(I,J)=TTBL(ITH+1,IP)
        T01(I,J)=TTBL(ITH,IP+1)
        T11(I,J)=TTBL(ITH+1,IP+1)
      ENDDO
      ENDDO
!--------------PARCEL TEMPERATURE AT 500MB----------------------------
      DO J=JSTA,JEND
      DO I=1,IM
        IF(TPSP(I,J).GE.H5E4)THEN
          PARTMP(I,J)=(T00(I,J)+(T10(I,J)-T00(I,J))*PP(I,J)    &
                               +(T01(I,J)-T00(I,J))*QQ(I,J)    &
              +(T00(I,J)-T10(I,J)-T01(I,J)+T11(I,J))*PP(I,J)*QQ(I,J))
        ELSE
          PARTMP(I,J)=TBT(I,J)*APEBT(I,J)*D8202
        ENDIF
!--------------LIFTED INDEX---------------------------------------------
        SLINDX(I,J)=T500(I,J)-PARTMP(I,J)
      ENDDO
      ENDDO
!       write(*,*) ' in otlift t500 partmp ',t500(1,1),partmp(1,1)
!       write(*,*) ' in otlift tbt ',tbt(1,1)
!
      RETURN
      END
