      SUBROUTINE MISCLN
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    MISCLN      POSTS MISCELLANEOUS FIELDS
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-20
!     
! ABSTRACT:
!     THIS ROUTINE HAS BECOME THE CATCH-ALL FOR MISCELLANEOUS
!     OUTPUT FIELDS POSTED BY THE ETA POST PROCESSOR.  
!     CURRENTLY THIS ROUTINE POSTS THE FOLLOWING FIELDS:
!        (1) TROPOPAUSE LEVEL Z,P, T, U, V, AND VERTICAL WIND SHEAR,
!        (2) MAX WIND LEVEL Z, P, U, AND V,
!        (3) FD LEVEL T, Q, U, AND V,
!        (4) FREEZING LEVEL Z AND RH,
!        (5) CONSTANT MASS (BOUNDARY) FIELDS,
!        (6) LFM LOOK-ALIKE FIELDS, AND
!        (7) NGM LOOK-ALIKE FIELDS.
!
!   .     
!     
! PROGRAM HISTORY LOG:
!   92-12-20  RUSS TREADON
!   93-06-19  RUSS TREADON - ADDED TYPE 2 CAPE POSTING.
!   94-11-07  MIKE BALDWIN - ADDED HELICITY POSTING.
!   96-03-26  MIKE BALDWIN - CHANGE ETA BOUNDARY LAYER LABELS FOR GRIB
!   96-11-19  MIKE BALDWIN - BACK OUT PREVIOUS CHANGE 
!   97-04-25  MIKE BALDWIN - CHANGE ETA BOUNDARY LAYER LABELS FOR GRIB
!   97-04-29  GEOFF MANIKIN - ADDED TROPOPAUSE HEIGHT AND
!                             MAX WIND LEVEL FIELDS
!   98-06-15  T BLACK       - CONVERSION FROM 1-D TO 2-D
!   98-07-17  MIKE BALDWIN - REMOVED LABL84
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   02-04-23  MIKE BALDWIN - WRF VERSION
!     
! USAGE:    CALL MISCLN
!   INPUT ARGUMENT LIST:
!
!   OUTPUT ARGUMENT LIST: 
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       TRPAUS  - COMPUTE TROPOPAUSE LEVEL FIELDS.
!       CALMXW  - COMPUTE MAX WIND LEVEL FIELDS.
!       SCLFLD  - SCALE ARRAY ELEMENTS BY CONSTANT.
!       GRIBIT  - OUTPUT FIELD TO GRIB FILE.
!       CALPOT  - CALCULATE POTENTIAL TEMPERATURE.
!       FDLVL   - COMPUTE FD LEVEL DATA (AGL OR MSL).
!       FRZLVL  - COMPUTE FREEZING LEVEL DATA.
!       BOUND   - BOUND ARRAY ELEMENTS BETWEEN MINIMUM AND MAXIMUM VALUES.
!       BNDLYR  - COMPUTE BOUNDARY LAYER FIELDS.
!       CALDWP  - CALCULATE DEWPOINT TEMPERATURE.
!       OTLFT   - COMPUTE LIFTED INDEX AT 500MB.
!       CALLCL  - COMPUTE LCL DATA.
!       LFMFLD  - COMPUTE LFM LOOK-ALIKE FIELDS.
!       NGMFLD  - COMPUTE NGM LOOK-ALIKE FIELDS.
!       CALTHTE - COMPUTE THETA-E.
!       CALHEL  - COMPUTE HELICITY AND STORM MOTION.
!
!     LIBRARY:
!       COMMON - RQSTFLD
!                CTLBLK
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
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
!     SET LOCAL PARAMETERS.  MAKE SURE NFD AND NBND AGREE
!     WITH THE VALUES SET IN SUBROUTINES FDLVL AND BNDLYR,
!     RESPECTIVELY.
      integer,PARAMETER :: NFD=15,NBND=6
      real,PARAMETER :: C2K=273.15
      real,parameter:: con_rd      =2.8705e+2 ! gas constant air    (J/kg/K)
      real,parameter:: con_rv      =4.6150e+2 ! gas constant H2O
      real,parameter:: con_eps     =con_rd/con_rv
      real,parameter:: con_epsm1   =con_rd/con_rv-1
!     
!     DECLARE VARIABLES.
!     
      LOGICAL NORTH
      LOGICAL FIELD1,FIELD2
      LOGICAL DONE(IM,JSTA:JEND),DONE1(IM,JSTA:JEND)
      INTEGER LVLBND(IM,JM,NBND),LB2(IM,JM),LPBL(IM,JM)
      REAL P1D(IM,JM),T1D(IM,JM),Q1D(IM,JM),U1D(IM,JM),V1D(IM,JM)
      REAL SHR1D(IM,JM),Z1D(IM,JM),RH1D(IM,JM)
      REAL OMGBND(IM,JM,NBND),PWTBND(IM,JM,NBND)
      REAL QCNVBND(IM,JM,NBND)
      REAL PBND(IM,JM,NBND),TBND(IM,JM,NBND),QBND(IM,JM,NBND)
      REAL UBND(IM,JM,NBND),VBND(IM,JM,NBND),RHBND(IM,JM,NBND)
      REAL WBND(IM,JM,NBND)
      REAL T78483(IM,JM),T89671(IM,JM),P78483(IM,JM),P89671(IM,JM)
      REAL PKU1
      REAL QM8510(IM,JM),RH4710(IM,JM),RH8498(IM,JM)
      REAL RH4796(IM,JM),RH1847(IM,JM),UST(IM,JM),VST(IM,JM)
      REAL RH3310(IM,JM),RH6610(IM,JM),RH3366(IM,JM),PW3310(IM,JM)
      REAL RH4410(IM,JM),RH7294(IM,JM),RH4472(IM,JM)
      REAL HTFD(NFD),T7D(IM,JM,NFD),Q7D(IM,JM,NFD),U7D(IM,JM,NFD),V6D(IM,JM,NFD) &
         ,P7D(IM,JM,NFD)      
      REAL PETABND(NBND),SIGBND(NBND),HELI(IM,JM)
      REAL EGRID1(IM,JM),EGRID2(IM,JM),EGRID3(IM,JM)
      REAL EGRID4(IM,JM),EGRID5(IM,JM)
      REAL GRID1(IM,JM),GRID2(IM,JM)
      REAL MAXWP(IM,JM),MAXWZ(IM,JM),MAXWU(IM,JM), MAXWV(IM,JM)    &
     &  ,MAXWT(IM,JM)      
      REAL GUST(IM,JM) 
!     
      integer I,J,L,ITYPE,ISVALUE,LBND,ILVL,IFD,ITYPEFDLVL(NFD)
      real DPBND,PKL1,FAC1,FAC2,PL,TL,QL,QSAT,RHL,TVRL,TVRBLO,     &
           ES1,ES2,QS1,QS2,RH1,RH2,ZSF,DEPTH
      real,external :: fpvsnew
!     
!     SET FD LEVEL HEIGHTS IN GEOPOTENTAL METERS.
      DATA HTFD  / 30.E0,50.E0,80.E0,100.E0,305.E0,457.E0,610.E0,   &
	   914.E0,1524.E0,1829.E0,2134.E0,2743.E0,3658.E0,4572.E0,6000.E0/
!     
!     SET MIDPOINT "SIGMA" VALUES FOR ETA BOUNDARY LAYERS.
      DATA SIGBND / 0.985,0.955,0.925,0.895,0.865,0.835 /
      DATA PETABND / 15.,45.,75.,105.,135.,165. /
!     
!****************************************************************************
!     START MISCLN HERE.
!     
!        HELICITY AND STORM MOTION.
       IF (IGET(162).GT.0.OR.IGET(163).GT.0.OR.IGET(164).GT.0) THEN
        IF (IGET(162).GT.0) THEN
          IF(LVLS(1,IGET(162)).GT.0)DEPTH=3000.0
          CALL CALHEL(DEPTH,UST,VST,HELI)
          IF (IGET(162).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=HELI(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            ID(10)   = 30
            ID(11)   = 0
            CALL GRIBIT(IGET(162),LVLS(1,IGET(162)),GRID1,IM,JM)
          ENDIF

          IF(LVLS(2,IGET(162)).GT.0)DEPTH=1000.0
          CALL CALHEL(DEPTH,UST,VST,HELI)
          IF (IGET(162).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=HELI(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            ID(10)   = 10
            ID(11)   = 0
            CALL GRIBIT(IGET(162),LVLS(1,IGET(162)),GRID1,IM,JM)
          ENDIF
         ENDIF



         IF (IGET(163).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=UST(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            ID(10)   = 60
            ID(11)   = 0 
            CALL GRIBIT(IGET(163),LVLS(1,IGET(163)),GRID1,IM,JM)
         ENDIF
         IF (IGET(164).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=VST(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            ID(10)   = 60
            ID(11)   = 0 
            CALL GRIBIT(IGET(164),LVLS(1,IGET(164)),GRID1,IM,JM)
         ENDIF
       ENDIF
!     
!
!
!     ***BLOCK 1:  TROPOPAUSE P, Z, T, U, V, AND WIND SHEAR.
!    
      IF ( (IGET(054).GT.0).OR.(IGET(055).GT.0).OR.       &
           (IGET(056).GT.0).OR.(IGET(057).GT.0).OR.       &
           (IGET(177).GT.0).OR.                           &
           (IGET(058).GT.0).OR.(IGET(108).GT.0) ) THEN
! Chuang: Use GFS algorithm per Iredell's and DiMego's decision on unification
          DO J=JSTA,JEND
           DO I=1,IM
! INPUT
            CALL TPAUSE(LM,PMID(I,J,1:LM),UH(I,J,1:LM)    & 
! INPUT
      	      ,VH(I,J,1:LM),T(I,J,1:LM),ZMID(I,J,1:LM)    &
! OUTPUT
              ,P1D(I,J),U1D(I,J),V1D(I,J),T1D(I,J)        &
! OUTPUT
              ,Z1D(I,J),SHR1D(I,J))	               ! OUTPUT
           END DO
	  END DO 	  
!
!        TROPOPAUSE PRESSURE.
         IF (IGET(054).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=P1D(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(054),LVLS(1,IGET(054)),GRID1,IM,JM)
         ENDIF

!        ICAO HEIGHT OF TROPOPAUSE
         IF (IGET(399).GT.0) THEN
	    CALL ICAOHEIGHT(P1D,  & !input
                         GRID1)   ! output  
!            print*,'sample TROPOPAUSE ICAO HEIGHTS',GRID1(im/2,(jsta+jend)/2)
            ID(1:25) = 0
            CALL GRIBIT(IGET(399),LVLS(1,IGET(399)),GRID1,IM,JM)
          ENDIF

!        TROPOPAUSE HEIGHT.
         IF (IGET(177).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=Z1D(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(177),LVLS(1,IGET(177)),GRID1,IM,JM)
         ENDIF
!
!        TROPOPAUSE TEMPERATURE.
         IF (IGET(055).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=T1D(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(055),LVLS(1,IGET(055)),GRID1,IM,JM)
         ENDIF
!
!        TROPOPAUSE POTENTIAL TEMPERATURE.
         IF (IGET(108).GT.0) THEN
            CALL CALPOT(P1D,T1D,EGRID1)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID1(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(108),LVLS(1,IGET(108)),GRID1,IM,JM)
         ENDIF
!     
!        TROPOPAUSE U WIND AND/OR V WIND.
         IF ((IGET(056).GT.0).OR.(IGET(057).GT.0)) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=U1D(I,J)
                 GRID2(I,J)=V1D(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            IF (IGET(056).GT.0) CALL GRIBIT(IGET(056),      &
                 LVLS(1,IGET(056)),GRID1,IM,JM)
            ID(1:25) = 0
            IF (IGET(057).GT.0) CALL GRIBIT(IGET(057),      &
                 LVLS(1,IGET(057)),GRID2,IM,JM)
         ENDIF
!
!        TROPOPAUSE WIND SHEAR.
         IF (IGET(058).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=SHR1D(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(058),LVLS(1,IGET(058)),GRID1,IM,JM)
         ENDIF
      ENDIF
!
!
!
!     ***BLOCK 2:  MAX WIND LEVEL  P, Z, U, AND V
!
!        MAX WIND LEVEL CALCULATIONS
         IF ((IGET(173).GT.0) .OR. (IGET(174).GT.0) .OR.    &
            (IGET(175).GT.0) .OR. (IGET(176).GT.0)) THEN
!            CALL CALMXW(MAXWP,MAXWZ,MAXWU,MAXWV,MAXWT)
! Chuang: Use GFS algorithm per Iredell's and DiMego's decision on unification
          DO J=JSTA,JEND
           DO I=1,IM
! INPUT
            CALL MXWIND(LM,PMID(I,J,1:LM),UH(I,J,1:LM)      &
! INPUT
      	      ,VH(I,J,1:LM),T(I,J,1:LM),ZMID(I,J,1:LM)      &
! OUTPUT
              ,MAXWP(I,J),MAXWU(I,J),MAXWV(I,J)             &
! OUTPUT
              ,MAXWT(I,J),MAXWZ(I,J))	               
           END DO
	  END DO 
         ENDIF
!        PRESSURE OF MAX WIND LEVEL
         IF (IGET(173).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=MAXWP(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(173),LVLS(1,IGET(173)),GRID1,IM,JM)
          ENDIF
!        ICAO HEIGHT OF MAX WIND LEVEL
         IF (IGET(398).GT.0) THEN
	    CALL ICAOHEIGHT(MAXWP,  & !input
                         GRID1)   ! output  
!            print*,'sample MAX WIND ICAO HEIGHTS',GRID1(im/2,(jsta+jend)/2)
            ID(1:25) = 0
            CALL GRIBIT(IGET(398),LVLS(1,IGET(398)),GRID1,IM,JM)
          ENDIF
!        HEIGHT OF MAX WIND LEVEL
         IF (IGET(174).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=MAXWZ(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(174),LVLS(1,IGET(174)),GRID1,IM,JM)
          ENDIF

!        MAX WIND LEVEL U WIND AND/OR V WIND.
         IF ((IGET(175).GT.0).OR.(IGET(176).GT.0)) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=MAXWU(I,J)
                 GRID2(I,J)=MAXWV(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            IF (IGET(175).GT.0) CALL GRIBIT(IGET(175),      &
                 LVLS(1,IGET(175)),GRID1,IM,JM)
            ID(1:25) = 0
            IF (IGET(176).GT.0) CALL GRIBIT(IGET(176),      &
                 LVLS(1,IGET(176)),GRID2,IM,JM)
         ENDIF
!        TEMPERATURE OF MAX WIND LEVEL
         IF (IGET(314).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=MAXWT(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(314),LVLS(1,IGET(314)),GRID1,IM,JM)
          ENDIF
!
!
!
!     ***BLOCK 3:  FD LEVEL T, Q, U, AND V.
!     
      IF ( (IGET(059).GT.0).OR.(IGET(060).GT.0).OR.        &
           (IGET(061).GT.0).OR.(IGET(451).GT.0) ) THEN
!
!     DETERMINE WHETHER TO DO MSL OR AGL FD LEVELS
!
         ITYPEFDLVL=1
         DO IFD = 1,NFD
           IF (IGET(059).GT.0) THEN
            IF (LVLS(IFD,IGET(059)).GT.1) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(060).GT.0) THEN
            IF (LVLS(IFD,IGET(060)).GT.1) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(061).GT.0) THEN
            IF (LVLS(IFD,IGET(061)).GT.1) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(451).GT.0) THEN
	    IF (LVLS(IFD,IGET(451)).GT.1) ITYPEFDLVL(IFD)=2
	   ENDIF
         ENDDO
	write(6,*) 'call FDLVL with ITYPEFDLVL: ', ITYPEFDLVL
         CALL FDLVL(NFD,ITYPEFDLVL,HTFD,T7D,Q7D,U7D,V6D,P7D)
!     
         DO 10 IFD = 1,NFD
            ID(1:25) = 0
            ISVALUE = NINT(HTFD(IFD))
            ID(11) = ISVALUE
            if(ITYPEFDLVL(IFD)==2)ID(9)=105
!
!           FD LEVEL TEMPERATURE.
            IF (IGET(059).GT.0) THEN
              IF (LVLS(IFD,IGET(059)).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=T7D(I,J,IFD)
               ENDDO
               ENDDO
               CALL GRIBIT(IGET(059),LVLS(IFD,IGET(059)),GRID1,IM,JM)
              ENDIF
            ENDIF
!
!           FD LEVEL SPEC HUMIDITY.
            IF (IGET(451).GT.0) THEN
	      IF (LVLS(IFD,IGET(451)).GT.0) THEN
	       DO J=JSTA,JEND
	       DO I=1,IM
	         GRID1(I,J)=Q7D(I,J,IFD)
               ENDDO
               ENDDO
	       CALL GRIBIT(IGET(451),LVLS(IFD,IGET(451)),GRID1,IM,JM)
	      ENDIF
	    ENDIF
!
!           FD LEVEL PRESSURE
            IF (IGET(482).GT.0) THEN
	      IF (LVLS(IFD,IGET(482)).GT.0) THEN
	       DO J=JSTA,JEND
	       DO I=1,IM
	         GRID1(I,J)=P7D(I,J,IFD)
               ENDDO
               ENDDO
	       CALL GRIBIT(IGET(482),LVLS(IFD,IGET(482)),GRID1,IM,JM)
	      ENDIF
	    ENDIF
!
!           FD LEVEL U WIND AND/OR V WIND.
            IF ((IGET(060).GT.0).OR.(IGET(061).GT.0)) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=U7D(I,J,IFD)
                 GRID2(I,J)=V6D(I,J,IFD)
               ENDDO
               ENDDO
               IF (IGET(060).GT.0) THEN
                 IF (LVLS(IFD,IGET(060)).GT.0) CALL GRIBIT(      &
                    IGET(060),LVLS(IFD,IGET(060)),GRID1,IM,JM)
               ENDIF
               IF (IGET(061).GT.0) THEN
                 IF (LVLS(IFD,IGET(061)).GT.0) CALL GRIBIT(      &
                    IGET(061),LVLS(IFD,IGET(061)),GRID2,IM,JM)
               ENDIF
            ENDIF
 10      CONTINUE
      ENDIF
!     
!
!
!     ***BLOCK 4:  FREEZING LEVEL Z AND RH.
!     
      IF ( (IGET(062).GT.0).OR.(IGET(063).GT.0) ) THEN
         CALL FRZLVL(Z1D,RH1D)
!
!        FREEZING LEVEL HEIGHT.
         IF (IGET(062).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=Z1D(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL BOUND (GRID1,D00,H99999)
            CALL GRIBIT(IGET(062),LVLS(1,IGET(062)),GRID1,IM,JM)
         ENDIF
!
!        FREEZING LEVEL RELATIVE HUMIDITY.
         IF (IGET(063).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH1D(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL SCLFLD(GRID1,H100,IM,JM)
            CALL BOUND(GRID1,H1,H100)
            CALL GRIBIT(IGET(063),LVLS(1,IGET(063)),GRID1,IM,JM)
         ENDIF
      ENDIF
      IF (IGET(165).GT.0 .OR. IGET(350).GT.0) THEN
         CALL FRZLVL2(Z1D,RH1D)
!
!        HIGHEST FREEZING LEVEL HEIGHT.
          IF (IGET(165).GT.0)THEN  
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=Z1D(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL BOUND (GRID1,D00,H99999)
            CALL GRIBIT(IGET(165),LVLS(1,IGET(165)),GRID1,IM,JM)
          END IF
	  
	  IF (IGET(350).GT.0)THEN  
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH1D(I,J)*100.
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL BOUND (GRID1,H1,H100)
            CALL GRIBIT(IGET(350),LVLS(1,IGET(350)),GRID1,IM,JM)
          END IF
      ENDIF
!     
!
!
!     ***BLOCK 5:  BOUNDARY LAYER FIELDS.
!     
      IF ( (IGET(067).GT.0).OR.(IGET(068).GT.0).OR.       &
           (IGET(069).GT.0).OR.(IGET(070).GT.0).OR.       &
           (IGET(071).GT.0).OR.(IGET(072).GT.0).OR.       &
           (IGET(073).GT.0).OR.(IGET(074).GT.0).OR.       &
           (IGET(088).GT.0).OR.(IGET(089).GT.0).OR.       &
           (IGET(090).GT.0).OR.(IGET(075).GT.0).OR.       &
           (IGET(109).GT.0).OR.(IGET(110).GT.0).OR.       &
           (IGET(031).GT.0).OR.(IGET(032).GT.0).OR.       &
           (IGET(107).GT.0).OR.(IGET(091).GT.0).OR.       &
           (IGET(092).GT.0).OR.(IGET(093).GT.0).OR.       &
           (IGET(094).GT.0).OR.(IGET(095).GT.0).OR.       &
           (IGET(096).GT.0).OR.(IGET(097).GT.0).OR.       &
           (IGET(098).GT.0).OR.(IGET(221).GT.0) ) THEN
!
!        COMPUTE ETA BOUNDARY LAYER FIELDS.
         CALL BNDLYR(PBND,TBND,QBND,RHBND,UBND,VBND,      &
              WBND,OMGBND,PWTBND,QCNVBND,LVLBND)

!     
!        LOOP OVER NBND BOUNDARY LAYERS.
         DO 20 LBND = 1,NBND
            ID(1:25) = 0
            ID(10)   = NINT(PETABND(LBND)+15.)
            ID(11)   = NINT(PETABND(LBND)-15.)
!     
!           BOUNDARY LAYER PRESSURE.
            IF (IGET(067).GT.0) THEN
              IF (LVLS(LBND,IGET(067)).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=PBND(I,J,LBND)
               ENDDO
               ENDDO
               CALL GRIBIT(IGET(067),LVLS(LBND,IGET(067)),GRID1,IM,JM)
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER TEMPERATURE.
            IF (IGET(068).GT.0) THEN
              IF (LVLS(LBND,IGET(068)).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=TBND(I,J,LBND)
               ENDDO
               ENDDO
               CALL GRIBIT(IGET(068),LVLS(LBND,IGET(068)),GRID1,IM,JM)
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER POTENTIAL TEMPERATURE.
            IF (IGET(069).GT.0) THEN
              IF (LVLS(LBND,IGET(069)).GT.0) THEN
               CALL CALPOT(PBND(1,1,LBND),TBND(1,1,LBND),EGRID1)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID1(I,J)
               ENDDO
               ENDDO
               CALL GRIBIT(IGET(069),LVLS(LBND,IGET(069)),GRID1,IM,JM)
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER RELATIVE HUMIDITY.
            IF (IGET(072).GT.0) THEN
              IF (LVLS(LBND,IGET(072)).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RHBND(I,J,LBND)
               ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               CALL GRIBIT(IGET(072),LVLS(LBND,IGET(072)),GRID1,IM,JM)
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER DEWPOINT TEMPERATURE.
            IF (IGET(070).GT.0) THEN
              IF (LVLS(LBND,IGET(070)).GT.0) THEN
               CALL CALDWP(PBND(1,1,LBND),QBND(1,1,LBND),EGRID1,    &
                    TBND(1,1,LBND))
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID1(I,J)
               ENDDO
               ENDDO
               CALL GRIBIT(IGET(070),LVLS(LBND,IGET(070)),GRID1,IM,JM)
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER SPECIFIC HUMIDITY.
            IF (IGET(071).GT.0) THEN
              IF (LVLS(LBND,IGET(071)).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=QBND(I,J,LBND)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,H1M12,H99999)
               CALL GRIBIT(IGET(071),LVLS(LBND,IGET(071)),GRID1,IM,JM)
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER MOISTURE CONVERGENCE.
            IF (IGET(088).GT.0) THEN
              IF (LVLS(LBND,IGET(088)).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=QCNVBND(I,J,LBND)
               ENDDO
               ENDDO
               CALL GRIBIT(IGET(088),LVLS(LBND,IGET(088)),GRID1,IM,JM)
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER U WIND AND/OR V WIND.
!
            FIELD1=.FALSE.
            FIELD2=.FALSE.
!
            IF(IGET(073).GT.0)THEN
              IF(LVLS(LBND,IGET(073)).GT.0)FIELD1=.TRUE.
            ENDIF
            IF(IGET(074).GT.0)THEN
              IF(LVLS(LBND,IGET(074)).GT.0)FIELD2=.TRUE.
            ENDIF
!
            IF(FIELD1.OR.FIELD2)THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=UBND(I,J,LBND)
                 GRID2(I,J)=VBND(I,J,LBND)
               ENDDO
               ENDDO
!
               IF (IGET(073).GT.0) THEN
                 IF (LVLS(LBND,IGET(073)).GT.0) CALL GRIBIT(IGET(073), &
                    LVLS(LBND,IGET(073)),GRID1,IM,JM)
               ENDIF
               IF (IGET(074).GT.0) THEN
                 IF (LVLS(LBND,IGET(074)).GT.0) CALL GRIBIT(IGET(074), &
                    LVLS(LBND,IGET(074)),GRID2,IM,JM)
               ENDIF
            ENDIF
!     
!           BOUNDARY LAYER OMEGA.
            IF (IGET(090).GT.0) THEN
              IF (LVLS(LBND,IGET(090)).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=OMGBND(I,J,LBND)
               ENDDO
               ENDDO
               CALL GRIBIT(IGET(090),LVLS(LBND,IGET(090)),GRID1,IM,JM)
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER PRECIPITBLE WATER.
            IF (IGET(089).GT.0) THEN
              IF (LVLS(LBND,IGET(089)).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=PWTBND(I,J,LBND)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               CALL GRIBIT(IGET(089),LVLS(LBND,IGET(089)),GRID1,IM,JM)
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER LIFTED INDEX.
            IF (IGET(075).GT.0) THEN
              IF (LVLS(LBND,IGET(075)).GT.0) THEN
               CALL OTLFT(PBND(1,1,LBND),TBND(1,1,LBND),    &
                    QBND(1,1,LBND),EGRID1)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID1(I,J)
               ENDDO
               ENDDO
               CALL GRIBIT(IGET(075),LVLS(LBND,IGET(075)),GRID1,IM,JM)
              ENDIF
            ENDIF
!
!        END OF ETA BOUNDARY LAYER LOOP.
 20      CONTINUE
!     
!        BEST LIFTED INDEX FROM BOUNDARY LAYER FIELDS.
!     
         IF (IGET(031).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
              EGRID1(I,J) = H99999
              EGRID2(I,J) = H99999
            ENDDO
            ENDDO
!
            DO 50 LBND = 1,NBND
               CALL OTLFT(PBND(1,1,LBND),TBND(1,1,LBND),      &
                    QBND(1,1,LBND),EGRID2)
               DO J=JSTA,JEND
               DO I=1,IM
                 EGRID1(I,J)=AMIN1(EGRID1(I,J),EGRID2(I,J))
               ENDDO
               ENDDO
 50         CONTINUE
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID1(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            ID(10)   = PETABND(NBND)+15.
            ID(11)   = PETABND(1)-15.
            CALL GRIBIT(IGET(031),LVLS(1,IGET(031)),GRID1,IM,JM)
         ENDIF
!     
!        BEST BOUNDARY LAYER CAPE AND CINS.
!     
         FIELD1=.FALSE.
         FIELD2=.FALSE.
!
         IF(IGET(032).GT.0)THEN
           IF(LVLS(2,IGET(032)).GT.0)FIELD1=.TRUE.
         ENDIF
         IF(IGET(107).GT.0)THEN
           IF(LVLS(2,IGET(107)).GT.0)FIELD2=.TRUE.
         ENDIF
!
         IF(FIELD1.OR.FIELD2)THEN
           ITYPE = 2
!
           DO J=JSTA,JEND
           DO I=1,IM
             EGRID1(I,J) = -H99999
             EGRID2(I,J) = -H99999
           ENDDO
           ENDDO
!
           DO 80 LBND = 1,NBND
           CALL CALTHTE(PBND(1,1,LBND),TBND(1,1,LBND),        &
                        QBND(1,1,LBND),EGRID1)
           DO J=JSTA,JEND
           DO I=1,IM
             IF (EGRID1(I,J).GT.EGRID2(I,J)) THEN
               EGRID2(I,J) = EGRID1(I,J)
	       LB2(I,J)  = LVLBND(I,J,LBND)
               P1D(I,J)  = PBND(I,J,LBND)
               T1D(I,J)  = TBND(I,J,LBND)
               Q1D(I,J)  = QBND(I,J,LBND)
             ENDIF
           ENDDO
           ENDDO
 80        CONTINUE
!
           DPBND=0.
           CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,LB2,EGRID1,   &
      	           EGRID2,EGRID3,EGRID4,EGRID5) 
!
           IF (IGET(032).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID1(I,J)
               ENDDO
               ENDDO
             CALL BOUND(GRID1,D00,H99999)
             ID(1:25) = 0
             ID(09)   = 116
             ID(10)   = PETABND(NBND)+15.
             ID(11)   = PETABND(1)-15.
             CALL GRIBIT(IGET(032),LVLS(1,IGET(032)),GRID1,IM,JM)
           ENDIF
!
           IF (IGET(107).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID2(I,J)
               ENDDO
               ENDDO
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = -1.*GRID1(I,J)
             ENDDO
             ENDDO
!
             CALL BOUND(GRID1,D00,H99999)
!
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = -1.*GRID1(I,J)
             ENDDO
             ENDDO
!
             ID(1:25) = 0
             ID(09)   = 116
             ID(10)   = PETABND(NBND)+15.
             ID(11)   = PETABND(1)-15.
             CALL GRIBIT(IGET(107),LVLS(1,IGET(107)),GRID1,IM,JM)
           ENDIF
         ENDIF
!

!    PBL HEIGHT 
         IF(IGET(221).GT.0) THEN
	   DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J)=PBLH(I,J)
           ENDDO
           ENDDO
	   ID(1:25) = 0
	   CALL GRIBIT(IGET(221),LVLS(1,IGET(221)),GRID1,IM,JM)
         END IF
!        BOUNDARY LAYER LIFTING CONDENSATION PRESSURE AND HEIGHT.
!        EGRID1 IS LCL PRESSURE.  EGRID2 IS LCL HEIGHT.
!
         IF ( (IGET(109).GT.0).OR.(IGET(110).GT.0) ) THEN
            CALL CALLCL(PBND(1,1,1),TBND(1,1,1),          &
                 QBND(1,1,1),EGRID1,EGRID2)
            IF (IGET(109).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID2(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(109),ILVL, GRID1,IM,JM)
            ENDIF
            IF (IGET(110).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID1(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(110),ILVL, GRID1,IM,JM)
            ENDIF
         ENDIF
!     
!        NGM BOUNDARY LAYER FIELDS.
!     
         IF ( (IGET(091).GT.0).OR.(IGET(092).GT.0).OR.      &
              (IGET(093).GT.0).OR.(IGET(094).GT.0).OR.      &
              (IGET(095).GT.0).OR.(IGET(095).GT.0).OR.      &
              (IGET(096).GT.0).OR.(IGET(097).GT.0).OR.      &
              (IGET(098).GT.0) ) THEN
!
!  COMPUTE SIGMA 0.89671 AND 0.78483 TEMPERATURES
!    INTERPOLATE LINEAR IN LOG P
            IF (IGET(097).GT.0.OR.IGET(098).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 P78483(I,J)=ALOG(PINT(I,J,NINT(LMH(I,J)))*0.78483)
                 P89671(I,J)=ALOG(PINT(I,J,NINT(LMH(I,J)))*0.89671)
               ENDDO
               ENDDO
!$omp  parallel do
!$omp& private(fac1,fac2,pkl1,pku1,t78483,t89671)
               DONE =.FALSE.
	       DONE1=.FALSE.
               DO L=2,LM
                DO J=JSTA,JEND
                DO I=1,IM                  
                  PKL1=0.5*(ALPINT(I,J,L)+ALPINT(I,J,L+1))
                  PKU1=0.5*(ALPINT(I,J,L)+ALPINT(I,J,L-1))
		  IF(I.EQ.1 .AND. J.EQ.1)PRINT*,'L,P89671,PKL1,PKU1= ', &
                    L,P89671(I,J), PKL1, PKU1   		  
                  IF(P78483(I,J).LT.PKL1.AND.P78483(I,J).GT.PKU1)THEN
                    FAC1=(PKL1-P78483(I,J))/(PKL1-PKU1)
                    FAC2=(P78483(I,J)-PKU1)/(PKL1-PKU1)
                    T78483(I,J)=T(I,J,L)*FAC2+T(I,J,L-1)*FAC1
		    DONE1(I,J)=.TRUE.
                  ENDIF
                  IF(P89671(I,J).LT.PKL1.AND.P89671(I,J).GT.PKU1)THEN
                    FAC1=(PKL1-P89671(I,J))/(PKL1-PKU1)
                    FAC2=(P89671(I,J)-PKU1)/(PKL1-PKU1)
                    T89671(I,J)=T(I,J,L)*FAC2+T(I,J,L-1)*FAC1
                    DONE(I,J)=.TRUE.
		    IF(I.EQ.1 .AND. J.EQ.1)PRINT*,'done(1,1)= ',done(1,1)
                  ENDIF
                ENDDO
                ENDDO
               ENDDO
!	       print*,'done(1,1)= ',done(1,1)
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(.NOT. DONE(I,J))THEN
		   PL=PINT(I,J,LM-1)
                   TL=0.5*(T(I,J,LM-2)+T(I,J,LM-1))
                   QL=0.5*(Q(I,J,LM-2)+Q(I,J,LM-1))
                   QSAT=PQ0/PL *EXP(A2*(TL-A3)/(TL-A4))
!
                   RHL=QL/QSAT
!
                   IF(RHL.GT.1.)THEN
                    RHL=1.
                    QL =RHL*QSAT
                   ENDIF
!
                   IF(RHL.LT.RHmin)THEN
                    RHL=RHmin
                    QL =RHL*QSAT
                   ENDIF
!
                   TVRL  =TL*(1.+0.608*QL)
                   TVRBLO=TVRL*(P89671(I,J)/PL)**RGAMOG
                   T89671(I,J)  =TVRBLO/(1.+0.608*QL)
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
		   PL=PINT(I,J,LM-1)
                   TL=0.5*(T(I,J,LM-2)+T(I,J,LM-1))
                   QL=0.5*(Q(I,J,LM-2)+Q(I,J,LM-1))
                   QSAT=PQ0/PL *EXP(A2*(TL-A3)/(TL-A4))
!
                   RHL=QL/QSAT
!
                   IF(RHL.GT.1.)THEN
                    RHL=1.
                    QL =RHL*QSAT
                   ENDIF
!
                   IF(RHL.LT.RHmin)THEN
                    RHL=RHmin
                    QL =RHL*QSAT
                   ENDIF
!
                   TVRL  =TL*(1.+0.608*QL)
                   TVRBLO=TVRL*(P78483(I,J)/PL)**RGAMOG
                   T78483(I,J)  =TVRBLO/(1.+0.608*QL)
!     
                 END IF
		 
               END DO
             END DO
!     
!           SIGMA 0.89671 TEMPERATURE
             IF (IGET(097).GT.0) THEN
               ID(1:25) = 0
               ISVALUE = 8967
               ID(11) = ISVALUE
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=T89671(I,J)
                 IF(T89671(I,J).GT.350.)PRINT*,'LARGE T89671 ',   &
                   I,J,T89671(I,J)
               ENDDO
               ENDDO
               CALL GRIBIT(IGET(097),LVLS(1,IGET(097)),GRID1,IM,JM)
             ENDIF
!     
!           SIGMA 0.78483 TEMPERATURE
             IF (IGET(098).GT.0) THEN
               ID(1:25) = 0
               ISVALUE = 7848
               ID(11) = ISVALUE
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=T78483(I,J)
               ENDDO
               ENDDO
               CALL GRIBIT(IGET(098),LVLS(1,IGET(098)),GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!           NGM SIGMA LAYER 0.98230 FIELDS.  THESE FIELDS ARE 
!           THE FIRST ETA LAYER BOUNDARY LAYER FIELDS. 
!     
!     
            IF ( (IGET(091).GT.0).OR.(IGET(092).GT.0).OR.      &
                 (IGET(093).GT.0).OR.(IGET(094).GT.0).OR.      &
                 (IGET(095).GT.0).OR.(IGET(095).GT.0).OR.      &
                 (IGET(096).GT.0) ) THEN
!     
               ID(1:25) = 0
               ISVALUE = 9823
               ID(11) = ISVALUE
!     
!              PRESSURE.
               IF (IGET(091).GT.0) THEN
                 DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J)=PBND(I,J,1)
                 ENDDO
                 ENDDO
                  CALL GRIBIT(IGET(091),LVLS(1,IGET(091)),GRID1,IM,JM)
               ENDIF
!     
!              TEMPERATURE.
               IF (IGET(092).GT.0) THEN
                 DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J)=TBND(I,J,1)
                 ENDDO
                 ENDDO
                  CALL GRIBIT(IGET(092),LVLS(1,IGET(092)),       &
                       GRID1,IM,JM)
               ENDIF
!     
!              SPECIFIC HUMIDITY.
               IF (IGET(093).GT.0) THEN
                 DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J)=QBND(I,J,1)
                 ENDDO
                 ENDDO
                  CALL BOUND(GRID1,H1M12,H99999)
                  CALL GRIBIT(IGET(093),LVLS(1,IGET(093)),GRID1,IM,JM)
               ENDIF
!     
!              RELATIVE HUMIDITY.
               IF (IGET(094).GT.0) THEN
                 DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J)=RHBND(I,J,1)
                 ENDDO
                 ENDDO
                  CALL SCLFLD(GRID1,H100,IM,JM)
                  CALL BOUND(GRID1,H1,H100)
                  CALL GRIBIT(IGET(094),LVLS(1,IGET(094)),GRID1,IM,JM)
               ENDIF
!     
!              U AND/OR V WIND.
               IF ((IGET(095).GT.0).OR.(IGET(096).GT.0)) THEN
                 DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J)=UBND(I,J,1)
                   GRID2(I,J)=VBND(I,J,1)
                 ENDDO
                 ENDDO
                  IF (IGET(095).GT.0) CALL GRIBIT(IGET(095),       &
                       LVLS(1,IGET(095)),GRID1,IM,JM)
                  IF (IGET(096).GT.0) CALL GRIBIT(IGET(096),       &
                       LVLS(1,IGET(096)),GRID2,IM,JM)
               ENDIF
            ENDIF
         ENDIF
!     
!     ENDIF FOR BOUNDARY LAYER BLOCK.
!
      ENDIF
!     
!
!
!     ***BLOCK 6:  MISCELLANEOUS LAYER MEAN LFM AND NGM FIELDS.
!     
      IF ( (IGET(066).GT.0).OR.(IGET(081).GT.0).OR.        &
           (IGET(082).GT.0).OR.(IGET(104).GT.0).OR.        &
           (IGET(099).GT.0).OR.(IGET(100).GT.0).OR.        &
           (IGET(101).GT.0).OR.(IGET(102).GT.0).OR.        &
           (IGET(103).GT.0) ) THEN
!     
!        LFM "MEAN" RELATIVE HUMIDITIES AND PRECIPITABLE WATER.
!     
         IF ( (IGET(066).GT.0).OR.(IGET(081).GT.0).OR.     &
              (IGET(082).GT.0).OR.(IGET(104).GT.0) ) THEN
            CALL LFMFLD(RH3310,RH6610,RH3366,PW3310)
            ID(1:25) = 0
!     
!           SIGMA 0.33-1.00 MEAN RELATIVE HUMIIDITY.
            IF (IGET(066).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH3310(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10) = 33
               ID(11) = 100
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               CALL GRIBIT(IGET(066),LVLS(1,IGET(066)),   &
                    GRID1,IM,JM)
            ENDIF
!     
!           SIGMA 0.66-1.00 MEAN RELATIVE HUMIIDITY.
            IF (IGET(081).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH6610(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10) = 67
               ID(11) = 100
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               CALL GRIBIT(IGET(081),LVLS(1,IGET(081)),GRID1,IM,JM)
            ENDIF
!     
!           SIGMA 0.33-0.66 MEAN RELATIVE HUMIIDITY.
            IF (IGET(082).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH3366(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10) = 33
               ID(11) = 67
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               CALL GRIBIT(IGET(082),LVLS(1,IGET(082)),GRID1,IM,JM)
            ENDIF
!     
!           SIGMA 0.33-1.00 PRECIPITABLE WATER.
            IF (IGET(104).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=PW3310(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10) = 33
               ID(11) = 100
               CALL BOUND(GRID1,D00,H99999)
               CALL GRIBIT(IGET(104),LVLS(1,IGET(104)),GRID1,IM,JM)
            ENDIF
         ENDIF
!     
!        VARIOUS LAYER MEAN NGM SIGMA FIELDS.
!     
         IF ( (IGET(099).GT.0).OR.(IGET(100).GT.0).OR.    &
              (IGET(101).GT.0).OR.(IGET(102).GT.0).OR.    &
              (IGET(103).GT.0) ) THEN
            CALL NGMFLD(RH4710,RH4796,RH1847,RH8498,QM8510)
!     
!           SIGMA 0.47191-1.00000 RELATIVE HUMIDITY.
            IF (IGET(099).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH4710(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10)   = 47
               ID(11)   = 100
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               CALL GRIBIT(IGET(099),LVLS(1,IGET(099)),GRID1,IM,JM)
            ENDIF
!     
!           SIGMA 0.47191-0.96470 RELATIVE HUMIDITY.
            IF (IGET(100).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH4796(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10)   = 47
               ID(11)   = 96
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               CALL GRIBIT(IGET(100),LVLS(1,IGET(100)),GRID1,IM,JM)
            ENDIF
!     
!           SIGMA 0.18019-0.47191 RELATIVE HUMIDITY.
            IF (IGET(101).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH1847(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10)   = 18
               ID(11)   = 47
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               CALL GRIBIT(IGET(101),LVLS(1,IGET(101)),GRID1,IM,JM)
            ENDIF
!     
!           SIGMA 0.84368-0.98230 RELATIVE HUMIDITY.
            IF (IGET(102).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH8498(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10)   = 84
               ID(11)   = 98
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               CALL GRIBIT(IGET(102),LVLS(1,IGET(102)),GRID1,IM,JM)
            ENDIF
!     
!           SIGMA 0.85000-1.00000 MOISTURE CONVERGENCE.
            IF (IGET(103).GT.0) THEN
!           CONVERT TO DIVERGENCE FOR GRIB
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=-1.0*QM8510(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10)   = 85
               ID(11)   = 100
               CALL GRIBIT(IGET(103),LVLS(1,IGET(103)),GRID1,IM,JM)
            ENDIF
         ENDIF
      ENDIF

      IF ( (IGET(318).GT.0).OR.(IGET(319).GT.0).OR.     &
           (IGET(320).GT.0))THEN
       CALL LFMFLD_GFS(RH4410,RH7294,RH4472,RH3310)
!     
!           SIGMA 0.44-1.00 MEAN RELATIVE HUMIIDITY.
            IF (IGET(318).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH4410(I,J)*100.
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10) = 44
               ID(11) = 100
               CALL BOUND(GRID1,D00,H100)
               CALL GRIBIT(IGET(318),LVLS(1,IGET(318)),  &
                    GRID1,IM,JM)
            ENDIF
!     
!           SIGMA 0.72-0.94 MEAN RELATIVE HUMIIDITY.
            IF (IGET(319).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH7294(I,J)*100.
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10) = 72
               ID(11) = 94
               CALL BOUND(GRID1,D00,H100)
               CALL GRIBIT(IGET(319),LVLS(1,IGET(319)),GRID1,IM,JM)
            ENDIF
!     
!           SIGMA 0.44-0.72 MEAN RELATIVE HUMIIDITY.
            IF (IGET(320).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH4472(I,J)*100.
               ENDDO
               ENDDO
               ID(1:25) = 0
               ID(10) = 44
               ID(11) = 72
               CALL BOUND(GRID1,D00,H100)
               CALL GRIBIT(IGET(320),LVLS(1,IGET(320)),GRID1,IM,JM)
            ENDIF
      END IF  	    	    
! GFS computes sigma=0.9950 T, THETA, U, V from lowest two model level fields 
         IF ( (IGET(321).GT.0).OR.(IGET(322).GT.0).OR.     &
              (IGET(323).GT.0).OR.(IGET(324).GT.0).OR.     &
              (IGET(325).GT.0).OR.(IGET(326).GT.0) ) THEN
           DO J=JSTA,JEND
	     DO I=1,IM
               EGRID2(I,J)=0.995*PINT(I,J,LM+1)
               EGRID1(I,J)=LOG(PMID(I,J,LM)/EGRID2(I,J))   &
        	   /LOG(PMID(I,J,LM)/PMID(I,J,LM-1))
	     END DO
	   END DO
! Temperature	   
	   IF (IGET(321).GT.0) THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=T(I,J,LM)+(T(I,J,LM-1)-T(I,J,LM)) &
      	       *EGRID1(I,J)
             ENDDO
             ENDDO
             ID(1:25) = 0
             ID(11) = 9950
             CALL GRIBIT(IGET(321),LVLS(1,IGET(321)),GRID1,IM,JM)
            ENDIF
! Potential Temperature	    
	    IF (IGET(322).GT.0) THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID2(I,J)=T(I,J,LM)+(T(I,J,LM-1)-T(I,J,LM))     &
      	       *EGRID1(I,J)
             ENDDO
             ENDDO
	     CALL CALPOT(EGRID2,GRID2,GRID1)
             ID(1:25) = 0
             ID(11) = 9950
             CALL GRIBIT(IGET(322),LVLS(1,IGET(322)),GRID1,IM,JM)
            ENDIF
! RH	    
	    IF (IGET(323).GT.0) THEN
             DO J=JSTA,JEND
             DO I=1,IM
	       ES1=FPVSNEW(T(I,J,LM))
	       ES1=MIN(ES1,PMID(I,J,LM))
	       QS1=CON_EPS*ES1/(PMID(I,J,LM)+CON_EPSM1*ES1)
	       RH1=Q(I,J,LM)/QS1
	       ES2=FPVSNEW(T(I,J,LM-1))
	       ES2=MIN(ES2,PMID(I,J,LM-1))
	       QS2=CON_EPS*ES2/(PMID(I,J,LM-1)+CON_EPSM1*ES2)
	       RH2=Q(I,J,LM-1)/QS2
               GRID1(I,J)=(RH1+(RH2-RH1)*EGRID1(I,J))*100.
             ENDDO
             ENDDO
             CALL BOUND(GRID1,D00,H100)
             ID(1:25) = 0
             ID(11) = 9950
             CALL GRIBIT(IGET(323),LVLS(1,IGET(323)),GRID1,IM,JM)
            ENDIF	    
! U	   
	   IF (IGET(324).GT.0) THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=UH(I,J,LM)+(UH(I,J,LM-1)-UH(I,J,LM))    &
      	       *EGRID1(I,J)
             ENDDO
             ENDDO
             ID(1:25) = 0
             ID(11) = 9950
             CALL GRIBIT(IGET(324),LVLS(1,IGET(324)),GRID1,IM,JM)
            ENDIF	    
! V	   
	   IF (IGET(325).GT.0) THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=VH(I,J,LM)+(VH(I,J,LM-1)-VH(I,J,LM))    &
           	          *EGRID1(I,J)
             ENDDO
             ENDDO
             ID(1:25) = 0
             ID(11) = 9950
             CALL GRIBIT(IGET(325),LVLS(1,IGET(325)),GRID1,IM,JM)
            ENDIF	       
! OMGA	   
	   IF (IGET(326).GT.0) THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=OMGA(I,J,LM)+(OMGA(I,J,LM-1)-OMGA(I,J,LM))  &
      	       *EGRID1(I,J)
             ENDDO
             ENDDO
             ID(1:25) = 0
             ID(11) = 9950
             CALL GRIBIT(IGET(326),LVLS(1,IGET(326)),GRID1,IM,JM)
            ENDIF
         END IF
	 
!       MIXED LAYER CAPE AND CINS
!
         FIELD1=.FALSE.
         FIELD2=.FALSE.
!
         IF(IGET(032).GT.0)THEN
           IF(LVLS(3,IGET(032)).GT.0)FIELD1=.TRUE.
         ENDIF
         IF(IGET(107).GT.0)THEN
           IF(LVLS(3,IGET(107)).GT.0)FIELD2=.TRUE.
         ENDIF
!
         IF(FIELD1.OR.FIELD2)THEN
           ITYPE = 2
!
           DO J=JSTA,JEND
           DO I=1,IM
             EGRID1(I,J) = -H99999
             EGRID2(I,J) = -H99999
           ENDDO
           ENDDO
            
            
           DO J=JSTA,JEND
           DO I=1,IM
               LB2(I,J)  = (LVLBND(I,J,1) + LVLBND(I,J,2) +           &
                            LVLBND(I,J,3))/3
               P1D(I,J)  = (PBND(I,J,1) + PBND(I,J,2) + PBND(I,J,3))/3
               T1D(I,J)  = (TBND(I,J,1) + TBND(I,J,2) + TBND(I,J,3))/3
               Q1D(I,J)  = (QBND(I,J,1) + QBND(I,J,2) + QBND(I,J,3))/3
           ENDDO
           ENDDO
!
           DPBND=0.
           CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,LB2,EGRID1,           &
                EGRID2,EGRID3,EGRID4,EGRID5)
 
           IF (IGET(032).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID1(I,J)
               ENDDO
               ENDDO

               CALL BOUND(GRID1,D00,H99999)
               ID(1:25) = 0
               ID(09)   = 116
               ID(10)   = PETABND(3)+15.
               ID(11)   = PETABND(1)-15.
	       CALL GRIBIT(IGET(32),LVLS(3,IGET(32)),GRID1,IM,JM)
           ENDIF
                                                                                               
           IF (IGET(107).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=-1.*EGRID2(I,J)
               ENDDO
               ENDDO	    
!
               CALL BOUND(GRID1,D00,H99999)
!
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = -1.*GRID1(I,J)
               ENDDO
               ENDDO
!
               ID(1:25) = 0
               ID(09)   = 116
               ID(10)   = PETABND(3)+15.
               ID(11)   = PETABND(1)-15.
	       CALL GRIBIT(IGET(107),LVLS(3,IGET(107)),          &
                    GRID1,IM,JM)
           ENDIF
         ENDIF
              
!        MIXED LAYER LIFTING CONDENSATION PRESSURE AND HEIGHT.
!        EGRID1 IS LCL PRESSURE.  EGRID2 IS LCL HEIGHT.
!
!         IF ( (IGET(109).GT.0).OR.(IGET(110).GT.0) ) THEN
!            CALL CALLCL(P1D,T1D,Q1D,EGRID1,EGRID2)
!            IF (IGET(109).GT.0) THEN
!	       DO J=JSTA,JEND
!               DO I=1,IM
!                 GRID1(I,J)=EGRID2(I,J)
!               ENDDO
!               ENDDO
!           
!               ID(1:25) = 0
!	       
!	       CALL GRIBIT(IGET(109),1,
!     X              GRID1,IM,JM)
!            ENDIF
!	    
!            IF (IGET(110).GT.0) THEN
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
!         ENDIF
!
!       MOST UNSTABLE CAPE-LOWEST 300 MB
!
         FIELD1=.FALSE.
         FIELD2=.FALSE.
!
         IF(IGET(032).GT.0)THEN
           IF(LVLS(4,IGET(032)).GT.0)FIELD1=.TRUE.
         ENDIF
              
         IF(IGET(107).GT.0)THEN
           IF(LVLS(4,IGET(107)).GT.0)FIELD2=.TRUE.
         ENDIF
!
         IF(FIELD1.OR.FIELD2)THEN
           ITYPE = 1
!
           DO J=JSTA,JEND
           DO I=1,IM
             EGRID1(I,J) = -H99999
             EGRID2(I,J) = -H99999
           ENDDO
           ENDDO
              
           DPBND=300.E2
           CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,LB2,EGRID1,     &
                   EGRID2,EGRID3,EGRID4,EGRID5)
!
           IF (IGET(032).GT.0) THEN
	       DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID1(I,J)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               ID(1:25) = 0
               ID(09)   = 116
               ID(10) = 255
               ID(11) = 0
	       CALL GRIBIT(IGET(32),4,GRID1,IM,JM)
           ENDIF
                
           IF (IGET(107).GT.0) THEN
	       DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=-1.0*EGRID2(I,J)
               ENDDO
               ENDDO

               CALL BOUND(GRID1,D00,H99999)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = -1.*GRID1(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               ID(09)   = 116
               ID(10) = 255
               ID(11) = 0
	       CALL GRIBIT(IGET(107),4,GRID1,IM,JM)
            ENDIF

!    EQUILLIBRIUM HEIGHT
           IF (IGET(443).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID4(I,J)
               ENDDO
               ENDDO
             ID(1:25) = 0
             CALL GRIBIT(IGET(443),LVLS(1,IGET(443)),GRID1,IM,JM)
           ENDIF
              
!    GENERAL THUNDER PARAMETER
        IF (IGET(458).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 IF (CPRATE(I,J) .GT. PTHRESH) THEN
                  GRID1(I,J)=EGRID5(I,J)
                 ELSE
                  GRID1(I,J)=0
                 ENDIF
               ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               ID(1:25) = 0
               CALL GRIBIT(IGET(458),1,GRID1,IM,JM)
           ENDIF

!      PRESSURE OF LEVEL FROM WHICH 300 MB MOST UNSTABLE CAPE
!             PARCEL WAS LIFTED
           IF (IGET(246).GT.0) THEN
	       DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID3(I,J)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               ID(1:25) = 0
               ID(02) = 129
               ID(09) = 116
               ID(10) = 255
               ID(11) = 0
	       CALL GRIBIT(IGET(246),1,GRID1,IM,JM)
           ENDIF
         ENDIF  

! CALCULATE LPBL
      IF (IGET(245).GT.0) THEN
       DO 101 J=JSTA,JEND
        DO 101 I=1,IM
         ZSF=ZINT(I,J,NINT(LMH(I,J))+1)
         DO L=NINT(LMH(I,J)),1,-1
          IF(ZINT(I,J,L) .GT. PBLH(I,J)+ZSF)THEN
           LPBL(I,J)=L+1
           IF(LPBL(I,J).GE.LP1) LPBL(I,J) = LM
           GO TO 101
          END IF
         END DO
 101   CONTINUE
       CALL CALGUST(LPBL,PBLH,GUST)
       DO J=JSTA,JEND
       DO I=1,IM
         if(GUST(I,J) .gt. 200. .and. gust(i,j).lt.spval)    &
      	 print*,'big gust at ',i,j
         GRID1(I,J)=GUST(I,J)
       ENDDO
       ENDDO      
       ID(1:25) = 0
       CALL GRIBIT(IGET(245),1,GRID1,IM,JM)      
      END IF
!    
 
!     END OF ROUTINE.
!     
      RETURN
      END
