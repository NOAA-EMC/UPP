      SUBROUTINE CALPW(PW,IDECID)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALPW       COMPUTES 
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-24       
!     
! ABSTRACT:  
!     THIS ROUTINE COMPUTES PRECIPITABLE WATER IN A COLUMN
!     EXTENDING FROM THE FIRST ATMOSPHERIC ETA LAYER TO THE
!     MODEL TOP.  THE DEFINITION USED IS
!                                 TOP
!            PRECIPITABLE WATER = SUM (Q+CLDW) DP*HTM/G
!                                 BOT
!     WHERE,
!        BOT IS THE FIRST ETA LAYER,
!        TOP IS THE MODEL TOP,
!        Q IS THE SPECIFIC HUMIDITY (KG/KG) IN THE LAYER
!        CLDW IS THE CLOUD WATER (KG/KG) IN THE LAYER
!        DP (Pa) IS THE LAYER THICKNESS.
!        HTM IS THE HEIGHT MASK AT THAT LAYER (=0 IF BELOW GROUND)
!        G IS THE GRAVITATIONAL CONSTANT
!     
! PROGRAM HISTORY LOG:
!   92-12-24  RUSS TREADON
!   96-03-04  MIKE BALDWIN - ADD CLOUD WATER AND SPEED UP CODE
!   98-06-15  T BLACK      - CONVERSION FROM 1-D TO 2-D
!   00-01-04  JIM TUCCILLO - MPI VERSION                 
!   02-06-19  MIKE BALDWIN - WRF VERSION 
!   04-12-30  H CHUANG      - UPDATE TO CALCULATE TOTAL COLUMN FOR OTHER
!                                     HYDROMETEORS                
!   11-12-14  SARAH LU     - UPDATE TO CALCULATE AEROSOL OPTICAL DEPTH
!     
! USAGE:    CALL CALPW(PW)
!   INPUT ARGUMENT LIST:
!     PW       - ARRAY OF PRECIPITABLE WATER.
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
!       COMMON   - LOOPS
!                  MASKS
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
!     
      use vrbls3d, only: q, qqw, qqi, qqr, qqs, cwm, qqg, t, rswtt, train, tcucn, mcvg,&
              pmid, o3, ext, pint, rlwtt
      use masks, only: htm
      use params_mod, only: tfrz, gi
      use ctlblk_mod, only: lm, jsta, jend, im, jm
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     
!     SET DENSITY OF WATER AT 1 ATMOSPHERE PRESSURE, 0C.
!     UNITS ARE KG/M**3.
      real,PARAMETER :: RHOWAT=1.E3
      real,parameter:: con_rd      =2.8705e+2 ! gas constant air    (J/kg/K)
      real,parameter:: con_rv      =4.6150e+2 ! gas constant H2O
      real,parameter:: con_eps     =con_rd/con_rv
      real,parameter:: con_epsm1   =con_rd/con_rv-1
!     
!     DECLARE VARIABLES.
!     
      integer,intent(in)  ::  IDECID
      real,dimension(IM,JM),intent(inout) :: PW
      INTEGER LLMH,I,J,L
      REAL ALPM,DZ,PM,PWSUM,RHOAIR,DP,ES
      real,external :: FPVSNEW
      REAL QDUM(IM,JM)
      REAL PWS(IM,JM),QS(IM,JM)
!
!***************************************************************
!     START CALPW HERE.
!
!     INITIALIZE PW TO 0.    
!     
      PW = 0.
      PWS = 0.
!     
!     OUTER LOOP OVER VERTICAL DIMENSION.
!     INNER LOOP OVER HORIZONTAL GRID.
!     
!$omp  parallel do
!$omp& private(dp)
      DO L = 1,LM
        IF (IDECID .LE. 1) THEN
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=Q(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID .EQ. 2) THEN
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=QQW(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID .EQ. 3) THEN
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=QQI(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID .EQ. 4) THEN
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=QQR(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID .EQ. 5) THEN
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=QQS(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID .EQ. 6) THEN
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=CWM(I,J,L)
            ENDDO
          ENDDO
! SRD
        ELSE IF (IDECID .EQ. 16) THEN
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=QQG(I,J,L)
            ENDDO
          ENDDO
! SRD
        ELSE IF (IDECID .EQ. 7) THEN
!-- Total supercooled liquid
          DO J=JSTA,JEND
            DO I=1,IM
              IF (T(I,J,L) .GE. TFRZ) THEN
                Qdum(I,J)=0.
              ELSE
                Qdum(I,J)=QQW(I,J,L)+QQR(I,J,L)
              ENDIF
            ENDDO
          ENDDO
        ELSE IF (IDECID .EQ. 8) THEN
!-- Total melting ice
          DO J=JSTA,JEND
            DO I=1,IM
              IF (T(I,J,L) .LE. TFRZ) THEN
                Qdum(I,J)=0.
              ELSE
                Qdum(I,J)=QQI(I,J,L)+QQS(I,J,L)
              ENDIF
            ENDDO
          ENDDO
        ELSE IF (IDECID .EQ. 9) THEN
! SHORT WAVE T TENDENCY
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=RSWTT(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID .EQ. 10) THEN
! LONG WAVE T TENDENCY
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=RLWTT(I,J,L)
            ENDDO
          ENDDO	  
        ELSE IF (IDECID .EQ. 11) THEN
! LATENT HEATING FROM GRID SCALE RAIN/EVAP
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=TRAIN(I,J,L)
            ENDDO
          ENDDO	  
        ELSE IF (IDECID .EQ. 12) THEN
! LATENT HEATING FROM CONVECTION
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=TCUCN(I,J,L)
            ENDDO
          ENDDO	  	  
        ELSE IF (IDECID .EQ. 13) THEN
! MOISTURE CONVERGENCE
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J)=MCVG(I,J,L)
            ENDDO
          ENDDO
! RH
	ELSE IF (IDECID .EQ. 14) THEN
          DO J=JSTA,JEND
            DO I=1,IM
	      Qdum(I,J)=Q(I,J,L)
	      ES=FPVSNEW(T(I,J,L))
              ES=MIN(ES,PMID(I,J,L))
              QS(I,J)=CON_EPS*ES/(PMID(I,J,L)+CON_EPSM1*ES)
	    ENDDO
	  END DO
! OZONE
	ELSE IF (IDECID .EQ. 15) THEN
          DO J=JSTA,JEND
            DO I=1,IM
	      Qdum(I,J)=O3(I,J,L)
	    ENDDO
	  END DO	  

! AEROSOL EXTINCTION (GOCART)
	ELSE IF (IDECID .EQ. 17) THEN
          DO J=JSTA,JEND
            DO I=1,IM
 	      Qdum(I,J)=EXT(I,J,L)
	    ENDDO
	  END DO	  

        ENDIF

        DO J=JSTA,JEND
          DO I=1,IM
            DP   =PINT(I,J,L+1)-PINT(I,J,L)
            PW(I,J)=PW(I,J)+Qdum(I,J)*DP*GI*HTM(I,J,L)
            IF (IDECID .EQ. 17) THEN
             PW(I,J)=PW(I,J)+Qdum(I,J)*MAX(DP,0.)*GI*HTM(I,J,L)
            ENDIF
	    IF (IDECID .EQ. 14) PWS(I,J)=PWS(I,J)                           & 
      	    +QS(I,J)*DP*GI*HTM(I,J,L)
          ENDDO
        ENDDO

      ENDDO

      IF (IDECID .EQ. 14)THEN
        DO J=JSTA,JEND
          DO I=1,IM
            PW(I,J)=max(0.,PW(I,J)/PWS(I,J)*100.) 
          ENDDO
        ENDDO
      END IF
!  convert ozone from kg/m2 to dobson units, which give the depth of the
!  ozone layer in 1e-5 m if brought to natural temperature and pressure.    
      IF (IDECID .EQ. 15)PW(:,:)=PW(:,:)/2.14e-5

!
!     END OF ROUTINE.
!     
      RETURN
      END
