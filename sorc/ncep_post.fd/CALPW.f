!> @file
!                .      .    .     
!> SUBPROGRAM:    CALPW       COMPUTES 
!!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-24       
!!     
!! ABSTRACT:  
!!     THIS ROUTINE COMPUTES PRECIPITABLE WATER IN A COLUMN
!!     EXTENDING FROM THE FIRST ATMOSPHERIC ETA LAYER TO THE
!!     MODEL TOP.  THE DEFINITION USED IS
!!                                 TOP
!!            PRECIPITABLE WATER = SUM (Q+CLDW) DP*HTM/G
!!                                 BOT
!!     WHERE,
!!        BOT IS THE FIRST ETA LAYER,
!!        TOP IS THE MODEL TOP,
!!        Q IS THE SPECIFIC HUMIDITY (KG/KG) IN THE LAYER
!!        CLDW IS THE CLOUD WATER (KG/KG) IN THE LAYER
!!        DP (Pa) IS THE LAYER THICKNESS.
!!        HTM IS THE HEIGHT MASK AT THAT LAYER (=0 IF BELOW GROUND)
!!        G IS THE GRAVITATIONAL CONSTANT
!!     
!! PROGRAM HISTORY LOG:
!!   92-12-24  RUSS TREADON
!!   96-03-04  MIKE BALDWIN - ADD CLOUD WATER AND SPEED UP CODE
!!   98-06-15  T BLACK      - CONVERSION FROM 1-D TO 2-D
!!   00-01-04  JIM TUCCILLO - MPI VERSION                 
!!   02-06-19  MIKE BALDWIN - WRF VERSION 
!!   04-12-30  H CHUANG      - UPDATE TO CALCULATE TOTAL COLUMN FOR OTHER
!!                                     HYDROMETEORS                
!!   14-11-12  SARAH LU     - UPDATE TO CALCULATE AEROSOL OPTICAL DEPTH
!!   15-07-02  SARAH LU     - UPDATE TO CALCULATE SCATTERING AEROSOL
!!                            OPTICAL DEPTH (18)
!!   15-07-04  SARAH LU     - CORRECT PW INTEGRATION FOR AOD (17)
!!   15-07-10  SARAH LU     - UPDATE TO CALCULATE ASYMETRY PARAMETER
!!   19-07-25  Li(Kate) Zhang - MERGE SARHA LU's update for FV3-Chem
!!   20-11-10  JESSE MENG   - USE UPP_PHYSICS MODULE
!!     
!! USAGE:    CALL CALPW(PW)
!!   INPUT ARGUMENT LIST:
!!     PW       - ARRAY OF PRECIPITABLE WATER.
!!
!!   OUTPUT ARGUMENT LIST: 
!!     NONE     
!!     
!!   OUTPUT FILES:
!!     NONE
!!     
!!   SUBPROGRAMS CALLED:
!!     UTILITIES:
!!       NONE
!!     LIBRARY:
!!       COMMON   - LOOPS
!!                  MASKS
!!     
!!   ATTRIBUTES:
!!     LANGUAGE: FORTRAN
!!     MACHINE : CRAY C-90
!!
      SUBROUTINE CALPW(PW,IDECID)

!     
      use vrbls3d,    only: q, qqw, qqi, qqr, qqs, cwm, qqg, t, rswtt,    &
                            train, tcucn, mcvg, pmid, o3, ext, pint, rlwtt, &
                            taod5503d,sca, asy
      use vrbls4d,    only: smoke
      use masks,      only: htm
      use params_mod, only: tfrz, gi
      use ctlblk_mod, only: lm, jsta, jend, im, spval
      use upp_physics, only: FPVSNEW
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
      real,dimension(IM,jsta:jend),intent(inout) :: PW
      INTEGER LLMH,I,J,L
      REAL ALPM,DZ,PM,PWSUM,RHOAIR,DP,ES
      REAL QDUM(IM,jsta:jend), PWS(IM,jsta:jend),QS(IM,jsta:jend)
!
!***************************************************************
!     START CALPW HERE.
!
!     INITIALIZE PW TO 0.    
!     
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=1,IM
          PW(i,j)  = 0.
          PWS(i,j) = 0.
        ENDDO
      ENDDO
!     
!     OUTER LOOP OVER VERTICAL DIMENSION.
!     INNER LOOP OVER HORIZONTAL GRID.
!     
!!$omp  parallel do private(i,j,l,es,dp)
      DO L = 1,LM
        IF (IDECID <= 1) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = Q(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 2) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = QQW(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 3) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = QQI(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 4) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = QQR(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 5) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = QQS(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 6) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = CWM(I,J,L)
            ENDDO
          ENDDO
! SRD
        ELSE IF (IDECID == 16) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = QQG(I,J,L)
            ENDDO
          ENDDO
! SRD
        ELSE IF (IDECID == 7) THEN
!-- Total supercooled liquid
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              IF (T(I,J,L) >= TFRZ) THEN
                Qdum(I,J) = 0.
              ELSE
                Qdum(I,J) = QQW(I,J,L) + QQR(I,J,L)
              ENDIF
            ENDDO
          ENDDO
        ELSE IF (IDECID == 8) THEN
!-- Total melting ice
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              IF (T(I,J,L) <= TFRZ) THEN
                Qdum(I,J) = 0.
              ELSE
                Qdum(I,J) = QQI(I,J,L) + QQS(I,J,L)
              ENDIF
            ENDDO
          ENDDO
        ELSE IF (IDECID == 9) THEN
! SHORT WAVE T TENDENCY
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = RSWTT(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 10) THEN
! LONG WAVE T TENDENCY
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = RLWTT(I,J,L)
            ENDDO
          ENDDO  
        ELSE IF (IDECID == 11) THEN
! LATENT HEATING FROM GRID SCALE RAIN/EVAP
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = TRAIN(I,J,L)
            ENDDO
          ENDDO  
        ELSE IF (IDECID == 12) THEN
! LATENT HEATING FROM CONVECTION
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = TCUCN(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 13) THEN
! MOISTURE CONVERGENCE
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = MCVG(I,J,L)
            ENDDO
          ENDDO
! RH
        ELSE IF (IDECID == 14) THEN
!$omp  parallel do private(i,j,es)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = Q(I,J,L)
              ES        = min(FPVSNEW(T(I,J,L)),PMID(I,J,L))
              QS(I,J)   = CON_EPS*ES/(PMID(I,J,L)+CON_EPSM1*ES)
            ENDDO
          END DO
! OZONE
        ELSE IF (IDECID == 15) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = O3(I,J,L)
            ENDDO
          END DO

! AEROSOL EXTINCTION (GOCART)
        ELSE IF (IDECID == 17) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = EXT(I,J,L)
            ENDDO
          END DO
!
! E. James - 8 Dec 2017
! FIRE SMOKE (tracer_1a FROM HRRR-SMOKE)
        ELSE IF (IDECID == 18) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = SMOKE(I,J,L,1)/1000000000.
            ENDDO
          END DO
!
! E. James - 8 Dec 2017
! HRRR-SMOKE AOD
        ELSE IF (IDECID == 19) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = TAOD5503D(I,J,L)
            ENDDO
          END DO
!LZhang -July 2019
! SCATTERING AEROSOL OPTICAL THICKNESS (GOCART V2)
        ELSE IF (IDECID == 20) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = SCA(I,J,L)
            ENDDO
          END DO

! ASYMMETRY PARAMETER (GOCART V2)
        ELSE IF (IDECID == 21) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              Qdum(I,J) = ASY(I,J,L)
            ENDDO
          END DO
        ENDIF
!
!$omp  parallel do private(i,j,dp)
        DO J=JSTA,JEND
          DO I=1,IM
            if(PINT(I,J,L+1) <spval .and. Qdum(I,J) < spval) then
             DP      = PINT(I,J,L+1) - PINT(I,J,L)
             PW(I,J) = PW(I,J) + Qdum(I,J)*DP*GI*HTM(I,J,L)
            IF (IDECID == 17 .or. IDECID == 20 .or. IDECID == 21) THEN
             PW(I,J) = PW(I,J) + Qdum(I,J)*MAX(DP,0.)*GI*HTM(I,J,L)
            ENDIF
            IF (IDECID == 19) THEN
             PW(I,J) = PW(I,J) + Qdum(I,J)
            ENDIF
            IF (IDECID == 14) PWS(I,J) = PWS(I,J) + QS(I,J)*DP*GI*HTM(I,J,L)
            else
             PW(I,J) = spval
             PWS(I,J) = spval
            endif
          ENDDO
        ENDDO
      ENDDO                 ! l loop

      
      IF (IDECID == 14)THEN
!$omp  parallel do private(i,j,dp)
        DO J=JSTA,JEND
          DO I=1,IM
            if( PW(I,J)<spval) then
            PW(I,J) = max(0.,PW(I,J)/PWS(I,J)*100.) 
            endif
          ENDDO
        ENDDO
      END IF
!  convert ozone from kg/m2 to dobson units, which give the depth of the
!  ozone layer in 1e-5 m if brought to natural temperature and pressure.    

      IF (IDECID == 15) then
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            if( PW(I,J)<spval) then
            PW(I,J) = PW(I,J) / 2.14e-5
            endif
          ENDDO
        ENDDO
      endif
!
!     END OF ROUTINE.
!     
      RETURN
      END
