!> @file
!> @brief Subroutine that computes precipitable water.
!>
!><pre>
!> This routine computes precipitable water in a column
!> extending from the first atmospheric ETA layer to the
!> model top. The definition used is
!>                      TOP
!> precipitable water = sum (Q+CLDW) DP*HTM/G
!>                      BOT
!> where,
!> BOT is the first ETA layer,
!> TOP is the model top,
!> Q is the specific humidity (kg/kg) in the layer
!> CLDW is the cloud water (kg/kg) in the layer
!> DP (Pa) is the layer thickness.
!> HTM is the height mask at that layer (=0 if below ground)
!> G is the gravitational constant.
!></pre>
!>     
!> @param[in] PW Array of precipitable water.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-24 | Russ Treadon   | Initial
!> 1996-03-04 | Mike Baldwin   | Add cloud water and speed up code
!> 1998-06-15 | T Black        | Convesion from 1-D to 2-D
!> 2000-01-04 | Jim Tuccillo   | MPI Version            
!> 2002-06-19 | Mike Baldwin   | WRF Version 
!> 2004-12-30 | H Chuang       | Update to calculate total column for other hydrometeors
!> 2014-11-12 | Sarah Lu       | Update tp calculate aerosol optical depth
!> 2015-07-02 | Sarah Lu       | Update to calculate scattering aerosal optical depth (18)
!> 2015-07-04 | Sarah Lu       | Correct PW integration for AOD (17)
!> 2015-07-10 | Sarah Lu       | Update to calculate asymetry parameter
!> 2019-07-25 | Li(Kate) Zhang | Merge Sarah Lu's update for FV3-Chem
!> 2020-11-10 | Jesse Meng     | Use UPP_PHYSICS Module
!> 2021-09-02 | Bo Cui         | Decompose UPP in X direction          
!> 2022-11-16 | Eric James     | Adding calculation of vertically integrated dust from RRFS
!> 2023-02-23 | Eric James     | Adding vertically integrated coarse PM from RRFS
!>     
!> @author Russ Treadon W/NP2 @date 1992-12-24
      SUBROUTINE CALPW(PW,IDECID)

!     
      use vrbls3d,    only: q, qqw, qqi, qqr, qqs, cwm, qqg, t, rswtt,    &
                            train, tcucn, mcvg, pmid, o3, ext, pint, rlwtt, &
                            taod5503d,sca, asy
      use vrbls4d,    only: smoke, fv3dust, coarsepm
      use masks,      only: htm
      use params_mod, only: tfrz, gi
      use ctlblk_mod, only: lm, jsta, jend, im, spval, ista, iend
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
      real,dimension(ista:iend,jsta:jend),intent(inout) :: PW
      INTEGER LLMH,I,J,L
      REAL ALPM,DZ,PM,PWSUM,RHOAIR,DP,ES
      REAL QDUM(ista:iend,jsta:jend), PWS(ista:iend,jsta:jend),QS(ista:iend,jsta:jend)
!
!***************************************************************
!     START CALPW HERE.
!
!     INITIALIZE PW TO 0.    
!     
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
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
            DO I=ISTA,IEND
              Qdum(I,J) = Q(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 2) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = QQW(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 3) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = QQI(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 4) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = QQR(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 5) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = QQS(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 6) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = CWM(I,J,L)
            ENDDO
          ENDDO
! SRD
        ELSE IF (IDECID == 16) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = QQG(I,J,L)
            ENDDO
          ENDDO
! SRD
        ELSE IF (IDECID == 7) THEN
!-- Total supercooled liquid
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
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
            DO I=ISTA,IEND
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
            DO I=ISTA,IEND
              Qdum(I,J) = RSWTT(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 10) THEN
! LONG WAVE T TENDENCY
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = RLWTT(I,J,L)
            ENDDO
          ENDDO  
        ELSE IF (IDECID == 11) THEN
! LATENT HEATING FROM GRID SCALE RAIN/EVAP
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = TRAIN(I,J,L)
            ENDDO
          ENDDO  
        ELSE IF (IDECID == 12) THEN
! LATENT HEATING FROM CONVECTION
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = TCUCN(I,J,L)
            ENDDO
          ENDDO
        ELSE IF (IDECID == 13) THEN
! MOISTURE CONVERGENCE
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = MCVG(I,J,L)
            ENDDO
          ENDDO
! RH
        ELSE IF (IDECID == 14) THEN
!$omp  parallel do private(i,j,es)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = Q(I,J,L)
              ES        = min(FPVSNEW(T(I,J,L)),PMID(I,J,L))
              QS(I,J)   = CON_EPS*ES/(PMID(I,J,L)+CON_EPSM1*ES)
            ENDDO
          END DO
! OZONE
        ELSE IF (IDECID == 15) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = O3(I,J,L)
            ENDDO
          END DO

! AEROSOL EXTINCTION (GOCART)
        ELSE IF (IDECID == 17) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = EXT(I,J,L)
            ENDDO
          END DO
!
! E. James - 8 Dec 2017
! FIRE SMOKE (tracer_1a FROM HRRR-SMOKE)
        ELSE IF (IDECID == 18) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = SMOKE(I,J,L,1)/(1E9)
            ENDDO
          END DO
!
! E. James - 8 Dec 2017
! HRRR-SMOKE AOD
        ELSE IF (IDECID == 19) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = TAOD5503D(I,J,L)
            ENDDO
          END DO
!LZhang -July 2019
! SCATTERING AEROSOL OPTICAL THICKNESS (GOCART V2)
        ELSE IF (IDECID == 20) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = SCA(I,J,L)
            ENDDO
          END DO

! ASYMMETRY PARAMETER (GOCART V2)
        ELSE IF (IDECID == 21) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = ASY(I,J,L)
            ENDDO
          END DO

! E. James - 14 Sep 2022
! DUST (from RRFS)
        ELSE IF (IDECID == 22) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = FV3DUST(I,J,L,1)/(1E9)
            ENDDO
          END DO

! E. James - 23 Feb 2023
! COARSEPM (from RRFS)
        ELSE IF (IDECID == 23) THEN
!$omp  parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              Qdum(I,J) = COARSEPM(I,J,L,1)/(1E9)
            ENDDO
          END DO
        ENDIF
!
!$omp  parallel do private(i,j,dp)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            if(PINT(I,J,L+1) <spval .and. Qdum(I,J) < spval) then
             DP      = PINT(I,J,L+1) - PINT(I,J,L)
            IF (IDECID == 19) THEN
             PW(I,J) = PW(I,J) + Qdum(I,J)
	    ELSE
	     PW(I,J) = PW(I,J) + Qdum(I,J)*MAX(DP,0.)*GI*HTM(I,J,L)
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
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
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
          DO I=ISTA,IEND
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
