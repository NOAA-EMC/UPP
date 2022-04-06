!> @file
!> @brief LFMFLD_GFS computes layer mean LFM fields.
!>
!> This routine computes three layer mean relative humidities
!> and a precipitable water field from ETA level data.  The
!> computed fields are intended to mimic similar fields com-
!> puted by the LFM.  The algorithm used here is fairly pri-
!> mative.
!> <pre>
!> In each column above a mass point on the ETA grid we set the following target pressures:
!>     Sigma layer 1.00 pressure:  Surface pressure
!>     Sigma layer 0.66 pressure:  0.50 * Surface pressure
!>     Sigma layer 0.33 pressure:  0.4356 * Surface pressure
!> </pre>
!> Given there pressures a surface up summation is made of
!> relative humidity and/or precipitable water between these
!> target pressures.  Each term in the summation is weighted
!> By the thickness of the ETA layer.  The final layer mean
!> is this sum normalized by the total depth of the layer. 
!> There is, obviously, no normalization for precipitable water.
!>
!> @param[out] RH3310 Sigma layer 0.33-1.00 mean relative humidity.
!> @param[out] RH6610 Sigma layer 0.66-1.00 mean relative humidity.
!> @param[out] RH3366 Sigma layer 0.33-0.66 mean relative humidity.
!> @param[out] PW3310 Sigma layer 0.33-1.00 precipitable water.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-22 | Russ Treadon | Initial
!> 1993-07-27 | Russ Treadon | Modified summation limits from 0.66*PSFC to 0.75*PSFC and 0.33*PSFC to 0.50*PSFC, where PSFC is the surfaces pressure.  The reason for this change was recognition that in the LFM 0.33 and 0.66 were measured from the surface to the tropopause not the top of the model.
!> 1993-09-13 | Russ Treadon | RH calculations were made internal to the routine.
!> 1996-03-04 | Mike Baldwin | Change PW CALC to include CLD WTR
!> 1998-06-16 | T Black      | Conversion from 1-D to 2-D
!> 1998-08-17 | Mike Baldwin | Compute RH over ice
!> 1998-12-22 | Mike Baldwin | Back out RH over ice
!> 2000-01-04 | Jim Tuccillo | MPI Version
!> 2002-04-24 | Mike Baldwin | WRF Version
!> 2006-11-06 | H CHUANG     | Modify to output GFS LFM fields which have different thickness as MESO and use DP rather than DZ
!> 2019-10-30 | Bo Cui       | Remove "GOTO" statement
!> 2020-11-10 | Jesse Meng   | Use UPP_PHYSICS Module
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE LFMFLD_GFS(RH4410,RH7294,RH4472,RH3310)

!     
!
      use vrbls3d, only: pint, q, t, pmid
      use masks, only: lmh
      use params_mod, only: d00
      use ctlblk_mod, only: jsta, jend, spval, im
      use upp_physics, only: FPVSNEW
!     
    implicit none
!
      real,PARAMETER :: RHOWAT=1.E3
      real,parameter:: con_rd      =2.8705e+2 ! gas constant air    (J/kg/K)
      real,parameter:: con_rv      =4.6150e+2 ! gas constant H2O
      real,parameter:: con_eps     =con_rd/con_rv
      real,parameter:: con_epsm1   =con_rd/con_rv-1
      real,parameter:: strh1=0.44,strh2=0.72,strh3=0.44,strh4=0.33 &
                      ,sbrh1=1.00,sbrh2=0.94,sbrh3=0.72,sbrh4=1.00
!     
!     DECLARE VARIABLES.
!     
      REAL ALPM, DZ, ES, PM, PWSUM, QM, QS
      REAL,dimension(IM,jsta:jend),intent(out) :: RH4410, RH7294, RH4472    &
                                                 ,RH3310    
!
      integer I,J,L,LLMH
      real P4410, P7294,P4472,P3310,Q4410,Q7294,Q4472,Q3310,QS4410, &
         QS7294,QS4472,QS3310,PS,P33,DP1,DP2,DP3,DP4

!***********************************************************************
!     START LFMFLD HERE
!     
!
!     LOOP OVER HORIZONTAL GRID.
!     
      DO 30 J=JSTA,JEND
      DO 30 I=1,IM
!     
!        ZERO VARIABLES.
         RH4410(I,J) = D00
         RH4472(I,J) = D00
         RH7294(I,J) = D00
         RH3310(I,J) = D00
         P4410       = D00
         P7294       = D00
         P4472       = D00
         P3310       = D00
         Q4410       = D00
         Q7294       = D00
         Q4472       = D00
         Q3310       = D00
         QS4410      = D00
         QS7294      = D00
         QS4472      = D00
         QS3310      = D00
!     
!        SET BOUNDS FOR PRESSURES AND SURFACE L.
 
         LLMH = NINT(LMH(I,J))
	 PS=PINT(I,J,LLMH+1)
	 P33  = 0.33*PS
!     
!        ACCULMULATE RELATIVE HUMIDITIES AND PRECIPITABLE WATER.
!
         DO 10 L = LLMH,1,-1
!     
!           GET P, Z, T, AND Q AT MIDPOINT OF ETA LAYER.
            
            DP1 = MAX(MIN(PINT(I,J,L+1),SBRH1*PS)      &
                     -MAX(PINT(I,J,L),STRH1*PS),0.)            
            DP2 = MAX(MIN(PINT(I,J,L+1),SBRH2*PS)      &
                     -MAX(PINT(I,J,L),STRH2*PS),0.)
            DP3 = MAX(MIN(PINT(I,J,L+1),SBRH3*PS)      &
                     -MAX(PINT(I,J,L),STRH3*PS),0.)
            DP4 = MAX(MIN(PINT(I,J,L+1),SBRH4*PS)      &
                     -MAX(PINT(I,J,L),STRH4*PS),0.)
     
            PM   = PINT(I,J,L)
            QM   = Q(I,J,L)
            QM   = MAX(QM,D00)
            ES   = min(FPVSNEW(T(I,J,L)),PMID(I,J,L))
            QS=CON_EPS*ES/(PMID(I,J,L)+CON_EPSM1*ES)
!
!
!           JUMP OUT OF THIS LOOP IF WE ARE ABOVE THE HIGHEST TARGET PRESSURE.
            IF (PM<=P33) exit
!     
!           0.44-1.00 RELATIVE HUMIDITY.
!            IF ((PM<=P10).AND.(PM>=P44)) THEN
               P4410     = P4410 + DP1
               Q4410     = Q4410 + QM*DP1
               QS4410    = QS4410+ QS*DP1
!            ENDIF
!     
!           0.33-1.00 RELATIVE HUMIDITY 
!            IF ((PM<=P10).AND.(PM>=P33)) THEN
               P3310      = P3310 + DP4
               Q3310     = Q3310 + QM*DP4
               QS3310    = QS3310+ QS*DP4
!            ENDIF
!     
!           0.44-0.72 RELATIVE HUMIDITY.
!            IF ((PM<=P66).AND.(PM>=P33)) THEN
               P4472     = P4472 + DP3
               Q4472     = Q4472 + QM*DP3
               QS4472    = QS4472+ QS*DP3
!            ENDIF
!           0.72-0.94 RELATIVE HUMIDITY.
!            IF ((PM<=P66).AND.(PM>=P33)) THEN
               P7294     = P7294 + DP2
               Q7294     = Q7294 + QM*DP2
               QS7294    = QS7294+ QS*DP2
!            ENDIF
!
 10      CONTINUE
!     
!        NORMALIZE TO GET MEAN RELATIVE HUMIDITIES.  AT
!        ONE TIME WE DIVIDED PRECIPITABLE WATER BY DENSITY
!        TO GET THE EQUIVALENT WATER DEPTH IN METERS.  NO MORE.
         IF (P4410>D00) THEN
            RH4410(I,J) = Q4410/QS4410
         ELSE
            RH4410(I,J) = SPVAL
         ENDIF
!     
         IF (P3310>D00) THEN
            RH3310(I,J) = Q3310/QS3310
         ELSE
            RH3310(I,J) = SPVAL
         ENDIF
!     
         IF (P4472>D00) THEN
            RH4472(I,J) = Q4472/QS4472
         ELSE
            RH4472(I,J) = SPVAL
         ENDIF

         IF (P7294>D00) THEN
            RH7294(I,J) = Q7294/QS7294
         ELSE
            RH7294(I,J) = SPVAL
         ENDIF
 30   CONTINUE
!     
!     END OF ROUTINE.
!     
      RETURN
      END
