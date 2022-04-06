!> @file
!> @brief lfmfld() computes layer mean LFM fields.
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
!> 2019-10-30 | Bo Cui       | Remove "GOTO" statement
!> 2020-11-10 | Jesse Meng   | Use UPP_PHYSICS Module
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE LFMFLD(RH3310,RH6610,RH3366,PW3310)

!     
!
      use vrbls3d, only: pint, alpint, zint, t, q, cwm
      use masks, only: lmh
      use params_mod, only: d00, d50, pq0, a2, a3, a4, h1, d01, gi
      use ctlblk_mod, only: jsta, jend, modelname, spval, im
      use physcons_post, only: con_rd, con_rv, con_eps, con_epsm1
      use upp_physics, only: FPVSNEW

      implicit none

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      real,PARAMETER :: RHOWAT=1.E3
!     
!     DECLARE VARIABLES.
!     
      REAL ALPM, DZ, ES, PM, PWSUM, QM, QS, TM, DP, RH
      REAL,dimension(IM,jsta:jend),intent(inout) :: RH3310, RH6610, RH3366
      REAL,dimension(IM,jsta:jend),intent(inout) :: PW3310
      real Z3310,Z6610,Z3366,P10,P33,P66
      integer I,J,L,LLMH
!
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
         RH3310(I,J) = D00
         PW3310(I,J) = D00
         RH6610(I,J) = D00
         RH3366(I,J) = D00
         Z3310     = D00
         Z6610     = D00
         Z3366     = D00
!     
!        SET BOUNDS FOR PRESSURES AND SURFACE L.
         P10  = PINT(I,J,NINT(LMH(I,J)))
         P66  = 0.75*P10
         P33  = 0.50*P10
         LLMH = NINT(LMH(I,J))
!     
!        ACCULMULATE RELATIVE HUMIDITIES AND PRECIPITABLE WATER.
!
         DO 10 L = LLMH,1,-1
!     
!           GET P, Z, T, AND Q AT MIDPOINT OF ETA LAYER.
            ALPM = D50*(ALPINT(I,J,L)+ALPINT(I,J,L+1))
            DZ   = ZINT(I,J,L)-ZINT(I,J,L+1)
            DP   = PINT(I,J,L+1)-PINT(I,J,L)
            PM   = EXP(ALPM)
            TM   = T(I,J,L)
            QM   = Q(I,J,L)
            QM   = AMAX1(QM,D00)
!
!            QS=PQ0/PM*EXP(A2*(TM-A3)/(TM-A4))
	    IF(MODELNAME == 'GFS')THEN
	      ES = min(FPVSNEW(TM),PM)
	      QS = CON_EPS*ES/(PM+CON_EPSM1*ES)
	    ELSE      
              QS=PQ0/PM*EXP(A2*(TM-A3)/(TM-A4))
	    END IF
            RH   = QM/QS
            IF (RH>H1) THEN
               RH = H1
               QM = RH*QS
            ENDIF
            IF (RH<D01) THEN
               RH = D01
               QM = RH*QS
            ENDIF
!
!           JUMP OUT OF THIS LOOP IF WE ARE ABOVE THE HIGHEST TARGET PRESSURE.
            IF (PM<=P33) exit     
!     
!           0.66-1.00 RELATIVE HUMIDITY.
            IF ((PM<=P10).AND.(PM>=P66)) THEN
               Z6610     = Z6610 + DZ
               RH6610(I,J) = RH6610(I,J) + RH*DZ
            ENDIF
!     
!           0.33-1.00 RELATIVE HUMIDITY AND PRECIPITABLE WATER.
            IF ((PM<=P10).AND.(PM>=P33)) THEN
               Z3310      = Z3310 + DZ
               RH3310(I,J)= RH3310(I,J)+RH*DZ
               PW3310(I,J)= PW3310(I,J)+(Q(I,J,L)+CWM(I,J,L))*DP*GI
            ENDIF
!     
!           0.33-0.66 RELATIVE HUMIDITY.
            IF ((PM<=P66).AND.(PM>=P33)) THEN
               Z3366     = Z3366 + DZ
               RH3366(I,J) = RH3366(I,J) + RH*DZ
            ENDIF
!
 10      CONTINUE
!     
!        NORMALIZE TO GET MEAN RELATIVE HUMIDITIES.  AT
!        ONE TIME WE DIVIDED PRECIPITABLE WATER BY DENSITY
!        TO GET THE EQUIVALENT WATER DEPTH IN METERS.  NO MORE.
         IF (Z6610>D00) THEN
            RH6610(I,J) = RH6610(I,J)/Z6610
         ELSE
            RH6610(I,J) = SPVAL
         ENDIF
!     
         IF (Z3310>D00) THEN
            RH3310(I,J) = RH3310(I,J)/Z3310
         ELSE
            RH3310(I,J) = SPVAL
         ENDIF
!     
         IF (Z3366>D00) THEN
            RH3366(I,J) = RH3366(I,J)/Z3366
         ELSE
            RH3366(I,J) = SPVAL
         ENDIF
 30   CONTINUE
!     
!     
!     END OF ROUTINE.
!     
      RETURN
      END
