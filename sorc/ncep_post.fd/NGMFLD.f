!> @file
!> @brief ngmfld() computes layer mean NGM fields
!>
!> This routine computes a handful of NGM layer mean 
!> fields.  This is done to provide a fully complete 
!> ETA NGM look-alike output file.  The sigma (layer)
!> fields computed bu this routine are tabulated below.
!><pre>
!>       Sigma (layer)         Field(s)
!>      ---------------     --------------
!>      0.47191-1.00000          RH
!>      0.47171-0.96470          RH
!>      0.18019-0.47191          RH
!>      0.84368-0.98230          RH
!>      0.85000-1.00000         MCONV
!> where
!>      RH    = Relative humidity
!>      MCONV = Moisture convergence
!></pre>
!> Layer means are a summation over ETA layers mapping into
!> The pressure range corresponding to the sigma range above.
!> The calculation of these bounding pressures is done at
!> each horizontal grid point based on the surface pressure.
!> Each term in the summation is weighted by the thickness of
!> the ETA layer.  The final layer mean is this sum normalized
!> by the total depth of the layer.
!>
!> @param[out] RH4710 Sigma layer 0.47-1.00 mean relative humidity.
!> @param[out] RH4796 Sigma layer 0.47-0.96 mean relative humidity.
!> @param[out] RH1847 Sigma layer 0.18-0.47 mean relative humidity.
!> @param[out] RH8498 Sigma layer 0.84-0.98 mean relative humidity.
!> @param[out] QM8510 Sigma layer 0.85-1.00 mean moisture convergence.
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
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE NGMFLD(RH4710,RH4796,RH1847,RH8498,QM8510)

!     
!     
!     INCLUDE PARAMETERS
      use vrbls3d,    only: q, uh, vh, pint, alpint, zint, t
      use masks,      only: lmh
      use params_mod, only: d00, d50, h1m12, pq0, a2, a3, a4, h1, d01, small
      use ctlblk_mod, only: jsta, jend, lm, jsta_2l, jend_2u, jsta_m2, jend_m2,&
                            spval, im
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      real,PARAMETER :: SIG100=1.00000, SIG98=0.98230, SIG96=0.96470
      real,PARAMETER :: SIG89 =0.89671, SIG85=0.85000, SIG84=0.84368
      real,PARAMETER :: SIG78 =0.78483, SIG47=0.47191, SIG18=0.18018
!     
!     DECLARE VARIABLES.
      LOGICAL GOT8510,GOT4710,GOT4796,GOT1847,GOT8498
      REAL,dimension(IM,jsta_2l:jend_2u),intent(out) :: QM8510,RH4710,RH8498, &
                                                        RH4796,RH1847
      REAL,dimension(im,jsta_2l:jend_2u) :: Z8510,Z4710,Z8498,Z4796,Z1847
      real,dimension(im,jsta_2l:jend_2u) ::  Q1D, U1D, V1D, QCNVG
!
      integer I,J,L
      real P100,P85,P98,P96,P84,P47,P18,ALPM,DE,PM,TM,QM,     &
           QMCVG,QS,RH,DZ
!********************************************************************
!     START NGMFLD HERE.
!     
!     INITIALIZE ARRAYS.
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=1,IM
           QM8510(I,J) = D00
           RH4710(I,J) = D00
           RH8498(I,J) = D00
           RH4796(I,J) = D00
           RH1847(I,J) = D00
           Z8510(I,J)  = D00
           Z8498(I,J)  = D00
           Z4710(I,J)  = D00
           Z4796(I,J)  = D00
           Z1847(I,J)  = D00
        ENDDO
      ENDDO
!     
!     LOOP OVER HORIZONTAL GRID.
!     
!!$omp  parallel do                                                 &
!     & private(dz,p100,p18,p47,p84,p85,                            &
!     &         p96,p98,pm,qdiv,qk,qkhn,qkhs,qkm1,qm,qm8510,        &
!     &         qmcvg,qs,qudx,qvdy,r2dx,r2dy,rh,rh1847,rh4710,      &
!     &         rh4796,rh8498,tm,tmt0,tmt15,z1847,z4710,z4796,      &
!     &         z8498,z8510,q1d,u1d,v1d,qcnvg)

     DO L=1,LM
!          COMPUTE MOISTURE CONVERGENCE
!$omp parallel do private(i,j)
       DO J=JSTA_2L,JEND_2U
         DO I=1,IM
           Q1D(I,J) = Q(I,J,L)
           U1D(I,J) = UH(I,J,L)
           V1D(I,J) = VH(I,J,L)
         ENDDO
       ENDDO
       CALL CALMCVG(Q1D,U1D,V1D,QCNVG)
!          COMPUTE MOISTURE CONVERGENCE
      DO J=JSTA_M2,JEND_M2
      DO I=2,IM-1
!
!        SET TARGET PRESSURES.
         
         P100  = PINT(I,J,NINT(LMH(I,J)))
         P98   = SIG98*P100
         P96   = SIG96*P100
         P85   = SIG85*P100
         P84   = SIG84*P100
         P47   = SIG47*P100
         P18   = SIG18*P100
!     
!     
!        COMPUTE LAYER MEAN FIELDS AT THE GIVEN K.
!
!          COMPUTE P, Z, T, AND Q AT THE MIDPOINT OF THE CURRENT ETA LAYER.
           ALPM = D50*(ALPINT(I,J,L)+ALPINT(I,J,L+1))
           DZ   = ZINT(I,J,L)-ZINT(I,J,L+1)
           PM   = EXP(ALPM)
           TM   = T(I,J,L)
           QM   = Q(I,J,L)
           QM   = AMAX1(QM,H1M12)
           QMCVG= QCNVG(I,J)
!
!     
!          COMPUTE RELATIVE HUMIDITY.
!
           QS=PQ0/PM*EXP(A2*(TM-A3)/(TM-A4))
!
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
!          SIGMA 0.85-1.00 MOISTURE CONVERGENCE.
           IF ((PM<=P100).AND.(PM>=P85)) THEN
              Z8510(I,J)  = Z8510(I,J) + DZ
              QM8510(I,J) = QM8510(I,J) + QMCVG*DZ
           ENDIF
!    
!          SIGMA 0.47-1.00 RELATIVE HUMIDITY.
           IF ((PM<=P100).AND.(PM>=P47)) THEN
              Z4710(I,J)  = Z4710(I,J) + DZ
              RH4710(I,J) = RH4710(I,J) + RH*DZ
           ENDIF
!
!          SIGMA 0.84-0.98 RELATIVE HUMIDITY.
           IF ((PM<=P98).AND.(PM>=P84)) THEN
              Z8498(I,J)  = Z8498(I,J) + DZ
              RH8498(I,J) = RH8498(I,J) + RH*DZ
           ENDIF
!     
!          SIGMA 0.47-0.96 RELATIVE HUMIDITY.
           IF ((PM<=P96).AND.(PM>=P47)) THEN
              Z4796(I,J)  = Z4796(I,J) + DZ
              RH4796(I,J) = RH4796(I,J) + RH*DZ
           ENDIF
!     
!          SIGMA 0.18-0.47 RELATIVE HUMIDITY.
           IF ((PM<=P47).AND.(PM>=P18)) THEN
              Z1847(I,J)  = Z1847(I,J) + DZ
              RH1847(I,J) = RH1847(I,J) + RH*DZ
           ENDIF
!
      ENDDO
      ENDDO
      ENDDO
!     
      DO J=JSTA_M2,JEND_M2
      DO I=2,IM-1
!        NORMALIZE TO GET LAYER MEAN VALUES.
         IF (Z8510(I,J)>0) THEN
            QM8510(I,J) = QM8510(I,J)/Z8510(I,J)
         ELSE
            QM8510(I,J) = SPVAL
         ENDIF
         IF (ABS(QM8510(I,J)-SPVAL)<SMALL)QM8510(I,J)=H1M12
!
         IF (Z4710(I,J)>0) THEN
            RH4710(I,J) = RH4710(I,J)/Z4710(I,J)
         ELSE
            RH4710(I,J) = SPVAL
         ENDIF
!
         IF (Z8498(I,J)>0) THEN
            RH8498(I,J) = RH8498(I,J)/Z8498(I,J)
         ELSE
            RH8498(I,J) = SPVAL
         ENDIF
!
         IF (Z4796(I,J)>0) THEN
            RH4796(I,J) = RH4796(I,J)/Z4796(I,J)
         ELSE
            RH4796(I,J) = SPVAL
         ENDIF
!
         IF (Z1847(I,J)>0) THEN
            RH1847(I,J) = RH1847(I,J)/Z1847(I,J)
         ELSE
            RH1847(I,J) = SPVAL
         ENDIF
      ENDDO
      ENDDO
!
!     
!     END OF ROUTINE.
!     
      RETURN
      END

