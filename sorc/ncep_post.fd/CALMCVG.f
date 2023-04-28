!> @file
!> @brief Subroutine that computes moisture convergence.
!>
!><pre>
!> Given specific humidity, Q, and the U-V wind components
!> This routine evaluates the vector operation, 
!>                  DEL DOT (Q*VEC)
!> where,
!>    DEL is the vector gradient operator,
!>    DOT is the standard dot product operator, and
!>    VEC is the vector wind.
!> Minus one times the resulting scalar field is the 
!> moisture convergence which is returned by this routine.
!></pre>
!>   
!> @param[in] Q1D real Specific humidity at P-points (kg/kg).
!> @param[in] U1D real U-wind component (m/s) at P-points.
!> @param[in] V1D real V-wind component (m/s) at P-points.
!> @param[out] QCNVG real Moisture convergence (1/s) at P-points.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1993-01-22 | Russ Treadon | Initial
!> 1998-06-08 | T Black      | Conversion From 1-D To 2-D
!> 2000-01-04 | Jim Tuccillo | MPI Version              
!> 2002-04-23 | Mike Baldwin | WRF C-Grid Version     
!> 2005-07-07 | Binbin Zhou  | Add RSM A Grid
!> 2006-04-25 | H Chuang     | Bug fixes to correctly compute MC at boundaries 
!> 2021-04-01 | J Meng       | Computation on defined points only
!> 2021-09-02 | B CUI        | REPLACE EXCH_F to EXCH             
!> 2021-09-30 | J MENG       | 2D DECOMPOSITION
!>     
!> @author Russ Treadon W/NP2 @date 1993-01-22
!-----------------------------------------------------------------------
!> @brief Subroutine that computes moisture convergence.
!>   
!> @param[in] Q1D real Specific humidity at P-points (kg/kg).
!> @param[in] U1D real U-wind component (m/s) at P-points.
!> @param[in] V1D real V-wind component (m/s) at P-points.
!> @param[out] QCNVG real Moisture convergence (1/s) at P-points.
!-----------------------------------------------------------------------
      SUBROUTINE CALMCVG(Q1D,U1D,V1D,QCNVG)

!
!     
!     
      use masks,        only: dx, dy, hbm2
      use params_mod,   only: d00, d25
      use ctlblk_mod,   only: jsta_2l, jend_2u, spval, jsta_m, jend_m,       &
                              jsta_m2, jend_m2, im, jm,                      &
                              ista_2l, iend_2u, ista_m, iend_m, ista_m2, iend_m2 
      use gridspec_mod, only: gridtype
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     DECLARE VARIABLES.
!     
      REAL,dimension(ista_2l:iend_2u,jsta_2l:jend_2u),intent(in)    ::  Q1D, U1D, V1D
      REAL,dimension(ista_2l:iend_2u,jsta_2l:jend_2u),intent(inout) ::  QCNVG

      REAL R2DY, R2DX
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u) ::  UWND, VWND, QV
      INTEGER IHE(JM),IHW(JM),IVE(JM),IVW(JM)
      integer I,J,ISTA2,IEND2
      real QVDY,QUDX
!     
!***************************************************************************
!     START CALMCVG HERE.
!
!     
!     INITIALIZE MOISTURE CONVERGENCE ARRAY.  LOAD TEMPORARY WIND ARRAYS.
!     
       CALL EXCH(Q1D)
       CALL EXCH(U1D)
       CALL EXCH(V1D)

!$omp  parallel do private(i,j)
      DO J=JSTA_2L,JEND_2U
!        DO I=1,IM
        DO I=ISTA_2L,IEND_2U
          IF(U1D(I,J)<SPVAL)THEN
          QCNVG(I,J) = 0.
          ELSE
          QCNVG(I,J) = SPVAL
          ENDIF
          UWND(I,J)  = U1D(I,J)
          VWND(I,J)  = V1D(I,J)
          IF (UWND(I,J) == SPVAL) UWND(I,J) = D00
          IF (VWND(I,J) == SPVAL) VWND(I,J) = D00
        ENDDO
      ENDDO
!
      IF(gridtype == 'A')THEN
!$omp  parallel do private(i,j,qudx,qvdy,r2dx,r2dy)
       DO J=JSTA_M,JEND_M
!         DO I=2,IM-1
         DO I=ISTA_M,IEND_M
           IF(Q1D(I,J+1)<SPVAL.AND.Q1D(I,J-1)<SPVAL.AND.          &
              Q1D(I+1,J)<SPVAL.AND.Q1D(I-1,J)<SPVAL.AND. &
              Q1D(I,J)<SPVAL) THEN
             R2DX   = 1./(2.*DX(I,J))   !MEB DX?
             R2DY   = 1./(2.*DY(I,J))   !MEB DY?  
             QUDX   = (Q1D(I+1,J)*UWND(I+1,J)-Q1D(I-1,J)*UWND(I-1,J))*R2DX
             QVDY   = (Q1D(I,J+1)*VWND(I,J+1)-Q1D(I,J-1)*VWND(I,J-1))*R2DY
             QCNVG(I,J) = -(QUDX + QVDY)
           ELSE
             QCNVG(I,J) = SPVAL
           ENDIF
         ENDDO
       ENDDO
      ELSE IF(gridtype == 'E')THEN

       DO J=JSTA_M,JEND_M
         IHE(J) = MOD(J+1,2)
         IHW(J) = IHE(J)-1
         IVE(J) = MOD(J,2)
         IVW(J) = IVE(J)-1 
       END DO
     
!$omp  parallel do private(i,j)
       DO J=JSTA_M,JEND_M
!         ISTA = 1+MOD(J+1,2)
!         IEND = IM-MOD(J,2)
!         DO I=ISTA,IEND
         DO I=ISTA_M,IEND_M
          IF(Q1D(I,J-1)<SPVAL.AND.Q1D(I+IVW(J),J)<SPVAL.AND.&
             Q1D(I+IVE(J),J)<SPVAL.AND.Q1D(I,J+1)<SPVAL) THEN
           QV(I,J) = D25*(Q1D(I,J-1)+Q1D(I+IVW(J),J)                   &
                         +Q1D(I+IVE(J),J)+Q1D(I,J+1))
          ELSE
           QV(I,J) = SPVAL
          ENDIF
         END DO
       END DO

       CALL EXCH(QV)
!      CALL EXCH(VWND)

!
!$omp  parallel do private(i,j,qudx,qvdy,r2dx,r2dy)
       DO J=JSTA_M2,JEND_M2
!         IEND = IM-1-MOD(J,2)
!         DO I=2,IEND
         DO I=ISTA_M,IEND_M-MOD(J,2)
          IF(QV(I+IHE(J),J)<SPVAL.AND.UWND(I+IHE(J),J)<SPVAL.AND.&
             QV(I+IHW(J),J)<SPVAL.AND.UWND(I+IHW(J),J)<SPVAL.AND.&
             QV(I,J)<SPVAL.AND.QV(I,J-1)<SPVAL.AND.QV(I,J+1)<SPVAL) THEN
           R2DX   = 1./(2.*DX(I,J))
           R2DY   = 1./(2.*DY(I,J))
           QUDX   = (QV(I+IHE(J),J)*UWND(I+IHE(J),J)                   &
                    -QV(I+IHW(J),J)*UWND(I+IHW(J),J))*R2DX
           QVDY   = (QV(I,J+1)*VWND(I,J+1)-QV(I,J-1)*VWND(I,J-1))*R2DY

           QCNVG(I,J) = -(QUDX + QVDY) * HBM2(I,J)
          ELSE
           QCNVG(I,J) = SPVAL
          ENDIF
         ENDDO
       ENDDO
      ELSE IF(gridtype=='B')THEN
     
!       CALL EXCH(UWND)
!
!$omp  parallel do private(i,j,qudx,qvdy,r2dx,r2dy)
       DO J=JSTA_M,JEND_M
!        DO I=2,IM-1
        DO I=ISTA_M,IEND_M
         IF(UWND(I,J)<SPVAL.AND.UWND(I,J-1)<SPVAL.AND.&
            UWND(I-1,J)<SPVAL.AND.UWND(I-1,J-1)<SPVAL.AND.&
            Q1D(I,J)<SPVAL.AND.Q1D(I+1,J)<SPVAL.AND.Q1D(I-1,J)<SPVAL.AND.&
            VWND(I,J)<SPVAL.AND.VWND(I-1,J)<SPVAL.AND.&
            VWND(I,J-1)<SPVAL.AND.VWND(I-1,J-1)<SPVAL.AND.&
            Q1D(I,J+1)<SPVAL.AND.Q1D(I,J-1)<SPVAL) THEN
          R2DX   = 1./DX(I,J)
          R2DY   = 1./DY(I,J)
          QUDX=(0.5*(UWND(I,J)+UWND(I,J-1))*0.5*(Q1D(I,J)+Q1D(I+1,J))        &
               -0.5*(UWND(I-1,J)+UWND(I-1,J-1))*0.5*(Q1D(I,J)+Q1D(I-1,J)))*R2DX
          QVDY=(0.5*(VWND(I,J)+VWND(I-1,J))*0.5*(Q1D(I,J)+Q1D(I,J+1))        &
               -0.5*(VWND(I,J-1)+VWND(I-1,J-1))*0.5*(Q1D(I,J)+Q1D(I,J-1)))*R2DY
  
          QCNVG(I,J) = -(QUDX + QVDY)
         ELSE
          QCNVG(I,J) = SPVAL
         ENDIF
!	  print*,'mcvg=',i,j,r2dx,r2dy,QCNVG(I,J)
        ENDDO
       ENDDO
      ENDIF
!meb not sure about the indexing for the c-grid
!
!     END OF ROUTINE.
!     
      RETURN
      END

