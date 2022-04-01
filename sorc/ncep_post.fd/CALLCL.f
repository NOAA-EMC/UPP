!> @file
!> @brief Subroutine that computes LCL heights and pressure.
!>
!> This routine computes the lifting condensation level
!> pressure and height in each column at mass points.
!> The height is above ground level. The equation used
!> to find the LCL pressure is from Boltan (1980, MWR)
!> and is the same as that used in subroutine CALCAPE.
!> 
!> This is a test version. Still to be resolved 
!> is the "best" parcel to lift.
!>
!> @param[in] P1D Array of parcel pressures (Pa).
!> @param[in] T1D Array of parcel temperatures (K).
!> @param[in] Q1D Array of parcel specific humidities (kg/kg).
!> @param[out] PLCL Parcel Pressure at LCL (Pa).
!> @param[out] ZLCL Parcel AGL height at LCL (m).
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1993-03-15 | Russ Treadon | Initial
!> 1998-06-16 | T Black      | Convesion from 1-D to 2-D
!> 2000-01-04 | Jim Tuccillo | MPI Version            
!> 2002-04-24 | Mike Baldwin | WRF Version            
!> 2019-10-30 | Bo Cui       | Remove "GOTO" Statement
!> 2021-07-28 | W Meng       | Restriction compuatation from undefined grids
!>
!> @author Russ Treadon W/NP2 @date 1993-03-15
      SUBROUTINE CALLCL(P1D,T1D,Q1D,PLCL,ZLCL)

!     
!     
      use vrbls3d, only: alpint, zint
      use vrbls2d, only: fis
      use masks, only: lmh
      use params_mod, only: eps, oneps, d01, h1m12, gi, d00
      use ctlblk_mod, only: jsta, jend, spval, jsta_m, jend_m, im
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
      real,PARAMETER :: D35=3.5, D4805=4.805,  H2840=2840.
      real,PARAMETER :: H55=55., D2845=0.2845, D28=0.28
!
!     DECLARE VARIABLES.
!     
      REAL,dimension(IM,jsta:jend), intent(in)    :: P1D,T1D,Q1D
      REAL,dimension(IM,jsta:jend), intent(inout) :: PLCL,ZLCL
      REAL TLCL(IM,jsta:jend)
      integer I,J,L,LLMH
      real DLPLCL,ZSFC,DZ,DALP,ALPLCL,RMX,EVP,ARG,RKAPA
!     
!**********************************************************************
!     START CALLCL HERE.
!     
!     LOAD OUTPUT ARRAYS WITH SPECIAL VALUE.
!     
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=1,IM
          PLCL(I,J) = SPVAL
          TLCL(I,J) = SPVAL
          ZLCL(I,J) = SPVAL
        ENDDO
      ENDDO
!     
!     COMPUTE PRESSURE, TEMPERATURE AND AGL HEIGHT AT LCL.
!
! Bo Cui 10/30/2019, remove "GOTO" statement

      DO 30 J=JSTA_M,JEND_M
      DO 30 I=2,IM-1
!     DO 30 I=1,IM
      IF(P1D(I,J)<spval.and.Q1D(I,J)<spval)THEN
      EVP       = P1D(I,J)*Q1D(I,J)/(EPS+ONEPS*Q1D(I,J))
      RMX       = EPS*EVP/(P1D(I,J)-EVP)
      RKAPA     = 1.0 / (D2845*(1.0-D28*RMX))
      ARG       = MAX(H1M12,EVP*D01)
      TLCL(I,J) = H55 + H2840 / (D35*LOG(T1D(I,J))-LOG(ARG)-D4805)
      PLCL(I,J) = P1D(I,J)*(TLCL(I,J)/T1D(I,J))**RKAPA
      ALPLCL    = LOG(PLCL(I,J))
      LLMH      = NINT(LMH(I,J))
      ZSFC      = FIS(I,J)*GI
!
      DO 20 L=LLMH,1,-1
      IF(ALPINT(I,J,L) < ALPLCL)THEN
        DLPLCL    = ALPLCL        - ALPINT(I,J,L+1)
        DALP      = ALPINT(I,J,L) - ALPINT(I,J,L+1)
        DZ        = ZINT(I,J,L)   - ZINT(I,J,L+1)
        ZLCL(I,J) = max(D00, ZINT(I,J,L+1) + DZ*DLPLCL/DALP - ZSFC)
        EXIT
      ENDIF
 20   CONTINUE
      ENDIF
 30   CONTINUE
!     
!     END OF ROUTINE.
!     
      RETURN
      END
