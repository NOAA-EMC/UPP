!> @file
!> @brief trpaus2() computes tropopause level fields.
!> 
!> This routine computes tropopause data.
!> Code is adapted from original routine trpaus().
!> 
!> At each mass point a downward search is made for the 
!> first occurrence of a local stability parameter
!> (based on the square of the Brunt-Vaisalla frequency)
!> exceeding the value of Rd / Cp - as used by ECMWF.
!> A maximum tropopause pressure of ~500mb is enforced.
!> Once the tropopause is located in a column, pressure,
!> temperature, u and v winds, and vertical wind shear
!> are computed.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-22 | Russ Treadon  | Initial
!> 1997-03-06 | Geoff Manikin | Changed criteria for determining the tropopause and added height
!> 1998-06-15 | T Black       | Conversion from 1-D TO 2-D
!> 2000-01-04 | Jim Tuccillo  | MPI Version
!> 2002-04-23 | Mike Baldwin  | WRF Version
!> 2019-10-30 | Bo Cui        | ReMOVE "GOTO" STATEMENT
!> 2021-09-13 | JESSE MENG    | 2D DECOMPOSITION
!> 2026-05-20 | C Hill        | Alternative, ECMWF-employed algorithm 
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
!------------------------------------------------------------------------------
!> @brief Computes tropopause data.
!> 
!> @param[out] PTROP Tropopause pressure.
!> @param[out] TTROP Tropopause temperature.
!> @param[out] ZTROP Tropopause height.
!> @param[out] UTROP Tropopause u wind component.
!> @param[out] VTROP Tropopause v wind component.
!> @param[out] SHTROP Vertical wind shear at tropopause.
!>
      SUBROUTINE TRPAUS2(PTROP,TTROP,ZTROP,UTROP,VTROP,SHTROP)

!     
!     
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE PARAMETERS.
!
       use vrbls3d,    only: pint, t, zint, pmid, zmid, uh, vh
       use masks,      only: lmh
       use physcons_post, only: CON_G, CON_RD, CON_ROCP
       use ctlblk_mod, only: jsta, jend, spval, im, jm, lm, ista, iend
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
!     PARAMETER TRPPRM SPECIFIES A LOCAL STABILITY VALUE
!     (IN ---) FOR IDENTIFYING THE TROPOPAUSE.  WE START 
!     LOOKING FOR THE TROPOPAUSE BEGINNING AT PRESSURE LEVEL
!     'PSTART' (IN PASCALS) AND DOWNWARD TO 'PEND'.
      real,PARAMETER :: PSTART=7.0E3, PEND=5.0E4
      real,PARAMETER :: N02=2.5E-4
!     
!     DECLARE VARIABLES.
!     
      REAL,dimension(IM,JM) :: PTROP,TTROP,ZTROP,UTROP,VTROP,SHTROP
      REAL TRPPRM(LM)
!
      integer I,J,LLMH,L
      real PM,DELT,DZ,DP,RSQDIF
!     
!*****************************************************************************
!     START TRPAUS HERE.
!     
!     LOOP OVER THE HORIZONTAL GRID.
!    
      DO J=JSTA,JEND
      DO I=ISTA,IEND
         PTROP(I,J)  = SPVAL
         TTROP(I,J)  = SPVAL
         ZTROP(I,J)  = SPVAL
         UTROP(I,J)  = SPVAL
         VTROP(I,J)  = SPVAL
         SHTROP(I,J) = SPVAL
      ENDDO
      ENDDO
!
!!!$omp parallel do private(i,j,delt,dz,dp,l,llmh,pm,pmd,rsqdif,trpprm)

       DO J=JSTA,JEND
        DO I=ISTA,IEND
!     
!        COMPUTE A STABILITY PARAMETER AT EACH ETA LAYER DESCENDING
!        FROM THE TOP LEVEL. THE FIRST ETA LAYER BELOW PRESSURE
!        LEVEL "PSTART" WHERE THE {PARAMETER > Rd / Cp} IS
!        LABELED THE TROPOPAUSE.
!
        LLMH=NINT(LMH(I,J))
!
        TRPLVL = .FALSE.
        loopL: DO L = 2,LLMH-1
        PM     = PINT(I,J,L)
        DELT   = T(I,J,L-1)-T(I,J,L)
        DP     = PMID(I,J,L-1)-PMID(I,J,L)
        TRPPRM(L) = ((PMID(I,J,L)/T(I,J,L))*(DELT/DP)) &
     &              +((CON_RD*T(I,J,L)/(CON_G**2.))*N02)
!
        IF ((PM>PSTART).AND.(PM<PEND)) THEN
         IF (TRPPRM(L)>CON_ROCP) THEN
          TRPLVL = .TRUE.
!
          PTROP(I,J)  = PMID(I,J,L)
          TTROP(I,J)  = T(I,J,L)
          ZTROP(I,J)  = ZMID(I,J,L)
!
          UTROP (I,J) = UH(I,J,L)
          VTROP (I,J) = VH(I,J,L)
          DZ        = ZINT(I,J,L)-ZINT(I,J,L+1)
          RSQDIF    = SQRT(((UH(I,J,L-1)-UH(I,J,L+1))*0.5)**2. &
     &                    +((VH(I,J,L-1)-VH(I,J,L+1))*0.5)**2.)
          SHTROP(I,J) = RSQDIF/DZ
         ENDIF
        ENDIF

        IF (TRPLVL == .TRUE.) EXIT loopL
        ENDDO loopL

       ENDDO !end I
      ENDDO !end J

!     
!     END OF ROUTINE.
!     
      RETURN
      END
