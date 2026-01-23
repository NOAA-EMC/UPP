!> @file
!> @brief SUBROUTINE THAT COMPUTES WIND GUST FROM TURBULENCE AND CONVECTIVE COMPONENTS
!>        FOR HAFS MODEL.
!> 
!> THIS SUBROUTINE IS BASED ON THE PYTHON CODE PROVIDED BY ANDREW HAZELTON (HAFS TEAM),
!> WHICH COMPUTES A NEW WIND GUST FACTOR BASED ON TURBULENCE AND CONVECTIVE COMPONTNETS.
!> 
!> @param[in] SPEED850 WIND SPEED AT 850 MB
!> @param[in] SPEED950 WIND SPEED AT 950 MB
!> @param[inout] GUST FACTOR FROM CONVECTIVE AND TURBULENT COMPONENTS
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2026-01-06 | Karina Asmar | Initial
!>   
!> @author Karina Asmar NCEP/EMC @date 2026-01-06

      SUBROUTINE CALGUSTCONV(SPEED850,SPEED950,GUSTCONV)
!     
!     
      use vrbls2d , only: u10,v10, ustar
      use ctlblk_mod, only: ista, iend, jsta, jend, ista_2l, iend_2u, jsta_2l, jend_2u, spval

      implicit none

      INCLUDE "mpif.h"
! 
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE PARAMETERS.
!
!     DECLARE VARIABLES.
!     
      REAL,intent(in)    :: SPEED850(ista_2l:iend_2u,jsta_2l:jend_2u)
      REAL,intent(in)    :: SPEED950(ista_2l:iend_2u,jsta_2l:jend_2u)
      REAL,intent(inout) :: GUSTCONV(ista_2l:iend_2u,jsta_2l:jend_2u)

      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u) :: WIND10, WSD, GUST1_NEW, GUST2_NEW, WSTT1_NEW, &
                                                            WSTT2_NEW, GF1_NEW, GF2_NEW

!
      integer I,J
!     
!     
!*****************************************************************************
!> CALGUSTCONV computes new gust wind factor for HAFS based on convective and
!> turbulent components.
!     START CALGUSTCONV HERE.
!     
!    
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          GUSTCONV(I,J) = SPVAL 
        ENDDO
      ENDDO

! CALCULATE GUST FACTORS
!$omp parallel do private(i,j)
     DO J=JSTA,JEND
       DO I=ISTA,IEND
       IF(U10(I,J)<SPVAL.AND.V10(I,J)<SPVAL.AND.SPEED850(I,J)<SPVAL.AND.SPEED950(I,J)<SPVAL) THEN
       ! CALCULATE WINDS AND WIND SPEED DIFFERENCE
       WIND10(I,J) = SQRT(U10(I,J)**2 + V10(I,J)**2)
       WSD(I,J) = SPEED850(I,J) - SPEED950(I,J)
       IF (WSD(I,J) < 0.0) WSD(I,J) = 0.0      ! MAKE SURE WIND SPEED IS >= ZERO

       ! CALCULATE NEW GUST FACTOR
       GUST1_NEW(I,J) = 3.0 * USTAR(I,J)  ! TURBULENT COMPONENT
       WSTT1_NEW(I,J) = WIND10(I,J) + GUST1_NEW(I,J)
       GF1_NEW(I,J) = WSTT1_NEW(I,J) / WIND10(I,J)
       GUST2_NEW(I,J) = 0.3 * WSD(I,J)  ! CONVECTIVE MIXING COMPONENT

       GUSTCONV(I,J) = (WSTT1_NEW(I,J) + GUST2_NEW(I,J))
       ENDIF
       ENDDO
     ENDDO
!
!     END OF ROUTINE.
!     
      RETURN
      END
