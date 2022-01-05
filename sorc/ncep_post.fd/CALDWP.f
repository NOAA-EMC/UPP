!> @file
!> @brief Subroutine related to dewpoint temperature.
!
!> Computes dewpoint from P, T, and Q
!>     
!> @param[in] P1D Pressure (Pa)
!> @param[in] Q1D Specific humidity (kg/kg)
!> @param[in] T1D Temperature (K)
!> @param[out] TDWP Dewpoint temperature (K)
!>
!> Program history
!> - 92-12-22  Russ Treadon
!> - 93-10-04  Russ Treadon - Added check to bound dewpoint
!>                            temperature to not exceed the
!>                            ambient temperature.
!> - 98-06-08  T BLACK      - Conversion from 1-D to 2-D
!> - 00-01-04  Jim Tuccillo - MPI version                
!> - 21-07-23  Wen Meng     - Retrict computation from undefined points
!>     
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE CALDWP(P1D,Q1D,TDWP,T1D)

!
!
!     SET PARAMETERS.
     use params_mod, only: eps, oneps, d001, h1m12
     use ctlblk_mod, only: jsta, jend, im, spval, ista, iend
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     DECLARE VARIABLES.
!     
     REAL,dimension(ista:iend,jsta:jend),intent(in)    ::  P1D,Q1D,T1D
     REAL,dimension(ista:iend,jsta:jend),intent(inout) ::  TDWP

     REAL EVP(ista:iend,jsta:jend)
     integer I,J
!     
!****************************************************************************
!     START CALDWP HERE.
!     
!     COMPUTE VAPOR PRESSURE.  CONVERT TO CENITBARS.
!     
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IF(P1D(I,j)<spval .and. Q1D(I,J)<spval) THEN
          EVP(I,J) = P1D(I,J)*Q1D(I,J)/(EPS+ONEPS*Q1D(I,J))
          EVP(I,J) = MAX(H1M12,EVP(I,J)*D001)
          ELSE
          EVP(I,J) = spval
          ENDIF
        ENDDO
      ENDDO
!     
!     COMPUTE DEWPOINT TEMPERATURE.
!     
      CALL DEWPOINT(EVP,TDWP)
!     
!     ENSURE DEWPOINT TEMPERATURE DOES NOT EXCEED AMBIENT TEMPERATURE.
!     
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          TDWP(I,J) = min(TDWP(I,J),T1D(I,J))
        ENDDO
      ENDDO
!
!     END OF ROUTINE.
!     
      RETURN
      END
