!> @file
!> @brief Subroutine that computes potential temperature.
!>
!> Given pressure and temperature this routine returns
!> the potential temperature.
!>
!> @param[in] P1D pressures (Pa).
!> @param[in] T1D temperatures (K).
!> @param[out] THETA potential temperatures (K).
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-24 | Russ Treadon | Initial
!> 1998-06-15 | T Black      | Convesion from 1-D to 2-D
!> 2000-01-04 | Jim Tuccillo | MPI Version            
!> 2002-04-24 | Mike Baldwin | WRF Version      
!> 2021-09-02 | Bo Cui       | Decompose UPP in X direction          
!>
!> @author Russ Treadon W/NP2 @date 1992-12-24
!-----------------------------------------------------------------------
!> @brief Subroutine that computes potential temperature.
!>
!> @param[in] P1D pressures (Pa).
!> @param[in] T1D temperatures (K).
!> @param[out] THETA potential temperatures (K).
!-----------------------------------------------------------------------
      SUBROUTINE CALPOT(P1D,T1D,THETA)

!     
      use ctlblk_mod, only: jsta, jend, spval, im,  ista, iend
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      implicit none
!     
!     SET REQUIRED CONSTANTS.
      real,PARAMETER :: CAPA=0.28589641,P1000=1000.E2
!
!     DECLARE VARIABLES.
!     
      real,dimension(ista:iend,jsta:jend),intent(in)    :: P1D,T1D
      real,dimension(ista:iend,jsta:jend),intent(inout) :: THETA

      integer I,J
!     
!**********************************************************************
!     START CALPOT HERE.
!     
!     COMPUTE THETA
!     
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IF(T1D(I,J) < SPVAL) THEN
!           IF(ABS(P1D(I,J)) > 1.0) THEN
            IF(P1D(I,J) > 1.0) THEN
              THETA(I,J) = T1D(I,J) * (P1000/P1D(I,J))**CAPA
            ELSE
              THETA(I,J) = 0.0
            ENDIF
          ELSE
            THETA(I,J) = SPVAL
          ENDIF
        ENDDO
      ENDDO
!     do j = 180, 185
!        print *, ' me, j, p1d,t1d,theta = ',
!    *   me, j, p1d(10,j),t1d(10,j),theta (10,j)
!     end do
!       stop
!     
!     END OF ROUTINE.
!     
      RETURN
      END
