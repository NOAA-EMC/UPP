!> @file
!> @brief bound() clips data in passed array.
!> 
!> @author Russ Treadon W/NP2 @date 1993-01-18
!> 
!> This routine bounds data in the passed array FLD (im x jm elements long) 
!> and clips data values such that on exiting the routine
!> @code
!>              FMIN <= FLD(I,J) <= FMAX
!> @endcode
!> for all points.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1993-01-18 | Russ Treadon | Initial
!> 1993-05-07 | Russ Treadon | Added DOCBLOC
!> 1998-05-29 | T Black      | Conversion from 1-D to 2-D
!> 2000-01-04 | Jim Tuccillo | MPI Version           
!> 2002-04-24 | Mike Baldwin | WRF Version
!> 2021-09002 | Bo Cui       | Decompose UPP in X direction
!>
!> @author Russ Treadon W/NP2 @date 1993-01-18
!---------------------------------------------------------------------------------------
!> @brief Clips data in passed array.
!> 
!> @param[in] FMIN Lower (inclusive) bound for data.
!> @param[in] FMAX Upper (inclusive) bound for data.
!> @param[out] FLD Array whose elements are bounded by [FMIN,FMAX].
!>
      SUBROUTINE BOUND(FLD,FMIN,FMAX)

!     
     use ctlblk_mod, only: jsta, jend, spval, im, jm, ista, iend
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     implicit none
!     
!     DECLARE VARIABLES.
      REAL,intent(in) :: FMAX, FMIN
      REAL,intent(inout) :: FLD(IM,JM)
      integer i,j
!     
!     
!**********************************************************************
!     START BOUND HERE.
!     
!     BOUND ARRAY.
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          if(fld(i,j) /= spval) then
            FLD(I,J) = min(FMAX, MAX(FMIN,FLD(I,J)))
          end if
        ENDDO
      ENDDO
!     
!     END OF ROUTINE.
!     
      RETURN
      END
