!> @file
!> @brief sclfld() scale array element by constant.
!>
!> @author Russ Treadon W/NP2 @date 1992-09-13

!> This routine multiples (scales) the first IMO*JMO
!> elements of array fld by the real scalar scale.
!> Array elements which equal a special value will
!> not be scaled by scale.  They will be left as is.
!> The special value, spval, is passed through common
!> block options.  It is set in include file options.
!>
!> @param[in] FLD Array whose elements are to be scaled.
!> @param[in] SCALE Constant by which to scale elements of fld.
!> @param[in] IMO,JMO Dimension of array fld.
!> @param[out] FLD Array whose elements have been scaled by scale.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-09-13 | Russ Treadon  | Initial
!> 2000-01-04 | Jim Tuccillo  | MPI Version
!>
!> @author Russ Treadon W/NP2 @date 1992-09-13
      SUBROUTINE SCLFLD(FLD,SCALE,IMO,JMO)
!

!
      use params_mod, only: small
      use ctlblk_mod, only: jsta, jend, spval
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!     
!     DECLARE VARIABLES.
!     
      integer,intent(in) :: IMO,JMO
      REAL,intent(in) ::  SCALE
      REAL,dimension(imo,jmo),intent(inout) :: FLD
      integer I,J
!     
!     
!***********************************************************************
!     START SCLFLD HERE
!     
!     MULTIPLY EACH ELEMENT OF FLD BY SCALE.
!     
!$omp  parallel do
      DO J=JSTA,JEND
      DO I=1,IMO
        IF(ABS(FLD(I,J)-SPVAL)>SMALL) FLD(I,J)=SCALE*FLD(I,J)
      ENDDO
      ENDDO
!     
!     END OF ROUTINE.
!     
      RETURN
      END
