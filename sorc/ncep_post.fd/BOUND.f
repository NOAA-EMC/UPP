!> @file
!
!>
!!                .      .    .     
!! SUBPROGRAM:    BOUND       CLIPS DATA IN PASSED ARRAY
!!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-01-18
!!     
!! ABSTRACT:  THIS ROUTINE BOUNDS DATA IN THE PASSED ARRAY 
!!   FLD (IMxJM ELEMENTS LONG) AND CLIPS DATA VALUES SUCH 
!!   THAT ON EXITING THE ROUTINE 
!!                FMIN <= FLD(I,J) <= FMAX
!!   FOR ALL POINTS.
!!
!!     
!! PROGRAM HISTORY LOG:
!!   93-01-18  RUSS TREADON
!!   93-05-07  RUSS TREADON - ADDED DOCBLOC
!!   98-05-29  BLACK - CONVERSION FROM 1-D TO 2-D
!!   00-01-04  JIM TUCCILLO - MPI VERSION        
!!   02-04-24  MIKE BALDWIN - WRF VERSION        
!!   21-09-02  Bo Cui       - Decompose UPP in X direction
!!     
!! USAGE:    CALL bound(fld,fmin,fmax)
!!   INPUT ARGUMENT LIST:
!!     FMIN     - LOWER (INCLUSIVE) BOUND FOR DATA.
!!     FMAX     - UPPER (INCLUSIVE) BOUND FOR DATA.
!!
!!   OUTPUT ARGUMENT LIST: 
!!     FLD      - ARRAY WHOSE ELEMENTS ARE BOUNDED BY [FMIN,FMAX].
!!     
!!   OUTPUT FILES:
!!     NONE
!!     
!!   SUBPROGRAMS CALLED:
!!     UTILITIES:
!!       NONE
!!     LIBRARY:
!!       NONE
!!
!!   ATTRIBUTES:
!!     LANGUAGE: FORTRAN 90
!!     MACHINE : CRAY C-90
!!
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
