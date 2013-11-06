      SUBROUTINE CALPOT(P1D,T1D,THETA)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALPOT      COMPUTES POTENTIAL TEMPERATURE
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-24
!     
! ABSTRACT: 
!     GIVEN PRESSURE AND TEMPERATURE THIS ROUTINE RETURNS
!     THE POTENTIAL TEMPERATURE.
!   .     
!     
! PROGRAM HISTORY LOG:
!   92-12-24  RUSS TREADON
!   98-06-15  T BLACK - CONVERSION FROM 1-D TO 2-D
!   00-01-04  JIM TUCCILLO - MPI VERSION            
!   02-04-24  MIKE BALDWIN - WRF VERSION         
!     
! USAGE:    CALL CALPOT(P1D,T1D,THETA)
!   INPUT ARGUMENT LIST:
!     P1D      - PRESSURE (PA)
!     T1D      - TEMPERATURE (K)
!
!   OUTPUT ARGUMENT LIST: 
!     THETA    - POTENTIAL TEMPERATURE (K)
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       NONE
!     LIBRARY:
!       NONE      
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN 90
!     MACHINE : CRAY C-90
!$$$  
!     
      use ctlblk_mod, only: jsta, jend, spval, im, jm
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      implicit none
!     
!     SET REQUIRED CONSTANTS.
      real,PARAMETER :: CAPA=0.28589641,P1000=1000.E2
!
!     DECLARE VARIABLES.
!     
      real,dimension(IM,JM),intent(in) :: P1D,T1D
      real,dimension(IM,JM),intent(inout) :: THETA

      REAL FAC
      integer I,J
!     
!**********************************************************************
!     START CALPOT HERE.
!     
!     COMPUTE THETA
!     
      DO J=JSTA,JEND
      DO I=1,IM
        IF(T1D(I,J).LT.SPVAL)THEN
         IF(ABS(P1D(I,J)).GT.1)THEN
          FAC=(P1000/P1D(I,J))**CAPA
          THETA(I,J)=FAC*T1D(I,J)
         ELSE
          THETA(I,J)=0.0
         ENDIF
        ELSE
         THETA(I,J)=SPVAL
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
