      SUBROUTINE CALDWP(P1D,Q1D,TDWP,T1D)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALDWP      COMPUTES 
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-22
!     
! ABSTRACT:  COMPUTES DEWPOINT FROM P, T, AND Q
!   .     
!     
! PROGRAM HISTORY LOG:
!   92-12-22  RUSS TREADON
!   93-10-04  RUSS TREADON - ADDED CHECK TO BOUND DEWPOINT
!                            TEMPERATURE TO NOT EXCEED THE
!                            AMBIENT TEMPERATURE.
!   98-06-08  T BLACK      - CONVERSION FROM 1-D TO 2-D
!   00-01-04  JIM TUCCILLO - MPI VERSION                
!     
! USAGE:    CALL CALDWP(P1D,Q1D,TDWP,T1D)
!   INPUT ARGUMENT LIST:
!     P1D      - PRESSURE (PA)
!     Q1D      - SPECIFIC HUMIDITY (KG/KG)
!     T1D      - TEMPERATURE (K)
!
!   OUTPUT ARGUMENT LIST: 
!     TDWP     - DEWPOINT TEMPERATURE (K)

!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       DEWPOINT - COMPUTES DEWPOINT GIVEN VAPOR PRESSURE.
!     LIBRARY:
!       NONE
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN 90
!     MACHINE : CRAY C-90
!$$$  
!
!
!     SET PARAMETERS.
     use params_mod, only: eps, oneps, d001, h1m12
     use ctlblk_mod, only: jsta, jend, im, jm
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     DECLARE VARIABLES.
!     
     REAL,dimension(IM,JM),intent(in) ::  P1D,Q1D,T1D
     REAL,dimension(IM,JM),intent(inout) ::  TDWP

     REAL EVP(IM,JM)
     integer I,J
!     
!****************************************************************************
!     START CALDWP HERE.
!     
!     COMPUTE VAPOR PRESSURE.  CONVERT TO CENITBARS.
!     
      DO J=JSTA,JEND
      DO I=1,IM
        EVP(I,J) = P1D(I,J)*Q1D(I,J)/(EPS+ONEPS*Q1D(I,J))
        EVP(I,J) = EVP(I,J)*D001
        EVP(I,J) = AMAX1(H1M12,EVP(I,J))
      ENDDO
      ENDDO
!     
!     COMPUTE DEWPOINT TEMPERATURE.
!     
      CALL DEWPOINT(EVP,TDWP)
!     
!     ENSURE DEWPOINT TEMPERATURE DOES NOT EXCEED AMBIENT TEMPERATURE.
!     
      DO J=JSTA,JEND
      DO I=1,IM
        IF (TDWP(I,J).GT.T1D(I,J)) TDWP(I,J) = T1D(I,J)
      ENDDO
      ENDDO
!
!     END OF ROUTINE.
!     
      RETURN
      END
