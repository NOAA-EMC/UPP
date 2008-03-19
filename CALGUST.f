      SUBROUTINE CALGUST(LPBL,ZPBL,GUST)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALGUST      COMPUTE MAX WIND LEVEL 
!   PRGRMMR: MANIKIN        ORG: W/NP2   DATE: 97-03-04       
!     
! ABSTRACT:  
!     THIS ROUTINE COMPUTES SURFACE WIND GUST BY MIXING
!  DOWN MOMENTUM FROM THE LEVEL AT THE HEIGHT OF THE PBL
!     
!     
! PROGRAM HISTORY LOG:
!   03-10-15 GEOFF MANIKIN
!   05-03-09 H CHUANG - WRF VERSION
!   05-07-07 BINBIN ZHOU - ADD RSM   
!   
! USAGE:    CALL CALGUST(GUST) 
!   INPUT ARGUMENT LIST:
!     NONE     
!
!   OUTPUT ARGUMENT LIST: 
!     GUST    - SPEED OF THE MAXIMUM SFC WIND GUST 
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       H2V     
!
!     LIBRARY:
!       COMMON   - 
!                  LOOPS
!                  OPTIONS
!                  MASKS
!                  INDX
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN 90
!     MACHINE : CRAY C-90
!$$$  
!     
!     
      use vrbls3d
      use vrbls2d 
      use params_mod
      use ctlblk_mod
! 
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE PARAMETERS.
!
!     DECLARE VARIABLES.
!     
      INTEGER,intent(in) ::  LPBL(IM,JM)
      REAL,intent(in) ::  ZPBL(IM,jsta_2l:jend_2u)
      REAL,intent(inout) :: GUST(IM,JM)

      integer I,J,IE,IW
      real ZSFC,DELWIND,USFC,VSFC,SFCWIND,WIND,U0,V0
!     
!     
!*****************************************************************************
!     START CALMXW HERE.
!     
!     LOOP OVER THE GRID.
!    
      DO J=JSTA,JEND
      DO I=1,IM
        GUST(I,J) = SPVAL 
      ENDDO
      ENDDO
!
!     ASSUME THAT U AND V HAVE UPDATED HALOS
!
!$omp  parallel do
!$omp& private(ie,iw,mxww,u0,v0,wind)
      DO 20 J=JSTA_M,JEND_M
      DO 20 I=2,IM-1
       L=LPBL(I,J) 
       IF(MODELNAME .EQ. 'NMM')THEN
        IE=I+MOD(J+1,2) 
        IW=I+MOD(J+1,2)-1
	
        USFC=D25*(U10(I,J-1)+U10(IW,J)+U10(IE,J)+U10(I,J+1)) 
        VSFC=D25*(V10(I,J-1)+V10(IW,J)+V10(IE,J)+V10(I,J+1))
        SFCWIND=SQRT(USFC**2 + VSFC**2)
        U0 = D25*(UH(I,J-1,L)+UH(IW,J,L)+UH(IE,J,L)+UH(I,J+1,L))
        V0 = D25*(VH(I,J-1,L)+VH(IW,J,L)+VH(IE,J,L)+VH(I,J+1,L))
        WIND=SQRT(U0**2 + V0**2)
        
       ELSE
        USFC=U10(I,J)
        VSFC=V10(I,J)
        SFCWIND=SQRT(USFC**2 + VSFC**2) 
        U0=UH(I,J,L)
        V0=VH(I,J,L)
        WIND=SQRT(U0**2 + V0**2)
       END IF
       DELWIND=WIND - SFCWIND
       ZSFC=FIS(I,J)*GI
       DELWIND=DELWIND*(1.0-AMIN1(0.5,ZPBL(I,J)/2000.))
       GUST(I,J)=SFCWIND+DELWIND
   10 CONTINUE
   20 CONTINUE

!     END OF ROUTINE.
!     
      RETURN
      END
