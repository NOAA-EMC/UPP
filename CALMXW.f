      SUBROUTINE CALMXW(MXWP,MXWZ,MXWU,MXWV,MXWT)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALMXW      COMPUTE MAX WIND LEVEL 
!   PRGRMMR: MANIKIN        ORG: W/NP2   DATE: 97-03-04       
!     
! ABSTRACT:  
!     THIS ROUTINE COMPUTES MAX WIND LEVEL.  AT EACH POINT,
!   IT FINDS THE MAX WIND ABOVE 500 MB AND DETERMINES THE
!   PRESSURE AND HEIGHT AT THAT LEVEL.
!     
!     
! PROGRAM HISTORY LOG:
!   97-03-04 GEOFF MANIKIN
!   98-06-15 T BLACK - CONVERSION FROM 1-D TO 2-D
!   00-01-02 JIM TUCCILLO - MPI VERSION
!   02-01-15  MIKE BALDWIN - WRF VERSION
!   05-02-24 H CHUANG - ADD WRF NMM COMPONENTS 
!   05-07-07 BINBIN ZHOU - ADD RSM 
!   
! USAGE:    CALL  CALMXW(MXWP,MXWZ,MXWU,MXWV)
!   INPUT ARGUMENT LIST:
!     NONE     
!
!   OUTPUT ARGUMENT LIST: 
!     MXWP    - PRESSURE LEVEL OF THE MAX WIND
!     MXWZ    - HEIGHT OF THE MAX WIND
!     MXWU    - U COMPONENT OF THE ACTUAL MAX WIND 
!     MXWV    - V COMPONENT OF THE ACTUAL MAX WIND
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
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
      use masks
      use params_mod
      use ctlblk_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
! 
!     DECLARE VARIABLES.
!     
      real,dimension(IM,JM),intent(inout) :: MXWP,MXWZ,MXWU,MXWV,MXWT
!
      REAL MXWW 
      INTEGER IHE(JM),IHW(JM) 
      integer I,J,LLMH,IE,IW,L
      real CRITP,U0,V0,WIND
!     
!     
!*****************************************************************************
!     START CALMXW HERE.
!     
!     LOOP OVER THE GRID.
!    
      CRITP=5.0E4
!
      DO J=JSTA,JEND
      DO I=1,IM
        MXWU(I,J) = SPVAL
        MXWV(I,J) = SPVAL
        MXWP(I,J) = SPVAL
        MXWZ(I,J) = SPVAL 
	MXWT(I,J) = SPVAL
      ENDDO
      ENDDO
!
!$omp  parallel do
!$omp& private(ie,iw,mxww,u0,v0,wind)
      IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN 
       DO 20 J=JSTA,JEND
       DO 20 I=1,IM
        MXWW  = -1000.
        LLMH=NINT(LMH(I,J))
!
        DO 10 L= LLMH-1,1,-1
         U0 = UH(I,J,L)
         V0 = VH(I,J,L)
         WIND = SQRT(U0**2 + V0**2)

!  MAX WIND LEVEL MUST BE ABOVE THE 500 MB 

         IF (WIND .GT. MXWW .and. PMID(I,J,L) .LT. CRITP) THEN
           MXWU(I,J) = U0
           MXWV(I,J) = V0
           MXWW = WIND 
           MXWP(I,J) = PMID(I,J,L)
           MXWZ(I,J) = ZMID(I,J,L)
	   MXWT(I,J) = T(I,J,L)
         ENDIF
   10   CONTINUE
   20  CONTINUE
      
      ELSE IF(MODELNAME .EQ. 'GFS')THEN
       DO 22 J=JSTA,JEND
       DO 22 I=1,IM
        MXWW  = -1000.
!
        DO 24 L= LM-1,1,-1
         U0 = UH(I,J,L)
         V0 = VH(I,J,L)
         WIND = SQRT(U0**2 + V0**2)

!  MAX WIND LEVEL MUST BE ABOVE THE 500 MB AND BELOW 100MB FOR GFS

         IF (WIND .GT. MXWW .and. PMID(I,J,L) .LT. CRITP          &
            .and. PMID(I,J,L) .GT. 10000.) THEN
           MXWU(I,J) = U0
           MXWV(I,J) = V0
           MXWW = WIND 
           MXWP(I,J) = PMID(I,J,L)
           MXWZ(I,J) = ZMID(I,J,L)
	   MXWT(I,J) = T(I,J,L)
         ENDIF
   24   CONTINUE
   22  CONTINUE

      ELSE IF(MODELNAME .EQ. 'NMM')THEN
       
       DO J=JSTA_M,JEND_M
        IHE(J)=MOD(J+1,2)
        IHW(J)=IHE(J)-1
       ENDDO
       
       DO 40 J=JSTA_M,JEND_M
       DO 40 I=2,IM-1
        IE=I+IHE(J)
        IW=I+IHW(J)
        MXWW  = -1000.
        LLMH=NINT(LMH(I,J))
!
        DO 30 L= LLMH-1,1,-1
         U0 = D25*(U(I,J-1,L)+U(IW,J,L)+U(IE,J,L)+U(I,J+1,L))
         V0 = D25*(V(I,J-1,L)+V(IW,J,L)+V(IE,J,L)+V(I,J+1,L))
         WIND = SQRT(U0**2 + V0**2)
	 
	 IF (WIND .GT. MXWW .and. PMID(I,J,L) .LT. CRITP) THEN
           MXWU(I,J) = U0
           MXWV(I,J) = V0
           MXWW = WIND
           MXWP(I,J) = PMID(I,J,L) 
           MXWZ(I,J)=ZMID(I,J,L)
	   MXWT(I,J) = T(I,J,L)
         END IF
   30   CONTINUE
   40  CONTINUE

      END IF
!     END OF ROUTINE.
!     
      RETURN
      END
