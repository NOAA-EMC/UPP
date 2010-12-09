      SUBROUTINE CALHEL(DEPTH,UST,VST,HELI)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALHEL       COMPUTES STORM RELATIVE HELICITY
!   PRGRMMR: BALDWIN         ORG: W/NP2      DATE: 94-08-22       
!     
! ABSTRACT:
!     THIS ROUTINE COMPUTES ESTIMATED STORM MOTION AND
!     STORM-RELATIVE ENVIRONMENTAL HELICITY.  
!     (DAVIES-JONES ET AL 1990) THE ALGORITHM PROCEEDS AS 
!     FOLLOWS.
!     
!     THE STORM MOTION COMPUTATION NO LONGER EMPLOYS THE DAVIES AND
!     JOHNS (1993) METHOD WHICH DEFINED STORM MOTION AS 30 DEGREES TO
!     THE RIGHT OF THE 0-6 KM MEAN WIND AT 75% OF THE SPEED FOR MEAN
!     SPEEDS LESS THAN 15 M/S AND 20 DEGREES TO THE RIGHT FOR SPEEDS
!     GREATER THAN 15 M/S.   INSTEAD, WE NOW USE THE DYNAMIC METHOD
!     (BUNKERS ET AL. 1998) WHICH HAS BEEN FOUND TO DO BETTER IN
!     CASES WITH 'NON-CLASSIC' HODOGRAPHS (SUCH AS NORTHWEST-FLOW
!     EVENTS) AND DO AS WELL OR BETTER THAN THE OLD METHOD IN MORE
!     CLASSIC SITUATIONS. 
!     
! PROGRAM HISTORY LOG:
!   94-08-22  MICHAEL BALDWIN
!   97-03-27  MICHAEL BALDWIN - SPEED UP CODE
!   98-06-15  T BLACK         - CONVERSION FROM 1-D TO 2-D
!   00-01-04  JIM TUCCILLO    - MPI VERSION
!   00-01-10  G MANIKIN       - CHANGED TO BUNKERS METHOD
!   02-05-22  G MANIKIN       - NOW ALLOW CHOICE OF COMPUTING
!                               HELICITY OVER TWO DIFFERENT
!                               (0-1 and 0-3 KM) DEPTHS
!   03-03-25  G MANIKIN       - MODIFIED CODE TO COMPUTE MEAN WINDS
!                               USING ARITHMETIC AVERAGES INSTEAD OF
!                               MASS WEIGHTING;  DIFFERENCES ARE MINOR
!                               BUT WANT TO BE CONSISTENT WITH THE
!                               BUNKERS METHOD
!   04-04-16  M PYLE          - MINIMAL MODIFICATIONS, BUT PUT INTO
!                                NMM WRFPOST CODE
!   05=02-25  H CHUANG        - ADD COMPUTATION FOR ARW A GRID
!   05-07-07  BINBIN ZHOU     - ADD RSM FOR A GRID  
!   
! USAGE:    CALHEL(UST,VST,HELI)
!   INPUT ARGUMENT LIST:
!     DPTH      - DEPTH IN METERS OVER WHICH HELICITY SHOULD BE COMPUTED;
!                 ALLOWS ONE TO DISTINGUISH 0-3 KM AND 0-1 KM VALUES
!
!   OUTPUT ARGUMENT LIST: 
!     UST      - ESTIMATED U COMPONENT (M/S) OF STORM MOTION.
!     VST      - ESTIMATED V COMPONENT (M/S) OF STORM MOTION.
!     HELI     - STORM-RELATIVE HELICITY (M**2/S**2)
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!
!     LIBRARY:
!       COMMON   - VRBLS
!                  LOOPS
!                  PHYS 
!                  EXTRA
!                  MASKS
!                  OPTIONS
!                  INDX
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN 90
!     MACHINE : IBM SP
!$$$  
!
      use vrbls3d
      use vrbls2d
      use masks
      use params_mod
      use lookup_mod,only :ITB,JTB,ITBQ,JTBQ
      use ctlblk_mod
      use gridspec_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
      real,PARAMETER :: P150=15000.0,P300=30000.0,S15=15.0
      real,PARAMETER :: D3000=3000.0,PI6=0.5235987756,PI9=0.34906585
      real,PARAMETER :: D5500=5500.0,D6000=6000.0,D7000=7000.0
      real,PARAMETER :: D500=500.0
!     
!     DECLARE VARIABLES
!     
      real,intent(in) ::  DEPTH
      REAL,dimension(IM,JM),intent(inout) ::  UST,VST,HELI
!
      REAL HTSFC(IM,JM)
!
      REAL UST6(IM,JM),VST6(IM,JM)
      REAL UST5(IM,JM),VST5(IM,JM)
      REAL UST1(IM,JM),VST1(IM,JM)
      INTEGER COUNT6(IM,JM),COUNT5(IM,JM),COUNT1(IM,JM)
	
      INTEGER IVE(JM),IVW(JM)
      integer I,J,IW,IE,JS,JN,JVN,JVS,L
      integer ISTART,ISTOP,JSTART,JSTOP
      real Z2,DZABV,UMEAN5,VMEAN5,UMEAN1,VMEAN1,UMEAN6,VMEAN6,      &
           USHR,VSHR,Z1,Z3,DZ,DZ1,DZ2,DU1,DU2,DV1,DV2
!     
!****************************************************************
!     START CALHEL HERE
!     
!     INITIALIZE ARRAYS.
!     
!$omp  parallel do
      DO J=JSTA,JEND
      DO I=1,IM
         UST(I,J)    = 0.0
         VST(I,J)    = 0.0
         HELI(I,J)   = 0.0
         UST1(I,J)   = 0.0
         VST1(I,J)   = 0.0
         UST5(I,J)   = 0.0
         VST5(I,J)   = 0.0
         UST6(I,J)   = 0.0
         VST6(I,J)   = 0.0
         COUNT6(I,J) = 0
         COUNT5(I,J) = 0
         COUNT1(I,J) = 0
      ENDDO
      ENDDO
      IF(gridtype=='E')THEN
        JVN=1
        JVS=-1
	do J=JSTA,JEND
	IVE(J)=MOD(J,2)
	IVW(J)=IVE(J)-1
	enddo
	ISTART=2
        ISTOP=IM-1
        JSTART=JSTA_M
        JSTOP=JEND_M
      ELSE IF(gridtype=='B')THEN
        JVN=1
        JVS=0
	do J=JSTA,JEND
	IVE(J)=1
	IVW(J)=0
	enddo
	ISTART=2
        ISTOP=IM-1
        JSTART=JSTA_M
        JSTOP=JEND_M
      ELSE
        JVN=0
        JVS=0
        do J=JSTA,JEND
        IVE(J)=0
        IVW(J)=0
        enddo
	ISTART=1
        ISTOP=IM
        JSTART=JSTA
        JSTOP=JEND 
      END IF 
!
!     LOOP OVER HORIZONTAL GRID.
!
!      CALL EXCH(RES(1,jsta_2l)
!      CALL EXCH(PD()

!      DO L = 1,LP1
!        CALL EXCH(ZINT(1,jsta_2l,L))
!      END DO
! 
!$omp  parallel do
!$omp& private(htsfc,ie,iw,pdslvk,pkl,psfck)
      IF(gridtype/='A')CALL EXCH(FIS(1:IM,JSTA_2L:JEND_2U))
      DO L = 1,LM
        IF(gridtype/='A')CALL EXCH(ZMID(1:IM,JSTA_2L:JEND_2U,L)) 
        DO J=JSTART,JSTOP
        DO I=ISTART,ISTOP
          IE=I+IVE(J)
          IW=I+IVW(J)
          JN=J+JVN 
          JS=J+JVS
!mp          PDSLVK=(PD(IW,J)*RES(IW,J)+PD(IE,J)*RES(IE,J)+
!mp     1           PD(I,J+1)*RES(I,J+1)+PD(I,J-1)*RES(I,J-1))*0.25
!mp          PSFCK=AETA(LMV(I,J))*PDSLVK+PT
          IF (gridtype=='B')THEN
	   HTSFC(I,J)=(FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(IE,JN))/4.0/G
!     
!     COMPUTE MASS WEIGHTED MEAN WIND IN THE 0-6 KM LAYER, THE
!  0-0.5 KM LAYER, AND THE 5.5-6 KM LAYER 
!
           Z2=0.25*(ZMID(IW,J,L)+ZMID(IE,J,L)+                       &
                  ZMID(I,JN,L)+ZMID(IE,JN,L))                       
	  ELSE
	   HTSFC(I,J)=(FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(I,JS))/4.0/G	  
!     
!     COMPUTE MASS WEIGHTED MEAN WIND IN THE 0-6 KM LAYER, THE
!  0-0.5 KM LAYER, AND THE 5.5-6 KM LAYER 
!
           Z2=0.25*(ZMID(IW,J,L)+ZMID(IE,J,L)+                       &
                  ZMID(I,JN,L)+ZMID(I,JS,L))
          END IF
          DZABV=Z2-HTSFC(I,J)
  
          IF (DZABV.LE.D6000 .AND. L.LE.NINT(LMV(I,J))) THEN
               UST6(I,J) = UST6(I,J) + UH(I,J,L) 
               VST6(I,J) = VST6(I,J) + VH(I,J,L)
               COUNT6(I,J) = COUNT6(I,J) + 1 
          ENDIF

          IF (DZABV.LT.D6000 .AND. DZABV.GE.D5500 .AND.              &
             L.LE.NINT(LMV(I,J))) THEN
               UST5(I,J) = UST5(I,J) + UH(I,J,L)
               VST5(I,J) = VST5(I,J) + VH(I,J,L)
               COUNT5(I,J) = COUNT5(I,J) + 1
          ENDIF         

          IF (DZABV.LT.D500 .AND. L.LE.NINT(LMV(I,J))) THEN
               UST1(I,J) = UST1(I,J) + UH(I,J,L)
               VST1(I,J) = VST1(I,J) + VH(I,J,L) 
               COUNT1(I,J) = COUNT1(I,J) + 1
          ENDIF

        ENDDO
        ENDDO
      ENDDO
!
! CASE WHERE THERE IS NO LEVEL WITH HEIGHT BETWEEN 5500 AND 6000
!
      DO J=JSTART,JSTOP
      DO I=ISTART,ISTOP
        IF (COUNT5(I,J) .EQ. 0) THEN
         DO L=LM,1,-1
          IE=I+IVE(J)
          IW=I+IVW(J)
          JN=J+JVN
          JS=J+JVS
	  IF (gridtype=='B')THEN
	   Z2=0.25*(ZMID(IW,J,L)+ZMID(IE,J,L)+                       &
                  ZMID(I,JN,L)+ZMID(IE,JN,L))                       
	  ELSE
	   Z2=0.25*(ZMID(IW,J,L)+ZMID(IE,J,L)+                       &
                  ZMID(I,JN,L)+ZMID(I,JS,L))
	  END IF	  

          DZABV=Z2-HTSFC(I,J)
          IF (DZABV.LT.D7000 .AND. DZABV.GE.D6000) THEN 
               UST5(I,J) = UST5(I,J) + UH(I,J,L)
               VST5(I,J) = VST5(I,J) + VH(I,J,L)
               COUNT5(I,J) = 1
               GOTO 30
          ENDIF
         ENDDO
        ENDIF
30    CONTINUE
      ENDDO
      ENDDO

!
!$omp  parallel do
!$omp& private(umean6,vmean6,umean5,vmean5,umean1,vmean1,ushr,vshr)

      DO J=JSTART,JSTOP
      DO I=ISTART,ISTOP
         IF (COUNT6(I,J).GT.0.0 .AND. COUNT1(I,J) .GT. 0.0         &
            .AND. COUNT5(I,J) .GT. 0.0) THEN
           UMEAN5 = UST5(I,J) / COUNT5(I,J)
           VMEAN5 = VST5(I,J) / COUNT5(I,J)
           UMEAN1 = UST1(I,J) / COUNT1(I,J)
           VMEAN1 = VST1(I,J) / COUNT1(I,J)
           UMEAN6 = UST6(I,J) / COUNT6(I,J)
           VMEAN6 = VST6(I,J) / COUNT6(I,J)
           
!
!      COMPUTE STORM MOTION VECTOR
!      IT IS DEFINED AS 7.5 M/S TO THE RIGHT OF THE 0-6 KM MEAN
!      WIND CONSTRAINED ALONG A LINE WHICH IS BOTH PERPENDICULAR
!      TO THE 0-6 KM MEAN VERTICAL WIND SHEAR VECTOR AND PASSES
!      THROUGH THE 0-6 KM MEAN WIND.  THE WIND SHEAR VECTOR IS
!      SET AS THE DIFFERENCE BETWEEN THE 5.5-6 KM WIND (THE HEAD
!      OF THE SHEAR VECTOR) AND THE 0-0.5 KM WIND (THE TAIL).
!      THIS IS FOR THE RIGHT-MOVING CASE;  WE IGNORE THE LEFT MOVER.

           USHR = UMEAN5 - UMEAN1
           VSHR = VMEAN5 - VMEAN1

           UST(I,J) = UMEAN6 + (7.5*VSHR/SQRT(USHR*USHR+VSHR*VSHR))
           VST(I,J) = VMEAN6 - (7.5*USHR/SQRT(USHR*USHR+VSHR*VSHR))
         ELSE
           UST(I,J) = 0.0
           VST(I,J) = 0.0
        ENDIF
      ENDDO
      ENDDO
!
!       COMPUTE STORM-RELATIVE HELICITY
!
!$omp  parallel do
!$omp& private(du1,du2,dv1,dv2,dz,dz1,dz2,dzabv,ie,iw,z1,z2,z3)
      DO L = 2,LM-1
        if(GRIDTYPE /= 'A')then
          call exch(ZINT(1,jsta_2l,L))
          call exch(ZINT(1,jsta_2l,L+1))
        end if
	DO J=JSTART,JSTOP
        DO I=ISTART,ISTOP
          IW=I+IVW(J)
          IE=I+IVE(J)
          JN=J+JVN
          JS=J+JVS
	  IF (gridtype=='B')THEN
	   Z2=0.25*(ZMID(IW,J,L)+ZMID(IE,J,L)+                       &
                  ZMID(I,JN,L)+ZMID(IE,JN,L))                       
	  ELSE
	   Z2=0.25*(ZMID(IW,J,L)+ZMID(IE,J,L)+                       &
                  ZMID(I,JN,L)+ZMID(I,JS,L))
	  END IF	    
          DZABV=Z2-HTSFC(I,J)
!
          IF(DZABV.LT.DEPTH.AND.L.LE.NINT(LMV(I,J)))THEN
            IF (gridtype=='B')THEN
	      Z1=0.25*(ZMID(IW,J,L+1)+ZMID(IE,J,L+1)+                       &
                  ZMID(I,JN,L+1)+ZMID(IE,JN,L+1))
	      Z3=0.25*(ZMID(IW,J,L-1)+ZMID(IE,J,L-1)+                       &
                  ZMID(I,JN,L-1)+ZMID(IE,JN,L-1))	      	  
	      DZ=0.25*((ZINT(IW,J,L)+ZINT(IE,J,L)+                 &
                      ZINT(I,JN,L)+ZINT(IE,JN,L))-                &
                     (ZINT(IW,J,L+1)+ZINT(IE,J,L+1)+             &
                      ZINT(I,JN,L+1)+ZINT(IE,JN,L+1)))	  	                         
	    ELSE
	      Z1=0.25*(ZMID(IW,J,L+1)+ZMID(IE,J,L+1)+                       &
                  ZMID(I,JN,L+1)+ZMID(I,JS,L+1))
	      Z3=0.25*(ZMID(IW,J,L-1)+ZMID(IE,J,L-1)+                       &
                  ZMID(I,JN,L-1)+ZMID(I,JS,L-1))	 
              DZ=0.25*((ZINT(IW,J,L)+ZINT(IE,J,L)+                 &
                      ZINT(I,JS,L)+ZINT(I,JN,L))-                &
                     (ZINT(IW,J,L+1)+ZINT(IE,J,L+1)+             &
                      ZINT(I,JS,L+1)+ZINT(I,JN,L+1)))
	    END IF	      
            DZ1=Z1-Z2
            DZ2=Z2-Z3
            DU1=UH(I,J,L+1)-UH(I,J,L)
            DU2=UH(I,J,L)-UH(I,J,L-1)
            DV1=VH(I,J,L+1)-VH(I,J,L)
            DV2=VH(I,J,L)-VH(I,J,L-1)
            HELI(I,J)=((VH(I,J,L)-VST(I,J))*                     &
                      (DZ2*(DU1/DZ1)+DZ1*(DU2/DZ2))              &
                      -(UH(I,J,L)-UST(I,J))*                     &
                      (DZ2*(DV1/DZ1)+DZ1*(DV2/DZ2)))             &
                      *DZ/(DZ1+DZ2)+HELI(I,J) 
           ENDIF
        ENDDO
        ENDDO
      ENDDO
!
!     END OF ROUTINE.
!
      RETURN
      END
