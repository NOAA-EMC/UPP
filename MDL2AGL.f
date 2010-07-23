      SUBROUTINE MDL2AGL
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    MDL2P       VERT INTRP OF MODEL LVLS TO AGL HEIGHT
!   PRGRMMR: CHUANG           ORG: W/NP22     DATE: 05-05-23       
!     
! ABSTRACT:
!     FOR MOST APPLICATIONS THIS ROUTINE IS THE WORKHORSE
!     OF THE POST PROCESSOR.  IN A NUTSHELL IT INTERPOLATES
!     DATA FROM MODEL TO AGL HEIGHT SURFACES. 
!   .     
!     
! PROGRAM HISTORY LOG:
!   05-09-20  H CHUANG AND B ZHOU - ADD WIND DIFFERENCES OVER 2000 FT
!     
! USAGE:    CALL MDL2P
!   INPUT ARGUMENT LIST:
!
!   OUTPUT ARGUMENT LIST: 
!     NONE       
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       SCLFLD   - SCALE ARRAY ELEMENTS BY CONSTANT.
!       CALPOT   - COMPUTE POTENTIAL TEMPERATURE.
!       CALRH    - COMPUTE RELATIVE HUMIDITY.
!       CALDWP   - COMPUTE DEWPOINT TEMPERATURE.
!       BOUND    - BOUND ARRAY ELEMENTS BETWEEN LOWER AND UPPER LIMITS.
!       CALMCVG  - COMPUTE MOISTURE CONVERGENCE.
!       CALVOR   - COMPUTE ABSOLUTE VORTICITY.
!       CALSTRM  - COMPUTE GEOSTROPHIC STREAMFUNCTION.
!
!     LIBRARY:
!       COMMON   - CTLBLK
!                  RQSTFLD
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN 90
!     MACHINE : IBM SP
!$$$  
!
!
      use vrbls3d
      use vrbls2d
      use masks
      use params_mod
      use ctlblk_mod
      use rqstfld_mod
      use gridspec_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
      INCLUDE "mpif.h"
!     
!     INCLUDE MODEL DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
!     GAMMA AND RGAMOG ARE USED IN THE EXTRAPOLATION OF VIRTUAL
!     TEMPERATURES BEYOND THE UPPER OF LOWER LIMITS OF DATA.
!     
      integer,PARAMETER :: LAGL=2,LAGL2=1
!
!     DECLARE VARIABLES.
!     
      LOGICAL IOOMG,IOALL
      REAL UAGL(IM,JM),VAGL(IM,JM)
      REAL GRID1(IM,JM),GRID2(IM,JM)
!
      INTEGER NL1X(IM,JSTA:JEND), IHE(JM),IHW(JM),MAXLL,MINLL
      INTEGER LXXX,IERR
      INTEGER ISTART,ISTOP,JSTART,JSTOP
!
!
!--- Definition of the following 2D (horizontal) dummy variables
!
!  C1D   - total condensate
!  QW1   - cloud water mixing ratio
!  QI1   - cloud ice mixing ratio
!  QR1   - rain mixing ratio
!  QS1   - snow mixing ratio
!  DBZ1  - radar reflectivity
!  DBZR1 - radar reflectivity from rain
!  DBZI1 - radar reflectivity from ice (snow + graupel + sleet)
!  DBZC1 - radar reflectivity from parameterized convection (bogused)
!
!      REAL C1D(IM,JM),QW1(IM,JM),QI1(IM,JM),QR1(IM,JM)
!     &,    QS1(IM,JM) ,DBZ1(IM,JM)
      REAL DBZ1(IM,JM),DBZR1(IM,JM),DBZI1(IM,JM),DBZC1(IM,JM)       &
      ,    ZAGL(LAGL),ZAGL2(LAGL2),ZAGL3(LAGL2)                     &
! CRA
      ,    PAGL(IM,JM),TAGL(IM,JM),QAGL(IM,JM),PV,RHO
! CRA

     real PAGLU,PAGLL,TAGLU,TAGLL,QAGLU,QAGLL

     integer I,J,L,II,JJ,LP,LL,LLMH,ie,iw,jn,js
     real UAGLL,UAGLU,VAGLL,VAGLU,FACT,ZDUM
!
!     
!******************************************************************************
!
!     START MDL2P. 
!     
!     SET TOTAL NUMBER OF POINTS ON OUTPUT GRID.
!
!---------------------------------------------------------------
      ZAGL(1)=4000.
      ZAGL(2)=1000.
      ZAGL2(1)=609.6  ! 2000 ft
!
!     *** PART I ***
!
!     VERTICAL INTERPOLATION OF EVERYTHING ELSE.  EXECUTE ONLY
!     IF THERE'S SOMETHING WE WANT.
!
      IF (IGET(253).GT.0 .OR. IGET(279).GT.0 .OR. IGET(280).GT.0 .OR.   &
     &    IGET(281).GT.0 ) THEN
!
!---------------------------------------------------------------------
!***
!***  BECAUSE SIGMA LAYERS DO NOT GO UNDERGROUND,  DO ALL
!***  INTERPOLATION ABOVE GROUND NOW.
!***
!
        DO 310 LP=1,LAGL
         IF(LVLS(LP,IGET(253)).GT.0 .OR.LVLS(LP,IGET(279)).GT.0 .OR.    &
     &      LVLS(LP,IGET(280)).GT.0 .OR.LVLS(LP,IGET(281)).GT.0 ) THEN 
!
          jj=float(jsta+jend)/2.0
          ii=float(im)/3.0
          DO J=JSTA,JEND
          DO I=1,IM

!
	   DBZ1(I,J)=SPVAL
	   DBZR1(I,J)=SPVAL
	   DBZI1(I,J)=SPVAL
	   DBZC1(I,J)=SPVAL
!
!***  LOCATE VERTICAL INDEX OF MODEL MIDLAYER JUST BELOW
!***  THE AGL LEVEL TO WHICH WE ARE INTERPOLATING.
!
           LLMH=NINT(LMH(I,J))
           NL1X(I,J)=LLMH+1
           DO L=LLMH,2,-1
            ZDUM=ZMID(I,J,L)-ZINT(I,J,LLMH+1)
            IF(ZDUM.GE.ZAGL(LP))THEN
             NL1X(I,J)=L+1
	     GO TO 30
            ENDIF
           ENDDO
   30      CONTINUE	   
!
!  IF THE AGL LEVEL IS BELOW THE LOWEST MODEL MIDLAYER
!  BUT STILL ABOVE THE LOWEST MODEL BOTTOM INTERFACE,
!  WE WILL NOT CONSIDER IT UNDERGROUND AND THE INTERPOLATION
!  WILL EXTRAPOLATE TO THAT POINT
!
           IF(NL1X(I,J).EQ.(LLMH+1) .AND. ZAGL(LP).GT.0.)THEN
            NL1X(I,J)=LM
           ENDIF
!
!        if(NL1X(I,J).EQ.LMP1)print*,'Debug: NL1X=LMP1 AT '
!     1 ,i,j,lp
         ENDDO
         ENDDO
!
!mptest        IF(NHOLD.EQ.0)GO TO 310
!
!$omp  parallel do
!$omp& private(nn,i,j,ll,fact,qsat,rhl)
!hc        DO 220 NN=1,NHOLD
!hc        I=IHOLD(NN)
!hc        J=JHOLD(NN)
!        DO 220 J=JSTA,JEND
         DO 220 J=JSTA,JEND
         DO 220 I=1,IM
          LL=NL1X(I,J)
!---------------------------------------------------------------------
!***  VERTICAL INTERPOLATION OF GEOPOTENTIAL, TEMPERATURE, SPECIFIC
!***  HUMIDITY, CLOUD WATER/ICE, OMEGA, WINDS, AND TKE.
!---------------------------------------------------------------------
!
!HC        IF(NL1X(I,J).LE.LM)THEN
          LLMH = NINT(LMH(I,J))
          IF(NL1X(I,J).LE.LLMH)THEN
!
!---------------------------------------------------------------------
!          INTERPOLATE LINEARLY IN LOG(P)
!***  EXTRAPOLATE ABOVE THE TOPMOST MIDLAYER OF THE MODEL
!***  INTERPOLATION BETWEEN NORMAL LOWER AND UPPER BOUNDS
!***  EXTRAPOLATE BELOW LOWEST MODEL MIDLAYER (BUT STILL ABOVE GROUND)
!---------------------------------------------------------------------
!
!          FACT=(ALSL(LP)-ALOG(PMID(I,J,LL)))/
!     &         (ALOG(PMID(I,J,LL))-ALOG(PMID(I,J,LL-1)))
           ZDUM=ZAGL(LP)+ZINT(I,J,NINT(LMH(I,J))+1)
           FACT=(ZDUM-ZMID(I,J,LL))/(ZMID(I,J,LL)-ZMID(I,J,LL-1))
!	  
	 DBZ1(I,J)=DBZ(I,J,LL)+(DBZ(I,J,LL)-DBZ(I,J,LL-1))*FACT
	 DBZR1(I,J)=DBZR(I,J,LL)+(DBZR(I,J,LL)-DBZR(I,J,LL-1))*FACT
	 DBZI1(I,J)=DBZI(I,J,LL)+(DBZI(I,J,LL)-DBZI(I,J,LL-1))*FACT
	 DBZC1(I,J)=DBZC(I,J,LL)+(DBZC(I,J,LL)-DBZC(I,J,LL-1))*FACT
!           IF(I.eq.ii.and.j.eq.jj)print*,'Debug AGL RADAR REF',
!     &     i,j,ll,zagl(lp),ZINT(I,J,NINT(LMH(I,J))+1)
!     &      ,ZMID(I,J,LL-1),ZMID(I,J,LL)
!     &     ,DBZ(I,J,LL-1),DBZ(I,J,LL),DBZ1(I,J)
!     &     ,DBZR(I,J,LL-1),DBZR(I,J,LL),DBZR1(I,J)
!     &     ,DBZI(I,J,LL-1),DBZI(I,J,LL),DBZI1(I,J)
!     &     ,DBZC(I,J,LL-1),DBZC(I,J,LL),DBZC1(I,J)
	   DBZ1(I,J)=AMAX1(DBZ1(I,J),DBZmin)
	   DBZR1(I,J)=AMAX1(DBZR1(I,J),DBZmin)
	   DBZI1(I,J)=AMAX1(DBZI1(I,J),DBZmin)
	   DBZC1(I,J)=AMAX1(DBZC1(I,J),DBZmin)
!
! FOR UNDERGROUND AGL LEVELS, ASSUME TEMPERATURE TO CHANGE 
! ADIABATICLY, RH TO BE THE SAME AS THE AVERAGE OF THE 2ND AND 3RD
! LAYERS FROM THE GOUND, WIND TO BE THE SAME AS THE LOWEST LEVEL ABOVE
! GOUND
          ELSE
	   DBZ1(I,J)=DBZmin
	   DBZR1(I,J)=DBZmin
	   DBZI1(I,J)=DBZmin
	   DBZC1(I,J)=DBZmin
          END IF
  220    CONTINUE
!
!     
!---------------------------------------------------------------------
!        *** PART II ***
!---------------------------------------------------------------------
!
!        OUTPUT SELECTED FIELDS.
!
!---------------------------------------------------------------------
!
!
!---  Radar Reflectivity
          IF((IGET(253).GT.0) )THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=DBZ1(I,J)
             ENDDO
             ENDDO
             ID(1:25)=0
             ID(02)=129
             ID(11) = NINT(ZAGL(LP))
             CALL GRIBIT(IGET(253),LP,GRID1,IM,JM)
          END IF    
!---  Radar reflectivity from rain
          IF((IGET(279).GT.0) )THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=DBZR1(I,J)
             ENDDO
             ENDDO
             ID(1:25)=0
             ID(02)=129
             ID(11) = NINT(ZAGL(LP))
             CALL GRIBIT(IGET(279),LP,GRID1,IM,JM)
          END IF    
!---  Radar reflectivity from all ice habits (snow + graupel + sleet, etc.)
          IF((IGET(280).GT.0) )THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=DBZI1(I,J)
             ENDDO
             ENDDO
             ID(1:25)=0
             ID(02)=129
             ID(11) = NINT(ZAGL(LP))
             CALL GRIBIT(IGET(280),LP,GRID1,IM,JM)
          END IF    
!---  Radar reflectivity from parameterized convection
          IF((IGET(281).GT.0) )THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=DBZC1(I,J)
             ENDDO
             ENDDO
             ID(1:25)=0
             ID(02)=129
             ID(11) = NINT(ZAGL(LP))
             CALL GRIBIT(IGET(281),LP,GRID1,IM,JM)
          END IF    
!          
         ENDIF ! FOR LEVEL
!     
!***  END OF MAIN VERTICAL LOOP
!     
  310   CONTINUE
!***  ENDIF FOR IF TEST SEEING IF WE WANT ANY OTHER VARIABLES
!
      ENDIF
! SRD
!---  Max Derived Radar Reflectivity
          IF((IGET(421).GT.0) )THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=REFD_MAX(I,J)
             ENDDO
             ENDDO
             ID(1:25)=0
             ID(02)=129
             ID(9)=105
             ID(11) = 10
             ID(20) = 2
             ID(19) = IFHR
             IF (IFHR.EQ.0) THEN
               ID(18) = 0
             ELSE
               ID(18) = IFHR - 1
             ENDIF
             CALL GRIBIT(IGET(421),LP,GRID1,IM,JM)
          END IF

!---  Max Updraft Helicity
          IF((IGET(420).GT.0) )THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=UP_HELI_MAX(I,J)
             ENDDO
             ENDDO
             ID(1:25)=0
             ID(02)=129
!             ID(11) = NINT(ZAGL(2))
             ID(9) = 106
             ID(10) = 50
             ID(11) = 20
             ID(20) = 2
             ID(19) = IFHR
             IF (IFHR.EQ.0) THEN
               ID(18) = 0
             ELSE
               ID(18) = IFHR - 1
             ENDIF
             CALL GRIBIT(IGET(420),LP,GRID1,IM,JM)
          END IF

!---  Max Column Integrated Graupel
          IF((IGET(429).GT.0) )THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=GRPL_MAX(I,J)
             ENDDO
             ENDDO
             ID(1:25)=0
             ID(02)=129
             ID(20) = 2
             ID(19) = IFHR
             IF (IFHR.EQ.0) THEN
               ID(18) = 0
             ELSE
               ID(18) = IFHR - 1
             ENDIF
             CALL GRIBIT(IGET(429),LP,GRID1,IM,JM)
          END IF
! SRD

!
      IF((IGET(259).GT.0) )THEN
!
!---------------------------------------------------------------------
!***
!***  BECAUSE SIGMA LAYERS DO NOT GO UNDERGROUND,  DO ALL
!***  INTERPOLATION ABOVE GROUND NOW.
!***
!
        DO 320 LP=1,LAGL2
         IF(LVLS(LP,IGET(259)).GT.0)THEN 
!
          jj=(jsta+jend)/2
          ii=(im)/2
          DO J=JSTA,JEND
          DO I=1,IM

!
	   UAGL(I,J)=SPVAL
	   VAGL(I,J)=SPVAL
!
!***  LOCATE VERTICAL INDEX OF MODEL MIDLAYER JUST BELOW
!***  THE AGL LEVEL TO WHICH WE ARE INTERPOLATING.
!
           LLMH=NINT(LMH(I,J))
           NL1X(I,J)=LLMH+1
           DO L=LLMH,2,-1
            ZDUM=ZMID(I,J,L)-ZINT(I,J,LLMH+1)
            IF(ZDUM.GE.ZAGL2(LP))THEN
             NL1X(I,J)=L+1
	     GO TO 40
            ENDIF
           ENDDO
   40      CONTINUE	   
!
!  IF THE AGL LEVEL IS BELOW THE LOWEST MODEL MIDLAYER
!  BUT STILL ABOVE THE LOWEST MODEL BOTTOM INTERFACE,
!  WE WILL NOT CONSIDER IT UNDERGROUND AND THE INTERPOLATION
!  WILL EXTRAPOLATE TO THAT POINT
!
           IF(NL1X(I,J).EQ.(LLMH+1) .AND. ZAGL2(LP).GT.0.)THEN
            NL1X(I,J)=LM
           ENDIF
!
!        if(NL1X(I,J).EQ.LMP1)print*,'Debug: NL1X=LMP1 AT '
!     1 ,i,j,lp
         ENDDO
         ENDDO
!
!mptest        IF(NHOLD.EQ.0)GO TO 310
!
!$omp  parallel do
!$omp& private(nn,i,j,ll,fact,qsat,rhl)
!hc        DO 220 NN=1,NHOLD
!hc        I=IHOLD(NN)
!hc        J=JHOLD(NN)
!        DO 220 J=JSTA,JEND
         DO J=JSTA,JEND
          IF(gridtype=='A')THEN
           IHW(J)=-1
           IHE(J)=1 
          ELSE IF(gridtype=='E')THEN
           IHW(J)=-MOD(J,2)
           IHE(J)=IHW(J)+1
          END IF	
         ENDDO
	 IF(global)then
	   ISTART=1
           ISTOP=IM
           JSTART=JSTA
           JSTOP=JEND
	 ELSE
	   ISTART=2
           ISTOP=IM-1
           JSTART=JSTA_M
           JSTOP=JEND_M
	 END IF    
	 IF(gridtype/='A')THEN 
!	  MAXLL=maxval(NL1X)
	  MINLL=minval(NL1X)
	  print*,'MINLL before all reduce= ',MINLL
	  CALL MPI_ALLREDUCE(MINLL,LXXX,1,MPI_INTEGER,MPI_MIN,MPI_COMM_COMP,IERR)
	  MINLL=LXXX
	  print*,'exchange wind in MDL2AGL from ',MINLL
	  DO LL=MINLL,LM
	   call exch(UH(1:IM,JSTA_2L:JEND_2U,LL))
	   call exch(VH(1:IM,JSTA_2L:JEND_2U,LL))
	  END DO
	 END IF   
         DO 230 J=JSTART,JSTOP
         DO 230 I=ISTART,ISTOP
          LL=NL1X(I,J)
!---------------------------------------------------------------------
!***  VERTICAL INTERPOLATION OF GEOPOTENTIAL, TEMPERATURE, SPECIFIC
!***  HUMIDITY, CLOUD WATER/ICE, OMEGA, WINDS, AND TKE.
!---------------------------------------------------------------------
!
!HC        IF(NL1X(I,J).LE.LM)THEN
          LLMH = NINT(LMH(I,J))
          IF(NL1X(I,J).LE.LLMH)THEN
!
!---------------------------------------------------------------------
!          INTERPOLATE LINEARLY IN LOG(P)
!***  EXTRAPOLATE ABOVE THE TOPMOST MIDLAYER OF THE MODEL
!***  INTERPOLATION BETWEEN NORMAL LOWER AND UPPER BOUNDS
!***  EXTRAPOLATE BELOW LOWEST MODEL MIDLAYER (BUT STILL ABOVE GROUND)
!---------------------------------------------------------------------
!
!          FACT=(ALSL(LP)-ALOG(PMID(I,J,LL)))/
!     &         (ALOG(PMID(I,J,LL))-ALOG(PMID(I,J,LL-1)))
           ZDUM=ZAGL2(LP)+ZINT(I,J,NINT(LMH(I,J))+1)
           FACT=(ZDUM-ZMID(I,J,LL))/(ZMID(I,J,LL)-ZMID(I,J,LL-1))
!	  
           IF(gridtype=='A')THEN
	    UAGLU=UH(I,J,LL-1)
            UAGLL=UH(I,J,LL)
     
            VAGLU=VH(I,J,LL-1)
            VAGLL=VH(I,J,LL)
	   ELSE IF(gridtype=='E')THEN
            UAGLU=(UH(I+IHE(J),J,LL-1)+UH(I+IHW(J),J,LL-1)+         &
     &	       UH(I,J-1,LL-1)+UH(I,J+1,LL-1))/4.0
            UAGLL=(UH(I+IHE(J),J,LL)+UH(I+IHW(J),J,LL)+             &
     &	       UH(I,J-1,LL)+UH(I,J+1,LL))/4.0
     
            VAGLU=(VH(I+IHE(J),J,LL-1)+VH(I+IHW(J),J,LL-1)+         &
     &	       VH(I,J-1,LL-1)+VH(I,J+1,LL-1))/4.0
            VAGLL=(VH(I+IHE(J),J,LL)+VH(I+IHW(J),J,LL)+             &
     &	       VH(I,J-1,LL)+VH(I,J+1,LL))/4.0
           ELSE IF(gridtype=='B')THEN
	    IE=I
            IW=I-1
            JN=J
            JS=J-1
            UAGLU=(UH(IE,J,LL-1)+UH(IW,J,LL-1)+         &
     &	       UH(IE,JS,LL-1)+UH(IW,JS,LL-1))/4.0
            UAGLL=(UH(IE,J,LL)+UH(IW,J,LL)+         &
     &	       UH(IE,JS,LL)+UH(IW,JS,LL))/4.0
     
            VAGLU=(VH(IE,J,LL-1)+VH(IW,J,LL-1)+         &
     &	       VH(IE,JS,LL-1)+VH(IW,JS,LL-1))/4.0
            VAGLL=(VH(IE,J,LL)+VH(IW,J,LL)+         &
     &	       VH(IE,JS,LL)+VH(IW,JS,LL))/4.0
           END IF
           UAGL(I,J)=UAGLL+(UAGLL-UAGLU)*FACT
	   VAGL(I,J)=VAGLL+(VAGLL-VAGLU)*FACT
           IF(I==II.AND.J==JJ)PRINT*,                                &
     &	   'DEBUG LLWS: I,J,NL1X,UU,UL,VU,VL,ZSFC,ZMIDU,ZMIDL,U,V= ' &
     &,     i,j,ll,UAGLU,UAGLL,VAGLU,VAGLL,ZINT(I,J,NINT(LMH(I,J))+1)&
     &,     ZMID(I,J,LL-1),ZMID(I,J,LL),UAGL(I,J),VAGL(I,J)          &
     &,     U10(I,J),V10(I,J)
!
! FOR UNDERGROUND AGL LEVELS, ASSUME TEMPERATURE TO CHANGE 
! ADIABATICLY, RH TO BE THE SAME AS THE AVERAGE OF THE 2ND AND 3RD
! LAYERS FROM THE GOUND, WIND TO BE THE SAME AS THE LOWEST LEVEL ABOVE
! GOUND
          ELSE
	   IF(gridtype=='A')THEN
            UAGL(I,J)=UH(I,J,NINT(LMV(I,J)))  
	    VAGL(I,J)=VH(I,J,NINT(LMV(I,J)))
	   ELSE IF(gridtype=='E')THEN
	    UAGL(I,J)=(UH(I+IHE(J),J,NINT(LMV(I+IHE(J),J)))          &
     &	      +UH(I+IHW(J),J,NINT(LMV(I+IHW(J),J)))+                 &
     &	       UH(I,J-1,NINT(LMV(I,J-1)))+UH(I,J+1,NINT(LMV(I,J+1))))/4.0   
	    VAGL(I,J)=(VH(I+IHE(J),J,NINT(LMV(I+IHE(J),J)))          &
     &	      +VH(I+IHW(J),J,NINT(LMV(I+IHW(J),J)))+                 &
     &	       VH(I,J-1,NINT(LMV(I,J-1)))+VH(I,J+1,NINT(LMV(I,J+1))))/4.0
           ELSE IF(gridtype=='B')THEN
	    IE=I
            IW=I-1
            JN=J
            JS=J-1
	    UAGL(I,J)=(UH(IE,J,NINT(LMV(IE,J)))          &
     &	      +UH(IW,J,NINT(LMV(IW,J)))+                 &
     &	       UH(IE,JS,NINT(LMV(IE,JS)))+UH(IW,JS,NINT(LMV(IW,JS))))/4.0   
	    VAGL(I,J)=(VH(IE,J,NINT(LMV(IE,J)))          &
     &	      +VH(IW,J,NINT(LMV(IW,J)))+                 &
     &	       VH(IE,JS,NINT(LMV(IE,JS)))+VH(IW,JS,NINT(LMV(IW,JS))))/4.0
           END IF
          END IF
  230    CONTINUE
!
!     
!---------------------------------------------------------------------
!        *** PART II ***
!---------------------------------------------------------------------
!
!        OUTPUT SELECTED FIELDS.
!
!---------------------------------------------------------------------
!
!
!---  Wind Shear (wind speed difference in knots between sfc and 2000 ft)

	     DO J=JSTA,JEND
             DO I=1,IM
	       IF(ABS(UAGL(I,J)-SPVAL).GT.SMALL .AND.               &
                  ABS(VAGL(I,J)-SPVAL).GT.SMALL)THEN        
                GRID1(I,J)=SQRT((UAGL(I,J)-U10(I,J))**2+            &
      	          (VAGL(I,J)-V10(I,J))**2)*1.943*ZAGL2(LP)/         &
                  (ZAGL2(LP)-10.)
               ELSE
	        GRID1(I,J)=SPVAL
	       END IF	 
             ENDDO
             ENDDO
             ID(1:25)=0
	     ID(10) = NINT(ZAGL2(LP))
             ID(11) = 0
             CALL GRIBIT(IGET(259),LP,GRID1,IM,JM)
!          
         ENDIF ! FOR LEVEL
!     
!***  END OF MAIN VERTICAL LOOP
!     
  320   CONTINUE
!***  ENDIF FOR IF TEST SEEING IF WE WANT ANY OTHER VARIABLES
!
      ENDIF
! CRA
      IF (IGET(411).GT.0 .OR. IGET(412).GT.0 .OR. IGET(413).GT.0) THEN
!
!---------------------------------------------------------------------
!***
!***  BECAUSE SIGMA LAYERS DO NOT GO UNDERGROUND,  DO ALL
!***  INTERPOLATION ABOVE GROUND NOW.
!***
!
        DO 330 LP=1,LAGL2
         IF(LVLS(LP,IGET(411)).GT.0 .OR. LVLS(LP,IGET(412)).GT.0     &
     &      .OR. LVLS(LP,IGET(413)).GT.0 ) THEN

!
          jj=float(jsta+jend)/2.0
          ii=float(im)/3.0
          DO J=JSTA_2L,JEND_2U
          DO I=1,IM

!
           PAGL(I,J)=SPVAL
           TAGL(I,J)=SPVAL
           QAGL(I,J)=SPVAL
           UAGL(I,J)=SPVAL
           VAGL(I,J)=SPVAL
!
!***  LOCATE VERTICAL INDEX OF MODEL MIDLAYER JUST BELOW
!***  THE AGL LEVEL TO WHICH WE ARE INTERPOLATING.
!
           LLMH=NINT(LMH(I,J))
           NL1X(I,J)=LLMH+1
           DO L=LLMH,2,-1
            ZDUM=ZMID(I,J,L)-ZINT(I,J,LLMH+1)
            IF(ZDUM.GE.ZAGL3(LP))THEN
             NL1X(I,J)=L+1
             GO TO 50
            ENDIF
           ENDDO
   50      CONTINUE
!
!  IF THE AGL LEVEL IS BELOW THE LOWEST MODEL MIDLAYER
!  BUT STILL ABOVE THE LOWEST MODEL BOTTOM INTERFACE,
!  WE WILL NOT CONSIDER IT UNDERGROUND AND THE INTERPOLATION
!  WILL EXTRAPOLATE TO THAT POINT
!
           IF(NL1X(I,J).EQ.(LLMH+1) .AND. ZAGL3(LP).GT.0.)THEN
            NL1X(I,J)=LM
           ENDIF
!
!        if(NL1X(I,J).EQ.LMP1)print*,'Debug: NL1X=LMP1 AT '
!     1 ,i,j,lp
         ENDDO
         ENDDO
!
!mptest        IF(NHOLD.EQ.0)GO TO 310
!
!$omp  parallel do
!$omp& private(nn,i,j,ll,fact,qsat,rhl)
!chc        DO 220 NN=1,NHOLD
!chc        I=IHOLD(NN)
!chc        J=JHOLD(NN)
!        DO 220 J=JSTA,JEND
         DO 240 J=JSTA_2L,JEND_2U
         DO 240 I=1,IM
          LL=NL1X(I,J)
!---------------------------------------------------------------------
!***  VERTICAL INTERPOLATION OF GEOPOTENTIAL, TEMPERATURE, SPECIFIC
!***  HUMIDITY, CLOUD WATER/ICE, OMEGA, WINDS, AND TKE.
!---------------------------------------------------------------------
!
!CHC        IF(NL1X(I,J).LE.LM)THEN
          LLMH = NINT(LMH(I,J))
          IF(NL1X(I,J).LE.LLMH)THEN
!
!---------------------------------------------------------------------
!          INTERPOLATE LINEARLY IN LOG(P)
!***  EXTRAPOLATE ABOVE THE TOPMOST MIDLAYER OF THE MODEL
!***  INTERPOLATION BETWEEN NORMAL LOWER AND UPPER BOUNDS
!***  EXTRAPOLATE BELOW LOWEST MODEL MIDLAYER (BUT STILL ABOVE GROUND)
!---------------------------------------------------------------------
!
!          FACT=(ALSL(LP)-ALOG(PMID(I,J,LL)))/
!     &         (ALOG(PMID(I,J,LL))-ALOG(PMID(I,J,LL-1)))
           ZDUM=ZAGL3(LP)+ZINT(I,J,NINT(LMH(I,J))+1)
           FACT=(ZDUM-ZMID(I,J,LL))                             &
              /(ZMID(I,J,LL)-ZMID(I,J,LL-1))
!
           PAGLU=ALOG(PMID(I,J,LL-1))
           PAGLL=ALOG(PMID(I,J,LL))

           TAGLU=T(I,J,LL-1)
           TAGLL=T(I,J,LL)

           QAGLU=Q(I,J,LL-1)
           QAGLL=Q(I,J,LL)

           UAGLU=UH(I,J,LL-1)
           UAGLL=UH(I,J,LL)

           VAGLU=VH(I,J,LL-1)
           VAGLL=VH(I,J,LL)

           PAGL(I,J)=EXP(PAGLL+(PAGLL-PAGLU)*FACT)
           TAGL(I,J)=TAGLL+(TAGLL-TAGLU)*FACT
           QAGL(I,J)=QAGLL+(QAGLL-TAGLU)*FACT
           UAGL(I,J)=UAGLL+(UAGLL-UAGLU)*FACT
           VAGL(I,J)=VAGLL+(VAGLL-VAGLU)*FACT
!
! FOR UNDERGROUND AGL LEVELS, ASSUME TEMPERATURE TO CHANGE
! ADIABATICLY, RH TO BE THE SAME AS THE AVERAGE OF THE 2ND AND 3RD
! LAYERS FROM THE GOUND, WIND TO BE THE SAME AS THE LOWEST LEVEL ABOVE
! GOUND
          ELSE
            PAGL(I,J)=PMID(I,J,NINT(LMV(I,J)))
            TAGL(I,J)=T(I,J,NINT(LMV(I,J)))
            QAGL(I,J)=Q(I,J,NINT(LMV(I,J)))
            UAGL(I,J)=UH(I,J,NINT(LMV(I,J)))
            VAGL(I,J)=VH(I,J,NINT(LMV(I,J)))
          END IF
  240 CONTINUE
!
!
!---------------------------------------------------------------------
!        *** PART II ***
!---------------------------------------------------------------------
!
!        OUTPUT SELECTED FIELDS.
!
!---------------------------------------------------------------------
!
!
!---  Wind Energy Potential -- 0.5 * moist air density * wind speed^3
          IF((IGET(411).GT.0) ) THEN
            DO J=JSTA,JEND
            DO I=1,IM
              QAGL(I,J)=QAGL(I,J)/1000.0
              PV=QAGL(I,J)*PAGL(I,J)/(EPS*(1-QAGL(I,J)) + QAGL(I,J))
              RHO=(1/TAGL(I,J))*(((PAGL(I,J)-PV)/RD) + PV/461.495)
              GRID1(I,J)=0.5*RHO*(SQRT(UAGL(I,J)**2+VAGL(I,J)**2))**3
            ENDDO
            ENDDO
            ID(1:25)=0
            ID(11) = NINT(ZAGL3(LP))
            CALL GRIBIT(IGET(411),LP,GRID1,IM,JM)
          ENDIF
!--- U Component of wind
          IF((IGET(412).GT.0) ) THEN
            DO J=JSTA,JEND
            DO I=1,IM
              GRID1(I,J)=UAGL(I,J)
            ENDDO
            ENDDO
            ID(1:25)=0
            ID(11) = NINT(ZAGL3(LP))
            CALL GRIBIT(IGET(412),LP,GRID1,IM,JM)
          ENDIF
!--- V Component of wind
          IF((IGET(413).GT.0) ) THEN
            DO J=JSTA,JEND
            DO I=1,IM
              GRID1(I,J)=VAGL(I,J)
            ENDDO
            ENDDO
            ID(1:25)=0
            ID(11) = NINT(ZAGL3(LP))
            CALL GRIBIT(IGET(413),LP,GRID1,IM,JM)
          ENDIF
!
         ENDIF ! FOR LEVEL

!
!***  END OF MAIN VERTICAL LOOP
!
  330   CONTINUE
!***  ENDIF FOR IF TEST SEEING IF WE WANT ANY OTHER VARIABLES
!
      ENDIF
! CRA
!
!     END OF ROUTINE.
!
      RETURN
      END
