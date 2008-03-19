      SUBROUTINE FRZLVL(ZFRZ,RHFRZ)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    FRZLVL      COMPUTES FRZING LVL Z AND RH
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-22       
!     
! ABSTRACT:
!     THIS ROUTINE COMPUTES THE FREEZING LEVEL HEIGHT AND RELATIVE
!     HUMIDITY AT THIS LEVEL FOR EACH MASS POINT ON THE ETA GRID.
!     THE COMPUTED FREEZING LEVEL HEIGHT IS THE MEAN SEA LEVEL
!     HEIGHT.  AT EACH MASS POINT WE MOVE UP FROM THE SURFACE TO 
!     FIND THE FIRST ETA LAYER WHERE THE TEMPERATURE IS LESS THAN
!     273.16K.  VERTICAL INTERPOLATION IN TEMPERATURE TO THE FREEZING
!     TEMPERATURE GIVES THE FREEZING LEVEL HEIGHT.  PRESSURE AND 
!     SPECIFIC HUMIDITY ARE INTERPOLATED TO THIS LEVEL AND ALONG WITH
!     THE TEMPERATURE PROVIDE THE FREEZING LEVEL RELATIVE HUMIDITY.
!     IF THE SURFACE (SKIN) TEMPERATURE IS BELOW FREEZING, THE ROUTINE
!     USES SURFACE BASED FIELDS TO COMPUTE THE RELATIVE HUMIDITY.
!     
!     NOTE THAT IN POSTING FREEZING LEVEL DATA THE LFM LOOK-ALIKE FILE
!     (IE, GRID 26), WE PACK 273.15K AS THE FREEZING TEMPERATURE.  ALL
!     OTHER OUTPUT GRIDS USE 273.16K
!   .     
!     
! PROGRAM HISTORY LOG:
!   92-12-22  RUSS TREADON
!   93-06-10  RUSS TREADON - CORRECTED FREEZING LEVEL HEIGHTS TO BE
!                            WITH REPSECT TO MEAN SEA LEVEL, NOT  
!                            ABOVE GROUND LEVEL.
!   98-06-15  T BLACK      - CONVERSION FROM 1-D TO 2-D
!   98-08-17  MIKE BALDWIN - COMPUTE RH OVER ICE IF NECESSARY
!   98-12-22  MIKE BALDWIN - BACK OUT RH OVER ICE
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   01-10-25  H CHUANG     - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-01-15  MIKE BALDWIN - WRF VERSION
!     
! USAGE:    CALL FRZLVL(ZFRZ,RHFRZ)
!   INPUT ARGUMENT LIST:
!     NONE     
!
!   OUTPUT ARGUMENT LIST: 
!     ZFRZ     - ABOVE GROUND LEVEL FREEZING HEIGHT.
!     RHFRZ    - RELATIVE HUMIDITY AT FREEZING LEVEL.
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       NONE
!     LIBRARY:
!       COMMON   - 
!                  LOOPS
!                  PVRBLS
!                  MASKS
!                  MAPOT
!                  POSTVAR
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
!     
!     
      use vrbls3d
      use vrbls2d
      use masks
      use params_mod
      use ctlblk_mod
!
      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL RHFRZ(IM,JM),ZFRZ(IM,JM)
      integer I,J,LLMH,L
      real HTSFC,PSFC,TSFC,QSFC,QSAT,RHSFC,DELZ,DELT,DELQ,DELALP,     &
           DELZP,ZL,DZABV,QFRZ,ALPL,ALPH,ALPFRZ,PFRZ,QSFRZ,RHZ,ZU,    &
           DZFR
!     
!*********************************************************************
!     START FRZLVL.
!
!
!     
!     LOOP OVER HORIZONTAL GRID.
!     
!$omp  parallel do
!$omp& private(alpfrz,alph,alpl,delalp,delq,delt,delz,
!$omp&         delzp,dzabv,dzfr,htsfc,l,llmh,psfc,qfrz,
!$omp&         qsat,qsfc,qsfrz,rhsfc,rhz,tsfc,
!$omp&         zl,zu)
      DO 20 J=JSTA,JEND
      DO 20 I=1,IM
         HTSFC    = FIS(I,J)*GI
         LLMH     = NINT(LMH(I,J))
         RHFRZ(I,J) = D00
         ZFRZ(I,J)  = HTSFC
         PSFC    = PINT(I,J,LLMH+1)
!     
!        CHECK IF FREEZING LEVEL IS AT THE GROUND.
!     
         IF(SM(I,J)/=SPVAL .AND. THZ0(I,J)/=SPVAL .AND.        &
      	   THS(I,J)/=SPVAL)THEN
          TSFC = (SM(I,J)*THZ0(I,J)+(1.-SM(I,J))*THS(I,J))     &
      	    *(PINT(I,J,NINT(LMH(I,J))+1)/P1000)**CAPA
         ELSE
! GFS analysis does not have flux file to retrieve TSFC from	 
	  TSFC=T(I,J,LM)+D0065*(ZMID(I,J,LM)-HTSFC)
	 END IF  
         IF (TSFC.LE.TFRZ) THEN
            ZFRZ(I,J) = HTSFC+(TSFC-TFRZ)/D0065
	    IF(SM(I,J)/=SPVAL .AND. QZ0(I,J)/=SPVAL .AND.      &
      	      QS(I,J)/=SPVAL)THEN
             QSFC    = SM(I,J)*QZ0(I,J)+(1.-SM(I,J))*QS(I,J)
! GFS does not output QS		   
            ELSE IF(QSHLTR(I,J)/=SPVAL)THEN
             QSFC=QSHLTR(I,J)
	    ELSE
	     QSFC=Q(I,J,LM) 
            END IF  
!
            QSAT=PQ0/PSFC*EXP(A2*(TSFC-A3)/(TSFC-A4))
!
            RHSFC   = QSFC/QSAT
            RHSFC   = AMAX1(0.01,RHSFC)
            RHSFC   = AMIN1(RHSFC,1.0)
            RHFRZ(I,J)= RHSFC
            GOTO 20
         ENDIF
!     
!        OTHERWISE, LOCATE THE FREEZING LEVEL ALOFT.
!
         DO 10 L = LLMH,1,-1
            IF (T(I,J,L).LE.TFRZ) THEN
               IF (L.LT.LLMH-1) THEN
                  DELZ = D50*(ZINT(I,J,L)-ZINT(I,J,L+2))
                  ZL   = D50*(ZINT(I,J,L+1)+ZINT(I,J,L+2))
                  DELT = T(I,J,L)-T(I,J,L+1)
                  ZFRZ(I,J) = ZL + (TFRZ-T(I,J,L+1))/DELT*DELZ
!     
                  DZABV = ZFRZ(I,J)-ZL
                  DELQ  = Q(I,J,L)-Q(I,J,L+1)
                  QFRZ  = Q(I,J,L+1) + DELQ/DELZ*DZABV
                  QFRZ  = AMAX1(0.0,QFRZ)
!     
!
                  ALPL   = ALPINT(I,J,L+2)
                  ALPH   = ALPINT(I,J,L)
                  DELALP = ALPH - ALPL
                  DELZP  = ZINT(I,J,L)-ZINT(I,J,L+2)
                  DZFR   = ZFRZ(I,J) - ZINT(I,J,L+2)
                  ALPFRZ = ALPL + DELALP/DELZP*DZFR
                  PFRZ   = EXP(ALPFRZ)
                  QSFRZ  = PQ0/PFRZ
!     
                  RHZ      = QFRZ/QSFRZ
                  RHZ      = AMAX1(0.01,RHZ)
                  RHZ      = AMIN1(RHZ,1.0)
                  RHFRZ(I,J) = RHZ
!     
               ELSE
                  ZU      = ZMID(I,J,L)
                  ZL      = HTSFC
                  DELZ    = ZU-ZL
		  IF(SM(I,J)/=SPVAL .AND. THZ0(I,J)/=SPVAL .AND.       &
      	            THS(I,J)/=SPVAL)THEN
                   TSFC    = (SM(I,J)*THZ0(I,J)+(1.-SM(I,J))*THS(I,J)) &
                     *(PINT(I,J,NINT(LMH(I,J))+1)/P1000)**CAPA
                  ELSE
! GFS analysis does not have flux file to retrieve TSFC from	 
	           TSFC=T(I,J,LM)+D0065*(ZMID(I,J,LM)-HTSFC)
	          END IF 		   
                  DELT    = T(I,J,L)-TSFC
                  ZFRZ(I,J) = ZL + (TFRZ-TSFC)/DELT*DELZ
!     
                  DZABV   = ZFRZ(I,J)-ZL
		  IF(SM(I,J)/=SPVAL .AND. QZ0(I,J)/=SPVAL .AND.       &
      	            QS(I,J)/=SPVAL)THEN
                   QSFC    = SM(I,J)*QZ0(I,J)+(1.-SM(I,J))*QS(I,J)
! GFS does not output QS		   
		  ELSE IF(QSHLTR(I,J)/=SPVAL)THEN
		   QSFC=QSHLTR(I,J)
		  ELSE
		   QSFC=Q(I,J,LM) 
		  END IF  
                  DELQ    = Q(I,J,L)-QSFC
                  QFRZ    = QSFC + DELQ/DELZ*DZABV
                  QFRZ    = AMAX1(0.0,QFRZ)
!     
                  ALPH    = ALPINT(I,J,L)
                  ALPL    = ALOG(PSFC)
                  DELALP  = ALPH-ALPL
                  ALPFRZ  = ALPL + DELALP/DELZ*DZABV
                  PFRZ    = EXP(ALPFRZ)
                  QSFRZ   = PQ0/PFRZ
!
                  RHZ     = QFRZ/QSFRZ
                  RHZ     = AMAX1(0.01,RHZ)
                  RHZ     = AMIN1(RHZ,1.0)
                  RHFRZ(I,J)= RHZ
               ENDIF
!     
!              BOUND FREEZING LEVEL RH.  FREEZING LEVEL HEIGHT IS
!              MEASURED WITH RESPECT TO MEAN SEA LEVEL.
!
               RHFRZ(I,J) = AMAX1(0.01,RHFRZ(I,J))
               RHFRZ(I,J) = AMIN1(RHFRZ(I,J),1.00)
               ZFRZ(I,J)  = AMAX1(0.0,ZFRZ(I,J))
               GOTO 20
            ENDIF
 10      CONTINUE
 20   CONTINUE
!     
!     END OF ROUTINE.
!     
      RETURN
      END
