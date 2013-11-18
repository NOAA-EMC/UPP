      SUBROUTINE CALLCL(P1D,T1D,Q1D,PLCL,ZLCL)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALLCL      COMPUTES LCL HEIGHTS AND PRESSURE
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-03-15
!     
! ABSTRACT:
!     THIS ROUTINE COMPUTES THE LIFTING CONDENSATION LEVEL 
!     PRESSURE AND HEIGHT IN EACH COLUMN AT MASS POINTS.
!     THE HEIGHT IS ABOVE GROUND LEVEL.  THE EQUATION USED
!     TO FIND THE LCL PRESSURE IS FROM BOLTAN (1980,MWR) 
!     AND IS THE SAME AS THAT USED IN SUBROUTINE CALCAPE.
!     
!     THIS ROUTINE IS A TEST VERSION.  STILL TO BE RESOLVED
!     IS THE "BEST" PARCEL TO LIFT.
!   .     
!     
! PROGRAM HISTORY LOG:
!   93-03-15  RUSS TREADON
!   98-06-16  T BLACK - CONVERSION FROM 1-D TO 2-D
!   00-01-04  JIM TUCCILLO - MPI VERSION            
!   02-04-24  MIKE BALDWIN - WRF VERSION            
!     
! USAGE:    CALL CALLCL(P1D,T1D,Q1D,PLCL,ZLCL)
!   INPUT ARGUMENT LIST:
!     P1D      - ARRAY OF PARCEL PRESSURES (PA)
!     T1D      - ARRAY OF PARCEL TEMPERATURES (K)
!     Q1D      - ARRAY OF PARCEL SPECIFIC HUMIDITIES (KG/KG)
!
!   OUTPUT ARGUMENT LIST: 
!     PLCL     - PARCEL PRESSURE AT LCL (PA)
!     ZLCL     - PARCEL AGL HEIGHT AT LCL (M)
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       NONE
!     LIBRARY:
!       COMMON   - LOOPS
!                  OPTIONS
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN 90
!     MACHINE : CRAY C-90
!$$$  
!     
!     
      use vrbls3d, only: alpint, zint
      use vrbls2d, only: fis
      use masks, only: lmh
      use params_mod, only: eps, oneps, d01, h1m12, gi, d00
      use ctlblk_mod, only: jsta, jend, spval, jsta_m, jend_m, im, jm
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
      real,PARAMETER :: D35=3.5,D4805=4.805,H2840=2840.,H55=55.
      real,PARAMETER :: D2845=0.2845,D28=0.28
!
!     DECLARE VARIABLES.
!     
      REAL,dimension(IM,JM),intent(in) :: P1D,T1D,Q1D
      REAL,dimension(IM,JM),intent(inout) ::  PLCL,ZLCL
      REAL TLCL(IM,JM)
      integer I,J,L,LLMH
      real DLPLCL,ZSFC,DZ,DALP,ALPLCL,CKAPA,RMX,EVP,DENOM,ARG,RKAPA
!     
!**********************************************************************
!     START CALLCL HERE.
!     
!     LOAD OUTPUT ARRAYS WITH SPECIAL VALUE.
!     
      DO J=JSTA,JEND
      DO I=1,IM
        PLCL(I,J)=SPVAL
        TLCL(I,J)=SPVAL
        ZLCL(I,J)=SPVAL
      ENDDO
      ENDDO

!     
!     COMPUTE PRESSURE, TEMPERATURE AND AGL HEIGHT AT LCL.
!
      DO 30 J=JSTA_M,JEND_M
      DO 30 I=2,IM-1
      EVP      =P1D(I,J)*Q1D(I,J)/(EPS+ONEPS*Q1D(I,J))
      RMX      =EPS*EVP/(P1D(I,J)-EVP)
      CKAPA    =D2845*(1.-D28*RMX)
      RKAPA    =1./CKAPA
      ARG      =EVP*D01
      ARG      =AMAX1(H1M12,ARG)
      DENOM    =D35*ALOG(T1D(I,J))-ALOG(ARG)-D4805
      TLCL(I,J)=H2840/DENOM+H55
      PLCL(I,J)=P1D(I,J)*(TLCL(I,J)/T1D(I,J))**RKAPA
      ALPLCL   =ALOG(PLCL(I,J))
      LLMH     =NINT(LMH(I,J))
!
      DO 20 L=LLMH,1,-1
      IF(ALPINT(I,J,L).LT.ALPLCL)THEN
        DLPLCL   =ALPLCL-ALPINT(I,J,L+1)
        DALP     =ALPINT(I,J,L)-ALPINT(I,J,L+1)
        DZ       =ZINT(I,J,L)-ZINT(I,J,L+1)
        ZLCL(I,J)=ZINT(I,J,L+1)+DZ*DLPLCL/DALP
        ZSFC     =FIS(I,J)*GI
        ZLCL(I,J)=ZLCL(I,J)-ZSFC
        ZLCL(I,J)=AMAX1(D00,ZLCL(I,J))
        GOTO 30
      ENDIF
 20   CONTINUE
 30   CONTINUE
!     
!     END OF ROUTINE.
!     
      RETURN
      END
