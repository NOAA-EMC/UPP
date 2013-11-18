      SUBROUTINE CALWXT_EXPLICIT_POST(LMH,THS,PMID,PREC,SR,F_RIMEF,IWX)
! 
!     FILE: CALWXT.f
!     WRITTEN: 24 AUGUST 2005, G MANIKIN and B FERRIER 
!
!     ROUTINE TO COMPUTE PRECIPITATION TYPE USING EXPLICIT FIELDS
!       FROM THE MODEL MICROPHYSICS

      use params_mod, only: p1000, capa
      use ctlblk_mod, only: jsta, jend, modelname, pthresh, im, jsta_2l, jend_2u,&
               lm, jm
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
!  LIST OF VARIABLES NEEDED
!    PARAMETERS:
!
!    INPUT:
      real,dimension(im,jsta_2l:jend_2u,lm),intent(in):: F_RimeF, pmid
      REAL,dimension(im,jsta_2l:jend_2u),intent(in):: LMH, PREC, THS, SR
      integer,dimension(im,jm),intent(inout) :: IWX
      integer I,J,LMHK
      real PSFC,TSKIN,SNOW
!
!     ALLOCATE LOCAL STORAGE
!
!$omp  parallel do
      DO J=JSTA,JEND
      DO I=1,IM
        IWX(I,J) = 0
      ENDDO
      ENDDO

!GSM  THE RSM IS CURRENTLY INCOMPATIBLE WITH THIS ROUTINE
!GSM   ACCORDING TO B FERRIER, THERE MAY BE A WAY TO WRITE
!GSM   A VERSION OF THIS ALGORITHM TO WORK WITH THE RSM
!GSM   MICROPHYSICS, BUT IT DOESN'T EXIST AT THIS TIME
 
      IF (MODELNAME .EQ. 'RSM') GOTO 810 
!
!$omp  parallel do
!$omp& private(lmhk,psfc,tskin)
      DO 800 J=JSTA,JEND
      DO 800 I=1,IM
      LMHK=LMH(I,J)
!
!   SKIP THIS POINT IF NO PRECIP THIS TIME STEP 
!
      IF (PREC(I,J).LE.PTHRESH) GOTO 800
!
!  A SNOW RATIO LESS THAN 0.5 ELIMINATES SNOW AND SLEET
!   USE THE SKIN TEMPERATURE TO DISTINGUISH RAIN FROM FREEZING RAIN
!   NOTE THAT 2-M TEMPERATURE MAY BE A BETTER CHOICE IF THE MODEL
!   HAS A COLD BIAS FOR SKIN TEMPERATURE
! 
      IF (SR(I,J).LT.0.5) THEN
!        SURFACE (SKIN) POTENTIAL TEMPERATURE AND TEMPERATURE.
         PSFC=PMID(I,J,LMHK)
         TSKIN=THS(I,J)*(PSFC/P1000)**CAPA 

       IF (TSKIN.LT.273.15) THEN
!          FREEZING RAIN = 4
           IWX(I,J)=IWX(I,J)+4
       ELSE
!          RAIN = 8
           IWX(I,J)=IWX(I,J)+8
       ENDIF
       ELSE
!  
!  DISTINGUISH SNOW FROM SLEET WITH THE RIME FACTOR
! 
        IF(F_RimeF(I,J,LMHK).GE.10) THEN
!          SLEET = 2
           IWX(I,J)=IWX(I,J)+2
        ELSE
           SNOW = 1
           IWX(I,J)=IWX(I,J)+1 
       ENDIF
      ENDIF
 800  CONTINUE
 810  RETURN 
      END
