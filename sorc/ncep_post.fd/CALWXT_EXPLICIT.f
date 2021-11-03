      SUBROUTINE CALWXT_EXPLICIT_POST(LMH,THS,PMID,PREC,SR,F_RIMEF,IWX)
! 
!     FILE: CALWXT.f
!     WRITTEN: 24 AUGUST 2005, G MANIKIN and B FERRIER 
!
!     ROUTINE TO COMPUTE PRECIPITATION TYPE USING EXPLICIT FIELDS
!       FROM THE MODEL MICROPHYSICS
!
!     PROGRAM HISTORY LOG:
!     21-10-31  JESSE MENG - 2D DECOMPOSITION

      use params_mod, only: p1000, capa
      use ctlblk_mod, only: jsta, jend, modelname, pthresh, im, jsta_2l,  &
                            jend_2u, lm, ista, iend, ista_2l, iend_2u
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
!  LIST OF VARIABLES NEEDED
!    PARAMETERS:
!
!    INPUT:
      real,dimension(ista_2l:iend_2u,jsta_2l:jend_2u,lm),intent(in)    :: F_RimeF, pmid
      REAL,dimension(ista_2l:iend_2u,jsta_2l:jend_2u),   intent(in)    :: LMH, PREC, THS, SR
      integer,dimension(ista:iend,jsta:jend),      intent(inout) :: IWX
      integer I,J,LMHK
      real PSFC,TSKIN,SNOW
!
!     ALLOCATE LOCAL STORAGE
!
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IWX(I,J) = 0
        ENDDO
      ENDDO

!
!$omp  parallel do private(j,i,lmhk,psfc,tskin)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          LMHK=LMH(I,J)
!
!   SKIP THIS POINT IF NO PRECIP THIS TIME STEP 
!
          IF (PREC(I,J) <= PTHRESH) cycle
!
!  A SNOW RATIO LESS THAN 0.5 ELIMINATES SNOW AND SLEET
!   USE THE SKIN TEMPERATURE TO DISTINGUISH RAIN FROM FREEZING RAIN
!   NOTE THAT 2-M TEMPERATURE MAY BE A BETTER CHOICE IF THE MODEL
!   HAS A COLD BIAS FOR SKIN TEMPERATURE
! 
          IF (SR(I,J) < 0.5) THEN
!        SURFACE (SKIN) POTENTIAL TEMPERATURE AND TEMPERATURE.
            PSFC  = PMID(I,J,LMHK)
            TSKIN = THS(I,J)*(PSFC/P1000)**CAPA 

            IF (TSKIN < 273.15) THEN
!          FREEZING RAIN = 4
              IWX(I,J) = IWX(I,J)+4
            ELSE
!          RAIN = 8
              IWX(I,J) = IWX(I,J)+8
            ENDIF
          ELSE
!  
!  DISTINGUISH SNOW FROM SLEET WITH THE RIME FACTOR
! 
            IF(F_RimeF(I,J,LMHK) >= 10) THEN
!          SLEET = 2
              IWX(I,J) = IWX(I,J)+2
            ELSE
!              SNOW = 1
              IWX(I,J) = IWX(I,J)+1 
            ENDIF
          ENDIF
        enddo
      enddo
!
      RETURN 
      END
