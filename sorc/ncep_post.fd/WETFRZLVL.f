!> @file
!> wetfrzlvl() computes level of 0 wet bulb.
!>
!> This routine computes the lowest height with a wet bulb
!> temperature of freezing for each mass point on the eta grid.  
!> The computed wet bulb zero height is the mean sea level
!> height.  At each mass point we move up from the surface to 
!> find the first eta layer where the tw is less than
!> 273.16K.  Vertical interpolation in temperature to the freezing
!> temperature gives the freezing level height.  Pressure and 
!> specific humidity are interpolated to this level and along with
!> the temperature provide the freezing level relative humidity.
!> If the surface (skin) temperature is below freezing, the routine
!> uses surface based fields to compute the relative humidity.
!>
!> @param[in] TWET Wet bulb temperatures.
!> @param[out] ZWET Above ground level height of level with 0 wet bulb.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2003-11-14 | Geoff Manikin | Initial
!> 2004-12-06 | Geoff Manikin | Corrected computation of SFC temperature
!> 2005-03-11 | H CHUANG      | WRF Version
!> 2021-07-26 | W Meng        | Restrict computation from undefined grids
!>
!> @author Geoff Manikin W/NP2 @date 2003-11-14
      SUBROUTINE WETFRZLVL(TWET,ZWET)

!     
!     
      use vrbls3d, only: pint, zint, t
      use vrbls2d, only:  fis, thz0, ths
      use masks, only: lmh, sm
      use params_mod, only: gi, p1000, capa, tfrz, d0065, d50
      use ctlblk_mod, only: jsta, jend, im, jsta_2l, jend_2u, lm, spval
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL,intent(in) :: TWET(IM,JSTA_2L:JEND_2U,LM)
      REAL,intent(out) :: ZWET(IM,jsta:jend)
!     
      integer I,J,LLMH,L
      real HTSFC,THSFC,PSFC,TSFC,DELZ,DELT,ZL,ZU
!*********************************************************************
!     START FRZLVL.
!
!     LOOP OVER HORIZONTAL GRID.
!     
!!$omp  parallel do
!!$omp& private(delt,delz,htsfc,l,llmh
!!$omp&         tsfc,zl,zu)
      DO J=JSTA,JEND
      DO I=1,IM
         IF(FIS(I,J)==spval)THEN
           ZWET(I,J)=spval
           CYCLE
         ENDIF
         HTSFC     = FIS(I,J)*GI
         LLMH      = NINT(LMH(I,J))
         ZWET(I,J) = HTSFC
!     
!        CHECK IF FREEZING LEVEL IS AT THE GROUND.
!        IF YES, ESTIMATE UNDERGROUND FREEZING LEVEL USING 6.5C/KM LAPSE RATE
!        AND ASSUME RH TO BE EQUAL TO RH AT SFC
!     
         THSFC = (SM(I,J)*THZ0(I,J)+(1.-SM(I,J))*THS(I,J))
         PSFC  = PINT(I,J,LLMH+1)
         TSFC  = THSFC*(PSFC/P1000)**CAPA

         IF (TSFC<=TFRZ) THEN
!            ZWET(I,J) = HTSFC
            ZWET(I,J) = HTSFC+(TSFC-TFRZ)/D0065
            CYCLE   
         ENDIF
!     
!        OTHERWISE, LOCATE THE FREEZING LEVEL ALOFT.
!
         loopL:DO L = LLMH,1,-1
            IF (TWET(I,J,L)<=TFRZ) THEN
               IF (L<LLMH-1) THEN
                  DELZ = D50*(ZINT(I,J,L)-ZINT(I,J,L+2))
                  ZL   = D50*(ZINT(I,J,L+1)+ZINT(I,J,L+2))
                  DELT = TWET(I,J,L)-TWET(I,J,L+1)
                  ZWET(I,J) = ZL + (TFRZ-TWET(I,J,L+1))/DELT*DELZ
               ELSE
                  ZU      = D50*(ZINT(I,J,L)+ZINT(I,J,L+1))
                  ZL      = HTSFC
                  DELZ    = ZU-ZL
                  TSFC    = SM(I,J)*THZ0(I,J)+(1.-SM(I,J))*THS(I,J)    &
                   *(PINT(I,J,NINT(LMH(I,J))+1)/P1000)**CAPA
                  DELT    = T(I,J,L)-TSFC
		  IF(DELT /= 0.)THEN  
                   ZWET(I,J) = ZL + (TFRZ-TSFC)/DELT*DELZ
		  ELSE
		   ZWET(I,J) = HTSFC+(TSFC-TWET(I,J,L))/D0065
		  END IF  
                  IF (ZWET(I,J) > ZU) THEN
                    ZWET(I,J)=ZU
                  ENDIF
                   IF ((-1*ZWET(I,J)) > ZU) THEN
                    ZWET(I,J)=ZU
                  endif
               ENDIF
               EXIT loopL 
            ENDIF
         ENDDO loopL

      ENDDO !end I
      ENDDO !end J
!     
!     END OF ROUTINE.
!     
      RETURN
      END
