!> @file
!> @brief Subroutine that computes FRZING LVL, Z and RH.
!>
!> This routine computes the isothermal level height and relative
!> humidity at this level for each mass point on the ETA grid.
!> The computed isothermal level height is the mean sea level
!> height. At each mass point we move up from the surface to  
!> find the last ETA layer where the temperature is less than
!> isotherm and the temp in the layer below is above isotherm.
!> Vertical interpolation in temperature to the isotherm
!> temperature gives the isothermal level height. Pressure and   
!> specific humidity are interpolated to this level and along with
!> the temperature provide the isothermal level relative humidity.
!> If the entire atmosphere is below isotherm, the routine
!> uses surface based fields to compute the relative humidity.
!> 
!> Note that in posting freezing level data the LFM look-alike file
!> (IE, GRID 26), we pack 273.15K as the freezing temperature. All 
!> other output grids use 273.16K.
!>
!> @param[in] isotherm isothermal value of height to be output.
!> @param[out] ZFRZ Above ground level/ZFL at isotherm height.
!> @param[out] RHFRZ Relative humidity at isotherm level.
!> @param[out] PFRZL pressure at isotherm level.
!> 
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-22 | Russ Treadon | Initial
!> 1993-06-05 | Russ Treadon | Corrected freezing level heights to be with respect to mean sea level, not above ground level
!> 1995-03-10 | Mike Baldwin | Get highest freezing level
!> 1998-06-15 | T Black      | Conversion from 1-D to 2-D
!> 1998-08-17 | Mike Baldwin | Compute RH over ice if necessary
!> 1998-12-22 | Mike Baldwin | Back out RH over ice
!> 2000-01-04 | Jim Tuccillo | MPI version           
!> 2001-10-25 | H Chuang     | Modified to process hybrid model output
!> 2002-01-15 | Mike Baldwin | WRF version
!> 2010-08-27 | T. Smirnova  | Added PFRZL to the output
!> 2016-01-21 | C. Alexander | Generalized function for any isotherm
!> 2019-10-30 | Bo Cui       | Remove "GOTO" statement
!> 2020-11-10 | Jesse Meng   | Use UPP_PHYSICS module
!> 2021-07-28 | W. Meng      | Restrict compuatation from undefined grids
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE FRZLVL2(ISOTHERM,ZFRZ,RHFRZ,PFRZL)

!     
      use vrbls3d, only: pint, t, zmid, pmid, q, zint, alpint
      use vrbls2d, only: fis, tshltr, pshltr, qz0, qs, qshltr
      use masks, only: lmh, sm
      use params_mod, only: gi, d00, capa, d0065, tfrz, pq0, a2, a3, a4, d50
      use ctlblk_mod, only: jsta, jend, spval, lm, modelname, im
      use physcons_post, only: con_rd, con_rv, con_eps, con_epsm1
      use upp_physics, only: FPVSNEW

      implicit none

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      implicit none
!     
!     DECLARE VARIABLES.
!
      REAL,PARAMETER::PUCAP=300.0E2
      real,intent(in)                   ::  ISOTHERM
      REAL,dimension(im,jsta:jend),intent(out) ::  RHFRZ, ZFRZ, PFRZL
!jw
      integer I,J,L,LICE,LLMH
      real HTSFC,PSFC,QSFC,RHSFC,QW,QSAT,DELZ,DELT,DELQ,DELALP,DELZP,  &
           ZL,ZU,DZABV,QFRZ,ALPL,ALPH,ALPFRZ,PFRZ,QSFRZ,RHZ,DZFR, &
           TSFC,es
!     
!*********************************************************************
!     START FRZLVL.
!
!     LOOP OVER HORIZONTAL GRID.
!     

      DO 20 J=JSTA,JEND
      DO 20 I=1,IM
         IF(FIS(I,J)<spval)THEN
         HTSFC    = FIS(I,J)*GI
         LLMH     = NINT(LMH(I,J))
         RHFRZ(I,J) = D00
         ZFRZ(I,J)  = HTSFC
         PSFC     = PINT(I,J,LLMH)
         PFRZL(I,J) = PSFC
!     
!        FIND THE HIGHEST LAYER WHERE THE TEMPERATURE
!        CHANGES FROM ABOVE TO BELOW ISOTHERM.
!     
!         TSFC = (SM(I,J)*THZ0(I,J)+(1.-SM(I,J))*THS(I,J))    &
!         	 *(PINT(I,J,NINT(LMH(I,J))+1)/P1000)**CAPA	 
         IF(TSHLTR(I,J)/=SPVAL .AND. PSHLTR(I,J)/=SPVAL)THEN
          TSFC=TSHLTR(I,J)*(PSHLTR(I,J)*1.E-5)**CAPA
         ELSE
! GFS analysis does not have flux file to retrieve TSFC from
          TSFC=T(I,J,LM)+D0065*(ZMID(I,J,LM)-HTSFC-2.0)
         END IF
         LICE=LLMH
! Per AWC's request, put a 300 mb cap for highest isothermal level so that it
! does not go into stratosphere
         DO L = LLMH-1,1,-1
            IF (PMID(I,J,L)>=PUCAP .AND. &
	    (T(I,J,L)<=ISOTHERM.AND.T(I,J,L+1)>ISOTHERM))LICE=L
         ENDDO
!     
!        CHECK IF ISOTHERM LEVEL IS AT THE GROUND.
!     
         IF (LICE==LLMH.AND.TSFC<=ISOTHERM) THEN
            ZFRZ(I,J) = HTSFC+2.0+(TSFC-ISOTHERM)/D0065
            QSFC    = SM(I,J)*QZ0(I,J)+(1.-SM(I,J))*QS(I,J)
            IF(QSHLTR(I,J)/=SPVAL)THEN
             PSFC=PSHLTR(I,J)
             QSFC=QSHLTR(I,J)
            ELSE
             QSFC=Q(I,J,LM)
             PSFC=PMID(I,J,LM)
            END IF
            PFRZL(I,J) = PSFC
!
            IF(MODELNAME == 'GFS' .OR. MODELNAME == 'RAPR')THEN
             ES=FPVSNEW(TSFC)
             ES=MIN(ES,PSFC)
             QSAT=CON_EPS*ES/(PSFC+CON_EPSM1*ES)
            ELSE
             QSAT=PQ0/PSFC  &
               *EXP(A2*(TSFC-A3)/(TSFC-A4))
            END IF
!
            RHSFC   = QSFC/QSAT
            RHSFC   = AMAX1(0.01,RHSFC)
            RHSFC   = AMIN1(RHSFC,1.0)
            RHFRZ(I,J)= RHSFC
!     
!        OTHERWISE, LOCATE THE ISOTHERM LEVEL ALOFT.
!
         ELSE IF (LICE<LLMH) THEN
                  L=LICE
                  DELZ = D50*(ZINT(I,J,L)-ZINT(I,J,L+2))
                  ZL   = D50*(ZINT(I,J,L+1)+ZINT(I,J,L+2))
                  DELT = T(I,J,L)-T(I,J,L+1)
                  ZFRZ(I,J) = ZL+(ISOTHERM-T(I,J,L+1))/DELT*DELZ
!     
                  DZABV = ZFRZ(I,J)-ZL
                  DELQ  = Q(I,J,L)-Q(I,J,L+1)
                  QFRZ  = Q(I,J,L+1) + DELQ/DELZ*DZABV
                  QFRZ  = AMAX1(0.0,QFRZ)
!     
                  ALPL   = ALPINT(I,J,L+2)
                  ALPH   = ALPINT(I,J,L)
                  DELALP = ALPH - ALPL
                  DELZP  = ZINT(I,J,L)-ZINT(I,J,L+2)
                  DZFR   = ZFRZ(I,J) - ZINT(I,J,L+2)
                  ALPFRZ = ALPL + DELALP/DELZP*DZFR
                  PFRZ   = EXP(ALPFRZ)
                  PFRZL(I,J) = PFRZ
                  IF(MODELNAME == 'GFS'.OR.MODELNAME == 'RAPR')THEN
                    ES=FPVSNEW(ISOTHERM)
                    ES=MIN(ES,PFRZ)
                    QSFRZ=CON_EPS*ES/(PFRZ+CON_EPSM1*ES)
                  ELSE
                    QSFRZ=PQ0/PFRZ  & 
                     *EXP(A2*(ISOTHERM-A3)/(ISOTHERM-A4))
                  END IF
!                  QSFRZ  = PQ0/PFRZ
!     
                  RHZ      = QFRZ/QSFRZ
                  RHZ      = AMAX1(0.01,RHZ)
                  RHZ      = AMIN1(RHZ,1.0)
                  RHFRZ(I,J) = RHZ
!     
         ELSE
                  L=LICE
                  ZU      = ZMID(I,J,L)
                  ZL      = HTSFC+2.0
                  DELZ    = ZU-ZL
                  IF(TSHLTR(I,J)/=SPVAL .AND. PSHLTR(I,J)/=SPVAL)THEN
                   TSFC=TSHLTR(I,J)*(PSHLTR(I,J)*1.E-5)**CAPA
                  ELSE
! GFS analysis does not have flux file to retrieve TSFC from
                   TSFC=T(I,J,LM)+D0065*(ZMID(I,J,LM)-HTSFC-2.0)
                  END IF
                  DELT    = T(I,J,L)-TSFC
                  ZFRZ(I,J) = ZL + (ISOTHERM-TSFC)/DELT*DELZ
!     
                  DZABV   = ZFRZ(I,J)-ZL
! GFS does not output QS
                  IF(QSHLTR(I,J)/=SPVAL)THEN
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
                  PFRZL(I,J) = PFRZ
                  IF(MODELNAME == 'GFS'.OR.MODELNAME == 'RAPR')THEN
                    ES=FPVSNEW(ISOTHERM)
                    ES=MIN(ES,PFRZ)
                    QSFRZ=CON_EPS*ES/(PFRZ+CON_EPSM1*ES)
                  ELSE
                    QSFRZ=PQ0/PFRZ  &
                     *EXP(A2*(ISOTHERM-A3)/(ISOTHERM-A4))
                  END IF
!
                  RHZ     = QFRZ/QSFRZ
                  RHZ     = AMAX1(0.01,RHZ)
                  RHZ     = AMIN1(RHZ,1.0)
                  RHFRZ(I,J)= RHZ
         ENDIF
!     
!              BOUND ISOTHERM LEVEL RH.  ISOTHERM LEVEL HEIGHT IS
!              MEASURED WITH RESPECT TO MEAN SEA LEVEL.
!
               RHFRZ(I,J) = AMAX1(0.01,RHFRZ(I,J))
               RHFRZ(I,J) = AMIN1(RHFRZ(I,J),1.00)
               ZFRZ(I,J)  = AMAX1(0.0,ZFRZ(I,J))
         ELSE
               RHFRZ(I,J) = spval
               ZFRZ(I,J)  = spval
         ENDIF
 20   CONTINUE
!     
!     END OF ROUTINE.
!     
      RETURN
      END
