!> @file
!> @brief Subroutine that computes FRZING LVL, Z and RH.
!>
!> This routine computes the freezing level height and relative
!> humidity at this level for each mass point on the ETA grid.
!> The computed freezing level height is the mean sea level
!> height. At each mass point we move up from the surface to  
!> find the first ETA layer where the temperature is less than
!> 273.16K. Vertical interpolation in temperature to the freezing 
!> temperature gives the freezing level height. Pressure and   
!> specific humidity are interpolated to this level and along with
!> the temperature provide the freezing level relative humidity.
!> If the surface (skin) temperature is below freezing, the routine
!> uses surface based fields to compute the relative humidity.
!> 
!> Note that in posting freezing level data the LFM look-alike file
!> (IE, GRID 26), we pack 273.15K as the freezing temperature. All 
!> other output grids use 273.16K.
!>
!> @param[out] ZFRZ Above ground level freezing height.
!> @param[out] RHFRZ Relative humidity at freezing level.
!> @param[out] PFRZL pressure at freezing level.
!> 
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-22 | Russ Treadon | Initial
!> 1993-06-05 | Russ Treadon | Corrected freezing level heights to be with respect to mean sea level, not above ground level
!> 1998-06-15 | T Black      | Conversion from 1-D to 2-D
!> 1998-08-17 | Mike Baldwin | Compute RH over ice if necessary
!> 1998-12-22 | Mike Baldwin | Back out RH over ice
!> 2000-01-04 | Jim Tuccillo | MPI version           
!> 2001-10-25 | H Chuang     | Modified to process hybrid model output
!> 2002-01-15 | Mike Baldwin | WRF version
!> 2010-08-27 | T. Smirnova  | Added PFRZL to the output
!> 2019-10-30 | Bo Cui       | Remove "GOTO" statement
!> 2020-11-10 | Jesse Meng   | Use UPP_PHYSICS module
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE FRZLVL(ZFRZ,RHFRZ,PFRZL)

!     
!     
      use vrbls3d, only: pint, t, zmid, q, pmid
      use vrbls2d, only: fis, tshltr, pshltr, qshltr
      use masks, only: lmh
      use params_mod, only: gi, d00, capa, d0065, tfrz, pq0, a2, a3, a4
      use ctlblk_mod, only: jsta, jend, spval, lm, modelname, im
      use physcons_post, only: con_rd, con_rv, con_eps, con_epsm1
      use upp_physics, only: FPVSNEW

      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL,dimension(im,jsta:jend) :: RHFRZ, ZFRZ, PFRZL
      integer I,J,LLMH,L
      real HTSFC,PSFC,TSFC,QSFC,QSAT,RHSFC,DELZ,DELT,DELQ,DELALP,     &
           DELZP,ZL,DZABV,QFRZ,ALPL,ALPH,ALPFRZ,PFRZ,QSFRZ,RHZ,ZU,    &
           DZFR,ES
!     
!*********************************************************************
!     START FRZLVL.
!
!
!     
!     LOOP OVER HORIZONTAL GRID.
!     
!!$omp  parallel do                                                   &
!    & private(i,j,alpfrz,alph,alpl,delalp,delq,delt,delz,            &
!    &         delzp,dzabv,dzfr,htsfc,l,llmh,psfc,qfrz,               &
!    &         qsat,qsfc,qsfrz,rhsfc,rhz,tsfc,                        &
!    &         zl,zu)

       DO 20 J=JSTA,JEND
       DO 20 I=1,IM
         HTSFC    = FIS(I,J)*GI
         LLMH     = NINT(LMH(I,J))
         RHFRZ(I,J) = D00
         ZFRZ(I,J)  = HTSFC
         PSFC    = PINT(I,J,LLMH+1)
         PFRZL(I,J) = PSFC
!     
!        CHECK IF FREEZING LEVEL IS AT THE GROUND.
!     
!         IF(SM(I,J)/=SPVAL .AND. THZ0(I,J)/=SPVAL .AND.        &
!      	   THS(I,J)/=SPVAL)THEN
!          TSFC = (SM(I,J)*THZ0(I,J)+(1.-SM(I,J))*THS(I,J))     &
!      	    *(PINT(I,J,NINT(LMH(I,J))+1)/P1000)**CAPA
! Per AWC's request, use 2m T instead of skin T so that freezing level
! would be above ground more often
         IF(TSHLTR(I,J)/=SPVAL .AND. PSHLTR(I,J)/=SPVAL)THEN
          TSFC=TSHLTR(I,J)*(PSHLTR(I,J)*1.E-5)**CAPA
         ELSE
! GFS analysis does not have flux file to retrieve TSFC from	 
	  TSFC=T(I,J,LM)+D0065*(ZMID(I,J,LM)-HTSFC-2.0)
	 END IF  
         IF (TSFC<=TFRZ) THEN
!            ZFRZ(I,J) = HTSFC+(TSFC-TFRZ)/D0065
	    ZFRZ(I,J) = HTSFC+2.0+(TSFC-TFRZ)/D0065
!	    IF(SM(I,J)/=SPVAL .AND. QZ0(I,J)/=SPVAL .AND.      &
!      	      QS(I,J)/=SPVAL)THEN
!             QSFC    = SM(I,J)*QZ0(I,J)+(1.-SM(I,J))*QS(I,J)
! GFS does not output QS		   
!            ELSE IF(QSHLTR(I,J)/=SPVAL)THEN
	    IF(QSHLTR(I,J)/=SPVAL)THEN
	     PSFC=PSHLTR(I,J)
             QSFC=QSHLTR(I,J)
	    ELSE
	     QSFC=Q(I,J,LM)
	     PSFC=PMID(I,J,LM)  
            END IF  
!
            IF(MODELNAME == 'GFS' .OR. MODELNAME == 'RAPR')THEN
	     ES=FPVSNEW(TSFC)
	     ES=MIN(ES,PSFC)
	     QSAT=CON_EPS*ES/(PSFC+CON_EPSM1*ES)
	    ELSE 
             QSAT=PQ0/PSFC*EXP(A2*(TSFC-A3)/(TSFC-A4))
	    END IF 
!
            RHSFC   = QSFC/QSAT
            RHSFC   = AMAX1(0.01,RHSFC)
            RHSFC   = AMIN1(RHSFC,1.0)
            RHFRZ(I,J)= RHSFC
            PFRZL(I,J)= PSFC
            CYCLE 
         ENDIF
!     
!        OTHERWISE, LOCATE THE FREEZING LEVEL ALOFT.
!
         DO 10 L = LLMH,1,-1
            IF (T(I,J,L)<=TFRZ) THEN
               IF (L<LLMH) THEN
                  DELZ = ZMID(I,J,L)-ZMID(I,J,L+1)
                  ZL   = ZMID(I,J,L+1)
                  DELT = T(I,J,L)-T(I,J,L+1)
                  ZFRZ(I,J) = ZL + (TFRZ-T(I,J,L+1))/DELT*DELZ
!     
                  DZABV = ZFRZ(I,J)-ZL
                  DELQ  = Q(I,J,L)-Q(I,J,L+1)
                  QFRZ  = Q(I,J,L+1) + DELQ/DELZ*DZABV
                  QFRZ  = AMAX1(0.0,QFRZ)
!     
!
                  ALPL   = ALOG(PMID(I,J,L+1))
                  ALPH   = ALOG(PMID(I,J,L))
                  ALPFRZ = ALPL + (ALPH-ALPL)/DELZ*DZABV
                  PFRZ   = EXP(ALPFRZ)
                  PFRZL(I,J)  = PFRZ
		  IF(MODELNAME == 'GFS' .OR.MODELNAME == 'RAPR')THEN
	            ES=FPVSNEW(TFRZ)
	            ES=MIN(ES,PFRZ)
	            QSFRZ=CON_EPS*ES/(PFRZ+CON_EPSM1*ES)
	          ELSE 
                    QSFRZ=PQ0/PFRZ    &
                     *EXP(A2*(TFRZ-A3)/(TFRZ-A4))
                  END IF
!     
                  RHZ      = QFRZ/QSFRZ
                  RHZ      = AMAX1(0.01,RHZ)
                  RHZ      = AMIN1(RHZ,1.0)
                  RHFRZ(I,J) = RHZ
!     
               ELSE
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
                  ZFRZ(I,J) = ZL + (TFRZ-TSFC)/DELT*DELZ
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
                  ALPH    = ALOG(PMID(I,J,L))
                  ALPL    = ALOG(PSFC)
                  DELALP  = ALPH-ALPL
                  ALPFRZ  = ALPL + DELALP/DELZ*DZABV
                  PFRZ    = EXP(ALPFRZ)
!
                  PFRZL(I,J)  = PFRZ
		  IF(MODELNAME == 'GFS'.OR.MODELNAME == 'RAPR')THEN
	            ES=FPVSNEW(TFRZ)
	            ES=MIN(ES,PFRZ)
	            QSFRZ=CON_EPS*ES/(PFRZ+CON_EPSM1*ES)
	          ELSE 
                    QSFRZ=PQ0/PFRZ   &
                     *EXP(A2*(TFRZ-A3)/(TFRZ-A4))
                  END IF
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
!               RHFRZ(I,J) = AMAX1(0.01,RHFRZ(I,J))
!               RHFRZ(I,J) = AMIN1(RHFRZ(I,J),1.00)
               ZFRZ(I,J)  = AMAX1(0.0,ZFRZ(I,J))
               EXIT             
            ENDIF
 10      CONTINUE
20   CONTINUE
!     
!     END OF ROUTINE.
!     
      RETURN
      END
