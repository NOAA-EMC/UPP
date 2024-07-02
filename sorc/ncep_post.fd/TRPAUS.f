!> @file
!> @brief trpaus() computes tropopause level fields.
!> 
!> This routine computes tropopause data.  At each mass
!> point a surface up search is made for the first 
!> occurrence of a three layer mean lapse rate less than
!> or equal to a critical lapse rate.  This critcal lapse
!> rate is 2 deg/km. This is in accord with the WMO
!> definition of a tropopause. A maximum tropopause
!> pressure of 500mb is enforced. Once the tropopause
!> is located in a column, pressure, temperature, u
!> and v winds, and vertical wind shear are computed.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-22 | Russ Treadon  | Initial
!> 1997-03-06 | Geoff Manikin | Changed criteria for determining the tropopause and added height
!> 1998-06-15 | T Black       | Conversion from 1-D TO 2-D
!> 2000-01-04 | Jim Tuccillo  | MPI Version
!> 2002-04-23 | Mike Baldwin  | WRF Version
!> 2019-10-30 | Bo Cui        | ReMOVE "GOTO" STATEMENT
!> 2021-09-13 | JESSE MENG    | 2D DECOMPOSITION
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
!------------------------------------------------------------------------------
!> @brief Computes tropopause data.
!> 
!> @param[out] PTROP Tropopause pressure.
!> @param[out] TTROP Tropopause temperature.
!> @param[out] ZTROP Tropopause height.
!> @param[out] UTROP Tropopause u wind component.
!> @param[out] VTROP Tropopause v wind component.
!> @param[out] SHTROP Vertical wind shear at tropopause.
!>
      SUBROUTINE TRPAUS(PTROP,TTROP,ZTROP,UTROP,VTROP,SHTROP)

!     
!     
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE PARAMETERS.
!
       use vrbls3d,    only: pint, t, zint, uh, vh
       use masks,      only: lmh
       use params_mod, only: d50
       use ctlblk_mod, only: jsta, jend, spval, im, jm, lm, &
                             ista, iend
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
!     PARAMTER CRTLAP SPECIFIES THE CRITICAL LAPSE RATE
!     (IN K/M) IDENTIFYING THE TROPOPAUSE.  WE START 
!     LOOKING FOR THE TROPOPAUSE ABOVE PRESSURE LEVEL
!     PSTART (IN PASALS).
      real,PARAMETER :: CRTLAP=0.002E0, PSTART=5.0E4
!     
!     DECLARE VARIABLES.
!     
      REAL,dimension(IM,JM),intent(out) :: PTROP,TTROP,ZTROP,UTROP,  &
           VTROP,SHTROP
      REAL TLAPSE(LM),DZ2(LM),DELT2(LM),TLAPSE2(LM)
!
      integer I,J,LL,LLMH,L
      real PM,DELT,DZ,RSQDIF
!     
!*****************************************************************************
!     START TRPAUS HERE.
!     
!     LOOP OVER THE HORIZONTAL GRID.
!    
      DO J=JSTA,JEND
      DO I=ISTA,IEND
         PTROP(I,J)  = SPVAL
         TTROP(I,J)  = SPVAL
         ZTROP(I,J)  = SPVAL
         UTROP(I,J)  = SPVAL
         VTROP(I,J)  = SPVAL
         SHTROP(I,J) = SPVAL
      ENDDO
      ENDDO
!
!!$omp  parallel do
!!$omp& private(delt,delt2,dz,dz2,ie,iw,l,llmh,pm,rsqdif,
!!$omp&         tlapse,tlapse2,u0,u0l,uh,uh0,ul,
!!$omp&         v0,v0l,vh,vh0)
       DO J=JSTA,JEND
       loopI:DO I=ISTA,IEND
!     
!        COMPUTE THE TEMPERATURE LAPSE RATE (-DT/DZ) BETWEEN ETA 
!        LAYERS MOVING UP FROM THE GROUND.  THE FIRST ETA LAYER
!        ABOVE PRESSURE "PSTART" IN WHICH THE LAPSE RATE IS LESS
!        THAN THE CRITCAL LAPSE RATE IS LABELED THE TROPOPAUSE.
!
        LLMH=NINT(LMH(I,J))
!
        loopL: DO L=LLMH-1,2,-1
        PM     = PINT(I,J,L)
        DELT   = T(I,J,L-1)-T(I,J,L)
        DZ     = D50*(ZINT(I,J,L-1)-ZINT(I,J,L+1))
        TLAPSE(L) = -DELT/DZ
!
        IF ((TLAPSE(L)<CRTLAP).AND.(PM<PSTART)) THEN 
          IF (L == 2 .AND. TLAPSE(L) < CRTLAP) GOTO 15
          DZ2(L+1) = 0.
!
          DO 17 LL=L,3,-1
          DZ2(LL) = 0.
          DELT2(LL) = 0.
          TLAPSE2(LL) = 0.
          DZ2(LL) = (2./3.)*(ZINT(I,J,LL-2)-ZINT(I,J,L+1))
          IF ((DZ2(LL) > 2000.) .AND.                    &
              (DZ2(LL+1) > 2000.)) GO TO 15
          DELT2(LL) = T(I,J,LL-2)-T(I,J,L)
          TLAPSE2(LL) = -DELT2(LL)/DZ2(LL)
!
          IF (TLAPSE2(LL) > CRTLAP) THEN
            CYCLE loopL
          ENDIF
!
   17     CONTINUE 
        ELSE
          CYCLE loopL
        ENDIF 
!
   15   PTROP(I,J)  = D50*(PINT(I,J,L)+PINT(I,J,L+1))
        TTROP(I,J)  = T(I,J,L)
        ZTROP(I,J)= 0.5*(ZINT(I,J,L)+ZINT(I,J,L+1))
!
        UTROP (I,J) = UH(I,J,L)
        VTROP (I,J) = VH(I,J,L)
        DZ        = ZINT(I,J,L)-ZINT(I,J,L+1)
        RSQDIF    = SQRT(((UH(I,J,L-1)-UH(I,J,L+1))*0.5)**2  &
     &                  +((VH(I,J,L-1)-VH(I,J,L+1))*0.5)**2)
        SHTROP(I,J) = RSQDIF/DZ
        CYCLE loopI

        ENDDO loopL

!X         WRITE(88,*)'REACHED TOP FOR K,P,TLAPSE:  ',K,PM,TLAPSE

        DZ       = D50*(ZINT(I,J,2)-ZINT(I,J,3))
        PTROP(I,J) = D50*(PINT(I,J,2)+PINT(I,J,3))
        TTROP(I,J) = T(I,J,2)
        ZTROP(I,J)= D50*(ZINT(I,J,2)+ZINT(I,J,3))
        UTROP (I,J) = UH(I,J,2)
        VTROP (I,J) = VH(I,J,2)
        RSQDIF    = SQRT(((UH(I,J,1)-UH(I,J,3))*0.5)**2    &
     &                  +((VH(I,J,1)-VH(I,J,3))*0.5)**2)
        SHTROP(I,J) = RSQDIF/DZ

!X        WRITE(82,1010)I,J,L,PTROP(I,J)*D01,TTROP(I,J),
!X     X       UTROP(I,J),VTROP(I,J),SHTROP(I,J)
!     
      ENDDO loopI !end I
      ENDDO !end J

!     
!     END OF ROUTINE.
!     
      RETURN
      END
