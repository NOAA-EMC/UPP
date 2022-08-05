      SUBROUTINE MIXLEN(EL0,EL)
!
!     CALCULATES LAYER-AVERAGED BLACKADAR'S MIXING LENGTH, AND PBL TOP
!     AS CPBLT*(ASYMPTOTIC EL); AND THEN EL, ACCOUNT TAKEN OF STABILITY,
!     PBL TOP AND VERTICAL GRID DISTANCE RESTRICTIONS (SEE BELOW)
!
!     SET FROM EXISTING CODES BY L. LOBOCKI, JUNE 5, 1992
!       MODIFIED BY FEDOR MESINGER, OCTOBER 13, NOVEMBER 19
!       MODIFIED BY JIM TUCCILLO FOR MPI IMPLEMENTATION
!   01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-06-19  MIKE BALDWIN - WRF VERSION
!   21-03-11  B Cui - change local arrays to dimension (im,jsta:jend)
!   21-09-30  J MENG - 2D DECOMPOSITION
!   
!
!     INPUT:
!     ------
!
!     ZINT (IM,jsta_2l:jend_2u,LP1) - ETA INTERFACES HEIGHT FIELD
!     T    (IM,jsta_2l:jend_2u,LM)  - TEMPERATURE
!     PMID (IM,jsta_2l:jend_2u,LM)  - PRESSURE IN LAYERS
!     Q2   (IM,jsta_2l:jend_2u,LM)  - TURBULENCE KINETIC ENERGY * 2
!     HGT  (IM,jsta_2l:jend_2u)     - SURFACE ELEVATION ARRAY
!     HTM  (IM,jsta_2l:jend_2u,LM)  - HEIGHT TOPOGRAPHY MASK ARRAY
!     EL0  (IM,JM)     - ARRAY OF ASYMPTOTIC VALUES FOR MIXING LENGTH
!
!     OUTPUT:
!     -------
!
!     EL   (IM,jsta_2l:jend_2u,LM) - FIELD OF RESULTING MASTER LENGTH SCALES
!
!
!     SCRATCH AREAS:
!     --------------
!
!     VKRMZ(IM,JM)
!
!     RELEVANT CONSTANTS:
!     -------------------
!
!     VON KARMAN CONSTANT:
      use vrbls3d, only: zint, pmid, t, q2
      use masks, only: lmh, htm
      use params_mod, only: EPSQ2, CAPA
      use ctlblk_mod, only: jsta, jend, jsta_m, jend_m, im, jm, jsta_2l, jend_2u,&
              lm, lm1, spval,&
              ista, iend, ista_m, iend_m, ista_2l, iend_2u

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      real,PARAMETER :: VKRM=0.4
!     CONSTANTS NEEDED FOR THE EL(BL,ST,ZI) SCHEME:
      real,PARAMETER :: FRG=4.*9.8,DRDRFF=0.54,CPBLT=10.,     &
        CSH=0.23*0.5, EPSN2=1.E-7
!
!     ------------------------------------------------------------------
!
      real,intent(in) :: el0(ista_2l:iend_2u,jsta_2l:jend_2u)
      real,intent(out) ::  EL(ista_2l:iend_2u,jsta_2l:jend_2u,LM)
      real HGT(ISTA:IEND,JSTA:JEND),APE(ISTA_M:IEND_M,JSTA_M:JEND_M,2)
!
      integer I,J,L
      real ZL,VKRMZ,ENSQ,Q2KL,ELST,ZIAG,ELVGD

!***********************************************************************
!
!$omp  parallel do
      DO L=1,LM
        DO J=JSTA,JEND
        DO I=ISTA,IEND
          EL(I,J,L)=0.
        ENDDO
        ENDDO
      ENDDO
        DO J=JSTA,JEND
        DO I=ISTA,IEND
          HGT(I,J)=ZINT(I,J,NINT(LMH(I,J))+1)
        ENDDO
        ENDDO
!
!---THE AVERAGE EL SCHEME---------------------------(FM, AUGUST 19 MEMO)
!   FIRST GET EL IN THE LAYERS
!
!$omp  parallel do private(i,j,l,vkrmz,zl)
      DO L=1,LM
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF(HGT(I,J)<spval)THEN
            ZL        = 0.5*(ZINT(I,J,L)+ZINT(I,J,L+1))
            VKRMZ     = (ZL-HGT(I,J))*VKRM
            EL(I,J,L) = EL0(I,J)*VKRMZ/(EL0(I,J)+VKRMZ)
            ELSE
            EL(I,J,L) = spval
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!***
!***  GET NOW THE INTERFACE EL BY TWO-POINT AVERAGING OF LAYER VALUES
!***
      DO L=1,LM1
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF(HGT(I,J)<spval)THEN
            EL(I,J,L) = 0.5*(EL(I,J,L)+EL(I,J,L+1))*HTM(I,J,L+1)
            ELSE
            EL(I,J,L) = spval
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IF(HGT(I,J)<spval)THEN
          EL(I,J,LM) = 0.0
          ELSE
          EL(I,J,LM) = spval
          ENDIF
        ENDDO
      ENDDO
!---STABILITY, PBL TOP, AND VERTICAL GRID DISTANCE RESTRICTIONS:--------
!   COMPUTE EL STABLE AND
!   * USE THE SMALLER OF EL BLACKADAR, EL STABLE IF WITHIN PBL;
!   * USE THE SMALLEST OF EL STABLE, ELVGD, AND VKRMZ IF ABOVE PBL
!       (ASSUME PBL TOP IS AT CPBLT*EL0(K));
!$omp  parallel do private(i,j)
      DO J=JSTA_M,JEND_M
        DO I=ISTA_M,IEND_M
          APE(I,J,1) = (1.E5/PMID(I,J,1))**CAPA
        ENDDO
      ENDDO
!
      DO L=1,LM1
!$omp  parallel do private(i,j,elst,elvgd,ensq,q2kl,ziag)
        DO J=JSTA_M,JEND_M
          DO I=ISTA_M,IEND_M
            IF(T(I,J,L)<spval)THEN
            APE(I,J,2) = (1.E5/PMID(I,J,L+1))**CAPA
            ENSQ = HTM(I,J,L+1)*                                     &
                   FRG*(T(I,J,L)*APE(I,J,1)-T(I,J,L+1)*APE(I,J,2))/  &
                  ((T(I,J,L)*APE(I,J,1)+T(I,J,L+1)*APE(I,J,2))*     &
                   (ZINT(I,J,L)-ZINT(I,J,L+2))+EPSN2)
            ENSQ = AMAX1(ENSQ,EPSN2)
            Q2KL = AMAX1(EPSQ2,Q2(I,J,L))
            ELST = DRDRFF*SQRT(Q2KL/ENSQ)
!WAS        ELST = DRDRFF*SQRT(Q2(I,J,L)/ENSQ)
            ZIAG = ZINT(I,J,L+1)-HGT(I,J)
!
            IF(ZIAG < CPBLT*EL0(I,J))THEN
              EL(I,J,L) = AMIN1(EL(I,J,L),ELST)
            ELSE
              ELVGD     = CSH*(ZINT(I,J,L)-ZINT(I,J,L+2))
              EL(I,J,L) = AMIN1(ELST,ELVGD,VKRM*ZIAG)
            ENDIF
            APE(I,J,1) = APE(I,J,2)
            ELSE
            EL(I,J,L) = spval
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!
      RETURN
      END

