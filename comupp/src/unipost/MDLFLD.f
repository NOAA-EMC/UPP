      SUBROUTINE MDLFLD
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    MDLFLD      SLP AND NATIVE LEVEL POSTING
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-21       
!     
! ABSTRACT:
!     THIS ROUTINE DOES SEVERAL THINGS.  IT IS THE FIRST 
!     ROUTINE CALLED BY POST PROCESSOR SUBROUTINE PROCESS 
!     WHICH SETS THE ORDER IN WHICH FIELDS ARE POSTED.  THE
!     NEGATIVE SPECIFIC HUMIDITY IS CLIPPED.
!     COMPUTE THE STANDARD NMC SEA LEVEL PRESSURE IF THIS OPTION
!     IS ACTIVATED.  FINALLY WE COMPUTE/POST REQUESTED FIELDS ON
!     MODEL LAYERS.
!
!   .     
!     
! PROGRAM HISTORY LOG:
!   92-12-21  RUSS TREADON
!   93-09-01  RUSS TREADON - ADDED ADDITIONAL OUTPUT FIELDS.
!   96-03-20  MIKE BALDWIN - ADDED CLOUD TOP TEMPS, CHANGE CLOUD WATER
!                            TO CONTAIN WATER ONLY
!   97-04-29  GEOFF MANIKIN - MOVED CLOUD TOP TEMPS TO CLDRAD
!   98-06-01  T BLACK - CONVERSION FROM 1-D TO 2-D
!   98-07-20  MIKE BALDWIN - REMOVED LABL84
!   98-08-18  T BLACK - REMOVED EXCESS SPACE IN EXTRA.com
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   01-10-22  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-01-15  MIKE BALDWIN - WRF VERSION
!   04-11-17  H CHUANG, B FERRIER, AND Y JIN - ADD HYDROMETEORS, 
!					VISIBILITY & RADAR REFLECTIVITY
!   05-07-07  B ZHOU ADD RSM MODEL A GRID     
!   05-08-18  B ZHOU ADD /VISB/ COMMON BLOCK TO PASS VISIBILITY TO
!                        AVIATION SUBROUTINE TO CALCULATE FLIGHT
!                        CONDITION RESTRICTION
!
! USAGE:    CALL MDLFLD
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
!       BOUND    - BOUND ARRAY ELEMENTS BETWEEN LOWER AND UPPER LIMITS.
!       SCLFLD   - SCALE ARRAY ELEMENTS BY SCALAR CONSTANT.
!       NGMSLP   - COMPUTE SLP USING STANDARD NMC REDUCTION METHOD.
!       CALPOT   - COMPUTE POTENTIAL TEMPERATURE.
!       CALRH    - COMPUTE RELATIVE HUMIDITY.
!       CALDWP   - COMPUTE DEWPOINT TEMPERATURE.
!       CALMCVG  - COMPUTE MOISTURE CONVERGENCE.
!       CALVOR   - COMPUTE ABSOLUTE VORTICITY.
!       CALSTRM  - COMPUTE GEOSTROPHIC STREAMFUNCTION.
!       CALMICT  - COMPUTES NEW CLOUD FIELDS AND RADAR REFLECTIVITY FACTOR
!     LIBRARY:
!       COMMON   - 
!                  RQSTFLD
!                  CTLBLK
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
!
      use vrbls3d
      use vrbls2d
      use masks
      use params_mod
      use pmicrph_mod
!      use paramr_mod
      use ctlblk_mod
      use rqstfld_mod
      use gridspec_mod
!     
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      REAL, PARAMETER :: CURATE=24.*1000., CTIM1=0., CTIM2=24.*3600.    &
     &, RAINCON=0.8333*1.1787E4, SNOCON=0.94*1.4594E5                   &
! specify in params now
!
!--- 88D reflectivity algorithm, Z = 300.*R**1.4 , R is rain rate in mm/h
!
     &, DBZmax=80., ZR_A=300., ZR_B=1.4
!
!--- Modification of Slingo (1987) to enhance convective cloudiness
!
      REAL CC(10), PPT(10)
      DATA CC / 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 /
      DATA PPT/  0., .14, .31, .70, 1.6, 3.4, 7.7, 17., 38., 85. /
      INTEGER ICBOT(IM,JM),ICTOP(IM,JM) 

!     
!     DECLARE VARIABLES.
!     
      LOGICAL NORTH,NEED(IM,JM)
      REAL EGRID1(IM,JM),EGRID2(IM,JM),EGRID3(IM,JM),EGRID4(IM,JM)  &
     &,     EL0(IM,JM)   &
     &,     P1D(IM,JM),T1D(IM,JM),Q1D(IM,JM),EGRID5(IM,JM)          &
     &,     C1D(IM,JM),FI1D(IM,JM),FR1D(IM,JM),FS1D(IM,JM)          &
     &,     QW1(IM,JM),QI1(IM,JM),QR1(IM,JM),QS1(IM,JM)    &
     &,     GRID1(IM,JM), GRID2(IM,JM)    &
     &,     CUREFL_S(IM,JM), CUREFL(IM,JM), CUREFL_I(IM,JM)    &
     &,     Zfrz(IM,JM), DBZ1(IM,JM),DBZR1(IM,JM),DBZI1(IM,JM)    &
     &,     DBZC1(IM,JM)
!
      REAL, ALLOCATABLE :: EL(:,:,:),RICHNO(:,:,:),PBLRI(:,:)     &
     &   ,PBLREGIME(:,:)      
!
      REAL QI(IM,JM),QINT(IM,JM)
      REAL TT(IM,JM),PPP(IM,JM),QV(IM,JM),QCD(IM,JM),QICE1(IM,JM)
      REAL QRAIN1(IM,JM),QSNO1(IM,JM), refl(im,jm),QG1(IM,JM)
      REAL RH(IM,JM)
      integer I,J,L,Lctop,LLMH,IICE,LL,II,JJ,IFINCR,ITHEAT,NC,NMOD
      real RDTPHS,CFRdum,PMOD,CC1,CC2,P1,P2,CUPRATE,FACR,RRNUM,   &
           RAINRATE,TERM1,TERM2,TERM3,QROLD,SNORATE,&
           DENS,DELZ,FCTR

      REAL rain,ronv,slor,snow,rhoqs,temp_c,sonv,slos             &
       ,graupel,rhoqg,gonv,slog
      real alpha, rhod, bb
      real ze_s, ze_r, ze_g, ze_max, ze_nc, ze_conv, ze_sum

      real ze_smax, ze_rmax,ze_gmax
      real ze_nc_1km, ze_nc_4km, dz

      REAL T700(IM,JM),TH700(IM,JM),SDUMMY(IM,2)
      REAL PSFC,TSFC,ZSFC,ZSL,TAUCR,GORD,DP,CONST               
      integer iz1km,iz4km

      real LAPSES, EXPo,EXPINV,TSFCNEW
      REAL GAM,GAMD,GAMS

      PARAMETER (ZSL=0.0)
      PARAMETER (TAUCR=RD*GI*290.66,CONST=0.005*G/RD)
      PARAMETER (GORD=G/RD,DP=60.E2)

        GAMS = 0.0065
        GAMD = 0.0100

        LAPSES = 0.0065
! deg K / meter
        EXPo = ROG*LAPSES
        EXPINV = 1./EXPo

!      REAL VIS(IM,JM)
!     

! ADD by B Zhou
!      COMMON /VISB/VIS
!
!     
!*****************************************************************************
!     START SUBROUTINE MDLFLD.
!
!     ALLOCATE LOCAL ARRAYS
!
      ALLOCATE(EL     (IM,JSTA_2L:JEND_2U,LM))     
      ALLOCATE(RICHNO (IM,JSTA_2L:JEND_2U,LM))     
      ALLOCATE(PBLRI  (IM,JSTA_2L:JEND_2U))     
!HC COMMENT OUT THE CALL BECAUSE MEMBRANE SLP REDUCTION IS
!HC NOW DONE IN MDL2P
!
!     OUTPUT SEA LEVEL PRESSURE IF REQUESTED.
!     FIRST, MESINGER'S SEA LEVEL PRESSURE.
!HC      IF (IGET(023).GT.0) THEN
!HC         DO J=JSTA,JEND
!HC         DO I=1,IM
!HC           GRID1(I,J)=PSLP(I,J)
!HC         ENDDO
!HC         ENDDO
!HC         ID(1:25) = 0
!HC         CALL GRIBIT(IGET(023),LVLS(1,IGET(023)),
!HC     X        GRID1,IM,JM)
!HC      ENDIF
!     
!     SECOND, STANDARD NGM SEA LEVEL PRESSURE.
      IF (IGET(105).GT.0) THEN
         CALL NGMSLP
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J)=SLP(I,J)
         ENDDO
         ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(105),LVLS(1,IGET(105)),GRID1,IM,JM)
      ENDIF
!
!--- Calculate convective cloud fractions following radiation in
!    NMM; used in subroutine CALRAD_WCLOUD for satellite radiances
!
      IF (MODELNAME=='NMM' .OR. imp_physics==5) THEN
        print*,'DTQ2 in MDLFLD= ',DTQ2
        RDTPHS=24.*3.6E6/DTQ2
        DO J=JSTA,JEND
          DO I=1,IM
          IF ((HBOT(I,J)-HTOP(I,J)) .LE. 1.0) THEN
            ICBOT(I,J)=0
            ICTOP(I,J)=0
            CNVCFR(I,J)=0.
          ELSE
            ICBOT(I,J)=NINT(HBOT(I,J))
            ICTOP(I,J)=NINT(HTOP(I,J))
            CFRdum=CC(1)
            PMOD=RDTPHS*CPRATE(I,J)       ! mm/day
            IF (PMOD .GT. PPT(1)) THEN
              DO NC=1,10
                IF(PMOD.GT.PPT(NC)) NMOD=NC
              ENDDO
              IF (NMOD .GE. 10) THEN
                CFRdum=CC(10)
              ELSE
                CC1=CC(NMOD)
                CC2=CC(NMOD+1)
                P1=PPT(NMOD)
                P2=PPT(NMOD+1)
                CFRdum=CC1+(CC2-CC1)*(PMOD-P1)/(P2-P1)
              ENDIF   !--- End IF (NMOD .GE. 10) ...
              CFRdum=MIN(H1, CFRdum)
            ENDIF     !--- End IF (PMOD .GT. PPT(1)) ...
!            CNVCFR(I,J)=100.*CFRdum
            CNVCFR(I,J)=CFRdum
          ENDIF       !--- End IF (HBOT(I,J)-HTOP(I,J) .LE. 1.0) ...
          ENDDO       !--- DO I=1,IM
        ENDDO         !--- DO J=JSTA,JEND
      ENDIF           !-- IF (MODELNAME=='NMM' .OR. imp_physics==5) THEN
!
!    Calculate convective radar reflectivity at the surface (CUREFL_S), 
!    and the decrease in reflectivity above the 0C level (CUREFL_I)
!
      IF(imp_physics.eq.5)THEN
       RDTPHS=3.6E6/DTQ2
       DO J=JSTA,JEND
        DO I=1,IM
          CUPRATE=RDTPHS*CPRATE(I,J)            !--- Cu precip rate, R (mm/h)
!          CUPRATE=CUPPT(I,J)*1000./TRDLW        !--- mm/h
          Zfrz(I,J)=ZMID(I,J,NINT(LMH(I,J)))  !-- Initialize to lowest model level
          DO L=1,NINT(LMH(I,J))               !-- Start from the top, work down
             IF (T(I,J,L) .GE. TFRZ) THEN
                Zfrz(I,J)=ZMID(I,J,L)         !-- Find highest level where T>0C
                EXIT
             ENDIF
          ENDDO       !--- DO L=1,NINT(LMH(I,J))
!          IF (CUPRATE .LE. 0. .OR. CUPPT(I,J).LE.0.) THEN
          IF (CUPRATE .LE. 0.) THEN ! bug fix, post doesn not use CUPPT 
             CUREFL_S(I,J)=0.
             CUREFL_I(I,J)=0.
          ELSE
             CUREFL_S(I,J)=ZR_A*CUPRATE**ZR_B   !--- Use Z=A*R**B
             Lctop=NINT(HTOP(I,J))              !--- Cu cld top level
  !
  !--- Assume convective reflectivity (Z, not dBZ) above 0C level decreases
  !    with height by two orders of magnitude (20 dBZ) from the 0C level up
  !    to cloud top.  If cloud top temperature is above 0C, assume 20 dBZ
  !    decrease occurs in the first 1 km above the 0C level.
  !
             CUREFL_I(I,J)=-2./MAX( 1000., ZMID(I,J,Lctop)-Zfrz(I,J) )
          ENDIF       !--- IF (CUPRATE .LE. 0. .OR. CUPPT(I,J).LE.0.) THEN
        ENDDO         !--- End DO I
       ENDDO    

!
!--- Calculate each hydrometeor category & GRID-SCALE cloud fraction
!    (Jin, Aug-Oct '01; Ferrier, Feb '02)
!
       DO L=1,LM
        DO J=JSTA,JEND
        DO I=1,IM
          P1D(I,J)=PMID(I,J,L)
          T1D(I,J)=T(I,J,L)
          Q1D(I,J)=Q(I,J,L)
          C1D(I,J)=CWM(I,J,L)
          FI1D(I,J)=F_ice(I,J,L)
          FR1D(I,J)=F_rain(I,J,L)
          FS1D(I,J)=MAX(H1, F_RimeF(I,J,L))
    !
    !--- Estimate radar reflectivity factor at level L
    !
          CUREFL(I,J)=0.
          IF (CUREFL_S(I,J) .GT. 0.) THEN
             FCTR=0.
             LLMH = NINT(LMH(I,J)) 
             Lctop=NINT(HTOP(I,J))              !--- Cu cld top level
             IF (L.GE.Lctop .AND. L.LE.LLMH) THEN
                DELZ=ZMID(I,J,L)-Zfrz(I,J)
                IF (DELZ .LE. 0.) THEN
                   FCTR=1.        !-- Below the highest freezing level
                ELSE
       !
       !--- Reduce convective radar reflectivity above freezing level
       !
                   FCTR=10.**(CUREFL_I(I,J)*DELZ)
                ENDIF             !-- End IF (DELZ .LE. 0.)
             ENDIF                !-- End IF (L.GE.HTOP(I,J) .OR. L.LE.LLMH)
             CUREFL(I,J)=FCTR*CUREFL_S(I,J)
          ENDIF                   !-- End IF (CUREFL_S(I,J) .GT. 0.)

        ENDDO         !-- End DO I loop
        ENDDO         !-- End DO J loop 
  !
  !--- Determine composition of condensate in terms of cloud water,
  !    rain, and ice (cloud ice & precipitation ice) following
  !    GSMDRIVE in the model; composition of cloud ice & precipitation
  !    ice (snow) follows algorithm in GSMCOLUMN; radar reflectivity
  !    is derived to be consistent with microphysical assumptions 
  !
        CALL CALMICT(P1D,T1D,Q1D,C1D,FI1D,FR1D,FS1D,CUREFL        &
     &               ,QW1,QI1,QR1,QS1,DBZ1,DBZR1,DBZI1,DBZC1)
        DO J=JSTA,JEND
        DO I=1,IM
          LLMH = NINT(LMH(I,J))
          IF (L .GT. LLMH) THEN
            QQW(I,J,L)=D00
            QQI(I,J,L)=D00
            QQR(I,J,L)=D00
            QQS(I,J,L)=D00
            CFR(I,J,L)=D00
            DBZ(I,J,L)=DBZmin
            DBZR(I,J,L)=DBZmin
            DBZI(I,J,L)=DBZmin
            DBZC(I,J,L)=DBZmin
          ELSE
            QQW(I,J,L)=MAX(D00, QW1(I,J))
            QQI(I,J,L)=MAX(D00, QI1(I,J))
            QQR(I,J,L)=MAX(D00, QR1(I,J))
            QQS(I,J,L)=MAX(D00, QS1(I,J))
            DBZ(I,J,L)=MAX(DBZmin, DBZ1(I,J))
            DBZR(I,J,L)=MAX(DBZmin, DBZR1(I,J))
            DBZI(I,J,L)=MAX(DBZmin, DBZI1(I,J))
            DBZC(I,J,L)=MAX(DBZmin, DBZC1(I,J))
          ENDIF       !-- End IF (L .GT. LMH(I,J)) ...
        ENDDO         !-- End DO I loop
        ENDDO         !-- End DO J loop
                                        
       ENDDO           !-- End DO L loop        

      ELSE
! compute radar reflectivity for non-ferrier's scheme      
        print*,'calculating radar ref for non-Ferrier scheme' 
! Determine IICE FLAG
        IF(IMP_PHYSICS.EQ.1 .OR. IMP_PHYSICS.EQ.3)THEN
          IICE=0
        ELSE
          IICE=1
        END IF
        PRINT*,'IICE= ',IICE

        IF(IMP_PHYSICS.NE.8) THEN
!tgs - non-Thompson schemes

        DO L=1,LM
         DO J=JSTA,JEND
          DO I=1,IM
            IF(T(I,J,L) .LT. 1.0E-3)print*,'ZERO T'    
            IF(T(I,J,L) .gt. 1.0E-3)                            &
     &       DENS=PMID(I,J,L)/(RD*T(I,J,L)*(Q(I,J,L)*D608+1.0))      ! DENSITY
! PATCH to set QQR, QQS, AND QQG to zeros if they are negative so that post won't abort
            IF(QQR(I,J,L).LT. 0.0)QQR(I,J,L)=0.0
            IF(QQS(I,J,L).LT. 0.0)QQS(I,J,L)=0.0    ! jkw
            IF (IICE.EQ.0) THEN
               IF (T(I,J,L) .GE. TFRZ) THEN
                  DBZ(I,J,L)=((QQR(I,J,L)*DENS)**1.75)*         &
     &               3.630803E-9 * 1.E18                  ! Z FOR RAIN
                  DBZR(I,J,L)=DBZ(I,J,L)
               ELSE
!mptest            DBZ(I,J,L)=((QQR(I,J,L)*DENS)**1.75)*  &
                  DBZ(I,J,L)=((QQS(I,J,L)*DENS)**1.75)*         &
     &               2.18500E-10 * 1.E18                  ! Z FOR SNOW
                  DBZI(I,J,L)=DBZ(I,J,L)
               ENDIF
            ELSEIF (IICE.EQ.1) THEN
	       IF(QQS(I,J,L).LT. 0.0)QQS(I,J,L)=0.0
	       IF(QQG(I,J,L).LT. 0.0)QQG(I,J,L)=0.0
               DBZR(I,J,L)=((QQR(I,J,L)*DENS)**1.75)*           &
     &               3.630803E-9 * 1.E18                  ! Z FOR RAIN
               DBZI(I,J,L)= DBZI(I,J,L)+((QQS(I,J,L)*DENS)**1.75)* &
     &               2.18500E-10 * 1.E18                  ! Z FOR SNOW
               IF (QQG(I,J,L) < SPVAL) &
                 DBZI(I,J,L)= DBZI(I,J,L)+((QQG(I,J,L)*DENS)**1.75)* &
     &               1.033267E-9 * 1.E18                  ! Z FOR GRAUP
               DBZ(I,J,L)=DBZR(I,J,L)+DBZI(I,J,L)
!               IF(L.EQ.27.and.QQR(I,J,L).gt.1.e-4)print*,
!     &'sample QQR DEN,DBZ= ',QQR(I,J,L),DENS,DBZ(I,J,L)
            ENDIF
            IF (DBZ(I,J,L).GT.0.) DBZ(I,J,L)=10.0*LOG10(DBZ(I,J,L)) ! DBZ
            IF (DBZR(I,J,L).GT.0.)DBZR(I,J,L)=10.0*LOG10(DBZR(I,J,L)) ! DBZ
            IF (DBZI(I,J,L).GT.0.)      &
     &         DBZI(I,J,L)=10.0*LOG10(DBZI(I,J,L)) ! DBZ
            LLMH = NINT(LMH(I,J))
            IF(L.GT.LLMH)THEN
             DBZ(I,J,L)=DBZmin
             DBZR(I,J,L)=DBZmin
             DBZI(I,J,L)=DBZmin
            ELSE
             DBZ(I,J,L)=MAX(DBZmin, DBZ(I,J,L))
             DBZR(I,J,L)=MAX(DBZmin, DBZR(I,J,L))
             DBZI(I,J,L)=MAX(DBZmin, DBZI(I,J,L))
            END IF
           ENDDO
          ENDDO
         ENDDO
!tgs
        ELSE
! for Thompson microphisics scheme (option 8), developed at GSD/ESRL
! 13 January 2009
      call paramr       ! compute constants for reflectivity algorithm

      bb = 0.           !  bright band effect - yes or no (0)

      alpha = 0.224 ! = (1000kg/m^3/917kg/m^3)**2)*(0.176/0.930)
!                      1000kg/m^3 is density of liquid water
!                       917kg/m^3 is density of solid ice
!                      0.176 = dielectric factor of ice
!                      0.930 = dielectric factor of liquid water

      ze_smax = -1.E30
      ze_rmax = -1.E30
      ze_gmax = -1.E30

         DO J=JSTA,JEND
          DO I=1,IM
        refl(i,j) = -10.
        ze_max = -10.

          iz1km = 0
          iz4km = 0

        DO L=1,LM
          LL=LM-L+1 
            IF(T(I,J,LL) .LT. 1.0E-3)print*,'ZERO T'    
            IF(T(I,J,LL) .gt. 1.0E-3)                            &
             RHOD=PMID(I,J,LL)/                                  & 
               (RD*T(I,J,LL)*(Q(I,J,LL)*D608+1.0))      ! DENSITY
             DZ=ZINT(i,j,ll)-ZINT(i,j,lm)
!      Particle size distributions and reflectivity
!      ---------------------------------------------
!       Much of this code borrowed from EXMOISG loop 20 to get particle size 
!       distributions 

!jmb--    Note that SLOR, SLOS and SLOG are inverse slopes!!!! Also,
!          RONV,SONV,GONV, M-P zero intercept values, normalized by
!          max allowable values. 

!         Have to set min values of hydrometeors (r1) large enough to avoid
!          underflow problems with log later on.  

!   -- rain
              ze_r = 1.e-35
              if (qqr(i,j,ll).lt.1.e-6) go to 124

              rain = max(r1,qqr(i,j,ll))
              ronv = (const1r*tanh((qr0 - rain)/delqr0) +        &
               const2r)/ron
              SLOR=(RHOd*RAIN/(TOPR*RONV))**0.25
              ze_r = 720.*ronv*ron*slor**7 ! Stoelinga Eq. 2, reflectivity

124         continue

!   -- snow
              ze_s = 1.e-35
              if (qqs(i,j,ll).lt.1.e-6) go to 125
              snow = max(r1,qqs(i,j,ll))
!             New SONV formulation based on Fig. 7, curve_3 of Houze et al 1979
              rhoqs=RHOd*snow
              temp_C = min(-0.001, T(i,j,ll)-273.15)
              sonv = (min(2.0E8, 2.0E6*exp(-0.12*temp_C)))/son
              slos=(rhoqs/(tops*sonv))**0.25
              ze_s = 720.*alpha*sonv*son*slos**7*(dsnow/drain)**2
!               From Stoelinga Eq. 5, reflectivity

!             For bright band, increase reflectivity by factor of 5.28,
!              which is ratio of dielectric factors for water/ice (.930/.176)
              IF (T(i,j,ll) .gt. 273.15)                         &
               ze_s = ze_s*(1. + 4.28*bb)

125         continue

!   -- graupel
              ze_g = 1.e-35
              if (qqg(i,j,ll).lt.1.e-6) go to 126
              graupel = max(r1,qqg(i,j,ll))
              rhoqg=RHOd*graupel
              gonv=1.
              gonv=const_ng1*(rhoqg**const_ng2)
              gonv = max(1.e4, min(gonv,gon))
              gonv = gonv/gon
              slog=(rhoqg/(topg*gonv))**0.25
              ze_g = 720.*alpha*gonv*gon*slog**7*(dgraupel/drain)**2
!               Stoelinga Eq. 5 applied to graupel

!             For bright band
              IF (t(i,j,ll) .gt. 273.15)                         &
               ze_g = ze_g*(1. + 4.28*bb)

126         continue

!   -- total grid scale
              ze_nc = ze_r + ze_s + ze_g

              if (iz1km.eq.0 .and. dz.gt.1000.) then
                 ze_nc_1km = ze_nc
                 iz1km = 1
              end if

              if (iz4km.eq.0 .and. dz.gt.4000.) then
                 ze_nc_4km = ze_nc
                 iz4km = 1
              end if

              ze_rmax = max(ze_r,ze_rmax)
              ze_smax = max(ze_s,ze_smax)
              ze_gmax = max(ze_g,ze_gmax)


!           Reflectivities are in units of m^6/m^3
!            convert to mm^6/m^3 and take log base 10 to get 
!            reflectivities in dbZe (decibels).
!            comp_refl_r(j,k) = 10.*LOG10(ze_r*1.E18)
!            comp_refl_s(j,k) = 10.*LOG10(ze_s*1.E18)
!            comp_refl_g(j,k) = 10.*LOG10(ze_g*1.E18)
!           comp_refl_nc(j,k) = 10.*LOG10(ze_nc*1.E18) 


!         Total composite reflectivity, including convection, in dbZe
          ze_sum = ze_nc*1.E18  ! + ze_conv
          ze_max = max(ze_max, ze_sum )

             DBZ(i,j,ll) = ze_sum
             DBZR(i,j,ll) = ze_r*1.E18
             DBZI(i,j,ll) = (ze_s+ze_g)*1.E18

           ENDDO
!         parameterized convection
!         -------------------------
!          conv_prate(i,j)  is convective pcpn rate, assumed in mm/h
!         ze_conv = 300.*conv_prate**1.4 ! Units: mm^6/m^3

          RDTPHS=3.6E6/DT
          CUPRATE=RDTPHS*CPRATE(I,J)            !--- Cu precip rate, R (mm/h)

!        ze_conv= max(0.1,300*(4.*CUPRATE)**1.4)
! -- switch to time-step conv precip in RR
        ze_conv= max(0.1,300*(CUPRATE)**1.4)

!  Combine max resolved reflectivity component
!    and sub-grid scale component
!         Total composite reflectivity, including convection, in dbZe
          ze_sum = ze_max  + ze_conv
          refl(i,j) = 10.*LOG10(ze_sum)
!          refl1km(i,j) = 10.*LOG10(ze_nc_1km*1.E18 + ze_conv)
!          refl4km(i,j) = 10.*LOG10(ze_nc_4km*1.E18 + ze_conv)

          ENDDO
         ENDDO

       ze_rmax = 10.*log10(ze_rmax*1.e18)
       ze_smax = 10.*log10(ze_smax*1.e18)
       ze_gmax = 10.*log10(ze_gmax*1.e18)

       write (6,*) 'dbze_max-r/s/g',ze_rmax,ze_smax,ze_gmax
      ENDIF     !tgs endif for Thompson scheme

      END IF
!     
!     OUTPUT/CALCULATE PRESSURE, OMEGA, POTENTIAL TEMPERATURE,
!     DEWPOINT TEMPERATURE, RELATIVE HUMIDITY, AND 
!     ABSOLUTE VORTICITY ON MDL SURFACES.
!     
!
      IF ( (IGET(001).GT.0).OR.(IGET(077).GT.0).OR.      &
           (IGET(002).GT.0).OR.(IGET(003).GT.0).OR.      &
           (IGET(004).GT.0).OR.(IGET(005).GT.0).OR.      &
           (IGET(006).GT.0).OR.(IGET(083).GT.0).OR.      &
           (IGET(007).GT.0).OR.(IGET(008).GT.0).OR.      &
           (IGET(009).GT.0).OR.(IGET(010).GT.0).OR.      &
           (IGET(084).GT.0).OR.(IGET(011).GT.0).OR.      &
           (IGET(041).GT.0).OR.(IGET(124).GT.0).OR.      &
           (IGET(078).GT.0).OR.(IGET(079).GT.0).OR.      &
           (IGET(125).GT.0).OR.(IGET(145).GT.0).OR.      &
           (IGET(140).GT.0).OR.(IGET(040).GT.0).OR.      &
           (IGET(181).GT.0).OR.(IGET(182).GT.0).OR.      &
           (IGET(199).GT.0).OR.(IGET(185).GT.0).OR.      &
           (IGET(186).GT.0).OR.(IGET(187).GT.0).OR.      &
           (IGET(250).GT.0).OR.(IGET(252).GT.0).OR.      &
           (IGET(276).GT.0).OR.(IGET(277).GT.0).OR.      &
           (IGET(278).GT.0).OR.(IGET(264).GT.0) )  THEN

      DO 190 L=1,LM

!           PRESSURE ON MDL SURFACES.
            IF (IGET(001).GT.0) THEN
             IF (LVLS(L,IGET(001)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=PMID(I,J,LL)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(001),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!
!---  CLOUD WATER on MDL SURFACE  (Jin, '01; Ferrier, Feb '02)
!
          IF (IGET(124) .GT. 0) THEN
            IF (LVLS(L,IGET(124)) .GT. 0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=QQW(I,J,LL)
                 if(GRID1(I,J)<1e-20) GRID1(I,J)=0.0
               ENDDO
               ENDDO	    
               ID(1:25) = 0
	       CALL GRIBIT(IGET(124),L,GRID1,IM,JM)
            ENDIF
          ENDIF 
!
!---  CLOUD ICE ON MDL SURFACE  (Jin, '01; Ferrier, Feb '02)
!
          IF (IGET(125) .GT. 0) THEN
            IF (LVLS(L,IGET(125)) .GT. 0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=QQI(I,J,LL)
                 if(GRID1(I,J)<1e-20) GRID1(I,J)=0.0
               ENDDO
               ENDDO	   	    
               ID(1:25) = 0
	       CALL GRIBIT(IGET(125),L,GRID1,IM,JM)
            ENDIF
          ENDIF
!
!---  RAIN ON MDL SURFACE  (Jin, '01; Ferrier, Feb '02)
!
          IF (IGET(181) .GT. 0) THEN
            IF (LVLS(L,IGET(181)) .GT. 0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=QQR(I,J,LL)
                 if(GRID1(I,J)<1e-20) GRID1(I,J)=0.0
               ENDDO
               ENDDO
               ID(1:25) = 0
	       CALL GRIBIT(IGET(181),L,GRID1,IM,JM)
            ENDIF
          ENDIF
!
!---  SNOW ON MDL SURFACE  (Jin, '01; Ferrier, Feb '02)
!
          IF (IGET(182) .GT. 0) THEN
            IF (LVLS(L,IGET(182)) .GT. 0)THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=QQS(I,J,LL)
                 if(GRID1(I,J)<1e-20) GRID1(I,J)=0.0
               ENDDO
               ENDDO	    
               ID(1:25) = 0
	       CALL GRIBIT(IGET(182),L,GRID1,IM,JM)
            ENDIF
          ENDIF
!
!---  GRAUPEL ON MDL SURFACE   --tgs
!
          IF (IGET(415) .GT. 0) THEN
            IF (LVLS(L,IGET(415)) .GT. 0)THEN
               LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
            if(QQG(I,J,LL).lt.1.e-12)QQG(I,J,LL)=0.     !tgs
                 GRID1(I,J)=QQG(I,J,LL)
               ENDDO
               ENDDO    
               ID(1:25) = 0
               CALL GRIBIT(IGET(415),L,GRID1,IM,JM)
            ENDIF
          ENDIF

!
!---  Total cloud fraction on MDL surfaces.  (Ferrier, Nov '04)
!
          IF (IGET(145) .GT. 0) THEN
            IF (LVLS(L,IGET(145)) .GT. 0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 IF(abs(CFR(I,J,LL)-SPVAL).GT.SMALL)     &
     &                 GRID1(I,J)=CFR(I,J,LL)*H100
               ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H100)
               ID(1:25) = 0
               CALL GRIBIT(IGET(145),L,GRID1,IM,JM)
            ENDIF
          ENDIF
!
!---  Equivalent radar reflectivity factor.  
!
          IF (IGET(250) .GT. 0) THEN
            IF (LVLS(L,IGET(250)) .GT. 0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=DBZ(I,J,LL)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,DBZmin,DBZmax)
               ID(1:25) = 0
	       ID(02)=129
               CALL GRIBIT(IGET(250),L,GRID1,IM,JM) 
            ENDIF
          ENDIF

!
!--- TOTAL CONDENSATE ON MDL SURFACE (CWM array; Ferrier, Feb '02)
!
          IF (IGET(199).GT.0) THEN
            IF (LVLS(L,IGET(199)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=CWM(I,J,LL)
               ENDDO
               ENDDO	     
               ID(1:25) = 0
               ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
	       CALL GRIBIT(IGET(199),L,GRID1,IM,JM)
            ENDIF
          ENDIF
!
!---  F_rain ON MDL SURFACE  (Jin, '01; Ferrier, Feb '02)
!
          IF (IGET(185).GT.0) THEN
            IF (LVLS(L,IGET(185)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=F_rain(I,J,LL)
               ENDDO
               ENDDO	    
               ID(1:25) = 0
               ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
	       CALL GRIBIT(IGET(185),L,GRID1,IM,JM)
            ENDIF
          ENDIF
!
!---  F_ice ON MDL SURFACE  (Jin, '01; Ferrier, Feb '02)
!
          IF (IGET(186).GT.0) THEN
            IF (LVLS(L,IGET(186)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=F_ice(I,J,LL)
               ENDDO
               ENDDO	    
               ID(1:25) = 0
               ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
	       CALL GRIBIT(IGET(186),L,GRID1,IM,JM)
            ENDIF
          ENDIF
!
!---  F_RimeF ON MDL SURFACE  (Jin, '01; Ferrier, Feb '02)
!
          IF (IGET(187).GT.0) THEN
            IF (LVLS(L,IGET(187)).GT.0) THEN
!--- Filter "rime factor" for non-zero precip rates and % frozen precip
              LL=LM-L+1
              DO J=JSTA,JEND
                DO I=1,IM
                 GRID1(I,J)=F_RimeF(I,J,LL)		 
                ENDDO
              ENDDO
              ID(1:25) = 0
              ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
	      CALL GRIBIT(IGET(187),L,GRID1,IM,JM)
            ENDIF
          ENDIF
!	  
!           HEIGHTS ON MDL SURFACES.
            IF (IGET(077).GT.0) THEN
             IF (LVLS(L,IGET(077)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=ZMID(I,J,LL)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(077),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!           TEMPERATURE ON MDL SURFACES.
            IF (IGET(002).GT.0) THEN
             IF (LVLS(L,IGET(002)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=T(I,J,LL)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(002),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!           POTENTIAL TEMPERATURE ON MDL SURFACES.
            IF (IGET(003).GT.0) THEN
             IF (LVLS(L,IGET(003)).GT.0) THEN
              LL=LM-L+1
!              IF(MODELNAME .EQ. 'NCAR' .OR. MODELNAME == 'RAPR' .OR. MODELNAME == 'RAPR')THEN
!               DO J=JSTA,JEND
!               DO I=1,IM
!                 GRID1(I,J)=TH(I,J,LL)
!               ENDDO
!               ENDDO
!              ELSE
               DO J=JSTA,JEND
               DO I=1,IM
                 P1D(I,J)=PMID(I,J,LL)
                 T1D(I,J)=T(I,J,LL)
               ENDDO
               ENDDO
               CALL CALPOT(P1D,T1D,EGRID3)
                                                                                
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID3(I,J)
               ENDDO
               ENDDO
!              END IF
              ID(1:25) = 0
              CALL GRIBIT(IGET(003),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!           RELATIVE HUMIDITY ON MDL SURFACES.
            IF (IGET(006).GT.0) THEN
             IF (LVLS(L,IGET(006)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 P1D(I,J)=PMID(I,J,LL)
                 T1D(I,J)=T(I,J,LL)
                 Q1D(I,J)=Q(I,J,LL)
               ENDDO
               ENDDO
	       IF(MODELNAME == 'GFS')THEN
	        CALL CALRH_GFS(P1D,T1D,Q1D,EGRID4)
	       ELSE IF (MODELNAME == 'RAPR')THEN
                CALL CALRH_GSD(P1D,T1D,Q1D,EGRID4)
               ELSE
                CALL CALRH(P1D,T1D,Q1D,EGRID4)
	       END IF               
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID4(I,J)*100.
		 EGRID2(I,J)=Q(I,J,LL)/EGRID4(I,J) ! Revert QS to compute cloud cover later
               ENDDO
               ENDDO
!               CALL BOUND(GRID1,H1,H100)
               ID(1:25) = 0
               CALL GRIBIT(IGET(006),L,GRID1,IM,JM)
             ENDIF
            ENDIF

!     
!           DEWPOINT ON MDL SURFACES.
            IF (IGET(004).GT.0) THEN
             IF (LVLS(L,IGET(004)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 P1D(I,J)=PMID(I,J,LL)
                 T1D(I,J)=T(I,J,LL)
                 Q1D(I,J)=Q(I,J,LL)
               ENDDO
               ENDDO
               CALL CALDWP(P1D,Q1D,EGRID3,T1D)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID3(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(004),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!           SPECIFIC HUMIDITY ON MDL SURFACES.
            IF (IGET(005).GT.0) THEN
             IF (LVLS(L,IGET(005)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=Q(I,J,LL)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,H1M12,H99999)
               ID(1:25) = 0
               CALL GRIBIT(IGET(005),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!           MOISTURE CONVERGENCE ON MDL SURFACES.
            IF (IGET(083).GT.0 .OR. IGET(295).GT.0) THEN
             IF (LVLS(L,IGET(083)).GT.0 .OR. IGET(295).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA_2L,JEND_2U
               DO I=1,IM
                 Q1D(I,J)=Q(I,J,LL)
                 EGRID1(I,J)=UH(I,J,LL)
                 EGRID2(I,J)=VH(I,J,LL)
               ENDDO
               ENDDO
               CALL CALMCVG(Q1D,EGRID1,EGRID2,EGRID3)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID3(I,J)
		 MCVG(I,J,LL)=EGRID3(I,J)
               ENDDO
               ENDDO
	       IF(IGET(083).GT.0 .AND. LVLS(L,IGET(083)).GT.0)THEN
                ID(1:25) = 0
                CALL GRIBIT(IGET(083),L,GRID1,IM,JM)
	       END IF	
             ENDIF
            ENDIF
!     
!           U AND/OR V WIND ON MDL SURFACES.
!MEB needs to be modified to do u at u-points and v at v-points
            IF (IGET(007).GT.0.OR.IGET(008).GT.0) THEN
             IF (LVLS(L,IGET(007)).GT.0.OR.LVLS(L,IGET(008)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=UH(I,J,LL)
                 GRID2(I,J)=VH(I,J,LL)
               ENDDO
               ENDDO
               ID(1:25) = 0
               IF (IGET(007).GT.0) CALL GRIBIT(IGET(007),L,GRID1,IM,JM)
               ID(1:25) = 0
               IF (IGET(008).GT.0) CALL GRIBIT(IGET(008),L,GRID2,IM,JM)
             ENDIF
            ENDIF
!     
!           OMEGA ON MDL SURFACES.
            IF (IGET(009).GT.0) THEN
             IF (LVLS(L,IGET(009)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=OMGA(I,J,LL)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(009),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!           W ON MDL SURFACES.
            IF (IGET(264).GT.0) THEN
             IF (LVLS(L,IGET(264)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=WH(I,J,LL)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(264),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!           ABSOLUTE VORTICITY ON MDL SURFACES.
            IF (IGET(010).GT.0) THEN
             IF (LVLS(L,IGET(010)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA_2L,JEND_2U
               DO I=1,IM
                 EGRID1(I,J)=UH(I,J,LL)
                 EGRID2(I,J)=VH(I,J,LL)
               ENDDO
               ENDDO
               CALL CALVOR(EGRID1,EGRID2,EGRID3)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID3(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(010),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!           GEOSTROPHIC STREAMFUNCTION ON MDL SURFACES.
            IF (IGET(084).GT.0) THEN
             IF (LVLS(L,IGET(084)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 EGRID1(I,J)=ZMID(I,J,LL)
               ENDDO
               ENDDO
               CALL CALSTRM(EGRID1,EGRID2)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EGRID2(I,J)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(084),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!           TURBULENT KINETIC ENERGY ON MDL SURFACES.
            IF (IGET(011).GT.0) THEN
             IF (LVLS(L,IGET(011)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=Q2(I,J,LL)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(011),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!    
!           CLOUD WATER CONTENT
!HC            IF (IGET(124).GT.0) THEN
!HC             IF (LVLS(L,IGET(124)).GT.0) THEN
!HC              DO J=JSTA,JEND
!HC              DO I=1,IM
!HC                IF(CWM(I,J,L).LT.0..AND.CWM(I,J,L).GT.-1.E-10)
!HC     1            CWM(I,J,L)=0.
!HC                 GRID1(I,J)=CWM(I,J,L)
!HC              ENDDO
!HC              ENDDO
!HC              ID(1:25) = 0
!HC              CALL GRIBIT(IGET(124),L,GRID1,IM,JM)
!HC             ENDIF
!HC            ENDIF
!     
!           CLOUD ICE CONTENT.
!commented out until QICE is brought into post
!           IF (IGET(125).GT.0) THEN
!            IF (LVLS(L,IGET(125)).GT.0) THEN
!              DO J=JSTA,JEND
!              DO I=1,IM
!                GRID1(I,J)=QICE(I,J,L)
!              ENDDO
!              ENDDO
!              ID(1:25) = 0
!              CALL GRIBIT(IGET(125),L,GRID1,IM,JM)
!            ENDIF
!           ENDIF
!     
!           CLOUD FRACTION
!     
!commented out until CFRC is brought into post
!           IF (IGET(145).GT.0) THEN
!            IF (LVLS(L,IGET(145)).GT.0) THEN
!              DO J=JSTA,JEND
!              DO I=1,IM
!                GRID1(I,J)=CFRC(I,J,L)
!              ENDDO
!              ENDDO
!              ID(1:25) = 0
!              CALL GRIBIT(IGET(145),L,GRID1,IM,JM)
!            ENDIF
!           ENDIF
!     
!           TEMPERATURE TENDENCY DUE TO RADIATIVE FLUX CONVERGENCE
!commented out until TTND is brought into post
           IF (IGET(140).GT.0) THEN
            IF (LVLS(L,IGET(140)).GT.0) THEN
	      LL=LM-L+1
              DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J)=TTND(I,J,LL)
              ENDDO
              ENDDO
                 ID(1:25) = 0
                 CALL GRIBIT(IGET(140),L,GRID1,IM,JM)
            ENDIF
           ENDIF
!     
!           TEMPERATURE TENDENCY DUE TO SHORT WAVE RADIATION.
!commented out until RSWTT is brought into post
           IF (IGET(040).GT.0) THEN
            IF (LVLS(L,IGET(040)).GT.0) THEN
	      LL=LM-L+1
              DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J)=RSWTT(I,J,LL)
              ENDDO
              ENDDO
                 ID(1:25) = 0
                 CALL GRIBIT(IGET(040),L,GRID1,IM,JM)
            ENDIF
           ENDIF
!     
!           TEMPERATURE TENDENCY DUE TO LONG WAVE RADIATION.
!commented out until RLWTT is brought into post
           IF (IGET(041).GT.0) THEN
            IF (LVLS(L,IGET(041)).GT.0) THEN
	      LL=LM-L+1
              DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J)=RLWTT(I,J,LL)
              ENDDO
              ENDDO
                 ID(1:25) = 0
                 CALL GRIBIT(IGET(041),L,GRID1,IM,JM)
            ENDIF
           ENDIF
!
!     
!        PROCESS NEXT MDL LEVEL.
!
!           LATENT HEATING FROM GRID SCALE RAIN/EVAP. (TIME AVE)
           IF (IGET(078).GT.0) THEN
            IF (LVLS(L,IGET(078)).GT.0) THEN
	       LL=LM-L+1 
               IF(AVRAIN.GT.0.)THEN 
                 RRNUM=1./AVRAIN
               ELSE
                 RRNUM=0.
               ENDIF
!$omp  parallel do
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=TRAIN(I,J,LL)*RRNUM
               ENDDO
               ENDDO
               ID(1:25) = 0
               ITHEAT     = INT(THEAT)
	       IF (ITHEAT .NE. 0) THEN
                IFINCR     = MOD(IFHR,ITHEAT)
	       ELSE
	        IFINCR=0
	       END IF		
               ID(19) = IFHR
	       IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
               ID(20) = 3
               IF (IFINCR.EQ.0) THEN
                  ID(18) = IFHR-ITHEAT
               ELSE
                  ID(18) = IFHR-IFINCR
               ENDIF
	       IF(IFMIN .GE. 1)ID(18)=ID(18)*60
               IF (ID(18).LT.0) ID(18) = 0
               CALL GRIBIT(IGET(078),L,GRID1,IM,JM)
            END IF
           ENDIF
!
!           LATENT HEATING FROM CONVECTION. (TIME AVE)
           IF (IGET(079).GT.0) THEN
            IF (LVLS(L,IGET(079)).GT.0) THEN
	       LL=LM-L+1 
               IF(AVCNVC.GT.0.)THEN
                 RRNUM=1./AVCNVC
               ELSE
                 RRNUM=0.
               ENDIF
!$omp  parallel do
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = TCUCN(I,J,LL)*RRNUM
               ENDDO
               ENDDO
               ID(1:25) = 0
               ITHEAT     = INT(THEAT)
	       IF (ITHEAT .NE. 0) THEN
                IFINCR     = MOD(IFHR,ITHEAT)
	       ELSE
	        IFINCR=0
	       END IF	
               ID(19) = IFHR
	       IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
               ID(20) = 3
               IF (IFINCR.EQ.0) THEN
                  ID(18) = IFHR-ITHEAT
               ELSE
                  ID(18) = IFHR-IFINCR
               ENDIF
	       IF(IFMIN .GE. 1)ID(18)=ID(18)*60
               IF (ID(18).LT.0) ID(18) = 0
               CALL GRIBIT(IGET(079),L,GRID1,IM,JM)
            END IF
           ENDIF
	   
!
!           OZONE
           IF (IGET(267).GT.0) THEN
            IF (LVLS(L,IGET(267)).GT.0) THEN
	       LL=LM-L+1 
               ID(1:25) = 0
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = O3(I,J,LL)
               ENDDO
               ENDDO
               CALL GRIBIT(IGET(267),L,GRID1,IM,JM)
            END IF
           ENDIF
	   
 190     CONTINUE
!
!     END OF MDL SURFACE OUTPUT BLOCK.
!
      ENDIF
!   VISIBILITY
!     IF (IGET(180).GT.0) THEN
!comment out until we get QICE, QSNOW brought into post
!MEB   RDTPHS= 1./(NPHS*DT)
!MEB modifying this Eta-specific code, assuming WRF physics will
!MEB explicitly predict vapor/water/ice/rain/snow
!MEB comments starting with MEB are lines associated with this
!MEB Eta-specific code
!            NEED TO CALCULATE RAIN WATER AND SNOW MIXING RATIOS
!      DO J=JSTA,JEND
!      DO I=1,IM
!MEB     IF (PREC(I,J).EQ.0) THEN
!MEB       QSNO(I,J)=0.
!MEB       QRAIN(I,J)=0.
!MEB     ELSE
!MEB       LLMH=LMH(I,J)
!MEB       SNORATE=SR(I,J)*PREC(I,J)*RDTPHS
!MEB       RAINRATE=(1-SR(I,J))*PREC(I,J)*RDTPHS
!MEB       TERM1=(T(I,J,LM)/PSLP(I,J))**0.4167
!MEB       TERM2=(T(I,J,LLMH)/PMID(I,J,LMH(I,J)))**0.5833
!MEB       TERM3=RAINRATE**0.8333
!MEB       QRAIN(I,J)=RAINCON*TERM1*TERM2*TERM3
!MEB       TERM4=(T(I,J,LM)/PSLP(I,J))**0.47
!MEB       TERM5=(T(I,J,LLMH)/PMID(I,J,LMH(I,J)))**0.53
!MEB       TERM6=SNORATE**0.94
!MEB       QSNO(I,J)=SNOCON*TERM4*TERM5*TERM6
!MEB     ENDIF
!        LLMH=NINT(LMH(I,J))
!        QRAIN1(I,J)=QRAIN(I,J,LLMH)
!        QSNO1(I,J)=QSNOW(I,J,LLMH)
!        TT(I,J)=T(I,J,LLMH)
!        QV(I,J)=Q(I,J,LLMH)
!        QCD(I,J)=CWM(I,J,LLMH)
!        QICE1(I,J)=QICE(I,J,LLMH)
!        PPP(I,J)=PMID(I,J,LLMH)
!      ENDDO
!      ENDDO
!      CALL CALVIS(QV,QCD,QRAIN1,QICE1,QSNO1,TT,PPP,VIS)
!              DO J=JSTA,JEND
!              DO I=1,IM
!                GRID1(I,J)=VIS(I,J)
!              ENDDO
!              ENDDO
!      ID(1:25) = 0
!      CALL GRIBIT(IGET(180),LVLS(1,IGET(180)),
!    X           GRID1,IM,JM)
!      ENDIF
!
!     INSTANTANEOUS CONVECTIVE PRECIPITATION RATE.
!
!      IF (IGET(249).GT.0) THEN
!         RDTPHS=1000./DTQ2
!         DO J=JSTA,JEND
!         DO I=1,IM
!           GRID1(I,J)=CPRATE(I,J)*RDTPHS
!           GRID1(I,J)=SPVAL
!         ENDDO
!         ENDDO
!         ID(1:25) = 0
!	 CALL GRIBIT(IGET(249),LM,GRID1,IM,JM)
!      ENDIF
!
!     COMPOSITE RADAR REFLECTIVITY (maximum dBZ in each column)
!
      IF (IGET(252).GT.0) THEN
        IF(IMP_PHYSICS.NE.8) THEN
         DO J=JSTA,JEND
            DO I=1,IM
               GRID1(I,J)=DBZmin
               DO L=1,NINT(LMH(I,J))
                  GRID1(I,J)=MAX( GRID1(I,J), DBZ(I,J,L) )
               ENDDO
            ENDDO
         ENDDO
         ELSE
!tgs - for Thompson scheme
         DO J=JSTA,JEND
            DO I=1,IM
               GRID1(I,J)=refl(i,j)
            ENDDO
         ENDDO
        ENDIF
         ID(1:25) = 0
	 ID(02)=129
         CALL GRIBIT(IGET(252),LM,GRID1,IM,JM)
      ENDIF
!
!--   COMPOSITE RADAR REFLECTIVITY FROM RAIN (maximum dBZ in each column due to rain)
!
      IF (IGET(276).GT.0) THEN
         DO J=JSTA,JEND
            DO I=1,IM
               GRID1(I,J)=DBZmin
               DO L=1,NINT(LMH(I,J))
                  GRID1(I,J)=MAX( GRID1(I,J), DBZR(I,J,L) )
               ENDDO
            ENDDO
         ENDDO
         ID(1:25) = 0
         ID(02)=129
         CALL GRIBIT(IGET(276),LM,GRID1,IM,JM)
      ENDIF
!
!--   COMPOSITE RADAR REFLECTIVITY FROM ICE
!     (maximum dBZ in each column due to all ice habits; snow + graupel + etc.)
!
      IF (IGET(277).GT.0) THEN
         DO J=JSTA,JEND
            DO I=1,IM
               GRID1(I,J)=DBZmin
               DO L=1,NINT(LMH(I,J))
                  GRID1(I,J)=MAX( GRID1(I,J), DBZI(I,J,L) )
               ENDDO
            ENDDO
         ENDDO
         ID(1:25) = 0
         ID(02)=129
         CALL GRIBIT(IGET(277),LM,GRID1,IM,JM)
      ENDIF
!
!--   COMPOSITE RADAR REFLECTIVITY FROM PARAMETERIZED CONVECTION
!     (maximum dBZ in each column due to parameterized convection, as bogused into
!      post assuming a constant reflectivity from the surface to the 0C level, 
!      and decreasing with height at higher levels)
!
      IF (IGET(278).GT.0) THEN
         DO J=JSTA,JEND
            DO I=1,IM
               GRID1(I,J)=DBZmin
               DO L=1,NINT(LMH(I,J))
                  GRID1(I,J)=MAX( GRID1(I,J), DBZC(I,J,L) )
               ENDDO
            ENDDO
         ENDDO
         ID(1:25) = 0
         ID(02)=129
         CALL GRIBIT(IGET(278),LM,GRID1,IM,JM)
      ENDIF
! SRD -- converted to kft
! J.Case, ENSCO Inc. (5/26/2008) -- Output Echo Tops (Highest HGT in meters
! of the 18-dBZ reflectivity on a model level)

      IF (IGET(426).GT.0) THEN
         DO J=JSTA,JEND
            DO I=1,IM
               GRID1(I,J)=0.0
               DO L=1,NINT(LMH(I,J))
                  IF (DBZ(I,J,L).GE.18.0) THEN
                     GRID1(I,J)=ZMID(I,J,L)*3.2808/1000.
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(426),LM,GRID1,IM,JM)
      ENDIF

! J.Case (end mods)
! SRD

!
!---   VISIBILITY
!
      IF (IGET(180).GT.0) THEN
        RDTPHS=1./DTQ2
  !
  !--- Needed values at 1st level above ground  (Jin, '01; Ferrier, Feb '02)
  !
        DO J=JSTA,JEND
          DO I=1,IM
            LLMH=NINT(LMH(I,J))
            Q1D(I,J)=Q(I,J,LLMH)
           if(Q1D(I,J).le.0.) Q1D(I,J)=0.         !tgs
            QW1(I,J)=QQW(I,J,LLMH)
            QR1(I,J)=QQR(I,J,LLMH)
            QI1(I,J)=QQI(I,J,LLMH)
            QS1(I,J)=QQS(I,J,LLMH)
            QG1(I,J)=QQG(I,J,LLMH)      !tgs
            T1D(I,J)=T(I,J,LLMH)
            P1D(I,J)=PMID(I,J,LLMH)
!HC Because instantanous convective precip rate is not yet available as wrf output,
!HC cuppt is used as a replacement for now  
!HC Only adding convective precipitation rate when using Ferrier's scheme
           IF(imp_physics.eq.5)THEN
            IF (CPRATE(I,J) .GT. 0.) THEN
!            IF (CUPPT(I,J) .GT. 0.) THEN
               RAINRATE=(1-SR(I,J))*CPRATE(I,J)*RDTPHS
!               RAINRATE=(1-SR(I,J))*CUPPT(I,J)/(TRDLW*3600.)
               TERM1=(T(I,J,LM)/PMID(I,J,LM))**0.4167
               TERM2=(T1D(I,J)/P1D(I,J))**0.5833
               TERM3=RAINRATE**0.8333
	       QROLD=1.2*QR1(I,J)
               QR1(I,J)=QR1(I,J)+RAINCON*TERM1*TERM2*TERM3
               IF (SR(I,J) .GT. 0.) THEN
                  SNORATE=SR(I,J)*CPRATE(I,J)*RDTPHS
!                  SNORATE=SR(I,J)*CUPPT(I,J)/(TRDLW*3600.)
                  TERM1=(T(I,J,LM)/PMID(I,J,LM))**0.47
                  TERM2=(T1D(I,J)/P1D(I,J))**0.53
                  TERM3=SNORATE**0.94
                  QS1(I,J)=QS1(I,J)+SNOCON*TERM1*TERM2*TERM3
               ENDIF
            ENDIF
	   END IF 
	   
	   IF(imp_physics.eq.99)THEN ! use rain rate for visibility
            IF (prec(i,j) < spval .and. prec(I,J) > 0.) THEN
!            IF (CUPPT(I,J) .GT. 0.) THEN
               RAINRATE=(1-SR(I,J))*PREC(I,J)*RDTPHS
!               RAINRATE=(1-SR(I,J))*CUPPT(I,J)/(TRDLW*3600.)
               TERM1=(T(I,J,LM)/PMID(I,J,LM))**0.4167
               TERM2=(T1D(I,J)/P1D(I,J))**0.5833
               TERM3=RAINRATE**0.8333
	       QROLD=1.2*QR1(I,J)
               QR1(I,J)=QR1(I,J)+RAINCON*TERM1*TERM2*TERM3
               IF (sr(i,j) < spval .and. SR(I,J) > 0.) THEN
                  SNORATE=SR(I,J)*PREC(I,J)*RDTPHS
!                  SNORATE=SR(I,J)*CUPPT(I,J)/(TRDLW*3600.)
                  TERM1=(T(I,J,LM)/PMID(I,J,LM))**0.47
                  TERM2=(T1D(I,J)/P1D(I,J))**0.53
                  TERM3=SNORATE**0.94
                  QS1(I,J)=QS1(I,J)+SNOCON*TERM1*TERM2*TERM3
               ENDIF
            ENDIF
	   END IF
	   
          ENDDO
        ENDDO
  !
  !-- Visibility using Warner-Stoelinga algorithm  (Jin, 01)
  !
        ii=im/2
        jj=(jsta+jend)/2
!        print*,'Debug: Visbility ',Q1D(ii,jj),QW1(ii,jj),QR1(ii,jj)
!     +,QI1(ii,jj) ,QS1(ii,jj),T1D(ii,jj),P1D(ii,jj)

        CALL CALVIS(Q1D,QW1,QR1,QI1,QS1,T1D,P1D,VIS)
        print*,'Debug: Visbility ',VIS(ii,jj)

!        print*,'Debug: Visbility ',Q1D(ii,jj),QW1(ii,jj),QR1(ii,jj),QI1(ii,jj)
!     +,QS1(ii,jj),T1D(ii,jj),P1D(ii,jj)
	DO J=JSTA,JEND
	DO I=1,IM
	  IF(abs(vis(i,j)).gt.24135.1)print*,'bad visbility'    &
       ,i,j,Q1D(i,j),QW1(i,j),QR1(i,j),QI1(i,j)                 &
       ,QS1(i,j),T1D(i,j),P1D(i,j),vis(i,j)	  
	  GRID1(I,J)=VIS(I,J)
	END DO
	END DO  
        ID(1:25) = 0
	CALL GRIBIT(IGET(180),LM,GRID1,IM,JM)
       ENDIF

!
! --- GSD VISIBILITY
!
      IF (IGET(410).GT.0) THEN
        CALL CALVIS_GSD(VIS)
        DO J=JSTA,JEND
        DO I=1,IM
          GRID1(I,J)=VIS(I,J)
        END DO
        END DO
        ID(1:25) = 0
        CALL GRIBIT(IGET(410),LM,GRID1,IM,JM)
       ENDIF

!
!     
!     ASYMPTOTIC AND FREE ATMOSPHERE MASTER LENGTH SCALE (EL), PLUS
!     GRADIENT RICHARDSON NUMBER.
!
      IF ( (IGET(111).GT.0) .OR. (IGET(146).GT.0) .OR.           &
           (IGET(147).GT.0) ) THEN
!     
!        COMPUTE ASYMPTOTIC MASTER LENGTH SCALE.
         CALL CLMAX(EL0,EGRID2,EGRID3,EGRID4,EGRID5)
!     
!        IF REQUESTED, POST ASYMPTOTIC MASTER LENGTH SCALE.
         IF (IGET(147).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EL0(I,J)
               ENDDO
               ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(147),LM,GRID1,IM,JM)
         ENDIF
!     
!        IF REQUESTED, POST FREE ATMOSPHERE MASTER LENGTH SCALE
!        AND/OR THE GRADIENT RICHARDSON NUMBER.    
!
         IF ( (IGET(111).GT.0) .OR. (IGET(146).GT.0) ) THEN
!     
!           COMPUTE FREE ATMOSPHERE MASTER LENGTH SCALE.
!$omp  parallel do
            DO L=1,LM
               DO J=JSTA,JEND
               DO I=1,IM
                 EL(I,J,L)=D00
               ENDDO
               ENDDO
            ENDDO

            IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM' .OR. MODELNAME == 'RAPR')THEN
!             CALL MIXLEN(EL0,EL)  
            ELSE IF(MODELNAME .EQ. 'NMM')THEN
              DO L=1,LM
               DO J=JSTA,JEND
               DO I=1,IM
                 EL(I,J,L)=EL_MYJ(I,J,L)  !NOW EL COMES OUT OF WRF NMM
               ENDDO
               ENDDO
              ENDDO
            END IF
!     
!           COMPUTE GRADIENT RICHARDSON NUMBER IF REQUESTED.
!     
            IF ( (IGET(111).GT.0) ) CALL CALRCH(EL,RICHNO)
!
!           LOOP OVER MDL LAYERS.
            DO 200 L = 1,LM
!     
!              POST MIXING LENGTH.
!
            IF (IGET(146).GT.0) THEN
             IF (LVLS(L,IGET(146)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=EL(I,J,LL)
               ENDDO
               ENDDO
                  ID(1:25) = 0
                  CALL GRIBIT(IGET(146),L,GRID1,IM,JM)
             ENDIF
            ENDIF
!     
!              POST GRADIENT RICHARDSON NUMBER.
!
            IF(L .LT. LM)THEN
             IF (IGET(111).GT.0) THEN
              IF (LVLS(L,IGET(111)).GT.0) THEN
	       LL=LM-L+1
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RICHNO(I,J,LL)
               ENDDO
               ENDDO
               ID(1:25) = 0
               CALL GRIBIT(IGET(111),L,GRID1,IM,JM)
              ENDIF
             ENDIF
            END IF
 200        CONTINUE
!
!
         ENDIF
      ENDIF
!     
!           COMPUTE PBL HEIGHT BASED ON RICHARDSON NUMBER
!     
            IF ( (IGET(289).GT.0) ) CALL CALPBL(PBLRI)

            IF (IGET(289).GT.0) THEN
                DO J=JSTA,JEND
                DO I=1,IM
                     GRID1(I,J)=PBLRI(I,J)
!                     PBLH(I,J)=PBLRI(I,J)
                ENDDO
                ENDDO
                ID(1:25) = 0
!		ID(02)=129
                CALL GRIBIT(IGET(289),LM,GRID1,IM,JM)
            ENDIF
!     
!           COMPUTE PBL REGIME BASED ON WRF version of BULK RICHARDSON NUMBER
!     

            IF (IGET(344).GT.0) THEN
	        allocate(PBLREGIME(im,jsta_2l:jend_2u))
	        CALL CALPBLREGIME(PBLREGIME)
		print*,'back from callpblregime'
                DO J=JSTA,JEND
                DO I=1,IM
                  GRID1(I,J)=PBLREGIME(I,J)
                ENDDO
                ENDDO
                ID(1:25) = 0
		ID(02)=129
                CALL GRIBIT(IGET(344),LM,GRID1,IM,JM)
		deallocate(PBLREGIME)
            ENDIF
!
!     RADAR ECHO TOP (highest 18.3 dBZ level in each column)
!
      IF(IGET(400).GT.0)THEN
        DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J)=SPVAL	       	      
            DO L=1,NINT(LMH(I,J))
	      IF(DBZ(I,J,L)>18.3)then
	        GRID1(I,J)=ZMID(I,J,L)
		go to 201
	      END IF  
            ENDDO
 201        CONTINUE 	       
!	       if(grid1(i,j)<0.)print*,'bad echo top',
!     +           i,j,grid1(i,j),dbz(i,j,1:lm)	       
          ENDDO
        ENDDO
        ID(1:25) = 0
	ID(02)=129
        CALL GRIBIT(IGET(400),LM,GRID1,IM,JM)
      ENDIF	    
!     

      DEALLOCATE(EL)
      DEALLOCATE(RICHNO)
      DEALLOCATE(PBLRI)
      print*,'getting out of MDLFLD'
!     
!     END OF ROUTINE.
!     
      RETURN
      END
