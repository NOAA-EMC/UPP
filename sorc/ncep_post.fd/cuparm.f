!> @file
!> @brief module: cuparm_mod defines variables used for cumulus parameterization
  module cuparm_mod
!
    implicit none
!
    real,parameter :: H1=1.E0, &    !< _____
      H1D5=1.5E0,       &     !< _____
      H2D5=2.5E0,       &     !< _____
      H3000=3000.E0,    &     !< _____
      H10E5=100000.E0,  &     !< _____
      D00=0.E0,         &     !< _____
      D125=.125E0,      &     !< _____
      D50=.5E0,         &     !< _____
      D608=.608E0,      &     !< _____
      G=9.8E0,          &     !< _____
      CP=1004.6E0,      &     !< _____
      CAPA=0.28589641E0, &    !< _____
      ROG=287.04/9.8,   &     !< _____
      ELWV=2.50E6,      &     !< _____
      ELIVW=2.72E6,     &     !< _____
      ROW=1.E3,         &     !< _____
      EPSQ=2.E-12,      &     !< _____
      A2=17.2693882E0,  &     !< _____
      A3=273.16E0,      &     !< _____
      A4=35.86E0,       &     !< _____
      T0=273.16E0,      &     !< _____
      T1=274.16E0,      &     !< _____
      PQ0=379.90516E0,  &     !< _____
      STRESH=1.10E0,    &     !< _____
      STABS=1.0E0,      &     !< _____
      STABD=.90E0,      &     !< _____
      STABFC=1.00E0,    &     !< _____
      DTTOP=0.0E0,      &     !< _____
!---VVVVV
      RHF=0.10,      &        !< _____
      EPSUP=1.00,    &        !< _____
      EPSDN=1.05,    &        !< _____
      EPSTH=0.0,     &        !< _____
      PBM=13000.,    &        !< _____
      PQM=20000.,    &        !< _____
      PNO=1000.,     &        !< _____
      PONE=2500.,    &        !< _____
      ZSH=2000.,     &        !< _____
      PFRZ=15000.,   &        !< _____
      PSHU=45000.,   &        !< _____

!    &, RHF=0.20,EPSUP=0.93,EPSDN=1.00,EPSTH=0.3
!    &, RHF=0.20,EPSUP=1.00,EPSDN=1.00,EPSTH=0.3
!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
!
!    AUGUST '91: SCHEME HAVING THE OPTION OF USING DIFFERENT FAST AND
!    SLOW PROFILES FOR SEA AND FOR LAND POINTS; AND ALSO THE "SEA" AND
!    THE "LAND" SCHEME EVERYWHERE.  OVER LAND PROFILES DEPART FROM THE
!    FAST (DRY) PROFILES ONLY FOR PRECIPITATION/TIME STEP >
!    A PRESCRIBED VALUE (CURRENTLY, IN THE VERSION #3 DONE WEDNESDAY
!    18 SEPTEMBER, 1/4 INCH/24 H).  USE OF VARIOUS SWITCHES AS FOLLOWS.
!
!       THE "OLD" ("HARD", =ZAVISA OCT.1990) LAND SCHEME WITH FIXED
!       LAND PROFILES IS RUN BY
!       * SETTING OCT90=.TRUE. IN THE FIRST EXECUTABLE LINE FOLLOWING
!            THESE COMMENTS (THIS REACTIVATES EFI=H1 OVER LAND IF
!            .NOT.UNIS, AND CLDEFI(K)=EFIMN AS THE LEFTOVER CLDEFI
!            VALUE AT SWAP POINTS);
!       * DEFINING FAST LAND PROFILES SAME AS FAST SEA PROFILES (OR BY
!            CHOOSING ANOTHER SET OF LAND PROFILES ZAVISA USED AT
!            EARLIER TIME);
!       * SETTING FSL=1.; AND
!       * DEFINING STEFI (STARTING EFI) EQUAL TO AVGEFI.
!            (THE LAST THREE POINTS ARE HANDLED BY SWITCHING AROUND THE
!            "CFM" COMMENTS AT TWO PLACES)
!
!       THE "OLD,OLD" (APPR. ORIGINAL BETTS) SCHEME IS RUN BY
!       * SPECIFYING UNIL=.TRUE.;
!       * SETTING FSL=1.;
!       * SETTING OCT90=.TRUE.
!            (WITH THESE SETTINGS FAST LAND PROFILES ONLY ARE USED).
!                                                                     FM
     FSS=.85E0,EFIMN=.20E0,EFMNT=.70E0,FCC=.50,FCP=H1-FCC,              &
!
!         IN THIS VERSION 3.5, OVER LAND AND FOR THE FAST PROFILES, DSPB
!         IS PRESCRIBED TO BE 25 PERCENT DRIER THAN THE FAST SEA VALUE
!         (IN ROUGH AGREEMENT WITH BINDER, QJ, IN PRESS) WHILE DSP0 AND
!         DSPT ARE EACH 20 PERCENT DRIER THAN THE CORRESPONDING FAST
!         SEA VALUES.                WITH FSL=.875 THIS MAKES THE
!         AVERAGE OF THE FAST AND THE SLOW LAND PROFILES SOMEWHAT DRIER
!         THAN THE OCT90 FIXED LAND PROFILES.                         FM
!
     DSPBFL=-4843.75E0,DSP0FL=-7050.00E0,DSPTFL=-2250.0E0,FSL=.850E0,   &
!***   ACTIVATE THE FOLLOWING LINE IF OCT90=.TRUE. (AND COMMENT OUT THE
!***   PRECEDING LINE):
!    DSPBFL=-3875.E0,DSP0FL=-5875.E0,DSPTFL=-1875.E0,FSL=1.0E0,         &
     DSPBFS=-3875.E0,DSP0FS=-5875.E0,DSPTFS=-1875.E0,                   &
     DSPBSL=DSPBFL*FSL,DSP0SL=DSP0FL*FSL,DSPTSL=DSPTFL*FSL,             &
     DSPBSS=DSPBFS*FSS,DSP0SS=DSP0FS*FSS,DSPTSS=DSPTFS*FSS,             &
!*** NEW CONVECTION SCHEME WITH CROSSING DSP PROFILES ******************
!+-  &, UNIS=.FALSE.,EFIMN=.71E0,EFMNT=.71,FCC=0.5,FCP=H1-FCC
!+-  &, DSPBL=-3875.E0,DSP0L=-5875.E0,DSPTL=-1875.E0
!+-  &, DSPBS=-2875.E0,DSP0S=-5125.E0,DSPTS=-4875.E0
!+-  &, DSPBF=-4375.E0,DSP0F=-4375.E0,DSPTF=-1000.E0
!*** BETTS CONVECTION SCHEME *******************************************
!    &, UNIS=.FALSE.,EFIMN=.9999E0,EFMNT=.9999E0,FCC=.50,FCP=H1-FCC
!    &, DSPBL=-3875.E0,DSP0L=-5875.E0,DSPTL=-1875.E0
!    &, DSPBF=-3875.E0,DSP0F=-5875.E0,DSPTF=-1875.E0
!    &, DSPBS=-3875.E0,DSP0S=-5875.E0,DSPTS=-1875.E0
!***********************************************************************
     TREL=3000.,                 &     !< _____
     EPSNTP=.0010E0,             &     !< _____
     EFIFC=5.0E0,                &     !< _____
     AVGEFI=(EFIMN+1.E0)*.5E0,   &     !< _____
     DSPC=-3000.E0,              &     !< _____
     EPSP=1.E-7,                 &     !< _____
     STEFI=1.E0,                 &     !< _____
!*** ACTIVATE THE FOLLOWING LINE AND COMMENT OUT THE PRECEDING LINE IF
!*** OCT90=.TRUE.
!    &, STEFI=AVGEFI                      &   !< _____
     SLOPBL=(DSPBFL-DSPBSL)/(H1-EFIMN),   &   !< _____
     SLOP0L=(DSP0FL-DSP0SL)/(H1-EFIMN),   &   !< _____
     SLOPTL=(DSPTFL-DSPTSL)/(H1-EFIMN),   &   !< _____
     SLOPBS=(DSPBFS-DSPBSS)/(H1-EFIMN),   &   !< _____
     SLOP0S=(DSP0FS-DSP0SS)/(H1-EFIMN),   &   !< _____
     SLOPTS=(DSPTFS-DSPTSS)/(H1-EFIMN),   &   !< _____
     SLOPE=(H1   -EFMNT)/(H1-EFIMN)       &   !< _____
   real, parameter ::         &
     A23M4L=A2*(A3-A4)*ELWV,  &  !< _____
     ELOCP=ELIVW/CP,          &  !< _____
     CPRLG=CP/(ROW*G*ELWV),   &  !< _____
     RCP=H1/CP                &  !< _____
   logical,parameter :: &
     UNIS=.FALSE.,      &    !< _____
     UNIL=.FALSE.,      &    !< _____
     OCT90=.FALSE.      &    !< _____
  end module cuparm_mod
