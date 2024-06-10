!> @file
!> @brief module: cuparm_mod defines variables used for cumulus parameterization
  module cuparm_mod

    implicit none

    real, parameter :: H1=1.E0,  &   !< 1.0
      H1D5=1.5E0,       &     !< 1.5
      H2D5=2.5E0,       &     !< 2.5
      H3000=3000.E0,    &     !< 3000
      H10E5=100000.E0,  &     !< 100,000
      D00=0.E0,         &     !< Decimal number 0
      D125=.125E0,      &     !< Decimal number 0.125
      D50=.5E0,         &     !< Decimal number 0.5
      D608=.608E0,      &     !< Decimal number 0.608
      G=9.8E0,          &     !< Acceleration due to gravity
      CP=1004.6E0,      &     !< Specific heat capacity of dry air at constant pressure (kJ/kg-K)
      CAPA=0.28589641E0, &    !< R/Cp (universal gas constant over specific heat capacity of dry air at constant pressure)
      ROG=287.04/9.8,   &     !< RD over G - Gas constant for dry air divided by acceleration due to gravity
      ELWV=2.50E6,      &     !< Latent heat of vaporization of water
      ELIVW=2.72E6,     &     !< Latent heat of vaporization of water in J/kg, used in calculations involving energy transfer during evaporation
      ROW=1.E3,         &     !< Density (rho) of water
      EPSQ=2.E-12,      &     !< Minimum q (specific humidity) for computing precipitation type ?
      A2=17.2693882E0,  &     !< Constant used to parameterize specific humidity at 2m in WRFPOST: qs=pq0/p*exp(a2*(t-a3)/(t-a4))
      A3=273.16E0,      &     !< Constant used to parameterize specific humidity at 2m in WRFPOST: qs=pq0/p*exp(a2*(t-a3)/(t-a4))
      A4=35.86E0,       &     !< Constant used to parameterize specific humidity at 2m in WRFPOST: qs=pq0/p*exp(a2*(t-a3)/(t-a4))
      T0=273.16E0,      &     !< Triple point of water (K)
      T1=274.16E0,      &     !< 1 degree above triple point of water (K)
      PQ0=379.90516E0,  &     !< Constant used to parameterize specific humidity at 2m in WRFPOST: qs=pq0/p*exp(a2*(t-a3)/(t-a4))
      STRESH=1.10E0,    &     !< No longer used/supported
      STABS=1.0E0,      &     !< No longer used/supported
      STABD=.90E0,      &     !< No longer used/supported
      STABFC=1.00E0,    &     !< No longer used/supported
      DTTOP=0.0E0,      &     !< No longer used/supported
!---VVVVV
      RHF=0.10,      &        !< No longer used/supported (Relative humidity factor)
      EPSUP=1.00,    &        !< No longer used/supported
      EPSDN=1.05,    &        !< No longer used/supported
      EPSTH=0.0,     &        !< No longer used/supported
      PBM=13000.,    &        !< _____ ?
      PQM=20000.,    &        !< No longer used/supported
      PNO=1000.,     &        !< No longer used/supported
      PONE=2500.,    &        !< _____ ?
      ZSH=2000.,     &        !< No longer used/supported
      PFRZ=15000.,   &        !< _____ ?
      PSHU=45000.,   &        !< No longer used/supported

!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
!
!    AUGUST '91: SCHEME HAVING THE OPTION OF USING DIFFERENT FAST AND
!    SLOW PROFILES FOR SEA AND FOR LAND POINTS; AND ALSO THE "SEA" AND
!    THE "LAND" SCHEME EVERYWHERE. OVER LAND PROFILES DEPART FROM THE
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
     FSS=.85E0, &             !< Factor to obtain “slow” profiles over the sea ?
     EFIMN=.20E0, &           !< Minimum convective cloud efficiency
     EFMNT=.70E0, &           !< _____ ?
     FCC=.50, &               !< _____ ?
     FCP=H1-FCC, &            !< _____ ?
!
!         IN THIS VERSION 3.5, OVER LAND AND FOR THE FAST PROFILES, DSPB
!         IS PRESCRIBED TO BE 25 PERCENT DRIER THAN THE FAST SEA VALUE
!         (IN ROUGH AGREEMENT WITH BINDER, QJ, IN PRESS) WHILE DSP0 AND
!         DSPT ARE EACH 20 PERCENT DRIER THAN THE CORRESPONDING FAST
!         SEA VALUES.                WITH FSL=.875 THIS MAKES THE
!         AVERAGE OF THE FAST AND THE SLOW LAND PROFILES SOMEWHAT DRIER
!         THAN THE OCT90 FIXED LAND PROFILES.                         FM
!
!> Deficit Saturation Pressure (DSP)
!> DSPB = cloud base
!> DSP0 = freezing level
!> DSPT = cloud top
!
     DSPBFL=-4843.75E0, &     !< Deficit Saturation Pressure (DSP) at cloud base - fast land ?
     DSP0FL=-7050.00E0, &     !< Deficit Saturation Pressure (DSP) at freezing level - fast land ?
     DSPTFL=-2250.0E0, &      !< Deficit Saturation Pressure (DSP) at cloud top - fast land ?
     FSL=.850E0,   &          !< Factor to obtain "slow" profiles over land ?
!***   ACTIVATE THE FOLLOWING LINE IF OCT90=.TRUE. (AND COMMENT OUT THE
!***   PRECEDING LINE):
!    DSPBFL=-3875.E0, &            !< Deficit Saturation Pressure (DSP) at cloud base - fast land ?
     DSP0FL=-5875.E0, &            !< Deficit Saturation Pressure (DSP) at freezing level - fast land ?
     DSPTFL=-1875.E0, &            !< Deficit Saturation Pressure (DSP) at cloud top - fast land ?
     FSL=1.0E0, &                  !< Factor to obtain "slow" profiles over land ?
     DSPBFS=-3875.E0, &            !< Deficit Saturation Pressure (DSP) at cloud base  - fast sea ?
     DSP0FS=-5875.E0, &            !< Deficit Saturation Pressure (DSP) at freezing level  - fast sea ?
     DSPTFS=-1875.E0,   &          !< Deficit Saturation Pressure (DSP) at cloud top  - fast sea ?
     DSPBSL=DSPBFL*FSL, &          !< Deficit Saturation Pressure (DSP) at cloud base - slow land ?
     DSP0SL=DSP0FL*FSL, &          !< Deficit Saturation Pressure (DSP) at freezing level - slow land ?
     DSPTSL=DSPTFL*FSL, &          !< Deficit Saturation Pressure (DSP) at cloud top - slow land ?
     DSPBSS=DSPBFS*FSS, &          !< Deficit Saturation Pressure (DSP) at cloud base - slow sea ?
     DSP0SS=DSP0FS*FSS, &          !< Deficit Saturation Pressure (DSP) at freezing level - slow sea ?
     DSPTSS=DSPTFS*FSS, &          !< Deficit Saturation Pressure (DSP) at cloud top - slow sea ?
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
     TREL=3000.,                 &     !< Relaxation time
     EPSNTP=.0010E0,             &     !< No longer used/supported
     EFIFC=5.0E0,                &     !< No longer used/supported
     AVGEFI=(EFIMN+1.E0)*.5E0,   &     !< Average cloud efficiency
     DSPC=-3000.E0,              &     !< No longer used/supported
     EPSP=1.E-7,                 &     !< No longer used/supported
     STEFI=1.E0,                 &     !< Starting efficiency (EFI)
!*** ACTIVATE THE FOLLOWING LINE AND COMMENT OUT THE PRECEDING LINE IF
!*** OCT90=.TRUE.
!    STEFI=AVGEFI,                        &   
     SLOPBL=(DSPBFL-DSPBSL)/(H1-EFIMN),   &   !< Slope - cloud base over land ?
     SLOP0L=(DSP0FL-DSP0SL)/(H1-EFIMN),   &   !< Slope - freezing level over land ?
     SLOPTL=(DSPTFL-DSPTSL)/(H1-EFIMN),   &   !< Slope - cloud top over land ?
     SLOPBS=(DSPBFS-DSPBSS)/(H1-EFIMN),   &   !< Slope - cloud base over sea ?
     SLOP0S=(DSP0FS-DSP0SS)/(H1-EFIMN),   &   !< Slope - freezing level over sea ?
     SLOPTS=(DSPTFS-DSPTSS)/(H1-EFIMN),   &   !< Slope - cloud top over sea ?
     SLOPE=(H1   -EFMNT)/(H1-EFIMN)           !< Slope ?
   real,parameter :: & 
     A23M4L=A2*(A3-A4)*ELWV,  &  !< Coefficient derived from A2, A3, A4, and ELWV
     ELOCP=ELIVW/CP,          &  !< ELIVW over CP
     CPRLG=CP/(ROW*G*ELWV),   &  !< Ratio of CP to product of water density, acceleration due to gravity, and latent heat of vaporization of water
     RCP=H1/CP                   !< Reciprocal of CP (1/CP)
   logical,parameter :: &
     UNIS=.FALSE.,      &    !< _____ ?
     UNIL=.FALSE.,      &    !< _____ ?
     OCT90=.FALSE.           !< Flag to activate old Zavisa Oct. 1990 land scheme with fixed land profiles
  end module cuparm_mod

