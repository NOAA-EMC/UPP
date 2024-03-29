!> @file
!> @brief module: VRBLS3D declares 3D variables that are used throughout
!the UPP code
!   01-10-22  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-04-17  BALDWIN  - MODIFIED TO INCLUDE ALL 3D ARRAYS
!   11-10-18  SARAH LU - MODIFIED TO INCLUDE AEROSOL OPTICAL PROPERTIES
!   11-12-15  SARAH LU - MODIFIED TO INCLUDE AEROSOL DIAG FIELDS
!   12-01-06  SARAH LU - MODIFIED TO INCLUDE AIR DENSITY AND LAYER THICKNESS
!   15-07-02  SARAH LU - MODIFIED TO INCLUDE SCATTERING AEROSOL OPTICAL THICKNESS
!   23-03-22  WM LEWIS - ADDED EFFECTIVE RADIUS ARRAYS
!   23-08-16  Yali Mao - Add CIT (Convectively-Induced Turbulence) for GTG4
      module vrbls3d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!

      real, allocatable :: UH(:,:,:) &    !< u-component wind
      ,VH(:,:,:) &         !< v-component wind
      ,WH(:,:,:) &         !< geometric vertical velocity
      ,U(:,:,:) &          !< _____
      ,V(:,:,:) &          !< _____
      ,T(:,:,:) &          !< temperature
      ,Q(:,:,:) &          !< specific humidity
      ,CWM(:,:,:) &        !< total condensate mixing ratio 
      ,Q2(:,:,:) &         !< turbulence kinetic energy _____???
      ,PMID(:,:,:) &       !< pressure _____ input pressure? mid layer pressure? PRESSURE IN LAYERS? 
      ,PMIDV(:,:,:) &      !< _____ 
      ,PINT(:,:,:) &       !< pressure on interfaces 
      ,ALPINT(:,:,:) &     !< _____ 
      ,ZMID(:,:,:) &       !< _____ 
      ,ZINT(:,:,:) &       !< ETA INTERFACES HEIGHT FIELD _____?
      ,OMGA(:,:,:) &       !< Omega - vertical velocity
      ,T_ADJ(:,:,:) &      !< _____ 
      ,F_ice(:,:,:) &      !< fraction of ice
      ,F_rain(:,:,:) &     !< fraction of rain
      ,F_RimeF(:,:,:) &    !< mass ratio of rimed ice
      ,QQW(:,:,:) &        !< cloud water mixing ratio
      ,QQI(:,:,:) &        !< ice mixing ratio
      ,QQR(:,:,:) &        !< rain mixing ratio
      ,QQS(:,:,:) &        !< snow mixing ratio
      ,QQG(:,:,:) &        !< graupel mixing ratio
      ,QQH(:,:,:) &        !< hail mixing ratio
      ,QQNW(:,:,:) &       !< cloud water number concentration 
      ,QQNI(:,:,:) &       !< ice number concentration 
      ,QQNR(:,:,:) &       !< rain number concentration 
      ,QC_BL(:,:,:) &      !< cloud water mixing ratio in PBL schemes _____ 
      ,QRIMEF(:,:,:) &     !< rime factor * ice mixing ratio _____
      ,CFR(:,:,:) &        !< Instantaneous 3d cloud fraction
      ,DBZ(:,:,:) &        !< Equivalent radar reflectivity factor _____
      ,DBZR(:,:,:) &       !< Equivalent radar reflectivity factor from rain
      ,DBZI(:,:,:) &       !< Equivalent radar reflectivity factor from ice (all forms)
      ,DBZC(:,:,:) &       !< Equivalent radar reflectivity factor from parameterized convection
      ,TTND(:,:,:) &       !< Temperature tendency due to radiative flux convergence _____
      ,RSWTT(:,:,:) &      !< Temperature tendency due to shortwave radiation
      ,RLWTT(:,:,:) &      !< Temperature tendency due to longwave radiation
      ,REF_10CM(:,:,:) &   !< Reflectivity
      ,EXCH_H(:,:,:) &     !< Exchange coefficient
      ,TRAIN(:,:,:) &      !< Temperature tendency due to latent heating from grid scale_____
      ,TCUCN(:,:,:) &      !< Temperature tendency due to latent heating from convection
      ,EL_PBL(:,:,:)  &    !< Mixing length _____
      ,MCVG(:,:,:) &       !< Moisture convergence
      ,EXTCOF55(:,:,:) &   !< Unified extinction ext550/Aerosol optical depth
      ,NLICE(:,:,:) &      !< Time-averaged number concentration of large ice
      ,CFR_RAW(:,:,:) &    !< Raw cloud fraction _____???
!! Wm Lewis: added
      ,NRAIN(:,:,:) &         !< Number concentration of rain drops
      ,EFFRI(:,:,:) &         !< Thompson scheme cloud ice effective radius
      ,EFFRL(:,:,:) &         !< Thompson scheme cloud water effective radius
      ,EFFRS(:,:,:) &         !< Thompson scheme snow effective radius
      ,radius_cloud(:,:,:) &  !< _____
      ,radius_ice(:,:,:) &    !< _____
      ,radius_snow(:,:,:) &   !< _____
! KRS Add HWRF fields     
      ,REFL_10CM(:,:,:) &     !< Reflectivity
! Add GFS fields     
      ,O3(:,:,:) &            !< Ozone mixing ratio
      ,O(:,:,:) &             !< _____
      ,O2(:,:,:) &            !< _____
! Add GFS D3D fields
      ,vdifftt(:,:,:)         & !< Vertical diffusion temperature tendency_____
      ,tcucns(:,:,:)          & !< _____
      ,vdiffmois(:,:,:)       & !< Vertical diffusion moisture_____
      ,dconvmois(:,:,:)       & !< _____
      ,sconvmois(:,:,:)       & !< _____
      ,nradtt(:,:,:)          & !< _____  
      ,o3vdiff(:,:,:)         & !< _____
      ,o3prod(:,:,:)          & !< _____
      ,o3tndy(:,:,:)          & !< _____
      ,mwpv(:,:,:)            & !< _____
      ,unknown(:,:,:)         & !< _____
      ,vdiffzacce(:,:,:)      & !< _____
      ,zgdrag(:,:,:)          & !< _____
      ,cnvctummixing(:,:,:)   & !< _____
      ,vdiffmacce(:,:,:)      & !< _____
      ,mgdrag(:,:,:)          & !< _____
      ,cnvctvmmixing(:,:,:)   & !< _____
      ,ncnvctcfrac(:,:,:)     & !< _____
      ,cnvctumflx(:,:,:)      & !< _____
      ,cnvctdmflx(:,:,:)      & !< _____  
      ,cnvctdetmflx(:,:,:)    & !< _____
      ,cnvctzgdrag(:,:,:)     & !< _____
      ,cnvctmgdrag(:,:,:)     & !< _____   
      ,QQNWFA(:,:,:)          & !< _____
      ,QQNIFA(:,:,:)          & !< _____
      ,TAOD5503D(:,:,:)       & !< _____
      ,AEXTC55(:,:,:)         & !< _____
!
! Add aerosol optical properties for GOCART (NGAC)
      ,ext(:,:,:)           & !< aerosol extinction coefficient _____?
      ,asy(:,:,:)           & !< asymmetry parameter _____?
      ,ssa(:,:,:)           & !< single-scattering albedo
      ,sca(:,:,:)           & !< aerosol scattering coefficient _____?
! Add aerosol diagnosis fields for GOCART (NGAC)
      ,duem(:,:,:)         & !< Dust emission fluxes
      ,dusd(:,:,:)         & !< Dust sedimentation fluxes
      ,dudp(:,:,:)         & !< Dust dry deposition fluxes _____ ??? In Unified model variables table, it's dupd instead of dudp... bug?
      ,duwt(:,:,:)         & !< Dust wet deposition fluxes
      ,dusv(:,:,:)         & !< Dust scavenging fluxes
      ,sssv(:,:,:)         & !< Seasalt scavenging fluxes
      ,suem(:,:,:)         & !< _____
      ,susd(:,:,:)         & !< _____
      ,sudp(:,:,:)         & !< _____
      ,suwt(:,:,:)         & !< _____
      ,ssem(:,:,:)         & !< Seasalt emission fluxes
      ,sssd(:,:,:)         & !< Seasalt emission/sedimentation
      ,ssdp(:,:,:)         & !< Seasalt dry deposition fluxes
      ,sswt(:,:,:)         & !< Seasalt wet deposition fluxes
      ,ocem(:,:,:)         & !< Organic carbon emission fluxes
      ,ocsd(:,:,:)         & !< Organic carbon sedimentation fluxes
      ,ocdp(:,:,:)         & !< Organic carbon dry deposition fluxes
      ,ocwt(:,:,:)         & !< Organic carbon large wet deposition fluxes
      ,ocsv(:,:,:)         & !< Organic carbon convective wet deposition fluxes
      ,bcsv(:,:,:)         & !< Black carbon convective wet deposition fluxes
      ,bcem(:,:,:)         & !< Black carbon emission fluxes
      ,bcsd(:,:,:)         & !< Black carbon sedimentation fluxes
      ,bcdp(:,:,:)         & !< Black carbon dry deposition fluxes
      ,bcwt(:,:,:)         & !< Black carbon large wet deposition fluxes
! Add air density and thickness for GOCART (NGAC)
      ,dpres(:,:,:) &         !< Layer thickness in pressure on hybrid levels
      ,rhomid(:,:,:) &        !< Air density rhomid

! Add NCAR GFIP ICING
      ,icing_gfip(:,:,:) &    !< _____
      ,icing_gfis(:,:,:) &    !< _____
! Add NCAR GTG turbulence
      ,catedr(:,:,:) &        !< Clean air turbulence (CAT) eddy dissipation parameter (EDR)
      ,mwt(:,:,:) &           !< Mountain wave turbulence
      ,gtg(:,:,:) &           !< Graphical turbulence guidance
      ,cit(:,:,:) &           !< Convectively-induced turbulence
! AQF
      ,avgozcon(:,:,:) &      !< Average ozone concentration
      ,avgpmtf(:,:,:)         !< Average particulate matter (fine)

      end module vrbls3d
