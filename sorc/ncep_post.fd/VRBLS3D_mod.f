!> @file
!> @brief VRBLS3D declares 3D variables that are used throughout the UPP code
!> 
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2001-10-22 | H CHUANG | MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!> 2002-04-17 | BALDWIN  | MODIFIED TO INCLUDE ALL 3D ARRAYS
!> 2011-10-18 | SARAH LU | MODIFIED TO INCLUDE AEROSOL OPTICAL PROPERTIES
!> 2011-12-15 | SARAH LU | MODIFIED TO INCLUDE AEROSOL DIAG FIELDS
!> 2012-01-06 | SARAH LU | MODIFIED TO INCLUDE AIR DENSITY AND LAYER THICKNESS
!> 2015-07-02 | SARAH LU | MODIFIED TO INCLUDE SCATTERING AEROSOL OPTICAL THICKNESS
!> 2023-03-22 | WM LEWIS | ADDED EFFECTIVE RADIUS ARRAYS
!> 2023-08-16 | Yali Mao | Add CIT (Convectively-Induced Turbulence) for GTG4
      module vrbls3d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!

      real, allocatable :: UH(:,:,:) &    !< U-component wind (U at P-points including 2 row halo)
      ,VH(:,:,:) &         !< V-component wind (V at P-points including 2 row halo)
      ,WH(:,:,:) &         !< Geometric vertical velocity
      ,U(:,:,:) &          !< U-component wind
      ,V(:,:,:) &          !< V-component wind
      ,T(:,:,:) &          !< Temperature
      ,Q(:,:,:) &          !< Specific humidity
      ,CWM(:,:,:) &        !< Total condensate mixing ratio 
      ,Q2(:,:,:) &         !< Turbulence kinetic energy
      ,PMID(:,:,:) &       !< Mid-layer pressure 
      ,PMIDV(:,:,:) &      !< Model midlayer for v-point just below the pressure level to which we are interpolating ? 
      ,PINT(:,:,:) &       !< Model layer interface pressure
      ,ALPINT(:,:,:) &     !< _____? 
      ,ZMID(:,:,:) &       !< Mid-layer height
      ,ZINT(:,:,:) &       !< Model layer interface height
      ,OMGA(:,:,:) &       !< Omega - vertical velocity
      ,T_ADJ(:,:,:) &      !< No longer used/supported; temperature adjustment factor?
      ,F_ice(:,:,:) &      !< Fraction of condensate in form of ice
      ,F_rain(:,:,:) &     !< Fraction of condensate in form of rain
      ,F_RimeF(:,:,:) &    !< Rime Factor - ratio of total ice growth to deposition growth
      ,QQW(:,:,:) &        !< cloud water mixing ratio
      ,QQI(:,:,:) &        !< ice mixing ratio
      ,QQR(:,:,:) &        !< rain mixing ratio
      ,QQS(:,:,:) &        !< snow mixing ratio
      ,QQG(:,:,:) &        !< graupel mixing ratio
      ,QQH(:,:,:) &        !< hail mixing ratio
      ,QQNW(:,:,:) &       !< cloud water number concentration 
      ,QQNI(:,:,:) &       !< ice number concentration 
      ,QQNR(:,:,:) &       !< rain number concentration 
      ,QC_BL(:,:,:) &      !< cloud water mixing ratio in PBL schemes 
      ,QRIMEF(:,:,:) &     !< rime factor * ice mixing ratio ?
      ,CFR(:,:,:) &        !< Instantaneous 3d cloud fraction
      ,DBZ(:,:,:) &        !< Radar reflectivity factor ?
      ,DBZR(:,:,:) &       !< Radar reflectivity factor from rain
      ,DBZI(:,:,:) &       !< Radar reflectivity factor from ice (all forms)
      ,DBZC(:,:,:) &       !< Radar reflectivity factor from parameterized convection
      ,TTND(:,:,:) &       !< Temperature tendency due to radiative flux convergence
      ,RSWTT(:,:,:) &      !< Temperature tendency due to shortwave radiation
      ,RLWTT(:,:,:) &      !< Temperature tendency due to longwave radiation
      ,REF_10CM(:,:,:) &   !< Reflectivity
      ,EXCH_H(:,:,:) &     !< Exchange coefficient
      ,TRAIN(:,:,:) &      !< Temperature tendency due to latent heating from grid scale
      ,TCUCN(:,:,:) &      !< Temperature tendency due to latent heating from convection
      ,EL_PBL(:,:,:)  &    !< Mixing length ?
      ,MCVG(:,:,:) &       !< Moisture convergence
      ,EXTCOF55(:,:,:) &   !< Unified extinction ext550/Aerosol optical depth
      ,NLICE(:,:,:) &      !< Time-averaged number concentration of large ice
      ,CFR_RAW(:,:,:) &    !< Cloud fraction (unprocessed)
!! Wm Lewis: added
      ,NRAIN(:,:,:) &         !< Number concentration of rain drops
      ,EFFRI(:,:,:) &         !< Thompson scheme cloud ice effective radius (for RRFS)
      ,EFFRL(:,:,:) &         !< Thompson scheme cloud water effective radius (for RRFS)
      ,EFFRS(:,:,:) &         !< Thompson scheme snow effective radius (for RRFS)
      ,radius_cloud(:,:,:) &  !< Effective radius, cloud drops
      ,radius_ice(:,:,:) &    !< Effective radius, cloud ice
      ,radius_snow(:,:,:) &   !< Effective radius, snow
! KRS Add HWRF fields     
      ,REFL_10CM(:,:,:) &     !< Reflectivity (for HWRF)
! Add GFS fields     
      ,O3(:,:,:) &            !< Ozone mixing ratio
      ,O(:,:,:) &             !< Atomic oxygen mixing ratio ?
      ,O2(:,:,:) &            !< Molecular oxygen mixing ratio ?
! Add GFS D3D fields
      ,vdifftt(:,:,:)         & !< Vertical diffusion temperature tendency ?
      ,tcucns(:,:,:)          & !< Temperature tendency due to latent heating from (shallow?) convection ?
      ,vdiffmois(:,:,:)       & !< Vertical diffusion moistening rate ?
      ,dconvmois(:,:,:)       & !< Deep convective moistening rate ?
      ,sconvmois(:,:,:)       & !< Shallow convective moistening rate ?
      ,nradtt(:,:,:)          & !< Net radiation temperature tendency ? 
      ,o3vdiff(:,:,:)         & !< Ozone vertical diffusion ?
      ,o3prod(:,:,:)          & !< Ozone production ?
      ,o3tndy(:,:,:)          & !< Ozone tendency ?
      ,mwpv(:,:,:)            & !< Mass-weighted potential vorticity ?
      ,unknown(:,:,:)         & !< _____?
      ,vdiffzacce(:,:,:)      & !< Vertical diffusion zonal acceleration ?
      ,zgdrag(:,:,:)          & !< Gravity wave drag zonal acceleration ?
      ,cnvctummixing(:,:,:)   & !< Convective zonal momentum mixing acceleration ?
      ,vdiffmacce(:,:,:)      & !< Vertical diffusion meridional acceleration ?
      ,mgdrag(:,:,:)          & !< Gravity wave drag meridional acceleration ?
      ,cnvctvmmixing(:,:,:)   & !< Convective meridional momentum mixing acceleration ?
      ,ncnvctcfrac(:,:,:)     & !< Non-convective cloud fraction (%) ?
      ,cnvctumflx(:,:,:)      & !< Convective updraft mass flux ?
      ,cnvctdmflx(:,:,:)      & !< Convective downdraft mass flux ?
      ,cnvctdetmflx(:,:,:)    & !< Convective detrainment mass flux ?
      ,cnvctzgdrag(:,:,:)     & !< Convective gravity wave drag zonal acceleration ?
      ,cnvctmgdrag(:,:,:)     & !< Convective gravity wave drag meridional acceleration ?
      ,QQNWFA(:,:,:)          & !< Water-friendly aerosol number concentration
      ,QQNIFA(:,:,:)          & !< Ice-friendly aerosol number concentration
      ,TAOD5503D(:,:,:)       & !< 3D aerosol optical depth at 550 nm
      ,AEXTC55(:,:,:)         & !< Aerosol extinction coefficient
!
! Add aerosol optical properties for GOCART (NGAC)
      ,ext(:,:,:)           & !< aerosol extinction coefficient
      ,asy(:,:,:)           & !< asymmetry parameter
      ,ssa(:,:,:)           & !< single-scattering albedo
      ,sca(:,:,:)           & !< aerosol scattering coefficient ?
! Add aerosol diagnosis fields for GOCART (NGAC)
      ,duem(:,:,:)         & !< Dust emission fluxes
      ,dusd(:,:,:)         & !< Dust sedimentation fluxes
      ,dudp(:,:,:)         & !< Dust dry deposition fluxes
      ,duwt(:,:,:)         & !< Dust wet deposition fluxes
      ,dusv(:,:,:)         & !< Dust scavenging fluxes
      ,sssv(:,:,:)         & !< Seasalt scavenging fluxes
      ,suem(:,:,:)         & !< Sulfate emission mass flux
      ,susd(:,:,:)         & !< Sulfate sedimentation mass flux
      ,sudp(:,:,:)         & !< Sulfate dry deposition mass flux
      ,suwt(:,:,:)         & !< Sulfate wet deposition mass flux
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
      ,icing_gfip(:,:,:) &    !< Global Forecast Icing Potential
      ,icing_gfis(:,:,:) &    !< Global Forecast Icing Severity
! Add NCAR GTG turbulence
      ,catedr(:,:,:) &        !< Clean air turbulence (CAT) eddy dissipation parameter (EDR)
      ,mwt(:,:,:) &           !< Mountain wave turbulence
      ,gtg(:,:,:) &           !< Graphical turbulence guidance
      ,cit(:,:,:) &           !< Convectively-induced turbulence
! AQF
      ,avgozcon(:,:,:) &      !< Average ozone concentration
      ,avgpmtf(:,:,:)         !< Average particulate matter (fine)

      end module vrbls3d
