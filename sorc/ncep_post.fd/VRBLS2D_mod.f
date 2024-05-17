!> @file
!> @brief module: VRBLS2D declares 2D variables that are used throughout the
!UPP code
      module vrbls2d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      real, allocatable :: &
      U10   (:,:) &        !< 10m u-wind component
      ,AKMS  (:,:) &       !< Surface exchange coefficient of momentum
      ,AKHS  (:,:) &       !< Surface exchange coefficient of heat/moisture
      ,THS   (:,:) &       !< Surface potential temperature
      ,QS(:,:) &           !< Specific humidity
      ,UZ0(:,:) &          !< U-wind component at roughness length (m/s)
      ,VZ0(:,:) &          !< V-wind component at roughness length (m/s)
      ,THZ0(:,:) &         !< Potential temperature at roughness length (K)
      ,QZ0(:,:) &          !< Specific humidity at roughness length (kg/kg)
      ,SNO   (:,:) &       !< Instantaneous snow water equivalent _____
      ,TSHLTR   (:,:) &    !< 2m temperature
      ,QSHLTR(:,:) &       !< 2m specific humidity
      ,MRSHLTR(:,:) &      !< Shelter mixing ratio
      ,V10(:,:) &          !< 10 m v-wind component
      ,ACPREC(:,:) &       !< Accumulated total precipitation
      ,CUPREC(:,:) &       !< Accumulated total convective precipitation
      ,ANCPRC(:,:) &       !< Accumulated total grid-scale precipitation
      ,CUPPT(:,:) &        !< Accumulated convective rain since last call to radiation
      ,SMSTAV(:,:) &       !< Soil moisture availability [%]      
      ,SSROFF(:,:) &       !< Storm surface runoff
      ,BGROFF(:,:) &       !< Underground/subsurface runoff
      ,VEGFRC(:,:) &       !< Vegetation fraction
      ,SHDMIN(:,:) &       !< Annual MIN vegetation fraction
      ,SHDMAX(:,:) &       !< Annual MAX vegetation fraction
      ,LAI(:,:) &          !< Leaf area index   
      ,ACSNOW(:,:) &       !< Accumulated snowfall
      ,ACSNOM(:,:) &       !< Accumulated snowmelt
      ,CMC(:,:) &          !< Canopy moisture content (m)
      ,SST(:,:) &          !< Sea surface temperature
      ,RSWIN(:,:) &        !< Incoming shortwave radiation at the surface
      ,RLWIN(:,:) &        !< Incoming longwave radiation at the surface
      ,RLWTOA(:,:) &       !< Outgoing longwave radiation at the top of the atmopshere 
      ,LWDNBC(:,:) &       !< Downward longwave radiation at the surface 
      ,LWUPBC(:,:) &       !< Upward longwave radiation at the surface
      ,TG(:,:) &           !< Ground temperature
      ,SFCSHX(:,:) &       !< Surface sensible heat flux
      ,PSLP(:,:) &         !< Reduced sea-level pressure array
      ,T700(:,:) &         !< Temperature at 700 hPa
      ,Z500(:,:) &         !< Geopotential height at 500 hPa
      ,Z700(:,:) &         !< Geopotential height at 700 hPa
      ,SFCLHX(:,:) &       !< Time averaged surface latent heat flux
      ,FIS(:,:) &          !< Geopotential height of the surface
      ,T500(:,:) &         !< Temperature at 500 hPa
      ,Z1000(:,:) &        !< Geopotential height at 1000 hPa
      ,SLP(:,:) &          !< Sea level pressure
      ,CFRACL(:,:) &       !< Radiation state variable - low cloud fraction
      ,CFRACM(:,:) &       !< Radiation state variable - middle cloud fraction
      ,CFRACH(:,:) &       !< Radiation state variable - high cloud fraction
      ,ACFRST(:,:) &       !< Radiation state variable - accumulated stratiform cloud fraction
      ,ACFRCV(:,:) &       !< Radiation state variable - accumulated convective cloud fraction
      ,NCFRST(:,:) &       !< Radiation state variable - times stratiform cloud >0 between rad calls
      ,NCFRCV(:,:) &       !< Radiation state variable - times convec cloud >0 between rad calls
      ,HBOT(:,:) &         !< Bottom of convection level
      ,HTOP(:,:) &         !< Top of convection level
      ,ASWIN(:,:) &        !< Time-averaged incoming shortwave radiation at the surface
      ,ALWIN(:,:) &        !< Time-averaged incoming longwave radiation at the surface
      ,ASWOUT(:,:) &       !< Time-averaged outgoing shortwave radiation at the surface
      ,ALWOUT(:,:) &       !< Time-averaged outgoing longwave radiation at the surface
      ,ASWTOA(:,:) &       !< Time-averaged outgoing shortwave radiation at the model top
      ,ALWTOA(:,:) &       !< Time-average outgoing longwave radiation at the model top
      ,CZEN(:,:) &         !< Cosine of solar zenith angle
      ,CZMEAN(:,:) &       !< Mean cosine of solar zenith angle
      ,SIGT4(:,:) &        !< Sigma of temperature (Stefan-Boltzmann * T**4)
      ,RSWOUT(:,:) &       !< Instantaneous outgoing shortwave radiation from the surface
      ,RADOT(:,:) &        !< Instantaneous outgoing longwave radiation from the surface _____ Radiative emission from surface
      ,SMSTOT(:,:) &       !< Total soil moisture
      ,PCTSNO(:,:) &       !< Snow cover percentage
      ,PSHLTR(:,:) &       !< Shelter-level pressure
      ,TH10(:,:) &         !< Potential temperature at 10 meters above the surface (anemometer level)
      ,Q10(:,:) &          !< Specific humidity at 10 meters above the surface (anemometer level)
      ,SR(:,:) &           !< Snow ratio
      ,PREC(:,:) &         !< Precipitation
      ,SUBSHX(:,:) &       !< Accumulated deep soil heat flux
      ,SNOPCX(:,:) &       !< Snow phase change heat flux
      ,SFCUVX(:,:) &       !< Total momentum flux _____ ??? Surface ultraviolet radiation 
      ,SFCEVP(:,:) &       !< Accumulated surface evaporation
      ,POTEVP(:,:) &       !< Potential evaporation
      ,Z0(:,:) &           !< Roughness length
      ,USTAR(:,:) &        !< Frictional velocity
      ,TWBS(:,:) &         !< Instantaneous surface sensible heat flux
      ,QWBS(:,:) &         !< Instantaneous surface latent heat flux
      ,SFCEXC(:,:) &       !< Surface exchange coefficient 
      ,GRNFLX(:,:) &       !< Instantaneous ground heat flux 
      ,SOILTB(:,:) &       !< Deep ground soil temperature
      ,F(:,:) &            !< ??? Coriolis-related - possibly (Coriolis * DT/2) or Coriolis sine latitude term
      ,ALBEDO(:,:) &       !< Surface albedo
      ,CLDFRA(:,:) &       !< Instantaneous 3D cloud fraction
      ,CPRATE(:,:) &       !< Convective precipitation rate
      ,CNVCFR(:,:) &       !< Convective cloud fraction
      ,PBLH(:,:) &         !< Planetary boundary layer (PBL) height
      ,PBLHGUST(:,:) &     !< Effective PBL height diagnosed from the theta-v profile, rather than the Ri profile
      ,HBOTD(:,:) &        !< Bottom of the deep convection layer
      ,HTOPD(:,:) &        !< Top of the deep convection layer
      ,HBOTS(:,:) &        !< Bottom of the shallow convection layer
      ,HTOPS(:,:) &        !< Top of the shallow convection layer
      ,CLDEFI(:,:) &       !< Convective cloud efficiency
      ,ALBASE(:,:) &       !< Base (snow-free) albedo
      ,SI(:,:) &           !< Snow depth in mm
      ,LSPA(:,:) &         !< Land surface precipitation accumulation
      ,RSWINC(:,:) &       !< Incoming clear-sky shortwave radiation at the surface (clear-sky equivalent of RSWIN)
      ,VIS(:,:) &          !< Visibility
      ,PD(:,:) &           !< Surface pressure minus PTOP
      ,MXSNAL(:,:) &       !< Maximum snow albedo
      ,MIXHT(:,:) &        !< Mixing height on surface
      ,SNONC(:,:) &        !< Accumulated total grid-scale snow/ice ? 
      ,EPSR(:,:) &         !< Radiative emissivity
      ,RSWTOA(:,:) &       !< Outgoing shortwave radiation flux at top of atmosphere
      ,TEQL(:,:) &         !< Equivalent temperature _____?
! Variables saved for input to IFI
      ,IFI_APCP(:,:) &     !< In-flight icing (IFI) total accumulated precipitation at surface
      ,CAPE(:,:) &         !< Convective available potential energy
      ,CIN(:,:) &          !< Convective inhibition
! HWRF additions
      ,MDLTAUX(:,:) &      !< Zonal (u-) component of the wind stress
      ,MDLTAUY(:,:) &      !< Meridional (v-) component of the wind stress
      ,CD10(:,:) &         !< Drag coefficient at 10 meters above the ground
      ,CH10(:,:)  &        !< Heat transfer coefficient at 10 meters above the ground
      ,ACSWUPT(:,:) &      !< Accumulated shortwave upwelling radiation flux at top
      ,SWDNT(:,:) &        !< Instantaneous shortwave downward radiation flux at top
      ,ACSWDNT(:,:) &      !< Accumulated shortwave downward radiation flux at top
! NAMB additions
      ,SNOAVG(:,:) &       !< Average snow cover
      ,PSFCAVG(:,:) &      !< Average surface pressure ?
      ,T10AVG(:,:) &       !< Time-averaged temperature at 10 meters above ground
      ,AKHSAVG(:,:) &      !< Time-averaged mass exchange coefficient or surface exchange coefficient for heat?
      ,AKMSAVG(:,:) &      !< Time-averaged wind exchange coefficient or surface exchange coefficient for momentum?
      ,T10M(:,:) &         !< Temperature at 10 meters above ground
      ,U10MAX(:,:) &       !< Maximum hourly zonal (u-) wind speed at 10 meters above ground level
      ,V10MAX(:,:) &       !< Maximum hourly meridional (v-) wind speed at 10 meters above ground level
      ,u10h(:,:) &         !< Hourly zonal (u-) wind speed at 10 meters above ground level
      ,v10h(:,:) &         !< Hourly meridional (v-) wind speed at 10 meters above ground level
      ,PRATE_MAX(:,:) &    !< Maximum precipitation rate in mm/h
      ,FPRATE_MAX(:,:) &   !< Maximum frozen precipitation rate in mm/h
! GSD addition
      ,WSPD10MAX(:,:) &       !< Maximum hourly wind speed at 10 meters above ground level
      ,W_UP_MAX(:,:) &        !< Maximum hourly updraft velocity
      ,W_DN_MAX(:,:) &        !< Maximum hourly downdraft velocity
      ,REFD_MAX(:,:) &        !< Maximum hourly 1 km above ground level reflectivity
      ,UP_HELI_MAX(:,:) &     !< Maximum hourly updraft helicity
      ,UP_HELI_MAX16(:,:) &   !< Maximum 1-6km hourly updraft helicity
      ,GRPL_MAX(:,:) &        !< Maximum column-integrated graupel
      ,QRMAX(:,:) &           !< Maximum rain water mixing ratio
      ,UP_HELI(:,:) &         !< Updraft helicity
      ,UP_HELI16(:,:) &       !< Updraft helicity in 1-6km layer
      ,LTG1_MAX(:,:) &        !< Maximum lightning threat index 1
      ,LTG2_MAX(:,:) &        !< Maximum lightning threat index 2
      ,LTG3_MAX(:,:) &        !< Maximum lightning threat index 3
      ,UP_HELI_MIN(:,:) &     !< Minimum updraft helicity
      ,UP_HELI_MIN16(:,:) &   !< Minimum updraft helicity over 1-6 km layer
      ,UP_HELI_MAX02(:,:) &   !< Maximum updraft helicity over 0-2 km layer
      ,UP_HELI_MIN02(:,:) &   !< Minimum updraft helicity over 0-2 km layer
      ,UP_HELI_MAX03(:,:) &   !< Maximum updraft helicity over 0-3 km layer
      ,UP_HELI_MIN03(:,:) &   !< Minimum updraft helicity over 0-3 km layer
      ,REL_VORT_MAX(:,:) &    !< Maximum relative vertical vorticity
      ,REL_VORT_MAX01(:,:) &  !< Maximum relative vertical vorticity at 0-1 km
      ,REL_VORT_MAXHY1(:,:) & !< Maximum relative vertical vorticity at hybrid level 1
      ,WSPD10UMAX(:,:) &      !< Maximum u-component of wind at 10 meters above ground level
      ,WSPD10VMAX(:,:) &      !< Maximum hourly v-component of wind at 10 meters above ground level
      ,REFDM10C_MAX(:,:) &    !< Maximum hourly -10C reflectivity
      ,HAIL_MAX2D(:,:) &      !< Maximum hail diameter in column
      ,HAIL_MAXK1(:,:) &      !< Maximum hail diameter at k=1
      ,HAIL_MAXHAILCAST(:,:) & !< Maximum hail diameter at surface from HAILCAST algorithm (HRRR and RRFS applications)
      ,NCI_LTG(:,:) &         !< Convective initiation lightning
      ,NCA_LTG(:,:) &         !< Convective activity lightning
      ,NCI_WQ(:,:) &          !< Convective Initiation Vertical Hydrometeor Flux
      ,NCA_WQ(:,:) &          !< Convective Activity Vertical Hydrometeor Flux
      ,NCI_REFD(:,:) &        !< Convective Initiation Reflectivity
      ,NCA_REFD(:,:) &        !< Convective Activity Reflectivity
      ,RAINC_BUCKET1(:,:) &   !< Accumulated cumulus precipitation over BUCKET_DT1 periods of time
      ,RAINNC_BUCKET1(:,:) &  !< Accumulated grid-scale precipitation over BUCKET_DT1 periods of time
      ,RAINC_BUCKET(:,:) &    !< Accumulated cumulus precipitation over BUCKET_DT periods of time
      ,RAINNC_BUCKET(:,:) &   !< Accumulated grid-scale precipitation over BUCKET_DT periods of time
      ,SNOW_BUCKET(:,:) &     !< Accumulated grid-scale snow over BUCKET_DT periods of time
      ,GRAUP_BUCKET(:,:) &    !< Accumulated grid-scale graupel over BUCKET_DT periods of time
      ,PCP_BUCKET(:,:) &      !< Bucket total precipitation over BUCKET_DT periods of time
      ,ACGRAUP(:,:) &         !< Accumulated graupel/sleet
      ,ACFRAIN(:,:) &         !< Accumulated freezing rain
      ,FRZRN_BUCKET(:,:) &    !< Freezing rain bucket
      ,SNOW_ACM(:,:) &        !< Accumulated snowfall
      ,SNOW_BKT(:,:) &        !< Snowfall bucket
      ,SNOW_BUCKET1(:,:) &    !< Accumulated grid-scale snow over BUCKET_DT1 periods of time
      ,GRAUP_BUCKET1(:,:) &   !< Accumulated grid-scale graupel over BUCKET_DT1 periods of time
      ,PCP_BUCKET1(:,:) &     !< Bucket total precipitation over BUCKET_DT1 periods of time
      ,SNOWNC(:,:) &          !< Accumulated total grid-scale snow/ice precipitation (per time step?)
      ,GRAUPELNC(:,:) &       !< Time step non-convective graupel in [m]
      ,TMAX(:,:) &            !< Maximum 2m temperature
      ,W_MEAN(:,:) &          !< Mean vertical velocity
      ,TSNOW(:,:) &           !< Snow temperature
      ,QVG(:,:) &             !< Water vapor mixing ratio at the surface
      ,QV2m(:,:) &            !< Water vapor mixing ratio at 2 meters above ground level
      ,QVl1(:,:) &            !< Water vapor mixing ratio at level 1
      ,REFC_10CM(:,:) &       !< Composite (10cm) radar reflectivity
      ,REF1KM_10CM(:,:) &     !< (10cm) radar reflectivity at 1 km
      ,REF4KM_10CM(:,:) &     !< (10cm) radar reflectivity at 4 km
      ,SWRADmean(:,:) &       !< Time-averaged incoming shortwave radiation
      ,U10mean(:,:) &         !< Time-averaged u-component wind speed at 10 meters above ground level
      ,V10mean(:,:) &         !< Time-averaged v-component wind speed at 10 meters above ground level
      ,SPDUV10mean(:,:) &     !< Time-averaged wind speed at 10 meters above ground level
      ,SWNORMmean(:,:) &      !< Time-averaged SWNORM (terrain-normal downwelling shortwave radiation)
      ,SNFDEN(:,:) &          !< Snowfall density
      ,SNDEPAC(:,:) &         !< Accumulated depth of snowfall
      ,SWDDNI(:,:) &          !< Instantaneous shortwave surface downward direct normal irradiance
      ,SWDDIF(:,:) &          !< Instantaneous shortwave surface downward diffuse irradiance
      ,SWDNBC(:,:) &          !< Shortwave surface downward clear-sky shortwave irradiance
      ,SWDDNIC(:,:) &         !< Clear-sky shortwave surface downward direct normal irradiance
      ,SWDDIFC(:,:) &         !< Clear-sky shortwave surface downward diffuse horizontal irradiance
      ,SWUPBC(:,:) &          !< Clear-sky surface upwelling shortwave flux
      ,SWUPT(:,:) &           !< Upward shortwave flux at top of atmosphere
      ,TAOD5502D(:,:) &       !< Total aerosol optical depth at 550 nm
      ,AERASY2D(:,:) &        !< Aerosol asymmetry parameter
      ,AERSSA2D(:,:) &        !< Aerosol single-scattering albedo
      ,MEAN_FRP(:,:) &        !< Instantaneous mean fire radiative power
      ,HWP(:,:) &             !< Hourly wildfire potential
      ,LWP(:,:) &             !< Liquid water path
      ,IWP(:,:) &             !< Ice water path
      ,XLAIXY(:,:) &          !< Leaf area index ?
      ,SMOKE_AVE(:,:) &       !< Hourly averaged smoke
      ,DUST_AVE(:,:)   &      !< Hourly averaged fine dust (PM 2.5)
      ,COARSEPM_AVE(:,:) &    !< Hourly averaged coarse dust (PM 10)
! add new fields for GFS
      ,SFCUX(:,:) &           !< Time-averaged zonal momentum flux
      ,SFCVX(:,:) &           !< Time-averaged meridional momentum flux
      ,SFCUXI(:,:) &          !< Instantaneous zonal momentum flux
      ,SFCVXI(:,:) &          !< Instantaneous meridional momentum flux
      ,AVGALBEDO(:,:) &       !< Mid-day average albedo
      ,AVGCPRATE(:,:) &       !< Convective precipitation in m per physics time step
      ,AVGPREC(:,:) &         !< Average precipitation rate in m per physics time step
      ,PTOP(:,:) &            !< Instantaneous convective cloud top pressure
      ,PBOT(:,:) &            !< Instantaneous convective cloud bottom pressure
      ,AVGCFRACH(:,:) &       !< Average high cloud fraction
      ,AVGCFRACM(:,:) &       !< Average mid cloud fraction
      ,AVGCFRACL(:,:) &       !< Average low cloud fraction
      ,AVGTCDC(:,:) &         !< Time-averaged column cloud fraction
      ,AUVBIN(:,:) &          !< Time-averaged incoming surface UV-B
      ,AUVBINC(:,:) &         !< Time-averaged incoming surface clear-sky UV-B
      ,ptopl(:,:) &           !< Time-averaged low cloud top pressure
      ,pbotl(:,:) &           !< Time-averaged low cloud bottom pressure
      ,Ttopl(:,:) &           !< Time-averaged low cloud top temperature 
      ,ptopm(:,:) &           !< Time-averaged middle cloud top pressure
      ,pbotm(:,:) &           !< Time-averaged middle cloud bottom pressure 
      ,Ttopm(:,:) &           !< Time-averaged middle cloud top temperature
      ,ptoph(:,:) &           !< Time-averaged high cloud top pressure
      ,pboth(:,:) &           !< Time-averaged high cloud bottom pressure 
      ,Ttoph(:,:) &           !< Time-averaged high cloud top temperature
      ,sfcugs(:,:) &          !< _____?
      ,sfcvgs(:,:) &          !< _____?
      ,PBLCFR(:,:) &          !< Boundary layer cloud cover
      ,cldwork(:,:) &         !< Cloud work function
      ,gtaux(:,:) &           !< Time-averaged zonal gravity wave stress
      ,gtauy(:,:) &           !< Time-averaged meridional gravity wave stress
      ,runoff(:,:) &          !< Accumulated total (base+surface) runoff
      ,maxtshltr(:,:) &       !< Shelter max temperature
      ,mintshltr(:,:) &       !< Shelter min temperature
      ,maxrhshltr(:,:) &      !< Shelter max relative humidity
      ,minrhshltr(:,:) &      !< Shelter min relative humidity
      ,dzice(:,:) &           !< Ice thickness
      ,maxqshltr(:,:) &       !< Shelter max specific humidity
      ,minqshltr(:,:) &       !< Shelter min specific humidity
      ,alwinc(:,:) &          !< Time-averaged surface clear-sky incoming longwave
      ,alwoutc(:,:) &         !< Time-averaged surface clear-sky outgoing longwave
      ,alwtoac(:,:) &         !< Time-averaged clear-sky outgoing longwave at top of atmosphere 
      ,aswinc(:,:) &          !< Time-averaged surface clear-sky incoming shortwave
      ,aswoutc(:,:) &         !< Time-averaged surface clear-sky outgoing shortwave
      ,aswtoac(:,:) &         !< Time-averaged clear-sky outgoing shortwave at top of atmosphere 
      ,aswintoa(:,:) &        !< Time-averaged model top incoming shortwave
      ,smcwlt(:,:) &          !< Wilting point
      ,suntime(:,:) &         !< Sunshine duration
      ,fieldcapa(:,:) &       !< Field capacity 
      ,avisbeamswin(:,:) &    !< Time-averaged surface visible beam downward solar flux
      ,avisdiffswin(:,:) &    !< Time-averaged surface visible diffuse downward solar flux
      ,airbeamswin(:,:) &     !< Time-averaged surface near ir beam downward solar flux
      ,airdiffswin(:,:) &     !< Time-averaged surface near ir diffuse downward solar flux
      ,snowfall(:,:) &        !< Total accumulated snowfall ?
      ,acond(:,:) &           !< Aerodynamic conductance on surface
      ,edir(:,:) &            !< Direct soil evaporation (W/m2)
      ,ecan(:,:) &            !< Canopy water evaporation (W/m2)
      ,etrans(:,:) &          !< Plant transpiration (W/m2)
      ,esnow(:,:) &           !< Snow sublimation (W/m2)
      ,avgedir(:,:) &         !< Direct soil evaporation (6-hr average?)
      ,avgecan(:,:) &         !< Accumulated evaporation of intercepted water (6-hr average?)
      ,avgetrans(:,:)&        !< Plant transpiration (6-hr average?)
      ,avgesnow(:,:) &        !< Snow sublimation (6-hr average?)
      ,avgpotevp(:,:) &       !< Time-averaged accumulated potential evaporation
      ,avgprec_cont(:,:) &    !< Average precipitation rate - continuous bucket
      ,avgcprate_cont(:,:) &  !< Convective precipitation - coninuous bucket
      ,ti(:,:) &              !< Sea ice skin temperature
      ,aod550(:,:) &          !< Instantaneous aerosol optical depth at 550 nm
      ,du_aod550(:,:) &       !< Instantaneous aerosol optical depth at 550 nm (dust)
      ,ss_aod550(:,:) &       !< Instantaneous aerosol optical depth at 550 nm (seasalt)
      ,su_aod550(:,:) &       !< Instantaneous aerosol optical depth at 550 nm (sulfates)
      ,bc_aod550(:,:) &       !< Instantaneous aerosol optical depth at 550 nm (organic carbon)
      ,oc_aod550(:,:) &       !< Instantaneous aerosol optical depth at 550 nm (black carbon)
      ,landfrac(:,:) &        !< Land fraction
      ,paha(:,:) &            !< Averaged precipitation advected heat flux
      ,pahi(:,:) &            !< Instantaneous precipitation advected heat flux
      ,tecan(:,:) &           !< Accumulated evaporation of intercepted water
      ,tetran(:,:) &          !< Accumulated plant transpiration
      ,tedir(:,:) &           !< Accumulated soil surface evaporation
      ,twa(:,:) &             !< Total water storage in aquifer
      ,fdnsst(:,:) &          !< Foundation temperature
      ,pwat(:,:)              !< Precipitable water
      integer, allocatable :: IVGTYP(:,:) &     !< Vegetation type
      ,ISLTYP(:,:) &          !< Soil type
      ,ISLOPE(:,:) &          !< Slope type
      ,IEQL(:,:)              !< EQ level (highest positively buoyant level ?)
      
! Add 2d aerosol diagnosis fields for GOCART (NEMS-GFS Aerosol Component [NGAC])
      real, allocatable ::                                                   &
       DUSMASS(:,:) &         !< Dust (PM10) surface mass concentration
      ,DUCMASS(:,:) &         !< Dust (PM10) column mass density
      ,DUSMASS25(:,:) &       !< Dust (PM25) surface mass concentration
      ,DUCMASS25(:,:) &       !< Dust (PM25) column mass density
      ,SUSMASS(:,:) &         !< Sulfate surface mass concentration - no longer used
      ,SUCMASS(:,:) &         !< Sulfate column mass density - no longer used
      ,SUSMASS25(:,:) &       !< Sulfate (PM25) surface mass concentration - no longer used
      ,SUCMASS25(:,:) &       !< Sulfate (PM25) column mass density - no longer used
      ,OCSMASS(:,:) &         !< Organic Carbon Surface Mass Concentration - no longer used
      ,OCCMASS(:,:) &         !< Organic Carbon Column Mass Density - no longer used
      ,OCSMASS25(:,:) &       !< Organic Carbon (PM25) Surface Mass Concentration - no longer used
      ,OCCMASS25(:,:) &       !< Organic Carbon (PM25) Column Mass Density - no longer used
      ,BCSMASS(:,:) &         !< Black Carbon Surface Mass Concentration - no longer used
      ,BCCMASS(:,:) &         !< Black Carbon Column Mass Density - no longer used
      ,BCSMASS25(:,:) &       !< Black Carbon (PM25) Surface Mass Concentration - no longer used
      ,BCCMASS25(:,:) &       !< Black Carbon (PM25) Column Mass Density - no longer used
      ,SSSMASS(:,:) &         !< Sea Salt Surface Mass Concentration - no longer used
      ,SSCMASS(:,:) &         !< Sea Salt Column Mass Density - no longer used
      ,SSSMASS25(:,:) &       !< Sea Salt (PM25) Surface Mass Concentration - no longer used
      ,SSCMASS25(:,:) &       !< Sea Salt (PM25) Column Mass Density - no longer used
      ,DUSTCB(:,:) &          !< GFS output dust in nemsio (GOCART)
      ,SSCB(:,:) &            !< GFS output sea salt in nemsio (GOCART)
      ,OCCB(:,:) &            !< GFS output organic carbon in nemsio (GOCART)
      ,BCCB(:,:) &            !< GFS output black carbon in nemsio (GOCART)
      ,SULFCB(:,:) &          !< GFS output sulfate in netcdf (GOCART)
      ,DUSTALLCB(:,:) &       !< GFS output dust in nemsio (GOCART)
      ,SSALLCB(:,:) &         !< GFS output sea salt in nemsio (GOCART)
      ,DUSTPM(:,:) &          !< PM25 dust
      ,SSPM(:,:) &            !< PM25 sea salt
      ,PP25CB(:,:) &          !< GFS output pp25 in nemsio (GOCART)
      ,DUSTPM10(:,:) &        !< Dust 10 micrometers mass density concentration on model surface
      ,PP10CB(:,:) &          !< GFS output pp10 in nemsio (GOCART)
      ,NO3CB(:,:) &           !< GFS output nitrate in netcdf (GOCART)
      ,NH4CB(:,:) &           !< GFS output NH4 in netcdf (GOCART)
      ,maod(:,:)              !< MIE AOD at 550nm !for FV3-Chem ___?

! Add new field for AQM
      real, allocatable :: aqm_aod550(:,:)   !< Air quality model (AQM) aerosol optical depth at 550nm
 
!
      end module vrbls2d
