!> @file
!> @brief module: VRBLS2D declares 2D variables that are used throughout the
!UPP code
      module vrbls2d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      real, allocatable :: &
      U10   (:,:) &        !< Max hourly u-component of 10m agl wind
      ,AKMS  (:,:) &       !< Time-averaged wind exchange coefficient
      ,AKHS  (:,:) &       !< Time-averaged mass exchange coefficient
      ,THS   (:,:) &       !< Surface potential temperature
      ,QS(:,:) &           !< Specific humidity 
      ,UZ0(:,:) &          !< Wind speed at the surface ???
      ,VZ0(:,:) &          !< Wind direction at the surface ???
      ,THZ0(:,:) &         !< Surface potential temperature ???
      ,QZ0(:,:) &          !< Specific humidity at the surface ???
      ,SNO   (:,:) &       !< Snow depth
      ,TSHLTR   (:,:) &    !< 2m temperature
      ,QSHLTR(:,:) &       !< 2m specific humidity
      ,MRSHLTR(:,:) &      !< ___? something to do with moisture?
      ,V10(:,:) &          !< 10 m v-wind component
      ,ACPREC(:,:) &       !< Accumulated precipitation over a certain period
      ,CUPREC(:,:) &       !< Cumulative precipitation over a certain period
      ,ANCPRC(:,:) &       !< A form of accumulated or cumulative precipitation ___?
      ,CUPPT(:,:) &        !< Cumulative precipitation ___?
      ,SMSTAV(:,:) &       !< Average soil moisture content       
      ,SSROFF(:,:) &       !< Surface water runoff ___?
      ,BGROFF(:,:) &       !< Groundwater runoff___?
      ,VEGFRC(:,:) &       !< Vegetation fraction___?
      ,SHDMIN(:,:) &       !< Minimum soil moisture content
      ,SHDMAX(:,:) &       !< Maximum soil moisture content
      ,LAI(:,:) &          !< Leaf area index   
      ,ACSNOW(:,:) &       !< Accumulated snowfall over a certain period
      ,ACSNOM(:,:) &       !< Accumulated snowmelt over a certain period
      ,CMC(:,:) &          !< Convective cloud cover____?
      ,SST(:,:) &          !< Sea surface temperature
      ,RSWIN(:,:) &        !< Incoming shortwave radiation at the surface
      ,RLWIN(:,:) &        !< Incoming longwave radiation at the surface
      ,RLWTOA(:,:) &       !< Outgoing longwave radiation at the top of the atmopshere 
      ,LWDNBC(:,:) &       !< Downward longwave radiation at the bottom of the atmosphere 
      ,LWUPBC(:,:) &       !< Upward longwave radiation at the bottom of the atmosphere
      ,TG(:,:) &           !< Ground temperature
      ,SFCSHX(:,:) &       !< Surface sensible heat flux
      ,PSLP(:,:) &         !< Sea-level pressure ___?
      ,T700(:,:) &         !< Temperature at 700 hPa
      ,Z500(:,:) &         !< Geopotential height at 500 hPa
      ,Z700(:,:) &         !< Geopotential height at 700 hPa
      ,SFCLHX(:,:) &       !< Time averaged surface latent heat flux
      ,FIS(:,:) &          !< ___?
      ,T500(:,:) &         !< Temperature at 500 hPa
      ,Z1000(:,:) &        !< Geopotential height at 1000 hPa
      ,SLP(:,:) &          !< Sea level pressure
      ,CFRACL(:,:) &       !< Low cloud fractions ___?
      ,CFRACM(:,:) &       !< Medium cloud fractions___?
      ,CFRACH(:,:) &       !< High cloud fractions___?
      ,ACFRST(:,:) &       !< Convective Stratiform cloud fraction ___?
      ,ACFRCV(:,:) &       !< Convective cloud fraction ___?
      ,NCFRST(:,:) &       !< None-convective stratiform cloud fraction___?
      ,NCFRCV(:,:) &       !< ___?
      ,HBOT(:,:) &         !< ___?
      ,HTOP(:,:) &         !< ___?
      ,ASWIN(:,:) &        !< Time averaged incoming surface shortwave
      ,ALWIN(:,:) &        !< Time averaged incoming surface longwave
      ,ASWOUT(:,:) &       !< Time averaged outgoing sfc shortwave
      ,ALWOUT(:,:) &       !< Time averaged outgoing surface longwave
      ,ASWTOA(:,:) &       !< Time averaged model top outgoing shortwave
      ,ALWTOA(:,:) &       !< Time average outgoing model top longwave
      ,CZEN(:,:) &         !< Zenithal cloud cover
      ,CZMEAN(:,:) &       !< Mean cloud cover at the zenith
      ,SIGT4(:,:) &        !< significant tornado paramter___?
      ,RSWOUT(:,:) &       !< Inst outgoing sfc shortwave
      ,RADOT(:,:) &        !< Instataneous outgoing surface longwave
      ,SMSTOT(:,:) &       !< Total soil moisture
      ,PCTSNO(:,:) &       !< Snow percentage cover
      ,PSHLTR(:,:) &       !< Surface potential temperature
      ,TH10(:,:) &         !< Temperature at 10 meters above the surface 
      ,Q10(:,:) &          !< Specific humidity at 10 meters above the surface
      ,SR(:,:) &           !< Surface runoff
      ,PREC(:,:) &         !< Precipitation
      ,SUBSHX(:,:) &       !< Sublimation ___?
      ,SNOPCX(:,:) &       !< Snow phase change
      ,SFCUVX(:,:) &       !< surface ultraviolet radiation exposure indices
      ,SFCEVP(:,:) &       !< Surface evaporation
      ,POTEVP(:,:) &       !< Potential evaporation
      ,Z0(:,:) &           !< Roughness length
      ,USTAR(:,:) &        !< Friction velocity
      ,TWBS(:,:) &         !< Wet bulb temperature
      ,QWBS(:,:) &         !< Wet bulb specific humidity
      ,SFCEXC(:,:) &       !< Surface exchange coefficient 
      ,GRNFLX(:,:) &       !< Ground heat flux 
      ,SOILTB(:,:) &       !< Soil temperature
      ,F(:,:) &            !< ___?
      ,ALBEDO(:,:) &       !< Surface albedo
      ,CLDFRA(:,:) &       !< Cloud fraction
      ,CPRATE(:,:) &       !< Convective precipitation rate
      ,CNVCFR(:,:) &       !< Convective cloud fraction
      ,PBLH(:,:) &         !< Planetary boundary layer height
      ,PBLHGUST(:,:) &     !< Planetary boundary layer height gust
      ,HBOTD(:,:) &        !< Depth of the bottom layer___?
      ,HTOPD(:,:) &        !< Depth of the top layer
      ,HBOTS(:,:) &        !< Surface-based layer height___?
      ,HTOPS(:,:) &        !< Surface-based layer height___?
      ,CLDEFI(:,:) &       !< Cloud ice fraction
      ,ALBASE(:,:) &       !< Base altitude of the aerosol layer
      ,SI(:,:) &           !< Snow ice fraction
      ,LSPA(:,:) &         !< Land surface pressure anomaly
      ,RSWINC(:,:) &       !< Incoming shortwave radiation at the surface____?
      ,VIS(:,:) &          !< Visibility
      ,PD(:,:) &           !< Pressure at the surface
      ,MXSNAL(:,:) &       !< Maximum snow albedo
      ,MIXHT(:,:) &        !< Mixed layer height 
      ,SNONC(:,:) &        !< Snow number concentration 
      ,EPSR(:,:) &         !< Error in shortwave radiation___?
      ,RSWTOA(:,:) &       !< Outgoing shortwave radiation at the top of the atmosphere___?
      ,TEQL(:,:) &         !< Equivalent temperature___?
! Variables saved for input to IFI
      ,IFI_APCP(:,:) &     !< Instantaneous forecast interval accumulated precipitation
      ,CAPE(:,:) &         !< Convective available potential energy
      ,CIN(:,:) &          !< Convective inhibition
! HWRF additions
      ,MDLTAUX(:,:) &      !< Zonal component of the wind stress
      ,MDLTAUY(:,:) &      !< Meridional component of the wind stress
      ,CD10(:,:) &         !< Drag coefficient at 10 meters height
      ,CH10(:,:)  &        !< Heat transfer coefficient at 10 meters height
      ,ACSWUPT(:,:) &      !< Accumulated shortwave upward radiation flux at the top of the atmosphere___?
      ,SWDNT(:,:) &        !< Shortwave downward radiation flux at the surface___?
      ,ACSWDNT(:,:) &      !< Accumulated shortwave downward radiation flux at the surface___?
! NAMB additions
      ,SNOAVG(:,:) &       !< Average snow depth
      ,PSFCAVG(:,:) &      !< Average surface pressure
      ,T10AVG(:,:) &       !< Average temperature at 10 meters height
      ,AKHSAVG(:,:) &      !< Time-averaged mass exchange coefficient__??
      ,AKMSAVG(:,:) &      !< Time-averaged wind exchange coefficient__??
      ,T10M(:,:) &         !< Temperature at 10 meters height
      ,U10MAX(:,:) &       !< Maximum hourly zonal wind speed at 10 meters height
      ,V10MAX(:,:) &       !< Maximum hourly meridional wind speed at 10 meters height
      ,u10h(:,:) &         !< Hourly zonal wind speed at 10 meters height
      ,v10h(:,:) &         !< Hourly meridional wind speed at 10 meters height
      ,PRATE_MAX(:,:) &    !< Maximum precipitation rate 
      ,FPRATE_MAX(:,:) &   !< Fractional coverage of maximum precipitation rate
! GSD addition
      ,WSPD10MAX(:,:) &       !< Maximum hourly wind speed at 10 meters heigh
      ,W_UP_MAX(:,:) &        !< Maximum upward vertical velocity
      ,W_DN_MAX(:,:) &        !< Maximum downward vertical velocity
      ,REFD_MAX(:,:) &        !< Maximum reflectivity factor
      ,UP_HELI_MAX(:,:) &     !< Maximum updraft helicity
      ,UP_HELI_MAX16(:,:) &   !< Maximum updraft helicity over 16 grid points
      ,GRPL_MAX(:,:) &        !< Maximum graupel mixing ratio
      ,QRMAX(:,:) &           !< Maximum rain mixing ratio
      ,UP_HELI(:,:) &         !< Updraft helicity
      ,UP_HELI16(:,:) &       !< Updraft helicity over 16 grid points
      ,LTG1_MAX(:,:) &        !< Maximum lightning flash rate in category 1___?
      ,LTG2_MAX(:,:) &        !< Maximum lightning flash rate in category 2___?
      ,LTG3_MAX(:,:) &        !< Maximum lightning flash rate in category 3___?
      ,UP_HELI_MIN(:,:) &     !< Minimum updraft helicity
      ,UP_HELI_MIN16(:,:) &   !< Minimum updraft helicity over 16 grid points
      ,UP_HELI_MAX02(:,:) &   !< Maximum updraft helicity over 0-2 km layer
      ,UP_HELI_MIN02(:,:) &   !< Minimum updraft helicity over 0-2 km layer
      ,UP_HELI_MAX03(:,:) &   !< Maximum updraft helicity over 0-3 km layer
      ,UP_HELI_MIN03(:,:) &   !< Minimum updraft helicity over 0-3 km layer
      ,REL_VORT_MAX(:,:) &    !< Maximum relative vorticity
      ,REL_VORT_MAX01(:,:) &  !< Maximum relative vorticity at 0-1 km layer
      ,REL_VORT_MAXHY1(:,:) & !< Maximum relative vorticity at hybrid 1 layer
      ,WSPD10UMAX(:,:) &      !< Maximum hourly u-component of wind at 10 meters height
      ,WSPD10VMAX(:,:) &      !< Maximum hourly v-component of wind at 10 meters height
      ,REFDM10C_MAX(:,:) &    !< Maximum reflectivity factor at model level 10 ___?
      ,HAIL_MAX2D(:,:) &      !< Maximum hail size in 2D
      ,HAIL_MAXK1(:,:) &      !< Maximum hail size at hybrid level 1
      ,HAIL_MAXHAILCAST(:,:) & !< Maximum hail size forecasted by HAILCAST
      ,NCI_LTG(:,:) &         !< Number of intra-cloud lightning flashes
      ,NCA_LTG(:,:) &         !< Number of cloud-to-air lightning flashes
      ,NCI_WQ(:,:) &          !< Number of intra-cloud water quantity
      ,NCA_WQ(:,:) &          !< Number of cloud-to-air water quantity
      ,NCI_REFD(:,:) &        !< Number of intra-cloud reflectivity factor
      ,NCA_REFD(:,:) &        !< Number of cloud-to-air reflectivity factor
      ,RAINC_BUCKET1(:,:) &   !< Rainc from bucket model 1___?
      ,RAINNC_BUCKET1(:,:) &  !< Rainnc from bucket model 1 ___?
      ,RAINC_BUCKET(:,:) &    !< Rainc from bucket mode___?
      ,RAINNC_BUCKET(:,:) &   !< Rainnc from bucket model___?
      ,SNOW_BUCKET(:,:) &     !< Snow from bucket model
      ,GRAUP_BUCKET(:,:) &    !< Graupel from bucket model
      ,PCP_BUCKET(:,:) &      !< Precipitation from bucket model
      ,ACGRAUP(:,:) &         !< Accumulated graupel
      ,ACFRAIN(:,:) &         !< Accumulated freezing rain
      ,FRZRN_BUCKET(:,:) &    !< Freezing rain from bucket model
      ,SNOW_ACM(:,:) &        !< Snow accumulation
      ,SNOW_BKT(:,:) &        !< Snow from bucke model
      ,SNOW_BUCKET1(:,:) &    !< Snow from bucket model 1
      ,GRAUP_BUCKET1(:,:) &   !< Graupel from bucket model 1
      ,PCP_BUCKET1(:,:) &     !< Precipitation from bucket model 1
      ,SNOWNC(:,:) &          !< Snow number concentration
      ,GRAUPELNC(:,:) &       !< Graupel number concentration
      ,TMAX(:,:) &            !< Maximum temperature
      ,W_MEAN(:,:) &          !< Mean vertical velocity
      ,TSNOW(:,:) &           !< Temperature at the time of snowfall
      ,QVG(:,:) &             !< Vapor growth rate
      ,QV2m(:,:) &            !< Specific humidity at 2 meters above ground level
      ,QVl1(:,:) &            !< Specific humidity at the first level
      ,REFC_10CM(:,:) &       !< Reflectivity at 10 cm
      ,REF1KM_10CM(:,:) &     !< Reflectivity at 1 km for 10 cm___?
      ,REF4KM_10CM(:,:) &     !< Reflectivity at 4 km for 10 cm___?
      ,SWRADmean(:,:) &       !< Mean shortwave radiation
      ,U10mean(:,:) &         !< Mean u-component of wind at 10 meters above ground level
      ,V10mean(:,:) &         !< Mean v-component of wind at 10 meters above ground level
      ,SPDUV10mean(:,:) &     !< Mean horizontal wind speed at 10 meters above ground level
      ,SWNORMmean(:,:) &      !< Mean normal shortwave radiation
      ,SNFDEN(:,:) &          !< Snowfall density
      ,SNDEPAC(:,:) &         !< Snow depth accumulation
      ,SWDDNI(:,:) &          !< Diffuse downwelling shortwave radiatio___?
      ,SWDDIF(:,:) &          !< Direct downwelling shortwave radiation___? 
      ,SWDNBC(:,:) &          !< Normal beam downward shortwave radiation flux at the surface__?
      ,SWDDNIC(:,:) &         !< Diffuse downward shortwave radiation flux at the surface___?
      ,SWDDIFC(:,:) &         !< Direct downward shortwave radiation flux at the surface___?
      ,SWUPBC(:,:) &          !< Normal beam upward shortwave radiation flux at the surface__?
      ,SWUPT(:,:) &           !< Diffuse upward shortwave radiation flux at the surface___?
      ,TAOD5502D(:,:) &       !< Aerosol optical depth at 550 nm
      ,AERASY2D(:,:) &        !< Aerosol asymmetry parameter
      ,AERSSA2D(:,:) &        !< Aerosol single-scattering albedo
      ,MEAN_FRP(:,:) &        !< Mean fire radiative power ___?
      ,HWP(:,:) &             !< Hourly wildfire potential
      ,LWP(:,:) &             !< Liquid water path ___?
      ,IWP(:,:) &             !< Ice water path___?
      ,XLAIXY(:,:) &          !< ___?
! add new fields for GFS
      ,SFCUX(:,:) &           !< Time averaged zonal momentum flux
      ,SFCVX(:,:) &           !< Time averaged meridional momentum flux
      ,SFCUXI(:,:) &          !< Instantaneous zonal momentum flux
      ,SFCVXI(:,:) &          !< Instantaneous meridional momentum flux
      ,AVGALBEDO(:,:) &       !< Mid day avg albedo
      ,AVGCPRATE(:,:) &       !< Convective precip - coninuous bucket
      ,AVGPREC(:,:) &         !< Average precip rate in m per physics time step
      ,PTOP(:,:) &            !< Instantaneous convective cloud top pressure
      ,PBOT(:,:) &            !< Instantaneous convective cloud bottom pressure
      ,AVGCFRACH(:,:) &       !< Average high cloud fraction
      ,AVGCFRACM(:,:) &       !< Average mid cloud fraction
      ,AVGCFRACL(:,:) &       !< Average low cloud fraction
      ,AVGTCDC(:,:) &         !< Time-averaged column cloud fraction
      ,AUVBIN(:,:) &          !< Time averaged incoming sfc uv-b
      ,AUVBINC(:,:) &         !< Time averaged incoming sfc clear sky uv-b
      ,ptopl(:,:) &           !< Time averaged low cloud top pressure
      ,pbotl(:,:) &           !< Time averaged low cloud bottom pressure
      ,Ttopl(:,:) &           !< Time averaged low cloud top temperature 
      ,ptopm(:,:) &           !< Time averaged middle cloud top pressure
      ,pbotm(:,:) &           !< Time averaged middle cloud bottom pressure 
      ,Ttopm(:,:) &           !< Time averaged middle cloud top temperature
      ,ptoph(:,:) &           !< Time averaged high cloud top pressure
      ,pboth(:,:) &           !< Time averaged high cloud top pressure 
      ,Ttoph(:,:) &           !< Time averaged high cloud top temperature
      ,sfcugs(:,:) &          !< ___?
      ,sfcvgs(:,:) &          !< ___?
      ,PBLCFR(:,:) &          !< Boundary layer cloud cover
      ,cldwork(:,:) &         !< Cloud work function
      ,gtaux(:,:) &           !< Time averaged zonal gravity wave stress
      ,gtauy(:,:) &           !< Time averaged meridional gravity wave stress
      ,runoff(:,:) &          !< Accumulated total (base+surface) runoff
      ,maxtshltr(:,:) &       !< Shelter max temperature
      ,mintshltr(:,:) &       !< Shelter min temperature
      ,maxrhshltr(:,:) &      !< Shelter max rh
      ,minrhshltr(:,:) &      !< Shelter min rh
      ,dzice(:,:) &           !< Ice thickness
      ,maxqshltr(:,:) &       !< Shelter max specific humidity
      ,minqshltr(:,:) &       !< Shelter min specific humidity
      ,alwinc(:,:) &          !< Time averaged surface clear sky incoming lw
      ,alwoutc(:,:) &         !< Time averaged surface clear sky outgoing lw
      ,alwtoac(:,:) &         !< Time averaged toa clear sky outgoing lw
      ,aswinc(:,:) &          !< Time averaged surface clear sky incoming sw
      ,aswoutc(:,:) &         !< Time averaged surface clear sky outgoing sw
      ,aswtoac(:,:) &         !< Time averaged toa clear sky outgoing sw
      ,aswintoa(:,:) &        !< Time averaged model top incoming shortwave
      ,smcwlt(:,:) &          !< Wilting point
      ,suntime(:,:) &         !< Sunshine duration
      ,fieldcapa(:,:) &       !< Field capacity 
      ,avisbeamswin(:,:) &    !< Time averaged surface visible beam downward solar flux
      ,avisdiffswin(:,:) &    !< Time averaged surface visible diffuse downward solar flux
      ,airbeamswin(:,:) &     !< Time averaged surface near ir beam downward solar flux
      ,airdiffswin(:,:) &     !< Time averaged surface near ir diffuse downward solar flux
      ,snowfall(:,:) &        !< Total accumulated snowfall___?
      ,acond(:,:) &           !< Aerodynamic conductance
      ,edir(:,:) &            !< Soil evaporation___?
      ,ecan(:,:) &            !< Accumulated evaporation of intercepted water__?
      ,etrans(:,:) &          !< Plant transpiration__?
      ,esnow(:,:) &           !< Snow sublimation___?
      ,avgedir(:,:) &         !< Direct soil evaporation__?
      ,avgecan(:,:) &         !< Canopy water evaporation
      ,avgetrans(:,:)&        !< Plant transpiration___?
      ,avgesnow(:,:) &        !< Snow sublimation
      ,avgpotevp(:,:) &       !< Time averaged accumulated potential evaporation
      ,avgprec_cont(:,:) &    !< Average precip rate - continuous bucket
      ,avgcprate_cont(:,:) &  !< Convective precip - coninuous bucket
      ,ti(:,:) &              !< Sea ice skin temperature
      ,aod550(:,:) &          !< Instantaneous aod550 optical depth
      ,du_aod550(:,:) &       !< Instantaneous aod550 optical depth (dust)
      ,ss_aod550(:,:) &       !< Instantaneous aod550 optical depth (seasalt)
      ,su_aod550(:,:) &       !< Instantaneous aod550 optical depth (sulfates)
      ,bc_aod550(:,:) &       !< Instantaneous aod550 optical depth (organic carbon)
      ,oc_aod550(:,:) &       !< Instantaneous aod550 optical depth (black carbon)
      ,landfrac(:,:) &        !< Land fraction
      ,paha(:,:) &            !< Averaged precipitation advected heat flux
      ,pahi(:,:) &            !< Instantaneous precipitation advected heat flux
      ,tecan(:,:) &           !< Accumulated evaporation of intercepted water
      ,tetran(:,:) &          !< Accumulated plant transpiration
      ,tedir(:,:) &           !< Accumulated soil surface evaporation
      ,twa(:,:) &             !< Total water storage in aquifer
      ,fdnsst(:,:) &          !< ___?
      ,pwat(:,:)              !< Precipitable water
      integer, allocatable :: IVGTYP(:,:) &     !< Vegetation type
      ,ISLTYP(:,:) &          !< Soil type
      ,ISLOPE(:,:) &          !< Slope type
      ,IEQL(:,:)              !< ___?
      
! Add 2d aerosol diagnosis fields for GOCART (NGAC)
      real, allocatable ::                                                   &
       DUSMASS(:,:) &         !< Mass of dust aerosols in the atmosphere___?
      ,DUCMASS(:,:) &         !< Mass of dust aerosols in the cloud___?
      ,DUSMASS25(:,:) &       !< Mass of dust aerosols with diameter greater than 2.5 micrometers in the atmosphere___?
      ,DUCMASS25(:,:) &       !< Mass of dust aerosols with diameter greater than 2.5 micrometers in the cloud
      ,SUSMASS(:,:) &         !< Mass of sulfate aerosols in the atmosphere
      ,SUCMASS(:,:) &         !< Mass of sulfate aerosols in the cloud
      ,SUSMASS25(:,:) &       !< Mass of sulfate aerosols with diameter greater than 2.5 micrometers in the atmosphere
      ,SUCMASS25(:,:) &       !< Mass of sulfate aerosols with diameter greater than 2.5 micrometers in the cloud
      ,OCSMASS(:,:) &         !< Mass of organic carbon aerosols in the atmosphere
      ,OCCMASS(:,:) &         !< Mass of organic carbon aerosols in the cloud
      ,OCSMASS25(:,:) &       !< Mass of organic carbon aerosols with diameter greater than 2.5 micrometers in the atmosphere
      ,OCCMASS25(:,:) &       !< Mass of organic carbon aerosols with diameter greater than 2.5 micrometers in the cloud
      ,BCSMASS(:,:) &         !< Mass of black carbon aerosols in the atmosphere
      ,BCCMASS(:,:) &         !< Mass of black carbon aerosols in the cloud
      ,BCSMASS25(:,:) &       !< Mass of black carbon aerosols with diameter greater than 2.5 micrometers in the atmosphere
      ,BCCMASS25(:,:) &       !< Mass of black carbon aerosols with diameter greater than 2.5 micrometers in the cloud
      ,SSSMASS(:,:) &         !< Mass of sea salt aerosols in the atmosphere
      ,SSCMASS(:,:) &         !< Mass of sea salt aerosols in the cloud
      ,SSSMASS25(:,:) &       !< Mass of sea salt aerosols with diameter greater than 2.5 micrometers in the atmosphere  
      ,SSCMASS25(:,:) &       !< Mass of sea salt aerosols with diameter greater than 2.5 micrometers in the cloud
      ,DUSTCB(:,:) &          !< Dust concentration in the atmosphere boundary layer
      ,SSCB(:,:) &            !< Sea salt concentration in the atmosphere boundary layer
      ,OCCB(:,:) &            !< Organic carbon concentration in the atmosphere boundary layer
      ,BCCB(:,:) &            !< Black carbon concentration in the atmosphere boundary layer
      ,SULFCB(:,:) &          !< Sulfate concentration in the atmosphere boundary layer
      ,DUSTALLCB(:,:) &       !< Total dust concentration in the atmosphere boundary layer
      ,SSALLCB(:,:) &         !< Total sea salt concentration in the atmosphere boundary layer
      ,DUSTPM(:,:) &          !< Dust particulate matter concentration
      ,SSPM(:,:) &            !< Sea salt particulate matter concentration 
      ,PP25CB(:,:) &          !< Particulate matter with diameter less than 2.5 micrometers concentration in the atmosphere boundary layer
      ,DUSTPM10(:,:) &        !< Dust particulate matter with diameter less than 10 micrometers concentration
      ,PP10CB(:,:) &          !< Particulate matter with diameter less than 10 micrometers concentration in the atmosphere boundary layer
      ,NO3CB(:,:) &           !< Nitrate concentration in the atmosphere boundary layer
      ,NH4CB(:,:) &           !< Ammonium concentration in the atmosphere boundary layer
      ,maod(:,:)              !< _____ !lzhang, add for FV3-Chem ___?

! Add new field for AQM
      real, allocatable :: aqm_aod550(:,:)   !< AQM aerosol optical depth
 
!
      end module vrbls2d
