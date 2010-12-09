       SUBROUTINE CALRAD_WCLOUD
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALRAD      
!   PRGRMMR: CHUANG        ORG: EMC      DATE: 07-01-17       
!     
! ABSTRACT:
!     THIS ROUTINE COMPUTES MODEL DERIVED BRIGHTNESS TEMPERATURE
!     USING CRTM. IT IS PATTERNED AFTER GSI SETUPRAD WITH TREADON'S HELP     
! PROGRAM HISTORY LOG:
!
! USAGE:    CALL MDLFLD
!   INPUT ARGUMENT LIST:
!     NONE
!   OUTPUT ARGUMENT LIST: 
!     NONE
!
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!
!     LIBRARY:
!     /nwprod/lib/sorc/crtm2
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : IBM
!$$$  
       use vrbls3d
       use vrbls2d
       use masks
       use soil
      
       use kinds, only: r_kind,r_single,i_kind
       use crtm_module, only: crtm_atmosphere_type,crtm_surface_type,crtm_geometryinfo_type, &
       crtm_allocate_surface,o3_id,wet_soil,crtm_forward,mass_mixing_ratio_units, &
       crtm_allocate_atmosphere,grass_scrub,grass_soil, meadow_grass,urban_concrete, &
       irrigated_low_vegetation,broadleaf_pine_forest,pine_forest,compacted_soil, &
       broadleaf_forest,broadleaf_brush,tundra,tilled_soil,scrub,scrub_soil,&
       crtm_options_type,crtm_destroy,crtm_init,SPECIFIC_AMOUNT_UNITS
       use crtm_rtsolution_define, only: crtm_rtsolution_type, crtm_allocate_rtsolution, &
       crtm_destroy_rtsolution
       use crtm_spccoeff, only: sc
       use crtm_atmosphere_define, only:h2o_id,crtm_assign_atmosphere, &
       crtm_destroy_atmosphere,volume_mixing_ratio_units,crtm_zero_atmosphere
       use crtm_surface_define, only: crtm_assign_surface,crtm_destroy_surface, &
       crtm_zero_surface
       use crtm_channelinfo_define, only: crtm_channelinfo_type
       use crtm_parameters, only: limit_exp,toa_pressure,max_n_layers
       use crtm_cloud_define, only:  water_cloud,ice_cloud,rain_cloud,snow_cloud
       use message_handler, only: success,warning, display_message

      
       use params_mod
       use rqstfld_mod
       use ctlblk_mod
!     
   implicit none

!     DECLARE VARIABLES.
!     
! Mapping land surface type of GFS to CRTM
!  Note: index 0 is water, and index 13 is ice. The two indices are not
!        used and just assigned to COMPACTED_SOIL.
        integer, parameter, dimension(0:13) :: gfs_to_crtm=(/COMPACTED_SOIL,     &
      &  BROADLEAF_FOREST, BROADLEAF_FOREST, BROADLEAF_PINE_FOREST, PINE_FOREST, &
      &  PINE_FOREST, BROADLEAF_BRUSH, SCRUB, SCRUB, SCRUB_SOIL, TUNDRA,         &
      &  COMPACTED_SOIL, TILLED_SOIL, COMPACTED_SOIL/)

! Mapping land surface type of NMM to CRTM
!  Note: index 16 is water, and index 24 is ice. The two indices are not
!        used and just assigned to COMPACTED_SOIL.
       integer, parameter, dimension(24) :: nmm_to_crtm=(/URBAN_CONCRETE,       &
      &   COMPACTED_SOIL, IRRIGATED_LOW_VEGETATION, GRASS_SOIL, MEADOW_GRASS,   &
      &   MEADOW_GRASS, MEADOW_GRASS, SCRUB, GRASS_SCRUB, MEADOW_GRASS,         &
      &   BROADLEAF_FOREST, PINE_FOREST, BROADLEAF_FOREST, PINE_FOREST,         &
      &   BROADLEAF_PINE_FOREST, COMPACTED_SOIL, WET_SOIL, WET_SOIL,            &
      &   IRRIGATED_LOW_VEGETATION, TUNDRA, TUNDRA, TUNDRA, TUNDRA,             &
      &   COMPACTED_SOIL/)

      integer, parameter:: ndat=100
! CRTM structure variable declarations.
      integer,parameter::  n_absorbers = 2
!      integer,parameter::  n_clouds = 4 
      integer,parameter::  n_aerosols = 0
      integer(i_kind),parameter:: n_sensors=1
      character(len=20),allocatable,dimension(:):: sensorlist ! CRTM satellite/sensor list

      integer(i_kind) sensorindex
      integer(i_kind) lunin,nobs,nchanl,nreal
      integer(i_kind) error_status,itype
      integer(i_kind) err1,err2,err3,err4
      integer(i_kind) i,j,k,msig
      integer(i_kind) lcbot,lctop   !bsf
      integer jdn
  
      real(r_kind),parameter:: r100=100.0_r_kind
      real,parameter:: ozsmall = 1.e-10 ! to convert to mass mixing ratio
      real(r_kind) tsfc 
      real(r_kind),dimension(4):: sfcpct
      real(r_kind) snodepth,snoeqv,vegcover
      real snofrac
      real(r_kind),dimension(im,jsta:jend):: tb1,tb2,tb3,tb4
      real,dimension(im,jm):: grid1,grid2
      real sun_zenith,sun_azimuth, dpovg
      real q_conv   !bsf
  
      character(10)::obstype
      character(20)::isis
  
      logical hirs2,msu,goessndr,hirs3,hirs4,hirs,amsua,amsub,airs,hsb  &
     &      ,goes_img,mhs
      logical avhrr,avhrr_navy,lextra,ssu
      logical ssmi,ssmis,amsre,amsre_low,amsre_mid,amsre_hig,change
      logical ssmis_las,ssmis_uas,ssmis_env,ssmis_img
      logical sea,mixed,land,ice,snow,toss
      logical micrim,microwave
!  logical,dimension(nobs):: luse

      type(crtm_atmosphere_type),dimension(1):: atmosphere
      type(crtm_surface_type),dimension(1) :: surface
      type(crtm_geometryinfo_type),dimension(1) :: geometryinfo
      type(crtm_options_type),dimension(1)      :: options

      type(crtm_rtsolution_type),allocatable,dimension(:,:):: rtsolution
      type(crtm_channelinfo_type),allocatable,dimension(:) :: channelinfo
!     
      integer ii,jj,n_clouds
      integer,external :: iw3jdn
!
!*****************************************************************************
!     START SUBROUTINE CALRAD.
      if (iget(327) > 0 .or. iget(328) > 0 .or. iget(329) > 0       &
     & .or. iget(330) > 0 ) then
! specify numbers of cloud species     
       if(MODELNAME == 'GFS')then
        n_clouds = 2 ! GFS uses Zhao scheme
       else
        n_clouds = 4
       end if  
! Initialize ozone to zeros for WRF NMM and ARW for now
       if (MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR')o3=0.0
! Compute solar zenith angle for GFS
       if (MODELNAME == 'GFS')then
        jdn=iw3jdn(idat(3),idat(1),idat(2))
	do j=jsta,jend
	 do i=1,im
	  call zensun(jdn,float(idat(4)),gdlat(i,j),gdlon(i,j)       &
     &	    ,pi,sun_zenith,sun_azimuth)
          czen(i,j)=cos(sun_zenith)
	 end do
	end do
!	ii=361
!	jj=278
	ii=1035
	jj=219
	if(jj>=jsta .and. jj<=jend)                                  &
     &    print*,'sample GFS zenith angle=',acos(czen(ii,jj))*rtd   
       end if	       
! Initialize CRTM.  Load satellite sensor array.
! The optional arguments Process_ID and Output_Process_ID limit
! generation of runtime informative output to mpi task
! Output_Process_ID (which here is set to be task 0)
       print*,'success in CALRAD= ',success
       allocate( channelinfo(n_sensors))
       allocate(sensorlist(n_sensors))
       sensorlist(1)='imgr_g12'

       error_status = crtm_init(channelinfo,SensorID=sensorlist,&
          Process_ID=0,Output_Process_ID=0 )

       if (error_status /= 0_i_kind)                                  &
     &   write(6,*)'ERROR*** crtm_init error_status=',error_status
  
!   lunin=1 ! will read data file in the future, only simulate GOES for now
!   open(lunin,file='obs_setup',form='unformatted') ! still need to find out filename
!   rewind lunin
  
! Loop over data types to process    
!  do is=1,ndat
!    read(lunin,end=125) obstype,isis,nreal,nchanl
! Test GOES for now.
       obstype='goes_img'
       isis='imgr_g12'
!       isis='amsua_n15'

! Initialize logical flags for satellite platform

       hirs2      = obstype == 'hirs2'
       hirs3      = obstype == 'hirs3'
       hirs4      = obstype == 'hirs4'
       hirs       = hirs2 .or. hirs3 .or. hirs4
       msu        = obstype == 'msu'
       ssu        = obstype == 'ssu'
       goessndr   = obstype == 'sndr'  .or. obstype == 'sndrd1' .or.    &
     &                obstype == 'sndrd2'.or. obstype == 'sndrd3' .or.  &
     &                obstype == 'sndrd4'
       amsua      = obstype == 'amsua'
       amsub      = obstype == 'amsub'
       mhs        = obstype == 'mhs'
       airs       = obstype == 'airs'
       hsb        = obstype == 'hsb'
       goes_img   = obstype == 'goes_img'
       avhrr      = obstype == 'avhrr'
       avhrr_navy = obstype == 'avhrr_navy'
       ssmi       = obstype == 'ssmi'
       amsre_low  = obstype == 'amsre_low'
       amsre_mid  = obstype == 'amsre_mid'
       amsre_hig  = obstype == 'amsre_hig'
       amsre      = amsre_low .or. amsre_mid .or. amsre_hig
       ssmis      = obstype == 'ssmis'
       ssmis_las  = obstype == 'ssmis_las'
       ssmis_uas  = obstype == 'ssmis_uas'
       ssmis_img  = obstype == 'ssmis_img'
       ssmis_env  = obstype == 'ssmis_env'

       ssmis=ssmis_las.or.ssmis_uas.or.ssmis_img.or.ssmis_env

       micrim=ssmi .or. ssmis .or. amsre   ! only used for MW-imager-QC and id_qc(ch)

       microwave=amsua .or. amsub  .or. mhs .or. msu .or. hsb .or.  &
     &             micrim 

! Determine specific sensor
   sensorindex = 0
   sensor_search: do j = 1, n_sensors
     if (channelinfo(j)%sensor_id == isis ) then
         sensorindex = j
         exit sensor_search
     endif
   end do sensor_search
   if (sensorindex == 0 ) then
      write(6,*)'SETUPRAD:  ***WARNING*** problem with sensorindex=',isis,&
           ' --> CAN NOT PROCESS isis=',isis,'   TERMINATE PROGRAM EXECUTION'
      stop
   endif

! Set CRTM to process given satellite/sensor
!       Error_Status = CRTM_Set_ChannelInfo(isis, ChannelInfo)
!       if (error_status /= success)                                 &
!     &    write(6,*)'ERROR*** crtm_set_channelinfo error_status=',  &
!     &    error_status,' for satsensor=',isis

! Allocate structures for radiative transfer
       print*,'channel number= ',channelinfo(sensorindex)%n_channels
       allocate(rtsolution  (channelinfo(sensorindex)%n_channels,1))
       err1=0; err2=0; err3=0; err4=0
       if(lm > max_n_layers)then
        write(6,*) 'CALRAD: lm > max_n_layers - '//                 &
     &	   'increase crtm max_n_layers ',                            &
     &     lm,max_n_layers
        stop 2
       end if
       err1 = crtm_allocate_atmosphere(lm,n_absorbers,n_clouds      &
     &       ,n_aerosols,atmosphere(1))
       err2 = crtm_allocate_surface(channelinfo(sensorindex)%n_channels,surface(1))
       err3 = crtm_allocate_rtsolution(lm,rtsolution)
       if (err1/=success) write(6,*)' ERROR** allocating atmosphere.err=',&
     &       err1
       if (err2/=success) write(6,*)' ***ERROR** allocating surface.err=',&
     &       err2
       if (err3/=success) write(6,*)' ***ERROR** allocating rtsolution.err=',&
     &       err3

       atmosphere(1)%n_layers = lm
!      atmosphere(1)%level_temperature_input = 0
       atmosphere(1)%absorber_id(1) = H2O_ID
       atmosphere(1)%absorber_id(2) = O3_ID
       atmosphere(1)%absorber_units(1) = MASS_MIXING_RATIO_UNITS
!       atmosphere(1)%absorber_units(2) = VOLUME_MIXING_RATIO_UNITS
       atmosphere(1)%absorber_units(2) = MASS_MIXING_RATIO_UNITS ! Use mass mixing ratio
       atmosphere(1)%level_pressure(0) = TOA_PRESSURE
! Ddefine Clouds
       if(MODELNAME == 'GFS')then
        atmosphere(1)%cloud(1)%n_layers = lm
        atmosphere(1)%cloud(1)%Type = WATER_CLOUD
        atmosphere(1)%cloud(2)%n_layers = lm
        atmosphere(1)%cloud(2)%Type = ICE_CLOUD
       else
        atmosphere(1)%cloud(1)%n_layers = lm
        atmosphere(1)%cloud(1)%Type = WATER_CLOUD
        atmosphere(1)%cloud(2)%n_layers = lm
        atmosphere(1)%cloud(2)%Type = ICE_CLOUD
        atmosphere(1)%cloud(3)%n_layers = lm
        atmosphere(1)%cloud(3)%Type = RAIN_CLOUD
        atmosphere(1)%cloud(4)%n_layers = lm
        atmosphere(1)%cloud(4)%Type = SNOW_CLOUD
       end if        

!    if(nchanl /= channelinfo(sensorindex)%n_channels) write(6,*)'***ERROR** nchanl,n_channels ', &
!           nchanl,channelinfo(sensorindex)%n_channels

! Load surface sensor data structure
       surface(1)%sensordata%n_channels = channelinfo(sensorindex)%n_channels
       surface(1)%sensordata%sensor_id  = channelinfo(sensorindex)%wmo_sensor_id

! Loop through all grid points
!      ii=361
!      jj=278
      do j=jsta,jend
       do i=1,im

!    Load geometry structure
!    geometryinfo(1)%sensor_zenith_angle = zasat*rtd  ! local zenith angle ???????
        geometryinfo(1)%sensor_zenith_angle = 0. ! 44.
        geometryinfo(1)%source_zenith_angle = acos(czen(i,j))*rtd ! solar zenith angle
!        geometryinfo(1)%source_zenith_angle = 0. 
        geometryinfo(1)%sensor_scan_angle   = 0. ! scan angle, assuming nadir
        if(abs(geometryinfo(1)%sensor_zenith_angle) > h1 ) then
          geometryinfo(1)%distance_ratio =                                 &
     &        abs( sin(geometryinfo(1)%sensor_scan_angle*dtr)/             &
     &        sin(geometryinfo(1)%sensor_zenith_angle*dtr) )
        endif
        if(i==ii.and.j==jj)print*,'sample geometry ',                   &
     &  geometryinfo(1)%sensor_zenith_angle                                &
     &  ,geometryinfo(1)%source_zenith_angle                               &
     &  ,czen(i,j)*rtd 
!  Set land/sea, snow, ice percentages and flags
        if (MODELNAME == 'GFS')then ! GFS uses 13 veg types
         itype=IVGTYP(I,J)
         itype = min(max(0,ivgtyp(i,j)),13)
!         IF(itype <= 0 .or. itype > 13)itype=7 !use scrub for ocean point
         if(sno(i,j)/=spval)then
          snoeqv=sno(i,j)
         else
          snoeqv=0.
         end if
         if(i==ii.and.j==jj)print*,'sno,itype,ivgtyp B cing snfrc = ',  &
     &      snoeqv,itype,IVGTYP(I,J)
         if(sm(i,j) > 0.1)then
          sfcpct(4)=0.
         else 
          call snfrac_gfs(i,j,SNOeqv,itype,sfcpct(4))
         end if
         if(i==ii.and.j==jj)print*,'sno,itype,ivgtyp,sfcpct(4) = ',     &
     &      snoeqv,itype,IVGTYP(I,J),sfcpct(4)
        else          
         itype=IVGTYP(I,J)
         IF(itype == 0)itype=8
         CALL SNFRAC (SNO(I,J),itype,sfcpct(4))
        end if 
!	CALL SNFRAC (SNO(I,J),IVGTYP(I,J),snofrac)
!	sfcpct(4)=snofrac
	if(sm(i,j) > 0.1)then ! water
!	 tsfc=sst(i,j)
         tsfc = ths(i,j)*(pint(i,j,nint(lmh(i,j))+1)/p1000)**capa
         vegcover=0.0
	 if(sfcpct(4) > 0.0_r_kind)then ! snow and water
          sfcpct(1) = 1.0_r_kind-sfcpct(4)
          sfcpct(2) = 0.0_r_kind
	  sfcpct(3) = 0.0_r_kind
	 else ! pure water
	  sfcpct(1) = 1.0_r_kind
          sfcpct(2) = 0.0_r_kind
	  sfcpct(3) = 0.0_r_kind
	 end if  
        else ! land and sea ice
	 tsfc = ths(i,j)*(pint(i,j,nint(lmh(i,j))+1)/p1000)**capa
         vegcover=vegfrc(i,j)
	 if(sice(i,j) > 0.1)then ! sea ice
	  if(sfcpct(4) > 0.0_r_kind)then ! sea ice and snow
	   sfcpct(3) = 1.0_r_kind-sfcpct(4)
	   sfcpct(1) = 0.0_r_kind
           sfcpct(2) = 0.0_r_kind
	  else ! pure sea ice
	   sfcpct(3)= 1.0_r_kind
	   sfcpct(1) = 0.0_r_kind
           sfcpct(2) = 0.0_r_kind
	  end if
	 else ! land
	  if(sfcpct(4) > 0.0_r_kind)then ! land and snow
	   sfcpct(2)= 1.0_r_kind-sfcpct(4)
           sfcpct(1) = 0.0_r_kind
           sfcpct(3) = 0.0_r_kind
	  else ! pure land
	   sfcpct(2)= 1.0_r_kind
           sfcpct(1) = 0.0_r_kind
           sfcpct(3) = 0.0_r_kind
	  end if  
         end if
	end if 
        if(si(i,j)/=spval)then 	
	 snodepth = si(i,j)
        else
         snodepth = 0.
        end if

        sea  = sfcpct(1)  >= 0.99_r_kind
        land = sfcpct(2)  >= 0.99_r_kind
        ice  = sfcpct(3)  >= 0.99_r_kind
        snow = sfcpct(4)  >= 0.99_r_kind
        mixed = .not. sea  .and. .not. ice .and.     &
     &             .not. land .and. .not. snow
        if((sfcpct(1)+sfcpct(2)+sfcpct(3)+sfcpct(4))       &
     &    >1._r_kind)print*,'ERROR sfcpct ',i,j,sfcpct(1)  &
     &    ,sfcpct(2),sfcpct(3),sfcpct(4)
!    Load surface structure

!    Define land characteristics

!    **NOTE:  The model surface type --> CRTM surface type
!             mapping below is specific to the versions NCEP
!             GFS and NNM as of September 2005
!    itype = ivgtyp(i,j)
       if (MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR') then
        itype = min(max(1,ivgtyp(i,j)),24)
        surface(1)%land_type = nmm_to_crtm(itype)
       else
        itype = min(max(0,ivgtyp(i,j)),13)
        surface(1)%land_type = gfs_to_crtm(itype)
       end if

       surface(1)%wind_speed            = sqrt(u10(i,j)*u10(i,j)   &
     &                                   +v10(i,j)*v10(i,j))    
       surface(1)%water_coverage        = sfcpct(1)
       surface(1)%land_coverage         = sfcpct(2)
       surface(1)%ice_coverage          = sfcpct(3)
       surface(1)%snow_coverage         = sfcpct(4)
       
       surface(1)%land_temperature      = tsfc
       surface(1)%snow_temperature      = min(tsfc,280._r_kind)
       surface(1)%water_temperature     = max(tsfc,270._r_kind)
       surface(1)%ice_temperature       = min(tsfc,280._r_kind)
       if(smstot(i,j)/=spval)then
        surface(1)%soil_moisture_content = smstot(i,j)/10. !convert to cgs !???
       else
        surface(1)%soil_moisture_content = 0.05 ! default crtm value
       end if		
       surface(1)%vegetation_fraction   = vegcover
!       surface(1)%vegetation_fraction   = vegfrc(i,j)
       surface(1)%soil_temperature      = 283.
!       surface(1)%soil_temperature      = stc(i,j,1)
       surface(1)%snow_depth            = snodepth ! in mm
       
       if(surface(1)%wind_speed<0. .or. surface(1)%wind_speed>200.)  &
     &      print*,'bad 10 m wind'
       if(surface(1)%water_coverage<0. .or. surface(1)%water_coverage>1.) &
     &      print*,'bad water coverage'
       if(surface(1)%land_coverage<0. .or. surface(1)%land_coverage>1.)  &
     &      print*,'bad land coverage'
       if(surface(1)%ice_coverage<0. .or. surface(1)%ice_coverage>1.)  &
     &      print*,'bad ice coverage'
       if(surface(1)%snow_coverage<0. .or. surface(1)%snow_coverage>1.)  &
     &      print*,'bad snow coverage'
       if(surface(1)%land_temperature<0. .or. surface(1)%land_temperature>350.)  &
     &      print*,'bad land T'
       if(surface(1)%soil_moisture_content<0. .or. surface(1)%soil_moisture_content>600.) &
     &      print*,'bad soil_moisture_content'
       if(surface(1)%vegetation_fraction<0. .or. surface(1)%vegetation_fraction>1.) &
     &      print*,'bad vegetation cover'
       if(surface(1)%snow_depth<0. .or.  surface(1)%snow_depth>10000.) &
     &      print*,'bad snow_depth'
       if(MODELNAME == 'GFS' .and. (itype<0 .or. itype>13)) &
     &      print*,'bad veg type'
       
       
       if(i==ii.and.j==jj)print*,'sample surface in CALRAD=', &
     &   i,j,surface(1)%wind_speed,surface(1)%water_coverage,       &
     &   surface(1)%land_coverage,surface(1)%ice_coverage,          &
     &   surface(1)%snow_coverage,surface(1)%land_temperature,      &
     &   surface(1)%snow_temperature,surface(1)%water_temperature,  &
     &   surface(1)%ice_temperature,surface(1)%vegetation_fraction, &
     &   surface(1)%soil_temperature,surface(1)%snow_depth,         &
     &   surface(1)%land_type,sm(i,j)

!       Load profiles into model layers

!       Load atmosphere profiles into RTM model layers
!       CRTM counts from top down just as post does
       if(i==ii.and.j==jj)print*,'TOA= ',atmosphere(1)%level_pressure(0)
       do k = 1,lm
        atmosphere(1)%level_pressure(k) = pint(i,j,k+1)/r100
        atmosphere(1)%pressure(k)       = pmid(i,j,k)/r100
        atmosphere(1)%temperature(k)    = t(i,j,k)
	atmosphere(1)%absorber(k,1)     = max(0.  &
     &                                 ,q(i,j,k)*h1000/(h1-q(i,j,k))) ! use mixing ratio like GSI
        atmosphere(1)%absorber(k,2)     = max(ozsmall,o3(i,j,k)*h1000) ! convert to g/kg
! fill in cloud mixing ratio later  
        if(atmosphere(1)%level_pressure(k)<0. .or. atmosphere(1)%level_pressure(k)>1060.) &
     &      print*,'bad atmosphere(1)%level_pressure'  &
     &      ,i,j,k,atmosphere(1)%level_pressure(k)     
        if(atmosphere(1)%pressure(k)<0. .or.   &
     &      atmosphere(1)%pressure(k)>1060.)  &
     &      print*,'bad atmosphere(1)%pressure'  &
     &      ,i,j,k,atmosphere(1)%pressure(k) 
        if(atmosphere(1)%temperature(k)<0. .or.   &
     &      atmosphere(1)%temperature(k)>400.)  &
     &      print*,'bad atmosphere(1)%temperature'
!        if(atmosphere(1)%absorber(k,1)<0. .or.   &
!     &      atmosphere(1)%absorber(k,1)>1.)  &
!     &      print*,'bad atmosphere water vapor'
!        if(atmosphere(1)%absorber(k,2)<0. .or.   &
!     &      atmosphere(1)%absorber(k,1)>1.)  &
!     &      print*,'bad atmosphere o3'
        
        if(i==ii.and.j==jj)print*,'sample atmosphere in CALRAD=',  &
     &	  i,j,k,atmosphere(1)%level_pressure(k),atmosphere(1)%pressure(k),  &
     &    atmosphere(1)%temperature(k),atmosphere(1)%absorber(k,1),  &
     &    atmosphere(1)%absorber(k,2)
! Specify clouds
        dpovg=(pint(i,j,k+1)-pint(i,j,k))/g !crtm uses column integrated field
        if(MODELNAME == 'GFS')then
         atmosphere(1)%cloud(1)%effective_radius(k) = 10.
	 atmosphere(1)%cloud(1)%water_content(k) = max(0.,qqw(i,j,k)*dpovg)
! GFS uses temperature and ice concentration dependency formulation to determine effetive radis for cloud ice
! since GFS does not output ice concentration yet, use default 50 um
	 atmosphere(1)%cloud(2)%effective_radius(k) = 50.	 
	 atmosphere(1)%cloud(2)%water_content(k) = max(0.,qqi(i,j,k)*dpovg)
         if(atmosphere(1)%cloud(1)%water_content(k)<0. .or.   &
     &      atmosphere(1)%cloud(1)%water_content(k)>1.)  &
     &      print*,'bad atmosphere cloud water'
         if(atmosphere(1)%cloud(2)%water_content(k)<0. .or.   &
     &      atmosphere(1)%cloud(2)%water_content(k)>1.)  &
     &      print*,'bad atmosphere cloud ice'
        else
	 atmosphere(1)%cloud(1)%effective_radius(k) = 10.
	 atmosphere(1)%cloud(1)%water_content(k) = max(0.,qqw(i,j,k)*dpovg)
	 atmosphere(1)%cloud(2)%effective_radius(k) = 25.
	 atmosphere(1)%cloud(2)%water_content(k) = max(0.,qqi(i,j,k)*dpovg)
	 atmosphere(1)%cloud(3)%effective_radius(k) = 200.
	 atmosphere(1)%cloud(3)%water_content(k) = max(0.,qqr(i,j,k)*dpovg)
	 atmosphere(1)%cloud(4)%effective_radius(k) = 250.
	 atmosphere(1)%cloud(4)%water_content(k) = max(0.,qqs(i,j,k)*dpovg)
        end if 
       end do

!bsf - start
!-- Add subgrid-scale convective clouds for WRF runs
       if (MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR') then
         lcbot=nint(hbot(i,j))
         lctop=nint(htop(i,j))
         if (lcbot-lctop > 1) then
!-- q_conv = assumed grid-averaged cloud water/ice condensate from convection (Cu)
!   In "params" Qconv=0.1e-3 and TRAD_ice=253.15; cnvcfr is the Cu cloud fraction
!   calculated as a function of Cu rain rate (Slingo, 1987) in subroutine MDLFLD
           q_conv = cnvcfr(i,j)*Qconv
           do k = lctop,lcbot
             dpovg=(pint(i,j,k+1)-pint(i,j,k))/g
             if (t(i,j,k) < TRAD_ice) then
	     atmosphere(1)%cloud(2)%water_content(k) =   &
     &              atmosphere(1)%cloud(2)%water_content(k) + dpovg*q_conv
             else
	     atmosphere(1)%cloud(1)%water_content(k) =   &
     &              atmosphere(1)%cloud(1)%water_content(k) + dpovg*q_conv
             endif
           end do   !-- do k = lctop,lcbot
         endif      !-- if (lcbot-lctop > 1) then
       endif        !-- if (MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR') then
!bsf - end

!     call crtm forward model
       error_status = crtm_forward(atmosphere,surface,                 &
     &         geometryinfo,channelinfo(sensorindex:sensorindex),&
     &         rtsolution)
       if (error_status /=0) then
        print*,'***ERROR*** during crtm_forward call ',  &
     &       error_status
       else 	 
        tb1(i,j)=rtsolution(1,1)%brightness_temperature
        tb2(i,j)=rtsolution(2,1)%brightness_temperature
        tb3(i,j)=rtsolution(3,1)%brightness_temperature	 
        tb4(i,j)=rtsolution(4,1)%brightness_temperature
	if(i==ii.and.j==jj)print*,'sample rtsolution in CALRAD=',  &
     &    rtsolution(1,1)%brightness_temperature,  &
     &    rtsolution(2,1)%brightness_temperature,  &
     &    rtsolution(3,1)%brightness_temperature,  &
     &    rtsolution(4,1)%brightness_temperature   
	if(i==ii.and.j==jj)print*,'sample TB in CALRAD=', &
     &    tb1(i,j),tb2(i,j),tb3(i,j),tb4(i,j)   
!        if(tb1(i,j) < 400. )  &
!     &        print*,'good tb1 ',i,j,tb1(i,j),gdlat(i,j),gdlon(i,j)
!        if(tb2(i,j) > 400.)print*,'bad tb2 ',i,j,tb2(i,j)
!        if(tb3(i,j) > 400.)print*,'bad tb3 ',i,j,tb3(i,j)
!        if(tb4(i,j) > 400.)print*,'bad tb4 ',i,j,tb4(i,j)
       end if  	 
       end do ! end loop for i
      end do ! end loop for j 
  
      error_status = crtm_destroy(channelinfo)
      if (error_status /= success) &
     &   print*,'ERROR*** crtm_destroy error_status=',error_status
  
      if( iget(327) > 0 ) then
       do j=jsta,jend
        do i=1,im
         grid1(i,j)=tb1(i,j)
!         if(grid1(i,j)>400. .or. grid1(i,j)/=grid1(i,j))grid1(i,j)=400.
        enddo
       enddo
       id(1:25) = 0
       id(02) = 129
       call gribit(iget(327),lvls(1,iget(327)), grid1,im,jm)
      endif
  
      if( iget(328) > 0 ) then !water vapor channel
       do j=jsta,jend
        do i=1,im
         grid1(i,j)=tb2(i,j)
        enddo
       enddo
       id(1:25) = 0
       id(02) = 129
       call gribit(iget(328),lvls(1,iget(328)), grid1,im,jm)
      endif 

      if( iget(376) > 0 ) then ! water vapor channel brightness counts
       do j=jsta,jend
        do i=1,im
! convert to brightness value for direct comparison with NESDID products
! Formulation taken from NESDIS web site
! http://www.oso.noaa.gov/goes/goes-calibration/gvar-conversion.htm
         if(tb2(i,j)>163. .and. tb2(i,j)<=242.)then
	  grid1(i,j)=NINT(418.-tb2(i,j))*1.0
	 else if(tb2(i,j)>242. .and. tb2(i,j)<=330.)then
	  grid1(i,j)=NINT(660.-2.0*tb2(i,j))*1.0
	 else
	  grid1(i,j)=0.0 
	 end if  
        enddo
       enddo
       id(1:25) = 0
       id(02) = 129
       call gribit(iget(376),lvls(1,iget(376)), grid1,im,jm)
      endif  

    
      if( iget(329) > 0 ) then ! IR channel
       do j=jsta,jend
        do i=1,im
         grid1(i,j)=tb3(i,j)
        enddo
       enddo
       id(1:25) = 0
       id(02) = 129
       call gribit(iget(329),lvls(1,iget(329)), grid1,im,jm)
      endif  
    
      if( iget(377) > 0 ) then ! IR channel brightness counts
       do j=jsta,jend
        do i=1,im
! convert to brightness value for direct comparison with NESDID products
! Formulation taken from NESDIS web site
! http://www.oso.noaa.gov/goes/goes-calibration/gvar-conversion.htm
         if(tb3(i,j)>163. .and. tb3(i,j)<=242.)then
	  grid1(i,j)=NINT(418.-tb3(i,j))*1.0
	 else if(tb3(i,j)>242. .and. tb3(i,j)<=330.)then
	  grid1(i,j)=NINT(660.-2.0*tb3(i,j))*1.0
	 else
	  grid1(i,j)=0.0 
	 end if  
        enddo
       enddo
       id(1:25) = 0
       id(02) = 129
       call gribit(iget(377),lvls(1,iget(377)), grid1,im,jm)
      endif  


      if( iget(330) > 0 ) then
       do j=jsta,jend
        do i=1,im
         grid1(i,j)=tb4(i,j)
!         if(grid1(i,j)>400.)grid1(i,j)=400.
        enddo
       enddo
       id(1:25) = 0
       id(02) = 129
       call gribit(iget(330),lvls(1,iget(330)), grid1,im,jm)
      endif
    
! Deallocate arrays
      err1=0; err2=0; err3=0; err4=0
      err1 = crtm_destroy_atmosphere(atmosphere(1))
      err2 = crtm_destroy_surface(surface(1))
      err3 = crtm_destroy_rtsolution(rtsolution)
      if (err1/=success) write(6,*)'ERROR** destroy atmosphere.err=',&
     &      err1
      if (err2/=success) write(6,*)'ERROR** destroy surface..err=', &
     &      err2
      if (err3/=success) write(6,*)'ERROR** destroy rtsolution..err=', &
     &      err3
      deallocate (rtsolution)
!     

      endif ! for all iget logical
      return
      end
