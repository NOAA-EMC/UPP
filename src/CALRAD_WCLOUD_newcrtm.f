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
  !   11-02-06 Jun WANG   - addgrib2 option 
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
  use vrbls3d, only: o3, pint, pmid, t, q, qqw, qqi, qqr, f_rimef, nlice, qqs, qqg
  use vrbls2d, only: czen, ivgtyp, sno, pctsno, ths, vegfrc, si, u10h, v10h, u10,&
       v10, smstot, hbot, htop, cnvcfr
  use masks, only: gdlat, gdlon, sm, lmh, sice
  use soil, only:
  use gridspec_mod, only: gridtype
  use cmassi_mod, only: TRAD_ice
  use kinds, only: r_kind,r_single,i_kind
  use crtm_module, only: crtm_atmosphere_type,crtm_surface_type,crtm_geometry_type, &
       crtm_surface_create,o3_id,co2_id,wet_soil,crtm_forward,mass_mixing_ratio_units, &
       crtm_atmosphere_create,grass_scrub,grass_soil, meadow_grass,urban_concrete, &
       irrigated_low_vegetation,broadleaf_pine_forest,pine_forest,compacted_soil, &
       broadleaf_forest,broadleaf_brush,tundra,tilled_soil,scrub,scrub_soil,&
       crtm_options_type,crtm_destroy,crtm_init,SPECIFIC_AMOUNT_UNITS, &
       success,crtm_options_destroy,crtm_options_create, crtm_options_associated, &
       invalid_land
  use crtm_rtsolution_define, only: crtm_rtsolution_type, crtm_rtsolution_create, &
       crtm_rtsolution_destroy, crtm_rtsolution_associated 
  use crtm_spccoeff, only: sc
  use crtm_atmosphere_define, only:h2o_id,crtm_atmosphere_associated, &
       crtm_atmosphere_destroy,volume_mixing_ratio_units,crtm_atmosphere_zero
  use crtm_surface_define, only: crtm_surface_destroy, crtm_surface_associated, &
       crtm_surface_zero
  use crtm_channelinfo_define, only: crtm_channelinfo_type
  use crtm_parameters, only: limit_exp,toa_pressure,max_n_layers,MAX_SENSOR_SCAN_ANGLE
  use crtm_cloud_define, only:  water_cloud,ice_cloud,rain_cloud,snow_cloud,graupel_cloud,hail_cloud
  use message_handler, only: success,warning, display_message

  use params_mod, only: pi, rtd, p1000, capa, h1000, h1, g, rd, d608, qconv
  use rqstfld_mod, only: iget, id, lvls, iavblfld
  use ctlblk_mod, only: modelname, ivegsrc, novegtype, imp_physics, lm, spval, icu_physics,&
              grib, cfld, fld_info, datapd, idat, im, jsta, jend, jm
!     
  implicit none

  !     DECLARE VARIABLES.
  !     
  ! Mapping land surface type of GFS to CRTM
  !  Note: index 0 is water, and index 13 is ice. The two indices are not
  !        used and just assigned to COMPACTED_SOIL.
  integer, parameter, dimension(0:13) :: gfs_to_crtm=(/COMPACTED_SOIL,     &
         BROADLEAF_FOREST, BROADLEAF_FOREST, BROADLEAF_PINE_FOREST, PINE_FOREST, &
         PINE_FOREST, BROADLEAF_BRUSH, SCRUB, SCRUB, SCRUB_SOIL, TUNDRA,         &
         COMPACTED_SOIL, TILLED_SOIL, COMPACTED_SOIL/)

  ! Mapping land surface type of NMM to CRTM
  !  Note: index 16 is water, and index 24 is ice. The two indices are not
  !        used and just assigned to COMPACTED_SOIL.
  !       integer, parameter, dimension(24) :: nmm_to_crtm=(/URBAN_CONCRETE,       &
  !      &   COMPACTED_SOIL, IRRIGATED_LOW_VEGETATION, GRASS_SOIL, MEADOW_GRASS,   &
  !      &   MEADOW_GRASS, MEADOW_GRASS, SCRUB, GRASS_SCRUB, MEADOW_GRASS,         &
  !      &   BROADLEAF_FOREST, PINE_FOREST, BROADLEAF_FOREST, PINE_FOREST,         &
  !      &   BROADLEAF_PINE_FOREST, COMPACTED_SOIL, WET_SOIL, WET_SOIL,            &
  !      &   IRRIGATED_LOW_VEGETATION, TUNDRA, TUNDRA, TUNDRA, TUNDRA,             &
  !      &   COMPACTED_SOIL/)

  integer, allocatable:: nmm_to_crtm(:)
  integer, parameter:: ndat=100
  ! CRTM structure variable declarations.
  integer,parameter::  n_absorbers = 2
  !      integer,parameter::  n_clouds = 4 
  integer,parameter::  n_aerosols = 0
  ! Add your sensors here
  integer(i_kind),parameter:: n_sensors=7
  character(len=20),parameter,dimension(1:n_sensors):: sensorlist= &
      (/'imgr_g12            ', &
        'imgr_g11            ', &
	'amsre_aqua          ', &
	'tmi_trmm            ', &
	'ssmi_f15            ', &
	'ssmis_f20           ', &
        'ssmis_f17           '/)
  character(len=10),parameter,dimension(1:n_sensors):: obslist=  &
      (/'goes_img  ', &
        'goes_img  ', &
	'amsre     ', &
	'tmi       ', &
	'ssmi      ', &
	'ssmis     ', &
        'ssmis     '/)
!
  integer(i_kind) sensorindex
  integer(i_kind) lunin,nobs,nchanl,nreal
  integer(i_kind) error_status,itype
  integer(i_kind) err1,err2,err3,err4
  integer(i_kind) i,j,k,msig
  integer(i_kind) lcbot,lctop   !bsf
  integer jdn,ichan,ixchan,igot
  integer isat
  
  real(r_kind),parameter:: r100=100.0_r_kind
  real,parameter:: ozsmall = 1.e-10 ! to convert to mass mixing ratio
  real(r_kind) tsfc 
  real(r_kind),dimension(4):: sfcpct
  real(r_kind) snodepth,vegcover
  real snoeqv
  real snofrac
  real(r_kind),dimension(im,jsta:jend):: tb1,tb2,tb3,tb4
  real(r_kind),allocatable :: tb(:,:,:)
  real,dimension(im,jm):: grid1,grid2
  real sun_zenith,sun_azimuth, dpovg
  real sat_zenith
  real q_conv   !bsf
  real,parameter:: constoz = 604229.0_r_kind 
  real sublat,sublon
  real RHO,RHOX
  character(10)::obstype
  character(20)::isis
  
  logical hirs2,msu,goessndr,hirs3,hirs4,hirs,amsua,amsub,airs,hsb  &
            ,goes_img,mhs
  logical avhrr,avhrr_navy,lextra,ssu
  logical ssmi,ssmis,amsre,amsre_low,amsre_mid,amsre_hig,change
  logical ssmis_las,ssmis_uas,ssmis_env,ssmis_img
  logical sea,mixed,land,ice,snow,toss
  logical micrim,microwave
  !  logical,dimension(nobs):: luse
  logical, parameter :: debugprint = .false.
  type(crtm_atmosphere_type),dimension(1):: atmosphere
  type(crtm_surface_type),dimension(1) :: surface
  type(crtm_geometry_type),dimension(1) :: geometryinfo
  type(crtm_options_type),dimension(1)      :: options

  type(crtm_rtsolution_type),allocatable,dimension(:,:):: rtsolution
  type(crtm_channelinfo_type),allocatable,dimension(:) :: channelinfo
!     
  integer ii,jj,n_clouds,n
  integer,external :: iw3jdn
  !
  !*****************************************************************************
  ! Mapping land surface type of NMM to CRTM
  !      allocate(nmm_to_crtm(novegtype) )
  if(MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR')then 
     if(ivegsrc==1)then  !IGBP veg type
        allocate(nmm_to_crtm(novegtype) )
        nmm_to_crtm=(/PINE_FOREST, BROADLEAF_FOREST, PINE_FOREST,       &
             BROADLEAF_FOREST,BROADLEAF_PINE_FOREST, SCRUB, SCRUB_SOIL, &
             BROADLEAF_BRUSH,BROADLEAF_BRUSH, SCRUB, BROADLEAF_BRUSH,   &
             TILLED_SOIL, URBAN_CONCRETE,TILLED_SOIL, INVALID_LAND,     &
             COMPACTED_SOIL, INVALID_LAND, TUNDRA,TUNDRA, TUNDRA/)
     else if(ivegsrc==0)then ! USGS veg type
        allocate(nmm_to_crtm(novegtype) )
        nmm_to_crtm=(/URBAN_CONCRETE,       &
             COMPACTED_SOIL, IRRIGATED_LOW_VEGETATION, GRASS_SOIL, MEADOW_GRASS,   &
             MEADOW_GRASS, MEADOW_GRASS, SCRUB, GRASS_SCRUB, MEADOW_GRASS,         &
             BROADLEAF_FOREST, PINE_FOREST, BROADLEAF_FOREST, PINE_FOREST,         &
             BROADLEAF_PINE_FOREST, COMPACTED_SOIL, WET_SOIL, WET_SOIL,            &
             IRRIGATED_LOW_VEGETATION, TUNDRA, TUNDRA, TUNDRA, TUNDRA,             &
             COMPACTED_SOIL/)
     else
        print*,'novegtype=',novegtype
        print*,'model veg type not supported by post in calling crtm ' 
        print*,'skipping generation of simulated radiance' 
        return
     end if 
  end if 
      
  !     START SUBROUTINE CALRAD.
  ifactive: if (iget(327) > 0 .or. iget(328) > 0 .or. iget(329) > 0       &
       .or. iget(330) > 0 .or. iget(446) > 0 .or. iget(447) > 0  & 
       .or. iget(448) > 0 .or. iget(449) > 0  .or. iget(456) > 0   &
       .or. iget(457) > 0 .or. iget(458) > 0 .or. iget(459) > 0  &
       .or. iget(460) > 0 .or. iget(461) > 0 .or. iget(462) > 0  &
       .or. iget(463) > 0 .or. iget(483) > 0 .or. iget(484) > 0  &
       .or. iget(485) > 0 .or. iget(486) > 0 .or. iget(488) > 0  &
       .or. iget(489) > 0 .or. iget(490) > 0 .or. iget(491) > 0  &
       .or. iget(492) > 0 .or. iget(493) > 0 .or. iget(494) > 0  &
       .or. iget(495) > 0 .or. iget(496) > 0 .or. iget(497) > 0  &
       .or. iget(498) > 0 .or. iget(499) > 0 .or. iget(800) > 0  &
       .or. iget(801) > 0 .or. iget(802) > 0 .or. iget(803) > 0  &
       .or. iget(804) > 0 .or. iget(805) > 0 .or. iget(806) > 0  &
       .or. iget(807) > 0) then
     ! specify numbers of cloud species    
     if(imp_physics==99)then ! Zhao Scheme
        n_clouds=2 ! GFS uses Zhao scheme
     else if(imp_physics==5 .or. imp_physics==85 .or. imp_physics==95)then
        n_clouds=6  ! change to 6 cloud types because microwave is sensitive to density
     else if(imp_physics==8 .or. imp_physics==6 .or. imp_physics==2)then
        n_clouds=5
     end if
 
     ! Initialize debug print gridpoint index to middle of tile:
     ii=im/2
     jj=(jsta+jend)/2

     ! Initialize ozone to zeros for WRF NMM and ARW for now
     if (MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR')o3=0.0
     ! Compute solar zenith angle for GFS
     if (MODELNAME /= 'NMM')then
        jdn=iw3jdn(idat(3),idat(1),idat(2))
	do j=jsta,jend
	   do i=1,im
	      call zensun(jdn,float(idat(4)),gdlat(i,j),gdlon(i,j)       &
      	                  ,pi,sun_zenith,sun_azimuth)
              czen(i,j)=cos(sun_zenith)
	   end do
	end do
	if(jj>=jsta .and. jj<=jend)                                  &
            print*,'sample GFS zenith angle=',acos(czen(ii,jj))*rtd   
     end if	       
     ! Initialize CRTM.  Load satellite sensor array.
     ! The optional arguments Process_ID and Output_Process_ID limit
     ! generation of runtime informative output to mpi task
     ! Output_Process_ID (which here is set to be task 0)
     print*,'success in CALRAD= ',success
     allocate( channelinfo(n_sensors))

     error_status = crtm_init(sensorlist,channelinfo,   &
          Process_ID=0,Output_Process_ID=0 )
     print*, 'channelinfo after init= ',channelinfo(1)%sensor_id, &
              channelinfo(2)%sensor_id
     if (error_status /= 0_i_kind)                                  &
         write(6,*)'ERROR*** crtm_init error_status=',error_status
  
     ! Discard all ssmis_f17 channels except the four we output:
     call select_channels(channelinfo(7),4,(/ 15,16,17,18 /))

     !   lunin=1 ! will read data file in the future, only simulate GOES for now
     !   open(lunin,file='obs_setup',form='unformatted') ! still need to find out filename
     !   rewind lunin
  
     ! Loop over data types to process    
     sensordo: do isat=1,n_sensors
        !    read(lunin,end=125) obstype,isis,nreal,nchanl
        ! Test GOES for now.
        !       obstype='goes_img'
        !       isis='imgr_g12'
        obstype=obslist(isat) 
        isis=trim(sensorlist(isat))
        !       print*,'obstype, isis= ',obstype,isis
        ! only call crtm forward model for this sensor if requested in the
        ! post control file
        sensor_avail: if((isis=='imgr_g12' .and. (iget(327) > 0 .or. iget(328) > 0 &
             .or. iget(329) > 0 .or. iget(330) > 0 .or. iget(456) > 0   &
             .or. iget(457) > 0 .or. iget(458) > 0 .or. iget(459) > 0 )) .OR. &
             (isis=='imgr_g11' .and. (iget(446) > 0 .or. iget(447) > 0 &
             .or. iget(448) > 0 .or. iget(449) > 0 .or. iget(460) > 0   &
             .or. iget(461) > 0 .or. iget(462) > 0 .or. iget(463) > 0)) .OR. &
             (isis=='amsre_aqua' .and. (iget(483) > 0 .or. iget(484) > 0  &
             .or. iget(485) > 0 .or. iget(486) > 0)) .OR. &
             (isis=='tmi_trmm' .and. (iget(488) > 0 .or. iget(489) > 0  &
             .or. iget(490) > 0 .or. iget(491) > 0)) .OR. &
             (isis=='ssmi_f15' .and. (iget(492) > 0 .or. iget(493) > 0  &
             .or. iget(494) > 0 .or. iget(495) > 0)) .OR. &
             (isis=='ssmis_f20' .and. (iget(496) > 0 .or. iget(497) > 0  &
             .or. iget(498) > 0 .or. iget(499) > 0)) .OR. &
             (isis=='ssmis_f17' .and. (iget(800) > 0 .or. iget(801) > 0  &
             .or. iget(802) > 0 .or. iget(803) > 0 .or. iget(804) > 0 &
             .or. iget(805) > 0 .or. iget(806) > 0 .or. iget(807) > 0)) )then
           print*,'obstype, isis= ',obstype,isis
           !       isis='amsua_n15'

           ! Initialize logical flags for satellite platform

           hirs2      = obstype == 'hirs2'
           hirs3      = obstype == 'hirs3'
           hirs4      = obstype == 'hirs4'
           hirs       = hirs2 .or. hirs3 .or. hirs4
           msu        = obstype == 'msu'
           ssu        = obstype == 'ssu'
           goessndr   = obstype == 'sndr'  .or. obstype == 'sndrd1' .or.    &
                          obstype == 'sndrd2'.or. obstype == 'sndrd3' .or.  &
                          obstype == 'sndrd4'
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

           ssmis=ssmis_las.or.ssmis_uas.or.ssmis_img.or.ssmis_env.or.ssmis

           micrim=ssmi .or. ssmis .or. amsre   ! only used for MW-imager-QC and id_qc(ch)

           microwave=amsua .or. amsub .or. mhs .or. msu .or. hsb .or. micrim

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
           allocate(tb(im,jsta:jend,channelinfo(sensorindex)%n_channels))
           err1=0; err2=0; err3=0; err4=0
           if(lm > max_n_layers)then
              write(6,*) 'CALRAD: lm > max_n_layers - '//                 &
      	                 'increase crtm max_n_layers ',                   & 
                         lm,max_n_layers
              stop 2
           end if
           CALL crtm_atmosphere_create(atmosphere(1),lm,n_absorbers,n_clouds &
                        ,n_aerosols)
           CALL crtm_surface_create(surface(1),channelinfo(sensorindex)%n_channels)
           CALL crtm_rtsolution_create(rtsolution,lm)
           if (.NOT.(crtm_atmosphere_associated(atmosphere(1)))) &
               write(6,*)' ***ERROR** creating atmosphere.'
           if (.NOT.(crtm_surface_associated(surface(1)))) &
               write(6,*)' ***ERROR** creating surface.'
           if (.NOT.(ANY(crtm_rtsolution_associated(rtsolution)))) &
               write(6,*)' ***ERROR** creating rtsolution.'

           atmosphere(1)%n_layers = lm
           !      atmosphere(1)%level_temperature_input = 0
                  atmosphere(1)%absorber_id(1) = H2O_ID
                  atmosphere(1)%absorber_id(2) = O3_ID
                  atmosphere(1)%absorber_units(1) = MASS_MIXING_RATIO_UNITS
                  atmosphere(1)%absorber_units(2) = VOLUME_MIXING_RATIO_UNITS
           !       atmosphere(1)%absorber_units(2) = MASS_MIXING_RATIO_UNITS ! Use mass mixing ratio
                  atmosphere(1)%level_pressure(0) = TOA_PRESSURE
           ! Define Clouds
           if(imp_physics==99)then
              atmosphere(1)%cloud(1)%n_layers = lm
              atmosphere(1)%cloud(1)%Type = WATER_CLOUD
              atmosphere(1)%cloud(2)%n_layers = lm
              atmosphere(1)%cloud(2)%Type = ICE_CLOUD
           else if(imp_physics==5 .or. imp_physics==85 .or. imp_physics==95)then
              atmosphere(1)%cloud(1)%n_layers = lm
              atmosphere(1)%cloud(1)%Type = WATER_CLOUD
              atmosphere(1)%cloud(2)%n_layers = lm
              atmosphere(1)%cloud(2)%Type = ICE_CLOUD
              atmosphere(1)%cloud(3)%n_layers = lm
              atmosphere(1)%cloud(3)%Type = RAIN_CLOUD
              atmosphere(1)%cloud(4)%n_layers = lm
              atmosphere(1)%cloud(4)%Type = SNOW_CLOUD
              atmosphere(1)%cloud(5)%n_layers = lm
              atmosphere(1)%cloud(5)%Type = GRAUPEL_CLOUD
      	      atmosphere(1)%cloud(6)%n_layers = lm
              atmosphere(1)%cloud(6)%Type = HAIL_CLOUD
           else if(imp_physics==8 .or. imp_physics==6 .or. imp_physics==2)then
              atmosphere(1)%cloud(1)%n_layers = lm
              atmosphere(1)%cloud(1)%Type = WATER_CLOUD
              atmosphere(1)%cloud(2)%n_layers = lm
              atmosphere(1)%cloud(2)%Type = ICE_CLOUD
              atmosphere(1)%cloud(3)%n_layers = lm
              atmosphere(1)%cloud(3)%Type = RAIN_CLOUD
              atmosphere(1)%cloud(4)%n_layers = lm
              atmosphere(1)%cloud(4)%Type = SNOW_CLOUD
              atmosphere(1)%cloud(5)%n_layers = lm
              atmosphere(1)%cloud(5)%Type = GRAUPEL_CLOUD
           end if        

           !    if(nchanl /= channelinfo(sensorindex)%n_channels) write(6,*)'***ERROR** nchanl,n_channels ', &
           !           nchanl,channelinfo(sensorindex)%n_channels

           ! Load surface sensor data structure
           surface(1)%sensordata%n_channels = channelinfo(sensorindex)%n_channels
           surface(1)%sensordata%sensor_id  = channelinfo(sensorindex)%sensor_id
           surface(1)%sensordata%WMO_sensor_id = channelinfo(sensorindex)%WMO_sensor_id
           surface(1)%sensordata%WMO_Satellite_id = channelinfo(sensorindex)%WMO_Satellite_id
           surface(1)%sensordata%sensor_channel = channelinfo(sensorindex)%sensor_channel

           ! Loop through all grid points
           !
           ! First compute nadir simulated radiance because it is what is being computed operationally
           nadir: if (iget(327) > 0 .or. iget(328) > 0 .or. iget(329) > 0       &
                .or. iget(330) > 0 .or. iget(446) > 0 .or. iget(447) > 0  & 
                .or. iget(448) > 0 .or. iget(449) > 0 .or. iget(483) > 0  & 
                .or. iget(484) > 0 .or. iget(485) > 0 .or. iget(486) > 0  &
                .or. iget(488) > 0 .or. iget(489) > 0 .or. iget(490) > 0  &
                .or. iget(491) > 0 .or. iget(492) > 0 .or. iget(493) > 0  &
                .or. iget(494) > 0 .or. iget(495) > 0 .or. iget(496) > 0  &
                .or. iget(497) > 0 .or. iget(498) > 0 .or. iget(499) > 0  &
                .or. iget(800) > 0 .or. iget(801) > 0 .or. iget(802) > 0  &
                .or. iget(803) > 0) then
              do j=jsta,jend
                 do i=1,im

                    !    Load geometry structure
                    !    geometryinfo(1)%sensor_zenith_angle = zasat*rtd  ! local zenith angle ???????
                    geometryinfo(1)%sensor_zenith_angle=0.
                    geometryinfo(1)%sensor_scan_angle=0.

                    !only call crtm if we have right satellite zenith angle
                    IF(geometryinfo(1)%sensor_zenith_angle <= MAX_SENSOR_SCAN_ANGLE &
                         .and. geometryinfo(1)%sensor_zenith_angle >= 0.0_r_kind)THEN
                       geometryinfo(1)%source_zenith_angle = acos(czen(i,j))*rtd ! solar zenith angle
                       geometryinfo(1)%sensor_scan_angle   = 0. ! scan angle, assuming nadir
                       if(i==ii.and.j==jj)print*,'sample geometry ',                   &
                               geometryinfo(1)%sensor_zenith_angle                     &
                              ,geometryinfo(1)%source_zenith_angle                     &
                              ,czen(i,j)*rtd 
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
                                                     snoeqv,itype,IVGTYP(I,J)
                          if(sm(i,j) > 0.1)then
                             sfcpct(4)=0.
                          else 
                             call snfrac_gfs(SNOeqv,IVGTYP(I,J),snofrac)
                             sfcpct(4)=snofrac
                          end if
                          if(i==ii.and.j==jj)print*,'sno,itype,ivgtyp,sfcpct(4) = ',     &
                                                     snoeqv,itype,IVGTYP(I,J),sfcpct(4)
                       else if(MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR')then
                          sfcpct(4)=pctsno(i,j)
                       else          
                          itype=IVGTYP(I,J)
                          IF(itype == 0)itype=8
                          CALL SNFRAC (SNO(I,J),IVGTYP(I,J),snofrac)
	                  sfcpct(4)=snofrac
                       end if 
                       !	CALL SNFRAC (SNO(I,J),IVGTYP(I,J),snofrac)
                       !	sfcpct(4)=snofrac
	               if(sm(i,j) > 0.1)then ! water
                          !  	 tsfc=sst(i,j)
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
                       ! Chuang: for igbp type 15 (snow/ice), the main type needs to be set to ice or snow
                       ! to prevent crtm forward model from failing	
                       if(ivegsrc==1 .and. itype==15 .and. sfcpct(4)<1.0_r_kind)then
                          if(debugprint)print*,'changing land type for veg type 15',i,j,itype,sfcpct(1:4)
	                  sfcpct(1)=0.0_r_kind
	                  sfcpct(2)=0.0_r_kind
	                  sfcpct(3)=0.0_r_kind
	                  sfcpct(4)=1.0_r_kind
                          !print*,'change main land type to snow for veg type 15 ',i,j
	               end if 

                       sea  = sfcpct(1)  >= 0.99_r_kind
                       land = sfcpct(2)  >= 0.99_r_kind
                       ice  = sfcpct(3)  >= 0.99_r_kind
                       snow = sfcpct(4)  >= 0.99_r_kind
                       mixed = .not. sea  .and. .not. ice .and.     &
                                .not. land .and. .not. snow
                       if((sfcpct(1)+sfcpct(2)+sfcpct(3)+sfcpct(4)) >1._r_kind) &
                           print*,'ERROR sfcpct ',i,j,sfcpct(1)  &
                              ,sfcpct(2),sfcpct(3),sfcpct(4)
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

                       if(gridtype=='B' .or. gridtype=='E')then
                          surface(1)%wind_speed         = sqrt(u10h(i,j)*u10h(i,j)   &
                                                              +v10h(i,j)*v10h(i,j))
                       else
                          surface(1)%wind_speed         = sqrt(u10(i,j)*u10(i,j)   &
                                                              +v10(i,j)*v10(i,j))
                       end if
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
                       ! Debug print
                       if(debugprint)then       
                          if(surface(1)%wind_speed<0. .or. surface(1)%wind_speed>200.)  &
                             print*,'bad 10 m wind'
                          if(surface(1)%water_coverage<0. .or. surface(1)%water_coverage>1.) &
                             print*,'bad water coverage'
                          if(surface(1)%land_coverage<0. .or. surface(1)%land_coverage>1.)  &
                             print*,'bad land coverage'
                          if(surface(1)%ice_coverage<0. .or. surface(1)%ice_coverage>1.)  &
                             print*,'bad ice coverage'
                          if(surface(1)%snow_coverage<0. .or. surface(1)%snow_coverage>1.)  &
                             print*,'bad snow coverage'
                          if(surface(1)%land_temperature<0. .or. surface(1)%land_temperature>350.)  &
                             print*,'bad land T'
                          if(surface(1)%soil_moisture_content<0. .or. surface(1)%soil_moisture_content>600.) &
                             print*,'bad soil_moisture_content'
                          if(surface(1)%vegetation_fraction<0. .or. surface(1)%vegetation_fraction>1.) &
                             print*,'bad vegetation cover'
                          if(surface(1)%snow_depth<0. .or.  surface(1)%snow_depth>10000.) &
                             print*,'bad snow_depth'
                          if(MODELNAME == 'GFS' .and. (itype<0 .or. itype>13)) &
                             print*,'bad veg type'
                       end if
       
                       if(i==ii.and.j==jj)print*,'sample surface in CALRAD=', &
                             i,j,surface(1)%wind_speed,surface(1)%water_coverage,       &
                             surface(1)%land_coverage,surface(1)%ice_coverage,          &
                             surface(1)%snow_coverage,surface(1)%land_temperature,      &
                             surface(1)%snow_temperature,surface(1)%water_temperature,  &
                             surface(1)%ice_temperature,surface(1)%vegetation_fraction, &
                             surface(1)%soil_temperature,surface(1)%snow_depth,         &
                             surface(1)%land_type,sm(i,j)

                       !       Load profiles into model layers

                       !       Load atmosphere profiles into RTM model layers
                       !       CRTM counts from top down just as post does
                       if(i==ii.and.j==jj)print*,'TOA= ',atmosphere(1)%level_pressure(0)
                       do k = 1,lm
                          atmosphere(1)%level_pressure(k) = pint(i,j,k+1)/r100
                          atmosphere(1)%pressure(k)       = pmid(i,j,k)/r100
                          atmosphere(1)%temperature(k)    = t(i,j,k)
                          atmosphere(1)%absorber(k,1)     = max(0.  &
                                          ,q(i,j,k)*h1000/(h1-q(i,j,k))) ! use mixing ratio like GSI
                          atmosphere(1)%absorber(k,2)     = max(ozsmall,o3(i,j,k)*constoz)
                          !        atmosphere(1)%absorber(k,2)     = max(ozsmall,o3(i,j,k)*h1000) ! convert to g/kg
                          ! fill in cloud mixing ratio later  
                          if(debugprint)then
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
                          end if
                          if(i==ii.and.j==jj)print*,'sample atmosphere in CALRAD=',  &
   	                        i,j,k,atmosphere(1)%level_pressure(k),atmosphere(1)%pressure(k),  &
                                atmosphere(1)%temperature(k),atmosphere(1)%absorber(k,1),  &
                                atmosphere(1)%absorber(k,2)
                          ! Specify clouds
                          dpovg=(pint(i,j,k+1)-pint(i,j,k))/g !crtm uses column integrated field
                          if(imp_physics==99)then
                             atmosphere(1)%cloud(1)%effective_radius(k) = 10.
	                        atmosphere(1)%cloud(1)%water_content(k) = max(0.,qqw(i,j,k)*dpovg)
                             ! GFS uses temperature and ice concentration dependency formulation to 
                             ! determine effetive radis for cloud ice
                             ! since GFS does not output ice concentration yet, use default 50 um
	                        atmosphere(1)%cloud(2)%effective_radius(k) = 50.	 
	                        atmosphere(1)%cloud(2)%water_content(k) = max(0.,qqi(i,j,k)*dpovg)
                             if(debugprint)then
                                if(atmosphere(1)%cloud(1)%water_content(k)<0. .or.   &
                                        atmosphere(1)%cloud(1)%water_content(k)>1.)  &
                                   print*,'bad atmosphere cloud water'
                                if(atmosphere(1)%cloud(2)%water_content(k)<0. .or.   &
                                        atmosphere(1)%cloud(2)%water_content(k)>1.)  &
                                   print*,'bad atmosphere cloud ice'
                             end if
                          else if(imp_physics==5 .or. imp_physics==85 .or. imp_physics==95)then
	                     atmosphere(1)%cloud(1)%effective_radius(k) = 10.
	                     atmosphere(1)%cloud(1)%water_content(k) = max(0.,qqw(i,j,k)*dpovg)
	                     atmosphere(1)%cloud(2)%effective_radius(k) = 25.
	                     atmosphere(1)%cloud(2)%water_content(k) = max(0.,qqi(i,j,k)*dpovg)
	                     atmosphere(1)%cloud(3)%effective_radius(k) = 200.
	                     atmosphere(1)%cloud(3)%water_content(k) = max(0.,qqr(i,j,k)*dpovg)
                             !	 atmosphere(1)%cloud(4)%effective_radius(k) = 250.
	                     RHO=pmid(i,j,k)/(RD*T(I,J,K)*(1.+D608*Q(I,J,K)))
	                     if(F_RimeF(i,j,k)<=5.0)then
	                        RHOX=100
	                        if(NLICE(I,J,K)>0.) &
	                           atmosphere(1)%cloud(4)%effective_radius(k) =   &
                                   1.0E6*1.5*(RHO*qqs(i,j,k)/(PI*RHOX*NLICE(I,J,K)))**(1./3.) !convert to microns
	                        atmosphere(1)%cloud(4)%water_content(k) = max(0.,qqs(i,j,k)*dpovg)
	                        atmosphere(1)%cloud(5)%effective_radius(k) = 0.
	                        atmosphere(1)%cloud(5)%water_content(k) =0.
	                        atmosphere(1)%cloud(6)%effective_radius(k) = 0.
	                        atmosphere(1)%cloud(6)%water_content(k) =0.
	                     else if(F_RimeF(i,j,k)<=20.0)then 
	                        atmosphere(1)%cloud(4)%effective_radius(k) = 0.
	                        atmosphere(1)%cloud(4)%water_content(k) =0.
	                        RHOX=400.
                                if(NLICE(I,J,K)>0.) &
	                           atmosphere(1)%cloud(5)%effective_radius(k) =   &
                                         1.0E6*1.5*(RHO*qqs(i,j,k)/(PI*RHOX*NLICE(I,J,K)))**(1./3.)
	                        atmosphere(1)%cloud(5)%water_content(k) =max(0.,qqs(i,j,k)*dpovg)
	                        atmosphere(1)%cloud(6)%effective_radius(k) = 0.
	                        atmosphere(1)%cloud(6)%water_content(k) =0.
	                     else
	                        atmosphere(1)%cloud(4)%effective_radius(k) = 0.
	                        atmosphere(1)%cloud(4)%water_content(k) =0.
	                        atmosphere(1)%cloud(5)%effective_radius(k) = 0.
	                        atmosphere(1)%cloud(5)%water_content(k) =0.
	                        RHOX=900.
                                if(NLICE(I,J,K)>0.) &
	                           atmosphere(1)%cloud(6)%effective_radius(k) =  &
                                       1.0E6*1.5*(RHO*qqs(i,j,k)/(PI*RHOX*NLICE(I,J,K)))**(1./3.)
	                           atmosphere(1)%cloud(6)%water_content(k) =max(0.,qqs(i,j,k)*dpovg)
	                     end if  
                             if(debugprint .and. i==im/2 .and. j==jsta)   &
                                print*,'sample precip ice radius= ',i,j,k, F_RimeF(i,j,k), &
	                        atmosphere(1)%cloud(4)%effective_radius(k), atmosphere(1)%cloud(4)%water_content(k), &
	                        atmosphere(1)%cloud(5)%effective_radius(k), atmosphere(1)%cloud(5)%water_content(k), &
	                        atmosphere(1)%cloud(6)%effective_radius(k), atmosphere(1)%cloud(6)%water_content(k)
	
                          else if(imp_physics==8 .or. imp_physics==6 .or. imp_physics==2)then
                             atmosphere(1)%cloud(1)%effective_radius(k) = 10.
                             atmosphere(1)%cloud(1)%water_content(k) = max(0.,qqw(i,j,k)*dpovg)
                             atmosphere(1)%cloud(2)%effective_radius(k) = 25.
                             atmosphere(1)%cloud(2)%water_content(k) = max(0.,qqi(i,j,k)*dpovg)
                             atmosphere(1)%cloud(3)%effective_radius(k) = 200.
                             atmosphere(1)%cloud(3)%water_content(k) = max(0.,qqr(i,j,k)*dpovg)
                             atmosphere(1)%cloud(4)%effective_radius(k) = 250.
                             atmosphere(1)%cloud(4)%water_content(k) = max(0.,qqs(i,j,k)*dpovg)
                             atmosphere(1)%cloud(5)%effective_radius(k) = 350.
                             atmosphere(1)%cloud(5)%water_content(k) = max(0.,qqg(i,j,k)*dpovg)
                          end if 
                       end do

                       !bsf - start
                       !-- Add subgrid-scale convective clouds for WRF runs
                       if(icu_physics==2) then
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
                                   atmosphere(1)%cloud(2)%water_content(k) + dpovg*q_conv
                                else
                                   atmosphere(1)%cloud(1)%water_content(k) =   &
                                   atmosphere(1)%cloud(1)%water_content(k) + dpovg*q_conv
                                endif
                             end do   !-- do k = lctop,lcbot
                          endif      !-- if (lcbot-lctop > 1) then
                       endif        !-- if (MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR') then
                       !bsf - end

                       !     call crtm forward model
                       error_status = crtm_forward(atmosphere,surface,                    &
                                      geometryinfo,channelinfo(sensorindex:sensorindex),  &
                                      rtsolution)
                       if (error_status /=0) then
                          print*,'***ERROR*** during crtm_forward call ', error_status
                          do n=1,channelinfo(sensorindex)%n_channels
                             tb(i,j,n)=spval
                          end do 
                          !         tb2(i,j)=spval
                          !         tb3(i,j)=spval
                          !         tb4(i,j)=spval
                       else 	 
                          do n=1,channelinfo(sensorindex)%n_channels
                             tb(i,j,n)=rtsolution(n,1)%brightness_temperature
                          end do 
                          !        tb1(i,j)=rtsolution(1,1)%brightness_temperature
                          !        tb2(i,j)=rtsolution(2,1)%brightness_temperature
                          !        tb3(i,j)=rtsolution(3,1)%brightness_temperature	 
                          !        tb4(i,j)=rtsolution(4,1)%brightness_temperature
                          if(i==ii.and.j==jj) then
                             do n=1,channelinfo(sensorindex)%n_channels
 3301                           format('Sample rtsolution(',I0,',',I0,') in CALRAD = ',F0.3)
                                print 3301,n,1,rtsolution(n,1)%brightness_temperature
                             enddo
                             do n=1,channelinfo(sensorindex)%n_channels
 3302                           format('Sample tb(',I0,',',I0,',',I0,') in CALRAD = ',F0.3)
                                print 3302,ii,jj,n,tb(ii,jj,n)
                             enddo
                          endif
                          !        if(tb1(i,j) < 400. )  &
                          !     &        print*,'good tb1 ',i,j,tb1(i,j),gdlat(i,j),gdlon(i,j)
                          !        if(tb2(i,j) > 400.)print*,'bad tb2 ',i,j,tb2(i,j)
                          !        if(tb3(i,j) > 400.)print*,'bad tb3 ',i,j,tb3(i,j)
                          !        if(tb4(i,j) > 400.)print*,'bad tb4 ',i,j,tb4(i,j)
                       end if
                    else
                       do n=1,channelinfo(sensorindex)%n_channels
                          tb(i,j,n)=spval
                       end do
                       !       tb1(i,j)=spval
                       !       tb2(i,j)=spval
                       !       tb3(i,j)=spval
                       !       tb4(i,j)=spval
                    END IF ! endif block for allowable satellite zenith angle 
                 end do ! end loop for i
              end do ! end loop for j 
  
              !      error_status = crtm_destroy(channelinfo)
              !      if (error_status /= success) &
              !     &   print*,'ERROR*** crtm_destroy error_status=',error_status
  
              if (isis=='imgr_g12')then  ! writing goes 12 to grib
                 do ixchan=1,4   ! write brightness temperatures
                    ichan=ixchan
                    igot=iget(326+ixchan)
                    if(igot>0) then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,ichan)
                           enddo
                        enddo
                        id(1:25) = 0
                        id(02) = 129
                        if (grib=="grib1") then
                           call gribit(igot,lvls(1,igot), grid1,im,jm)
                        else
                           cfld=cfld+1
                           fld_info(cfld)%ifld=IAVBLFLD(igot)
                           datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                        endif
                    endif
                 enddo

                 do ixchan=1,2  ! write brightness counts for only two channels
                   ichan=1+ixchan
                   igot=iget(375+ixchan)
                   if(igot>0) then
                     do j=jsta,jend
                       do i=1,im
                         ! convert to brightness value for direct comparison with NESDID products
                         ! Formulation taken from NESDIS web site
                         ! http://www.oso.noaa.gov/goes/goes-calibration/gvar-conversion.htm
                         if(tb(i,j,ichan)>163. .and. tb(i,j,ichan)<=242.)then
                           grid1(i,j)=NINT(418.-tb(i,j,ichan))*1.0
                         else if(tb(i,j,ichan)>242. .and. tb(i,j,ichan)<=330.)then
                           grid1(i,j)=NINT(660.-2.0*tb(i,j,ichan))*1.0
                         else
                           grid1(i,j)=0.0
                         end if
                       enddo
                     enddo
                     id(1:25) = 0
                     id(02) = 129
                     if (grib=="grib1") then
                       call gribit(igot,lvls(1,igot), grid1,im,jm)
                     else
                       cfld=cfld+1
                       fld_info(cfld)%ifld=IAVBLFLD(igot)
                       datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                     endif
                   endif
                 enddo
              end if  ! end of outputting goes 12 

              if (isis=='imgr_g11')then  ! writing goes 11 to grib
                 do ixchan=1,4
                    ichan=ixchan
                    igot=445+ixchan
                    if(igot>0) then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,ichan)
                          enddo
                       enddo
                       id(1:25) = 0
                       id(02) = 130
                       if (grib=="grib1") then
                          call gribit(igot,lvls(1,igot), grid1,im,jm)
                       else
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif ! IGOT
                 enddo
              end if  ! end of outputting goes 11
      
              if (isis=='amsre_aqua')then  ! writing amsre to grib (37 & 89 GHz)
                 do ixchan=1,4
                    ichan=8+ixchan
                    igot=iget(482+ixchan)
                    if(igot>0) then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,ichan)
                          enddo
                       enddo
                       id(1:25) = 0
                       id(02) = 133
                       if (grib=="grib1") then
                          call gribit(igot,lvls(1,igot), grid1,im,jm)
                       else
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              end if  ! end of outputting amsre
      
              if (isis=='tmi_trmm')then  ! writing trmm to grib (37 & 85.5 GHz)
                 do ixchan=1,4
                    ichan=5+ixchan
                    igot=iget(487+ixchan)
                    if(igot>0) then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,ichan)
                          enddo
                       enddo
                       id(1:25) = 0
                       id(02) = 133
                       if (grib=="grib1") then
                          call gribit(igot,lvls(1,igot), grid1,im,jm)
                       else
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              end if  ! end of outputting trmm
      
              if (isis=='ssmi_f15')then  ! writing ssmi to grib (37 & 85 GHz)
                 do ixchan=1,4
                    ichan=3+ixchan
                    igot=iget(491+ixchan)
                    if(igot>0) then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,ichan)
                          enddo
                       enddo
                       id(1:25) = 0
                       id(02) = 133
                       if (grib=="grib1") then
                          call gribit(igot,lvls(1,igot), grid1,im,jm)
                       else
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              end if  ! end of outputting ssmi
      
              if (isis=='ssmis_f20')then  ! writing ssmi to grib (37 & 91 GHz)
                 do ixchan=1,4
                    ichan=14+ixchan
                    igot=iget(495+ixchan)
                    if(igot>0) then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,ichan)
                          enddo
                       enddo
                       id(1:25) = 0
                       id(02) = 133
                       if (grib=="grib1") then
                          call gribit(igot,lvls(1,igot), grid1,im,jm)
                       else
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              end if  ! end of outputting ssmis

              if(isis=='ssmis_f17') then ! writing f17 ssmis to grib
                 do ixchan=1,4
                    ichan=ixchan ! using select_channels, we discard channels 1-14 and 19-24
                    igot=iget(799+ixchan) ! iget(800) ... iget(803)
                    if(igot > 0) then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,ichan)
                          enddo
                       enddo
                       id(1:25) = 0
                       id(02) = 2
                       id(09) = 112
                       id(10) = 167
                       id(11) = ichan+14
                       if (grib=="grib1") then
                          call gribit(igot,lvls(1,igot), grid1,im,jm)
                       else
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif  ! IGOT
                 enddo
              endif ! end of outputting f17 ssmis

           end if nadir ! end if for computing nadir simulated radiance    

           ! rerun crtm again for users who wish to factor in satellite zenith angle      
           nonnadir: if (iget(456) > 0 .or. iget(457) > 0 .or. iget(458) > 0  &
                   .or. iget(459) > 0 .or. iget(460) > 0 .or. iget(461) > 0         &
                   .or. iget(462) > 0 .or. iget(463) > 0 .or. iget(804) > 0         &
                   .or. iget(805) > 0 .or. iget(806) > 0 .or. iget(807) > 0) then
              do j=jsta,jend
                 do i=1,im

                    !    Load geometry structure
                    !    geometryinfo(1)%sensor_zenith_angle = zasat*rtd  ! local zenith angle ???????
                    ! compute satellite zenith angle
                    if(isis=='imgr_g12')then
                       sublat=0.0
                       sublon=-75.0
                    else if(isis=='imgr_g11')then
                       sublat=0.0
                       sublon=-135.0
                    end if

                    if(isis=='ssmis_f17') then
                       ! Use a constant zenith angle of 53 degrees for F17 SSMIS:
                       sat_zenith=53.0
                    else
                       ! For other imagers (GOES-11 and 12), calculate based on satellite location:
                       call GEO_ZENITH_ANGLE(i,j,gdlat(i,j),gdlon(i,j)  &
                            ,sublat,sublon,sat_zenith)
                    endif

                    geometryinfo(1)%sensor_zenith_angle=sat_zenith
	            geometryinfo(1)%sensor_scan_angle=sat_zenith

                    if(i==ii .and. j==jj) then
                       print *,'zenith info: zenith=',sat_zenith,' scan=',sat_zenith, &
                             ' MAX_SENSOR_SCAN_ANGLE=',MAX_SENSOR_SCAN_ANGLE
                    endif

                    !        geometryinfo(1)%sensor_zenith_angle = 0. ! 44.
                    !only call crtm if we have right satellite zenith angle
                    IF(geometryinfo(1)%sensor_zenith_angle <= MAX_SENSOR_SCAN_ANGLE &
                         .and. geometryinfo(1)%sensor_zenith_angle >= 0.0_r_kind)THEN
                       geometryinfo(1)%source_zenith_angle = acos(czen(i,j))*rtd ! solar zenith angle
                       geometryinfo(1)%sensor_scan_angle   = 0. ! scan angle, assuming nadir
                       if(i==ii.and.j==jj)print*,'sample geometry ',                   &
                          geometryinfo(1)%sensor_zenith_angle                          &
                          ,geometryinfo(1)%source_zenith_angle                         &
                          ,czen(i,j)*rtd 
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
                                             snoeqv,itype,IVGTYP(I,J)
                          if(sm(i,j) > 0.1)then
                             sfcpct(4)=0.
                          else 
	                     call snfrac_gfs(SNOeqv,IVGTYP(I,J),snofrac)
	                     sfcpct(4)=snofrac
                          end if
                          if(i==ii.and.j==jj)print*,'sno,itype,ivgtyp,sfcpct(4) = ',     &
                                             snoeqv,itype,IVGTYP(I,J),sfcpct(4)
                       else if(MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR')then
                          sfcpct(4)=pctsno(i,j)
                       else          
                          itype=IVGTYP(I,J)
                          IF(itype == 0)itype=8
                          CALL SNFRAC (SNO(I,J),IVGTYP(I,J),snofrac)
	                  sfcpct(4)=snofrac
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
                               .not. land .and. .not. snow
                       if((sfcpct(1)+sfcpct(2)+sfcpct(3)+sfcpct(4)) >1._r_kind)  &
                          print*,'ERROR sfcpct ',i,j,sfcpct(1)  &
                                 ,sfcpct(2),sfcpct(3),sfcpct(4)
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

                       if(gridtype=='B' .or. gridtype=='E')then
                          surface(1)%wind_speed            = sqrt(u10h(i,j)*u10h(i,j)   &
                                                                  +v10h(i,j)*v10h(i,j))
                       else
                          surface(1)%wind_speed            = sqrt(u10(i,j)*u10(i,j)   &
                                                                  +v10(i,j)*v10(i,j))
                       end if

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
                       ! Debug print
                       if(debugprint)then       
                          if(surface(1)%wind_speed<0. .or. surface(1)%wind_speed>200.)  &
                             print*,'bad 10 m wind'
                          if(surface(1)%water_coverage<0. .or. surface(1)%water_coverage>1.) &
                             print*,'bad water coverage'
                          if(surface(1)%land_coverage<0. .or. surface(1)%land_coverage>1.)  &
                             print*,'bad land coverage'
                          if(surface(1)%ice_coverage<0. .or. surface(1)%ice_coverage>1.)  &
                             print*,'bad ice coverage'
                          if(surface(1)%snow_coverage<0. .or. surface(1)%snow_coverage>1.)  &
                             print*,'bad snow coverage'
                          if(surface(1)%land_temperature<0. .or. surface(1)%land_temperature>350.)  &
                             print*,'bad land T'
                          if(surface(1)%soil_moisture_content<0. .or. surface(1)%soil_moisture_content>600.) &
                             print*,'bad soil_moisture_content'
                          if(surface(1)%vegetation_fraction<0. .or. surface(1)%vegetation_fraction>1.) &
                             print*,'bad vegetation cover'
                          if(surface(1)%snow_depth<0. .or.  surface(1)%snow_depth>10000.) &
                             print*,'bad snow_depth'
                          if(MODELNAME == 'GFS' .and. (itype<0 .or. itype>13)) &
                             print*,'bad veg type'
                       end if
       
                       if(i==ii.and.j==jj)print*,'sample surface in CALRAD=',           &
                             i,j,surface(1)%wind_speed,surface(1)%water_coverage,       &
                             surface(1)%land_coverage,surface(1)%ice_coverage,          &
                             surface(1)%snow_coverage,surface(1)%land_temperature,      &
                             surface(1)%snow_temperature,surface(1)%water_temperature,  &
                             surface(1)%ice_temperature,surface(1)%vegetation_fraction, &
                             surface(1)%soil_temperature,surface(1)%snow_depth,         &
                             surface(1)%land_type,sm(i,j)

                       !       Load profiles into model layers

                       !       Load atmosphere profiles into RTM model layers
                       !       CRTM counts from top down just as post does
                       if(i==ii.and.j==jj)print*,'TOA= ',atmosphere(1)%level_pressure(0)
                       do k = 1,lm
                          atmosphere(1)%level_pressure(k) = pint(i,j,k+1)/r100
                          atmosphere(1)%pressure(k)       = pmid(i,j,k)/r100
                          atmosphere(1)%temperature(k)    = t(i,j,k)
	                  atmosphere(1)%absorber(k,1)     = max(0.  &
                                                              ,q(i,j,k)*h1000/(h1-q(i,j,k))) ! use mixing ratio like GSI
                          atmosphere(1)%absorber(k,2)     = max(ozsmall,o3(i,j,k)*constoz)
                          !        atmosphere(1)%absorber(k,2)     = max(ozsmall,o3(i,j,k)*h1000) ! convert to g/kg
                          ! fill in cloud mixing ratio later  
                          if(debugprint)then
                             if(atmosphere(1)%level_pressure(k)<0. .or. atmosphere(1)%level_pressure(k)>1060.) &
                                print*,'bad atmosphere(1)%level_pressure'  &
                                       ,i,j,k,atmosphere(1)%level_pressure(k)     
                             if(atmosphere(1)%pressure(k)<0. .or.   &
                                              atmosphere(1)%pressure(k)>1060.)  &
                                print*,'bad atmosphere(1)%pressure'  &
                                       ,i,j,k,atmosphere(1)%pressure(k) 
                             if(atmosphere(1)%temperature(k)<0. .or.   &
                                              atmosphere(1)%temperature(k)>400.)  &
                                print*,'bad atmosphere(1)%temperature'
                             !        if(atmosphere(1)%absorber(k,1)<0. .or.   &
                             !     &      atmosphere(1)%absorber(k,1)>1.)  &
                             !     &      print*,'bad atmosphere water vapor'
                             !        if(atmosphere(1)%absorber(k,2)<0. .or.   &
                             !     &      atmosphere(1)%absorber(k,1)>1.)  &
                             !     &      print*,'bad atmosphere o3'
                          end if
                          if(i==ii.and.j==jj)print*,'sample atmosphere in CALRAD=',  &
      	                           i,j,k,atmosphere(1)%level_pressure(k),atmosphere(1)%pressure(k),  &
                                   atmosphere(1)%temperature(k),atmosphere(1)%absorber(k,1),  &
                                   atmosphere(1)%absorber(k,2)
                          ! Specify clouds
                          dpovg=(pint(i,j,k+1)-pint(i,j,k))/g !crtm uses column integrated field
                          if(imp_physics==99)then
                             atmosphere(1)%cloud(1)%effective_radius(k) = 10.
	                     atmosphere(1)%cloud(1)%water_content(k) = max(0.,qqw(i,j,k)*dpovg)
                             ! GFS uses temperature and ice concentration dependency formulation to determine 
                             ! effetive radis for cloud ice since GFS does not output ice concentration yet, 
                             !use default 50 um
	                     atmosphere(1)%cloud(2)%effective_radius(k) = 50.	 
                             atmosphere(1)%cloud(2)%water_content(k) = max(0.,qqi(i,j,k)*dpovg)
                             if(debugprint)then
                                if(atmosphere(1)%cloud(1)%water_content(k)<0. .or.   &
                                        atmosphere(1)%cloud(1)%water_content(k)>1.)  &
                                   print*,'bad atmosphere cloud water'
                                if(atmosphere(1)%cloud(2)%water_content(k)<0. .or.   &
                                        atmosphere(1)%cloud(2)%water_content(k)>1.)  &
                                   print*,'bad atmosphere cloud ice'
                             end if
                          else if(imp_physics==5 .or. imp_physics==85 .or. imp_physics==95)then
                             atmosphere(1)%cloud(1)%effective_radius(k) = 10.
                             atmosphere(1)%cloud(1)%water_content(k) = max(0.,qqw(i,j,k)*dpovg)
                             atmosphere(1)%cloud(2)%effective_radius(k) = 25.
                             atmosphere(1)%cloud(2)%water_content(k) = max(0.,qqi(i,j,k)*dpovg)
                             atmosphere(1)%cloud(3)%effective_radius(k) = 200.
                             atmosphere(1)%cloud(3)%water_content(k) = max(0.,qqr(i,j,k)*dpovg)
                             RHO=pmid(i,j,k)/(RD*T(I,J,K)*(1.+D608*Q(I,J,K)))
                             if(F_RimeF(i,j,k)<=5.0)then
                                RHOX=100
                                if(NLICE(I,J,K)>0.)                                      &
                                   atmosphere(1)%cloud(4)%effective_radius(k) =          &
                                      1.0E6*1.5*(RHO*qqs(i,j,k)/(PI*RHOX*NLICE(I,J,K)))**(1./3.) !convert to microns
                                atmosphere(1)%cloud(4)%water_content(k) = max(0.,qqs(i,j,k)*dpovg)
                                atmosphere(1)%cloud(5)%effective_radius(k) = 0.
                                atmosphere(1)%cloud(5)%water_content(k) =0.
                                atmosphere(1)%cloud(6)%effective_radius(k) = 0.
                                atmosphere(1)%cloud(6)%water_content(k) =0.
                             else if(F_RimeF(i,j,k)<=20.0)then
                                atmosphere(1)%cloud(4)%effective_radius(k) = 0.
                                atmosphere(1)%cloud(4)%water_content(k) =0.
                                RHOX=400.
                                if(NLICE(I,J,K)>0.)                                       &
                                   atmosphere(1)%cloud(5)%effective_radius(k) =           &
                                      1.0E6*1.5*(RHO*qqs(i,j,k)/(PI*RHOX*NLICE(I,J,K)))**(1./3.)
                                atmosphere(1)%cloud(5)%water_content(k) =max(0.,qqs(i,j,k)*dpovg)
                                atmosphere(1)%cloud(6)%effective_radius(k) = 0.
                                atmosphere(1)%cloud(6)%water_content(k) =0.
                             else
                                atmosphere(1)%cloud(4)%effective_radius(k) = 0.
                                atmosphere(1)%cloud(4)%water_content(k) =0.
                                atmosphere(1)%cloud(5)%effective_radius(k) = 0.
                                atmosphere(1)%cloud(5)%water_content(k) =0.
                                RHOX=900.
                                if(NLICE(I,J,K)>0.)                                     &
                                   atmosphere(1)%cloud(6)%effective_radius(k) =         &
                                      1.0E6*1.5*(RHO*qqs(i,j,k)/(PI*RHOX*NLICE(I,J,K)))**(1./3.)
                                atmosphere(1)%cloud(6)%water_content(k) =max(0.,qqs(i,j,k)*dpovg)
                             end if
                             if(debugprint .and. i==im/2 .and. j==jsta)                     &
                                print*,'sample precip ice radius= ',i,j,k, F_RimeF(i,j,k),  &
                                       atmosphere(1)%cloud(4)%effective_radius(k),          &
                                       atmosphere(1)%cloud(4)%water_content(k),             &
                                       atmosphere(1)%cloud(5)%effective_radius(k),          &
                                       atmosphere(1)%cloud(5)%water_content(k),             &
                                       atmosphere(1)%cloud(6)%effective_radius(k),          &
                                       atmosphere(1)%cloud(6)%water_content(k)

                          else if(imp_physics==8 .or. imp_physics==6 .or. imp_physics==2)then
                             atmosphere(1)%cloud(1)%effective_radius(k) = 10.
                             atmosphere(1)%cloud(1)%water_content(k) = max(0.,qqw(i,j,k)*dpovg)
                             atmosphere(1)%cloud(2)%effective_radius(k) = 25.
                             atmosphere(1)%cloud(2)%water_content(k) = max(0.,qqi(i,j,k)*dpovg)
                             atmosphere(1)%cloud(3)%effective_radius(k) = 200.
                             atmosphere(1)%cloud(3)%water_content(k) = max(0.,qqr(i,j,k)*dpovg)
                             atmosphere(1)%cloud(4)%effective_radius(k) = 250.
                             atmosphere(1)%cloud(4)%water_content(k) = max(0.,qqs(i,j,k)*dpovg)
                             atmosphere(1)%cloud(5)%effective_radius(k) = 350.
                             atmosphere(1)%cloud(5)%water_content(k) = max(0.,qqg(i,j,k)*dpovg)
                          end if 
                       end do

                       !bsf - start
                       !-- Add subgrid-scale convective clouds for WRF runs
                       if(icu_physics==2) then
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
                                             atmosphere(1)%cloud(2)%water_content(k) + dpovg*q_conv
                                else
	                           atmosphere(1)%cloud(1)%water_content(k) =   &
                                              atmosphere(1)%cloud(1)%water_content(k) + dpovg*q_conv
                                endif
                             end do   !-- do k = lctop,lcbot
                          endif      !-- if (lcbot-lctop > 1) then
                       endif        !-- if (MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR') then
                       !bsf - end

                       !     call crtm forward model
                       error_status = crtm_forward(atmosphere,surface,                    &
                                      geometryinfo,channelinfo(sensorindex:sensorindex),  &
                                      rtsolution)
                       if (error_status /=0) then
                          print*,'***ERROR*** during crtm_forward call ',  &
                                 error_status
                          do n=1,channelinfo(sensorindex)%n_channels
                             tb(i,j,n)=spval
                          end do
                       else
                          do n=1,channelinfo(sensorindex)%n_channels
                             tb(i,j,n)=rtsolution(n,1)%brightness_temperature
                          end do
                          if(i==ii.and.j==jj) then
                             do n=1,channelinfo(sensorindex)%n_channels
 3303                           format('Sample rtsolution(',I0,',',I0,') in CALRAD = ',F0.3)
                                print 3303,n,1,rtsolution(n,1)%brightness_temperature
                             enddo
                             do n=1,channelinfo(sensorindex)%n_channels
 3304                           format('Sample tb(',I0,',',I0,',',I0,') in CALRAD = ',F0.3)
                                print 3304,i,j,n,tb(i,j,n)
                             enddo
                          endif
                          !        if(tb1(i,j) < 400. )  &
                          !     &        print*,'good tb1 ',i,j,tb1(i,j),gdlat(i,j),gdlon(i,j)
                          !        if(tb2(i,j) > 400.)print*,'bad tb2 ',i,j,tb2(i,j)
                          !        if(tb3(i,j) > 400.)print*,'bad tb3 ',i,j,tb3(i,j)
                          !        if(tb4(i,j) > 400.)print*,'bad tb4 ',i,j,tb4(i,j)
                       end if
                    else
                       do n=1,channelinfo(sensorindex)%n_channels
                          tb(i,j,n)=spval
                       end do
                    END IF ! endif block for allowable satellite zenith angle 
                 end do ! end loop for i
              end do ! end loop for j 
  
               !      error_status = crtm_destroy(channelinfo)
               !      if (error_status /= success) &
               !     &   print*,'ERROR*** crtm_destroy error_status=',error_status

              if(isis=='ssmis_f17') then ! writing f17 ssmis to grib
                 do ixchan=1,4
                    !ichan=14+ixchan  ! channel number
                    ichan=ixchan ! using select_channels, we discard channels 1-14 and 19-24
                    igot=iget(803+ixchan) ! iget(804) ... iget(807)
                    if(igot > 0) then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,ichan)
                          enddo
                       enddo
                       id(1:25) = 0
                       id(02) = 2
                       id(09) = 112
                       id(10) = 117
                       !id(11) = ichan
                       id(11) = ichan+14
                       if(grib=="grib1") then
                          call gribit(igot,lvls(1,igot), grid1,im,jm)
                       else if(grib=="grib2" )then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              endif
  
              if (isis=='imgr_g12')then  ! writing goes 12 to grib
                 do ixchan=1,4
                    ichan=ixchan
                    igot=iget(455+ixchan)
                    if(igot>0) then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,ichan)
                          enddo
                       enddo
                       id(1:25) = 0
                       id(02) = 129
                       if(grib=="grib1") then
                          call gribit(igot,lvls(1,igot), grid1,im,jm)
                       else if(grib=="grib2" )then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              end if  ! end of outputting goes 12 

              if (isis=='imgr_g11')then  ! writing goes 11 to grib
                 do ixchan=1,4
                    ichan=ixchan 
                    igot=iget(459+ixchan)
                    if(igot>0) then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,ichan)
                          enddo
                       enddo
                       id(1:25) = 0
                       id(02) = 130
                       if(grib=="grib1") then
                          call gribit(igot,lvls(1,igot), grid1,im,jm)
                       else if(grib=="grib2" )then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              end if  ! end of outputting goes 11
           end if nonnadir  ! end if for computing simulated radiance with zenith angle correction
      
           ! Deallocate arrays
           CALL crtm_atmosphere_destroy(atmosphere(1))
           CALL crtm_surface_destroy(surface(1))
           CALL crtm_rtsolution_destroy(rtsolution)
           if (crtm_atmosphere_associated(atmosphere(1))) &
              write(6,*)' ***ERROR** destroying atmosphere.'
           if (crtm_surface_associated(surface(1))) &
              write(6,*)' ***ERROR** destroying surface.'
           if (ANY(crtm_rtsolution_associated(rtsolution))) &
              write(6,*)' ***ERROR** destroying rtsolution.'
           deallocate(tb)
           deallocate (rtsolution)
       
!     
        end if  sensor_avail  ! end of if block for only calling crtm when sepcific sensor is requested in the control file
     end do  sensordo ! end looping for different satellite

     error_status = crtm_destroy(channelinfo)
     if (error_status /= success) &
         print*,'ERROR*** crtm_destroy error_status=',error_status
     deallocate(channelinfo) 
  endif ifactive ! for all iget logical
  return
end SUBROUTINE CALRAD_WCLOUD
