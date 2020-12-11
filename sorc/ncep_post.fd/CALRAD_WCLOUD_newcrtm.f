!> @file
!
!> THIS ROUTINE COMPUTES MODEL DERIVED BRIGHTNESS TEMPERATURE
!! USING CRTM. IT IS PATTERNED AFTER GSI SETUPRAD WITH TREADON'S HELP     
!!     
!! PROGRAM HISTORY LOG:
!! -  11-02-06 Jun WANG   - addgrib2 option 
!! -  14-12-09 WM LEWIS ADDED:
!!            FUNCTION EFFR TO COMPUTE EFFECTIVE PARTICLE RADII 
!!            CHANNEL SELECTION USING LVLS FROM WRF_CNTRL.PARM
!! -  19-04-01 Sharon NEBUDA - Added output option for GOES-16 & GOES-17 ABI IR Channels 7-16
!! -  20-04-09 Tracy Hertneky - Added Himawari-8 AHI CH7-CH16
!!
!!   OUTPUT FILES:
!!     NONE
!!     
!!   SUBPROGRAMS CALLED:
!!     UTILITIES:
!!
!!     LIBRARY:
!!     /nwprod/lib/sorc/crtm2
!!
!! @author CHUANG @date 07-01-17       
!!     
      SUBROUTINE CALRAD_WCLOUD

  use vrbls3d, only: o3, pint, pmid, t, q, qqw, qqi, qqr, f_rimef, nlice, nrain, qqs, qqg, &
                     qqnr, qqni, qqnw
  use vrbls2d, only: czen, ivgtyp, sno, pctsno, ths, vegfrc, si, u10h, v10h, u10,&
       v10, smstot, hbot, htop, cnvcfr
  use masks, only: gdlat, gdlon, sm, lmh, sice
  use soil, only:
  use gridspec_mod, only: gridtype
  use cmassi_mod, only: TRAD_ice
  use kinds, only: r_kind,r_single,r_double,i_kind
  use crtm_module, only: crtm_atmosphere_type,crtm_surface_type,crtm_geometry_type, &
       crtm_surface_create,o3_id,co2_id,crtm_forward,mass_mixing_ratio_units, &
       crtm_atmosphere_create, &
       crtm_options_type,crtm_destroy,crtm_init,SPECIFIC_AMOUNT_UNITS, &
       success,crtm_options_destroy,crtm_options_create, crtm_options_associated
       
  use crtm_rtsolution_define, only: crtm_rtsolution_type, crtm_rtsolution_create, &
       crtm_rtsolution_destroy, crtm_rtsolution_associated 
  use crtm_spccoeff, only: sc
  use crtm_atmosphere_define, only:h2o_id,crtm_atmosphere_associated, &
       crtm_atmosphere_destroy,volume_mixing_ratio_units,crtm_atmosphere_zero
  use crtm_surface_define, only: crtm_surface_destroy, crtm_surface_associated, &
       crtm_surface_zero
  use crtm_channelinfo_define, only: crtm_channelinfo_type
!  use crtm_channelinfo_define, only: crtm_channelinfo_type, AntCorr_type
  use crtm_parameters, only: limit_exp,toa_pressure,max_n_layers,MAX_SENSOR_SCAN_ANGLE
  use crtm_cloud_define, only:  water_cloud,ice_cloud,rain_cloud,snow_cloud,graupel_cloud,hail_cloud
  use message_handler, only: success,warning, display_message

  use params_mod, only: pi, rtd, p1000, capa, h1000, h1, g, rd, d608, qconv
  use rqstfld_mod, only: iget, id, lvls, iavblfld
  use ctlblk_mod, only: modelname, ivegsrc, novegtype, imp_physics, lm, spval, icu_physics,&
              grib, cfld, fld_info, datapd, idat, im, jsta, jend, jm, me
!     
  implicit none

  !     DECLARE VARIABLES.
  !     
  ! Mapping land surface type of GFS to CRTM
  !  Note: index 0 is water, and index 13 is ice. The two indices are not
  !        used and just assigned to COMPACTED_SOIL.
  !integer, parameter, dimension(0:13) :: gfs_to_crtm=(/COMPACTED_SOIL,     &
  !       BROADLEAF_FOREST, BROADLEAF_FOREST, BROADLEAF_PINE_FOREST, PINE_FOREST, &
  !       PINE_FOREST, BROADLEAF_BRUSH, SCRUB, SCRUB, SCRUB_SOIL, TUNDRA,         &
  !       COMPACTED_SOIL, TILLED_SOIL, COMPACTED_SOIL/)

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

  ! For land, the land types
  !INTEGER, PARAMETER :: N_VALID_LAND_TYPES = 20
  INTEGER, PARAMETER :: INVALID_LAND             =  0
  INTEGER, PARAMETER :: COMPACTED_SOIL           =  1
  INTEGER, PARAMETER :: TILLED_SOIL              =  2
  INTEGER, PARAMETER :: SAND                     =  3
  INTEGER, PARAMETER :: ROCK                     =  4
  INTEGER, PARAMETER :: IRRIGATED_LOW_VEGETATION =  5
  INTEGER, PARAMETER :: MEADOW_GRASS             =  6
  INTEGER, PARAMETER :: SCRUB                    =  7
  INTEGER, PARAMETER :: BROADLEAF_FOREST         =  8
  INTEGER, PARAMETER :: PINE_FOREST              =  9
  INTEGER, PARAMETER :: TUNDRA                   = 10
  INTEGER, PARAMETER :: GRASS_SOIL               = 11
  INTEGER, PARAMETER :: BROADLEAF_PINE_FOREST    = 12
  INTEGER, PARAMETER :: GRASS_SCRUB              = 13
  INTEGER, PARAMETER :: SOIL_GRASS_SCRUB         = 14
  INTEGER, PARAMETER :: URBAN_CONCRETE           = 15
  INTEGER, PARAMETER :: PINE_BRUSH               = 16
  INTEGER, PARAMETER :: BROADLEAF_BRUSH          = 17
  INTEGER, PARAMETER :: WET_SOIL                 = 18
  INTEGER, PARAMETER :: SCRUB_SOIL               = 19
  INTEGER, PARAMETER :: BROADLEAF70_PINE30       = 20
 

  integer, allocatable:: model_to_crtm(:)
  integer, parameter:: ndat=100
  ! CRTM structure variable declarations.
  integer,parameter::  n_absorbers = 2
  !      integer,parameter::  n_clouds = 4 
  integer,parameter::  n_aerosols = 0
  ! Add your sensors here
  integer(i_kind),parameter:: n_sensors=22
  character(len=20),parameter,dimension(1:n_sensors):: sensorlist= &
      (/'imgr_g15            ', &
        'imgr_g13            ', &
        'imgr_g12            ', &
        'imgr_g11            ', &
        'amsre_aqua          ', &
        'tmi_trmm            ', &
        'ssmi_f13            ', &
        'ssmi_f14            ', &
        'ssmi_f15            ', &
        'ssmis_f16           ', &
        'ssmis_f17           ', &
        'ssmis_f18           ', &
        'ssmis_f19           ', &
        'ssmis_f20           ', &
        'seviri_m10          ', &
        'imgr_mt2            ', &
        'imgr_mt1r           ', &
        'imgr_insat3d        ', &
        'abi_gr              ', &
        'abi_g16             ', &
        'abi_g17             ', &
        'ahi_himawari8       '/)
  character(len=13),parameter,dimension(1:n_sensors):: obslist=  &
      (/'goes_img     ', &
        'goes_img     ', &
        'goes_img     ', &
        'goes_img     ', &
        'amsre        ', &
        'tmi          ', &
        'ssmi         ', &
        'ssmi         ', &
        'ssmi         ', &
        'ssmis        ', &
        'ssmis        ', &
        'ssmis        ', &
        'ssmis        ', &
        'ssmis        ', &
        'seviri       ', &
        'imgr_mt2     ', &
        'imgr_mt1r    ', &
        'imgr_insat3d ', &
        'abi          ', &
        'abi          ', &
        'abi          ', &
        'ahi_himawari8'/)
  character(len=20),dimension(1:n_sensors):: sensorlist_local
!
  integer(i_kind) sensorindex
  integer(i_kind) lunin,nobs,nchanl,nreal
  integer(i_kind) error_status,itype
  integer(i_kind) err1,err2,err3,err4
  integer(i_kind) i,j,k,msig
  integer(i_kind) lcbot,lctop   !bsf
  integer jdn,ichan,ixchan,igot
  integer isat
 
! Wm Lewis: added 
  real :: EFFR 
 
  real(r_kind),parameter:: r100=100.0_r_kind
  real,parameter:: ozsmall = 1.e-10 ! to convert to mass mixing ratio
  real(r_kind) tsfc 
  real(r_double),dimension(4):: sfcpct
  real(r_kind) snodepth,vegcover
  real snoeqv
  real snofrac
  real(r_kind),dimension(im,jsta:jend):: tb1,tb2,tb3,tb4
  real(r_kind),allocatable :: tb(:,:,:)
  real,dimension(im,jm):: grid1
  real sun_zenith,sun_azimuth, dpovg, sun_zenith_rad
  real sat_zenith
  real q_conv   !bsf
  real,parameter:: constoz = 604229.0_r_kind 
  real sublat,sublon
  real RHO,RHOX
  character(13)::obstype
  character(20)::isis
  character(20)::isis_local

  logical hirs2,msu,goessndr,hirs3,hirs4,hirs,amsua,amsub,airs,hsb  &
            ,goes_img,abi,seviri, mhs,insat3d
  logical avhrr,avhrr_navy,lextra,ssu
  logical ssmi,ssmis,amsre,amsre_low,amsre_mid,amsre_hig,change
  logical ssmis_las,ssmis_uas,ssmis_env,ssmis_img
  logical sea,mixed,land,ice,snow,toss
  logical micrim,microwave
  logical post_abig16, post_abig17, post_abigr ! if true, user requested at least one abi channel
  logical fix_abig16, fix_abig17   ! if true, abi_g16, abi_g17 fix files are available
  logical post_ahi8 ! if true, user requested at least on ahi channel (himawari8)
  !  logical,dimension(nobs):: luse
  logical, parameter :: debugprint = .false.
  type(crtm_atmosphere_type),dimension(1):: atmosphere
  type(crtm_surface_type),dimension(1) :: surface
  type(crtm_geometry_type),dimension(1) :: geometryinfo
  type(crtm_options_type),dimension(1)      :: options

  type(crtm_rtsolution_type),allocatable,dimension(:,:):: rtsolution
  type(crtm_channelinfo_type),allocatable,dimension(:) :: channelinfo
!     
  integer ii,jj,n_clouds,n,nc
  integer,external :: iw3jdn
  !

  !*****************************************************************************
  ! This code and sensorlist_local, isis_local can be modified/removed when the 
  ! linked CRTM version is updated with fix files abi_g16 & abi_g17
   fix_abig16 = .False.
   fix_abig17 = .False.
   do n=1, n_sensors
     sensorlist_local(n) = sensorlist(n)
     if (sensorlist(n) == 'abi_g16') then  ! check if fix file is available
       inquire(file='abi_g16.SpcCoeff.bin',exist=fix_abig16)
       if (.not.fix_abig16) sensorlist_local(n) = 'abi_gr              '
     endif
     if (sensorlist(n) == 'abi_g17') then
       inquire(file='abi_g17.SpcCoeff.bin',exist=fix_abig17)
       if (.not.fix_abig17) sensorlist_local(n) = 'abi_gr              '
     endif
   enddo
  

  ! Mapping land surface type of NMM to CRTM
  !if(MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR')then 
   if(ivegsrc==1)then  !IGBP veg type
      allocate(model_to_crtm(novegtype) )
      model_to_crtm=(/PINE_FOREST, BROADLEAF_FOREST, PINE_FOREST,       &
           BROADLEAF_FOREST,BROADLEAF_PINE_FOREST, SCRUB, SCRUB_SOIL, &
           BROADLEAF_BRUSH,BROADLEAF_BRUSH, SCRUB, BROADLEAF_BRUSH,   &
           TILLED_SOIL, URBAN_CONCRETE,TILLED_SOIL, INVALID_LAND,     &
           COMPACTED_SOIL, INVALID_LAND, TUNDRA,TUNDRA, TUNDRA/)
   else if(ivegsrc==0)then ! USGS veg type
      allocate(model_to_crtm(novegtype) )
      model_to_crtm=(/URBAN_CONCRETE,       &
           COMPACTED_SOIL, IRRIGATED_LOW_VEGETATION, GRASS_SOIL, MEADOW_GRASS,   &
           MEADOW_GRASS, MEADOW_GRASS, SCRUB, GRASS_SCRUB, MEADOW_GRASS,         &
           BROADLEAF_FOREST, PINE_FOREST, BROADLEAF_FOREST, PINE_FOREST,         &
           BROADLEAF_PINE_FOREST, COMPACTED_SOIL, WET_SOIL, WET_SOIL,            &
           IRRIGATED_LOW_VEGETATION, TUNDRA, TUNDRA, TUNDRA, TUNDRA,             &
           COMPACTED_SOIL/)
   else if(ivegsrc==2)then ! old GFS veg type
      allocate(model_to_crtm(0:novegtype) )
      model_to_crtm=(/COMPACTED_SOIL,     &
         BROADLEAF_FOREST, BROADLEAF_FOREST, BROADLEAF_PINE_FOREST, &
         PINE_FOREST, PINE_FOREST, BROADLEAF_BRUSH, SCRUB, SCRUB, SCRUB_SOIL, &
         TUNDRA, COMPACTED_SOIL, TILLED_SOIL, COMPACTED_SOIL/)
   else
      print*,'novegtype=',novegtype
      print*,'model veg type not supported by post in calling crtm ' 
      print*,'skipping generation of simulated radiance' 
      return
   end if 
  !end if 

  !10 channels, easier to set a logical
  post_abig16=.false.
  do n = 927, 927+9  ! 927 set in RQSTFLD.f
    if (iget(n) > 0) post_abig16=.true.
  enddo
  post_abig17=.false.
  do n = 937, 937+9  ! 937 set in RQSTFLD.f
    if (iget(n) > 0) post_abig17=.true.
  enddo
  post_abigr=.false.
  do n = 958, 958+9  ! 958 set in RQSTFLD.f
    if (iget(n) > 0) post_abigr=.true.
  enddo
  post_ahi8=.false.
  do n = 969, 969+9  ! 969 set in RQSTFLD.f
    if (iget(n) > 0) post_ahi8=.true.
  enddo


  !     DO NOT FORGET TO ADD YOUR NEW IGET HERE (IF YOU'VE ADDED ONE)      
  !     START SUBROUTINE CALRAD.
  ifactive: if (iget(327) > 0 .or. iget(328) > 0 .or. iget(329) > 0       &
       .or. iget(330) > 0 .or. iget(446) > 0 .or. iget(447) > 0  & 
       .or. iget(448) > 0 .or. iget(449) > 0 .or. iget(456) > 0  &
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
       .or. iget(807) > 0 .or. iget(809) > 0                     &
       .or. iget(810) > 0 .or. iget(811) > 0 .or. iget(812) > 0  &
       .or. iget(813) > 0 .or. iget(814) > 0 .or. iget(815) > 0  &
       .or. iget(816) > 0 .or. iget(817) > 0 .or. iget(818) > 0  &
       .or. iget(819) > 0 .or. iget(820) > 0 .or. iget(821) > 0  &
       .or. iget(822) > 0 .or. iget(823) > 0 .or. iget(824) > 0  &
       .or. iget(825) > 0 .or. iget(826) > 0 .or. iget(827) > 0  &
       .or. iget(828) > 0 .or. iget(829) > 0 .or. iget(830) > 0  &
       .or. iget(831) > 0 .or. iget(832) > 0 .or. iget(833) > 0  &
       .or. iget(834) > 0 .or. iget(835) > 0 .or. iget(836) > 0  &
       .or. iget(837) > 0 .or. iget(838) > 0 .or. iget(839) > 0  &
       .or. iget(840) > 0 .or. iget(841) > 0 .or. iget(842) > 0  &
       .or. iget(843) > 0 .or. iget(844) > 0 .or. iget(845) > 0  &
       .or. iget(846) > 0 .or. iget(847) > 0 .or. iget(848) > 0  &
       .or. iget(849) > 0 .or. iget(850) > 0 .or. iget(851) > 0  &
       .or. iget(852) > 0 .or. iget(856) > 0 .or. iget(857) > 0  &
       .or. iget(860) > 0 .or. iget(861) > 0  &
       .or. iget(862) > 0 .or. iget(863) > 0 .or. iget(864) > 0  &
       .or. iget(865) > 0 .or. iget(866) > 0 .or. iget(867) > 0  &
       .or. iget(868) > 0 .or. iget(869) > 0 .or. iget(870) > 0  &
       .or. iget(871) > 0 .or. iget(872) > 0 .or. iget(873) > 0  &
       .or. iget(874) > 0 .or. iget(875) > 0 .or. iget(876) > 0  &
       .or. iget(877) > 0 .or. iget(878) > 0 .or. iget(879) > 0  &
       .or. iget(880) > 0 .or. iget(881) > 0 .or. iget(882) > 0  &
       .or. post_ahi8 & 
       .or. post_abig16 .or. post_abig17 .or. post_abigr ) then

     ! specify numbers of cloud species    
     ! Thompson==8, Ferrier==5,95, WSM6==6, Lin==2
     if(imp_physics==99 .or. imp_physics==98)then ! Zhao Scheme
        n_clouds=2 ! GFS uses Zhao scheme
     else if(imp_physics==5 .or. imp_physics==85 .or. imp_physics==95)then
        n_clouds=6  ! change to 6 cloud types because microwave is sensitive to density
     else if(imp_physics==8 .or. imp_physics==6 .or. imp_physics==2 &
        .or. imp_physics==28 .or. imp_physics==11)then
        n_clouds=5
     else
        n_clouds=0
        print*,'Warning: number of cloud species (n_clouds) being set to zero for imp_physics=',imp_physics
     end if

     ! Initialize debug print gridpoint index to middle of tile:
     ii=im/2
     jj=(jsta+jend)/2

     ! Initialize ozone to zeros for WRF NMM and ARW for now
     if (MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR' &
      )o3=0.0
     ! Compute solar zenith angle for GFS, ARW now computes czen in INITPOST
!     if (MODELNAME == 'GFS')then
        jdn=iw3jdn(idat(3),idat(1),idat(2))
	do j=jsta,jend
	   do i=1,im
	      call zensun(jdn,float(idat(4)),gdlat(i,j),gdlon(i,j)       &
      	                  ,pi,sun_zenith,sun_azimuth)
              sun_zenith_rad=sun_zenith/rtd              
              czen(i,j)=cos(sun_zenith_rad)
	   end do
	end do
        if(jj>=jsta .and. jj<=jend.and.debugprint)                   &
            print*,'sample GFS zenith angle=',acos(czen(ii,jj))*rtd   
!     end if	       
     ! Initialize CRTM.  Load satellite sensor array.
     ! The optional arguments Process_ID and Output_Process_ID limit
     ! generation of runtime informative output to mpi task
     ! Output_Process_ID (which here is set to be task 0)
     if(me==0)print*,'success in CALRAD= ',success
     allocate( channelinfo(n_sensors))

     error_status = crtm_init(sensorlist_local,channelinfo,   &
          Process_ID=0,Output_Process_ID=0 )
     if(me==0)print*, 'channelinfo after init= ',channelinfo(1)%sensor_id, &
              channelinfo(2)%sensor_id
     if (error_status /= 0_i_kind)                                  &
         write(6,*)'ERROR*** crtm_init error_status=',error_status

     ! Restrict channel list to those which are selected for simulation 
     ! in the LVLS filed of wrf_cntrl.parm (does not currently apply 
     ! to all sensors / channels).

     ! GOES-13
     if(iget(868)>0)then
     call select_channels_L(channelinfo(2),4,(/ 1,2,3,4 /),lvls(1:4,iget(868)),iget(868))
     endif
     ! GOES-15
     if(iget(872)>0)then
     call select_channels_L(channelinfo(1),4,(/ 1,2,3,4 /),lvls(1:4,iget(872)),iget(872))
     endif
     ! GOES-16 
     if(post_abig16)then
       nchanl=0
       do n = 927, 927+9  ! 927 set in RQSTFLD.f
         if (iget(n) > 0) then
           nchanl = nchanl+1
         endif
       enddo
       if (nchanl > 0 .and. nchanl <10) then 
         do n = 927, 927+9  ! 927 set in RQSTFLD.f
           if (iget(n) == 0) channelinfo(19)%Process_Channel(n-927+1)=.False.  !  turn off channel processing
         enddo
       endif
     endif
     ! GOES-17 
     if(post_abig17)then
       nchanl=0
       do n = 937, 937+9  ! 937 set in RQSTFLD.f
         if (iget(n) > 0) then
           nchanl = nchanl+1
         endif
       enddo
       if (nchanl > 0 .and. nchanl <10) then 
         do n = 937, 937+9  ! 927 set in RQSTFLD.f
           if (iget(n) == 0) channelinfo(20)%Process_Channel(n-937+1)=.False.  !  turn off channel processing
         enddo
       endif
     endif
     ! GOES-R for NADIR output 
     if(post_abigr)then
       nchanl=0
       do n = 958, 958+9  ! 958 set in RQSTFLD.f
         if (iget(n) > 0) then
           nchanl = nchanl+1
         endif
       enddo
       if (nchanl > 0 .and. nchanl <10) then 
         do n = 958, 958+9  ! 958 set in RQSTFLD.f
           if (iget(n) == 0) channelinfo(20)%Process_Channel(n-958+1)=.False.  !  turn off channel processing
         enddo
       endif
     endif

     ! Himawari-8 AHI infrared 
     if(post_ahi8)then
       nchanl=0
       do n = 969, 969+9  ! 969 set in RQSTFLD.f
         if (iget(n) > 0) then
           nchanl = nchanl+1
         endif
       enddo
       if (nchanl > 0 .and. nchanl <10) then 
         do n = 969, 969+9  ! 969 set in RQSTFLD.f
           if (iget(n) == 0) channelinfo(22)%Process_Channel(n-969+1)=.False.  !  turn off channel processing
         enddo
       endif
     endif

     ! SSMI, F13-F15 (19H,19V,??H,37H,37V,85H,85V)
     if(iget(800)>0)then
     call select_channels_L(channelinfo(7),7,(/ 1,2,3,4,5,6,7 /),lvls(1:7,iget(800)),iget(800))
     endif
     if(iget(806)>0)then
     call select_channels_L(channelinfo(8),7,(/ 1,2,3,4,5,6,7 /),lvls(1:7,iget(806)),iget(806))
     endif
     if(iget(812)>0)then
     call select_channels_L(channelinfo(9),7,(/ 1,2,3,4,5,6,7 /),lvls(1:7,iget(812)),iget(812))
     endif
     ! SSMIS, F16-F20 (183H,19H,19V,37H,37V,91H,91V)
     if(iget(818)>0)then
     call select_channels_L(channelinfo(10),24,(/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 /),lvls(1:24,iget(818)),iget(818))
     endif
     if(iget(825)>0)then
     call select_channels_L(channelinfo(11),24,(/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 /),lvls(1:24,iget(825)),iget(825))
     endif
     if(iget(832)>0)then
     call select_channels_L(channelinfo(12),24,(/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 /),lvls(1:24,iget(832)),iget(832))
     endif
     if(iget(839)>0)then
     call select_channels_L(channelinfo(13),24,(/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 /),lvls(1:24,iget(839)),iget(839))
     endif
     if(iget(846)>0)then
     call select_channels_L(channelinfo(14),24,(/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 /),lvls(1:24,iget(846)),iget(846))
     endif
     ! SEVIRI
     if(iget(876)>0)then
     call select_channels_L(channelinfo(15),8,(/ 1,2,3,4,5,6,7,8 /),lvls(1:8,iget(876)),iget(876))
     endif
     ! MT2
     if(iget(860)>0)then
     call select_channels_L(channelinfo(16),4,(/ 1,2,3,4 /),lvls(1:4,iget(860)),iget(860))
     endif
     ! MT1R
     if(iget(864)>0)then
     call select_channels_L(channelinfo(17),4,(/ 1,2,3,4 /),lvls(1:4,iget(864)),iget(864))
     endif
     ! INSAT 3D (Kalpana)
     if(iget(865)>0)then
     call select_channels_L(channelinfo(18),4,(/ 1,2,3,4 /),lvls(1:4,iget(865)),iget(865))
     endif
     ! Himiwari-8 AHI infrared
     if(iget(969)>0)then
     call select_channels_L(channelinfo(19),10,(/1,2,3,4,5,6,7,8,9,10/),lvls(1:10,iget(969)),iget(969))
     endif


     ! Loop over data types to process    
     sensordo: do isat=1,n_sensors

        if(me==0)print*,'n_sensor,obstype,isis',isat,obslist(isat),sensorlist(isat)

        obstype=obslist(isat) 
        isis=trim(sensorlist(isat))

        sensor_avail: if( &
             (isis=='imgr_g12' .and. (iget(327) > 0 .or. iget(328) > 0 &
             .or. iget(329) > 0 .or. iget(330) > 0 .or. iget(456) > 0   &
             .or. iget(457) > 0 .or. iget(458) > 0 .or. iget(459) > 0 )) .OR. &
             (isis=='imgr_g11' .and. (iget(446) > 0 .or. iget(447) > 0 &
             .or. iget(448) > 0 .or. iget(449) > 0 .or. iget(460) > 0   &
             .or. iget(461) > 0 .or. iget(462) > 0 .or. iget(463) > 0)) .OR. &
             (isis=='amsre_aqua' .and. (iget(483) > 0 .or. iget(484) > 0  &
             .or. iget(485) > 0 .or. iget(486) > 0)) .OR. &
             (isis=='tmi_trmm' .and. (iget(488) > 0 .or. iget(489) > 0  &
             .or. iget(490) > 0 .or. iget(491) > 0)) .OR. &
             (isis=='ssmi_f13' .and. iget(800) > 0 ) .OR. &
             (isis=='ssmi_f14' .and. iget(806) > 0 ) .OR. &
             (isis=='ssmi_f15' .and. iget(812) > 0 ) .OR. &
             (isis=='ssmis_f16' .and. iget(818) > 0) .OR. &
             (isis=='ssmis_f17' .and. iget(825) > 0) .OR. &
             (isis=='ssmis_f18' .and. iget(832) > 0) .OR. &
             (isis=='ssmis_f19' .and. iget(839) > 0) .OR. &
             (isis=='ssmis_f20' .and. iget(846) > 0) .OR. &
             (isis=='imgr_mt2' .and. iget(860)>0) .OR. &
             (isis=='imgr_mt1r' .and. iget(864)>0) .OR. &
             (isis=='imgr_insat3d' .and. iget(865)>0) .OR. &
             (isis=='imgr_g13' .and. iget(868)>0) .OR. &
             (isis=='imgr_g15' .and. iget(872)>0) .OR. &
             (isis=='abi_g16'  .and. post_abig16) .OR. &
             (isis=='abi_g17'  .and. post_abig17) .OR. &
             (isis=='abi_gr'   .and. post_abigr) .OR. &
             (isis=='seviri_m10' .and. iget(876)>0) .OR. &
             (isis=='ahi_himawari8' .and. post_ahi8) )then
           if(me==0)print*,'obstype, isis= ',obstype,isis
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
           abi        = obstype == 'abi'
           seviri     = obstype == 'seviri'
           insat3d    = obstype == 'imgr_insat3d'
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
           ! check sensor list
           sensorindex = 0
           sensor_search: do j = 1, n_sensors
              isis_local = isis  ! allows abi_g16 & abi_g17 output using abi_gr fix files
              if (isis=='abi_g16' .and. .not.fix_abig16) then
                isis_local='abi_gr              '
              endif
              if (isis=='abi_g17' .and. .not.fix_abig17) then
                isis_local='abi_gr              '
              endif
              if (channelinfo(j)%sensor_id == isis_local ) then
                 sensorindex = j
                 exit sensor_search
              endif
           end do sensor_search
           if (sensorindex == 0 ) then
              write(6,*)'SETUPRAD:  ***WARNING*** problem with sensorindex=',isis,&
                   ' --> CAN NOT PROCESS isis=',isis,'   TERMINATE PROGRAM EXECUTION'
              stop 19
           endif

!          set Satellite IDs for F19 and F20 to valid values since CRTM will not
!          simulate an instrument w/o a WMO ID:
           if(isis=='ssmis_f19')channelinfo(sensorindex)%WMO_Satellite_Id=287
           if(isis=='ssmis_f20')channelinfo(sensorindex)%WMO_Satellite_Id=289
!          quiet verbose output warning messages
           if(isis=='abi_g16')channelinfo(sensorindex)%WMO_Satellite_Id=270
           if(isis=='abi_g16')channelinfo(sensorindex)%WMO_Sensor_Id=617
           if(isis=='abi_g17')channelinfo(sensorindex)%WMO_Satellite_Id=271
           if(isis=='abi_g17')channelinfo(sensorindex)%WMO_Sensor_Id=617
           if(isis=='abi_gr')channelinfo(sensorindex)%WMO_Satellite_Id=270
           if(isis=='abi_gr')channelinfo(sensorindex)%WMO_Sensor_Id=617

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
           if(imp_physics==99 .or. imp_physics==98)then
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
           else if(imp_physics==8 .or. imp_physics==6 .or. imp_physics==2 .or. imp_physics==28 &
            .or. imp_physics==11)then
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

           ! run crtm for nadir instruments / channels
           nadir: if ( (isis=='imgr_g12' .and. (iget(327)>0 .or. &
                       iget(328)>0 .or. iget(329)>0 .or. iget(330)>0)) .or. &
                       (isis=='imgr_g11' .and. (iget(446)>0 .or. &
                       iget(447)>0 .or. iget(448)>0 .or. iget(449)>0)) .or. &
                       (isis=='amsre_aqua' .and. (iget(483) > 0 .or. iget(484) > 0  &
                       .or. iget(485) > 0 .or. iget(486) > 0)) .OR. &
                       (isis=='tmi_trmm' .and. (iget(488) > 0 .or. iget(489) > 0  &
                       .or. iget(490) > 0 .or. iget(491) > 0)) .OR. &
                        (isis=='abi_gr'  .and. post_abigr) )then

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
                       if(i==ii.and.j==jj.and.debugprint)print*,'sample geometry ',                   &
                               geometryinfo(1)%sensor_zenith_angle                     &
                              ,geometryinfo(1)%source_zenith_angle                     &
                              ,czen(i,j)*rtd 
                       !  Set land/sea, snow, ice percentages and flags
                       if(MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR')then
                          sfcpct(4)=pctsno(i,j)
                       else if(ivegsrc==1)then
                          itype=IVGTYP(I,J)
                          IF(itype == 0)itype=8
                          if(sno(i,j)<spval)then
                             snoeqv=sno(i,j)
                          else
                             snoeqv=0.
                          end if
                          CALL SNFRAC (SNO(I,J),IVGTYP(I,J),snofrac)
                          sfcpct(4)=snofrac
                       else if(ivegsrc==2)then
                          itype=IVGTYP(I,J)
                          itype = min(max(0,ivgtyp(i,j)),13)
                          if(sno(i,j)<spval)then
                             snoeqv=sno(i,j)
                          else
                             snoeqv=0.
                          end if
                          if(i==ii.and.j==jj.and.debugprint)print*,'sno,itype,ivgtyp B cing snfrc = ',  &
                                                     snoeqv,itype,IVGTYP(I,J)
                          if(sm(i,j) > 0.1)then
                             sfcpct(4)=0.
                          else
                             call snfrac_gfs(SNOeqv,IVGTYP(I,J),snofrac)
                             sfcpct(4)=snofrac
                          end if
                       end if

!                       if (MODELNAME == 'GFS')then ! GFS uses 13 veg types
!                          itype=IVGTYP(I,J)
!                          itype = min(max(0,ivgtyp(i,j)),13)
!                          !         IF(itype <= 0 .or. itype > 13)itype=7 !use scrub for ocean point
!                          if(sno(i,j)/=spval)then
!                             snoeqv=sno(i,j)
!                          else
!                             snoeqv=0.
!                          end if
!                          if(i==ii.and.j==jj)print*,'sno,itype,ivgtyp B cing snfrc = ',  &
!                                                     snoeqv,itype,IVGTYP(I,J)
!                          if(sm(i,j) > 0.1)then
!                             sfcpct(4)=0.
!                          else 
!                             call snfrac_gfs(SNOeqv,IVGTYP(I,J),snofrac)
!                             sfcpct(4)=snofrac
!                          end if
!                          if(i==ii.and.j==jj)print*,'sno,itype,ivgtyp,sfcpct(4) = ',     &
!                                                     snoeqv,itype,IVGTYP(I,J),sfcpct(4)
!                       else if(MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR')then
!                          sfcpct(4)=pctsno(i,j)
!                       else          
!                          itype=IVGTYP(I,J)
!                          IF(itype == 0)itype=8
!                          CALL SNFRAC (SNO(I,J),IVGTYP(I,J),snofrac)
!	                  sfcpct(4)=snofrac
!                       end if 
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
                       if(ivegsrc==0)then
                         itype = min(max(0,ivgtyp(i,j)),novegtype)
                       else
                         itype = min(max(1,ivgtyp(i,j)),novegtype)
                       end if
                       surface(1)%land_type = model_to_crtm(itype)
                       
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
                       end if
                       if(i==ii.and.j==jj.and.debugprint)print*,'sample surface in CALRAD=', &
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
                       if(i==ii.and.j==jj.and.debugprint)print*,'TOA= ',atmosphere(1)%level_pressure(0)
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
                          if(i==ii.and.j==jj.and.debugprint)print*,'sample atmosphere in CALRAD=',  &
   	                        i,j,k,atmosphere(1)%level_pressure(k),atmosphere(1)%pressure(k),  &
                                atmosphere(1)%temperature(k),atmosphere(1)%absorber(k,1),  &
                                atmosphere(1)%absorber(k,2)
                          ! Specify clouds
                          dpovg=(pint(i,j,k+1)-pint(i,j,k))/g !crtm uses column integrated field
                          if(imp_physics==99 .or. imp_physics==98)then
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
                             atmosphere(1)%cloud(2)%effective_radius(k) = 75.
                             atmosphere(1)%cloud(2)%water_content(k) = max(0.,qqi(i,j,k)*dpovg)
                             RHOX=1000.
                             RHO=pmid(i,j,k)/(RD*T(I,J,K)*(1.+D608*Q(I,J,K)))
                             atmosphere(1)%cloud(3)%effective_radius(k) = 0.
                             if(NRAIN(I,J,K)>0.) &
                                   atmosphere(1)%cloud(3)%effective_radius(k) =   &
                                   1.0E6*1.5*(RHO*qqr(i,j,k)/(PI*RHOX*NRAIN(I,J,K)))**(1./3.)
                             atmosphere(1)%cloud(3)%water_content(k) = max(0.,qqr(i,j,k)*dpovg)
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

                          else if(imp_physics==8 .or. imp_physics==6 .or. imp_physics==2 .or. & 
                                  imp_physics==28 .or. imp_physics==11)then
                             atmosphere(1)%cloud(1)%water_content(k)=max(0.,qqw(i,j,k)*dpovg)
                             atmosphere(1)%cloud(2)%water_content(k)=max(0.,qqi(i,j,k)*dpovg)
                             atmosphere(1)%cloud(3)%water_content(k)=max(0.,qqr(i,j,k)*dpovg)
                             atmosphere(1)%cloud(4)%water_content(k)=max(0.,qqs(i,j,k)*dpovg)
                             atmosphere(1)%cloud(5)%water_content(k)=max(0.,qqg(i,j,k)*dpovg)
                             atmosphere(1)%cloud(1)%effective_radius(k)=effr(pmid(i,j,k),t(i,j,k), &
                             q(i,j,k),qqw(i,j,k),qqi(i,j,k),qqr(i,j,k),f_rimef(i,j,k),nlice(i,j,k), &
                             nrain(i,j,k),qqs(i,j,k),qqg(i,j,k),qqnr(i,j,k),qqni(i,j,k),qqnw(i,j,k),imp_physics,'C')
                             atmosphere(1)%cloud(2)%effective_radius(k)=effr(pmid(i,j,k),t(i,j,k), &
                             q(i,j,k),qqw(i,j,k),qqi(i,j,k),qqr(i,j,k),f_rimef(i,j,k),nlice(i,j,k), &
                             nrain(i,j,k),qqs(i,j,k),qqg(i,j,k),qqnr(i,j,k),qqni(i,j,k),qqnw(i,j,k),imp_physics,'I')
                             atmosphere(1)%cloud(3)%effective_radius(k)=effr(pmid(i,j,k),t(i,j,k), &
                             q(i,j,k),qqw(i,j,k),qqi(i,j,k),qqr(i,j,k),f_rimef(i,j,k),nlice(i,j,k), &
                             nrain(i,j,k),qqs(i,j,k),qqg(i,j,k),qqnr(i,j,k),qqni(i,j,k),qqnw(i,j,k),imp_physics,'R')
                             atmosphere(1)%cloud(4)%effective_radius(k)=effr(pmid(i,j,k),t(i,j,k), &
                             q(i,j,k),qqw(i,j,k),qqi(i,j,k),qqr(i,j,k),f_rimef(i,j,k),nlice(i,j,k), &
                             nrain(i,j,k),qqs(i,j,k),qqg(i,j,k),qqnr(i,j,k),qqni(i,j,k),qqnw(i,j,k),imp_physics,'S')
                             atmosphere(1)%cloud(5)%effective_radius(k)=effr(pmid(i,j,k),t(i,j,k), &
                             q(i,j,k),qqw(i,j,k),qqi(i,j,k),qqr(i,j,k),f_rimef(i,j,k),nlice(i,j,k), &
                             nrain(i,j,k),qqs(i,j,k),qqg(i,j,k),qqnr(i,j,k),qqni(i,j,k),qqnw(i,j,k),imp_physics,'G')
                          end if 
                       end do
!Meng 09/2018 modify two model layer having identical pressure
                       do k = 1,lm-1
                          if (abs(atmosphere(1)%pressure(k)-atmosphere(1)%pressure(k+1)) &
                              < 0.005) then
                             atmosphere(1)%pressure(k)=atmosphere(1)%pressure(k)+0.005
                          endif
                       enddo

                       !bsf - start
                       !-- Add subgrid-scale convective clouds for WRF runs
                       if(icu_physics==2) then
                          lcbot=nint(hbot(i,j))
                          lctop=nint(htop(i,j))
                          if (lcbot-lctop > 1) then
                             !-- q_conv = assumed grid-averaged cloud water/ice condensate from conimp_physics,vection (Cu)
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
                       if (grib=="grib2") then
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
                             grid1(i,j) = tb(i,j,ichan)
                          enddo
                       enddo
                       if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              end if  ! end of outputting trmm

              if (isis=='imgr_g11')then  ! writing goes 11 to grib
                 do ixchan=1,4
                    ichan=ixchan
                    igot=445+ixchan
                    if(igot>0) then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j) = tb(i,j,ichan)
                          enddo
                       enddo
                       if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif ! IGOT
                 enddo
              end if  ! end of outputting goes 11

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
                        if (grib=="grib2") then
                           cfld=cfld+1
                           fld_info(cfld)%ifld=IAVBLFLD(igot)
                           datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                        endif
                    endif
                 enddo
               endif ! end of outputting goes 12
              if (isis=='abi_gr')then  ! writing goes-r nadir to grib2
                 nc=0
                 do ixchan=1,10
                   ichan=ixchan
                   igot=iget(957+ixchan)
                   if(igot>0)then
                    do j=jsta,jend
                     do i=1,im
                      grid1(i,j)=tb(i,j,ichan)
                     enddo
                    enddo
                    if(grib=="grib2" )then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(igot)
                     datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                    endif
                   endif
                 enddo ! channel loop
              end if  ! end of outputting goes-r nadir


           end if nadir ! end if for computing nadir simulated radiance    

           ! run crtm for non-nadir instruments / channels 

           nonnadir: if((isis=='ssmi_f13' .and. iget(800) > 0 ) .OR. &
                        (isis=='ssmi_f14' .and. iget(806) > 0 ) .OR. &
                        (isis=='ssmi_f15' .and. iget(812) > 0 ) .OR. &
                        (isis=='ssmis_f16' .and. iget(818) > 0) .OR. &
                        (isis=='ssmis_f17' .and. iget(825) > 0) .OR. &
                        (isis=='ssmis_f18' .and. iget(832) > 0) .OR. &
                        (isis=='ssmis_f19' .and. iget(839) > 0) .OR. &
                        (isis=='ssmis_f20' .and. iget(846) > 0) .OR. &
                        (isis=='imgr_mt2' .and. iget(860)>0) .OR. &
                        (isis=='imgr_mt1r' .and. iget(864)>0) .OR. &
                        (isis=='imgr_insat3d' .and. iget(865)>0) .OR. &
                        (isis=='imgr_g13' .and. iget(868)>0) .OR. &
                        (isis=='imgr_g15' .and. iget(872)>0) .OR. &
                        (isis=='abi_g16'  .and. post_abig16) .OR. &
                        (isis=='abi_g17'  .and. post_abig17) .OR. &
                        (isis=='seviri_m10' .and. iget(876)>0) .OR. &
                        (isis=='ahi_himawari8' .and. post_ahi8) .OR. &
                        (isis=='imgr_g12' .and. (iget(456)>0 .or. &
                        iget(457)>0 .or. iget(458)>0 .or. iget(459)>0)) .or. &
                        (isis=='imgr_g11' .and. (iget(460)>0 .or. &
                        iget(461)>0 .or. iget(462)>0 .or. iget(463)>0)))then

              do j=jsta,jend
                 do i=1,im
                    !    Load geometry structure
                    !    geometryinfo(1)%sensor_zenith_angle = zasat*rtd  ! local zenith angle ???????
                    ! compute satellite zenith angle
                    if(isis=='imgr_g12')then
                       sublat=0.0
                       sublon=-75.0
                    else if(isis=='seviri_m10')then
                       sublat=0.0
                       sublon=0.0
                    else if(isis=='imgr_g13')then
                       sublat=0.0
                       sublon=-75.0
                    else if(isis=='imgr_g15')then
                       sublat=0.0
                       sublon=-135.0
                    else if(isis=='abi_g16')then  ! positions should be controlled by runtime setting or fix file
                       sublat=0.0
                       sublon=-75.2
                    else if(isis=='abi_g17')then
                       sublat=0.0
                       sublon=-137.2
                    else if(isis=='imgr_g11')then
                       sublat=0.0
                       sublon=-135.0
                    else if(isis=='imgr_mt2') then
                       sublat=0.0
                       sublon=145.0
                    else if(isis=='imgr_mt1r') then
                       sublat=0.0
                       sublon=140.0
                    else if(isis=='imgr_insat3d') then
                       sublat=0.0
                       sublon=74.0
                    else if(isis=='ahi_himawari8') then
                       sublat=0.0
                       sublon=140.7
                    end if

!                   use zenith angle = 53.1 for SSMI and SSMIS:
                    if(isis=='ssmis_f16'.or.isis=='ssmis_f17'.or.isis=='ssmis_f18'.or. &
                       isis=='ssmis_f19'.or.isis=='ssmis_f20'.or.isis=='ssmi_f13'.or. &
                       isis=='ssmi_f14'.or.isis=='ssmi_f15')then
                       sat_zenith=53.1
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

                       if(MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR')then
                          sfcpct(4)=pctsno(i,j)
                       else if(ivegsrc==1)then
                          itype=IVGTYP(I,J)
                          IF(itype == 0)itype=8
                          if(sno(i,j)<spval)then
                             snoeqv=sno(i,j)
                          else
                             snoeqv=0.
                          end if
                          CALL SNFRAC (SNO(I,J),IVGTYP(I,J),snofrac)
                          sfcpct(4)=snofrac
                       else if(ivegsrc==2)then
                          itype=IVGTYP(I,J)
                          itype = min(max(0,ivgtyp(i,j)),13)
                          if(sno(i,j)<spval)then
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
                       !DTC added based on nadir section
                       ! Chuang: for igbp type 15 (snow/ice), the main type
                       ! needs to be set to ice or snow
                       ! to prevent crtm forward model from failing
                       if(novegtype==20 .and. itype==15 .and.sfcpct(4)<1.0_r_kind)then
                          if(debugprint)print*,'changing land type for veg type 15',i,j,itype,sfcpct(1:4)
                         sfcpct(1)=0.0_r_kind
                         sfcpct(2)=0.0_r_kind
                         sfcpct(3)=0.0_r_kind
                         sfcpct(4)=1.0_r_kind
                          !print*,'change main land type to snow for veg type 15
                          !',i,j
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
                       if(ivegsrc==0)then
                         itype = min(max(0,ivgtyp(i,j)),novegtype)
                       else
                         itype = min(max(1,ivgtyp(i,j)),novegtype)
                       end if
                       surface(1)%land_type = model_to_crtm(itype)

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
                          if(imp_physics==99 .or. imp_physics==98)then
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
                             atmosphere(1)%cloud(2)%effective_radius(k) = 75.
                             atmosphere(1)%cloud(2)%water_content(k) = max(0.,qqi(i,j,k)*dpovg)
                             RHOX=1000.
                             RHO=pmid(i,j,k)/(RD*T(I,J,K)*(1.+D608*Q(I,J,K)))
                             atmosphere(1)%cloud(3)%effective_radius(k) = 0.
                             if(NRAIN(I,J,K)>0.) &
                                   atmosphere(1)%cloud(3)%effective_radius(k) =   &
                                   1.0E6*1.5*(RHO*qqr(i,j,k)/(PI*RHOX*NRAIN(I,J,K)))**(1./3.)
                             atmosphere(1)%cloud(3)%water_content(k) = max(0.,qqr(i,j,k)*dpovg)
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

                          else if(imp_physics==8 .or. imp_physics==6 .or. imp_physics==2 .or. &
                                  imp_physics==28 .or.  imp_physics==11)then
                             atmosphere(1)%cloud(1)%water_content(k)=max(0.,qqw(i,j,k)*dpovg)
                             atmosphere(1)%cloud(2)%water_content(k)=max(0.,qqi(i,j,k)*dpovg)
                             atmosphere(1)%cloud(3)%water_content(k)=max(0.,qqr(i,j,k)*dpovg)
                             atmosphere(1)%cloud(4)%water_content(k)=max(0.,qqs(i,j,k)*dpovg)
                             atmosphere(1)%cloud(5)%water_content(k)=max(0.,qqg(i,j,k)*dpovg)
                             atmosphere(1)%cloud(1)%effective_radius(k)=effr(pmid(i,j,k),t(i,j,k), &
                             q(i,j,k),qqw(i,j,k),qqi(i,j,k),qqr(i,j,k),f_rimef(i,j,k),nlice(i,j,k), &
                             nrain(i,j,k),qqs(i,j,k),qqg(i,j,k),qqnr(i,j,k),qqni(i,j,k),qqnw(i,j,k),imp_physics,'C')
                             atmosphere(1)%cloud(2)%effective_radius(k)=effr(pmid(i,j,k),t(i,j,k), &
                             q(i,j,k),qqw(i,j,k),qqi(i,j,k),qqr(i,j,k),f_rimef(i,j,k),nlice(i,j,k), &
                             nrain(i,j,k),qqs(i,j,k),qqg(i,j,k),qqnr(i,j,k),qqni(i,j,k),qqnw(i,j,k),imp_physics,'I')
                             atmosphere(1)%cloud(3)%effective_radius(k)=effr(pmid(i,j,k),t(i,j,k), &
                             q(i,j,k),qqw(i,j,k),qqi(i,j,k),qqr(i,j,k),f_rimef(i,j,k),nlice(i,j,k), &
                             nrain(i,j,k),qqs(i,j,k),qqg(i,j,k),qqnr(i,j,k),qqni(i,j,k),qqnw(i,j,k),imp_physics,'R')
                             atmosphere(1)%cloud(4)%effective_radius(k)=effr(pmid(i,j,k),t(i,j,k), &
                             q(i,j,k),qqw(i,j,k),qqi(i,j,k),qqr(i,j,k),f_rimef(i,j,k),nlice(i,j,k), &
                             nrain(i,j,k),qqs(i,j,k),qqg(i,j,k),qqnr(i,j,k),qqni(i,j,k),qqnw(i,j,k),imp_physics,'S')
                             atmosphere(1)%cloud(5)%effective_radius(k)=effr(pmid(i,j,k),t(i,j,k), &
                             q(i,j,k),qqw(i,j,k),qqi(i,j,k),qqr(i,j,k),f_rimef(i,j,k),nlice(i,j,k), &
                             nrain(i,j,k),qqs(i,j,k),qqg(i,j,k),qqnr(i,j,k),qqni(i,j,k),qqnw(i,j,k),imp_physics,'G')
                          end if 
                       end do
!Meng 09/2018 modify two model layer having identical pressure
                       do k = 1,lm-1
                          if (abs(atmosphere(1)%pressure(k)-atmosphere(1)%pressure(k+1)) &
                              < 0.005) then
                             atmosphere(1)%pressure(k)=atmosphere(1)%pressure(k)+0.005
                          endif
                       enddo
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

              if (isis=='ssmi_f13')then  ! writing ssmi to grib (37 & 85 GHz)
              nc=0
              do ixchan=1,7
                ichan=ixchan
                igot=iget(800)
                if(igot>0) then
                if(lvls(ixchan,igot).eq.1)then
                  nc=nc+1
                  do j=jsta,jend
                    do i=1,im
                      grid1(i,j)=tb(i,j,nc)
                    enddo
                  enddo
                  if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                  endif
                 endif
                endif
              enddo
              end if  ! end of outputting ssmi f13
              if (isis=='ssmi_f14')then  ! writing ssmi to grib (19,37 & 85 GHz)
              nc=0
              do ixchan=1,7
                ichan=ixchan
                igot=iget(806)
                if(igot>0) then
                if(lvls(ixchan,igot).eq.1)then
                  nc=nc+1
                  do j=jsta,jend
                    do i=1,im
                      grid1(i,j)=tb(i,j,nc)
                    enddo
                  enddo
                  if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                  endif
                 endif
                endif
              enddo
              end if  ! end of outputting ssmi f14
              if (isis=='ssmi_f15')then  ! writing ssmi to grib (19,37 & 85 GHz)
              nc=0
              do ixchan=1,7
                ichan=ixchan
                igot=iget(812)
                if(igot>0) then
                if(lvls(ixchan,igot).eq.1)then
                  nc=nc+1
                  do j=jsta,jend
                    do i=1,im
                      grid1(i,j)=tb(i,j,nc)
                    enddo
                  enddo
                  if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                  endif
                 endif
                endif
              enddo
              end if  ! end of outputting ssmi f15
              if (isis=='ssmis_f16')then  ! writing ssmis to grib (183,19,37 & 85GHz)
              nc=0
              do ixchan=1,24
                ichan=ixchan
                igot=iget(818)
                if(igot>0) then
                print*,'ixchan,lvls=',ixchan,lvls(ixchan,igot)
                if(lvls(ixchan,igot).eq.1)then
                  nc=nc+1
                  do j=jsta,jend
                    do i=1,im
                      grid1(i,j)=tb(i,j,nc)
                    enddo
                  enddo
                  if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                  endif
                 endif
                endif
              enddo
              end if  ! end of outputting ssmis f16
              if (isis=='ssmis_f17')then  ! writing ssmis to grib (183,19,37 &85GHz)
              nc=0
              do ixchan=1,24
                ichan=ixchan
                igot=iget(825)
                if(igot>0) then
                if(lvls(ixchan,igot).eq.1)then
                  nc=nc+1
                  do j=jsta,jend
                    do i=1,im
                      grid1(i,j)=tb(i,j,nc)
                    enddo
                  enddo
                  if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                  endif
                 endif
                endif
              enddo
              end if  ! end of outputting ssmis f17
              if (isis=='ssmis_f18')then  ! writing ssmis to grib (183,19,37 &85GHz)
              nc=0
              do ixchan=1,24
                ichan=ixchan
                igot=iget(832)
                if(igot>0) then
                if(lvls(ixchan,igot).eq.1)then
                  nc=nc+1
                  do j=jsta,jend
                    do i=1,im
                      grid1(i,j)=tb(i,j,nc)
                    enddo
                  enddo
                  if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                  endif
                 endif
                endif
              enddo
              end if  ! end of outputting ssmis f18
              if (isis=='ssmis_f19')then  ! writing ssmis to grib (183,19,37 &85GHz)
              nc=0
              do ixchan=1,24
                ichan=ixchan
                igot=iget(839)
                if(igot>0) then
                if(lvls(ixchan,igot).eq.1)then
                  nc=nc+1
                  do j=jsta,jend
                    do i=1,im
                      grid1(i,j)=tb(i,j,nc)
                    enddo
                  enddo
                  if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                  endif
                 endif
                endif
              enddo
              end if  ! end of outputting ssmis f19
              if (isis=='ssmis_f20')then  ! writing ssmis to grib (183,19,37 &85GHz)
              nc=0
              do ixchan=1,24
                ichan=ixchan
                igot=iget(846)
                if(igot>0) then
                if(lvls(ixchan,igot).eq.1)then
                  nc=nc+1
                  do j=jsta,jend
                    do i=1,im
                      grid1(i,j)=tb(i,j,nc)
                    enddo
                  enddo
                  if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                  endif
                 endif
                endif
              enddo
              end if  ! end of outputting ssmis f20
              if(isis=='imgr_mt2') then ! writing MTSAT-2 to grib
                 nc=0
                 do ichan=1,4
                    igot=iget(860)
                      if(lvls(ichan,igot).eq.1)then
                       nc=nc+1
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,nc)
                          enddo
                       enddo
                       if(grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              endif
              if(isis=='imgr_mt1r') then ! writing MTSAT-1r to grib
                 nc=0
                 do ichan=1,4
                    igot=iget(864) 
                      if(lvls(ichan,igot).eq.1)then
                       nc=nc+1
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,nc)
                          enddo
                       enddo
                       if(grib=="grib2" )then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              endif
              if_insat3d: if(isis=='imgr_insat3d') then ! writing MTSAT-1r to grib
                 nc=0
                 do ichan=1,4
                    igot=iget(865) 
                      if(lvls(ichan,igot).eq.1)then
                       nc=nc+1
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,nc)
                          enddo
                       enddo
                       if(grib=="grib2" )then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              endif if_insat3d
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
                       if(grib=="grib2" )then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                enddo
              endif ! end of outputting goes 11
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
                       if(grib=="grib2" )then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              end if  ! end of outputting goes 12
              if (isis=='seviri_m10')then  ! writing msg/severi 10
                 nc=0
                 do ixchan=1,8
                   ichan=ixchan
                   igot=iget(876)
                   if(igot>0) then
                   if(lvls(ixchan,igot).eq.1)then
                    nc=nc+1
                    do j=jsta,jend
                     do i=1,im
                      grid1(i,j)=tb(i,j,nc)
                     enddo
                    enddo
                    if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                    endif
                   endif
                   endif
                 enddo
              end if  ! end of outputting msg/seviri 10
              if (isis=='imgr_g13')then  ! writing goes 13 to grib
                 nc=0
                 do ixchan=1,4
                   ichan=ixchan
                   igot=iget(868)
                   if(igot>0) then
                   if(lvls(ixchan,igot).eq.1)then
                    nc=nc+1
                    do j=jsta,jend
                     do i=1,im
                      grid1(i,j)=tb(i,j,nc)
                     enddo
                    enddo
                    if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                    endif
                   endif
                   endif
                 enddo
              end if  ! end of outputting goes 13
              if (isis=='imgr_g15')then  ! writing goes 15 to grib
                 nc=0
                 do ixchan=1,4
                   ichan=ixchan
                   igot=iget(872)
                   if(igot>0) then
                   if(lvls(ixchan,igot).eq.1)then
                    nc=nc+1
                    do j=jsta,jend
                     do i=1,im
                      grid1(i,j)=tb(i,j,nc)
                     enddo
                    enddo
                    if (grib=="grib2") then
                          cfld=cfld+1
                          fld_info(cfld)%ifld=IAVBLFLD(igot)
                          datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                    endif
                   endif
                   endif
                 enddo
              end if  ! end of outputting goes 15
              if (isis=='abi_g16')then  ! writing goes 16 to grib
                 nc=0
                 do ixchan=1,10
                   ichan=ixchan
                   igot=iget(926+ixchan)
                   if(igot>0)then
                    do j=jsta,jend
                     do i=1,im
                      grid1(i,j)=tb(i,j,ichan)
                     enddo
                    enddo
                    if(grib=="grib2" )then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(igot)
                     datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                    endif
                   endif
                 enddo ! channel loop
              end if  ! end of outputting goes 16
              if (isis=='abi_g17')then  ! writing goes 16 to grib
                 nc=0
                 do ixchan=1,10
                   ichan=ixchan
                   igot=iget(936+ixchan)
                   if(igot>0)then
                    do j=jsta,jend
                     do i=1,im
                      grid1(i,j)=tb(i,j,ichan)
                     enddo
                    enddo
                    if(grib=="grib2" )then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(igot)
                     datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                    endif
                   endif
                 enddo ! channel loop
              end if  ! end of outputting goes 17
              if(isis=='ahi_himawari8') then ! writing Himawari-8 AHI to grib
                 do ichan=1,10
                    igot=iget(968+ichan)
                      if(igot>0)then
                       do j=jsta,jend
                          do i=1,im
                             grid1(i,j)=tb(i,j,ichan)
                          enddo
                       enddo
                       if(grib=="grib2" )then
                        cfld=cfld+1
                        fld_info(cfld)%ifld=IAVBLFLD(igot)
                        datapd(1:im,1:jend-jsta+1,cfld)=grid1(1:im,jsta:jend)
                       endif
                    endif
                 enddo
              endif ! end of outputting himawari-8 ahi

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

REAL FUNCTION EFFR(pmid,t,q,qqw,qqi,qqr,f_rimef, nlice, nrain, &
                   qqs,qqg,qqnr,qqni,qqnw,mp_opt,species)

!       JASON OTKIN AND WILLIAM LEWIS
!       09 DECEMBER 2014
!       Greg Thompson, 20200924

  use params_mod, only: pi, rd, d608, rg

        implicit none

        real :: pmid,t,q,qqw,qqi,qqr,qqs,qqg,f_rimef,nlice,nrain
        real :: qqnr,qqni,qqnw
        character(LEN=1) :: species

        integer                         :: n,count,count1,mp_opt
        real :: rho, ncc, rhox
        real :: n0_s, n0_r, n0_g
        real :: lambdar, lambdas, lambdag

!-------------------------------------------------------------------------------
!  GAMMA FUNCTION & RELATED VARIABLES
!-------------------------------------------------------------------------------

        real :: gamma
        real :: gamma_crg, gamma_s
!       real :: gamma_i

        real :: WGAMMA, GAMMLN

        real    :: rc,am_c,bm_c,cce(3,15),ccg(3,15),ocg1(15),ocg2(15)
        integer :: nu_c

        real, dimension(0:15), parameter:: g_ratio = (/6,24,60,120,210, &
     &              336,504,720,990,1320,1716,2184,2730,3360,4080,4896/)

        real    :: rr, mu_r, am_r, bm_r, cre(3), crg(3), ore1, org1, org2
        real    :: mvd_r, ron_sl, ron_r0, ron_c0, ron_c1, ron_c2, obmr

        real    :: ri, mu_i, am_i, bm_i, cie(3), cig(3), oig1, oig2, obmi

        real    :: rs, am_s, oams, cse(3)
        real    :: loga, a, b, tc0, smob, smo2, smoc
        REAL, PARAMETER:: mu_s = 0.6357
        REAL, PARAMETER:: Kap0 = 490.6
        REAL, PARAMETER:: Kap1 = 17.46
        REAL, PARAMETER:: Lam0 = 20.78
        REAL, PARAMETER:: Lam1 = 3.29

!-------------------------------------------------------------------------------
!  MINIMUM/MAXIMUM CONSTANTS FOR ALL SCHEMES
!-------------------------------------------------------------------------------

        real, parameter :: eps=0.622, beta_crg=3., beta_i=2.,beta_s=2.4

        real, parameter :: min_qc=1.e-7, min_qr=1.e-7, min_qi=1.e-8,min_qs=1.e-8, min_qg=1.e-7
        real, parameter :: min_c=2.e-6,  min_r=20.e-6, min_i=4.e-6,min_s=20.e-6, min_g=20.e-6
        real, parameter :: max_c=1.e-2,  max_r=1.e-2,  max_i=1.e-3,max_s=2.e-2,  max_g=5.e-0

        real    :: am_g, bm_g, mu_g
        real    :: cgg(3), cge(3), oge1, obmg, ogg1, ogg2

        real    :: ygra1, zans1, rg2
        double precision :: no_exp, no_min, lm_exp, lamg, lamc, lamr, lami, lams

!-------------------------------------------------------------------------------
!  WSM6-SPECIFIC ARRAYS
!-------------------------------------------------------------------------------

        real :: wsm6_nci, xmi, xmitemp

!-------------------------------------------------------------------------------
!  CONSTANTS FOR WSM6 MICROPHYSICS SCHEME
!-------------------------------------------------------------------------------

        real, parameter :: wsm6_cnp=3.e8, wsm6_rhor=1000.
        real, parameter :: wsm6_rhos=100., wsm6_rhog=500.
        real, parameter :: wsm6_dimax=500.e-6, wsm6_dicon=11.9
        real, parameter :: wsm6_alpha=.12, wsm6_t0c=273.15
        real, parameter :: wsm6_n0s=2.e6, wsm6_n0smax=1.e11
        real, parameter :: wsm6_n0r=8.e6, wsm6_n0g=4.e6
        real, parameter :: wsm6_qmin=1.e-15

!-------------------------------------------------------------------------------
!  CONSTANTS FOR LIN MICROPHYSICS SCHEME
!-------------------------------------------------------------------------------

        real, parameter :: lin_rhoi=100., lin_rhor=1000., lin_rhos=100.
        real, parameter :: lin_rhog=400., lin_cnp=3.e8
        real, parameter :: lin_n0r=8.e6,  lin_n0s=3.e6,   lin_n0g=4.e6

!-------------------------------------------------------------------------------
!  CONSTANTS FOR NEW THOMPSON MICROPHYSICS SCHEME (FOR WRF VERSIONS 3.1 AND UP)
!-------------------------------------------------------------------------------

        real, parameter :: nthom_rhor=1000., nthom_rhos=100.
!       WM LEWIS updated rhog to 500 from 400
        real, parameter :: nthom_rhog=500.,  nthom_rhoi=890.
        real, parameter :: nthom_gon_min=1.e4, nthom_gon_max=3.e6
        real, parameter :: nthom_nt_c=100.e6

        real, parameter :: nthom_min_nci=5.e2
        real, parameter :: nthom_min_ncr=1.e-6

        real, parameter :: nthom_bm_s=2.0               !this is important

        real :: nci2, ncc2, ncr2

        real, dimension(10), parameter :: &
        nthom_sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
                      0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
        real, dimension(10), parameter :: &
        nthom_sb = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
                      0.060366,  0.000079,  0.000594, 0.0,      -0.003577/)

!-------------------------------------------------------------------------------
!  CONSTANTS FOR GFDL MICROPHYSICS SCHEME - which is Lin for precip clouds
!-------------------------------------------------------------------------------

        real, parameter :: gfdl_rhoi=100., gfdl_rhor=1000., gfdl_rhos=100.
        real, parameter :: gfdl_rhog=400., gfdl_cnp=3.e8
        real, parameter :: gfdl_tice = 273.16

        real, parameter :: gfdl_qmin = 1.0e-5, gfdl_ccn = 1.0e8, gfdl_beta = 1.22
        real, parameter :: gfdl_gammar = 17.837789, gfdl_gammas = 8.2850630, gfdl_gammag = 11.631769
        real, parameter :: gfdl_alphar = 0.8, gfdl_alphas = 0.25, gfdl_alphag = 0.5
        real, parameter :: gfdl_n0r=8.e6, gfdl_n0s=3.e6, gfdl_n0g=4.e6

        real, parameter :: gfdl_rewmin = 5.0,  gfdl_rewmax = 10.0
        real, parameter :: gfdl_reimin = 10.0, gfdl_reimax = 150.0
        real, parameter :: gfdl_rermin = 0.0,  gfdl_rermax = 10000.0
        real, parameter :: gfdl_resmin = 0.0,  gfdl_resmax = 10000.0
        real, parameter :: gfdl_regmin = 0.0,  gfdl_regmax = 10000.0



        if(mp_opt.eq.6) then                        !WSM6 SCHEME

          n0_r = wsm6_n0r
          n0_g = wsm6_n0g
          n0_s = wsm6_n0s

        elseif(mp_opt.eq.2)then                     !LIN SCHEME

          n0_r = lin_n0r
          n0_g = lin_n0g
          n0_s = lin_n0s
   
        endif

        gamma_crg = 6.0      ! gamma(1.0 + beta_crg)
        gamma_s = 2.981134   ! gamma(1.0 + beta_s)
!       gamma_i = 2.0        ! gamma(1.0 + beta_i)

!------------------------------------------------------------------------------
!  SET DIAMETER ARRAYS TO ZERO, COMPUTE DENSITY
!------------------------------------------------------------------------------

        effr=0. 

        rho=pmid/(rd*t*(1.+D608*q))


 if(mp_opt.eq.6)then

     SELECT CASE(species)

     CASE("C")

     if ( qqw.gt.min_qc ) then !cloud diameter: assume constant # concentration
       effr = 1.0E6*(( 6. * rho * qqw ) / &
       (pi * wsm6_rhor * wsm6_cnp))**(1/3.)

     endif

     CASE("R")

     if ( qqr.gt.min_qr ) then !rain diameter: assume gamma distribution
       effr = 1.0E6*( ( 6. * rho * qqr ) / &
       ( pi * wsm6_rhor * n0_r * gamma_crg ) ) ** (1/(1+beta_crg ) )
     endif

     CASE("G")

     if ( qqg.gt.min_qg ) then !graupel diameter: assume gamma distribution
       effr = 1.0E6*( ( 6. * rho * qqg ) / &
       ( pi * wsm6_rhog * n0_g * gamma_crg ) ) ** (1/(1+beta_crg ) )
     endif

     CASE("S")

     if ( qqs.gt.min_qs ) then !snow diameter: assume gamma distribution
       effr = 1.0E6*( ( 6. * rho * qqs ) / &
       ( pi * wsm6_rhos * n0_s * gamma_s   ) ) ** ( 1/(1+beta_s) )
     endif

!  ICE DIAMETER: CALCULATED USING METHOD OUTLINED IN WRF BROWSER.  Refer to
!  phys/module_mp_wsm6.F (Vice:fallout of ice crystal).

     CASE("I")

     if ( qqi.gt.min_qi ) then !ice diameter
!       wsm6_nci = min(max(5.38e7*(rho*max(qqi,wsm6_qmin)),1.e3),1.e6)
!       xmi = rho * qqi / wsm6_nci
!       effr = 1.0E6*min( sqrt(xmi), wsm6_dimax)
!!      from wsm6, HWRF ver 3.6:
!!          temp = (den(i,k)*max(qci(i,k,2),qmin))
!!          temp = sqrt(sqrt(temp*temp*temp))
!!          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
!!      diameter  = max(min(dicon * sqrt(xmi),dimax), 1.e-25)
       xmitemp=rho*max(qqi,wsm6_qmin)
       xmitemp=sqrt(sqrt(xmitemp*xmitemp*xmitemp))
       xmi= min(max(5.38e7*xmitemp,1.e3),1.e6)
       effr = 1.0E6*max(min(wsm6_dicon * sqrt(xmi),wsm6_dimax), 1.e-25)
     endif

     END SELECT 

 elseif(mp_opt.eq.2)then

     SELECT CASE(species)

     CASE("C")

     if ( qqw > min_qc ) then !cloud diameter: assume constant # concentration
       effr = 1.0E6*(( 6. * rho * qqw ) / &
       (pi * lin_rhor * lin_cnp))**(1/3.)
     endif

     CASE("R")

     if ( qqr > min_qr ) then !rain diameter: assume gamma distribution
       effr = 1.0E6*( ( 6. * rho * qqr ) / &
       ( pi * lin_rhor * n0_r * gamma_crg ) ) ** (1/(1+beta_crg ) )
     endif

     CASE("I")

     if ( qqi > min_qi ) then !ice diameter: assume constant # concentrtion
       effr = 1.0E6*( ( 6. * rho * qqi ) / &
       ( pi * lin_rhoi * lin_cnp ) ) ** ( 1/3.)
     endif

     CASE("S")

     if ( qqs > min_qs ) then !snow diameter: assume gamma distribution
       effr = 1.0E6*( ( 6. * rho * qqs ) / &
       ( pi * lin_rhos * n0_s * gamma_s   ) ) ** ( 1/(1+beta_s) )
     endif

     CASE("G")

     if ( qqg > min_qg ) then !graupel diameter: assume gamma distribution
       effr = 1.0E6*( ( 6. * rho * qqg ) / &
       ( pi * lin_rhog * n0_g * gamma_crg ) ) ** (1/(1+beta_crg ) )
     endif

     END SELECT

 elseif(mp_opt.eq.8 .or. mp_opt.eq.28)then

!  rain section

          bm_r   = 3.0
          mu_r   = 0.0
          obmr   = 1.0 / bm_r
          am_r   = pi * nthom_rhor / 6.0

          cre(1) = bm_r + 1.
          cre(2) = mu_r + 1.
          cre(3) = bm_r + mu_r + 1.

          crg(1) = WGAMMA(cre(1))
          crg(2) = WGAMMA(cre(2))
          crg(3) = WGAMMA(cre(3))

          ore1   = 1. / cre(1)
          org1   = 1. / crg(1)
          org2   = 1. / crg(2)

!  cloud section

          bm_c   = bm_r

          do n = 1, 15
             cce(1,n) = n + 1.             ! Substitute variable value of mu_c
             cce(2,n) = bm_c + n + 1.      ! Substitute variable value of mu_c

             ccg(1,n) = WGAMMA(cce(1,n))
             ccg(2,n) = WGAMMA(cce(2,n))

             ocg1(n)   = 1./ccg(1,n)
             ocg2(n)   = 1./ccg(2,n)
          enddo

!  ice section

          am_i   = pi * nthom_rhoi / 6.0
          bm_i   = 3.0
          mu_i   = 0.

          cie(1) = mu_i + 1.
          cie(2) = bm_i + mu_i + 1.

          cig(1) = WGAMMA(cie(1))
          cig(2) = WGAMMA(cie(2))

          oig1   = 1./cig(1)
          oig2   = 1./cig(2)
          obmi   = 1./bm_i

!  snow section

          am_s   = 0.069

          oams   = 1./am_s

          cse(1) = nthom_bm_s + 1.

!  graupel section
          bm_g   = 3.0
          mu_g   = 0.0
          obmg   = 1.0 / bm_g
          am_g   = pi * nthom_rhog / 6.0

          cge(1) = bm_g + 1.
          cge(2) = mu_g + 1.
          cge(3) = bm_g + mu_g + 1.

          cgg(1) = WGAMMA(cge(1))
          cgg(2) = WGAMMA(cge(2))
          cgg(3) = WGAMMA(cge(3))

          oge1   = 1. / cge(1)
          ogg1   = 1. / cgg(1)
          ogg2   = 1. / cgg(2)

!CLOUD DIAMETER CALCULATIONS

     SELECT CASE (species)

     CASE("C")

            if(qqw .ge. min_qc) then

              rc = MAX(1.E-12, qqw * rho)

              if (mp_opt.eq.8) then
                 ncc2 = nthom_nt_c
              elseif (mp_opt.eq.28) then
                 ncc2 = MAX(1.E-6, qqnw * rho)
              endif

              if (ncc2 .lt. 10.e6) then
                nu_c = 15
              else
                nu_c   = min (15, NINT(1000.e6/ncc2) + 2)
              endif

              lamc = (ncc2/rc)**obmr * (am_r*g_ratio(nu_c))**obmr

              effr = 1.0E6*MAX(4.01E-6, MIN(SNGL(1.0D0*DBLE(3.+nu_c)/lamc),50.E-6))

!           old UPP
!             effr = 2.*10.

            endif

!RAIN DIAMETER CALCULATIONS

     CASE("R")

            if( qqr > min_qr) then

              rr = MAX(1.E-12, qqr * rho)
              ncr2 = MAX(1.E-6, qqnr * rho)
              lamr = (ncr2/rr)**obmr * (am_r*crg(3)*org2)**obmr

              effr = 1.0E6*MAX(50.01E-6, MIN(SNGL(1.0D0*DBLE(3.+mu_r)/lamr),1999.E-6))

!             old UPP
!              effr=2.*200.

!              print*,'effr_rain=',effr/2.

            endif

!ICE DIAMETER CACLULATIONS

     CASE("I")

            if(qqi .ge. min_qi) then

              ri = MAX(1.E-12, qqi * rho)
              nci2 = MAX(1.E-6, qqni * rho)

              lami = (nci2/ri)**obmi * (am_i*cig(2)*oig1)**obmi

              effr = 1.0E6*MAX(10.01E-6, MIN(SNGL(1.0D0*DBLE(3.+mu_i)/lami),250.E-6))

!             old UPP
!               effr=2.*25.

            endif

!SNOW DIAMETER CALCULATIONS

     CASE("S")

            rs = MAX(1.E-12, qqs * rho)

            if(qqs .ge. min_qs) then

              tc0  = min(-0.1, t-273.15)
              smob = rs*oams

              if (nthom_bm_s.gt.(2.0-1.e-3) .and. nthom_bm_s.lt.(2.0+1.e-3))then
                  smo2 = smob
              else
                  loga = nthom_sa(1) + nthom_sa(2)*tc0 + nthom_sa(3)*nthom_bm_s+               &
                         nthom_sa(4)*tc0*nthom_bm_s + nthom_sa(5)*tc0*tc0 +&
                         nthom_sa(6)*nthom_bm_s*nthom_bm_s +nthom_sa(7)*tc0*tc0*nthom_bm_s +   &
                         nthom_sa(8)*tc0*nthom_bm_s*nthom_bm_s +nthom_sa(9)*tc0*tc0*tc0 +      &
                         nthom_sa(10)*nthom_bm_s*nthom_bm_s*nthom_bm_s

                  a    = 10.0**loga

                  b    = nthom_sb(1) + nthom_sb(2)*tc0 + nthom_sb(3)*nthom_bm_s+               &
                         nthom_sb(4)*tc0*nthom_bm_s + nthom_sb(5)*tc0*tc0 +&
                         nthom_sb(6)*nthom_bm_s*nthom_bm_s +nthom_sb(7)*tc0*tc0*nthom_bm_s +   &
                         nthom_sb(8)*tc0*nthom_bm_s*nthom_bm_s +nthom_sb(9)*tc0*tc0*tc0 +      &
                         nthom_sb(10)*nthom_bm_s*nthom_bm_s*nthom_bm_s
                  smo2 = (smob/a)**(1./b)
              endif

              !Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
              loga      = nthom_sa(1) + nthom_sa(2)*tc0 + nthom_sa(3)*cse(1) +&
                          nthom_sa(4)*tc0*cse(1) + nthom_sa(5)*tc0*tc0 +&
                          nthom_sa(6)*cse(1)*cse(1) + nthom_sa(7)*tc0*tc0*cse(1)+      &
                          nthom_sa(8)*tc0*cse(1)*cse(1) +nthom_sa(9)*tc0*tc0*tc0 +     &
                          nthom_sa(10)*cse(1)*cse(1)*cse(1)

              a       = 10.0**loga

              b       = nthom_sb(1)+ nthom_sb(2)*tc0 + nthom_sb(3)*cse(1) +&
                        nthom_sb(4)*tc0*cse(1) + nthom_sb(5)*tc0*tc0 +&
                        nthom_sb(6)*cse(1)*cse(1) + nthom_sb(7)*tc0*tc0*cse(1) +&
                        nthom_sb(8)*tc0*cse(1)*cse(1) + nthom_sb(9)*tc0*tc0*tc0 +       &
                        nthom_sb(10)*cse(1)*cse(1)*cse(1)

              smoc      = a * smo2**b

              effr = 1.0E6*MAX(50.E-6, MIN(smoc/smob, 1999.E-6))

!              print*,'snow effr=',effr

!             changing snow effr recovers "old" UPP Thompson almost exactly;
!             i.e. the snow effr is the source of the cold discprepancy.

!             old UPP
!              effr=2.*250.

            endif

     CASE("G")

            if(qqg .ge. min_qg) then

                rg2 = MAX(1.E-12, qqg * rho)

                ygra1 = alog10(max(1.E-9, rg2))

                zans1 = 3. + 2./7. * (ygra1+7.)
                zans1 = MAX(2., MIN(zans1, 7.))

                no_exp = 10.**(zans1)

                lm_exp = (no_exp*am_g*cgg(1)/rg2)**oge1

                lamg = lm_exp * (cgg(3)*ogg2*ogg1)**obmg

                effr= 1.0E6*MAX(99.E-6, MIN(SNGL((3.0+mu_g)/lamg), 9999.E-6))

!           old UPP
!            effr=350.

            endif

     END SELECT

  elseif(mp_opt.eq.11)then ! GFDL 

     SELECT CASE(species)

     CASE("C")

! cloud water (martin et al., 1994)
     if (qqw > min_qc) then
       effr = exp (1.0 / 3.0 * log ((3. * qqw ) / (4. * pi * gfdl_rhor * gfdl_ccn))) * 1.0e6
       effr = max (gfdl_rewmin, min (gfdl_rewmax, effr))
       effr = effr*2.  ! because need diameter here, converted to radius at exit
     end if

     CASE("I")

! cloud ice (heymsfield and mcfarquhar, 1996)
     if (qqi > min_qi) then
       if ((t-gfdl_tice) .lt. - 50) then
         effr = gfdl_beta / 9.917 * exp ((1 - 0.891) * log (1.0e3 * qqi)) * 1.0e3
       elseif ((t-gfdl_tice) .lt. - 40.) then
         effr = gfdl_beta / 9.337 * exp ((1 - 0.920) * log (1.0e3 * qqi)) * 1.0e3
       elseif ((t-gfdl_tice) .lt. - 30.) then
         effr = gfdl_beta / 9.208 * exp ((1 - 0.945) * log (1.0e3 * qqi)) * 1.0e3
       else
         effr = gfdl_beta / 9.387 * exp ((1 - 0.969) * log (1.0e3 * qqi)) * 1.0e3
       endif
       effr = max (gfdl_reimin, min (gfdl_reimax, effr))
       effr = effr*2.  ! because need diameter here, converted to radius at exit
     end if

     CASE("R")

     if ( qqr > min_qr ) then !rain diameter: assume gamma distribution
       lambdar = exp (0.25 * log (pi * gfdl_rhor * gfdl_n0r / qqr))
       effr  =  0.5*exp (log (gfdl_gammar / 6.) / gfdl_alphar) / lambdar * 1.0e6
       effr = max (gfdl_rermin, min (gfdl_rermax, effr))
       effr = effr*2.  ! because need diameter here, converted to radius at exit
     endif


     CASE("S")

     if ( qqs > min_qs ) then !snow diameter: assume gamma distribution
       lambdas = exp (0.25 * log (pi * gfdl_rhos * gfdl_n0s / qqs))
       effr = 0.5 * exp (log (gfdl_gammas / 6.) / gfdl_alphas) / lambdas * 1.0e6
       effr = max (gfdl_resmin, min (gfdl_resmax, effr))
       effr = effr*2.  ! because need diameter here, converted to radius at exit
     endif

     CASE("G")

     if ( qqg > min_qg ) then !graupel diameter: assume gamma distribution
       lambdag = exp (0.25 * log (pi * gfdl_rhog * gfdl_n0g / qqg))
       effr = 0.5 * exp (log (gfdl_gammag / 6.) / gfdl_alphag) / lambdag * 1.0e6
       effr = max (gfdl_regmin, min (gfdl_regmax, effr))
       effr = effr*2.  ! because need diameter here, converted to radius at exit
     endif

     END SELECT


  elseif(mp_opt.eq.5.or.mp_opt.eq.85.or.mp_opt.eq.95)then

     SELECT CASE (species)

     CASE("C")

      effr=2.*10.

     CASE("I")

      effr=2.*25.

     CASE("R")

      if( qqr > min_qr) then
      rhox=1000.
      effr=2.*1.0E6*1.5*(rho*qqr/(pi*rhox*nrain))**(1./3.)

!      old UPP    
!      effr=2.*200.
!      effr=min(200.,effr)
!      print*,'effr_rain=',effr/2.
      endif

     CASE("S")

      if(F_RimeF<=5.0)then
        RHOX=100.
          if(NLICE>0.) then
            effr  = 2.*1.0E6*1.5*(RHO*qqs/(PI*RHOX*NLICE))**(1./3.)
          endif
      endif

     CASE("G")

      if(F_RimeF>5.0.and.F_RimeF<=20.0)then
        RHOX=400.
          if(NLICE>0.) then
            effr  = 2.*1.0E6*1.5*(RHO*qqs/(PI*RHOX*NLICE))**(1./3.)
          endif
      endif

     CASE("H")

      if(F_RimeF>20.0)then
        RHOX=900.
          if(NLICE>0.) then
            effr  = 2.*1.0E6*1.5*(RHO*qqs/(PI*RHOX*NLICE))**(1./3.)
          endif
      endif

     END SELECT


  endif

!-----------------------------------------
! DIAMETER -> RADIUS
!-----------------------------------------

  effr = 0.5*effr


end function EFFR

      REAL FUNCTION GAMMLN(XX)
!     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
      IMPLICIT NONE
      REAL, INTENT(IN):: XX
      DOUBLE PRECISION, PARAMETER:: STP = 2.5066282746310005D0
      DOUBLE PRECISION, DIMENSION(6), PARAMETER:: &
               COF = (/76.18009172947146D0, -86.50532032941677D0, &
                       24.01409824083091D0, -1.231739572450155D0, &
                      .1208650973866179D-2, -.5395239384953D-5/)
      DOUBLE PRECISION:: SER,TMP,X,Y
      INTEGER:: J

      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*LOG(TMP)-TMP
      SER=1.000000000190015D0
      DO 11 J=1,6
        Y=Y+1.D0
        SER=SER+COF(J)/Y
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER/X)
      END FUNCTION GAMMLN

      REAL FUNCTION WGAMMA(y)

      IMPLICIT NONE
      REAL, INTENT(IN):: y

        real    :: GAMMLN

      WGAMMA = EXP(GAMMLN(y))

      END FUNCTION WGAMMA

