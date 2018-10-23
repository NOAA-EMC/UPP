module grib2_all_tables_module
!
!----------------module documentation block------------------------------------
!  PURPOSE: Define the variables in the derived data type discipline 
!  provided in Section 0
!
!  HISTORY: 12/04/2009 - Creation.                 V. Krishna Kumar
!                                                  NCEP Central Operations
!                                                  Systems Integration Branch
!  revision log:
!   2012/01/25   Jun Wang    Add template 4.44 and 4.48
!   2012/02/20   Jun Wang    Add complex packing
!   2014/07/08   Boi Vuong   Corrected Scaled value of second fixed surfaces 
!                            in template 4.8 and Added generating process model HRRR
!   2015/01/09   Boi Vuong   Added template 4.1, 411 and 4.12 and update code tables
!                            routines: get_g2_typeof ensfcst, get_g2_typeofderivefcst
!------------------------------------------------------------------------------
implicit none
integer, parameter :: MAXSUBCEN=25
integer, parameter :: MAXVER=30
integer, parameter :: MAXLOCVER=20
integer, parameter :: MAXREFTIME=15,MAXPRODSTATUS=17
integer, parameter :: MAXTYPEOFDATA=20
integer, parameter :: MAXTYPEOFGENPROC=30
integer, parameter :: MAXTYPEOFENSFCST=30
integer, parameter :: MAXTYPEOFDERIVEFCST=30
integer, parameter :: MAXFIXEDSURFACETYPES=200
integer, parameter :: MAXUNITOFTIMERANGE=30
integer, parameter :: MAXGENPROC=150
integer, parameter :: MAXORIGINCENTERS=500
integer, parameter :: MAXTYPEOFORIGFIELDVAL=15
integer, parameter :: MAXTYPEOFCOMPRESSION=5
integer, parameter :: MAXTYPEOFPACKINGMETHOD=15
integer, parameter :: MAXSTATPROCESSTYPES=50
integer, parameter :: MAXTYPEOFTIMEINTVLS=15
integer, parameter :: MAXTYPEOFAEROSOL=115
integer, parameter :: MAXTYPEOFINTVLS=13
integer, parameter :: MAXORDOFSPTDIFF=3

type subcenters
     character(len=20) :: subcenkey
     integer :: subcenval
end type subcenters

type(subcenters), dimension(MAXSUBCEN) :: tablec

data tablec(1) /subcenters('ncep_reanl',1)/
data tablec(2) /subcenters('ncep_ensem',2)/
data tablec(3) /subcenters('ncep_nco',3)/
data tablec(4) /subcenters('ncep_emc',4)/
data tablec(5) /subcenters('ncep_hpc',5)/
data tablec(6) /subcenters('ncep_opc',6)/
data tablec(7) /subcenters('ncep_cpc',7)/
data tablec(8) /subcenters('ncep_awc',8)/
data tablec(9) /subcenters('ncep_spc',9)/
data tablec(10) /subcenters('ncep_tpc',10)/
data tablec(11) /subcenters('nws_tdl',11)/
data tablec(12) /subcenters('nesdis_ora',12)/
data tablec(13) /subcenters('faa',13)/
data tablec(14) /subcenters('nws_mdl',14)/
data tablec(15) /subcenters('narr',15)/
data tablec(16) /subcenters('sec',16)/

type version_no
     character(len=20) :: verskey
     integer :: versval
end type version_no

type(version_no),dimension(MAXVER) :: table1_0

data table1_0(1) /version_no('expt',0)/
data table1_0(2) /version_no('v2001',1)/
data table1_0(3) /version_no('v2003',2)/
data table1_0(4) /version_no('v2005',3)/
data table1_0(5) /version_no('v2007',4)/
data table1_0(6) /version_no('v2009',5)/
data table1_0(7) /version_no('v2010',6)/
data table1_0(8) /version_no('v052011',7)/
data table1_0(9) /version_no('v112011',8)/
data table1_0(10) /version_no('v052012',9)/
data table1_0(11) /version_no('v112012',10)/
data table1_0(12) /version_no('v052013',11)/
data table1_0(13) /version_no('v112013',12)/
data table1_0(14) /version_no('v052014',13)/
data table1_0(15) /version_no('v112014',14)/
data table1_0(16) /version_no('v052015',15)/
data table1_0(17) /version_no('v112015',16)/
data table1_0(18) /version_no('v052016',17)/
data table1_0(19) /version_no('v112016',18)/
data table1_0(20) /version_no('v052017',19)/
data table1_0(21) /version_no('v112017',20)/
data table1_0(22) /version_no('v052018',21)/
data table1_0(23) /version_no('v112018',22)/
data table1_0(24) /version_no('v052019',23)/
data table1_0(25) /version_no('v112019',24)/
data table1_0(26) /version_no('v052020',25)/
data table1_0(27) /version_no('v112020',26)/
data table1_0(28) /version_no('v052021',27)/
data table1_0(29) /version_no('v112021',28)/
data table1_0(30) /version_no('preoper',29)/
!  
!
type local_table_vers_no
     character(len=20) :: locverskey
     integer :: locversval
end type local_table_vers_no

type(local_table_vers_no),dimension(MAXLOCVER) :: table1_1
!
data table1_1(1) /local_table_vers_no('local_tab_no',0)/
data table1_1(2) /local_table_vers_no('local_tab_yes1',1)/
data table1_1(3) /local_table_vers_no('local_tab_yes2',2)/
data table1_1(4) /local_table_vers_no('local_tab_yes3',3)/
data table1_1(5) /local_table_vers_no('local_tab_yes4',4)/
data table1_1(6) /local_table_vers_no('local_tab_yes5',5)/
data table1_1(7) /local_table_vers_no('missing',255)/
! 
!
type sigreftime
     character(len=20) :: sigrefkey
     integer :: sigrefval
end type sigreftime

! Declare a variable of type discipline

type(sigreftime),dimension(MAXREFTIME) :: table1_2

data table1_2(1) /sigreftime('anal',0)/
data table1_2(2) /sigreftime('fcst',1)/
data table1_2(3) /sigreftime('vfcst',2)/
data table1_2(4) /sigreftime('obstime',3)/
data table1_2(5) /sigreftime('missing',255)/
!
!
type prod_status
     character(len=20) :: prodstatuskey
     integer :: prodstatusval
end type prod_status

! Declare a variable of type discipline

type(prod_status),dimension(MAXPRODSTATUS) :: table1_3

data table1_3(1) /prod_status('oper',0)/
data table1_3(2) /prod_status('oper_test',1)/
data table1_3(3) /prod_status('res',2)/
data table1_3(4) /prod_status('reanal',3)/
data table1_3(5) /prod_status('tigge',4)/
data table1_3(6) /prod_status('tigge_test',5)/
data table1_3(7) /prod_status('missing',255)/
data table1_3(8) /prod_status('s2s_oper',6)/
data table1_3(9) /prod_status('s2s_test',7)/
!
!
type type_of_data
     character(len=20) :: typeofdatakey
     integer :: typeofdataval
end type type_of_data
!
type(type_of_data),dimension(MAXTYPEOFDATA) :: table1_4

data table1_4(1) /type_of_data('anal',0)/
data table1_4(2) /type_of_data('fcst',1)/
data table1_4(3) /type_of_data('anal_fcst',2)/
data table1_4(4) /type_of_data('con_fcst',3)/
data table1_4(5) /type_of_data('per_fcst',4)/
data table1_4(6) /type_of_data('con_per_fcst',5)/
data table1_4(7) /type_of_data('proc_sat_obs',6)/
data table1_4(8) /type_of_data('proc_rad_obs',7)/
data table1_4(9) /type_of_data('event_prob',8)/
data table1_4(10) /type_of_data('missing',255)/
!
!
type type_of_gen_proc
     character(len=30) :: typeofgenprockey
     integer :: typeofgenprocval
end type type_of_gen_proc
!
type(type_of_gen_proc),dimension(MAXTYPEOFGENPROC) :: table4_3

data table4_3(1) /type_of_gen_proc('anal',0)/
data table4_3(2) /type_of_gen_proc('init',1)/
data table4_3(3) /type_of_gen_proc('fcst',2)/
data table4_3(4) /type_of_gen_proc('bias_corr_fcst',3)/
data table4_3(5) /type_of_gen_proc('ens_fcst',4)/
data table4_3(6) /type_of_gen_proc('prob_fcst',5)/
data table4_3(7) /type_of_gen_proc('fcst_err',6)/
data table4_3(8) /type_of_gen_proc('anal_err',7)/
data table4_3(9) /type_of_gen_proc('obs',8)/
data table4_3(10) /type_of_gen_proc('clim',9)/
data table4_3(11) /type_of_gen_proc('prob_wt_fcst',10)/
data table4_3(12) /type_of_gen_proc('fcst_con_ind',192)/
data table4_3(13) /type_of_gen_proc('bias_corr_ens_fcst',11)/
data table4_3(14) /type_of_gen_proc('missing',255)/
data table4_3(15) /type_of_gen_proc('post_proc_anal',12)/
data table4_3(16) /type_of_gen_proc('post_proc_fcst',13)/
data table4_3(17) /type_of_gen_proc('nowcast',14)/
data table4_3(18) /type_of_gen_proc('hincsast',15)/
!
!
type unit_of_time_range
     character(len=30) :: unitoftimerangekey
     integer :: unitoftimerangeval
end type unit_of_time_range
!
type(unit_of_time_range),dimension(MAXUNITOFTIMERANGE) :: table4_4

data table4_4(1) /unit_of_time_range('minute',0)/
data table4_4(2) /unit_of_time_range('hour',1)/
data table4_4(3) /unit_of_time_range('day',2)/
data table4_4(4) /unit_of_time_range('month',3)/
data table4_4(5) /unit_of_time_range('year',4)/
data table4_4(6) /unit_of_time_range('decade',5)/
data table4_4(7) /unit_of_time_range('normal_30y',6)/
data table4_4(8) /unit_of_time_range('century',7)/
data table4_4(9) /unit_of_time_range('3hours',10)/
data table4_4(10) /unit_of_time_range('6hours',11)/
data table4_4(11) /unit_of_time_range('12hours',12)/
data table4_4(12) /unit_of_time_range('second',13)/
data table4_4(13) /unit_of_time_range('missing',255)/
!
!
type fixed_surface_types
     character(len=80) :: fixedsurfacetypeskey
     integer :: fixedsurfacetypesval
end type fixed_surface_types
!
type(fixed_surface_types),dimension(MAXFIXEDSURFACETYPES) :: table4_5

data table4_5(1) /fixed_surface_types('surface',1)/
data table4_5(2) /fixed_surface_types('cloud_base',2)/
data table4_5(3) /fixed_surface_types('cloud_top',3)/
data table4_5(4) /fixed_surface_types('0C_isotherm',4)/
data table4_5(5) /fixed_surface_types('lvl_of_adiab_cond_from_sfc',5)/
data table4_5(6) /fixed_surface_types('max_wind',6)/
data table4_5(7) /fixed_surface_types('tropopause',7)/
data table4_5(8) /fixed_surface_types('top_of_atmos',8)/
data table4_5(9) /fixed_surface_types('sea_bottom',9)/
data table4_5(10) /fixed_surface_types('entire_atmos',10)/
data table4_5(11) /fixed_surface_types('cb_base',11)/
data table4_5(12) /fixed_surface_types('cb_top',12)/
data table4_5(13) /fixed_surface_types('isothermal',20)/
data table4_5(14) /fixed_surface_types('isobaric_sfc',100)/
data table4_5(15) /fixed_surface_types('mean_sea_lvl',101)/
data table4_5(16) /fixed_surface_types('spec_alt_above_mean_sea_lvl',102)/
data table4_5(17) /fixed_surface_types('spec_hgt_lvl_above_grnd',103)/
data table4_5(18) /fixed_surface_types('sigma_lvl',104)/
data table4_5(19) /fixed_surface_types('hybrid_lvl',105)/
data table4_5(20) /fixed_surface_types('depth_bel_land_sfc',106)/
data table4_5(21) /fixed_surface_types('isentropic_lvl',107)/
data table4_5(22) /fixed_surface_types('spec_pres_above_grnd',108)/
data table4_5(23) /fixed_surface_types('pot_vort_sfc',109)/
data table4_5(24) /fixed_surface_types('eta_lvl',111)/
data table4_5(25) /fixed_surface_types('mixed_lyr_depth',117)/
data table4_5(26) /fixed_surface_types('depth_below_sea_lvl',160)/
data table4_5(27) /fixed_surface_types('entire_atmos_single_lyr',200)/
data table4_5(28) /fixed_surface_types('entire_ocean_single_lyr',201)/
data table4_5(29) /fixed_surface_types('hghst_trop_frz_lvl',204)/
data table4_5(30) /fixed_surface_types('grid_scale_cloud_bot_lvl',206)/
data table4_5(31) /fixed_surface_types('grid_scale_cloud_top_lvl',207)/
data table4_5(32) /fixed_surface_types('bound_lyr_cloud_bot_lvl',209)/
data table4_5(33) /fixed_surface_types('bound_lyr_cloud_top_lvl',210)/
data table4_5(34) /fixed_surface_types('bound_lyr_cloud_lyr',211)/
data table4_5(35) /fixed_surface_types('low_cloud_bot_lvl',212)/
data table4_5(36) /fixed_surface_types('low_cloud_top_lvl',213)/
data table4_5(37) /fixed_surface_types('low_cloud_lyr',214)/
data table4_5(38) /fixed_surface_types('cloud_ceilng',215)/
data table4_5(39) /fixed_surface_types('planetary_bound_lyr',220)/
data table4_5(40) /fixed_surface_types('lyr_betwn_2hybrid_lvls',221)/
data table4_5(41) /fixed_surface_types('mid_cloud_bot_lvl',222)/
data table4_5(42) /fixed_surface_types('mid_cloud_top_lvl',223)/
data table4_5(43) /fixed_surface_types('mid_cloud_lyr',224)/
data table4_5(44) /fixed_surface_types('high_cloud_bot_lvl',232)/
data table4_5(45) /fixed_surface_types('high_cloud_top_lvl',233)/
data table4_5(46) /fixed_surface_types('high_cloud_lyr',234)/
data table4_5(47) /fixed_surface_types('ocean_isotherm_lvl',235)/
data table4_5(48) /fixed_surface_types('lyr_betwn_2depths_below_ocean_sfc',236)/
data table4_5(49) /fixed_surface_types('bot_of_ocean_mix_lyr',237)/
data table4_5(50) /fixed_surface_types('bot_of_ocean_isoth_lyr',238)/
data table4_5(51) /fixed_surface_types('lyr_ocean_sfc_26c_ocean_isotherm_lvl',239)/
data table4_5(52) /fixed_surface_types('ocean_mix_lyr',240)/
data table4_5(53) /fixed_surface_types('ordered_seq_of_data',241)/
data table4_5(54) /fixed_surface_types('convective_cloud_bot_lvl',242)/
data table4_5(55) /fixed_surface_types('convective_cloud_top_lvl',243)/
data table4_5(56) /fixed_surface_types('convective_cloud_lyr',244)/
data table4_5(57) /fixed_surface_types('lwst_lvl_of_wet_bulb_zero',245)/
data table4_5(58) /fixed_surface_types('max_equiv_pot_temp_lvl',246)/
data table4_5(59) /fixed_surface_types('equil_lvl',247)/
data table4_5(60) /fixed_surface_types('shall_convective_cloud_bot_lvl',248)/
data table4_5(61) /fixed_surface_types('shall_convective_cloud_top_lvl',249)/
data table4_5(62) /fixed_surface_types('deep_convective_cloud_bot_lvl',251)/
data table4_5(63) /fixed_surface_types('deep_convective_cloud_top_lvl',252)/
data table4_5(64) /fixed_surface_types('lwst_bot_lvl_of_supercooled_liq_water_lyr',253)/
data table4_5(65) /fixed_surface_types('hghst_top_lvl_of_supercooled_liq_water_lyr',254)/
data table4_5(66) /fixed_surface_types('missing',255)/
data table4_5(67) /fixed_surface_types('hybrid_height_lvl',118)/
data table4_5(68) /fixed_surface_types('hybrid_pres_lvl',119)/
data table4_5(69) /fixed_surface_types('gen_vertical_height_coor',150)/
data table4_5(70) /fixed_surface_types('depth_below_water_lvl',161)/
data table4_5(71) /fixed_surface_types('lake_or_river_bottom',162)/
data table4_5(72) /fixed_surface_types('bottom_of_sediment_layer',163)/
data table4_5(73) /fixed_surface_types('bottom_of_therm_sediment_layer',164)/
data table4_5(74) /fixed_surface_types('bottom_of_sediment_layer_thermal_wave',165)/
data table4_5(75) /fixed_surface_types('maxing_layer',166)/
data table4_5(76) /fixed_surface_types('bottom_root_zone',167)/
!
!
type type_of_ens_fcst
     character(len=50) :: typeofensfcstkey
     integer :: typeofensfcstval
end type type_of_ens_fcst
!
type(type_of_ens_fcst),dimension(MAXTYPEOFENSFCST) :: table4_6

data table4_6(1) /type_of_ens_fcst('unpert_hi_res_ctrl_fcst',0)/
data table4_6(2) /type_of_ens_fcst('unpert_lo_res_ctrl_fcst',1)/
data table4_6(3) /type_of_ens_fcst('neg_pert_fcst',2)/
data table4_6(4) /type_of_ens_fcst('pos_pert_fcst',3)/
data table4_6(5) /type_of_ens_fcst('multi_model_fcst',4)/
data table4_6(6) /type_of_ens_fcst('pert_ens_member',192)/
!
!
type type_of_derive_fcst
     character(len=50) :: typeofderivefcstkey
     integer :: typeofderivefcstval
end type type_of_derive_fcst
!
type(type_of_derive_fcst),dimension(MAXTYPEOFDERIVEFCST) :: table4_7

data table4_7(1) /type_of_derive_fcst('unweighted_mean_all_mem',0)/
data table4_7(2) /type_of_derive_fcst('weighted_mean_all_mem',1)/
data table4_7(3) /type_of_derive_fcst('std_devn_res_cluster_mean',2)/
data table4_7(4) /type_of_derive_fcst('std_devn_res_cluster_mean_norm',3)/
data table4_7(5) /type_of_derive_fcst('spread_all_mem',4)/
data table4_7(6) /type_of_derive_fcst('large_anomaly_index',5)/
data table4_7(7) /type_of_derive_fcst('unweighted_mean_cluster',6)/
data table4_7(8) /type_of_derive_fcst('interquartile_range',7)/
data table4_7(9) /type_of_derive_fcst('min_all_ens_mem',8)/
data table4_7(10) /type_of_derive_fcst('max_all_ens_mem',9)/
data table4_7(11) /type_of_derive_fcst('unweighted_mode_all_mem',192)/
data table4_7(12) /type_of_derive_fcst('percentile_val_10',193)/
data table4_7(13) /type_of_derive_fcst('percentile_val_50',194)/
data table4_7(14) /type_of_derive_fcst('percentile_val_90',195)/
data table4_7(15) /type_of_derive_fcst('stat_decide_mem',196)/
data table4_7(16) /type_of_derive_fcst('clim_percentile',197)/
!
!
type statistical_processing_types
     character(len=80) :: statprocesstypeskey
     integer :: statprocesstypesval
end type statistical_processing_types
!
type(statistical_processing_types),dimension(MAXSTATPROCESSTYPES) :: table4_10

data table4_10(1) /statistical_processing_types('AVE',0)/
data table4_10(2) /statistical_processing_types('ACM',1)/
data table4_10(3) /statistical_processing_types('MAX',2)/
data table4_10(4) /statistical_processing_types('MIN',3)/
data table4_10(5) /statistical_processing_types('diff_end-beg',4)/
data table4_10(6) /statistical_processing_types('rms',5)/
data table4_10(7) /statistical_processing_types('std_devn',6)/
data table4_10(8) /statistical_processing_types('covariance',7)/
data table4_10(9) /statistical_processing_types('diff_beg-end',8)/
data table4_10(10) /statistical_processing_types('ratio',9)/
data table4_10(11) /statistical_processing_types('std_anomaly',10)/
data table4_10(12) /statistical_processing_types('clim_mean',192)/
data table4_10(13) /statistical_processing_types('ave_n_fcsts',193)/
data table4_10(14) /statistical_processing_types('ave_n_unin_anal',194)/
data table4_10(15) /statistical_processing_types('ave_fcst_acc_24',195)/
data table4_10(16) /statistical_processing_types('ave_succ_fcst_acc',196)/
data table4_10(17) /statistical_processing_types('ave_fcst_ave_24',197)/
data table4_10(18) /statistical_processing_types('ave_succ_fcst_ave',198)/
data table4_10(19) /statistical_processing_types('clim_ave_n_anal',199)/
data table4_10(20) /statistical_processing_types('clim_ave_n_fcst',200)/
data table4_10(21) /statistical_processing_types('clim_rms_diff',201)/
data table4_10(22) /statistical_processing_types('clim_std_n_fcst',202)/
data table4_10(23) /statistical_processing_types('clim_std_n_anal',203)/
data table4_10(24) /statistical_processing_types('ave_fcst_acc_6',204)/
data table4_10(25) /statistical_processing_types('ave_fcst_ave_6',205)/
data table4_10(26) /statistical_processing_types('ave_fcst_acc_12',206)/
data table4_10(27) /statistical_processing_types('ave_fcst_ave_12',207)/
data table4_10(28) /statistical_processing_types('missing',255)/
data table4_10(29) /statistical_processing_types('summation',11)/
!
!
type type_of_time_intervals
     character(len=80) :: typeoftimeintervalskey
     integer :: typeoftimeintervalsval
end type type_of_time_intervals
!
type(type_of_time_intervals),dimension(MAXTYPEOFTIMEINTVLS) :: table4_11

data table4_11(1) /type_of_time_intervals('reserved',0)/
data table4_11(2) /type_of_time_intervals('same_fcst_time_start_time_fcst_inc',1)/
data table4_11(3) /type_of_time_intervals('same_start_time_fcst_fcst_time_inc',2)/
data table4_11(4) /type_of_time_intervals('start_time_fcst_inc_fcst_time_dec',3)/
data table4_11(5) /type_of_time_intervals('start_time_fcst_dec_fcst_time_inc',4)/
data table4_11(6) /type_of_time_intervals('fltng_time_betwn_fcst_time_end_time_intvl',5)/
data table4_11(7) /type_of_time_intervals('local1',192)/
data table4_11(8) /type_of_time_intervals('local2',193)/
data table4_11(9) /type_of_time_intervals('local3',194)/
data table4_11(10) /type_of_time_intervals('missing',255)/
!
!
type type_of_intervals
     character(len=80) :: typeofintervalskey
     integer :: typeofintervalsval
end type type_of_intervals
!
type(type_of_intervals),dimension(MAXTYPEOFINTVLS) :: table4_91

data table4_91(1) /type_of_intervals('smaller_than_first_limit',0)/
data table4_91(2) /type_of_intervals('greater_than_second_limit',1)/
data table4_91(3) /type_of_intervals('between_first_second_limit_noincl2ndlmt',2)/
data table4_91(4) /type_of_intervals('greater_than_first_limit',3)/
data table4_91(5) /type_of_intervals('smaller_than_second_limit',4)/
data table4_91(6) /type_of_intervals('smaller_or_equal_first_limit',5)/
data table4_91(7) /type_of_intervals('greater_or_equal_second_limit',6)/
data table4_91(8) /type_of_intervals('between_first_second_limit',7)/
data table4_91(9) /type_of_intervals('greater_or_equal_first_limit',8)/
data table4_91(10) /type_of_intervals('smaller_or_equal_second_limit',9)/
data table4_91(11) /type_of_intervals('between_first_second_limit_noincl1stlmt',10)/
data table4_91(12) /type_of_intervals('equall_to_first_limit',11)/
data table4_91(13) /type_of_intervals('missing',255)/
!
!
type type_of_aerosol
     character(len=80) :: typeofaerosolkey
     integer :: typeofaerosolval
end type type_of_aerosol
!
type(type_of_aerosol),dimension(MAXTYPEOFAEROSOL) :: table4_233

data table4_233(1) /type_of_aerosol('ozone',0)/
data table4_233(2) /type_of_aerosol('water_vapor',1)/
data table4_233(3) /type_of_aerosol('methane',2)/
data table4_233(4) /type_of_aerosol('carbon_dioxide',3)/
data table4_233(5) /type_of_aerosol('carbon_monoxide',4)/
data table4_233(6) /type_of_aerosol('nitrogen_dioxide',5)/
data table4_233(7) /type_of_aerosol('nitrous_oxide',6)/
data table4_233(8) /type_of_aerosol('formaldehyde',7)/
data table4_233(9) /type_of_aerosol('sulphur_dioxide',8)/
data table4_233(10) /type_of_aerosol('ammonia',9)/
data table4_233(11) /type_of_aerosol('ammonium',10)/
data table4_233(12) /type_of_aerosol('nitrogen_monoxide',11)/
data table4_233(13) /type_of_aerosol('atomic_oxygen',12)/
data table4_233(14) /type_of_aerosol('nitrate_radical',13)/
data table4_233(15) /type_of_aerosol('hydroperoxyl_radical',14)/
data table4_233(16) /type_of_aerosol('dinitrogen_pentoxide',15)/
data table4_233(17) /type_of_aerosol('nitrous_acid',16)/
data table4_233(18) /type_of_aerosol('nitric_acid',17)/
data table4_233(19) /type_of_aerosol('peroxynitric_acid',18)/
data table4_233(20) /type_of_aerosol('hydrogen_peroxide',19)/
data table4_233(21) /type_of_aerosol('molecular_hydrogen',20)/
data table4_233(22) /type_of_aerosol('atomic_nitrogen',21)/
data table4_233(23) /type_of_aerosol('sulphate',22)/
data table4_233(24) /type_of_aerosol('radon',23)/
data table4_233(25) /type_of_aerosol('elemental_mercury',24)/
data table4_233(26) /type_of_aerosol('divalent_mercury',25)/
data table4_233(27) /type_of_aerosol('atomic_chlorine',26)/
data table4_233(28) /type_of_aerosol('chlorine_monoxide',27)/
data table4_233(29) /type_of_aerosol('dichlorine_peroxide',28)/
data table4_233(30) /type_of_aerosol('hypochlorous_acid',29)/
data table4_233(31) /type_of_aerosol('chlorine_nitrate',30)/
data table4_233(32) /type_of_aerosol('chlorine_dioxide',31)/
data table4_233(33) /type_of_aerosol('atomic_bromide',32)/
data table4_233(34) /type_of_aerosol('bromine_monoxide',33)/
data table4_233(35) /type_of_aerosol('bromine_chloride',34)/
data table4_233(36) /type_of_aerosol('hydrogen_bromide',35)/
data table4_233(37) /type_of_aerosol('hypobromous_acid',36)/
data table4_233(38) /type_of_aerosol('bromine_nitrate',37)/
data table4_233(39) /type_of_aerosol('hydroxyl_radical',10000)/
data table4_233(40) /type_of_aerosol('methyl_peroxy_radical',10001)/
data table4_233(41) /type_of_aerosol('methyl_hydroperoxide',10002)/
data table4_233(42) /type_of_aerosol('methanol',10004)/
data table4_233(43) /type_of_aerosol('formic_acid',10005)/
data table4_233(44) /type_of_aerosol('hydrogen_cyanide',10006)/
data table4_233(45) /type_of_aerosol('aceto_nitrile',10007)/
data table4_233(46) /type_of_aerosol('ethane',10008)/
data table4_233(47) /type_of_aerosol('ethene',10009)/
data table4_233(48) /type_of_aerosol('ethyne',10010)/
data table4_233(49) /type_of_aerosol('ethanol',10011)/
data table4_233(50) /type_of_aerosol('acetic_acid',10012)/
data table4_233(51) /type_of_aerosol('peroxyacetyl_nitrate',10013)/
data table4_233(52) /type_of_aerosol('propane',10014)/
data table4_233(53) /type_of_aerosol('propene',10015)/
data table4_233(54) /type_of_aerosol('butanes',10016)/
data table4_233(55) /type_of_aerosol('isoprene',10017)/
data table4_233(56) /type_of_aerosol('alpha_pinene',10018)/
data table4_233(57) /type_of_aerosol('beta_pinene',10019)/
data table4_233(58) /type_of_aerosol('limonene',10020)/
data table4_233(59) /type_of_aerosol('benzene',10021)/
data table4_233(60) /type_of_aerosol('toluene',10022)/
data table4_233(61) /type_of_aerosol('xylene',10023)/
data table4_233(62) /type_of_aerosol('dumethyl_sulphide',10500)/
data table4_233(63) /type_of_aerosol('hydrogen_chloride',20001)/
data table4_233(64) /type_of_aerosol('cfc-11',20002)/
data table4_233(65) /type_of_aerosol('cfc-12',20003)/
data table4_233(66) /type_of_aerosol('cfc-113',20004)/
data table4_233(67) /type_of_aerosol('cfc-113a',20005)/
data table4_233(68) /type_of_aerosol('cfc-114',20006)/
data table4_233(69) /type_of_aerosol('cfc-115',20007)/
data table4_233(70) /type_of_aerosol('hcfc-22',20008)/
data table4_233(71) /type_of_aerosol('hcfc-141b',20009)/
data table4_233(72) /type_of_aerosol('hcfc-142b',20010)/
data table4_233(73) /type_of_aerosol('halon-1202',20011)/
data table4_233(74) /type_of_aerosol('halon-1211',20012)/
data table4_233(75) /type_of_aerosol('halon-1301',20013)/
data table4_233(76) /type_of_aerosol('halon-2402',20014)/
data table4_233(77) /type_of_aerosol('methyl_chloride',20015)/
data table4_233(78) /type_of_aerosol('carbon_tetrachloride',20016)/
data table4_233(79) /type_of_aerosol('hcc-140a',20017)/
data table4_233(80) /type_of_aerosol('methyl_bromide',20018)/
data table4_233(81) /type_of_aerosol('hexachlorocyclohexane',20019)/
data table4_233(82) /type_of_aerosol('alpha_hexachlorocyclohexane',20020)/
data table4_233(83) /type_of_aerosol('hexachlorobiphenyl',20021)/
data table4_233(84) /type_of_aerosol('radioactive_pollutant',30000)/
data table4_233(85) /type_of_aerosol('hox_radical',60000)/
data table4_233(86) /type_of_aerosol('total_inorganic_org_peroxy_radicals',60001)/
data table4_233(87) /type_of_aerosol('passive_ozone',60002)/
data table4_233(88) /type_of_aerosol('nox_nitrogen',60003)/
data table4_233(89) /type_of_aerosol('all_nitrogen_oxides',60004)/
data table4_233(90) /type_of_aerosol('total_inorganic_chlorine',60005)/
data table4_233(91) /type_of_aerosol('total_inorganic_bromine',60006)/
data table4_233(92) /type_of_aerosol('total_inorganic_chlorine_noHclClono2Clox',60007)/
data table4_233(93) /type_of_aerosol('total_inorganic_bromine_noHbrBrono2Brox',60008)/
data table4_233(94) /type_of_aerosol('lumped_alkanes',60009)/
data table4_233(95) /type_of_aerosol('lumped_alkenes',60010)/
data table4_233(96) /type_of_aerosol('lumped_aromatic_comp',60011)/
data table4_233(97) /type_of_aerosol('lumped_terpenes',60012)/
data table4_233(98) /type_of_aerosol('non_methane_volatile_org_comp_carbon)',60013)/
data table4_233(99) /type_of_aerosol('anthropogenic_non_methane_voiatile_org_comp_carbon',60014)/
data table4_233(100) /type_of_aerosol('biogenic_non_methane_volatile_org_comp_carbon',60015)/
data table4_233(101) /type_of_aerosol('lumped_oxygenated_hydrocarbon',60016)/
data table4_233(102) /type_of_aerosol('total_aerosol',62000)/
data table4_233(103) /type_of_aerosol('dust_dry',62001)/
data table4_233(104) /type_of_aerosol('water_in_ambient',62002)/
data table4_233(105) /type_of_aerosol('ammonium_dry',62003)/
data table4_233(106) /type_of_aerosol('nitrate_dry',62004)/
data table4_233(107) /type_of_aerosol('nitric_acid_trihydrate',62005)/
data table4_233(108) /type_of_aerosol('sulphate_dry',62006)/
data table4_233(109) /type_of_aerosol('mercury_dry',62007)/
data table4_233(110) /type_of_aerosol('sea_salt_dry',62008)/
data table4_233(111) /type_of_aerosol('black_carbon_dry',62009)/
data table4_233(112) /type_of_aerosol('particulate_org_matter_dry',62010)/
data table4_233(113) /type_of_aerosol('primary_particulate_org_matter_dry',62011)/
data table4_233(114) /type_of_aerosol('secondary_particulate_org_matter_dry',62012)/
data table4_233(115) /type_of_aerosol('missing',65535)/
!
!
type type_of_orig_field_vals
     character(len=50) :: typeoforigfieldvalskey
     integer :: typeoforigfieldvals
end type type_of_orig_field_vals
!
type(type_of_orig_field_vals), dimension(MAXTYPEOFORIGFIELDVAL) ::table5_1
     data table5_1(1) /type_of_orig_field_vals('fltng_pnt',0)/
     data table5_1(2) /type_of_orig_field_vals('integer',1)/
     data table5_1(3) /type_of_orig_field_vals('local1',192)/
     data table5_1(4) /type_of_orig_field_vals('local2',193)/
     data table5_1(5) /type_of_orig_field_vals('local3',194)/
     data table5_1(6) /type_of_orig_field_vals('local4',195)/
     data table5_1(7) /type_of_orig_field_vals('local5',196)/
     data table5_1(8) /type_of_orig_field_vals('local6',197)/
     data table5_1(9) /type_of_orig_field_vals('local7',198)/
     data table5_1(10) /type_of_orig_field_vals('local8',199)/
     data table5_1(11) /type_of_orig_field_vals('local9',200)/
     data table5_1(12) /type_of_orig_field_vals('local10',201)/
     data table5_1(13) /type_of_orig_field_vals('missing',255)/
!
!
type order_of_sptdiff_vals
     character(len=50) :: ordofsptdiffkey
     integer :: ordofsptdiffvals
end type order_of_sptdiff_vals
!
type(order_of_sptdiff_vals), dimension(MAXORDOFSPTDIFF) :: table5_6
     data table5_6(1) /order_of_sptdiff_vals('1st_ord_sptdiff',1)/
     data table5_6(2) /order_of_sptdiff_vals('2nd_ord_sptdiff',2)/
     data table5_6(3) /order_of_sptdiff_vals('missing',255)/
!
!
type type_of_compression
     character(len=50) :: typeofcompressionkey
     integer :: typeofcompressionvals
end type type_of_compression
!
type(type_of_compression), dimension(MAXTYPEOFCOMPRESSION) :: table5_40
     data table5_40(1) /type_of_compression('lossless',0)/
     data table5_40(2) /type_of_compression('lossy',1)/
     data table5_40(3) /type_of_compression('missing',255)/
!
!
type type_of_packingmethod
     character(len=50) :: packingmethodkey
     integer :: packingmethodvals
end type type_of_packingmethod
!
type(type_of_packingmethod), dimension(MAXTYPEOFPACKINGMETHOD) :: table5_0
     data table5_0(1) /type_of_packingmethod('simple_packing',0)/
     data table5_0(2) /type_of_packingmethod('maxtric_simple_packing',1)/
     data table5_0(3) /type_of_packingmethod('complex_packing',2)/
     data table5_0(4) /type_of_packingmethod('complex_packing_spatial_diff',3)/
     data table5_0(5) /type_of_packingmethod('ieee_floating_point',4)/
     data table5_0(6) /type_of_packingmethod('jpeg2000',40)/
     data table5_0(7) /type_of_packingmethod('png',41)/
     data table5_0(8) /type_of_packingmethod('spectral_simple_packing',50)/
     data table5_0(9) /type_of_packingmethod('spectral_complex_packing',51)/
     data table5_0(10) /type_of_packingmethod('simple_packing_log_preprcs',61)/
     data table5_0(11) /type_of_packingmethod('run_length_packing_lvl_val',200)/
!
!
type origin_centers
     character(len=50) :: origincenterskey
     integer :: origincentersval
end type origin_centers
!
type(origin_centers),dimension(MAXORIGINCENTERS) :: on388_table0

        data on388_table0(1) /origin_centers('melbourne1',1)/
        data on388_table0(2) /origin_centers('melbourne2',2)/
        data on388_table0(3) /origin_centers('melbourne3',3)/
        data on388_table0(4) /origin_centers('moscow1',4)/
        data on388_table0(5) /origin_centers('moscow2',5)/
        data on388_table0(6) /origin_centers('moscow3',6)/
        data on388_table0(7) /origin_centers('nws_ncep',7)/
        data on388_table0(8) /origin_centers('nws_nwstg',8)/
        data on388_table0(9) /origin_centers('nws_other',9)/
        data on388_table0(10) /origin_centers('cairo1',10)/
        data on388_table0(11) /origin_centers('cairo2',11)/
        data on388_table0(12) /origin_centers('dakar1',12)/
        data on388_table0(13) /origin_centers('dakar2',13)/
        data on388_table0(14) /origin_centers('nairobi1',14)/
        data on388_table0(15) /origin_centers('nairobi2',15)/
        data on388_table0(16) /origin_centers('casablanca',16)/
        data on388_table0(17) /origin_centers('tunis',17)/
        data on388_table0(18) /origin_centers('tunis_casablanca1',18)/
        data on388_table0(19) /origin_centers('tunis-casablanca2',19)/
        data on388_table0(20) /origin_centers('las_palmas',20)/
        data on388_table0(21) /origin_centers('algiers',21)/
        data on388_table0(22) /origin_centers('acmad',22)/
        data on388_table0(23) /origin_centers('mozambique',23)/
        data on388_table0(24) /origin_centers('pretoria',24)/
        data on388_table0(25) /origin_centers('la_reunion',25)/
        data on388_table0(26) /origin_centers('khabarovsk1',26)/
        data on388_table0(27) /origin_centers('khabarovsk2',27)/
        data on388_table0(28) /origin_centers('new_delhi1',28)/
        data on388_table0(29) /origin_centers('new_delhi2',29)/
        data on388_table0(30) /origin_centers('novosibirsk1',30)/
        data on388_table0(31) /origin_centers('novosibirsk2',31)/
        data on388_table0(32) /origin_centers('tashkent',32)/
        data on388_table0(33) /origin_centers('jeddah',33)/
        data on388_table0(34) /origin_centers('jma_tokyo1',34)/
        data on388_table0(35) /origin_centers('jma_tokyo2',35)/
        data on388_table0(36) /origin_centers('bankok',36)/
        data on388_table0(37) /origin_centers('ulan_bator',37)/
        data on388_table0(38) /origin_centers('beijing1',38)/
        data on388_table0(39) /origin_centers('beijing2',39)/
        data on388_table0(40) /origin_centers('seoul',40)/
        data on388_table0(41) /origin_centers('buenos_aires1',41)/
        data on388_table0(42) /origin_centers('buenos_aires2',42)/
        data on388_table0(43) /origin_centers('brasilia1',43)/
        data on388_table0(44) /origin_centers('brasilia2',44)/
        data on388_table0(45) /origin_centers('santiago',45)/
        data on388_table0(46) /origin_centers('brazilian_inpe',46)/
        data on388_table0(47) /origin_centers('columbia',47)/
        data on388_table0(48) /origin_centers('ecuador',48)/
        data on388_table0(49) /origin_centers('peru',49)/
        data on388_table0(50) /origin_centers('venezuela',50)/
        data on388_table0(51) /origin_centers('miami',51)/
        data on388_table0(52) /origin_centers('tpc_nhc',52)/
        data on388_table0(53) /origin_centers('cms_montreal1',53)/
        data on388_table0(54) /origin_centers('cms_montreal2',54)/
        data on388_table0(55) /origin_centers('san_francisco',55)/
        data on388_table0(56) /origin_centers('arinc_center',56)/
        data on388_table0(57) /origin_centers('usaf_gwc',57)/
        data on388_table0(58) /origin_centers('us_navy_fnoc',58)/
        data on388_table0(59) /origin_centers('noaa_fsl_boulder',59)/
        data on388_table0(60) /origin_centers('ncar_boulder',60)/
        data on388_table0(61) /origin_centers('service_argos',61)/
        data on388_table0(62) /origin_centers('us_naval_ocean_off',62)/
        data on388_table0(63) /origin_centers('honolulu',64)/
        data on388_table0(64) /origin_centers('darwin1',65)/
        data on388_table0(65) /origin_centers('darwin2',66)/
        data on388_table0(66) /origin_centers('melbourne4',67)/
        data on388_table0(67) /origin_centers('wellington1',69)/
        data on388_table0(68) /origin_centers('wellington2',70)/
        data on388_table0(69) /origin_centers('nadi',71)/
        data on388_table0(70) /origin_centers('singapore',72)/
        data on388_table0(71) /origin_centers('malaysia',73)/
        data on388_table0(72) /origin_centers('ukmo_exeter1',74)/
        data on388_table0(73) /origin_centers('ukmo_exeter2',75)/
        data on388_table0(74) /origin_centers('moscow4',76)/
        data on388_table0(75) /origin_centers('offenbach1',78)/
        data on388_table0(76) /origin_centers('offenbach2',79)/
        data on388_table0(77) /origin_centers('rome1',80)/
        data on388_table0(78) /origin_centers('rome2',81)/
        data on388_table0(79) /origin_centers('norrkoping1',82)/
        data on388_table0(80) /origin_centers('norrkoping2',83)/
        data on388_table0(81) /origin_centers('french_weather_toulouse1',84)/
        data on388_table0(82) /origin_centers('french_weather_toulouse2',85)/
        data on388_table0(83) /origin_centers('helsinki',86)/
        data on388_table0(84) /origin_centers('belgrade',87)/
        data on388_table0(85) /origin_centers('oslo',88)/
        data on388_table0(86) /origin_centers('prague',89)/
        data on388_table0(87) /origin_centers('episkopi',90)/
        data on388_table0(88) /origin_centers('ankara',91)/
        data on388_table0(89) /origin_centers('frankfurt_main',92)/
        data on388_table0(90) /origin_centers('london',93)/
        data on388_table0(91) /origin_centers('copenhagen',94)/
        data on388_table0(92) /origin_centers('rota',95)/
        data on388_table0(93) /origin_centers('athens',96)/
        data on388_table0(94) /origin_centers('esa',97)/
        data on388_table0(95) /origin_centers('ecmwf',98)/
        data on388_table0(96) /origin_centers('de_bilt_netherlands',99)/
        data on388_table0(97) /origin_centers('brazzaville',100)/
        data on388_table0(98) /origin_centers('abidjan',101)/
        data on388_table0(99) /origin_centers('libyan_arab_jamahiriya',102)/
        data on388_table0(100) /origin_centers('madagascar',103)/
        data on388_table0(101) /origin_centers('mauritius',104)/
        data on388_table0(102) /origin_centers('niger',105)/
        data on388_table0(103) /origin_centers('seychelles',106)/
        data on388_table0(104) /origin_centers('uganda',107)/
        data on388_table0(105) /origin_centers('tanzania',108)/
        data on388_table0(106) /origin_centers('zimbabwe',109)/
        data on388_table0(107) /origin_centers('hong_kong',110)/
        data on388_table0(108) /origin_centers('afghanistan',111)/
        data on388_table0(109) /origin_centers('bahrain',112)/
        data on388_table0(110) /origin_centers('bangladesh',113)/
        data on388_table0(111) /origin_centers('bhutan',114)/
        data on388_table0(112) /origin_centers('cambodia',115)/
        data on388_table0(113) /origin_centers('dprk',116)/
        data on388_table0(114) /origin_centers('iran',117)/
        data on388_table0(115) /origin_centers('iraq',118)/
        data on388_table0(116) /origin_centers('kazakhstan',119)/
        data on388_table0(117) /origin_centers('kuwait',120)/
        data on388_table0(118) /origin_centers('kyrgyz_republic',121)/
        data on388_table0(119) /origin_centers('lao_pdr',122)/
        data on388_table0(120) /origin_centers('macao_china',123)/
        data on388_table0(121) /origin_centers('maldives',124)/
        data on388_table0(122) /origin_centers('myanmar',125)/
        data on388_table0(123) /origin_centers('nepal',126)/
        data on388_table0(124) /origin_centers('oman',127)/
        data on388_table0(125) /origin_centers('pakistan',128)/
        data on388_table0(126) /origin_centers('qatar',129)/
        data on388_table0(127) /origin_centers('yemen',130)/
        data on388_table0(128) /origin_centers('sri_lanka',131)/
        data on388_table0(129) /origin_centers('tajikistan',132)/
        data on388_table0(130) /origin_centers('turkmenistan',133)/
        data on388_table0(131) /origin_centers('uae',134)/
        data on388_table0(132) /origin_centers('uzbekistan',135)/
        data on388_table0(133) /origin_centers('viet_nam ',136)/
        data on388_table0(134) /origin_centers('bolivia',140)/
        data on388_table0(135) /origin_centers('guyana',141)/
        data on388_table0(136) /origin_centers('paraguay',142)/
        data on388_table0(137) /origin_centers('suriname',143)/
        data on388_table0(138) /origin_centers('uruguay',144)/
        data on388_table0(139) /origin_centers('french_guyana',145)/
        data on388_table0(140) /origin_centers('brazilian_navy_hydro_center',146)/
        data on388_table0(141) /origin_centers('antigua_barbuda',150)/
        data on388_table0(142) /origin_centers('bahamas',151)/
        data on388_table0(143) /origin_centers('barbados',152)/
        data on388_table0(144) /origin_centers('belize',153)/
        data on388_table0(145) /origin_centers('british_caribbean_terr_center',154)/
        data on388_table0(146) /origin_centers('san_jose',155)/
        data on388_table0(147) /origin_centers('cuba',156)/
        data on388_table0(148) /origin_centers('dominica',157)/
        data on388_table0(149) /origin_centers('dominican_republic',158)/
        data on388_table0(150) /origin_centers('el_salvador',159)/
        data on388_table0(151) /origin_centers('us_noaa_nesdis',160)/
        data on388_table0(152) /origin_centers('us_noaa_oar',161)/
        data on388_table0(153) /origin_centers('guatemala',162)/
        data on388_table0(154) /origin_centers('haiti',163)/
        data on388_table0(155) /origin_centers('honduras',164)/
        data on388_table0(156) /origin_centers('jamaica',165)/
        data on388_table0(157) /origin_centers('mexico',166)/
        data on388_table0(158) /origin_centers('netherlands_antilles_aruba',167)/
        data on388_table0(159) /origin_centers('nicaragua',168)/
        data on388_table0(160) /origin_centers('panama',169)/
        data on388_table0(161) /origin_centers('saint_lucia',170)/
        data on388_table0(162) /origin_centers('trinidad_tobago',171)/
        data on388_table0(163) /origin_centers('french_departments',172)/
        data on388_table0(164) /origin_centers('cook_islands',190)/
        data on388_table0(165) /origin_centers('french_polynesia',191)/
        data on388_table0(166) /origin_centers('tonga',192)/
        data on388_table0(167) /origin_centers('vanuatu',193)/
        data on388_table0(168) /origin_centers('brunei',194)/
        data on388_table0(169) /origin_centers('indonesia',195)/
        data on388_table0(170) /origin_centers('kiribati',196)/
        data on388_table0(171) /origin_centers('federated_states_micronesia',197)/
        data on388_table0(172) /origin_centers('new_caledonia',198)/
        data on388_table0(173) /origin_centers('niue',199)/
        data on388_table0(174) /origin_centers('papua_new_guinea',200)/
        data on388_table0(175) /origin_centers('philippines',201)/
        data on388_table0(176) /origin_centers('samoa',202)/
        data on388_table0(177) /origin_centers('solomon_islands',203)/
        data on388_table0(178) /origin_centers('frascati',210)/
        data on388_table0(179) /origin_centers('lanion',211)/
        data on388_table0(180) /origin_centers('lisboa',212)/
        data on388_table0(181) /origin_centers('reykjavik',213)/
        data on388_table0(182) /origin_centers('madrid',214)/
        data on388_table0(183) /origin_centers('zurich',215)/
        data on388_table0(184) /origin_centers('service_argos_toulouse_fr',216)/
        data on388_table0(185) /origin_centers('bratislava',217)/
        data on388_table0(186) /origin_centers('budapest',218)/
        data on388_table0(187) /origin_centers('ljubljana',219)/
        data on388_table0(188) /origin_centers('warsaw',220)/
        data on388_table0(189) /origin_centers('zagreb',221)/
        data on388_table0(190) /origin_centers('albania',222)/
        data on388_table0(191) /origin_centers('armenia',223)/
        data on388_table0(192) /origin_centers('austria',224)/
        data on388_table0(193) /origin_centers('azerbaijan',225)/
        data on388_table0(194) /origin_centers('belarus',226)/
        data on388_table0(195) /origin_centers('belgium',227)/
        data on388_table0(196) /origin_centers('bosnia_herzegovina',228)/
        data on388_table0(197) /origin_centers('bulgaria',229)/
        data on388_table0(198) /origin_centers('cyprus',230)/
        data on388_table0(199) /origin_centers('estonia',231)/
        data on388_table0(200) /origin_centers('georgia',232)/
        data on388_table0(201) /origin_centers('dublin',233)/
        data on388_table0(202) /origin_centers('israel',234)/
        data on388_table0(203) /origin_centers('jordan',235)/
        data on388_table0(204) /origin_centers('latvia',236)/
        data on388_table0(205) /origin_centers('lebanon',237)/
        data on388_table0(206) /origin_centers('lithuania',238)/
        data on388_table0(207) /origin_centers('luxembourg',239)/
        data on388_table0(208) /origin_centers('malta',240)/
        data on388_table0(209) /origin_centers('monaco',241)/
        data on388_table0(210) /origin_centers('romania',242)/
        data on388_table0(211) /origin_centers('syrian_arab_republic',243)/
        data on388_table0(212) /origin_centers('macedonia',244)/
        data on388_table0(213) /origin_centers('ukraine',245)/
        data on388_table0(214) /origin_centers('republic_moldova',246)/
        data on388_table0(215) /origin_centers('eumetsat_op_cen',254)/
        data on388_table0(216) /origin_centers('missing',255)/
!
!
type gen_proc
     character(len=30) :: genprockey
     integer :: genprocval
end type gen_proc
!
type(gen_proc),dimension(MAXGENPROC) :: on388_tablea

data on388_tablea(1) /gen_proc('res',0)/
data on388_tablea(2) /gen_proc('uvim',2)/
data on388_tablea(3) /gen_proc('ncep_arl_tdm',3)/
data on388_tablea(4) /gen_proc('ncep_arl_smoke',4)/
data on388_tablea(5) /gen_proc('sat_der_prec_temp',5)/
data on388_tablea(6) /gen_proc('gwind_wave_mod',10)/
data on388_tablea(7) /gen_proc('multi_grid_wave_mod',11)/
data on388_tablea(8) /gen_proc('prob_st_surg',12)/
data on388_tablea(9) /gen_proc('lfm_anal',19)/
data on388_tablea(10) /gen_proc('snow_cov_anal',25)/
data on388_tablea(11) /gen_proc('for_gen_field',30)/
data on388_tablea(12) /gen_proc('val_add_post_proc_field',31)/
data on388_tablea(13) /gen_proc('ngm',39)/
data on388_tablea(14) /gen_proc('goi_gfs',42)/
data on388_tablea(15) /gen_proc('goi_fnl',43)/
data on388_tablea(16) /gen_proc('ssta',44)/
data on388_tablea(17) /gen_proc('coast_ocm',45)/
data on388_tablea(18) /gen_proc('hycom_glob',46)/
data on388_tablea(19) /gen_proc('hycom_npb',47)/
data on388_tablea(20) /gen_proc('hycom_nab',48)/
data on388_tablea(21) /gen_proc('ozone_anal_tiros',49)/
data on388_tablea(22) /gen_proc('ozone_anal_nimbus',52)/
data on388_tablea(23) /gen_proc('lfm_fofm',53)/
data on388_tablea(24) /gen_proc('roi',64)/
data on388_tablea(25) /gen_proc('t80l18gfs',68)/
data on388_tablea(26) /gen_proc('t80l18mrf',69)/
data on388_tablea(27) /gen_proc('qlm',70)/
data on388_tablea(28) /gen_proc('fogfm_opc',73)/
data on388_tablea(29) /gen_proc('gulf_of_mex_wind_wave',74)/
data on388_tablea(30) /gen_proc('gulf_of_alas_wind_wave',75)/
data on388_tablea(31) /gen_proc('bias_corr_mrf',76)/
data on388_tablea(32) /gen_proc('t126l28gfs',77)/
data on388_tablea(33) /gen_proc('t126l28mrf',78)/
data on388_tablea(34) /gen_proc('backup_from_prev_run',79)/
data on388_tablea(35) /gen_proc('t62l28mrf',80)/
data on388_tablea(36) /gen_proc('anal_gfs',81)/
data on388_tablea(37) /gen_proc('anal_gdas',82)/
data on388_tablea(38) /gen_proc('meso_nam12km',84)/
data on388_tablea(39) /gen_proc('ruc_fsl_isen_60km_40n',86)/
data on388_tablea(40) /gen_proc('cac_ensem_fcsts_spect',87)/
data on388_tablea(41) /gen_proc('nww3_owm',88)/
data on388_tablea(42) /gen_proc('nmm_8km',89)/
data on388_tablea(43) /gen_proc('t62l28extmrf',90)/
data on388_tablea(44) /gen_proc('t62l28extgfs',91)/
data on388_tablea(45) /gen_proc('t62l28mrffnl',92)/
data on388_tablea(46) /gen_proc('t62l28gdasmrf',93)/
data on388_tablea(47) /gen_proc('t170l42mrf',94)/
data on388_tablea(48) /gen_proc('t126l42mrf',95)/
data on388_tablea(49) /gen_proc('gfs_avn',96)/
data on388_tablea(50) /gen_proc('cfs_t62l64_l40mom3',98)/
data on388_tablea(51) /gen_proc('misc_test_id',99)/
data on388_tablea(52) /gen_proc('ruc_sanal_60km_40n',100)/
data on388_tablea(53) /gen_proc('ruc_sanal_40km_40n',101)/
data on388_tablea(54) /gen_proc('ruc_fsl_isen_20km_40n',105)/
data on388_tablea(55) /gen_proc('gefs',107)/
data on388_tablea(56) /gen_proc('lamp',108)/
data on388_tablea(57) /gen_proc('rtma',109)/
data on388_tablea(58) /gen_proc('nam_15k',110)/
data on388_tablea(59) /gen_proc('nam_gen_sref',111)/
data on388_tablea(60) /gen_proc('wrf_nmm_ncep',112)/
data on388_tablea(61) /gen_proc('prod_ncep_sref',113)/
data on388_tablea(62) /gen_proc('naefs_prod_ncep_cmc',114)/
data on388_tablea(63) /gen_proc('down_scal_gfs_nam_ext',115)/
data on388_tablea(64) /gen_proc('wrf_em_ncar_arwrf',116)/
data on388_tablea(65) /gen_proc('ice_conc_anal',120)/
data on388_tablea(66) /gen_proc('wna_reg_wav_mod',121)/
data on388_tablea(67) /gen_proc('alas_wat_reg_wav_mod',122)/
data on388_tablea(68) /gen_proc('na_hurr_wav_mod',123)/
data on388_tablea(69) /gen_proc('enp_reg_wav_mod',124)/
data on388_tablea(70) /gen_proc('np_hurr_wav_mod',125)/
data on388_tablea(71) /gen_proc('sea_ice_fcst_mod',126)/
data on388_tablea(72) /gen_proc('lake_ice_fcst_mod',127)/
data on388_tablea(73) /gen_proc('glob_oce_fcst_mod',128)/
data on388_tablea(74) /gen_proc('godas',129)/
data on388_tablea(75) /gen_proc('merge_fields_ruc_nam_gfs',130)/
data on388_tablea(76) /gen_proc('great_lakes_wave_mod',131)/
data on388_tablea(77) /gen_proc('narr',140)/
data on388_tablea(78) /gen_proc('ldafs',141)/
data on388_tablea(79) /gen_proc('nwsrfs',150)/
data on388_tablea(80) /gen_proc('nwsffgs',151)/
data on388_tablea(81) /gen_proc('wsr_88d_s2_prec_anal',152)/
data on388_tablea(82) /gen_proc('wsr_88d_s3_prec_anal',153)/
data on388_tablea(83) /gen_proc('qpf_ncep',180)/
data on388_tablea(84) /gen_proc('rfcqpf_ncep',181)/
data on388_tablea(85) /gen_proc('rfcqpe_ncep',182)/
data on388_tablea(86) /gen_proc('ndfd_ncep_hpc',183)/
data on388_tablea(87) /gen_proc('ncwd_ncep_awc',190)/
data on388_tablea(88) /gen_proc('cipap_ncep_awc',191)/
data on388_tablea(89) /gen_proc('anal_ncep_awc',192)/
data on388_tablea(90) /gen_proc('fcst_ncep_awc',193)/
data on388_tablea(91) /gen_proc('cdas2',195)/
data on388_tablea(92) /gen_proc('cdas2_regen',196)/
data on388_tablea(93) /gen_proc('cdas',197)/
data on388_tablea(94) /gen_proc('cdas_regen',198)/
data on388_tablea(95) /gen_proc('cfsr_t382l64_l40mom4',199)/
data on388_tablea(96) /gen_proc('cpc_man_fcst',200)/
data on388_tablea(97) /gen_proc('cpc_auto_prod',201)/
data on388_tablea(98) /gen_proc('epa_usne',210)/
data on388_tablea(99) /gen_proc('epa_use',211)/
data on388_tablea(100) /gen_proc('spc_man_fcst',215)/
data on388_tablea(101) /gen_proc('ncep_opc_auto_prod',220)/
data on388_tablea(102) /gen_proc('missing',255)/
data on388_tablea(103) /gen_proc('ngac',117)/
data on388_tablea(104) /gen_proc('hrrr',83)/
data on388_tablea(105) /gen_proc('ncep_arl_dust',6)/
data on388_tablea(106) /gen_proc('hrricane_mult_wave',13)/
data on388_tablea(107) /gen_proc('extratropical_storm_surge',14)/
data on388_tablea(108) /gen_proc('nearshore_wave_prediction',15)/

contains
!
!
     subroutine get_g2_subcenters(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_subcenters
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 subcenters 
!   value for a given short key name based on Table C 
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_subcenters(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for subcenter 
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 subcenter value from table c 
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXSUBCEN
        if (trim(tablec(n)%subcenkey).eq.trim(key)) then
            value=tablec(n)%subcenval
            return
        endif
     enddo
           print *,'get_g2_subcenters key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_subcenters
!
!
     subroutine get_g2_versionno(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_versionno
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 version
!   number for a given short key name based on Table 1.0
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_versionno(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for version number
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 version number value from table 1.0
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXVER
        if (trim(table1_0(n)%verskey).eq.trim(key)) then
            value=table1_0(n)%versval
            return
        endif
     enddo
           print *,'get_g2_versionno key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_versionno
!
!
     subroutine get_g2_loctabversno(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_loctabversno
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 local table version
!   number for a given short key name based on Table 1.1
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_loctabversno(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for local table version number
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 local table version number value from table 1.1
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXLOCVER
        if (trim(table1_1(n)%locverskey).eq.trim(key)) then
            value=table1_1(n)%locversval
            return
        endif
     enddo
           print *,'get_g2_loctabversno key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_loctabversno
!
!
     subroutine get_g2_sigreftime(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_sigreftime
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 significant
!   reference time value for a given short key name based on Table 1.2
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_sigreftime(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for significant reference time 
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 significant value from table 1.2 
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
!     integer, parameter :: MAXREFTIME=15
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXREFTIME
        if (trim(table1_2(n)%sigrefkey).eq.trim(key)) then
            value=table1_2(n)%sigrefval
            return
        endif
     enddo
           print *,'get_g2_sigreftime key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_sigreftime
!
!
     subroutine get_g2_prodstatus(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_prodstatus
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 production 
!   status of data value for a given short key name based on Table 1.3
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_prodstatus(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for production status of data
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 significant value from table 1.3
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXPRODSTATUS
        if (trim(table1_3(n)%prodstatuskey).eq.trim(key)) then
            value=table1_3(n)%prodstatusval
            return
        endif
     enddo
           print *,'get_g2_prodstatus key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_prodstatus
!
!
     subroutine get_g2_typeofdata(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_typeofdata
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 type of
!   data value for a given short key name based on Table 1.4
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_typeofdata(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for production status of data
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 type of data value from table 1.4
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXTYPEOFDATA
        if (trim(table1_4(n)%typeofdatakey).eq.trim(key)) then
            value=table1_4(n)%typeofdataval
            return
        endif
     enddo
           print *,'get_g2_typeofdata key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_typeofdata
!
!
     subroutine get_g2_typeofgenproc(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_typeofgenproc
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 - Type of Generating
!   Process value for a given short key name based on Table 4.3 of Section 4, Octet 12
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_typeofgenproc(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for type of generating process from Table 4.3 
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 value from Table 4.3 
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXTYPEOFGENPROC
        if (trim(table4_3(n)%typeofgenprockey).eq.trim(key)) then
            value=table4_3(n)%typeofgenprocval
            return
        endif
     enddo
           print *,'get_g2_typeofgenproc key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_typeofgenproc
!
!
     subroutine get_g2_unitoftimerange(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_unitoftimerange
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 - Indicator of unit of time 
!   range value for a given short key name based on Table 4.4 of Section 4, Octet 18
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_unitoftimerange(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key  - GRIB2 character short key for indicator of unit of time range from Table 4.4
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 value from Table 4.4
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXUNITOFTIMERANGE
        if (trim(table4_4(n)%unitoftimerangekey).eq.trim(key)) then
            value=table4_4(n)%unitoftimerangeval
            return
        endif
     enddo

     value=255
           print *,'get_g2_unitoftimerange key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_unitoftimerange
!
!
     subroutine get_g2_fixedsurfacetypes(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_fixedsurfacetypes
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 - Fixed Surface Types and Units 
!   value for a given short key name based on Table 4.5 of Section 4, Octets 23 and 29
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_fixedsurfacetypes(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key  - GRIB2 character short key for fixed surface types from Table 4.5
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 value from Table 4.5
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXFIXEDSURFACETYPES
        if (trim(table4_5(n)%fixedsurfacetypeskey).eq.trim(key)) then
            value=table4_5(n)%fixedsurfacetypesval
            return
        endif
     enddo

     value=table4_5(66)%fixedsurfacetypesval
           print *,'get_g2_fixedsurfacetypes key ', trim(key), value,  &
                   ' not found.'
           return
     end subroutine get_g2_fixedsurfacetypes
!
!
     subroutine get_g2_statprocesstypes(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_statprocesstypes
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 - Type of statistica
! processing
!   value for a given short key name based on Table 4.10 of Section 4 Octets 47 (template 8)
!   60 (temp 9), 48 (temp 10), 50 (temp 11), 49 (temp 12), 81 (temp 13), 77 (temp 14), 27 (temp 1001),
!   25 (temp 1002) and 39 (temp 1101)
!
! PROGRAM HISTORY LOG:
! 2009-04-02  V. Krishna Kumar
!
! USAGE:    CALL get_g2_statprocesstypes(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key  - GRIB2 character short key for type of statistical processing from Table 4.10
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 value from Table 4.10
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*),intent(in) :: key
!    integer,intent(out) :: value,ierr
     integer :: value,ierr
     integer :: n
!
     do n=1,MAXSTATPROCESSTYPES
        if (trim(table4_10(n)%statprocesstypeskey).eq.key) then
            value=table4_10(n)%statprocesstypesval
            return
        endif
     enddo
           print *,'get_g2_statprocesstypes key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_statprocesstypes
!
!
     subroutine get_g2_typeoftimeintervals(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_typeoftimeintervals
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-04-03
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 - Type of time intervals
!   value for a given short key name based on Table 4.11 of Section 4 Octets 48 (template 8)
!   61 (temp 9), 49 (temp 10), 51 (temp 11), 50 (temp 12), 82 (temp 13), 78 (tem p 14), 28 (temp 1001),
!   and 40 (temp 1101)
!
! PROGRAM HISTORY LOG:
! 2009-04-03  V. Krishna Kumar
!
! USAGE:    CALL get_g2_typeoftimeintervals(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key  - GRIB2 character short key for type of statistical processing from Table 4.11
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 value from Table 4.11
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXTYPEOFTIMEINTVLS
        if (trim(table4_11(n)%typeoftimeintervalskey).eq.key) then
            value=table4_11(n)%typeoftimeintervalsval
            return
        endif
     enddo
           print *,'get_g2_typeoftimeintervals key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_typeoftimeintervals
!
!
     subroutine get_g2_typeofintervals(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_typeofintervals
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-04-03
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 - Type of intervals
!   value for a given short key name based on Table 4.91 of Section 4 Octets 14 (template 44)
!   14 (temp 45), 14 (temp 46), 15 (temp 47),14 and 25 (temp 48)
!
! PROGRAM HISTORY LOG:
! 2012-01-25  Jun Wang : set type of intervals for table 4.91
!
! USAGE:    CALL get_g2_typeofintervals(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key  - GRIB2 character short key for type of intervals from Table 4.91
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 value from Table 4.91
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXTYPEOFINTVLS
        if (trim(table4_91(n)%typeofintervalskey).eq.trim(key)) then
            value=table4_91(n)%typeofintervalsval
            return
        endif
     enddo

     if(trim(key).eq.'') then
       value=255
       return
     endif
           print *,'get_g2_typeofintervals key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_typeofintervals
!
!
     subroutine get_g2_typeofaerosol(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_typeofaerosol
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-04-03
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 - Type of aerosol
!   value for a given short key name based on Table 4.233 of Section 4 Octets 12-13 
!   (template 44), 12-13 (temp 45), 12-13 (temp 46), 13-14 (temp 47),
!   12-13 (temp 48)
!
! PROGRAM HISTORY LOG:
! 2012-01-25  Jun Wang : set type of aerosol for table 4.233
!
! USAGE:    CALL get_g2_typeofaerosol(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key  - GRIB2 character short key for type of aerosol from Table 4.233
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 value from Table 4.233
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXTYPEOFAEROSOL
        if (trim(table4_233(n)%typeofaerosolkey).eq.trim(key)) then
            value=table4_233(n)%typeofaerosolval
            return
        endif
     enddo

     if (trim(key).eq.'') then
       value=65535
       return
     endif
           print *,'get_g2_typeofaerosol key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_typeofaerosol
!
!
     subroutine get_g2_on388origincenters(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_on388origincenters
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB1 - PDS Ocet5            
!   GRIB2 - Section 1, Octet 6-7 National/International Originating Centers
!   value for a given short key name based on ON388 - Table 0 
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_on388origincenters(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB1 character short key for originating center based on ON388 - Table 0 
!
!   OUTPUT ARGUMENT LIST:
!     value  - corresponding GRIB1-PDS Octet 5/GRIB2-Section 1, Octets 6-7 value from ON388 - Table 0 
!     ierr   - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXORIGINCENTERS
        if (trim(on388_table0(n)%origincenterskey).eq.trim(key)) then
            value=on388_table0(n)%origincentersval
            return
        endif
     enddo
           print *,'get_g2_on388origincenters key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_on388origincenters
!
!
     subroutine get_g2_on388genproc(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_on388genproc
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2009-12-10
!
! ABSTRACT: This subroutine returns the corresponding GRIB1 - PDS Ocet6
!   data value (Generating process or model) from originating center 7 (USNWS NCEP)
!   for a given short key name based on ON388 - Table A
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_on388genproc(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB1 character short key for model based on ON388 - Table A
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB1 - PDS Octet 6 value from ON388 - Table A
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXGENPROC
        if (trim(on388_tablea(n)%genprockey).eq.trim(key)) then
            value=on388_tablea(n)%genprocval
            return
        endif
     enddo
           print *,'get_g2_on388genproc key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_on388genproc
!
!
     subroutine get_g2_typeoforigfieldvals(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_typeoforigfieldvals
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-03-08
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 - Table 5.1
!   Type of Original Field Values for a given short key name based on GRIB2 - Table 5.1
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_typeoforigfieldvals(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for type of original field values based on Table 5.1
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 - Table 5.1 value
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXTYPEOFORIGFIELDVAL
        if (trim(table5_1(n)%typeoforigfieldvalskey).eq.trim(key)) then
            value=table5_1(n)%typeoforigfieldvals
            return
        endif
     enddo

           print *,'get_g2_typeoforigfieldvals key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_typeoforigfieldvals
!
!
     subroutine get_g2_ordofspcdiffvals(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_grpspltmthdvals
!   PRGMMR: J.Wang                   ORG: NCEP/EMC    DATE: 2012-02-20
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 - Table 5.6
!   Order of spatial differencing for a given short key name based on GRIB2 - Table 5.6
!   default is 1st order spatial differencing
!
! PROGRAM HISTORY LOG:
! 2012-02-20  J.Wang
!
! USAGE:    CALL get_g2_ordofspcdiff(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for Order of spatial differencing based on Table 5.6
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 - Table 5.6 value
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXORDOFSPTDIFF
        if (trim(table5_6(n)%ordofsptdiffkey).eq.trim(key)) then
            value=table5_6(n)%ordofsptdiffvals
            return
        endif
     enddo
           print *,'get_g2_ordofsptdiffvals key ', key,   &
                   ' not found.'
           value=1
           return
     end subroutine get_g2_ordofspcdiffvals
!
!
     subroutine get_g2_typeofcompression(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_typeofcompression
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-03-08
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 - Table 5.40
!   Type of compression for a given short key name based on GRIB2 - Table 5.40
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
!
! USAGE:    CALL get_g2_typeofcompression(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for type of compression based on Table 5.40
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 - Table 5.40 value
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXTYPEOFCOMPRESSION
        if (trim(table5_40(n)%typeofcompressionkey).eq.trim(key)) then
            value=table5_40(n)%typeofcompressionvals
            return
        endif
     enddo
           print *,'get_g2_typeofcompression key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_typeofcompression
!
!
     subroutine get_g2_sec5packingmethod(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_sec5tmplnum
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-03-08
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 template number
!   for a given short key name based on GRIB2 - 
!
! PROGRAM HISTORY LOG:
! 2009-12-10  V. Krishna Kumar
! 2010-03-15  Jun Wang : get section5 template number
!
! USAGE:    CALL get_g2_sec5tmplnum(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for packing method based on Table 5.0
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 - Table 5.0 value
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXTYPEOFPACKINGMETHOD
        if (trim(table5_0(n)%packingmethodkey).eq.trim(key)) then
            value=table5_0(n)%packingmethodvals
            return
        endif
     enddo
     print *,'get_g2_sec5packingmethod key ', key,   &
             ' not found.'
     return
     end subroutine get_g2_sec5packingmethod
!
!
     subroutine g2sec0(idisc,listsec0)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec0
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-03-01
!
! ABSTRACT: This subroutine returns the section 0 list for a given discipline
!   value 
!
! PROGRAM HISTORY LOG:
! 2010-03-01  V. Krishna Kumar
!
! USAGE:    CALL g2sec0(idisc,listsec0)
!   INPUT ARGUMENT LIST:
!     idisc  - GRIB2 Discipline (From Table 0.0)
!
!   OUTPUT ARRAY:
!     listsec0(1)  - GRIB2 Discipline (From Table 0.0)
!     listsec0(2)  - Edition number - 2 for GRIB2
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     integer :: idisc
     integer :: listsec0(2)
!
     listsec0(1) = idisc
     listsec0(2) = 2      ! Edition number - 2 for GRIB2
     end subroutine g2sec0
!
!
     subroutine g2sec1(origin_key,subcen_key,vers_key,lvers_key,sigreftime_key,refyear_val, &
                       refmon_val,refday_val,refhour_val,refmin_val,refsec_val,prodstatus_key, &
                       typeofdata_key,listsec1)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec1
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-03-01
!
! ABSTRACT: This subroutine returns the section 1 list for given keys
!    
! PROGRAM HISTORY LOG:
! 2010-03-01  V. Krishna Kumar
!
! USAGE:    CALL g2sec1(origin_key,subcen_key,vers_key,lvers_key,sigreftime_key,refyear_val, &
!                       refmon_val,refday_val,refhour_val,refmin_val,refsec_val,prodstatus_key, &
!                       typeofdata_key,listsec1)
!   INPUT ARGUMENT LIST:
!     origin_key - Identification of originating/generating center (See Table 0 {GRIB1})
!     subcen_key - Identification of originating/generating subcenter (See Table C)
!     vers_key - GRIB master tables version number (currently 2) (See Table 1.0) (See note 1 below)
!     lvers_key - Version number of GRIB local tables used to augment Master Tables (see Table 1.1)
!     sigreftime_key - Significance of reference time (See Table 1.2) 
!     refyear_val - Year (4 digits)
!     refmon_val - Month
!     refday_val - Day
!     refhour_val - Hour
!     refmin_val - Minute
!     refsec_val - Second
!     prodstatus_key - Production Status of Processed data in the GRIB message (See Table 1.3)
!     typeofdata_key - Type of processed data in this GRIB message (See Table 1.4)
!
!   OUTPUT ARRAY:
!     listsec1  - GRIB2 Section 1 Identification Section values array 
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
!    integer,intent(inout) :: listsec1(13)
     integer :: listsec1(13)
     integer :: refyear_val,refmon_val,refday_val,refhour_val,refmin_val,refsec_val
     character(len=*) :: origin_key,subcen_key,vers_key,lvers_key,  &
                          sigreftime_key,prodstatus_key,typeofdata_key
!
     integer(4) :: value,ierr
!
     call get_g2_on388origincenters(origin_key,value,ierr)
     listsec1(1) = value
!
     call get_g2_subcenters(subcen_key,value,ierr) 
     listsec1(2) = value
!
     call get_g2_versionno(vers_key,value,ierr)
     listsec1(3) = value
!    
     call get_g2_loctabversno(lvers_key,value,ierr)
     listsec1(4) = value
!
     call get_g2_sigreftime(sigreftime_key,value,ierr)
     listsec1(5) = value
!
! Set the time yyyy,mm,dd,hh,min,sec
!
     listsec1(6) = refyear_val
     listsec1(7) = refmon_val
     listsec1(8) = refday_val
     listsec1(9) = refhour_val
     listsec1(10) = refmin_val
     listsec1(11) = refsec_val
!
     call get_g2_prodstatus(prodstatus_key,value,ierr)
     listsec1(12) = value
!
     call get_g2_typeofdata(typeofdata_key,value,ierr)
     listsec1(13) = value
!
     end subroutine g2sec1
!
!
     subroutine g2sec4_temp0(icatg,iparm,typ_gen_proc_key,                         &
                             gen_proc_or_mod_key,hrs_obs_cutoff,min_obs_cutoff,    &
                             unit_of_time_key,fcst_time,lvl_type1,scale_fac1,      &
                             scaled_val1,lvl_type2,scale_fac2,scaled_val2,         &
                             ipdstmpl0)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec4_temp0
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-03-01
!
! ABSTRACT: This subroutine returns the Grib2 Section 4 Template 4.0 list for given keys
!           PDT 4.0 - Analysis or forecast at a horizontal level or in a
!                     horizontal layer at a point in time.
!
! PROGRAM HISTORY LOG:
! 2010-03-01  V. Krishna Kumar
!
! USAGE:    CALL g2sec4_temp0(icatg,iparm,typ_gen_proc_key,gen_proc_or_mod_key,  
!                             hrs_obs_cutoff,min_obs_cutoff,unit_of_time_key,  
!                             fcst_time,lvl_type1,scale_fac1,scaled_val1,lvl_type2,  
!                             scale_fac2,scaled_val2,ipdstmpl0)
!   INPUT ARGUMENT LIST:
!      icatg - Parameter category (see Code table 4.1)
!      iparm - Parameter number (see Code table 4.2)
!      typ_gen_proc_key - Type of generating process (see Code table 4.3)
!      bckgnd_gen_proc_id - Background generating process identifier (defined by originating centre)
!      gen_proc_or_mod_key - Analysis or forecast generating process identified (see Code ON388 Table A)
!      hrs_obs_cutoff - Hours of observational data cutoff after reference time (see Note)
!      min_obs_cutoff - Minutes of observational data cutoff after reference time (see Note)
!      unit_of_time_key - Indicator of unit of time range (see Code table 4.4)
!      fcst_time - Forecast time in units defined by octet 18
!      lvl_type1 - Type of first fixed surface (see Code table 4.5)
!      scale_fac1 - Scale factor of first fixed surface
!      scaled_val1 - Scaled value of first fixed surface
!      lvl_type2 - Type of second fixed surfaced (see Code table 4.5)
!      scale_fac2 - Scale factor of second fixed surface
!      scaled_val2 - Scaled value of second fixed surfaces
!
!   OUTPUT ARRAY:
!      ipdstmpl0  - GRIB2 PDS Template 4.0 listing
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
     integer(4),intent(in) :: icatg,iparm,hrs_obs_cutoff,min_obs_cutoff,         &
                   fcst_time,scale_fac1,scaled_val1,scale_fac2,scaled_val2
!     integer(4),intent(inout)  :: bckgnd_gen_proc_id    ! defined by the center
!
     character(len=*),intent(in) :: typ_gen_proc_key,gen_proc_or_mod_key,       &
                          unit_of_time_key,lvl_type1,lvl_type2 
!
     integer(4),intent(inout)  :: ipdstmpl0(15)
!
!local vars
     integer(4) :: value,ierr
     integer(4) :: bckgnd_gen_proc_id    ! defined by the center
!
     bckgnd_gen_proc_id=0    ! defined by the center
!
     ipdstmpl0(1) = icatg
     ipdstmpl0(2) = iparm
!
     call get_g2_typeofgenproc(typ_gen_proc_key,value,ierr)
     ipdstmpl0(3) = value
!
     ipdstmpl0(4) = bckgnd_gen_proc_id
!
     call get_g2_on388genproc(gen_proc_or_mod_key,value,ierr)
     ipdstmpl0(5) = value
!
     ipdstmpl0(6) = hrs_obs_cutoff
     ipdstmpl0(7) = min_obs_cutoff
!
     call get_g2_unitoftimerange(unit_of_time_key,value,ierr)
     ipdstmpl0(8) = value
     ipdstmpl0(9) = fcst_time
!
     call get_g2_fixedsurfacetypes(lvl_type1,value,ierr)
     ipdstmpl0(10) = value
     ipdstmpl0(11) = scale_fac1
     ipdstmpl0(12) = scaled_val1
!
     call get_g2_fixedsurfacetypes(lvl_type2,value,ierr)
     ipdstmpl0(13) = value
!
     ipdstmpl0(14) = scale_fac2
     ipdstmpl0(15) = scaled_val2
!
     end subroutine g2sec4_temp0
!
!
     subroutine g2sec4_temp8(icatg,iparm,typ_gen_proc_key,gen_proc_or_mod_key,     &
                             hrs_obs_cutoff,min_obs_cutoff,unit_of_time_key,       &
                             fcst_time,lvl_type1,scale_fac1,scaled_val1,lvl_type2, &
                             scale_fac2,scaled_val2,year_intvl,                    &
                             mon_intvl,day_intvl,hour_intvl,min_intvl,sec_intvl,   &
                             num_time_range,stat_miss_val,type_of_stat_proc,       &
                             type_of_time_inc,stat_unit_time_key,                  &
                             leng_time_range_stat,stat_unit_time_key_succ,         &
                             time_inc_betwn_succ_fld,ipdstmpl8)
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec4_temp8
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-03-01
!
! ABSTRACT: This subroutine returns the Grib2 Section 4 Template 4.8 list for given keys
!           PDT 4.8 - Average, accumulation, extreme values or other statistically 
!                     processed values at a horizontal level or in a horizontal layer
!                     in a continuous or non-continuous time interval.
!
! PROGRAM HISTORY LOG:
! 2010-03-01  V. Krishna Kumar
! 2010-04-20  Jun Wang
!
! USAGE:    CALL g2sec4_temp8(icatg,iparm,typ_gen_proc_key,gen_proc_or_mod_key, &
!                             hrs_obs_cutoff,min_obs_cutoff,unit_of_time_key, &
!                             fcst_time,lvl_type1,scale_fac1,scaled_val1,lvl_type2, &
!                             scale_fac2,scaled_val2,year_intvl, &
!                             mon_intvl,day_intvl,hour_intvl,min_intvl,sec_intvl, &
!                             num_time_range,stat_miss_val,type_of_stat_proc, &
!                             type_of_time_inc,stat_unit_time_key, &
!                             leng_time_range_stat,stat_unit_time_key_succ, &
!                             time_inc_betwn_succ_fld,ipdstmpl8)
!   INPUT ARGUMENT LIST:
!
!	icatg - Parameter category (see Code Table 4.1)
!	iparm - Parameter number (see Code Table 4.2)
!	typ_gen_proc_key - Type of generating process (see Code Table 4.3)
!	bckgnd_gen_proc_id - Background generating process identifier (defined by originating centre)
!	gen_proc_or_mod_key - Analysis or forecast generating process identified (see Code ON388 Table A)
!	hrs_obs_cutoff - Hours after reference time data cutoff (see Note 1)
!	min_obs_cutoff - Minutes after reference time data cutoff
!	unit_of_time_key - Indicator of unit of time range (see Code Table 4.4)
!	fcst_time - Forecast time in units defined by octet 18 (see Note 2)
!	lvl_type1 - Type of first fixed surface (see Code Table 4.5)
!	scale_fac1 - Scale factor of first fixed surface
!	scaled_val1 - Scaled value of first fixed surface
!	lvl_type2 - Type of second fixed surfaced (see Code Table 4.5)
!	scale_fac2 - Scale factor of second fixed surface
!	scaled_val2 - Scaled value of second fixed surfaces
!       year_intvl - Year Time of end of overall time interval
!       mon_intvl - Month Time of end of overall time interval
!       day_intvl - Day Time of end of overall time interval
!       hour_intvl - Hour Time of end of overall time interval
!       min_intvl - Minute Time of end of overall time interval
!       sec_intvl - Second Time of end of overall time interval
!       num_time_range - n number of time ranges specifications describing
!                        the time intervals used to calculate the
!                        statistically-processed field
!       stat_miss_val - Total number of data values missing in statistical process
!                       Specification of the outermost (or only) time range over
!                       which statistical processing is done
!       type_of_stat_proc - Statistical process used to calculate the processed
!                           field from the field at each time increment during the
!                           time range (see Code Table 4.10)
!       type_of_time_inc - Type of time increment between successive fields
!                          used in the statistical processing (see Code Table 4.11)
!       stat_unit_time_key - Indicator of unit of time for time range over which
!                            statistical processing is done (see Code Table 4.4)
!       leng_time_range_stat - Length of the time range over which statistical processing
!                              is done, in units defined by the previous octet
!       stat_unit_time_key_succ - Indicator of unit of time for the increment between the
!                                 successive fields used (see Code table 4.4)
!       time_inc_betwn_succ_fld - Time increment between successive fields,
!                                  in units defined by the previous octet (see Notes 3 & 4)
!
!   OUTPUT ARRAY:
!      ipdstmpl8  - GRIB2 PDS Template 4.8 listing
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
     integer(4),intent(in) :: icatg,iparm,hrs_obs_cutoff,min_obs_cutoff,fcst_time, &
                              scale_fac1,scaled_val1,scale_fac2,scaled_val2 
     integer(4),intent(in) :: year_intvl,mon_intvl,day_intvl,hour_intvl,min_intvl, &
                              sec_intvl,num_time_range,stat_miss_val, &
                              leng_time_range_stat,time_inc_betwn_succ_fld
!
     character(len=*),intent(in) :: typ_gen_proc_key,gen_proc_or_mod_key, &
                          unit_of_time_key,lvl_type1,lvl_type2, &
                          type_of_stat_proc,type_of_time_inc, &
                          stat_unit_time_key,stat_unit_time_key_succ
!
     integer(4)               :: bckgnd_gen_proc_id    ! defined by the center
!
     integer(4),intent(inout) :: ipdstmpl8(29)         ! currently works only for n=1
                                                       ! later on, this will be generalized
!
!-- local vars
     integer(4) :: value,ierr
!
     bckgnd_gen_proc_id=0
!
     ipdstmpl8(1) = icatg
     ipdstmpl8(2) = iparm
!
     call get_g2_typeofgenproc(typ_gen_proc_key,value,ierr)
     ipdstmpl8(3) = value
!
     ipdstmpl8(4) = bckgnd_gen_proc_id
!
     call get_g2_on388genproc(gen_proc_or_mod_key,value,ierr)
     ipdstmpl8(5) = value
!
     ipdstmpl8(6) = hrs_obs_cutoff
     ipdstmpl8(7) = min_obs_cutoff
!
     call get_g2_unitoftimerange(unit_of_time_key,value,ierr)
     ipdstmpl8(8) = value
     ipdstmpl8(9) = fcst_time
!
     call get_g2_fixedsurfacetypes(lvl_type1,value,ierr)
     ipdstmpl8(10) = value
     ipdstmpl8(11) = scale_fac1
     ipdstmpl8(12) = scaled_val1
!
     call get_g2_fixedsurfacetypes(lvl_type2,value,ierr)
     ipdstmpl8(13) = value
!
     ipdstmpl8(14) = scale_fac2
     ipdstmpl8(15) = scaled_val2
     ipdstmpl8(16) = year_intvl
     ipdstmpl8(17) = mon_intvl
     ipdstmpl8(18) = day_intvl
     ipdstmpl8(19) = hour_intvl
     ipdstmpl8(20) = min_intvl
     ipdstmpl8(21) = sec_intvl
!
     ipdstmpl8(22) = num_time_range ! choose n=1 for this case
     ipdstmpl8(23) = stat_miss_val  ! choose 0 for this case
!
     call get_g2_statprocesstypes(type_of_stat_proc,value,ierr)
     ipdstmpl8(24) = value  ! types_of_stat_proc='accumulation'
!
     call get_g2_typeoftimeintervals(type_of_time_inc,value,ierr)
     ipdstmpl8(25) = value  ! type_of_time_inc='same_start_time_fcst_fcst_time_inc'
                            ! value = 2 (Successive times processed have same start
                            !       time of forecast, forecast time is incremented)
!
     call get_g2_unitoftimerange(stat_unit_time_key,value,ierr)
     ipdstmpl8(26) = value  ! stat_unit_time_key='hour'
                            ! value = 1
     ipdstmpl8(27) = leng_time_range_stat  ! value = 6
!
     call get_g2_unitoftimerange(stat_unit_time_key_succ,value,ierr)
                            ! stat_unit_time_key_succ='missing'
     ipdstmpl8(28) = value  ! value = 255
!
     ipdstmpl8(29) = time_inc_betwn_succ_fld   ! value = 0
!
     end subroutine g2sec4_temp8
!
!
     subroutine g2sec4_temp44(icatg,iparm,aer_type,typ_intvl_size,                 &
                             scale_fac1_size,scale_val1_size,scale_fac2_size,      &
                             scale_val2_size,typ_gen_proc_key,                     &
                             gen_proc_or_mod_key,hrs_obs_cutoff,min_obs_cutoff,    &
                             unit_of_time_key,fcst_time,lvl_type1,scale_fac1,      &
                             scaled_val1,lvl_type2,scale_fac2,scaled_val2,         &
                             ipdstmpl44)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec4_temp44
!   PRGMMR: J. WANG                  ORG: NCEP/EMC  DATE: 2012-01-25
!
! ABSTRACT: This subroutine returns the Grib2 Section 4 Template 4.44 list for given keys
!           PDT 4.44 - Analysis or forecast at a horizontal level or in a 
!                      horizontal layer at a point in time for aerosol
!
! PROGRAM HISTORY LOG:
! 2012-01-25  Jun Wang        generate pds template 4.44
!
! USAGE:    CALL g2sec4_temp44(icatg,iparm,aer_type,typ_intvl_size,scale_fac1_size,
!                             scale_val1_size,scale_fac2_size,scale_val2_size,
!                             typ_gen_proc_key,gen_proc_or_mod_key,  
!                             hrs_obs_cutoff,min_obs_cutoff,unit_of_time_key,  
!                             fcst_time,lvl_type1,scale_fac1,scaled_val1,lvl_type2,  
!                             scale_fac2,scaled_val2,ipdstmpl44)
!   INPUT ARGUMENT LIST:
!      icatg - Parameter category (see Code table 4.1)
!      iparm - Parameter number (see Code table 4.2)
!      aer_type - Aetosol type (see Code table 4.233)
!      typ_intvl_size - Type of interval for first and second size (see Code table 4.91)
!      scale_fac1_size - Scale factor of first size
!      scale_val1_size - Scale value of first size in meters
!      scale_fac2_size - Scale factor of second size
!      scale_val2_size - Scale value of second size in meters
!      typ_gen_proc_key - Type of generating process (see Code table 4.3)
!      bckgnd_gen_proc_id - Background generating process identifier (defined by originating centre)
!      gen_proc_or_mod_key - Analysis or forecast generating process identified (see Code ON388 Table A)
!      hrs_obs_cutoff - Hours of observational data cutoff after reference time (see Note)
!      min_obs_cutoff - Minutes of observational data cutoff after reference time (see Note)
!      unit_of_time_key - Indicator of unit of time range (see Code table 4.4)
!      fcst_time - Forecast time in units defined by octet 18
!      lvl_type1 - Type of first fixed surface (see Code table 4.5)
!      scale_fac1 - Scale factor of first fixed surface
!      scaled_val1 - Scaled value of first fixed surface
!      lvl_type2 - Type of second fixed surfaced (see Code table 4.5)
!      scale_fac2 - Scale factor of second fixed surface
!      scaled_val2 - Scaled value of second fixed surfaces
!
!   OUTPUT ARRAY:
!      ipdstmpl44  - GRIB2 PDS Template 4.44 listing
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
     integer(4),intent(in) :: icatg,iparm,hrs_obs_cutoff,min_obs_cutoff,         &
                   scale_fac1_size,scale_fac2_size,                              &
                   fcst_time,scale_fac1,scaled_val1, scale_fac2,scaled_val2
     real,intent(in) :: scale_val1_size,scale_val2_size
!
     character(len=*),intent(in) :: aer_type,typ_intvl_size,typ_gen_proc_key,    &
                   gen_proc_or_mod_key,unit_of_time_key,lvl_type1,lvl_type2 
!
     integer(4),intent(inout)  :: ipdstmpl44(21)
!
!local vars
     integer(4) :: value,ierr
     integer(4) :: bckgnd_gen_proc_id    ! defined by the center
!
     bckgnd_gen_proc_id=0    ! defined by the center
!
     ipdstmpl44(1) = icatg
     ipdstmpl44(2) = iparm
!
     call get_g2_typeofaerosol(aer_type,value,ierr)
     ipdstmpl44(3) = value
!
     call get_g2_typeofintervals(typ_intvl_size,value,ierr)
     ipdstmpl44(4) = value
     ipdstmpl44(5) = scale_fac1_size
     ipdstmpl44(6) = scale_val1_size
     ipdstmpl44(7) = scale_fac2_size
     ipdstmpl44(8) = scale_val2_size
!
     call get_g2_typeofgenproc(typ_gen_proc_key,value,ierr)
     ipdstmpl44(9) = value
!
     ipdstmpl44(10) = bckgnd_gen_proc_id
!
     call get_g2_on388genproc(gen_proc_or_mod_key,value,ierr)
     ipdstmpl44(11) = value
!
     ipdstmpl44(12) = hrs_obs_cutoff
     ipdstmpl44(13) = min_obs_cutoff
!
     call get_g2_unitoftimerange(unit_of_time_key,value,ierr)
     ipdstmpl44(14) = value
     ipdstmpl44(15) = fcst_time
!
     call get_g2_fixedsurfacetypes(lvl_type1,value,ierr)
     ipdstmpl44(16) = value
     ipdstmpl44(17) = scale_fac1
     ipdstmpl44(18) = scaled_val1
!
     call get_g2_fixedsurfacetypes(lvl_type2,value,ierr)
     ipdstmpl44(19) = value
!
     ipdstmpl44(20) = scale_fac2
     ipdstmpl44(21) = scaled_val2
!
     end subroutine g2sec4_temp44
!
!
     subroutine g2sec4_temp48(icatg,iparm,aer_type,typ_intvl_size,                 &
                             scale_fac1_size,scale_val1_size,scale_fac2_size,      &
                             scale_val2_size,typ_intvl_wavelength,                 &
                             scale_fac1_wavelength,scale_val1_wavelength,          &
                             scale_fac2_wavelength,scale_val2_wavelength,          &
                             typ_gen_proc_key, gen_proc_or_mod_key,                &
                             hrs_obs_cutoff,min_obs_cutoff,                        &
                             unit_of_time_key,fcst_time,lvl_type1,scale_fac1,      &
                             scaled_val1,lvl_type2,scale_fac2,scaled_val2,         &
                             ipdstmpl48)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec4_temp48
!   PRGMMR: J. WANG                  ORG: NCEP/EMC  DATE: 2012-01-25
!
! ABSTRACT: This subroutine returns the Grib2 Section 4 Template 4.0 list for given keys
!           PDT 4.48 - Analysis or forecast at a horizontal level or in a
!                      horizontal layer at a point in time for aerosol. 
!
! PROGRAM HISTORY LOG:
! 2012-01-25  Jun Wang        generate pds template 4.48
!
! USAGE:    CALL g2sec4_temp48(icatg,iparm,aer_type,typ_intvl_size,scale_fac1_size,
!                             scale_val1_size,scale_fac2_size,scale_val2_size,
!                             typ_intvl_wavelength,scale_val1_wavelength,
!                             scale_val1_wavelength,scale_fac2_wavelength,
!                             scale_val2_wavelength,
!                             typ_gen_proc_key,gen_proc_or_mod_key,  
!                             hrs_obs_cutoff,min_obs_cutoff,unit_of_time_key,  
!                             fcst_time,lvl_type1,scale_fac1,scaled_val1,lvl_type2,  
!                             scale_fac2,scaled_val2,ipdstmpl0)
!   INPUT ARGUMENT LIST:
!      icatg - Parameter category (see Code table 4.1)
!      iparm - Parameter number (see Code table 4.2)
!      aer_type - Aetosol type (see Code table 4.233)
!      typ_intvl_size - Type of interval for first and second size (see Code table 4.91)
!      scale_fac1_size - Scale factor of first size
!      scale_val1_size - Scale value of first size in meters
!      scale_fac2_size - Scale factor of second size
!      scale_val2_size - Scale value of second size in meters
!      typ_intvl_wavelength - Type of interval for first and second wavelength (see Code table 4.91)
!      scale_fac1_wavelength - Scale factor of first wavelength
!      scale_val1_wavelength - Scale value of first wavelength in meters
!      scale_fac2_wavelength - Scale factor of second wavelength
!      scale_val2_wavelength - Scale value of second wavelength in meters
!      typ_gen_proc_key - Type of generating process (see Code table 4.3)
!      bckgnd_gen_proc_id - Background generating process identifier (defined by originating centre)
!      gen_proc_or_mod_key - Analysis or forecast generating process identified (see Code ON388 Table A)
!      hrs_obs_cutoff - Hours of observational data cutoff after reference time (see Note)
!      min_obs_cutoff - Minutes of observational data cutoff after reference time (see Note)
!      unit_of_time_key - Indicator of unit of time range (see Code table 4.4)
!      fcst_time - Forecast time in units defined by octet 18
!      lvl_type1 - Type of first fixed surface (see Code table 4.5)
!      scale_fac1 - Scale factor of first fixed surface
!      scaled_val1 - Scaled value of first fixed surface
!      lvl_type2 - Type of second fixed surfaced (see Code table 4.5)
!      scale_fac2 - Scale factor of second fixed surface
!      scaled_val2 - Scaled value of second fixed surfaces
!
!   OUTPUT ARRAY:
!      ipdstmpl48  - GRIB2 PDS Template 4.48 listing
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
     integer(4),intent(in) :: icatg,iparm,hrs_obs_cutoff,min_obs_cutoff,         &
                   scale_fac1_size,scale_fac2_size, scale_fac1_wavelength,       &
                   scale_fac2_wavelength,                                        &
                   fcst_time,scale_fac1,scaled_val1,                             &
                   scale_fac2,scaled_val2
     real,intent(in) :: scale_val1_size,scale_val2_size,scale_val1_wavelength,   &
                   scale_val2_wavelength
!
     character(len=*),intent(in) :: aer_type,typ_intvl_size,                     &
                   typ_intvl_wavelength,typ_gen_proc_key,                        &
                   gen_proc_or_mod_key,unit_of_time_key,lvl_type1,lvl_type2 
!
     integer(4),intent(inout)  :: ipdstmpl48(26)
!
!local vars
     integer(4) :: value,ierr
     integer(4) :: bckgnd_gen_proc_id    ! defined by the center
!
     bckgnd_gen_proc_id=0    ! defined by the center
!
     ipdstmpl48(1) = icatg
     ipdstmpl48(2) = iparm
!
     call get_g2_typeofaerosol(aer_type,value,ierr)
     ipdstmpl48(3) = value
!
     call get_g2_typeofintervals(typ_intvl_size,value,ierr)
     ipdstmpl48(4) = value
     ipdstmpl48(5) = scale_fac1_size
     ipdstmpl48(6) = nint(scale_val1_size)
     ipdstmpl48(7) = scale_fac2_size
     ipdstmpl48(8) = nint(scale_val2_size)
!
     call get_g2_typeofintervals(typ_intvl_wavelength,value,ierr)
     ipdstmpl48(9) = value
     ipdstmpl48(10) = scale_fac1_wavelength
     ipdstmpl48(11) = nint(scale_val1_wavelength)
     ipdstmpl48(12) = scale_fac2_wavelength
     ipdstmpl48(13) = nint(scale_val2_wavelength)
!
     call get_g2_typeofgenproc(typ_gen_proc_key,value,ierr)
     ipdstmpl48(14) = value
!
     ipdstmpl48(15) = bckgnd_gen_proc_id
!
     call get_g2_on388genproc(gen_proc_or_mod_key,value,ierr)
     ipdstmpl48(16) = value
!
     ipdstmpl48(17) = hrs_obs_cutoff
     ipdstmpl48(18) = min_obs_cutoff
!
     call get_g2_unitoftimerange(unit_of_time_key,value,ierr)
     ipdstmpl48(19) = value
     ipdstmpl48(20) = fcst_time
!
     call get_g2_fixedsurfacetypes(lvl_type1,value,ierr)
     ipdstmpl48(21) = value
     ipdstmpl48(22) = scale_fac1
     ipdstmpl48(23) = scaled_val1
!
     call get_g2_fixedsurfacetypes(lvl_type2,value,ierr)
     ipdstmpl48(24) = value
!
     ipdstmpl48(25) = scale_fac2
     ipdstmpl48(26) = scaled_val2
!
     end subroutine g2sec4_temp48
!
!
     subroutine get_g2_typeofensfcst(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_typeofensfcst
!   PRGMMR: Boi Vuong         ORG: W/SIB     DATE: 2015-01-09
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 type of
!    ensemble forecast value for a given short key name based on Table 4.6
!
! PROGRAM HISTORY LOG:
! 2015-01-09  Boi vuong
!
! USAGE:    CALL get_g2_typeofensfcst(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for type of ensemble forecast
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 type of ensemble forecast value from table 4.6
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXTYPEOFENSFCST
        if (trim(table4_6(n)%typeofensfcstkey).eq.trim(key)) then
            value=table4_6(n)%typeofensfcstval
            return
        endif
     enddo
           print *,'get_g2_typeofensfcst key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_typeofensfcst
!
!
     subroutine get_g2_typeofderivefcst(key,value,ierr)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    get_g2_typeofderivefcst
!   PRGMMR: Boi Vuong         ORG: W/SIB     DATE: 2015-01-09
!
! ABSTRACT: This subroutine returns the corresponding GRIB2 type of
!    derive forecast value for a given short key name based on Table 4.7
!
! PROGRAM HISTORY LOG:
! 2015-01-09  Boi vuong
!
! USAGE:    CALL get_g2_typeofderivefcst(key,value,ierr)
!   INPUT ARGUMENT LIST:
!     key      - GRIB2 character short key for type of derive forecast
!
!   OUTPUT ARGUMENT LIST:
!     value    - corresponding GRIB2 type of derive forecast value from table 4.7
!     ierr     - error messages
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     character(len=*) :: key
     integer :: value,n,ierr
!
     do n=1,MAXTYPEOFDERIVEFCST
        if (trim(table4_7(n)%typeofderivefcstkey).eq.trim(key)) then
            value=table4_7(n)%typeofderivefcstval
            return
        endif
     enddo
           print *,'get_g2_typeofderivefcst key ', key,   &
                   ' not found.'
           return
     end subroutine get_g2_typeofderivefcst
!
!
     subroutine g2sec5_temp0(dec_scale_fac,bin_scale_fac,tlnumbits,ifield5)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec5_temp0
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-03-01
!
! ABSTRACT: This subroutine returns the section 5 list array for a given decimal 
!   scale factor (D) and type of original field values (Table 5.1) value from
!   GRIB2 - GRID Template 5.0 Grid point data - simple packing
!
! PROGRAM HISTORY LOG:
! 2010-03-01  V. Krishna Kumar
! 2012-02-21  J. Wang   add binary scale factor and number of bits in argument list
!
! USAGE:    CALL g2sec5_temp0(dec_scale_fac,bin_scale_fac,tlnumbits,ifield5)
!   INPUT ARGUMENT LIST:
!     dec_scale_fac  - Decimal scale factor (E)
!     bin_scale_fac  - binary scale factor (D)
!     tlnumbits      - Number of bits used 
!
!   OUTPUT ARRAY:
!     ifield5  - GRIB2 - GRID Template 5.0  listing
!
!********************************************************************************
! ifield5(1): reference value(R) (IEEE 32-bit floating-point value)             *
! ifield5(2): binary scale factor (E)                                           *
! ifield5(3): decimal scale factor (D)                                          *
! ifield5(4): number of bits used for each packed value for simple packing      *
!             or for each group reference value for complex packing or          *
!             spatial differencing                                              *
! ifield5(5): type of original field values (See Code Table 5.1)                *
!********************************************************************************
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     integer(4),intent(in) :: bin_scale_fac,dec_scale_fac,tlnumbits
     integer(4),intent(out) :: ifield5(5)
!     character(len=50) :: type_of_field
     integer(4) :: value,ierr
!
     ifield5(1) = 0 ! Any value. Will be later overwritten
     ifield5(2) = bin_scale_fac     
     ifield5(3) = dec_scale_fac
     ifield5(4) = tlnumbits  
     ifield5(5) = 0             !g2 lib only 0 
!
     end subroutine g2sec5_temp0
!
!
     subroutine g2sec5_temp2(dec_scale_fac,bin_scale_fac,ifield5)
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec5_temp2
!   PRGMMR: J.Wang                   ORG: W/EMC    DATE: 2012-02-20
!
! ABSTRACT: This subroutine returns the section 5 list array with a given binary,
!    and decimal scale factor from GRIB2 - GRID Template 5.2 Grid point data - 
!    complex packing 
!
! PROGRAM HISTORY LOG:
! 2012-02-21  J. Wang
!
! USAGE:    CALL g2sec5_temp2(dec_scale_fac,bin_scale_fac,ifield5)
!   INPUT ARGUMENT LIST:
!     bin_scale_fac  - binary scale factor (E)
!     dec_scale_fac  - Decimal scale factor (D)
!
!   OUTPUT ARRAY:
!     ifield5  - GRIB2 - GRID Template 5.2  listing
!
!********************************************************************************
! ifield5(1): reference value(R) (IEEE 32-bit floating-point value)             *
! ifield5(2): binary scale factor (E)                                           *
! ifield5(3): decimal scale factor (D)                                          *
! ifield5(4): number of bits used for each packed value for simple packing      *
!             or for each group reference value for complex packing or          *
!             spatial differencing                                              *
!********************************************************************************
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     integer(4),intent(inout)  :: ifield5(16)
     integer(4),intent(in) :: dec_scale_fac,bin_scale_fac
!
     integer(4) :: value,ierr
!
     ifield5=0
     ifield5(1) = 0 ! Any value. Will be later overwritten
     ifield5(2) = bin_scale_fac
     ifield5(3) = dec_scale_fac
!
     end subroutine g2sec5_temp2
!
!
     subroutine g2sec5_temp3(dec_scale_fac,bin_scale_fac,order_of_sptdiff,   &
       ifield5)
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec5_temp3
!   PRGMMR: J.Wang                   ORG: W/EMC    DATE: 2012-02-20
!
! ABSTRACT: This subroutine returns the section 5 list array with a given binary,
!    and decimal scale factor from GRIB2 - GRID Template 5.3 Grid point data - 
!    complex packing with spatial difference
!
! PROGRAM HISTORY LOG:
! 2012-02-21  J. Wang
!
! USAGE:    CALL g2sec5_temp3(dec_scale_fac,bin_scale_fac,order_of_sptdiff,ifield5)
!   INPUT ARGUMENT LIST:
!     bin_scale_fac  - binary scale factor (E)
!     dec_scale_fac  - Decimal scale factor (D)
!     order_of_sptdiff - Order of spatial difference
!
!   OUTPUT ARRAY:
!     ifield5  - GRIB2 - GRID Template 5.3  listing
!
!********************************************************************************
! ifield5(1): reference value(R) (IEEE 32-bit floating-point value)             *
! ifield5(2): binary scale factor (E)                                           *
! ifield5(3): decimal scale factor (D)                                          *
! ifield5(4): number of bits used for each packed value for simple packing      *
!             or for each group reference value for complex packing or          *
!             spatial differencing                                              *
!********************************************************************************
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     integer(4),intent(in) :: dec_scale_fac,bin_scale_fac
     character(*),intent(in) :: order_of_sptdiff
     integer(4),intent(out) :: ifield5(18)
!
     integer(4) :: value,ierr
!
     ifield5=0
     ifield5(1) = 0 ! Any value. Will be later overwritten
     ifield5(2) = bin_scale_fac
     ifield5(3) = dec_scale_fac
!
     call get_g2_ordofspcdiffvals(order_of_sptdiff,value,ierr)
     ifield5(17) = value
!
     end subroutine g2sec5_temp3
!
!
     subroutine g2sec5_temp40(dec_scale_fac,bin_scale_fac,tlnumbits,                   &
                               type_of_compression,ifield5)
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec5_temp40
!   PRGMMR: V. Krishna Kumar         ORG: W/NP12    DATE: 2010-03-01
!
! ABSTRACT: This subroutine returns the section 5 list array for a given decimal
!   scale factor (D),type of original field value (Table 5.40) and type of compression used 
!   from GRIB2 - GRID Template 5.40 Grid point data - JPEG 2000 Code Stream Format  
!
! PROGRAM HISTORY LOG:
! 2010-03-01  V. Krishna Kumar
! 2010-04-07  Jun Wang    add total number of bits
!
! USAGE:    CALL g2sec5_temp40(bin_scale_fac,dec_scale_fac,tlnumbit,type_of_field,type_of_comp,ifield5)
!   INPUT ARGUMENT LIST:
!     dec_scale_fac  - Decimal scale factor (D)
!     bin_scale_fac  - binary scale factor (B)
!     tlnumbits      - total number of bits
!     type_of_field - Type of original field values (see Code Table 5.40)
!     type_of_comp - Type of original field values (see Code Table 5.40)
!
!   OUTPUT ARRAY:
!     ifield5  - GRIB2 - GRID Template 5.40  listing
!
!********************************************************************************
! ifield5(1): reference value(R) (IEEE 32-bit floating-point value)             *
! ifield5(2): binary scale factor (E)                                           *
! ifield5(3): decimal scale factor (D)                                          *
! ifield5(4): number of bits required to hold the resulting scaled and          *
!             reference data values (i.e. The depth of the grayscale image.)    *
!             (see Note 2)                                                      *
! ifield5(5): type of original field values (See Code Table 5.1)                *
! ifield5(6): type of compression used (See Code Table 5.40)                    *
! ifield5(7): target compression ration, M:1 (with respect to the bit-depth     *
!             specified in octet 20), when octet 22 indicates Lossy Compression.* 
!                                         Otherwise, set to missing (see Note 3)*
!********************************************************************************
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
     integer(4),intent(in) :: bin_scale_fac,dec_scale_fac,tlnumbits
     character(*),intent(in) :: type_of_compression
     integer(4),intent(inout) :: ifield5(7)
!
!--- local variable
     integer(4) :: value,ierr
     integer,parameter :: MAX_NUMBIT=16
     integer ibm
     integer,allocatable   :: mg(:)
!
     ifield5(1) = 0 ! Any value. Will be later overwritten
     ifield5(2) = bin_scale_fac
     ifield5(3) = dec_scale_fac
     ifield5(4) = tlnumbits 
     ifield5(5) = 0                  !g2lib assumes original data were reals 
!
!     call get_g2_typeoforigfieldvals(type_of_field,value,ierr)
!     ifield5(5) = value
!
     call get_g2_typeofcompression(type_of_compression,value,ierr)
     ifield5(6) = value
!
     ifield5(7) = 255 
     end subroutine g2sec5_temp40

!=======================================================================
end module grib2_all_tables_module
