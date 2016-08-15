     subroutine g2sec4_temp11(icatg,iparm,typ_gen_proc_key,gen_proc_or_mod_key,    &
                             hrs_obs_cutoff,min_obs_cutoff,unit_of_time_key,       &
                             fcst_time,lvl_type1,scale_fac1,scaled_val1,lvl_type2, &
                             scale_fac2,scaled_val2,type_ens_fcst_key,             &
                             perturb_num,num_fcst_ens,year_intvl,                  &
                             mon_intvl,day_intvl,hour_intvl,min_intvl,sec_intvl,   &
                             num_time_range,stat_miss_val,type_of_stat_proc,       &
                             type_of_time_inc,stat_unit_time_key,                  &
                             leng_time_range_stat,stat_unit_time_key_succ,         &
                             time_inc_betwn_succ_fld,ipdstmpl11)
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec4_temp11
!   PRGMMR: Boi Vuong         ORG: W/SIB    DATE: 2015-01-09
!
! ABSTRACT: This subroutine returns the Grib2 Section 4 Template 4.11 list for given keys
!           PDT 4.11 - Individual ensemble forecast, control and perturbed, at a
!                      horizontal level or in a horizontal layer, in a continuous
!                      or non-continuous time interval
!
! PROGRAM HISTORY LOG:
! 2015-01-09  Boi Vuong
!
! USAGE:    CALL g2sec4_temp11(icatg,iparm,typ_gen_proc_key,gen_proc_or_mod_key,    &
!                             hrs_obs_cutoff,min_obs_cutoff,unit_of_time_key,       &
!                             fcst_time,lvl_type1,scale_fac1,scaled_val1,lvl_type2, &
!                             scale_fac2,scaled_val2,type_ens_fcst_key,             &
!                             perturb_num,num_fcst_ens,year_intvl,                  &
!                             mon_intvl,day_intvl,hour_intvl,min_intvl,sec_intvl,   &
!                             num_time_range,stat_miss_val,type_of_stat_proc,       &
!                             type_of_time_inc,stat_unit_time_key,                  &
!                             leng_time_range_stat,stat_unit_time_key_succ,         &
!                             time_inc_betwn_succ_fld,ipdstmpl11)
!   INPUT ARGUMENT LIST:
!
!       icatg - Parameter category (see Code Table 4.1)
!       iparm - Parameter number (see Code Table 4.2)
!       typ_gen_proc_key - Type of generating process (see Code Table 4.3)
!       bckgnd_gen_proc_id - Background generating process identifier (defined by originating centre)
!       gen_proc_or_mod_key - Analysis or forecast generating process identified (see Code ON388 Table A)
!       hrs_obs_cutoff - Hours after reference time data cutoff (see Note 1)
!       min_obs_cutoff - Minutes after reference time data cutoff
!       unit_of_time_key - Indicator of unit of time range (see Code Table 4.4)
!       fcst_time - Forecast time in units defined by octet 18 (see Note 2)
!       lvl_type1 - Type of first fixed surface (see Code Table 4.5)
!       scale_fac1 - Scale factor of first fixed surface
!       scaled_val1 - Scaled value of first fixed surface
!       lvl_type2 - Type of second fixed surfaced (see Code Table 4.5)
!       scale_fac2 - Scale factor of second fixed surface
!       scaled_val2 - Scaled value of second fixed surfaces
!       type_ens_fcst_key - Type of ensemble forecast (see Code table 4.6)
!       perturb_num - Perturbation ensemble number
!       num_fcst_ens - number of forecasts in ensemble
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
!      ipdstmpl11  - GRIB2 PDS Template 4.11 listing
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!

     use grib2_all_tables_module

     integer(4),intent(in) :: icatg,iparm,hrs_obs_cutoff,min_obs_cutoff,fcst_time, &
                              scale_fac1,scaled_val1,scale_fac2,scaled_val2
     integer(4),intent(in) :: year_intvl,mon_intvl,day_intvl,hour_intvl,min_intvl, &
                              sec_intvl,num_time_range,stat_miss_val, &
                              leng_time_range_stat,time_inc_betwn_succ_fld
     integer(4),intent(in)  :: perturb_num, num_fcst_ens
!
     character(len=*),intent(in) :: typ_gen_proc_key,gen_proc_or_mod_key, &
                          unit_of_time_key,lvl_type1,lvl_type2,           &
                          type_of_stat_proc,type_of_time_inc,             &
                          stat_unit_time_key,stat_unit_time_key_succ,     &
                          type_ens_fcst_key
!
     integer(4)               :: bckgnd_gen_proc_id    ! defined by the center
!
     integer(4),intent(inout) :: ipdstmpl11(32)        ! currently works only for n=1
                                                       ! later on, this will be generalized
!
!-- local vars
     integer(4) :: value,ierr
!
     bckgnd_gen_proc_id=0
!
     ipdstmpl11(1) = icatg
     ipdstmpl11(2) = iparm
!
     call get_g2_typeofgenproc(typ_gen_proc_key,value,ierr)
     ipdstmpl11(3) = value
!
     ipdstmpl11(4) = bckgnd_gen_proc_id
!
     call get_g2_on388genproc(gen_proc_or_mod_key,value,ierr)
     ipdstmpl11(5) = value
!
     ipdstmpl11(6) = hrs_obs_cutoff
     ipdstmpl11(7) = min_obs_cutoff
!
     call get_g2_unitoftimerange(unit_of_time_key,value,ierr)
     ipdstmpl11(8) = value
     ipdstmpl11(9) = fcst_time
!
     call get_g2_fixedsurfacetypes(lvl_type1,value,ierr)
     ipdstmpl11(10) = value
     ipdstmpl11(11) = scale_fac1
     ipdstmpl11(12) = scaled_val1
!
     call get_g2_fixedsurfacetypes(lvl_type2,value,ierr)
     ipdstmpl11(13) = value
!
     ipdstmpl11(14) = scale_fac2
     ipdstmpl11(15) = scaled_val2
!
     call get_g2_typeofensfcst(type_ens_fcst_key,value,ierr)
     ipdstmpl11(16) = value
!
     ipdstmpl11(17) = perturb_num
     ipdstmpl11(18) = num_fcst_ens
!
     ipdstmpl11(19) = year_intvl
     ipdstmpl11(20) = mon_intvl
     ipdstmpl11(21) = day_intvl
     ipdstmpl11(22) = hour_intvl
     ipdstmpl11(23) = min_intvl
     ipdstmpl11(24) = sec_intvl
!
     ipdstmpl11(25) = num_time_range ! choose n=1 for this case
     ipdstmpl11(26) = stat_miss_val  ! choose 0 for this case
!
     call get_g2_statprocesstypes(type_of_stat_proc,value,ierr)
     ipdstmpl11(27) = value  ! types_of_stat_proc='accumulation'
!
     call get_g2_typeoftimeintervals(type_of_time_inc,value,ierr)
     ipdstmpl11(28) = value  ! type_of_time_inc='same_start_time_fcst_fcst_time_inc'
                             ! value = 2 (Successive times processed have same start
                             !       time of forecast, forecast time is incremented)
!
     call get_g2_unitoftimerange(stat_unit_time_key,value,ierr)
     ipdstmpl11(29) = value  ! stat_unit_time_key='hour'
                             ! value = 1
     ipdstmpl11(30) = leng_time_range_stat  ! value = 6
!
     call get_g2_unitoftimerange(stat_unit_time_key_succ,value,ierr)
                             ! stat_unit_time_key_succ='missing'
     ipdstmpl11(31) = value  ! value = 255
!
     ipdstmpl11(32) = time_inc_betwn_succ_fld   ! value = 0
!
     return
     end
