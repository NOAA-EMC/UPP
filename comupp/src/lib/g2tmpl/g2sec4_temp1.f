     subroutine g2sec4_temp1(icatg,iparm,typ_gen_proc_key,                         &
                             gen_proc_or_mod_key,hrs_obs_cutoff,min_obs_cutoff,    &
                             unit_of_time_key,fcst_time,lvl_type1,scale_fac1,      &
                             scaled_val1,lvl_type2,scale_fac2,scaled_val2,         &
                             type_ens_fcst_key,perturb_num,num_fcst_ens,           &
                             ipdstmpl1)
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    g2sec4_temp1
!   PRGMMR: Boi Vuong         ORG: W/SIB    DATE: 2015-01-09
!
! ABSTRACT: This subroutine returns the Grib2 Section 4 Template 4.1 list for given keys
!           PDT 4.1 - Individual ensemble forecast, control and perturbed, at a
!                     horizontal level or in a horizontal layer at a point in time.
!
! PROGRAM HISTORY LOG:
! 2015-01-09  Boi Vuong
!
! USAGE:    CALL g2sec4_temp1(icatg,iparm,typ_gen_proc_key,gen_proc_or_mod_key,
!                             hrs_obs_cutoff,min_obs_cutoff,unit_of_time_key,
!                             fcst_time,lvl_type1,scale_fac1,scaled_val1,lvl_type2,
!                             scale_fac2,scaled_val2,type_ens_fcst_key,perturb_num,
!                             num_fcst_ens,ipdstmpl1)
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
!      type_ens_fcst_key - Type of ensemble forecast (see Code table 4.6)
!      perturb_num - Perturbation ensemble number
!      num_fcst_ens - number of forecasts in ensemble
!
!   OUTPUT ARRAY:
!      ipdstmpl1  - GRIB2 PDS Template 4.1 listing
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!

     use grib2_all_tables_module

     integer(4),intent(in) :: icatg,iparm,hrs_obs_cutoff,min_obs_cutoff,       &
                   fcst_time,scale_fac1,scaled_val1,scale_fac2,scaled_val2
     integer(4),intent(in)  :: perturb_num, num_fcst_ens
!
     character(len=*),intent(in) :: typ_gen_proc_key,gen_proc_or_mod_key,      &
                          unit_of_time_key,lvl_type1,lvl_type2,                &
                          type_ens_fcst_key
!
     integer(4),intent(inout)  :: ipdstmpl1(18)
!
!local vars
     integer(4) :: value,ierr
     integer(4) :: bckgnd_gen_proc_id    ! defined by the center
!
     bckgnd_gen_proc_id=0    ! defined by the center
!
     ipdstmpl1(1) = icatg
     ipdstmpl1(2) = iparm
!
     call get_g2_typeofgenproc(typ_gen_proc_key,value,ierr)
     ipdstmpl1(3) = value
!
     ipdstmpl1(4) = bckgnd_gen_proc_id
!
     call get_g2_on388genproc(gen_proc_or_mod_key,value,ierr)
     ipdstmpl1(5) = value
!
     ipdstmpl1(6) = hrs_obs_cutoff
     ipdstmpl1(7) = min_obs_cutoff
!
     call get_g2_unitoftimerange(unit_of_time_key,value,ierr)
     ipdstmpl1(8) = value
     ipdstmpl1(9) = fcst_time
!
     call get_g2_fixedsurfacetypes(lvl_type1,value,ierr)
     ipdstmpl1(10) = value
     ipdstmpl1(11) = scale_fac1
     ipdstmpl1(12) = scaled_val1
!
     call get_g2_fixedsurfacetypes(lvl_type2,value,ierr)
     ipdstmpl1(13) = value
!
     ipdstmpl1(14) = scale_fac2
     ipdstmpl1(15) = scaled_val2
!
     call get_g2_typeofensfcst(type_ens_fcst_key,value,ierr)
     ipdstmpl1(16) = value
!
     ipdstmpl1(17) = perturb_num
     ipdstmpl1(18) = num_fcst_ens
!
     return
     end
