module xml_data_post_t
   use READ_XML_PRIMITIVES
   use XMLPARSE
   implicit none
   integer, private :: lurep_
   logical, private :: strict_

type param_t
   integer                                         :: post_avblfldidx=-9999
   character(len=80)                                :: shortname=''
   character(len=300)                                :: longname=''
   integer                                         :: mass_windpoint=1
   character(len=30)                                :: pdstmpl='tmpl4_0'
   character(len=30)                                :: pname=''
   character(len=10)                                :: table_info=''
   character(len=80)                                :: stats_proc=''
   character(len=80)                                :: fixed_sfc1_type=''
   integer, dimension(:), pointer                  :: scale_fact_fixed_sfc1 => null()
   real, dimension(:), pointer                     :: level => null()
   character(len=80)                                :: fixed_sfc2_type=''
   integer, dimension(:), pointer                  :: scale_fact_fixed_sfc2 => null()
   real, dimension(:), pointer                     :: level2 => null()
   character(len=80)                                :: aerosol_type=''
   character(len=80)                                :: typ_intvl_size=''
   integer                                         :: scale_fact_1st_size=0
   real                                            :: scale_val_1st_size=0.0
   integer                                         :: scale_fact_2nd_size=0
   real                                            :: scale_val_2nd_size=0.0
   character(len=80)                                :: typ_intvl_wvlen=''
   integer                                         :: scale_fact_1st_wvlen=0
   real                                            :: scale_val_1st_wvlen=0.0
   integer                                         :: scale_fact_2nd_wvlen=0
   real                                            :: scale_val_2nd_wvlen=0.0
   real, dimension(:), pointer                     :: scale => null()
   integer                                         :: stat_miss_val=0
   integer                                         :: leng_time_range_prev=0
   integer                                         :: time_inc_betwn_succ_fld=0
   character(len=80)                                :: type_of_time_inc=''
   character(len=20)                                :: stat_unit_time_key_succ=''
   character(len=20)                                :: bit_map_flag=''
end type param_t

type paramset_t
   character(len=6)                                :: datset=''
   integer                                         :: grid_num=255
   character(len=20)                                :: sub_center=''
   character(len=20)                                :: version_no=''
   character(len=20)                                :: local_table_vers_no=''
   character(len=20)                                :: sigreftime=''
   character(len=20)                                :: prod_status=''
   character(len=20)                                :: data_type=''
   character(len=20)                                :: gen_proc_type=''
   character(len=30)                                :: time_range_unit=''
   character(len=50)                                :: orig_center=''
   character(len=30)                                :: gen_proc=''
   character(len=50)                                :: packing_method=''
   character(len=30)                                :: order_of_sptdiff='1st_ord_sptdiff'
   character(len=20)                                :: field_datatype=''
   character(len=30)                                :: comprs_type=''
   type(param_t), dimension(:), pointer            :: param => null()
end type paramset_t

type post_avblfld_t
   type(param_t), dimension(:), pointer            :: param => null()
end type post_avblfld_t
   type(paramset_t), dimension(:), pointer         :: paramset => null()
   type(post_avblfld_t),save                       :: post_avblflds
contains
subroutine read_xml_type_param_t_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(param_t), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(param_t), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_param_t( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_param_t_array

subroutine read_xml_type_param_t( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(param_t), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_post_avblfldidx
   logical                                         :: has_shortname
   logical                                         :: has_longname
   logical                                         :: has_mass_windpoint
   logical                                         :: has_pdstmpl
   logical                                         :: has_pname
   logical                                         :: has_table_info
   logical                                         :: has_stats_proc
   logical                                         :: has_fixed_sfc1_type
   logical                                         :: has_scale_fact_fixed_sfc1
   logical                                         :: has_level
   logical                                         :: has_fixed_sfc2_type
   logical                                         :: has_scale_fact_fixed_sfc2
   logical                                         :: has_level2
   logical                                         :: has_aerosol_type
   logical                                         :: has_typ_intvl_size
   logical                                         :: has_scale_fact_1st_size
   logical                                         :: has_scale_val_1st_size
   logical                                         :: has_scale_fact_2nd_size
   logical                                         :: has_scale_val_2nd_size
   logical                                         :: has_typ_intvl_wvlen
   logical                                         :: has_scale_fact_1st_wvlen
   logical                                         :: has_scale_val_1st_wvlen
   logical                                         :: has_scale_fact_2nd_wvlen
   logical                                         :: has_scale_val_2nd_wvlen
   logical                                         :: has_scale
   logical                                         :: has_stat_miss_val
   logical                                         :: has_leng_time_range_prev
   logical                                         :: has_time_inc_betwn_succ_fld
   logical                                         :: has_type_of_time_inc
   logical                                         :: has_stat_unit_time_key_succ
   logical                                         :: has_bit_map_flag
   has_post_avblfldidx                  = .false.
   has_shortname                        = .false.
   has_longname                         = .false.
   has_mass_windpoint                   = .false.
   has_pdstmpl                          = .false.
   has_pname                            = .false.
   has_table_info                       = .false.
   has_stats_proc                       = .false.
   has_fixed_sfc1_type                  = .false.
   has_scale_fact_fixed_sfc1            = .false.
   allocate(dvar%scale_fact_fixed_sfc1(0))
   has_level                            = .false.
   allocate(dvar%level(0))
   has_fixed_sfc2_type                  = .false.
   has_scale_fact_fixed_sfc2            = .false.
   allocate(dvar%scale_fact_fixed_sfc2(0))
   has_level2                           = .false.
   allocate(dvar%level2(0))
   has_aerosol_type                     = .false.
   has_typ_intvl_size                   = .false.
   has_scale_fact_1st_size              = .false.
   has_scale_val_1st_size               = .false.
   has_scale_fact_2nd_size              = .false.
   has_scale_val_2nd_size               = .false.
   has_typ_intvl_wvlen                  = .false.
   has_scale_fact_1st_wvlen             = .false.
   has_scale_val_1st_wvlen              = .false.
   has_scale_fact_2nd_wvlen             = .false.
   has_scale_val_2nd_wvlen              = .false.
   has_scale                            = .false.
   allocate(dvar%scale(0))
   has_stat_miss_val                    = .false.
   has_leng_time_range_prev             = .false.
   has_time_inc_betwn_succ_fld          = .false.
   has_type_of_time_inc                 = .false.
   has_stat_unit_time_key_succ          = .false.
   has_bit_map_flag                     = .false.
   call init_xml_type_param_t(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata .ne. 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ .lt. noatt_ .and. noatt_ .gt. 1 ) then
         att_      = att_ + 1
         if ( att_ .le. noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag .eq. starttag ) then
         exit
      endif
      if ( endtag .and. noattribs .eq. 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('post_avblfldidx')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%post_avblfldidx, has_post_avblfldidx )
      case('shortname')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%shortname, has_shortname )
      case('longname')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%longname, has_longname )
      case('mass_windpoint')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%mass_windpoint, has_mass_windpoint )
      case('pdstmpl')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%pdstmpl, has_pdstmpl )
      case('pname')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%pname, has_pname )
      case('table_info')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%table_info, has_table_info )
      case('stats_proc')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%stats_proc, has_stats_proc )
      case('fixed_sfc1_type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%fixed_sfc1_type, has_fixed_sfc1_type )
      case('scale_fact_fixed_sfc1')
         call read_xml_integer_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scale_fact_fixed_sfc1, has_scale_fact_fixed_sfc1 )
      case('level')
         call read_xml_real_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%level, has_level )
      case('fixed_sfc2_type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%fixed_sfc2_type, has_fixed_sfc2_type )
      case('scale_fact_fixed_sfc2')
         call read_xml_integer_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scale_fact_fixed_sfc2, has_scale_fact_fixed_sfc2 )
      case('level2')
         call read_xml_real_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%level2, has_level2 )
      case('aerosol_type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%aerosol_type, has_aerosol_type )
      case('typ_intvl_size')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%typ_intvl_size, has_typ_intvl_size )
      case('scale_fact_1st_size')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scale_fact_1st_size, has_scale_fact_1st_size )
      case('scale_val_1st_size')
         call read_xml_real( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scale_val_1st_size, has_scale_val_1st_size )
      case('scale_fact_2nd_size')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scale_fact_2nd_size, has_scale_fact_2nd_size )
      case('scale_val_2nd_size')
         call read_xml_real( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scale_val_2nd_size, has_scale_val_2nd_size )
      case('typ_intvl_wvlen')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%typ_intvl_wvlen, has_typ_intvl_wvlen )
      case('scale_fact_1st_wvlen')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scale_fact_1st_wvlen, has_scale_fact_1st_wvlen )
      case('scale_val_1st_wvlen')
         call read_xml_real( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scale_val_1st_wvlen, has_scale_val_1st_wvlen )
      case('scale_fact_2nd_wvlen')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scale_fact_2nd_wvlen, has_scale_fact_2nd_wvlen )
      case('scale_val_2nd_wvlen')
         call read_xml_real( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scale_val_2nd_wvlen, has_scale_val_2nd_wvlen )
      case('scale')
         call read_xml_real_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scale, has_scale )
      case('stat_miss_val')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%stat_miss_val, has_stat_miss_val )
      case('leng_time_range_prev')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%leng_time_range_prev, has_leng_time_range_prev )
      case('time_inc_betwn_succ_fld')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%time_inc_betwn_succ_fld, has_time_inc_betwn_succ_fld )
      case('type_of_time_inc')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%type_of_time_inc, has_type_of_time_inc )
      case('stat_unit_time_key_succ')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%stat_unit_time_key_succ, has_stat_unit_time_key_succ )
      case('bit_map_flag')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%bit_map_flag, has_bit_map_flag )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_post_avblfldidx ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on post_avblfldidx')
   endif
   if ( .not. has_shortname ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on shortname')
   endif
   if ( .not. has_longname ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on longname')
   endif
   if ( .not. has_mass_windpoint ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on mass_windpoint')
   endif
   if ( .not. has_pdstmpl ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on pdstmpl')
   endif
   if ( .not. has_pname ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on pname')
   endif
   if ( .not. has_table_info ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on table_info')
   endif
   if ( .not. has_stats_proc ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on stats_proc')
   endif
   if ( .not. has_fixed_sfc1_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on fixed_sfc1_type')
   endif
   if ( .not. has_scale_fact_fixed_sfc1 ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scale_fact_fixed_sfc1')
   endif
   if ( .not. has_fixed_sfc2_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on fixed_sfc2_type')
   endif
   if ( .not. has_scale_fact_fixed_sfc2 ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scale_fact_fixed_sfc2')
   endif
   if ( .not. has_level2 ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on level2')
   endif
   if ( .not. has_aerosol_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on aerosol_type')
   endif
   if ( .not. has_typ_intvl_size ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on typ_intvl_size')
   endif
   if ( .not. has_scale_fact_1st_size ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scale_fact_1st_size')
   endif
   if ( .not. has_scale_val_1st_size ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scale_val_1st_size')
   endif
   if ( .not. has_scale_fact_2nd_size ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scale_fact_2nd_size')
   endif
   if ( .not. has_scale_val_2nd_size ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scale_val_2nd_size')
   endif
   if ( .not. has_typ_intvl_wvlen ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on typ_intvl_wvlen')
   endif
   if ( .not. has_scale_fact_1st_wvlen ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scale_fact_1st_wvlen')
   endif
   if ( .not. has_scale_val_1st_wvlen ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scale_val_1st_wvlen')
   endif
   if ( .not. has_scale_fact_2nd_wvlen ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scale_fact_2nd_wvlen')
   endif
   if ( .not. has_scale_val_2nd_wvlen ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scale_val_2nd_wvlen')
   endif
   if ( .not. has_scale ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scale')
   endif
   if ( .not. has_stat_miss_val ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on stat_miss_val')
   endif
   if ( .not. has_leng_time_range_prev ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on leng_time_range_prev')
   endif
   if ( .not. has_time_inc_betwn_succ_fld ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on time_inc_betwn_succ_fld')
   endif
   if ( .not. has_type_of_time_inc ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on type_of_time_inc')
   endif
   if ( .not. has_stat_unit_time_key_succ ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on stat_unit_time_key_succ')
   endif
   if ( .not. has_bit_map_flag ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on bit_map_flag')
   endif
end subroutine read_xml_type_param_t
subroutine init_xml_type_param_t_array( dvar )
   type(param_t), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_param_t_array
subroutine init_xml_type_param_t(dvar)
   type(param_t) :: dvar
   dvar%level = 0.
end subroutine init_xml_type_param_t
subroutine read_xml_type_paramset_t_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(paramset_t), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(paramset_t), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_paramset_t( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_paramset_t_array

subroutine read_xml_type_paramset_t( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(paramset_t), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_datset
   logical                                         :: has_grid_num
   logical                                         :: has_sub_center
   logical                                         :: has_version_no
   logical                                         :: has_local_table_vers_no
   logical                                         :: has_sigreftime
   logical                                         :: has_prod_status
   logical                                         :: has_data_type
   logical                                         :: has_gen_proc_type
   logical                                         :: has_time_range_unit
   logical                                         :: has_orig_center
   logical                                         :: has_gen_proc
   logical                                         :: has_packing_method
   logical                                         :: has_order_of_sptdiff
   logical                                         :: has_field_datatype
   logical                                         :: has_comprs_type
   logical                                         :: has_param
   has_datset                           = .false.
   has_grid_num                         = .false.
   has_sub_center                       = .false.
   has_version_no                       = .false.
   has_local_table_vers_no              = .false.
   has_sigreftime                       = .false.
   has_prod_status                      = .false.
   has_data_type                        = .false.
   has_gen_proc_type                    = .false.
   has_time_range_unit                  = .false.
   has_orig_center                      = .false.
   has_gen_proc                         = .false.
   has_packing_method                   = .false.
   has_order_of_sptdiff                 = .false.
   has_field_datatype                   = .false.
   has_comprs_type                      = .false.
   has_param                            = .false.
   allocate(dvar%param(0))
   call init_xml_type_paramset_t(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata .ne. 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ .lt. noatt_ .and. noatt_ .gt. 1 ) then
         att_      = att_ + 1
         if ( att_ .le. noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag .eq. starttag ) then
         exit
      endif
      if ( endtag .and. noattribs .eq. 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('datset')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%datset, has_datset )
      case('grid_num')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%grid_num, has_grid_num )
      case('sub_center')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%sub_center, has_sub_center )
      case('version_no')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%version_no, has_version_no )
      case('local_table_vers_no')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%local_table_vers_no, has_local_table_vers_no )
      case('sigreftime')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%sigreftime, has_sigreftime )
      case('prod_status')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%prod_status, has_prod_status )
      case('data_type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%data_type, has_data_type )
      case('gen_proc_type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%gen_proc_type, has_gen_proc_type )
      case('time_range_unit')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%time_range_unit, has_time_range_unit )
      case('orig_center')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%orig_center, has_orig_center )
      case('gen_proc')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%gen_proc, has_gen_proc )
      case('packing_method')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%packing_method, has_packing_method )
      case('order_of_sptdiff')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%order_of_sptdiff, has_order_of_sptdiff )
      case('field_datatype')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%field_datatype, has_field_datatype )
      case('comprs_type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%comprs_type, has_comprs_type )
      case('param')
         call read_xml_type_param_t_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%param, has_param )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_datset ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on datset')
   endif
   if ( .not. has_grid_num ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on grid_num')
   endif
   if ( .not. has_sub_center ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on sub_center')
   endif
   if ( .not. has_version_no ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on version_no')
   endif
   if ( .not. has_local_table_vers_no ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on local_table_vers_no')
   endif
   if ( .not. has_sigreftime ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on sigreftime')
   endif
   if ( .not. has_prod_status ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on prod_status')
   endif
   if ( .not. has_data_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on data_type')
   endif
   if ( .not. has_gen_proc_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on gen_proc_type')
   endif
   if ( .not. has_time_range_unit ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on time_range_unit')
   endif
   if ( .not. has_orig_center ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on orig_center')
   endif
   if ( .not. has_gen_proc ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on gen_proc')
   endif
   if ( .not. has_packing_method ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on packing_method')
   endif
   if ( .not. has_order_of_sptdiff ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on order_of_sptdiff')
   endif
   if ( .not. has_field_datatype ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on field_datatype')
   endif
   if ( .not. has_comprs_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on comprs_type')
   endif
   if ( .not. has_param ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on param')
   endif
end subroutine read_xml_type_paramset_t
subroutine init_xml_type_paramset_t_array( dvar )
   type(paramset_t), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_paramset_t_array
subroutine init_xml_type_paramset_t(dvar)
   type(paramset_t) :: dvar
end subroutine init_xml_type_paramset_t
subroutine read_xml_type_post_avblfld_t_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(post_avblfld_t), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(post_avblfld_t), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_post_avblfld_t( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_post_avblfld_t_array

subroutine read_xml_type_post_avblfld_t( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(post_avblfld_t), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_param
   has_param                            = .false.
   allocate(dvar%param(0))
   call init_xml_type_post_avblfld_t(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata .ne. 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ .lt. noatt_ .and. noatt_ .gt. 1 ) then
         att_      = att_ + 1
         if ( att_ .le. noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag .eq. starttag ) then
         exit
      endif
      if ( endtag .and. noattribs .eq. 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('param')
         call read_xml_type_param_t_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%param, has_param )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_param ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on param')
   endif
end subroutine read_xml_type_post_avblfld_t
subroutine init_xml_type_post_avblfld_t_array( dvar )
   type(post_avblfld_t), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_post_avblfld_t_array
subroutine init_xml_type_post_avblfld_t(dvar)
   type(post_avblfld_t) :: dvar
end subroutine init_xml_type_post_avblfld_t
subroutine read_xml_file_post_t(fname, lurep, errout)
   character(len=*), intent(in)           :: fname
   integer, intent(in), optional          :: lurep
   logical, intent(out), optional         :: errout

   type(XML_PARSE)                        :: info
   logical                                :: error
   character(len=80)                      :: tag
   character(len=80)                      :: starttag
   logical                                :: endtag
   character(len=80), dimension(1:2,1:20) :: attribs
   integer                                :: noattribs
   character(len=200), dimension(1:100)   :: data
   integer                                :: nodata
   logical                                         :: has_paramset
   logical                                         :: has_post_avblflds
   has_paramset                         = .false.
   allocate(paramset(0))
   has_post_avblflds                    = .false.

   call init_xml_file_post_t
   call xml_open( info, fname, .true. )
   call xml_options( info, report_errors=.true., ignore_whitespace=.true.)
   lurep_ = 0
   if ( present(lurep) ) then
      lurep_ = lurep
      call xml_options( info, report_lun=lurep )
   endif
   do
      call xml_get( info, starttag, endtag, attribs, noattribs, &
         data, nodata)
      if ( starttag .ne. '!--' ) exit
   enddo
   if ( starttag .ne. "postxml" ) then
      call xml_report_errors( info, &
         'XML-file should have root element "postxml"')
      error = .true.
      call xml_close(info)
      return
   endif
   strict_ = .true.
   error = .false.
   do
      call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
      if ( xml_error(info) ) then
         write(lurep_,*) 'Error reading input file!'
         error = .true.
         return
      endif
      if ( endtag .and. tag .eq. starttag ) then
         exit
      endif
      if ( endtag .and. noattribs .eq. 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('paramset')
         call read_xml_type_paramset_t_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            paramset, has_paramset )
      case('post_avblflds')
         call read_xml_type_post_avblfld_t( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            post_avblflds, has_post_avblflds )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_paramset ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on paramset')
   endif
   if ( .not. has_post_avblflds ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on post_avblflds')
   endif
   if ( present(errout) ) errout = error
end subroutine
subroutine init_xml_file_post_t

end subroutine

end module
