#ifndef USE_IFI

! IFI stubs

module upp_ifi_mod
  implicit none
  logical, parameter :: have_libIFI = .false.
  
contains !==============================================================

  subroutine set_ifi_dims()
    use CTLBLK_mod, only: ifi_nflight, ifi_flight_levels
    implicit none

    ! Bogus fill value for flight levels to prevent a crash
    ifi_nflight = 0
    allocate(ifi_flight_levels(1))
    ifi_flight_levels = 10000
  end subroutine set_ifi_dims

  ! --------------------------------------------------------------------

  subroutine run_ifi()
    write(0,*) 'WARNING: This post cannot produce IFI icing products because it was not compiled with IFI.'
  end subroutine run_ifi
end module upp_ifi_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Actual IFI code.
module upp_ifi_mod
  use ifi_type_mod
  use ifi_mod
  use iso_c_binding, only: c_bool, c_int64_t, c_ptr, c_double
  implicit none
  private

  public run_ifi, set_ifi_dims, have_libIFI

  logical, parameter :: have_libIFI = .true.
  real, parameter :: feet2meters = 0.3048

  type(IFIConfig) :: ifi_config
  logical :: have_read_ifi_config = .false.

contains !==============================================================

  subroutine ifi_check(status,error_message)
    use mpi, only: MPI_Abort, MPI_COMM_WORLD
    implicit none
    integer(c_int64_t), intent(in) :: status
    character(*), intent(in) :: error_message
    integer :: ierr

    ! Exit program if status is non-zero
    if(status/=0) then
      write(0,'("IFI Failed: ",A)') trim(error_message)
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine ifi_check

  ! --------------------------------------------------------------------

  subroutine read_ifi_config()
    implicit none
    ! Only read the config if we haven't already:
    if(.not.have_read_ifi_config) then
      call ifi_check(ifi_config%init(),'read IFI config')
      have_read_ifi_config = .true.
    endif
  end subroutine read_ifi_config

  ! --------------------------------------------------------------------

  subroutine set_ifi_dims()
    use CTLBLK_mod, only: ifi_nflight, ifi_flight_levels
    implicit none

    ! Warning: do not deallocate config_flight_levels_feet.
    ! The IFI C++ library manages that array.
    integer(kind=c_int64_t), pointer :: config_flight_levels_feet(:)

    ! Make sure the config was read in
    call read_ifi_config()

    ! Get the flight levels
    call ifi_check(ifi_config%get_flight_levels_feet(config_flight_levels_feet), &
         'cannot get flight levels in feet from IFI config')

    ! Convert from integer to real:
    ifi_nflight = size(config_flight_levels_feet)
    allocate(ifi_flight_levels(ifi_nflight))
    ifi_flight_levels = config_flight_levels_feet
  end subroutine set_ifi_dims

  ! --------------------------------------------------------------------

  subroutine send_data(vars,name,ient)
    use ctlblk_mod, only: spval, jsta, jend, lm, im, cfld, datapd, fld_info, ifi_flight_levels, jsta_2l, jend_2u
    use rqstfld_mod, only: iget, iavblfld, lvlsxml, lvls
    implicit none

    class(IFIData) :: vars
    integer, intent(in) :: ient
    character(*), intent(in) :: name

    ! Locals
    integer(c_int64_t) :: &
         ims,ime,jms,jme,kms,kme, ids,ide,jds,jde,kds,kde, ips,ipe,jps,jpe,kps,kpe
    real(kind=ifi_real_t), pointer :: data(:)
    logical(c_bool) :: missing_value_is_set
    real(kind=ifi_real_t) :: missing_value
    integer :: i,j,k,nj,jpad,ilen,jlen

    ! WARNING: do not deallocate the data pointer. It is managed by the IFI C++ library.
    data => vars%get_data(trim(name),missing_value_is_set,missing_value, &
           ims,ime,jms,jme,kms,kme, ids,ide,jds,jde,kds,kde, ips,ipe,jps,jpe,kps,kpe)
    if(.not.missing_value_is_set) then
      missing_value = MISSING
    endif

    ! Get dimensions and do some sanity checks
    nj = jpe-jps+1
    jpad = jps-jms
    if(nj/=jend-jsta+1 .or. jpad/=jsta-jsta_2l) then
      call ifi_check(-1,'Internal error: IFI output j bounds do not match UPP j bounds')
    end if
    ilen=ide-ids+1
    jlen=jde-jds+1
    
    ! Go level-by-level writing grib2 if requested
    do k=1,lm
      if(LVLS(k,IGET(ient))>0) then
        cfld = cfld+1
        fld_info(cfld)%ifld = IAVBLFLD(IGET(ient))
        fld_info(cfld)%lvl  = ifi_flight_levels(k)*feet2meters
        !$OMP PARALLEL DO PRIVATE(i,j,k) COLLAPSE(2)
        do j=1,nj
          do i=1,im
            datapd(i,j,cfld) = data((k*jlen+j+jpad)*ilen+i)
            if(datapd(i,j,cfld)==missing_value) then
              datapd(i,j,cfld) = spval
            endif
          enddo
        enddo
      endif
    enddo
  end subroutine send_data
  
  ! --------------------------------------------------------------------

  subroutine run_ifi()
    use ctlblk_mod, only: spval, lm, lp1, im, jsta_2l,jend_2u, ITPREC, IFHR, IFMIN, grib, jm, jsta,jend
    use vrbls3d, only: &
         zmid, & ! = IFI "HGT"
         zint, & ! zint(:,:,LM+1) = IFI "HGT_surface"
         QQI, & ! = IFI "CIMIXR"
         QQG, & ! = IFI "GRLE"
         QQW, & ! = IFI "CLMR"
         pmid, & ! = IFI "PRES"
         QQR, & ! "RWMR"
         q, & ! "SPFH"
         t, & ! "TMP"
         QQS, & ! "SNMR"
         wh ! "VVEL"
    use vrbls2d, only : &
         ACPCP, & ! IFI "APCP_surface" over ITPREC bucket time
         CAPE, CIN
    use rqstfld_mod, only: IGET
    use iso_c_binding, only: c_bool, c_int64_t
    implicit none

    type(IFIData) :: hybr_vars, cat_vars
    type(IFIConfig) :: config
    type(IFIAlgo) :: algo

    real(c_double) :: fhr
    
    if(grib/='grib2' .or. (IGET(1100)<=0 .and. IGET(1101)<=0 .and. IGET(1102)<=0)) then
      return ! nothing to do
    endif

    print '(A)','Running libIFI to get icing products.'

    ! Read config and initialize input data structures:

    call read_ifi_config()
    call ifi_check(hybr_vars%init(int(1,c_int64_t),int(im,c_int64_t),int(1,c_int64_t),int(jm,c_int64_t),&
         int(1,c_int64_t),int(lm,c_int64_t),int( 1,c_int64_t),int(im,c_int64_t),int(jsta,c_int64_t),&
         int(jend,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t)),&
         'could not initialize IFI input data structures')

    ! Copy 2D vars to IFI internal structures:
    
    call ifi_check(hybr_vars%add_ij_var('topography',int(1,c_int64_t),int(im,c_int64_t),&
         int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t), &
         zint(:,:,LP1)),'could not send topography to IFI')
    if(ITPREC>1e-5) then
      call ifi_check(hybr_vars%add_ij_var('ACPCP',int(1,c_int64_t),int(im,c_int64_t),&
      int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t), &
           ACPCP/ITPREC),'could not send ACPCP to IFI')
    else
      call ifi_check(hybr_vars%add_ij_var('ACPCP',int(1,c_int64_t),int(im,c_int64_t),&
           int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t), &
           ACPCP),'could not send ACPCP to IFI')
    endif
    call ifi_check(hybr_vars%add_ij_var('CAPE',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),CAPE), &
         'could not send CAPE to IFI')
    call ifi_check(hybr_vars%add_ij_var('CIN',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),CAPE), &
         'could not send CIN to IFI')

    ! Copy 3D vars to IFI internal structures:
    
    call ifi_check(hybr_vars%add_ijk_var('HGT',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t),zmid), &
         'could not send HGT (zmid) to IFI')
    call ifi_check(hybr_vars%add_ijk_var('CIMIXR',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t),QQI), &
         'could not send CIMIXR (QQI) to IFI')
    call ifi_check(hybr_vars%add_ijk_var('CLMR',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t),QQW), &
         'could not send CLWR (QQW) to IFI')
    call ifi_check(hybr_vars%add_ijk_var('GRLE',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t),QQG), &
         'could not send GRLE (QQG) to IFI')
    call ifi_check(hybr_vars%add_ijk_var('PRES',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t),pmid), &
         'could not send PRES (pmid) to IFI')
    call ifi_check(hybr_vars%add_ijk_var('RWMR',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t),QQR), &
         'could not send RWMR (QQR) to IFI')
    call ifi_check(hybr_vars%add_ijk_var('SPFH',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t),Q), &
         'could not send SPFH (Q) to IFI')
    call ifi_check(hybr_vars%add_ijk_var('TMP',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t),T), &
         'could not send TMP (T) to IFI')
    call ifi_check(hybr_vars%add_ijk_var('SNMR',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t),QQS), &
         'could not send SNMR (QQS) to IFI')
    call ifi_check(hybr_vars%add_ijk_var('VVEL',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t),wh), &
         'could not send VVEL (wh) to IFI')

    ! Initialize the IFI algorithm

    fhr=IFHR+IFMIN/60.0
    call ifi_check(algo%init(config,fhr,hybr_vars), &
         'could not initialize IFI algorithm')

    ! Change n>=spval to Config.MISSING in all input fields
    
    call ifi_check(algo%cap_to_missing(spval),&
         'could not change post spval to Config.MISSING in input to IFI')

    ! Run the IFI algorithm
    
    call ifi_check(algo%calc_exner(),'calc_exner() failed')
    call ifi_check(algo%hybrid_to_pressure(),'hybrid_to_pressure() failed')
    call ifi_check(algo%discard_hybrid_level_vars(),'discard_hybrid_level_vars() failed')
    call ifi_check(algo%smooth_fields(),'smooth_fields() failed')
    call ifi_check(algo%imprecise_pressure_level_vars(),'imprecise_pressure_level_vars() failed')
    call ifi_check(algo%derive_fields(),'derive_fields() failed')
    call ifi_check(algo%imprecise_derived_vars(),'imprecise_derived_vars() failed')
    call ifi_check(algo%run_fip_algo(),'run_fip_algo() failed')
    call ifi_check(algo%imprecise_fip_algo_vars(),'imprecise_fip_algo_vars() failed')
    call ifi_check(algo%discard_pressure_level_vars(),'discard_pressure_level_vars() failed')
    call ifi_check(algo%discard_derived_vars(),'discard_derived_vars() failed')
    call ifi_check(algo%pressure_to_flight(),'pressure_to_flight() failed')
    call ifi_check(algo%discard_fip_algo_vars(),'discard_fip_algo_vars() failed')
    call ifi_check(algo%make_icing_category(),'make_icing_category() failed')

    ! Get the final output fields:

    call ifi_check(algo%get_cat_vars(cat_vars),'get_cat_vars()')
    call send_data(cat_vars,'ICE_PROB',1100)
    call send_data(cat_vars,'SLD',1101)
    call send_data(cat_vars,'ICE_SEV_CAT',1102)
    
    ! When this subroutine ends, a Fortran-2003-compliant compiler
    ! will free all memory IFI uses, by calling the destructors (final
    ! routines) for cat_vars, hybr_vars, algo, and config.

  end subroutine run_ifi

end module upp_ifi_mod

#endif
