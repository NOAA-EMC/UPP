#ifndef USE_IFI

! IFI stubs

module upp_ifi_mod
  use iso_c_binding, only: c_float
  implicit none

  private
  public run_ifi, set_ifi_dims, ifi_real_t

  integer, parameter :: ifi_real_t = c_float
  
contains !==============================================================

  subroutine set_ifi_dims()
    use CTLBLK_mod, only: ifi_nflight, ifi_flight_levels
    implicit none
    integer :: i

    ! Bogus fill value for flight levels to prevent a crash
    ifi_nflight = 60
    allocate(ifi_flight_levels(ifi_nflight))
    do i=1,ifi_nflight
      ifi_flight_levels(i) = 500*i
    enddo
  end subroutine set_ifi_dims

  ! --------------------------------------------------------------------

  subroutine run_ifi()
    ! Fill any requested IFI fields with missing data.
    call send_missing_data(1100)  ! ICE_PROB = missing
    call send_missing_data(1101)  ! SLD = missing
    call send_missing_data(1102)  ! ICE_SEV_CAT = missing
  end subroutine run_ifi

  ! --------------------------------------------------------------------

  subroutine send_missing_data(ient)
    use CTLBLK_mod, only: ifi_nflight, ifi_flight_levels
    use ctlblk_mod, only: spval, jsta, jend, lm, im, cfld, datapd, fld_info, ifi_flight_levels, jsta_2l, jend_2u
    use rqstfld_mod, only: iget, iavblfld, lvlsxml, lvls
    implicit none

    integer, intent(in) :: ient

    ! Locals
    integer :: i,j,k, ista,iend

    logical, save :: wrote_message = .false. ! guarded by an OMP CRITICAL below

    ! FIXME: Once the post is i-decomposed, replace these.
    ista=1
    iend=im

    ! Go level-by-level writing grib2 if requested
    do k=1,ifi_nflight
      if(LVLS(k,IGET(ient))>0) then

        if(.not.wrote_message) then
          !$OMP CRITICAL
          if(.not.wrote_message) then
            write(0,'(A)') 'This post cannot produce IFI icing products because it was not compiled with libIFI.'
            wrote_message = .true.
          endif
          !$OMP END CRITICAL
        endif

        cfld = cfld+1
        fld_info(cfld)%ifld = IAVBLFLD(IGET(ient))
        fld_info(cfld)%lvl  = k

        !$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)
        do j=jsta,jend
          do i=ista,iend
            datapd(i-ista+1,j-jsta+1,cfld) = spval
          enddo
        enddo
      endif
    enddo
  end subroutine send_missing_data
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
  public run_ifi, set_ifi_dims, ifi_real_t

  real, parameter :: feet2meters = 0.3048

  type(IFIConfig) :: ifi_config
  logical :: have_read_ifi_config = .false.
  logical, parameter :: write_ifi_netcdf_debug_files = .true.
  
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

  subroutine ifi_minmax_callback(minv,maxv) bind(C)
    use ctlblk_mod, only: mpi_comm_comp, me
    use mpi
    implicit none
    real(kind=ifi_real_t) :: minv,maxv
    real(kind=8) :: sendbuf(1),recvbuf(1)
    integer :: ierr

    sendbuf(1)=minv
    call MPI_Allreduce(sendbuf,recvbuf,1,MPI_REAL8,MPI_MIN,mpi_comm_comp,ierr)
    minv=recvbuf(1)

    sendbuf(1)=maxv
    call MPI_Allreduce(sendbuf,recvbuf,1,MPI_REAL8,MPI_MAX,mpi_comm_comp,ierr)
    maxv=recvbuf(1)
  end subroutine ifi_minmax_callback

  ! --------------------------------------------------------------------

  subroutine ifi_halo_callback(d,ims1,ime1,jms1,jme1, &
                                 ids1,ide1,jds1,jde1, &
                                 ips1,ipe1,jps1,jpe1) BIND(C)
    use ctlblk_mod, only: jsta,jend,im,jsta_2l,jend_2u
    implicit none
    real(kind=ifi_real_t) :: d(*)
    integer(kind=c_int64_t), value :: &
         ims1,ime1,jms1,jme1, &
         ids1,ide1,jds1,jde1, &
         ips1,ipe1,jps1,jpe1
    integer(kind=c_int64_t) :: &
         ims,ime,jms,jme, &
         ids,ide,jds,jde, &
         ips,ipe,jps,jpe

    ! Convert from C++ to Fortran indexing:
    ims=ims1+1 ; ime=ime1+1 ; jms=jms1+1 ; jme=jme1+1
    ids=ids1+1 ; ide=ide1+1 ; jds=jds1+1 ; jde=jde1+1
    ips=ips1+1 ; ipe=ipe1+1 ; jps=jps1+1 ; jpe=jpe1+1

    if(jms/=jsta_2l .or. jme/=jend_2u .or. jsta/=jps .or. jend/=jpe &
         .or. ips/=1 .or. ipe/=im .or. ims/=1 .or. ime/=im) then
94    format('   ',A,' = ',I0,', ',I0)
95    format('   ',A,' = ',I0,', ',I0,', ',I0)
      print '(A)','Warning: IFI halo callback j bounds do not match UPP j bounds'
      print 94,'jps,jsta',jps,jsta
      print 94,'jpe,jend',jpe,jend
      print 94,'jms,jsta_2l',jms,jsta_2l
      print 94,'jme,jend_2u',jme,jend_2u
      print 95,'ips,ims,1',ips,ims,1
      print 95,'ipe,ime,im',ipe,ime,im
      !call ifi_check(-1,'Internal error: IFI halo callback j bounds do not match UPP j bounds')
    end if
    call EXCH_IFI(d)
  end subroutine ifi_halo_callback

  ! --------------------------------------------------------------------

  subroutine set_ifi_dims()
    use CTLBLK_mod, only: ifi_nflight, ifi_flight_levels
    implicit none

    ! Warning: do not deallocate config_flight_levels_feet.
    ! The IFI C++ library manages that array.
    integer(kind=c_int64_t), pointer :: config_flight_levels_feet(:)
    integer :: i

    ! Make sure the config was read in
    call read_ifi_config()

    ! Get the flight levels
    call ifi_check(ifi_config%get_flight_levels_feet(config_flight_levels_feet), &
         'cannot get flight levels in feet from IFI config')

    ! Convert from integer to real:
    ifi_nflight = size(config_flight_levels_feet)
    allocate(ifi_flight_levels(ifi_nflight))
    ifi_flight_levels = config_flight_levels_feet

!     print '(A)','IFI flight levels:'
! 38  format('  ifi_flight_level[',I0,'] = ',F15.5)
!     do i=1,ifi_nflight
!       print 38,i,ifi_flight_levels(i)
!     enddo

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
    integer :: i,j,k,nj,jpad,ilen,jlen,kstartm1,jstartm1,iloc,ndata

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
94    format('   ',A,' = ',I0,', ',I0)
83    format('Warning: ',A,': IFI output j bounds do not match UPP j bounds')
      print 83,trim(name)
      print 94,'jps,jsta',jps,jsta
      print 94,'jpe,jend',jpe,jend
      print 94,'jms,jsta_2l',jms,jsta_2l
      print 94,'jme,jend_2u',jme,jend_2u
      print 94,'jpad,jsta-jsta_2l',jpad,jsta-jsta_2l
      print 94,'nj,jend-jsta+1',nj,jend-jsta+1
      !call ifi_check(-1,'Internal error: IFI output j bounds do not match UPP j bounds')
    end if
    ilen=ime-ims+1
    jlen=jme-jms+1
    ndata=ilen*jlen*(kme-kms+1)
    
    ! Go level-by-level writing grib2 if requested
    do k=kps,kpe
      kstartm1=ilen*jlen*(k-kms)
      if(LVLS(k,IGET(ient))>0) then
        cfld = cfld+1
        fld_info(cfld)%ifld = IAVBLFLD(IGET(ient))
        fld_info(cfld)%lvl  = k ! ifi_flight_levels(k)*feet2meters
        !$OMP PARALLEL DO PRIVATE(i,j,iloc,jstartm1) COLLAPSE(2)
        do j=jps,jpe
          do i=ips,ipe
            jstartm1 = kstartm1 + (j-jms)*ilen
            iloc = jstartm1+(i-ims)+1
            if(iloc<1 .or. iloc>ndata) then
              call ifi_check(-1,'Internal error: out of bounds in send_data.')
            endif
            if(data(iloc)==missing_value) then
              datapd(i-ips+1,j-jps+1,cfld) = spval
            else
              datapd(i-ips+1,j-jps+1,cfld) = data(iloc)
            endif
          enddo
        enddo
      endif
    enddo
  end subroutine send_data

  ! --------------------------------------------------------------------

  subroutine gather_ifi(local_buf_c,global_buf_c,receiving_rank) bind(C)
      use iso_c_binding, only: c_float, c_double,c_int64_t,c_f_pointer
      use mpi
      use ctlblk_mod, only: ME, MPI_COMM_COMP, JSTA, JEND, IM, JM, ICNT, IDSP
      implicit none

      integer(kind=c_int64_t), value :: receiving_rank
      type(c_ptr), value :: global_buf_c, local_buf_c
      real(kind=ifi_real_t), pointer :: global_buf(:,:), local_buf(:,:)

      integer :: type, iret

      call c_f_pointer(local_buf_c,local_buf,(/ im, jend-jsta+1 /))
      if(me==receiving_rank) then
        call c_f_pointer(global_buf_c,global_buf,(/ im, jm /))
      else
        call c_f_pointer(global_buf_c,global_buf,(/ 1,1 /))
      endif

      call gather_for_write(local_buf,global_buf,receiving_rank)
      
      ! if(ifi_real_t==c_double) then
      !   type=MPI_REAL8
      ! else
      !   type=MPI_REAL4
      ! endif

      ! call MPI_Gatherv(local_buf,icnt(me),type,global_buf,icnt,idsp,type,receiving_rank,mpi_comm_comp,iret)

  end subroutine gather_ifi

  ! --------------------------------------------------------------------

  subroutine gather_for_write(local_buf,global_buf,receiving_rank) bind(C)
      use iso_c_binding, only: c_float, c_double,c_int64_t,c_f_pointer
      use mpi
      use ctlblk_mod, only: ME, MPI_COMM_COMP, JSTA, JEND, IM, JM
      implicit none

      integer(kind=c_int64_t), value :: receiving_rank
      real(kind=ifi_real_t) :: global_buf(:,:), local_buf(:,:)

      integer :: type, iret, comm_size, count, i
      integer, allocatable :: counts(:),displ(:)
      
      if(ifi_real_t==c_double) then
        type=MPI_REAL8
      else
        type=MPI_REAL4
      endif

      call MPI_Comm_size(MPI_COMM_COMP,comm_size,iret)
      allocate(counts(comm_size))
      allocate(displ(comm_size))

      count = im*(jend-jsta+1)
      call MPI_Allgather(count,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_COMP,iret)

      displ(1)=0
      do i=2,comm_size
        displ(i) = displ(i-1) + counts(i-1)
      enddo

      call MPI_Gatherv(local_buf,count,type,global_buf,counts,displ,type,receiving_rank,mpi_comm_comp,iret)

      deallocate(counts,displ)

  end subroutine gather_for_write

  ! --------------------------------------------------------------------

  subroutine scatter_ifi(local_buf_c,global_buf_c,sending_rank)bind(C)
      use iso_c_binding, only: c_float, c_double, c_int64_t, c_f_pointer
      use mpi
      use ctlblk_mod, only: ME, MPI_COMM_COMP, ICNT, IDSP, JSTA, JEND,  &
     &                      IM, JM
      implicit none

      integer(kind=c_int64_t), value :: sending_rank
      type(c_ptr), value :: global_buf_c, local_buf_c
      real(kind=ifi_real_t), pointer :: global_buf(:,:), local_buf(:,:)

      integer :: rank, type, iret

      call c_f_pointer(local_buf_c,local_buf,(/ im, jend-jsta+1 /))
      if(me==sending_rank) then
        call c_f_pointer(global_buf_c,global_buf,(/ im, jm /))
      else
        call c_f_pointer(global_buf_c,global_buf,(/ 1,1 /))
      endif

      if(ifi_real_t==c_double) then
        type=MPI_REAL8
      else
        type=MPI_REAL4
      endif

      call mpi_scatterv(global_buf,icnt,idsp,type,local_buf,icnt(me),   &
     &                  type,sending_rank,MPI_COMM_COMP,iret)
      
  end subroutine scatter_ifi
  
  ! --------------------------------------------------------------------

  subroutine run_ifi()
    use ctlblk_mod, only: spval, lm, lp1, im, jsta_2l,jend_2u, ITPREC, IFHR, IFMIN, grib, jm, jsta,jend, me, num_procs, mpi_comm_comp
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
         OMGA ! "VVEL"
    use vrbls2d, only : &
         APCP, & ! IFI "APCP_surface" over ITPREC bucket time
         CAPE, CIN
    use rqstfld_mod, only: IGET
    use iso_c_binding, only: c_bool, c_int64_t

    implicit none

    INTERFACE ! implemented in EXCH.f
      SUBROUTINE EXCH_IFI(A) BIND(C)
        use ifi_type_mod, only: ifi_real_t
        implicit none
        real(kind=ifi_real_t) :: a(*)
      END SUBROUTINE EXCH_IFI
    END INTERFACE

    type(IFIData) :: hybr_vars, pres_vars, derived_vars, fip_algo_vars, flight_vars, cat_vars
    type(IFIAlgo) :: algo
    character(88) :: outfile

    real(c_double) :: fcst_lead_sec

    ! FIXME: Once the post is i-decomposed, replace the appropriate 1..im in this routine.
    
    if(grib/='grib2' .or. (IGET(1100)<=0 .and. IGET(1101)<=0 .and. IGET(1102)<=0)) then
      return ! nothing to do
    endif

    if(me==0) then
      print '(A)','Running libIFI to get icing products.'
      call ifi_print_copyright()
    end if


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
      call ifi_check(hybr_vars%add_ij_var('APCP_surface',int(1,c_int64_t),int(im,c_int64_t),&
      int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t), &
           APCP/ITPREC),'could not send APCP to IFI')
    else
      call ifi_check(hybr_vars%add_ij_var('APCP_surface',int(1,c_int64_t),int(im,c_int64_t),&
           int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t), &
           APCP),'could not send APCP to IFI')
    endif
    call ifi_check(hybr_vars%add_ij_var('CAPE_surface',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),CAPE), &
         'could not send CAPE to IFI')
    call ifi_check(hybr_vars%add_ij_var('CIN_surface',int(1,c_int64_t),int(im,c_int64_t),int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),CIN), &
         'could not send CIN to IFI')

    ! Copy 3D vars to IFI internal structures, inverting K dimension
    
    call add_ijk_var('HGT','zmid',zmid)
    call add_ijk_var('CIMIXR','QQI',QQI)
    call add_ijk_var('CLMR','QQW',QQW)
    call add_ijk_var('GRLE','QQG',QQG)
    call add_ijk_var('PRES','pmid',pmid)
    call add_ijk_var('RWMR','QQR',QQR)
    call add_ijk_var('SPFH','Q',Q)
    call add_ijk_var('TMP','T',T)
    call add_ijk_var('SNMR','QQS',QQS)
    call add_ijk_var('VVEL','OMGA',OMGA)

308 format(A,'_',I0,'.nc')

    fcst_lead_sec=IFHR*3600.0+IFMIN*60.0

    if(write_ifi_netcdf_debug_files) then
      write(outfile,308) 'hybr_vars',me
      call write_fip_output(hybr_vars,fcst_lead_sec,trim(outfile),.true.,'z0','z1')
    endif
    
    ! Initialize the IFI algorithm

    call ifi_check(algo%init(ifi_config,fcst_lead_sec,hybr_vars,ME,NUM_PROCS,MPI_COMM_COMP), &
         'could not initialize IFI algorithm')
    call ifi_check(algo%set_halo_callback(ifi_halo_callback), &
         'send halo callback function to IFI algorithm')
    call ifi_check(algo%set_minmax_callback(ifi_minmax_callback), &
         'send minmax callback function to IFI algorithm')
    call ifi_check(algo%set_gather_callback(gather_ifi), &
         'send gather callback function to IFI algorithm')
    call ifi_check(algo%set_scatter_callback(scatter_ifi), &
         'send scatter callback function to IFI algorithm')

    ! Change n>=spval to Config.MISSING in all input fields
    
    call ifi_check(algo%cap_to_missing(9e9),&
         'could not change post spval to Config.MISSING in input to IFI')

    ! Run the IFI algorithm
    
    call ifi_check(algo%calc_exner(),'calc_exner() failed')
    call ifi_check(algo%hybrid_to_pressure(),'hybrid_to_pressure() failed')

    call ifi_check(algo%get_pres_vars(pres_vars),'get_pres_vars()')
    if(write_ifi_netcdf_debug_files) then
      write(outfile,308) 'pres_vars',me
      call write_fip_output(pres_vars,fcst_lead_sec,trim(outfile),.true.,'z1','z0')
    endif

    call ifi_check(algo%discard_hybrid_level_vars(),'discard_hybrid_level_vars() failed')
    call ifi_check(algo%smooth_fields(),'smooth_fields() failed')

    if(write_ifi_netcdf_debug_files) then
      write(outfile,308) 'smooth_vars',me
      call write_fip_output(pres_vars,fcst_lead_sec,trim(outfile),.true.,'z1','z0')
    endif

    !call ifi_check(algo%imprecise_pressure_level_vars(),'imprecise_pressure_level_vars() failed')
    call ifi_check(algo%derive_fields(),'derive_fields() failed')

    call ifi_check(algo%get_derived_vars(derived_vars),'get_derived_vars()')
    if(write_ifi_netcdf_debug_files) then
      write(outfile,308) 'derived_vars',me
      call write_fip_output(derived_vars,fcst_lead_sec,trim(outfile),.false.,'z1','z0')
    endif

    !call ifi_check(algo%imprecise_derived_vars(),'imprecise_derived_vars() failed')
    call ifi_check(algo%run_fip_algo(),'run_fip_algo() failed')

    call ifi_check(algo%get_fip_algo_vars(pres_vars),'get_fip_algo_vars()')
    if(write_ifi_netcdf_debug_files) then
      write(outfile,308) 'fip_algo_vars',me
      call write_fip_output(pres_vars,fcst_lead_sec,trim(outfile),.false.,'z1','z0')
    endif

    !call ifi_check(algo%imprecise_fip_algo_vars(),'imprecise_fip_algo_vars() failed')
    call ifi_check(algo%discard_pressure_level_vars(),'discard_pressure_level_vars() failed')
    call ifi_check(algo%discard_derived_vars(),'discard_derived_vars() failed')
    call ifi_check(algo%pressure_to_flight(),'pressure_to_flight() failed')

    call ifi_check(algo%get_flight_vars(flight_vars),'get_flight_vars()')
    if(write_ifi_netcdf_debug_files) then
      write(outfile,308) 'flight_vars',me
      call write_fip_output(flight_vars,fcst_lead_sec,trim(outfile),.false.,'z1','z0')
    endif

    call ifi_check(algo%discard_fip_algo_vars(),'discard_fip_algo_vars() failed')
    call ifi_check(algo%make_icing_category(),'make_icing_category() failed')

    ! Get the final output fields:

    call ifi_check(algo%get_cat_vars(cat_vars),'get_cat_vars()')

    if(write_ifi_netcdf_debug_files) then
      write(outfile,308) 'cat_vars',me
      call write_fip_output(cat_vars,fcst_lead_sec,trim(outfile),.false.,'z1','z0')
    endif
    
    call send_data(cat_vars,'ICE_PROB',1100)
    call send_data(cat_vars,'SLD',1101)
    call send_data(cat_vars,'ICE_SEV_CAT',1102)
    
    ! When this subroutine ends, a Fortran-2003-compliant compiler
    ! will free all memory IFI uses, by calling the destructors (final
    ! routines) for cat_vars, hybr_vars, algo, and config.

  contains
    
    subroutine add_ijk_var(ifi_name,upp_name,upp_var)
      implicit none
      character(len=*), intent(in) :: ifi_name,upp_name
      real, intent(in) :: upp_var(IM,JSTA_2L:JEND_2U,LM)
      real(kind=ifi_real_t) :: ifi_var(IM,JSTA_2L:JEND_2U,LM)
      integer i,j,k

!$OMP PARALLEL DO COLLAPSE(2)
      do k=1,lm
        do j=jsta_2l,jend_2u
          do i=1,im
            ifi_var(i,j,k) = upp_var(i,j,lm-k+1)
          enddo
        enddo
      enddo
      
      call ifi_check(hybr_vars%add_ijk_var(ifi_name,int(1,c_int64_t),int(im,c_int64_t),&
           int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),&
           int(1,c_int64_t),int(lm,c_int64_t),ifi_var), &
           'could not send '//ifi_name//' ('//upp_name//') to IFI')
    end subroutine add_ijk_var

  end subroutine run_ifi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! FIXME: DELETE THIS TEST CODE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_range(var,count,missing_value_is_set,missing_value,min_not_miss,max_not_miss,all_missing)
    USE ieee_arithmetic
    use mpi
    use ctlblk_mod, only: mpi_comm_comp
    implicit none
    real(kind=ifi_real_t), intent(in) :: var(*)
    real(kind=ifi_real_t), intent(in) :: missing_value
    real(kind=ifi_real_t), intent(out) :: min_not_miss,max_not_miss
    integer, intent(in) :: count
    logical, intent(out) :: all_missing
    logical(c_bool), intent(in) :: missing_value_is_set

    real(kind=ifi_real_t) :: minv,maxv,epsilon, global_minv,global_maxv
    logical :: first
    integer :: i,type,iret
    real(kind=ifi_real_t), parameter :: zero = 0

    if(missing_value_is_set) then
      epsilon = abs(missing_value)*1e-4
      if(.not. epsilon+1>epsilon) then
        epsilon=1
      endif
    endif

    ! Initialize to out-of-bounds values so they'll stay that way if no valid values are found:
    minv = ieee_value(zero,ieee_positive_inf)
    maxv = ieee_value(zero,ieee_negative_inf)
    
    ! Find the min and max values in the array:
    first=.true.
    if(missing_value_is_set) then
      !$OMP PARALLEL DO REDUCTION(min:minv) REDUCTION(max:maxv)
      do i=1,count
        if(.not. abs(var(i)-missing_value)>epsilon) then
          if(first) then
            first=.false.
            minv=var(i)
            maxv=var(i)
          else
            minv=min(minv,var(i))
            maxv=max(maxv,var(i))
          end if
        endif
      end do
    else
      !$OMP PARALLEL DO REDUCTION(min:minv) REDUCTION(max:maxv)
      do i=1,count
        if(var(i)+1>var(i)) then
          if(first) then
            first=.false.
            minv=var(i)
            maxv=var(i)
          else
            minv=min(minv,var(i))
            maxv=max(maxv,var(i))
          end if
        endif
      end do
    endif

    if(ifi_real_t==c_double) then
      type=MPI_REAL8
    else
      type=MPI_REAL4
    endif

    !call MPI_Allreduce(minv,global_minv,1,type,MPI_MIN,mpi_comm_comp,iret)
    !call MPI_Allreduce(maxv,global_maxv,1,type,MPI_MAX,mpi_comm_comp,iret)
    global_minv=minv
    global_maxv=maxv
    
    ! If min or max are inf, -inf, or NaN, assume all values are missing.
    all_missing = (.not.(global_minv==global_minv) .or. .not.(global_maxv==global_maxv) &
         .or. .not.(global_minv+1>global_minv) .or. .not.(global_maxv+1>global_maxv))

    if(all_missing) then
      min_not_miss=-1
      max_not_miss=1
    else
      min_not_miss=global_minv
      max_not_miss=global_maxv
    endif
  end subroutine find_range

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nc_check(code,file,message)
    use netcdf
    implicit none
    integer, intent(in) :: code
    character(*), intent(in) :: file, message
    integer :: i
    if(code/=0) then
      write(0,93) trim(file),trim(message),trim(nf90_strerror(code)),code
      stop 1
    endif
    93 format(A,': ',A,': ',A,'(',I0,')')
  end subroutine nc_check

  subroutine write_fip_output(ifi_data,fcst_lead_sec,output_file,rename,z2dname,z3dname)
    use ifi_mod
    use iso_c_binding
    use netcdf
    use ctlblk_mod, only: me
    implicit none
    integer, parameter :: maxname=80
    integer, parameter :: maxvars=999
    logical, intent(in) :: rename
    double precision, intent(in) :: fcst_lead_sec
    character(len=*), intent(in) :: z2dname,z3dname
    type var_info
      character(len=maxname) :: varname, outname
      integer :: varid, ndims, dims(4), dimids(4), count
      real(ifi_real_t), pointer :: data(:)
      integer :: ims,ime,jms,jme,kms,kme
      integer :: ids,ide,jds,jde,kds,kde
      integer :: ips,ipe,jps,jpe,kps,kpe
      logical :: should_dealloc
      real(kind=ifi_real_t) :: missing_value
      real(kind=ifi_real_t) :: min_not_miss,max_not_miss
      logical(c_bool) :: missing_value_is_set
      logical :: all_missing
    end type var_info
    character(len=*), intent(in) :: output_file

    type(IFIData) :: ifi_data
    integer(c_int64_t) :: ids,ide, jds,jde, kds,kde
    integer(c_int64_t) :: ips,ipe, jps,jpe, kps,kpe
    integer(c_int64_t) :: ims,ime, jms,jme, kms,kme
    logical(c_bool) :: missing_value_is_set
    real(kind=ifi_real_t) :: missing_value
    integer :: count

    character(len=200) :: varname,nextname

    integer :: dimids(4), ncid, ivar, id, nvars, nx0, ny0, nz0, ntime, dims(4)
    type(var_info),target :: var_data(maxvars)

    call ifi_check(ifi_data%get_dims(ids,ide, jds,jde, kds,kde, ips,ipe, jps,jpe, kps,kpe),"get_dims")

    nx0=ide-ids+1
    ny0=jde-jds+1
    nz0=kde-kds+1
    ntime=1

    dims = (/nx0,ny0,nz0,ntime/)

    varname=' ' ! special value indicating "start at the beginning"

    ivar=0
    do while(ivar<maxvars)
      write(0,*) 'bottom of loop rank ',me
      nextname='internal_error_nextname_was_not_updated'
      call ifi_check(ifi_data%next_varname(varname,nextname),"next_varname")
      varname=nextname

      if(len_trim(nextname)<1) then
        exit
      else
        ivar=ivar+1
      endif
      
      if(ivar==1 .and. me==0) then
59      format(A,': write to file')
        print 59,output_file
        
        call nc_check(nf90_create(output_file,ior(NF90_NOCLOBBER,NF90_64BIT_OFFSET),ncid),&
             output_file,'nf90_create')

        call nc_check(nf90_put_att(ncid,NF90_GLOBAL,"fcst_lead_sec",fcst_lead_sec),&
             output_file,"nf90_put_att fcst_lead_sec")
        
        call nc_check(nf90_def_dim(ncid,"x0",nx0,dimids(1)),output_file,"nf90_def_dim x0")
        call nc_check(nf90_def_dim(ncid,"y0",ny0,dimids(2)),output_file,"nf90_def_dim y0")
        call nc_check(nf90_def_dim(ncid,z3dname,nz0,dimids(3)),output_file,"nf90_def_dim "//z3dname)
        call nc_check(nf90_def_dim(ncid,"time",ntime,dimids(4)),output_file,"nf90_def_dim time")
      endif

      var_data(ivar)%varname = varname
      var_data(ivar)%outname = varname
      var_data(ivar)%should_dealloc = .false.
      write(0,*) 'before read_var ',me
      var_data(ivar)%data => read_var(var_data(ivar),trim(varname),var_data(ivar)%should_dealloc)
      write(0,*) 'after read_var ',me
      if(me==0) then
        write(0,*) 'rank 0 defines var'
        call set_dims(var_data(ivar),dims,dimids)
        var_data(ivar)%varid = def_var(var_data(ivar),rename)
      endif
      write(0,*) 'bottom of loop rank ',me
    enddo

    if(me/=0) then
      ! Ranks that are not writing are done now
      return
    endif

    nvars=ivar

    if(nvars<1) then
      write(0,"(A,A)") output_file,': no data to write!'
      return
    endif

    if(me==0) then
      call nc_check(nf90_enddef(ncid),output_file,'nf90_enddef')
    endif

    write(0,*) 'before write loop ',me
    do ivar=1,nvars
24    format('  put var ',A)
      print 24,trim(var_data(ivar)%varname)
      if(me==0) then
        call find_range(var_data(ivar)%data,var_data(ivar)%count, &
             var_data(ivar)%missing_value_is_set,var_data(ivar)%missing_value,&
             var_data(ivar)%min_not_miss,var_data(ivar)%max_not_miss, &
             var_data(ivar)%all_missing)
        call write_var(var_data(ivar))
      endif
      if(var_data(ivar)%should_dealloc) then
        deallocate(var_data(ivar)%data)
      endif
      nullify(var_data(ivar)%data)
    enddo

    if(me==0) then
      call nc_check(nf90_close(ncid),output_file,"nf90_close")
    endif
  contains

    subroutine set_dims(var,dims,dimids)
      implicit none
      type(var_info), intent(inout), target :: var
      integer, intent(in) :: dims(:), dimids(:)
      character(len=:), pointer :: varname

      varname=>var%varname(1:len_trim(var%varname))

      if(me/=0) then
        ! Ranks that are not writing are done now
        return
      endif

      ! These must also be in IFITest.cc IFITest::write_netcdf

      if(varname=='x' .or. varname=='x0') then
        var%dims=(/nx0,1,1,1/)
        var%dimids=(/dimids(1),-1,-1,-1/)
        var%ndims=1
      else if(varname=='y' .or. varname=='y0') then
        var%dims=(/ny0,1,1,1/)
        var%dimids=(/dimids(2),-1,-1,-1/)
        var%ndims=1
      else if(varname=='z' .or. varname=='z0' .or.  varname=='z1' .or. &
           varname=='pressure_levels' .or. varname=='exner_levels') then
        var%dims=(/nz0,1,1,1/)
        var%dimids=(/dimids(3),-1,-1,-1/)
        var%ndims=1
      else if(varname=='latitude' .or. varname=='longitude') then
        var%dims=(/nx0,ny0,1,1/)
        var%dimids=(/dimids(1),dimids(2),-1,-1/)
        var%ndims=2
      else if(varname=='time') then
        var%dims=(/ntime,1,1,1/)
        var%dimids=(/dimids(4),-1,-1,-1/)
        var%ndims=1
      else if(kme<=kms) then
        var%dims=(/nx0,ny0,ntime,1/)
        var%dimids=(/dimids(1),dimids(2),dimids(4),-1/)
        var%ndims=3
      else
        var%dims=(/nx0,ny0,nz0,ntime/)
        var%dimids=dimids
        var%ndims=4
      endif

      var%count=var%dims(1)*var%dims(2)*var%dims(3)*var%dims(4)
    end subroutine set_dims

    subroutine write_var(var)
      implicit none
      type(var_info), intent(inout) :: var

      integer :: ones(var%ndims), dims(var%ndims)
      real(ifi_real_t) :: put(var%count)
      integer :: i,j,k,n,ilen,jlen,m,imlen,jmlen

      ones = 1
      dims = var%dims(1:var%ndims)

      if(me/=0) then
        ! Ranks that are not writing are done now
        return
      endif
      
      call nc_check(nf90_put_var(ncid=ncid,varid=var%varid,values=var%data, &
                                 start=ones,count=dims), &
                    output_file,"nf90_put_var "//trim(var%outname))

    end subroutine write_var

    integer function def_var(var,rename)
      use iso_c_binding, only: c_float
      implicit none
      logical :: rename
      type(var_info), intent(inout) :: var

      integer :: varid, xtype, dimids(var%ndims)
      character(len=100) :: outname

      if(rename) then
        var%outname=var%varname
        select case(trim(var%varname))
          case('CIMIXR')
            var%outname = 'ICMR'
          case('SPFH')
            var%outname = 'MIXR'
          case('GRLE')
            var%outname = 'GRMR'
          case('CAPE_surface')
            var%outname = 'CAPE'
          case('CIN_surface')
            var%outname = 'CIN'
          case('CLMR')
            var%outname = 'CLWMR'
          case('APCP_surface')
            var%outname = 'APCP1Hr'
        end select
      endif

      var%ims=ims ; var%ime=ime ; var%jms=jms ; var%jme=jme ; var%kms=kms ; var%kme=kme
      var%ids=ids ; var%ide=ide ; var%jds=jds ; var%jde=jde ; var%kds=kds ; var%kde=kde
      var%ips=ips ; var%ipe=ipe ; var%jps=jps ; var%jpe=jpe ; var%kps=kps ; var%kpe=kpe

      dimids = var%dimids(1:var%ndims)

      if(ifi_real_t==c_float) then
        xtype = NF90_FLOAT
      else
        xtype = NF90_DOUBLE
      endif

      call nc_check(nf90_def_var(ncid,trim(var%outname),xtype,dimids,def_var), &
           output_file,"nf90_def_var "//trim(var%outname))

      call nc_check(nf90_put_att(ncid,def_var,"min_value",var%min_not_miss), &
           output_file,"nf90_put_att "//trim(var%outname)//" min_value")

      call nc_check(nf90_put_att(ncid,def_var,"max_value",var%max_not_miss), &
           output_file,"nf90_put_att "//trim(var%outname)//" max_value")

      if(var%missing_value_is_set) then
        call nc_check(nf90_put_att(ncid,def_var,"_FillValue",-9999.0), &
             output_file,"nf90_put_att "//trim(var%outname)//" _FillValue")
      end if

    end function def_var

    function read_var(var,varname,should_dealloc)
      implicit none
      type(var_info), intent(inout) :: var
      real(kind=ifi_real_t), pointer :: read_var(:)
      real(kind=ifi_real_t), pointer :: local_data_1D(:)
      character(len=*), intent(in) :: varname
      logical, intent(out) :: should_dealloc

      real(kind=ifi_real_t), allocatable :: local_data(:,:),global_data(:,:)
      real(kind=ifi_real_t), pointer :: global_data_1D(:)

      integer :: nxny_local,nxny_global,nz,count,i,j,k
      integer :: local_ilen,local_jlen,local_klen,local_index
      integer :: global_ilen,global_jlen,global_klen,global_index
      type(c_ptr) :: global_cptr,local_cptr

      local_data_1D => ifi_data%get_data(trim(varname),missing_value_is_set,missing_value, &
           ims,ime,jms,jme,kms,kme, ids,ide,jds,jde,kds,kde, ips,ipe,jps,jpe,kps,kpe)

      global_ilen=ide-ids+1
      global_jlen=jde-jds+1
      global_klen=kde-kds+1
      local_ilen=ipe-ips+1
      local_jlen=jpe-jps+1
      local_klen=kpe-kps+1

      if(.not.associated(local_data_1D)) then
        write(0,38) trim(varname)
38      format("IFI did not produce ",A," variable.")
        stop 19
      end if

      count=global_ilen*global_jlen*global_klen
      if(count<=0) then
        write(0,39) trim(varname),count
39      format("IFI variable ",A," had no data (size=",I0,")")
        stop 20
      endif

      if(ids==ide .or. jds==jde) then
        ! This is a 1D variable so we write it as is from rank 0
        read_var => local_data_1D
        should_dealloc = .false.
        return
      endif

40    format('var ',A,' has bad ',A,'. Got: ',I0,' but expected: ',I0)

      if(global_ilen /= nx0) then
        write(0,40) trim(varname),'global_ilen',global_ilen,nx0
        stop 31
      endif

      if(global_jlen /= ny0) then
        write(0,40) trim(varname),'global_jlen',global_jlen,ny0
        stop 32
      endif

      if(global_klen/=1 .and. global_klen /= nz0) then
        write(0,40) trim(varname),'global_klen',global_klen,nz0
        stop 33
      endif


      allocate(local_data(local_ilen,local_jlen))
      if(me==0) then
        allocate(global_data(global_ilen,global_jlen))
        allocate(global_data_1D(global_ilen*global_jlen*global_klen))
      else
        allocate(global_data(1,1))
        allocate(global_data_1D(1))
      endif

      do k=kds,kde
        do j=jps,jpe
          do i=ips,ipe
            local_data(i-ips+1,j-jps+1) = local_data_1D( 1 + (i-ims) + ((j-jms) + (k-kms)*(jme-jms+1))*(ime-ims+1) )
          enddo
        enddo
        call gather_for_write(local_data,global_data,0)
        if(me==0) then
          do j=1,global_jlen
            do i=1,global_ilen
              global_data_1D(1 + (i-1) + ((j-1) + (k-kds)*global_jlen)*global_ilen) = &
                   global_data(i,j)
            enddo
          enddo
        endif
      enddo

      ! Ranks that are not writing data are done.
      if(me/=0) then
        ! This rank does not have the global data, and will not use the data anyway,
        ! but it wants an array, so we'll send it the local data that is managed
        ! internally by libIFI.
        read_var => local_data_1D
        should_dealloc = .false.
        deallocate(global_data_1D)
      else
        read_var => global_data_1D
        should_dealloc = .true.
      endif
      
      deallocate(global_data)
      deallocate(local_data)
    end function read_var
  end subroutine write_fip_output


end module upp_ifi_mod


#endif
