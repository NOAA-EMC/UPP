!> @file
!> @brief module: upp_ifi_mod declares variables used in the in-flight icing
!(IFI) routines 
module upp_ifi_mod

#ifdef USE_IFI
  ! Use the libIFI API modules.
  use ifi_type_mod
  use ifi_mod
#endif

  use iso_c_binding, only: c_bool, c_int64_t, c_ptr, c_double, c_float
  implicit none

  private

  public run_ifi, set_ifi_dims, ifi_real_t, write_ifi_debug_files, &
       first_supported_ifi_hour, last_supported_ifi_hour

  logical :: write_ifi_debug_files = .false.
  real, parameter :: first_supported_ifi_hour = 1
  real, parameter :: last_supported_ifi_hour = 21

#ifndef USE_IFI
  ! Stubs need their own version of the ifi_real_t since it is defined in ifi_type_mod
  integer, parameter :: ifi_real_t = c_float

#else
  ! Actual IFI code needs a few more module-scope variables.

  real, parameter :: feet2meters = 0.3048

  type(IFIConfig) :: ifi_config
  logical :: have_read_ifi_config = .false.

  ! Communication-related variables for gathers and
  ! scatters. Gathering is done in the i direction first, all to one
  ! rank in each row. Then that rank in each row gathers to the
  ! destination rank. Scatter follows the same process, in reverse.

  integer :: ifi_mpi_real_kind=-999
  integer :: row_comm = -999 ! communicator for i-direction messages
  integer :: col_comm = -999 ! communicator for j-direction messages
  integer :: row_comm_rank=-999, col_comm_rank=-999 ! Rank in the i & j comm communicators
  integer :: row_comm_size=-999, col_comm_size=-999 ! Size of the i & j comm communicators
  integer :: grid_rank ! rank within mpi_comm_comp

  integer, allocatable :: ista_grid(:) ! ista for each rank in mpi_comm_comp in i direction
  integer, allocatable :: jsta_grid(:) ! jsta for each rank in mpi_comm_comp in j direction

  real, allocatable :: rearrange(:,:) ! for rearranging data on local processor
  integer, allocatable :: row_comm_count(:) ! scatterv&gatherv receive counts for i-direction collectives
  integer, allocatable :: row_comm_displ(:) ! scatterv&gatherv displacements for i-direction collectives
  integer, allocatable :: row_ista(:) ! i starts for rearranging data (icomm_size)
  integer, allocatable :: row_iend(:) ! i ends for rearranging data (icomm_size)
  real, allocatable :: rearrange_row1d(:) ! for rearranging data for an entire row (roots of icomm communicators)

  real, allocatable :: rearrange_row2d(:,:) ! for rearranging data for an entire row (roots of icomm communicators)
  integer, allocatable :: col_comm_count(:) ! scatterv&gatherv receive counts for j-direction collectives
  integer, allocatable :: col_comm_displ(:) ! scatterv&gatherv displacements for j-direction collectives
  integer, allocatable :: col_jsta(:) ! j starts for rearranging data (jcomm_size)
  integer, allocatable :: col_jend(:) ! j ends for rearranging data (jcomm_size)

  integer, allocatable :: isize_grid(:) ! iend-ista+1 for each gridpoint
#endif
  
contains !==============================================================

  subroutine send_missing_data(ient)
    ! Fills IFI fields in index ient (of IGET) with missing values.
    !
    ! This subroutine is both in the stub and actual IFI code.
    ! The stubs always send missing data, while the actual code
    ! sends missing data for unsupported hours.
    use CTLBLK_mod, only: ifi_nflight, ifi_flight_levels
    use ctlblk_mod, only: spval, ista, iend, jsta, jend, lm, im, cfld, datapd, fld_info, ifi_flight_levels, jsta_2l, jend_2u,me
    use rqstfld_mod, only: iget, iavblfld, lvlsxml, lvls
    implicit none

    integer, intent(in) :: ient

    ! Locals
    integer :: i,j,k

    logical, save :: wrote_message = .false. ! guarded by an OMP CRITICAL below

    if(size(IGET)<ient) then
      return
    endif

    if(IGET(ient)<1 .or. IGET(ient)>size(LVLS,2)) then
      return
    endif

    ! Go level-by-level writing grib2 if requested
    do k=1,ifi_nflight
      if(k>size(LVLS,1)) then
        write(0,*) 'ERROR: LVLS is too small to handle ',ifi_nflight,' vertical levels'
        return
      endif
      if(LVLS(k,IGET(ient))>0) then
        if(.not.wrote_message) then
          !$OMP CRITICAL
          if(.not.wrote_message.and.me==0) then
            write(*,*) 'This post cannot produce IFI icing products because it was not compiled with libIFI.'
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

#ifndef USE_IFI
! IFI stubs

  subroutine set_ifi_dims()
    use CTLBLK_mod, only: ifi_nflight, ifi_flight_levels
    implicit none
    integer :: i

    ! Bogus fill value for flight levels to prevent a crash
    ifi_nflight = 60
    if(allocated(ifi_flight_levels)) then
      deallocate(ifi_flight_levels)
    endif
    allocate(ifi_flight_levels(ifi_nflight))
    do i=1,ifi_nflight
      ifi_flight_levels(i) = 500*i
    enddo
  end subroutine set_ifi_dims

  ! --------------------------------------------------------------------

  subroutine run_ifi()
    ! Fill any requested IFI fields with missing data.
    call send_missing_data(1007)  ! ICE_PROB = missing
    call send_missing_data(1008)  ! SLD = missing
    call send_missing_data(1009)  ! ICE_SEV_CAT = missing
    call send_missing_data(1010)  ! ICE_SEV_CAT = missing
  end subroutine run_ifi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Actual IFI code

  subroutine make_communicators
    use ctlblk_mod, only: jsta, jend, ista, iend, mpi_comm_comp, im, jm
    use mpi
    use iso_c_binding, only: c_sizeof
    implicit none

    integer, allocatable :: key(:) ! ensure ranks are in a consistent order
    integer :: size, ierr, local_count, accum, irank, i
    real(kind=ifi_real_t) :: testreal

    if(allocated(jsta_grid)) then
      return ! Already made communicators
    endif

    ! Use the right MPI type for the ifi real kind
    if(c_sizeof(testreal)==4) then
      ifi_mpi_real_kind=MPI_REAL4
    else
      ifi_mpi_real_kind=MPI_REAL8
    endif

    ! Get the correct size of the communicator before we try to split it.
    call MPI_Comm_size(mpi_comm_comp,size,ierr)
    call MPI_Comm_rank(mpi_comm_comp,grid_rank,ierr)

! 929 format('Rank ',I0,' ista=',I0,' jsta=',I0)
!     print 929,grid_rank,ista,jsta

    if(allocated(ista_grid)) then
      deallocate(ista_grid)
      deallocate(jsta_grid)
      deallocate(rearrange)
      deallocate(rearrange_row1d)
      deallocate(rearrange_row2d)
      deallocate(row_ista)
      deallocate(row_iend)
      deallocate(row_comm_count)
      deallocate(row_comm_displ)
      deallocate(col_jsta)
      deallocate(col_jend)
      deallocate(col_comm_count)
      deallocate(col_comm_displ)
    endif

    ! Get the start locations on every rank. We need this for the key,
    ! and for determining root ranks.
    allocate(ista_grid(size))
    allocate(jsta_grid(size))
    ista_grid=-1
    jsta_grid=-1
    call MPI_Allgather(ista,1,MPI_INTEGER,ista_grid,1,MPI_INTEGER,mpi_comm_comp,ierr)
    call MPI_Allgather(jsta,1,MPI_INTEGER,jsta_grid,1,MPI_INTEGER,mpi_comm_comp,ierr)

    ! Ensure ranks are arranged in a consistent order.
    allocate(key(size))
    do i=1,size
      if(ista_grid(i)==-1) then
        write(0,*) 'invalid ista_grid ',ista_grid(i)
        call mpi_abort(mpi_comm_world,1,ierr)
      endif
      if(jsta_grid(i)==-1) then
        write(0,*) 'invalid jsta_grid ',jsta_grid(i)
        call mpi_abort(mpi_comm_world,1,ierr)
      endif
      key(i) = ista_grid(i) + im*(jsta_grid(i)-1)
    enddo

    ! Create a communicator for row gather/scatter (all ranks that share a jsta)
    if(all(jsta_grid==jsta_grid(1))) then
      row_comm = mpi_comm_comp
    else
      call mpi_comm_split(mpi_comm_comp,jsta_grid(grid_rank+1),key(grid_rank+1),row_comm,ierr)
    endif
    if(row_comm==MPI_COMM_NULL .or. row_comm==0) then
      write(0,*) 'MPI_Comm_split gave MPI_COMM_NULL for row_comm'
      call mpi_abort(mpi_comm_world,1,ierr)
    endif
    call mpi_comm_rank(row_comm,row_comm_rank,ierr)
    call mpi_comm_size(row_comm,row_comm_size,ierr)

    ! Create a communicator for column gather/scatter (all ranks that share an ista)
    if(all(ista_grid==ista_grid(1))) then
      col_comm = mpi_comm_comp
    else
      call mpi_comm_split(mpi_comm_comp,ista_grid(grid_rank+1),key(grid_rank+1),col_comm,ierr)
    endif
    if(col_comm==MPI_COMM_NULL .or. col_comm==0) then
      write(0,*) 'MPI_Comm_split gave MPI_COMM_NULL for col_comm'
      call mpi_abort(mpi_comm_world,1,ierr)
    endif
    call mpi_comm_rank(col_comm,col_comm_rank,ierr)
    call mpi_comm_size(col_comm,col_comm_size,ierr)

    ! Done with the key.
    deallocate(key)

    ! Allocate the arrays used to rearrange data.
    allocate(rearrange(ista:iend,jsta:jend))
    allocate(rearrange_row1d(im*(jend-jsta+1)))
    allocate(rearrange_row2d(1:im,jsta:jend))

    ! Calculate information for row MPI_Gatherv and MPI_Scatterv calls.
    allocate(row_ista(row_comm_size))
    call MPI_Allgather(ista,1,MPI_INTEGER,row_ista,1,MPI_INTEGER,row_comm,ierr)
    allocate(row_iend(row_comm_size))
    call MPI_Allgather(iend,1,MPI_INTEGER,row_iend,1,MPI_INTEGER,row_comm,ierr)
    allocate(row_comm_count(row_comm_size))
    allocate(row_comm_displ(row_comm_size))
    call get_count_and_displ(row_comm,row_comm_size,(iend-ista+1)*(jend-jsta+1), &
         row_comm_count,row_comm_displ)

    ! Calculate information for j-direction MPI_Gatherv and MPI_Scatterv calls.
    allocate(col_jsta(col_comm_size))
    call MPI_Allgather(jsta,1,MPI_INTEGER,col_jsta,1,MPI_INTEGER,col_comm,ierr)
    allocate(col_jend(col_comm_size))
    call MPI_Allgather(jend,1,MPI_INTEGER,col_jend,1,MPI_INTEGER,col_comm,ierr)
    allocate(col_comm_count(col_comm_size))
    allocate(col_comm_displ(col_comm_size))
    call get_count_and_displ(col_comm,col_comm_size,im*(jend-jsta+1), &
         col_comm_count,col_comm_displ)

  end subroutine make_communicators

  ! --------------------------------------------------------------------

  subroutine get_count_and_displ(comm,comm_size,here,counts,displacements)
    use mpi
    implicit none
    integer, intent(in) :: comm,comm_size,here
    integer, intent(inout) :: counts(comm_size),displacements(comm_size)
    integer :: accum, ierr, i

    ! Get counts from all ranks
    call MPI_Allgather(here,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)

    ! Displacements are a cumulative sum starting at 0 for rank 0
    accum=0
    do i=1,comm_size
      displacements(i) = accum
      accum = accum+counts(i)
    end do
  end subroutine get_count_and_displ

  ! --------------------------------------------------------------------

  subroutine find_comm_roots(global_rank,row_root,col_root)
    ! Find the roots of the row and column communicators for a gather
    ! or scatter with the specified rank as the location of the global
    ! array.
    use mpi
    implicit none

    integer, intent(in) :: global_rank
    integer, intent(inout) :: row_root,col_root
    integer :: r, destination_ista, destination_jsta, ierr

    ! print *,'find roots for rank ',global_rank

    destination_ista=ista_grid(global_rank+1)
    destination_jsta=jsta_grid(global_rank+1)

    row_root=-1
    do r=1,row_comm_size
      if(row_ista(r)==destination_ista) then
        row_root=r-1
        exit
      endif
    enddo

    col_root=-1
    do r=1,col_comm_size
      if(col_jsta(r)==destination_jsta) then
        col_root=r-1
        exit
      endif
    enddo

! 201 format('root is row_root=',I0,' col_root=',I0,' for rank ',I0)
!     print 201,row_root,col_root,global_rank
    
    if(col_root==-1) then
      write(0,'(A,I0)') 'ABORT: Could not find column root rank for rank ',global_rank
      call mpi_abort(mpi_comm_world,1,ierr)
    endif
    if(row_root==-1) then
      write(0,'(A,I0)') 'ABORT: Could not find row root rank for rank ',global_rank
      call mpi_abort(mpi_comm_world,1,ierr)
    endif
  end subroutine find_comm_roots

  ! --------------------------------------------------------------------

  subroutine global_gather(local,global,global_rank,ista_local,iend_local,jsta_local,jend_local)
    ! Gather from local arrays on all ranks to global array on specified rank
    use mpi
    use ctlblk_mod, only: jsta, jend, ista, iend, im, jm
    use iso_c_binding, only: c_sizeof
    implicit none

    real(kind=ifi_real_t), intent(in) :: local(ista_local:iend_local,jsta_local:jend_local)
    real(kind=ifi_real_t), intent(out) :: global(im,jm)
    integer, intent(in) :: global_rank,ista_local,iend_local,jsta_local,jend_local

    integer :: i,j,r,inindex,row_root,col_root, destination_ista, destination_jsta,ni_rank,nj_rank,idxstart, ierr

    if(col_comm==0 .or. col_comm==MPI_COMM_NULL) then
      write(0,*) 'somehow, col_comm became invalid ',col_comm
      call mpi_abort(mpi_comm_world,1,ierr)
    endif

    if(row_comm==0 .or. row_comm==MPI_COMM_NULL) then
      write(0,*) 'somehow, row_comm became invalid ',row_comm
      call mpi_abort(mpi_comm_world,1,ierr)
    endif

    ! Find out who is the root of the i and j direction communications.
    call find_comm_roots(global_rank,row_root,col_root)

    ! Store the data in a contiguous local array
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j)
    do j=jsta,jend
      do i=ista,iend
        rearrange(i,j) = local(i,j)
      enddo
    enddo

    ! Gather across all ranks in row
    call MPI_Gatherv(rearrange,(jend-jsta+1)*(iend-ista+1),ifi_mpi_real_kind, &
         rearrange_row1d,row_comm_count,row_comm_displ,ifi_mpi_real_kind,row_root,row_comm,ierr)

    if(row_comm_rank==row_root) then
      ! Rearrange from i-j-rank dimensions to i-j dimensions.
      idxstart=1
      do r=1,row_comm_size

        ni_rank=row_iend(r)-row_ista(r)+1
        nj_rank=jend-jsta+1
        ! print *,'r,ista,iend=',r,row_ista(r),row_iend(r)

        !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(inindex)
        do j=0,nj_rank-1
          do i=0,ni_rank-1
            inindex = idxstart + i + ni_rank*j
            rearrange_row2d(i+row_ista(r),j+jsta) = rearrange_row1d(inindex)
          enddo
        enddo

        idxstart = idxstart + ni_rank*nj_rank
      enddo

      ! Gather along columns
      ! print *,'global gatherv to rank ',grid_rank
      call MPI_Gatherv(rearrange_row2d,im*(jend-jsta+1),ifi_mpi_real_kind, &
           global,col_comm_count,col_comm_displ,ifi_mpi_real_kind,col_root,col_comm,ierr)
    endif

    ! print *,'global_gather success'
  end subroutine global_gather

  ! --------------------------------------------------------------------

  subroutine global_scatter(local,global,global_rank,ista_local,iend_local,jsta_local,jend_local)
    ! Scatter from global array at specified rank to local arrays on all ranks.
    ! NOTE: Does NOT update halo regions!
    use mpi
    use ctlblk_mod, only: jsta, jend, ista, iend, im, jm
    use iso_c_binding, only: c_sizeof
    implicit none

    real(kind=ifi_real_t), intent(out) :: local(ista_local:iend_local,jsta_local:jend_local)
    real(kind=ifi_real_t), intent(in) :: global(im,jm)
    integer, intent(in) :: global_rank,ista_local,iend_local,jsta_local,jend_local

    integer :: i,j,r,outindex,idxstart,ni_rank,nj_rank,col_root,ierr,row_root

    if(col_comm==0 .or. col_comm==MPI_COMM_NULL) then
      write(0,*) 'somehow, col_comm became invalid ',col_comm
      call mpi_abort(mpi_comm_world,1,ierr)
    endif

    if(row_comm==0 .or. row_comm==MPI_COMM_NULL) then
      write(0,*) 'somehow, row_comm became invalid ',row_comm
      call mpi_abort(mpi_comm_world,1,ierr)
    endif

    ! Find out who is the root of the i and j direction communications.
    call find_comm_roots(global_rank,row_root,col_root)

    if(row_comm_rank==row_root) then
      ! Distribute in J direction.
      ! print *,'global scatterv from rank ',grid_rank
      call MPI_Scatterv(global, col_comm_count, col_comm_displ, ifi_mpi_real_kind, &
           rearrange_row2d ,im*(jend-jsta+1), ifi_mpi_real_kind, col_root, col_comm, ierr)
      ! Rearrange from i-j dimensions to i-j-rank dimensions
      idxstart=1
      do r=1,row_comm_size

        ni_rank=row_iend(r)-row_ista(r)+1
        nj_rank=jend-jsta+1

        !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(outindex)
        do j=0,nj_rank-1
          do i=0,ni_rank-1
            outindex = idxstart + i + ni_rank*j
            rearrange_row1d(outindex) = rearrange_row2d(i+row_ista(r),j+jsta)
          enddo
        enddo

        idxstart = idxstart + ni_rank*nj_rank
      enddo
    endif

    ! Distribute across row
    call MPI_Scatterv(rearrange_row1d,row_comm_count,row_comm_displ,ifi_mpi_real_kind, &
         rearrange,(iend-ista+1)*(jend-jsta+1),ifi_mpi_real_kind,row_root,row_comm,ierr)

    ! Copy back to contiguous local array.
    ! NOTE: Does NOT update the halo!
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j)
    do j=jsta,jend
      do i=ista,iend
        local(i,j) = rearrange(i,j)
      enddo
    enddo

    ! print *,'global_scatter success'

  end subroutine global_scatter

  ! --------------------------------------------------------------------

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
    integer :: i

    ! Make sure the config was read in
    call read_ifi_config()

    ! Get the flight levels
    call ifi_check(ifi_config%get_flight_levels_feet(config_flight_levels_feet), &
         'cannot get flight levels in feet from IFI config')

    ! Convert from integer to real:
    ifi_nflight = size(config_flight_levels_feet)
    if(allocated(ifi_flight_levels)) then
      deallocate(ifi_flight_levels)
    endif
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
    use ctlblk_mod, only: spval, jsta, jend, lm, cfld, datapd, fld_info, ifi_flight_levels, jsta_2l, jend_2u, ista, iend, ista_2l, iend_2u
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
    integer :: i,j,k,nj_local,jpad,ilen,jlen,kstartm1,jstartm1,iloc,ndata,ji,ni_local,ipad

    if(.not. IGET(ient)>0) then
      return
    endif

    ! WARNING: do not deallocate the data pointer. It is managed by the IFI C++ library.
    data => vars%get_data(trim(name),missing_value_is_set,missing_value, &
           ims,ime,jms,jme,kms,kme, ids,ide,jds,jde,kds,kde, ips,ipe,jps,jpe,kps,kpe)
    if(.not.missing_value_is_set) then
      missing_value = MISSING
    endif

    ! Get dimensions and do some sanity checks
    nj_local = jpe-jps+1
    ni_local = ipe-ips+1
    jpad = jps-jms
    ipad = ips-ims
    if(nj_local/=jend-jsta+1 .or. jpad/=jsta-jsta_2l .or. ni_local/=iend-ista+1 .or. ipad/=ista-ista_2l) then
94    format('   ',A,' = ',I0,', ',I0)
83    format('Warning: ',A,': IFI output j bounds do not match UPP j bounds')
      print 83,trim(name)
      print 94,'jps,jsta',jps,jsta
      print 94,'jpe,jend',jpe,jend
      print 94,'jms,jsta_2l',jms,jsta_2l
      print 94,'jme,jend_2u',jme,jend_2u
      print 94,'jpad,jsta-jsta_2l',jpad,jsta-jsta_2l
      print 94,'nj_local,jend-jsta+1',nj_local,jend-jsta+1
      print 94,'ips,ista',ips,ista
      print 94,'ipe,iend',ipe,iend
      print 94,'ims,ista_2l',ims,ista_2l
      print 94,'ime,iend_2u',ime,iend_2u
      print 94,'ipad,ista-ista_2l',ipad,ista-ista_2l
      print 94,'ni_local,iend-ista+1',ni_local,iend-ista+1
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
            if(iloc<1 .or. iloc>size(data)) then
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
      use ctlblk_mod, only: ME, MPI_COMM_COMP, JSTA, JEND, IM, JM, ICNT, IDSP, ISTA, IEND
      implicit none

      integer(kind=c_int64_t), value :: receiving_rank
      type(c_ptr), value :: global_buf_c, local_buf_c
      real(kind=ifi_real_t), pointer :: global_buf(:,:), local_buf(:,:)

      integer :: type, iret

      call c_f_pointer(local_buf_c,local_buf,(/ iend-ista+1, jend-jsta+1 /))
      if(me==receiving_rank) then
        call c_f_pointer(global_buf_c,global_buf,(/ im, jm /))
      else
        call c_f_pointer(global_buf_c,global_buf,(/ 1,1 /))
      endif

      call gather_for_write(local_buf,global_buf,receiving_rank)
      
  end subroutine gather_ifi

  ! --------------------------------------------------------------------

  subroutine gather_for_write(local_buf,global_buf,receiving_rank_c) bind(C)
      use iso_c_binding, only: c_float, c_double,c_int64_t,c_f_pointer
      use mpi
      use ctlblk_mod, only: ME, MPI_COMM_COMP, JSTA, JEND, IM, JM, ISTA, ISTA, IEND
      implicit none

      integer(kind=c_int64_t), value :: receiving_rank_c
      real(kind=ifi_real_t) :: global_buf(:,:), local_buf(:,:)
      integer :: receiving_rank
      receiving_rank = receiving_rank_c
      call global_gather(local_buf,global_buf,receiving_rank,ista,iend,jsta,jend)
  end subroutine gather_for_write

  ! --------------------------------------------------------------------

  subroutine scatter_ifi(local_buf_c,global_buf_c,sending_rank_c)bind(C)
      use iso_c_binding, only: c_float, c_double, c_int64_t, c_f_pointer
      use mpi
      use ctlblk_mod, only: JSTA, JEND, IM, JM, ISTA, IEND, ME
      implicit none

      integer(kind=c_int64_t), value :: sending_rank_c
      type(c_ptr), value :: global_buf_c, local_buf_c
      real(kind=ifi_real_t), pointer :: global_buf(:,:), local_buf(:,:)
      integer :: sending_rank

      sending_rank=sending_rank_c

      call c_f_pointer(local_buf_c,local_buf,(/ iend-ista+1, jend-jsta+1 /))
      if(me==sending_rank) then
        call c_f_pointer(global_buf_c,global_buf,(/ im, jm /))
      else
        call c_f_pointer(global_buf_c,global_buf,(/ 1,1 /))
      endif

      call global_scatter(local_buf,global_buf,sending_rank,ista,iend,jsta,jend)
  end subroutine scatter_ifi

  ! --------------------------------------------------------------------

  subroutine run_ifi()
    use ctlblk_mod, only: spval, lm, lp1, im, jsta_2l,jend_2u, ITPREC, IFHR, IFMIN, grib, &
         jm, jsta,jend, me, num_procs, mpi_comm_comp, ista, iend, ista_2l, iend_2u, dtq2
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
         IFI_APCP, & ! IFI "APCP_surface" over ITPREC bucket time
         CAPE, CIN, &
         AVGPREC_CONT ! Backup plan for points where IFI_APCP is invalid or 0
    use rqstfld_mod, only: IGET
    use iso_c_binding, only: c_bool, c_int64_t

    implicit none

    INTERFACE ! implemented in EXCH_c_float.f
      SUBROUTINE EXCH_c_float(A) BIND(C)
        use ifi_type_mod, only: ifi_real_t
        implicit none
        real(kind=ifi_real_t) :: a(*)
      END SUBROUTINE EXCH_c_float
    END INTERFACE

    type(IFIData) :: hybr_vars, pres_vars, derived_vars, fip_algo_vars, flight_vars, cat_vars
    type(IFIAlgo) :: algo
    character(88) :: outfile
    real :: to_hourly
    real(c_double) :: fcst_lead_sec
    integer :: i, j

    if(IFHR < first_supported_ifi_hour .or. IFHR > last_supported_ifi_hour) then
1838  format('You requested IFI fields for hour ',I0)
      print 1838, IFHR
1839  format('IFI fields are only available from hours ',I0,' through ',I0,'.')
      print 1839, first_supported_ifi_hour, last_supported_ifi_hour
1840  format(A)
      print 1840, 'Will output missing data for IFI fields for this time.'
      call send_missing_data(1007)
      call send_missing_data(1008)
      call send_missing_data(1009)
      call send_missing_data(1010)
      return
    endif

    ! FIXME: Once the post is i-decomposed, replace the appropriate 1..im in this routine.
    
    if(grib/='grib2' .or. (IGET(1007)<=0 .and. IGET(1008)<=0 .and. IGET(1009)<=0 .and. IGET(1010)<=0)) then
      if(me==0) then
        print '(A)','IFI fields were not requested; skipping IFI'
      endif
      return ! nothing to do
    endif

    if(me==0) then
      print '(A)','Running libIFI to get icing products.'
      call ifi_print_copyright()

      if(write_ifi_debug_files) then
        print '(A)','Will output IFI debug files.'
      else
        print '(A)','Will NOT output IFI debug files.'
      endif
    endif

    if(write_ifi_debug_files) then
      call make_communicators()
    endif

    ! Read config and initialize input data structures:

    call read_ifi_config()
    call ifi_check(hybr_vars%init(int(1,c_int64_t),int(im,c_int64_t),int(1,c_int64_t),int(jm,c_int64_t),&
         int(1,c_int64_t),int(lm,c_int64_t),int( ista,c_int64_t),int(iend,c_int64_t),int(jsta,c_int64_t),&
         int(jend,c_int64_t),int(1,c_int64_t),int(lm,c_int64_t)),&
         'could not initialize IFI input data structures')

    ! Copy 2D vars to IFI internal structures:
    
    call ifi_check(hybr_vars%add_ij_var('topography',int(ista_2l,c_int64_t),int(iend_2u,c_int64_t),&
                                                     int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t), &
         zint(:,:,LP1)),'could not send topography to IFI')

    to_hourly=3600*1000.0/dtq2
!$OMP PARALLEL DO COLLAPSE(2)
    do j=jsta,jend
       do i=ista,iend
          if(ifi_apcp(i,j)>9e9 .or. ifi_apcp(i,j)<-9e9 .or. ifi_apcp(i,j)==0) then
             ifi_apcp(i,j) = avgprec_cont(i,j)*to_hourly
          else if(ITPREC>1e-5) then
             ifi_apcp(i,j) = ifi_apcp(i,j)/itprec
          endif
       enddo
    enddo

    call ifi_check(hybr_vars%add_ij_var('APCP_surface',int(ista_2l,c_int64_t),int(iend_2u,c_int64_t),&
                                                       int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t), &
         IFI_APCP),'could not send IFI_APCP to IFI')
    call ifi_check(hybr_vars%add_ij_var('CAPE_surface',int(ista_2l,c_int64_t),int(iend_2u,c_int64_t), &
                                                       int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),CAPE), &
         'could not send CAPE to IFI')
    call ifi_check(hybr_vars%add_ij_var('CIN_surface',int(ista_2l,c_int64_t),int(iend_2u,c_int64_t), &
                                                      int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),CIN), &
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

    if(write_ifi_debug_files) then
      write(outfile,308) 'hybr_vars',me
      call write_fip_output(hybr_vars,fcst_lead_sec,trim(outfile),.true.,'z0','z1')
    endif
   
    ! Initialize the IFI algorithm

    call ifi_check(algo%init(ifi_config,fcst_lead_sec,hybr_vars,ME,NUM_PROCS,MPI_COMM_COMP), &
         'could not initialize IFI algorithm')

    ! Run the IFI algorithm
    
    call ifi_check(algo%calc_exner(),'calc_exner() failed')
    call ifi_check(algo%hybrid_to_pressure(),'hybrid_to_pressure() failed')

    call ifi_check(algo%get_pres_vars(pres_vars),'get_pres_vars()')
    if(write_ifi_debug_files) then
      write(outfile,308) 'pres_vars',me
      call write_fip_output(pres_vars,fcst_lead_sec,trim(outfile),.true.,'z1','z0')
    endif

    call ifi_check(algo%discard_hybrid_level_vars(),'discard_hybrid_level_vars() failed')

    call ifi_check(algo%derive_fields(),'derive_fields() failed')

    call ifi_check(algo%get_derived_vars(derived_vars),'get_derived_vars()')
    if(write_ifi_debug_files) then
      write(outfile,308) 'derived_vars',me
      call write_fip_output(derived_vars,fcst_lead_sec,trim(outfile),.false.,'z1','z0')
    endif

    call ifi_check(algo%run_fip_algo(),'run_fip_algo() failed')

    call ifi_check(algo%get_fip_algo_vars(pres_vars),'get_fip_algo_vars()')
    if(write_ifi_debug_files) then
      write(outfile,308) 'fip_algo_vars',me
      call write_fip_output(pres_vars,fcst_lead_sec,trim(outfile),.false.,'z1','z0')
    endif

    call ifi_check(algo%discard_pressure_level_vars(),'discard_pressure_level_vars() failed')
    call ifi_check(algo%discard_derived_vars(),'discard_derived_vars() failed')
    call ifi_check(algo%pressure_to_flight(),'pressure_to_flight() failed')

    call ifi_check(algo%get_flight_vars(flight_vars),'get_flight_vars()')
    if(write_ifi_debug_files) then
      write(outfile,308) 'flight_vars',me
      call write_fip_output(flight_vars,fcst_lead_sec,trim(outfile),.false.,'z1','z0')
    endif

    call ifi_check(algo%discard_fip_algo_vars(),'discard_fip_algo_vars() failed')
    call ifi_check(algo%make_icing_category(),'make_icing_category() failed')

    ! Get the final output fields:

    call ifi_check(algo%get_cat_vars(cat_vars),'get_cat_vars()')

    if(write_ifi_debug_files) then
      write(outfile,308) 'cat_vars',me
      call write_fip_output(cat_vars,fcst_lead_sec,trim(outfile),.false.,'z1','z0')
    endif
    
    call send_data(cat_vars,'ICE_PROB',1007)
    call send_data(cat_vars,'SLD',1008)
    call send_data(cat_vars,'ICE_SEV_CAT',1009)
    call send_data(cat_vars,'WMO_ICE_SEV_CAT',1010)
    
    ! When this subroutine ends, a Fortran-2003-compliant compiler
    ! will free all memory IFI uses, by calling the destructors (final
    ! routines) for cat_vars, hybr_vars, algo, and config.

  contains
    
    subroutine add_ijk_var(ifi_name,upp_name,upp_var)
      implicit none
      character(len=*), intent(in) :: ifi_name,upp_name
      real, intent(in) :: upp_var(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U,LM)
      real(kind=ifi_real_t) :: ifi_var(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U,LM)
      integer i,j,k

!$OMP PARALLEL DO COLLAPSE(2)
      do k=1,lm
        do j=jsta_2l,jend_2u
          do i=ista_2l,iend_2u
            ifi_var(i,j,k) = upp_var(i,j,lm-k+1)
          enddo
        enddo
      enddo
      
      call ifi_check(hybr_vars%add_ijk_var(ifi_name,int(ista_2l,c_int64_t),int(iend_2u,c_int64_t),&
           int(jsta_2l,c_int64_t),int(jend_2u,c_int64_t),&
           int(1,c_int64_t),int(lm,c_int64_t),ifi_var), &
           'could not send '//ifi_name//' ('//upp_name//') to IFI')
    end subroutine add_ijk_var

  end subroutine run_ifi

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
    integer :: i,type,iret
    real(kind=ifi_real_t), parameter :: zero = 0

    print *,'find range begin'

    if(missing_value_is_set) then
      epsilon = abs(missing_value)*1e-4
      if(.not. epsilon+1>epsilon) then
        epsilon=1
      endif
    endif

    ! Initialize to out-of-bounds values so they'll stay that way if no valid values are found:
    minv = 1e30
    maxv = -1e30
    
    ! Find the min and max values in the array:
    if(missing_value_is_set) then
      !$OMP PARALLEL DO REDUCTION(min:minv) REDUCTION(max:maxv)
      do i=1,count
        if(.not. abs(var(i)-missing_value)>epsilon .and. var(i)<9e9 .and. var(i)>-9e9) then
           minv=min(minv,var(i))
           maxv=max(maxv,var(i))
        endif
      end do
    else
      !$OMP PARALLEL DO REDUCTION(min:minv) REDUCTION(max:maxv)
      do i=1,count
        if(var(i)+1>var(i)) then
          if(var(i)<9e9 .and. var(i)>-9e9) then
             minv=min(minv,var(i))
             maxv=max(maxv,var(i))
          endif
        endif
      end do
    endif

    if(ifi_real_t==c_double) then
      type=MPI_REAL8
    else
      type=MPI_REAL4
    endif

    ! call MPI_Allreduce(minv,global_minv,1,type,MPI_MIN,mpi_comm_comp,iret)
    ! call MPI_Allreduce(maxv,global_maxv,1,type,MPI_MAX,mpi_comm_comp,iret)
    global_minv=minv
    global_maxv=maxv
    
    ! If min or max are inf, -inf, or NaN, assume all values are missing.
    all_missing = global_minv<-9e9 .or. global_minv>9e9 .or. global_maxv<-9e9 .or. global_maxv>9e9

    print *,'range min,max = ',global_minv,global_maxv

    if(all_missing) then
      min_not_miss=-1
      max_not_miss=1
    else
      min_not_miss=global_minv
      max_not_miss=global_maxv
    endif

    print *,'find range end'
  end subroutine find_range

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nc_check(code,file,message)
    use netcdf
    use mpi
    implicit none
    integer, intent(in) :: code
    character(*), intent(in) :: file, message
    integer :: i, ierr
    if(code/=0) then
      write(0,93) trim(file),trim(message),trim(nf90_strerror(code)),code
      call mpi_abort(mpi_comm_world,1,ierr)
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
    integer :: count

    character(len=200) :: varname,nextname

    integer :: dimids(4), ncid, ivar, id, nvars, nx0, ny0, nz0, ntime, dims(4)
    type(var_info),target :: var_data(maxvars)

    if(me==0) then
      write(*,'(A,A)') trim(output_file),': writing IFI debug data'
    endif

    call ifi_check(ifi_data%get_dims(ids,ide, jds,jde, kds,kde, ips,ipe, jps,jpe, kps,kpe),"get_dims")

    nx0=ide-ids+1
    ny0=jde-jds+1
    nz0=kde-kds+1
    ntime=1

    dims = (/nx0,ny0,nz0,ntime/)

    varname=' ' ! special value indicating "start at the beginning"

    ivar=0
    do while(ivar<maxvars)
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
      var_data(ivar)%data => read_var(var_data(ivar),trim(varname),var_data(ivar)%should_dealloc)
      if(me==0) then
        call set_dims(var_data(ivar),dims,dimids)
        var_data(ivar)%varid = def_var(var_data(ivar),rename)
      endif
    enddo

    if(me/=0) then
      ! Ranks that are not writing are done now
      return
    endif

    nvars=ivar

    if(nvars<1) then
      write(*,"(A,A)") output_file,': no data to write!'
      return
    endif

    call nc_check(nf90_enddef(ncid),output_file,'nf90_enddef')

    write(*,*) 'before write loop ',me
    do ivar=1,nvars
24     format('  put var ',A)
       print 24,trim(var_data(ivar)%varname)

       call find_range(var_data(ivar)%data,var_data(ivar)%count, &
            var_data(ivar)%missing_value_is_set,var_data(ivar)%missing_value,&
            var_data(ivar)%min_not_miss,var_data(ivar)%max_not_miss, &
            var_data(ivar)%all_missing)
       call write_var(var_data(ivar))
      if(var_data(ivar)%should_dealloc) then
        deallocate(var_data(ivar)%data)
      endif
      nullify(var_data(ivar)%data)
    enddo

    call nc_check(nf90_close(ncid),output_file,"nf90_close")
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

      call nc_check(nf90_put_att(ncid,var%varid,"min_value",var%min_not_miss), &
           output_file,"nf90_put_att "//trim(var%outname)//" min_value")

      call nc_check(nf90_put_att(ncid,var%varid,"max_value",var%max_not_miss), &
           output_file,"nf90_put_att "//trim(var%outname)//" max_value")

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
      use mpi
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
      integer :: global_ilen,global_jlen,global_klen,global_index, ierr
      type(c_ptr) :: global_cptr,local_cptr

      local_data_1D => ifi_data%get_data(trim(varname),var%missing_value_is_set,var%missing_value, &
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
        call mpi_abort(mpi_comm_world,1,ierr)
      end if

      count=global_ilen*global_jlen*global_klen
      if(count<=0) then
        write(0,39) trim(varname),count
39      format("IFI variable ",A," had no data (size=",I0,")")
        call mpi_abort(mpi_comm_world,1,ierr)
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
        call mpi_abort(mpi_comm_world,1,ierr)
      endif

      if(global_jlen /= ny0) then
        write(0,40) trim(varname),'global_jlen',global_jlen,ny0
        call mpi_abort(mpi_comm_world,1,ierr)
      endif

      if(global_klen/=1 .and. global_klen /= nz0) then
        write(0,40) trim(varname),'global_klen',global_klen,nz0
        call mpi_abort(mpi_comm_world,1,ierr)
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

#endif

end module upp_ifi_mod
