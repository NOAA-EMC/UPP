!> @file
!> @brief Subroutines that exchange one halo row.
!>
!> These routines are to exchange one halo row.
!> 
!> @param[in] A Array to have halos exchanged.
!> @param[out] A Array with halos exchanged.
!>
!> ### Program history log:
!> Date       | Programmer          | Comments
!> -----------|---------------------|----------
!> 2000-01-06 | Jim Tuccillo        | Initial
!> 2021-06-01 | George Vandenberghe | 2D decomposition             
!>
!> @note The 1st line is an inlined compiler directive that turns off -qcheck
!> during compilation, even if it's specified as a compiler option in the
!> makefile (Tuccillo, personal communication;  Ferrier, Feb '02).
!>
!> @author Jim Tuccillo IBM @date 2000-01-06
      SUBROUTINE EXCH_c_float(A)

      use ctlblk_mod, only: num_procs, jend, iup, jsta, idn, mpi_comm_comp, im,&
          icoords,ibcoords,bufs,ibufs,me,numx, &  
          jsta_2l, jend_2u,ileft,iright,ista_2l,iend_2u,ista,iend,jm,modelname
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use iso_c_binding, only: c_sizeof, c_float
      use mpi
      implicit none
!
      real(kind=c_float),intent(inout) :: a ( ista_2l:iend_2u,jsta_2l:jend_2u )
      real(kind=c_float), allocatable :: coll(:), colr(:)
      integer, allocatable :: icoll(:), icolr(:)
      integer status(MPI_STATUS_SIZE)
      integer ierr, jstam1, jendp1,j
      integer size,ubound,lbound
      integer msglenl, msglenr
      integer i,ii,jj, ibl,ibu,jbl,jbu,icc,jcc 
      integer iwest,ieast
      integer ifirst
      integer mpi_kind

      logical, parameter :: checkcoords = .false.
      
      data ifirst/0/
      allocate(coll(jm))
      allocate(colr(jm))
      allocate(icolr(jm))
      allocate(icoll(jm)) 
      ibl=max(ista-1,1)
      ibu=min(im,iend+1)
      jbu=min(jm,jend+1)
      jbl=max(jsta-1,1)
!

!     write(0,*) 'mype=',me,'num_procs=',num_procs,'im=',im,'jsta_2l=', &
!             jsta_2l,'jend_2u=',jend_2u,'jend=',jend,'iup=',iup,'jsta=', &
!             jsta,'idn=',idn
      if ( num_procs <= 1 ) return
!
!  for global model apply cyclic boundary condition

      IF(MODELNAME == 'GFS') then
        if(ifirst .le.  0 .and. me .eq. 0) print *,'  CYCLIC BC APPLIED'
        if(ileft .eq. MPI_PROC_NULL)  iwest=1         ! get eastern bc from western boundary of full domain
        if(iright .eq. MPI_PROC_NULL)  ieast=1        ! get western bc from eastern boundary of full domain
        if(ileft .eq. MPI_PROC_NULL)  ileft=me+(numx-1)
        if(iright .eq. MPI_PROC_NULL)  iright=(me-numx) +1 
      endif

      jstam1 = max(jsta_2l,jsta-1)                        ! Moorthi

!  send last row to iup's first row+  and receive first row-  from idn's last row

      call mpi_sendrecv(a(ista,jend),iend-ista+1,MPI_REAL4,iup,1,             &
     &                  a(ista,jstam1),iend-ista+1,MPI_REAL4,idn,1,           &
     &                  MPI_COMM_COMP,status,ierr)

      if ( ierr /= 0 ) then
        print *, ' problem with first  sendrecv in exch, ierr = ',ierr
        stop 6661
      endif

      if (checkcoords) then
        if(ifirst .le. 0) then !IFIRST ONLY
          call mpi_sendrecv(ibcoords(ista,jend),iend-ista+1,MPI_INTEGER,iup,1,            &
         &                  ibcoords(ista,jstam1),iend-ista+1,MPI_INTEGER,idn,1,          &
         &                  MPI_COMM_COMP,status,ierr)
          if ( ierr /= 0 ) then
            print *, ' problem with second sendrecv in exch, ierr = ',ierr
            stop 7661
          endif
          do i=ista,iend
            ii=ibcoords(i,jstam1)/10000
            jj=ibcoords(i,jstam1)-(ii*10000)
            if(ii .ne. i .or. jj .ne. jstam1 ) print *,' GWVX JEXCH CHECK FAIL ',ii,jj,ibcoords(i,jstam1),i
          end do
        endif !IFIRST
      endif !checkcoords

!  build the I columns to send and receive

      msglenl=jend-jsta+1
      msglenr=jend-jsta+1
      if(iright .lt. 0) msglenr=1
      if(ileft .lt. 0) msglenl=1

      do j=jsta,jend
       coll(j)=a(ista,j)
      end do

      call mpi_barrier(mpi_comm_comp,ierr)

!  send first col    to  ileft  last  col+  and receive last  col+ from ileft first col 

      call mpi_sendrecv(coll(jsta),msglenl    ,MPI_REAL4,ileft,1,           &
    &                  colr(jsta),msglenr    ,MPI_REAL4,iright,1,           &
    &                  MPI_COMM_COMP,status,ierr)

      if ( ierr /= 0 ) then
        print *, ' problem with third  sendrecv in exch, ierr = ',ierr
        stop 6662
      endif

      if(ifirst .le. 0) then ! IFIRST ONLY
         call mpi_sendrecv(icoll(jsta),msglenl    ,MPI_INTEGER,ileft,1,      & 
        &                  icolr(jsta),msglenr    ,MPI_INTEGER,iright,1,     & 
        &                  MPI_COMM_COMP,status,ierr)
        if ( ierr /= 0 ) then
          print *, ' problem with fourth sendrecv in exch, ierr = ',ierr
          stop 7662
        endif
      endif !IFIRST

      if(iright .ge. 0) then
        do j=jsta,jend
          a(iend+1,j)=colr(j)
          if(checkcoords) then
            if(ifirst .le. 0) then !IFIRST ONLY
              ibcoords(iend+1,j)=icolr(j) 
               ii=ibcoords(iend+1,j)/10000
               jj=ibcoords( iend+1,j)-(ii*10000)
               if( j .ne. jj .or. ii .ne. iend+1 .and. ii .ne. im .and. ii .ne. 1) &
               write(0,921) j,iend+1,ii,jj,ibcoords(iend+1,j),'IEXCH COORD FAIL j,iend+1,ii,jj,ibcoord '
            endif !IFIRST
          endif !checkcoords
        end do
      endif  ! for iright

 921  format(5i10,a50)

!     print *,'mype=',me,'in EXCH, after first mpi_sendrecv'

      if ( ierr /= 0 ) then
        print *, ' problem with fifth sendrecv in exch, ierr = ',ierr
        stop 6663
      end if
      jendp1 = min(jend+1,jend_2u)                          ! Moorthi

!GWV.  change from full im row exchange to iend-ista+1 subrow exchange,

      do j=jsta,jend
        colr(j)=a(iend,j)
      end do

!  send first row to idown's last row+  and receive last row+  from iup's first row

      call mpi_sendrecv(a(ista,jsta),iend-ista+1,MPI_REAL4,idn,1,             &
     &                  a(ista,jendp1),iend-ista+1,MPI_REAL4,iup,1,           &
     &                  MPI_COMM_COMP,status,ierr)
      if ( ierr /= 0 ) then
        print *, ' problem with sixth  sendrecv in exch, ierr = ',ierr
        stop 6664
      endif

      if (checkcoords) then
        if (ifirst .le. 0) then
          call mpi_sendrecv(ibcoords(ista,jsta),iend-ista+1,MPI_INTEGER,idn,1, &
     &                  ibcoords(ista,jendp1),iend-ista+1,MPI_INTEGER,iup,1,   &
     &                  MPI_COMM_COMP,status,ierr)
          if ( ierr /= 0 ) then
            print *, ' problem with seventh sendrecv in exch, ierr = ',ierr
            stop 7664
          endif
        endif ! IFIRST
      endif ! checkcoords

!  send last col to iright first col- and receive first col- from ileft last col 

      call mpi_sendrecv(colr(jsta),msglenr    ,MPI_REAL4,iright,1 ,         &
    &                  coll(jsta),msglenl    ,MPI_REAL4,ileft ,1,           &
    &                  MPI_COMM_COMP,status,ierr)

      if ( ierr /= 0 ) then
         print *, ' problem with eighth sendrecv in exch, ierr = ',ierr
         stop 6665
      endif

      if (ifirst .le. 0) then
        call mpi_sendrecv(icolr(jsta),msglenr    ,MPI_integer,iright,1 ,    &
        &                  icoll(jsta),msglenl    ,MPI_integer,ileft ,1,    &
        &                  MPI_COMM_COMP,status,ierr)
        if ( ierr /= 0 ) then
          print *, ' problem with ninth  sendrecv in exch, ierr = ',ierr
          stop 7665
        endif
      endif !IFIRST

      if(ileft .ge. 0) then
        do j=jsta,jend
          a(ista-1,j)=coll(j)
          if(checkcoords) then
            if(ifirst .le. 0) then
              ibcoords(ista-1,j)=icoll(j) 
              ii=ibcoords(ista-1,j)/10000
              jj=ibcoords( ista-1,j)-(ii*10000)
              if( j .ne. jj .or. ii .ne. ista-1 .and. ii .ne. im .and. ii .ne. 1) &
                write(0,921) j,ista-1,ii,jj,ibcoords(ista-1,j),'EXCH COORD FAIL j,ista-1,ii,jj,ibcoord '
            endif !IFIRST
          endif !checkcoords
        end do
      endif

!  interior check

      if(checkcoords) then
        if(ifirst .le. 0) then
          do j=jsta,jend
            do i=ista,iend
              ii=ibcoords(i,j)/10000
              jj=ibcoords( i,j)-(ii*10000)
              if(ii .ne. i .or. jj .ne. j) write(0,151) 'INFAILED IJ ',i,j,ibcoords(i,j),ibl,jbl,ibu,jbu
            end do
          end do
        endif !IFIRST
      endif !checkcoords

!!   corner points.   After the exchanges above, corner points are replicated in
!    neighbour halos so we can get them from the neighbors rather than
!    calculating more corner neighbor numbers  
! A(ista-1,jsta-1) is in the ileft     a(iend,jsta-1) location 
! A(ista-1,jend+1) is in the ileft     a(iend,jend+1) location 
! A(iend+1,jsta-1) is in the iright     a(ista,jsta-1) location 
! A(iend+1,jend+1) is in the iright    a(ista,jend+1) location 
!GWVx      ibl=max(ista-1,1)
!GWVx      ibu=min(im,iend+1)

      ibl=max(ista-1,1)
      ibu=min(im,iend+1)
      if(modelname == 'GFS') then
        ibl=max(ista-1,0)
        ibu=min(im+1,iend+1)
      endif

      jbu=min(jm,jend+1)
      jbl=max(jsta-1,1)

      call mpi_sendrecv(a(iend,jbl   ),1,    MPI_REAL4,iright,1 ,            &
    &                  a(ibl   ,jbl   ),1,   MPI_REAL4,ileft ,1,           &
    &                  MPI_COMM_COMP,status,ierr)
      if ( ierr /= 0 ) then
        print *, ' problem with tenth  sendrecv in exch, ierr = ',ierr
        stop 6771
      endif

      call mpi_sendrecv(a(iend,jbu   ),1,    MPI_REAL4,iright,1 ,            &
    &                  a(ibl   ,jbu   ),1,   MPI_REAL4,ileft ,1,           &
    &                  MPI_COMM_COMP,status,ierr)

      if ( ierr /= 0 ) then
        print *, ' problem with  eleventh sendrecv in exch, ierr = ',ierr
        stop 6772
      endif

      call mpi_sendrecv(a(ista,jbl   ),1,    MPI_REAL4,ileft ,1,            &
    &                  a(ibu   ,jbl   ),1,   MPI_REAL4,iright,1,           &
    &                  MPI_COMM_COMP,status,ierr)

      if ( ierr /= 0 ) then
        print *, ' problem with twelft sendrecv in exch, ierr = ',ierr
        stop 6773
      endif

      call mpi_sendrecv(a(ista,jbu   ),1,    MPI_REAL4,ileft ,1 ,            &
    &                  a(ibu   ,jbu   ),1,   MPI_REAL4,iright,1,           &
    &                  MPI_COMM_COMP,status,ierr)

      if ( ierr /= 0 ) then
        print *, ' problem with thirteenth  sendrecv in exch, ierr = ',ierr
        stop 6774
      endif

 139  format(a20,5(i10,i6,i6,'<>'))

      if(checkcoords) then
        if(ifirst .le. 0) then
          call mpi_sendrecv(ibcoords(iend,jbl   ),1    ,MPI_INTEGER,iright,1 ,          &
        &                  ibcoords(ibl   ,jbl   ),1   ,MPI_INTEGER,ileft ,1,           &
        &                  MPI_COMM_COMP,status,ierr)

          call mpi_sendrecv(ibcoords(iend,jbu   ),1    ,MPI_INTEGER,iright,1,           &
        &                  ibcoords(ibl   ,jbu   ),1   ,MPI_INTEGER,ileft ,1,           &
        &                  MPI_COMM_COMP,status,ierr)
          call mpi_sendrecv(ibcoords(ista,jbl   ),1    ,MPI_INTEGER,ileft ,1,           &
        &                  ibcoords(ibu   ,jbl   ),1   ,MPI_INTEGER,iright,1,           &
        &                  MPI_COMM_COMP,status,ierr)
          call mpi_sendrecv(ibcoords(ista,jbu   ),1    ,MPI_INTEGER,ileft ,1 ,          &
        &                  ibcoords(ibu   ,jbu   ),1   ,MPI_INTEGER,iright,1,           &
                           MPI_COMM_COMP,status,ierr)

!    corner check for coordnates

          icc=ibl
          jcc=jbl
          ii=ibcoords(icc,jcc)/10000
          jj=ibcoords(icc,jcc)-(ii*10000)

          if(ii .ne. icc .and. icc .ne. 0) write(0,151) ' CORNER FAILI ilb  ll ',icc,jcc,ibcoords(icc,jcc),ii,jj
          if( jj .ne. jcc)  write(0,151) ' CORNER FAILJ ilb  ll ',icc,jcc,ibcoords(icc,jcc),ii,jj

          icc=ibu
          jcc=jbl
          ii=ibcoords(icc,jcc)/10000
          jj=ibcoords(icc,jcc)-(ii*10000)
          if(ii .ne. icc .and. icc .ne. im+1 ) write(0,151) ' CORNER FAILI ilb ul  ',icc,jcc,ibcoords(icc,jcc),ii,jj
          if( jj .ne. jcc  ) write(0,151) ' CORNER FAILJ ilb ul  ',icc,jcc,ibcoords(icc,jcc),ii,jj

          icc=ibu
          jcc=jbu
          ii=ibcoords(icc,jcc)/10000
          jj=ibcoords(icc,jcc)-(ii*10000)
          if(ii .ne. icc  .and. icc .ne. im+1) write(0,151) ' CORNER FAILI ilb uu  ',icc,jcc,ibcoords(icc,jcc),ii,jj
          if( jj .ne. jcc  ) write(0,151) ' CORNER FAILJ ilb ul  ',icc,jcc,ibcoords(icc,jcc),ii,jj

          icc=ibl
          jcc=jbu
          ii=ibcoords(icc,jcc)/10000.
          jj=ibcoords(icc,jcc)-(ii*10000)
          if(ii .ne. icc  .and. icc .ne. 0 ) write(0,151) ' CORNER FAILI ilb lu  ',icc,jcc,ibcoords(icc,jcc),ii,jj
          if( jj .ne. jcc  ) write(0,151) ' CORNER FAILJ ilb ul  ',icc,jcc,ibcoords(icc,jcc),ii,jj

!         if(ileft .ge. 0) then
!119  format(' GWX LEFT EXCHANGE ileft,me,ibcoords(ista-1,jend+1),ibcoords(ista-1,jend-1),ista-1,jend-1,jend+1', &
!     10i10)                                                                           
!         endif

!        if(iright .ge. 0) then
!!         write(0,129) iright,me,ibcoords(ista+1,jend+1),ibcoords(ista+1,jend-1),ista-1,jend-1,jend+1 !GWVX
!129  format(' GWX RIGHT  EXCHANGE iright,me,ibcoords(ista+1,jend+1),ibcoords(ista-1,jend+1),ista-1,jend-1,jend+1', &
!     10i10)                                                                           
!        endif

!  interior check

          do j=jsta,jend
            do i=ista,iend
              ii=ibcoords(i,j)/10000
              jj=ibcoords( i,j)-(ii*10000)
              if(ii .ne. i .or. jj .ne. j) write(0,151) 'GWVX FAILED IJ ',i,j,ibcoords(i,j),ibl,jbl,ibu,jbu
            end do
          end do 

 151  format(a70,10i10)

! bounds check
! first check top and bottom halo rows

          j=jbu 
          do i=ista,iend
            ii=ibcoords(i,j)/10000
            jj=ibcoords( i,j)-(ii*10000)
            if(ii .ne. i .or. jj .ne. j) write(0,151) 'GWVX FAILEDI JBU IJ ',i,j,ibcoords(i,j),ibl,jbl,ibu,jbu
          end do

          j=jbl
          do i=ista,iend
            ii=ibcoords(i,j)/10000
            jj=ibcoords( i,j)-(ii*10000)
            if(ii .ne. i .or. jj .ne. j) write(0,151) 'GWVX FAILEDI JBL IJ ',i,j,ibcoords(i,j),ibl,jbl,ibu,jbu
          end do

! second and last, check left and right halo columns

          i=ibl
          do j=jsta,jend
            ii=ibcoords(i,j)/10000
            jj=ibcoords( i,j)-(ii*10000)
            if(ii .ne. i .and. ii .ne. im  .or. jj .ne. j) write(0,151) 'GWVX FAILED IBL IJ ',ii,i,j,ibcoords(i,j),ibl,jbl,ibu,jbu
          end do

          i=ibu
          do j=jsta,jend
            ii=ibcoords(i,j)/10000
            jj=ibcoords( i,j)-(ii*10000)
            if(ii .ne. i .and. ii .ne. 1  .or. jj .ne. j) write(0,151) 'GWVX FAILED IBU ii i j ibcoords ibl,jbl,ibu,jbu',ii,i,j,ibcoords(i,j),ibl,jbl,ibu,jbu
          end do

          if(me .eq. 0) write(0,*) '  IFIRST CHECK'

        endif ! IFIRST
      endif !checkcoords

! end halo checks 
      if ( ierr /= 0 ) then
         print *, ' problem with second sendrecv in exch, ierr = ',ierr
         stop
      end if
      call mpi_barrier(mpi_comm_comp,ierr)
      ifirst=ifirst+1
      end
