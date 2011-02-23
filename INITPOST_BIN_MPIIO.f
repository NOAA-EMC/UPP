      SUBROUTINE INITPOST_BIN_MPIIO
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    INITPOST    INITIALIZE POST FOR RUN
!   PRGRMMR: RUSS TREADON    ORG: W/NP2      DATE: 93-11-10
!     
! ABSTRACT:  THIS ROUTINE INITIALIZES CONSTANTS AND
!   VARIABLES AT THE START OF AN ETA MODEL OR POST 
!   PROCESSOR RUN.
!
!   THIS ROUTINE ASSUMES THAT INTEGERS AND REALS ARE THE SAME SIZE
!   .     
!     
! PROGRAM HISTORY LOG:
!   93-11-10  RUSS TREADON - ADDED DOCBLOC
!   98-05-29  BLACK - CONVERSION OF POST CODE FROM 1-D TO 2-D
!   99-01 20  TUCCILLO - MPI VERSION
!   01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-06-19  MIKE BALDWIN - WRF VERSION
!   02-08-15  H CHUANG - UNIT CORRECTION AND GENERALIZE PROJECTION OPTIONS
!   02-10-31  H CHUANG - MODIFY TO READ WRF BINARY OUTPUT
!     
! USAGE:    CALL INIT
!   INPUT ARGUMENT LIST:
!     NONE     
!
!   OUTPUT ARGUMENT LIST: 
!     NONE
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       NONE
!     LIBRARY:
!       COMMON   - CTLBLK
!                  LOOKUP
!                  SOILDEPTH
!
!    
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
      use vrbls3d
      use vrbls2d
      use soil
      use masks
      use kinds, only             : i_llong
      use params_mod
      use lookup_mod
      use ctlblk_mod
      use gridspec_mod
      use wrf_io_flags_mod
!
      implicit none
!
!     INCLUDE/SET PARAMETERS.
      INCLUDE "mpif.h"
!
      character(len=31) :: VarName
      integer :: Status
      character startdate*80,SysDepInfo*80
! 
!     NOTE: SOME INTEGER VARIABLES ARE READ INTO DUMMY ( A REAL ). THIS IS OK
!     AS LONG AS REALS AND INTEGERS ARE THE SAME SIZE.
!
!     ALSO, EXTRACT IS CALLED WITH DUMMY ( A REAL ) EVEN WHEN THE NUMBERS ARE
!     INTEGERS - THIS IS OK AS LONG AS INTEGERS AND REALS ARE THE SAME SIZE.

      INTEGER IDATE(8),JDATE(8)
!
!     DECLARE VARIABLES.
!     
      REAL SLDPTH2(NSOIL)
      REAL RINC(5)
      REAL DUMMY ( IM, JM )
      REAL DUMMY2 ( IM, JM ),MSFT(IM,JM)
      INTEGER IDUMMY ( IM, JM )
      REAL, ALLOCATABLE::  DUM3D_IKJ ( :,:,: ), DUM3D_IKJ2(:,:,:)
      real, allocatable::  pvapor(:,:)
      real, allocatable::  pvapor_orig(:,:)
      REAL, ALLOCATABLE :: thv(:,:,:)      

      character*132, allocatable :: datestr_all(:)
      character*132, allocatable :: varname_all(:)
      integer, allocatable       :: domainend_all(:,:)
      integer, allocatable       :: start_block(:)
      integer, allocatable       :: end_block(:)
      integer, allocatable       :: start_byte(:)
      integer, allocatable       :: end_byte(:)
      integer(kind=i_llong), allocatable           :: file_offset(:)
      integer this_offset, this_length

      character*80     :: titlestring

      real :: dumcst
      real :: tmp
      real :: pvapornew, qmean, tsph, tlmh, rho, dz
      integer :: itmp, iret1, iret2, irtn, imn, iyear, iday, istatus, ioutcount, i, j, l, je, js, iunit
      integer :: ll,jj,ii,n, iret, jev, index, nrecs, igdout, ierr
      integer :: nsrfc,nrdlw,nrdsw,nheat,nclod
      integer :: k1,k
      real :: ZSF,ZPBLTOP
!
!***********************************************************************
!     START INIT HERE.
!
      WRITE(6,*)'INITPOST:  ENTER INITPOST'
!     
      gridtype='A'
!     
!     STEP 1.  READ MODEL OUTPUT FILE
!
!
!***
!
! LMH always = LM for sigma-type vert coord
! LMV always = LM for sigma-type vert coord

       do j = jsta_2l, jend_2u
        do i = 1, im
            LMV ( i, j ) = lm
            LMH ( i, j ) = lm
        end do
       end do


! HTM VTM all 1 for sigma-type vert coord

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            HTM ( i, j, l ) = 1.0
            VTM ( i, j, l ) = 1.0
        end do
       end do
      end do
!
!      DateStr = '2002-03-05_18:00:00'
!  how do I get the filename?
         call ext_int_ioinit(SysDepInfo,Status)
          print*,'called ioinit', Status
         call ext_int_open_for_read( trim(fileName), 0, 0, " ",        &
     &  DataHandle, Status)
          print*,'called open for read', Status
       if ( Status /= 0 ) then
         print*,'error opening ',fileName, ' Status = ', Status ; stop
       endif
! get date/time info
!  this routine will get the next time from the file, not using it
      print *,'DateStr before calling ext_int_get_next_time=',DateStr
!      call ext_int_get_next_time(DataHandle, DateStr, Status)
      print *,'DateStri,Status,DataHandle = ',DateStr,Status,DataHandle

!  The end j row is going to be jend_2u for all variables except for V.
      JS=JSTA_2L
      JE=JEND_2U
      IF (JEND_2U.EQ.JM) THEN
       JEV=JEND_2U+1
      ELSE
       JEV=JEND_2U
      ENDIF
!
      call ext_int_get_dom_ti_char(DataHandle,'TITLE',titlestring, status )
        print*,'TITLE= ',trim(titlestring)
!
! Getting start time
      call ext_int_get_dom_ti_char(DataHandle,'START_DATE',startdate,   &
        status )
        print*,'startdate= ',startdate
      jdate=0
      idate=0
      read(startdate,15)iyear,imn,iday,ihrst,imin      
 15   format(i4,1x,i2,1x,i2,1x,i2,1x,i2)
      print*,'start yr mo day hr min=',iyear,imn,iday,ihrst,imin
      print*,'processing yr mo day hr min=',idat(3),idat(1),idat(2),     &
        idat(4),idat(5)
      idate(1)=iyear
      idate(2)=imn
      idate(3)=iday
      idate(5)=ihrst
      idate(6)=imin
      SDAT(1)=imn
      SDAT(2)=iday
      SDAT(3)=iyear
      jdate(1)=idat(3)
      jdate(2)=idat(1)
      jdate(3)=idat(2)
      jdate(5)=idat(4)
      jdate(6)=idat(5)
      CALL W3DIFDAT(JDATE,IDATE,0,RINC)
      ifhr=nint(rinc(2)+rinc(1)*24.)
      ifmin=nint(rinc(3))
      print*,' in INITPOST ifhr ifmin fileName=',ifhr,ifmin,fileName
!  OK, since all of the variables are dimensioned/allocated to be
!  the same size, this means we have to be careful int getVariable
!  to not try to get too much data.  For example, 
!  DUM3D is dimensioned IM+1,JM+1,LM+1 but there might actually
!  only be im,jm,lm points of data available for a particular variable.  
! get metadata
        call ext_int_get_dom_ti_real(DataHandle,'DX',tmp                 &
     & ,1,ioutcount,istatus)
        dxval=nint(tmp)
        write(6,*) 'dxval= ', dxval
        call ext_int_get_dom_ti_real(DataHandle,'DY',tmp                 &
     & ,1,ioutcount,istatus)
        dyval=nint(tmp)
        write(6,*) 'dyval= ', dyval

        call ext_int_get_dom_ti_integer(DataHandle,'MP_PHYSICS'                 &
     & ,itmp,1,ioutcount,istatus)
        imp_physics=itmp
        print*,'MP_PHYSICS= ',imp_physics,istatus

        call ext_int_get_dom_ti_integer(DataHandle,'SF_SURFACE_PHYSICS' &
        ,itmp,1,ioutcount,istatus)
        iSF_SURFACE_PHYSICS=itmp
        print*,'SF_SURFACE_PHYSICS= ',iSF_SURFACE_PHYSICS

        call ext_int_get_dom_ti_integer(DataHandle,'CU_PHYSICS'                 &
     & ,itmp,1,ioutcount,istatus)
        icu_physics=itmp
        print*,'CU_PHYSICS= ',icu_physics

        call ext_int_get_dom_ti_real(DataHandle,'DT',tmp                 &
     &    ,1,ioutcount,istatus)
        if(istatus==0)DT=tmp
        print*,'DT= ',DT,istatus

!need to get period of time for precipitation buckets
!      call ext_int_get_dom_ti_real(DataHandle,'PREC_ACC_DT'  &
!        ,tmp,1,ioutcount,istatus)
!      prec_acc_dt=tmp
!      print*,'PREC_ACC_DT= ',prec_acc_dt

        call ext_int_get_dom_ti_real(DataHandle,'CEN_LAT',tmp                 &
     & ,1,ioutcount,istatus)
        cenlat=nint(1000.*tmp)
        write(6,*) 'cenlat= ', cenlat
        call ext_int_get_dom_ti_real(DataHandle,'CEN_LON',tmp                 &
     & ,1,ioutcount,istatus)
        cenlon=nint(1000.*tmp)
        write(6,*) 'cenlon= ', cenlon

!        if(maptype.ne.6)then
        call ext_int_get_dom_ti_real(DataHandle,'TRUELAT1',tmp                 &
     & ,1,ioutcount,istatus)
        truelat1=nint(1000.*tmp)
        write(6,*) 'truelat1= ', truelat1
        call ext_int_get_dom_ti_real(DataHandle,'TRUELAT2',tmp                 &
     & ,1,ioutcount,istatus)
        truelat2=nint(1000.*tmp)
        write(6,*) 'truelat2= ', truelat2
!        endif
	call ext_int_get_dom_ti_real(DataHandle,'STAND_LON',tmp                 &
     & ,1,ioutcount,istatus)
        STANDLON=nint(1000.*tmp)
        write(6,*) 'STANDLON= ', STANDLON

        call ext_int_get_dom_ti_integer(DataHandle,'MAP_PROJ',itmp            &
     & ,1,ioutcount,istatus)
        maptype=itmp
        write(6,*) 'maptype is ', maptype

! get 3-D variables
! closing wrf io api
      call ext_int_ioclose ( DataHandle, Status )
! start calling mpi io
      iunit=33
      call count_recs_wrf_binary_file(iunit, fileName, nrecs)
      print*,'- FILE CONTAINS ',nrecs, ' RECORDS'
      allocate (datestr_all(nrecs))
      allocate (varname_all(nrecs))
      allocate (domainend_all(3,nrecs))
      allocate (start_block(nrecs))
      allocate (end_block(nrecs))
      allocate (start_byte(nrecs))
      allocate (end_byte(nrecs))
      allocate (file_offset(nrecs))

      call inventory_wrf_binary_file(iunit, filename, nrecs,                 &
     &                datestr_all,varname_all,domainend_all,                 &
     &      start_block,end_block,start_byte,end_byte,file_offset)

      close(iunit)

      call mpi_file_open(mpi_comm_world, filename                 &
     & , mpi_mode_rdonly,mpi_info_null, iunit, ierr)
      if (ierr /= 0) then
       print*,"Error opening file with mpi io"
       stop
      end if

      print*,'im,jm,lm= ',im,jm,lm

! assign SLDPTH to be the same as eta
! jkw comment out because Pleim Xiu only has 2 layers
!         SLDPTH(1)=0.10
!         SLDPTH(2)=0.3
!         SLDPTH(3)=0.6
!         SLDPTH(4)=1.0
! Initialize soil depth to some bogus value 
! to alert user if not found in wrfout file
      do I=1,NSOIL
       SLDPTH(I) = 0.0
      end do

! or get SLDPTH from wrf output

	index=1
      IF(iSF_SURFACE_PHYSICS==3)then ! RUC LSM
       VarName='ZS'
       call retrieve_index(index,VarName,varname_all,nrecs,iret)
       call mpi_file_read_at(iunit,file_offset(index+1)   &
       ,SLLEVEL,NSOIL,mpi_real4,mpi_status_ignore, ierr)
       if (iret /= 0 .or. ierr /= 0) then
         print*,"Error reading ", VarName,"Assigned missing values"
         SLLEVEL=SPVAL
       end if
      ELSE
       VarName='DZS'
       call retrieve_index(index,VarName,varname_all,nrecs,iret)
	write(0,*) 'iret from retrieve_index: ', iret
        call mpi_file_read_at(iunit,file_offset(index+1)                 &
     & ,SLDPTH2,NSOIL,mpi_real4                                          &
     & , mpi_status_ignore, ierr)

        if (iret /= 0 .or. ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SLDPTH2=SPVAL
        end if
      END IF 
! if SLDPTH in wrf output is non-zero, then use it
      DUMCST=0.0
      DO N=1,NSOIL
       DUMCST=DUMCST+SLDPTH2(N)
      END DO
      IF(ABS(DUMCST-0.).GT.1.0E-2)THEN
       DO N=1,NSOIL
        SLDPTH(N)=SLDPTH2(N)
        print*, 'N, SLDPTH(N): ', N, SLDPTH(N)
       END DO
      END IF


!      DO N=1,NSOIL
!       IF(SLDPTH2(N) .LT. SPVAL) SLDPTH(N)=SLDPTH2(N)
!      END DO

      print*,'SLDPTH= ',(SLDPTH(N),N=1,NSOIL)

!-------------------------------------------------------

      VarName='U'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)

	allocate(DUM3D_IKJ(IM+1,LM,JM))

        call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,((im+1)*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0 .and. ierr .eq. 0) then
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im+1
            u ( i, j, l ) = dum3d_IKJ ( i, LM-L+1, j )
        end do
       end do
!  fill up UH which is U at P-points including 2 row halo
       do j = jsta_2l, jend_2u
        do i = 1, im
            UH (I,J,L) = (dum3d_IKJ(I,LM-L+1,J)+ &
     &                    dum3d_IKJ(I+1,LM-L+1,J))*0.5
        end do
       end do
      end do

      ii=im/2
      jj=(jsta+jend)/2
      ll=lm
      write(*,*) 'U,   Level, Maximum,   Minimum   single '
      DO l=1,lm
         write(*,*) l,maxval(UH (:,:,l)),minval(UH(:,:,l)),UH (ii,jj,l)
      ENDDO
	
	else
	U=SPVAL
	UH=SPVAL
	endif

	deallocate(dum3d_ikj)

!----------------------

      VarName='V'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)

	allocate(DUM3D_IKJ(IM,LM,JM+1))

       call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*(jm+1)*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0 .and. ierr .eq. 0) then
      do l = 1, lm
       do j = jsta_2l, jev-1
        do i = 1, im
            v ( i, j, l ) = dum3d_ikj ( I, LM-L+1, J )
        end do
       end do
!  fill up VH which is V at P-points including 2 row halo
       do j = jsta_2l, jend_2u
        do i = 1, im
          VH(I,J,L) = (dum3d_ikj(I,LM-L+1,J)+ &
     &                 dum3d_ikj(I,LM-L+1,J+1))*0.5
        end do
       end do
      end do
      print*,'finish reading V'
      write(*,*) 'V,   Level, Maximum,   Minimum   single '
      DO l=1,lm
         write(*,*) l,maxval(VH (:,:,l)),minval(VH(:,:,l)),VH (ii,jj,l)
      ENDDO
	else
	V=SPVAL
	VH=SPVAL
	endif

	deallocate(dum3d_ikj)

!------------------------

      VarName='W'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)

	allocate(dum3d_ikj(IM,LM+1,JM))

        call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*(lm+1)),mpi_real4 &
     & , mpi_status_ignore, ierr)


!  fill up WH which is W at P-points including 2 row halo
	if (iret .eq. 0 .and. ierr .eq. 0) then
      DO L=1,LM
        DO I=1,IM
         DO J=JSTA_2L,JEND_2U
          WH(I,J,L) = (DUM3D_IKJ(I,LM-L+1,J)+ &
     &                 DUM3D_IKJ(I,LM-L+2,J))*0.5
         ENDDO
        ENDDO
      ENDDO
      print*,'finish reading W'
      write(*,*) 'W,   Level, Maximum,   Minimum   single '
      DO l=1,lm
         write(*,*) l,maxval(WH (:,:,l)),minval(WH(:,:,l)),WH (ii,jj,l)
      ENDDO
	else
	WH=SPVAL
	endif

	deallocate(dum3d_ikj)

!------------------------
!
! reading geopotential

      VarName='PH'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)

	allocate(dum3d_ikj2(IM,LM+1,JM))

        call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ2,(im*jm*(lm+1)),mpi_real4 &
     & , mpi_status_ignore, ierr)
	
      VarName='PHB'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)

	allocate(dum3d_ikj(IM,LM+1,JM))
        call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*(lm+1)),mpi_real4 &
     & , mpi_status_ignore, ierr)

      print*,'finish reading geopotential'
! ph/phb are geopotential z=(ph+phb)/9.801
      DO L=1,LM+1
        DO I=1,IM
         DO J=JS,JE
          ZINT(I,J,L)=(DUM3D_IKJ(I,LM-L+2,J)+ &
     &                 DUM3D_IKJ2(I,LM-L+2,J))/G
         ENDDO
        ENDDO
      ENDDO
      DO L=1,LM
       DO I=1,IM
        DO J=JS,JE
         ZMID(I,J,L)=(ZINT(I,J,L+1)+ZINT(I,J,L))*0.5  ! ave of z
        ENDDO
       ENDDO
      ENDDO
!      print*,'ZMID at ',ii,jj,ll,' = ',ZMID(ii,jj,ll)      
!      print*,'ZINT at ',ii,jj,ll+1,' = ',ZINT(ii,jj,ll+1)
      write(*,*) 'PH,   Level, Maximum,   Minimum   single '
      DO l=1,lm
         write(*,*) l,maxval(ZMID(:,:,l)),minval(ZMID(:,:,l)), &
     &                       ZMID(ii,jj,l)
      ENDDO

	deallocate(dum3d_ikj, dum3d_ikj2)

!------------------------

!
! reading potential temperature
      VarName='T'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
	
	ALLOCATE(DUM3D_IKJ(IM,LM,JM))

       call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)


	if (iret .eq. 0  .and. ierr .eq. 0) then
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            t ( i, j, l ) = dum3d_ikj ( i, LM-L+1, j ) + 300.
!MEB  this is theta the 300 is my guess at what T0 is
        end do
       end do
      end do
      print*,'finish reading T'
      write(*,*) 'TH,   Level, Maximum,   Minimum   single '
      DO l=1,lm
         write(*,*) l,maxval(T(:,:,l)),minval(T(:,:,l)),T(ii,jj,l)
      ENDDO
	else
	T=SPVAL
	endif

	deallocate(dum3d_ikj)

!------------------------
!
! reading sfc pressure
      VarName='MU'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
       call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)
	if (iret .ne. 0 .or. ierr .ne. 0) then
	DUMMY=SPVAL
	endif

!------------------------------------------------------------------

      VarName='MUB'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
       call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY2,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)
	if (iret .eq. 0 .and. ierr .eq. 0) then
      DO I=1,IM
            DO J=JS,JE
                 PINT (I,J,LM+1) = DUMMY(I,J)+DUMMY2(I,J) !PSFC-PT
            ENDDO
      ENDDO
      endif
!
      write(*,*) IM,JS,JE
      write(*,*) 'PSFC-PT,   Level, Maximum,   Minimum   single '
      write(*,*) l,maxval(PINT(:,:,LM+1)),minval(PINT(:,:,LM+1)), &
     &                       PINT(ii,jj,LM+1)

!-----------------------------------

      VarName='P'
      call retrieve_index(index,VarName,varname_all,nrecs,iret1)

	ALLOCATE(DUM3D_IKJ2(IM,LM,JM))

      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ2,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

!-----------------------------------

      VarName='PB'
      call retrieve_index(index,VarName,varname_all,nrecs,iret2)

	ALLOCATE(DUM3D_IKJ(IM,LM,JM))
      ALLOCATE ( thv(IM,JSTA_2L:JEND_2U,LM) )
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret1 .eq. 0 .and. iret2 .eq. 0 .and. ierr .eq. 0) then
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            PMID(I,J,L)=DUM3D_IKJ(I,LM-L+1,J)+DUM3D_IKJ2(I,LM-L+1,J)
!            thv ( i, j, l ) = T(I,J,L)*(Q(I,J,L)*0.608+1.)
! now that I have P, convert theta to t
            t ( i, j, l ) = T(I,J,L)*(PMID(I,J,L)*1.E-5)**CAPA

        end do
       end do
      end do
	 else
	write(0,*) 'MISSING PB or P...setting PMID and T to SPVAL'
	T=SPVAL
        PMID=SPVAL
	endif
       write(*,*) lm,jsta_2l, jend_2u,im
      write(*,*) 'P,   Level, Maximum,   Minimum   single '
      DO l=1,lm
         write(*,*) l,maxval(PMID(:,:,l)),minval(PMID(:,:,l)), &
     &                       PMID(ii,jj,l)
      ENDDO
      write(*,*) 'T,   Level, Maximum,   Minimum   single '
      DO l=1,lm
         write(*,*) l,maxval(t(:,:,l)),minval(t(:,:,l)), &
     &                       t(ii,jj,l)
      ENDDO

      DO L=2,LM
         DO I=1,IM
            DO J=JSTA_2L,JEND_2U
!              PINT(I,J,L)=EXP((ALOG(PMID(I,J,L-1))+
!     &                 ALOG(PMID(I,J,L)))*0.5)  ! ave of ln p
              PINT(I,J,L)=(PMID(I,J,L-1)+PMID(I,J,L))*0.5
              ALPINT(I,J,L)=ALOG(PINT(I,J,L))

	if (I .eq. ii .and. J .eq. jj) then
	write(6,*) 'L, pint(I,J,L): ', L, pint(I,J,L)
	endif

            ENDDO
         ENDDO
      END DO

!      print*,'PINT at ',ii,jj,ll,' = ',pint(ii,jj,ll)
!      print*,'T at ',ii,jj,ll,' = ',t(ii,jj,ll)
      
	deallocate(dum3d_ikj,dum3d_ikj2)

!-------------------------------------------

      VarName='SR'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)
	if (iret .eq. 0 .and. ierr .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            SR ( i, j ) = dummy ( i, j )
        end do
       end do      
	write(6,*) 'maxval SR: ', maxval(SR)
        else
	write(0,*) 'MISSING SR set to ZERO'
	SR=0.
	endif
!-------------------------------------------

	VarName='SMSTAV'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0 .and. ierr .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
	SMSTAV( i, j )= dummy ( i, j )
        enddo
       enddo
        else
         write(0,*) VarName ,'not found....set to SPVAL'
        SMSTAV=SPVAL
        SMSTAV=0.
        endif

!-------------------------------------------

	VarName='SMSTOT'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0 .and. ierr .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
	SMSTOT( i, j )= dummy ( i, j )
        enddo
       enddo
        else
         write(0,*) VarName ,'not found....set to SPVAL'
        SMSTOT=SPVAL
        SMSTOT=0.
        endif


!-------------------------------------------
! reading 2m theta
      VarName='TH2'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)
	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            TSHLTR ( i, j ) = dummy ( i, j )
        end do
       end do
        write(*,*) ' TH2'
        write(*,*) maxval(TSHLTR),minval(TSHLTR),TSHLTR(ii,jj)

        else
	TSHLTR=SPVAL
	endif

!-------------------------------------------

! reading 10 m wind
      VarName='U10'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY2,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            U10 ( i, j ) = dummy2( i, j )
        end do
       end do
!       print*,'U10 at ',ii,jj,' = ',U10(ii,jj)
        write(*,*) ' U10'
        write(*,*) maxval(U10),minval(U10),U10(ii,jj)
	else
	U10=SPVAL
	endif

!-------------------------------------------

      VarName='V10'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY2,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jev-1
        do i = 1, im
            V10 ( i, j ) = dummy2( i, j )
        end do
       end do
!       print*,'V10 at ',ii,jj,' = ',V10(ii,jj)
        write(*,*) ' V10'
        write(*,*) maxval(V10),minval(V10),V10(ii,jj)
	else
	V10=spval
	ENDIF

!--------------------------------------------------

       do j = jsta_2l, jend_2u
        do i = 1, im
            TH10 ( i, j ) = SPVAL
	    Q10 ( i, j ) = SPVAL
        end do
       end do
   
!--------------------------------------------------

! reading water vapor mixing ratio
      VarName='QVAPOR'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)

	allocate(DUM3D_IKJ(IM,LM,JM))
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
!HC CONVERT MIXING RATIO TO SPECIFIC HUMIDITY
            q ( i, j, l ) = dum3d_ikj(i,lm-l+1,j)/ &
                           (1.0+dum3d_ikj(i,lm-l+1,j))
! now that I have T,q,P  compute omega from wh
            omga(I,J,L) = -WH(I,J,L)*pmid(i,j,l)*G/ &
                              (RD*t(i,j,l)*(1.+D608*q(i,j,l)))
            thv ( i, j, l ) = T(I,J,L)*(1.E5/PMID(I,J,L))**CAPA  &
             *(Q(I,J,L)*0.608+1.)
        end do
       end do
      end do
      write(*,*) 'Q,   Level, Maximum,   Minimum   single '
      DO l=1,lm
         write(*,*) l,maxval(q(:,:,l)),minval(q(:,:,l)),q(ii,jj,l)
      ENDDO

	else
	Q=SPVAL
	OMGA=SPVAL
	endif

	deallocate(dum3d_ikj)

!      print*,'Q at ',ii,jj,ll,' = ',Q(ii,jj,ll)
!------------------------------------------------------

!!!!!!!!!!!!!
! Pyle's fixes for ARW SLP

        allocate(pvapor(IM,jsta_2l:jend_2u))
        allocate(pvapor_orig(IM,jsta_2l:jend_2u))
        DO J=jsta,jend
        DO I=1,IM


        pvapor(I,J)=0.
       do L=1,LM
       dz=ZINT(I,J,L)-ZINT(I,J,L+1)
       rho=PMID(I,J,L)/(RD*T(I,J,L))


        if (L .le. LM-1) then
        QMEAN=0.5*(Q(I,J,L)+Q(I,J,L+1))
        else
        QMEAN=Q(I,J,L)
        endif


       pvapor(I,J)=pvapor(I,J)+G*rho*dz*QMEAN
       enddo


! test elim
!       pvapor(I,J)=0.


        pvapor_orig(I,J)=pvapor(I,J)


      ENDDO
      ENDDO

      do L=1,405
        call exch(pvapor(1,jsta_2l))
        do J=JSTA_M,JEND_M
        do I=2,IM-1

        pvapornew=AD05*(4.*(pvapor(I-1,J)+pvapor(I+1,J) &
     &                  +pvapor(I,J-1)+pvapor(I,J+1)) &
     &                  +pvapor(I-1,J-1)+pvapor(I+1,J-1) &
     &                  +pvapor(I-1,J+1)+pvapor(I+1,J+1)) &
     &                  -CFT0*pvapor(I,J)

        pvapor(I,J)=pvapornew

        enddo
        enddo
        enddo   ! iteration loop

! southern boundary
        if (JS .eq. 1) then
        J=1
        do I=2,IM-1
        pvapor(I,J)=pvapor_orig(I,J)+(pvapor(I,J+1)-pvapor_orig(I,J+1))
        enddo
        endif

! northern boundary

        if (JE .eq. JM) then
        J=JM
        do I=2,IM-1
        pvapor(I,J)=pvapor_orig(I,J)+(pvapor(I,J-1)-pvapor_orig(I,J-1))
        enddo
        endif

! western boundary
        I=1
        do J=JS,JE
        pvapor(I,J)=pvapor_orig(I,J)+(pvapor(I+1,J)-pvapor_orig(I+1,J))
        enddo

! eastern boundary
        I=IM
        do J=JS,JE
        pvapor(I,J)=pvapor_orig(I,J)+(pvapor(I-1,J)-pvapor_orig(I-1,J))
        enddo

!      DO J=Jsta,jend
!      DO I=1,IM
!              PINT(I,J,LM+1)=PINT(I,J,LM+1)+PVAPOR(I,J)
!      ENDDO
!      ENDDO

        write(6,*) 'surface pvapor field (post-smooth)'

!        deallocate(pvapor)
        deallocate(pvapor_orig)

! reading 2 m mixing ratio
      VarName='Q2'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
!HC CONVERT FROM MIXING RATIO TO SPECIFIC HUMIDITY
          IF(MODELNAME == 'RAPR')THEN
!tgs - for RR set it equal to 1st level
            QSHLTR ( i, j ) =  q ( i, j, lm )
          ELSE
            QSHLTR ( i, j ) = dummy ( i, j )/(1.0+dummy ( i, j ))
          END IF 
        end do
       end do
        else
         write(0,*) 'Q2 not found....QSHLTR set to SPVAL'
          QSHLTR=SPVAL
        endif



!!!!!!!!!!!!!
! reading cloud water mixing ratio
! Brad comment out the output of individual species for Ferrier's scheme within 
! ARW in Registry file

      qqw=0.
      qqr=0.
      qqs=0.
      qqi=0.
      qqg=0. 
      cwm=0.

      allocate(DUM3D_IKJ(IM,LM,JM))

      if(imp_physics.ne.5 .and.     &
         imp_physics.ne.85 .and.    &
         imp_physics.ne.0)then 

      VarName='QCLOUD'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
! partition cloud water and ice for WSM3 
	    if(imp_physics.eq.3)then 
             if(t(i,j,l) .ge. TFRZ)then  
              qqw ( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
	     else
	      qqi  ( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
	     end if
            else ! bug fix provided by J CASE
             qqw ( i, j, l ) = dum3d_ikj ( i, LM-L+1, j ) 
	    end if  	     
        end do
       end do
      end do
!      print*,'qqw at ',ii,jj,ll,' = ',qqw(ii,jj,ll)
	else
	QQW=SPVAL
      end if
	endif

      if(imp_physics.ne.5 .and.        &
         imp_physics.ne.85 .and.       &
         imp_physics.ne.0)then
      VarName='QRAIN'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_ikj,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
! partition rain and snow for WSM3 	
          if(imp_physics .eq. 3)then
	    if(t(i,j,l) .ge. TFRZ)then  
             qqr ( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
	    else
	     qqs ( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
	    end if
           else
            qqr ( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )  
	   end if 
        end do
       end do
      end do
	else
	QQR=SPVAL
      end if
	endif
!      print*,'qqr at ',ii,jj,ll,' = ',qqr(ii,jj,ll)

!---------------------------------------

      if(imp_physics.ne.1 .and. imp_physics.ne.3 .and.   &
         imp_physics.ne.85 .and. imp_physics.ne.5 .and.  &
         imp_physics.ne.0)then
      VarName='QICE'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            qqi ( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
        end do
       end do
      end do
	else
	QQI=SPVAL
	endif
      end if

!---------------------------------------
      
      if(imp_physics.ne.1 .and. imp_physics.ne.3 .and.    &
         imp_physics.ne.85 .and. imp_physics.ne.5 .and.   &
         imp_physics.ne.0)then
      VarName='QSNOW'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)
	if (iret .eq. 0) then
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            qqs ( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
        end do
       end do
      end do
!      print*,'qqs at ',ii,jj,ll,' = ',qqs(ii,jj,ll)
	else
	QQS=SPVAL
	endif
      end if

!---------------------------------------

      if(imp_physics.eq.2 .or. imp_physics.eq.6 .or. imp_physics.eq.8)then
      VarName='QGRAUP'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
!HC            cwm ( i, j, l ) = dum3d ( i, j, l )
!HC CONVERT MIXING RATIO TO SPECIFIC HUMIDITY
            qqg ( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
        end do
       end do
      end do
!      print*,'qqg at ',ii,jj,ll,' = ',qqg(ii,jj,ll)
	else
	QQG=SPVAL
	endif
      end if

!---------------------------------------
      
      if((imp_physics.ne.5) .and. (imp_physics.ne.85))then
!HC SUM UP ALL CONDENSATE FOR CWM
       do l = 1, lm
        do j = jsta_2l, jend_2u
         do i = 1, im
          IF(QQR(I,J,L).LT.SPVAL)THEN
           CWM(I,J,L)=QQR(I,J,L)
          END IF
          IF(QQI(I,J,L).LT.SPVAL)THEN
           CWM(I,J,L)=CWM(I,J,L)+QQI(I,J,L)
          END IF
          IF(QQW(I,J,L).LT.SPVAL)THEN
           CWM(I,J,L)=CWM(I,J,L)+QQW(I,J,L)
          END IF
          IF(QQS(I,J,L).LT.SPVAL)THEN
           CWM(I,J,L)=CWM(I,J,L)+QQS(I,J,L)
          END IF
          IF(QQG(I,J,L).LT.SPVAL)THEN
           CWM(I,J,L)=CWM(I,J,L)+QQG(I,J,L)
          END IF 
         end do
        end do
       end do
      else
       VarName='CWM'

       call retrieve_index(index,VarName,varname_all,nrecs,iret)
       call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

       if (iret .eq. 0) then
         do l = 1, lm
          do j = jsta_2l, jend_2u
           do i = 1, im
            CWM ( i, j, l ) = dum3d_ikj ( i, l, j )
! Chuang
!            CWM ( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
           end do
          end do
         end do 
	else
	 CWM=SPVAL
	endif
        print*,'sample CWM= ',l,cwm(im/2,jsta,1:lm)
        VarName='F_ICE_PHY'
        call retrieve_index(index,VarName,varname_all,nrecs,iret)
        call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

       if (iret .eq. 0) then
         do l = 1, lm
          do j = jsta_2l, jend_2u
           do i = 1, im
            F_ICE( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
           end do
          end do
         end do
        else
         F_ICE=SPVAL
        endif

        VarName='F_RAIN_PHY'
        call retrieve_index(index,VarName,varname_all,nrecs,iret)
        call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

       if (iret .eq. 0) then
         do l = 1, lm
          do j = jsta_2l, jend_2u
           do i = 1, im
            F_RAIN( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
           end do
          end do
         end do
        else
         F_RAIN=SPVAL
        endif

        VarName='F_RIMEF_PHY'
        call retrieve_index(index,VarName,varname_all,nrecs,iret)
        call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
     & , mpi_status_ignore, ierr)

       if (iret .eq. 0) then
         do l = 1, lm
          do j = jsta_2l, jend_2u
           do i = 1, im
            F_RIMEF( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
           end do
          end do
         end do
        else
         F_RIMEF=SPVAL
        endif

      endif
! after all cloud fields

      IF(MODELNAME == 'RAPR')THEN
       VarName='TKE_PBL'
      ELSE
       VarName='TKE'
      END IF
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
        call mpi_file_read_at(iunit,file_offset(index+1) &
      ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
      , mpi_status_ignore, ierr)

      if (iret .eq. 0) then
        do l = 1, lm
         do j = jsta_2l, jend_2u
          do i = 1, im
           q2( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
          end do
         end do
        end do
       else
        q2=spval
       end if 

      deallocate(dum3d_ikj)

      VarName='HTOP'
      IF(ICU_PHYSICS .EQ. 3 .or. ICU_PHYSICS .EQ. 5) VarName='CUTOP'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

      if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            HTOP( i, j ) = float(LM)-dummy(i,j)+1.0
        end do
       end do
      else
        HTOP=SPVAL
      endif 

      VarName='HBOT'
      IF(ICU_PHYSICS .EQ. 3 .or. ICU_PHYSICS .EQ. 5) VarName='CUBOT'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

      if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            HBOT( i, j ) = float(LM)-dummy(i,j)+1.0
        end do
       end do
      else
        HBOT=SPVAL
      endif

      VarName='CUPPT'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

      if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            CUPPT( i, j ) = dummy(i,j)
        end do
       end do
      else
        CUPPT=SPVAL
      endif

!-----------------------------------
   
!
! reading soil temperature
      VarName='TSLB'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)

	allocate(DUM3D_IKJ(IM,NSOIL,JM))

      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*nsoil),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
      do l = 1, nsoil
       write(0,*) 'do l ',l
       do j = jsta_2l, jend_2u
        do i = 1, im
            stc ( i, j, l ) = dum3d_ikj ( i, l, j )
! flip soil layer again because wrf soil variable vertical indexing
! is the same with eta and vertical indexing was flipped for both
! atmospheric and soil layers within getVariable
!            stc ( i, j, l ) = dum3d_ikj ( i, nsoil-l+1, j)
        end do
       end do
      end do
      print*,'STC at ',ii,jj,N,' = ',stc(ii,jj,1),stc(ii,jj,2) &
     &,stc(ii,jj,3),stc(ii,jj,4)

	else
	STC=SPVAL
	endif
!----------------------------------------------

!
!
! reading soil moisture
      VarName='SMOIS'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)

      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*nsoil),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then

      do l = 1, nsoil
       do j = jsta_2l, jend_2u
        do i = 1, im
!            smc ( i, j, l ) = dum3d_ikj ( i, nsoil-l+1, j )
            smc ( i, j, l ) = dum3d_ikj ( i, l, j )
        end do
       end do
      end do
!      print*,'SMC at ',ii,jj,N,' = ',smc(ii,jj,1),smc(ii,jj,2)
!     &,smc(ii,jj,3),smc(ii,jj,4),smc(ii,jj,5)

	else
	SMC=SPVAL
	endif

!-----------------------------------------------

      VarName='SH2O'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUM3D_IKJ,(im*jm*nsoil),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       
      do l = 1, nsoil
       do j = jsta_2l, jend_2u
        do i = 1, im
!           sh2o ( i, j, l ) = dum3d_ikj ( i, nsoil-l+1, j )
            sh2o ( i, j, l ) = dum3d_ikj ( i, l, j )
        end do
       end do
      end do 
	else
	SH2O=SPVAL
	endif
        deallocate(dum3d_ikj) 

! bitmask out high, middle, and low cloud cover
       do j = jsta_2l, jend_2u
        do i = 1, im
            CFRACH ( i, j ) = SPVAL/100.
            CFRACL ( i, j ) = SPVAL/100.
            CFRACM ( i, j ) = SPVAL/100.
        end do
       end do

! WRF ARW outputs 3D cloud cover now
      allocate(DUM3D_IKJ(im,lm,jm))
      VarName='CLDFRA'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
      ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
      , mpi_status_ignore, ierr)

      if (iret .eq. 0) then
        do l = 1, lm
         do j = jsta_2l, jend_2u
          do i = 1, im
           CFR( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
          end do
         end do
        end do
       else
        cfr=spval
       end if

! Read in extinction coefficient for aerosol at 550 nm
      VarName='EXTCOF55'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
      ,DUM3D_IKJ,(im*jm*lm),mpi_real4 &
      , mpi_status_ignore, ierr)

      if (iret .eq. 0) then
        do l = 1, lm
         do j = jsta_2l, jend_2u
          do i = 1, im
           extcof55( i, j, l ) = dum3d_ikj ( i, LM-L+1, j )
          end do
         end do
        end do
       else
        extcof55=spval
       end if


      deallocate(dum3d_ikj)
!-------------------------------------------------------

      VarName='SEAICE'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
      do j = jsta_2l, jend_2u
        do i = 1, im
            SICE( i, j ) = dummy ( i, j )
        end do
       end do
	else
	SICE=SPVAL
	endif

! --------------------------------------------------

! reading SURFACE RUNOFF 
      VarName='SFROFF'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            SSROFF ( i, j ) = dummy ( i, j )
        end do
       end do
	else
	SSROFF=SPVAL
	endif

!----------------------------------------------------------------

!
! reading UNDERGROUND RUNOFF
      VarName='UDROFF'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            BGROFF ( i, j ) = dummy ( i, j )
        end do
       end do
	else
	BGROFF=SPVAL
	endif

! reading SFC EVAPORATION
      VarName='SFCEVP'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            SFCEVP( i, j ) = dummy ( i, j )
        end do
       end do
        else
        SFCEVP=SPVAL
        endif

! reading SFC EXCHANGE COEFF
      VarName='SFCEXC'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            SFCEXC( i, j ) = dummy ( i, j )
        end do
       end do
        else
        SFCEXC=SPVAL
        endif
!--------------------------------------------------------------

! reading VEGETATION TYPE 
      VarName='IVGTYP'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,IDUMMY,(im*jm),mpi_integer4 &
     & , mpi_status_ignore, ierr)

	if (IRET .EQ. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            IVGTYP ( i, j ) = idummy ( i, j ) 
        end do
       end do 
	else
	IVGTYP=SPVAL
	endif

!------------------------------------------------------------------
       
      VarName='ISLTYP' 

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,IDUMMY,(im*jm),mpi_integer4 &
     & , mpi_status_ignore, ierr)
	if (iret .eq. 0) then
      do j = jsta_2l, jend_2u
        do i = 1, im
            ISLTYP ( i, j ) = idummy ( i, j ) 
        end do
       end do
       print*,'MAX ISLTYP=', maxval(idummy)
	else
	ISLTYP=SPVAL
	endif

      VarName='ISLOPE'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,IDUMMY,(im*jm),mpi_integer4 &
     & , mpi_status_ignore, ierr)
        if (iret .eq. 0) then
      do j = jsta_2l, jend_2u
        do i = 1, im
            ISLOPE( i, j ) = idummy ( i, j )
        end do
       end do
        else
        ISLTYP=SPVAL
        endif

!-------------------------------------------------------
       
      VarName='VEGFRA'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY2,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)
	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            VEGFRC ( i, j ) = dummy2 ( i, j )
        end do
       end do
!      print*,'VEGFRC at ',ii,jj,' = ',VEGFRC(ii,jj) 
	else
	VEGFRC=SPVAL
        endif
       
!-------------------------------------------------------

      VarName='ACSNOW'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY2,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)
        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            ACSNOW ( i, j ) = dummy2 ( i, j )
        end do
       end do
        else
        ACSNOW=SPVAL
        endif

!-------------------------------------------------------

      VarName='ACSNOM'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY2,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)
        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            ACSNOM ( i, j ) = dummy2 ( i, j )
        end do
       end do
        else
        ACSNOM=SPVAL
        endif

!-------------------------------------------------------
      
      VarName='GRDFLX'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
      do j = jsta_2l, jend_2u
        do i = 1, im
            GRNFLX(I,J) = dummy ( i, j )
        end do
       end do    
	else
	GRNFLX=SPVAL
	endif

!-------------------------------------------------------------
 
      VarName='SNOW'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)
     
	if (iret .eq. 0) then
      do j = jsta_2l, jend_2u
        do i = 1, im
            SNO ( i, j ) = dummy ( i, j )
        end do
       end do
	else
	SNO=SPVAL
	endif

!-------------------------------------------------------------
     
      VarName='CANWAT'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            CMC ( i, j ) = dummy ( i, j )
        end do
       end do
	else
	CMC=SPVAL
	endif

!-------------------------------------------------------------

      VarName='SST'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            SST ( i, j ) = dummy ( i, j )
        end do
       end do
!      print*,'SST at ',ii,jj,' = ',sst(ii,jj)      
	else
	SST=SPVAL
	endif

!-------------------------------------------------------------

      VarName='THZ0'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            THZ0 ( i, j ) = dummy ( i, j )
        end do
       end do
!      print*,'SST at ',ii,jj,' = ',sst(ii,jj)
        else
        THZ0=SPVAL
        endif

!-------------------------------------------------------------

      VarName='QZ0'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            QZ0 ( i, j ) = dummy ( i, j )
        end do
       end do
!      print*,'SST at ',ii,jj,' = ',sst(ii,jj)
        else
        QZ0=SPVAL
        endif

!-------------------------------------------------------------

      VarName='UZ0'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            UZ0 ( i, j ) = dummy ( i, j )
        end do
       end do
!      print*,'SST at ',ii,jj,' = ',sst(ii,jj)
        else
        UZ0=SPVAL
        endif
!-------------------------------------------------------------

      VarName='VZ0'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            VZ0 ( i, j ) = dummy ( i, j )
        end do
       end do
!      print*,'SST at ',ii,jj,' = ',sst(ii,jj)
        else
        VZ0=SPVAL
        endif
!-------------------------------------------------------------

      VarName='QSFC'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            QS ( i, j ) = dummy ( i, j )
        end do
       end do
!      print*,'SST at ',ii,jj,' = ',sst(ii,jj)
        else
        QS=SPVAL
        endif

!-------------------------------------------------------------

      VarName='Z0'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            Z0 ( i, j ) = dummy ( i, j )
        end do
       end do
!      print*,'SST at ',ii,jj,' = ',sst(ii,jj)
        else
        Z0=SPVAL
        endif
!-------------------------------------------------------
!    
      VarName='MAPFAC_M'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4 &
     & , mpi_status_ignore, ierr)

	if (iret .eq. 0) then
      do j = jsta_2l, jend_2u
        do i = 1, im
            MSFT ( i, j ) = dummy ( i, j ) 
        end do
       end do
	else
	MSFT=SPVAL
	endif
       

!------------------------------------------------------

! reading terrain height
      VarName='HGT'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            FIS ( i, j ) = dummy ( i, j ) * G
!            if(i.eq.80.and.j.eq.42)print*,'Debug: sample fis,zint='
!     1,dummy( i, j ),zint(i,j,lm+1)
        end do
       end do
!       print*,'FIS at ',ii,jj,ll,' = ',FIS(ii,jj)

	else
	FIS=SPVAL
	endif

!-----------------------------------------------------
        write(6,*) 'past getting of HGT'

      VarName='TSK'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)
	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
!HC            THS ( i, j ) = dummy ( i, j ) ! this is WRONG (should be theta)
!HC CONVERT SKIN TEMPERATURE TO SKIN POTENTIAL TEMPERATURE
! CHC: deriving outgoing longwave fluxes by assuming emmissitivity=1
            THS ( i, j ) = dummy ( i, j ) ! wait to convert to theta later 
            RADOT ( i, j ) = DUMMY(i,j)**4.0/STBOL    
        end do
       end do
!       print*,'THS at ',ii,jj,' = ',THS(ii,jj)

	else
	THS=SPVAL
	RADOT=SPVAL
	endif
!-----------------------------------------------
!
! reading p_top
      VarName='P_TOP'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,PT,1,mpi_real4, mpi_status_ignore, ierr)

	if (iret .eq. 0) then
      print*,'P_TOP = ',PT
      DO I=1,IM
            DO J=JS,JE
                 PINT (I,J,LM+1) = PINT (I,J,LM+1)+PT
                 PINT (I,J,LM+1) = PINT (I,J,LM+1)+PVAPOR(I,J)
                 THS ( i, j ) = THS ( i, j )  &
     &                       *(P1000/PINT(I,J,NINT(LMH(I,J))+1))**CAPA
                 PINT (I,J,1) = PT
                 ALPINT(I,J,LM+1)=ALOG(PINT(I,J,LM+1))
                 ALPINT(I,J,1)=ALOG(PINT(I,J,1))
            ENDDO
         ENDDO
!      print*,'PSFC at ',ii,jj,' = ',PINT (ii,jj,lm+1)
!      print*,'THS at ',ii,jj,' = ',THS(ii,jj)
	else
	PINT=SPVAL
	THS=SPVAL
	ALPINT=SPVAL
	endif

	deallocate(PVAPOR)


!C
!C RAINC is "ACCUMULATED TOTAL CUMULUS PRECIPITATION" 
!C RAINNC is "ACCUMULATED TOTAL GRID SCALE PRECIPITATION"

	write(6,*) 'getting RAINC'
      VarName='RAINC'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            CUPREC ( i, j ) = dummy ( i, j ) * 0.001
        end do
       end do
!       print*,'CUPREC at ',ii,jj,' = ',CUPREC(ii,jj)
	else
	CUPREC=SPVAL
	endif

!-------------------------------------------------

      write(6,*) 'getting RAINNC'
      VarName='RAINNC'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            ANCPRC ( i, j ) = dummy ( i, j )* 0.001
	    ACPREC ( i, j ) = ANCPRC(I,J)+CUPREC(I,J)
        end do
       end do
!       print*,'ANCPRC at ',ii,jj,' = ',ANCPRC(ii,jj)
	write(6,*) 'past getting RAINNC'
	else
	ACPREC=SPVAL
	endif

!-- RAINC_bucket is "ACCUMULATED CUMULUS PRECIPITATION OVER BUCKET_DT PERIODS OF TIME"
        write(6,*) 'getting PREC_ACC_C, [mm]'
!      VarName='RAINC_BUCKET'
      VarName='PREC_ACC_C'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            rainc_bucket ( i, j ) = dummy ( i, j )
        end do
       end do
        else
        rainc_bucket=SPVAL
        endif

!-- RAINNC_bucket  is "ACCUMULATED GRID SCALE  PRECIPITATION OVER BUCKET_DT PERIODS OF TIME"
        write(6,*) 'getting PREC_ACC_NC, [mm]'
!      VarName='RAINNC_BUCKET'
      VarName='PREC_ACC_NC'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            rainnc_bucket ( i, j ) = dummy ( i, j )
        end do
       end do
        else
        rainnc_bucket=SPVAL
        endif

       do j = jsta_2l, jend_2u
        do i = 1, im
            PCP_BUCKET(I,J)=rainc_bucket(I,J)+rainnc_bucket(I,J)
        end do
       end do

!-------------------------------------------------

      VarName='RAINCV'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

	if (iret .eq. 0 .and. ierr .eq. 0) then
         do j = jsta_2l, jend_2u
          do i = 1, im
             CPRATE ( i, j ) = dummy ( i, j )* 0.001
          enddo
	 enddo
        else
	 CPRATE=0.
	write(6,*) 'NO RAINCV field...CPRATE set to ZERO'
	endif

!-------------------------------------------------

      VarName='RAINNCV'
      DUMMY2=0.
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
	if (iret .eq. 0) then
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY2,(im*jm),mpi_real4, mpi_status_ignore, ierr)
	endif

         do j = jsta_2l, jend_2u
          do i = 1, im
             PREC ( i, j ) = (dummy(i,j) + dummy2(I,J))* 0.001
          enddo
	 enddo

!-------------------------------------------------

      VarName='SNOWNCV'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
        if (iret .eq. 0) then
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)
        endif

         do j = jsta_2l, jend_2u
          do i = 1, im
!-- SNOW is in [m] per time sep
             snownc ( i, j ) = dummy(i,j) * 0.001 
          enddo
         enddo

!-- RAINNC_bucket  is "ACCUMULATED GRID SCALE  PRECIPITATION OVER BUCKET_DT PERIODS OF TIME"
      VarName='SNOW_ACC_NC'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
        if (iret .eq. 0) then
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)
        endif

         do j = jsta_2l, jend_2u
          do i = 1, im
             snow_bucket ( i, j ) = dummy(i,j) 
          enddo
         enddo

      VarName='GRAUPELNCV' 
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
        if (iret .eq. 0) then
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)
        endif
        do j = jsta_2l, jend_2u
          do i = 1, im
             graupelnc ( i, j ) = dummy(i,j)
          enddo
         enddo
!-------------------------------------------------

! RSWIN can be output normally in SURFCE
       do j = jsta_2l, jend_2u
        do i = 1, im
             CZEN ( i, j ) = 1.0 
             CZMEAN ( i, j ) = CZEN ( i, j )
        end do
       end do       
       
      VarName='GLW'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            RLWIN ( i, j ) = dummy ( i, j )
        end do
       end do
	else
	RLWIN=SPVAL
	endif


! ncar wrf does not output sigt4 so make sig4=sigma*tlmh**4
       do j = jsta_2l, jend_2u
        do i = 1, im
             TLMH=T(I,J,NINT(LMH(I,J)))
             SIGT4 ( i, j ) =  5.67E-8*TLMH*TLMH*TLMH*TLMH
        end do
       end do
! Top of the atmosphere outgoing LW radiation
      VarName='OLR'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            RLWTOA ( i, j ) = dummy ( i, j )
        end do
       end do
        else
        RLWTOA=SPVAL
        endif      
! NCAR WRF does not output accumulated fluxes so set the bitmap of these fluxes to 0
      do j = jsta_2l, jend_2u
        do i = 1, im
!	   RLWTOA(I,J)=SPVAL
	   RSWINC(I,J)=SPVAL
           ASWIN(I,J)=SPVAL  
	   ASWOUT(I,J)=SPVAL
	   ALWIN(I,J)=SPVAL
	   ALWOUT(I,J)=SPVAL
	   ALWTOA(I,J)=SPVAL
	   ASWTOA(I,J)=SPVAL
	   ARDLW=1.0
	   ARDSW=1.0
	   NRDLW=1
	   NRDSW=1
        end do
       end do

      VarName='XLAT'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)   &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            GDLAT ( i, j ) = dummy ( i, j )
! compute F = 2*omg*sin(xlat)
            f(i,j) = 1.454441e-4*sin(gdlat(i,j)*DTR)
        end do
       end do
! pos north
      print*,'read past GDLAT'
!      print*,'GDLAT at ',ii,jj,' = ',GDLAT(ii,jj)
	else
	GDLAT=SPVAL
	F=SPVAL
	endif

!-------------------------------------------------------------

      VarName='XLONG'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)   &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            GDLON ( i, j ) = dummy ( i, j )
        end do
       end do
!       print*,'GDLON at ',ii,jj,' = ',GDLON(ii,jj)
	else
	GDLON=SPVAL
	endif

!------------------------------------

      VarName='SWDOWN'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)   &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, iret2)

      VarName='ALBEDO'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)     &
     & ,DUMMY2,(im*jm),mpi_real4, mpi_status_ignore, ierr)

	if(iret .eq. 0 .and. iret2 .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            ALBEDO ( i, j ) = dummy2 ( i, j )
            RSWIN ( i, j ) = dummy ( i, j )
            RSWOUT ( i, j ) = RSWIN ( i, j ) * ALBEDO ( i, j )
        end do
       end do
	else
	ALBEDO=SPVAL
	RSWIN=SPVAL
	RSWOUT=SPVAL
	endif
!----------------------------------

      VarName='TMN'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)  &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

	if (IRET .EQ. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            TG ( i, j ) = dummy ( i, j )
            SOILTB ( i, j ) = dummy ( i, j )
        end do
       end do
	else
	TG=SPVAL
	SOILTB=SPVAL
	endif

!------------------------------------------------

!
! XLAND 1 land 2 sea
      VarName='XLAND'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
	do i = 1, im
            SM ( i, j ) = dummy ( i, j ) - 1.0
        end do
       end do
	else
	write(0,*) 'no XLAND...setting to SPVAL'
	SM=SPVAL
	endif
       
!----------------------------------------

      VarName='UST'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            USTAR ( i, j ) = dummy ( i, j ) 
        end do
       end do 
	else
	USTAR=SPVAL
	endif

!----------------------------------------

      VarName='AKHS'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            AKHS ( i, j ) = dummy ( i, j )
        end do
       end do
        else
        AKHS=SPVAL
        endif

!----------------------------------------

      VarName='AKMS'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)

        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            AKMS ( i, j ) = dummy ( i, j )
        end do
       end do
        else
        AKMS=SPVAL
        endif

!--------------------------------------------------
      IF(MODELNAME /= 'RAPR')THEN      
       VarName='PBLH'
       print*,'start reading PBLH'
       call retrieve_index(index,VarName,varname_all,nrecs,iret)
       call mpi_file_read_at(iunit,file_offset(index+1) &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)
       if (iret .eq. 0) then
        do j = jsta_2l, jend_2u
         do i = 1, im
            PBLH ( i, j ) = dummy ( i, j ) 
         end do
        end do
       else
	PBLH=SPVAL
       endif
      ELSE
! PBL depth from GSD
       do j = jsta_2l, jend_2u
        do i = 1, im
!   Is there any mixed layer at all?
          if (thv(i,j,lm-1) .lt. thv(i,j,lm)) then
            ZSF=ZINT(I,J,NINT(LMH(I,J))+1)
!   Calculate k1 level as first above PBL top
            do 34 k=3,LM
              k1 = k
! - give theta-v at the sfc a 0.5K boost in
!         the PBL height definition
              if (thv(i,j,lm-k+1).gt.thv(i,j,lm)  &
                     +0.5) go to 341
 34         continue
 341        continue
           zpbltop = zmid(i,j,lm-k1+1) +  &
                   (thv(i,j,lm)+0.5-thv(i,j,lm-k1+1))  &
                 * (zmid(i,j,lm-k1+2)-zmid(i,j,lm-k1+1))  &
                 / (thv(i,j,lm-k1+2) - thv(i,j,lm-k1+1))

            PBLH ( i, j ) = zpbltop - zsf
          else
            PBLH ( i, j ) = 0.
          endif
        end do
       end do
       ENDIF

       deallocate(thv)
       print*,'done reading or deriving PBLH'
!-------------------------------------------------------------------

!
      VarName='HFX'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)    &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)
	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            TWBS( i, j ) = dummy ( i, j )
        end do
       end do
	else
	TWBS=SPVAL
	endif

!-------------------------------------------------------------------
      IF(iSF_SURFACE_PHYSICS.NE.3) then       
       VarName='LH'
      ELSE
       VarName='QFX'
      END IF   

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)   &
     & ,DUMMY,(im*jm),mpi_real4, mpi_status_ignore, ierr)
	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            QWBS(I,J) = dummy ( i, j )
        end do
       end do
	else
	QWBS=SPVAL
	endif

! NCAR WRF does not output accumulated fluxes so bitmask out these fields
      do j = jsta_2l, jend_2u
        do i = 1, im
           SFCSHX(I,J)=SPVAL  
	   SFCLHX(I,J)=SPVAL
	   SUBSHX(I,J)=SPVAL
	   SNOPCX(I,J)=SPVAL
	   SFCUVX(I,J)=SPVAL
	   POTEVP(I,J)=SPVAL
	   NCFRCV(I,J)=SPVAL
	   NCFRST(I,J)=SPVAL
	   ASRFC=1.0
	   NSRFC=1
        end do
       end do

      VarName='SNOWC'

      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)                  &
     & ,DUMMY,(im*jm),mpi_real4,mpi_status_ignore, ierr)
	if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            PCTSNO( i, j ) = dummy ( i, j )
        end do
       end do
	else
	PCTSNO=SPVAL
	endif

! SRD
! get 2-d variables

      VarName='WSPD10MAX'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)                  &
     & ,DUMMY,(im*jm),mpi_real4,mpi_status_ignore, ierr)
        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
            WSPD10MAX( i, j ) = dummy ( i, j )
        end do
       end do
        else
            WSPD10MAX=spval
        end if

      VarName='W_UP_MAX'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)                  &
     & ,DUMMY,(im*jm),mpi_real4,mpi_status_ignore, ierr)
        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
           W_UP_MAX( i, j ) = dummy ( i, j )
        end do
       end do
        else
          W_UP_MAX=spval
        end if

      VarName='W_DN_MAX'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)                  &
     & ,DUMMY,(im*jm),mpi_real4,mpi_status_ignore, ierr)
        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
           W_DN_MAX( i, j ) = dummy ( i, j )
        end do
       end do
        else
          W_DN_MAX=spval
        end if

      VarName='W_MEAN'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)                  &
     & ,DUMMY,(im*jm),mpi_real4,mpi_status_ignore, ierr)
        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
           W_MEAN( i, j ) = dummy ( i, j )
        end do
       end do
        else
          W_MEAN=spval
        end if

      VarName='REFD_MAX'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)                  &
     & ,DUMMY,(im*jm),mpi_real4,mpi_status_ignore, ierr)
        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
           REFD_MAX( i, j ) = dummy ( i, j )
        end do
       end do
        else
          REFD_MAX=spval
        end if

      VarName='UP_HELI_MAX'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)                  &
     & ,DUMMY,(im*jm),mpi_real4,mpi_status_ignore, ierr)
        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
           UP_HELI_MAX( i, j ) = dummy ( i, j )
        end do
       end do
        else
          UP_HELI_MAX=spval
        end if

      VarName='GRPL_MAX'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      call mpi_file_read_at(iunit,file_offset(index+1)                  &
     & ,DUMMY,(im*jm),mpi_real4,mpi_status_ignore, ierr)
        if (iret .eq. 0) then
       do j = jsta_2l, jend_2u
        do i = 1, im
           GRPL_MAX( i, j ) = dummy ( i, j )
        end do
       end do
        else
          GRPL_MAX=spval
        end if

! pos east
       call collect_loc(gdlat,dummy)
       if(me.eq.0)then
        latstart=nint(dummy(1,1)*1000.)
        latlast=nint(dummy(im,jm)*1000.)
       end if
       write(6,*) 'laststart,latlast B calling bcast= ', latstart,latlast
       call mpi_bcast(latstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(latlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*) 'laststart,latlast A calling bcast= ', latstart,latlast
       call collect_loc(gdlon,dummy)
       if(me.eq.0)then
        lonstart=nint(dummy(1,1)*1000.)
        lonlast=nint(dummy(im,jm)*1000.)
       end if
       write(6,*)'lonstart,lonlast B calling bcast= ', lonstart,lonlast
       call mpi_bcast(lonstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(lonlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*)'lonstart,lonlast A calling bcast= ', lonstart,lonlast
!

!!
!! 
!!
        write(6,*) 'filename in INITPOST=', filename

!	status=nf_open(filename,NF_NOWRITE,ncid)
!	        write(6,*) 'returned ncid= ', ncid
!        status=nf_get_att_real(ncid,varid,'DX',tmp)
!	dxval=int(tmp)
!        status=nf_get_att_real(ncid,varid,'DY',tmp)
!	dyval=int(tmp)
!        status=nf_get_att_real(ncid,varid,'CEN_LAT',tmp)
!	cenlat=int(1000.*tmp)
!        status=nf_get_att_real(ncid,varid,'CEN_LON',tmp)
!	cenlon=int(1000.*tmp)
!        status=nf_get_att_real(ncid,varid,'TRUELAT1',tmp)
!	truelat1=int(1000.*tmp)
!        status=nf_get_att_real(ncid,varid,'TRUELAT2',tmp)
!	truelat2=int(1000.*tmp)
!        status=nf_get_att_real(ncid,varid,'MAP_PROJ',tmp)
!        maptype=int(tmp)
!	status=nf_close(ncid)

!	dxval=30000.
! 	dyval=30000.
!
!        write(6,*) 'dxval= ', dxval
!        write(6,*) 'dyval= ', dyval
!        write(6,*) 'cenlat= ', cenlat
!        write(6,*) 'cenlon= ', cenlon
!        write(6,*) 'truelat1= ', truelat1
!        write(6,*) 'truelat2= ', truelat2
!        write(6,*) 'maptype is ', maptype
!

!MEB not sure how to get these 
       do j = jsta_2l, jend_2u
        do i = 1, im
          if(msft(i,j)>small .and. msft(i,j)/=spval)then
            DX ( i, j ) = dxval/MSFT(I,J)
            DY ( i, j ) = dyval/MSFT(I,J)
          endif
        end do
       end do

! Convert DXVAL and DYVAL for ARW rotated
! latlon from meters to radian
        if(maptype==6)then
         dxval=(DXVAL * 360.)/(ERAD*2.*pi)*1000.
         dyval=(DYVAL * 360.)/(ERAD*2.*pi)*1000.
         print*,'dx and dy for arw rotated latlon= ', &
         dxval,dyval
        end if
!MEB not sure how to get these 

! close up shop
      call mpi_file_close(iunit, ierr)

!tgs Define smoothing flag for isobaric output
              IF(MODELNAME == 'RAPR')THEN
                SMFLAG=.TRUE.
              ELSE
                SMFLAG=.FALSE.
              ENDIF

! generate look up table for lifted parcel calculations

      THL=210.
      PLQ=70000.

      CALL TABLE(PTBL,TTBL,PT,RDQ,RDTH,RDP,RDTHE,PL,THL,QS0,SQS,STHE,THE0)

      CALL TABLEQ(TTBLQ,RDPQ,RDTHEQ,PLQ,THL,STHEQ,THE0Q)


!     
!     
      IF(ME.EQ.0)THEN
        WRITE(6,*)'  SPL (POSTED PRESSURE LEVELS) BELOW: '
        WRITE(6,51) (SPL(L),L=1,LSM)
   50   FORMAT(14(F4.1,1X))
   51   FORMAT(8(F8.1,1X))
      ENDIF
!     
!     COMPUTE DERIVED TIME STEPPING CONSTANTS.
!
!MEB need to get DT
!      DT = 120. !MEB need to get DT
      NPHS = 1  !CHUANG SET IT TO 1 BECAUSE ALL THE INST PRECIP ARE ACCUMULATED 1 TIME STEP
      DTQ2 = DT * NPHS  !MEB need to get physics DT
      TSPH = 3600./DT   !MEB need to get DT
!MEB need to get DT

! Randomly specify accumulation period because WRF EM does not
! output accumulation fluxes yet and accumulated fluxes are bit
! masked out

      TSRFC=1.0
      TRDLW=1.0
      TRDSW=1.0
      THEAT=1.0
      TCLOD=1.0
      TPREC=float(ifhr)  ! WRF EM does not empty precip buket at all
      print*,'TSRFC TRDLW TRDSW= ',TSRFC, TRDLW, TRDSW
!how am i going to get this information?
!      NPREC  = INT(TPREC *TSPH+D50)
!      NHEAT  = INT(THEAT *TSPH+D50)
!      NCLOD  = INT(TCLOD *TSPH+D50)
!      NRDSW  = INT(TRDSW *TSPH+D50)
!      NRDLW  = INT(TRDLW *TSPH+D50)
!      NSRFC  = INT(TSRFC *TSPH+D50)
!how am i going to get this information?
!     
!     IF(ME.EQ.0)THEN
!       WRITE(6,*)' '
!       WRITE(6,*)'DERIVED TIME STEPPING CONSTANTS'
!       WRITE(6,*)' NPREC,NHEAT,NSRFC :  ',NPREC,NHEAT,NSRFC
!       WRITE(6,*)' NCLOD,NRDSW,NRDLW :  ',NCLOD,NRDSW,NRDLW
!     ENDIF
!
!     COMPUTE DERIVED MAP OUTPUT CONSTANTS.
      DO L = 1,LSM
         ALSL(L) = ALOG(SPL(L))
      END DO
!
!HC WRITE IGDS OUT FOR WEIGHTMAKER TO READ IN AS KGDSIN
        if(me.eq.0)then
        print*,'writing out igds'
        igdout=110
!        open(igdout,file='griddef.out',form='unformatted'
!     +  ,status='unknown')
        if(maptype == 1)THEN  ! Lambert conformal
          print*,'IGDS= ',im,jm,LATSTART,LONSTART,8  &
          ,STANDLON,DXVAL,DYVAL,0,64,TRUELAT2,TRUELAT1
          WRITE(igdout)3
          WRITE(6,*)'igd(1)=',3
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
!          WRITE(igdout)CENLON
          WRITE(igdout)STANDLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)TRUELAT2
          WRITE(igdout)TRUELAT1
          WRITE(igdout)255
        ELSE IF(MAPTYPE .EQ. 2)THEN  !Polar stereographic
          WRITE(igdout)5
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
          WRITE(igdout)CENLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)TRUELAT2  !Assume projection at +-90
          WRITE(igdout)TRUELAT1
          WRITE(igdout)255
        ELSE IF(MAPTYPE .EQ. 3)THEN  !Mercator
          WRITE(igdout)1
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
          WRITE(igdout)latlast
          WRITE(igdout)lonlast
          WRITE(igdout)TRUELAT1
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)255
        END IF
        end if

	print*, 'end of INITPOST_BIN_MPIIO'
!     
!
      RETURN
      END SUBROUTINE INITPOST_BIN_MPIIO
