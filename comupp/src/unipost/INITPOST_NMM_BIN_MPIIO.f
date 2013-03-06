      SUBROUTINE INITPOST_NMM_BIN_MPIIO
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
!   05-12-05  H CHUANG - ADD CAPABILITY TO OUTPUT OFF-HOUR FORECAST WHICH HAS
!               NO INPACTS ON ON-HOUR FORECAST
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
      use wrf_io_flags_mod
      use params_mod
      use lookup_mod
      use ctlblk_mod
      use gridspec_mod
      use module_io_int_idx, only: io_int_index, io_int_loc, r_info
      use initpost_nmm_bin_mpiio_read, only: fetch_data

!
!     INCLUDE/SET PARAMETERS.
!     
      INCLUDE "mpif.h"

      character(len=31) :: VarName
      integer :: Status
      character startdate*19,SysDepInfo*80,cgar*1
      character startdate2(19)*4
! 
!     NOTE: SOME INTEGER VARIABLES ARE READ INTO DUMMY ( A REAL ). THIS IS OK
!     AS LONG AS REALS AND INTEGERS ARE THE SAME SIZE.
!
!     ALSO, EXTRACT IS CALLED WITH DUMMY ( A REAL ) EVEN WHEN THE NUMBERS ARE
!     INTEGERS - THIS IS OK AS LONG AS INTEGERS AND REALS ARE THE SAME SIZE.
      LOGICAL RUNB,SINGLRST,SUBPOST,NEST,HYDRO
      LOGICAL IOOMG,IOALL
      CHARACTER*32 LABEL
      CHARACTER*40 CONTRL,FILALL,FILMST,FILTMP,FILTKE,FILUNV                  &
         , FILCLD,FILRAD,FILSFC
      CHARACTER*4 RESTHR
      CHARACTER FNAME*80,ENVAR*50,BLANK*4
      INTEGER IDATB(3),IDATE(8),JDATE(8)
!
!     DECLARE VARIABLES.
!
      INTEGER, ALLOCATABLE, DIMENSION(:,:)    :: IBUF
      REAL,    ALLOCATABLE, DIMENSION(:)      :: SLDPTH2, RINC, ETA1, ETA2
      REAL,    ALLOCATABLE, DIMENSION(:,:)    :: BUF, DUMMY
      REAL,    ALLOCATABLE, DIMENSION(:,:,:)  :: FI, &
                                                 BUFSOIL, BUF3D, BUF3D2, BUF3DX
!jw
      integer ii,jj,js,je,jev,iyear,imn,iday,itmp,ioutcount,istatus,   &
              nsrfc,nrdlw,nrdsw,nheat,nclod,                           &
              iunit,nrecs,I,J,L

      character*80        :: titlestring

      type(r_info), pointer         :: r(:) => NULL()
      integer(kind=mpi_offset_kind) :: pos
      integer                       :: n

!
      DATA BLANK/'    '/
!
!***********************************************************************
!     START INIT HERE.
!
      WRITE(6,*)'INITPOST:  ENTER INITPOST'

      ! Allocate the local arrays
      allocate(SLDPTH2(NSOIL), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate SLDPTH2'
       stop
      end if
      allocate(RINC(5), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate RINC'
       stop
      end if
      allocate(ETA1(LM+1), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate ETA1'
       stop
      end if
      allocate(ETA2(LM+1), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate ETA2'
       stop
      end if
      allocate(DUMMY(IM, JM), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate DUMMY'
       stop
      end if
      allocate(FI(IM,JM,2), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate FI'
       stop
      end if
      allocate(ibuf(im,jsta_2l:jend_2u), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate ibuf'
       stop
      end if
      allocate(buf(im,jsta_2l:jend_2u), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate buf'
       stop
      end if
      allocate(bufsoil(im,nsoil,jsta_2l:jend_2u), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate bufsoil'
       stop
      end if
      allocate(buf3d(im,jm,lm), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate buf3d'
       stop
      end if
      allocate(buf3d2(im,jm,lp1), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate buf3d2'
       stop
      end if
      allocate(buf3dx(im,lm,jm), stat=ierr)
      if (ierr /= 0) then
       write(6,*)'Error unable to allocate buf3dx'
       stop
      end if
!     
!     
!     STEP 1.  READ MODEL OUTPUT FILE
!
!***
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

!  The end j row is going to be jend_2u for all variables except for V.
      JS=JSTA_2L
      JE=JEND_2U
      IF (JEND_2U.EQ.JM) THEN
       JEV=JEND_2U+1
      ELSE
       JEV=JEND_2U
      ENDIF

      ! Get an index of the file
      call io_int_index(filename, r, ierr)
      if (ierr /= 0) then
       print*,'Error obtinaing index of: ', trim(filename)
       stop
      end if

      ! MPI Open the file
      call mpi_file_open(mpi_comm_world, filename,                      &
                         mpi_mode_rdonly, mpi_info_null, iunit, ierr)
      if (ierr /= 0) then
       print*,"Error opening file with mpi io"
       stop
      end if

      call fetch_data(iunit, r, 'TITLE', dst=titlestring, ierr=ierr)

!  OK, since all of the variables are dimensioned/allocated to be
!  the same size, this means we have to be careful int getVariable
!  to not try to get too much data.  For example, 
!  DUM3D is dimensioned IM+1,JM+1,LM+1 but there might actually
!  only be im,jm,lm points of data available for a particular variable.  
! get metadata

      call fetch_data(iunit, r, 'MP_PHYSICS', dst=imp_physics, ierr=ierr)
      if (ierr /= 0) then
         imp_physics=5        ! assume ferrier if nothing specified
      endif

! Initializes constants for Ferrier microphysics
      if(imp_physics==5 .or. imp_physics==85 .or. imp_physics==95)then
       CALL MICROINIT(imp_physics)
      end if

      call fetch_data(iunit, r,'CU_PHYSICS', dst=icu_physics, ierr=ierr)
      if (ierr /= 0) then
         icu_physics=4        ! assume SAS if nothing specified
      endif
      if(icu_physics==84) icu_physics=4  ! HWRF SAS = SAS
      print*,'CU_PHYSICS= ',icu_physics

      call fetch_data(iunit,r,'SF_SURFACE_PHYSICS',dst=isf_physics,ierr=ierr)
      print*,'SF_PHYSICS= ',isf_physics

      call fetch_data(iunit, r, 'START_DATE', dst=startdate, ierr=ierr)
      if (ierr /= 0) then
        print*,"Error reading START_DATE using MPIIO"
      else
        print*,'START_DATE from MPIIO READ= ', startdate
      end if

      jdate=0
      idate=0
      read(startdate,15)iyear,imn,iday,ihrst,imin       
 15   format(i4,1x,i2,1x,i2,1x,i2,1x,i2)
      print*,'start yr mo day hr min =',iyear,imn,iday,ihrst,imin
      print*,'processing yr mo day hr min='                             &
         ,idat(3),idat(1),idat(2),idat(4),idat(5)
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

! Getting tstart
      tstart=0.
      call fetch_data(iunit, r, 'TSTART', dst=tstart, ierr=ierr)
      print*,'tstart= ',tstart

      IF(tstart .GT. 1.0E-2)THEN
       ifhr=ifhr+NINT(tstart)
       rinc=0
       idate=0
       rinc(2)=-1.0*ifhr
       call w3movdat(rinc,jdate,idate)
       SDAT(1)=idate(2)
       SDAT(2)=idate(3)
       SDAT(3)=idate(1)
       IHRST=idate(5)       
       print*,'new forecast hours for restrt run= ',ifhr
       print*,'new start yr mo day hr min =',sdat(3),sdat(1)               &
             ,sdat(2),ihrst,imin
      END IF

      RESTRT=.TRUE.  ! set RESTRT as default

      call fetch_data(iunit, r, 'DX', dst=garb, ierr=ierr)
      print*,'DX from MPIIO READ= ',garb
      dxval=nint(garb*1000.) ! E-grid dlamda in degree
      write(6,*) 'dxval= ', dxval

      call fetch_data(iunit, r, 'DY', dst=garb, ierr=ierr)
      print*,'DY from MPIIO READ= ',garb
      dyval=nint(garb*1000.) ! E-grid dlamda in degree
      write(6,*) 'dyval= ', dyval

      call fetch_data(iunit, r, 'DT', dst=dt, ierr=ierr)
      write(6,*) 'DT= ', DT

      call fetch_data(iunit, r, 'CEN_LAT', dst=garb, ierr=ierr)
      print*,'CEN_LAT from MPIIO READ= ',garb
      cenlat=nint(garb*1000.)
      write(6,*) 'cenlat= ', cenlat

      call fetch_data(iunit, r, 'CEN_LON', dst=garb, ierr=ierr)
      print*,'CEN_LON from MPIIO READ= ',garb
      cenlon=nint(garb*1000.)
      write(6,*) 'cenlon= ', cenlon

      ! Does HWRF (NMM) use TRUELAT1 and TRUELAT2?
      call fetch_data(iunit, r, 'TRUELAT1', dst=garb, ierr=ierr)
      print*,'TRUELAT1 from MPIIO READ= ',garb
      TRUELAT1=nint(garb*1000.)
      write(6,*) 'truelat1= ', TRUELAT1

      call fetch_data(iunit, r, 'TRUELAT2', dst=garb, ierr=ierr)
      print*,'TRUELAT2 from MPIIO READ= ',garb
      TRUELAT2=nint(garb*1000.)
      write(6,*) 'truelat2= ', TRUELAT2

      call fetch_data(iunit, r, 'MAP_PROJ', dst=maptype, ierr=ierr)
      write(6,*) 'maptype is ', maptype

      call fetch_data(iunit, r, 'GRIDTYPE', dst=gridtype, ierr=ierr)
      write(6,*) 'gridtype is ', gridtype

      VarName='HBM2'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        HBM2=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1) 
        call fetch_data(iunit, r, VarName, pos, n, HBM2, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          HBM2=SPVAL
        end if
      end if

      VarName='SM'      
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SM=SPVAL
      else        
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, sm, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SM=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
             if (j.eq.jm/2 .and. mod(i,10).eq.0)                    &   
                print*,'sample SM= ',i,j,sm(i,j)
           enddo
          enddo 
        end if 
      end if

      VarName='SICE'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SICE=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, SICE, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SICE=SPVAL
        end if
      end if
      

      VarName='PD'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        PD=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, PD, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          PD=SPVAL
        end if
      end if

      VarName='FIS'
      call io_int_loc(VarName, r, pos, n, iret)
      
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        FIS=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, FIS, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          FIS=SPVAL
        end if
      end if

      VarName='T'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        T=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          T=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
!            T ( i, j, l ) = buf3d ( i, ll, j )
             T ( i, j, l ) = buf3d ( i, j, ll )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample T= ',    &
               i,j,l,T ( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if	 
	  

      VarName='Q'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        Q=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          Q=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             Q ( i, j, l ) = buf3d ( i, j, ll )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample Q= ',    &
               i,j,l,Q ( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if
      
      print*,'finish reading mixing ratio'
      ii=im/2
      jj=(jsta+jend)/2
!      print*,'Q at ',ii,jj,ll,' = ',Q(ii,jj,ll)

      VarName='U'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        U=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          U=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             U ( i, j, l ) = buf3d ( i, j, ll )
	     UH( i, j, l ) = U( i, j, l )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample U= ',    &
               i,j,l,U ( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if

      VarName='V'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        V=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          V=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             V ( i, j, l ) = buf3d ( i, j, ll )
	     VH( i, j, l ) = V( i, j, l )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample V= ',   &
               i,j,l,V ( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after V'
      
      varname='DX_NMM'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        DX=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, dx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          DX=SPVAL
        end if
      end if

      varname='ETA1'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ETA1=SPVAL
      else
        call fetch_data(iunit, r, VarName, pos, n, ETA1, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ETA1=SPVAL
        end if
      end if

      varname='ETA2'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ETA2=SPVAL
      else
        call fetch_data(iunit, r, VarName, pos, n, ETA2, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ETA2=SPVAL
        end if
      end if
      
      open(75,file='ETAPROFILE.txt',form='formatted',                    &
              status='unknown')
        DO L=1,lm+1 
	 IF(L .EQ. 1)THEN
	  write(75,1020)L, 0., 0.
	 ELSE 
	  write(75,1020)L, ETA1(lm+2-l), ETA2(lm+2-l)
	 END IF     
!         print*,'L, ETA1, ETA2= ',L, ETA1(l), ETA2(l)
        END DO
 1020   format(I3,2E17.10)	
	close (75)

      varname='PDTOP'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        PDTOP=SPVAL
      else
        ! n = 1
        call fetch_data(iunit, r, VarName, pos, n, pdtop, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          PDTOP=SPVAL
        end if
      end if

        varname='PT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        PT=SPVAL
      else
        ! n = 1
        call fetch_data(iunit, r, VarName, pos, n, pt, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          PT=SPVAL
        end if
      end if
      
      print*,'PT, PDTOP= ',PT,PDTOP
	
      varname='PBLH'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        PBLH=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, pblh, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          PBLH=SPVAL
        end if
      end if

     varname='MIXHT' !PLee (3/07)
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        MIXHT=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
        n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, mixht, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          MIXHT=SPVAL
        end if
      end if

      varname='USTAR'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        USTAR=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, ustar, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          USTAR=SPVAL
        end if
      end if

      varname='Z0'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        Z0=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, z0, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          Z0=SPVAL
        end if
      end if
      
      varname='THS'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        THS=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, ths, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          THS=SPVAL
        end if
      end if
	
      VarName='QS'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        QS=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, qs, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          QS=SPVAL
        end if
      end if

      varname='TWBS'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        TWBS=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, twbs, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          TWBS=SPVAL
        end if
      end if

      varname='QWBS'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        QWBS=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, qwbs, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          QWBS=SPVAL
        end if
      end if

      varname='PREC' ! instantaneous precip rate?
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        PREC=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, prec, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          PREC=SPVAL
        end if
      end if

      varname='ACPREC' ! accum total precip
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ACPREC=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, ACPREC, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ACPREC=SPVAL
        end if
      end if
      
      varname='CUPREC' ! accum cumulus precip
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CUPREC=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, cuprec, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CUPREC=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
	     ANCPRC(I,J)=ACPREC(I,J)-CUPREC(I,J)
           enddo
          enddo
        end if
      end if
      write(0,*)' after CUPREC'

      varname='LSPA'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        LSPA=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, lspa, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          LSPA=SPVAL
        end if
      end if

      varname='SNO'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SNO=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, sno, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SNO=SPVAL
        end if
      end if
     
      varname='SI'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SI=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, si, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SI=SPVAL
        end if
      end if
      
      varname='CLDEFI'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CLDEFI=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, cldefi, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CLDEFI=SPVAL
        end if
      end if

      varname='TH10'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        TH10=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, th10, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          TH10=SPVAL
        end if
      end if	
       
      varname='Q10'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        Q10=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, q10, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          Q10=SPVAL
        end if
      end if

      varname='PSHLTR'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        PSHLTR=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, pshltr, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          PSHLTR=SPVAL
        end if
      end if

      varname='TSHLTR'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        TSHLTR=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, tshltr, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          TSHLTR=SPVAL
        end if
      end if

      varname='QSHLTR'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        QSHLTR=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, qshltr, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          QSHLTR=SPVAL
	end if  
      end if
      write(0,*)' after QSHLTR'
      
      VarName='Q2'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        Q2=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          Q2=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             Q2( i, j, l ) = buf3d ( i, j, ll )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample Q2= ',   &
               i,j,l,Q2( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after Q2'

      varname='AKHS_OUT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        AKHS=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, akhs, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          AKHS=SPVAL
        end if
      end if	

      varname='AKMS_OUT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        AKMS=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, akms, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          AKMS=SPVAL
        end if
      end if		
	
      varname='ALBASE'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ALBASE=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, albase, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ALBASE=SPVAL
        end if
      end if	
	
      varname='ALBEDO'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ALBEDO=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, albedo, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ALBEDO=SPVAL
        end if
      end if	

      varname='CZEN'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CZEN=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, czen, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CZEN=SPVAL
        end if
      end if

      varname='CZMEAN'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CZMEAN=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, czmean, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CZMEAN=SPVAL
        end if
      end if	
       print*,'max CZMEAN= ',maxval(czmean) 

      varname='GLAT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        GDLAT=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, buf, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          GDLAT=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
             F(I,J)=1.454441e-4*sin(buf(I,J))   ! 2*omeg*sin(phi)
             GDLAT(I,J)=buf(I,J)*RTD
	     
           enddo
          enddo
        end if
      end if
      
      varname='GLON'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        GDLON=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, buf, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          GDLON=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
             GDLON(I,J)=buf(I,J)*RTD
	     if(i.eq.409.and.j.eq.835)print*,'GDLAT GDLON in INITPOST='  &
      	     ,i,j,GDLAT(I,J),GDLON(I,J)
           enddo
          enddo
        end if
      end if
      
       if(jsta.le.594.and.jend.ge.594)print*,'gdlon(120,594)= ',       &
       gdlon(120,594)

      varname='MXSNAL'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        MXSNAL=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, mxsnal, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          MXSNAL=SPVAL
        endif
      end if	
      write(0,*)' after MXSNAL'
	
      varname='RADOT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        RADOT=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, radot, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          RADOT=SPVAL
        end if
      end if
      
      varname='SIGT4'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SIGT4=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, sigt4, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SIGT4=SPVAL
        end if
      end if
       
      varname='TGROUND'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        TG=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, tg, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          TG=SPVAL
        end if
      end if

      varname='CWM'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CWM=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CWM=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             CWM( i, j, l ) = buf3d ( i, j, ll )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample CWM= ',   &
               i,j,l,CWM( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if     
      write(0,*)' after CWM'

      varname='F_ICE'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        F_ice=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3dx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          F_ice=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             F_ice( i, j, l ) = buf3dx ( i, ll, j )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample F_ice= ', &
               i,j,l,F_ice( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if	
      write(0,*)' after F_ICE'

      varname='F_RAIN'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        F_rain=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3dx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          F_rain=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             F_rain( i, j, l ) = buf3dx ( i, ll, j )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample F_rain= ',&
               i,j,l,F_rain( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after F_RAIN'

      varname='F_RIMEF'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        F_RimeF=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3dx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          F_RimeF=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             F_RimeF( i, j, l ) = buf3dx ( i, ll, j )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,                &
               'sample F_RimeF= ',i,j,l,F_RimeF( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after F_RimeF'

       varname='CLDFRA'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CFR=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CFR=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             CFR( i, j, l ) = buf3d ( i, j, ll )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample CFR= ', &
               i,j,l,CFR( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after CLDFRA'

      varname='SR'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SR=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, sr, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SR=SPVAL
        end if
      end if	

      varname='CFRACH'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CFRACH=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, cfrach, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CFRACH=SPVAL
        end if
      end if

      varname='CFRACL'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CFRACL=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, cfracl, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CFRACL=SPVAL
        end if
      end if

      varname='CFRACM'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CFRACM=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, cfracm, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CFRACM=SPVAL
        end if
      end if
      write(6,*) 'maxval CFRACM: ', maxval(CFRACM)

      varname='ISLOPE'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ISLOPE=NINT(SPVAL)
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, islope, ierr)
         if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ISLOPE=NINT(SPVAL)
        end if
      end if	
	
!	varname='SOILTB'
!	write(6,*) 'call getVariableB for : ', VarName
!      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,
!     &  IM,1,JM,1,IM,JS,JE,1)

! either assign SLDPTH to be the same as eta (which is original
! setup in WRF LSM) or extract thickness of soil layers from wrf
! output

! assign SLDPTH to be the same as eta
! jkw comment out because Pleim Xiu only has 2 layers
! jkw         SLDPTH(1)=0.10
! jkw         SLDPTH(2)=0.3
! jkw         SLDPTH(3)=0.6
! jkw         SLDPTH(4)=1.0
! Initialize soil depth to some bogus value
! to alert user if not found in wrfout file
       do I=1,NSOIL
        SLDPTH(I) = 0.0
       end do

      if (isf_PHYSICS == 3) then
! get SLDPTH from wrf output
        VarName='SLDPTH'
      call io_int_loc(VarName, r, pos, n, iret)
        if (iret /= 0) then
          print*,VarName," not found in file-Assigned missing values"
          SLDPTH2=SPVAL
        else
          n = NSOIL
          call fetch_data(iunit, r, VarName, pos, n, SLDPTH2, ierr)
          if (ierr /= 0) then
            print*,"Error reading ", VarName,"Assigned missing values"
            SLDPTH2=SPVAL
          end if
        end if

        DUMCST=0.0
        DO N=1,NSOIL
          DUMCST=DUMCST+SLDPTH2(N)
        END DO
        IF(ABS(DUMCST-0.).GT.1.0E-2)THEN
          DO N=1,NSOIL
            SLLEVEL(N)=SLDPTH2(N)
          END DO
        END IF
        print*,'SLLEVEL ',(SLLEVEL(N),N=1,NSOIL)

      else ! isf_PHYSICS /= 3
        VarName='DZSOIL'
      call io_int_loc(VarName, r, pos, n, iret)
        if (iret /= 0) then
          print*,VarName," not found in file-Assigned missing values"
          SLDPTH2=SPVAL
        else
          n = NSOIL
          call fetch_data(iunit, r, VarName, pos, n, SLDPTH2, ierr)
          if (ierr /= 0) then
            print*,"Error reading ", VarName,"Assigned missing values"
            SLDPTH2=SPVAL
          end if
        end if ! if (iret /= 0)

        DUMCST=0.0
        DO N=1,NSOIL
          DUMCST=DUMCST+SLDPTH2(N)
        END DO
        IF(ABS(DUMCST-0.).GT.1.0E-2)THEN
          DO N=1,NSOIL
            SLDPTH(N)=SLDPTH2(N)
          END DO
        END IF
        print*,'SLDPTH= ',(SLDPTH(N),N=1,NSOIL)
      end if   ! if (isf_PHYSICS==3)

      VarName='CMC'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CMC=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, cmc, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CMC=SPVAL
        end if
      end if
      
      varname='GRNFLX'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        GRNFLX=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, grnflx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          GRNFLX=SPVAL
        end if
      end if

      varname='PCTSNO'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        PCTSNO=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, pctsno, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          PCTSNO=SPVAL
        end if
      end if	
	
      varname='SOILTB'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SOILTB=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, soiltb, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SOILTB=SPVAL
        end if
      end if

      varname='VEGFRC'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        VEGFRC=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, vegfrc, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          VEGFRC=SPVAL
        end if
      end if

      VarName='SH2O'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SH2O=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im*nsoil
	n=im*(jend_2u-jsta_2l+1)*nsoil
        call fetch_data(iunit, r, VarName, pos, n, bufsoil, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SH2O=SPVAL
        else
	  do l = 1, nsoil
           do j = jsta_2l, jend_2u
            do i = 1, im
             SH2O(I,J,L)=bufSOIL(I,L,J)
	    enddo 
           enddo
          enddo
        end if
      end if
      write(0,*)' after SH2O'

      VarName='SMC'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SMC=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im*nsoil
	n=im*(jend_2u-jsta_2l+1)*nsoil
        call fetch_data(iunit, r, VarName, pos, n, bufsoil, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SMC=SPVAL
        else
	  do l = 1, nsoil
           do j = jsta_2l, jend_2u
            do i = 1, im
             SMC(I,J,L)=bufSOIL(I,L,J)
	    enddo 
           enddo
          enddo
        end if
      end if

      print*,'SMC at ',ii,jj,N,' = ',smc(ii,jj,1),smc(ii,jj,2)         &
      ,smc(ii,jj,3),smc(ii,jj,4)

      VarName='STC'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        STC=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im*nsoil
	n=im*(jend_2u-jsta_2l+1)*nsoil
        call fetch_data(iunit, r, VarName, pos, n, bufsoil, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          STC=SPVAL
        else
	  do l = 1, nsoil
           do j = jsta_2l, jend_2u
            do i = 1, im
             STC(I,J,L)=bufSOIL(I,L,J)
	    enddo 
           enddo
          enddo
        end if
      end if
    
      if(jj.ge.jsta.and.jj.le.jend)                                    &
        print*,'STC at ',ii,jj,' = ',stc(ii,jj,1),stc(ii,jj,2)         &
      ,stc(ii,jj,3),stc(ii,jj,4)
      write(0,*)' after STC'

      VarName='PINT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        PINT=SPVAL
      else
	n=im*jm*lp1
        call fetch_data(iunit, r, VarName, pos, n, buf3d2, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          PINT=SPVAL
        else
	  do l = 1, lp1
	   ll=lp1-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             PINT( i, j, l ) = buf3d2 ( i, j, ll )
!      if(i==1.and.j==250)then
!        write(0,*)' l=',l,' iin PINT=',pint(i,j,l)
!      endif
             ALPINT(I,J,L)=ALOG(PINT(I,J,L))     
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'PINT= ',       &
               i,j,l,PINT ( i, j, l )
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after PINT'

      do l = 2, lm+1
       do j = jsta_2l, jend_2u
        do i = 1, im
!         PMID ( i, j, l-1 ) = EXP((ALOG(PINT(I,J,L-1))+
!     &               ALOG(PINT(I,J,L)))*0.5)
         PMID ( i, j, l-1 ) = (PINT(I,J,L-1)+                              &
                     PINT(I,J,L))*0.5 ! representative of what model does
      if(i==1.and.j==250.and.l==2)then
        write(0,*)' pmid=',pmid(i,j,l-1)                                   &
      ,           ' pint=',pint(i,j,l-1),pint(i,j,l)
      endif
         if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'PMID= ',              &
               i,j,l-1,PMID ( i, j, l-1 )
        end do
       end do
      end do 
      write(0,*)' after PMID'

      do l = 1, lm
       do j = jsta, jend
        do i = 1, im-MOD(J,2) 
	 IF(J .EQ. 1 .AND. I .LT. IM)THEN   !SOUTHERN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J,L)+PMID(I+1,J,L))
         ELSE IF(J.EQ.JM .AND. I.LT.IM)THEN   !NORTHERN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J,L)+PMID(I+1,J,L))
         ELSE IF(I .EQ. 1 .AND. MOD(J,2) .EQ. 0) THEN   !WESTERN EVEN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J-1,L)+PMID(I,J+1,L))
	 ELSE IF(I .EQ. IM .AND. MOD(J,2) .EQ. 0                             &  
      	 .AND. J .LT. JM) THEN   !EASTERN EVEN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J-1,L)+PMID(I,J+1,L))  
         ELSE IF (MOD(J,2) .LT. 1) THEN
           PMIDV(I,J,L)=0.25*(PMID(I,J,L)+PMID(I-1,J,L)                       &
             +PMID(I,J+1,L)+PMID(I,J-1,L))
         ELSE
           PMIDV(I,J,L)=0.25*(PMID(I,J,L)+PMID(I+1,J,L)                       &
             +PMID(I,J+1,L)+PMID(I,J-1,L))
         END IF  
        end do
       end do
      end do
      write(0,*)' after PMIDV'


!!!!! COMPUTE Z
       do j = jsta_2l, jend_2u
        do i = 1, im
            ZINT(I,J,LM+1)=FIS(I,J)/G
	if (I .eq. 1 .and. J .eq. jsta_2l) then
                   write(6,*) 'G,ZINT: ', G,ZINT(I,J,LM+1)
	endif
            FI(I,J,1)=FIS(I,J)
        end do
       end do
      write(0,*)' after FI'

! SECOND, INTEGRATE HEIGHT HYDROSTATICLY
      DO L=LM,1,-1
       do j = jsta_2l, jend_2u
        do i = 1, im
         FI(I,J,2)=HTM(I,J,L)*T(I,J,L)*(Q(I,J,L)*D608+1.0)*RD*                &
                   (ALPINT(I,J,L+1)-ALPINT(I,J,L))+FI(I,J,1)
         ZINT(I,J,L)=FI(I,J,2)/G
      if(i==1.and.j==250.and.l<=2)then
        write(0,*)' zint=',zint(i,j,l),' fi=',fi(i,j,2)                       &
      ,           fi(i,j,1)
        write(0,*)' t=',t(i,j,l),' q=',q(i,j,l)
        write(0,*)' alpint=',alpint(i,j,l+1),alpint(i,j,l)
      endif
         if(i.eq.ii.and.j.eq.jj)                                              &
        print*,'L,sample HTM,T,Q,ALPINT(L+1),ALPINT(l),ZINT= '                &
        ,l,HTM(I,J,L),T(I,J,L),Q(I,J,L),ALPINT(I,J,L+1),                      &
        ALPINT(I,J,L),ZINT(I,J,L)
         FI(I,J,1)=FI(I,J,2)
        ENDDO
       ENDDO
      END DO
      print*,'finish deriving geopotential in nmm'
      write(0,*)' after ZINT lm=',lm,' js=',js,' je=',je,' im=',im
      write(0,*)' zmid lbounds=',lbound(zmid),' ubounds=',ubound(zmid)
      write(0,*)' zint lbounds=',lbound(zint),' ubounds=',ubound(zint)
      write(0,*)' pmid lbounds=',lbound(pmid),' ubounds=',ubound(pmid)
      write(0,*)' pint lbounds=',lbound(pint),' ubounds=',ubound(pint)
!
      DO L=1,LM
!      write(0,*)' zmid l=',l
       DO J=JS,JE
!      write(0,*)' zmid j=',j
        DO I=1,IM
!      write(0,*)' zmid i=',i
!         ZMID(I,J,L)=(ZINT(I,J,L+1)+ZINT(I,J,L))*0.5  ! ave of z
!      write(0,*)' pmid=',pmid(i,j,l)
!      write(0,*)' pint=',pint(i,j,l),pint(i,j,l+1)
!      write(0,*)' zint=',zint(i,j,l),zint(i,j,l+1)
         FACT=(ALOG(PMID(I,J,L))-ALOG(PINT(I,J,L)))/                      &
               (ALOG(PINT(I,J,L+1))-ALOG(PINT(I,J,L)))	 
         ZMID(I,J,L)=ZINT(I,J,L)+(ZINT(I,J,L+1)-ZINT(I,J,L))*FACT
        ENDDO
       ENDDO
      ENDDO

      write(0,*)' before W'
      VarName='W'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        WH=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d2, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          WH=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             WH( i, j, l ) = buf3d2 ( i, j, ll )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'WH= ',               &
               i,j,l,WH ( i, j, l )
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after W'

      VarName='ACFRCV'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ACFRCV=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, acfrcv, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ACFRCV=SPVAL
        end if
      end if
      write(0,*)' after ACFRCV'
      
      write(6,*) 'MAX ACFRCV: ', maxval(ACFRCV)

      VarName='ACFRST'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ACFRST=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, acfrst, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ACFRST=SPVAL
        end if
      end if
      write(6,*) 'max ACFRST ', maxval(ACFRST)
      write(0,*)' after ACFRST'

!insert-mp
      VarName='SSROFF'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SSROFF=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, ssroff, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SSROFF=SPVAL
        end if
      end if
      write(0,*)' after SSROFF'

! reading UNDERGROUND RUNOFF
      VarName='BGROFF'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        BGROFF=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, bgroff, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          BGROFF=SPVAL
        end if
      end if
      write(0,*)' after BGROFF'
      
      VarName='RLWIN'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        RLWIN=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, rlwin, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          RLWIN=SPVAL
        end if
      end if
      write(0,*)' after RLWIN'

      VarName='RLWTOA'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        RLWTOA=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, rlwtoa, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          RLWTOA=SPVAL
        end if
      end if
      write(0,*)' after RLWTOA'

      VarName='ALWIN'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ALWIN=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, alwin, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ALWIN=SPVAL
        end if
      end if
      write(0,*)' after ALWIN'
      
      VarName='ALWOUT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ALWOUT=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, alwout, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ALWOUT=SPVAL
        end if
      end if

      VarName='ALWTOA'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ALWTOA=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, alwtoa, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ALWTOA=SPVAL
        end if
      end if

      VarName='RSWIN'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        RSWIN=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, rswin, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          RSWIN=SPVAL
        end if
      end if
      
      VarName='RSWINC'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        RSWINC=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, rswinc, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          RSWINC=SPVAL
        end if
      end if
      
       print*,'max RSWINC= ',maxval(RSWINC)

      VarName='RSWOUT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        RSWOUT=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, rswout, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          RSWOUT=SPVAL
        end if
      end if

      VarName='ASWIN'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ASWIN=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, aswin, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ASWIN=SPVAL
        end if
      end if
      
      VarName='ASWOUT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ASWOUT=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, aswout, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ASWOUT=SPVAL
        end if
      end if

      VarName='ASWTOA'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ASWTOA=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, aswtoa, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ASWTOA=SPVAL
        end if
      end if
      write(0,*)' after ASWTOA'

      VarName='SFCSHX'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SFCSHX=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, sfcshx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SFCSHX=SPVAL
        end if
      end if
      
      VarName='SFCLHX'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SFCLHX=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, sfclhx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SFCLHX=SPVAL
        end if
      end if
      
      VarName='SUBSHX'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SUBSHX=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, subshx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SUBSHX=SPVAL
        end if
      end if

      VarName='SNOPCX'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SNOPCX=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, snopcx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SNOPCX=SPVAL
        end if
      end if
      write(0,*)' after SNOPCX'
	
      VarName='SFCUVX'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SFCUVX=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, sfcuvx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SFCUVX=SPVAL
        end if
      end if

      VarName='POTEVP'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        POTEVP=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, potevp, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          POTEVP=SPVAL
        end if
      end if
      write(0,*)' after POTEVP'

      varname='RLWTT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        RLWTT=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          RLWTT=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             RLWTT( i, j, l ) = buf3d ( i, j, ll )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample RLWTT= ', &
                  i,j,l,RLWTT( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after RLWTT'

      varname='RSWTT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        RSWTT=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          RSWTT=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             RSWTT( i, j, l ) = buf3d ( i, j, ll )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample RSWTT= ', &
                  i,j,l,RSWTT( i, j, l )
             ttnd ( i, j, l ) = rswtt(i,j,l) + rlwtt(i,j,l)	     
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after RSWTT'

      varname='TCUCN'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        TCUCN=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          TCUCN=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             TCUCN( i, j, l ) = buf3d ( i, j, ll )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample TCUCN= ', &
                  i,j,l,TCUCN( i, j, l )	 
             if(l.eq.lm.and.ABS(TCUCN( i, j, l )).gt.1.0e-4)              &
          print*,'nonzero TCUCN',i,j,l,TCUCN( i, j, l )    
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after TCUCN'
	
      varname='TRAIN'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        TRAIN=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3d, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          TRAIN=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             TRAIN( i, j, l ) = buf3d ( i, j, ll )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample TRAIN= ',  &
                  i,j,l,TRAIN( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after TRAIN'

      VarName='NCFRCV'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        NCFRCV=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, ibuf, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          NCFRCV=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
             NCFRCV(I,J)=FLOAT(ibuf(I,J))
           enddo
          enddo
        end if
      end if
      write(0,*)' after NCFRCV'

      VarName='NCFRST'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        NCFRST=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, ibuf, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          NCFRST=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
             NCFRST(I,J)=FLOAT(IBUF(I,J))
           enddo
          enddo
        end if
      end if
      
! set default to not empty buket
      NSRFC=0
      NRDLW=0
      NRDSW=0
      NHEAT=0
      NCLOD=0
      NPREC=0

      VarName='NPHS0'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        NPHS=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, NPHS, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          NPHS=NINT(SPVAL)
	end if  
      end if
      write(6,*) 'NPHS= ', NPHS

      VarName='NPREC'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        NPREC=NINT(SPVAL)
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, NPREC, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          NPREC=NINT(SPVAL)
	end if  
      end if
      write(6,*) 'NPREC= ', NPREC

      VarName='NCLOD'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        NCLOD=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, NCLOD, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          NCLOD=SPVAL
        end if
      end if
      write(6,*) 'NCLOD= ', NCLOD
      
      VarName='NHEAT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        NHEAT=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, NHEAT, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          NHEAT=SPVAL
        end if
      end if
      write(6,*) 'NHEAT= ', NHEAT      

      VarName='NRDLW'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        NRDLW=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, NRDLW, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          NRDLW=SPVAL
        end if
      end if
      write(6,*) 'NRDLW= ', NRDLW

      VarName='NRDSW'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        NRDSW=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, NRDSW, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          NRDSW=SPVAL
        end if
      end if
	write(6,*) 'NRDSW= ', NRDSW

      VarName='NSRFC'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        NSRFC=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, NSRFC, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          NSRFC=SPVAL
        end if
      end if
	write(6,*) 'NSRFC= ', NSRFC

      VarName='AVRAIN'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        AVRAIN=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, AVRAIN, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          AVRAIN=SPVAL
        end if
      end if
      write(6,*) 'AVRAIN= ', AVRAIN
      write(0,*)' after AVRAIN'

      VarName='AVCNVC'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        AVCNVC=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, AVCNVC, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          AVCNVC=SPVAL
        end if
      end if
      write(6,*) 'AVCNVC= ', AVCNVC

      VarName='ARDLW'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ARDLW=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, ARDLW, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ARDLW=SPVAL
        end if
      end if
      write(6,*) 'ARDLW= ', ARDLW

      VarName='ARDSW'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ARDSW=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, ARDSW, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ARDSW=SPVAL
        end if
      end if
	write(6,*) 'ARDSW= ', ARDSW
	
      VarName='ASRFC'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ASRFC=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, ASRFC, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ASRFC=SPVAL
        end if
      end if
	write(6,*) 'ASRFC= ', ASRFC	

      VarName='APHTIM'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        APHTIM=SPVAL
      else
        n = 1
        call fetch_data(iunit, r, VarName, pos, n, APHTIM, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          APHTIM=SPVAL
        end if
      end if

! reading 10 m wind
      VarName='U10'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        U10=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, u10, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          U10=SPVAL
        end if
      end if
      if(jj.ge.jsta.and.jj.le.jend)                                     &
               print*,'U10 at ',ii,jj,' = ',U10(ii,jj)
      write(0,*)' after U10'

      VarName='V10'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        V10=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, v10, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          V10=SPVAL
        end if
      end if
      if(jj.ge.jsta.and.jj.le.jend)                                     &
           print*,'V10 at ',ii,jj,' = ',V10(ii,jj)
!
!
! reading SMSTAV
      VarName='SMSTAV'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SMSTAV=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, smstav, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SMSTAV=SPVAL
        end if
      end if

      VarName='SMSTOT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SMSTOT=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, smstot, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SMSTOT=SPVAL
        end if
      end if
! reading VEGETATION TYPE 
      VarName='IVGTYP'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        IVGTYP=NINT(SPVAL)
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, ivgtyp, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          IVGTYP=NINT(SPVAL)
        end if
      end if	

      VarName='ISLTYP' 
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ISLTYP=NINT(SPVAL)
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, isltyp, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ISLTYP=NINT(SPVAL)
        end if
      end if

      VarName='SFCEVP'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SFCEVP=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, sfcevp, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SFCEVP=SPVAL
        end if
      end if

      VarName='SFCEXC'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SFCEXC=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, sfcexc, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SFCEXC=SPVAL
        end if
      end if

      VarName='ACSNOW'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ACSNOW=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, acsnow, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ACSNOW=SPVAL
        end if
      end if
       
      VarName='ACSNOM'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        ACSNOM=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, acsnom, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          ACSNOM=SPVAL
        end if
      end if

      VarName='SST'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        SST=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, sst, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          SST=SPVAL
        end if
      end if
      if(jj.ge.jsta.and.jj.le.jend)                                      &
            print*,'SST at ',ii,jj,' = ',sst(ii,jj)      
      write(0,*)' after SST'

! ADDED TAUX AND TAUY in POST --------------- zhan's doing
      VarName='TAUX'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        MDLTAUX=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
        n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, mdltaux, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          MDLTAUX=SPVAL
        end if
      end if
      if(jj.ge.jsta.and.jj.le.jend)                &
        print*,'MDLTAUX at ',ii,jj,' = ',mdltaux(ii,jj)
      write(0,*)' after MDLTAUX'

      VarName='TAUY'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        MDLTAUY=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
        n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, mdltauy, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          MDLTAUY=SPVAL
        end if
      end if
      if(jj.ge.jsta.and.jj.le.jend)                 &
        print*,'MDLTAUY at ',ii,jj,' = ',mdltauy(ii,jj)
      write(0,*)' after MDLTAUY'
! zhang's dong ends

      VarName='EL_PBL'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        EL_PBL=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3dx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          EL_PBL=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             EL_PBL( i, j, l ) = buf3dx ( i, j ,ll)
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample EL= ', &
                  i,j,l,EL_PBL( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after EL_PBL'

      VarName='EXCH_H'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        EXCH_H=SPVAL
      else
	n=im*jm*lm
        call fetch_data(iunit, r, VarName, pos, n, buf3dx, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          EXCH_H=SPVAL
        else
	  do l = 1, lm
	   ll=lm-l+1
           do j = jsta_2l, jend_2u
            do i = 1, im
             EXCH_H( i, j, l ) = buf3dx ( i, j, ll )
	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample EXCH= ', &
                  i,j,l,EXCH_H( i, j, l )	     
            end do
           end do
          end do 
	end if 
      end if
      write(0,*)' after EXCH_H'

      VarName='THZ0'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        THZ0=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, thz0, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          THZ0=SPVAL
        end if
      end if
      print*,'THZ0 at ',ii,jj,' = ',THZ0(ii,jj)

      VarName='QZ0'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        QZ0=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, qz0, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          QZ0=SPVAL
        end if
      end if
      print*,'QZ0 at ',ii,jj,' = ',QZ0(ii,jj)

      VarName='UZ0'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        UZ0=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, uz0, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          UZ0=SPVAL
	end if  
      end if
      if(jj.ge.jsta.and.jj.le.jend)                                     &
           print*,'UZ0 at ',ii,jj,' = ',UZ0(ii,jj)

      VarName='VZ0'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        VZ0=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, vz0, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          VZ0=SPVAL
        end if
      end if
      if(jj.ge.jsta.and.jj.le.jend) print*,'VZ0 at ',ii,jj,'=',VZ0(ii,jj)


!
! Very confusing story ...
!
! Retrieve htop and hbot => They are named CNVTOP, CNVBOT in the model and
! with HBOTS,HTOPS (shallow conv) and HBOTD,HTOPD (deep conv) represent
! the 3 sets of convective cloud base/top arrays tied to the frequency
! that history files are written.
!
! IN THE *MODEL*, arrays HBOT,HTOP are similar to CNVTOP,CNVBOT but are
! used in radiation and are tied to the frequency of radiation updates.
!
! For historical reasons model arrays CNVTOP,CNVBOT are renamed HBOT,HTOP
! and manipulated throughout the post. 

! retrieve htop and hbot
!      VarName='HTOP'
      VarName='CNVTOP'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        HTOP=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, buf, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          HTOP=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
             HTOP ( i, j ) = float(LM)-buf(i,j)+1.0
             HTOP ( i, j ) = max(1.0,min(HTOP(I,J),float(LM)))
           enddo
          enddo
        end if
      end if
       print*,'maxval HTOP: ', maxval(HTOP)
      write(0,*)' after HTOP'

!      VarName='HBOT'
      VarName='CNVBOT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        HBOT=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, buf, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          HBOT=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
             HBOT ( i, j ) = float(LM)-buf(i,j)+1.0
             HBOT ( i, j ) = max(1.0,min(HBOT(I,J),float(LM)))
           enddo
          enddo
        end if
      end if      
       print*,'maxval HBOT: ', maxval(HBOT)
      write(0,*)' after HBOT'

      VarName='HTOPD'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        HTOPD=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, buf, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          HTOPD=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
             HTOPD ( i, j ) = float(LM)-buf(i,j)+1.0
           enddo
          enddo
        end if
      end if
       print*,'maxval HTOPD: ', maxval(HTOPD)

      VarName='HBOTD'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        HBOTD=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, buf, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          HBOTD=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
             HBOTD( i, j ) = float(LM)-buf(i,j)+1.0
           enddo
          enddo
        end if
      end if
       print*,'maxval HBOTD: ', maxval(HBOTD)

      VarName='HTOPS'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        HTOPS=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, buf, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          HTOPS=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
             HTOPS( i, j ) = float(LM)-buf(i,j)+1.0
           enddo
          enddo
        end if
      end if
       print*,'maxval HTOPS: ', maxval(HTOPS)
                                                                                 
      VarName='HBOTS'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        HBOTS=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, buf, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          HBOTS=SPVAL
        else
          do j = jsta_2l, jend_2u
           do i = 1, im
             HBOTS ( i, j ) = float(LM)-buf(i,j)+1.0
           enddo
          enddo
        end if
      end if
       print*,'maxval HBOTS: ', maxval(HBOTS)
      write(0,*)' after HBOTS'

      VarName='CUPPT'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CUPPT=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, cuppt, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CUPPT=SPVAL
        end if
      end if
       print*,'maxval CUPPT: ', maxval(CUPPT)

      VarName='CPRATE'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        CPRATE=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, cprate, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          CPRATE=SPVAL
        end if
      end if
       print*,'maxval CPRATE: ', maxval(CPRATE)

      VarName='HBM2'
      call io_int_loc(VarName, r, pos, n, iret)
      if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        HBM2=SPVAL
      else
        pos=pos+(jsta_2l-1)*4*im
	n=im*(jend_2u-jsta_2l+1)
        call fetch_data(iunit, r, VarName, pos, n, hbm2, ierr)
        if (ierr /= 0) then
          print*,"Error reading ", VarName,"Assigned missing values"
          HBM2=SPVAL
        end if
      end if

!!!! DONE GETTING

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            IF(ABS(T(I,J,L)).GT.1.0E-3)                                &
              OMGA(I,J,L) = -WH(I,J,L)*PMID(I,J,L)*G/                   &
                       (RD*T(I,J,L)*(1.+D608*Q(I,J,L)))

        end do
       end do
      end do
      write(0,*)' after OMGA'

! pos east
      call collect_loc(gdlat,dummy)
      if(me.eq.0)then
        latstart=nint(dummy(1,1)*1000.)
        latlast=nint(dummy(im,jm)*1000.)
! temporary patch for nmm wrf for moving nest. gopal's doing
! jkw changed if statement as per MP's suggestion
! jkw        if(mod(im,2).ne.0) then
! chuang: test
        icen=(im+1)/2
        jcen=(jm+1)/2
        if(mod(im,2).ne.0)then !per Pyle, jm is always odd
         if(mod(jm+1,4).ne.0)then
          cenlat=nint(dummy(icen,jcen)*1000.)
         else
          cenlat=nint(0.5*(dummy(icen-1,jcen)+dummy(icen,jcen))*1000.)
         end if
        else
         if(mod(jm+1,4).ne.0)then
          cenlat=nint(0.5*(dummy(icen,jcen)+dummy(icen+1,jcen))*1000.)
         else
          cenlat=nint(dummy(icen,jcen)*1000.)
         end if
        end if

!        if(mod(im,2).eq.0) then
!           icen=(im+1)/2
!           jcen=(jm+1)/2
!           cenlat=nint(dummy(icen,jcen)*1000.)
!        else
!           icen=im/2
!           jcen=(jm+1)/2
!           cenlat=nint(0.5*(dummy(icen,jcen)+dummy(icen+1,jcen))*1000.)
!        end if
        
      end if
      write(6,*) 'laststart,latlast,cenlat B calling bcast= ', &
     &    latstart,latlast,cenlat
      call mpi_bcast(latstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
      call mpi_bcast(latlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
      call mpi_bcast(cenlat,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
      write(6,*) 'laststart,latlast,cenlat A calling bcast= ', &
     &    latstart,latlast,cenlat

      call collect_loc(gdlon,dummy)
      if(me.eq.0)then
        lonstart=nint(dummy(1,1)*1000.)
        lonlast=nint(dummy(im,jm)*1000.)
! temporary patch for nmm wrf for moving nest. gopal's doing
!lrb changed if statement as per MP's suggestion
!lrb        if(mod(im,2).ne.0) then
!Chuang: test
        icen=(im+1)/2
        jcen=(jm+1)/2
        if(mod(im,2).ne.0)then !per Pyle, jm is always odd
         if(mod(jm+1,4).ne.0)then
          cenlon=nint(dummy(icen,jcen)*1000.)
         else
          cenlon=nint(0.5*(dummy(icen-1,jcen)+dummy(icen,jcen))*1000.)
         end if
        else
         if(mod(jm+1,4).ne.0)then
          cenlon=nint(0.5*(dummy(icen,jcen)+dummy(icen+1,jcen))*1000.)
         else
          cenlon=nint(dummy(icen,jcen)*1000.)
         end if
        end if

!        if(mod(im,2).eq.0) then
!           icen=(im+1)/2
!           jcen=(jm+1)/2
!           cenlon=nint(dummy(icen,jcen)*1000.)
!        else
!           icen=im/2
!           jcen=(jm+1)/2
!           cenlon=nint(0.5*(dummy(icen,jcen)+dummy(icen+1,jcen))*1000.)
!        end if
       end if

       write(6,*)'lonstart,lonlast,cenlon B calling bcast= ', &
     &      lonstart,lonlast,cenlon
       call mpi_bcast(lonstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(lonlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(cenlon,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*)'lonstart,lonlast,cenlon A calling bcast= ', &
     &      lonstart,lonlast,cenlon
!
        write(6,*) 'filename in INITPOST=', filename


!MEB not sure how to get these 
       do j = jsta_2l, jend_2u
        do i = 1, im
!            DX ( i, j ) = dxval  !MEB ???
!            DY ( i, j ) = dyval*DTR*ERAD  

!!!!!!!!!!!!!!!!!!!!! DY ????

            DY ( i, j ) =   0.001*ERAD*DYVAL*DTR  ! like A*DPH
        end do
       end do
!MEB not sure how to get these 
      write(0,*)' after DX DY'

! close up shop
!      call ext_int_ioclose ( DataHandle, Status )

! generate look up table for lifted parcel calculations

      THL=210.
      PLQ=70000.

      CALL TABLE(PTBL,TTBL,PT,                                       &
                RDQ,RDTH,RDP,RDTHE,PL,THL,QS0,SQS,STHE,THE0)

      CALL TABLEQ(TTBLQ,RDPQ,RDTHEQ,PLQ,THL,STHEQ,THE0Q)
      write(0,*)' after TABLEQ'


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
!      NPHS = 4  !MEB need to get physics DT
      DTQ2 = DT * NPHS  !MEB need to get physics DT
      TSPH = 3600./DT   !MEB need to get DT

      TSRFC=float(NSRFC)/TSPH
      IF(NSRFC.EQ.0)TSRFC=float(ifhr)  !in case buket does not get emptied
      TRDLW=float(NRDLW)/TSPH
      IF(NRDLW.EQ.0)TRDLW=float(ifhr)  !in case buket does not get emptied
      TRDSW=float(NRDSW)/TSPH
      IF(NRDSW.EQ.0)TRDSW=float(ifhr)  !in case buket does not get emptied
      THEAT=float(NHEAT)/TSPH
      IF(NHEAT.EQ.0)THEAT=float(ifhr)  !in case buket does not get emptied
      TCLOD=float(NCLOD)/TSPH
      IF(NCLOD.EQ.0)TCLOD=float(ifhr)  !in case buket does not get emptied
      TPREC=float(NPREC)/TSPH
      IF(NPREC.EQ.0)TPREC=float(ifhr)  !in case buket does not get emptied
!       TPREC=float(ifhr)
      print*,'TSRFC TRDLW TRDSW= ',TSRFC, TRDLW, TRDSW
!MEB need to get DT

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
      write(0,*)' after ALSL'
!
!HC WRITE IGDS OUT FOR WEIGHTMAKER TO READ IN AS KGDSIN
        if(me.eq.0)then
        print*,'writing out igds'
        igdout=110
!        open(igdout,file='griddef.out',form='unformatted'
!     +  ,status='unknown')
        if(maptype .eq. 1)THEN  ! Lambert conformal
          WRITE(igdout)3
          WRITE(6,*)'igd(1)=',3
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
        ELSE IF(MAPTYPE.EQ.0 .OR. MAPTYPE.EQ.203)THEN  !A STAGGERED E-GRID
          WRITE(igdout)203
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)136
          WRITE(igdout)CENLAT
          WRITE(igdout)CENLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)64
          WRITE(igdout)0
          WRITE(igdout)0
          WRITE(igdout)0
        END IF
        end if
      write(0,*)' after writes'

      call mpi_file_close(iunit,ierr)

      ! Deallocate the local arrays
      deallocate(SLDPTH2)
      deallocate(RINC)
      deallocate(ETA1)
      deallocate(ETA2)
      deallocate(DUMMY)
      deallocate(FI)
      deallocate(ibuf)
      deallocate(buf)
      deallocate(bufsoil)
      deallocate(buf3d)
      deallocate(buf3d2)
      deallocate(buf3dx)

      RETURN
      END
