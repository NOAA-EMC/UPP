      SUBROUTINE INITPOST_GFS(NREC,iunit,gfile)
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
      use gfsio_module
      type(gfsio_gfile) :: gfile
!
!     INCLUDE/SET PARAMETERS.
!     
      include 'wrf_io_flags.h'
!      INCLUDE "parmeta"
      INCLUDE "params"
!      INCLUDE "parmsoil"
      INCLUDE "parm.tbl"
      INCLUDE "mpif.h"
      integer,parameter:: MAXPTS=1000000 ! max im*jm points
! This version of INITPOST shows how to initialize, open, read from, and
! close a NetCDF dataset. In order to change it to read an internal (binary)
! dataset, do a global replacement of _ncd_ with _int_. 

!     character(len=32) :: fileName
!     character(len=19) :: DateStr ! 2002-03-05_03:00:00
!      integer :: DataHandle
      character(len=20) :: VarName
      character(len=20) :: VcoordName
      integer :: Status
      character startdate*19,SysDepInfo*80,cgar*1
      character startdate2(19)*4
! 
!     NOTE: SOME INTEGER VARIABLES ARE READ INTO DUMMY ( A REAL ). THIS IS OK
!     AS LONG AS REALS AND INTEGERS ARE THE SAME SIZE.
!
!     ALSO, EXTRACT IS CALLED WITH DUMMY ( A REAL ) EVEN WHEN THE NUMBERS ARE
!     INTEGERS - THIS IS OK AS LONG AS INTEGERS AND REALS ARE THE SAME SIZE.
      LOGICAL RUN,RUNB,RESTRT,SINGLRST
     1,       SIGMA,SUBPOST,NEST,HYDRO
      LOGICAL IOOMG,IOALL
      CHARACTER*32 LABEL
      CHARACTER*40 CONTRL,FILALL,FILMST,FILTMP,FILTKE,FILUNV
     &, FILCLD,FILRAD,FILSFC
      CHARACTER*4 RESTHR
      CHARACTER FNAME*80,ENVAR*50,BLANK*4
      INTEGER IDATB(3),IDATE(8),JDATE(8)
      INTEGER JPDS(200),JGDS(200),KPDS(200),KGDS(200)
      LOGICAL*1 LB(IM_JM)
      INTEGER IRET
      REAL BUFF(IM_JM)
!     
!     INCLUDE COMMON BLOCKS.
!
      INCLUDE "LOOKUP.comm"
      INCLUDE "CTLBLK.comm"
!     INCLUDE "SOILDEPTH.comm"
      INCLUDE "GRIDSPEC.comm"
      INCLUDE "RQSTFLD.comm"
!
!     DECLARE VARIABLES.
!     
      REAL fhour
      REAL SLDPTH2(NSOIL)
      REAL RINC(5)
      REAL ETA1(LM), ETA2(LM)
      REAL DUM1D (LM+1)
      REAL DUMMY ( IM, JM )
      REAL DUMMY2 ( IM, JM )
      REAL FI(IM,JM,2)
      INTEGER IDUMMY ( IM, JM )
!      REAL DUM3D ( IM, LM, JM )
!      REAL DUM3D2 ( IM, LM+1, JM ),DUMSOIL ( IM, NSOIL, JM )
!mp
	INTEGER CENLAT,CENLON,TRUELAT1,TRUELAT2
     	INTEGER LATSTART,LONSTART,LATLAST,LONLAST
	INTEGER DXVAL, DYVAL

      character*8,allocatable:: recname(:)
      character*16,allocatable  :: reclevtyp(:)
      integer,allocatable:: reclev(:)
      integer ierr
   
      integer ibuf(im,jsta_2l:jend_2u)
      real buf(im,jsta_2l:jend_2u),bufsoil(im,nsoil,jsta_2l:jend_2u)
     +  ,buf3d(im,jsta_2l:jend_2u,lm),buf3d2(im,lp1,jsta_2l:jend_2u)
!
      DATA BLANK/'    '/
!
!***********************************************************************
!     START INIT HERE.
!
      WRITE(6,*)'INITPOST:  ENTER INITPOST'
!     
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
!  how do I get the filename? 
!      fileName = '/ptmp/wx20mb/wrfout_01_030500'
!      DateStr = '2002-03-05_18:00:00'
!  how do I get the filename?
!         call ext_int_ioinit(SysDepInfo,Status)
!          print*,'called ioinit', Status
!         call ext_int_open_for_read( trim(fileName), 0, 0, " ", 
!     &  DataHandle, Status)
!          print*,'called open for read', Status
!       if ( Status /= 0 ) then
!         print*,'error opening ',fileName, ' Status = ', Status ; stop
!       endif
! get date/time info
!  this routine will get the next time from the file, not using it
!      print *,'DateStr before calling ext_int_get_next_time=',DateStr
!      call ext_int_get_next_time(DataHandle, DateStr, Status)
!      print *,'DateStri,Status,DataHandle = ',DateStr,Status,DataHandle

!  The end j row is going to be jend_2u for all variables except for V.
      JS=JSTA_2L
      JE=JEND_2U
      IF (JEND_2U.EQ.JM) THEN
       JEV=JEND_2U+1
      ELSE
       JEV=JEND_2U
      ENDIF
! get start date
      if (me == 0)then
       print*,'nrec=',nrec
       allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
       call gfsio_getfilehead(gfile,iret=iret 
     1   ,idate=idate(1:4),fhour=fhour,recname=recname
     2   ,reclevtyp=reclevtyp,reclev=reclev)
       if(iret/=0)print*,'error getting idate,fhour'
       print *,'printing an inventory of GFS Grib file'
       print *,'recname=',(trim(recname(i)),i=1,nrec)
       print *,'reclevtyp=',(trim(reclevtyp(i)),i=1,nrec)
       print *,'reclev=',(reclev(i),i=1,nrec)
       deallocate(recname,reclevtyp, reclev)
! can't get idate and fhour, specify them for now
!       idate(4)=2006
!       idate(2)=9  
!       idate(3)=16
!       idate(1)=0
!       fhour=6.0
        print*,'idate before broadcast = ',(idate(i),i=1,4)
      end if
      call mpi_bcast(idate(1),4,MPI_INTEGER,0,mpi_comm_comp,iret)
      call mpi_bcast(fhour,1,MPI_REAL,0,mpi_comm_comp,iret)
      print*,'idate after broadcast = ',(idate(i),i=1,4)
      print*,'fhour = ',fhour
      iyear=idate(4)+2000 ! ask Jun if need to add 2000
      imn=idate(3) ! ask Jun 
      iday=idate(2) ! ask Jun
      ihrst=idate(1)
      jdate=0
      idate=0 
!      read(startdate,15)iyear,imn,iday,ihrst,imin       
 15   format(i4,1x,i2,1x,i2,1x,i2,1x,i2)
      print*,'start yr mo day hr min =',iyear,imn,iday,ihrst,imin
      print*,'processing yr mo day hr min='
     +,idat(3),idat(1),idat(2),idat(4),idat(5)
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
!      CALL W3DIFDAT(JDATE,IDATE,2,RINC)
!      ifhr=nint(rinc(2))
      CALL W3DIFDAT(JDATE,IDATE,0,RINC)
      ifhr=nint(rinc(2)+rinc(1)*24.)
      ifmin=nint(rinc(3))
!      if(ifhr /= nint(fhour))print*,'find wrong Grib file';stop
      print*,' in INITPOST ifhr ifmin fileName=',ifhr,ifmin,fileName
      
! Getting tstart
      tstart=0.
!      VarName='TSTART'
!      call retrieve_index(index,VarName,varname_all,nrecs,iret)
!      if (iret /= 0) then
!        print*,VarName," not found in file"
!      else
!        call mpi_file_read_at(iunit,file_offset(index)+5*4
!     + ,garb,1,mpi_real4
!     + , mpi_status_ignore, ierr)
!        if (ierr /= 0) then
!          print*,"Error reading ", VarName," using MPIIO"
!        else
!          print*,VarName, ' from MPIIO READ= ',garb
!	  tstart=garb
!        end if	
!      end if
      print*,'tstart= ',tstart
      
! Getiing restart
      
      RESTRT=.TRUE.  ! set RESTRT as default
!      call ext_int_get_dom_ti_integer(DataHandle,'RESTARTBIN',itmp
!     + ,1,ioutcount,istatus)
      
!      IF(itmp .LT. 1)THEN
!        RESTRT=.FALSE.
!      ELSE
!        RESTRT=.TRUE.
!      END IF
     
!      print*,'status for getting RESTARTBIN= ',istatus
     
!      print*,'Is this a restrt run? ',RESTRT
            
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
       print*,'new start yr mo day hr min =',sdat(3),sdat(1)
     +       ,sdat(2),ihrst,imin
      END IF 
      
      imp_physics=99 !set GFS mp physics to 99 for Zhao scheme
      print*,'MP_PHYSICS= ',imp_physics

! waiting to retrieve lat lon infor from raw GFS output

!      VarName='DX'


!      VarName='DY'

! GFS does not need DT to compute accumulated fields, set it to one
!      VarName='DT'
      DT=1
! GFS does not need truelat
!      VarName='TRUELAT1'

!      VarName='TRUELAT2'

! Specigy maptype=4 for Gaussian grid
      maptype=4
      write(6,*) 'maptype is ', maptype	  
! HBM2 is most likely not in Grib message, set them to ones
      HBM2=1.0
! sample print point
      ii=im/2
      jj=jm/2

! try to get kgds from flux grib file and then convert to igds that is used by GRIBIT.f

      if(me == 0)then       
       jpds=-1.0
       jgds=-1.0
                                                                                       
       call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,ierr)
       call R63W72(KPDS,KGDS,JPDS,IGDS(1:18))
!       print*,'IGDS for GFS= ',(IGDS(I),I=1,18)
      end if
      call mpi_bcast(igds(1),18,MPI_INTEGER,0,mpi_comm_comp,iret)      
      print*,'IGDS for GFS= ',(IGDS(I),I=1,18)

! start retrieving data using gfsio, first land/sea mask

!      VarName='land'
!      VcoordName='sfc'
!      l=1
      
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!       if (iret /= 0) then
!        print*,VarName," not found in file-Assigned missing values"
!        dummy=spval
!       else
!
!        do j = 1, jm
!           do i = 1, im
!             dummy(I,J)=1.0 - dummy(I,J) ! convert Grib message to 2D
!             if (j.eq.jm/2 .and. mod(i,10).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)
!     
!           enddo
!          enddo
!       end if
!      end if	
!
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,sm(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! start retrieving data using getgb, first land/sea mask
      Index=50
      VarName=avbl(index)
      jpds=-1
      jgds=-1
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        print*,'jpds(5),jpds(6)= ',jpds(5),jpds(6)
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,ierr)
       if (ierr /= 0) then
        print*,VarName," not found in file-Assigned missing values"
     +  ,"error code= ",ierr
        dummy=spval
       else
                                                                                       
        do j = 1, jm
           do i = 1, im
             dummy(I,J)=1.0 - dummy(I,J) ! convert to sea mask 
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)
                                                                                       
           enddo
          enddo
       end if
      end if
        
      print*,'scattering sea mask'                                                                               
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,sm(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,ierr)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',VarName, ' after scatter= ',i,j,sm(i,j)
        end do
       end do
!      itmp=icnt(1)
!      call mpi_scatter(dummy(1,1),itmp,mpi_real
!     + ,sm(1,jsta),itmp,mpi_real,0,MPI_COMM_COMP,ierr)
       print*,'error code from scattering sm= ',ierr
!        if (abs(ierr-0).gt.1)print*,'Error scattering array';stop
!        print*,'done scattering sea mask'
! sea ice mask using GFSIO
!      VarName='icec'
!      VcoordName='sfc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,sice(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! sea ice mask using getgb
      Index=51
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,ierr)
       if (ierr /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       end if
      end if
                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,sice(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,ierr)
       print*,'error code from scattering sice= ',ierr
!      if (ierr /= 0)print*,'Error scattering array';stop

! GFS does not output PD  
      pd=spval

       con_g       =9.80665e+0! gravity
! Terrain height * G   using gfsio 
      VarName='hgt'
      VcoordName='sfc'
      l=1
      if(me == 0)then
        print*,'retrieving sfc height from using gfsio'
        call gfsio_readrecv(gfile,trim(VarName)
     + ,trim(VcoordName),l,dummy,iret=iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
	else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)*con_G
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',trim(VarName), ' = ',i,j,dummy(i,j)     
           enddo
          enddo 
        end if
      end if	
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,fis(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! Terrain height * G   using getgb
!      Index=25
!      VarName=avbl(index)
!      jpds(5)=iq(index)
!      jpds(6)=is(index)
!      if(me == 0)then
!        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
!     +    ,k,kpds,kgds,lb,dummy,ierr)
!       if (ierr /= 0) then
!        print*,VarName," not found in file-Assigned missing values"
!        dummy=spval
!       else
                                                                                       
!        do j = 1, jm
!           do i = 1, im
!             dummy(I,J)=dummy(I,J)*G ! convert to gpm
!             if (j.eq.jm/2 .and. mod(i,10).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)
!                                                                                       
!           enddo
!          enddo
!       end if
!      end if
                                                                                       
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,fis(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,ierr)
      
! model level T
      print*,'start retrieving GFS T using gfsio'
      VarName='tmp'
      VcoordName='layer'
      do l=1,lm	
       if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
         print*,VarName," not found at level ",l,
     +	 " - Assigned missing values", "iret= ",iret
         dummy=spval
        end if
       end if	
       ll=lm-l+1
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,t(1,jsta,ll),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,ll,t(i,j,ll)
        end do
       end do
!       if (iret /= 0)print*,'Error scattering array';stop
      End do ! do loop for l

! model level q      
      VarName='spfh'
      VcoordName='layer'
      do l=1,lm	
       if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
	 print*,VarName," not found at level ",l,
     +	 " - Assigned missing values"
         dummy=spval
        end if
       end if	
       ll=lm-l+1
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,q(1,jsta,ll),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,ll,q(i,j,ll)
        end do
       end do
!       if (iret /= 0)print*,'Error scattering array';stop
      End do ! do loop for l
	
! model level u      
      VarName='ugrd'
      VcoordName='layer'
      do l=1,lm	
       if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
	 print*,VarName," not found at level ",l,
     +	 " - Assigned missing values"
         dummy=spval
        end if
       end if	
       ll=lm-l+1
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,uh(1,jsta,ll),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,ll,uh(i,j,ll)
        end do
       end do
!       if (iret /= 0)print*,'Error scattering array';stop
      End do ! do loop for l      
      
! model level v      
      VarName='vgrd'
      VcoordName='layer'
      do l=1,lm	
       if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
	 print*,VarName," not found at level ",l,
     +	 " - Assigned missing values"
         dummy=spval
        end if
       end if	
       ll=lm-l+1
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,vh(1,jsta,ll),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,ll,vh(i,j,ll)
        end do
       end do
!       if (iret /= 0)print*,'Error scattering array';stop
      End do ! do loop for l  
      
! model level pressure      
      VarName='pres'
      VcoordName='layer'
      do l=1,lm	
       if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
	 print*,VarName," not found at level ",l,
     +	 " - Assigned missing values"
         dummy=spval
        end if
       end if	
       ll=lm-l+1
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,pmid(1,jsta,ll),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,ll,pmid(i,j,ll)
        end do
       end do
!       if (iret /= 0)print*,'Error scattering array';stop
      End do ! do loop for l     
      
! GFS is on A grid and does not need PMIDV        

! Surface pressure  using gfsio 
      VarName='pres'
      VcoordName='sfc'
      l=1
      if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
      end if	
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,pint(1,jsta,lp1),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
      do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample PSFC after scatter= '
     +   ,i,j,pint(i,j,lp1)
        end do
       end do
!      if (iret /= 0)print*,'Error scattering array';stop

! Surface pressure using getgb
!      Index=24
!      VarName=avbl(index)
!      jpds(5)=iq(index)
!      jpds(6)=is(index)
!      if(me == 0)then
!        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
!     +    ,k,kpds,kgds,lb,dummy,ierr)
!       if (ierr /= 0) then
!        print*,VarName," not found in file-Assigned missing values"
!        dummy=spval
!       end if
!      end if                                                                                       
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,pint(1,jsta,lp1),icnt(me),mpi_real,0,MPI_COMM_COMP,ierr)
     
      do j=jsta,jend
	do i=1,im
	  ALPINT(I,J,LP1)=ALOG(PINT(I,J,LP1))     
!	    if(i==ii .and. j==jj)print*,'sample PSFC'
!     +        ,i,j,pint(i,j,lp1)	    
	end do
      end do       
      
! dp     
      VarName='dpres'
      VcoordName='layer'
      do l=1,lm	
       if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
	 print*,VarName," not found at level ",l,
     +	 " - Assigned missing values"
         dummy=spval
        end if
       end if	
       ll=lm-l+1
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,buf3d(1,jsta,ll),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!       if (iret /= 0)print*,'Error scattering array';stop
      End do ! do loop for l           
      
! construct interface pressure from psfc and dp from bottom up
!      do l=lm,1,-1
!        do j=jsta,jend
!	  do i=1,im
!	    pint(i,j,l)=pint(i,j,l+1)-buf3d(i,j,l)
!	    if(pint(i,j,l)< 0.)pint(i,j,l)=0. ! set model top to 0
!	    if(l > 1)ALPINT(I,J,L)=ALOG(PINT(I,J,L))     
!	    if(i==ii .and. j==jj)print*,'sample pint,pmid'
!     +        ,i,j,l,pint(i,j,l),pmid(i,j,l)	    
!	  end do
!	end do
!      end do	    

! construct interface pressure from model top (which is zero) and dp from top down
! PDTOP
      pdtop=spval
      pt=0.
      print*,'PT, PDTOP= ',PT,PDTOP
      do j=jsta,jend
        do i=1,im
	  pint(i,j,1)=PT
	end do
      end do	  
      do l=2,lm
        do j=jsta,jend
	  do i=1,im	    
	    pint(i,j,l)=pint(i,j,l-1)+buf3d(i,j,l-1)
	    ALPINT(I,J,L)=ALOG(PINT(I,J,L))     
	    if(i==ii .and. j==jj)print*,'sample pint,pmid'
     +        ,i,j,l,pint(i,j,l),pmid(i,j,l)	    
	  end do
	end do
      end do
      
!!!!! COMPUTE Z, GFS integrates Z on mid-layer instead
!!! use GFS contants to see if height becomes more aggreable to GFS pressure grib file
       con_rv      =4.6150e+2 ! gas constant H2O
       con_rd      =2.8705e+2 ! gas constant air
       con_fvirt   =con_rv/con_rd-1.

       do j = jsta, jend
        do i = 1, im
            ZINT(I,J,LM+1)=FIS(I,J)/con_G
	    FI(I,J,1)=FIS(I,J)+T(I,J,LM)
     +	    *(Q(I,J,LM)*con_fvirt+1.0)*con_rd ! using GFS consts	  
!     +	    *(Q(I,J,LM)*D608+1.0)*RD
     +      *(ALPINT(I,J,LM+1)-ALOG(PMID(I,J,LM)))
            ZMID(I,J,LM)=FI(I,J,1)/con_G
!            FI(I,J,1)=FIS(I,J)
        end do
       end do

! SECOND, INTEGRATE HEIGHT HYDROSTATICLY, GFS integrate height on mid-layer
      DO L=LM-1,1,-1
       do j = jsta, jend
        do i = 1, im
         FI(I,J,2)=0.5*(T(I,J,L)*(Q(I,J,L)*con_fvirt+1.0)
     1	           +T(I,J,L+1)*(Q(I,J,L+1)*con_fvirt+1.0))*con_rd*
!         FI(I,J,2)=0.5*(T(I,J,L)*(Q(I,J,L)*D608+1.0)     	 
!     1	           +T(I,J,L+1)*(Q(I,J,L+1)*D608+1.0))*RD*
     2             (ALOG(PMID(I,J,L+1))-ALOG(PMID(I,J,L)))
     3             +FI(I,J,1)
         ZMID(I,J,L)=FI(I,J,2)/con_G
	 
         if(i.eq.ii.and.j.eq.jj)
     1  print*,'L,sample T,Q,ALPMID(L+1),ALPMID(L),ZMID= '
     2  ,l,T(I,J,L),Q(I,J,L),ALOG(PMID(I,J,L+1)),
     3  ALOG(PMID(I,J,L)),ZMID(I,J,L)
     
         FI(I,J,1)=FI(I,J,2)
        ENDDO
       ENDDO
      END DO

      DO L=LM,2,-1  ! omit computing model top height because it's infinity
       DO J=JSTA,JEND
        DO I=1,IM
!         ZMID(I,J,L)=(ZINT(I,J,L+1)+ZINT(I,J,L))*0.5  ! ave of z
         FACT=(ALPINT(I,J,L)-ALOG(PMID(I,J,L)))/
     &         (ALOG(PMID(I,J,L-1))-ALOG(PMID(I,J,L)))          
         ZINT(I,J,L)=ZMID(I,J,L)+(ZMID(I,J,L-1)-ZMID(I,J,L))
     &                 *FACT 
         if(i.eq.ii.and.j.eq.jj)
     1     print*,'L ZINT= ',l,zint(i,j,l)	 	 
        ENDDO
       ENDDO
      ENDDO

! WRF way of integrating height on interfaces       
!      DO L=LM,1,-1
!       do j = jsta, jend
!        do i = 1, im
!         FI(I,J,2)=HTM(I,J,L)*T(I,J,L)*(Q(I,J,L)*D608+1.0)*RD*
!     1             (ALPINT(I,J,L+1)-ALPINT(I,J,L))+FI(I,J,1)
!         ZINT(I,J,L)=FI(I,J,2)/G
!         if(i.eq.ii.and.j.eq.jj)
!     1  print*,'L,sample HTM,T,Q,ALPINT(L+1),ALPINT(l),ZINT= '
!     2  ,l,HTM(I,J,L),T(I,J,L),Q(I,J,L),ALPINT(I,J,L+1),
!     3  ALPINT(I,J,L),ZINT(I,J,L)
!         FI(I,J,1)=FI(I,J,2)
!        ENDDO
!       ENDDO
!      END DO
!
!      DO L=1,LM
!       DO I=1,IM
!        DO J=JS,JE
!         FACT=(ALOG(PMID(I,J,L))-ALOG(PINT(I,J,L)))/
!     &         (ALOG(PINT(I,J,L+1))-ALOG(PINT(I,J,L)))	 
!         ZMID(I,J,L)=ZINT(I,J,L)+(ZINT(I,J,L+1)-ZINT(I,J,L))
!     &       *FACT
!        ENDDO
!       ENDDO
!      ENDDO

! GFS does not output W, will try to derive it
!      VarName='W'      

! ozone mixing ratio     
      VarName='o3mr'
      VcoordName='layer'
      do l=1,lm	
       if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
	 print*,VarName," not found at level ",l,
     +	 " - Assigned missing values"
         dummy=spval
        end if
       end if	
       ll=lm-l+1
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,o3(1,jsta,ll),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,ll,o3(i,j,ll)
        end do
       end do
!       if (iret /= 0)print*,'Error scattering array';stop
      End do ! do loop for l          

! cloud water and ice mixing ratio  for zhao scheme  
! need to look up old eta post to derive cloud water/ice from cwm 
! Zhao scheme does not produce suspended rain and snow
      qqr=0.
      qqs=0.
      qqi=0.
      VarName='clwmr'
      VcoordName='layer'
      do l=1,lm	
       if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
	 print*,VarName," not found at level ",l,
     +	 " - Assigned missing values"
         dummy=spval
        end if
       end if	
       ll=lm-l+1
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,cwm(1,jsta,ll),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
!	 if(t(i,j,ll) < (TFRZ-15.) )then ! dividing cloud water from ice
!	  qqi(i,j,ll)=cwm(i,j,ll)
!	 else 
	  qqw(i,j,ll)=cwm(i,j,ll)
!	 end if 
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,ll,cwm(i,j,ll)
        end do
       end do
!       if (iret /= 0)print*,'Error scattering array';stop
      End do ! do loop for l 
      
! waiting to retrieve lat lon from GFS to compute DX(i,j)      
!      varname='DX_NMM'      
! PBL height using gfsio
!      VarName='hpbl'
!      VcoordName='sfc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,pblh(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop 

! PBL height using getgb
      Index=221
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,ierr)
       if (ierr /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,pblh(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,ierr)

! frictional velocity using gfssio
!      VarName='uustar'
!      VcoordName='sfc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,ustar(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop           

! frictional velocity using getgb
      Index=45
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,ustar(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)

! roughness length using gfsio
!      VarName='sfcr'
!      VcoordName='sfc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,z0(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! roughness length using getgb
      Index=44
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,z0(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)

! surface potential T  using gfsio  
!      VarName='tmp'
!      VcoordName='sfc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!	else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)*(p1000/pint(i,j,lm+1))**CAPA ! convert to THS
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo 
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,ths(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! surface potential T  using getgb
      Index=26
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval 	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,ths(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
      
      do j=jsta,jend
       do i=1,im
        ths(i,j)=ths(i,j)*(p1000/pint(i,j,lp1))
     +	        **CAPA ! convert to THS
        if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' psfc = ',i,j,ths(i,j)
     +   ,pint(i,j,lp1)
       end do
      end do 

! GFS does not have surface specific humidity
      QS=SPVAL           

! GFS does not have inst sensible heat flux
      twbs=SPVAL   
      
! GFS does not have inst latent heat flux
      qwbs=SPVAL    
          
!  GFS does not have time step and physics time step, make up ones since they
! are not really used anyway
      NPHS=1
      DT=1.
      DTQ2 = DT * NPHS  !MEB need to get physics DT
      TSPH = 3600./DT   !MEB need to get DT
! All GFS time-averaged quantities are in 6 hour bucket
      TPREC=6.0

! convective precip in m per physics time step using gfsio
!      VarName='cprat'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!	else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)*dtq2/1000. ! convert to m
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo  
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , avgcprate(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! convective precip in m per physics time step using getgb
      Index=272
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)*dtq2/1000. ! convert to m
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,avgcprate(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)      
      
      do j=jsta,jend
        do i=1,im
	  cprate(i,j)=avgcprate(i,j) !set avg cprate to inst one to derive other fields
	end do
       end do	  
       print*,'maxval CPRATE: ', maxval(CPRATE)

! precip rate in m per physics time step using gfsio  
!      VarName='prate'
!      VcoordName='sfc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!	else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)*dtq2/1000. ! convert to m 
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo 
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,avgprec(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! precip rate in m per physics time step using getgb
      Index=271
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)*dtq2/1000. ! convert to m
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,avgprec(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret) 
      
      do j=jsta,jend
        do i=1,im
	  prec(i,j)=avgprec(i,j) !set avg cprate to inst one to derive other fields
	end do
       end do	

! GFS does not have accumulated total, gridscale, and convective precip, will use inst precip to derive in SURFCE.f

      
! GFS does not have similated precip
      lspa=spval  

! inst snow water eqivalent using gfsio                    
!      VarName='weasd'
!      VcoordName='sfc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,sno(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! inst snow water eqivalent using getgb
      Index=119
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,sno(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)


! snow depth in mm using gfsio                    
!      VarName='snod'
!      VcoordName='sfc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)*1000. ! convert to mm
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,si(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! snow depth in mm using getgb
      Index=224
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)*1000. ! convert to mm
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,si(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret) 

! GFS does not have convective cloud efficiency
      CLDEFI=SPVAL
      
! GFS does not have 10 m theta
      TH10=SPVAL

! GFS does not have 10 m humidity
      Q10=SPVAL

! 2m T using gfsio                    
!      VarName='tmp'
!      VcoordName='2m above gnc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,tshltr(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! 2m T using getgb
      Index=106
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=2
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,tshltr(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)      

! GFS does not have 2m pres, estimate it, also convert t to theta 
      Do j=jsta,jend
        Do i=1,im
          PSHLTR(I,J)=pint(I,J,lm+1)*EXP(-0.068283/tshltr(i,j))
          tshltr(i,j)= tshltr(i,j)*(p1000/PSHLTR(I,J))**CAPA ! convert to theta
          if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample 2m T and P after scatter= '
     +   ,i,j,tshltr(i,j),pshltr(i,j)
        end do
      end do

! 2m specific humidity using gfsio                    
!      VarName='spfh'
!      VcoordName='2m above gnc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,qshltr(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! 2m specific humidity using getgb
      Index=112
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,qshltr(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
     
! GFS does not have TKE because it uses MRF scheme
      Q2=SPVAL
      
! GFS does not have surface exchange coeff
 
! GFS does not have snow free albedo
      ALBASE=SPVAL
	
! mid day avg albedo in fraction using gfsio                   
!      VarName='albdo'
!      VcoordName='sfc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)/100. ! convert to fraction
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + ,avgalbedo(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! mid day avg albedo in fraction using getgb      
      Index=266
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)/100. ! convert to mm
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,avgalbedo(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)	
     
! time averaged column cloud fraction
      Index=144
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)/100. ! convert to mm
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,avgtcdc(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)     	
	
! GFS probably does not use zenith angle
      Czen=spval
      CZMEAN=SPVAL      

! waiting to read in lat lon from GFS soon
!      varname='GLAT'
!      call retrieve_index(index,VarName,varname_all,nrecs,iret)
!      if (iret /= 0) then
!        print*,VarName," not found in file-Assigned missing values"
!        GDLAT=SPVAL
!      else
!        this_offset=file_offset(index+1)+(jsta_2l-1)*4*im
!	this_length=im*(jend_2u-jsta_2l+1)
!        call mpi_file_read_at(iunit,this_offset
!     + ,buf,this_length,mpi_real4
!     + , mpi_status_ignore, ierr)
!        if (ierr /= 0) then
!          print*,"Error reading ", VarName,"Assigned missing values"
!          GDLAT=SPVAL
!        else
!          do j = jsta_2l, jend_2u
!           do i = 1, im
!             F(I,J)=1.454441e-4*sin(buf(I,J))   ! 2*omeg*sin(phi)
!             GDLAT(I,J)=buf(I,J)*RTD
	     
!           enddo
!          enddo
!        end if
!      end if
      
!      varname='GLON'
!      call retrieve_index(index,VarName,varname_all,nrecs,iret)
!      if (iret /= 0) then
!        print*,VarName," not found in file-Assigned missing values"
!        GDLON=SPVAL
!      else
!        this_offset=file_offset(index+1)+(jsta_2l-1)*4*im
!	this_length=im*(jend_2u-jsta_2l+1)
!        call mpi_file_read_at(iunit,this_offset
!     + ,buf,this_length,mpi_real4
!     + , mpi_status_ignore, ierr)
!        if (ierr /= 0) then
!          print*,"Error reading ", VarName,"Assigned missing values"
!          GDLON=SPVAL
!        else
!          do j = jsta_2l, jend_2u
!           do i = 1, im
!             GDLON(I,J)=buf(I,J)*RTD
!	     if(i.eq.409.and.j.eq.835)print*,'GDLAT GDLON in INITPOST='
!     +	     ,i,j,GDLAT(I,J),GDLON(I,J)
!           enddo
!          enddo
!        end if
!      end if
      
!       if(jsta.le.594.and.jend.ge.594)print*,'gdlon(120,594)= ',
!     + gdlon(120,594)

! maximum snow albedo in fraction using gfsio                 
!      VarName='mxsalb'
!      VcoordName='sfc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)/100. ! convert to fraction
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , mxsnal(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! maximum snow albedo in fraction using getgb
      Index=227
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)/100. ! convert to fraction
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,mxsnal(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
     
! GFS does not have inst surface outgoing longwave	
      radot=spval

! GFS probably does not use sigt4, set it to sig*t^4
      Do j=jsta,jend
        Do i=1,im
          TLMH=T(I,J,LM)
          Sigt4(I,j)= 5.67E-8*TLMH*TLMH*TLMH*TLMH
        End do
      End do

! TG is not used, skip it for now

! will retrive f_ice when GFS switches to Ferrier scheme
!      varname='F_ICE'
!      call retrieve_index(index,VarName,varname_all,nrecs,iret)
!      if (iret /= 0) then
!        print*,VarName," not found in file-Assigned missing values"
!        F_ice=SPVAL
!     else
!       this_offset=file_offset(index+1)+(jsta_2l-1)*4*im*lm
!	this_length=im*(jend_2u-jsta_2l+1)*lm
!        call mpi_file_read_at(iunit,this_offset
!     + ,buf3d,this_length,mpi_real4
!     + , mpi_status_ignore, ierr)        
!        if (ierr /= 0) then
!          print*,"Error reading ", VarName,"Assigned missing values"
!          F_ice=SPVAL
!        else
!	  do l = 1, lm
!	   ll=lm-l+1
!           do j = jsta_2l, jend_2u
!            do i = 1, im
!             F_ice( i, j, l ) = buf3d ( i, ll, j )
!	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample F_ice= ',
!     +         i,j,l,F_ice( i, j, l )	     
!            end do
!           end do
!          end do 
!	end if 
!      end if	

!      varname='F_RAIN'
!      call retrieve_index(index,VarName,varname_all,nrecs,iret)
!      if (iret /= 0) then
!        print*,VarName," not found in file-Assigned missing values"
!        F_rain=SPVAL
!      else
!        this_offset=file_offset(index+1)+(jsta_2l-1)*4*im*lm
!	this_length=im*(jend_2u-jsta_2l+1)*lm
!        call mpi_file_read_at(iunit,this_offset
!     + ,buf3d,this_length,mpi_real4
!     + , mpi_status_ignore, ierr)      
!        if (ierr /= 0) then
!          print*,"Error reading ", VarName,"Assigned missing values"
!          F_rain=SPVAL
!        else
!	  do l = 1, lm
!	   ll=lm-l+1
!           do j = jsta_2l, jend_2u
!            do i = 1, im
!             F_rain( i, j, l ) = buf3d ( i, ll, j )
!	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample F_rain= ',
!     +         i,j,l,F_rain( i, j, l )	     
!            end do
!           end do
!          end do 
!	end if 
!      end if

!      varname='F_RIMEF'
!      call retrieve_index(index,VarName,varname_all,nrecs,iret)
!      if (iret /= 0) then
!        print*,VarName," not found in file-Assigned missing values"
!        F_RimeF=SPVAL
!      else
!       this_offset=file_offset(index+1)+(jsta_2l-1)*4*im*lm
!	this_length=im*(jend_2u-jsta_2l+1)*lm
!        call mpi_file_read_at(iunit,this_offset
!     + ,buf3d,this_length,mpi_real4
!     + , mpi_status_ignore, ierr)
!        if (ierr /= 0) then
!          print*,"Error reading ", VarName,"Assigned missing values"
!          F_RimeF=SPVAL
!        else
!	  do l = 1, lm
!	   ll=lm-l+1
!           do j = jsta_2l, jend_2u
!            do i = 1, im
!             F_RimeF( i, j, l ) = buf3d ( i, ll, j )
!	     if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,
!     +         'sample F_RimeF= ',i,j,l,F_RimeF( i, j, l )	     
!            end do
!           end do
!          end do 
!	end if 
!      end if

! GFS does not have model level cloud fraction
      CFR=SPVAL

! ask murthy if there is snow rate in GFS
!      varname='SR'
!      call retrieve_index(index,VarName,varname_all,nrecs,iret)
!      if (iret /= 0) then
!        print*,VarName," not found in file-Assigned missing values"
!        SR=SPVAL
!      else
!        this_offset=file_offset(index+1)+(jsta_2l-1)*4*im
!	this_length=im*(jend_2u-jsta_2l+1)
!        call mpi_file_read_at(iunit,this_offset
!     + ,sr,this_length,mpi_real4
!     + , mpi_status_ignore, ierr)
!        if (ierr /= 0) then
!          print*,"Error reading ", VarName,"Assigned missing values"
!          SR=SPVAL
!        end if
!      end if	

! GFS could fields are averaged over 6 hour bucket
      TCLOD=6.0

! GFS does not have inst cloud fraction for high, middle, and low cloud
      cfrach=spval
      cfracl=spval
      cfracm=spval

! ave high cloud fraction using gfsio
!      VarName='tcdc'
!      VcoordName='high cld lay' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)/100. ! convert to fraction
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , avgcfrach(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop 

! ave high cloud fraction using getgb
      Index=302
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)/100. ! convert to fraction
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,avgcfrach(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)

! low cloud fraction using gfsio
!      VarName='tcdc'
!      VcoordName='low cld lay' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)/100. ! convert to fraction
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , avgcfracl(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop      

! ave low cloud fraction using getgb
      Index=300
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)/100. ! convert to fraction
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,avgcfracl(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)

! mid cloud fraction using gfsio
!      VarName='tcdc'
!      VcoordName='mid cld lay' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)/100. ! convert to fraction
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , avgcfracm(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop      
!      write(6,*) 'maxval CFRACM: ', maxval(CFRACM)
      
! ave middle cloud fraction using getgb
      Index=301
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)/100. ! convert to fraction
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,avgcfracm(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)      

! inst convective cloud fraction using getgb
      Index=301
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=71
      jpds(6)=244
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,cnvcfr(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)

! slope type using gfsio
!      VarName='sltyp'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , islope(1,jsta),icnt(me),mpi_integer,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop 

! slope type using getgb
      Index=223
      VarName=avbl(index)
      jpds=-1.0
       jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,islope(1,jsta),icnt(me),mpi_integer,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop 
      
      do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,islope(i,j)
        end do
       end do      


! plant canopy sfc wtr in m using gfsio
!      VarName='cnwat'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)/1000. ! convert from kg*m^2 to m
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , cmc(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop      

! plant canopy sfc wtr in m using getgb
      Index=118
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0
      if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
       if (iret /= 0) then
        print*,VarName," not found in file-Assigned missing values"
        dummy=spval
       else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)/1000. ! convert from kg*m^2 to m
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 	
       end if
      end if                                                                                       
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,cmc(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)

! GFS does not have inst ground heat flux
      grnflx=spval    

! GFS does not have snow cover yet
!      VarName='gflux'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , pctsno(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! asuume tg3 in GFS is the same as soiltb in wrf nmm. It's in sfc file, will
! be able to read it when it merges to gfs io
      VarName='tg3'
      VcoordName='sfc' 
      l=1
      if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
      end if	
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + , soiltb(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      	
! vegetation fraction in fraction. It's in sfc file, will
! be able to read it when it merges to gfs io
      VarName='veg'
      VcoordName='sfc' 
      l=1
      if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
      end if	
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + , vegfrc(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! liquid volumetric soil mpisture in fraction using gfsio
!      VarName='soill'
!      VcoordName='soil layer' 
!      do l=1,nsoil
!       if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!       end if	
!       call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , sh2o(1,jsta,l),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!       if (iret /= 0)print*,'Error scattering array';stop
!      end do 
! GFS doesn not yet output soil layer thickness, assign SLDPTH to be the same as eta

         SLDPTH(1)=0.10
         SLDPTH(2)=0.3
         SLDPTH(3)=0.6
         SLDPTH(4)=1.0
	 
! liquid volumetric soil mpisture in fraction using getgb
      Index=225
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      do l=1,nsoil
       if(l == 1)then
        jpds(7)=nint(sldpth(1)*100.)
       else 
        soilayert=0 
        do n=1,l-1
         soilayert=soilayert+sldpth(n)*100.
	end do
	soilayerb=soilayert+sldpth(l)*100. 
        jpds(7)=nint(soilayert*256.+soilayerb) 
       end if
!       print*,'soil jpds7= ',jpds(7)	                                                                                  
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if
                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,sh2o(1,jsta,l),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,l,sh2o(i,j,l)
        end do
       end do
       
      End do ! do loop for l     

! volumetric soil moisture using gfsio
!      VarName='soilw'
!      VcoordName='soil layer' 
!      do l=1,nsoil
!       if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!       end if	
!       call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , smc(1,jsta,l),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!       if (iret /= 0)print*,'Error scattering array';stop
!      end do 

! volumetric soil moisture using getgb
      Index=117
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      do l=1,nsoil
       if(l == 1)then
        jpds(7)=nint(sldpth(1)*100.)
       else 
        soilayert=0 
        do n=1,l-1
         soilayert=soilayert+sldpth(n)*100.
	end do
	soilayerb=soilayert+sldpth(l)*100. 
        jpds(7)=nint(soilayert*256.+soilayerb) 
       end if                                                                             
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if
                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,smc(1,jsta,l),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
     
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,l,smc(i,j,l)
        end do
       end do 
     
      End do ! do loop for l


! soil temperature using gfsio
!      VarName='tmp'
!      VcoordName='soil layer' 
!      do l=1,nsoil
!       if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!       end if	
!       call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , stc(1,jsta,l),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!       if (iret /= 0)print*,'Error scattering array';stop
!      end do 

! soil temperature using getgb
      Index=116
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=11  ! GFS used 11 for soil T instead of 85
      jpds(6)=is(index)
      do l=1,nsoil
       if(l == 1)then
        jpds(7)=nint(sldpth(1)*100.)
       else 
        soilayert=0 
        do n=1,l-1
         soilayert=soilayert+sldpth(n)*100.
	end do
	soilayerb=soilayert+sldpth(l)*100. 
        jpds(7)=nint(soilayert*256.+soilayerb) 
       end if                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if
                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,stc(1,jsta,l),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,l,stc(i,j,l)
        end do
       end do 
      End do ! do loop for l

     
! GFS does not output time averaged convective and strat cloud fraction, set acfrcv to spval, ncfrcv to 1
      acfrcv=spval
      ncfrcv=1.0
! GFS does not output time averaged cloud fraction, set acfrst to spval, ncfrst to 1
      acfrst=spval
      ncfrst=1.0

! GFS does not have storm runoff
      ssroff=spval

! GFS does not have UNDERGROUND RUNOFF
      bgroff=spval

! GFS incoming sfc longwave has been averaged over 6 hr bucket, set ARDLW to 1
      ardlw=1.0
      trdlw=6.0

! GFS does not have inst incoming sfc longwave
      rlwin=spval
       
! GFS does not have inst model top outgoing longwave
      rlwtoa=spval

! time averaged incoming sfc longwave using gfsio
!      VarName='dlwrf'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , alwin(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! time averaged incoming sfc longwave using getgb
      Index=127
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,alwin(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,alwin(i,j)
        end do
       end do
   
      
! time averaged outgoing sfc longwave using gfsio
!      VarName='ulwrf'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , alwout(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! time averaged outgoing sfc longwave using getgb
      Index=129
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,alwout(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,alwout(i,j)
        end do
       end do      

! time averaged outgoing model top longwave using gfsio
!      VarName='ulwrf'
!      VcoordName='nom. top' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , alwtoa(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! time averaged outgoing model top longwave using getgb
      Index=131
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,alwtoa(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
      
! GFS does not have inst incoming sfc shortwave
      rswin=spval 

! GFS does not have inst incoming clear sky sfc shortwave
      rswinc=spval      

! GFS does not have inst outgoing sfc shortwave
      rswout=spval
           
! GFS incoming sfc longwave has been averaged, set ARDLW to 1
      ardsw=1.0
      trdsw=6.0

! time averaged incoming sfc shortwave using gfsio
!      VarName='dswrf'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , aswin(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! time averaged incoming sfc shortwave using getgb
      Index=126
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,aswin(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret) 
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,aswin(i,j)
        end do
       end do     

! time averaged incoming sfc uv-b using getgb
      Index=298
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,auvbin(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret) 
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,auvbin(i,j)
        end do
       end do  
       
! time averaged incoming sfc clear sky uv-b using getgb
      Index=297
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,auvbinc(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret) 
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,auvbinc(i,j)
        end do
       end do          
      
! time averaged outgoing sfc shortwave using gfsio
!      VarName='uswrf'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , aswout(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! time averaged outgoing sfc shortwave using getgb
      Index=128
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,aswout(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,aswout(i,j)
        end do
       end do
      
! time averaged model top outgoing sfc shortwave
!      VarName='uswrf'
!      VcoordName='nom. top' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , aswtoa(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! time averaged model top outgoing sfc shortwave      
      Index=130
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,aswtoa(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,aswtoa(i,j)
        end do
       end do      

! time averaged surface sensible heat flux, multiplied by -1 because wrf model flux
! has reversed sign convention using gfsio
!      VarName='shtfl'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!	else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)*-1.0 
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo 
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , sfcshx(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! time averaged surface sensible heat flux, multiplied by -1 because wrf model flux
! has reversed sign convention using getgb
      Index=43
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
	else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)*(-1.0) 
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,sfcshx(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)


! GFS surface flux has been averaged, set  ASRFC to 1 
      asrfc=1.0  
      tsrfc=6.0

! time averaged surface latent heat flux, multiplied by -1 because wrf model flux
! has reversed sign vonvention using gfsio
!      VarName='lhtfl'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!	else
!         do j = 1, jm
!           do i = 1, im
!             dummy(I,J)= dummy(i,j)*-1.0 
!             if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
!           enddo
!          enddo 
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , sfclhx(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop  

! time averaged surface latent heat flux, multiplied by -1 because wrf model flux
! has reversed sign vonvention using getgb
      Index=42
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
	else
         do j = 1, jm
           do i = 1, im
             dummy(I,J)= dummy(i,j)*(-1.0) 
             if (j.eq.jm/2 .and. mod(i,50).eq.0)
     + print*,'sample ',VarName, ' = ',i,j,dummy(i,j)     
           enddo
          enddo 
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,sfclhx(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
      
! time averaged ground heat flux using gfsio
!      VarName='gflux'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , subshx(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! time averaged ground heat flux using getgb      
      Index=135
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,subshx(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)      

! GFS does not have snow phase change heat flux
      snopcx=spval

! time averaged zonal momentum flux using gfsio
!      VarName='uflx'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , sfcux(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! time averaged zonal momentum flux using getgb      
      Index=269
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,sfcux(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,sfcux(i,j)
        end do
       end do
      
! time averaged meridional momentum flux using gfsio
!      VarName='vflx'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , sfcvx(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! time averaged meridional momentum flux using getgb
      Index=270
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,sfcvx(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
     
! GFS does not use total momentum flux
      sfcuvx=spval

! time averaged accumulated potential evaporation
      Index=242
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,potevp(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
      

! GFS does not have temperature tendency due to long wave radiation
      rlwtt=spval
      
! GFS does not have temperature tendency due to long wave radiation
      rswtt=spval
      
! GFS does not have temperature tendency due to latent heating from convection
      tcucn=spval

! set avrain to 1
      avrain=1.0
      avcnvc=1.0
      theat=6.0 ! just in case GFS decides to output T tendency   
      
! GFS does not have temperature tendency due to latent heating from grid scale
      train=spval

! set theat to 6.0 for now even though it won't be used
      theat=6.0	

! 10 m u using gfsio
!      VarName='ugrd'
!      VcoordName='10 m above gr' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , u10(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! 10 m u using getgb
      Index=64
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=10                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,u10(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
       do j=jsta,jend
        do i=1,im
         if (j.eq.jm/2 .and. mod(i,50).eq.0)
     +   print*,'sample ',trim(VarName), ' after scatter= '
     +   ,i,j,u10(i,j)
        end do
       end do
            
! 10 m v using gfsio
!      VarName='vgrd'
!      VcoordName='10 m above gr' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , v10(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! 10 m v using getgb
      Index=65
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=10                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,v10(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
      
! GFS does not have soil moisture availability 
      smstav=spval

! GFS does not have total soil moisture 
      smstot=spval
      
! vegetation type, it's in GFS surface file, hopefully will merge into gfsio soon 
      VarName='vtype'
      VcoordName='sfc' 
      l=1
      if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
      end if	
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + , ivgtyp(1,jsta),icnt(me),mpi_integer,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! soil type, it's in GFS surface file, hopefully will merge into gfsio soon
      VarName='stype'
      VcoordName='sfc' 
      l=1
      if(me == 0)then
        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
     +	,l,dummy,iret=iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
      end if	
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + , isltyp(1,jsta),icnt(me),mpi_integer,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! GFS does not have accumulated surface evaporation
      sfcevp=spval

! GFS does not have surface exchange coeefficient
      sfcexc=spval
      
! GFS does not have averaged accumulated snow
      acsnow=spval

! GFS does not have snow melt
      acsnom=spval
       
! GFS does not have sst????
      sst=spval

! GFS does not have mixing length
      el_myj=spval      

! GFS does not output exchange coefficient
      exch_h=spval
      
! GFS does not have THZ0, use THS to substitute
      do j=jsta,jend
        do i=1,im
	  thz0(i,j)=ths(i,j)
	end do
      end do      	   
      if(jj.ge.jsta.and.jj.le.jend)
     + print*,'THZ0 at ',ii,jj,' = ',THZ0(ii,jj)

! GFS does not output humidity at roughness length
      qz0=spval
      
! GFS does not output u at roughness length
      uz0=spval
      
! GFS does not output humidity at roughness length
      vz0=spval      

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

! retrieve inst convective cloud top, GFS has cloud top pressure instead of index,
! will need to modify CLDRAD.f to use pressure directly instead of index
!      VarName='pres'
!      VcoordName='convect-cld top' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , ptop(1,jsta),icnt(me),mpi_integer,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! retrieve inst convective cloud top using getgb
      Index=189
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,ptop(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
        print*,'maxval PTOP: ', maxval(PTOP)

      htop=0.0	
      do j=jsta,jend
        do i=1,im
	  do l=1,lm
	    if(ptop(i,j) <= 0.0)ptop(i,j)=spval
	    if(ptop(i,j) <= pmid(i,j,l))then
	     htop(i,j)=l
	     if(i==ii .and. j==jj)print*,'sample ptop,pmid pmid-1,pint= ',
     +	        ptop(i,j),pmid(i,j,l),pmid(i,j,l-1),pint(i,j,l),htop(i,j)
             go to 48
	    end if
	  end do
 48	  continue   
        end do
       end do

! retrieve inst convective cloud bottom, GFS has cloud top pressure instead of index,
! will need to modify CLDRAD.f to use pressure directly instead of index
!      VarName='pres'
!      VcoordName='convect-cld bot' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecv(gfile,trim(VarName),trim(VcoordName)
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real
!     + , pbot(1,jsta),icnt(me),mpi_integer,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! retrieve inst convective cloud bottom using getgb
      Index=188
      VarName=avbl(index)
      jpds=-1.0
      jgds=-1.0
      jpds(5)=iq(index)
      jpds(6)=is(index)
      jpds(7)=0                                                                                    
       if(me == 0)then
        call getgb(iunit,0,im_jm,0,jpds,jgds,kf
     +    ,k,kpds,kgds,lb,dummy,iret)
        if (iret /= 0) then
         print*,VarName," not found in file-Assigned missing values"
         dummy=spval
        end if
       end if                                                                                    
       call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real
     + ,pbot(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
        print*,'maxval PBOT: ', maxval(PBOT)
      
      hbot=0.0 
      do j=jsta,jend
        do i=1,im
	  if(pbot(i,j) <= 0.0)pbot(i,j)=spval
	  do l=lm,1,-1
	    if(pbot(i,j) >= pmid(i,j,l))then
	     hbot(i,j)=l
	     if(i==ii .and. j==jj)print*,'sample pbot,pmid= ',
     +	        pbot(i,j),pmid(i,j,l),hbot(i,j)
             go to 49
	    end if
	  end do
 49	  continue   
        end do
       end do	    
! GFS does not have deep convective cloud top and bottom fields
      HTOPD=SPVAL
      HBOTD=SPVAL   
      HTOPS=SPVAL
      HBOTS=SPVAL 
      CUPPT=SPVAL 

!!!! DONE GETTING

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            IF(ABS(T(I,J,L)).GT.1.0E-3)
     &        OMGA(I,J,L) = -WH(I,J,L)*PMID(I,J,L)*G/
     &                 (RD*T(I,J,L)*(1.+D608*Q(I,J,L)))

        end do
       end do
      end do

    
! pos east
       call collect_loc(gdlat,dummy)
       if(me.eq.0)then
        latstart=nint(dummy(1,1)*1000.)
        latlast=nint(dummy(im,jm)*1000.)
       end if
       write(6,*) 'laststart,latlast B calling bcast= '
     +, latstart,latlast
       call mpi_bcast(latstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(latlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*) 'laststart,latlast A calling bcast= '
     +, latstart,latlast
       call collect_loc(gdlon,dummy)
       if(me.eq.0)then
        lonstart=nint(dummy(1,1)*1000.)
        lonlast=nint(dummy(im,jm)*1000.)
       end if
       write(6,*)'lonstart,lonlast B calling bcast= '
     +, lonstart,lonlast
       call mpi_bcast(lonstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(lonlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*)'lonstart,lonlast A calling bcast= '
     +, lonstart,lonlast
!
!        ncdump -h

!!
!! 
!!
        write(6,*) 'filename in INITPOST=', filename,' is'

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
!            DX ( i, j ) = dxval  !MEB ???
!            DY ( i, j ) = dyval*DTR*ERAD  

!!!!!!!!!!!!!!!!!!!!! DY ????

            DY ( i, j ) =   0.001*ERAD*DYVAL*DTR  ! like A*DPH
        end do
       end do
!MEB not sure how to get these 

! close up shop
!      call ext_int_ioclose ( DataHandle, Status )

! generate look up table for lifted parcel calculations

      THL=210.
      PLQ=70000.

      CALL TABLE(PTBL,TTBL,PT,
     &          RDQ,RDTH,RDP,RDTHE,PL,THL,QS0,SQS,STHE,THE0)

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
!      NPHS = 4  !MEB need to get physics DT
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
!     
!

      RETURN
      END


