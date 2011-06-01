      SUBROUTINE INITPOST_NMM_BIN
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
      use wrf_io_flags_mod
      use params_mod
      use lookup_mod
      use ctlblk_mod
      use gridspec_mod
! Used for data inventory
      use module_internal_header_util ! open for data inventory

      implicit none
!     
      INCLUDE "mpif.h"
! This version of INITPOST shows how to initialize, open, read from, and
! close a NetCDF dataset. In order to change it to read an internal (binary)
! dataset, do a global replacement of _ncd_ with _int_. 

      character(len=31) :: VarName
      integer :: Status
      character startdate*132,SysDepInfo*80
! 
!     NOTE: SOME INTEGER VARIABLES ARE READ INTO DUMMY ( A REAL ). THIS IS OK
!     AS LONG AS REALS AND INTEGERS ARE THE SAME SIZE.
!
!     ALSO, EXTRACT IS CALLED WITH DUMMY ( A REAL ) EVEN WHEN THE NUMBERS ARE
!     INTEGERS - THIS IS OK AS LONG AS INTEGERS AND REALS ARE THE SAME SIZE.
      LOGICAL RUNB,SINGLRST,SUBPOST,NEST,HYDRO
      LOGICAL IOOMG,IOALL
      CHARACTER*32 LABEL
      CHARACTER*40 CONTRL,FILALL,FILMST,FILTMP,FILTKE,FILUNV   &
     &, FILCLD,FILRAD,FILSFC
      CHARACTER*4 RESTHR
      CHARACTER FNAME*80,ENVAR*50,BLANK*4
      INTEGER IDATB(3),IDATE(8),JDATE(8)
!     
!     DECLARE VARIABLES.
!     
      REAL SLDPTH2(NSOIL)
      REAL RINC(5)
      REAL ETA1(LM), ETA2(LM)
      REAL DUM1D (LM+1)
      REAL DUMMY ( IM, JM )
      REAL DUMMY2 ( IM, JM )
      REAL FI(IM,JM,2)
      INTEGER IDUMMY ( IM, JM )
      REAL DUM3D ( IM+1, JM+1, LM+1 )
      REAL DUM3D2 ( IM+1, JM+1, LM+1 )
!NCAR (jw)
      integer ii,jj,js,je,jev,iyear,imn,iday,itmp,ioutcount,istatus,    &
              nsrfc,nrdlw,nrdsw,nheat,nclod,                            &
              I,J,L,LL,N,LONEND,LATEND,IMM,INAV,IRTN,                   &
              IFDX,IFDY,IGDOUT,ICEN,JCEN
      real TSPH,fact,dumcst,tstart,tmp
      integer isf_physics, isf_fphysics

!
! Declarations for  :
!  Comments and code provided for use with copygb - R.Rozumalski - NWS
      INTEGER idxave, dlat, latnm, latsm, nlat, lonem, lonwm, dlon, nlon

      DATA BLANK/'    '/

!mhu for variables survey
      REAL, DIMENSION(6)     :: tmp_array
      INTEGER, DIMENSION(6)     :: itmp_array
      INTEGER, DIMENSION(512)     :: hdrbuf
      INTEGER  hdrbufsize, itypesize,code
      CHARACTER*79 locElement,dumstr,Element
      INTEGER  rtypesize, IALL
      LOGICAL keepgoing

      INTEGER                        :: locDataHandle
      CHARACTER(len=79)              :: locDateStr
      CHARACTER*(79)                 :: locVarName
      integer                        :: locFieldType
      integer                        :: locComm
      integer                        :: locIOComm
      integer                        :: locDomainDesc
      character*132                  :: locMemoryOrder
      character*132                  :: locStagger
      character*132 , dimension (3)  :: locDimNames
      integer ,dimension(3)          :: locDomainStart, locDomainEnd
      integer ,dimension(3)          :: locMemoryStart, locMemoryEnd
      integer ,dimension(3)          :: locPatchStart,  locPatchEnd

      REAL, DIMENSION(10)    :: Field

      character*80     :: titlestring, version
      integer          :: istart,iend,istat

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
!     fileName = '/ptmp/wx20mb/wrfout_01_030500'
!     DateStr = '2002-03-05_18:00:00'
!  how do I get the filename?
      call ext_int_ioinit(SysDepInfo,Status)
      print*,'called ioinit', Status
      call ext_int_open_for_read( trim(fileName), 0, 0, " ",         & 
           DataHandle, Status)
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
! ===data inventory
      if(2==1) then
       CALL wrf_sizeof_integer( itypesize )
       CALL wrf_sizeof_real   ( rtypesize )
       hdrbufsize = 0
       ioutcount=0
       locDataHandle=0
       itmp_array = 0
       keepgoing = .true.
       Element='GGGGGGGTITLE'
       VarName='GGG'
       IALL=0
       locDateStr="20080111"
       locVarName='AAA'
       Field=0.0
       locFieldType=0
       locDomainStart=1
!CCCC       locDomainEnd=10
       locMemoryStart=1
       locMemoryEnd=10
       locPatchStart=1
       locPatchEnd=10
       DO WHILE ( keepgoing )
        IALL=IALL+1
!!        write(*,*) 'Here =',IALL
        if(IALL == 300) keepgoing = .false.
        READ( unit=DataHandle , iostat = istat ) hdrbuf
        IF ( istat .EQ. 0 ) THEN
         code = hdrbuf(2)
!!         write(*,*) 'Code=', code
         IF ( code .EQ. 220 ) THEN
           CALL int_get_ti_header_char( hdrbuf, hdrbufsize, itypesize,  &
                 locDataHandle, locElement, dumstr, startdate, code )
           IF ( TRIM(locElement) .EQ. TRIM(Element) ) THEN
              keepgoing = .false. ;  Status = 0
           ENDIF
!!            write(*,*)  'startdate=', startdate
!!            write(*,*)  hdrbufsize,'=', trim(locElement)
             write(*,*) IALL, code, trim(locElement)
         ELSEIF ( code .EQ. 140 ) THEN
            CALL int_get_ti_header_real(hdrbuf,hdrbufsize,itypesize,    &
              rtypesize,locDataHandle,locElement,                       &
              tmp_array,ioutcount,code )
            IF ( TRIM(locElement) .EQ. TRIM(Element) ) THEN
              keepgoing = .false. ;  Status = 0
            ENDIF
!!           write(*,*)  ioutcount,'=',tmp_array
!!           write(*,*)  hdrbufsize,'=', trim(locElement)
             write(*,*) IALL, code, trim(locElement)
         ELSEIF ( code .EQ. 180 ) THEN
            CALL int_get_ti_header_integer(hdrbuf,hdrbufsize,itypesize, &
                   rtypesize,locDataHandle,locElement,                  &
                   itmp_array,ioutcount,code )
            IF ( TRIM(locElement) .EQ. TRIM(Element) ) THEN
              keepgoing = .false. ;  Status = 0
            ENDIF
!!             write(*,*)  ioutcount,'=',itmp_array
!!            write(*,*)  hdrbufsize,'=', trim(locElement)
             write(*,*) IALL, code, trim(locElement)
          ELSEIF(code .EQ. 530) then
           CALL int_get_write_field_header(hdrbuf,hdrbufsize,itypesize, &
                   rtypesize, locDataHandle,locDateStr, locVarName,     &
                   tmp_array, locFieldType , locComm , locIOComm,       &
                   locDomainDesc,locMemoryOrder,locStagger,locDimNames, &
                   locDomainStart,locDomainEnd, locMemoryStart,         &
                   locMemoryEnd, locPatchStart,locPatchEnd)
              IF ( TRIM(locVarName) .EQ. TRIM(VarName) ) THEN
                 write( *,*) TRIM(locVarName)
              endif
!!              write( *,*) locFieldType
!!              write(*,*)  hdrbufsize,'=', trim(locVarName)
              write(*,*) IALL, code, trim(locVarName), locFieldType
              READ( unit=DataHandle , iostat = istat )
          ENDIF
        ENDIF
       ENDDO
      stop 123
      endif
!di END OF WHAT WAS Cdi 

      print*,'HERE'

!CCC ==============================

      call ext_int_get_dom_ti_char(DataHandle,'TITLE',titlestring, status)
        print*,'TITLE= ',trim(titlestring)

!
! Getting start time
      print*,'startdate into ext_int_get_dom_ti_char ',startdate
      call ext_int_get_dom_ti_char(DataHandle,'START_DATE',startdate, status)
      print*,'startdate= ',startdate
      jdate=0
      idate=0
      read(startdate,15)iyear,imn,iday,ihrst,imin       
 15   format(i4,1x,i2,1x,i2,1x,i2,1x,i2)
      print*,'start yr mo day hr min =',iyear,imn,iday,ihrst,imin
      print*,'processing yr mo day hr min='                            &
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
!      CALL W3DIFDAT(JDATE,IDATE,2,RINC)
!      ifhr=nint(rinc(2))
      CALL W3DIFDAT(JDATE,IDATE,0,RINC)
      ifhr=nint(rinc(2)+rinc(1)*24.)
      ifmin=nint(rinc(3))
      print*,' in INITPOST ifhr ifmin fileName=',ifhr,ifmin,fileName

! Getting tstart
      if(2==1) then   ! comment out by ming hu
        tstart=0.
        call ext_int_get_dom_ti_real(DataHandle,'TSTART',tmp              &  
          ,1,ioutcount,istatus)
        tstart=tmp 
        print*,'status for getting TSTART= ',istatus 
        IF( abs(istatus-0) .GT. 1)THEN
         PRINT*,'you do not have tstart in your WRF output,'//            &
           'many grid navigation will be read in incorrectly, STOPPING'       
         STOP   
        END IF
        print*,'status for getting TSTART= ',istatus     
        print*,'TSTART= ',TSTART 
      
! Getting restart
      
        RESTRT=.TRUE.  ! set RESTRT as default
!        call ext_int_get_dom_ti_integer(DataHandle,'RESTARTBIN',itmp
!       + ,1,ioutcount,istatus)
      
!        IF(itmp .LT. 1)THEN
!          RESTRT=.FALSE.
!        ELSE
!          RESTRT=.TRUE.
!        END IF
     
!        print*,'status for getting RESTARTBIN= ',istatus
!        print*,'Is this a restrt run? ',RESTRT
            
        IF(RESTRT)THEN
         ifhr=ifhr+NINT(tstart)
         print*,'new forecast hours for restrt run= ',ifhr
        END IF
      endif ! mhu

! FROM UPP - but TSTART commented out before
!      IF(tstart .GT. 1.0E-2)THEN
!       ifhr=ifhr+NINT(tstart)
!       rinc=0
!       idate=0
!       rinc(2)=-1.0*ifhr
!       call w3movdat(rinc,jdate,idate)
!       SDAT(1)=idate(2)
!       SDAT(2)=idate(3)
!       SDAT(3)=idate(1)
!       IHRST=idate(5)       
!       print*,'new forecast hours for restrt run= ',ifhr
!       print*,'new start yr mo day hr min =',sdat(3),sdat(1)            &  
!             ,sdat(2),ihrst,imin
!      END IF
      
      
!  OK, since all of the variables are dimensioned/allocated to be
!  the same size, this means we have to be careful int getVariable
!  to not try to get too much data.  For example, 
!  DUM3D is dimensioned IM+1,JM+1,LM+1 but there might actually
!  only be im,jm,lm points of data available for a particular variable.  
! get metadata

      call ext_int_get_dom_ti_integer(DataHandle,'MP_PHYSICS'           &
       ,itmp,1,ioutcount,istatus)
      imp_physics=itmp
      if (imp_physics .eq. 85) imp_physics = 5   ! HWRF
      print*,'MP_PHYSICS= ',imp_physics

      call ext_int_get_dom_ti_integer(DataHandle,'SF_SURFACE_PHYSICS'   &
       ,itmp,1,ioutcount,istatus)
      isf_physics=itmp
      print*,'SF_PHYSICS= ',isf_fphysics

      call ext_int_get_dom_ti_real(DataHandle,'DX',tmp                &
          ,1,ioutcount,istatus)
      dxval=nint(tmp*1000.) ! E-grid dlamda in degree
      write(6,*) 'dxval= ', dxval
      call ext_int_get_dom_ti_real(DataHandle,'DY',tmp                &
          ,1,ioutcount,istatus)
      dyval=nint(1000.*tmp)
      write(6,*) 'dyval= ', dyval
      call ext_int_get_dom_ti_real(DataHandle,'DT',tmp                &
          ,1,ioutcount,istatus)
      DT=tmp
      write(6,*) 'DT= ', DT

      call ext_int_get_dom_ti_real(DataHandle,'CEN_LAT',tmp           &
          ,1,ioutcount,istatus)
!      tmp=50.0
      cenlat=nint(1000.*tmp)
      write(6,*) 'cenlat= ', cenlat
      call ext_int_get_dom_ti_real(DataHandle,'CEN_LON',tmp           &
          ,1,ioutcount,istatus)
!      tmp=-111.0
      cenlon=nint(1000.*tmp)
      write(6,*) 'cenlon= ', cenlon
!
! These values are not set correctly 1.e+20 and cause the assignment
!   to int to overflow - tms
!tms      call ext_int_get_dom_ti_real(DataHandle,'TRUELAT1',tmp          &
!tms         ,1,ioutcount,istatus)
!tms      if (istatus .EQ. 0) then
!tms        truelat1=nint(1000.*tmp)
!tms      else
!tms        truelat1=1000
!tms      endif
!tms      write(6,*) 'truelat1= ', truelat1
!tms      call ext_int_get_dom_ti_real(DataHandle,'TRUELAT2',tmp          &
!tms          ,1,ioutcount,istatus)
!tms      truelat2=nint(1000.*tmp)
!tms      write(6,*) 'truelat2= ', truelat2
      call ext_int_get_dom_ti_integer(DataHandle,'MAP_PROJ',itmp      &
          ,1,ioutcount,istatus)
      maptype=itmp
!     maptype=203
      write(6,*) 'maptype is ', maptype

! get 3-D variables
      print*,'im,jm,lm= ',im,jm,lm
!
!tms
! remove this TOYVAR not in this data file
! Come back later and figure out how you should align this
!      VarName='TOYVAR'
!  write(6,*) 'call getVariableB for : ', VarName
!      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,
!     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

!tms  to skip some variables
      DO i=1, 9  !skip a bunch of stuff to align on LU_INDEX
         READ( unit=DataHandle) hdrbuf
         code = hdrbuf(2)
         WRITE(6,*) 'skip variable: ', code
         IF(code .EQ. 530) READ( unit=DataHandle)
      END DO

      VarName='LU_INDEX'
      write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,     &
        IM,1,JM,1,IM,JS,JE,1)
     
! to skip some variables
      DO i=1, 20
         READ( unit=DataHandle) hdrbuf
         code = hdrbuf(2)
         WRITE(6,*) 'skip variable: ', code
         IF(code .EQ. 530) READ( unit=DataHandle)
      END DO

!mhu      VarName='LMH'
!mhu  write(6,*) 'call getIVariable for : ', VarName
!mhu      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,
!mhu     &   IM,1,JM,1,IM,JS,JE,1)

!mhu      VarName='LMV'
!mhu  write(6,*) 'call getIVariable for : ', VarName
!mhu      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,
!mhu     &   IM,1,JM,1,IM,JS,JE,1)

      VarName='HBM2'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          HBM2(I,J)=DUMMY2(I,J)
        enddo
      enddo
  
      VarName='HBM3'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      VarName='VBM2'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      VarName='VBM3'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      VarName='SM'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          SM(I,J)=DUMMY2(I,J)
          if (j.eq.jm/2 .and. mod(i,10).eq.0) print*,'sample SM= ',i,j,sm(i,j)
        enddo
      enddo

      VarName='SICE'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
          SICE(I,J)=DUMMY2(I,J)
        enddo
      enddo

!mhu      VarName='HTM'
!mhu  write(6,*) 'call getVariableB for : ', VarName
!mhu      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,
!mhu     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

!mhu      VarName='VTM'
!mhu  write(6,*) 'call getVariableB for : ', VarName
!mhu      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,
!mhu     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      VarName='PD'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          PD(I,J)=DUMMY2(I,J)
        enddo
      enddo      

!       do j = jsta_2l, jend_2u
!        do i = 1, im
!  PD(I,J)=DUMMY2(I,J)
!        enddo
!       enddo

      VarName='FIS'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          FIS(I,J)=DUMMY2(I,J)
        enddo
      enddo

      VarName='RES'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,     &
              IM,1,JM,1,IM,JS,JE,1)

      VarName='T'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
              IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            T ( i, j, l ) = dum3d ( i, j, l )
	         if(i.eq.im/2.and.j.eq.(jsta+jend)/2)     &
              print*,'sample T= ', i,j,l,T ( i, j, l )
        end do
       end do
      end do

      VarName='Q'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
             IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            q ( i, j, l ) = dum3d ( i, j, l )
!            if(l.eq.1. and. q(i,j,l).lt.1.0e-7)
!     &      print*,'1st level Q < 1.0e-7 in INTPOST'
!     & ,i,j,l,q(i,j,l)
        end do
       end do
      end do
      print*,'finish reading mixing ratio'
      ii=im/2
      jj=(jsta+jend)/2
!      print*,'Q at ',ii,jj,ll,' = ',Q(ii,jj,ll)

      VarName='U'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
              IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            u ( i, j, l ) = dum3d ( i, j, l )
            UH( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do 
      ii=im/2
      jj=(jsta+jend)/2
      ll=lm
!      print*,'U at ',ii,jj,ll,' = ',UH (ii,jj,ll)


      VarName='V'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
              IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            v ( i, j, l ) = dum3d ( i, j, l )
            VH( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do
      print*,'finish reading V'
!      print*,'VH at ',ii,jj,ll,' = ',VH (ii,jj,ll)

      varname='DX_NMM'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
       do i = 1, im
         DX(I,J)=DUMMY2(I,J)
       enddo
      enddo

! to skip some variables
      DO i=1, 2
         READ( unit=DataHandle) hdrbuf
         code = hdrbuf(2)
         WRITE(6,*) 'skip variable: ', code
         IF(code .EQ. 530) READ( unit=DataHandle)
      END DO

      varname='ETA1'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,ETA1,       &   
              LM,1,1,1,LM,1,1,1)

! to skip some variables
      DO i=1, 2
         READ( unit=DataHandle) hdrbuf
         code = hdrbuf(2)
         WRITE(6,*) 'skip variable: ', code
         IF(code .EQ. 530) READ( unit=DataHandle)
      END DO

      varname='ETA2'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,ETA2,       &
              LM,1,1,1,LM,1,1,1)

      open(75,file='ETAPROFILE.txt',form='formatted', status='unknown')
      DO L=1,lm+1 
        IF(L .EQ. 1)THEN
          write(75,1020)L, 0., 0.
        ELSE 
          write(75,1020)L, ETA1(lm+2-l), ETA2(lm+2-l)
        END IF     
!         print*,'L, ETA1, ETA2= ',L, ETA1(l), ETA2(l)
      END DO
 1020 format(I3,2E17.10)  
      close (75)

! to skip some variables
      DO i=1, 2
         READ( unit=DataHandle) hdrbuf
         code = hdrbuf(2)
         WRITE(6,*) 'skip variable: ', code
         IF(code .EQ. 530) READ( unit=DataHandle)
      END DO
  
      varname='PDTOP'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              1,1,1,1,1,1,1,1)
	   write(6,*) 'PDTOP= ', DUMMY2(1,1)
      PDTOP=DUMMY2(1,1)

      varname='PT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'PT ', DUMMY2(1,1)
      PT=DUMMY2(1,1)

      varname='PBLH'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          PBLH(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='MIXHT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      varname='USTAR'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,     &
              IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          USTAR(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='Z0'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          Z0(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='THS'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          THS(I,J)=DUMMY2(I,J)
        enddo
      enddo

      VarName='QS'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
              IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          QS ( i, j ) = dummy ( i, j )
        end do
      end do
      print*,'QS at ',ii,jj,' = ',QS(ii,jj)

      varname='TWBS'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          TWBS(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='QWBS'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          QWBS(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='TAUX' 
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)
	
      varname='TAUY' 
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,     &
              IM,1,JM,1,IM,JS,JE,1)

      varname='PREC' ! instantaneous precip rate?
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,     &
              IM,1,JM,1,IM,JS,JE,1)
  
       do j = jsta_2l, jend_2u
         do i = 1, im
           PREC(I,J)=DUMMY2(I,J)
         enddo
       enddo

      varname='APREC' ! instantaneous precip rate?
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      varname='ACPREC' ! accum total precip
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          ACPREC(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='CUPREC' ! accum cumulus precip
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          CUPREC(I,J)=DUMMY2(I,J)
          ANCPRC(I,J)=ACPREC(I,J)-CUPREC(I,J)
        enddo
      enddo

! hoping to read instantanous convective precip rate soon, initialize it to spval
! for now

!      do j = jsta_2l, jend_2u
!       do i = 1, im
!        CPRATE(I,J)=SPVAL
!       enddo
!      enddo

      varname='LSPA'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
           LSPA(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='SNO'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          SNO(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='SI'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          SI(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='CLDEFI'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
           CLDEFI(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='TH10'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          TH10(I,J)=DUMMY2(I,J)
        enddo
      enddo
       
      varname='Q10'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          Q10(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='PSHLTR'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          PSHLTR(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='TSHLTR'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          TSHLTR(I,J)=DUMMY2(I,J)
          if(i.eq.409.and.j.eq.835)print*,'2m T and P in INITPOST=',    &
              i,j,TSHLTR(i,j),PSHLTR(I,J)
        enddo
      enddo

      varname='QSHLTR'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
              IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          QSHLTR(I,J) = DUMMY2(I,J) 
        enddo
      enddo

      VarName='Q2'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            q2 ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do
      print*,'finish reading TKE'
      print*,'Q2 at ',ii,jj,ll,' = ',Q2(ii,jj,ll)

      varname='AKHS_OUT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
         AKHS(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='AKMS_OUT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
         AKMS(I,J)=DUMMY2(I,J)
        enddo
      enddo
       
!      VarName='T_ADJ'
!  write(6,*) 'call getVariableB for : ', VarName
!      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
!       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!      do l = 1, lm
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            T_ADJ ( i, j, l ) = dum3d ( i, j, l )
!        end do
!       end do
!      end do

      varname='ALBASE'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          ALBASE(I,J)=DUMMY2(I,J)
        enddo
      enddo
  
      varname='ALBEDO'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          ALBEDO(I,J)=DUMMY2(I,J)
        enddo
      enddo

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

      varname='CNVBOT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
     
      do j = jsta_2l, jend_2u
        do i = 1, im
          HBOT ( i, j ) = float(LM)-dummy2(i,j)+1.0
        enddo
      enddo

      varname='CNVTOP'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
     
      do j = jsta_2l, jend_2u
        do i = 1, im
          HTOP ( i, j ) = float(LM)-dummy2(i,j)+1.0
        enddo
      enddo

      varname='CZEN'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          CZEN(I,J)=DUMMY2(I,J)
        enddo
      enddo
      print*,'max CZEN=',maxval(czen)

      varname='CZMEAN'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          CZMEAN(I,J)=DUMMY2(I,J)
        enddo
      enddo
      print*,'max CZMEAN= ',maxval(czmean) 

      varname='EPSR'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      varname='GLAT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,     &
              IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          F(I,J)=1.454441e-4*sin(DUMMY2(I,J))   ! 2*omeg*sin(phi)
!    GDLAT(I,J)=DUMMY2(I,J)*(180./acos(-1.))
          GDLAT(I,J)=DUMMY2(I,J)*RTD
        enddo
       enddo
       if(jsta.le.594.and.jend.ge.594)print*,'gdlat(120,594)= ',         &
          gdlat(120,594) 

      varname='GLON'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
!	  GDLON(I,J)=(DUMMY2(I,J)*(180./acos(-1.)))
          GDLON(I,J)=DUMMY2(I,J)*RTD
          if(i.eq.409.and.j.eq.835)print*,'GDLAT GDLON in INITPOST=',    &
      	     i,j,GDLAT(I,J),GDLON(I,J)
        enddo
      enddo
      if(jsta.le.594.and.jend.ge.594)print*,'gdlon(120,594)= ',         &
        gdlon(120,594)

      varname='MXSNAL'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          MXSNAL(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='RADOT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          RADOT(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='SIGT4'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
           SIGT4(I,J)=DUMMY2(I,J)
        enddo
      enddo
       
      varname='TGROUND'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
  
      do j = jsta_2l, jend_2u
        do i = 1, im
           TG(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='CWM'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
          cwm ( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do

      varname='RRW'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,     &
              IM,1,JM,1,IM,JS,JE,1)
  
      varname='F_ICE'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
                        IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
    
      do l = 1, lm
        do j = jsta_2l, jend_2u
          do i = 1, im
            F_ice( i, j, l ) = dum3d( i, j, l )
      if (I .eq. im/2 .and. J .eq. jend_2u/2) then
      write(0,*) 'I, J, L, F_ICE: ', I,J,L, F_ice(I,J,L)
      endif
          end do
        end do
      end do 
    
      varname='F_RAIN'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
                        IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
    
      do l = 1, lm
        do j = jsta_2l, jend_2u
          do i = 1, im
            F_rain( i, j, l ) = dum3d( i, j, l )
          end do
        end do
      end do  
    
      varname='F_RIMEF'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
                        IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
    
      do l = 1, lm
        do j = jsta_2l, jend_2u
          do i = 1, im
            F_RimeF( i, j, l ) = dum3d( i, j, l )
          end do
        end do
      end do  

      varname='CLDFRA'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
           CFR( i, j, l ) = dum3d( i, j, l )
        end do
       end do
      end do

      varname='SR'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
      write(6,*) 'maxval SR: ', maxval(DUMMY2)
  
      do j = jsta_2l, jend_2u
        do i = 1, im
           SR(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='CFRACH'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
      write(6,*) 'maxval CFRACH: ', maxval(DUMMY2)
  
      do j = jsta_2l, jend_2u
        do i = 1, im
           CFRACH(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='CFRACL'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
      write(6,*) 'maxval CFRACL: ', maxval(DUMMY2)

      do j = jsta_2l, jend_2u
        do i = 1, im
           CFRACL(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='CFRACM'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
           CFRACM(I,J)=DUMMY2(I,J)
        enddo
      enddo  
      write(6,*) 'maxval CFRACM: ', maxval(DUMMY2)

      varname='ISLOPE'
      write(6,*) 'call getVariableB for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
          ISLOPE(I,J)=IDUMMY(I,J)
        enddo
      enddo

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

      VarName='DZSOIL'
      call getVariableB(fileName,DateStr,DataHandle,VarName,SLDPTH2,    &
              1,1,1,NSOIL,1,1,1,NSOIL)

      IF(isf_PHYSICS /= 3)then 
       DUMCST=0.0
       DO N=1,NSOIL
          DUMCST=DUMCST+SLDPTH2(N)
       END DO
       IF(ABS(DUMCST-0.).GT.1.0E-2)THEN
         DO N=1,NSOIL
           SLDPTH(N)=SLDPTH2(N)
         END DO
       END IF
      ENDIF


!  varname='SOILTB'
!  write(6,*) 'call getVariableB for : ', VarName
!      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,
!     &  IM,1,JM,1,IM,JS,JE,1)


! or get SLDPTH from wrf output
      VarName='SLDPTH'
      call getVariableB(fileName,DateStr,DataHandle,VarName,SLDPTH2,  &
     & 1,1,1,NSOIL,1,1,1,NSOIL)

      IF(isf_SURFACE_PHYSICS == 3)then  ! RUC surface layer
! if SLDPTH in wrf output is non-zero, then use it
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

      ENDIF

      VarName='CMC'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
! flip soil layer again because wrf soil variable vertical indexing
! is the same with eta and vertical indexing was flipped for both
! atmospheric and soil layers within getVariable
          cmc ( i, j ) = dummy ( i, j)
        end do
      end do
      print*,'CMC at ',ii,jj,' = ',cmc(ii,jj)

      varname='GRNFLX'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
  		    GRNFLX(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='PCTSNO'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          PCTSNO(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='SOILTB'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          SOILTB(I,J)=DUMMY2(I,J)
        enddo
      enddo

      varname='VEGFRC'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

      do j = jsta_2l, jend_2u
        do i = 1, im
          VEGFRC(I,J)=DUMMY2(I,J)
        enddo
      enddo

      VarName='SH2O'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,NSOIL)

      do l = 1, nsoil
       do j = jsta_2l, jend_2u
        do i = 1, im
          sh2o ( i, j, l ) = dum3d ( i, j, nsoil-l+1)
        end do
       end do
      end do

      VarName='SMC'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,NSOIL)
      do l = 1, nsoil
       do j = jsta_2l, jend_2u
        do i = 1, im
! flip soil layer again because wrf soil variable vertical indexing
! is the same with eta and vertical indexing was flipped for both
! atmospheric and soil layers within getVariable
            smc ( i, j, l ) = dum3d ( i, j, nsoil-l+1)
        end do
       end do
      end do
      print*,'SMC at ',ii,jj,N,' = ',smc(ii,jj,1),smc(ii,jj,2)         &
         ,smc(ii,jj,3),smc(ii,jj,4)

      VarName='STC'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,NSOIL)
      do l = 1, nsoil
       do j = jsta_2l, jend_2u
        do i = 1, im
!           stc ( i, j, l ) = dum3d ( i, j, l )
! flip soil layer again because wrf soil variable vertical indexing
! is the same with eta and vertical indexing was flipped for both
! atmospheric and soil layers within getVariable
            stc ( i, j, l ) = dum3d ( i, j, nsoil-l+1)
        end do
       end do
      end do
      print*,'STC at ',ii,jj,N,' = ',stc(ii,jj,1),stc(ii,jj,2)         &
        ,stc(ii,jj,3),stc(ii,jj,4)

! to skip some variables
      DO i=1, 16
         READ( unit=DataHandle) hdrbuf
         code = hdrbuf(2)
         WRITE(6,*) 'skip variable: ', code
         IF(code .EQ. 530) READ( unit=DataHandle)
      END DO

      VarName='PINT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM+1)
      do l = 1, lm+1
       do j = jsta_2l, jend_2u
        do i = 1, im
           PINT ( i, j, l ) = dum3d ( i, j, l ) 
           ALPINT(I,J,L)=ALOG(PINT(I,J,L))     
           if (L .ge. 2) then
!            PMID ( i, j, l-1 ) = EXP((ALOG(PINT(I,J,L-1))+
!     &               ALOG(PINT(I,J,L)))*0.5)
             PMID ( i, j, l-1 ) = (PINT(I,J,L-1)+                       &
                     PINT(I,J,L))*0.5 ! representative of what model does
           endif
        end do
       end do
      end do

      do l = 1, lm
       do j = jsta, jend
        do i = 1, im-MOD(J,2) 
         IF(J .EQ. 1 .AND. I .LT. IM)THEN   !SOUTHERN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J,L)+PMID(I+1,J,L))
         ELSE IF(J.EQ.JM .AND. I.LT.IM)THEN   !NORTHERN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J,L)+PMID(I+1,J,L))
         ELSE IF(I .EQ. 1 .AND. MOD(J,2) .EQ. 0) THEN   !WESTERN EVEN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J-1,L)+PMID(I,J+1,L))
         ELSE IF(I .EQ. IM .AND. MOD(J,2) .EQ. 0                        &
      	 .AND. J .LT. JM) THEN   !EASTERN EVEN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J-1,L)+PMID(I,J+1,L))  
         ELSE IF (MOD(J,2) .LT. 1) THEN
           PMIDV(I,J,L)=0.25*(PMID(I,J,L)+PMID(I-1,J,L)                 &
             +PMID(I,J+1,L)+PMID(I,J-1,L))
         ELSE
           PMIDV(I,J,L)=0.25*(PMID(I,J,L)+PMID(I+1,J,L)                 &
             +PMID(I,J+1,L)+PMID(I,J-1,L))
         END IF  
        end do
       end do
      end do


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

! SECOND, INTEGRATE HEIGHT HYDROSTATICLY
      DO L=LM,1,-1
       do j = jsta_2l, jend_2u
        do i = 1, im
         FI(I,J,2)=HTM(I,J,L)*T(I,J,L)*(Q(I,J,L)*D608+1.0)*RD*          &   
                   (ALPINT(I,J,L+1)-ALPINT(I,J,L))+FI(I,J,1)
         ZINT(I,J,L)=FI(I,J,2)/G
         if(i.eq.ii.and.j.eq.jj)                                        &
           print*,'L,sample HTM,T,Q,ALPINT(L+1),ALPINT(l),ZINT= ',      &
           l,HTM(I,J,L),T(I,J,L),Q(I,J,L),ALPINT(I,J,L+1),              &
           ALPINT(I,J,L),ZINT(I,J,L)
         FI(I,J,1)=FI(I,J,2)
        ENDDO
       ENDDO
      END DO
      print*,'finish deriving geopotential in nmm'
!
      DO L=1,LM
       DO I=1,IM
        DO J=JS,JE
!         ZMID(I,J,L)=(ZINT(I,J,L+1)+ZINT(I,J,L))*0.5  ! ave of z   
         FACT=(ALOG(PMID(I,J,L))-ALOG(PINT(I,J,L)))/                   &
               (ALOG(PINT(I,J,L+1))-ALOG(PINT(I,J,L)))	 
         ZMID(I,J,L)=ZINT(I,J,L)+(ZINT(I,J,L+1)-ZINT(I,J,L))           &
             *FACT
        ENDDO
       ENDDO
      ENDDO

      VarName='W'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM+1)
!      do l = 1, lm+1
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            w ( i, j, l ) = dum3d ( i, j, l )
!        end do
!       end do
!      end do

!!! SEE IF NEEDED

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
!            wh ( i, j, l ) = 0.5*(w(i,j,l)+w(i,j,l+1))
!         if(lm.eq.1 .and. dum3d(i,j,l) .gt.1.0e-10)print*,
!     +   'nonzero w at 1st layer',i,j,dum3d(i,j,1)
!         if(dum3d( i, j, l+1 ) .ne. 0.)print*,'nonzero w in INITPOST'
         if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'WH= ',            &
               i,j,l,WH ( i, j, l )
         wh ( i, j, l ) = dum3d ( i, j, l+1 )
        end do
       end do
      end do

      print*,'W at ',ii,jj,ll,' = ',WH (ii,jj,ll)

      VarName='ACFRCV'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)

       do j = jsta_2l, jend_2u
        do i = 1, im
            ACFRCV ( i, j ) = dummy ( i, j )
        end do
       end do
        write(6,*) 'MAX ACFRCV: ', maxval(DUMMY)

      VarName='ACFRST'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)  

       do j = jsta_2l, jend_2u
        do i = 1, im
            ACFRST ( i, j ) = dummy ( i, j )
        end do
       end do
        write(6,*) 'max ACFRST ', maxval(DUMMY)

!insert-mp
      VarName='SSROFF'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SSROFF ( i, j ) = dummy ( i, j )
        end do
       end do

! reading UNDERGROUND RUNOFF
      VarName='BGROFF'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            BGROFF ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='RLWIN'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            RLWIN ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='RLWTOA'
        write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            RLWTOA ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='ALWIN'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ALWIN ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='ALWOUT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ALWOUT ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='ALWTOA'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ALWTOA( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='RSWIN'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            RSWIN( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='RSWINC'
        write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            RSWINC( i, j ) = dummy ( i, j )
        end do
       end do
       print*,'max RSWINC= ',maxval(RSWINC)

      VarName='RSWOUT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            RSWOUT( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='ASWIN'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ASWIN( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='ASWOUT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ASWOUT( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='ASWTOA'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ASWTOA( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='SFCSHX'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SFCSHX( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='SFCLHX'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SFCLHX( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='SUBSHX'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SUBSHX( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='SNOPCX'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SNOPCX( i, j ) = dummy ( i, j )
        end do
       end do
  
      VarName='SFCUVX'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SFCUVX( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='POTEVP'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            POTEVP( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='POTFLX'
        write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)

!mhu        varname='RLWTT'
!mhu       write(6,*) 'call getVariableB for : ', VarName
!mhu      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,
!mhu     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!mhu      do l = 1, lm
!mhu       do j = jsta_2l, jend_2u
!mhu        do i = 1, im
!mhu            rlwtt( i, j, l ) = dum3d ( i, j, l )
!mhu            if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample RLWTT= ',
!mhu     +         i,j,l,RLWTT( i, j, l )
!mhu        end do
!mhu       end do
!mhu      end do 

!mhu        varname='RSWTT'
!mhu        write(6,*) 'call getVariableB for : ', VarName
!mhu      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,
!mhu     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!mhu      do l = 1, lm
!mhu       do j = jsta_2l, jend_2u
!mhu        do i = 1, im
!mhu            rswtt ( i, j, l ) = dum3d ( i, j, l )
!mhu            if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample RSWTT= ',
!mhu     +         i,j,l,RSWTT( i, j, l )
!mhu        end do
!mhu       end do
!mhu      end do 

!mhu      do l = 1, lm
!mhu       do j = jsta_2l, jend_2u
!mhu        do i = 1, im
!mhu            ttnd ( i, j, l ) = rswtt(i,j,l) + rlwtt(i,j,l)
!mhu        end do
!mhu       end do
!mhu      end do
!mhu
!mhu        varname='TCUCN'
!mhu        write(6,*) 'call getVariableB for : ', VarName
!mhu      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,
!mhu     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!mhu      do l=1,lm
!mhu       do j = jsta_2l, jend_2u
!mhu        do i = 1, im
!mhu            tcucn ( i, j, l) = dum3d ( i, j, l )
!mhu            if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample TCUCN= ',
!mhu     +         i,j,l,TCUCN( i, j, l )
!mhu        end do
!mhu       end do
!mhu      end do
!mhu        varname='TRAIN'
!mhu        write(6,*) 'call getVariableB for : ', VarName
!mhu      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,
!mhu     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!mhu      do l=1,lm
!mhu       do j = jsta_2l, jend_2u
!mhu        do i = 1, im
!mhu           train ( i, j, l ) = dum3d ( i, j, l )
!mhu           if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample TRAIN= ',
!mhu   +         i,j,l,TRAIN( i, j, l )
!mhu        end do
!mhu       end do
!mhu      end do

      VarName='TLMIN'
        write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,     &
              IM,1,JM,1,IM,JS,JE,1)
      VarName='TLMAX'
        write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,     &
              IM,1,JM,1,IM,JS,JE,1)
      VarName='T02_MIN'
        write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,     &
              IM,1,JM,1,IM,JS,JE,1)
      VarName='T02_MAX'
        write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,     &
              IM,1,JM,1,IM,JS,JE,1)
      VarName='RH02_MIN'
        write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,     &
              IM,1,JM,1,IM,JS,JE,1)
      VarName='RH02_MAX'
        write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,     &
              IM,1,JM,1,IM,JS,JE,1)

      VarName='NCFRCV'
        write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,      &
        IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ncfrcv ( i, j ) = float(idummy ( i, j ))
        end do
       end do

      VarName='NCFRST'
        write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,      &
        IM,1,JM,1,IM,JS,JE,1)  
       do j = jsta_2l, jend_2u
        do i = 1, im
            ncfrst ( i, j ) = float(idummy ( i, j ))
        end do
       end do

! set default to not empty buket
      NSRFC=0
      NRDLW=0
      NRDSW=0
      NHEAT=0
      NCLOD=0
      NPREC=0

      VarName='NPHS0'
      write(6,*) 'call getVariableB for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NPHS,      &  
              1,1,1,1,1,1,1,1)
      write(6,*) 'NPHS= ', NPHS

      VarName='NCNVC0'
      write(6,*) 'call getVariableB for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,    &
              1,1,1,1,1,1,1,1)

      VarName='NPREC'
      write(6,*) 'call getVariableB for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NPREC,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'NPREC= ', NPREC

      VarName='NCLOD'
      write(6,*) 'call getVariableB for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NCLOD,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'NCLOD= ', NCLOD
      
      VarName='NHEAT'
      write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NHEAT,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'NHEAT= ', NHEAT      

      VarName='NRDLW'
      write(6,*) 'call getVariableB for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NRDLW,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'NRDLW= ', NRDLW

      VarName='NRDSW'
      write(6,*) 'call getVariableB for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NRDSW,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'NRDSW= ', NRDSW

      VarName='NSRFC'
      write(6,*) 'call getVariableB for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,NSRFC,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'NSRFC= ', NSRFC

      VarName='AVRAIN'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,AVRAIN,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'AVRAIN= ', AVRAIN

      VarName='AVCNVC'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,AVCNVC,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'AVCNVC= ', AVCNVC

! We are overwriting the previously read variable here - is that what
!  you want to do?  
      VarName='ACUTIM'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,AVCNVC,     &
              1,1,1,1,1,1,1,1)
      write(6,*) 'ACUTIM= ',AVCNVC 

      VarName='ARDLW'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,ARDLW,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'ARDLW= ', ARDLW

      VarName='ARDSW'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,ARDSW,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'ARDSW= ', ARDSW
  
      VarName='ASRFC'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,ASRFC,      &
              1,1,1,1,1,1,1,1)
      write(6,*) 'ASRFC= ', ASRFC	

      VarName='APHTIM'
        write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,APHTIM,      &
              1,1,1,1,1,1,1,1)
        write(6,*) 'APHTIM= ', APHTIM


!mp - end nmm-core specific vars

!!!!!!!!!!!!!!!!!
!  STOP
!!!!!!!!!!!!!!!!!
!
!

!
! reading TKE
!      VarName='TKE_PBL'
!      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
!       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!      do l = 1, lm
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            q2 ( i, j, l ) = dum3d ( i, j, l )
!        end do
!       end do
!      end do
!      print*,'TKE at ',ii,jj,ll,' = ',q2(ii,jj,ll)
!
      VarName='LANDMASK'
        write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)

! retrieve hydrometeo fields directly for non-Ferrier
      VarName='QVAPOR'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
              IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
!           q ( i, j, l ) = dum3d ( i, j, l )
!           if(l.eq.1)print*,'Debug: I,J,Q= ',i,j,q( i, j, l )
!CHC CONVERT MIXING RATIO TO SPECIFIC HUMIDITY
           q ( i, j, l ) = dum3d ( i, j, l )/(1.0+dum3d ( i, j, l ))
        end do
       end do
      end do
      print*,'finish reading specific humidity'
      if(jj.ge. jsta .and. jj.le.jend)print*,'sample Q= ',Q(ii,jj,ll)
      qqw=spval
      qqr=spval
      qqs=spval
      qqi=spval
      qqg=spval 
      cwm=spval
      
      VarName='QCLOUD'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
              IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
! partition cloud water and ice for WSM3 
          if(imp_physics.eq.3)then 
           if(t(i,j,l) .ge. TFRZ)then  
            qqw ( i, j, l ) = dum3d ( i, j, l )
           else
            qqi  ( i, j, l ) = dum3d ( i, j, l )
           end if
          else ! bug fix provided by J CASE
           qqw ( i, j, l ) = dum3d ( i, j, l )
          end if 
          cwm(i,j,l)=dum3d(i,j,l)        
        end do
       end do
      end do

      if(jj.ge. jsta .and. jj.le.jend)     &
        print*,'sample qqw= ',Qqw(ii,jj,ll)
       
      VarName='QRAIN'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
              IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
! partition rain and snow for WSM3   
         if(imp_physics .eq. 3)then
          if(t(i,j,l) .ge. TFRZ)then  
           qqr ( i, j, l ) = dum3d ( i, j, l )
          else
           qqs ( i, j, l ) = dum3d ( i, j, l )
          end if
         else ! bug fix provided by J CASE
          qqr ( i, j, l ) = dum3d ( i, j, l )  
         end if
         cwm(i,j,l)=cwm(i,j,l)+dum3d(i,j,l)
        end do
       end do
      end do
       
      if(jj.ge. jsta .and. jj.le.jend)print*,'sample qqr= ',Qqr(ii,jj,ll) 

!tms      VarName='QICE'
!tms      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
!tms              IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!tms      do l = 1, lm
!tms       do j = jsta_2l, jend_2u
!tms        do i = 1, im
!tms          qqi ( i, j, l ) = dum3d ( i, j, l )
!tms!          qqi ( i, j, l ) = 0
!tms           cwm(i,j,l)=cwm(i,j,l)+dum3d(i,j,l)
!tms        end do
!tms       end do
!tms      end do
!tms
!tms      if(jj.ge. jsta .and. jj.le.jend)print*,'sample qqi= ',Qqi(ii,jj,ll)
      
      VarName='QSNOW'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
              IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
          qqs ( i, j, l ) = dum3d ( i, j, l )
          cwm(i,j,l)=cwm(i,j,l)+dum3d(i,j,l)
        end do
       end do
      end do

      if(jj.ge. jsta .and. jj.le.jend)print*,'sample qqs= ',Qqs(ii,jj,ll)
       
!tms - QGRAUP not found in data remove it here
!tms      VarName='QGRAUP'
!tms      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,
!tms     &    IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!tms      do l = 1, lm
!tms       do j = jsta_2l, jend_2u
!tms        do i = 1, im
!tms          qqg ( i, j, l ) = dum3d ( i, j, l )
!tms!          qqg ( i, j, l ) = 0
!tms          cwm(i,j,l)=cwm(i,j,l)+dum3d(i,j,l)
!tms        end do
!tms       end do
!tms      end do

!tms      if(jj.ge. jsta .and. jj.le.jend)print*,'sample qqg= '
!tms     +,Qqg(ii,jj,ll)
          print*, 'QGRAUP not in data'

! end of retrieving hydrometeo 

      VarName='SMOIS'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)

! reading sfc pressure
      VarName='PSFC'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
        do j = jsta_2l, jend_2u
        do i = 1, im
!          IF(abs(PT+PDTOP+PD(I,J)-DUMMY(I,J)).GT.100.)
!     &    print*,'incon pt',PT+PDTOP+PD(I,J),DUMMY(I,J)
        end do
       end do     
!
! reading 2m theta
      VarName='TH2'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)

      
!! already in TSHLTR
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            TSHLTR ( i, j ) = dummy ( i, j )
!        end do
!       end do
!
! reading 10 m wind
      VarName='U10'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            U10 ( i, j ) = dummy ( i, j )
        end do
       end do
       print*,'U10 at ',ii,jj,' = ',U10(ii,jj)

      VarName='V10'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            V10 ( i, j ) = dummy ( i, j )
        end do
       end do
       print*,'V10 at ',ii,jj,' = ',V10(ii,jj)
!
! to skip some variables
      DO i=1, 1
         READ( unit=DataHandle) hdrbuf
         code = hdrbuf(2)
         WRITE(6,*) 'skip variable: ', code
         IF(code .EQ. 530) READ( unit=DataHandle)
      END DO

!
! reading SMSTAV
      VarName='SMSTAV'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SMSTAV ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='SMSTOT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SMSTOT ( i, j ) = dummy ( i, j )
        end do
       end do
!
! reading SURFACE RUNOFF 
      VarName='SFROFF'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
!
! reading UNDERGROUND RUNOFF
      VarName='UDROFF'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)

! reading VEGETATION TYPE 
      VarName='IVGTYP'
      write(6,*) 'call getVariableB for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY       &
        ,IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
            IVGTYP( i, j ) = idummy ( i, j )
        end do
       end do
      print*,'IVGTYP at ',ii,jj,' = ',IDUMMY(ii,jj) 

      VarName='ISLTYP' 
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY       &
        ,IM,1,JM,1,IM,JS,JE,1)
      do j = jsta_2l, jend_2u
        do i = 1, im
            ISLTYP( i, j ) = idummy ( i, j )
        end do
       end do
      print*,'ISLTYP at ',ii,jj,' = ',IDUMMY(ii,jj)

      VarName='VEGFRA'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
      print*,'VEGFRA at ',ii,jj,' = ',DUMMY(ii,jj) 

      VarName='SFCEVP'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SFCEVP( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='GRDFLX'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)

      VarName='SFCEXC'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SFCEXC( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='ACSNOW'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ACSNOW ( i, j ) = dummy ( i, j )
        end do
       end do
       
      VarName='ACSNOM'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            ACSNOM ( i, j ) = dummy ( i, j )
        end do
       end do

      VarName='SNOW'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)

      VarName='CANWAT'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            CMC ( i, j ) = dummy ( i, j )
!        end do
!       end do

      VarName='SST'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            SST ( i, j ) = dummy ( i, j )
        end do
       end do
      print*,'SST at ',ii,jj,' = ',sst(ii,jj)      

      VarName='WEASD'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
!      do j = jsta_2l, jend_2u
!        do i = 1, im
!            SNO ( i, j ) = dummy ( i, j )
!        end do
!       end do

! reading TKE
!mhu      VarName='TKE_PBL'
!mhu      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,
!mhu     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!     do l = 1, lm
!      do j = jsta_2l, jend_2u 
!       do i = 1, im
!           q2 ( i, j, l ) = dum3d ( i, j, l )
!       end do
!      end do
!     end do
!     print*,'TKE at ',ii,jj,ll,' = ',q2(ii,jj,ll)

!mhu      VarName='EL_PBL'
!mhu      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,
!mhu   &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!mhu      do l = 1, lm
!mhu       do j = jsta_2l, jend_2u
!mhu        do i = 1, im
!mhu            EL_PBL( i, j, l ) = dum3d ( i, j, l )
!mhu            if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample EL= ',
!mhu     +         i,j,l,EL_PBL( i, j, l )
!mhu        end do
!mhu       end do
!mhu      end do

!mhu      VarName='EXCH_H'
!mhu      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,
!mhu     &  IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!mhu      do l = 1, lm
!mhu       do j = jsta_2l, jend_2u
!mhu        do i = 1, im
!mhu            EXCH_H( i, j, l ) = dum3d ( i, j, l )
!mhu        end do
!mhu       end do
!mhu      end do

! to skip some variables (ZNT & NOAHRES)
      DO i=1, 2
         READ( unit=DataHandle) hdrbuf
         code = hdrbuf(2)
         WRITE(6,*) 'skip variable: ', code
         IF(code .EQ. 530) READ( unit=DataHandle)
      END DO


      VarName='THZ0'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            THZ0 ( i, j ) = dummy ( i, j )
        end do
       end do
      print*,'THZ0 at ',ii,jj,' = ',THZ0(ii,jj)

      VarName='QZ0'
      write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            QZ0 ( i, j ) = dummy ( i, j )
        end do
       end do
      print*,'QZ0 at ',ii,jj,' = ',QZ0(ii,jj)

      VarName='UZ0'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            UZ0 ( i, j ) = dummy ( i, j )
        end do
       end do
      print*,'UZ0 at ',ii,jj,' = ',UZ0(ii,jj)

      VarName='VZ0'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            VZ0 ( i, j ) = dummy ( i, j )
        end do
       end do
      print*,'VZ0 at ',ii,jj,' = ',VZ0(ii,jj)

      VarName='QSFC'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)

! retrieve htop and hbot
      VarName='HTOP'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1) 
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            HTOP ( i, j ) = float(LM)-dummy(i,j)+1.0 
!        end do
!       end do
!       print*,'maxval HTOP: ', maxval(HTOP)

      VarName='HBOT'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            HBOT ( i, j ) = float(LM)-dummy(i,j)+1.0
!            if(HBOT(i,j).gt.0.1 .and. (HTOP(i,j).gt.hbot(i,j)))
!     &      print*,'htop gt hbot ',i,j,htop(i,j),hbot(i,j)
!        end do
!       end do
       print*,'maxval HBOT: ', maxval(HBOT)

      VarName='HTOPD'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            HTOPD( i, j ) = float(LM)-dummy(i,j)+1.0
        end do
       end do
       print*,'maxval HTOPD: ', maxval(HTOPD)

      VarName='HBOTD'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            HBOTD( i, j ) = float(LM)-dummy(i,j)+1.0
        end do
       end do
       print*,'maxval HBOTD: ', maxval(HBOTD)

      VarName='HTOPS'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            HTOPS( i, j ) = float(LM)-dummy(i,j)+1.0
        end do
       end do
       print*,'maxval HTOPS: ', maxval(HTOPS)
                                                                                 
      VarName='HBOTS'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            HBOTS( i, j ) = float(LM)-dummy(i,j)+1.0
        end do
       end do
       print*,'maxval HBOTS: ', maxval(HBOTS)

!      VarName='SR'
!      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
!       IM,1,JM,1,IM,JS,JE,1)
!       do j = jsta_2l, jend_2u
!        do i = 1, im
!            SR ( i, j ) = dummy(i,j)
!        end do
!       end do
!       print*,'maxval SR: ', maxval(SR)

      VarName='CUPPT'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            CUPPT ( i, j ) = dummy ( i, j )
        end do
       end do
       print*,'maxval CUPPT: ', maxval(DUMMY)

      VarName='CPRATE'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)
       do j = jsta_2l, jend_2u
        do i = 1, im
            CPRATE ( i, j ) = dummy ( i, j )
        end do
       end do
       print*,'maxval CPRATE: ', maxval(DUMMY)

      VarName='SNOWH'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
       IM,1,JM,1,IM,JS,JE,1)

!!!! DONE GETTING

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            IF(ABS(T(I,J,L)).GT.1.0E-3)                                 &  
              OMGA(I,J,L) = -WH(I,J,L)*PMID(I,J,L)*G/                   &
                       (RD*T(I,J,L)*(1.+D608*Q(I,J,L)))

        end do
       end do
      end do

    
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
       write(6,*) 'laststart,latlast B calling bcast=',latstart,latlast
       call mpi_bcast(latstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(latlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(cenlat,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*) 'laststart,latlast A calling bcast= ',latstart,latlast
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
       write(6,*)'lonstart,lonlast B calling bcast=',lonstart,lonlast
       call mpi_bcast(lonstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(lonlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(cenlon,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*)'lonstart,lonlast A calling bcast= ',lonstart,lonlast
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
!  status=nf_close(ncid)

!  dxval=30000.
!   dyval=30000.
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
      call ext_int_ioclose ( DataHandle, Status )

! generate look up table for lifted parcel calculations

      THL=210.
      PLQ=70000.

      CALL TABLE(PTBL,TTBL,PT,                                         &  
                RDQ,RDTH,RDP,RDTHE,PL,THL,QS0,SQS,STHE,THE0)

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
!
!  Comments and code provided for use with copygb - R.Rozumalski - NWS
!
      IF (me.eq.0) THEN

         inav=10

         TRUELAT1 = CENLAT
         TRUELAT2 = CENLAT

         IFDX = NINT (dxval*107.)
         IFDY = NINT (dyval*110.)

         open(inav,file='copygb_gridnav.txt',form='formatted',          &
                   status='unknown')

         print *, ' MAPTYPE  :',maptype
         print *, ' IM       :',IM*2-1
         print *, ' JM       :',JM
         print *, ' LATSTART :',LATSTART
         print *, ' LONSTART :',LONSTART
         print *, ' CENLAT   :',CENLAT
         print *, ' CENLON   :',CENLON
         print *, ' TRUELAT2 :',TRUELAT2
         print *, ' TRUELAT1 :',TRUELAT1
         print *, ' DX       :',IFDX*0.001
         print *, ' DY       :',IFDY*0.001

         IF(MAPTYPE.EQ.0 .OR. MAPTYPE.EQ.203)THEN  !A STAGGERED E-GRID

            IMM = 2*IM-1
            IDXAVE = ( IFDY + IFDX ) * 0.5

            ! If the Center Latitude of the domain is located within 15 degrees
            ! of the equator then use a a regular Lat/Lon navigation for the
            ! remapped grid in copygb; otherwise, use a Lambert conformal.  Make
            ! sure to specify the correct pole for the S. Hemisphere (LCC).
            !
            IF ( abs(CENLAT).GT.15000) THEN
               write(6,*)'  Copygb LCC Navigation Information'
               IF (CENLAT .GT.0) THEN ! Northern Hemisphere
                  write(6,1000)    IMM,JM,LATSTART,LONSTART,CENLON,     &
                                   IFDX,IFDY,CENLAT,CENLAT
                  write(inav,1000) IMM,JM,LATSTART,LONSTART,CENLON,     &
                                   IFDX,IFDY,CENLAT,CENLAT
               ELSE  ! Southern Hemisphere
                  write(6,1001)    IMM,JM,LATSTART,LONSTART,CENLON,     &
                                   IFDX,IFDY,CENLAT,CENLAT
                  write(inav,1001) IMM,JM,LATSTART,LONSTART,CENLON,     &
                                   IFDX,IFDY,CENLAT,CENLAT
               END IF
            ELSE
               dlat = (latnm-latsm)/(JM-1)
               nlat = INT (dlat)

               if (lonem .lt. 0) lonem = 360000. + lonem
               if (lonwm .lt. 0) lonwm = 360000. + lonwm

               dlon = lonem-lonwm
               if (dlon .lt. 0.) dlon = dlon + 360000.
               dlon = (dlon)/(IMM-1)
               nlon = INT (dlon)

               write(6,*)'  Copygb Lat/Lon Navigation Information'
               write(6,2000)    IMM,JM,latsm,lonwm,latnm,lonem,nlon,nlat
               write(inav,2000) IMM,JM,latsm,lonwm,latnm,lonem,nlon,nlat
            ENDIF

 1000     format('255 3 ',2(I3,x),I6,x,I7,x,'8 ',I7,x,2(I6,x),'0 64',   &
                 2(x,I6))
 1001     format('255 3 ',2(I3,x),I6,x,I7,x,'8 ',I7,x,2(I6,x),'128 64', &
                 2(x,I6),' -90000 0')
 2000     format('255 0 ',2(I3,x),2(I7,x),'8 ',2(I7,x),2(I7,x),'64') 

         END IF ! IF (MAPTYPE  ...

      END IF ! IF (me.eq.0)

!  End of R.Rozumalski modifications


!CHC WRITE IGDS OUT FOR WEIGHTMAKER TO READ IN AS KGDSIN
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
      IFDX = NINT (dxval*107.)
      IFDY = NINT (dyval*110.)
      inav = 10
      open(inav,file='copygb_gridnav.txt',form='formatted',        &
               status='unknown')
      IMM = 2*IM-1
      write(10,1005) IMM,JM,LATSTART,LONSTART,CENLON,              &
                     IFDX,IFDY,CENLAT,CENLAT

1005  format('255 3 ',2(I3,x),I6,x,I7,x,'8 ',I7,x,2(I6,x),         &
             '0 64', 2(x,I6))
      close (inav)

! following for hurricane wrf post
      open(inav,file='copygb_hwrf.txt',form='formatted',           &
                status='unknown')
      LATEND=LATSTART+(JM-1)*dyval
      LONEND=LONSTART+(IMM-1)*dxval
      write(10,1010) IMM,JM,LATSTART,LONSTART,LATEND,LONEND,       &
                     dxval,dyval

1010 format('255 0 ',2(I3,x),I6,x,I7,x,'136 ',I6,x,I7,x,           &
             2(I6,x),'64')

! jkw use '64' for LLC rather than '0' which would be ULC
! jkw     +           2(I6,x),'0')

      close (inav)
!
      RETURN
      END
