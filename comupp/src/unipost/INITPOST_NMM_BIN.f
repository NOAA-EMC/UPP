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
!
!     INCLUDE/SET PARAMETERS.
!     
      INCLUDE "mpif.h"
! This version of INITPOST shows how to initialize, open, read from, and
! close a NetCDF dataset. In order to change it to read an internal (binary)
! dataset, do a global replacement of _ncd_ with _int_. 

      character(len=31) :: VarName
      integer :: Status
      character startdate*19,SysDepInfo*80
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
!jw
      integer ii,jj,js,je,jev,iyear,imn,iday,itmp,ioutcount,istatus,    &
              nsrfc,nrdlw,nrdsw,nheat,nclod
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
! Getting start time
      call ext_int_get_dom_ti_char(DataHandle,'START_DATE',startdate,  &
        status )
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
       print*,'new start yr mo day hr min =',sdat(3),sdat(1)            &  
             ,sdat(2),ihrst,imin
      END IF
      
      
!  OK, since all of the variables are dimensioned/allocated to be
!  the same size, this means we have to be careful int getVariable
!  to not try to get too much data.  For example, 
!  DUM3D is dimensioned IM+1,JM+1,LM+1 but there might actually
!  only be im,jm,lm points of data available for a particular variable.  
! get metadata
! NMM does not output mp_physics yet so hard-wire it to Ferrier scheme
        call ext_int_get_dom_ti_integer(DataHandle,'MP_PHYSICS'         &
          ,itmp,1,ioutcount,istatus)
        imp_physics=itmp
        print*,'MP_PHYSICS= ',imp_physics

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
!        tmp=50.0
        cenlat=nint(1000.*tmp)
        write(6,*) 'cenlat= ', cenlat
        call ext_int_get_dom_ti_real(DataHandle,'CEN_LON',tmp           &
          ,1,ioutcount,istatus)
!        tmp=-111.0
        cenlon=nint(1000.*tmp)
        write(6,*) 'cenlon= ', cenlon
        call ext_int_get_dom_ti_real(DataHandle,'TRUELAT1',tmp          &
          ,1,ioutcount,istatus)
        truelat1=nint(1000.*tmp)
        write(6,*) 'truelat1= ', truelat1
        call ext_int_get_dom_ti_real(DataHandle,'TRUELAT2',tmp          &
          ,1,ioutcount,istatus)
        truelat2=nint(1000.*tmp)
        write(6,*) 'truelat2= ', truelat2
        call ext_int_get_dom_ti_integer(DataHandle,'MAP_PROJ',itmp      &
          ,1,ioutcount,istatus)
        maptype=itmp
!       maptype=203
        write(6,*) 'maptype is ', maptype
        gridtype = 'E'

! get 3-D variables
      print*,'im,jm,lm= ',im,jm,lm
!
      VarName='LU_INDEX'
	write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,     &
        IM,1,JM,1,IM,JS,JE,1)
     
      VarName='LMH'
	write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,     &
        IM,1,JM,1,IM,JS,JE,1)

      VarName='LMV'
	write(6,*) 'call getIVariable for : ', VarName
      call getIVariable(fileName,DateStr,DataHandle,VarName,IDUMMY,     &
        IM,1,JM,1,IM,JS,JE,1)

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
        if (j.eq.jm/2 .and. mod(i,10).eq.0)                               &
           print*,'sample SM= ',i,j,sm(i,j)
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

      VarName='HTM'
	write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      VarName='VTM'
	write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      VarName='PD'
	write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
        IM,1,JM,1,IM,JS,JE,1)

       do j = jsta_2l, jend_2u
        do i = 1, im
	PD(I,J)=DUMMY2(I,J)
       enddo
       enddo      

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
	    if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample T= ',    &
               i,j,l,T ( i, j, l )
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
!     &	    print*,'1st level Q < 1.0e-7 in INTPOST'
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

       varname='ETA1'
	write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,ETA1,       &   
        LM,1,1,1,LM,1,1,1)

	varname='ETA2'
	write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,ETA2,       &
        LM,1,1,1,LM,1,1,1)

        open(75,file='ETAPROFILE.txt',form='formatted',                 &
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

	varname='USTAR'
	write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
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

	varname='PREC' ! instantaneous precip rate?
	write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
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
	 if(i.eq.409.and.j.eq.835)print*,'2m T and P in INITPOST='       &
      	     ,i,j,TSHLTR(i,j),PSHLTR(I,J)
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
!	write(6,*) 'call getVariableB for : ', VarName
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

	varname='GLAT'
	write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY2,      &
       IM,1,JM,1,IM,JS,JE,1)

       do j = jsta_2l, jend_2u
        do i = 1, im
	 F(I,J)=1.454441e-4*sin(DUMMY2(I,J))   ! 2*omeg*sin(phi)
!	  GDLAT(I,J)=DUMMY2(I,J)*(180./acos(-1.))
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
	 if(i.eq.409.and.j.eq.835)print*,'GDLAT GDLON in INITPOST='      &
      	     ,i,j,GDLAT(I,J),GDLON(I,J)
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

	varname='F_ICE'
	write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            F_ice( i, j, l ) = dum3d( i, j, l )
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

!	varname='SOILTB'
!	write(6,*) 'call getVariableB for : ', VarName
!      call getVariableB(fileName,DateStr,DataHandle,VarName,DUMMY,      &
!       IM,1,JM,1,IM,JS,JE,1)

! either assign SLDPTH to be the same as eta (which is original
! setup in WRF LSM) or extract thickness of soil layers from wrf
! output

! assign SLDPTH to be the same as eta

         SLDPTH(1)=0.10
         SLDPTH(2)=0.3
         SLDPTH(3)=0.6
         SLDPTH(4)=1.0

! or get SLDPTH from wrf output
      VarName='SLDPTH'
      call getVariableB(fileName,DateStr,DataHandle,VarName,SLDPTH2,  &
     & 1,1,1,NSOIL,1,1,1,NSOIL)
! if SLDPTH in wrf output is non-zero, then use it
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
!            stc ( i, j, l ) = dum3d ( i, j, l )
! flip soil layer again because wrf soil variable vertical indexing
! is the same with eta and vertical indexing was flipped for both
! atmospheric and soil layers within getVariable
            stc ( i, j, l ) = dum3d ( i, j, nsoil-l+1)
        end do
       end do
      end do
      print*,'STC at ',ii,jj,N,' = ',stc(ii,jj,1),stc(ii,jj,2)         &
        ,stc(ii,jj,3),stc(ii,jj,4)

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
         PMID ( i, j, l-1 ) = (PINT(I,J,L-1)+                           &  
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
        print*,'L,sample HTM,T,Q,ALPINT(L+1),ALPINT(l),ZINT= '          &
          ,l,HTM(I,J,L),T(I,J,L),Q(I,J,L),ALPINT(I,J,L+1),              &
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
!     +	 'nonzero w at 1st layer',i,j,dum3d(i,j,1)
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

        varname='RLWTT'
       write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            rlwtt( i, j, l ) = dum3d ( i, j, l )
            if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample RLWTT= ',&
               i,j,l,RLWTT( i, j, l )
        end do
       end do
      end do 

        varname='RSWTT'
        write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            rswtt ( i, j, l ) = dum3d ( i, j, l )
            if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample RSWTT= ',&
               i,j,l,RSWTT( i, j, l )
        end do
       end do
      end do 

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            ttnd ( i, j, l ) = rswtt(i,j,l) + rlwtt(i,j,l)
        end do
       end do
      end do

        varname='TCUCN'
        write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l=1,lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            tcucn ( i, j, l) = dum3d ( i, j, l )
            if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample TCUCN= ',&
               i,j,l,TCUCN( i, j, l )
        end do
       end do
      end do
        varname='TRAIN'
        write(6,*) 'call getVariableB for : ', VarName
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l=1,lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            train ( i, j, l ) = dum3d ( i, j, l )
            if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample TRAIN= ',&
               i,j,l,TRAIN( i, j, l )
        end do
       end do
      end do

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
!	STOP
!!!!!!!!!!!!!!!!!
!
!

!
! reading TKE
!      VarName='TKE_MYJ'
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
!     &	  print*,'incon pt',PT+PDTOP+PD(I,J),DUMMY(I,J)
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
      VarName='TKE_MYJ'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
!     do l = 1, lm
!      do j = jsta_2l, jend_2u 
!       do i = 1, im
!           q2 ( i, j, l ) = dum3d ( i, j, l )
!       end do
!      end do
!     end do
!     print*,'TKE at ',ii,jj,ll,' = ',q2(ii,jj,ll)

      VarName='EL_MYJ'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            EL_PBL( i, j, l ) = dum3d ( i, j, l )
            if(i.eq.im/2.and.j.eq.(jsta+jend)/2)print*,'sample EL= ',  &
               i,j,l,EL_PBL( i, j, l )
        end do
       end do
      end do

      VarName='EXCH_H'
      call getVariableB(fileName,DateStr,DataHandle,VarName,DUM3D,      &
       IM+1,1,JM+1,LM+1,IM,JS,JE,LM)
      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            EXCH_H( i, j, l ) = dum3d ( i, j, l )
        end do
       end do
      end do

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
       end if
       write(6,*) 'laststart,latlast B calling bcast=',latstart,latlast
       call mpi_bcast(latstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(latlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*) 'laststart,latlast A calling bcast=',latstart,latlast
       call collect_loc(gdlon,dummy)
       if(me.eq.0)then
        lonstart=nint(dummy(1,1)*1000.)
        lonlast=nint(dummy(im,jm)*1000.)
       end if
       write(6,*)'lonstart,lonlast B calling bcast=',lonstart,lonlast
       call mpi_bcast(lonstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(lonlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*)'lonstart,lonlast A calling bcast=',lonstart,lonlast
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
