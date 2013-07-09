       SUBROUTINE INITPOST_GFS_NEMS(NREC,iostatusFlux,iostatusD3D,   &
                                   iostatusAER,nfile,ffile,rfile)
!       SUBROUTINE INITPOST_GFS_NEMS(NREC,iostatusFlux,iostatusD3D,nfile,ffile)

!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    INITPOST    INITIALIZE POST FOR RUN
!   PRGRMMR: Hui-Ya Chuang    DATE: 2007-03-01
!     
! ABSTRACT:  THIS ROUTINE INITIALIZES CONSTANTS AND
!   VARIABLES AT THE START OF AN ETA MODEL OR POST 
!   PROCESSOR RUN.
!
! REVISION HISTORY
!   2011-02-07 Jun Wang    add grib2 option
!   2011-12-14 Sarah Lu    add aer option
!   2012-01-07 Sarah Lu    compute air density
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
      use vrbls4d, only: dust
      use vrbls3d, only: t, q, uh, vh, pmid, pint, alpint, dpres, zint, zmid, o3,&
              qqr, qqs, cwm, qqi, qqw, omga, rhomid, q2, cfr, rlwtt, rswtt, tcucn,&
              tcucns, train, el_pbl, exch_h, vdifftt, vdiffmois, dconvmois, nradtt,&
              o3vdiff, o3prod, o3tndy, mwpv, unknown, vdiffzacce, zgdrag,cnvctummixing,&
              vdiffmacce, mgdrag, cnvctvmmixing, ncnvctcfrac, cnvctumflx, cnvctdmflx,&
              cnvctzgdrag, sconvmois, cnvctmgdrag, cnvctdetmflx, duwt, duem, dusd, dudp
      use vrbls2d, only: f, pd, fis, pblh, ustar, z0, ths, qs, twbs, qwbs, avgcprate,&
              cprate, avgprec, prec, lspa, sno, si, cldefi, th10, q10, tshltr, pshltr,&
              tshltr, albase, avgalbedo, avgtcdc, czen, czmean, mxsnal, radot, sigt4,&
              cfrach, cfracl, cfracm, avgcfrach, qshltr, avgcfracl, avgcfracm, cnvcfr,&
              islope, cmc, grnflx, vegfrc, acfrcv, ncfrcv, acfrst, ncfrst, ssroff,&
              bgroff, rlwin, rlwtoa, cldwork, alwin, alwout, alwtoa, rswin, rswinc,&
              rswout, aswin, auvbin, auvbinc, aswout, aswtoa, sfcshx, sfclhx, subshx,&
              snopcx, sfcux, sfcvx, sfcuvx, gtaux, gtauy, potevp, u10, v10, smstav,&
              smstot, ivgtyp, isltyp, sfcevp, sfcexc, acsnow, acsnom, sst, thz0, qz0,&
              uz0, vz0, ptop, htop, pbot, hbot, ptopl, pbotl, ttopl, ptopm, pbotm, ttopm,&
              ptoph, pboth, pblcfr, ttoph, runoff, maxtshltr, mintshltr, maxrhshltr,&
              minrhshltr, dzice, smcwlt, suntime, fieldcapa, htopd, hbotd, htops, hbots,&
              cuppt, dusmass, ducmass, dusmass25, ducmass25
      use soil, only: sldpth, sh2o, smc, stc
      use masks, only: lmv, lmh, htm, vtm, gdlat, gdlon, dx, dy, hbm2, sm, sice
      use kinds, only: i_llong
      use nemsio_module, only: nemsio_gfile, nemsio_getfilehead, nemsio_getheadvar, nemsio_close
      use physcons, only: con_g, con_fvirt, con_rd, con_eps, con_epsm1
      use params_mod, only: erad, dtr, tfrz, h1, d608, rd, p1000, capa
      use lookup_mod, only: thl, plq, ptbl, ttbl, rdq, rdth, rdp, rdthe, pl, qs0, sqs, sthe,&
              ttblq, rdpq, rdtheq, stheq, the0q, the0
      use ctlblk_mod, only: me, mpi_comm_comp, icnt, idsp, jsta, jend, ihrst, idat, sdat, ifhr,&
              ifmin, filename, tprec, tclod, trdlw, trdsw, tsrfc, tmaxmin, td3d, restrt, sdat,&
              jend_m, imin, imp_physics, dt, spval, pdtop, pt, qmin, nbin_du, nphs, dtq2, ardlw,&
              ardsw, asrfc, avrain, avcnvc, theat, gdsdegr, spl, lsm, alsl, im, jm, im_jm, lm,&
              jsta_2l, jend_2u, nsoil, lp1, icu_physics, ivegsrc, novegtype
      use gridspec_mod, only: maptype, gridtype, latstart, latlast, lonstart, lonlast, cenlon,&
              dxval, dyval, truelat2, truelat1, psmapf, cenlat
      use rqstfld_mod, only: igds, avbl, iq, is
      use wrf_io_flags_mod, only:
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      implicit none
!
      type(nemsio_gfile),intent(inout) :: nfile,ffile,rfile
!
!     INCLUDE/SET PARAMETERS.
!     
      INCLUDE "mpif.h"
      integer,parameter:: MAXPTS=1000000 ! max im*jm points
!
!     real,parameter:: con_g       =9.80665e+0! gravity
!     real,parameter:: con_rv      =4.6150e+2 ! gas constant H2O
!     real,parameter:: con_rd      =2.8705e+2 ! gas constant air
!     real,parameter:: con_fvirt   =con_rv/con_rd-1.
!     real,parameter:: con_eps     =con_rd/con_rv
!     real,parameter:: con_epsm1   =con_rd/con_rv-1
!
!      real,external::FPVSNEW
! This version of INITPOST shows how to initialize, open, read from, and
! close a NetCDF dataset. In order to change it to read an internal (binary)
! dataset, do a global replacement of _ncd_ with _int_. 

      integer,intent(in) :: NREC,iostatusFlux,iostatusD3D,iostatusAER
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
      LOGICAL RUNB,SINGLRST,SUBPOST,NEST,HYDRO
      LOGICAL IOOMG,IOALL
      logical, parameter :: debugprint = .true.
      CHARACTER*32 LABEL
      CHARACTER*40 CONTRL,FILALL,FILMST,FILTMP,FILTKE,FILUNV            &  
         , FILCLD,FILRAD,FILSFC
      CHARACTER*4 RESTHR
      CHARACTER FNAME*80,ENVAR*50
      INTEGER IDATE(8),JDATE(8)
      INTEGER JPDS(200),JGDS(200),KPDS(200),KGDS(200)
      LOGICAL*1 LB(IM,JM)
      INTEGER IRET
      REAL BUFF(IM_JM)
!     
!     INCLUDE COMMON BLOCKS.
!
!     DECLARE VARIABLES.
!     
!      REAL fhour
      integer nfhour ! forecast hour from nems io file
      REAL RINC(5)
      REAL ETA1(LM), ETA2(LM)
      REAL DUM1D (LM+1)
      REAL DUMMY ( IM, JM )
      REAL DUMMY2 ( IM, JM )
      REAL FI(IM,JM,2)
      INTEGER IDUMMY ( IM, JM )
!jw
      integer ii,jj,js,je,iyear,imn,iday,itmp,ioutcount,istatus, &
              I,J,L,ll,k,kf,irtn,igdout,n,Index,nframe, &
	      impf,jmpf,nframed2,iunitd3d
      real TSTART,TLMH,TSPH,ES, FACT,soilayert,soilayerb,zhour,dum
      real, external :: fpvsnew

      character*8,allocatable:: recname(:)
      character*16,allocatable  :: reclevtyp(:)
      integer,allocatable:: reclev(:)
      real, allocatable:: glat1d(:),glon1d(:),qstl(:)
      integer ierr,idum
   
      integer ibuf(im,jsta_2l:jend_2u)
      real buf(im,jsta_2l:jend_2u),bufsoil(im,nsoil,jsta_2l:jend_2u)   &
        ,buf3d(im,jsta_2l:jend_2u,lm),buf3d2(im,lp1,jsta_2l:jend_2u)
      real LAT
!      REAL,  PARAMETER    :: QMIN = 1.E-15
      real                :: TV

!
!      DATA BLANK/'    '/
!
!***********************************************************************
!     START INIT HERE.
!
      WRITE(6,*)'INITPOST:  ENTER INITPOST_GFS_NEMS'
      WRITE(6,*)'me=',me,'LMV=',size(LMV,1),size(LMV,2),'LMH=', &
           size(LMH,1),size(LMH,2),'jsta_2l=',jsta_2l,'jend_2u=', &
          jend_2u,'im=',im
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
! get start date
      if (me == 0)then
       print*,'nrec=',nrec
       allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
       allocate(glat1d(im*jm),glon1d(im*jm))
       call nemsio_getfilehead(nfile,iret=iret                           &  
         ,idate=idate(1:7),nfhour=nfhour,recname=recname                  &
         ,reclevtyp=reclevtyp,reclev=reclev,lat=glat1d               &
         ,lon=glon1d,nframe=nframe)
       if(iret/=0)print*,'error getting idate,nfhour'
        print *,'latstar1=',glat1d(1),glat1d(im*jm)
!       print *,'printing an inventory of GFS nemsio file'
!       do i=1,nrec
!        print *,'recname=',(trim(recname(i)))
!        print *,'reclevtyp=',(trim(reclevtyp(i)))
!        print *,'reclev=',(reclev(i))
!       end do
!       deallocate (recname,reclevtyp,reclev)
       
!       call nemsio_getfilehead(ffile,nrec=idum)
!       print*,'nrec for flux file = ',idum
!       allocate(recname(idum),reclevtyp(idum),reclev(idum))
!       call nemsio_getfilehead(ffile,iret=iret,                       &  
!         recname=recname,reclevtyp=reclevtyp,reclev=reclev)
!       do i=1,idum
!        print *,'recname=',(trim(recname(i)))
!        print *,'reclevtyp=',(trim(reclevtyp(i)))
!        print *,'reclev=',(reclev(i))
!       end do
	 	
       do j=1,jm
         do i=1,im
	   dummy(i,j)  = glat1d((j-1)*im+i)
	   dummy2(i,j) = glon1d((j-1)*im+i)
	 end do
       end do	   
       deallocate(recname,reclevtyp,reclev,glat1d,glon1d)
! can't get idate and fhour, specify them for now
!       idate(4)=2006
!       idate(2)=9  
!       idate(3)=16
!       idate(1)=0
!       fhour=6.0
        print*,'idate before broadcast = ',(idate(i),i=1,7)
      end if
      call mpi_bcast(idate(1),7,MPI_INTEGER,0,mpi_comm_comp,iret)
      call mpi_bcast(nfhour,1,MPI_INTEGER,0,mpi_comm_comp,iret)
      call mpi_bcast(nframe,1,MPI_INTEGER,0,mpi_comm_comp,iret)
      print*,'idate after broadcast = ',(idate(i),i=1,4)
      print*,'nfhour = ',nfhour
      
! sample print point
      ii=im/2
      jj=jm/2
      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real                   &
       ,gdlat(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,ierr)
      call mpi_scatterv(dummy2(1,1),icnt,idsp,mpi_real                  &
       ,gdlon(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,ierr)
      
      print *,'before call EXCH,mype=',me,'max(gdlat)=',maxval(gdlat),'max(gdlon)=', &
        maxval(gdlon)
      CALL EXCH(gdlat(1,JSTA_2L))
      print *,'after call EXCH,mype=',me

      do j = jsta, jend_m
        do i = 1, im-1
          DX ( i, j ) = ERAD*COS(GDLAT(I,J)*DTR)                        &
      	    *(GDLON(I+1,J)-GDLON(I,J))*DTR  
          DY ( i, j ) =  ERAD*(GDLAT(I,J)-GDLAT(I,J+1))*DTR  ! like A*DPH
!	  F(I,J)=1.454441e-4*sin(gdlat(i,j)*DTR)   ! 2*omeg*sin(phi)
	  IF(i==ii.and.j==jend)print*,'sample LATLON, DY, DY='           &
            ,i,j,GDLAT(I,J),GDLON(I,J),DX(I,J),DY(I,J)
        end do
      end do
      
      do j=jsta,jend
        do i=1,im
	  F(I,J)=1.454441e-4*sin(gdlat(i,j)*DTR)   ! 2*omeg*sin(phi)
	end do
      end do
      
      impf=im
      jmpf=jm
      print*,'impf,jmpf,nframe= ',impf,jmpf,nframe 	   
      
!MEB not sure how to get these      
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

      
!      iyear=idate(4)+2000 ! older gfsio only has 2 digit year
      iyear = idate(1)
      imn   = idate(2) ! ask Jun 
      iday  = idate(3) ! ask Jun
      ihrst = idate(4)
      imin  = idate(5)
      jdate = 0
      idate = 0 
!
!      read(startdate,15)iyear,imn,iday,ihrst,imin       
 15   format(i4,1x,i2,1x,i2,1x,i2,1x,i2)
      print*,'start yr mo day hr min =',iyear,imn,iday,ihrst,imin
      print*,'processing yr mo day hr min='                            &
        ,idat(3),idat(1),idat(2),idat(4),idat(5)
!
      idate(1) = iyear
      idate(2) = imn
      idate(3) = iday
      idate(5) = ihrst
      idate(6) = imin
      SDAT(1)  = imn
      SDAT(2)  = iday
      SDAT(3)  = iyear
      jdate(1) = idat(3)
      jdate(2) = idat(1)
      jdate(3) = idat(2)
      jdate(5) = idat(4)
      jdate(6) = idat(5)
!
      print *,' idate=',idate
      print *,' jdate=',jdate
!      CALL W3DIFDAT(JDATE,IDATE,2,RINC)
!      ifhr=nint(rinc(2))
!
      CALL W3DIFDAT(JDATE,IDATE,0,RINC)
!
      print *,' rinc=',rinc
      ifhr=nint(rinc(2)+rinc(1)*24.)
      print *,' ifhr=',ifhr
      ifmin=nint(rinc(3))
!      if(ifhr /= nint(fhour))print*,'find wrong Grib file';stop
      print*,' in INITPOST ifhr ifmin fileName=',ifhr,ifmin,fileName
      
! GFS has the same accumulation bucket for precipitation and fluxes and it is written to header
! the header has the start hour information so post uses it to recontruct bucket
      if(me==0)then
       call nemsio_getheadvar(ffile,'zhour',zhour,iret=iret)
       if(iret==0)then
        tprec=1.0*ifhr-zhour
	tclod=tprec
	trdlw=tprec
	trdsw=tprec
	tsrfc=tprec
	tmaxmin=tprec
	td3d=tprec
	print*,'tprec from flux file header= ',tprec
       else 
        print*,'Error reading accumulation bucket from flux file header '
	print*,'will try to read bucket from env variable FHZER' 
        CALL GETENV('FHZER',ENVAR)
        read(ENVAR, '(I2)')idum
        tprec=idum*1.0
	tclod=tprec
	trdlw=tprec
	trdsw=tprec
	tsrfc=tprec
	tmaxmin=tprec
	td3d=tprec
        print*,'TPREC from FHZER= ',tprec
       end if
      end if
      
      call mpi_bcast(tprec,1,MPI_REAL,0,mpi_comm_comp,iret)
      call mpi_bcast(tclod,1,MPI_REAL,0,mpi_comm_comp,iret)
      call mpi_bcast(trdlw,1,MPI_REAL,0,mpi_comm_comp,iret)
      call mpi_bcast(trdsw,1,MPI_REAL,0,mpi_comm_comp,iret)
      call mpi_bcast(tsrfc,1,MPI_REAL,0,mpi_comm_comp,iret)
      call mpi_bcast(tmaxmin,1,MPI_REAL,0,mpi_comm_comp,iret)
      call mpi_bcast(td3d,1,MPI_REAL,0,mpi_comm_comp,iret)
      
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
       ifhr    = ifhr+NINT(tstart)
       rinc    = 0
       idate   = 0
       rinc(2) = -1.0*ifhr
       call w3movdat(rinc,jdate,idate)
       SDAT(1) = idate(2)
       SDAT(2) = idate(3)
       SDAT(3) = idate(1)
       IHRST   = idate(5)       
       print*,'new forecast hours for restrt run= ',ifhr
       print*,'new start yr mo day hr min =',sdat(3),sdat(1)           &   
             ,sdat(2),ihrst,imin
      END IF 
      
      imp_physics=99 !set GFS mp physics to 99 for Zhao scheme
      print*,'MP_PHYSICS= ',imp_physics
      
! Initializes constants for Ferrier microphysics       
      if(imp_physics==5 .or. imp_physics==85 .or. imp_physics==95)then
       CALL MICROINIT(imp_physics)
      end if      

! IVEGSRC=1 for IGBP, 0 for USGS, 2 for UMD
      VarName='IVEGSRC'
      if(me == 0)then
        call nemsio_getheadvar(nfile,trim(VarName),IVEGSRC,iret)
        if (iret /= 0) then
          print*,VarName," not found in file-Assigned 2 for UMD as default"
          IVEGSRC=1
        end if
      end if
      call mpi_bcast(IVEGSRC,1,MPI_INTEGER,0,mpi_comm_comp,iret)
      print*,'IVEGSRC= ',IVEGSRC

! set novegtype based on vegetation classification
      if(ivegsrc==2)then
       novegtype=13
      else if(ivegsrc==1)then
       novegtype=20
      else if(ivegsrc==0)then
       novegtype=24
      end if
      print*,'novegtype= ',novegtype

      VarName='CU_PHYSICS'
      if(me == 0)then
        call nemsio_getheadvar(nfile,trim(VarName),iCU_PHYSICS,iret)
        if (iret /= 0) then
          print*,VarName," not found in file-Assigned 4 for SAS as default"
          iCU_PHYSICS=4
        end if
      end if
      call mpi_bcast(iCU_PHYSICS,1,MPI_INTEGER,0,mpi_comm_comp,iret)
      print*,'CU_PHYSICS= ',iCU_PHYSICS
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
!      maptype=4
!      write(6,*) 'maptype is ', maptype	  
! HBM2 is most likely not in Grib message, set them to ones
      HBM2=1.0

! try to get kgds from flux grib file and then convert to igds that is used by GRIBIT.f
! flux files are now nemsio files so comment the following lines out
!      if(me == 0)then       
!       jpds=-1.0
!       jgds=-1.0
!       igds=0                                                                                
!       call getgb(iunit,0,im_jm,0,jpds,jgds,kf                          &  
!          ,k,kpds,kgds,lb,dummy,ierr)
!       if(ierr == 0)then
!        call R63W72(KPDS,KGDS,JPDS,IGDS(1:18))
!       print*,'in INITPOST_GFS,IGDS for GFS= ',(IGDS(I),I=1,18)
!       end if
!      end if
!      call mpi_bcast(igds(1),18,MPI_INTEGER,0,mpi_comm_comp,iret)      
!      print*,'IGDS for GFS= ',(IGDS(I),I=1,18)
      
! Specigy grid type
!      if(iostatusFlux==0)then
      if(IGDS(4)/=0)then
       maptype=IGDS(3)
      else if((im/2+1)==jm)then
       maptype=0 !latlon grid
      else
       maptype=4 ! default gaussian grid
      end if
      gridtype='A'

      write(6,*) 'maptype and gridtype is ', maptype,gridtype      
      
! start retrieving data using gfsio, first land/sea mask

!      VarName='land'
!      VcoordName='sfc'
!      l=1
      
!      if(me == 0)then
!        call gfsio_readrecvw34(gfile,trim(VarName),trim(VcoordName)
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
      VarName='land'  
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,impf,jmpf,nframe,sm)
      where(sm /= spval)sm=1.0-sm ! convert to sea mask
      if(debugprint)print*,'sample ',VarName,' = ',sm(im/2,(jsta+jend)/2)
      

! sea ice mask using getgb
      
      VarName='icec'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,sice)
      if(debugprint)print*,'sample ',VarName,' = ',sice(im/2,(jsta+jend)/2)
      
!      where(sice /=spval .and. sice >=1.0)sm=0.0 !sea ice has sea mask=0
! GFS flux files have land points with non-zero sea ice, per Iredell, these
! points have sea ice changed to zero, i.e., trust land mask more than sea ice
      where(sm/=spval .and. sm==0.0)sice=0.0 !specify sea ice=0 at land

! GFS does not output PD  
      pd=spval

! Terrain height * G   using nemsio 
      VarName='hgt'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,fis)
      where(fis /= spval)fis=fis*con_G
      if(debugprint)print*,'sample ',VarName,' = ',fis(im/2,(jsta+jend)/2)
      
! model level T
      print*,'start retrieving GFS T using nemsio'
      VarName='tmp'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,t(1,jsta_2l,ll))
        if(debugprint)print*,'sample ',ll,VarName,' = ',ll,t(im/2,(jsta+jend)/2,ll)      
      end do ! do loop for l

! model level q      
      VarName='spfh'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,q(1,jsta_2l,ll))
        if(debugprint)print*,'sample ',ll,VarName,' = ',ll,q(im/2,(jsta+jend)/2,ll)      
      end do ! do loop for l
	
! model level u      
      VarName='ugrd'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,uh(1,jsta_2l,ll))
        if(debugprint)print*,'sample ',ll,VarName,' = ',ll,uh(im/2,(jsta+jend)/2,ll)      
      end do ! do loop for l
      
! model level v      
      VarName='vgrd'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,vh(1,jsta_2l,ll))
        if(debugprint)print*,'sample ',ll,VarName,' = ',ll,vh(im/2,(jsta+jend)/2,ll)      
      end do ! do loop for l
      
! model level pressure      
      VarName='pres'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,pmid(1,jsta_2l,ll))
        if(debugprint)print*,'sample ',ll,VarName,' = ',ll,pmid(im/2,(jsta+jend)/2,ll)      
      end do ! do loop for l
      
! GFS is on A grid and does not need PMIDV        

! Surface pressure  using nemsio 
      VarName='pres'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,pint(1,jsta_2l,lp1))
      if(debugprint)print*,'sample surface pressure = ',pint(im/2,(jsta+jend)/2,lp1)
     
      do j=jsta,jend
	do i=1,im
	  ALPINT(I,J,LP1)=ALOG(PINT(I,J,LP1))     
!	    if(i==ii .and. j==jj)print*,'sample PSFC'
!     +        ,i,j,pint(i,j,lp1)	    
	end do
      end do       
      
! dp     
      VarName='dpres'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,dpres(1,jsta_2l,ll))
!      ,l,im,jm,nframe,buf3d(1,jsta_2l,ll))
!        if(debugprint)print*,'sample ',ll,VarName,' = ',ll,pmid(im/2,(jsta+jend)/2,ll)      
      end do ! do loop for l

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
!	    pint(i,j,l)=pint(i,j,l-1)+buf3d(i,j,l-1)
	    pint(i,j,l)=pint(i,j,l-1)+dpres(i,j,l-1)
	    ALPINT(I,J,L)=ALOG(PINT(I,J,L))     
	    if(i==ii .and. j==jj)print*,'sample pint,pmid'                     &   
              ,i,j,l,pint(i,j,l),pmid(i,j,l)	    
	  end do
	end do
      end do
      
!!!!! COMPUTE Z, GFS integrates Z on mid-layer instead
!!! use GFS contants to see if height becomes more aggreable to GFS pressure grib file

       do j = jsta, jend
        do i = 1, im
            ZINT(I,J,LM+1)=FIS(I,J)/con_G
	    FI(I,J,1)=FIS(I,J)+T(I,J,LM)                                &  
      	    *(Q(I,J,LM)*con_fvirt+1.0)*con_rd                           &
! using GFS consts	  
            *(ALPINT(I,J,LM+1)-ALOG(PMID(I,J,LM)))                      
            ZMID(I,J,LM)=FI(I,J,1)/con_G
!            FI(I,J,1)=FIS(I,J)
        end do
       end do

! SECOND, INTEGRATE HEIGHT HYDROSTATICLY, GFS integrate height on mid-layer
      DO L=LM-1,1,-1
       do j = jsta, jend
        do i = 1, im
         FI(I,J,2)=0.5*(T(I,J,L)*(Q(I,J,L)*con_fvirt+1.0)               &   
      	           +T(I,J,L+1)*(Q(I,J,L+1)*con_fvirt+1.0))*con_rd*      &
!         FI(I,J,2)=0.5*(T(I,J,L)*(Q(I,J,L)*D608+1.0)     	 
!     1	           +T(I,J,L+1)*(Q(I,J,L+1)*D608+1.0))*RD*
                   (ALOG(PMID(I,J,L+1))-ALOG(PMID(I,J,L)))              &
                   +FI(I,J,1)
         ZMID(I,J,L)=FI(I,J,2)/con_G
	 
         if(i.eq.ii.and.j.eq.jj)                                        &
        print*,'L,sample T,Q,ALPMID(L+1),ALPMID(L),ZMID= '              &
        ,l,T(I,J,L),Q(I,J,L),ALOG(PMID(I,J,L+1)),                       &
        ALOG(PMID(I,J,L)),ZMID(I,J,L)
     
         FI(I,J,1)=FI(I,J,2)
        ENDDO
       ENDDO
      END DO

      DO L=LM,2,-1  ! omit computing model top height because it's infinity
       DO J=JSTA,JEND
        DO I=1,IM
!         ZMID(I,J,L)=(ZINT(I,J,L+1)+ZINT(I,J,L))*0.5  ! ave of z
         FACT=(ALPINT(I,J,L)-ALOG(PMID(I,J,L)))/                        &  
               (ALOG(PMID(I,J,L-1))-ALOG(PMID(I,J,L)))          
         ZINT(I,J,L)=ZMID(I,J,L)+(ZMID(I,J,L-1)-ZMID(I,J,L))            &
                       *FACT 
         if(i.eq.ii.and.j.eq.jj) print*,'L ZINT= ',l,zint(i,j,l)	 	 
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

! ozone mixing ratio     
      VarName='o3mr'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,o3(1,jsta_2l,ll))
        if(debugprint)print*,'sample ',ll,VarName,' = ',ll,o3(im/2,(jsta+jend)/2,ll)      
      end do ! do loop for l      

! cloud water and ice mixing ratio  for zhao scheme  
! need to look up old eta post to derive cloud water/ice from cwm 
! Zhao scheme does not produce suspended rain and snow
      qqr=0.
      qqs=0.
!      qqi=0.
      VarName='clwmr'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,cwm(1,jsta_2l,ll))
        if(debugprint)print*,'sample ',ll,VarName,' = ',ll,cwm(im/2,(jsta+jend)/2,ll)      

        do j=jsta,jend
         do i=1,im
	  if(t(i,j,ll) < (TFRZ-15.) )then ! dividing cloud water from ice
	   qqi(i,j,ll)=cwm(i,j,ll)
	  else 
	   qqw(i,j,ll)=cwm(i,j,ll)
	  end if 
!         if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     +   print*,'sample ',trim(VarName), ' after scatter= '
!     +   ,i,j,ll,cwm(i,j,ll)
         end do
        end do
!       if (iret /= 0)print*,'Error scattering array';stop
      end do ! do loop for l 
      
! GFS does output omeg now
      VarName='vvel'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,omga(1,jsta_2l,ll))
        if(debugprint)print*,'sample l ',VarName,' = ',ll,omga(im/2,(jsta+jend)/2,ll)      
      end do ! do loop for l
      
! GFS output dust in nemsio (GOCART)
      DUST = SPVAL
      VarName='du001'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,dust(1,jsta_2l,ll,1))
        if(debugprint)print*,'sample l ',VarName,' = ',ll,dust(im/2,(jsta+jend)/2,ll,1)      
      end do ! do loop for l      
      
      VarName='du002'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,dust(1,jsta_2l,ll,2))
        if(debugprint)print*,'sample l ',VarName,' = ',ll,dust(im/2,(jsta+jend)/2,ll,2)      
      end do ! do loop for l 
      
      VarName='du003'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,dust(1,jsta_2l,ll,3))
        if(debugprint)print*,'sample l ',VarName,' = ',ll,dust(im/2,(jsta+jend)/2,ll,3)      
      end do ! do loop for l 
      
      VarName='du004'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,dust(1,jsta_2l,ll,4))
        if(debugprint)print*,'sample l ',VarName,' = ',ll,dust(im/2,(jsta+jend)/2,ll,4)      
      end do ! do loop for l 
      
      VarName='du005'
      VcoordName='mid layer'
      do l=1,lm	
        ll=lm-l+1
        call getnemsandscatter(me,nfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,dust(1,jsta_2l,ll,5))
        if(debugprint)print*,'sample l ',VarName,' = ',ll,dust(im/2,(jsta+jend)/2,ll,5)      
      end do ! do loop for l 

! -- compute air density RHOMID and remove negative tracer values
      do l=1,lm
        do j=jsta,jend
          do i=1,im
           TV=T(I,J,L)*(H1+D608*MAX(Q(I,J,L),QMIN))
           RHOMID(I,J,L)=PMID(I,J,L)/(RD*TV)
           do n = 1,  NBIN_DU
             IF ( dust(i,j,l,n) .LT. SPVAL) THEN
               DUST(i,j,l,n) = MAX(DUST(i,j,l,n), 0.0)    
             ENDIF
           enddo
          end do
        end do
      end do
!

! PBL height using nemsio
      VarName='hpbl'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,pblh)
      if(debugprint)print*,'sample ',VarName,' = ',pblh(im/2,(jsta+jend)/2)

! frictional velocity using nemsio
      VarName='fricv'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,ustar)
      if(debugprint)print*,'sample ',VarName,' = ',ustar(im/2,(jsta+jend)/2)

! roughness length using getgb
      VarName='sfcr'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,z0)
      if(debugprint)print*,'sample ',VarName,' = ',z0(im/2,(jsta+jend)/2)

! surface potential T  using getgb
      VarName='tmp'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,ths)
      where(ths/=spval)ths=ths*(p1000/pint(:,:,lp1))**CAPA ! convert to THS
      if(debugprint)print*,'sample ',VarName,' = ',ths(im/2,(jsta+jend)/2)

! GFS does not have surface specific humidity
      QS=SPVAL           

! GFS does not have inst sensible heat flux
      twbs=SPVAL   
      
! GFS does not have inst latent heat flux
      qwbs=SPVAL    
          
!  GFS does not have time step and physics time step, make up ones since they
! are not really used anyway
      NPHS=2.
      DT=80.
      DTQ2 = DT * NPHS  !MEB need to get physics DT
      TSPH = 3600./DT   !MEB need to get DT
! All GFS time-averaged quantities are in 6 hour bucket
!      TPREC=6.0

! convective precip in m per physics time step using gfsio
!      VarName='cprat'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecvw34(gfile,trim(VarName),trim(VcoordName)    &
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
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real                  &
!     + , avgcprate(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! convective precip in m per physics time step using getgb
      VarName='cprat'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,avgcprate)
      where(avgcprate /= spval)avgcprate=avgcprate*dtq2/1000. ! convert to m
      if(debugprint)print*,'sample ',VarName,' = ',avgcprate(im/2,(jsta+jend)/2)
      
      cprate=avgcprate   
!      print*,'maxval CPRATE: ', maxval(CPRATE)

! precip rate in m per physics time step using getgb
      VarName='prate'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,avgprec)
      where(avgprec /= spval)avgprec=avgprec*dtq2/1000. ! convert to m
      if(debugprint)print*,'sample ',VarName,' = ',avgprec(im/2,(jsta+jend)/2)
      
      prec=avgprec !set avg cprate to inst one to derive other fields

! GFS does not have accumulated total, gridscale, and convective precip, will use inst precip to derive in SURFCE.f

      
! GFS does not have similated precip
      lspa=spval  

! inst snow water eqivalent using nemsio
      VarName='weasd'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,sno)
      if(debugprint)print*,'sample ',VarName,' = ',sno(im/2,(jsta+jend)/2)

! snow depth in mm using nemsio
      VarName='snod'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,si)
      where(si /= spval)si=si*1000. ! convert to mm
      if(debugprint)print*,'sample ',VarName,' = ',si(im/2,(jsta+jend)/2)

! GFS does not have convective cloud efficiency
      CLDEFI=SPVAL
      
! GFS does not have 10 m theta
      TH10=SPVAL

! GFS does not have 10 m humidity
      Q10=SPVAL
      
! 2m T using nemsio
      VarName='tmp'
      VcoordName='2 m above gnd'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,tshltr)
      if(debugprint)print*,'sample ',VarName,' = ',tshltr(im/2,(jsta+jend)/2)

! GFS does not have 2m pres, estimate it, also convert t to theta 
      Do j=jsta,jend
        Do i=1,im
          PSHLTR(I,J)=pint(I,J,lm+1)*EXP(-0.068283/tshltr(i,j))
          tshltr(i,j)= tshltr(i,j)*(p1000/PSHLTR(I,J))**CAPA ! convert to theta
!          if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     +   print*,'sample 2m T and P after scatter= '
!     +   ,i,j,tshltr(i,j),pshltr(i,j)
        end do
      end do

! 2m specific humidity using gfsio                    
!      VarName='spfh'
!      VcoordName='2m above gnc'
!      l=1
!      if(me == 0)then
!        call gfsio_readrecvw34(gfile,trim(VarName),trim(VcoordName)    &
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real                  &
!     + ,qshltr(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop

! 2m specific humidity using nemsio
      VarName='spfh'
      VcoordName='2 m above gnd'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,qshltr)
      if(debugprint)print*,'sample ',VarName,' = ',qshltr(im/2,(jsta+jend)/2)
      
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
!        call gfsio_readrecvw34(gfile,trim(VarName),trim(VcoordName)    &
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
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real                  &
!     + ,avgalbedo(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! mid day avg albedo in fraction using nemsio
      VarName='albdo'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,avgalbedo)
      where(avgalbedo /= spval)avgalbedo=avgalbedo/100. ! convert to fraction
      if(debugprint)print*,'sample ',VarName,' = ',avgalbedo(im/2,(jsta+jend)/2)
     
! time averaged column cloud fractionusing nemsio
      VarName='tcdc'
      VcoordName='atmos col'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,avgtcdc)
      where(avgtcdc /= spval)avgtcdc=avgtcdc/100. ! convert to fraction
      if(debugprint)print*,'sample ',VarName,' = ',avgtcdc(im/2,(jsta+jend)/2)
	
! GFS probably does not use zenith angle
      Czen=spval
      CZMEAN=SPVAL      

! maximum snow albedo in fraction using nemsio
      VarName='mxsalb'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,mxsnal)
      where(mxsnal /= spval)mxsnal=mxsnal/100. ! convert to fraction
      if(debugprint)print*,'sample ',VarName,' = ',mxsnal(im/2,(jsta+jend)/2)
     
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

! GFS does not have model level cloud fraction -> derive cloud fraction
!      CFR=SPVAL
      allocate(qstl(lm))
      print*,'start deriving cloud fraction'
!      do j=jsta,jend
!        do i=1,im
!	  do l=1,lm
!	    if(i==im/2.and.j==jsta)print*,'sample T=',t(i,j,l)
!	    es=fpvsnew(t(i,j,l))
!	    if(i==im/2.and.j==jsta)print*,'sample ES=',es
!	    es=min(es,pmid(i,j,l))
!	    if(i==im/2.and.j==jsta)print*,'sample ES=',es
!	    qstl(l)=con_eps*es/(pmid(i,j,l)+con_epsm1*es) !saturation q for GFS
!          end do
!          call progcld1                                               
!...................................

!  ---  inputs:
!     &       ( pmid(i,j,1:lm)/100.,pint(i,j,1:lm+1)/100.,
!     &         t(i,j,1:lm),q(i,j,1:lm),qstl,cwm(i,j,1:lm),                         
!     &         gdlat(i,j),gdlon(i,j),                                  
!     &         1, lm, lm+1, 0,                                
!  ---  outputs:
!     &         cfr(i,j,1:lm)                                      
!     &        )
!          do l=1,lm
!	    cfr(i,j,l)=cldtot(l)
!	  end do   
!        end do
!      end do     
       
      do j=jsta,jend
        do i=1,im
	  do k = 1,lm
!	    if(i==im/2.and.j==jsta)print*,'sample T=',t(i,j,k)
	    es=fpvsnew(t(i,j,k))
!	    if(i==im/2.and.j==jsta)print*,'sample ES=',es
	    es=min(es,pmid(i,j,k))
!	    if(i==im/2.and.j==jsta)print*,'sample ES=',es
             if(pmid(i,j,k)>1.0)                                        &   
      	    qstl(k)=con_eps*es/(pmid(i,j,k)+con_epsm1*es) !saturation q for GFS
!          if(i==im/2.and.j==jsta)print*,'sample qstl=',k,qstl(k)  
          end do  
          call progcld1                                                 &
!...................................

!  ---  inputs:
             ( pmid(i,j,1:lm)/100.,t(i,j,1:lm),                         &
               q(i,j,1:lm),qstl(1:lm),cwm(i,j,1:lm),                    &    
               gdlat(i,j),gdlon(i,j),                                   &
               1, lm, 0,                                                &
!  ---  outputs:
               cfr(i,j,1:lm)                                            &
              )
!          do l=1,lm
!	    cfr(i,j,l)=cldtot(l)
!	  end do
        end do
      end do            
      deallocate(qstl)

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

! GFS does not have inst cloud fraction for high, middle, and low cloud
      cfrach=spval
      cfracl=spval
      cfracm=spval

! ave high cloud fraction using nemsio
      VarName='tcdc'
      VcoordName='high cld lay'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,avgcfrach)
      where(avgcfrach /= spval)avgcfrach=avgcfrach/100. ! convert to fraction
      if(debugprint)print*,'sample ',VarName,' = ',avgcfrach(im/2,(jsta+jend)/2)

! ave low cloud fraction using nemsio
      VarName='tcdc'
      VcoordName='low cld lay'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,avgcfracl)
      where(avgcfracl /= spval)avgcfracl=avgcfracl/100. ! convert to fraction
      if(debugprint)print*,'sample ',VarName,' = ',avgcfracl(im/2,(jsta+jend)/2)
      
! ave middle cloud fraction using nemsio
      VarName='tcdc'
      VcoordName='mid cld lay'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,avgcfracm)
      where(avgcfracm /= spval)avgcfracm=avgcfracm/100. ! convert to fraction
      if(debugprint)print*,'sample ',VarName,' = ',avgcfracm(im/2,(jsta+jend)/2)
      
! inst convective cloud fraction using nemsio
      VarName='tcdc'
      VcoordName='convect-cld laye'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,cnvcfr)
      where(cnvcfr /= spval)cnvcfr=cnvcfr/100. ! convert to fraction
      if(debugprint)print*,'sample ',VarName,' = ',cnvcfr(im/2,(jsta+jend)/2)
      
! slope type using nemsio
      VarName='sltyp'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,buf)
      where(buf /= spval)islope=nint(buf) 
      if(debugprint)print*,'sample ',VarName,' = ',islope(im/2,(jsta+jend)/2)

! plant canopy sfc wtr in m using nemsio
      VarName='cnwat'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,cmc)
      where(cmc /= spval)cmc=cmc/1000. ! convert from kg*m^2 to m
      if(debugprint)print*,'sample ',VarName,' = ',cmc(im/2,(jsta+jend)/2)
      
! GFS does not have inst ground heat flux
      grnflx=spval    

! GFS does not have snow cover yet
!      VarName='gflux'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecvw34(gfile,trim(VarName),trim(VcoordName)    &
!     +	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy,icnt,idsp,mpi_real                  &
!     + , pctsno(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      
! asuume tg3 in GFS is the same as soiltb in wrf nmm. It's in sfc file, will
! be able to read it when it merges to gfs io
! soiltb is not being put out, comment it out
!      VarName='tg3'
!      VcoordName='sfc' 
!      l=1
!      if(me == 0)then
!        call gfsio_readrecvw34(gfile,trim(VarName),trim(VcoordName)    &
!      	,l,dummy,iret=iret)
!        if (iret /= 0) then
!         print*,VarName," not found in file-Assigned missing values"
!         dummy=spval
!        end if
!      end if	
!      call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real                  &
!       , soiltb(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,iret)
!      if (iret /= 0)print*,'Error scattering array';stop
      	
! vegetation fraction in fraction. using nemsio
      VarName='veg'
      VcoordName='sfc'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,vegfrc)
      where(vegfrc /= spval)
       vegfrc=vegfrc/100. ! convert to fraction
      elsewhere (vegfrc == spval)
       vegfrc=0. ! set to zero to be reasonable input for crtm
      end where
      if(debugprint)print*,'sample ',VarName,' = ',vegfrc(im/2,(jsta+jend)/2)
      
! GFS doesn not yet output soil layer thickness, assign SLDPTH to be the same as nam

         SLDPTH(1)=0.10
         SLDPTH(2)=0.3
         SLDPTH(3)=0.6
         SLDPTH(4)=1.0
	 
! liquid volumetric soil mpisture in fraction using nemsio
      VarName='soill'
      VcoordName='0-10 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,sh2o(1,jsta_2l,1))
      if(debugprint)print*,'sample l',VarName,' = ',1,sh2o(im/2,(jsta+jend)/2,1)
      
      VarName='soill'
      VcoordName='10-40 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,sh2o(1,jsta_2l,2))
      if(debugprint)print*,'sample l',VarName,' = ',1,sh2o(im/2,(jsta+jend)/2,2)
      
      VarName='soill'
      VcoordName='40-100 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,sh2o(1,jsta_2l,3))
      if(debugprint)print*,'sample l',VarName,' = ',1,sh2o(im/2,(jsta+jend)/2,3)
      
      VarName='soill'
      VcoordName='100-200 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,sh2o(1,jsta_2l,4))
      if(debugprint)print*,'sample l',VarName,' = ',1,sh2o(im/2,(jsta+jend)/2,4)
      
! volumetric soil moisture using nemsio
      VarName='soilw'
      VcoordName='0-10 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,smc(1,jsta_2l,1))
      if(debugprint)print*,'sample l',VarName,' = ',1,smc(im/2,(jsta+jend)/2,1)
      
      VarName='soilw'
      VcoordName='10-40 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,smc(1,jsta_2l,2))
      if(debugprint)print*,'sample l',VarName,' = ',1,smc(im/2,(jsta+jend)/2,2)
      
      VarName='soilw'
      VcoordName='40-100 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,smc(1,jsta_2l,3))
      if(debugprint)print*,'sample l',VarName,' = ',1,smc(im/2,(jsta+jend)/2,3)
      
      VarName='soilw'
      VcoordName='100-200 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,smc(1,jsta_2l,4))
      if(debugprint)print*,'sample l',VarName,' = ',1,smc(im/2,(jsta+jend)/2,4)

! soil temperature using nemsio
      VarName='tmp'
      VcoordName='0-10 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,stc(1,jsta_2l,1))
      if(debugprint)print*,'sample l','stc',' = ',1,stc(im/2,(jsta+jend)/2,1)
      
      VarName='tmp'
      VcoordName='10-40 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,stc(1,jsta_2l,2))
      if(debugprint)print*,'sample stc = ',1,stc(im/2,(jsta+jend)/2,2)
      
      VarName='tmp'
      VcoordName='40-100 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,stc(1,jsta_2l,3))
      if(debugprint)print*,'sample stc = ',1,stc(im/2,(jsta+jend)/2,3)
      
      VarName='tmp'
      VcoordName='100-200 cm down'
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,stc(1,jsta_2l,4))
      if(debugprint)print*,'sample stc = ',1,stc(im/2,(jsta+jend)/2,4)

     
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
!      trdlw=6.0

! GFS does not have inst incoming sfc longwave
      rlwin=spval
       
! GFS does not have inst model top outgoing longwave
      rlwtoa=spval

! time averaged incoming sfc longwave using nemsio
      VarName='dlwrf'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,alwin)
                                                            
! time averaged outgoing sfc longwave using gfsio
      VarName='ulwrf'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,alwout)
      where(alwout /= spval) alwout=-alwout ! CLDRAD puts a minus sign before gribbing
      if(debugprint)print*,'sample l',VarName,' = ',1,alwout(im/2,(jsta+jend)/2)
                                                                                          
! time averaged outgoing model top longwave using gfsio
      VarName='ulwrf'
      VcoordName='nom. top' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,alwtoa)
      if(debugprint)print*,'sample l',VarName,' = ',1,alwtoa(im/2,(jsta+jend)/2)                                                 
      
! GFS does not have inst incoming sfc shortwave
      rswin=spval 

! GFS does not have inst incoming clear sky sfc shortwave
      rswinc=spval      

! GFS does not have inst outgoing sfc shortwave
      rswout=spval
           
! GFS incoming sfc longwave has been averaged, set ARDLW to 1
      ardsw=1.0
!      trdsw=6.0

! time averaged incoming sfc shortwave using gfsio
      VarName='dswrf'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,aswin)
      if(debugprint)print*,'sample l',VarName,' = ',1,aswin(im/2,(jsta+jend)/2)

! time averaged incoming sfc uv-b using getgb
      VarName='duvb'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,auvbin)
      if(debugprint)print*,'sample l',VarName,' = ',1,auvbin(im/2,(jsta+jend)/2)
       
! time averaged incoming sfc clear sky uv-b using getgb
      VarName='cduvb'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,auvbinc)
      if(debugprint)print*,'sample l',VarName,' = ',1,auvbinc(im/2,(jsta+jend)/2)
      
! time averaged outgoing sfc shortwave using gfsio
      VarName='uswrf'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,aswout)
      where(aswout /= spval) aswout=-aswout ! CLDRAD puts a minus sign before gribbing 
      if(debugprint)print*,'sample l',VarName,' = ',1,aswout(im/2,(jsta+jend)/2)
      
! time averaged model top outgoing shortwave
      VarName='uswrf'
      VcoordName='nom. top' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,aswtoa)
      if(debugprint)print*,'sample l',VarName,' = ',1,aswtoa(im/2,(jsta+jend)/2)
                                                                              
! time averaged surface sensible heat flux, multiplied by -1 because wrf model flux
! has reversed sign convention using gfsio
      VarName='shtfl'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,sfcshx)
      where (sfcshx /= spval)sfcshx=-sfcshx
      if(debugprint)print*,'sample l',VarName,' = ',1,sfcshx(im/2,(jsta+jend)/2)

! GFS surface flux has been averaged, set  ASRFC to 1 
      asrfc=1.0  
!      tsrfc=6.0

! time averaged surface latent heat flux, multiplied by -1 because wrf model flux
! has reversed sign vonvention using gfsio
      VarName='lhtfl'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,sfclhx)
      where (sfclhx /= spval)sfclhx=-sfclhx
      if(debugprint)print*,'sample l',VarName,' = ',1,sfclhx(im/2,(jsta+jend)/2)
                                                                                                
! time averaged ground heat flux using nemsio
      VarName='gflux'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,subshx)
      if(debugprint)print*,'sample l',VarName,' = ',1,subshx(im/2,(jsta+jend)/2)                                                

! GFS does not have snow phase change heat flux
      snopcx=spval

! time averaged zonal momentum flux using gfsio
      VarName='uflx'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,sfcux)
      if(debugprint)print*,'sample l',VarName,' = ',1,sfcux(im/2,(jsta+jend)/2)
      
! time averaged meridional momentum flux using nemsio
      VarName='vflx'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,sfcvx)
      if(debugprint)print*,'sample l',VarName,' = ',1,sfcvx(im/2,(jsta+jend)/2)
     
! GFS does not use total momentum flux
      sfcuvx=spval

! time averaged zonal gravity wave stress using nemsio
      VarName='u-gwd'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,gtaux)
      if(debugprint)print*,'sample l',VarName,' = ',1,gtaux(im/2,(jsta+jend)/2)
                                                                                          
! time averaged meridional gravity wave stress using getgb
      VarName='v-gwd'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,gtauy)
      if(debugprint)print*,'sample l',VarName,' = ',1,gtauy(im/2,(jsta+jend)/2)
                                                     
! time averaged accumulated potential evaporation
      VarName='pevpr'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,potevp)
      if(debugprint)print*,'sample l',VarName,' = ',1,potevp(im/2,(jsta+jend)/2)                                                 

! GFS does not have temperature tendency due to long wave radiation
      rlwtt=spval
      
! GFS does not have temperature tendency due to long wave radiation
      rswtt=spval
      
! GFS does not have temperature tendency due to latent heating from convection
      tcucn=spval
      tcucns=spval

! set avrain to 1
      avrain=1.0
      avcnvc=1.0
      theat=6.0 ! just in case GFS decides to output T tendency   
      
! GFS does not have temperature tendency due to latent heating from grid scale
      train=spval

! 10 m u using nemsio
      VarName='ugrd'
      VcoordName='10 m above gnd' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,u10)
      if(debugprint)print*,'sample l',VarName,' = ',1,u10(im/2,(jsta+jend)/2)
            
! 10 m v using gfsio
      VarName='vgrd'
      VcoordName='10 m above gnd' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,v10)
      if(debugprint)print*,'sample l',VarName,' = ',1,v10(im/2,(jsta+jend)/2)
      
! GFS does not have soil moisture availability 
      smstav=spval

! GFS does not have total soil moisture 
      smstot=spval
      
! vegetation type, it's in GFS surface file, hopefully will merge into gfsio soon 
      VarName='vgtyp'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,buf)
      where (buf /= spval)
       ivgtyp=nint(buf)
      elsewhere
       ivgtyp=0 !need to feed reasonable value to crtm
      end where 
      if(debugprint)print*,'sample l',VarName,' = ',1,ivgtyp(im/2,(jsta+jend)/2)
      
! soil type, it's in GFS surface file, hopefully will merge into gfsio soon
      VarName='sotyp'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,buf)
      where (buf /= spval)
       isltyp=nint(buf)
      elsewhere
       isltyp=0 !need to feed reasonable value to crtm
      end where 
      if(debugprint)print*,'sample l',VarName,' = ',1,isltyp(im/2,(jsta+jend)/2)
      
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
      EL_PBL=spval      

! GFS does not output exchange coefficient
      exch_h=spval
      
! GFS does not have THZ0, use THS to substitute
      thz0=ths
      if(debugprint)print*,'sample l',VarName,' = ',1,thz0(im/2,(jsta+jend)/2)

! GFS does not output humidity at roughness length
      qz0=spval
      
! GFS does not output u at roughness length
      uz0=spval
      
! GFS does not output humidity at roughness length
      vz0=spval      

! retrieve inst convective cloud top, GFS has cloud top pressure instead of index,
! will need to modify CLDRAD.f to use pressure directly instead of index
      VarName='pres'
      VcoordName='convect-cld top' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,ptop)
      if(debugprint)print*,'sample l',VarName,' = ',1,ptop(im/2,(jsta+jend)/2)
      
      htop=spval	
      do j=jsta,jend
        do i=1,im
	  if(ptop(i,j) <= 0.0)ptop(i,j)=spval
	  if(ptop(i,j) < spval)then
	   do l=1,lm
	    if(ptop(i,j) <= pmid(i,j,l))then
	     htop(i,j)=l
	     if(i==ii .and. j==jj)print*,'sample ptop,pmid pmid-1,pint= ',   &
      	        ptop(i,j),pmid(i,j,l),pmid(i,j,l-1),pint(i,j,l),htop(i,j)
             exit
	    end if
	   end do
	  end if 
        end do
       end do

! retrieve inst convective cloud bottom, GFS has cloud top pressure instead of index,
! will need to modify CLDRAD.f to use pressure directly instead of index
      VarName='pres'
      VcoordName='convect-cld bot' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,pbot)
      if(debugprint)print*,'sample l',VarName,VcoordName,' = ',1,pbot(im/2,(jsta+jend)/2)
      
      hbot=spval 
      do j=jsta,jend
        do i=1,im
	  if(pbot(i,j) <= 0.0)pbot(i,j)=spval
!	  if(.not.lb(i,j))print*,'false bitmask for pbot at '
!     +	    ,i,j,pbot(i,j)
          if(pbot(i,j) .lt. spval)then
	   do l=lm,1,-1
	    if(pbot(i,j) >= pmid(i,j,l))then
	     hbot(i,j)=l
	     if(i==ii .and. j==jj)print*,'sample pbot,pmid= ',    &
      	        pbot(i,j),pmid(i,j,l),hbot(i,j)
             exit
	    end if
	   end do
	  end if 
        end do
       end do	    

! retrieve time averaged low cloud top pressure using nemsio
      VarName='pres'
      VcoordName='low cld top' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,ptopl)
      if(debugprint)print*,'sample l',VarName,' = ',1,ptopl(im/2,(jsta+jend)/2)                                                                         

! retrieve time averaged low cloud bottom pressure using nemsio
      VarName='pres'
      VcoordName='low cld bot' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,pbotl)
      if(debugprint)print*,'sample l',VarName,' = ',1,pbotl(im/2,(jsta+jend)/2)
     
! retrieve time averaged low cloud top temperature using nemsio
      VarName='tmp'
      VcoordName='low cld top' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,Ttopl)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,Ttopl(im/2,(jsta+jend)/2)
                                                                                               
! retrieve time averaged middle cloud top pressure using nemsio
      VarName='pres'
      VcoordName='mid cld top' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,ptopm)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,ptopm(im/2,(jsta+jend)/2)
                                                             
! retrieve time averaged middle cloud bottom pressure using  nemsio
      VarName='pres'
      VcoordName='mid cld bot' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,pbotm)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,pbotm(im/2,(jsta+jend)/2)
      
! retrieve time averaged middle cloud top temperature using nemsio
      VarName='tmp'
      VcoordName='mid cld top' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,Ttopm)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,Ttopm(im/2,(jsta+jend)/2)
      
! retrieve time averaged high cloud top pressure using nemsio *********
      VarName='pres'
      VcoordName='high cld top' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,ptoph)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,ptoph(im/2,(jsta+jend)/2)
     
! retrieve time averaged high cloud bottom pressure using  nemsio
      VarName='pres'
      VcoordName='high cld bot' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,pboth)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,pboth(im/2,(jsta+jend)/2)
                                                                                               
! retrieve time averaged high cloud top temperature using nemsio
      VarName='tmp'
      VcoordName='high cld top' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,Ttoph)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,Ttoph(im/2,(jsta+jend)/2)
      
! retrieve boundary layer cloud cover using nemsio
      VarName='tcdc'
      VcoordName='bndary-layer cld' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,pblcfr)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,pblcfr(im/2,(jsta+jend)/2)
      where (pblcfr /= spval)pblcfr=pblcfr/100. ! convert to fraction
        
! retrieve cloud work function using nemsio
      VarName='cwork'
      VcoordName='atmos col' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,cldwork)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,cldwork(im/2,(jsta+jend)/2)
      
! retrieve water runoff using nemsio
      VarName='watr'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,runoff)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,runoff(im/2,(jsta+jend)/2)
      
! retrieve shelter max temperature using nemsio
      VarName='tmax'
      VcoordName='2 m above gnd' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,maxtshltr)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,maxtshltr(im/2,(jsta+jend)/2)     

! retrieve shelter max temperature using nemsio
      VarName='tmin'
      VcoordName='2 m above gnd' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,mintshltr)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,mintshltr(im/2,(jsta+jend)/2)
 
      MAXRHSHLTR=SPVAL
      MINRHSHLTR=SPVAL
      
! retrieve ice thickness using nemsio
      VarName='icetk'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,dzice)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,dzice(im/2,(jsta+jend)/2)

! retrieve wilting point using nemsio
      VarName='wilt'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,smcwlt)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,smcwlt(im/2,(jsta+jend)/2)
      
! retrieve sunshine duration using nemsio
      VarName='sunsd'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,suntime)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,suntime(im/2,(jsta+jend)/2)

! retrieve field capacity using nemsio
      VarName='fldcp'
      VcoordName='sfc' 
      l=1
      call getnemsandscatter(me,ffile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,fieldcapa)
      if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
      1,fieldcapa(im/2,(jsta+jend)/2)
      
! GFS does not have deep convective cloud top and bottom fields
      HTOPD=SPVAL
      HBOTD=SPVAL   
      HTOPS=SPVAL
      HBOTS=SPVAL 
      CUPPT=SPVAL 

!
!!!! DONE GETTING
! Will derive isobaric OMEGA from continuity equation later. 
!      OMGA=SPVAL
! retrieve d3d fields if it's listed
      print*,'iostatus for d3d file= ',iostatusD3D
      if(iostatusD3D==0)then ! start reading d3d file
! retrieve longwave tendency using getgb
        Index=41
        VarName=avbl(index)
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=is(index)
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,rlwtt(1,jsta_2l,ll))
        end do
		
! retrieve shortwave tendency using getgb
        Index=40
        VarName=avbl(index)
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=is(index)
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,rswtt(1,jsta_2l,ll))
        end do	
        
! retrieve vertical diffusion tendency using getgb
        Index=356
        VarName='VDIFF TNDY'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,vdifftt(1,jsta_2l,ll))
        end do	
	
! retrieve deep convective tendency using getgb
        Index=79
        VarName=avbl(index)
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=is(index)
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,tcucn(1,jsta_2l,ll))
        end do
	
! retrieve shallow convective tendency using getgb
        Index=358
        VarName='S CNVCT TNDY'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,tcucns(1,jsta_2l,ll))
        end do
	
! retrieve grid scale latent heat tendency using getgb
        Index=78
        VarName=avbl(index)
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=is(index)
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,train(1,jsta_2l,ll))
        end do
	
! retrieve vertical diffusion moistening using getgb
        Index=360
        VarName='Vertical diffusion moistening'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,vdiffmois(1,jsta_2l,ll))
        end do					

! retrieve deep convection moistening using getgb
        Index=361
        VarName='deep convection moistening'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,dconvmois(1,jsta_2l,ll))
        end do
	
! retrieve shallow convection moistening using getgb
        Index=362
        VarName='shallow convection moistening'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,sconvmois(1,jsta_2l,ll))
        end do	
	
! retrieve non-radiation tendency using getgb
        Index=363
        VarName='non-radiation tendency'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,nradtt(1,jsta_2l,ll))
        end do
	
! retrieve Vertical diffusion of ozone using getgb
        Index=364
        VarName='Vertical diffusion of ozone'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,o3vdiff(1,jsta_2l,ll))
        end do
	
! retrieve ozone production using getgb
        Index=365
        VarName='Ozone production'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,o3prod(1,jsta_2l,ll))
        end do
	
! retrieve ozone tendency using getgb
        Index=366
        VarName='Ozone tendency'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
	do l=1,lm 
	 jpds(7)=l
	 ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l       & 
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName                &
           ,jpds,jgds,kpds,o3tndy(1,jsta_2l,ll))
        end do	
	
! retrieve mass weighted PV using getgb
        Index=367
        VarName='Mass weighted PV'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l  &
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName           &
           ,jpds,jgds,kpds,mwpv(1,jsta_2l,ll))
        end do

! retrieve OZONE TNDY using getgb
        Index=368
        VarName='?'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l   &
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName            &
           ,jpds,jgds,kpds,unknown(1,jsta_2l,ll))
        end do

! retrieve vertical diffusion zonal acceleration
        Index=369
        VarName='VDIFF Z ACCE'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l  &
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName           &
           ,jpds,jgds,kpds,vdiffzacce(1,jsta_2l,ll))
        end do

! retrieve gravity drag zonal acceleration
        Index=370
        VarName='G DRAG Z ACCE'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l   &
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName            &
           ,jpds,jgds,kpds,zgdrag(1,jsta_2l,ll))
        end do

! retrieve convective U momemtum mixing
        Index=371
        VarName='CNVCT U M MIX'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l    &
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName             &
           ,jpds,jgds,kpds,cnvctummixing(1,jsta_2l,ll))
        end do

! retrieve vertical diffusion meridional acceleration
        Index=372
        VarName='VDIFF M ACCE'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l    &
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName             &
           ,jpds,jgds,kpds,vdiffmacce(1,jsta_2l,ll))
        end do

! retrieve gravity drag meridional acceleration
        Index=373
        VarName='G DRAG M ACCE'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l    &
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName             &
           ,jpds,jgds,kpds,mgdrag(1,jsta_2l,ll))
        end do

! retrieve convective V momemtum mixing
        Index=374
        VarName='CNVCT V M MIX'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l    &
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName             &
           ,jpds,jgds,kpds,cnvctvmmixing(1,jsta_2l,ll))
        end do

! retrieve nonconvective cloud fraction
        Index=375
        VarName='N CNVCT CLD FRA'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l    &
           ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName             &
           ,jpds,jgds,kpds,ncnvctcfrac(1,jsta_2l,ll))
        end do

! retrieve convective upward mass flux
        Index=391
        VarName='CNVCT U M FLX'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l &
          ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName          &
          ,jpds,jgds,kpds,cnvctumflx(1,jsta_2l,ll))
        end do
	
! retrieve convective downward mass flux
        Index=392
        VarName='CNVCT D M FLX'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l &
          ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName          &
          ,jpds,jgds,kpds,cnvctdmflx(1,jsta_2l,ll))
        end do	
	     
! retrieve nonconvective detraintment flux
        Index=393
        VarName='CNVCT DET M FLX'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l &
          ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName          &
          ,jpds,jgds,kpds,cnvctdetmflx(1,jsta_2l,ll))
        end do	     

! retrieve cnvct gravity drag zonal acceleration
        Index=394
        VarName='CNVCT G DRAG Z ACCE'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l &
          ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName          &
          ,jpds,jgds,kpds,cnvctzgdrag(1,jsta_2l,ll))
        end do

! retrieve cnvct gravity drag meridional acceleration
        Index=395
        VarName='CNVCT G DRAG M ACCE'
        jpds=-1.0
        jgds=-1.0
        jpds(5)=iq(index)
        jpds(6)=109
        do l=1,lm
         jpds(7)=l
         ll=lm-l+1 !flip 3d fields to count from top down
         call getgbandscatter(me,iunitd3d,im,jm,im_jm,jsta,jsta_2l &
          ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName          &
          ,jpds,jgds,kpds,cnvctmgdrag(1,jsta_2l,ll))
        end do
     
        call baclose(iunitd3d,status)
	print*,'done reading D3D fields'            
      end if ! end of d3d file read
      print *,'after d3d files reading,mype=',me

! Retrieve aer fields if it's listed (GOCART)
      print *, 'iostatus for aer file=', iostatusAER
      if(iostatusAER==0)then ! start reading aer file

! retrieve dust emission fluxes
      do K = 1, nbin_du
       if ( K == 1) VarName='DUEM001'
       if ( K == 2) VarName='DUEM002'
       if ( K == 3) VarName='DUEM003'
       if ( K == 4) VarName='DUEM004'
       if ( K == 5) VarName='DUEM005'
       VcoordName='atmos col'
       l=1
       call getnemsandscatter(me,rfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,duem(1,jsta_2l,K))
      if(debugprint)print*,'sample ',VarName,' = ',duem(im/2,(jsta+jend)/2,k)
      enddo

! retrieve dust sedimentation fluxes
      do K = 1, nbin_du
       if ( K == 1) VarName='DUSD001'
       if ( K == 2) VarName='DUSD002'
       if ( K == 3) VarName='DUSD003'
       if ( K == 4) VarName='DUSD004'
       if ( K == 5) VarName='DUSD005'
       VcoordName='atmos col'
       l=1
       call getnemsandscatter(me,rfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,dusd(1,jsta_2l,K))
      if(debugprint)print*,'sample ',VarName,' = ',dusd(im/2,(jsta+jend)/2,k)
      enddo

! retrieve dust dry deposition fluxes
      do K = 1, nbin_du
       if ( K == 1) VarName='DUDP001'
       if ( K == 2) VarName='DUDP002'
       if ( K == 3) VarName='DUDP003'
       if ( K == 4) VarName='DUDP004'
       if ( K == 5) VarName='DUDP005'
       VcoordName='atmos col'
       l=1
       call getnemsandscatter(me,rfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,dudp(1,jsta_2l,K))
        print *,'dudp,ck=',maxval(dudp(1:im,jsta:jend,k)), &
             minval(dudp(1:im,jsta:jend,k))
      if(debugprint)print*,'sample ',VarName,' = ',dudp(im/2,(jsta+jend)/2,k)
      enddo

! retrieve dust wet deposition fluxes
      do K = 1, nbin_du
       if ( K == 1) VarName='DUWT001'
       if ( K == 2) VarName='DUWT002'
       if ( K == 3) VarName='DUWT003'
       if ( K == 4) VarName='DUWT004'
       if ( K == 5) VarName='DUWT005'
       VcoordName='atmos col'
       l=1
       call getnemsandscatter(me,rfile,im,jm,jsta,jsta_2l &
       ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
       ,l,im,jm,nframe,duwt(1,jsta_2l,K))
      if(debugprint)print*,'sample ',VarName,' = ',duwt(im/2,(jsta+jend)/2,k)
      enddo

! retrieve sfc mass concentration
      VarName='DUSMASS'
      VcoordName='atmos col'
      l=1
      call getnemsandscatter(me,rfile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,dusmass)
      if(debugprint)print*,'sample ',VarName,' = ',dusmass(im/2,(jsta+jend)/2)

! retrieve col mass density
      VarName='DUCMASS'
      VcoordName='atmos col'
      l=1
      call getnemsandscatter(me,rfile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,ducmass)
      if(debugprint)print*,'sample ',VarName,' = ',ducmass(im/2,(jsta+jend)/2)

! retrieve sfc mass concentration (pm2.5)
      VarName='DUSMASS25'
      VcoordName='atmos col'
      l=1
      call getnemsandscatter(me,rfile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,dusmass25)
      if(debugprint)print*,'sample ',VarName,' = ',dusmass25(im/2,(jsta+jend)/2)

! retrieve col mass density (pm2.5)
      VarName='DUCMASS25'
      VcoordName='atmos col'
      l=1
      call getnemsandscatter(me,rfile,im,jm,jsta,jsta_2l &
      ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,VcoordName &
      ,l,im,jm,nframe,ducmass25)
      if(debugprint)print*,'sample ',VarName,' = ',ducmass25(im/2,(jsta+jend)/2)

      end if ! end of aer file read
      print *,'after aer files reading,mype=',me

! pos east
       call collect_loc(gdlat,dummy)
       if(me.eq.0)then
        latstart=nint(dummy(1,1)*gdsdegr)
        latlast=nint(dummy(im,jm)*gdsdegr)
	print*,'laststart,latlast B bcast= ',latstart,latlast,'gdsdegr=',gdsdegr,&
          'dummy(1,1)=',dummy(1,1),dummy(im,jm),'gdlat=',gdlat(1,1)
       end if
       call mpi_bcast(latstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(latlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*) 'laststart,latlast,me A calling bcast=',latstart,latlast,me
       call collect_loc(gdlon,dummy)
       if(me.eq.0)then
        lonstart=nint(dummy(1,1)*gdsdegr)
        lonlast=nint(dummy(im,jm)*gdsdegr)
       end if
       call mpi_bcast(lonstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(lonlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*)'lonstart,lonlast A calling bcast=',lonstart,lonlast
!
!        ncdump -h
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

! close up shop
!      call ext_int_ioclose ( DataHandle, Status )

! generate look up table for lifted parcel calculations

      THL=210.
      PLQ=70000.

      CALL TABLE(PTBL,TTBL,PT,                                     &  
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
!       TPREC=float(ifhr)
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
        !  Note: The calculation of the map scale factor at the standard
        !        lat/lon and the PSMAPF
        ! Get map factor at 60 degrees (N or S) for PS projection, which will
        ! be needed to correctly define the DX and DY values in the GRIB GDS
          if (TRUELAT1 .LT. 0.) THEN
            LAT = -60.
          else
            LAT = 60.
          end if

          CALL MSFPS (LAT,TRUELAT1*0.001,PSMAPF)

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
! close all files
        call nemsio_close(nfile,iret=status)
	call nemsio_close(ffile,iret=status)
	call nemsio_close(rfile,iret=status)
!	call baclose(iunit,status)
	
      RETURN
      END


