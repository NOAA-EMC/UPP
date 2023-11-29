!> @file
!> @brief wrfpost() drives the external wrf post processor.
!> @return wrfpost 
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-24 | Russ Treadon              | Coded etapost as stand alone code
!> 1998-05-29 | Black                     | Conversion of post code from 1-D to 2-D
!> 1900-02-04 | Jim Tuccillo              | Parallel version via MPI
!> 2001-02-15 | Jim Tuccillo              | Many common blocks replaced with modules to support fortran "allocate"s for the exact size of the arrays needed based on the number of mpi tasks. This was done to reduce the address space that the loader sees. These changes were necessary for running larger domains such as 12 kms
!> 2001-06-15 | JIM Tuccillo              | Added asyncronous I/O capability. if there are more than one mpi task, the io will be done aynchronously by the last MPI task.
!> 2002-06-17 | Mike Baldwin              | Convert etapost to wrfpost. Include wrf I/O api for input of model data. Modify code to deal with C-grid data. Streamline output to a call of one subroutine instead of three. Replace common blocks with a limited number of modules.
!> 2004-01-01 | H Chuang                  | Added nmm io module and binary options
!> 2005-07-08 | Binbin Zhou               | Added RSM model
!> 2005-12-05 | H Chuang                  | Added capability to output off-hour forecast which has no impacts on on-hour forecast
!> 2006-02-20 | Chuang, Black, and Rogers | Finalized complete list of NAM operational products from WRF
!> 2006-02-27 | H Chuang                  | Modified to post multiple forecast hours in one execution
!> 2006-03-03 | H Chuang                  | Added parrish's mpi binary io to read binary WRF file as random asscess so that variables in WRF output don't have to be read in in specific order 
!> 2011-02-06 | J Wang                    | Add grib2 option
!> 2011-12-14 | Sarah Lu                  | Add the option to read ngac aer file 
!> 2012-01-28 | J WANG                    | Use post available fields in xml file for grib2
!> 2013-06-25 | S Moorthi                 | Add gocart_on logical option to save memory
!> 2013-10-03 | J Wang                    |Add option for po to be pascal, and add gocart_on,d3d_on and popascal to namelist
!> 2020-03-25 | J Meng                    | Remove grib1
!> 2021-06-20 | W Meng                    | Remove reading grib1 and gfsio lib
!> 2021-07-07 | J MENG                    |2D DECOMPOSITION
!> 2021-10-22 | KaYee Wong                | Created formal fortran namelist for itag
!> 2021-11-03 | Tracy Hertneky            | Removed SIGIO option
!> 2022-01-14 | W Meng                    | Remove interfaces INITPOST_GS_NEMS, INITPOST_NEMS_MPIIO, INITPOST_NMM and INITPOST_GFS_NETCDF
!> 2022-03-15 | W Meng                    | Unify FV3 based interfaces
!> 2022-09-22 | L Zhang                   | Add option of nasa_on to process ufs-aerosols
!> 2022-11-08 | K Wang                    | Replace aqfamaq_on with aqf_on
!> 2023-01-24 | Sam Trahan                | write_ifi_debug_files flag for IFI debug capability
!> 2023-03-21 | Jesse Meng                | Add slrutah_on option to use U Utah SLR
!> 2023-04-04 | Li(Kate Zhang)  |Add namelist optoin for CCPP-Chem (UFS-Chem) 
!         and 2D diag. output (d2d_chem) for GEFS-Aerosols and CCPP-Chem model.
!> 2023-05-20 | Rahul Mahajan             | Bug fix for fileNameFlat as namelist configurable
!> 2023-08-16 | Yali Mao                  | Add gtg_on logical option
!> 2023-11-29 | Eric James                | Add method_blsn logical option
!> @author Mike Bladwin NSSL/SPC @date 2002-06-18
!---------------------------------------------------------------------
!> @return wrfpost
!---------------------------------------------------------------------
      PROGRAM WRFPOST

!
!
!============================================================================================================
!
!     This is an MPI code. All array indexing is with respect to the global indices. Loop indices 
!     look as follows for N MPI tasks.
!
!
!
!  Original                                            New
!  Index                                               Index
!
!   JM ----------------------------------------------- JEND
! JM-1 -                                             - JEND_M
! JM-2 -               MPI TASK N-1                  - JEND_M2
!      -                                             -
!      -                                             -
!      ----------------------------------------------- JSTA, JSTA_M, JSTA_M2
!      ----------------------------------------------- JEND, JEND_M, JEND_M2
!      -                                             -
!      -               MPI TASK N-2                  -
!      -                                             -
!      -                                             -
!      ----------------------------------------------- JSTA, JSTA_M, JSTA_M2
!
!                           .
!                           .
!                           .
!
!      ----------------------------------------------- JEND, JEND_M, JEND_M2
!      -                                             -
!      -               MPI TASK 1                    -
!      -                                             -
!      -                                             -
!      ----------------------------------------------- JSTA, JSTA_M, JSTA_M2
!      ----------------------------------------------- JEND, JEND_M, JEND_M2
!      -                                             - 
!      -               MPI TASK 0                    - 
!    3 -                                             - JSTA_M2
!    2 -                                             - JSTA_M
!    1 ----------------------------------------------- JSTA
!
!     1                                              IM               
!
!
!     Jim Tuccillo
!     Jan 2000
!
!     README - Jim Tuccillo Feb 2001
! 
!     Many common blocks have been replaced by modules to support Fortran
!     "allocate" commands. Many of the 3-D arrays are now allocated to be the
!     exact size required based on the number of MPI tasks. The dimensioning will be 
!        x ( im,jsta_2l:jend_2u,lm)
!     Most 2-D arrays continue to be dimensioned (im,jm). This is fine but please be aware 
!     that the EXCH routine for arrays dimensioned (im,jm) is different than arrays dimensioned
!     (im,jsta_2l:jend_2u). Also, be careful about passing any arrays dimensioned
!     (im,jst_2l:jend_2u,lm). See examples in the code as to the correct calling sequence and
!     EXCH routine to use.
!
!
!     ASYNCHRONOUS I/O HAS BEEN ADDED. THE LAST MPI TASK DOES THE I/O. IF THERE IS
!     ONLY ONE MPI TASK THN TASK ) DOES THE I/O.
!     THE CODE HAS GOTTEN A LITTLE KLUDGY. BASICLY, IM, IMX and IMOUT MUST BE EQUAL
!     AND REPRESENT THE VALUE USED IN THE MODEL. THE SAME HOLDS FOR JM, JMX and JMOUT.
!
!     Jim Tuccillo June 2001
!
!
!===========================================================================================
!
      use netcdf
      use nemsio_module, only: nemsio_getheadvar, nemsio_gfile, nemsio_init, nemsio_open, &
                               nemsio_getfilehead,nemsio_close
      use CTLBLK_mod,    only: filenameaer, me, num_procs, num_servers, mpi_comm_comp, datestr,      &
              mpi_comm_inter, filename, ioform, grib, idat, filenameflux, filenamed3d, gdsdegr,      &
              spldef, modelname, ihrst, lsmdef,vtimeunits, tprec, pthresh, datahandle, im, jm, lm,   &
              lp1, lm1, im_jm, isf_surface_physics, nsoil, spl, lsmp1, global, imp_physics,          &
              ista, iend, ista_m, iend_m, ista_2l, iend_2u,                                          &
              jsta, jend, jsta_m, jend_m, jsta_2l, jend_2u, novegtype, icount_calmict, npset, datapd,&
              lsm, fld_info, etafld2_tim, eta2p_tim, mdl2sigma_tim, cldrad_tim, miscln_tim,          &
              mdl2agl_tim, mdl2std_tim, mdl2thandpv_tim, calrad_wcloud_tim,nasa_on,gccpp_on,         &
              fixed_tim, time_output, imin, surfce2_tim, komax, ivegsrc, d3d_on, gocart_on,rdaod,    &
              readxml_tim, spval, fullmodelname, submodelname, hyb_sigp, filenameflat, aqf_on,numx,  &
              run_ifi_tim, slrutah_on, d2d_chem, gtg_on, method_blsn
      use grib2_module,   only: gribit2,num_pset,nrecout,first_grbtbl,grib_info_finalize
      use upp_ifi_mod, only: write_ifi_debug_files
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      implicit none
!
      type(nemsio_gfile) :: nfile,ffile,rfile
      INCLUDE "mpif.h"
!
!     DECLARE VARIABLES.
!     
!     SET HEADER WRITER FLAGS TO TRUE.
!
!temporary vars
!
      real(kind=8) :: time_initpost=0.,INITPOST_tim=0.,btim,bbtim
      real            rinc(5), untcnvt
      integer      :: status=0,iostatusD3D=0,iostatusFlux=0
      integer i,j,iii,l,k,ierr,nrec,ist,lusig,idrt,ncid3d,ncid2d,varid
      integer      :: PRNTSEC,iim,jjm,llm,ioutcount,itmp,iret,iunit,        &
                      iunitd3d,iyear,imn,iday,LCNTRL,ieof
      integer      :: iostatusAER
      logical      :: popascal
!
      integer      :: kpo,kth,kpv
      real,dimension(komax) :: po,th,pv
      namelist/nampgb/kpo,po,kth,th,kpv,pv,fileNameAER,d3d_on,gocart_on,gccpp_on, nasa_on,gtg_on,method_blsn,popascal &
                     ,hyb_sigp,rdaod,d2d_chem, aqf_on,slrutah_on, vtimeunits,numx,write_ifi_debug_files
      integer      :: itag_ierr
      namelist/model_inputs/fileName,IOFORM,grib,DateStr,MODELNAME,SUBMODELNAME &
                     ,fileNameFlux,fileNameFlat

      character startdate*19,SysDepInfo*80,IOWRFNAME*3,post_fname*255
      character cgar*1,cdum*4,line*10
!
!------------------------------------------------------------------------------
!     START HERE
!
      call start()
!
!     INITIALIZE MPI
      
      CALL SETUP_SERVERS(ME,                         &
     &                   NUM_PROCS,                  &
     &                   NUM_SERVERS,                &
     &                   MPI_COMM_COMP,              &
     &                   MPI_COMM_INTER)
!
!     ME IS THE RANK
!     NUM_PROCS IS THE NUMBER OF TASKS DOING POSTING
!     NUM_SERVERS IS ONE IF THERE ARE MORE THAN ONE TOTAL MPI TASKS, OTHERWISE ZERO
!     MPI_COMM_COMP IS THE INTRACOMMUNICATOR
!     MPI_COMM_INTER IS THE INTERCOMMUNICATOR FOR COMMUNCATION BETWEEN TASK 0 OF THE
!        TASKS DOING THE POSTING AND THE I/O SERVER
!
!
!     IF WE HAVE MORE THAN 1 MPI TASK THEN WE WILL FIRE UP THE IO SERVER
!     THE LAST TASK ( IN THE CONTEXT OF MPI_COMM_WORLD ) IS THE I/O SERVER
!
      if (me == 0) CALL W3TAGB('nems     ',0000,0000,0000,'np23   ')

      if ( me >= num_procs ) then
!
         call server
!
      else
        spval = 9.9e10
!
!**************************************************************************
!KaYee: Read itag in Fortran Namelist format
!Set default 
       SUBMODELNAME='NONE'
!Set control file name
       fileNameFlat='postxconfig-NT.txt'
       numx=1
!open namelist
       open(5,file='itag')
       read(5,nml=model_inputs,iostat=itag_ierr,err=888)
888    if (itag_ierr /= 0) then
       print*,'Incorrect namelist variable(s) found in the itag file,stopping.'
       stop
       endif
       if (me == 0) write(6, model_inputs)
       
!       if(MODELNAME == 'NMM')then
!        read(5,1114) VTIMEUNITS
! 1114   format(a4)
!        if (me==0) print*,'VALID TIME UNITS = ', VTIMEUNITS
!       endif
!
 303  format('MODELNAME="',A,'" SUBMODELNAME="',A,'"')

        if(me==0) write(*,*)'MODELNAME: ', MODELNAME, SUBMODELNAME

! assume for now that the first date in the stdin file is the start date
        read(DateStr,300) iyear,imn,iday,ihrst,imin
        if (me==0) write(*,*) 'in WRFPOST iyear,imn,iday,ihrst,imin',                &
                    iyear,imn,iday,ihrst,imin
 300    format(i4,1x,i2,1x,i2,1x,i2,1x,i2)

        IDAT(1) = imn
        IDAT(2) = iday
        IDAT(3) = iyear
        IDAT(4) = ihrst
        IDAT(5) = imin

 111    format(a256)
 112    format(a19)
 113    format(a20)
 114    format(a8)
 120    format(a5)
 121    format(a4)

!KaYee: Read in GFS/FV3 runs in Fortran Namelist Format.
      if(grib=='grib2') then
        gdsdegr = 1.d6
      endif
! 
! set default for kpo, kth, th, kpv, pv     
        kpo = 0
        po  = 0
        kth = 6
        th  = (/310.,320.,350.,450.,550.,650.,(0.,k=kth+1,komax)/) ! isentropic level to output
        kpv = 8
        pv  = (/0.5,-0.5,1.0,-1.0,1.5,-1.5,2.0,-2.0,(0.,k=kpv+1,komax)/)

        hyb_sigp    = .true.
        d3d_on      = .false.
        gocart_on   = .false.
        gccpp_on    = .false.
        nasa_on     = .false.
        aqf_on      = .false.
        slrutah_on  = .false.
        gtg_on   = .false.
        method_blsn = .true.
        popascal    = .false.
        fileNameAER = ''
        rdaod       = .false.
        d2d_chem     = .false.
        vtimeunits  = ''

        read(5,nampgb,iostat=iret,end=119)
 119    continue
        if (me == 0) write(6, nampgb)
       if(mod(num_procs,numx)/=0) then
         if (me==0) then
           print*,'total proces, num_procs=', num_procs 
           print*,'number of subdomain in x direction, numx=', numx 
           print*,'remainder of num_procs/numx = ', mod(num_procs,numx)
           print*,'Warning!!! the remainder of num_procs/numx is not 0, reset numx=1 &
     &             in this run or you adjust numx in the itag file to restart'
         endif
!        stop 9999
         numx=1
         if(me == 0) print*,'Warning!!!  Reset numx as 1, numx=',numx
       endif
       if(numx>num_procs/2) then
         if (me==0) then
           print*,'total proces, num_procs=', num_procs
           print*,'number of subdomain in x direction, numx=', numx
           print*,'Warning!!! numx cannot exceed num_procs/2, reset numx=1 in this run'
           print*,'or you adjust numx in the itag file to restart'
         endif
         numx=1
         if(me == 0) print*,'Warning!!!  Reset numx as 1, numx=',numx
       endif     
        if(me == 0) then
          print*,'komax,iret for nampgb= ',komax,iret 
          print*,'komax,kpo,kth,th,kpv,pv,fileNameAER,nasa_on,popascal= ',komax,kpo        &
     &           ,kth,th(1:kth),kpv,pv(1:kpv),trim(fileNameAER),nasa_on,popascal
          print*,'NUM_PROCS=',NUM_PROCS
          print*,'numx= ',numx
        endif

        IF(TRIM(IOFORM) /= 'netcdfpara' .AND. TRIM(IOFORM) /= 'netcdf' ) THEN
          numx=1
          if(me == 0) print*,'2D decomposition only supports netcdfpara IO.'
          if(me == 0) print*,'Reset numx= ',numx
        ENDIF

        IF(MODELNAME /= 'FV3R' .AND. MODELNAME /= 'GFS') THEN
          numx=1
          if(me == 0) print*,'2D decomposition only supports GFS and FV3R.'
          if(me == 0) print*,'Reset numx= ',numx
        ENDIF

! set up pressure level from POSTGPVARS or DEFAULT
        if(kpo == 0) then
! use default pressure levels
          if(me == 0) then
            print*,'using default pressure levels,spldef=',(spldef(l),l=1,lsmdef)
          endif
          lsm = lsmdef
          do l=1,lsm
            spl(l) = spldef(l)
          end do
        else
! use POSTGPVARS
          if(me == 0) then
            print*,'using pressure levels from POSTGPVARS'
          endif
          lsm = kpo
          if( .not. popascal ) then
            untcnvt = 100.
          else
            untcnvt = 1.
          endif
          if(po(lsm) < po(1))then ! post logic assumes asscending
            do l=1,lsm
              spl(l) = po(lsm-l+1)*untcnvt 
            end do
          else
            do l=1,lsm
              spl(l) = po(l)*untcnvt
            end do
          end if
        end if
        LSMP1 = LSM+1
      
 116    continue

! set PTHRESH for different models
        if(MODELNAME == 'NMM')then
          PTHRESH = 0.000004
        else
          PTHRESH = 0.000001
        end if  
!Chuang: add dynamical allocation
        if(TRIM(IOFORM) == 'netcdf' .OR. TRIM(IOFORM) == 'netcdfpara') THEN
         IF(MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR' .OR. MODELNAME == 'NMM') THEN
          call ext_ncd_ioinit(SysDepInfo,Status)
          call ext_ncd_open_for_read( trim(fileName), 0, 0, " ",          &
            DataHandle, Status)
          if ( Status /= 0 ) then
            print*,'error opening ',fileName, ' Status = ', Status ; stop
          endif
          call ext_ncd_get_dom_ti_integer(DataHandle                      &
            ,'WEST-EAST_GRID_DIMENSION',iim,1,ioutcount, status )
          im = iim-1
          call ext_ncd_get_dom_ti_integer(DataHandle                      &
             ,'SOUTH-NORTH_GRID_DIMENSION',jjm,1,ioutcount, status )
          jm = jjm-1
          call ext_ncd_get_dom_ti_integer(DataHandle                      &
            ,'BOTTOM-TOP_GRID_DIMENSION',llm,1,ioutcount, status )
          lm    = llm-1
          LP1   = LM+1
          LM1   = LM-1
          IM_JM = IM*JM
       
! Read and set global value for surface physics scheme
          call ext_ncd_get_dom_ti_integer(DataHandle                      &
            ,'SF_SURFACE_PHYSICS',itmp,1,ioutcount, status )
          iSF_SURFACE_PHYSICS = itmp
! set NSOIL to 4 as default for NOAH but change if using other
! SFC scheme
          NSOIL = 4
          IF(itmp      == 1) then !thermal diffusion scheme
            NSOIL = 5
          ELSE IF(itmp == 3) then ! RUC LSM
            NSOIL = 9
          ELSE IF(itmp == 7) then ! Pleim Xu
            NSOIL = 2
          END IF

          call ext_ncd_ioclose ( DataHandle, Status )
         ELSE
! use parallel netcdf lib directly to read FV3 output in netCDF
          spval = 9.99e20
          Status = nf90_open(trim(fileName),IOR(NF90_NOWRITE,NF90_MPIIO), &
                   ncid3d,comm=mpi_comm_world,info=mpi_info_null)
          if ( Status /= 0 ) then
            print*,'error opening ',fileName, ' Status = ', Status 
            stop
          endif
          Status = nf90_open(trim(fileNameFlux),IOR(NF90_NOWRITE,NF90_MPIIO), &
                   ncid2d,comm=mpi_comm_world,info=mpi_info_null)
          if ( Status /= 0 ) then
            print*,'error opening ',fileNameFlux, ' Status = ', Status
            stop
          endif
! read in LSM index and nsoil here
          Status=nf90_get_att(ncid2d,nf90_global,'landsfcmdl', iSF_SURFACE_PHYSICS)
          if(Status/=0)then
            print*,'landsfcmdl not found; assigning to 2'
            iSF_SURFACE_PHYSICS=2 !set LSM physics to 2 for NOAH
          endif
          if(iSF_SURFACE_PHYSICS<2)then
            iSF_SURFACE_PHYSICS=2 !set LSM physics to 2 for NOAH
          endif
          Status=nf90_get_att(ncid2d,nf90_global,'nsoil', NSOIL)
          if(Status/=0)then
            print*,'nsoil not found; assigning to 4'
            NSOIL=4 !set nsoil to 4 for NOAH
          endif
! read imp_physics
          Status=nf90_get_att(ncid2d,nf90_global,'imp_physics',imp_physics)
          if(Status/=0)then
            print*,'imp_physics not found; assigning to GFDL 11'
            imp_physics=11
          endif
! get dimesions
          Status = nf90_inq_dimid(ncid3d,'grid_xt',varid)
          if ( Status /= 0 ) then
           print*,Status,varid
           STOP 1
          end if
          Status = nf90_inquire_dimension(ncid3d,varid,len=im)
          if ( Status /= 0 ) then
           print*,Status
           STOP 1      
          end if   
          Status = nf90_inq_dimid(ncid3d,'grid_yt',varid)
          if ( Status /= 0 ) then
           print*,Status,varid
           STOP 1
          end if
          Status = nf90_inquire_dimension(ncid3d,varid,len=jm)
          if ( Status /= 0 ) then
           print*,Status
           STOP 1
          end if
          Status = nf90_inq_dimid(ncid3d,'pfull',varid)
          if ( Status /= 0 ) then
           print*,Status,varid
           STOP 1
          end if
          Status = nf90_inquire_dimension(ncid3d,varid,len=lm)
          if ( Status /= 0 ) then
           print*,Status
           STOP 1
          end if
          LP1   = LM+1
          LM1   = LM-1
          IM_JM = IM*JM
! set NSOIL to 4 as default for NOAH but change if using other
! SFC scheme
!          NSOIL = 4
         END IF 

        ELSE IF(TRIM(IOFORM) == 'binary'       .OR.                       &
                TRIM(IOFORM) == 'binarympiio' ) THEN
          print*,'WRF Binary format is no longer supported'
          STOP 9996
! NEMSIO format
        ELSE IF(TRIM(IOFORM) == 'binarynemsio' .or.                        &
          TRIM(IOFORM) == 'binarynemsiompiio' )THEN
      
          spval = 9.99e20
          IF(ME == 0)THEN
            call nemsio_init(iret=status)
            call nemsio_open(nfile,trim(filename),'read',iret=status)
            if ( Status /= 0 ) then
              print*,'error opening ',fileName, ' Status = ', Status ; stop
            endif
!---
            call nemsio_getfilehead(nfile,iret=status,nrec=nrec            &
                ,dimx=im,dimy=jm,dimz=lm,nsoil=nsoil)
            if ( Status /= 0 ) then
              print*,'error finding model dimensions '; stop
            endif
            call nemsio_getheadvar(nfile,'global',global,iret)
            if (iret /= 0)then 
              print*,"global not found in file-Assigned false"
              global = .FALSE.
            end if
            IF(MODELNAME == 'GFS') global = .TRUE.
! global NMMB has i=1 overlap with i=im so post will cut off i=im	     
            if(global .and. MODELNAME == 'NMM') im = im-1

          END IF

          CALL mpi_bcast(im,   1,MPI_INTEGER,0, mpi_comm_comp,status) 
          call mpi_bcast(jm,   1,MPI_INTEGER,0, mpi_comm_comp,status)
          call mpi_bcast(lm,   1,MPI_INTEGER,0, mpi_comm_comp,status)
          call mpi_bcast(nsoil,1,MPI_INTEGER,0, mpi_comm_comp,status)

          call mpi_bcast(global,1,MPI_LOGICAL,0,mpi_comm_comp,status)
          LP1   = LM+1
          LM1   = LM-1
          IM_JM = IM*JM

! opening GFS flux file
          IF(MODELNAME == 'GFS') THEN
!	    iunit=33
            call nemsio_open(ffile,trim(fileNameFlux),'read',iret=iostatusFlux)
            if ( iostatusFlux /= 0 ) then
              print*,'error opening ',fileNameFlux, ' Status = ', iostatusFlux
            endif
            iostatusD3D = -1
            iunitd3d    = -1
!
! opening GFS aer file
            call nemsio_open(rfile,trim(fileNameAER),'read',iret=iostatusAER)
            if ( iostatusAER /= 0  .and.  me == 0) then
              print*,'error opening AER ',fileNameAER, ' Status = ', iostatusAER
            endif
!
!           print*,'iostatusD3D in WRFPOST= ',iostatusD3D

          END IF 

        ELSE
          PRINT*,'UNKNOWN MODEL OUTPUT FORMAT, STOPPING'
          STOP 9999
        END IF  


        CALL MPI_FIRST()
        CALL ALLOCATE_ALL()
     
!
!       INITIALIZE POST COMMON BLOCKS 
!
        LCNTRL = 14
        REWIND(LCNTRL)

! EXP. initialize netcdf here instead
        bbtim = mpi_wtime()
        btim = mpi_wtime()
! set default novegtype
        if(MODELNAME == 'GFS')THEN
          novegtype = 13 
          ivegsrc   = 2
        else if(MODELNAME=='NMM' .and. TRIM(IOFORM)=='binarynemsio')then
          novegtype = 20
          ivegsrc   = 1
        else if(MODELNAME=='RAPR')then
          novegtype = 20
          ivegsrc   = 1
        else ! USGS
          novegtype = 24
          ivegsrc   = 0
        end if
      
! Reading model output for different models and IO format     
 
        IF(TRIM(IOFORM) == 'netcdf' .OR. TRIM(IOFORM) == 'netcdfpara') THEN
          IF(MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR') THEN
            CALL INITPOST
          ELSE IF (MODELNAME == 'FV3R' .OR. MODELNAME == 'GFS') THEN
! use parallel netcdf library to read output directly
            CALL INITPOST_NETCDF(ncid2d,ncid3d)
          ELSE
            PRINT*,'POST does not have netcdf option for model,',MODELNAME,' STOPPING,'
            STOP 9998
          END IF
        ELSE IF(TRIM(IOFORM) == 'binarympiio') THEN 
          IF(MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR' .OR. MODELNAME == 'NMM') THEN
            print*,'WRF BINARY IO FORMAT IS NO LONGER SUPPORTED, STOPPING'
            STOP 9996
          ELSE IF(MODELNAME == 'RSM') THEN                            
            print*,'MPI BINARY IO IS NOT YET INSTALLED FOR RSM, STOPPING'
            STOP 9997
          ELSE
            PRINT*,'POST does not have mpiio option for this model, STOPPING'
            STOP 9998
          END IF
        ELSE IF(TRIM(IOFORM) == 'binarynemsio') THEN 
          IF(MODELNAME == 'NMM') THEN
            CALL INITPOST_NEMS(NREC,nfile)
          ELSE
            PRINT*,'POST does not have nemsio option for model,',MODELNAME,' STOPPING,'
            STOP 9998

          END IF
       
        ELSE IF(TRIM(IOFORM) == 'binarynemsiompiio')THEN
          IF(MODELNAME == 'GFS') THEN
! close nemsio file for serial read
            call nemsio_close(nfile,iret=status)
            call nemsio_close(ffile,iret=status)
            call nemsio_close(rfile,iret=status)
            CALL INITPOST_GFS_NEMS_MPIIO(iostatusAER)
          ELSE
            PRINT*,'POST does not have nemsio mpi option for model,',MODELNAME, &
            'STOPPING,'
            STOP 9999

          END IF 

        ELSE
          PRINT*,'UNKNOWN MODEL OUTPUT FORMAT, STOPPING'
          STOP 9999
        END IF 
        INITPOST_tim  = INITPOST_tim +(mpi_wtime() - btim)
        IF(ME == 0)THEN
          WRITE(6,*)'WRFPOST:  INITIALIZED POST COMMON BLOCKS'
        ENDIF
!
!       IF GRIB2 read out post aviable fields xml file and post control file
!
        if(grib == "grib2") then
          btim=mpi_wtime()
          call READ_xml()
          READxml_tim = READxml_tim + (mpi_wtime() - btim)
        endif
! 
!     LOOP OVER THE OUTPUT GRID(S).  FIELD(S) AND  OUTPUT GRID(S) ARE SPECIFIED
!     IN THE CONTROL FILE.  WE PROCESS ONE GRID AND ITS FIELDS AT A TIME.
!     THAT'S WHAT THIS LOOP DOES.
!     
        icount_calmict = 0
        first_grbtbl   = .true.
        npset          = 0
!10   CONTINUE
!     
!        READ CONTROL FILE DIRECTING WHICH FIELDS ON WHICH
!        LEVELS AND TO WHICH GRID TO INTERPOLATE DATA TO.
!        VARIABLE IEOF/=0 WHEN THERE ARE NO MORE GRIDS TO PROCESS.
!
!                      --------    grib1 processing  ---------------
!                                 ------------------
!        if (grib == "grib1") then !DO NOT REVERT TO GRIB1. GRIB1 NOT SUPPORTED ANYMORE
!          IEOF = 0
!          do while (ieof == 0)
!            CALL READCNTRL(kth,IEOF)
!            IF(ME == 0)THEN
!              WRITE(6,*)'POST:  RETURN FROM READCNTRL.  ', 'IEOF=',IEOF
!            ENDIF
!
!           PROCESS SELECTED FIELDS.  FOR EACH SELECTED FIELD/LEVEL
!           WE GO THROUGH THE FOLLOWING STEPS:
!             (1) COMPUTE FIELD IF NEED BE
!             (2) WRITE FIELD TO OUTPUT FILE IN GRIB.
!
!            if (ieof == 0) then
!              CALL PROCESS(kth,kpv,th(1:kth),pv(1:kpv),iostatusD3D)
!              IF(ME == 0)THEN
!                WRITE(6,*)' '
!                WRITE(6,*)'WRFPOST:  PREPARE TO PROCESS NEXT GRID'
!              ENDIF
!            endif
!
!           PROCESS NEXT GRID.
!
!          enddo
!                      --------    grib2 processing  ---------------
!                                 ------------------
!        elseif (grib == "grib2") then
        if (me==0) write(*,*) ' in WRFPOST OUTFORM= ',grib
        if (me==0) write(*,*) '  GRIB1 IS NOT SUPPORTED ANYMORE'    
        if (grib == "grib2") then
          do while (npset < num_pset)
            npset = npset+1
            if (me==0) write(*,*)' in WRFPOST npset=',npset,' num_pset=',num_pset
            CALL SET_OUTFLDS(kth,th,kpv,pv)
            if (me==0) write(*,*)' in WRFPOST size datapd',size(datapd) 
            if(allocated(datapd)) deallocate(datapd)
!Jesse x-decomposition
!           allocate(datapd(im,1:jend-jsta+1,nrecout+100))
            allocate(datapd(1:iend-ista+1,1:jend-jsta+1,nrecout+100))
!$omp parallel do private(i,j,k)
            do k=1,nrecout+100
              do j=1,jend+1-jsta
!Jesse x-decomposition
!               do i=1,im
                do i =1,iend+1-ista
                  datapd(i,j,k) = 0.
                enddo
              enddo
            enddo
            call get_postfilename(post_fname)
            if (me==0) write(*,*)'post_fname=',trim(post_fname)
            if (me==0) write(*,*)'get_postfilename,post_fname=',trim(post_fname), &
                      'npset=',npset, 'num_pset=',num_pset,            &
                      'iSF_SURFACE_PHYSICS=',iSF_SURFACE_PHYSICS
!     
!           PROCESS SELECTED FIELDS.  FOR EACH SELECTED FIELD/LEVEL
!           WE GO THROUGH THE FOLLOWING STEPS:
!             (1) COMPUTE FIELD IF NEED BE
!             (2) WRITE FIELD TO OUTPUT FILE IN GRIB.
!
            CALL PROCESS(kth,kpv,th(1:kth),pv(1:kpv),iostatusD3D)
            IF(ME == 0) WRITE(6,*)'WRFPOST:  PREPARE TO PROCESS NEXT GRID'
!
!           write(*,*)'enter gribit2 before mpi_barrier'
            call mpi_barrier(mpi_comm_comp,ierr)

!           if(me==0)call w3tage('bf grb2  ')
!           write(*,*)'enter gribit2 after mpi barrier'
            call gribit2(post_fname)
            deallocate(datapd)
            deallocate(fld_info)
!
!           PROCESS NEXT GRID.
!
          enddo

        endif
!     
!-------
        call grib_info_finalize()
!
        IF(ME == 0) THEN
          WRITE(6,*)' '
          WRITE(6,*)'ALL GRIDS PROCESSED.'
          WRITE(6,*)' '
        ENDIF
!
        call DE_ALLOCATE

!       GO TO 98
 1000   CONTINUE
!exp      call ext_ncd_ioclose ( DataHandle, Status )
!
        IF(ME == 0) THEN
         print*, 'INITPOST_tim = ',  INITPOST_tim
         print*, 'MDLFLD_tim = ',  ETAFLD2_tim
         print*, 'MDL2P_tim =  ',ETA2P_tim 
         print*, 'MDL2SIGMA_tim =  ',MDL2SIGMA_tim 
         print*, 'MDL2AGL_tim =  ',MDL2AGL_tim 
         print*, 'SURFCE_tim =  ',SURFCE2_tim
         print*, 'CLDRAD_tim =  ',CLDRAD_tim 
         print*, 'MISCLN_tim = ',MISCLN_tim
         print*, 'MDL2STD_tim =  ',MDL2STD_tim    
         print*, 'FIXED_tim = ',FIXED_tim
         print*, 'MDL2THANDPV_tim =  ',MDL2THANDPV_tim
         print*, 'CALRAD_WCLOUD_tim = ',CALRAD_WCLOUD_tim    
         print*, 'RUN_IFI_tim = ',RUN_IFI_tim
         print*, 'Total time = ',(mpi_wtime() - bbtim)
         print*, 'Time for OUTPUT = ',time_output
         print*, 'Time for READxml = ',READxml_tim
        endif
!     
!       END OF PROGRAM.
!
!
!       MPI_LAST WILL SHUTDOWN THE IO SERVER, IF IT EXISTS
!
        CALL MPI_LAST
!
!
      end if
!
!
!
!      call summary()
      if (me == 0) CALL W3TAGE('UNIFIED_POST')
      CALL MPI_FINALIZE(IERR)


      STOP 0

      END

