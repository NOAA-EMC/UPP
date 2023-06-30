!> @file
!> @brief initpost_netcdf() initializes post for run.
!>
!> @author Hui-Ya Chuang @date 2016-03-04
!>
!> This routine initializes constants and
!> variables at the start of GFS model or post
!> processor run.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2007-03-01 | Hui-Ya Chuang | Initial. Start from INITPOST_GFS_NEMS_MPIIO.f
!> 2021-03-11 | Bo Cui        | Change local arrays to dimension (im,jsta:jend)
!> 2021-10-26 | Jesse Meng    | 2D DECOMPOSITION
!> 2022-02-07 | Wen Meng      | Changes for parallel netcdf read
!> 2022-03-15 | Wen Meng      | Unify regional and global interfaces
!> 2022-03-22 | Wen Meng      | Read PWAT from model
!> 2022-04-08 | Bo Cui        | 2D decomposition for unified fv3 read interfaces
!> 2022-06-05 | Hui-Ya Chuang | Modify dx/dy computation for RRFS domain over north pole
!> 2022-07-10 | Wen Meng      | Output lat/lon on four coner points of rotated lat-lon grids in text file.
!> 2022-07-18 | Wen Meng      | Read instant top of atmos ULWRF from model
!> 2022-09-18 | Li(Kate) Zhang| Add aerosol fileds for GEFS-Aerosols (gocart_on) and UFS-Aerosols(nasa_on) model
!> 2022-10-28 | Eric James    | Modifications to allow passing through soil moisture availability field from RUC LSM for RRFS
!> 2022-11-08 | Kai Wang      | Read time averaged PM2.5 and O3 concentration from model
!> 2022-11-08 | Wen Meng      | Remove instant PM2.5 calculation
!> 2022-11-16 | Eric James    | Read smoke, dust, biomass burning, and hourly wildfire potential from RRFS
!> 2022-12-07 | Wen Meng      | Read AOD from AQM model
!> 2022-12-23 | Eric Aligo    | Read six winter weather diagnostics from model
!> 2023-01-30 | Sam Trahan    | Read cldfra or cldfra_bl, whichever is available
!> 2023-02-23 | Eric James    | Read coarse PM and aodtot from RRFS
!> 2023-03-02 | Sam Trahan    | Read lightning threat index fields
!> 2023-03-22 | WM Lewis      | Read RRFS effective radii (EFFRL, EFFRI, EFFRS)
!> 2023-04-04 |Li(Kate Zhang)  |Add namelist optoin for CCPP-Chem(UFS-Chem) 
!         and 2D diag. output (d2d_chem) for GEFS-Aerosols and CCPP-Chem model.
!> 2023-04-17 | Eric James    | Read in unified ext550 extinction (and remove aodtot) for RRFS
!> 2023-04-21 | Eric James    | Read in / calculate some fields needed for GSL p-type diagnosis for RRFS
!> 2023-05-31 | Wen Meng      | Bug fix in qrmax initialization
!> 2023-06-14 | Wen Meng      ! Bug fix of reading seaswtc and modification of sndepac calculation
!>
!> @author Hui-Ya Chuang @date 2016-03-04
!----------------------------------------------------------------------
!> @brief INITPOST_NETCDF() This routine initializes constants and
!> variables at the start of GFS model or post processor run. 
!> 
!> @param[in] ncid2d integer netCDF ID of physics model output file.
!> @param[in] ncid3d integer netCDF ID of dynamics model output file.
!----------------------------------------------------------------------
      SUBROUTINE INITPOST_NETCDF(ncid2d,ncid3d)


      use netcdf
      use vrbls4d, only: dust, SALT, SUSO, SOOT, WASO, smoke, fv3dust, coarsepm,                &
              no3,nh4, PP25, PP10 
      use vrbls3d, only: t, q, uh, vh, pmid, pint, alpint, dpres, zint, zmid, o3,               &
              qqr, qqs, cwm, qqi, qqw, omga, rhomid, q2, cfr, rlwtt, rswtt, tcucn,              &
              tcucns, train, el_pbl, exch_h, vdifftt, vdiffmois, dconvmois, nradtt,             &
              o3vdiff, o3prod, o3tndy, mwpv, unknown, vdiffzacce, zgdrag,cnvctummixing,         &
              vdiffmacce, mgdrag, cnvctvmmixing, ncnvctcfrac, cnvctumflx, cnvctdmflx,           &
              cnvctzgdrag, sconvmois, cnvctmgdrag, cnvctdetmflx, duwt, duem, dusd, dudp,        &
              dusv,ssem,sssd,ssdp,sswt,sssv,bcem,bcsd,bcdp,bcwt,bcsv,ocem,ocsd,ocdp,ocwt,ocsv, &
              wh, qqg, ref_10cm, qqnifa, qqnwfa, avgpmtf, avgozcon, aextc55, taod5503d,         &
              effri, effrl, effrs

      use vrbls2d, only: f, pd, fis, pblh, ustar, z0, ths, qs, twbs, qwbs, avgcprate,           &
              cprate, avgprec, prec, lspa, sno, sndepac, si, cldefi, th10, q10, tshltr, pshltr, &
              tshltr, albase, avgalbedo, avgtcdc, czen, czmean, mxsnal, landfrac, radot, sigt4, &
              cfrach, cfracl, cfracm, avgcfrach, qshltr, avgcfracl, avgcfracm, cnvcfr,          &
              islope, cmc, grnflx, vegfrc, acfrcv, ncfrcv, acfrst, ncfrst, ssroff,              &
              bgroff, rlwin, rlwtoa, cldwork, alwin, alwout, alwtoa, rswin, rswinc,             &
              rswout, aswin, auvbin, auvbinc, aswout, aswtoa, sfcshx, sfclhx, subshx,           &
              snopcx, sfcux, sfcvx, sfcuxi, sfcvxi, sfcuvx, gtaux, gtauy, potevp, u10, v10, smstav,&
              smstot, ivgtyp, isltyp, sfcevp, sfcexc, acsnow, acsnom, sst, thz0, qz0,           &
              uz0, vz0, ptop, htop, pbot, hbot, ptopl, pbotl, ttopl, ptopm, pbotm, ttopm,       &
              ptoph, pboth, pblcfr, ttoph, runoff, tecan, tetran, tedir, twa, maxtshltr,        &
              mintshltr, maxrhshltr, fdnsst, acgraup, graup_bucket, acfrain, frzrn_bucket,      &
              snow_acm, snow_bkt, snownc, graupelnc, qrmax,                                     &
              minrhshltr, dzice, smcwlt, suntime, fieldcapa, htopd, hbotd, htops, hbots,        &
              cuppt, dusmass, ducmass, dusmass25, ducmass25, aswintoa,rel_vort_maxhy1,          &
              maxqshltr, minqshltr, acond, sr, u10h, v10h,refd_max, w_up_max, w_dn_max,         &
              up_heli_max,up_heli_min,up_heli_max03,up_heli_min03,rel_vort_max01,u10max, v10max,  &
              avgedir,avgecan,paha,pahi,avgetrans,avgesnow,avgprec_cont,avgcprate_cont,rel_vort_max, &
              avisbeamswin,avisdiffswin,airbeamswin,airdiffswin,refdm10c_max,wspd10max, &
              alwoutc,alwtoac,aswoutc,aswtoac,alwinc,aswinc,avgpotevp,snoavg, &
              ti,aod550,du_aod550,ss_aod550,su_aod550,oc_aod550,bc_aod550,prate_max,maod,dustpm10, &
              dustcb,bccb,occb,sulfcb,sscb,dustallcb,ssallcb,dustpm,sspm,pp25cb,pp10cb,no3cb,nh4cb,&
              pwat, ebb, hwp, aqm_aod550, ltg1_max,ltg2_max,ltg3_max
      use soil,  only: sldpth, sllevel, sh2o, smc, stc
      use masks, only: lmv, lmh, htm, vtm, gdlat, gdlon, dx, dy, hbm2, sm, sice
      use physcons_post, only: grav => con_g, fv => con_fvirt, rgas => con_rd,                     &
                            eps => con_eps, epsm1 => con_epsm1
      use params_mod, only: erad, dtr, tfrz, h1, d608, rd, p1000, capa,pi
      use lookup_mod, only: thl, plq, ptbl, ttbl, rdq, rdth, rdp, rdthe, pl, qs0, sqs, sthe,    &
                            ttblq, rdpq, rdtheq, stheq, the0q, the0
      use ctlblk_mod, only: me, mpi_comm_comp, icnt, idsp, jsta, jend, ihrst, idat, sdat, ifhr, &
              ifmin, filename, tprec, tclod, trdlw, trdsw, tsrfc, tmaxmin, td3d, restrt, sdat,  &
              jend_m, imin, imp_physics, dt, spval, pdtop, pt, qmin, nbin_du, nphs, dtq2, ardlw,&
              ardsw, asrfc, avrain, avcnvc, theat, gdsdegr, spl, lsm, alsl, im, jm, im_jm, lm,  &
              jsta_2l, jend_2u, nsoil, lp1, icu_physics, ivegsrc, novegtype, nbin_ss, nbin_bc,  &
              nbin_oc, nbin_su, nbin_no3, nbin_nh4, gocart_on,gccpp_on, nasa_on,pt_tbl,hyb_sigp,&
              filenameFlux, fileNameAER,                                               &
              iSF_SURFACE_PHYSICS,rdaod, d2d_chem, modelname, aqf_on,                         &
              ista, iend, ista_2l, iend_2u,iend_m
      use gridspec_mod, only: maptype, gridtype, latstart, latlast, lonstart, lonlast, cenlon,  &
              dxval, dyval, truelat2, truelat1, psmapf, cenlat,lonstartv, lonlastv, cenlonv,    &
              latstartv, latlastv,cenlatv,latstart_r,latlast_r,lonstart_r,lonlast_r, STANDLON,  &
              latse,lonse,latnw,lonnw
      use upp_physics, only: fpvsnew
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      implicit none
!
!      type(nemsio_gfile) :: nfile,ffile,rfile
      integer,parameter          :: nvar2d=48
!      character(nemsio_charkind) :: name2d(nvar2d)
      integer                    :: nvar3d, numDims
!      character(nemsio_charkind), allocatable :: name3din(:), name3dout(:)
!      character(nemsio_charkind) :: varname,levtype
!
!     INCLUDE/SET PARAMETERS.
!     
      INCLUDE "mpif.h"

!     integer,parameter:: MAXPTS=1000000 ! max im*jm points
!
!     real,parameter:: con_g       =9.80665e+0! gravity
!     real,parameter:: con_rv      =4.6150e+2 ! gas constant H2O
!     real,parameter:: con_rd      =2.8705e+2 ! gas constant air
!     real,parameter:: con_fvirt   =con_rv/con_rd-1.
!     real,parameter:: con_eps     =con_rd/con_rv
!     real,parameter:: con_epsm1   =con_rd/con_rv-1
!
! This version of INITPOST shows how to initialize, open, read from, and
! close a NetCDF dataset. In order to change it to read an internal (binary)
! dataset, do a global replacement of _ncd_ with _int_. 

      real, parameter    :: gravi = 1.0/grav
      character(len=20)  :: VarName, VcoordName
      integer            :: Status, fldsize, fldst, recn, recn_vvel
      character             startdate*19,SysDepInfo*80,cgar*1
      character             startdate2(19)*4, flatlon*40
      logical            :: read_lonlat=.true.
! 
!     NOTE: SOME INTEGER VARIABLES ARE READ INTO DUMMY ( A REAL ). THIS IS OK
!     AS LONG AS REALS AND INTEGERS ARE THE SAME SIZE.
!
!     ALSO, EXTRACT IS CALLED WITH DUMMY ( A REAL ) EVEN WHEN THE NUMBERS ARE
!     INTEGERS - THIS IS OK AS LONG AS INTEGERS AND REALS ARE THE SAME SIZE.
      LOGICAL RUNB,SINGLRST,SUBPOST,NEST,HYDRO,IOOMG,IOALL
      logical, parameter :: debugprint = .false., zerout = .false.
!     logical, parameter :: debugprint = .true.,  zerout = .false.
      logical :: convert_rad_to_deg=.false.
      CHARACTER*32 varcharval 
!      CHARACTER*40 CONTRL,FILALL,FILMST,FILTMP,FILTKE,FILUNV,FILCLD,FILRAD,FILSFC
      CHARACTER*4  RESTHR
      CHARACTER    FNAME*255,ENVAR*50
      INTEGER      IDATE(8),JDATE(8),JPDS(200),JGDS(200),KPDS(200),KGDS(200)
!     LOGICAL*1    LB(IM,JM)
!     
!     INCLUDE COMMON BLOCKS.
!
!     DECLARE VARIABLES.
!     
!      REAL fhour
      integer nfhour ! forecast hour from nems io file
      integer fhzero !bucket
      real dtp !physics time step
      real dz
      REAL RINC(5)

      REAL DUMMY(IM,JM)
!jw
      integer ii,jj,js,je,iyear,imn,iday,itmp,ioutcount,istatus,       &
              I,J,L,ll,k,kf,irtn,igdout,n,Index,nframe,                &
              nframed2,iunitd3d,ierr,idum,iret,nrec,idrt
      integer ncid3d,ncid2d,varid,nhcas,varid_bl,iret_bl
      real    TSTART,TLMH,TSPH,ES,FACT,soilayert,soilayerb,zhour,dum,  &
              tvll,pmll,tv, tx1, tx2

      character*20,allocatable :: recname(:)
      integer,     allocatable :: reclev(:), kmsk(:,:)
      real,        allocatable :: glat1d(:), glon1d(:), qstl(:)
      real,        allocatable :: wrk1(:,:), wrk2(:,:)
      real,        allocatable :: p2d(:,:),  t2d(:,:),  q2d(:,:),      &
                                  qs2d(:,:), cw2d(:,:), cfr2d(:,:)
      real, dimension(lm+1)    :: ak5, bk5
      real*8, allocatable      :: pm2d(:,:), pi2d(:,:)
      real,   allocatable      :: tmp(:)
      real                     :: buf(ista_2l:iend_2u,jsta_2l:jend_2u)
      real                     :: buf2(ista_2l:iend_2u,jsta_2l:jend_2u)
      real                     :: buf3d(ista_2l:iend_2u,jsta_2l:jend_2u,lm)
      real                     :: chem_2d(ista_2l:iend_2u,jsta_2l:jend_2u)
      real                     :: chemT(ista_2l:iend_2u,jsta_2l:jend_2u,lm)
      real                     :: dt1(ista_2l:iend_2u,jsta_2l:jend_2u,lm)
      real                     :: dt2(ista_2l:iend_2u,jsta_2l:jend_2u,lm)
      real                     :: dt3(ista_2l:iend_2u,jsta_2l:jend_2u,lm)
      real                     :: dt4(ista_2l:iend_2u,jsta_2l:jend_2u,lm)
      real                     :: dt5(ista_2l:iend_2u,jsta_2l:jend_2u,lm)

!     real buf(ista_2l:iend_2u,jsta_2l:jend_2u),bufsoil(im,nsoil,jsta_2l:jend_2u)   &
!         ,buf3d(ista_2l:iend_2u,jsta_2l:jend_2u,lm),buf3d2(im,lp1,jsta_2l:jend_2u)

      real LAT
      integer isa, jsa, latghf, jtem, idvc, idsl, nvcoord, ip1, nn, npass

      integer, parameter    :: npass2=5, npass3=30
      real, parameter       :: third=1.0/3.0
      INTEGER, DIMENSION(2) :: ij4min, ij4max
      REAL                  :: omgmin, omgmax
      real, allocatable :: d2d(:,:), u2d(:,:), v2d(:,:), omga2d(:,:)
      REAL, ALLOCATABLE :: ps2d(:,:),psx2d(:,:),psy2d(:,:)
      real, allocatable :: div3d(:,:,:)
      real(kind=4),allocatable :: vcrd(:,:)
      real                     :: dum_const 
      real, allocatable :: ext550(:,:,:)

      if (modelname == 'FV3R') then
         allocate(ext550(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      endif

!***********************************************************************
!     START INIT HERE.
!
      if(me==0)then
      WRITE(6,*)'INITPOST:  ENTER INITPOST_NETCDF'
      WRITE(6,*)'me=',me,  &
           'jsta_2l=',jsta_2l,'jend_2u=', &
           jend_2u,'im=',im, &
           'ista_2l=',ista_2l,'iend_2u=',iend_2u, &
           'ista=',ista,'iend=',iend, &
           'iend_m=',iend_m
      endif
!     
      isa = (ista+iend) / 2
      jsa = (jsta+jend) / 2

!$omp parallel do private(i,j)
      do j = jsta_2l, jend_2u
        do i= ista_2l, iend_2u
          buf(i,j) = spval
        enddo
      enddo

      Status=nf90_get_att(ncid3d,nf90_global,'ak',ak5)
      if(Status/=0)then
       print*,'ak not found; assigning missing value'
       ak5=spval
      else
       if(me==0)print*,'ak5= ',ak5
      end if 
      Status=nf90_get_att(ncid3d,nf90_global,'idrt',idrt)
      if(Status/=0)then
       if(me==0)print*,'idrt not in netcdf file,reading grid'
       Status=nf90_get_att(ncid3d,nf90_global,'grid',varcharval)
       if(Status/=0)then
        if(me==0)print*,'idrt and grid not in netcdf file, set default to latlon'
        idrt=0
        MAPTYPE=0
       else
        if(trim(varcharval)=='rotated_latlon')then
         MAPTYPE=207
         idrt=207
         Status=nf90_get_att(ncid3d,nf90_global,'cen_lon',dum_const)
         if(Status/=0)then
          print*,'cen_lon not found; assigning missing value'
          cenlon=spval
         else
          if(dum_const<0.)then
           cenlon=nint((dum_const+360.)*gdsdegr)
          else
           cenlon=dum_const*gdsdegr
          end if 
         end if
         Status=nf90_get_att(ncid3d,nf90_global,'cen_lat',dum_const)
         if(Status/=0)then
          print*,'cen_lat not found; assigning missing value'
          cenlat=spval
         else
          cenlat=dum_const*gdsdegr
         end if

         Status=nf90_get_att(ncid3d,nf90_global,'lon1',dum_const)
         if(Status/=0)then
          print*,'lonstart_r not found; assigning missing value'
          lonstart_r=spval
         else
          if(dum_const<0.)then
           lonstart_r=nint((dum_const+360.)*gdsdegr)
          else
           lonstart_r=dum_const*gdsdegr
          end if
         end if
         Status=nf90_get_att(ncid3d,nf90_global,'lat1',dum_const)
         if(Status/=0)then
          print*,'latstart_r not found; assigning missing value'
          latstart_r=spval
         else
          latstart_r=dum_const*gdsdegr
         end if

         Status=nf90_get_att(ncid3d,nf90_global,'lon2',dum_const)
         if(Status/=0)then
          print*,'lonlast_r not found; assigning missing value'
          lonlast_r=spval
         else
          if(dum_const<0.)then
           lonlast_r=nint((dum_const+360.)*gdsdegr)
          else
           lonlast_r=dum_const*gdsdegr
          end if
         end if
         Status=nf90_get_att(ncid3d,nf90_global,'lat2',dum_const)
         if(Status/=0)then
          print*,'latlast_r not found; assigning missing value'
          latlast_r=spval
         else
          latlast_r=dum_const*gdsdegr
         end if

         Status=nf90_get_att(ncid3d,nf90_global,'dlon',dum_const)
         if(Status/=0)then
          print*,'dlmd not found; assigning missing value'
          dxval=spval
         else
          dxval=dum_const*gdsdegr
         end if
         Status=nf90_get_att(ncid3d,nf90_global,'dlat',dum_const)
         if(Status/=0)then
          print*,'dphd not found; assigning missing value'
          dyval=spval
         else
          dyval=dum_const*gdsdegr
         end if

!         print*,'lonstart,latstart,cenlon,cenlat,dyval,dxval', &
!         lonstart,latstart,cenlon,cenlat,dyval,dxval

! Jili Dong add support for regular lat lon (2019/03/22) start
        else if(trim(varcharval)=='latlon')then
         MAPTYPE=0
         idrt=0

         Status=nf90_get_att(ncid3d,nf90_global,'lon1',dum_const)
         if(Status/=0)then
          print*,'lonstart not found; assigning missing value'
          lonstart=spval
         else
          if(dum_const<0.)then
           lonstart=nint((dum_const+360.)*gdsdegr)
          else
           lonstart=dum_const*gdsdegr
          end if
         end if
         Status=nf90_get_att(ncid3d,nf90_global,'lat1',dum_const)
         if(Status/=0)then
          print*,'latstart not found; assigning missing value'
          latstart=spval
         else
          latstart=dum_const*gdsdegr
         end if

         Status=nf90_get_att(ncid3d,nf90_global,'lon2',dum_const)
         if(Status/=0)then
          print*,'lonlast not found; assigning missing value'
          lonlast=spval
         else
          if(dum_const<0.)then
           lonlast=nint((dum_const+360.)*gdsdegr)
          else
           lonlast=dum_const*gdsdegr
          end if
         end if
         Status=nf90_get_att(ncid3d,nf90_global,'lat2',dum_const)
         if(Status/=0)then
          print*,'latlast not found; assigning missing value'
          latlast=spval
         else
          latlast=dum_const*gdsdegr
         end if

         Status=nf90_get_att(ncid3d,nf90_global,'dlon',dum_const)
         if(Status/=0)then
          print*,'dlmd not found; assigning missing value'
          dxval=spval
         else
          dxval=dum_const*gdsdegr
         end if
         Status=nf90_get_att(ncid3d,nf90_global,'dlat',dum_const)
         if(Status/=0)then
          print*,'dphd not found; assigning missing value'
          dyval=spval
         else
          dyval=dum_const*gdsdegr
         end if

!         print*,'lonstart,latstart,dyval,dxval', &
!         lonstart,lonlast,latstart,latlast,dyval,dxval

! Jili Dong add support for regular lat lon (2019/03/22) end 
 
        ELSE IF (trim(varcharval)=='lambert_conformal')then

          MAPTYPE=1
          idrt=1
          Status=nf90_get_att(ncid3d,nf90_global,'cen_lon',dum_const)
          if(Status/=0)then
            print*,'cen_lon not found; assigning missing value'
            cenlon=spval
          else
            if(dum_const<0.)then
              cenlon=nint((dum_const+360.)*gdsdegr)
            else
              cenlon=dum_const*gdsdegr
            end if
          end if
          Status=nf90_get_att(ncid3d,nf90_global,'cen_lat',dum_const)
          if(Status/=0)then
            print*,'cen_lat not found; assigning missing value'
            cenlat=spval
          else
            cenlat=dum_const*gdsdegr
          end if

          Status=nf90_get_att(ncid3d,nf90_global,'lon1',dum_const)
          if(Status/=0)then
            print*,'lonstart not found; assigning missing value'
            lonstart=spval
          else
            if(dum_const<0.)then
              lonstart=nint((dum_const+360.)*gdsdegr)
            else
              lonstart=dum_const*gdsdegr
            end if
          end if
          Status=nf90_get_att(ncid3d,nf90_global,'lat1',dum_const)
          if(Status/=0)then
            print*,'latstart not found; assigning missing value'
            latstart=spval
          else
            latstart=dum_const*gdsdegr
          end if

          Status=nf90_get_att(ncid3d,nf90_global,'stdlat1',dum_const)
          if(Status/=0)then
            print*,'stdlat1 not found; assigning missing value'
            truelat1=spval
          else
            truelat1=dum_const*gdsdegr
          end if
          Status=nf90_get_att(ncid3d,nf90_global,'stdlat2',dum_const)
          if(Status/=0)then
            print*,'stdlat2 not found; assigning missing value'
            truelat2=spval
          else
            truelat2=dum_const*gdsdegr
          end if

          Status=nf90_get_att(ncid3d,nf90_global,'dx',dum_const)
          if(Status/=0)then
            print*,'dx not found; assigning missing value'
            dxval=spval
          else
            dxval=dum_const*1.E3
          end if
          Status=nf90_get_att(ncid3d,nf90_global,'dy',dum_const)
          if(Status/=0)then
            print*,'dphd not found; assigning missing value'
            dyval=spval
          else
            dyval=dum_const*1.E3
          end if

          STANDLON = cenlon
!          print*,'lonstart,latstart,cenlon,cenlat,truelat1,truelat2, &
!                  stadlon,dyval,dxval', &
!          lonstart,latstart,cenlon,cenlat,truelat1,truelat2,standlon,dyval,dxval

        else if(trim(varcharval)=='gaussian')then
         MAPTYPE=4
         idrt=4
        else ! setting default maptype
         MAPTYPE=0
         idrt=0
        end if
       end if
      end if
      if(me==0)print*,'idrt MAPTYPE= ',idrt,MAPTYPE
!     STEP 1.  READ MODEL OUTPUT FILE
!
!
!***
!
! LMH and LMV  always = LM for sigma-type vert coord

!$omp parallel do private(i,j)
      do j = jsta_2l, jend_2u
        do i = ista_2l, iend_2u
          LMV(i,j) = lm
          LMH(i,j) = lm
        end do
      end do

! HTM VTM all 1 for sigma-type vert coord

!$omp parallel do private(i,j,l)
      do l = 1, lm
        do j = jsta_2l, jend_2u
          do i = ista_2l, iend_2u
            HTM (i,j,l) = 1.0
            VTM (i,j,l) = 1.0
          end do
        end do
      end do

      Status=nf90_get_att(ncid3d,nf90_global,'nhcas',nhcas)
      if(Status/=0)then
      if(me==0) print*,'nhcas not in netcdf file, set default to nonhydro'
       nhcas=0
      end if
      if(me==0)print*,'nhcas= ',nhcas
      if (nhcas == 0 ) then  !non-hydrostatic case
       nrec=18
       allocate (recname(nrec))
       recname=[character(len=20) :: 'ugrd','vgrd','spfh','tmp','o3mr', &
                                     'presnh','dzdt', 'clwmr','dpres',  &
                                     'delz','icmr','rwmr',              &
                                     'snmr','grle','smoke','dust',      &
                                     'coarsepm','ext550']
      else
       nrec=8
       allocate (recname(nrec))
       recname=[character(len=20) :: 'ugrd','vgrd','tmp','spfh','o3mr', &
                                     'hypres', 'clwmr','dpres']
      endif

!     write(*,*)'nrec=',nrec
      !allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
      allocate(glat1d(jm),glon1d(im))

! hardwire idate for now
!      idate=(/2017,08,07,00,0,0,0,0/)
! get cycle start time
      Status=nf90_inq_varid(ncid3d,'time',varid)
      if(Status/=0)then
       print*,'time not in netcdf file, stopping'
       STOP 1
      else
       Status=nf90_get_att(ncid3d,varid,'units',varcharval)
       if(Status/=0)then
         if(me==0)print*,'time unit not available'
       else
         if(me==0)print*,'time unit read from netcdf file= ',varcharval
! assume use hours as unit
!       idate_loc=index(varcharval,'since')+6
         read(varcharval,101)idate(1),idate(2),idate(3),idate(4),idate(5)
       end if
!       Status=nf90_inquire_dimension(ncid3d,varid,len=ntimes)
!       allocate(fhours(ntimes))
!       status = nf90_inq_varid(ncid3d,varid,fhours)
!       Status=nf90_get_var(ncid3d,varid,nfhour,start=(/1/), &
!              count=(/1/))
!       if(Status/=0)then
!        print*,'forecast hour not in netcdf file, stopping'
!        STOP 1
!       end if
      end if
 101  format(T13,i4,1x,i2,1x,i2,1x,i2,1x,i2)
      !print*,'idate= ',idate(1:5)

! Jili Dong check output format for coordinate reading
      Status=nf90_inq_varid(ncid3d,'grid_xt',varid)
      Status=nf90_inquire_variable(ncid3d,varid,ndims = numDims)
      if(numDims==1.and.modelname=="FV3R") then
        read_lonlat=.true.
      else
        read_lonlat=.false.
      end if
      

! Jili Dong add support for new write component output
! get longitude 
      if (read_lonlat) then
        Status=nf90_inq_varid(ncid3d,'lon',varid)
        Status=nf90_inquire_variable(ncid3d,varid,ndims = numDims)
        if(debugprint)print*,'number of dim for gdlon ',numDims
      else
        Status=nf90_inq_varid(ncid3d,'grid_xt',varid)
        Status=nf90_inquire_variable(ncid3d,varid,ndims = numDims)
        if(debugprint)print*,'number of dim for gdlon ',numDims
      end if
      if(numDims==1)then
        Status=nf90_get_var(ncid3d,varid,glon1d)  
        do j=jsta,jend
          do i=ista,iend
            gdlon(i,j) = real(glon1d(i),kind=4)
          end do
        end do
        lonstart = nint(glon1d(1)*gdsdegr)
        lonlast  = nint(glon1d(im)*gdsdegr)

! Jili Dong add support for regular lat lon (2019/03/22) start
       if (MAPTYPE == 0) then
        if(lonstart<0.)then
         lonstart=lonstart+360.*gdsdegr
        end if
        if(lonlast<0.)then
         lonlast=lonlast+360.*gdsdegr
        end if
       end if
! Jili Dong add support for regular lat lon (2019/03/22) end 

      else if(numDims==2)then
        Status=nf90_get_var(ncid3d,varid,dummy)
        if(maxval(abs(dummy))<2.0*pi)convert_rad_to_deg=.true. 
        if(convert_rad_to_deg)then
         do j=jsta,jend
          do i=ista,iend
            gdlon(i,j) = real(dummy(i,j),kind=4)*180./pi
          end do
         end do
        else
         do j=jsta,jend
          do i=ista,iend
            gdlon(i,j) = real(dummy(i,j),kind=4)
          end do
         end do
        end if
        if(convert_rad_to_deg)then
         lonstart = nint(dummy(1,1)*gdsdegr)*180./pi
         lonlast  = nint(dummy(im,jm)*gdsdegr)*180./pi
         lonse    = nint(dummy(im,1)*gdsdegr)*180./pi
         lonnw    = nint(dummy(1,jm)*gdsdegr)*180./pi
        else
         lonstart = nint(dummy(1,1)*gdsdegr)
         lonlast  = nint(dummy(im,jm)*gdsdegr)
         lonse    = nint(dummy(im,1)*gdsdegr)
         lonnw    = nint(dummy(1,jm)*gdsdegr)
        end if

! Jili Dong add support for regular lat lon (2019/03/22) start
       if (MAPTYPE == 0) then
        if(lonstart<0.)then
         lonstart=lonstart+360.*gdsdegr
        end if
        if(lonlast<0.)then
         lonlast=lonlast+360.*gdsdegr
        end if
       end if
! Jili Dong add support for regular lat lon (2019/03/22) end 

      end if
!      print*,'lonstart,lonlast ',lonstart,lonlast 
! Jili Dong add support for new write component output
! get latitude
      if (read_lonlat) then
        Status=nf90_inq_varid(ncid3d,'lat',varid)
        Status=nf90_inquire_variable(ncid3d,varid,ndims = numDims)
        if(debugprint)print*,'number of dim for gdlat ',numDims
      else
        Status=nf90_inq_varid(ncid3d,'grid_yt',varid)
        Status=nf90_inquire_variable(ncid3d,varid,ndims = numDims)
        if(debugprint)print*,'number of dim for gdlat ',numDims
      end if
      if(numDims==1)then
        Status=nf90_get_var(ncid3d,varid,glat1d)
        do j=jsta,jend
          do i=ista,iend
            gdlat(i,j) = real(glat1d(j),kind=4)
          end do
        end do
        latstart = nint(glat1d(1)*gdsdegr)
        latlast  = nint(glat1d(jm)*gdsdegr)
      else if(numDims==2)then
        Status=nf90_get_var(ncid3d,varid,dummy)
        if(maxval(abs(dummy))<pi)convert_rad_to_deg=.true.
        if(convert_rad_to_deg)then
         do j=jsta,jend
          do i=ista,iend
            gdlat(i,j) = real(dummy(i,j),kind=4)*180./pi
          end do
         end do
        else
         do j=jsta,jend
          do i=ista,iend
            gdlat(i,j) = real(dummy(i,j),kind=4)
          end do
         end do
        end if
        if(convert_rad_to_deg)then
         latstart = nint(dummy(1,1)*gdsdegr)*180./pi
         latlast  = nint(dummy(im,jm)*gdsdegr)*180./pi
         latse    = nint(dummy(im,1)*gdsdegr)*180./pi
         latnw    = nint(dummy(1,jm)*gdsdegr)*180./pi
        else
         latstart = nint(dummy(1,1)*gdsdegr)
         latlast  = nint(dummy(im,jm)*gdsdegr)
         latse    = nint(dummy(im,1)*gdsdegr)
         latnw    = nint(dummy(1,jm)*gdsdegr)
        end if
      end if
      !print*,'laststart,latlast = ',latstart,latlast
      if(debugprint)print*,'me sample gdlon gdlat= ' &
     ,me,gdlon(isa,jsa),gdlat(isa,jsa)

! Specify grid staggering type
      gridtype = 'A'
      maptype=idrt
      if (me == 0) print *, 'maptype and gridtype is ', &
      maptype,gridtype
 
      if(gridtype == 'A')then
        lonstartv=lonstart
        lonlastv=lonlast
        latstartv=latstart
        latlastv=latlast
        cenlatv=cenlat
        cenlonv=cenlon
      end if

      if(debugprint)then
        if (me == 0)then
          do i=1,nrec
            print *,'recname=',trim(recname(i))
          end do
        end if
      end if


      deallocate(glat1d,glon1d)

!      print*,'idate = ',(idate(i),i=1,7)
!      print*,'nfhour = ',nfhour
      
! sample print point
      ii = im/2
      jj = jm/2
      
!      print *,me,'max(gdlat)=', maxval(gdlat),  &
!                 'max(gdlon)=', maxval(gdlon)
      CALL EXCH(gdlat(ISTA_2L,JSTA_2L))
      CALL EXCH(gdlon(ISTA_2L,JSTA_2L))
!      print *,'after call EXCH,me=',me

!$omp parallel do private(i,j,ip1)
      do j = jsta, jend_m
        do i = ista, iend_m
          ip1 = i + 1
!          if (ip1 > im) ip1 = ip1 - im
          if(MAPTYPE==207)then
            DX(i,j) = erad*dxval*dtr/gdsdegr
          else
            DX(i,j) = ERAD*COS(GDLAT(I,J)*DTR) *(GDLON(IP1,J)-GDLON(I,J))*DTR
          endif
          if(MAPTYPE==207)then
            DY(i,j)= erad*dyval*dtr/gdsdegr
          else
            DY(i,j) = ERAD*(GDLAT(I,J+1)-GDLAT(I,J))*DTR  ! like A*DPH
          endif
!	  F(I,J)=1.454441e-4*sin(gdlat(i,j)*DTR)         ! 2*omeg*sin(phi)
!     if (i == ii .and. j == jj) print*,'sample LATLON, DY, DY='    &
!           ,i,j,GDLAT(I,J),GDLON(I,J),DX(I,J),DY(I,J)
        end do
      end do
      if(debugprint)print*,'me sample dx dy= ' &
     ,me,dx(isa,jsa),dy(isa,jsa)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          F(I,J) = 1.454441e-4*sin(gdlat(i,j)*DTR)   ! 2*omeg*sin(phi)
        end do
      end do
      
      iyear = idate(1)
      imn   = idate(2)  
      iday  = idate(3) 
      ihrst = idate(4)
      imin  = idate(5)
      jdate = 0
      idate = 0 
!
      if(me==0)then
      print*,'start yr mo day hr min =',iyear,imn,iday,ihrst,imin
      print*,'processing yr mo day hr min='                            &
             ,idat(3),idat(1),idat(2),idat(4),idat(5)
      endif
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
      !print *,' idate=',idate
      !print *,' jdate=',jdate
!
      CALL W3DIFDAT(JDATE,IDATE,0,RINC)
!
!      print *,' rinc=',rinc
      ifhr = nint(rinc(2)+rinc(1)*24.)
      !print *,' ifhr=',ifhr
      ifmin = nint(rinc(3))
!      if(ifhr /= nint(fhour))print*,'find wrong Grib file';stop
!      print*,' in INITPOST ifhr ifmin fileName=',ifhr,ifmin,fileName
      
! Getting tstart
      tstart = 0.
      !print*,'tstart= ',tstart
      
! Getiing restart
      
      RESTRT = .TRUE.  ! set RESTRT as default
            
      IF(tstart > 1.0E-2)THEN
        ifhr    = ifhr+NINT(tstart)
        rinc    = 0
        idate   = 0
        rinc(2) = -1.0*ifhr
        call w3movdat(rinc,jdate,idate)
        SDAT(1) = idate(2)
        SDAT(2) = idate(3)
        SDAT(3) = idate(1)
        IHRST   = idate(5)       
        !print*,'new forecast hours for restrt run= ',ifhr
        !print*,'new start yr mo day hr min =',sdat(3),sdat(1)           &
        !       ,sdat(2),ihrst,imin
      END IF 
      
! GFS does not need DT to compute accumulated fields, set it to one
!      VarName='DT'
      DT   = 1

      HBM2 = 1.0

! start reading 3d netcdf output
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(1),uh(ista_2l,jsta_2l,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(2),vh(ista_2l,jsta_2l,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(3),q(ista_2l,jsta_2l,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(4),t(ista_2l,jsta_2l,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(5),o3(ista_2l,jsta_2l,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(7),wh(ista_2l,jsta_2l,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(8),qqw(ista_2l,jsta_2l,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(9),dpres(ista_2l,jsta_2l,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(10),buf3d(ista_2l,jsta_2l,1),lm)
       do l=1,lm
       do j=jsta,jend
         do i=ista,iend
            cwm(i,j,l)=spval
! dong add missing value
           if (wh(i,j,l) < spval) then
            omga(i,j,l)=(-1.)*wh(i,j,l)*dpres(i,j,l)/abs(buf3d(i,j,l))
           else
            omga(i,j,l) = spval
           end if
!           if(t(i,j,l)>1000.)print*,'bad T ',t(i,j,l)
         enddo
       enddo
       enddo
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(11),qqi(ista_2l,jsta_2l,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(12),qqr(ista_2l,jsta_2l,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(13),qqs(ista_2l,jsta_2l,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(14),qqg(ista_2l,jsta_2l,1),lm)
! read for regional FV3
       if (modelname == 'FV3R') then
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(15),smoke(ista_2l,jsta_2l,1,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(16),fv3dust(ista_2l,jsta_2l,1,1),lm)
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(17),coarsepm(ista_2l,jsta_2l,1,1),lm)
       call read_netcdf_3d_para(ncid2d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,recname(18),ext550(ista_2l,jsta_2l,1),lm)
       endif

! Compute max QRAIN in the column to be used later in precip type computation
       do j = jsta_2l, jend_2u
        do i = ista_2l, iend_2u
           qrmax(i,j)=0.
        end do
       end do

! calculate CWM from FV3 output
       do l=1,lm
       do j=jsta,jend
         do i=ista,iend
            qrmax(i,j)=max(qrmax(i,j),qqr(i,j,l))
            cwm(i,j,l)=qqg(i,j,l)+qqs(i,j,l)+qqr(i,j,l)+qqi(i,j,l)+qqw(i,j,l)
         enddo
       enddo
       if(debugprint)print*,'sample l,t,q,u,v,w= ',isa,jsa,l &
       ,t(isa,jsa,l),q(isa,jsa,l),uh(isa,jsa,l),vh(isa,jsa,l) &
       ,wh(isa,jsa,l)
       if(debugprint)print*,'sample l cwm for FV3',l, &
          cwm(isa,jsa,l)
      end do 

! instantaneous 3D cloud fraction
      if ( imp_physics==11) then !GFDL MP
        VarName='cld_amt'
        call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
        spval,VarName,cfr(ista_2l,jsta_2l,1),lm)
      else

        iret_bl = nf90_inq_varid(ncid2d,'cldfra_bl',varid_bl)
        iret = nf90_inq_varid(ncid2d,'cldfra',varid)

        if(iret_bl==NF90_NOERR .and. iret==NF90_NOERR) then
          write(*,*) 'WARNING: BOTH cldfra_bl AND cldfra ARE AVAILABLE. USING cldfra.'
          VarName='cldfra'
        else if(iret_bl==NF90_NOERR) then
          VarName='cldfra_bl'
        else if(iret==NF90_NOERR) then
          VarName='cldfra'
        else
          VarName='nope'
        endif
          
        if(VarName /= 'nope') then
          call read_netcdf_3d_para(ncid2d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
               spval,VarName,cfr(ista_2l,jsta_2l,1),lm)
        endif
      endif
!      do l=1,lm
!       if(debugprint)print*,'sample ',VarName,'isa,jsa,l =' &
!          ,cfr(isa,jsa,l),isa,jsa,l
!      enddo

!     WL add cieffr for Thompson scheme cloud ice effective radius
      VarName='cieffr'
      call read_netcdf_3d_para(ncid2d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,effri(ista_2l,jsta_2l,1),lm)

!     WL add cleffr for Thompson scheme cloud water effective radius
      VarName='cleffr'
      call read_netcdf_3d_para(ncid2d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,effrl(ista_2l,jsta_2l,1),lm)

!     WL add cseffr for Thompson scheme snow effective radius
      VarName='cseffr'
      call read_netcdf_3d_para(ncid2d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,effrs(ista_2l,jsta_2l,1),lm)

!=====================================
! For AQF Hourly average field PM2.5
!=====================================

      if (aqf_on) then

       ! *********** VarName need to be in lower case ************
       ! === It will cause problem if not use the lower case =====
       ! *********************************************************

       !-- rename input o3_ave and pm25_ave to NCO grib2 name OZCON and PMTF

       VarName='o3_ave'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,VarName,avgozcon(ista_2l,jsta_2l,1),lm)

       VarName='pm25_ave'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,VarName,avgpmtf(ista_2l,jsta_2l,1),lm)

       VarName='aod'
       call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
       spval,VarName,aqm_aod550(ista_2l,jsta_2l))

      endif     ! -- aqf_on
!============================

! read for regional FV3
      if (modelname == 'FV3R') then
! max hourly updraft velocity
      VarName='upvvelmax'
      call read_netcdf_2d_para(ncid3d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,w_up_max(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' = ',w_up_max(isa,jsa)
! max hourly downdraft velocity
      VarName='dnvvelmax'
      call read_netcdf_2d_para(ncid3d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,w_dn_max(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' = ',w_dn_max(isa,jsa)
! max hourly updraft helicity
      VarName='uhmax25'
      call read_netcdf_2d_para(ncid3d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,up_heli_max(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' = ',up_heli_max(isa,jsa)
! min hourly updraft helicity
      VarName='uhmin25'
      call read_netcdf_2d_para(ncid3d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,up_heli_min(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' = ',up_heli_min(isa,jsa)
! max hourly 0-3km updraft helicity
      VarName='uhmax03'
      call read_netcdf_2d_para(ncid3d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,up_heli_max03(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' = ',up_heli_max03(isa,jsa)
! min hourly 0-3km updraft helicity
      VarName='uhmin03'
      call read_netcdf_2d_para(ncid3d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,up_heli_min03(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' = ',up_heli_min03(isa,jsa)

! max 0-1km relative vorticity max 
      VarName='maxvort01'
      call read_netcdf_2d_para(ncid3d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,rel_vort_max01(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' = ',rel_vort_max01(isa,jsa)
! max 0-2km relative vorticity max
      VarName='maxvort02'
      call read_netcdf_2d_para(ncid3d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,rel_vort_max(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' =',rel_vort_max(isa,jsa)
! max hybrid lev 1 relative vorticity max
      VarName='maxvorthy1'
      call read_netcdf_2d_para(ncid3d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,rel_vort_maxhy1(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' =',rel_vort_maxhy1(isa,jsa)
! biomass burning emissions
      VarName='ebb_smoke_hr'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ebb(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' =',ebb(isa,jsa)
! hourly wildfire potential
      VarName='hwp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,hwp(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' =',hwp(isa,jsa)
      endif

! lightning threat index 1
      VarName='ltg1_max'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ltg1_max(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' =',ltg1_max(isa,jsa)

! lightning threat index 2
      VarName='ltg2_max'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ltg2_max(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' =',ltg2_max(isa,jsa)

! lightning threat index 3
      VarName='ltg3_max'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ltg3_max(ista_2l,jsta_2l))
     if(debugprint)print*,'sample ',VarName,' =',ltg3_max(isa,jsa)

! surface pressure
      VarName='pressfc'
      call read_netcdf_2d_para(ncid3d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,pint(ista_2l,jsta_2l,lp1))
      do j=jsta,jend
        do i=ista,iend
!          if(pint(i,j,lp1)>1.0E6 .or. pint(ista_2l,jsta_2l,lp1)<50000.) &
!           print*,'bad psfc ',i,j,pint(i,j,lp1)
        end do
      end do
      if(debugprint)print*,'sample ',VarName,' =',pint(isa,jsa,lp1)

      pt = ak5(1)

      do j=jsta,jend
        do i=ista,iend
          pint(i,j,1)= pt
        end do
      end do

      do l=2,lp1
        do j=jsta,jend
          do i=ista,iend
            if (dpres(i,j,l-1)<spval .and. pint(i,j,l-1)<spval) then
              pint(i,j,l)= pint(i,j,l-1) + dpres(i,j,l-1)
            else
              pint(i,j,l)=spval
            endif
          enddo
        enddo
!        if (me == 0) print*,'sample model pint,pmid' ,ii,jj,l &
!          ,pint(ii,jj,l),pmid(ii,jj,l)
      end do

!compute pmid from averaged two layer pint
      do l=lm,1,-1
        do j=jsta,jend
          do i=ista,iend
            if (pint(i,j,l)<spval .and. pint(i,j,l+1)<spval) then
              pmid(i,j,l)=0.5*(pint(i,j,l)+pint(i,j,l+1))
            else
              pmid(i,j,l)=spval
            endif
          enddo
        enddo
      enddo

! surface height from FV3 
! dong set missing value for zint
!      zint=spval
      VarName='hgtsfc'
      call read_netcdf_2d_para(ncid3d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,zint(ista_2l,jsta_2l,lp1))
      if(debugprint)print*,'sample ',VarName,' =',zint(isa,jsa,lp1)
      do j=jsta,jend
        do i=ista,iend
          if (zint(i,j,lp1) /= spval) then
            fis(i,j)      = zint(i,j,lp1) * grav
          else
            fis(i,j)      = spval
          endif
        enddo
      enddo

      do l=lm,1,-1
        do j=jsta,jend
          do i=ista,iend
            if(zint(i,j,l+1)/=spval .and. buf3d(i,j,l)/=spval)then
!make sure delz is positive
             zint(i,j,l)=zint(i,j,l+1)+abs(buf3d(i,j,l))
!             if(zint(i,j,l)>1.0E6)print*,'bad H ',i,j,l,zint(i,j,l)
            else
             zint(i,j,l)=spval
            end if
          end do
        end do
        if(debugprint)print*,'sample zint= ',isa,jsa,l,zint(isa,jsa,l)
      end do

      do l=lp1,1,-1
        do j=jsta,jend
          do i=ista,iend
            alpint(i,j,l)=log(pint(i,j,l))
          end do
        end do
      end do

      do l=lm,1,-1
        do j=jsta,jend
          do i=ista,iend
            if(zint(i,j,l+1)/=spval .and. zint(i,j,l)/=spval &
            .and. pmid(i,j,l)/=spval)then
             zmid(i,j,l)=zint(i,j,l+1)+(zint(i,j,l)-zint(i,j,l+1))* &
                    (log(pmid(i,j,l))-alpint(i,j,l+1))/ &
                    (alpint(i,j,l)-alpint(i,j,l+1))
             if(zmid(i,j,l)>1.0E6)print*,'bad Hmid ',i,j,l,zmid(i,j,l)
            else
             zmid(i,j,l)=spval
            endif
          end do
        end do
      end do

      
      if (gocart_on .or.gccpp_on .or. nasa_on) then

! GFS output dust in nemsio (GOCART)
        dustcb=0.0
        dustallcb=0.0
!       DUST = SPVAL

        VarName='dust1'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt1(ista_2l,jsta_2l,1),lm)
        VarName='dust2'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt2(ista_2l,jsta_2l,1),lm)
        VarName='dust3'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt3(ista_2l,jsta_2l,1),lm)
        VarName='dust4'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt4(ista_2l,jsta_2l,1),lm)
        VarName='dust5'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt5(ista_2l,jsta_2l,1),lm)

      
        do l=1,lm

!$omp parallel do private(i,j)
          do j=jsta,jend
          do i=ista,iend
          dust(i,j,l,1)=dt1(i,j,l)
          dust(i,j,l,2)=dt2(i,j,l)
          dust(i,j,l,3)=dt3(i,j,l)
          dust(i,j,l,4)=dt4(i,j,l)
          dust(i,j,l,5)=dt5(i,j,l)
          

           dustcb(i,j)=dustcb(i,j)+&
           (dust(i,j,l,1)+0.38*dust(i,j,l,2))* &
           dpres(i,j,l)/grav


           dustallcb(i,j)=dustallcb(i,j)+ &
           (dust(i,j,l,1)+dust(i,j,l,2)+ &
           dust(i,j,l,3)+0.74*dust(i,j,l,4))* &
           dpres(i,j,l)/grav
           enddo
           enddo
        end do ! do loop for l

 
! GFS output sea salt in nemsio (GOCART)
        sscb=0.0
        ssallcb=0.0
        
        VarName='seas1'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt1(ista_2l,jsta_2l,1),lm)

        VarName='seas2'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt2(ista_2l,jsta_2l,1),lm)

        VarName='seas3'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt3(ista_2l,jsta_2l,1),lm)

        VarName='seas4'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt4(ista_2l,jsta_2l,1),lm)

        VarName='seas5'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt5(ista_2l,jsta_2l,1),lm)



        do l=1,lm
!$omp parallel do private(i,j)
          do j=jsta,jend
          do i=ista,iend
          salt(i,j,l,1)=dt1(i,j,l)
          salt(i,j,l,2)=dt2(i,j,l)
          salt(i,j,l,3)=dt3(i,j,l)
          salt(i,j,l,4)=dt4(i,j,l)
          salt(i,j,l,5)=dt5(i,j,l)

            sscb(i,j)=sscb(i,j)+ &
         (salt(i,j,l,1)+salt(i,j,l,2)+0.83*salt(i,j,l,3))*  &
           dpres(i,j,l)/grav


          ssallcb(i,j)=ssallcb(i,j)+ &
         (salt(i,j,l,1)+salt(i,j,l,2)+salt(i,j,l,3)+salt(i,j,l,4))* &
           dpres(i,j,l)/grav
           enddo
           enddo
        end do ! do loop for l
! GFS output black carbon in nemsio (GOCART)
        bccb=0.0


        VarName='bc1'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt1(ista_2l,jsta_2l,1),lm)

        VarName='bc2'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt2(ista_2l,jsta_2l,1),lm)

        do l=1,lm
!$omp parallel do private(i,j)
          do j=jsta,jend
          do i=ista,iend

          soot(i,j,l,1)=dt1(i,j,l)
          soot(i,j,l,2)=dt2(i,j,l)

            bccb(i,j)=bccb(i,j)+ &
        (soot(i,j,l,1)+soot(i,j,l,2))* &
           dpres(i,j,l)/grav
           enddo
           enddo
        end do ! do loop for l

        occb=0.0
! GFS output organic carbon in nemsio (GOCART)

        VarName='oc1'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt1(ista_2l,jsta_2l,1),lm)

        VarName='oc2'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt2(ista_2l,jsta_2l,1),lm)
        do l=1,lm
!$omp parallel do private(i,j)
          do j=jsta,jend
          do i=ista,iend
          waso(i,j,l,1)=dt1(i,j,l)
          waso(i,j,l,2)=dt2(i,j,l)

            occb(i,j)=occb(i,j)+ &
        (waso(i,j,l,1)+waso(i,j,l,2))* &
           dpres(i,j,l)/grav
           enddo
           enddo
        end do ! do loop for l

! GFS output sulfate in netcdf (GOCART)
        sulfcb=0.0

!       SUSO = SPVAL
        if (gocart_on .or. gccpp_on) then
        VarName='sulf'
        endif

        if (nasa_on) then
        VarName='so4'
        endif
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt1(ista_2l,jsta_2l,1),lm)
        
         do l=1,lm
!$omp parallel do private(i,j)
          do j=jsta,jend
          do i=ista,iend
          suso(i,j,l,1)=dt1(i,j,l)

            sulfcb(i,j)=sulfcb(i,j)+ &
        suso(i,j,l,1)* &
           dpres(i,j,l)/grav
           enddo
           enddo
        end do ! do loop for l

        if (nasa_on) then
! GFS output nitrate in netcdf (GOCART)
        no3cb=0.0
        VarName='no3an1'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt1(ista_2l,jsta_2l,1),lm)

        VarName='no3an2'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt2(ista_2l,jsta_2l,1),lm)

        VarName='no3an3'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt3(ista_2l,jsta_2l,1),lm)

         do l=1,lm
!$omp parallel do private(i,j)
          do j=jsta,jend
          do i=ista,iend
          no3(i,j,l,1)=dt1(i,j,l)
          no3(i,j,l,2)=dt2(i,j,l)
          no3(i,j,l,3)=dt3(i,j,l)

            no3cb(i,j)=no3cb(i,j)+ &
        (no3(i,j,l,1)+no3(i,j,l,2)+no3(i,j,l,3))* &
           dpres(i,j,l)/grav
           enddo
           enddo
        end do ! do loop for l

! GFS output NH4 in netcdf (GOCART)
        nh4cb=0.0
        VarName='nh4a'
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt1(ista_2l,jsta_2l,1),lm)

         do l=1,lm
!$omp parallel do private(i,j)
          do j=jsta,jend
          do i=ista,iend
          nh4(i,j,l,1)=dt1(i,j,l)

            nh4cb(i,j)=nh4cb(i,j)+ &
        nh4(i,j,l,1)* &
           dpres(i,j,l)/grav
           enddo
           enddo
        end do ! do loop for l


        endif !nasa_on

! GFS output pp25 in nemsio (GOCART)
        pp25cb=0.0

        if (gocart_on .or. gccpp_on) then
        VarName='pp25'
        endif

        if (nasa_on) then
        VarName='pm25'
        endif
       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt1(ista_2l,jsta_2l,1),lm)


! GFS output pp10 in nemsio (GOCART)
        pp10cb=0.0
        if (gocart_on .or. gccpp_on) then
        VarName='pp10'
        endif

        if (nasa_on) then
        VarName='pm10'
        endif

       call read_netcdf_3d_para(ncid3d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
       spval,VarName,dt2(ista_2l,jsta_2l,1),lm)

        do l=1,lm
!$omp parallel do private(i,j)
          do j=jsta,jend
          do i=ista,iend
          pp25(i,j,l,1)=dt1(i,j,l)
          pp10(i,j,l,1)=dt2(i,j,l)


            pp25cb(i,j)=pp25cb(i,j)+ &
        pp25(i,j,l,1)* dpres(i,j,l)/grav

            pp10cb(i,j)=pp10cb(i,j)+ &
        pp10(i,j,l,1)* dpres(i,j,l)/grav
           enddo
           enddo
        end do ! do loop for l
! -- compute air density RHOMID and remove negative tracer values
        do l=1,lm
!$omp parallel do private(i,j,tv)
          do j=jsta,jend
            do i=ista,iend

              TV = T(I,J,L) * (H1+D608*MAX(Q(I,J,L),QMIN))
              RHOMID(I,J,L) = PMID(I,J,L) / (RD*TV)
              do n = 1,  NBIN_DU
                IF ( dust(i,j,l,n) < SPVAL) THEN
                  DUST(i,j,l,n) = MAX(DUST(i,j,l,n), 0.0)
                ENDIF
              enddo
              do n = 1,  NBIN_SS
                IF ( salt(i,j,l,n) < SPVAL) THEN
                  SALT(i,j,l,n) = MAX(SALT(i,j,l,n), 0.0)
                ENDIF
              enddo
              do n = 1,  NBIN_OC
                IF ( waso(i,j,l,n) < SPVAL) THEN
                  WASO(i,j,l,n) = MAX(WASO(i,j,l,n), 0.0)
                ENDIF
              enddo
              do n = 1,  NBIN_BC
                IF ( soot(i,j,l,n) < SPVAL) THEN
                  SOOT(i,j,l,n) = MAX(SOOT(i,j,l,n), 0.0)
                ENDIF
              enddo
              do n = 1,  NBIN_SU
                IF ( suso(i,j,l,n) < SPVAL) THEN
                  SUSO(i,j,l,n) = MAX(SUSO(i,j,l,n), 0.0)
                ENDIF
              enddo
          if (nasa_on) then
              do n = 1,  NBIN_NO3
                IF ( no3(i,j,l,n) < SPVAL) THEN
                  no3(i,j,l,n) = MAX(no3(i,j,l,n), 0.0)
                ENDIF
              enddo
              do n = 1,  NBIN_NH4
                IF ( nh4(i,j,l,n) < SPVAL) THEN
                  nh4(i,j,l,n) = MAX(nh4(i,j,l,n), 0.0)
                ENDIF
              enddo
          endif !nasa_on
            end do
          end do
        end do
             l=lm
!$omp parallel do private(i,j)
          do j=jsta,jend
            do i=ista,iend
            dustcb(i,j) = MAX(dustcb(i,j), 0.0)
            dustallcb(i,j) = MAX(dustallcb(i,j), 0.0)
            sscb(i,j) = MAX(sscb(i,j), 0.0)
            ssallcb(i,j) = MAX(ssallcb(i,j), 0.0)
            bccb(i,j) = MAX(bccb(i,j), 0.0)
            occb(i,j) = MAX(occb(i,j), 0.0)
            sulfcb(i,j) = MAX(sulfcb(i,j), 0.0)
            if (nasa_on) then
            no3cb(i,j) = MAX(no3cb(i,j), 0.0)
            nh4cb(i,j) = MAX(nh4cb(i,j), 0.0)
            endif
            pp25cb(i,j) = MAX(pp25cb(i,j), 0.0)
            pp10cb(i,j) = MAX(pp10cb(i,j), 0.0)

!      PM25 dust and seasalt      
       dustpm(i,j)=(dust(i,j,l,1)+0.38*dust(i,j,l,2))*RHOMID(i,j,l) !ug/m3
       dustpm10(i,j)=(dust(i,j,l,1)+dust(i,j,l,2)+dust(i,j,l,3)+ &
       0.74*dust(i,j,l,4))*RHOMID(i,j,l) !ug/m3
       sspm(i,j)=(salt(i,j,l,1)+salt(i,j,l,2)+ &
       0.83*salt(i,j,l,3))*RHOMID(i,j,l)  !ug/m3 

       if (gocart_on .or. gccpp_on) then
!      PM10 concentration
       dusmass(i,j)=(dust(i,j,l,1)+dust(i,j,l,2)+dust(i,j,l,3)+ &
       0.74*dust(i,j,l,4)+salt(i,j,l,1)+salt(i,j,l,2)+salt(i,j,l,3)+ &
       salt(i,j,l,4) + soot(i,j,l,1)+soot(i,j,l,2)+waso(i,j,l,1)+ &
       waso(i,j,l,2) +suso(i,j,l,1)+pp25(i,j,l,1)+pp10(i,j,l,1)) &
       *RHOMID(i,j,l)  !ug/m3
!      PM25 concentration       
       dusmass25(i,j)=(dust(i,j,l,1)+0.38*dust(i,j,l,2)+ &
       salt(i,j,l,1)+salt(i,j,l,2)+0.83*salt(i,j,l,3) + &
       soot(i,j,l,1)+soot(i,j,l,2)+waso(i,j,l,1)+ &
       waso(i,j,l,2) +suso(i,j,l,1)+pp25(i,j,l,1))*RHOMID(i,j,l)  !ug/m3

!      PM10 column
        ducmass(i,j)=dustallcb(i,j)+ssallcb(i,j)+bccb(i,j)+ &
         occb(i,j)+sulfcb(i,j)+pp25cb(i,j)+pp10cb(i,j)
!      PM25 column
        ducmass25(i,j)=dustcb(i,j)+sscb(i,j)+bccb(i,j)+occb(i,j) &
         +sulfcb(i,j)+pp25cb(i,j)
       endif !gocart_on

       if (nasa_on) then
!      PM10 concentration
       dusmass(i,j)=pp10(i,j,l,1)*RHOMID(i,j,l)  !ug/m3
!      PM25 concentration       
       dusmass25(i,j)=pp25(i,j,l,1)*RHOMID(i,j,l)  !ug/m3

!      PM10 column
        ducmass(i,j)=pp10cb(i,j)
!      PM25 column
        ducmass25(i,j)=pp25cb(i,j)
       endif !nasa_on


            end do
          end do
       

      endif                     ! endif for gocart_on & nasa_on

!         ',ll,waso(isa,jsa,ll,2)
!         ',ll,waso(isa,jsa,ll,2)


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!

! done with 3d file, close it for now
      Status=nf90_close(ncid3d)
      deallocate(recname)

! IVEGSRC=1 for IGBP, 0 for USGS, 2 for UMD
      VarName='IVEGSRC'
      Status=nf90_get_att(ncid2d,nf90_global,'IVEGSRC',IVEGSRC)
      if (Status /= 0) then
       if(me==0)print*,VarName,' not found-Assigned 1 for IGBP as default'
       IVEGSRC=1
      end if
      if (me == 0) print*,'IVEGSRC= ',IVEGSRC

! set novegtype based on vegetation classification
      if(ivegsrc==2)then
       novegtype=13
      else if(ivegsrc==1)then
       novegtype=20
      else if(ivegsrc==0)then
       novegtype=24
      end if
      if (me == 0) print*,'novegtype= ',novegtype

      Status=nf90_get_att(ncid2d,nf90_global,'fhzero',fhzero)
      if (Status /= 0) then
       print*,'fhzero not found-Assigned 3 hours as default'
       fhzero=3
      end if
      if (me == 0) print*,'fhzero= ',fhzero
!
      Status=nf90_get_att(ncid2d,nf90_global,'dtp',dtp)
      if (Status /= 0) then
       print*,'dtp not found-Assigned 90s as default'
       dtp=90
      end if
      if (me == 0) print*,'dtp= ',dtp
! Initializes constants for Ferrier microphysics
      if(imp_physics==5 .or. imp_physics==85 .or. imp_physics==95) then
        CALL MICROINIT(imp_physics)
      end if

        tprec   = float(fhzero)
        if(ifhr>240)tprec=12.
        tclod   = tprec
        trdlw   = tprec
        trdsw   = tprec
        tsrfc   = tprec
        tmaxmin = tprec
        td3d    = tprec
        !print*,'tprec = ',tprec


      VarName='refl_10cm'
      call read_netcdf_3d_para(ncid2d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,REF_10CM(ista_2l,jsta_2l,1),lm)
!      do l=1,lm
!       if(debugprint)print*,'sample ',VarName,'isa,jsa,l =' &
!          ,REF_10CM(isa,jsa,l),isa,jsa,l
!      enddo

! turbulence kinetic energy (QKE = 2*TKE)
      VarName='qke'
      call read_netcdf_3d_para(ncid2d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,q2(ista_2l,jsta_2l,1),lm)
      do l=1,lm
      do j=jsta,jend
      do i=ista,iend
        q2(i,j,l)=q2(i,j,l)/2.0
      enddo
      enddo
      enddo

! ice-friendly aerosol number concentration
      VarName='nifa'
      call read_netcdf_3d_para(ncid2d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,qqnifa(ista_2l,jsta_2l,1),lm)

! water-friendly aerosol number concentration
      VarName='nwfa'
      call read_netcdf_3d_para(ncid2d,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,qqnwfa(ista_2l,jsta_2l,1),lm)

      VarName='land' 
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sm)
      if(debugprint)print*,'sample ',VarName,' =',sm((ista+iend)/2,(jsta+jend)/2)

!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= spval) sm(i,j) = 1.0 - sm(i,j)
        enddo
      enddo

! sea ice mask 

      VarName    = 'icec'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sice)
     if(debugprint)print*,'sample ',VarName,' = ',sice(isa,jsa)

!      where(sice /=spval .and. sice >=1.0)sm=0.0 !sea ice has sea
!      mask=0
! GFS flux files have land points with non-zero sea ice, per Iredell,
! these
! points have sea ice changed to zero, i.e., trust land mask more than
! sea ice
!     where(sm/=spval .and. sm==0.0)sice=0.0 !specify sea ice=0 at land

!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= spval .and. sm(i,j) == 0.0) sice(i,j) = 0.0
        enddo
      enddo


! PBL height using nemsio
      VarName    = 'hpbl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,pblh)
      if(debugprint)print*,'sample ',VarName,' = ',pblh(isa,jsa)

! frictional velocity using nemsio
      VarName='fricv'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ustar)
!     if(debugprint)print*,'sample ',VarName,' = ',ustar(isa,jsa)

! roughness length using getgb
      VarName='sfcr'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,z0)
!     if(debugprint)print*,'sample ',VarName,' = ',z0(isa,jsa)

! sfc exchange coeff
      VarName='sfexc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,SFCEXC)

! accumulated snowfall
      VarName='tsnowp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,SNOW_ACM)
! snowfall bucket
      VarName='tsnowpb'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,SNOW_BKT)

! accumulated graupel/sleet
      VarName='frozr'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,acgraup)

! graupel/sleet bucket
      VarName='frozrb'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,graup_bucket)

! accumulated freezing rain
      VarName='frzr'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,acfrain)

! freezing rain bucket
      VarName='frzrb'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,frzrn_bucket)

! time step snow (in m)
      VarName='snow'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,snownc)

! time step graupel (in m)
      VarName='graupel'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,graupelnc)

! aerodynamic conductance
      VarName='acond'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,acond)
      if(debugprint)print*,'sample ',VarName,' = ',acond(isa,jsa)
! mid day avg albedo
      VarName='albdo_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgalbedo)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (avgalbedo(i,j) /= spval) avgalbedo(i,j) = avgalbedo(i,j) * 0.01
        enddo
      enddo
      if(debugprint)print*,'sample ',VarName,' = ',avgalbedo(isa,jsa)

! surface potential T  using getgb
      VarName='tmpsfc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ths)

!     where(ths/=spval)ths=ths*(p1000/pint(:,:,lp1))**CAPA ! convert to THS

!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (ths(i,j) /= spval) then
!    write(*,*)' i=',i,' j=',j,' ths=',ths(i,j),' pint=',pint(i,j,lp1)
            ths(i,j) = ths(i,j) * (p1000/pint(i,j,lp1))**capa
          endif
          QS(i,j)    = SPVAL ! GFS does not have surface specific humidity
          twbs(i,j)  = SPVAL ! GFS does not have inst sensible heat flux
          qwbs(i,j)  = SPVAL ! GFS does not have inst latent heat flux
!assign sst
          if (sm(i,j) /= 0.0 .and. ths(i,j) < spval ) then
            if (sice(i,j) >= 0.15) then
              sst(i,j) = 271.4
            else
              sst(i,j) = ths(i,j) * (pint(i,j,lp1)/p1000)**capa
            endif
          else
              sst(i,j) = spval
          endif
        enddo
      enddo
     if(debugprint)print*,'sample ',VarName,' = ',ths(isa,jsa)

! foundation temperature
      VarName='tref'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,fdnsst)
      if(debugprint)print*,'sample ',VarName,' = ',fdnsst(isa,jsa)

          
!  GFS does not have time step and physics time step, make up ones since they
! are not really used anyway
!      NPHS=1.
!      DT=90.
!      DTQ2 = DT * NPHS  !MEB need to get physics DT
      DTQ2 = DTP !MEB need to get physics DT
      NPHS=1
      DT = DTQ2/NPHS   !MEB need to get DT
      TSPH = 3600./DT

! convective precip in m per physics time step using getgb
! read 3 hour bucket
      VarName='cpratb_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgcprate)
!     where(avgcprate /= spval)avgcprate=avgcprate*dtq2/1000. ! convert to m
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (avgcprate(i,j) /= spval) avgcprate(i,j) = avgcprate(i,j) * (dtq2*0.001)
        enddo
      enddo
!     if(debugprint)print*,'sample ',VarName,' = ',avgcprate(isa,jsa)
      
!      print*,'maxval CPRATE: ', maxval(CPRATE)

! read continuous bucket
      VarName='cprat_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgcprate_cont)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (avgcprate_cont(i,j) /= spval) avgcprate_cont(i,j) = &
            avgcprate_cont(i,j) * (dtq2*0.001)
        enddo
      enddo
!     if(debugprint)print*,'sample ',VarName,' = ',avgcprate_cont(isa,jsa)

!      print*,'maxval CPRATE: ', maxval(CPRATE)

! precip rate in m per physics time step using getgb
      VarName='prateb_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgprec)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if(avgprec(i,j) /= spval)avgprec(i,j)=avgprec(i,j)*(dtq2*0.001)
        enddo
      enddo

     if(debugprint)print*,'sample ',VarName,' = ',avgprec(isa,jsa)
      
!      prec = avgprec !set avg cprate to inst one to derive other fields

      VarName='prate_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgprec_cont)
!     where(avgprec /= spval)avgprec=avgprec*dtq2/1000. ! convert to m
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (avgprec_cont(i,j) /=spval)avgprec_cont(i,j)=avgprec_cont(i,j) &
               * (dtq2*0.001)
        enddo
      enddo

     if(debugprint)print*,'sample ',VarName,' = ',avgprec_cont(isa,jsa)
! precip rate in m per physics time step
      VarName='tprcp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,prec)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (prec(i,j) /= spval) prec(i,j)=prec(i,j)* (dtq2*0.001) &
              * 1000. / dtp
        enddo
      enddo

! convective precip rate in m per physics time step
      VarName='cnvprcp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,cprate)
!set cprate as 0.
      do j=jsta,jend
        do i=ista,iend
           if (cprate(i,j) /= spval) then
            cprate(i,j) = max(0.,cprate(i,j)) * (dtq2*0.001) &
                 * 1000. / dtp
          else
            cprate(i,j) = 0.
          endif
        enddo
      enddo

! GFS does not have accumulated total, gridscale, and convective precip, will use inst precip to derive in SURFCE.f

! max hourly surface precipitation rate
      VarName='pratemax'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,prate_max)
     if(debugprint)print*,'sample ',VarName,' = ',prate_max(isa,jsa)
! max hourly 1-km agl reflectivity
      VarName='refdmax'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,refd_max)
     if(debugprint)print*,'sample ',VarName,' = ',refd_max(isa,jsa)
! max hourly -10C reflectivity
      VarName='refdmax263k'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,refdm10c_max)
     if(debugprint)print*,'sample ',VarName,' = ',refdm10c_max(isa,jsa)

! max hourly u comp of 10m agl wind
      VarName='u10max'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,u10max)
     if(debugprint)print*,'sample ',VarName,' = ',u10max(isa,jsa)
! max hourly v comp of 10m agl wind
      VarName='v10max'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,v10max)
     if(debugprint)print*,'sample ',VarName,' = ',v10max(isa,jsa)
! max hourly 10m agl wind speed
      VarName='spd10max'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,wspd10max)
     if(debugprint)print*,'sample ',VarName,' = ',wspd10max(isa,jsa)

! inst snow water eqivalent using nemsio
      VarName='weasd'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sno)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j)==0.) sno(i,j) = spval
        enddo
      enddo
     if(debugprint)print*,'sample ',VarName,' = ',sno(isa,jsa)

! ave snow cover 
      VarName='snowc_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,snoavg)
! snow cover is multipled by 100 in SURFCE before writing it out
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j)==1.0 .and. sice(i,j)==0.) snoavg(i,j)=spval
          if(snoavg(i,j)/=spval)snoavg(i,j)=snoavg(i,j)/100.
        end do
      end do

! snow depth in mm using nemsio
      VarName='snod'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,si)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j)==1.0 .and. sice(i,j)==0.) si(i,j)=spval
          if (si(i,j) /= spval) si(i,j) = si(i,j) * 1000.0
          CLDEFI(i,j) = SPVAL ! GFS does not have convective cloud efficiency
          lspa(i,j)   = spval ! GFS does not have similated precip
          TH10(i,j)   = SPVAL ! GFS does not have 10 m theta
          TH10(i,j)   = SPVAL ! GFS does not have 10 m theta
          Q10(i,j)    = SPVAL ! GFS does not have 10 m humidity
          ALBASE(i,j) = SPVAL ! GFS does not have snow free albedo
        enddo
      enddo
     if(debugprint)print*,'sample ',VarName,' = ',si(isa,jsa)

! 2m T using nemsio
      VarName='tmp2m'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,tshltr)
     if(debugprint)print*,'sample ',VarName,' = ',tshltr(isa,jsa)

! GFS does not have 2m pres, estimate it, also convert t to theta 
      Do j=jsta,jend
        Do i=ista,iend
          PSHLTR(I,J)=pint(I,J,lm+1)*EXP(-0.068283/tshltr(i,j))
          tshltr(i,j)= tshltr(i,j)*(p1000/PSHLTR(I,J))**CAPA ! convert to theta
!          if (j == jm/2 .and. mod(i,50) == 0)
!     +   print*,'sample 2m T and P after scatter= '
!     +   ,i,j,tshltr(i,j),pshltr(i,j)
        end do
      end do

! 2m specific humidity using nemsio
      VarName='spfh2m'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,qshltr)
     if(debugprint)print*,'sample ',VarName,' = ',qshltr(isa,jsa)
      
! time averaged column cloud fractionusing nemsio
      VarName='tcdc_aveclm'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgtcdc)
!     where(avgtcdc /= spval)avgtcdc=avgtcdc/100. ! convert to fraction
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (avgtcdc(i,j) /= spval) avgtcdc(i,j) = avgtcdc(i,j) * 0.01
        enddo
      enddo
     if(debugprint)print*,'sample ',VarName,' = ',avgtcdc(isa,jsa)

! GFS probably does not use zenith angle
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
          Czen(i,j)   = spval
          CZMEAN(i,j) = SPVAL      
        enddo
      enddo

! maximum snow albedo in fraction using nemsio
      VarName='snoalb'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,mxsnal)

! land fraction
      VarName='lfrac'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,landfrac)
     
! GFS probably does not use sigt4, set it to sig*t^4
!$omp parallel do private(i,j,tlmh)
      Do j=jsta,jend
        Do i=ista,iend
          TLMH = T(I,J,LM) * T(I,J,LM)
          Sigt4(i,j) = 5.67E-8 * TLMH * TLMH
        End do
      End do

! TG is not used, skip it for now

! GFS does not have inst cloud fraction for high, middle, and low cloud
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
          cfrach(i,j) = spval
          cfracl(i,j) = spval
          cfracm(i,j) = spval
        enddo
      enddo

! ave high cloud fraction using nemsio
      VarName='tcdc_avehcl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgcfrach)
!     where(avgcfrach /= spval)avgcfrach=avgcfrach/100. ! convert to fraction
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (avgcfrach(i,j) /= spval) avgcfrach(i,j) = avgcfrach(i,j) * 0.01
        enddo
      enddo
     if(debugprint)print*,'sample ',VarName,' = ',avgcfrach(isa,jsa)

! ave low cloud fraction using nemsio
      VarName='tcdc_avelcl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgcfracl)
!     where(avgcfracl /= spval)avgcfracl=avgcfracl/100. ! convert to fraction
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (avgcfracl(i,j) /= spval) avgcfracl(i,j) = avgcfracl(i,j) * 0.01
        enddo
      enddo
     if(debugprint)print*,'sample ',VarName,' = ',avgcfracl(isa,jsa)
      
! ave middle cloud fraction using nemsio
      VarName='tcdc_avemcl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgcfracm)
!     where(avgcfracm /= spval)avgcfracm=avgcfracm/100. ! convert to fraction
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (avgcfracm(i,j) /= spval) avgcfracm(i,j) = avgcfracm(i,j) * 0.01
        enddo
      enddo
     if(debugprint)print*,'sample ',VarName,' = ',avgcfracm(isa,jsa)
      
! inst convective cloud fraction using nemsio
      VarName='tcdccnvcl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,cnvcfr)
!     where(cnvcfr /= spval)cnvcfr=cnvcfr/100. ! convert to fraction
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (cnvcfr(i,j) /= spval) cnvcfr (i,j)= cnvcfr(i,j) * 0.01
        enddo
      enddo
!     if(debugprint)print*,'sample ',VarName,' = ',cnvcfr(isa,jsa)
      
! slope type using nemsio
      VarName='sltyp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,buf)
!$omp parallel do private(i,j)
      do j = jsta_2l, jend_2u
        do i=ista,iend
          if (buf(i,j) < spval) then
             islope(i,j) = nint(buf(i,j))
          else
             islope(i,j) = 0
          endif
        enddo
      enddo
!     if(debugprint)print*,'sample ',VarName,' = ',islope(isa,jsa)

! plant canopy sfc wtr in m 
      VarName='cnwat'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,cmc)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (cmc(i,j) /= spval) cmc(i,j) = cmc(i,j) * 0.001
          if (sm(i,j) /= 0.0) cmc(i,j) = spval
        enddo
      enddo
!     if(debugprint)print*,'sample ',VarName,' = ',cmc(isa,jsa)
      
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
          grnflx(i,j) = spval ! GFS does not have inst ground heat flux
        enddo
      enddo

! frozen precip fraction using nemsio
      VarName='cpofp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sr)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if(sr(i,j) /= spval) then
!set range within (0,1)
            sr(i,j)=min(1.,max(0.,sr(i,j)))
          endif
        enddo
      enddo

! sea ice skin temperature
      VarName='tisfc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ti)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sice(i,j) == spval .or. sice(i,j) == 0.) ti(i,j)=spval
        enddo
      enddo

! vegetation fraction in fraction. using nemsio
      VarName='veg'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,vegfrc)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (vegfrc(i,j) /= spval) then
            vegfrc(i,j) = vegfrc(i,j) * 0.01
          else
            vegfrc(i,j) = 0.0
          endif
        enddo
      enddo
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) vegfrc(i,j) = spval
        enddo
      enddo
!     if(debugprint)print*,'sample ',VarName,' = ',vegfrc(isa,jsa)
      
! GFS doesn not yet output soil layer thickness, assign SLDPTH to be the same as nam

         SLDPTH(1) = 0.10
         SLDPTH(2) = 0.3
         SLDPTH(3) = 0.6
         SLDPTH(4) = 1.0

! Eric James, 1 Oct 2021: Because FV3 does not have 1d var "zs", used to
! assign soil depths for RUC LSM, hard wire 9 soil depths here
! so they aren't missing.

       IF (NSOIL==9) THEN
         SLLEVEL(1) = 0.0
         SLLEVEL(2) = 0.01
         SLLEVEL(3) = 0.04
         SLLEVEL(4) = 0.1
         SLLEVEL(5) = 0.3
         SLLEVEL(6) = 0.6
         SLLEVEL(7) = 1.0
         SLLEVEL(8) = 1.6
         SLLEVEL(9) = 3.0
       END IF
 
! liquid volumetric soil mpisture in fraction using nemsio
      VarName='soill1'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sh2o(ista_2l,jsta_2l,1))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) sh2o(i,j,1) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,sh2o(isa,jsa,1)

      VarName='soill2'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sh2o(ista_2l,jsta_2l,2))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) sh2o(i,j,2) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,sh2o(isa,jsa,2)

      VarName='soill3'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sh2o(ista_2l,jsta_2l,3))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) sh2o(i,j,3) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,sh2o(isa,jsa,3)

      VarName='soill4' 
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sh2o(ista_2l,jsta_2l,4))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) sh2o(i,j,4) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,sh2o(isa,jsa,4)

! volumetric soil moisture using nemsio
      VarName='soilw1'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,smc(ista_2l,jsta_2l,1))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) smc(i,j,1) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,smc(isa,jsa,1)
      
      VarName='soilw2'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,smc(ista_2l,jsta_2l,2))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) smc(i,j,2) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,smc(isa,jsa,2)
      
      VarName='soilw3'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,smc(ista_2l,jsta_2l,3))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) smc(i,j,3) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,smc(isa,jsa,3)
      
      VarName='soilw4'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,smc(ista_2l,jsta_2l,4))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) smc(i,j,4) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,smc(isa,jsa,4)

      IF (NSOIL==9) THEN

      VarName='soilw5'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,smc(ista_2l,jsta_2l,5))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) smc(i,j,5) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,smc(isa,jsa,5)

      VarName='soilw6'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,smc(ista_2l,jsta_2l,6))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) smc(i,j,6) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,smc(isa,jsa,6)

      VarName='soilw7'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,smc(ista_2l,jsta_2l,7))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) smc(i,j,7) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,smc(isa,jsa,7)

      VarName='soilw8'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,smc(ista_2l,jsta_2l,8))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) smc(i,j,8) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,smc(isa,jsa,8)

      VarName='soilw9'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,smc(ista_2l,jsta_2l,9))
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) smc(i,j,9) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l',VarName,' = ',1,smc(isa,jsa,9)

      END IF

! soil temperature using nemsio
      VarName='soilt1'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,stc(ista_2l,jsta_2l,1))
!     mask open water areas, combine with sea ice tmp
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) stc(i,j,1) = spval
          !if (sm(i,j) /= 0.0) stc(i,j,1) = spval
        enddo
      enddo
     if(debugprint)print*,'sample l','stc',' = ',1,stc(isa,jsa,1)
      
      VarName='soilt2'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,stc(ista_2l,jsta_2l,2))
!     mask open water areas, combine with sea ice tmp
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) stc(i,j,2) = spval
          !if (sm(i,j) /= 0.0) stc(i,j,2) = spval
        enddo
      enddo
     if(debugprint)print*,'sample stc = ',1,stc(isa,jsa,2)
      
      VarName='soilt3'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,stc(ista_2l,jsta_2l,3))
!     mask open water areas, combine with sea ice tmp
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) stc(i,j,3) = spval 
          !if (sm(i,j) /= 0.0) stc(i,j,3) = spval
        enddo
      enddo
     if(debugprint)print*,'sample stc = ',1,stc(isa,jsa,3)
      
      VarName='soilt4'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,stc(ista_2l,jsta_2l,4))
!     mask open water areas, combine with sea ice tmp
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) stc(i,j,4) = spval
          !if (sm(i,j) /= 0.0) stc(i,j,4) = spval
        enddo
      enddo
     if(debugprint)print*,'sample stc = ',1,stc(isa,jsa,4)

      IF (NSOIL==9) THEN

      VarName='soilt5'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,stc(ista_2l,jsta_2l,5))
!     mask open water areas, combine with sea ice tmp
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) stc(i,j,5) = spval
          !if (sm(i,j) /= 0.0) stc(i,j,5) = spval
        enddo
      enddo
     if(debugprint)print*,'sample stc = ',1,stc(isa,jsa,5)

      VarName='soilt6'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,stc(ista_2l,jsta_2l,6))
!     mask open water areas, combine with sea ice tmp
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) stc(i,j,6) = spval
          !if (sm(i,j) /= 0.0) stc(i,j,6) = spval
        enddo
      enddo
     if(debugprint)print*,'sample stc = ',1,stc(isa,jsa,6)

      VarName='soilt7'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,stc(ista_2l,jsta_2l,7))
!     mask open water areas, combine with sea ice tmp
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) stc(i,j,7) = spval
          !if (sm(i,j) /= 0.0) stc(i,j,7) = spval         enddo
        enddo
      enddo
     if(debugprint)print*,'sample stc = ',1,stc(isa,jsa,7)

      VarName='soilt8'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,stc(ista_2l,jsta_2l,8))
!     mask open water areas, combine with sea ice tmp
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) stc(i,j,8) = spval
          !if (sm(i,j) /= 0.0) stc(i,j,8) = spval
        enddo
      enddo
     if(debugprint)print*,'sample stc = ',1,stc(isa,jsa,8)

      VarName='soilt9'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,stc(ista_2l,jsta_2l,9))
!     mask open water areas, combine with sea ice tmp
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) stc(i,j,9) = spval
          !if (sm(i,j) /= 0.0) stc(i,j,9) = spval
        enddo
      enddo
     if(debugprint)print*,'sample stc = ',1,stc(isa,jsa,9)

      END IF
!
! E. James - 27 Sep 2022: this is for RRFS, adding smoke and dust
! extinction; it needs to be after ZINT is defined.
!
      if (modelname == 'FV3R') then
       do l = 1, lm
        do j = jsta_2l, jend_2u
         do i = ista_2l, iend_2u
          if(ext550(i,j,l)<spval)then
            taod5503d ( i, j, l) = ext550 ( i, j, l )
            dz = ZINT( i, j, l ) - ZINT( i, j, l+1 )
            aextc55 ( i, j, l ) = taod5503d ( i, j, l ) / dz
          endif
         if(debugprint.and.i==im/2.and.j==(jsta+jend)/2)print*,'sample taod5503d= ',   &
           i,j,l,taod5503d ( i, j, l )
         if(debugprint.and.i==im/2.and.j==(jsta+jend)/2)print*,'sample dz= ',          &
           dz
         if(debugprint.and.i==im/2.and.j==(jsta+jend)/2)print*,'sample AEXTC55= ',     &
           i,j,l,aextc55 ( i, j, l )
         end do
        end do
       end do
       deallocate(ext550)
      end if

!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          acfrcv(i,j) = spval ! GFS does not output time averaged convective and strat cloud fraction, set acfrcv to spval, ncfrcv to 1
          ncfrcv(i,j) = 1.0
          acfrst(i,j) = spval ! GFS does not output time averaged cloud fraction, set acfrst to spval, ncfrst to 1
          ncfrst(i,j) = 1.0
          bgroff(i,j) = spval ! GFS does not have UNDERGROUND RUNOFF
        enddo
      enddo
!     trdlw(i,j)  = 6.0
      ardlw = 1.0 ! GFS incoming sfc longwave has been averaged over 6 hr bucket, set ARDLW to 1

! time averaged incoming sfc longwave
      VarName='dlwrf_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,alwin)

! inst incoming sfc longwave
      VarName='dlwrf'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,rlwin)
                                                            
! time averaged outgoing sfc longwave
      VarName='ulwrf_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,alwout)

! inst outgoing sfc longwave 
      VarName='ulwrf'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,radot)

!     where(alwout /= spval) alwout=-alwout ! CLDRAD puts a minus sign before gribbing
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (alwout(i,j) /= spval) alwout(i,j) = -alwout(i,j)
        enddo
      enddo
!     if(debugprint)print*,'sample l',VarName,' = ',1,alwout(isa,jsa)

! time averaged outgoing model top longwave using gfsio
      VarName='ulwrf_avetoa'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,alwtoa)
!     if(debugprint)print*,'sample l',VarName,' = ',1,alwtoa(isa,jsa)

! instant outgoing model top longwave
      VarName='ulwrf_toa'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,rlwtoa)
!     if(debugprint)print*,'sample l',VarName,' = ',1,rlwtoa(isa,jsa)
      
! GFS incoming sfc longwave has been averaged, set ARDLW to 1
      ardsw=1.0
!     trdsw=6.0

! time averaged incoming sfc shortwave 
      VarName='dswrf_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,aswin)
!     if(debugprint)print*,'sample l',VarName,' = ',1,aswin(isa,jsa)

! inst incoming sfc shortwave 
      VarName='dswrf'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,rswin)

! inst incoming clear sky sfc shortwave
! FV3 do not output instant incoming clear sky sfc shortwave
      !$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
          rswinc(i,j) = spval 
        enddo
      enddo

! time averaged incoming sfc uv-b using getgb
      VarName='duvb_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,auvbin)
!     if(debugprint)print*,'sample l',VarName,' = ',1,auvbin(isa,jsa)
       
! time averaged incoming sfc clear sky uv-b using getgb
      VarName='cduvb_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,auvbinc)
!     if(debugprint)print*,'sample l',VarName,' = ',1,auvbinc(isa,jsa)
      
! time averaged outgoing sfc shortwave using gfsio
      VarName='uswrf_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,aswout)
!     where(aswout /= spval) aswout=-aswout ! CLDRAD puts a minus sign before gribbing 
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (aswout(i,j) /= spval) aswout(i,j) = -aswout(i,j)
        enddo
      enddo
!     if(debugprint)print*,'sample l',VarName,' = ',1,aswout(isa,jsa)

! inst outgoing sfc shortwave using gfsio
      VarName='uswrf'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,rswout)

! time averaged model top incoming shortwave
      VarName='dswrf_avetoa'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,aswintoa)
!     if(debugprint)print*,'sample l',VarName,' = ',1,aswintoa(isa,jsa)

! time averaged model top outgoing shortwave
      VarName='uswrf_avetoa'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,aswtoa)
!     if(debugprint)print*,'sample l',VarName,' = ',1,aswtoa(isa,jsa)

! time averaged surface sensible heat flux, multiplied by -1 because wrf model flux
! has reversed sign convention using gfsio
      VarName='shtfl_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sfcshx)
!     where (sfcshx /= spval)sfcshx=-sfcshx
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sfcshx(i,j) /= spval) sfcshx(i,j) = -sfcshx(i,j)
        enddo
      enddo
!     if(debugprint)print*,'sample l',VarName,' = ',1,sfcshx(isa,jsa)

! inst surface sensible heat flux
      VarName='shtfl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,twbs)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (twbs(i,j) /= spval) twbs(i,j) = -twbs(i,j)
        enddo
      enddo

! GFS surface flux has been averaged, set  ASRFC to 1 
      asrfc=1.0  
!      tsrfc=6.0

! time averaged surface latent heat flux, multiplied by -1 because wrf model flux
! has reversed sign vonvention using gfsio
      VarName='lhtfl_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sfclhx)
!     where (sfclhx /= spval)sfclhx=-sfclhx
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sfclhx(i,j) /= spval) sfclhx(i,j) = -sfclhx(i,j)
        enddo
      enddo
!     if(debugprint)print*,'sample l',VarName,' = ',1,sfclhx(isa,jsa)

! inst surface latent heat flux
      VarName='lhtfl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,qwbs)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (qwbs(i,j) /= spval) qwbs(i,j) = -qwbs(i,j)
        enddo
      enddo

      if(me==0)print*,'rdaod= ',rdaod
! inst aod550 optical depth
      if(rdaod) then
      VarName='aod550'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,aod550)

      VarName='du_aod550'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,du_aod550)

      VarName='ss_aod550'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ss_aod550)

      VarName='su_aod550'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,su_aod550)

      VarName='oc_aod550'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,oc_aod550)

      VarName='bc_aod550'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,bc_aod550)
      end if

! time averaged ground heat flux using nemsio
      VarName='gflux_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,subshx)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) subshx(i,j) = spval
        enddo
      enddo
!     if(debugprint)print*,'sample l',VarName,' = ',1,subshx(isa,jsa)

! inst ground heat flux using nemsio
      VarName='gflux'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,grnflx)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) grnflx(i,j) = spval
        enddo
      enddo

! time averaged zonal momentum flux using gfsio
      VarName='uflx_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sfcux)
!     if(debugprint)print*,'sample l',VarName,' = ',1,sfcux(isa,jsa)
      
! time averaged meridional momentum flux using nemsio
      VarName='vflx_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sfcvx)
!     if(debugprint)print*,'sample l',VarName,' = ',1,sfcvx(isa,jsa)

! dong read in inst surface flux 
! inst zonal momentum flux using gfsio
      VarName='uflx'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sfcuxi)
     if(debugprint)print*,'sample l',VarName,' = ',1,sfcuxi(isa,jsa)

! inst meridional momentum flux using nemsio
      VarName='vflx'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,sfcvxi)
     if(debugprint)print*,'sample l',VarName,' = ',1,sfcvxi(isa,jsa)

     
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
          sfcuvx(i,j) = spval ! GFS does not use total momentum flux
        enddo
      enddo

! time averaged zonal gravity wave stress using nemsio
      VarName='u-gwd_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,gtaux)
!     if(debugprint)print*,'sample l',VarName,' = ',1,gtaux(isa,jsa)

! time averaged meridional gravity wave stress using getgb
      VarName='v-gwd_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,gtauy)
!     if(debugprint)print*,'sample l',VarName,' = ',1,gtauy(isa,jsa)
                                                     
! time averaged accumulated potential evaporation
      VarName='pevpr_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgpotevp)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) avgpotevp(i,j) = spval
        enddo
      enddo
!     if(debugprint)print*,'sample l',VarName,' = ',1,potevp(isa,jsa)

! inst potential evaporation
      VarName='pevpr'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,potevp)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) potevp(i,j) = spval
        enddo
      enddo

      do l=1,lm
!$omp parallel do private(i,j)
        do j=jsta_2l,jend_2u
          do i=ista_2l,iend_2u
! GFS does not have temperature tendency due to long wave radiation
            rlwtt(i,j,l)  = spval
! GFS does not have temperature tendency due to short wave radiation
            rswtt(i,j,l)  = spval
! GFS does not have temperature tendency due to latent heating from convection
            tcucn(i,j,l)  = spval
            tcucns(i,j,l) = spval
! GFS does not have temperature tendency due to latent heating from grid scale
            train(i,j,l)  = spval
          enddo
        enddo
      enddo

! set avrain to 1
      avrain=1.0
      avcnvc=1.0
      theat=6.0 ! just in case GFS decides to output T tendency   
      
! 10 m u using nemsio
      VarName='ugrd10m'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,u10)

      do j=jsta,jend
        do i=ista,iend
          u10h(i,j)=u10(i,j)
        end do
      end do
!     if(debugprint)print*,'sample l',VarName,' = ',1,u10(isa,jsa)
            
! 10 m v using gfsio
      VarName='vgrd10m'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,v10)

      do j=jsta,jend
        do i=ista,iend
          v10h(i,j)=v10(i,j)
        end do
      end do
!     if(debugprint)print*,'sample l',VarName,' = ',1,v10(isa,jsa)
      
! vegetation type, it's in GFS surface file, hopefully will merge into gfsio soon 
      VarName='vtype'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,buf)
!     where (buf /= spval)
!      ivgtyp=nint(buf)
!     elsewhere
!      ivgtyp=0 !need to feed reasonable value to crtm
!     end where 
!$omp parallel do private(i,j)
      do j = jsta_2l, jend_2u
        do i=ista,iend
          if (buf(i,j) < spval) then
            ivgtyp(i,j) = nint(buf(i,j))
          else
            ivgtyp(i,j) = 0
          endif
        enddo
      enddo
!     if(debugprint)print*,'sample l',VarName,' = ',1,ivgtyp(isa,jsa)
      
! soil type, it's in GFS surface file, hopefully will merge into gfsio soon
      VarName='sotyp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,buf)
      l=1
!$omp parallel do private(i,j)
      do j = jsta_2l, jend_2u
        do i=ista,iend
          if (buf(i,j) < spval) then
            isltyp(i,j) = nint(buf(i,j))
          else
            isltyp(i,j) = 0 !need to feed reasonable value to crtm
          endif
        enddo
      enddo
!     if(debugprint)print*,'sample l',VarName,' = ',1,isltyp(isa,jsa)
      
      VarName='wetness'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,buf)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          smstav(i,j) = buf(i,j)
        enddo
      enddo
      VarName='snacc_land'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,buf)
      VarName='snacc_ice'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,buf2)
!$omp parallel do private(i,j)
      do j = jsta_2l, jend_2u
        do i=ista,iend
          if(buf(i,j)<spval) then
            sndepac(i,j) = buf(i,j)
          elseif(buf2(i,j)<spval) then
            sndepac(i,j) = buf2(i,j)
          else
            sndepac(i,j) = spval
          endif
        enddo
      enddo
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
!          smstav(i,j) = spval    ! GFS does not have soil moisture availability
!          smstot(i,j) = spval    ! GFS does not have total soil moisture
          sfcevp(i,j) = spval    ! GFS does not have accumulated surface evaporation
          acsnom(i,j) = spval    ! GFS does not have snow melt
!          sst(i,j)    = spval    ! GFS does not have sst????
          thz0(i,j)   = ths(i,j) ! GFS does not have THZ0, use THS to substitute
          qz0(i,j)    = spval    ! GFS does not output humidity at roughness length
          uz0(i,j)    = spval    ! GFS does not output u at roughness length
          vz0(i,j)    = spval    ! GFS does not output humidity at roughness length
        enddo
      enddo
      do l=1,lm
!$omp parallel do private(i,j)
        do j=jsta_2l,jend_2u
          do i=ista_2l,iend_2u
            EL_PBL(i,j,l) = spval    ! GFS does not have mixing length
            exch_h(i,j,l) = spval    ! GFS does not output exchange coefficient
          enddo
        enddo
      enddo
!     if(debugprint)print*,'sample l',VarName,' = ',1,thz0(isa,jsa)

! retrieve inst convective cloud top, GFS has cloud top pressure instead of index,
! will need to modify CLDRAD.f to use pressure directly instead of index
!      VarName='pres'
!      VcoordName='convect-cld top' 
!      l=1
!     if(debugprint)print*,'sample l',VarName,' = ',1,ptop(isa,jsa)
      VarName='prescnvclt'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ptop)

      
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          htop(i,j) = spval
          if(ptop(i,j) <= 0.0) ptop(i,j) = spval
        enddo
      enddo
      do j=jsta,jend
        do i=ista,iend
          if(ptop(i,j) < spval)then
            do l=1,lm
              if(ptop(i,j) <= pmid(i,j,l))then
                htop(i,j) = l
!                if(i==ii .and. j==jj)print*,'sample ptop,pmid pmid-1,pint= ',   &
!                ptop(i,j),pmid(i,j,l),pmid(i,j,l-1),pint(i,j,l),htop(i,j)
                 exit
              end if
            end do
          end if 
        end do
      end do

! retrieve inst convective cloud bottom, GFS has cloud top pressure instead of index,
! will need to modify CLDRAD.f to use pressure directly instead of index
      VarName='prescnvclb'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,pbot)
!     if(debugprint)print*,'sample l',VarName,VcoordName,' = ',1,pbot(isa,jsa)
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          hbot(i,j) = spval
          if(pbot(i,j) <= 0.0) pbot(i,j) = spval
        enddo
      enddo
      do j=jsta,jend
        do i=ista,iend
!	  if(.not.lb(i,j))print*,'false bitmask for pbot at '
!     +	    ,i,j,pbot(i,j)
          if(pbot(i,j) < spval)then
            do l=lm,1,-1
              if(pbot(i,j) >= pmid(i,j,l)) then
                hbot(i,j) = l
!                if(i==ii .and. j==jj)print*,'sample pbot,pmid= ',    &
!                                pbot(i,j),pmid(i,j,l),hbot(i,j)
                exit
              end if
            end do
          end if 
        end do
      end do
      if(debugprint)print*,'sample hbot = ',hbot(isa,jsa)
! retrieve time averaged low cloud top pressure using nemsio
      VarName='pres_avelct'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ptopl)
!     if(debugprint)print*,'sample l',VarName,' = ',1,ptopl(isa,jsa)

! retrieve time averaged low cloud bottom pressure using nemsio
      VarName='pres_avelcb'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,pbotl)
!     if(debugprint)print*,'sample l',VarName,' = ',1,pbotl(isa,jsa)
     
! retrieve time averaged low cloud top temperature using nemsio
      VarName='tmp_avelct'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,Ttopl)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ',1,Ttopl(isa,jsa)

! retrieve time averaged middle cloud top pressure using nemsio
      VarName='pres_avemct'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ptopm)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ',1,ptopm(isa,jsa)
                                                             
! retrieve time averaged middle cloud bottom pressure using  nemsio
      VarName='pres_avemcb'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,pbotm)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ',1,pbotm(isa,jsa)
      
! retrieve time averaged middle cloud top temperature using nemsio
      VarName='tmp_avemct'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,Ttopm)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ',1,Ttopm(isa,jsa)
      
! retrieve time averaged high cloud top pressure using nemsio *********
      VarName='pres_avehct'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,ptoph)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ',1,ptoph(isa,jsa)
     
! retrieve time averaged high cloud bottom pressure using  nemsio
      VarName='pres_avehcb'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,pboth)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ',1,pboth(isa,jsa)

! retrieve time averaged high cloud top temperature using nemsio
      VarName='tmp_avehct'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,Ttoph)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ',1,Ttoph(isa,jsa)
      
! retrieve boundary layer cloud cover using nemsio
      VarName='tcdc_avebndcl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,pblcfr)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ', 1,pblcfr(isa,jsa)
!     where (pblcfr /= spval)pblcfr=pblcfr/100. ! convert to fraction
!$omp parallel do private(i,j)
      do j = jsta_2l, jend_2u
        do i=ista,iend
          if (pblcfr(i,j) < spval) pblcfr(i,j) = pblcfr(i,j) * 0.01
        enddo
      enddo
        
! retrieve cloud work function 
      VarName='cwork_aveclm'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,cldwork)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ', 1,cldwork(isa,jsa)
      
! accumulated total (base+surface) runoff
      VarName='watr_acc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,runoff)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) runoff(i,j) = spval
        enddo
      enddo

! total water storage in aquifer
      VarName='wa_acc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,twa)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) twa(i,j) = spval
        enddo
      enddo
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ', 1,runoff(isa,jsa)

! accumulated evaporation of intercepted water
      VarName='ecan_acc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,tecan)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) tecan(i,j) = spval
        enddo
      enddo

! accumulated plant transpiration
      VarName='etran_acc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,tetran)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) tetran(i,j) = spval
        enddo
      enddo

! accumulated soil surface evaporation
      VarName='edir_acc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,tedir)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) tedir(i,j) = spval
        enddo
      enddo

! retrieve shelter max temperature using nemsio
      VarName='t02max'
      if(modelname=='GFS') VarName='tmax_max2m'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,maxtshltr)

! retrieve shelter min temperature using nemsio
      VarName='t02min'
      if(modelname=='GFS') VarName='tmin_min2m'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,mintshltr)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
!     1,mintshltr((ista+iend)/2,(jsta+jend)/2)

! retrieve shelter max RH
      VarName='rh02max'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,maxrhshltr)

! retrieve shelter min temperature using nemsio
      VarName='rh02min'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,minrhshltr)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ', &
!     1,mintshltr((ista+iend)/2,(jsta+jend)/2)

! retrieve shelter max specific humidity using nemsio
      VarName='spfhmax_max2m'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,maxqshltr)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ',
!     1,maxqshltr(isa,jsa)

! retrieve shelter min temperature using nemsio
      VarName='spfhmin_min2m'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,minqshltr)
 
! retrieve ice thickness using nemsio
      VarName='icetk'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,dzice)
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ', 1,dzice(isa,jsa)

! retrieve wilting point using nemsio
      VarName='wilt'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,smcwlt)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) smcwlt(i,j) = spval
        enddo
      enddo
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ', 1,smcwlt(isa,jsa)
      
! retrieve sunshine duration using nemsio
      VarName='sunsd_acc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,suntime)

! retrieve field capacity using nemsio
      VarName='fldcp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,fieldcapa)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) fieldcapa(i,j) = spval
        enddo
      enddo
!     if(debugprint)print*,'sample l',VcoordName,VarName,' = ', 1,fieldcapa(isa,jsa)

! retrieve time averaged surface visible beam downward solar flux
      VarName='vbdsf_ave'
       call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avisbeamswin)
      VcoordName='sfc'
      l=1

! retrieve time averaged surface visible diffuse downward solar flux
      VarName='vddsf_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avisdiffswin)

! retrieve time averaged surface near IR beam downward solar flux
      VarName='nbdsf_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,airbeamswin)

! retrieve time averaged surface near IR diffuse downward solar flux
      VarName='nddsf_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,airdiffswin)

! retrieve time averaged surface clear sky outgoing LW
      VarName='csulf'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,alwoutc)

! retrieve time averaged TOA clear sky outgoing LW
      VarName='csulftoa'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,alwtoac)

! retrieve time averaged surface clear sky outgoing SW
      VarName='csusf'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,aswoutc)

! retrieve time averaged TOA clear sky outgoing LW
      VarName='csusftoa'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,aswtoac)

! retrieve time averaged surface clear sky incoming LW
      VarName='csdlf'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,alwinc)

! retrieve time averaged surface clear sky incoming SW
      VarName='csdsf'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,aswinc)

! retrieve storm runoff using nemsio
      VarName='ssrun_acc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,SSROFF)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) ssroff(i,j) = spval
        enddo
      enddo

! retrieve direct soil evaporation
      VarName='evbs_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgedir)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) avgedir(i,j) = spval
        enddo
      enddo

! retrieve CANOPY WATER EVAP 
      VarName='evcw_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgecan)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) avgecan(i,j) = spval
        enddo
      enddo

! retrieve AVERAGED PRECIP ADVECTED HEAT FLUX
      VarName='pah_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,paha)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) paha(i,j) = spval
        enddo
      enddo

! retrieve instantaneous PRECIP ADVECTED HEAT FLUX
      VarName='pahi'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,pahi)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) pahi(i,j) = spval
        enddo
      enddo

! retrieve PLANT TRANSPIRATION 
      VarName='trans_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgetrans)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) avgetrans(i,j) = spval
        enddo
      enddo

! retrieve snow sublimation
      VarName='sbsno_ave'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,avgesnow)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j)==1.0 .and. sice(i,j)==0.) avgesnow(i,j)=spval
        enddo
      enddo

! retrive total soil moisture
      VarName='soilm'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,smstot)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) smstot(i,j) = spval
        enddo
      enddo

! retrieve snow phase change heat flux
      VarName='snohf'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,snopcx)
!     mask water areas
!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          if (sm(i,j) /= 0.0) snopcx(i,j) = spval
        enddo
      enddo

! retrieve pwater
      VarName='pwat'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
      spval,VarName,pwat)
      
! GFS does not have deep convective cloud top and bottom fields

!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          HTOPD(i,j) = SPVAL
          HBOTD(i,j) = SPVAL   
          HTOPS(i,j) = SPVAL
          HBOTS(i,j) = SPVAL 
          CUPPT(i,j) = SPVAL 
        enddo
      enddo


      if ((gocart_on .or. gccpp_on) .and. d2d_chem) then


! retrieve dust emission fluxes
       do K = 1, nbin_du
       if ( K == 1) VarName='duem001'
       if ( K == 2) VarName='duem002'
       if ( K == 3) VarName='duem003'
       if ( K == 4) VarName='duem004'
       if ( K == 5) VarName='duem005'

      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)

      duem(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve dust sedimentation fluxes
      do K = 1, nbin_du
       if ( K == 1) VarName='dust1sd'
       if ( K == 2) VarName='dust2sd'
       if ( K == 3) VarName='dust3sd'
       if ( K == 4) VarName='dust4sd'
       if ( K == 5) VarName='dust5sd'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      dusd(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve dust dry deposition fluxes
      do K = 1, nbin_du
       if ( K == 1) VarName='dust1dp'
       if ( K == 2) VarName='dust2dp'
       if ( K == 3) VarName='dust3dp'
       if ( K == 4) VarName='dust4dp'
       if ( K == 5) VarName='dust5dp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      dudp(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve dust wet deposition fluxes
      do K = 1, nbin_du
       if ( K == 1) VarName='dust1wtl'
       if ( K == 2) VarName='dust2wtl'
       if ( K == 3) VarName='dust3wtl'
       if ( K == 4) VarName='dust4wtl'
       if ( K == 5) VarName='dust5wtl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      duwt(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve dust scavenging fluxes
      do K = 1, nbin_du
       if ( K == 1) VarName='dust1wtc'
       if ( K == 2) VarName='dust2wtc'
       if ( K == 3) VarName='dust3wtc'
       if ( K == 4) VarName='dust4wtc'
       if ( K == 5) VarName='dust5wtc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      dusv(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve seasalt emission fluxes
      do K = 1, nbin_ss
       if ( K == 1) VarName='ssem001'
       if ( K == 2) VarName='ssem002'
       if ( K == 3) VarName='ssem003'
       if ( K == 4) VarName='ssem004'
       if ( K == 5) VarName='ssem005'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      ssem(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve seasalt emission fluxes
      do K = 1, nbin_ss
       if ( K == 1) VarName='seas1sd'
       if ( K == 2) VarName='seas2sd'
       if ( K == 3) VarName='seas3sd'
       if ( K == 4) VarName='seas4sd'
       if ( K == 5) VarName='seas5sd'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      sssd(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve seasalt dry deposition fluxes
      do K = 1, nbin_ss
       if ( K == 1) VarName='seas1dp'
       if ( K == 2) VarName='seas2dp'
       if ( K == 3) VarName='seas3dp'
       if ( K == 4) VarName='seas4dp'
       if ( K == 5) VarName='seas5dp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      ssdp(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve seasalt wet deposition fluxes
      do K = 1, nbin_ss
       if ( K == 1) VarName='seas1wtl'
       if ( K == 2) VarName='seas2wtl'
       if ( K == 3) VarName='seas3wtl'
       if ( K == 4) VarName='seas4wtl'
       if ( K == 5) VarName='seas5wtl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      sswt(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve seasalt scavenging fluxes
      do K = 1, nbin_ss
       if ( K == 1) VarName='seas1wtc'
       if ( K == 2) VarName='seas2wtc'
       if ( K == 3) VarName='seas3wtc'
       if ( K == 4) VarName='seas4wtc'
       if ( K == 5) VarName='seas5wtc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      sssv(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve bc emission fluxes
      do K = 1, nbin_bc
       if ( K == 1) VarName='bceman'
       if ( K == 2) VarName='bcembb'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      bcem(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve bc sedimentation fluxes
      do K = 1, nbin_bc
       if ( K == 1) VarName='bc1sd'
       if ( K == 2) VarName='bc2sd'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      bcsd(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve bc dry deposition fluxes
      do K = 1, nbin_bc
       if ( K == 1) VarName='bc1dp'
       if ( K == 2) VarName='bc2dp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      bcdp(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve bc large wet deposition fluxes
      do K = 1, nbin_bc
       if ( K == 1) VarName='bc1wtl'
       if ( K == 2) VarName='bc2wtl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      bcwt(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve bc convective wet deposition fluxes
      do K = 1, nbin_bc
       if ( K == 1) VarName='bc1wtc'
       if ( K == 2) VarName='bc2wtc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      bcsv(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve oc emission fluxes
      do K = 1, nbin_oc
       if ( K == 1) VarName='oceman'
       if ( K == 2) VarName='ocembb'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      ocem(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve oc sedimentation fluxes
      do K = 1, nbin_oc
       if ( K == 1) VarName='oc1sd'
       if ( K == 2) VarName='oc2sd'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      ocsd(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve oc dry deposition fluxes
      do K = 1, nbin_oc
       if ( K == 1) VarName='oc1dp'
       if ( K == 2) VarName='oc2dp'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      ocdp(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve oc large wet deposition fluxes
      do K = 1, nbin_oc
       if ( K == 1) VarName='oc1wtl'
       if ( K == 2) VarName='oc2wtl'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      ocwt(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve oc convective wet deposition fluxes
      do K = 1, nbin_oc
       if ( K == 1) VarName='oc1wtc'
       if ( K == 2) VarName='oc2wtc'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      ocsv(1:im,jsta_2l:jend_2u,K)=chem_2d(1:im,jsta_2l:jend_2u)
       enddo

! retrieve MIE AOD
       VarName='maod'
      call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u,&
      spval,VarName,chem_2d)
      maod(1:im,jsta_2l:jend_2u)=chem_2d(1:im,jsta_2l:jend_2u)

       endif ! gocart_on
! done with flux file, close it for now
      Status=nf90_close(ncid2d)
!      deallocate(tmp,recname,reclevtyp,reclev)

! pos east
!       call collect_loc(gdlat,dummy)
!       if(me == 0)then
!        latstart = nint(dummy(1,1)*gdsdegr)
!        latlast  = nint(dummy(im,jm)*gdsdegr)
!        print*,'laststart,latlast B bcast= ',latstart,latlast,'gdsdegr=',gdsdegr,&
!          'dummy(1,1)=',dummy(1,1),dummy(im,jm),'gdlat=',gdlat(1,1)
!       end if
!       call mpi_bcast(latstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
!       call mpi_bcast(latlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
!       write(6,*) 'laststart,latlast,me A calling bcast=',latstart,latlast,me
!       call collect_loc(gdlon,dummy)
!       if(me == 0)then
!        lonstart = nint(dummy(1,1)*gdsdegr)
!        lonlast  = nint(dummy(im,jm)*gdsdegr)
!       end if
!       call mpi_bcast(lonstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
!       call mpi_bcast(lonlast, 1,MPI_INTEGER,0,mpi_comm_comp,irtn)

!       write(6,*)'lonstart,lonlast A calling bcast=',lonstart,lonlast
!

! generate look up table for lifted parcel calculations

      THL    = 210.
      PLQ    = 70000.
      pt_TBL = 10000.          ! this is for 100 hPa added by Moorthi

      CALL TABLE(PTBL,TTBL,PT_TBL,                                     &
                 RDQ,RDTH,RDP,RDTHE,PL,THL,QS0,SQS,STHE,THE0)

      CALL TABLEQ(TTBLQ,RDPQ,RDTHEQ,PLQ,THL,STHEQ,THE0Q)

!     
!     
      IF(ME == 0)THEN
        WRITE(6,*)'  SPL (POSTED PRESSURE LEVELS) BELOW: '
        WRITE(6,51) (SPL(L),L=1,LSM)
   50   FORMAT(14(F4.1,1X))
   51   FORMAT(8(F8.1,1X))
      ENDIF
!     
!$omp parallel do private(l)
      DO L = 1,LSM
         ALSL(L) = LOG(SPL(L))
      END DO
!
!HC WRITE IGDS OUT FOR WEIGHTMAKER TO READ IN AS KGDSIN
      if(me == 0)then
        print*,'writing out igds'
        igdout = 110
!        open(igdout,file='griddef.out',form='unformatted'
!     +  ,status='unknown')
        if(maptype == 1)THEN  ! Lambert conformal
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
        ELSE IF(MAPTYPE  ==  2)THEN  !Polar stereographic
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
          if (TRUELAT1 < 0.) THEN
            LAT = -60.
          else
            LAT = 60.
          end if

          CALL MSFPS (LAT,TRUELAT1*0.001,PSMAPF)

        ELSE IF(MAPTYPE == 3) THEN  !Mercator
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
        ELSE IF(MAPTYPE == 0 .OR. MAPTYPE == 203)THEN  !A STAGGERED E-GRID
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
         ELSE IF(MAPTYPE == 207)THEN !Rotated lat-lon grid
           write(flatlon,1001)ifhr
           open(112,file=trim(flatlon),form='formatted', &
             status='unknown')
           write(112,1002)LATSTART/1000,LONSTART/1000,&
             LATSE/1000,LONSE/1000,LATNW/1000,LONNW/1000,&
             LATLAST/1000,LONLAST/1000
     1001 format('latlons_corners.txt.f',I3.3)
     1002 format(4(I6,I7,X))
           close(112)
        END IF
      end if
!     
!

      RETURN
      END

!----------------------------------------------------------------------
!> @brief read_netcdf_3d_para() reads dynamics variables from UFS model output. 
!> 
!> @param[in] ncid integer netCDF ID.
!> @param[in] im integer Full longitude domain.
!> @param[in] jm integer Full latitude domain.
!> @param[in] ista integer Start longitude latitude on a task subdomain.
!> @param[in] ista_2l integer Start longitude -2 of the subdomain.
!> @param[in] iend integer End longitude on a task subdomain.
!> @param[in] iend_2u integer End longitude +2 of the subdomain.
!> @param[in] jsta integer Start latitude on a task subdomain.
!> @param[in] jsta_2l integer Start latitude -2 of the subdomain.
!> @param[in] jend integer End latitude on a task subdomain.
!> @param[in] jend_2u integer End latitude +2 of the subdomain.
!> @param[in] spval real Missing value defined in UPP.
!> @param[in] varname character Variable name in netCDF file.
!> @param[out] buf real Variable values.
!> @param[in] lm integer Model levels.
!----------------------------------------------------------------------
      subroutine read_netcdf_3d_para(ncid,im,jm,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
                 spval,varname,buf,lm)

      use netcdf
      use ctlblk_mod, only : me
      use params_mod, only : small
      implicit none
      INCLUDE "mpif.h"

      character(len=20),intent(in) :: varname
      real,intent(in)    :: spval
      integer,intent(in) :: ncid,im,jm,lm,jsta_2l,jend_2u,jsta,jend
      integer,intent(in) :: ista_2l,iend_2u,ista,iend
      real,intent(out)   :: buf(ista_2l:iend_2u,jsta_2l:jend_2u,lm)
      integer            :: varid,iret,ii,jj,i,j,l,kk
      integer            :: start(3), count(3), stride(3)
      real,parameter     :: spval_netcdf=9.99e+20
      real               :: fill_value

      iret = nf90_inq_varid(ncid,trim(varname),varid)
      if (iret /= 0) then
        if (me == 0) print*,VarName," not found -Assigned missing values"
!$omp parallel do private(i,j,l)
          do l=1,lm
            do j=jsta,jend
              do i=ista,iend
                buf(i,j,l)=spval
              enddo
            enddo
          enddo
      else
        iret = nf90_get_att(ncid,varid,"_FillValue",fill_value)
        if (iret /= 0) fill_value = spval_netcdf
        start = (/ista,jsta,1/)
        ii=iend-ista+1
        jj=jend-jsta+1
        count = (/ii,jj,lm/)
        iret = nf90_get_var(ncid,varid,buf(ista:iend,jsta:jend,1:lm),start=start,count=count)
        if (iret /= 0) then
          print*," iret /=0, Error in reading varid "
        endif
        do l=1,lm
          do j=jsta,jend
            do i=ista,iend
              if(abs(buf(i,j,l)-fill_value)<small)buf(i,j,l)=spval
            end do
          end do
        end do
      endif

      end subroutine read_netcdf_3d_para

!----------------------------------------------------------------------
!> @brief read_netcdf_2d_para() reads physics variables from UFS model output. 
!> 
!> @param[in] ncid integer netCDF ID.
!> @param[in] ista integer Start longitude latitude on a task subdomain.
!> @param[in] ista_2l integer Start longitude -2 of the subdomain.
!> @param[in] iend integer End longitude on a task subdomain.
!> @param[in] iend_2u integer End longitude +2 of the subdomain.
!> @param[in] jsta integer Start latitude on a task subdomain.
!> @param[in] jsta_2l integer Start latitude -2 of the subdomain.
!> @param[in] jend integer End latitude on a task subdomain.
!> @param[in] jend_2u integer End latitude +2 of the subdomain.
!> @param[in] spval real Missing value defined in UPP.
!> @param[in] varname character Variable name in netCDF file.
!> @param[out] buf real Variable values.
!----------------------------------------------------------------------

      subroutine read_netcdf_2d_para(ncid,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
                 spval,VarName,buf)

      use netcdf
      use ctlblk_mod, only : me
      use params_mod, only : small
      implicit none
      INCLUDE "mpif.h"

      character(len=20),intent(in) :: VarName
      real,intent(in)    :: spval
      integer,intent(in) :: ncid,jsta_2l,jend_2u,jsta,jend,ista_2l,iend_2u,ista,iend
      real,intent(out)   :: buf(ista_2l:iend_2u,jsta_2l:jend_2u)
      integer            :: varid,iret,ii,jj,i,j,l,kk
      integer            :: start(2), count(2)
      real,parameter     :: spval_netcdf=9.99e+20
      real               :: fill_value

      iret = nf90_inq_varid(ncid,trim(varname),varid)
      if (iret /= 0) then
        if (me==0) print*,VarName," not found -Assigned missing values"
!$omp parallel do private(i,j)
        do j=jsta,jend
          do i=ista,iend
            buf(i,j)=spval
          enddo
        enddo
      else
        iret = nf90_get_att(ncid,varid,"_FillValue",fill_value)
        if (iret /= 0) fill_value = spval_netcdf
        start = (/ista,jsta/)
        ii=iend-ista+1
        jj=jend-jsta+1
        count = (/ii,jj/)
        iret = nf90_get_var(ncid,varid,buf(ista:iend,jsta:jend),start=start,count=count)
        if (iret /= 0) then
          print*," iret /=0, Error in reading varid "
        endif
        do j=jsta,jend
          do i=ista,iend
            if(abs(buf(i,j)-fill_value)<small)buf(i,j)=spval
          end do
        end do
      endif

      end subroutine read_netcdf_2d_para
