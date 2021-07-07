!> @file
!
!> SUBPROGRAM:    MPI_FIRST   SET UP MESSGAE PASSING INFO
!!   PRGRMMR: TUCCILLO        ORG: IBM
!!
!! ABSTRACT:
!!     SETS UP MESSAGE PASSING INFO
!!
!! PROGRAM HISTORY LOG:
!!   14-12-01   WM LEWIS: ADDED ADDNL VARIABLES FOR SAT OUTPUT
!!   00-01-06  TUCCILLO - ORIGINAL
!!   01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!!   02-06-19  MIKE BALDWIN - WRF VERSION
!!   11-12-16  SARAH LU - MODIFIED TO INITIALIZE AEROSOL FIELDS
!!   12-01-07  SARAH LU - MODIFIED TO INITIALIZE AIR DENSITY/LAYER THICKNESS
!!
!! USAGE:    CALL MPI_FIRST
!!   INPUT ARGUMENT LIST:
!!
!!   OUTPUT ARGUMENT LIST:
!!
!!   OUTPUT FILES:
!!     STDOUT  - RUN TIME STANDARD OUT.
!!
!!   SUBPROGRAMS CALLED:
!!       PARA_RANGE
!!     UTILITIES:
!!       NONE
!!     LIBRARY:
!!       COMMON - CTLBLK.comm
!!
!!   ATTRIBUTES:
!!     LANGUAGE: FORTRAN
!!     MACHINE : IBM RS/6000 SP
!!
!@PROCESS NOEXTCHK
      SUBROUTINE MPI_FIRST()

!
      use vrbls4d, only: dust, salt, soot, waso, suso, pp25, pp10
      use vrbls3d, only: u, v, t, q, uh, vh, wh, pmid, pmidv, pint, alpint, zmid,      &
              zint, q2, omga, t_adj, ttnd, rswtt, rlwtt, exch_h, train, tcucn,         &
              el_pbl, cwm, f_ice, f_rain, f_rimef, qqw, qqi, qqr, qqs,qqg, qqni, qqnr, &
              extcof55, cfr, dbz, dbzr, dbzi, dbzc, mcvg, nlice, nrain, o3, vdifftt,   &
              tcucns, vdiffmois, dconvmois, sconvmois, nradtt, o3vdiff, o3prod,        &
              o3tndy, mwpv, unknown, vdiffzacce, zgdrag, cnvctummixing, vdiffmacce,    &
              mgdrag, cnvctvmmixing, ncnvctcfrac, cnvctumflx, cnvctdmflx, cnvctdetmflx,&
              cnvctzgdrag, cnvctmgdrag, icing_gfip, asy, ssa, duem, dusd, dudp,        &
              duwt, suem, susd, sudp, suwt, ocem, ocsd, ocdp, ocwt, bcem, bcsd,        &
              bcdp, bcwt, ssem, sssd, ssdp, sswt, ext, dpres, rhomid
      use vrbls2d, only: wspd10max, w_up_max, w_dn_max, w_mean, refd_max, up_heli_max, &
              prate_max, fprate_max, swupt,                                            &
              up_heli_max16, grpl_max, up_heli, up_heli16, ltg1_max, ltg2_max,         &
              up_heli_min, up_heli_min16, up_heli_max02, up_heli_min02, up_heli_max03, &
              up_heli_min03, rel_vort_max, rel_vort_max01, wspd10umax, wspd10vmax,     &
              refdm10c_max, hail_max2d, hail_maxk1, ltg3_max,rel_vort_maxhy1,          &
              nci_ltg, nca_ltg, nci_wq, nca_wq, nci_refd,                              &
              u10, v10, tshltr, qshltr, mrshltr, smstav, ssroff, bgroff,               &
              nca_refd, vegfrc, acsnow, acsnom, cmc, sst, qz0, thz0, uz0, vz0, qs, ths,&
              sno, snonc, snoavg, psfcavg, t10m, t10avg, akmsavg, akhsavg, u10max,     &
              v10max, u10h, v10h, akms, akhs, cuprec, acprec, ancprc, cuppt,           &
              rainc_bucket, rainnc_bucket, pcp_bucket, snow_bucket, qrmax, tmax,       &
              snownc, graupelnc, tsnow, qvg, qv2m, rswin, rlwin, rlwtoa, tg, sfcshx,   &
              fis, t500, cfracl, cfracm, cfrach, acfrst, acfrcv, hbot, potevp,         &
              sfclhx, htop, aswin, alwin, aswout, alwout, aswtoa, alwtoa, czen, czmean,&
              sigt4, rswout, radot, ncfrst, ncfrcv, smstot, pctsno, pshltr, th10,      &
              q10, sr, prec, subshx, snopcx, sfcuvx, sfcevp, z0, ustar, pblh, mixht,   &
              twbs, qwbs, sfcexc, grnflx, soiltb, z1000, slp, pslp, f, albedo, albase, &
              cldfra, cprate, cnvcfr, ivgtyp, hbotd, htopd, hbots, isltyp, htops,      &
              cldefi, islope, si, lspa, rswinc, vis, pd, mxsnal, epsr, sfcux,          &
              sfcvx, sfcuxi, sfcvxi, avgalbedo, avgcprate, avgprec, ptop, pbot, avgcfrach, avgcfracm,  &
              avgcfracl, avgtcdc, auvbin, auvbinc, ptopl, pbotl, ttopl, ptopm,         &
              pbotm, ttopm, ptoph, pboth, ttoph, sfcugs, sfcvgs, pblcfr, cldwork,      &
              gtaux, gtauy, mdltaux, mdltauy, runoff, maxtshltr, mintshltr,            &
              maxrhshltr, minrhshltr, dzice, alwinc, alwoutc, alwtoac, aswinc,         &
              aswoutc,aswtoac, aswintoa, smcwlt, suntime, fieldcapa, avisbeamswin,     &
              avisdiffswin, airbeamswin, airdiffswin, snowfall, dusmass, ducmass,      &
              dusmass25, susmass, sucmass, susmass25, sucmass25, ocsmass, occmass,     &
              ocsmass25, occmass25, bcsmass, bccmass, bcsmass25, bccmass25,            &
              sssmass, sscmass, sssmass25, sscmass25, ducmass25,                       &
              dustcb, sscb, bccb, occb, sulfcb, dustallcb, ssallcb,dustpm,sspm, pp25cb,&
              pp10cb, ti
      use soil, only:  smc, stc, sh2o, sldpth, rtdpth, sllevel
      use masks, only: htm, vtm, hbm2, sm, sice, lmh, gdlat, gdlon, dx, dy, lmv
      use ctlblk_mod, only: me, num_procs, jm, jsta, jend, jsta_m, jsta_m2,ista,iend ,          &
              jend_m, jend_m2, iup, idn, icnt, im, idsp, jsta_2l, jend_2u,            &
              jvend_2u, lm, lp1, jsta_2l, jend_2u, nsoil, nbin_du, nbin_ss,           &
              nbin_bc, nbin_oc, nbin_su,                                               &
              ISTA_M,IEND_M,ISTA_M2,IEND_M2,                                          &
              ISTA_2L, IEND_2U,IVEND_2U             

!
!     use params_mod
!- - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - 
      implicit none
     
!
      include 'mpif.h'
!
      integer ierr,i,jsx,jex,isx,iex
      integer isumm
      integer numx !number of subdomain in x direction
!
       isumm=0
       numx=2

      if ( me == 0 ) then
!        print *, ' NUM_PROCS = ',num_procs
      end if

      if ( num_procs > 1024 ) then
         print *, ' too many MPI tasks, max is 1024, stopping'
         call mpi_abort(MPI_COMM_WORLD,1,ierr)
         stop
      end if
!
!     error check
!
      if ( num_procs > JM/2 ) then
         print *, ' too many MPI tasks, max is ',jm/2,' stopping'
         call mpi_abort(MPI_COMM_WORLD,1,ierr)
         stop
      end if
!
!     global loop ranges
!
!      call para_range(1,jm,num_procs,me,jsta,jend)
! GWVX temporary documentation
!  para_range2 supports a 2D decomposition.  The rest of the post
!  supports 1D still and the call here is the special case where each
!  processor gets all of the longitudes in the latitude 1D subdomain
!  jsta:jend.  The X decomposition will be specified by the third
!  argument (currently 1) and the Y decoposition will be specified by
!  the fourth argument (currently all of the ranks)   When X is
!  subdivided the third and fourth arguments will have to be integral
!  factors of num_procs and on 5/27/21 I am still working out a general
!  way to do this if the user doesn't select the factors
     ! call para_range2(im,jm,1,num_procs,me,ista,iend,jsta,jend)
      call para_range2(im,jm,numx,num_procs/numx,me,ista,iend,jsta,jend)
      jsta_m  = jsta
      jsta_m2 = jsta
      jend_m  = jend
      jend_m2 = jend
      if ( me == 0 ) then
         jsta_m  = 2
         jsta_m2 = 3
      end if
      if ( me == num_procs - 1 ) then
         jend_m  = jm - 1
         jend_m2 = jm - 2
      end if
!
!     neighbors
!
      iup = me + 1
      idn = me - 1
      if ( me == 0 ) then
         idn = MPI_PROC_NULL
      end if
      if ( me == num_procs - 1 ) then
         iup = MPI_PROC_NULL
      end if
!
!     print *, ' ME, NUM_PROCS = ',me,num_procs
!     print *, ' ME, JSTA, JSTA_M, JSTA_M2 = ',me,jsta,jsta_m,jsta_m2
!     print *, ' ME, JEND, JEND_M, JEND_M2 = ',me,jend,jend_m,jend_m2
!     print *, ' ME, IUP, IDN = ',me,iup,idn
!
!     counts, disps for gatherv and scatterv
!
      do i = 0, num_procs - 1
!         call para_range(1,jm,num_procs,i,jsx,jex) 
         call para_range2(im,jm,numx,num_procs/numx,i,isx,iex,jsx,jex) 
!         icnt(i) = (jex-jsx+1)*im
         icnt(i) = (jex-jsx+1)*(iex-isx+1)
         
!         idsp(i) = (jsx-1)*im
          idsp(i)=isumm
          isumm=isumm+(jex-jsx+1)*(iex-isx+1)
         if ( me == 0 ) then
           print *, ' i, icnt(i),idsp(i) = ',i,icnt(i),      &
            idsp(i)
         end if
      end do
!
!     extraction limits -- set to two rows    
!
      jsta_2l = max(jsta - 2,  1 )
      jend_2u = min(jend + 2, jm )
! special for c-grid v
      jvend_2u = min(jend + 2, jm+1 )
! special for c-grid v
      ista_2l=max(ista-2,1)
      iend_2u=min(iend+2,im)
      ivend_2u = min(iend + 2, im+1 )
!     print *, ' me, jvend_2u = ',me,jvend_2u
!
!     allocate arrays
!
!
!     FROM VRBLS3D
!
      print *, ' me, jsta_2l, jend_2u = ',me,jsta_2l, jend_2u,  &
               'jvend_2u=',jvend_2u,'im=',im,'jm=',jm,'lm=',lm, &
               'lp1=',lp1
      write(*,'(A,5I10)') 'MPI_FIRST me,jsta,jend,ista,iend,=',me,jsta,jend,ista,iend

      end
