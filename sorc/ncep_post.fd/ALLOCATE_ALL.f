!> @file
!
!> SET UP MESSGAE PASSING INFO
!! @author TUCCILLO ORG: IBM
!!
!! PROGRAM HISTORY LOG:
!! -  00-01-06  TUCCILLO - ORIGINAL
!! -  01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!! -  02-06-19  MIKE BALDWIN - WRF VERSION
!! -  11-12-16  SARAH LU - MODIFIED TO INITIALIZE AEROSOL FIELDS
!! -  12-01-07  SARAH LU - MODIFIED TO INITIALIZE AIR DENSITY/LAYER THICKNESS
!! -  15-07-04  SARAH LU - MODIFIED TO INITIALIZE SCA
!! -  15-07-21  Jun Wang - Add scavenging for DU, SS, OC, BC, remove 
!!                        SU diagnostic fields
!! -  19-07-24  Li(Kate) Zhang - Merge and update NGAC UPP for FV3-Chem
!! -  19-11-23  Wen Meng - Add sea ice skin T
!! -  20-11-06  Jesse Meng - Add UPP_MATH module variables
!! -  21-04-06  Wen Meng - Initializing all allocated arrays
!! -  21-04-16  Wen Meng - Initializing aextc55 and extc55 as 0. These
!!                      two arrays are involved in GSL visibility computation.
!!
!!   OUTPUT FILES:
!!   - STDOUT  - RUN TIME STANDARD OUT.
!!
!!   SUBPROGRAMS CALLED:
!!     - para_range()
!!   LIBRARY:
!!     - COMMON - CTLBLK.comm
!!
      SUBROUTINE ALLOCATE_ALL()
!
      use vrbls4d
      use vrbls3d
      use vrbls2d
      use soil
      use masks
      use upp_math, only: ddvdx, ddudy, uuavg
!
      !use params_mod
      use ctlblk_mod
!- - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - 
      implicit none
!
      include 'mpif.h'
!
      integer ierr,jsx,jex
      integer i,j,l,k
! Allocate arrays
      allocate(u(im,jsta_2l:jend_2u,lm))
      allocate(v(im,jsta_2l:jvend_2u,lm))
      allocate(t(im,jsta_2l:jend_2u,lm))
! CHUANG ADD POTENTIAL TEMP BECAUSE WRF OUTPUT THETA
!      allocate(th(im,jsta_2l:jend_2u,lm))   
      allocate(q(im,jsta_2l:jend_2u,lm))
!      allocate(w(im,jsta_2l:jend_2u,lp1))
      allocate(uh(im,jsta_2l:jend_2u,lm))
      allocate(vh(im,jsta_2l:jend_2u,lm))
      allocate(wh(im,jsta_2l:jend_2u,lm))
      allocate(pmid(im,jsta_2l:jend_2u,lm))
      allocate(pmidv(im,jsta_2l:jend_2u,lm))
      allocate(pint(im,jsta_2l:jend_2u,lp1))
      allocate(alpint(im,jsta_2l:jend_2u,lp1))
      allocate(zmid(im,jsta_2l:jend_2u,lm))
      allocate(zint(im,jsta_2l:jend_2u,lp1))
!      allocate(rainw(im,jsta_2l:jend_2u,lm))
      allocate(q2(im,jsta_2l:jend_2u,lm))
      allocate(omga(im,jsta_2l:jend_2u,lm))
      allocate(dpres(im,jsta_2l:jend_2u,lm))
      allocate(T_ADJ(im,jsta_2l:jend_2u,lm))
      allocate(ttnd(im,jsta_2l:jend_2u,lm))
      allocate(rswtt(im,jsta_2l:jend_2u,lm))
      allocate(rlwtt(im,jsta_2l:jend_2u,lm))
      allocate(exch_h(im,jsta_2l:jend_2u,lm)) 
      allocate(train(im,jsta_2l:jend_2u,lm))
      allocate(tcucn(im,jsta_2l:jend_2u,lm))
      allocate(EL_PBL(im,jsta_2l:jend_2u,lm))

!Initialization
!$omp parallel do private(i,j,l)
      do l=1,lm
        do j=jsta_2l,jend_2u
          do i=1,lm
            u(i,j,l)=0.
            v(i,j,l)=0.
            t(i,j,l)=spval
            q(i,j,l)=spval
            uh(i,j,l)=spval
            vh(i,j,l)=spval
            wh(i,j,l)=spval
            pmid(i,j,l)=spval
            pmidv(i,j,l)=spval
            zmid(i,j,l)=spval
            q2(i,j,l)=spval
            omga(i,j,l)=spval
            dpres(i,j,l)=spval
            T_ADJ(i,j,l)=spval 
            ttnd(i,j,l)=spval 
            rswtt(i,j,l)=spval 
            rlwtt(i,j,l)=spval 
            exch_h(i,j,l)=spval 
            train(i,j,l)=spval 
            tcucn(i,j,l)=spval 
            EL_PBL(i,j,l)=spval 
          enddo
        enddo
      enddo
!$omp parallel do private(i,j,l)
      do l=1,lp1
        do j=jsta_2l,jend_2u
          do i=1,lm
            pint(i,j,l)=spval
            alpint(i,j,l)=spval
            zint(i,j,l)=spval
          enddo
        enddo
      enddo

!     MP FIELD   
      allocate(cwm(im,jsta_2l:jend_2u,lm))
      allocate(F_ice(im,jsta_2l:jend_2u,lm))
      allocate(F_rain(im,jsta_2l:jend_2u,lm))
      allocate(F_RimeF(im,jsta_2l:jend_2u,lm))
      allocate(QQW(im,jsta_2l:jend_2u,lm))
      allocate(QRIMEF(im,jsta_2l:jend_2u,lm))
      allocate(QQI(im,jsta_2l:jend_2u,lm))
      allocate(QQR(im,jsta_2l:jend_2u,lm))
      allocate(QQS(im,jsta_2l:jend_2u,lm))
      allocate(QQG(im,jsta_2l:jend_2u,lm))
      allocate(QQNW(im,jsta_2l:jend_2u,lm))
      allocate(QQNI(im,jsta_2l:jend_2u,lm))
      allocate(QQNR(im,jsta_2l:jend_2u,lm))
      allocate(QQNWFA(im,jsta_2l:jend_2u,lm))
      allocate(QQNIFA(im,jsta_2l:jend_2u,lm))
      allocate(TAOD5503D(im,jsta_2l:jend_2u,lm))
      allocate(AEXTC55(im,jsta_2l:jend_2u,lm))
      allocate(EXTCOF55(im,jsta_2l:jend_2u,lm))
      allocate(QC_BL(im,jsta_2l:jend_2u,lm))
      allocate(CFR(im,jsta_2l:jend_2u,lm))
      allocate(CFR_RAW(im,jsta_2l:jend_2u,lm))
      allocate(DBZ(im,jsta_2l:jend_2u,lm))
      allocate(DBZR(im,jsta_2l:jend_2u,lm))
      allocate(DBZI(im,jsta_2l:jend_2u,lm))
      allocate(DBZC(im,jsta_2l:jend_2u,lm))
      allocate(mcvg(im,jsta_2l:jend_2u,lm))
      allocate(NLICE(im,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
      do l=1,lm
      do j=jsta_2l,jend_2u
          do i=1,lm 
            cwm(i,j,l)=spval
            F_ice(i,j,l)=spval
            F_rain(i,j,l)=spval
            F_RimeF(i,j,l)=spval
            QQW(i,j,l)=spval
            QRIMEF(i,j,l)=spval
            QQI(i,j,l)=spval
            QQR(i,j,l)=spval
            QQS(i,j,l)=spval
            QQG(i,j,l)=spval
            QQNW(i,j,l)=spval
            QQNI(i,j,l)=spval
            QQNR(i,j,l)=spval
            QQNWFA(i,j,l)=spval
            QQNIFA(i,j,l)=spval
            TAOD5503D(i,j,l)=spval
            AEXTC55(i,j,l)=0.
            EXTCOF55(i,j,l)=0.
            QC_BL(i,j,l)=spval
            CFR(i,j,l)=spval
            CFR_RAW(i,j,l)=spval
            DBZ(i,j,l)=spval
            DBZR(i,j,l)=spval
            DBZI(i,j,l)=spval
            DBZC(i,j,l)=spval
            mcvg(i,j,l)=spval
            NLICE(i,j,l)=spval
          enddo
        enddo
      enddo
!     Wm Lewis: added 
      allocate(NRAIN(im,jsta_2l:jend_2u,lm))
      allocate(radius_cloud(im,jsta_2l:jend_2u,lm))
      allocate(radius_ice(im,jsta_2l:jend_2u,lm))
      allocate(radius_snow(im,jsta_2l:jend_2u,lm))
! KRS: HWRF Addition for thompson reflectivity
! or non-ferrier physics. wrf-derived
      allocate(REFL_10CM(im,jsta_2l:jend_2u,lm))
!GFS FIELD
      allocate(o3(im,jsta_2l:jend_2u,lm))
      allocate(o(im,jsta_2l:jend_2u,lm))
      allocate(o2(im,jsta_2l:jend_2u,lm))
      allocate(tcucns(im,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
      do l=1,lm
      do j=jsta_2l,jend_2u
          do i=1,lm
            NRAIN(i,j,l)=spval
            radius_cloud(i,j,l)=spval
            radius_ice(i,j,l)=spval
            radius_snow(i,j,l)=spval
            REFL_10CM(i,j,l)=spval
            o3(i,j,l)=spval
            o(i,j,l)=spval
            o2(i,j,l)=spval
            tcucns(i,j,l)=spval
          enddo
        enddo
      enddo
! Add GFS d3d fields
      if (me == 0) print *,' d3d_on=',d3d_on
      if (d3d_on) then
        allocate(vdifftt(im,jsta_2l:jend_2u,lm))
!       allocate(tcucns(im,jsta_2l:jend_2u,lm))
        allocate(vdiffmois(im,jsta_2l:jend_2u,lm))
        allocate(dconvmois(im,jsta_2l:jend_2u,lm))
        allocate(sconvmois(im,jsta_2l:jend_2u,lm))
        allocate(nradtt(im,jsta_2l:jend_2u,lm))
        allocate(o3vdiff(im,jsta_2l:jend_2u,lm))
        allocate(o3prod(im,jsta_2l:jend_2u,lm))
        allocate(o3tndy(im,jsta_2l:jend_2u,lm))
        allocate(mwpv(im,jsta_2l:jend_2u,lm))
        allocate(unknown(im,jsta_2l:jend_2u,lm))
        allocate(vdiffzacce(im,jsta_2l:jend_2u,lm))
        allocate(zgdrag(im,jsta_2l:jend_2u,lm))
        allocate(cnvctummixing(im,jsta_2l:jend_2u,lm))
        allocate(vdiffmacce(im,jsta_2l:jend_2u,lm))
        allocate(mgdrag(im,jsta_2l:jend_2u,lm))
        allocate(cnvctvmmixing(im,jsta_2l:jend_2u,lm))
        allocate(ncnvctcfrac(im,jsta_2l:jend_2u,lm))
        allocate(cnvctumflx(im,jsta_2l:jend_2u,lm))
        allocate(cnvctdmflx(im,jsta_2l:jend_2u,lm))
        allocate(cnvctdetmflx(im,jsta_2l:jend_2u,lm))
        allocate(cnvctzgdrag(im,jsta_2l:jend_2u,lm))
        allocate(cnvctmgdrag(im,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
        do l=1,lm
          do j=jsta_2l,jend_2u
            do i=1,im
              vdifftt(i,j,l)=spval
              vdiffmois(i,j,l)=spval
              dconvmois(i,j,l)=spval
              sconvmois(i,j,l)=spval
              nradtt(i,j,l)=spval
              o3vdiff(i,j,l)=spval
              o3prod(i,j,l)=spval
              o3tndy(i,j,l)=spval
              mwpv(i,j,l)=spval
              unknown(i,j,l)=spval
              vdiffzacce(i,j,l)=spval
              zgdrag(i,j,l)=spval
              cnvctummixing(i,j,l)=spval
              vdiffmacce(i,j,l)=spval
              mgdrag(i,j,l)=spval
              cnvctvmmixing(i,j,l)=spval
              ncnvctcfrac(i,j,l)=spval
              cnvctumflx(i,j,l)=spval
              cnvctdmflx(i,j,l)=spval
              cnvctdetmflx(i,j,l)=spval
              cnvctzgdrag(i,j,l)=spval
              cnvctmgdrag(i,j,l)=spval
            enddo
          enddo
        enddo 
      endif
!
      allocate(htm(im,jsta_2l:jend_2u,lm))
      allocate(vtm(im,jsta_2l:jend_2u,lm))
! add GFIP ICING
      allocate(icing_gfip(im,jsta_2l:jend_2u,lm))        
      allocate(icing_gfis(im,jsta_2l:jend_2u,lm))        
!
! add GTG turbulence
      allocate(catedr(im,jsta_2l:jend_2u,lm))
      allocate(mwt(im,jsta_2l:jend_2u,lm))
      allocate(gtg(im,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
      do l=1,lm
        do j=jsta_2l,jend_2u
          do i=1,im
            htm(i,j,l)=spval
            vtm(i,j,l)=spval
            icing_gfip(i,j,l)=spval
            icing_gfis(i,j,l)=spval
            catedr(i,j,l)=spval
            mwt(i,j,l)=spval
            gtg(i,j,l)=spval
          enddo
        enddo
      enddo
!
!     FROM SOIL
!
      allocate(smc(im,jsta_2l:jend_2u,nsoil))
      allocate(stc(im,jsta_2l:jend_2u,nsoil))
      allocate(sh2o(im,jsta_2l:jend_2u,nsoil))
      allocate(SLDPTH(NSOIL))
      allocate(RTDPTH(NSOIL))
      allocate(SLLEVEL(NSOIL))
!Initialization
!$omp parallel do private(i,j,l)
      do l=1,nsoil
        do j=jsta_2l,jend_2u
          do i=1,im
            smc(i,j,l)=spval
            stc(i,j,l)=spval
            sh2o(i,j,l)=spval
          enddo
        enddo
      enddo
!$omp parallel do private(i)
      do i=1,NSOIL
        SLDPTH(i)=spval
        RTDPTH(i)=spval
        SLLEVEL(i)=spval
      enddo
!
!     FROM VRBLS2D
!
! SRD
      allocate(wspd10max(im,jsta_2l:jend_2u))
      allocate(w_up_max(im,jsta_2l:jend_2u))
      allocate(w_dn_max(im,jsta_2l:jend_2u))
      allocate(w_mean(im,jsta_2l:jend_2u))
      allocate(refd_max(im,jsta_2l:jend_2u))
      allocate(prate_max(im,jsta_2l:jend_2u))
      allocate(fprate_max(im,jsta_2l:jend_2u))
      allocate(up_heli_max(im,jsta_2l:jend_2u))
      allocate(up_heli_max16(im,jsta_2l:jend_2u))
      allocate(up_heli_min(im,jsta_2l:jend_2u))
      allocate(up_heli_min16(im,jsta_2l:jend_2u))
      allocate(up_heli_max02(im,jsta_2l:jend_2u))
      allocate(up_heli_min02(im,jsta_2l:jend_2u))
      allocate(up_heli_max03(im,jsta_2l:jend_2u))
      allocate(up_heli_min03(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          wspd10max(i,j)=spval
          w_up_max(i,j)=spval
          w_dn_max(i,j)=spval
          w_mean(i,j)=spval
          refd_max(i,j)=spval
          prate_max(i,j)=spval
          fprate_max(i,j)=spval
          up_heli_max(i,j)=spval
          up_heli_max16(i,j)=spval
          up_heli_min(i,j)=spval
          up_heli_min16(i,j)=spval
          up_heli_max02(i,j)=spval
          up_heli_min02(i,j)=spval
          up_heli_max03(i,j)=spval
          up_heli_min03(i,j)=spval
        enddo
      enddo
      allocate(rel_vort_max(im,jsta_2l:jend_2u))
      allocate(rel_vort_max01(im,jsta_2l:jend_2u))
      allocate(rel_vort_maxhy1(im,jsta_2l:jend_2u))
      allocate(wspd10umax(im,jsta_2l:jend_2u))
      allocate(wspd10vmax(im,jsta_2l:jend_2u))
      allocate(refdm10c_max(im,jsta_2l:jend_2u))
      allocate(hail_max2d(im,jsta_2l:jend_2u))
      allocate(hail_maxk1(im,jsta_2l:jend_2u))
      allocate(hail_maxhailcast(im,jsta_2l:jend_2u))
      allocate(grpl_max(im,jsta_2l:jend_2u))
      allocate(up_heli(im,jsta_2l:jend_2u))
      allocate(up_heli16(im,jsta_2l:jend_2u))
      allocate(ltg1_max(im,jsta_2l:jend_2u))
      allocate(ltg2_max(im,jsta_2l:jend_2u))
      allocate(ltg3_max(im,jsta_2l:jend_2u))
      allocate(nci_ltg(im,jsta_2l:jend_2u))
      allocate(nca_ltg(im,jsta_2l:jend_2u))
      allocate(nci_wq(im,jsta_2l:jend_2u))
      allocate(nca_wq(im,jsta_2l:jend_2u))
      allocate(nci_refd(im,jsta_2l:jend_2u))
      allocate(nca_refd(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          rel_vort_max(i,j)=spval
          rel_vort_max01(i,j)=spval
          rel_vort_maxhy1(i,j)=spval
          wspd10umax(i,j)=spval
          wspd10vmax(i,j)=spval
          refdm10c_max(i,j)=spval
          hail_max2d(i,j)=spval
          hail_maxk1(i,j)=spval
          hail_maxhailcast(i,j)=spval
          grpl_max(i,j)=spval
          up_heli(i,j)=spval
          up_heli16(i,j)=spval
          ltg1_max(i,j)=spval
          ltg2_max(i,j)=spval
          ltg3_max(i,j)=spval
          nci_ltg(i,j)=spval
          nca_ltg(i,j)=spval
          nci_wq(i,j)=spval
          nca_wq(i,j)=spval
          nci_refd(i,j)=spval
          nca_refd(i,j)=spval
        enddo
      enddo
! SRD
! CRA
      allocate(REF_10CM(im,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
      do l=1,lm
        do j=jsta_2l,jend_2u
          do i=1,im
            REF_10CM(i,j,l)=spval
          enddo
        enddo
      enddo
      allocate(REFC_10CM(im,jsta_2l:jend_2u))
      allocate(REF1KM_10CM(im,jsta_2l:jend_2u))
      allocate(REF4KM_10CM(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          REFC_10CM(i,j)=spval
          REF1KM_10CM(i,j)=spval
          REF4KM_10CM(i,j)=spval
        enddo
      enddo
! CRA
      allocate(u10(im,jsta_2l:jend_2u))
      allocate(v10(im,jsta_2l:jend_2u))
      allocate(tshltr(im,jsta_2l:jend_2u))
      allocate(qshltr(im,jsta_2l:jend_2u))
      allocate(mrshltr(im,jsta_2l:jend_2u))
      allocate(smstav(im,jsta_2l:jend_2u))
      allocate(ssroff(im,jsta_2l:jend_2u))
      allocate(bgroff(im,jsta_2l:jend_2u))
      allocate(vegfrc(im,jsta_2l:jend_2u))
      allocate(shdmin(im,jsta_2l:jend_2u))
      allocate(shdmax(im,jsta_2l:jend_2u))
      allocate(lai(im,jsta_2l:jend_2u))
      allocate(acsnow(im,jsta_2l:jend_2u))
      allocate(acgraup(im,jsta_2l:jend_2u))
      allocate(acfrain(im,jsta_2l:jend_2u))
      allocate(acsnom(im,jsta_2l:jend_2u))
      allocate(cmc(im,jsta_2l:jend_2u))
      allocate(sst(im,jsta_2l:jend_2u))
      allocate(qz0(im,jsta_2l:jend_2u))
      allocate(thz0(im,jsta_2l:jend_2u))
      allocate(uz0(im,jsta_2l:jend_2u))
      allocate(vz0(im,jsta_2l:jend_2u))
      allocate(qs(im,jsta_2l:jend_2u))
      allocate(ths(im,jsta_2l:jend_2u))
      allocate(sno(im,jsta_2l:jend_2u))
      allocate(snonc(im,jsta_2l:jend_2u))
      allocate(ti(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          u10(i,j)=spval
          v10(i,j)=spval
          tshltr(i,j)=spval
          qshltr(i,j)=spval
          mrshltr(i,j)=spval
          smstav(i,j)=spval
          ssroff(i,j)=spval
          bgroff(i,j)=spval
          vegfrc(i,j)=spval
          shdmin(i,j)=spval
          shdmax(i,j)=spval
          lai(i,j)=spval
          acsnow(i,j)=spval
          acgraup(i,j)=spval
          acfrain(i,j)=spval
          acsnom(i,j)=spval
          cmc(i,j)=spval
          sst(i,j)=spval
          qz0(i,j)=spval
          thz0(i,j)=spval
          uz0(i,j)=spval
          vz0(i,j)=spval
          qs(i,j)=spval
          ths(i,j)=spval
          sno(i,j)=spval
          snonc(i,j)=spval
          ti(i,j)=spval
        enddo
      enddo
! Time-averaged fileds
      allocate(u10mean(im,jsta_2l:jend_2u))
      allocate(v10mean(im,jsta_2l:jend_2u))
      allocate(spduv10mean(im,jsta_2l:jend_2u))
      allocate(swradmean(im,jsta_2l:jend_2u))
      allocate(swnormmean(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          u10mean(i,j)=spval
          v10mean(i,j)=spval
          spduv10mean(i,j)=spval
          swradmean(i,j)=spval
          swnormmean(i,j)=spval
        enddo
      enddo
!NAMstart
      allocate(snoavg(im,jsta_2l:jend_2u))
      allocate(psfcavg(im,jsta_2l:jend_2u))
      allocate(t10m(im,jsta_2l:jend_2u))
      allocate(t10avg(im,jsta_2l:jend_2u))
      allocate(akmsavg(im,jsta_2l:jend_2u))
      allocate(akhsavg(im,jsta_2l:jend_2u))
      allocate(u10max(im,jsta_2l:jend_2u))
      allocate(v10max(im,jsta_2l:jend_2u))
      allocate(u10h(im,jsta_2l:jend_2u))
      allocate(v10h(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          snoavg(i,j)=spval
          psfcavg(i,j)=spval
          t10m(i,j)=spval
          t10avg(i,j)=spval
          akmsavg(i,j)=spval
          akhsavg(i,j)=spval
          u10max(i,j)=spval
          v10max(i,j)=spval
          u10h(i,j)=spval
          v10h(i,j)=spval
        enddo
      enddo
!NAMend
      allocate(akms(im,jsta_2l:jend_2u))
      allocate(akhs(im,jsta_2l:jend_2u))
      allocate(cuprec(im,jsta_2l:jend_2u))
      allocate(acprec(im,jsta_2l:jend_2u))
      allocate(ancprc(im,jsta_2l:jend_2u))
      allocate(cuppt(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          akms(i,j)=spval
          akhs(i,j)=spval
          cuprec(i,j)=spval
          acprec(i,j)=spval
          ancprc(i,j)=spval
          cuppt(i,j)=spval
        enddo
      enddo
! GSDstart
      allocate(rainc_bucket(im,jsta_2l:jend_2u))
      allocate(rainc_bucket1(im,jsta_2l:jend_2u))
      allocate(rainnc_bucket(im,jsta_2l:jend_2u))
      allocate(rainnc_bucket1(im,jsta_2l:jend_2u))
      allocate(pcp_bucket(im,jsta_2l:jend_2u))
      allocate(pcp_bucket1(im,jsta_2l:jend_2u))
      allocate(snow_bucket(im,jsta_2l:jend_2u))
      allocate(snow_bucket1(im,jsta_2l:jend_2u))
      allocate(graup_bucket(im,jsta_2l:jend_2u))
      allocate(graup_bucket1(im,jsta_2l:jend_2u))
      allocate(qrmax(im,jsta_2l:jend_2u))
      allocate(tmax(im,jsta_2l:jend_2u))
      allocate(snownc(im,jsta_2l:jend_2u))
      allocate(graupelnc(im,jsta_2l:jend_2u))
      allocate(tsnow(im,jsta_2l:jend_2u))
      allocate(qvg(im,jsta_2l:jend_2u))
      allocate(qv2m(im,jsta_2l:jend_2u))
      allocate(qvl1(im,jsta_2l:jend_2u))
      allocate(snfden(im,jsta_2l:jend_2u))
      allocate(sndepac(im,jsta_2l:jend_2u))
      allocate(int_smoke(im,jsta_2l:jend_2u))
      allocate(mean_frp(im,jsta_2l:jend_2u))
      allocate(int_aod(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          rainc_bucket(i,j)=spval
          rainc_bucket1(i,j)=spval
          rainnc_bucket(i,j)=spval
          rainnc_bucket1(i,j)=spval
          pcp_bucket(i,j)=spval
          pcp_bucket1(i,j)=spval
          snow_bucket(i,j)=spval
          snow_bucket1(i,j)=spval
          graup_bucket(i,j)=spval
          graup_bucket1(i,j)=spval
          qrmax(i,j)=spval
          tmax(i,j)=spval
          snownc(i,j)=spval
          graupelnc(i,j)=spval
          tsnow(i,j)=spval
          qvg(i,j)=spval
          qv2m(i,j)=spval
          qvl1(i,j)=spval
          snfden(i,j)=spval
          sndepac(i,j)=spval
          int_smoke(i,j)=spval
          mean_frp(i,j)=spval
          int_aod(i,j)=spval
        enddo
      enddo
      allocate(smoke(im,jsta_2l:jend_2u,lm,nbin_sm))
!$omp parallel do private(i,j,l,k)
      do k=1,nbin_sm
        do l=1,lm
          do j=jsta_2l,jend_2u
            do i=1,im
              smoke(i,j,l,k)=spval
            enddo
          enddo
        enddo
      enddo
! GSDend
      allocate(rswin(im,jsta_2l:jend_2u))
      allocate(swddni(im,jsta_2l:jend_2u))
      allocate(swddif(im,jsta_2l:jend_2u))
      allocate(swdnbc(im,jsta_2l:jend_2u))
      allocate(swddnic(im,jsta_2l:jend_2u))
      allocate(swddifc(im,jsta_2l:jend_2u))
      allocate(swupbc(im,jsta_2l:jend_2u))
      allocate(swupt(im,jsta_2l:jend_2u))
      allocate(taod5502d(im,jsta_2l:jend_2u))
      allocate(aerasy2d(im,jsta_2l:jend_2u))
      allocate(aerssa2d(im,jsta_2l:jend_2u))
      allocate(lwp(im,jsta_2l:jend_2u))
      allocate(iwp(im,jsta_2l:jend_2u))
      allocate(rlwin(im,jsta_2l:jend_2u))
      allocate(lwdnbc(im,jsta_2l:jend_2u))
      allocate(lwupbc(im,jsta_2l:jend_2u))
      allocate(rlwtoa(im,jsta_2l:jend_2u))
      allocate(rswtoa(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          rswin(i,j)=spval
          swddni(i,j)=spval
          swddif(i,j)=spval
          swdnbc(i,j)=spval
          swddnic(i,j)=spval
          swddifc(i,j)=spval
          swupbc(i,j)=spval
          swupt(i,j)=spval
          taod5502d(i,j)=spval
          aerasy2d(i,j)=spval
          aerssa2d(i,j)=spval
          lwp(i,j)=spval
          iwp(i,j)=spval
          rlwin(i,j)=spval
          lwdnbc(i,j)=spval
          lwupbc(i,j)=spval
          rlwtoa(i,j)=spval
          rswtoa(i,j)=spval
        enddo
      enddo
      allocate(tg(im,jsta_2l:jend_2u))
      allocate(sfcshx(im,jsta_2l:jend_2u))
      allocate(sfclhx(im,jsta_2l:jend_2u))
      allocate(fis(im,jsta_2l:jend_2u))
      allocate(t500(im,jsta_2l:jend_2u))
      allocate(t700(im,jsta_2l:jend_2u))
      allocate(z500(im,jsta_2l:jend_2u))
      allocate(z700(im,jsta_2l:jend_2u))
      allocate(teql(im,jsta_2l:jend_2u))
      allocate(ieql(im,jsta_2l:jend_2u))
      allocate(cfracl(im,jsta_2l:jend_2u))
      allocate(cfracm(im,jsta_2l:jend_2u))
      allocate(cfrach(im,jsta_2l:jend_2u))
      allocate(acfrst(im,jsta_2l:jend_2u))
      allocate(acfrcv(im,jsta_2l:jend_2u))
      allocate(hbot(im,jsta_2l:jend_2u))
      allocate(htop(im,jsta_2l:jend_2u))
      allocate(aswin(im,jsta_2l:jend_2u))
      allocate(alwin(im,jsta_2l:jend_2u))
      allocate(aswout(im,jsta_2l:jend_2u))
      allocate(alwout(im,jsta_2l:jend_2u))
      allocate(aswtoa(im,jsta_2l:jend_2u))
      allocate(alwtoa(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          tg(i,j)=spval
          sfcshx(i,j)=spval
          sfclhx(i,j)=spval
          fis(i,j)=spval
          t500(i,j)=spval
          t700(i,j)=spval
          z700(i,j)=spval
          teql(i,j)=spval
          ieql(i,j)=spval
          cfracl(i,j)=spval
          cfracm(i,j)=spval
          cfrach(i,j)=spval
          acfrst(i,j)=spval
          acfrcv(i,j)=spval
          hbot(i,j)=spval
          htop(i,j)=spval
          aswin(i,j)=spval
          alwin(i,j)=spval
          aswout(i,j)=spval
          alwout(i,j)=spval
          aswtoa(i,j)=spval
          alwtoa(i,j)=spval
        enddo
      enddo
      allocate(czen(im,jsta_2l:jend_2u))
      allocate(czmean(im,jsta_2l:jend_2u))
      allocate(sigt4(im,jsta_2l:jend_2u))
      allocate(rswout(im,jsta_2l:jend_2u))
      allocate(radot(im,jsta_2l:jend_2u))
      allocate(ncfrst(im,jsta_2l:jend_2u))  ! real
      allocate(ncfrcv(im,jsta_2l:jend_2u))  ! real
      allocate(smstot(im,jsta_2l:jend_2u))
      allocate(pctsno(im,jsta_2l:jend_2u))
      allocate(pshltr(im,jsta_2l:jend_2u))
      allocate(th10(im,jsta_2l:jend_2u))
      allocate(q10(im,jsta_2l:jend_2u))
      allocate(sr(im,jsta_2l:jend_2u))
      allocate(prec(im,jsta_2l:jend_2u))
      allocate(subshx(im,jsta_2l:jend_2u))
      allocate(snopcx(im,jsta_2l:jend_2u))
      allocate(sfcuvx(im,jsta_2l:jend_2u))
      allocate(sfcevp(im,jsta_2l:jend_2u))
      allocate(potevp(im,jsta_2l:jend_2u))
      allocate(z0(im,jsta_2l:jend_2u))
      allocate(ustar(im,jsta_2l:jend_2u))
      allocate(pblh(im,jsta_2l:jend_2u))
      allocate(pblhgust(im,jsta_2l:jend_2u))
      allocate(mixht(im,jsta_2l:jend_2u))
      allocate(twbs(im,jsta_2l:jend_2u))
      allocate(qwbs(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          czen(i,j)=spval
          czmean(i,j)=spval
          sigt4(i,j)=spval
          rswout(i,j)=spval
          radot(i,j)=spval
          ncfrst(i,j)=spval
          ncfrcv(i,j)=spval
          smstot(i,j)=spval
          pctsno(i,j)=spval
          pshltr(i,j)=spval
          th10(i,j)=spval
          q10(i,j)=spval
          sr(i,j)=spval
          prec(i,j)=spval
          subshx(i,j)=spval
          snopcx(i,j)=spval
          sfcuvx(i,j)=spval
          sfcevp(i,j)=spval
          potevp(i,j)=spval
          z0(i,j)=spval
          ustar(i,j)=spval
          pblh(i,j)=spval
          pblhgust(i,j)=spval
          mixht(i,j)=spval
          twbs(i,j)=spval
          qwbs(i,j)=spval
        enddo
      enddo
      allocate(sfcexc(im,jsta_2l:jend_2u))
      allocate(grnflx(im,jsta_2l:jend_2u))
      allocate(soiltb(im,jsta_2l:jend_2u))
      allocate(z1000(im,jsta_2l:jend_2u))
      allocate(slp(im,jsta_2l:jend_2u))
      allocate(pslp(im,jsta_2l:jend_2u))
      allocate(f(im,jsta_2l:jend_2u))
      allocate(albedo(im,jsta_2l:jend_2u))
      allocate(albase(im,jsta_2l:jend_2u))
      allocate(cldfra(im,jsta_2l:jend_2u))
      allocate(cprate(im,jsta_2l:jend_2u))
      allocate(cnvcfr(im,jsta_2l:jend_2u))
      allocate(ivgtyp(im,jsta_2l:jend_2u))
      allocate(isltyp(im,jsta_2l:jend_2u))
      allocate(hbotd(im,jsta_2l:jend_2u))
      allocate(htopd(im,jsta_2l:jend_2u))
      allocate(hbots(im,jsta_2l:jend_2u))
      allocate(htops(im,jsta_2l:jend_2u))
      allocate(cldefi(im,jsta_2l:jend_2u))
      allocate(islope(im,jsta_2l:jend_2u))
      allocate(si(im,jsta_2l:jend_2u))
      allocate(lspa(im,jsta_2l:jend_2u))
      allocate(rswinc(im,jsta_2l:jend_2u))
      allocate(vis(im,jsta_2l:jend_2u))
      allocate(pd(im,jsta_2l:jend_2u))
      allocate(mxsnal(im,jsta_2l:jend_2u))
      allocate(epsr(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          sfcexc(i,j)=spval
          grnflx(i,j)=spval
          soiltb(i,j)=spval
          z1000(i,j)=spval
          slp(i,j)=spval
          pslp(i,j)=spval
          f(i,j)=spval
          albedo(i,j)=spval
          albase(i,j)=spval
          cldfra(i,j)=spval
          cprate(i,j)=spval
          cnvcfr(i,j)=spval
          ivgtyp(i,j)=spval
          isltyp(i,j)=spval
          hbotd(i,j)=spval
          htopd(i,j)=spval
          hbots(i,j)=spval
          htops(i,j)=spval
          cldefi(i,j)=spval
          islope(i,j)=spval
          si(i,j)=spval
          lspa(i,j)=spval
          rswinc(i,j)=spval
          vis(i,j)=spval
          pd(i,j)=spval
          mxsnal(i,j)=spval
          epsr(i,j)=spval
        enddo
      enddo
! add GFS fields
      allocate(sfcux(im,jsta_2l:jend_2u))
      allocate(sfcvx(im,jsta_2l:jend_2u))
      allocate(sfcuxi(im,jsta_2l:jend_2u))
      allocate(sfcvxi(im,jsta_2l:jend_2u))
      allocate(avgalbedo(im,jsta_2l:jend_2u))
      allocate(avgcprate(im,jsta_2l:jend_2u))
      allocate(avgprec(im,jsta_2l:jend_2u))
      allocate(avgprec_cont(im,jsta_2l:jend_2u))
      allocate(avgcprate_cont(im,jsta_2l:jend_2u))
      allocate(ptop(im,jsta_2l:jend_2u))
      allocate(pbot(im,jsta_2l:jend_2u))
      allocate(avgcfrach(im,jsta_2l:jend_2u))
      allocate(avgcfracm(im,jsta_2l:jend_2u))
      allocate(avgcfracl(im,jsta_2l:jend_2u))
      allocate(avgtcdc(im,jsta_2l:jend_2u))
      allocate(auvbin(im,jsta_2l:jend_2u))
      allocate(auvbinc(im,jsta_2l:jend_2u))
      allocate(ptopl(im,jsta_2l:jend_2u))
      allocate(pbotl(im,jsta_2l:jend_2u))
      allocate(Ttopl(im,jsta_2l:jend_2u))
      allocate(ptopm(im,jsta_2l:jend_2u))
      allocate(pbotm(im,jsta_2l:jend_2u))
      allocate(Ttopm(im,jsta_2l:jend_2u))
      allocate(ptoph(im,jsta_2l:jend_2u))
      allocate(pboth(im,jsta_2l:jend_2u))
      allocate(Ttoph(im,jsta_2l:jend_2u))
      allocate(sfcugs(im,jsta_2l:jend_2u))
      allocate(sfcvgs(im,jsta_2l:jend_2u))
      allocate(pblcfr(im,jsta_2l:jend_2u))
      allocate(cldwork(im,jsta_2l:jend_2u))
      allocate(gtaux(im,jsta_2l:jend_2u))
      allocate(gtauy(im,jsta_2l:jend_2u))
      allocate(cd10(im,jsta_2l:jend_2u))
      allocate(ch10(im,jsta_2l:jend_2u))
      allocate(mdltaux(im,jsta_2l:jend_2u))
      allocate(mdltauy(im,jsta_2l:jend_2u))
      allocate(runoff(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          sfcux(i,j)=spval
          sfcvx(i,j)=spval
          sfcuxi(i,j)=spval
          sfcvxi(i,j)=spval
          avgalbedo(i,j)=spval
          avgcprate(i,j)=spval
          avgprec(i,j)=spval
          avgprec_cont(i,j)=spval
          avgcprate_cont(i,j)=spval
          ptop(i,j)=spval
          pbot(i,j)=spval
          avgcfrach(i,j)=spval
          avgcfracm(i,j)=spval
          avgcfracl(i,j)=spval
          avgtcdc(i,j)=spval
          auvbin(i,j)=spval
          auvbinc(i,j)=spval
          ptopl(i,j)=spval
          pbotl(i,j)=spval
          Ttopl(i,j)=spval
          ptopm(i,j)=spval
          pbotm(i,j)=spval
          Ttopm(i,j)=spval
          ptoph(i,j)=spval
          pboth(i,j)=spval
          Ttoph(i,j)=spval
          sfcugs(i,j)=spval
          sfcvgs(i,j)=spval
          pblcfr(i,j)=spval
          cldwork(i,j)=spval
          gtaux(i,j)=spval
          gtauy(i,j)=spval
          cd10(i,j)=spval
          ch10(i,j)=spval
          mdltaux(i,j)=spval
          mdltauy(i,j)=spval
          runoff(i,j)=spval
        enddo
      enddo
      allocate(maxtshltr(im,jsta_2l:jend_2u))
      allocate(mintshltr(im,jsta_2l:jend_2u))
      allocate(maxrhshltr(im,jsta_2l:jend_2u))
      allocate(minrhshltr(im,jsta_2l:jend_2u))
      allocate(maxqshltr(im,jsta_2l:jend_2u))
      allocate(minqshltr(im,jsta_2l:jend_2u))
      allocate(dzice(im,jsta_2l:jend_2u))
      allocate(alwinc(im,jsta_2l:jend_2u))
      allocate(alwoutc(im,jsta_2l:jend_2u))
      allocate(alwtoac(im,jsta_2l:jend_2u))
      allocate(aswinc(im,jsta_2l:jend_2u))
      allocate(aswoutc(im,jsta_2l:jend_2u))
      allocate(aswtoac(im,jsta_2l:jend_2u))
      allocate(aswintoa(im,jsta_2l:jend_2u))
      allocate(smcwlt(im,jsta_2l:jend_2u))
      allocate(suntime(im,jsta_2l:jend_2u))
      allocate(fieldcapa(im,jsta_2l:jend_2u))
      allocate(avisbeamswin(im,jsta_2l:jend_2u))
      allocate(avisdiffswin(im,jsta_2l:jend_2u))
      allocate(airbeamswin(im,jsta_2l:jend_2u))
      allocate(airdiffswin(im,jsta_2l:jend_2u))
      allocate(snowfall(im,jsta_2l:jend_2u))
      allocate(acond(im,jsta_2l:jend_2u))
      allocate(edir(im,jsta_2l:jend_2u))
      allocate(ecan(im,jsta_2l:jend_2u))
      allocate(etrans(im,jsta_2l:jend_2u))
      allocate(esnow(im,jsta_2l:jend_2u))
      allocate(avgedir(im,jsta_2l:jend_2u))
      allocate(avgecan(im,jsta_2l:jend_2u))
      allocate(avgetrans(im,jsta_2l:jend_2u))
      allocate(avgesnow(im,jsta_2l:jend_2u))
      allocate(avgpotevp(im,jsta_2l:jend_2u))
      allocate(aod550(im,jsta_2l:jend_2u))
      allocate(du_aod550(im,jsta_2l:jend_2u))
      allocate(ss_aod550(im,jsta_2l:jend_2u))
      allocate(su_aod550(im,jsta_2l:jend_2u))
      allocate(oc_aod550(im,jsta_2l:jend_2u))
      allocate(bc_aod550(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          maxtshltr(i,j)=spval
          mintshltr(i,j)=spval
          maxrhshltr(i,j)=spval
          minrhshltr(i,j)=spval
          maxqshltr(i,j)=spval
          minqshltr(i,j)=spval
          dzice(i,j)=spval
          alwinc(i,j)=spval
          alwoutc(i,j)=spval
          alwtoac(i,j)=spval
          aswinc(i,j)=spval
          aswoutc(i,j)=spval
          aswtoac(i,j)=spval
          aswintoa(i,j)=spval
          smcwlt(i,j)=spval
          suntime(i,j)=spval
          fieldcapa(i,j)=spval
          avisbeamswin(i,j)=spval
          avisdiffswin(i,j)=spval
          airbeamswin(i,j)=spval
          airdiffswin(i,j)=spval
          snowfall(i,j)=spval
          acond(i,j)=spval
          edir(i,j)=spval
          ecan(i,j)=spval
          etrans(i,j)=spval
          esnow(i,j)=spval
          avgedir(i,j)=spval
          avgecan(i,j)=spval
          avgetrans(i,j)=spval
          avgesnow(i,j)=spval
          avgpotevp(i,j)=spval
          aod550(i,j)=spval
          du_aod550(i,j)=spval
          ss_aod550(i,j)=spval
          su_aod550(i,j)=spval
          oc_aod550(i,j)=spval
          bc_aod550(i,j)=spval
        enddo
      enddo
!
!     FROM MASKS
!
      allocate(hbm2(im,jsta_2l:jend_2u))
      allocate(sm(im,jsta_2l:jend_2u))
      allocate(sice(im,jsta_2l:jend_2u))
      allocate(lmh(im,jsta_2l:jend_2u))  ! real
      allocate(lmv(im,jsta_2l:jend_2u))  ! real
      allocate(gdlat(im,jsta_2l:jend_2u))
      allocate(gdlon(im,jsta_2l:jend_2u))
      allocate(dx(im,jsta_2l:jend_2u))
      allocate(dy(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,im
          hbm2(i,j)=spval
          sm(i,j)=spval
          sice(i,j)=spval
          lmh(i,j)=spval
          lmv(i,j)=spval
          gdlat(i,j)=spval
          gdlon(i,j)=spval
          dx(i,j)=spval
          dy(i,j)=spval
        enddo
      enddo

      if (me == 0) print *,' gocart_on=',gocart_on
      if (gocart_on) then
!  
! Add GOCART fields
! vrbls4d
        allocate(dust(im,jsta_2l:jend_2u,lm,nbin_du))
        allocate(salt(im,jsta_2l:jend_2u,lm,nbin_ss))
        allocate(soot(im,jsta_2l:jend_2u,lm,nbin_bc))
        allocate(waso(im,jsta_2l:jend_2u,lm,nbin_oc))
        allocate(suso(im,jsta_2l:jend_2u,lm,nbin_su))
        allocate(pp25(im,jsta_2l:jend_2u,lm,nbin_su))
        allocate(pp10(im,jsta_2l:jend_2u,lm,nbin_su))
!Initialization
!$omp parallel do private(i,j,l,k)
        do k=1,nbin_du
          do l=1,lm
            do j=jsta_2l,jend_2u
              do i=1,im
                dust(i,j,l,k)=spval
              enddo
            enddo
          enddo
        enddo
!$omp parallel do private(i,j,l,k)
        do k=1,nbin_ss
          do l=1,lm
            do j=jsta_2l,jend_2u
              do i=1,im
                salt(i,j,l,k)=spval
              enddo
            enddo
          enddo
        enddo
!$omp parallel do private(i,j,l,k)
        do k=1,nbin_bc
          do l=1,lm
            do j=jsta_2l,jend_2u
              do i=1,im
                soot(i,j,l,k)=spval
              enddo
            enddo
          enddo
        enddo
!$omp parallel do private(i,j,l,k)
        do k=1,nbin_oc
          do l=1,lm
            do j=jsta_2l,jend_2u
              do i=1,im
                waso(i,j,l,k)=spval
              enddo
            enddo
          enddo
        enddo
!$omp parallel do private(i,j,l,k)
        do k=1,nbin_su
          do l=1,lm
            do j=jsta_2l,jend_2u
              do i=1,im
                suso(i,j,l,k)=spval
                pp25(i,j,l,k)=spval
                pp10(i,j,l,k)=spval
              enddo
            enddo
          enddo
        enddo
! vrbls3d
        allocate(ext(im,jsta_2l:jend_2u,lm))
        allocate(asy(im,jsta_2l:jend_2u,lm))
        allocate(ssa(im,jsta_2l:jend_2u,lm))
        allocate(sca(im,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j)
        do l=1,lm
          do j=jsta_2l,jend_2u
            do i=1,im
              ext(i,j,l)=spval
              asy(i,j,l)=spval
              ssa(i,j,l)=spval
              sca(i,j,l)=spval
            enddo
          enddo
        enddo
        allocate(duem(im,jsta_2l:jend_2u,nbin_du))
        allocate(dusd(im,jsta_2l:jend_2u,nbin_du))
        allocate(dudp(im,jsta_2l:jend_2u,nbin_du))
        allocate(duwt(im,jsta_2l:jend_2u,nbin_du))
        allocate(dusv(im,jsta_2l:jend_2u,nbin_du))
        allocate(suem(im,jsta_2l:jend_2u,nbin_su))
        allocate(susd(im,jsta_2l:jend_2u,nbin_su))
        allocate(sudp(im,jsta_2l:jend_2u,nbin_su))
        allocate(suwt(im,jsta_2l:jend_2u,nbin_su))
        allocate(ocem(im,jsta_2l:jend_2u,nbin_oc))
        allocate(ocsd(im,jsta_2l:jend_2u,nbin_oc))
        allocate(ocdp(im,jsta_2l:jend_2u,nbin_oc))
        allocate(ocwt(im,jsta_2l:jend_2u,nbin_oc))
        allocate(ocsv(im,jsta_2l:jend_2u,nbin_oc))
        allocate(bcem(im,jsta_2l:jend_2u,nbin_bc))
        allocate(bcsd(im,jsta_2l:jend_2u,nbin_bc))
        allocate(bcdp(im,jsta_2l:jend_2u,nbin_bc))
        allocate(bcwt(im,jsta_2l:jend_2u,nbin_bc))
        allocate(bcsv(im,jsta_2l:jend_2u,nbin_bc))
        allocate(ssem(im,jsta_2l:jend_2u,nbin_ss))
        allocate(sssd(im,jsta_2l:jend_2u,nbin_ss))
        allocate(ssdp(im,jsta_2l:jend_2u,nbin_ss))
        allocate(sswt(im,jsta_2l:jend_2u,nbin_ss))
        allocate(sssv(im,jsta_2l:jend_2u,nbin_ss))
!Initialization
!$omp parallel do private(i,j,l)
        do l=1,nbin_du
          do j=jsta_2l,jend_2u
            do i=1,im
              duem(i,j,l)=spval
              dusd(i,j,l)=spval
              dudp(i,j,l)=spval
              duwt(i,j,l)=spval
              dusv(i,j,l)=spval
            enddo
          enddo
        enddo

        do l=1,nbin_su
          do j=jsta_2l,jend_2u
            do i=1,im
              suem(i,j,l)=spval
              susd(i,j,l)=spval
              sudp(i,j,l)=spval
              suwt(i,j,l)=spval
            enddo
          enddo
        enddo

        do l=1,nbin_oc
          do j=jsta_2l,jend_2u
            do i=1,im
              ocem(i,j,l)=spval
              ocsd(i,j,l)=spval
              ocdp(i,j,l)=spval
              ocwt(i,j,l)=spval
              ocsv(i,j,l)=spval
            enddo
          enddo
        enddo

        do l=1,nbin_bc
          do j=jsta_2l,jend_2u
            do i=1,im
              bcem(i,j,l)=spval
              bcsd(i,j,l)=spval
              bcdp(i,j,l)=spval
              bcwt(i,j,l)=spval
              bcsv(i,j,l)=spval
            enddo
          enddo
        enddo

        do l=1,nbin_ss
          do j=jsta_2l,jend_2u
            do i=1,im
              ssem(i,j,l)=spval
              sssd(i,j,l)=spval
              ssdp(i,j,l)=spval
              sswt(i,j,l)=spval
              sssv(i,j,l)=spval
            enddo
          enddo
        enddo
        allocate(rhomid(im,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
        do l=1,lm
          do j=jsta_2l,jend_2u
            do i=1,im
              rhomid(i,j,l)=spval
            enddo
          enddo
        enddo
! vrbls2d
        allocate(dusmass(im,jsta_2l:jend_2u))
        allocate(ducmass(im,jsta_2l:jend_2u))
        allocate(dusmass25(im,jsta_2l:jend_2u))
        allocate(ducmass25(im,jsta_2l:jend_2u))
        allocate(susmass(im,jsta_2l:jend_2u))
        allocate(sucmass(im,jsta_2l:jend_2u))
        allocate(susmass25(im,jsta_2l:jend_2u))
        allocate(sucmass25(im,jsta_2l:jend_2u))
        allocate(ocsmass(im,jsta_2l:jend_2u))
        allocate(occmass(im,jsta_2l:jend_2u))
        allocate(ocsmass25(im,jsta_2l:jend_2u))
        allocate(occmass25(im,jsta_2l:jend_2u))
        allocate(bcsmass(im,jsta_2l:jend_2u))
        allocate(bccmass(im,jsta_2l:jend_2u))
        allocate(bcsmass25(im,jsta_2l:jend_2u))
        allocate(bccmass25(im,jsta_2l:jend_2u))
        allocate(sssmass(im,jsta_2l:jend_2u))
        allocate(sscmass(im,jsta_2l:jend_2u))
        allocate(sssmass25(im,jsta_2l:jend_2u))
        allocate(sscmass25(im,jsta_2l:jend_2u))
        allocate(dustcb(im,jsta_2l:jend_2u))
        allocate(occb(im,jsta_2l:jend_2u))
        allocate(bccb(im,jsta_2l:jend_2u))
        allocate(sulfcb(im,jsta_2l:jend_2u))
        allocate(pp25cb(im,jsta_2l:jend_2u))
        allocate(pp10cb(im,jsta_2l:jend_2u))
        allocate(sscb(im,jsta_2l:jend_2u))
        allocate(dustallcb(im,jsta_2l:jend_2u))
        allocate(ssallcb(im,jsta_2l:jend_2u))
        allocate(dustpm(im,jsta_2l:jend_2u))
        allocate(sspm(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
       do j=jsta_2l,jend_2u
         do i=1,lm
           dusmass(i,j)=spval
           ducmass(i,j)=spval
           dusmass25(i,j)=spval
           ducmass25(i,j)=spval
           susmass(i,j)=spval
           sucmass(i,j)=spval
           susmass25(i,j)=spval
           sucmass25(i,j)=spval
           ocsmass(i,j)=spval
           occmass(i,j)=spval
           ocsmass25(i,j)=spval
           occmass25(i,j)=spval
           bcsmass(i,j)=spval
           bccmass(i,j)=spval
           bcsmass25(i,j)=spval
           bccmass25(i,j)=spval
           sssmass(i,j)=spval
           sscmass(i,j)=spval
           sssmass25(i,j)=spval
           sscmass25(i,j)=spval
           dustcb(i,j)=spval
           occb(i,j)=spval
           bccb(i,j)=spval
           sulfcb(i,j)=spval
           pp25cb(i,j)=spval
           pp10cb(i,j)=spval
           sscb(i,j)=spval
           dustallcb(i,j)=spval
           ssallcb(i,j)=spval
           dustpm(i,j)=spval
           sspm(i,j)=spval
         enddo
       enddo
      endif
! HWRF RRTMG output 
      allocate(acswupt(im,jsta_2l:jend_2u))
      allocate(swdnt(im,jsta_2l:jend_2u))
      allocate(acswdnt(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,lm
          acswupt(i,j)=spval
          swdnt(i,j)=spval
          acswdnt(i,j)=spval
        enddo
      enddo

! UPP_MATH MODULE DIFFERENTIAL EQUATIONS
      allocate(ddvdx(im,jsta_2l:jend_2u))
      allocate(ddudy(im,jsta_2l:jend_2u))
      allocate(uuavg(im,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=1,lm
          ddvdx(i,j)=spval
          ddudy(i,j)=spval
          uuavg(i,j)=spval
        enddo
      enddo
! 
      end
