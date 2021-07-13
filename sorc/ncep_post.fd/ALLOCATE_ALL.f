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
      allocate(u(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(v(ista_2l:iend_2u,jsta_2l:jvend_2u,lm))
      allocate(t(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
! CHUANG ADD POTENTIAL TEMP BECAUSE WRF OUTPUT THETA
!      allocate(th(ista_2l:iend_2u,jsta_2l:jend_2u,lm))   
      allocate(q(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
!      allocate(w(ista_2l:iend_2u,jsta_2l:jend_2u,lp1))
      allocate(uh(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(vh(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(wh(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(pmid(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(pmidv(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(pint(ista_2l:iend_2u,jsta_2l:jend_2u,lp1))
      allocate(alpint(ista_2l:iend_2u,jsta_2l:jend_2u,lp1))
      allocate(zmid(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(zint(ista_2l:iend_2u,jsta_2l:jend_2u,lp1))
!      allocate(rainw(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(q2(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(omga(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(dpres(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(T_ADJ(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(ttnd(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(rswtt(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(rlwtt(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(exch_h(ista_2l:iend_2u,jsta_2l:jend_2u,lm)) 
      allocate(train(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(tcucn(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(EL_PBL(ista_2l:iend_2u,jsta_2l:jend_2u,lm))

!Initialization
!$omp parallel do private(i,j,l)
      do l=1,lm
        do j=jsta_2l,jend_2u
          do i=ista_2l,iend_2u
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
          do i=ista_2l,iend_2u
            pint(i,j,l)=spval
            alpint(i,j,l)=spval
            zint(i,j,l)=spval
          enddo
        enddo
      enddo

!     MP FIELD   
      allocate(cwm(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(F_ice(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(F_rain(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(F_RimeF(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QQW(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QRIMEF(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QQI(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QQR(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QQS(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QQG(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QQNW(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QQNI(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QQNR(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QQNWFA(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QQNIFA(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(TAOD5503D(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(AEXTC55(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(EXTCOF55(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(QC_BL(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(CFR(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(CFR_RAW(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(DBZ(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(DBZR(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(DBZI(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(DBZC(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(mcvg(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(NLICE(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
      do l=1,lm
      do j=jsta_2l,jend_2u
          do i=ista_2l,iend_2u 
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
      allocate(NRAIN(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(radius_cloud(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(radius_ice(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(radius_snow(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
! KRS: HWRF Addition for thompson reflectivity
! or non-ferrier physics. wrf-derived
      allocate(REFL_10CM(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
!GFS FIELD
      allocate(o3(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(o(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(o2(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(tcucns(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
      do l=1,lm
      do j=jsta_2l,jend_2u
          do i=ista_2l,iend_2u
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
        allocate(vdifftt(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
!       allocate(tcucns(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(vdiffmois(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(dconvmois(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(sconvmois(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(nradtt(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(o3vdiff(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(o3prod(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(o3tndy(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(mwpv(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(unknown(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(vdiffzacce(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(zgdrag(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(cnvctummixing(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(vdiffmacce(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(mgdrag(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(cnvctvmmixing(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(ncnvctcfrac(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(cnvctumflx(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(cnvctdmflx(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(cnvctdetmflx(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(cnvctzgdrag(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(cnvctmgdrag(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
        do l=1,lm
          do j=jsta_2l,jend_2u
            do i=ista_2l,iend_2u
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
      allocate(htm(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(vtm(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
! add GFIP ICING
      allocate(icing_gfip(ista_2l:iend_2u,jsta_2l:jend_2u,lm))        
      allocate(icing_gfis(ista_2l:iend_2u,jsta_2l:jend_2u,lm))        
!
! add GTG turbulence
      allocate(catedr(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(mwt(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
      allocate(gtg(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
      do l=1,lm
        do j=jsta_2l,jend_2u
          do i=ista_2l,iend_2u
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
      allocate(smc(ista_2l:iend_2u,jsta_2l:jend_2u,nsoil))
      allocate(stc(ista_2l:iend_2u,jsta_2l:jend_2u,nsoil))
      allocate(sh2o(ista_2l:iend_2u,jsta_2l:jend_2u,nsoil))
      allocate(SLDPTH(NSOIL))
      allocate(RTDPTH(NSOIL))
      allocate(SLLEVEL(NSOIL))
!Initialization
!$omp parallel do private(i,j,l)
      do l=1,nsoil
        do j=jsta_2l,jend_2u
          do i=ista_2l,iend_2u
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
      allocate(wspd10max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(w_up_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(w_dn_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(w_mean(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(refd_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(prate_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(fprate_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(up_heli_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(up_heli_max16(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(up_heli_min(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(up_heli_min16(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(up_heli_max02(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(up_heli_min02(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(up_heli_max03(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(up_heli_min03(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
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
      allocate(rel_vort_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(rel_vort_max01(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(rel_vort_maxhy1(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(wspd10umax(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(wspd10vmax(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(refdm10c_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(hail_max2d(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(hail_maxk1(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(hail_maxhailcast(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(grpl_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(up_heli(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(up_heli16(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ltg1_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ltg2_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ltg3_max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(nci_ltg(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(nca_ltg(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(nci_wq(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(nca_wq(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(nci_refd(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(nca_refd(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
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
      allocate(REF_10CM(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
      do l=1,lm
        do j=jsta_2l,jend_2u
          do i=ista_2l,iend_2u
            REF_10CM(i,j,l)=spval
          enddo
        enddo
      enddo
      allocate(REFC_10CM(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(REF1KM_10CM(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(REF4KM_10CM(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
          REFC_10CM(i,j)=spval
          REF1KM_10CM(i,j)=spval
          REF4KM_10CM(i,j)=spval
        enddo
      enddo
! CRA
      allocate(u10(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(v10(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(tshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(qshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(mrshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(smstav(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ssroff(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(bgroff(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(vegfrc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(shdmin(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(shdmax(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(lai(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(acsnow(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(acgraup(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(acfrain(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(acsnom(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cmc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sst(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(qz0(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(thz0(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(uz0(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(vz0(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(qs(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ths(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sno(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(snonc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ti(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
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
      allocate(u10mean(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(v10mean(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(spduv10mean(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(swradmean(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(swnormmean(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
          u10mean(i,j)=spval
          v10mean(i,j)=spval
          spduv10mean(i,j)=spval
          swradmean(i,j)=spval
          swnormmean(i,j)=spval
        enddo
      enddo
!NAMstart
      allocate(snoavg(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(psfcavg(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(t10m(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(t10avg(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(akmsavg(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(akhsavg(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(u10max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(v10max(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(u10h(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(v10h(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
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
      allocate(akms(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(akhs(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cuprec(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(acprec(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ancprc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cuppt(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
          akms(i,j)=spval
          akhs(i,j)=spval
          cuprec(i,j)=spval
          acprec(i,j)=spval
          ancprc(i,j)=spval
          cuppt(i,j)=spval
        enddo
      enddo
! GSDstart
      allocate(rainc_bucket(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(rainc_bucket1(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(rainnc_bucket(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(rainnc_bucket1(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pcp_bucket(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pcp_bucket1(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(snow_bucket(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(snow_bucket1(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(graup_bucket(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(graup_bucket1(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(qrmax(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(tmax(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(snownc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(graupelnc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(tsnow(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(qvg(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(qv2m(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(qvl1(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(snfden(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sndepac(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(int_smoke(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(mean_frp(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(int_aod(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
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
      allocate(smoke(ista_2l:iend_2u,jsta_2l:jend_2u,lm,nbin_sm))
!$omp parallel do private(i,j,l,k)
      do k=1,nbin_sm
        do l=1,lm
          do j=jsta_2l,jend_2u
            do i=ista_2l,iend_2u
              smoke(i,j,l,k)=spval
            enddo
          enddo
        enddo
      enddo
! GSDend
      allocate(rswin(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(swddni(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(swddif(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(swdnbc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(swddnic(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(swddifc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(swupbc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(swupt(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(taod5502d(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(aerasy2d(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(aerssa2d(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(lwp(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(iwp(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(rlwin(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(lwdnbc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(lwupbc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(rlwtoa(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(rswtoa(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
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
      allocate(tg(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sfcshx(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sfclhx(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(fis(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(t500(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(t700(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(z500(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(z700(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(teql(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cfracl(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cfracm(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cfrach(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(acfrst(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(acfrcv(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(hbot(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(htop(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(aswin(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(alwin(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(aswout(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(alwout(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(aswtoa(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(alwtoa(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
          tg(i,j)=spval
          sfcshx(i,j)=spval
          sfclhx(i,j)=spval
          fis(i,j)=spval
          t500(i,j)=spval
          t700(i,j)=spval
          z700(i,j)=spval
          teql(i,j)=spval
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
      allocate(czen(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(czmean(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sigt4(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(rswout(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(radot(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ncfrst(ista_2l:iend_2u,jsta_2l:jend_2u))  ! real
      allocate(ncfrcv(ista_2l:iend_2u,jsta_2l:jend_2u))  ! real
      allocate(smstot(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pctsno(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(th10(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(q10(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(prec(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(subshx(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(snopcx(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sfcuvx(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sfcevp(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(potevp(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(z0(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ustar(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pblh(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pblhgust(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(mixht(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(twbs(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(qwbs(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
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
      allocate(sfcexc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(grnflx(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(soiltb(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(z1000(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(slp(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pslp(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(f(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(albedo(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(albase(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cldfra(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cprate(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cnvcfr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ivgtyp(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(isltyp(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(hbotd(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(htopd(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(hbots(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(htops(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cldefi(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(islope(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(si(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(lspa(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(rswinc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(vis(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pd(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(mxsnal(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(epsr(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
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
      allocate(sfcux(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sfcvx(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sfcuxi(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sfcvxi(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgalbedo(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgcprate(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgprec(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgprec_cont(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgcprate_cont(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ptop(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pbot(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgcfrach(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgcfracm(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgcfracl(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgtcdc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(auvbin(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(auvbinc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ptopl(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pbotl(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(Ttopl(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ptopm(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pbotm(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(Ttopm(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ptoph(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pboth(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(Ttoph(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sfcugs(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sfcvgs(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(pblcfr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cldwork(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(gtaux(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(gtauy(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(cd10(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ch10(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(mdltaux(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(mdltauy(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(runoff(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
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
      allocate(maxtshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(mintshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(maxrhshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(minrhshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(maxqshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(minqshltr(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(dzice(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(alwinc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(alwoutc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(alwtoac(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(aswinc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(aswoutc(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(aswtoac(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(aswintoa(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(smcwlt(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(suntime(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(fieldcapa(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avisbeamswin(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avisdiffswin(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(airbeamswin(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(airdiffswin(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(snowfall(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(acond(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(edir(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ecan(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(etrans(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(esnow(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgedir(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgecan(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgetrans(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgesnow(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(avgpotevp(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(aod550(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(du_aod550(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ss_aod550(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(su_aod550(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(oc_aod550(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(bc_aod550(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
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
      allocate(hbm2(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sm(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(sice(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(lmh(ista_2l:iend_2u,jsta_2l:jend_2u))  ! real
      allocate(lmv(ista_2l:iend_2u,jsta_2l:jend_2u))  ! real
      allocate(gdlat(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(gdlon(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(dx(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(dy(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
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
        allocate(dust(ista_2l:iend_2u,jsta_2l:jend_2u,lm,nbin_du))
        allocate(salt(ista_2l:iend_2u,jsta_2l:jend_2u,lm,nbin_ss))
        allocate(soot(ista_2l:iend_2u,jsta_2l:jend_2u,lm,nbin_bc))
        allocate(waso(ista_2l:iend_2u,jsta_2l:jend_2u,lm,nbin_oc))
        allocate(suso(ista_2l:iend_2u,jsta_2l:jend_2u,lm,nbin_su))
        allocate(pp25(ista_2l:iend_2u,jsta_2l:jend_2u,lm,nbin_su))
        allocate(pp10(ista_2l:iend_2u,jsta_2l:jend_2u,lm,nbin_su))
!Initialization
!$omp parallel do private(i,j,l,k)
        do k=1,nbin_du
          do l=1,lm
            do j=jsta_2l,jend_2u
              do i=ista_2l,iend_2u
                dust(i,j,l,k)=spval
              enddo
            enddo
          enddo
        enddo
!$omp parallel do private(i,j,l,k)
        do k=1,nbin_ss
          do l=1,lm
            do j=jsta_2l,jend_2u
              do i=ista_2l,iend_2u
                salt(i,j,l,k)=spval
              enddo
            enddo
          enddo
        enddo
!$omp parallel do private(i,j,l,k)
        do k=1,nbin_bc
          do l=1,lm
            do j=jsta_2l,jend_2u
              do i=ista_2l,iend_2u
                soot(i,j,l,k)=spval
              enddo
            enddo
          enddo
        enddo
!$omp parallel do private(i,j,l,k)
        do k=1,nbin_oc
          do l=1,lm
            do j=jsta_2l,jend_2u
              do i=ista_2l,iend_2u
                waso(i,j,l,k)=spval
              enddo
            enddo
          enddo
        enddo
!$omp parallel do private(i,j,l,k)
        do k=1,nbin_su
          do l=1,lm
            do j=jsta_2l,jend_2u
              do i=ista_2l,iend_2u
                suso(i,j,l,k)=spval
                pp25(i,j,l,k)=spval
                pp10(i,j,l,k)=spval
              enddo
            enddo
          enddo
        enddo
! vrbls3d
        allocate(ext(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(asy(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(ssa(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
        allocate(sca(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j)
        do l=1,lm
          do j=jsta_2l,jend_2u
            do i=ista_2l,iend_2u
              ext(i,j,l)=spval
              asy(i,j,l)=spval
              ssa(i,j,l)=spval
              sca(i,j,l)=spval
            enddo
          enddo
        enddo
        allocate(duem(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_du))
        allocate(dusd(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_du))
        allocate(dudp(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_du))
        allocate(duwt(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_du))
        allocate(dusv(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_du))
        allocate(suem(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_su))
        allocate(susd(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_su))
        allocate(sudp(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_su))
        allocate(suwt(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_su))
        allocate(ocem(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_oc))
        allocate(ocsd(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_oc))
        allocate(ocdp(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_oc))
        allocate(ocwt(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_oc))
        allocate(ocsv(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_oc))
        allocate(bcem(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_bc))
        allocate(bcsd(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_bc))
        allocate(bcdp(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_bc))
        allocate(bcwt(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_bc))
        allocate(bcsv(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_bc))
        allocate(ssem(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_ss))
        allocate(sssd(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_ss))
        allocate(ssdp(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_ss))
        allocate(sswt(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_ss))
        allocate(sssv(ista_2l:iend_2u,jsta_2l:jend_2u,nbin_ss))
!Initialization
!$omp parallel do private(i,j,l)
        do l=1,nbin_du
          do j=jsta_2l,jend_2u
            do i=ista_2l,iend_2u
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
            do i=ista_2l,iend_2u
              suem(i,j,l)=spval
              susd(i,j,l)=spval
              sudp(i,j,l)=spval
              suwt(i,j,l)=spval
            enddo
          enddo
        enddo

        do l=1,nbin_oc
          do j=jsta_2l,jend_2u
            do i=ista_2l,iend_2u
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
            do i=ista_2l,iend_2u
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
            do i=ista_2l,iend_2u
              ssem(i,j,l)=spval
              sssd(i,j,l)=spval
              ssdp(i,j,l)=spval
              sswt(i,j,l)=spval
              sssv(i,j,l)=spval
            enddo
          enddo
        enddo
        allocate(rhomid(ista_2l:iend_2u,jsta_2l:jend_2u,lm))
!Initialization
!$omp parallel do private(i,j,l)
        do l=1,lm
          do j=jsta_2l,jend_2u
            do i=ista_2l,iend_2u
              rhomid(i,j,l)=spval
            enddo
          enddo
        enddo
! vrbls2d
        allocate(dusmass(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(ducmass(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(dusmass25(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(ducmass25(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(susmass(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(sucmass(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(susmass25(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(sucmass25(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(ocsmass(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(occmass(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(ocsmass25(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(occmass25(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(bcsmass(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(bccmass(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(bcsmass25(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(bccmass25(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(sssmass(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(sscmass(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(sssmass25(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(sscmass25(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(dustcb(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(occb(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(bccb(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(sulfcb(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(pp25cb(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(pp10cb(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(sscb(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(dustallcb(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(ssallcb(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(dustpm(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(sspm(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
       do j=jsta_2l,jend_2u
         do i=ista_2l,iend_2u
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
      allocate(acswupt(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(swdnt(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(acswdnt(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
          acswupt(i,j)=spval
          swdnt(i,j)=spval
          acswdnt(i,j)=spval
        enddo
      enddo

! UPP_MATH MODULE DIFFERENTIAL EQUATIONS
      allocate(ddvdx(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ddudy(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(uuavg(ista_2l:iend_2u,jsta_2l:jend_2u))
!Initialization
!$omp parallel do private(i,j)
      do j=jsta_2l,jend_2u
        do i=ista_2l,iend_2u
          ddvdx(i,j)=spval
          ddudy(i,j)=spval
          uuavg(i,j)=spval
        enddo
      enddo
! 
      end
