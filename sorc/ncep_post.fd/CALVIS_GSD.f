!**********************************************************************c
      SUBROUTINE CALVIS_GSD(CZEN,VIS)

! SUBPROGRAM:    CALVIS      CALCULATE HORIZONTAL VISIBILITY   
!
!   PRGMMR:  BENJAMIN, STAN ORG: NOAA/FSL       DATE: 99-09-07
!
! ABSTRACT:
!
!   Started with Stoelinga-Warner algorithm for hydrometeors only.
!    Added coefficients for graupel.
!    Added algorithm for clear-air RH-based visibility.
!
!   This routine computes horizontal visibility (in km) at the
!   surface or lowest model layer, from qc, qr, qi, qs, and qg.  
!   qv--water vapor mixing ratio (kg/kg)
!   qc--cloud water mixing ratio (kg/kg)
!   qr--rain water mixing ratio  (kg/kg)
!   qi--cloud ice mixing ratio   (kg/kg)
!   qs--snow mixing ratio        (kg/kg)
!   qg--graupel mixing ratio     (kg/kg)
!   u/v - u/v wind components    (m/s)  
!   tb --            temperature (k)
!   pp--pressure                 (Pa)
!   rhb-- relative humidity      (0-100%)
!   aextc55--aerosol extinction coefficient (m**-1)
!
!
!   Independent of the above definitions, the scheme can use different
!   assumptions of the state of hydrometeors:
!        meth='r': Uses the four mixing ratios qrain, qsnow, qclw,
!           and qclice
!
!   The routine uses the following
!   expressions for extinction coefficient, beta (in km**-1),
!   with C being the mass concentration (in g/m**3):
!
!      cloud water:  beta = 144.7 * C ** (0.8800)
!      rain water:   beta =  2.24 * C ** (0.7500)
!      cloud ice:    beta = 327.8 * C ** (1.0000)
!      snow:         beta = 10.36 * C ** (0.7776)
!      graupel:      beta =  8.0  * C ** (0.7500)
!
!   These expressions were obtained from the following sources:
!
!      for cloud water: from Kunkel (1984)
!      for rainwater: from M-P dist'n, with No=8e6 m**-4 and
!         rho_w=1000 kg/m**3
!      for cloud ice: assume randomly oriented plates which follow
!         mass-diameter relationship from Rutledge and Hobbs (1983)
!      for snow: from Stallabrass (1985), assuming beta = -ln(.02)/vis
!      for graupel: guestimate by John Brown and Stan Benjamin,
!         similar to snow, but a smaller extinction coef seemed
!         reasonable.  27 Aug 99
!
!   The extinction coefficient for each water species present is
!   calculated, and then all applicable betas are summed to yield
!   a single beta. Then the following relationship is used to
!   determine visibility (in km), where epsilon is the threshhold
!   of contrast, usually taken to be .02:
!
!      vis = -ln(epsilon)/beta      [found in Kunkel (1984)]
!
!   The 'aextc55' field is 3-D and is derived from the 'aod_3d' field
!   by dividing by dz, the vertical thickness of that model level in 'm'.
!   This can be handled as a 2-D field if needed to save resources.
!
! HISTORY
! PROGRAM HISTORY LOG:
!    99-05-                        Version from Eta model and from 
!                                    Mark Stoelinga and Tom Warner
!    99-09-07      S. Benjamin     Modified for MM5 microphysics variables
!                                    include graupel mixing ratio
!    99-09         S. Benjamin     Added algorithm for RH-based clear-air
!                                    visibility
!    00-08         S. Benjamin     Added mods for base of 60km instead of 90km,
!                                    max RH from lowest 2 levels instead of
!                                    lev 2 only, max hydrometeor mix ratio
!                                    from lowest 5 levs instead of lev 1 only
!                                  Based on Schwartz stats and Smirnova et al
!                                    paper, and on METAR verif started this week
!    Dec 03        S. Benjamin     - updates
!                              - day/night distinction for vis constants
!                                  from Roy Rasmussen
!                              - low-level wind shear term 
!                                  - recommended by Evan Kuchera
!   2015-17        S. Benjamin, T. Smirnova - modifications for RH-based clear-air vis
!   2017-12        R. Ahmadov, Steve Albers - addition for attenuation from aerosols
!                              (not related to water vapor or RH at this point)
!   2021-05        Wen Meng        Unify CONST1 and VISRH. 
!   2021-05        Wen Meng  - Add checking for undefined points invloved in computation
!   2021-08        Wen Meng  - Restrict divided by 0.
!   2021-10        Jesse Meng - 2D DECOMPOSITION
!   2023-11        Tim Corrie, Eric James - addition of attenuation for blowing snow
!                           
!------------------------------------------------------------------
!

      use vrbls2d, only: sno, si, ustar
      use vrbls3d, only: qqw, qqi, qqs, qqr, qqg, t, pmid, q, u, v, extcof55, aextc55
      use params_mod, only: h1, d608, rd, g
      use ctlblk_mod, only: jm, im, jsta_2l, jend_2u, lm, modelname, spval, method_blsn,&
                                    ista_2l, iend_2u

      implicit none

      integer :: j, i, k, ll
      integer :: method
      real :: tx, pol, esx, es, e
      REAL VIS(ista_2l:iend_2u,jsta_2l:jend_2u)
      REAL RHB(ista_2l:iend_2u,jsta_2l:jend_2u,LM)
      REAL CZEN(ista_2l:iend_2u,jsta_2l:jend_2u)

      real :: z, ustar_t, u_p, lamda, r_bar, alpha
      real :: rho_sno
      real :: z_r, Q_s, C_r, c_z, c_alpha, vis_blsn, BETABLSN

      real celkel,tice,coeflc,coeflp,coeffc,coeffp,coeffg
      real exponlc,exponlp,exponfc,exponfp,exponfg,const1
      real rhoice,rhowat,qrain,qsnow,qgraupel,qclw,qclice,tv,rhoair,  &
        vovermd,conclc,conclp,concfc,concfp,concfg,betav

      real coeffp_dry, coeffp_wet, shear_fac, temp_fac
      real coef_snow, shear

      real coefrh,qrh,visrh
      real rhmax,shear5_cnt, shear8_cnt
      real shear5_cnt_lowvis, shear8_cnt_lowvis
      real shear4_cnt, shear4_cnt_lowvis
      integer night_cnt, lowsun_cnt

      real visrh10_cnt, vis1km_cnt, visrh_lower
      real vis3km_cnt
      real vis5km_cnt
      real vis_min, visrh_min
      real vis_night, zen_fac
!------------------------------------------------------------------

! Method used for clear-air visibility with extension for aerosols
      method = 3 
!                RH-only method (1), 
!                Aerosol method (2), 
!                Smoke added to RH method for clear air  (3)
!                   3 - option to add reducted visibility from smoke-based aerosols.
     
      CELKEL     = 273.15
      TICE       = CELKEL-10.
      COEFLC     = 144.7
      COEFLP     =   2.24
      COEFFC     = 327.8
      COEFFP     =  10.36

! - Initialize counters

      shear4_cnt = 0 
      shear5_cnt = 0 
      shear8_cnt = 0
      shear4_cnt_lowvis = 0
      shear5_cnt_lowvis = 0 
      shear8_cnt_lowvis = 0
      night_cnt = 0 
      lowsun_cnt = 0
      visrh10_cnt = 0 
      visrh_lower = 0
      vis1km_cnt = 0 
      vis3km_cnt = 0
      vis5km_cnt = 0

! - snow-based vis attenuation - coefficient values from Roy Rasmussen - Dec 2003
!      COEFFP_dry =  17.7
!      COEFFP_wet =   4.18
! - modified number - Stan B. - Dec 2007
!     after quick talks with John Brown and Ismail Gultepe
      COEFFP_dry =  10.0
      COEFFP_wet =   6.0

!     COEFFg     =   8.0
! - values from Roy Rasmussen - Dec 2003
!    Rasmussen et al. 2003,  J. App. Meteor.
!    Snow Nowcasting Using a Real-Time Correlation of Radar Reflectivity with Snow Gauge Accumulation
      COEFFg     =   4.0

      EXPONLC    =   0.8800
      EXPONLP    =   0.7500
      EXPONFC    =   1.0000
!     EXPONFP    =   0.7776
! - new value from Roy Rasmussen - Dec 2003  
      EXPONFP    =   1.0

      EXPONFg    =   0.75  
!     CONST1=-LOG(.02)
!     CONST1=3.912
      CONST1= 3.000

! visibility with respect to RH is
!   calculated from optical depth linearly
!   related to RH as follows:

!    vis = 60 exp (-2.5 * (RH-15)/80)
!       changed on 3/14/01

!  coefficient of 3 gives visibility of 5 km
!    at 95% RH

! Total visibility is minimum of vis-rh  (developed by Benjamin, Brown, Smirnova)
!   and vis-hydrometeors from Stoelinga/Warner

      RHOICE=917.
      RHOWAT=1000.

      vis_min = 1.e6
      visrh_min = 1.e6
 
      DO J=jsta_2l,jend_2u
      DO I=ista_2l,iend_2u
        VIS(I,J)=spval
! -checking undedined points
        if(T(I,J,LM)<spval .and. U(I,J,LM)<spval .and. V(I,J,LM)<spval &
           .and. PMID(I,J,LM)<spval) then 
!  - take max hydrometeor mixing ratios in lowest 3 levels (lowest 13 hPa, 100m with RAP/HRRR
      qrain = 0.
      qsnow = 0.
      qgraupel = 0.
      qclw = 0.
      qclice = 0.

          do k = 1,3
               LL=LM-k+1
            if(QQW(I,J,ll)<spval)qclw = max(qclw, QQW(I,J,ll) )
            if(QQI(I,J,ll)<spval)qclice = max(qclice, QQI(I,J,ll) )
            if(QQS(I,J,ll)<spval)qsnow = max(qsnow, QQS(I,J,ll) )
            if(QQR(I,J,ll)<spval)qrain = max(qrain, QQR(I,J,ll) )
            if(QQG(I,J,ll)<spval)qgraupel  = max(qgraupel, QQG(I,J,ll) )
! - compute relative humidity
        Tx=T(I,J,LL)-273.15
        RHB(I,J,LL)=0.
        POL = 0.99999683       + TX*(-0.90826951E-02 +    &
           TX*(0.78736169E-04   + TX*(-0.61117958E-06 +   &
           TX*(0.43884187E-08   + TX*(-0.29883885E-10 +   &
           TX*(0.21874425E-12   + TX*(-0.17892321E-14 +   &
           TX*(0.11112018E-16   + TX*(-0.30994571E-19)))))))))
        if(abs(POL) > 0.) THEN
        esx = 6.1078/POL**8

          ES = esx
          E = PMID(I,J,LL)/100.*Q(I,J,LL)/(0.62197+Q(I,J,LL)*0.37803)
          RHB(I,J,LL) = 100.*AMIN1(1.,E/ES)
       ENDIF

          enddo

!  - take max RH of levels 1 and 2 near the sfc
          rhmax = max (rhb(i,j,lm),rhb(i,j,lm-1))
          qrh = max(0.0,min(0.8,(rhmax/100.-0.15)))

!tgs 23 feb 2017 - increase of base value to 90 km to reduce attenuation
!                  from RH for clear-air visibility.  (i.e., increase clear-air vis overall)
       visrh = 90. * exp(-2.5*qrh)

!  -- add term to increase RH vis term for
!     low-level wind shear increasing from 4 to 6 ms-1
!     (using Evan Kuchera's paper as a guideline)

! -- calculate term for shear between levels 1 and 4
!   (about 25 hPa for HRRR/RAP)
         shear = sqrt( (u(i,j,lm-3)-u(i,j,lm))**2         &
                     +(v(i,j,lm-3)-v(i,j,lm))**2  )

        shear_fac = min(1.,max(0.,(shear-4.)/2.) )
        if (visrh<10.) visrh = visrh + (10.-visrh)*    &
           shear_fac

        if (shear>4.) shear4_cnt = shear4_cnt +1
        if (shear>5.) shear5_cnt = shear5_cnt +1
        if (shear>6.) shear8_cnt = shear8_cnt +1

        if (shear>4..and.visrh<10)                  &
          shear4_cnt_lowvis = shear4_cnt_lowvis +1
        if (shear>5..and.visrh<10)                  &
          shear5_cnt_lowvis = shear5_cnt_lowvis +1
        if (shear>6..and.visrh<10)                  &
          shear8_cnt_lowvis = shear8_cnt_lowvis +1

        if (visrh<10.) visrh10_cnt = visrh10_cnt+1
        if (czen(i,j)<0.) night_cnt = night_cnt + 1
        if (czen(i,j)<0.1) lowsun_cnt = lowsun_cnt + 1

        TV=T(I,J,lm)*(H1+D608*Q(I,J,lm))

        RHOAIR=PMID(I,J,lm)/(RD*TV)

          VOVERMD=(1.+Q(I,J,lm))/RHOAIR+(QCLW+QRAIN)/RHOWAT+    &
!          VOVERMD=(1.+Q(I,J,1))/RHOAIR+(QCLW+QRAIN)/RHOWAT+
                  (qgraupel+QCLICE+QSNOW)/RHOICE
          CONCLC=QCLW/VOVERMD*1000.
          CONCLP=QRAIN/VOVERMD*1000.
          CONCFC=QCLICE/VOVERMD*1000.
          CONCFP=QSNOW/VOVERMD*1000.
          CONCFg=Qgraupel/VOVERMD*1000.

          temp_fac = min(1.,max((t(i,j,lm)-271.15),0.) )

         coef_snow = coeffp_dry*(1.-temp_fac)                   &
                   + coeffp_wet* temp_fac    

          if (t(i,j,lm)< 270. .and. temp_fac==1.)          &
             write (6,*) 'Problem w/ temp_fac - calvis'

! Key calculation of attenuation from blowing snow -- updated 5 August 2022 by Tim Corrie
! Framework is from Letcher et al (2021)

        ustar_t = 0.2
        u_p = 2.8*ustar_t
        lamda = 0.45
        z = 2.0
        alpha = 15.0
        r_bar = 0.0002

!        print *, i,j


        if (si(i,j)<spval .and. si(i,j) .ge. 1.0) then
            z_r = 1.6*(ustar(i,j)**2./(2.*g))
            Q_s = max((0.68/ustar(i,j))*(RHOAIR/g)*(ustar(i,j)**2.-ustar_t**2.),0.0)
            C_r = (Q_s/u_p)*(lamda*g/ustar(i,j)**2.)*exp(-lamda*z_r*g/ustar(i,j)**2.)
            c_z = max(C_r * exp(-1.55*((0.05628*ustar(i,j))**-0.544 - z**-0.544)),1e-15)
            c_alpha = alpha/(alpha+2) !simplified version of (6) in Letcher et al (2021)    
            rho_sno = sno(i,j)/(si(i,j)/1.0e3)
            rho_sno = rho_sno*2. + 10.*max(0.,rho_sno-0.15)
            vis_blsn = (5.217*rho_sno*r_bar**1.011)/(1.82*c_z*c_alpha)
            BETABLSN = 3.912/(vis_blsn/1000.0)
            ! print to ensure quality
            !print *, "z_r", z_r
            !print *, "Q_s", Q_s
            !print *, "C_r", C_r
            !print *, "c_z", c_z
            !print *, "c_alpha", c_alpha
            !print *, "sno/SWE", sno(i,j)
            !print *, "si/SNOD", si(i,j)/1.0e3
            !print *, "rho_sno", rho_sno
            !print *, "vis_blsn", vis_blsn
            !print *, "BETABLSN", BETABLSN
        else
            BETABLSN = 0
            !print *, "BETABLSN", BETABLSN
        end if

! Key calculation of attenuation from each hydrometeor type (cloud, snow, graupel, rain, ice)
        BETAV=COEFFC*CONCFC**EXPONFC                            &
             + coef_SNOW*CONCFP**EXPONFP                        &
             + COEFLC*CONCLC**EXPONLC + COEFLP*CONCLP**EXPONLP    &
             + coeffg*concfg**exponfg  +1.E-10

! Addition of attenuation from aerosols if option selected
        if(method == 2 .or. method == 3)then ! aerosol method
            BETAV = BETAV + aextc55(i,j,lm)*1000.
            if(method_blsn) then ! BLSN method, updated 8 August 2022 by Tim Corrie
                BETAV = BETAV + BETABLSN
            endif
        endif

!  Calculation of visibility based on hydrometeor and aerosols.  (RH effect not yet included.)
        VIS(I,J)=MIN(90.,CONST1/(BETAV+extcof55(i,j,lm)))      ! max of 90km

        if (vis(i,j)<vis_min) vis_min = vis(i,j)
        if (visrh<visrh_min) visrh_min = visrh

        if (visrh<vis(i,j)) visrh_lower = visrh_lower + 1


! -- Dec 2003 - Roy Rasmussen (NCAR) expression for night vs. day vis
!   1.609 factor is number of km in mile.
       vis_night = 1.69 * ((vis(i,j)/1.609)**0.86) * 1.609

       zen_fac = min(0.1,max(czen(i,j),0.))/ 0.1
       vis(i,j) = zen_fac * vis(i,j) + (1.-zen_fac)*vis_night

       if(method == 1 .or. method == 3)then ! RH method (if lower vis)
         vis(i,j) = min(vis(i,j),visrh)
       endif

        if (vis(i,j)<1.) vis1km_cnt = vis1km_cnt + 1
        if (vis(i,j)<3.) vis3km_cnt = vis3km_cnt + 1
        if (vis(i,j)<5.) vis5km_cnt = vis5km_cnt + 1
! convert vis from km to [m]
        vis(i,j) = vis(i,j) * 1000.

        endif !end checking undefined points
      ENDDO
      ENDDO

!      write (6,*)
!      write (6,*) ' Visibility diagnostics follow: ------------'
!      write (6,*) ' -------------------------------------------'
!      write (6,*)                                                   &
!       '                                   any vis  /  vis < 10 km '
!      write (6,*)'No. of grid pts with shear (lev4-1) > 4m/s',      &
!          shear4_cnt, shear4_cnt_lowvis
!      write (6,*)'No. of grid pts with shear (lev4-1) > 5m/s',      &
!          shear5_cnt, shear5_cnt_lowvis
!      write (6,*)'No. of grid pts with shear (lev4-1) > 6m/s',      &
!          shear8_cnt, shear8_cnt_lowvis
!      write (6,*)
!      write (6,*)'No. of grid pts with vis-RH < 10 km',             &
!          visrh10_cnt
!      write (6,*)'No. of grid pts with vis    <  1 km',             &
!          vis1km_cnt
!      write (6,*)'No. of grid pts with vis    <  3 km',             &
!          vis3km_cnt
!      write (6,*)'No. of grid pts with vis    <  5 km',             &
!         vis5km_cnt
!      write (6,*)
!      write (6,*)'Min vis-hydrometeor, vis-RH', vis_min, visrh_min
!
!      write (6,*)'No. of grid pts with visRH < vis(hydrometeor)',   & 
!          visrh_lower
!      write (6,*)'% grid pts with night/cos(zen) < 0.1',            &
!          float(night_cnt)/float(IM*JM),float(lowsun_cnt)/          &
!          float(IM*JM)
!      write (6,*)
!
      RETURN
      END
