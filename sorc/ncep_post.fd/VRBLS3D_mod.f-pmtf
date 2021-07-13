!   01-10-22  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-04-17  BALDWIN  - MODIFIED TO INCLUDE ALL 3D ARRAYS
!   11-10-18  SARAH LU - MODIFIED TO INCLUDE AEROSOL OPTICAL PROPERTIES
!   11-12-15  SARAH LU - MODIFIED TO INCLUDE AEROSOL DIAG FIELDS
!   12-01-06  SARAH LU - MODIFIED TO INCLUDE AIR DENSITY AND LAYER THICKNESS
!   15-07-02  SARAH LU - MODIFIED TO INCLUDE SCATTERING AEROSOL OPTICAL THICKNESS
      module vrbls3d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      real, allocatable :: UH(:,:,:),VH(:,:,:),WH(:,:,:)             &
      ,U(:,:,:),V(:,:,:),T(:,:,:),Q(:,:,:)                           &
      ,CWM(:,:,:),Q2(:,:,:),PMID(:,:,:),PMIDV(:,:,:)                 &
      ,PINT(:,:,:),ALPINT(:,:,:),ZMID(:,:,:)                         &
      ,ZINT(:,:,:),OMGA(:,:,:)                                       &
      ,T_ADJ(:,:,:)                                                  &
      ,F_ice(:,:,:),F_rain(:,:,:),F_RimeF(:,:,:)                     &
      ,QQW(:,:,:), QQI(:,:,:), QQR(:,:,:), QQS(:,:,:), QQG(:,:,:)    &
      ,QQNW(:,:,:), QQNI(:,:,:),QQNR(:,:,:),QC_BL(:,:,:), QRIMEF(:,:,:) &
      ,CFR(:,:,:), DBZ(:,:,:), DBZR(:,:,:), DBZI(:,:,:), DBZC(:,:,:) &
      ,TTND(:,:,:),RSWTT(:,:,:),RLWTT(:,:,:), REF_10CM(:,:,:)        &
      ,EXCH_H(:,:,:),TRAIN(:,:,:),TCUCN(:,:,:),EL_PBL(:,:,:)         &
      ,MCVG(:,:,:),EXTCOF55(:,:,:),NLICE(:,:,:),CFR_RAW(:,:,:)       &
!! Wm Lewis: added
      ,NRAIN(:,:,:)                                                  &
      ,radius_cloud(:,:,:),radius_ice(:,:,:),radius_snow(:,:,:)      &
! KRS Add HWRF fields     
      ,REFL_10CM(:,:,:)             &
! Add GFS fields     
      ,O3(:,:,:),O(:,:,:),O2(:,:,:)              &
! Add GFS D3D fields
      ,vdifftt(:,:,:)         &
      ,tcucns(:,:,:)          &
      ,vdiffmois(:,:,:)       &
      ,dconvmois(:,:,:)       &
      ,sconvmois(:,:,:)       &
      ,nradtt(:,:,:)          &  
      ,o3vdiff(:,:,:)         &
      ,o3prod(:,:,:)          &
      ,o3tndy(:,:,:)          &
      ,mwpv(:,:,:)            &
      ,unknown(:,:,:)         &
      ,vdiffzacce(:,:,:)      &
      ,zgdrag(:,:,:)          &
      ,cnvctummixing(:,:,:)   &
      ,vdiffmacce(:,:,:)      &
      ,mgdrag(:,:,:)          &
      ,cnvctvmmixing(:,:,:)   &
      ,ncnvctcfrac(:,:,:)     &
      ,cnvctumflx(:,:,:)      &
      ,cnvctdmflx(:,:,:)      &  
      ,cnvctdetmflx(:,:,:)    &
      ,cnvctzgdrag(:,:,:)     &
      ,cnvctmgdrag(:,:,:)     &   
      ,QQNWFA(:,:,:)          &
      ,QQNIFA(:,:,:)          &
      ,TAOD5503D(:,:,:)       &
      ,AEXTC55(:,:,:)         &
!
! Add aerosol optical properties for GOCART (NGAC)
      ,ext(:,:,:), asy(:,:,:)           &
      ,ssa(:,:,:), sca(:,:,:)           &
! Add aerosol diagnosis fields for GOCART (NGAC)
      ,duem(:,:,:), dusd(:,:,:)         &
      ,dudp(:,:,:), duwt(:,:,:)         &
      ,dusv(:,:,:), sssv(:,:,:)         &
      ,suem(:,:,:), susd(:,:,:)         &
      ,sudp(:,:,:), suwt(:,:,:)         &
      ,ssem(:,:,:), sssd(:,:,:)         &
      ,ssdp(:,:,:), sswt(:,:,:)         &
      ,ocem(:,:,:), ocsd(:,:,:)         &
      ,ocdp(:,:,:), ocwt(:,:,:)         &
      ,ocsv(:,:,:), bcsv(:,:,:)         &
      ,bcem(:,:,:), bcsd(:,:,:)         &
      ,bcdp(:,:,:), bcwt(:,:,:)         &
! Add air density and thickness for GOCART (NGAC)
      ,dpres(:,:,:),rhomid(:,:,:)       &  

! Add NCAR GFIP ICING
      ,icing_gfip(:,:,:),icing_gfis(:,:,:) &
! Add NCAR GTG turbulence
      ,catedr(:,:,:),mwt(:,:,:),gtg(:,:,:) &

! AQF
      ,aacd(:,:,:), aalj(:,:,:)                    &
      ,aalk1j(:,:,:), aalk2j(:,:,:)                &
      ,abnz1j(:,:,:), abnz2j(:,:,:), abnz3j(:,:,:) &
      ,acaj(:,:,:), acet(:,:,:)                    &
      ,acli(:,:,:), aclj(:,:,:), aclk(:,:,:)       &
      ,acors(:,:,:), acro_primary(:,:,:)           &
      ,acrolein(:,:,:), aeci(:,:,:)                &
      ,aecj(:,:,:), afej(:,:,:)                    &
      ,aglyj(:,:,:)                                &
      ,ah2oi(:,:,:), ah2oj(:,:,:), ah2ok(:,:,:)    &
      ,ah3opi(:,:,:), ah3opj(:,:,:), ah3opk(:,:,:) &
      ,aiso1j(:,:,:), aiso2j(:,:,:), aiso3j(:,:,:) &
      ,aivpo1j(:,:,:), akj(:,:,:)                  &
      ,ald2(:,:,:), ald2_primary(:,:,:)            &
      ,aldx(:,:,:)                                 &
      ,alvoo1i(:,:,:), alvoo1j(:,:,:)              &
      ,alvoo2i(:,:,:), alvoo2j(:,:,:)              &
      ,alvpo1i(:,:,:), alvpo1j(:,:,:)              &
      ,amgj(:,:,:), amnj(:,:,:)                    &
      ,amgk(:,:,:), akk(:,:,:), acak(:,:,:)        &
      ,anai(:,:,:), anaj(:,:,:), anak(:,:,:)       &
      ,anh4i(:,:,:), anh4j(:,:,:), anh4k(:,:,:)    &
      ,ano3i(:,:,:), ano3j(:,:,:), ano3k(:,:,:)    &
      ,aolgaj(:,:,:), aolgbj(:,:,:), aorgcj(:,:,:) &
      ,aomi(:,:,:), aomj(:,:,:)                    &
      ,aothri(:,:,:), aothrj(:,:,:)                &
      ,apah1j(:,:,:), apah2j(:,:,:), apah3j(:,:,:) &
      ,apomi(:,:,:), apomj(:,:,:)                  &
      ,apcsoj(:,:,:), aseacat(:,:,:), asij(:,:,:)  &
      ,aso4i(:,:,:), aso4j(:,:,:), aso4k(:,:,:)    &
      ,asoil(:,:,:), asqtj(:,:,:)                  &
      ,asomi(:,:,:), asomj(:,:,:)                  &
      ,asvoo1i(:,:,:), asvoo1j(:,:,:)              &
      ,asvoo2i(:,:,:), asvoo2j(:,:,:)              &
      ,asvoo3j(:,:,:)                              &
      ,asvpo1i(:,:,:), asvpo1j(:,:,:)              &
      ,asvpo2i(:,:,:), asvpo2j(:,:,:)              &
      ,asvpo3j(:,:,:)                              &
      ,atij(:,:,:)                                 &
      ,atol1j(:,:,:), atol2j(:,:,:), atol3j(:,:,:) &
      ,atoti(:,:,:), atotj(:,:,:), atotk(:,:,:)    &
      ,atrp1j(:,:,:), atrp2j(:,:,:)                &
      ,axyl1j(:,:,:), axyl2j(:,:,:), axyl3j(:,:,:) &
      ,o3mr(:,:,:), ozcon(:,:,:)                   &
      ,pm25ac(:,:,:), pm25at(:,:,:), pm25co(:,:,:) &
      ,pm25hp(:,:,:), pm25cl(:,:,:), pm25ec(:,:,:) &
      ,pm25na(:,:,:), pm25mg(:,:,:), pm25k(:,:,:)  &
      ,pm25ca(:,:,:), pm25nh4(:,:,:)               &
      ,pm25no3(:,:,:), pm25oc(:,:,:), pm25om(:,:,:)&
      ,pm25soil(:,:,:), pm25so4(:,:,:)             &
      ,pmtf(:,:,:), pm25unspec1(:,:,:)

      end module vrbls3d
