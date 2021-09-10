      module vrbls2d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      real, allocatable ::                                                   &
      U10   (:,:),AKMS  (:,:),AKHS  (:,:),THS   (:,:),QS(:,:)                &
      ,UZ0(:,:),VZ0(:,:),THZ0(:,:),QZ0(:,:)                                  &
      ,SNO   (:,:),TSHLTR   (:,:),QSHLTR(:,:), MRSHLTR(:,:)                  &
      ,V10(:,:),ACPREC(:,:),CUPREC(:,:),ANCPRC(:,:),CUPPT(:,:)               &
      ,SMSTAV(:,:),SSROFF(:,:),BGROFF(:,:),VEGFRC(:,:)                       &
      ,SHDMIN(:,:),SHDMAX(:,:),LAI(:,:)                                      &
      ,ACSNOW(:,:),ACSNOM(:,:),CMC(:,:),SST(:,:)                             &
      ,RSWIN(:,:),RLWIN(:,:),RLWTOA(:,:)                                     &
      ,LWDNBC(:,:),LWUPBC(:,:)                                               &
      ,TG(:,:),SFCSHX(:,:),PSLP(:,:),T700(:,:),Z500(:,:),Z700(:,:)           &
      ,SFCLHX(:,:),FIS(:,:),T500(:,:),Z1000(:,:),SLP(:,:)                    &
      ,CFRACL(:,:),CFRACM(:,:),CFRACH(:,:),ACFRST(:,:)                       &
      ,ACFRCV(:,:),NCFRST(:,:),NCFRCV(:,:),HBOT(:,:)                         &
      ,HTOP(:,:),ASWIN(:,:),ALWIN(:,:),ASWOUT(:,:)                           &
      ,ALWOUT(:,:),ASWTOA(:,:),ALWTOA(:,:),CZEN(:,:)                         &
      ,CZMEAN(:,:),SIGT4(:,:),RSWOUT(:,:),RADOT(:,:)                         &
      ,SMSTOT(:,:),PCTSNO(:,:),PSHLTR(:,:),TH10(:,:)                         &
      ,Q10(:,:),SR(:,:),PREC(:,:),SUBSHX(:,:)                                &
      ,SNOPCX(:,:),SFCUVX(:,:),SFCEVP(:,:),POTEVP(:,:)                       &
      ,Z0(:,:),USTAR(:,:),TWBS(:,:),QWBS(:,:)                                &
      ,SFCEXC(:,:),GRNFLX(:,:),SOILTB(:,:),F(:,:)                            &
      ,ALBEDO(:,:),CLDFRA(:,:),CPRATE(:,:),CNVCFR(:,:)                       &
      ,PBLH(:,:),PBLHGUST(:,:),HBOTD(:,:),HTOPD(:,:),HBOTS(:,:),HTOPS(:,:)   &
      ,CLDEFI(:,:),ALBASE(:,:),SI(:,:),LSPA(:,:)                             &
      ,RSWINC(:,:),VIS(:,:),PD(:,:),MXSNAL(:,:),MIXHT(:,:)                   &
      ,SNONC(:,:),EPSR(:,:),RSWTOA(:,:),TEQL(:,:)                            &
! HWRF additions
      ,MDLTAUX(:,:),MDLTAUY(:,:),CD10(:,:),CH10(:,:)  &
      ,ACSWUPT(:,:),SWDNT(:,:),ACSWDNT(:,:) &
! NAMB additions
      ,SNOAVG(:,:),PSFCAVG(:,:),T10AVG(:,:),AKHSAVG(:,:),AKMSAVG(:,:)        &
      ,T10M(:,:),U10MAX(:,:),V10MAX(:,:),u10h(:,:),v10h(:,:)                 &
      ,PRATE_MAX(:,:),FPRATE_MAX(:,:)                                        &
! GSD addition
      ,WSPD10MAX(:,:),W_UP_MAX(:,:),W_DN_MAX(:,:),REFD_MAX(:,:)              &
      ,UP_HELI_MAX(:,:),UP_HELI_MAX16(:,:),GRPL_MAX(:,:),QRMAX(:,:)          &
      ,UP_HELI(:,:),UP_HELI16(:,:),LTG1_MAX(:,:),LTG2_MAX(:,:),LTG3_MAX(:,:) &
      ,UP_HELI_MIN(:,:),UP_HELI_MIN16(:,:)                                   &
      ,UP_HELI_MAX02(:,:),UP_HELI_MIN02(:,:)                                 &
      ,UP_HELI_MAX03(:,:),UP_HELI_MIN03(:,:)                                 &
      ,REL_VORT_MAX(:,:),REL_VORT_MAX01(:,:),REL_VORT_MAXHY1(:,:)            &
      ,WSPD10UMAX(:,:),WSPD10VMAX(:,:)                                       &
      ,REFDM10C_MAX(:,:),HAIL_MAX2D(:,:),HAIL_MAXK1(:,:)                     &
      ,HAIL_MAXHAILCAST(:,:)                                                 &
      ,NCI_LTG(:,:),NCA_LTG(:,:),NCI_WQ(:,:),NCA_WQ(:,:)                     &
      ,NCI_REFD(:,:),NCA_REFD(:,:),RAINC_BUCKET1(:,:),RAINNC_BUCKET1(:,:)   &
      ,RAINC_BUCKET(:,:),RAINNC_BUCKET(:,:),SNOW_BUCKET(:,:)                 &
      ,GRAUP_BUCKET(:,:),PCP_BUCKET(:,:),ACGRAUP(:,:),ACFRAIN(:,:)           &
      ,SNOW_BUCKET1(:,:),GRAUP_BUCKET1(:,:),PCP_BUCKET1(:,:)                 &
      ,SNOWNC(:,:),GRAUPELNC(:,:),TMAX(:,:),W_MEAN(:,:)                      &
      ,TSNOW(:,:),QVG(:,:),QV2m(:,:),QVl1(:,:)                               &
      ,REFC_10CM(:,:), REF1KM_10CM(:,:), REF4KM_10CM(:,:)                    &
      ,SWRADmean(:,:),U10mean(:,:),V10mean(:,:),SPDUV10mean(:,:)             &
      ,SWNORMmean(:,:),SNFDEN(:,:),SNDEPAC(:,:),SWDDNI(:,:),SWDDIF(:,:)      &
      ,SWDNBC(:,:),SWDDNIC(:,:),SWDDIFC(:,:), SWUPBC(:,:), SWUPT(:,:)        &
      ,TAOD5502D(:,:),AERASY2D(:,:),AERSSA2D(:,:),MEAN_FRP(:,:)              &
      ,LWP(:,:),IWP(:,:)                                                     &
      ,INT_SMOKE(:,:),INT_AOD(:,:)                                           &
! add new fields for GFS
      ,SFCUX(:,:),SFCVX(:,:),SFCUXI(:,:), SFCVXI(:,:),AVGALBEDO(:,:),AVGCPRATE(:,:)                   &
      ,AVGPREC(:,:),PTOP(:,:),PBOT(:,:),AVGCFRACH(:,:)                       &
      ,AVGCFRACM(:,:),AVGCFRACL(:,:),AVGTCDC(:,:)                            &
      ,AUVBIN(:,:),AUVBINC(:,:)                                              &
      ,ptopl(:,:),pbotl(:,:),Ttopl(:,:)                                      &
      ,ptopm(:,:),pbotm(:,:),Ttopm(:,:)                                      &
      ,ptoph(:,:),pboth(:,:),Ttoph(:,:)                                      &
      ,sfcugs(:,:),sfcvgs(:,:),PBLCFR(:,:)                                   &
      ,cldwork(:,:),gtaux(:,:),gtauy(:,:),runoff(:,:)                        &
      ,maxtshltr(:,:),mintshltr(:,:),maxrhshltr(:,:)                         &
      ,minrhshltr(:,:),dzice(:,:),maxqshltr(:,:),minqshltr(:,:)              &
      ,alwinc(:,:),alwoutc(:,:),alwtoac(:,:)                                 &
      ,aswinc(:,:),aswoutc(:,:),aswtoac(:,:),aswintoa(:,:)                   &
      ,smcwlt(:,:),suntime(:,:),fieldcapa(:,:)                               &
      ,avisbeamswin(:,:),avisdiffswin(:,:),airbeamswin(:,:)                  &
      ,airdiffswin(:,:),snowfall(:,:),acond(:,:),edir(:,:),ecan(:,:) &
      ,etrans(:,:),esnow(:,:),avgedir(:,:),avgecan(:,:),avgetrans(:,:)&
      ,avgesnow(:,:),avgpotevp(:,:),avgprec_cont(:,:),avgcprate_cont(:,:)&
      ,ti(:,:),aod550(:,:),du_aod550(:,:),ss_aod550(:,:),su_aod550(:,:)      &
      ,bc_aod550(:,:),oc_aod550(:,:)
      integer, allocatable :: IVGTYP(:,:),ISLTYP(:,:),ISLOPE(:,:) &
      ,IEQL(:,:)
! Add 2d aerosol diagnosis fields for GOCART (NGAC)
      real, allocatable ::                                                   &
       DUSMASS(:,:),DUCMASS(:,:),DUSMASS25(:,:),DUCMASS25(:,:)               &
      ,SUSMASS(:,:),SUCMASS(:,:),SUSMASS25(:,:),SUCMASS25(:,:)               &
      ,OCSMASS(:,:),OCCMASS(:,:),OCSMASS25(:,:),OCCMASS25(:,:)               &
      ,BCSMASS(:,:),BCCMASS(:,:),BCSMASS25(:,:),BCCMASS25(:,:)               &
      ,SSSMASS(:,:),SSCMASS(:,:),SSSMASS25(:,:),SSCMASS25(:,:)               &
      ,DUSTCB(:,:),SSCB(:,:),OCCB(:,:),BCCB(:,:),SULFCB(:,:)                 &
      ,DUSTALLCB(:,:),SSALLCB(:,:),DUSTPM(:,:),SSPM(:,:),PP25CB(:,:)         &
      ,PP10CB(:,:)!lzhang, add for FV3-Chem
 
!
      end module vrbls2d
