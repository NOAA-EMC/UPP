!   01-10-22  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-04-17  BALDWIN  - MODIFIED TO INCLUDE ALL 3D ARRAYS
      module vrbls3d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      real, allocatable :: UH(:,:,:),VH(:,:,:),WH(:,:,:) &
      ,U(:,:,:),V(:,:,:),T(:,:,:),Q(:,:,:) &
      ,CWM(:,:,:),Q2(:,:,:),PMID(:,:,:),PMIDV(:,:,:) &
      ,PINT(:,:,:),ALPINT(:,:,:),ZMID(:,:,:) &
      ,ZINT(:,:,:),OMGA(:,:,:) &
      ,T_ADJ(:,:,:) &
      ,F_ice(:,:,:),F_rain(:,:,:),F_RimeF(:,:,:) &
      ,QQW(:,:,:), QQI(:,:,:), QQR(:,:,:), QQS(:,:,:), QQG(:,:,:) &
      ,CFR(:,:,:), DBZ(:,:,:), DBZR(:,:,:), DBZI(:,:,:), DBZC(:,:,:) &
      ,TTND(:,:,:),RSWTT(:,:,:),RLWTT(:,:,:) &
      ,EXCH_H(:,:,:),TRAIN(:,:,:),TCUCN(:,:,:),EL_PBL(:,:,:) &
      ,MCVG(:,:,:),EXTCOF55(:,:,:),NLICE(:,:,:) &
! Add GFS fields     
      ,O3(:,:,:)             &
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
!
! Add NCAR GFIP ICING
      ,icing_gfip(:,:,:)
!
      end module vrbls3d
