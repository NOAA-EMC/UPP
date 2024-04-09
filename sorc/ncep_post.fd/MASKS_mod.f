!> @file
!> @brief module: masks declares values used in masks
      module masks
!
      implicit none
!
      REAL, ALLOCATABLE :: HBM2(:,:) &         !<
      ,SM(:,:) &         !< Soil moisture___?
      ,SICE(:,:) &       !< Sea ice mask
      ,GDLAT(:,:) &      !< Grid latitude___?
      ,GDLON(:,:) &      !< Grid longitude___?
      ,LMH(:,:) &        !< Land model height___?
      ,LMV(:,:) &        !< Land model vertical___?  
      ,HTM (:,:,:) &     !<  Horizontal temperature momentum___?
      ,VTM (:,:,:) &     !< Virtual temperature___?
      ,DX(:,:) &         !< Grid spacing in the x-direction___?
      ,DY(:,:) &         !< Grid spacing in the y-direction___?
!
      end module masks
