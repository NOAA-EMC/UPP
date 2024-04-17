!> @file
!> @brief module: masks declares values used in masks
      module masks
!
      implicit none
!
      REAL, ALLOCATABLE :: HBM2(:,:) &         !<
      ,SM(:,:) &         !< Soil moisture
      ,SICE(:,:) &       !< Sea ice mask
      ,GDLAT(:,:) &      !< Grid latitude
      ,GDLON(:,:) &      !< Grid longitude
      ,LMH(:,:) &        !< Land model height
      ,LMV(:,:) &        !< Land model vertical___?  
      ,HTM (:,:,:) &     !< Horizontal temperature momentum___?
      ,VTM (:,:,:) &     !< Virtual temperature___?
      ,DX(:,:) &         !< Grid spacing in the x-direction
      ,DY(:,:) &         !< Grid spacing in the y-direction
!
      end module masks
