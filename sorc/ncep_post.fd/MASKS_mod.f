!> @file
!> @brief module: masks declares values used in masks
      module masks
!
      implicit none
!
      REAL, ALLOCATABLE :: HBM2(:,:) &         !< _____
      ,SM(:,:) &         !< Land-sea mask ?
      ,SICE(:,:) &       !< Sea ice mask
      ,GDLAT(:,:) &      !< Grid latitude
      ,GDLON(:,:) &      !< Grid longitude
      ,LMH(:,:) &        !< Mass point at ETA surface mask ?
      ,LMV(:,:) &        !< Velocity point at ETA surface mask ?  
      ,HTM (:,:,:) &     !< Height topography mask array
      ,VTM (:,:,:) &     !< _____
      ,DX(:,:) &         !< Grid spacing in the x-direction ?
      ,DY(:,:)           !< Grid spacing in the y-direction ?
!
      end module masks
