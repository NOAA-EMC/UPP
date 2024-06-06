!> @file
!> @brief module: masks declares values used in masks
      module masks
!
      implicit none
!
      REAL, ALLOCATABLE :: HBM2(:,:) &         !< _____
      ,SM(:,:) &         !< Sea mask
      ,SICE(:,:) &       !< _____ Sea ice mask
      ,GDLAT(:,:) &      !< Grid latitude
      ,GDLON(:,:) &      !< Grid longitude
      ,LMH(:,:) &        !< Topography indexes array
      ,LMV(:,:) &        !< _____ Topography indexes array - vertical ?  
      ,HTM (:,:,:) &     !< Height topography mask array
      ,VTM (:,:,:) &     !< _____
      ,DX(:,:) &         !< _____ Grid spacing in the x-direction
      ,DY(:,:)           !< _____ Grid spacing in the y-direction
!
      end module masks
