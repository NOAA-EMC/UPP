!> @file
!> @brief module: soil declares soil-related variables
      module soil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      real, allocatable :: STC(:,:,:) &  !< Soil temperature 
       ,SMC(:,:,:) &     !< Volumetric soil moisture
       ,SH2O(:,:,:) &    !< Liquid volumetric soil moisture
       ,SLDPTH(:) &      !< Thickness of each soil layer
       ,RTDPTH(:) &      !< Depth of the bottom of the root zone
       ,SLLEVEL(:)       !< Soil level
      end module soil
