!> @file
!> @brief module: soil declares soil-related variables
      module soil
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      real, allocatable :: STC(:,:,:),SMC(:,:,:),SH2O(:,:,:)      &
             ,SLDPTH(:),RTDPTH(:),SLLEVEL(:)
      end module soil
