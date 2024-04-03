!> @file
!> @brief module: masks declares values used in masks
      module masks
!
      implicit none
!
      REAL, ALLOCATABLE ::  HBM2(:,:),SM(:,:),SICE(:,:)     &
     &,GDLAT(:,:),GDLON(:,:),LMH(:,:),LMV(:,:)              &
     &,HTM (:,:,:),VTM (:,:,:),DX(:,:),DY(:,:)
!
      end module masks
