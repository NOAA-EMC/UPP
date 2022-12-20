!   01-10-22  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-04-17  BALDWIN  - MODIFIED TO INCLUDE ALL 3D ARRAYS
!   11-10-18  SARAH LU - MODIFIED TO INCLUDE GOCART AEROSOLS
!   22-09-18  Li(Kate) Zhang - MODIFIED TO INCLUDE new NASA GOCART AEROSOLS of NO3 and NH4
      module vrbls4d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
      real, allocatable :: DUST(:,:,:,:)        ! dust
      real, allocatable :: SALT(:,:,:,:)        ! sea salt
      real, allocatable :: SOOT(:,:,:,:)        ! black carbon
      real, allocatable :: WASO(:,:,:,:)        ! organic carbon
      real, allocatable :: SUSO(:,:,:,:)        ! sulfate
      real, allocatable :: NO3(:,:,:,:)         ! no3
      real, allocatable :: NH4(:,:,:,:)         ! nh4
      real, allocatable :: SMOKE(:,:,:,:)
      real, allocatable :: FV3DUST(:,:,:,:)
      real, allocatable :: PP25(:,:,:,:)        ! PP25
      real, allocatable :: PP10(:,:,:,:)        ! PP10
!
      end module vrbls4d
