!> @file
!> @brief GRIDSPEC_mod assigns values to variables that define the model grid.
!----------------------------------------------------------------------------
     module GRIDSPEC_mod
!
!      COMMON /GRIDSPEC/
!     & DXVAL,DYVAL,CENLAT,CENLON,TRUELAT1,TRUELAT2,LATSTART,LONSTART
!     &,MAPTYPE,LATLAST,LONLAST,STANDLON
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none 
!
       integer DXVAL       !< grid cell size in x direction (can be degrees/meters)
       integer DYVAL       !< grid cell size in y direction (can be degrees/meters)
       integer CENLAT      !< cenlat: center latitude of grid
       integer CENLON      !< cenlon: center longitude of grid
       integer TRUELAT1    !< first latitude from the pole at which the secant cone cuts the sphere
       integer TRUELAT2    !< second latitude from the pole at which the secant cone cuts the sphere
       integer LATSTART    !< latitude of first grid point (lower left corner latitude)
       integer LONSTART    !< longitude of first grid point (lower left corner longitude)
       integer LATLAST     !< latitude of last grid point (upper right corner latitude)
       integer LONLAST     !< longitude of last grid point (upper right corner longitude)
       integer LATSTART_R  !< latitude of first grid point (lower left corner latitude)
       integer LONSTART_R  !< longitude of first grid point (lower left corner longitude)
       integer LATLAST_R   !< latitude of last grid point (upper right corner latitude)
       integer LONLAST_R   !< longitude of last grid point (upper right corner longitude)
       integer latnw       !< upper left corner latitude
       integer lonnw       !< upper left corner longitude
       integer latse       !< lower right corner latitude
       integer lonse       !< lower right corner longitude
       integer MAPTYPE     !< grid projection
       integer STANDLON    !< longitude of meridian parallel to y-axis (hardcoded as cenlon)
       integer latstartv   !< latitude of first grid point (lower left corner latitude)
       integer cenlatv     !< center latitude of grid
       integer lonstartv   !< longitude of first grid point (lower left corner longitude)
       integer cenlonv     !< center longitude of grid
       integer latlastv    !< latitude of last grid point (upper right corner latitude)
       integer lonlastv    !< longitude of last grid point (upper right corner longitude)
       real    PSMAPF      !< map scale factor
       character(len=1) gridtype !< type of grid staggering as in Arakawa grids (Arakawa-A through Arakawa-E)
!
     end module GRIDSPEC_mod
