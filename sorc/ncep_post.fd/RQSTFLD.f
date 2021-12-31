    module RQSTFLD_mod
!--------------------------------------------------------------------
! revision history
!   2011-02-06 Jun Wang add grib2 option
!   2011-10-18 Sarah Lu add GOCART aerosol fields
!   2011-12-18 Sarah Lu add GOCART aerosol optical properties, mass
!                           fields, and production/removal fluxes
!   2011-12-29 Sarah Lu add GOCART AOD at multiple channels
!   2012-01-06 Sarah Lu add GOCART SS, OC, BC, SU aerosols
!   2012-01-07 Sarah Lu add air density and dpres
!   2012-01-27 Sarah Lu use index 601-700 for GOCART
!   2012-01-30 Jun Wang add post available fields from xml file for grib2
!   2012-05-07 Tricia Slovacek use index 800-899 for satellite
!                       use index 900-949 for HWRF
!   2014-12-09 William Lewis added MSG/SEVIRI imager, 
!                      GOES-13 and GOES-15 imagers,
!                      and completed SSMI and SSMIS (F13-F20)
!   2019-04-01 Sharon Nebuda added GOES-16 GOES-17 ABI IR Channels 7-16
!   2019-04-22 Wen Meng increased model MXLVL to 500
!   2019-05-08 Wen Meng added continuous accumulated precipitations(417, 418,
!                       419).
!   2019-09-03 Jesse Meng added cape related variables for hrrr (950-957)
!   2021-08-31 Lin Zhu added ssmis-f17 channels 15-18 grib2 output
!   2021-11-03 Tracy Hertneky removed all grib1 settings
!--------------------------------------------------------------------

      implicit none
!
!     increase MXFLD each time you add a new field
      INTEGER, PARAMETER :: MXFLD=1200
      INTEGER, PARAMETER :: MXLVL=500
      
      CHARACTER*6 DATSET      
!
      LOGICAL RITEHD
!
      integer :: NFLD,IGET(MXFLD),                                &
                 LVLS(MXLVL,MXFLD),                               &
                 IDENT(MXFLD),IAVBLFLD(MXFLD),                    &
                 ID(25)
      integer :: num_post_afld
      integer,allocatable :: LVLSXML(:,:)

   end module RQSTFLD_mod
