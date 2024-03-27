!> @file
!> @brief module: svptbl declares variables related to saturation vapor pressure
!tables
   module svptbl_mod
!---------------------------------------------------------------------

    implicit none
!
    integer,PARAMETER :: NX=7501
    real C1XPVS0,C2XPVS0
    real C1XPVS,C2XPVS
    real TBPVS(NX),TBPVS0(NX)
!
   end module svptbl_mod
