!> @file
!> @brief SVPTBL declares variables related to saturation vapor pressure tables
   module svptbl_mod
!---------------------------------------------------------------------

    implicit none
!
    integer,PARAMETER :: NX=7501(:,:)         !< Table length
    real :: C1XPVS0(:,:)                      !< Coefficient 1 for saturation vapor pressure in TBPVS0
    real :: C2XPVS0(:,:)                      !< Coefficient 2 for saturation vapor pressure in TBPVS0
    real :: C1XPVS(:,:)                       !< Coefficient 1 for saturation vapor pressure in TBPVS
    real :: C2XPVS(:,:)                       !< Coefficient 2 for saturation vapor pressure in TBPVS
    real :: TBPVS(NX)(:,:)                    !< Table of saturation vapor pressure values
    real :: TBPVS0(NX)(:,:)                   !< Table of saturation vapor pressure values
!
   end module svptbl_mod
