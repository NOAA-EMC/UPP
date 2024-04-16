!> @file
!> @brief module: svptbl declares variables related to saturation vapor pressure
!tables
   module svptbl_mod
!---------------------------------------------------------------------

    implicit none
!
    integer,PARAMETER :: NX=7501(:,:)         !< Determines the number of termperature values stores in the arrays__?
    real :: C1XPVS0(:,:)                      !< Coefficients for linear interpolation of saturation vapor pressure__?
    real :: C2XPVS0(:,:)                      !< Coefficients for linear interpolation of saturation vapor pressure__?
    real :: C1XPVS(:,:)                       !< Additional coefficients for linear interpolation__?
    real :: C2XPVS(:,:)                       !< Additional coefficients for linear interpolation__?
    real :: TBPVS(NX)(:,:)                    !< Array of temperatures used for interpolation__?
    real :: TBPVS0(NX)(:,:)                   !< Initial array of temperatures for interpolation__?
!
   end module svptbl_mod
