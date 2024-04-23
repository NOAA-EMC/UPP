!> @file
!> @brief module: CMASSI_mod defines variables related to mass and precipitation
  module CMASSI_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   implicit none
!
!-----------------------------------------------------------------------
      REAL, PARAMETER :: DMImin=.05e-3, &    !< Minimum mean mass of precipitation ice particles
      DMImax=1.e-3,       &                  !< Maximum mean mass of precipitation ice particles 
      XMImin=1.e6*DMImin, &                  !< Minimum mean mass of precipitation ice particles multiplied by 1.0e6
      XMImax=1.e6*DMImax                     !< Maximum mean mass of precipitation ice particles multiplied by 1.0e6
      
      INTEGER, PARAMETER :: MDImin=XMImin, &
      MDImax=XMImax

!-----------------------------------------------------------------------

      REAL MASSI(MDImin:MDImax)        !< Mean mass of precipitation ice particles as functions of their mean size (in microns)

!--- Mean rain drop diameters vary from 50 microns to 1000 microns
! DMRmax definition is moved to microinit and has different values depending on imp_physics
      
      REAL, PARAMETER :: DMRmin=.05E-3, &          !< Minimum mean rain drop diameter 
      DelDMR=1.E-6,           &                    !< Delta mean rain drop diameter _____? 
      XMRmin=1.E6*DMRmin,     &                    !< Minimum mean rain drop diameter multiplied by 1.0E6
      N0r0=8.E6,              &                    !< Assumed intercept (m**-4) of rain drops if drop diameters are between 0.2 and 1.0 mm
      N0rmin=1.e4                                  !< Minimum intercept (m**-4) for rain drops

      REAL DMRmax    &                             !< Maximum mean rain drop diameter
      ,XMRmax                                      !< Maximum mean rain drop diameter multiplied by 1.0E6
      
      INTEGER, PARAMETER :: MDRmin=XMRmin          !< Minimum mean rain drop diameter multiplied by 1.0E6
      INTEGER MDRmax                               !< Maximum mean rain drop diameter multiplied by 1.0E6

! 
! --- Various rain lookup tables
! 

      REAL RQR_DRmin   &                       !< Minimum value of rain drop diameter lookup table
      ,RQR_DRmax       &                       !< Maximum value of rain drop diameter lookup table
      ,CN0r0           &                       !< Coefficient related to rain drop diameter
      ,CN0r_DMRmin     &                       !< Minimum value of rain drop diameter in a lookup table for a specific coefficient_
      ,CN0r_DMRmax                             !< Maximum value of rain drop diameter in a lookup table for a specific coefficient

! --- Other important parameters
!     (NLImax, FLARGE2 are used for the older version of the microphysics)
! 
      REAL T_ICE       &                       !< Temperature (C) threshold at which all remaining liquid water is glaciated to ice
      ,NLImax          &                       !< Maximum number of liquid particles in ice
      ,FLARGE2         &                       !< A parameter related to large-scale phenomena or processes___?
      ,TRAD_ice                                !< Traditional ice-related parameter

  end module  CMASSI_mod
