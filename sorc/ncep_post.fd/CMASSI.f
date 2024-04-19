!> @file
!> @brief module: CMASSI_mod defines variables related to mass and precipitation
  module CMASSI_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   implicit none
!
!-----------------------------------------------------------------------
!      Real Parameters.
      real, parameter :: DMImin = 0.05e-3          !< Minimum mean mass of precipitation ice particles
      real, parameter :: DMImax = 1.0e-3           !< Maximum mean mass of precipitation ice particles 
      real, parameter :: XMImin = 1.0e6 * DMImin   !< Minimum mean mass of precipitation ice particles multiplied by 1.0e6
      real, parameter :: XMImax = 1.0e6 * DMImax   !< Maximum mean mass of precipitation ice particles multiplied by 1.0e6
      real MASSI(MDImin:MDImax)                    !< Mean mass of precipitation ice particles as functions of their mean size 
      real, parameter :: DMRmin = 0.05E-3          !< Minimum mean rain drop diameter 
      real, parameter :: DelDMR = 1.0E-6           !< Delta mean rain drop diameter 
      real, parameter :: XMRmin = 1.0E6 * DMRmin   !< Minimum mean rain drop diameter multiplied by 1.0E6
      real, parameter :: N0r0 = 8.0E6              !< Number concentration of raindrops with a diameter greater than a certain threshhold
      real, parameter :: N0rmin = 1.0E4            !< Minimum value of concetraion of rain drops 
      real DMRmax                                  !< Maximum mean rain drop diameter 
      real XMRmax                                  !< Maximum mean rain drop diameter multiplied by 1.0E6
      real RQR_DRmin                               !< Minimum value of rain drop diameter lookup table
      real RQR_DRmax                               !< Maximum value of rain drop diameter lookup table
      real CN0r0                                   !< Coefficient related to rain drop diameter
      real CN0r_DMRmin                             !< Minimum value of rain drop diameter in a lookup table for a specific coefficient_
      real CN0r_DMRmax                             !< Maximum value of rain drop diameter in a lookup table for a specific coefficient
      real T_ICE                                   !< Temperature of ice
      real NLImax                                  !< Maximum number of liquid particles in ice
      real FLARGE2                                 !< A parameter related to large-scale phenomena or processes___?
      real TRAD_ice                                !< Traditional ice-related parameter

!      Integer Parameters.
      integer, parameter :: MDImin = XMImin        !< Minimum mean diameter 
      integer, parameter :: MDImax = XMImax        !< Maximum mean diameter 
      integer MDRmin = XMRmin                      !< Minimum mean rain drop diameter multiplied by 1.0E6
      integer MDRmax                               !< Maximum mean rain drop diameter multiplied by 1.0E6

!
  end module  CMASSI_mod
