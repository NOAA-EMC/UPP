!> @file
!> @brief CMASSI_mod defines variables related to mass and precipitation
!> See CCPP Ferrier-Aligo microphysics modules for more information
!> @defgroup CMASSI
!> Defines variables related to mass and precipitation
  module CTLBLK_mod

  module CMASSI_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   implicit none
!
!-----------------------------------------------------------------------
      REAL, PARAMETER :: DMImin=.05e-3, &    !< Minimum mean mass of precipitation ice particles.
      DMImax=1.e-3,       &                  !< Maximum mean mass of precipitation ice particles.
      XMImin=1.e6*DMImin, &                  !< Minimum mean mass of precipitation ice particles (in microns).
      XMImax=1.e6*DMImax                     !< Maximum mean mass of precipitation ice particles (in microns).
      
      INTEGER, PARAMETER :: MDImin=XMImin, & !< Minimum mean diameter of precipitation ice particles.
      MDImax=XMImax                          !< Maximum mean diameter of precipitation ice particles.

!-----------------------------------------------------------------------

      REAL MASSI(MDImin:MDImax)        !< Mean mass of precipitation ice particles as functions of their mean size (in microns).

!--- Mean rain drop diameters vary from 50 microns to 1000 microns
!> DMRmax definition is moved to microinit and has different values depending on imp_physics
      
      REAL, PARAMETER :: DMRmin=.05E-3, &          !< Minimum mean rain drop diameter (0.05 mm).
      DelDMR=1.E-6,           &                    !< One-micron interval (Lookup tables store solutions at 1 micron intervals [DelDMR] of mean rain drop diameter.).
      XMRmin=1.E6*DMRmin,     &                    !< Minimum mean rain drop diameter (in microns).
      N0r0=8.E6,              &                    !< Assumed intercept (m**-4) of rain drops if drop diameters are between 0.2 and 1.0 mm.
      N0rmin=1.e4                                  !< Minimum intercept (m**-4) for rain drops.

      REAL DMRmax    &                             !< Maximum mean rain drop diameter.
      ,XMRmax                                      !< Maximum mean rain drop diameter.
      
      INTEGER, PARAMETER :: MDRmin=XMRmin          !< Minimum mean rain drop diameter (in microns).
      INTEGER MDRmax                               !< Maximum mean rain drop diameter (in microns).

! 
!--- Various rain lookup tables
! 

      REAL RQR_DRmin   &                       !< Rain content (kg/m**3) for mean drop diameter of .05 mm.
      ,RQR_DRmax       &                       !< Rain content (kg/m**3) for mean drop diameter of 1.0 mm.
      ,CN0r0           &                       !< Constant derived from N0r0.
      ,CN0r_DMRmin     &                       !< Minimum (starting) value for rain lookup tables for mean rain drop diameters.
      ,CN0r_DMRmax                             !< Maximum (ending) value for rain lookup tables for mean rain drop diameters.

!--- Other important parameters
!     (NLImax, FLARGE2 are used for the older version of the microphysics)
! 
      REAL T_ICE       &                       !< Temperature (C) threshold at which all remaining liquid water is glaciated to ice.
      ,NLImax          &                       !< Maximum number concentrations (m**-3) of large ice (snow/graupel/sleet).
      ,FLARGE2         &                       !< Set in MICROINIT.F (no longer used).
      ,TRAD_ice                                !< Defined as 0.5*T_ICE+TFRZ, or 253.15K, in other routines. Possibly refers to thermal radiation of ice or ice nucleation temperature ? 


  end module  CMASSI_mod
