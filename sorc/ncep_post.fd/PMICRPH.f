!> @file
!> @brief PMICRPH_mod assigns constants related to microphysics 
!> (computed in paramr.f)
!-----------------------------------------------------------------
     module PMICRPH_mod

! -----  Constants related to microphysics
!          -- computed in paramr.f

!      REAL ABER1(31),ABER2(31)
! LOOKUP TABLE FOR A1 AND A2 IN BERGERON PROCESS

      implicit none

      REAL PI        !< Mathematical contant pi (3.14159...)
      REAL RON       !< slope intercept for rain (original M-P value)
      REAL SON       !< slope intercept for snow
      REAL GON       !< slope intercept for graupel (original M-P value)
      REAL BR        !< exponent for rain in fall speed (V(D)=A*D^B)
      REAL BS        !< exponent for snow in fall speed (V(D)=A*D^B)
      REAL BG        !< exponent for graupel in fall speed (V(D)=A*D^B)
      REAL DRAIN     !< density of rain
      REAL DSNOW     !< density of snow
      REAL DGRAUPEL  !< density of graupel
      REAL RON2      !< slope intercept for rain (modified)
      REAL DIACE_min !< minimum mass of ice
      REAL drain2    !< square of rain density
      REAL dsnow2    !< square of snow density
      REAL TOPR      !< top of slope (numerator Marshall-Palmer slope parameter) for rain
      REAL TOPS      !< top of slope (numerator Marshall-Palmer slope parameter) for snow
      REAL TOPG      !< top of slope (numerator Marshall-Palmer slope parameter) for graupel
      REAL ARAIN     !< A in fall speed for rain
      REAL ASNOW     !< A in fall speed for snow
      REAL AGRAUPEL  !< A in fall speed for graupel
      REAL TNO       !< constant in Cooper and Fletcher curves
      REAL ATO       !< constant in Cooper Fletcher curves (not used)
      REAL XSMAX     !< autoconversion to snow
      REAL BERC1     !< constant for Bergeron process (not used)
      REAL BP        !< B Prime (B') constant for computing freezing rate of cloud droplets
      REAL AP        !< A Prime (A') constant for computing freezing rate of cloud droplets
      REAL CNP       !< constant for computing cloud drop shape parameter
      REAL FRD1      !< Related to freezing of rain droplets (Lin, et al., 45)
      REAL FRA1      !< Related to freezing of rain droplets (Lin, et al., 45)
      REAL EFIS      !< collection efficiency of cloud ice by snow
      REAL EFIR      !< collection efficiency of cloud ice by rain
      REAL EFSR      !< collection efficiency of snow by rain
      REAL EFCS      !< collection efficiency of cloud water by snow
      REAL EFGI      !< collection efficiency of cloud ice by graupel
      REAL EFGC      !< collection efficiency of cloud water by graupel
      REAL EFGR      !< collection efficiency of graupel by rain
      REAL EFGS      !< collection efficiency of graupel by snow
      REAL EFCR      !< collection efficiency of cloud water by rain
      REAL ACRIS     !< Related to collection of cloud ice by snow
      REAL BACRIS    !< Related to collection of cloud ice by snow
      REAL CIR       !< collection of cloud ice by rain
      REAL CIRF      !< rate at which rain is frozen by collision with cloud ice
      REAL cpiacr0   !< constant for PIACR (not used)
      REAL cpiacr1   !< constant for PIACR (not used)
      REAL cpiacr2   !< constant for PIACR (not used)
      REAL cpiacr3   !< constant for PIACR (not used)
      REAL FRAIN     !< mean fall speed of rain
      REAL FSNOW     !< mean fall speed of snow
      REAL FGRAUPEL  !< mean fall speed of graupel
      REAL CSR       !< collection of snow by rain
      REAL CRS       !< collection of rain by snow
      REAL ACRCS     !< Related to collection of cloud water by snow using old particle size distribution
      REAL BACRCS    !< Related to collection of cloud water by snow using old particle size distribution
      REAL RMC       !< constant - no longer used/supported
      REAL ACRLS     !< Related to loss of snow due to collision with cloud water 
      REAL BACLS     !< Related to loss of snow due to collision with cloud water 
      REAL ACRCG     !< Related to collection of cloud water by graupel using old particle size distribution
      REAL BACRCG    !< Related to collection of cloud water by graupel using old particle size distribution
      REAL ACRIG     !< Related to collection of cloud ice by graupel
      REAL BACRIG    !< Related to collection of cloud ice by graupel
      REAL CRG       !< collection of rain by graupel
      REAL CSG       !< collection of snow by graupel
      REAL DEPG1     !< Depositional growth of graupel
      REAL DEPG2     !< Depositional growth of graupel
      REAL DEPG3     !< Depositional growth of graupel
      REAL DEPG4     !< Depositional growth of graupel
      REAL DEPS1     !< Depositional growth of snow
      REAL DEPS2     !< Depositional growth of snow
      REAL DEPS3     !< Depositional growth of snow
      REAL DEPS4     !< Depositional growth of snow
      REAL ACRCR     !< Related to collection of cloud water by rain
      REAL BACRCR    !< Related to collection of cloud water by rain
      REAL DEPR1     !< Depositional growth of rain
      REAL DEPR2     !< Depositional growth of rain
      REAL DEPR3     !< Depositional growth of rain
      REAL DEPR4     !< Depositional growth of rain
      REAL PSM1      !< Related to melting of snow
      REAL PSM2      !< Related to melting of snow
      REAL PSM3      !< Related to melting of snow
      REAL PSM4      !< Related to melting of snow
      REAL PGM1      !< Related to melting of graupel
      REAL PGM2      !< Related to melting of graupel
      REAL PGM3      !< Related to melting of graupel
      REAL PGM4      !< Related to melting of graupel
      REAL CW        !< constant for enhanced melting of graupel by rain and cloud water
      REAL HGFR      !< constant for homogeneous freezing of cloud droplets
      REAL XM01      !< constant used to calculate the minimum mass of ice
      REAL CNP1      !< Not used/no longer supported
      REAL DICE      !< density of ice
      REAL C1        !< aggregation of cloud ice
      REAL ALPHA1    !< constant used to calculate collection of snow by rain
      REAL BETA1     !< constant used to calculate collection of snow by rain
      REAL GAMMA3    !< constant used to calculate collection of snow by rain
!jmb--removed INT0 frm the real declaration since declared integer blo
      REAL CONST1A   !< constant for variable ‘son’ (slope intercept for snow)
      REAL CONST1B   !< constant for variable ‘son’ (slope intercept for snow)
      REAL XM0S      !< minimum mass of snow
      REAL XR0S      !< smallest size of snow
      REAL XM0G      !< minimum mass of graupel
      REAL ACRCS_new    !< Related to collection of cloud water by snow using new particle size distribution for snow (Roy R, Jul 99)
      REAL BACRCS_new   !< Related to collection of cloud water by snow using new particle size distribution for snow (Roy R, Jul 99)
      REAL ACRCG_new    !< Related to collection of cloud water by graupel using new particle size distribution for graupel (Roy R, Jul 99)
      REAL BACRCG_new   !< Related to collection of cloud water by graupel using new particle size distribution for graupel (Roy R, Jul 99)
      REAL const_ns1 !< constant for variable ‘son’ (slope intercept for snow)
      REAL const_ns2 !< constant for variable ‘son’ (slope intercept for snow)
      REAL const_ng1 !< constant for variable ‘gon’ (slope intercept for graupel)
      REAL const_ng2 !< constant for variable ‘gon’ (slope intercept for graupel)
      REAL xr0g      !< smallest size of graupel
      REAL r1        !< minimum value for hydrometeor mixing ratios
      REAL slor_r1   !< inverse slope value when rain mixing ratio is small 
      REAL slos_r1   !< inverse slope value when snow mixing ratio is small 
      REAL slog_r1   !< inverse slope value when graupel mixing ratio is small 
      REAL rho_not   !< Standard density (p/RT, values from ICAO standard atmosphere) used in computing density correction to fall speeds
      REAL qck1      !< Constant - no longer used/supported
      REAL qcth      !< Constant - no longer used/supported
      REAL ron_min   !< minimum allowed value for ron
      REAL qr0       !< center value of rain mixing ratio for transition from M-P slope-intercept for drizzle formed by a collision-coalescence process to M-P slope-intercept for traditional rain
      REAL delqr0    !< governs the sharpness of qr0 transition: small delt_qr0 makes the transition sharper
      REAL const1r   !< constant for variable ‘ron’ (slope intercept for rain)
      REAL const2r   !< constant for variable ‘ron’ (slope intercept for rain)
      REAL xnu       !< cloud drop shape parameter

      INTEGER INT0   !< constant for Bergeron process

     end module PMICRPH_mod
