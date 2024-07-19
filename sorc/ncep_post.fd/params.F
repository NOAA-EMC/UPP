!> @file
!> @brief params_mod defines physical constants and smoothing coefficients
module params_mod

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
! last update: 09/06/2006
!
!     SMOOTHING PARAMETERS.
    integer, parameter :: KSMUD=4      !< No longer used/supported
    integer, parameter :: KSLPD=1      !< No longer used/supported
    integer, parameter :: KSSLP=2      !< No longer used/supported
    integer, parameter :: KSFI=2       !< No longer used/supported
    integer, parameter :: KST=2        !< No longer used/supported
    integer, parameter :: KSUV=2       !< No longer used/supported
    integer, parameter :: KSRH=2       !< No longer used/supported
    integer, parameter :: KSAV=2       !< No longer used/supported
    integer, parameter :: KSLI=2       !< No longer used/supported
    integer, parameter :: KSLP=2       !< No longer used/supported
    integer, parameter :: NSUMD=1      !< No longer used/supported
                                                                                   
    real, parameter :: SMTHA=-1.8E-4   !< No longer used/supported
    real, parameter :: SMTHB=16.0      !< No longer used/supported
!
! CONSTANTS.
    real, parameter :: CM1=2937.4               !< Constant used in Magnus Formula to calculate saturation water vapor pressure
    real, parameter :: CM2=4.9283               !< Constant used in Magnus Formula to calculate saturation water vapor pressure
    real, parameter :: CM3=23.5518              !< Constant used in Magnus Formula to calculate saturation water vapor pressure
    real, parameter :: D00=0.0                  !< Decimal form of 0
    real, parameter :: D1=0.1                   !< Decimal number - one tenth
    real, parameter :: D01=0.01                 !< Decimal number - one hundreth
    real, parameter :: D001=0.001               !< Decimal number - one thousandth
    real, parameter :: D50=0.50                 !< Decimal number 0.5
    real, parameter :: H999=999.                !< 999
    real, parameter :: H1=1.0                   !< 1.0
    real, parameter :: H2=2.0                   !< 2.0
    real, parameter :: H4=4.0                   !< 4.0
    real, parameter :: H99999=99999.            !< 99999
    real, parameter :: H1M12=1.E-12             !< 1 * 10^(-12)
    real, parameter :: H1000=1000.              !< 1000
    real, parameter :: H1M5=1.E-5               !< 1/100000 == 0.00001
    real, parameter :: D125=.125                !< Decimal number 0.125
    real, parameter :: D25=0.25                 !< Decimal number 0.25
    real, parameter :: H100=100.                !< 100
    real, parameter :: H10000=10000.            !< 10000
    real, parameter :: H10E5=1.E5               !< 100000
    real, parameter :: CAPA=0.28589641E0        !< R/Cp (universal gas constant over specific heat capacity of dry air at constant pressure)
    real, parameter :: D115=0.115               !< Decimal number 0.125
    real, parameter :: D608=0.608               !< Decimal number 0.608
    real, parameter :: D05=0.05                 !< Decimal number 0.05
    real, parameter :: D0065=0.0065             !< Decimal number 0.0065
    real, parameter :: H1E2=100.                !< No longer used/supported
    real, parameter :: H1E4=1.E4                !< No longer used/supported
    real, parameter :: H1M2=1.E-2               !< No longer used/supported
    real, parameter :: G=9.81                   !< Acceleration due to gravity
    real, parameter :: RD=287.04                !< Gas constant for dry air
    real, parameter :: ROG=RD/G                 !< RD over G
    real, parameter :: GI=1./G                  !< Inverse of acceleration due to gravity
    real, parameter :: OVERRC=1.10              !< Overrelaxation coefficient
    real, parameter :: AD05=OVERRC*D05          !< Overrelaxation coefficient * 0.05
    real, parameter :: CFT0=OVERRC-H1           !< A constant calculated by subtracting H1 from OVERRC
    real, parameter :: PI=3.141592653589793     !< Pi
    real, parameter :: DTR=3.1415926/180.       !< Degrees to radians conversion
    real, parameter :: RTD=1./DTR               !< Radians to degrees conversion
    real, parameter :: ERAD=6.371E6             !< Earth's radius in meters
    real, parameter :: CP=1004.6                !< Specific heat capacity of dry air at constant pressure (kJ/kg-K)
    real, parameter :: RCAPA=1./CAPA            !< Reciprocal of CAPA
    real, parameter :: P1000=1000.E2            !< 100,000
    real, parameter :: DP10M=110                !< No longer used/supported
    real, parameter :: ELWV=2.50E6              !< Latent heat of vaporization of water in J/kg
    real, parameter :: CPDR=CP/RD               !< CP over RD
    real, parameter :: RG=1./GI                 !< Reciprocal of GI (equal to G)
    real, parameter :: ELDR=ELWV/RD             !< No longer used/supported
    real, parameter :: BETA=.00367              !< Thermal expansion coefficient of air
    real, parameter :: BTG=BETA*G               !< Beta times acceleration due to gravity
    real, parameter :: FMW=18.015               !< Mean molecular mass of wet air
    real, parameter :: FMD=28.964               !< Mean molecular mass of dry air
    real, parameter :: EPS=FMW/FMD              !< Ratio of dry air to water vapor gas constants
    real, parameter :: ONEPS=1.-EPS             !< _____ Represents the complement of the ratio of molecular weights of water vapor to dry air ?
    real, parameter :: TFRZ=273.15              !< Freezing temperature of water
    real, parameter :: RDOCP=RD/CP              !< RD over CP
    real, parameter :: LHEAT=2.5E6              !< Latent heat
    real, parameter :: LHTOCP=LHEAT/CP          !< Latent heat over CP
    real, parameter :: EPSQ2=0.02               !< Minimum value of turbulent kinetic energy (m2 s^-2) ?
    real, parameter :: PQ0=379.90516            !< Constant used to parameterize specific humidity at 2m in WRFPOST: qs=pq0/p*exp(a2*(t-a3)/(t-a4))
    real, parameter :: A2=17.2693882            !< Constant used to parameterize specific humidity at 2m in WRFPOST: qs=pq0/p*exp(a2*(t-a3)/(t-a4))
    real, parameter :: A3=273.16                !< Constant used to parameterize specific humidity at 2m in WRFPOST: qs=pq0/p*exp(a2*(t-a3)/(t-a4))
    real, parameter :: A4=35.86                 !< Constant used to parameterize specific humidity at 2m in WRFPOST: qs=pq0/p*exp(a2*(t-a3)/(t-a4))
    real, parameter :: EPSQ=1.E-12              !< Minimum q (specific humidity) for computing precipitation type
#ifdef GSDRAPR
    real, parameter :: QCLDmin=1.E-6            !< Minimum cloud mixing ratio - was 1.E-6  
#else
    real, parameter :: QCLDmin=1.E-5            !< Minimum cloud mixing ratio - was 1.E-6  
#endif
    real, parameter :: CLFRmin=D1               !< Minimum cloud fraction
    real, parameter :: NLImin=1.E3              !< Minimum number concentrations (m**-3) of large ice (snow/graupel/sleet) 
! move definition of NLImax and T_ICE to MICROINIT 2012012018
!    real, parameter :: NLImax=5.E3
!    real, parameter :: T_ICE=-30.
!    real, parameter :: TRAD_ice=0.5*T_ICE+TFRZ
    real, parameter :: Qconv=0.1E-3       !< Heat transfer coefficient of convection _____?
!
    real, parameter :: STBOL=1./5.67E-8   !< 1 over Stefan-Boltzmann constant (W m^-2 deg^-4)
    real, parameter :: DBZmin=-20.        !< Minimum detectable radar reflectivity
    real, parameter :: abscoef=8.0        !< Cloud water absorption coefficient
    real, parameter :: abscoefi=5.0       !< Cloud ice absorption coefficient
    real, parameter :: XLAI=4.0           !< Leaf area index
    real, parameter :: RHmin=1.0E-6       !< minimal RH bound
    real, parameter :: SMALL=1.E-6        !< _____?
!    real, parameter :: PTHRESH=0.000000
    real, parameter :: ELIVW=2.72E6       !< Latent heat of vaporization of water in J/kg, used in calculations involving energy transfer during evaporation
    real, parameter :: ELOCP=ELIVW/CP     !< ELIVW over CP
!
    real, parameter :: GAMMA=6.5E-3       !< Gamma - constant for lapse rate in troposphere
    real, parameter :: RGAMOG=RD*GAMMA/G  !< Gas constant for dry air times gamma over acceleration due to gravity

end module params_mod
