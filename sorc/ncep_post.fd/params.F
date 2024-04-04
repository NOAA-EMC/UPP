!> @file
!> @brief module: params_mod defines physical constants and smoothing
!coefficients
module params_mod

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
! last update: 09/06/2006
!
!     SMOOTHING PARAMETERS.
    integer, parameter :: KSMUD=4      !< _____
    integer, parameter :: KSLPD=1      !< _____
    integer, parameter :: KSSLP=2      !< _____
    integer, parameter :: KSFI=2       !< _____
    integer, parameter :: KST=2        !< _____
    integer, parameter :: KSUV=2       !< _____
    integer, parameter :: KSRH=2       !< _____
    integer, parameter :: KSAV=2       !< _____
    integer, parameter :: KSLI=2       !< _____
    integer, parameter :: KSLP=2       !< _____
    integer, parameter :: NSUMD=1      !< _____
                                                                                   
    real, parameter :: SMTHA=-1.8E-4   !< _____
    real, parameter :: SMTHB=16.0      !< _____
!
! CONSTANTS.
    real, parameter :: CM1=2937.4               !< _____
    real, parameter :: CM2=4.9283               !< _____
    real, parameter :: CM3=23.5518              !< _____
    real, parameter :: D00=0.0                  !< _____
    real, parameter :: D1=0.1                   !< _____
    real, parameter :: D01=0.01                 !< _____
    real, parameter :: D001=0.001               !< _____
    real, parameter :: D50=0.50                 !< _____
    real, parameter :: H999=999.                !< _____
    real, parameter :: H1=1.0                   !< _____
    real, parameter :: H2=2.0                   !< _____
    real, parameter :: H4=4.0                   !< _____
    real, parameter :: H99999=99999.            !< _____
    real, parameter :: H1M12=1.E-12             !< _____
    real, parameter :: H1000=1000.              !< _____
    real, parameter :: H1M5=1.E-5               !< _____
    real, parameter :: D125=.125                !< _____
    real, parameter :: D25=0.25                 !< _____
    real, parameter :: H100=100.                !< _____
    real, parameter :: H10000=10000.            !< _____
    real, parameter :: H10E5=1.E5               !< _____
    real, parameter :: CAPA=0.28589641E0        !< _____
    real, parameter :: D115=0.115               !< _____
    real, parameter :: D608=0.608               !< _____
    real, parameter :: D05=0.05                 !< _____
    real, parameter :: D0065=0.0065             !< _____
    real, parameter :: H1E2=100.                !< _____
    real, parameter :: H1E4=1.E4                !< _____
    real, parameter :: H1M2=1.E-2               !< _____
    real, parameter :: G=9.81                   !< _____
    real, parameter :: RD=287.04                !< _____
    real, parameter :: ROG=RD/G                 !< _____
    real, parameter :: GI=1./G                  !< _____
    real, parameter :: OVERRC=1.10              !< _____
    real, parameter :: AD05=OVERRC*D05          !< _____
    real, parameter :: CFT0=OVERRC-H1           !< _____
    real, parameter :: PI=3.141592653589793     !< _____
    real, parameter :: DTR=3.1415926/180.       !< _____
    real, parameter :: RTD=1./DTR               !< _____
    real, parameter :: ERAD=6.371E6             !< _____
    real, parameter :: CP=1004.6                !< _____
    real, parameter :: RCAPA=1./CAPA            !< _____
    real, parameter :: P1000=1000.E2,DP10M=110  !< _____
    real, parameter :: ELWV=2.50E6              !< _____
    real, parameter :: CPDR=CP/RD               !< _____
    real, parameter :: RG=1./GI                 !< _____
    real, parameter :: ELDR=ELWV/RD             !< _____
    real, parameter :: BETA=.00367              !< _____
    real, parameter :: BTG=BETA*G               !< _____
    real, parameter :: FMW=18.015               !< _____
    real, parameter :: FMD=28.964               !< _____
    real, parameter :: EPS=FMW/FMD              !< _____
    real, parameter :: ONEPS=1.-EPS             !< _____
    real, parameter :: TFRZ=273.15              !< _____
    real, parameter :: RDOCP=RD/CP              !< _____
    real, parameter :: LHEAT=2.5E6              !< _____
    real, parameter :: LHTOCP=LHEAT/CP          !< _____
    real, parameter :: EPSQ2=0.02               !< _____
    real, parameter :: PQ0=379.90516            !< _____
    real, parameter :: A2=17.2693882            !< _____
    real, parameter :: A3=273.16                !< _____
    real, parameter :: A4=35.86                 !< _____
    real, parameter :: EPSQ=1.E-12              !< _____
#ifdef GSDRAPR
    real, parameter :: QCLDmin=1.E-6  ! was 1.E-6  !< _____
#else
    real, parameter :: QCLDmin=1.E-5  ! was 1.E-6  !< _____
#endif
    real, parameter :: CLFRmin=D1               !< _____
    real, parameter :: NLImin=1.E3              !< _____
! move definition of NLImax and T_ICE to MICROINIT 2012012018
!    real, parameter :: NLImax=5.E3
!    real, parameter :: T_ICE=-30.
!    real, parameter :: TRAD_ice=0.5*T_ICE+TFRZ
    real, parameter :: Qconv=0.1E-3             !< _____
!
    real, parameter :: STBOL=1./5.67E-8   !< _____
    real, parameter :: DBZmin=-20.        !< _____
    real, parameter :: abscoef=8.0        !< _____
    real, parameter :: abscoefi=5.0       !< _____
    real, parameter :: XLAI=4.0           !< Leaf area index
    real, parameter :: RHmin=1.0E-6       !< minimal RH bound
    real, parameter :: SMALL=1.E-6        !< _____
!    real, parameter :: PTHRESH=0.000000
    real, parameter :: ELIVW=2.72E6       !< _____
    real, parameter :: ELOCP=ELIVW/CP     !< _____
!
    real, parameter :: GAMMA=6.5E-3       !< _____
    real, parameter :: RGAMOG=RD*GAMMA/G  !< _____

end module params_mod
