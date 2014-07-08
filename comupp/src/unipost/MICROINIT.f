      SUBROUTINE MICROINIT(imp_physics)
!
!-- ABSTRACT:
!     Initializes arrays for new cloud microphysics
!
!-- Program History Log:
!     02-02-08  B. Ferrier
!     04-11-19 H CHUANG - WRF VERSION
!
!-- Input argument list:
!     None
!
!-- Output argument list:
!     None
!
!-- Subprograms called:
!     Function FPVS
!
!-- Common blocks:
!     CMASSI
!     RMASS_TABLES
!     MAPOT
!     CRHgrd
!
!-- Attributes:
!     Language: FORTRAN 90
!     Machine : IBM SP
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use params_mod, only: tfrz, pi
      use cmassi_mod, only: dmrmax, t_ice, NLImax1, NLImax2, flarge2, xmrmax, &
          mdrmax, mdrmin, trad_ice, massi, rqr_drmin, n0r0, rqr_drmax, &
          cn0r0, cn0r_dmrmin, cn0r_dmrmax, dmrmin, CLImax, &
          QLImax1,QLImax2
!      use gridspec_mod
      use rhgrd_mod, only: rhgrd

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      REAL, PARAMETER :: RHOL=1000.
      real ax,C_N0r0
      integer i
      integer, intent(in):: imp_physics
      real, allocatable:: MASSR(:)
!
!------------------------ START EXECUTION ------------------------
!
!---  READ IN MASSI FROM LOOKUP TABLES 
!
      if(imp_physics==5)then
       DMRmax=1.E-3
       T_ICE=-40.
       NLImax1=10.E3
       NLImax2=1.0E3
!-- Note that QLImax1 must be < QLImax2 or else the results will be bad
       QLImax1=0.9E-3
       QLImax2=2.5E-3
       FLARGE2=0.07     !-- Set but no longer used
!
!-- Adjusts NLImax from a value of NLImax1 for ice contents <= QLImax1
!   to NLImax2 for ice contents >= QLImax2
!
       IF (QLImax1 >= QLImax2) THEN
          WRITE(0,*) 'QLImax1 must be < QLImax2 ... BAD RESULTS FOLLOW!!'
       ENDIF
       CLImax=ALOG(NLImax2/NLImax1)/(QLImax2-QLImax1)
       WRITE(0,*) 'CLImax=',CLImax
      else if(imp_physics==85)then
       DMRmax=.45E-3
       T_ICE=-40.
       NLImax1=20.E3
       NLImax2=20.E3
       CLImax=0.
       FLARGE2=0.2
      else  !-- Should be imp_physics==95
       DMRmax=.45E-3
       T_ICE=-40.  
       NLImax1=5.E3
       NLImax2=5.E3
       CLImax=0.
       FLARGE2=0.03
      end if 
      XMRmax=1.E6*DMRmax 
      MDRmax=XMRmax
      allocate(MASSR(MDRmin:MDRmax))
      TRAD_ice=0.5*T_ICE+TFRZ
      
      OPEN (UNIT=1,FILE="eta_micro_lookup.dat",convert='big_endian',FORM="UNFORMATTED")
      DO I=1,3
        READ(1)
      ENDDO
      READ(1) MASSR
      DO I=1,5
        READ(1)
      ENDDO
      READ(1) MASSI
      CLOSE(1)
      RQR_DRmin=N0r0*MASSR(MDRmin)    ! Rain content for mean drop diameter of .05 mm
      RQR_DRmax=N0r0*MASSR(MDRmax)    ! Rain content for mean drop diameter of .45 mm
!      PI=ACOS(-1.) ! defined in params now
      C_N0r0=PI*RHOL*N0r0
      CN0r0=1.E6/C_N0r0**.25
      CN0r_DMRmin=1./(PI*RHOL*DMRmin**4)
      CN0r_DMRmax=1./(PI*RHOL*DMRmax**4)
      print *,'MICROINIT: MDRmin, MASSR(MDRmin)=',MDRmin,MASSR(MDRmin)
      print *,'MICROINIT: MDRmax, MASSR(MDRmax)=',MDRmax,MASSR(MDRmax)
!      print *,  'ETA2P:MASSI(50)= ', MASSI(50)
!      print *,  'ETA2P:MASSI(450)= ', MASSI(450)
!      print *,  'ETA2P:MASSI(1000)= ', MASSI(1000)
!
!--- Initialize saturation vapor pressure lookup tables (functions FPVS, FPVS0)
!
      CALL GPVS
!
!--- Initialize RHgrd, grid-scale RH for onset of condensation. 
!    See GSMCONST in Eta model for algorithm with grid-size dependence.
!
!      AX=111.*(DPHD**2+DLMD**2)**.5
!      AX=111.*(DYVAL/1000.**2+DXVAL/1000.**2)**.5
!      AX=MIN(100., MAX(5., AX) )
!      RHgrd=0.90+.08*((100.-AX)/95.)**.5
      RHgrd=1.
      deallocate(MASSR)
!--- 
      RETURN
      END
