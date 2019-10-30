module physcons_post
  use machine_post,only:kind_phys
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
!  Physical constants as set in NMC handbook from Smithsonian tables.
!  Physical constants are given to 5 places.
!  1990/04/30: g and rd are made consistent with NWS usage.
!  2001/10/22: g made consistent with SI usage.
!  Math constants
!  real(kind=kind_phys),parameter:: con_pi      =3.141593e+0 ! pi
  real(kind=kind_phys),parameter:: con_pi     =4.0d0*atan(1.0d0)
  real(kind=kind_phys),parameter:: con_sqrt2   =1.414214e+0_kind_phys ! square root of 2
  real(kind=kind_phys),parameter:: con_sqrt3   =1.732051e+0_kind_phys ! square root of 3
!  Primary constants
  real(kind=kind_phys),parameter:: con_rerth   =6.3712e+6_kind_phys ! radius of earth     (m)
  real(kind=kind_phys),parameter:: con_g       =9.80665e+0_kind_phys! gravity             (m/s2)
  real(kind=kind_phys),parameter:: con_omega   =7.2921e-5_kind_phys ! ang vel of earth    (1/s)
  real(kind=kind_phys),parameter:: con_rd      =2.8705e+2_kind_phys ! gas constant air    (J/kg/K)
  real(kind=kind_phys),parameter:: con_rv      =4.6150e+2_kind_phys ! gas constant H2O    (J/kg/K)
  real(kind=kind_phys),parameter:: con_cp      =1.0046e+3_kind_phys ! spec heat air @p    (J/kg/K)
  real(kind=kind_phys),parameter:: con_cv      =7.1760e+2_kind_phys ! spec heat air @v    (J/kg/K)
  real(kind=kind_phys),parameter:: con_cvap    =1.8460e+3_kind_phys ! spec heat H2O gas   (J/kg/K)
  real(kind=kind_phys),parameter:: con_cliq    =4.1855e+3_kind_phys ! spec heat H2O liq   (J/kg/K)
  real(kind=kind_phys),parameter:: con_csol    =2.1060e+3_kind_phys ! spec heat H2O ice   (J/kg/K)
  real(kind=kind_phys),parameter:: con_hvap    =2.5000e+6_kind_phys ! lat heat H2O cond   (J/kg)
  real(kind=kind_phys),parameter:: con_hfus    =3.3358e+5_kind_phys ! lat heat H2O fusion (J/kg)
  real(kind=kind_phys),parameter:: con_psat    =6.1078e+2_kind_phys ! pres at H2O 3pt     (Pa)  
  real(kind=kind_phys),parameter:: con_sbc     =5.6730e-8_kind_phys ! stefan-boltzmann    (W/m2/K4)
  real(kind=kind_phys),parameter:: con_solr    =1.3533e+3_kind_phys ! solar constant      (W/m2)
  real(kind=kind_phys),parameter:: con_t0c     =2.7315e+2_kind_phys ! temp at 0C          (K)
  real(kind=kind_phys),parameter:: con_ttp     =2.7316e+2_kind_phys ! temp at H2O 3pt     (K)
!  Secondary constants
  real(kind=kind_phys),parameter:: con_rocp    =con_rd/con_cp
  real(kind=kind_phys),parameter:: con_cpor    =con_cp/con_rd
  real(kind=kind_phys),parameter:: con_rog     =con_rd/con_g
  real(kind=kind_phys),parameter:: con_fvirt   =con_rv/con_rd-1.
  real(kind=kind_phys),parameter:: con_eps     =con_rd/con_rv
  real(kind=kind_phys),parameter:: con_epsm1   =con_rd/con_rv-1.
  real(kind=kind_phys),parameter:: con_dldt    =con_cvap-con_cliq
  real(kind=kind_phys),parameter:: con_xpona   =-con_dldt/con_rv
  real(kind=kind_phys),parameter:: con_xponb   =-con_dldt/con_rv+con_hvap/(con_rv*con_ttp)
end module
