!> @file
!>
!> @brief upp_physics is a collection of UPP subroutines for physics variables calculation.
!> @author Jesse Meng @date 2020-05-20
!>
!> calcape() computes CAPE/CINS and other storm related variables.
!>
!> calcape2() computes additional storm related variables.
!>
!> caldiv() computes divergence.
!>
!> calgradps() computes gardients of a scalar field PS or LNPS. 
!>
!> calrh(), calrh_nam(), calrh_gfs(), calrh_gsd() compute RH using various algorithms.
!>
!> The NAM v4.1.18 algorithm (calrh_nam()) is selected as default for 
!> NMMB and FV3GFS, FV3GEFS, and FV3R for the UPP 2020 unification.
!>
!> calrh_pw() algorithm use at GSD for RUC and Rapid Refresh.
!>
!> calslr_roebber() computes snow solid-liquid-ratio slr using the Roebber algorithm.
!>      
!> calslr_uutah() computes snow solid-liquid-ratio slr using the UUtah Steenburgh algorithm.
!>   
!> calvor() computes absolute vorticity.   
!>      
!> fpvsnew() computes saturation vapor pressure.
!>
!> tvirtual() computes virtual temperature.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2020-05-20 | Jesse Meng | Initial
!> 2022-07-11 | Jesse Meng | CALSLR_ROEBBER
!> 2023-02-14 | Jesse Meng | CALSLR_UUTAH     
!> 2023-03-22 | Sam Trahan | Fix out-of-bounds access by not calling BOUND
!>
!> @author Jesse Meng @date 2020-05-20
  module upp_physics

  implicit none

  private

  public :: CALCAPE, CALCAPE2
  public :: CALDIV
  public :: CALGRADPS

  public :: CALRH
  public :: CALRH_GFS, CALRH_GSD, CALRH_NAM
  public :: CALRH_PW
  public :: CALSLR_ROEBBER, CALSLR_UUTAH
  public :: CALVOR

  public :: FPVSNEW
  public :: TVIRTUAL

  contains
!
!-------------------------------------------------------------------------------------
!> CALRH() computes relative humidity
!>
!> @param[in] P1 real Pressure (Pa)
!> @param[in] T1 real Temperature (K)
!> @param[inout] Q1 real Specific humidity (kg/kg)
!> @note P1/T1/Q1 refer to pressure/temperature/specific humidity at the selected
!> /requested level. For example, if the user wants relative humidity at 500mb, 
!> then T/P/Q at 500mb is input into the call to output RH@500mb
!> @param[out] RH real Relative humidity (decimal form)
!
      SUBROUTINE CALRH(P1,T1,Q1,RH)

      use ctlblk_mod, only: ista, iend, jsta, jend, MODELNAME
      implicit none

      REAL,dimension(ista:iend,jsta:jend),intent(in)    :: P1,T1
      REAL,dimension(ista:iend,jsta:jend),intent(inout) :: Q1
      REAL,dimension(ista:iend,jsta:jend),intent(out)   :: RH

      IF(MODELNAME == 'RAPR')THEN
         CALL CALRH_GSD(P1,T1,Q1,RH)
      ELSE
         CALL CALRH_NAM(P1,T1,Q1,RH)
      END IF

      END SUBROUTINE CALRH
!
!-------------------------------------------------------------------------------------
!
!> calrh_nam() computes relative humidity.
!>
!> This routine computes relative humidity given pressure, 
!> temperature, specific humidity. an upper and lower bound
!> of 100 and 1 percent relative humidity is enforced.  When
!> these bounds are applied the passed specific humidity 
!> array is adjusted as necessary to produce the set relative
!> humidity.
!>
!> @param[in] P1 Pressure (pa)
!> @param[in] T1 Temperature (K)
!> @param[inout] Q1 Specific humidity (kg/kg)
!> @param[out] RH Relative humidity (decimal form)
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> ????-??-?? | DENNIS DEAVEN | Initial
!> 1992-12-22 | Russ Treadon  | Modified as described above
!> 1998-06-08 | T Black       | Conversion from 1-D to 2-D
!> 1998-08-18 | Mike Baldwin  | Modify to compute RH over ice as in model
!> 1998-12-16 | Geoff Manikin | undo RH computation over ice
!> 2000-01-04 | Jim Tuccillo  | MPI Version
!> 2002-06-11 | Mike Baldwin  | WRF Version
!> 2006-03-19 | Wen Meng      | Modify top pressure to 1 pa
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
     SUBROUTINE CALRH_NAM(P1,T1,Q1,RH)
     use params_mod, only: PQ0, a2, a3, a4, rhmin
     use ctlblk_mod, only: ista, iend, jsta, jend, spval
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     SET PARAMETER.
!
!     DECLARE VARIABLES.
!     
      REAL,dimension(ista:iend,jsta:jend),intent(in)    :: P1,T1
      REAL,dimension(ista:iend,jsta:jend),intent(inout) :: Q1
      REAL,dimension(ista:iend,jsta:jend),intent(out)   :: RH
      REAL QC
      integer I,J
!***************************************************************
!
!     START CALRH.
!
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IF (T1(I,J) < spval) THEN
            IF (ABS(P1(I,J)) >= 1) THEN
              QC = PQ0/P1(I,J)*EXP(A2*(T1(I,J)-A3)/(T1(I,J)-A4))
!
              RH(I,J) = Q1(I,J)/QC
!
!   BOUNDS CHECK
!
              IF (RH(I,J) > 1.0) THEN
                RH(I,J) = 1.0
                Q1(I,J) = RH(I,J)*QC
              ENDIF
              IF (RH(I,J) < RHmin) THEN  !use smaller RH limit for stratosphere
                RH(I,J) = RHmin
                Q1(I,J) = RH(I,J)*QC
              ENDIF
!
            ENDIF
          ELSE
            RH(I,J) = spval
          ENDIF
        ENDDO
      ENDDO
!
!
!      END SUBROUTINE CALRH
      END SUBROUTINE CALRH_NAM
!
!-------------------------------------------------------------------------------------
!
!> calrh_gfs() computes relative humidity.
!>
!> This routine computes relative humidity given pressure, 
!> temperature, specific humidity. an upper and lower bound
!> of 100 and 1 percent relative humidity is enforced.  When
!> these bounds are applied the passed specific humidity 
!> array is adjusted as necessary to produce the set relative
!> humidity.
!>
!> @param[in] P1 Pressure (pa)
!> @param[in] T1 Temperature (K)
!> @param[inout] Q1 Specific humidity (kg/kg)
!> @param[out] RH Relative humidity (decimal form)
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> ????-??-?? | DENNIS DEAVEN | Initial
!> 1992-12-22 | Russ Treadon  | Modified as described above
!> 1998-06-08 | T Black       | Conversion from 1-D to 2-D
!> 1998-08-18 | Mike Baldwin  | Modify to compute RH over ice as in model
!> 1998-12-16 | Geoff Manikin | undo RH computation over ice
!> 2000-01-04 | Jim Tuccillo  | MPI Version
!> 2002-06-11 | Mike Baldwin  | WRF Version
!> 2013-08-13 | S. Moorthi    | Threading
!> 2006-03-19 | Wen Meng      | Modify top pressure to 1 pa
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE CALRH_GFS(P1,T1,Q1,RH)
      use params_mod, only: rhmin
      use ctlblk_mod, only: ista, iend, jsta, jend, spval
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      real,parameter:: con_rd      =2.8705e+2 ! gas constant air    (J/kg/K)
      real,parameter:: con_rv      =4.6150e+2 ! gas constant H2O 
      real,parameter:: con_eps     =con_rd/con_rv
      real,parameter:: con_epsm1   =con_rd/con_rv-1
!      real,external::FPVSNEW

!      INTERFACE
!        ELEMENTAL FUNCTION FPVSNEW (t)
!          REAL  FPVSNEW
!          REAL, INTENT(IN) :: t
!        END FUNCTION FPVSNEW
!      END INTERFACE
!
      REAL,dimension(ista:iend,jsta:jend),intent(in)   :: P1,T1
      REAL,dimension(ista:iend,jsta:jend),intent(inout):: Q1,RH
      REAL ES,QC
      integer :: I,J
!***************************************************************
!
!     START CALRH.
!
!$omp parallel do private(i,j,es,qc)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IF (T1(I,J) < spval .AND. P1(I,J) < spval.AND.Q1(I,J)/=spval) THEN
!           IF (ABS(P1(I,J)) > 1.0) THEN
!            IF (P1(I,J) > 1.0) THEN
            IF (P1(I,J) >= 1.0) THEN
              ES = MIN(FPVSNEW(T1(I,J)),P1(I,J))
              QC = CON_EPS*ES/(P1(I,J)+CON_EPSM1*ES)

!             QC=PQ0/P1(I,J)*EXP(A2*(T1(I,J)-A3)/(T1(I,J)-A4))

              RH(I,J) = min(1.0,max(Q1(I,J)/QC,rhmin))
              q1(i,j) = rh(i,j)*qc

!   BOUNDS CHECK
!
!             IF (RH(I,J) > 1.0) THEN
!               RH(I,J) = 1.0
!               Q1(I,J) = RH(I,J)*QC
!             ELSEIF (RH(I,J) < RHmin) THEN  !use smaller RH limit for stratosphere
!               RH(I,J) = RHmin
!               Q1(I,J) = RH(I,J)*QC
!             ENDIF

            ENDIF
          ELSE
            RH(I,J) = spval
          ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE CALRH_GFS
!
!-------------------------------------------------------------------------------------
!> CALRH_GSD() Compute RH with the NOAA GSL (formerly NOAA GSD) algorithm used for RUC and Rapid Refresh
!>
!> @param P1 real Pressure (Pa)
!> @param T1 real Temperature (K)
!> @param Q1 real Specific humidity (kg/kg)
!> @note P1/T1/Q1 refer to pressure/temperature/specific humidity at the selected
!> /requested level. For example, if the user wants relative humidity at 500mb, 
!> then T/P/Q at 500mb is input into the call to output RH@500mb
!> @param RHB real Relative humidity (decimal form)
!> 
!      
      SUBROUTINE CALRH_GSD(P1,T1,Q1,RHB)
!

      use ctlblk_mod, only: ista, iend, jsta, jend, spval

      implicit none

      integer :: j, i
      real :: tx, pol, esx, es, e
      real, dimension(ista:iend,jsta:jend) :: P1, T1, Q1, RHB


      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IF (T1(I,J) < spval .AND. P1(I,J) < spval .AND. Q1(I,J) < spval) THEN
! - compute relative humidity
          Tx=T1(I,J)-273.15
          POL = 0.99999683       + TX*(-0.90826951E-02 +    &
             TX*(0.78736169E-04   + TX*(-0.61117958E-06 +   &
             TX*(0.43884187E-08   + TX*(-0.29883885E-10 +   &
             TX*(0.21874425E-12   + TX*(-0.17892321E-14 +   &
             TX*(0.11112018E-16   + TX*(-0.30994571E-19)))))))))
          esx = 6.1078/POL**8

          ES = esx
          E = P1(I,J)/100.*Q1(I,J)/(0.62197+Q1(I,J)*0.37803)
          RHB(I,J) = MIN(1.,E/ES)
         ELSE
          RHB(I,J) = spval
         ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE CALRH_GSD
!
!-------------------------------------------------------------------------------------
!
!> CALRH_PW() algorithm used at GSL for RUC and Rapid Refresh.
!>
!> @param RHPW real Relative humidity with respect to precipitable water (entire atmosphere)
!> 

      SUBROUTINE CALRH_PW(RHPW)
!
!------------------------------------------------------------------
!

      use vrbls3d, only: q, pmid, t
      use params_mod, only: g
      use ctlblk_mod, only: lm, ista, iend, jsta, jend, spval
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none

      real,PARAMETER :: svp1=6.1153,svp2=17.67,svp3=29.65

      REAL, dimension(ista:iend,jsta:jend):: PW, PW_SAT, RHPW
      REAL deltp,sh,qv,temp,es,qs,qv_sat
      integer i,j,l,k,ka,kb

      pw     = 0.
      pw_sat = 0.
      rhpw   = 0.

      DO L=1,LM
        k=lm-l+1
       DO J=JSTA,JEND
        DO I=ISTA,IEND
! -- use specific humidity for PW calculation
         if(t(i,j,k)<spval.and.q(i,j,k)<spval) then
           sh = q(i,j,k)
           qv = sh/(1.-sh)
           KA = MAX(1,K-1)
           KB = MIN(LM,K+1)

!   assumes that P is in mb at this point - be careful!
           DELTP = 0.5*(PMID(I,J,KB)-PMID(I,J,KA))
           PW(I,J) = PW(I,J) + sh *DELTP/G

!Csgb -- Add more for RH w.r.t. PW-sat

          temp = T(I,J,K)
! --- use saturation mixing ratio w.r.t. water here
!       for this check.
          es = svp1*exp(SVP2*(Temp-273.15)/(Temp-SVP3))
! -- get saturation specific humidity (w.r.t. total air)
          qs = 0.62198*es/(pmid(i,j,k)*1.e-2-0.37802*es)
! -- get saturation mixing ratio (w.r.t. dry air)
          qv_sat = qs/(1.-qs)

          pw_sat(i,j) = pw_sat(i,j) + max(sh,Qs)*DELTP/G

!        if (i==120 .and. j==120 )                        &
!          write (6,*)'pw-sat', temp, sh, qs, pmid(i,j,kb)    &
!          ,pmid(i,j,ka),pw(i,j),pw_sat(i,j)

!sgb - This IS RH w.r.t. PW-sat.
           RHPW (i,j) = min(1.,PW(i,j) / pw_sat(i,j)) * 100.
          else
           RHPW (i,j) = spval
          endif
        ENDDO
       ENDDO
      ENDDO

      END SUBROUTINE CALRH_PW
!
!-------------------------------------------------------------------------------------
!
      elemental function fpvsnew(t)
!> fpvsnew() computes saturation vapor pressure.
!>
!> Compute saturation vapor pressure from the temperature.
!> A linear interpolation is done between values in a lookup table
!> computed in gpvs. See documentation for fpvsx for details.
!> Input values outside table range are reset to table extrema.
!> The interpolation accuracy is almost 6 decimal places.
!> On the Cray, fpvs is about 4 times faster than exact calculation.
!> This function should be expanded inline in the calling routine.
!>
!> @param[in] t Real(krealfp) Temperature in Kelvin.
!> @return fpvsnew Real(krealfp) Saturation vapor pressure in Pascals.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1991-05-07 | Iredell | Initial. Made into inlinable function
!> 1994-12-30 | Iredell | Expand table
!> 1999-03-01 | Iredell | F90 module
!> 2001-02-26 | Iredell | Ice phase
!>
!> @author N Phillips w/NMC2X2 @date 1982-12-30
      implicit none
      integer,parameter:: nxpvs=7501
      real,parameter:: con_ttp     =2.7316e+2 ! temp at H2O 3pt
      real,parameter:: con_psat    =6.1078e+2 ! pres at H2O 3pt
      real,parameter:: con_cvap    =1.8460e+3 ! spec heat H2O gas   (J/kg/K)
      real,parameter:: con_cliq    =4.1855e+3 ! spec heat H2O liq
      real,parameter:: con_hvap    =2.5000e+6 ! lat heat H2O cond
      real,parameter:: con_rv      =4.6150e+2 ! gas constant H2O
      real,parameter:: con_csol    =2.1060e+3 ! spec heat H2O ice
      real,parameter:: con_hfus    =3.3358e+5 ! lat heat H2O fusion
      real,parameter:: tliq=con_ttp
      real,parameter:: tice=con_ttp-20.0
      real,parameter:: dldtl=con_cvap-con_cliq
      real,parameter:: heatl=con_hvap
      real,parameter:: xponal=-dldtl/con_rv
      real,parameter:: xponbl=-dldtl/con_rv+heatl/(con_rv*con_ttp)
      real,parameter:: dldti=con_cvap-con_csol
      real,parameter:: heati=con_hvap+con_hfus
      real,parameter:: xponai=-dldti/con_rv
      real,parameter:: xponbi=-dldti/con_rv+heati/(con_rv*con_ttp)
      real tr,w,pvl,pvi
      real fpvsnew
      real,intent(in):: t
      integer jx
      real  xj,x,tbpvs(nxpvs),xp1
      real xmin,xmax,xinc,c2xpvs,c1xpvs
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      xmin=180.0
      xmax=330.0
      xinc=(xmax-xmin)/(nxpvs-1)
!   c1xpvs=1.-xmin/xinc
      c2xpvs=1./xinc
      c1xpvs=1.-xmin*c2xpvs
!    xj=min(max(c1xpvs+c2xpvs*t,1.0),real(nxpvs,krealfp))
      xj=min(max(c1xpvs+c2xpvs*t,1.0),float(nxpvs))
      jx=min(xj,float(nxpvs)-1.0)
      x=xmin+(jx-1)*xinc
      
      tr=con_ttp/x
      if(x>=tliq) then
        tbpvs(jx)=con_psat*(tr**xponal)*exp(xponbl*(1.-tr))
      elseif(x<tice) then
        tbpvs(jx)=con_psat*(tr**xponai)*exp(xponbi*(1.-tr))
      else
        w=(t-tice)/(tliq-tice)
        pvl=con_psat*(tr**xponal)*exp(xponbl*(1.-tr))
        pvi=con_psat*(tr**xponai)*exp(xponbi*(1.-tr))
        tbpvs(jx)=w*pvl+(1.-w)*pvi
      endif
      
      xp1=xmin+(jx-1+1)*xinc      
     
      tr=con_ttp/xp1
      if(xp1>=tliq) then
        tbpvs(jx+1)=con_psat*(tr**xponal)*exp(xponbl*(1.-tr))
      elseif(xp1<tice) then
        tbpvs(jx+1)=con_psat*(tr**xponai)*exp(xponbi*(1.-tr))
      else
        w=(t-tice)/(tliq-tice)
        pvl=con_psat*(tr**xponal)*exp(xponbl*(1.-tr))
        pvi=con_psat*(tr**xponai)*exp(xponbi*(1.-tr))
        tbpvs(jx+1)=w*pvl+(1.-w)*pvi
      endif
      
      fpvsnew=tbpvs(jx)+(xj-jx)*(tbpvs(jx+1)-tbpvs(jx))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end function fpvsnew
!
!-------------------------------------------------------------------------------------
!> calcape() computes CAPE and CINS.
!>
!> This routine computes CAPE and CINS given temperature,
!> pressure, and specific humidty.  In "storm and cloud 
!> dynamics" (1989, academic press) cotton and anthes define
!> CAPE (equation 9.16, p501) as
!>
!> @code
!>                  EL
!>        CAPE =  SUM G * LN(THETAP/THETAA) DZ 
!>                 LCL
!>     
!>     Where,
!>      EL    = Equilibrium level,
!>     LCL    = Lifting condenstation level,
!>       G    = Gravitational acceleration,
!>     THETAP = Lifted parcel potential temperature,
!>     THETAA = Ambient potential temperature.
!> @endcode
!>     
!>     Note that the integrand ln(THETAP/THETAA) approximately
!>     equals (THETAP-THETAA)/THETAA.  This ratio is often used
!>     in the definition of CAPE/CINS.
!>     
!>     Two types of CAPE/CINS can be computed by this routine.  The
!>     summation process is the same For both cases.  What differs
!>     is the definition of the parcel to lift.  FOR ITYPE=1 the
!>     parcel with the warmest THETA-E in A DPBND pascal layer above
!>     the model surface is lifted.  the arrays P1D, T1D, and Q1D
!>     are not used.  For itype=2 the arrays P1D, T1D, and Q1D
!>     define the parcel to lift in each column.  Both types of
!>     CAPE/CINS may be computed in a single execution of the post
!>     processor.
!>     
!>     This algorithm proceeds as follows.
!>     For each column, 
!>        (1)  Initialize running CAPE and CINS SUM TO 0.0
!>        (2)  Compute temperature and pressure at the LCL using
!>             look up table (PTBL).  Use either parcel that gives
!>             max THETAE in lowest DPBND above ground (ITYPE=1)
!>             or given parcel from t1D,Q1D,...(ITYPE=2).
!>        (3)  Compute the temp of a parcel lifted from the LCL.
!>             We know that the parcel's
!>             equivalent potential temperature (THESP) remains
!>             constant through this process.  we can
!>             compute tpar using this knowledge using look
!>             up table (subroutine TTBLEX).
!>        (4)  Find the equilibrium level.  This is defined as the
!>             highest positively buoyant layer.
!>             (If there is no positively buoyant layer, CAPE/CINS
!>              will be zero)
!>        (5)  Compute CAPE/CINS.  
!>             (A) Compute THETAP.  We know TPAR and P.
!>             (B) Compute THETAA.  We know T and P.  
!>        (6)  Add G*(THETAP-THETAA)*DZ to the running CAPE or CINS sum.
!>             (A) If THETAP > THETAA, add to the CAPE sum.
!>             (B) If THETAP < THETAA, add to the CINS sum.
!>        (7)  Are we at equilibrium level? 
!>             (A) If yes, stop the summation.
!>             (b) if no, contiunue the summation.
!>        (8)  Enforce limits on CAPE and CINS (i.e. no negative CAPE)
!>
!> @param[in] ITYPE INTEGER Flag specifying how parcel to lift is identified.  See comments above.
!> @param[in] DPBND Depth over which one searches for most unstable parcel.
!> @param[in] P1D Array of pressure of parcels to lift.
!> @param[in] T1D Array of temperature of parcels to lift.
!> @param[in] Q1D Array of specific humidity of parcels to lift.
!> @param[in] L1D Array of model level of parcels to lift.
!> @param[out] CAPE Convective available potential energy (J/kg).
!> @param[out] CINS Convective inhibition (J/kg).
!> @param[out] PPARC Pressure level of parcel lifted when one searches over a particular depth to compute CAPE/CIN.
!> @param[inout] ZEQL Equivalent level height.
!> @param THUND Thunder parameter.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1993-02-10 | Russ Treadon  | Initial
!> 1993-06-19 | Russ Treadon  | Generalized routine to allow for type 2 CAPE/CINS calculations
!> 1994-09-23 | Mike Baldwin  | Modified to use look up tables instead of complicated equations
!> 1994-10-13 | Mike Baldwin  | Modified to continue CAPE/CINS calc up to at highest buoyant layer
!> 1998-06-12 | T Black       | Conversion from 1-D TO 2-D
!> 1998-08-18 | T Black       | Compute APE internally
!> 2000-01-04 | Jim Tuccillo  | MPI Version              
!> 2002-01-15 | Mike Baldwin  | WRF Version
!> 2003-08-24 | G Manikin     | Added level of parcel being lifted as output from the routine and added the depth over which one searches for the most unstable parcel as input
!> 2010-09-09 | G Manikin     | Changed computation to use virtual temp added eq lvl hght and thunder parameter    
!> 2015-??-?? | S Moorthi     | Optimization and threading
!> 2021-07-28 | W Meng        | Restrict computation from undefined grids
!> 2021-09-01 | E Colon       | Equivalent level height index for RTMA
!>
!> @author Russ Treadon W/NP2 @date 1993-02-10
      SUBROUTINE CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,L1D,CAPE,    &  
                         CINS,PPARC,ZEQL,THUND)
      use vrbls3d,    only: pmid, t, q, zint
      use vrbls2d,    only: teql,ieql
      use masks,      only: lmh
      use params_mod, only: d00, h1m12, h99999, h10e5, capa, elocp, eps,  &
                            oneps, g
      use lookup_mod, only: thl, rdth, jtb, qs0, sqs, rdq, itb, ptbl,     &
                            plq, ttbl, pl, rdp, the0, sthe, rdthe, ttblq, &
                            itbq, jtbq, rdpq, the0q, stheq, rdtheq
      use ctlblk_mod, only: jsta_2l, jend_2u, lm, jsta, jend, im, me, spval, &
                            ista_2l, iend_2u, ista, iend
!     
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     INCLUDE/SET PARAMETERS.  CONSTANTS ARE FROM BOLTON (MWR, 1980).
      real,PARAMETER :: ISMTHP=2,ISMTHT=2,ISMTHQ=2
!     
!     DECLARE VARIABLES.
!
      integer,intent(in) :: ITYPE
      real,intent(in)    :: DPBND
      integer, dimension(ista:iend,Jsta:jend),intent(in)    :: L1D
      real,    dimension(ista:iend,Jsta:jend),intent(in)    :: P1D,T1D
      real,    dimension(ista:iend,jsta:jend),intent(inout) :: Q1D,CAPE,CINS,PPARC,ZEQL
!     
      integer, dimension(ista:iend,jsta:jend) :: IPTB, ITHTB, PARCEL, KLRES, KHRES, LCL, IDX
!     
      real,    dimension(ista:iend,jsta:jend) :: THESP, PSP, CAPE20, QQ, PP, THUND  
      REAL, ALLOCATABLE :: TPAR(:,:,:)

      LOGICAL THUNDER(ista:iend,jsta:jend), NEEDTHUN 
      real PSFCK,PKL,TBTK,QBTK,APEBTK,TTHBTK,TTHK,APESPK,TPSPK,        &
           BQS00K,SQS00K,BQS10K,SQS10K,BQK,SQK,TQK,PRESK,GDZKL,THETAP, &
           THETAA,P00K,P10K,P01K,P11K,TTHESK,ESATP,QSATP,TVP,TV
!      real,external :: fpvsnew
      integer I,J,L,KNUML,KNUMH,LBEG,LEND,IQ, KB,ITTBK

!     integer I,J,L,KNUML,KNUMH,LBEG,LEND,IQ,IT,LMHK, KB,ITTBK
!     
!**************************************************************
!     START CALCAPE HERE.
!     
      ALLOCATE(TPAR(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U,LM))
!
!     COMPUTE CAPE/CINS
!
!        WHICH IS: THE SUM FROM THE LCL TO THE EQ LEVEL OF
!             G * (LN(THETAP) - LN(THETAA)) * DZ
!
!             (POSITIVE AREA FOR CAPE, NEGATIVE FOR CINS)
!
!        WHERE:
!             THETAP IS THE PARCEL THETA
!             THETAA IS THE AMBIENT THETA
!             DZ IS THE THICKNESS OF THE LAYER
!
!         USING LCL AS LEVEL DIRECTLY BELOW SATURATION POINT
!         AND EQ LEVEL IS THE HIGHEST POSITIVELY BUOYANT LEVEL.
!  
!         IEQL = EQ LEVEL
!         P_thetaemax - real  pressure of theta-e max parcel (Pa)
!
!     INITIALIZE CAPE AND CINS ARRAYS
! 
!$omp  parallel do
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          CAPE(I,J)    = D00
          CAPE20(I,J)  = D00
          CINS(I,J)    = D00
          LCL(I,J)     = 0
          THESP(I,J)   = D00
          IEQL(I,J)    = LM
          PARCEL(I,J)  = LM
          PSP(I,J)     = D00
          PPARC(I,J)   = D00
          THUNDER(I,J) = .TRUE.
        ENDDO
      ENDDO
!
!$omp  parallel do
      DO L=1,LM
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            TPAR(I,J,L) = D00
          ENDDO
        ENDDO
      ENDDO
!     
!     TYPE 2 CAPE/CINS:
!     NOTE THAT FOR TYPE 1 CAPE/CINS ARRAYS P1D, T1D, Q1D 
!     ARE DUMMY ARRAYS.
!     
      IF (ITYPE == 2) THEN
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            Q1D(I,J) = MIN(MAX(H1M12,Q1D(I,J)),H99999)
          ENDDO
        ENDDO
      ENDIF
!-------FOR ITYPE=1--FIND MAXIMUM THETA E LAYER IN LOWEST DPBND ABOVE GROUND-------
!-------FOR ITYPE=2--FIND THETA E LAYER OF GIVEN T1D, Q1D, P1D---------------------
!--------------TRIAL MAXIMUM BUOYANCY LEVEL VARIABLES-------------------

      DO KB=1,LM
!hc     IF (ITYPE==2.AND.KB>1) cycle
        IF (ITYPE == 1 .OR. (ITYPE == 2 .AND. KB == 1)) THEN

!$omp  parallel do private(i,j,apebtk,apespk,bqk,bqs00k,bqs10k,iq,ittbk,    &
!$omp &         p00k,p01k,p10k,p11k,pkl,psfck,qbtk,sqk,sqs00k,              &
!$omp &         sqs10k,tbtk,tpspk,tqk,tthbtk,tthesk,tthk)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              PSFCK  = PMID(I,J,NINT(LMH(I,J)))
              PKL    = PMID(I,J,KB)
              IF(PSFCK<spval.and.PKL<spval)THEN

!hc           IF (ITYPE==1.AND.(PKL<PSFCK-DPBND.OR.PKL>PSFCK)) cycle
              IF (ITYPE ==2 .OR.                                                &
                 (ITYPE == 1 .AND. (PKL >= PSFCK-DPBND .AND. PKL <= PSFCK)))THEN
                IF (ITYPE == 1) THEN
                  TBTK   = T(I,J,KB)
                  QBTK   = max(0.0, Q(I,J,KB))
                  APEBTK = (H10E5/PKL)**CAPA
                ELSE
                  PKL    = P1D(I,J)
                  TBTK   = T1D(I,J)
                  QBTK   = max(0.0, Q1D(I,J))
                  APEBTK = (H10E5/PKL)**CAPA
                ENDIF

!----------Breogan Gomez - 2009-02-06
! To prevent QBTK to be less than 0 which leads to a unrealistic value of PRESK
!  and a floating invalid.

!               if(QBTK < 0) QBTK = 0

!--------------SCALING POTENTIAL TEMPERATURE & TABLE INDEX--------------
                TTHBTK  =  TBTK*APEBTK
                TTHK    = (TTHBTK-THL)*RDTH
                QQ(I,J) = TTHK - AINT(TTHK)
                ITTBK   = INT(TTHK) + 1
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
                IF(ITTBK < 1)   THEN
                  ITTBK   = 1
                  QQ(I,J) = D00
                ENDIF
                IF(ITTBK >= JTB) THEN
                  ITTBK   = JTB-1
                  QQ(I,J) = D00
                ENDIF
!--------------BASE AND SCALING FACTOR FOR SPEC. HUMIDITY---------------
                BQS00K = QS0(ITTBK)
                SQS00K = SQS(ITTBK)
                BQS10K = QS0(ITTBK+1)
                SQS10K = SQS(ITTBK+1)
!--------------SCALING SPEC. HUMIDITY & TABLE INDEX---------------------
                BQK     = (BQS10K-BQS00K)*QQ(I,J) + BQS00K
                SQK     = (SQS10K-SQS00K)*QQ(I,J) + SQS00K
                TQK     = (QBTK-BQK)/SQK*RDQ
                PP(I,J) = TQK-AINT(TQK)
                IQ      = INT(TQK)+1
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
                IF(IQ < 1)    THEN
                  IQ      = 1
                  PP(I,J) = D00
                ENDIF
                IF(IQ >= ITB)  THEN
                  IQ      = ITB-1
                  PP(I,J) = D00
                ENDIF
!--------------SATURATION PRESSURE AT FOUR SURROUNDING TABLE PTS.-------
                P00K = PTBL(IQ  ,ITTBK  )
                P10K = PTBL(IQ+1,ITTBK  )
                P01K = PTBL(IQ  ,ITTBK+1)
                P11K = PTBL(IQ+1,ITTBK+1)
!--------------SATURATION POINT VARIABLES AT THE BOTTOM-----------------
                TPSPK = P00K + (P10K-P00K)*PP(I,J) + (P01K-P00K)*QQ(I,J)  &
                             + (P00K-P10K-P01K+P11K)*PP(I,J)*QQ(I,J)

!!from WPP::tgs          APESPK=(H10E5/TPSPK)**CAPA
                if (TPSPK > 1.0e-3) then
                  APESPK = (max(0.,H10E5/ TPSPK))**CAPA
                else
                  APESPK = 0.0
                endif

                TTHESK = TTHBTK * EXP(ELOCP*QBTK*APESPK/TTHBTK)
!--------------CHECK FOR MAXIMUM THETA E--------------------------------
                IF(TTHESK > THESP(I,J)) THEN
                  PSP  (I,J)  = TPSPK
                  THESP(I,J)  = TTHESK
                  PARCEL(I,J) = KB
                ENDIF
              END IF 
              ENDIF !end PSFCK<spval.and.PKL<spval
            ENDDO  ! I  loop
          ENDDO    ! J  loop
        END IF
      ENDDO        ! KB loop

!----FIND THE PRESSURE OF THE PARCEL THAT WAS LIFTED
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            PPARC(I,J) = PMID(I,J,PARCEL(I,J))
          ENDDO
        ENDDO
!
!-----CHOOSE LAYER DIRECTLY BELOW PSP AS LCL AND------------------------
!-----ENSURE THAT THE LCL IS ABOVE GROUND.------------------------------
!-------(IN SOME RARE CASES FOR ITYPE=2, IT IS NOT)---------------------
      DO L=1,LM
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF (PMID(I,J,L) < PSP(I,J))    LCL(I,J) = L+1
          ENDDO
        ENDDO
      ENDDO
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IF (LCL(I,J) > NINT(LMH(I,J))) LCL(I,J) = NINT(LMH(I,J))
          IF (ITYPE  > 2) THEN
            IF (T(I,J,LCL(I,J)) < 263.15) THEN
              THUNDER(I,J) = .FALSE.
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!---------FIND TEMP OF PARCEL LIFTED ALONG MOIST ADIABAT (TPAR)---------
!-----------------------------------------------------------------------

      DO L=LM,1,-1
!--------------SCALING PRESSURE & TT TABLE INDEX------------------------
        KNUML = 0
        KNUMH = 0
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            KLRES(I,J) = 0
            KHRES(I,J) = 0
            IF(L <= LCL(I,J)) THEN
              IF(PMID(I,J,L) < PLQ)THEN
                KNUML = KNUML + 1
                KLRES(I,J) = 1
              ELSE
                KNUMH = KNUMH + 1
                KHRES(I,J) = 1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
!***
!***  COMPUTE PARCEL TEMPERATURE ALONG MOIST ADIABAT FOR PRESSURE<PLQ
!**
        IF(KNUML > 0) THEN
          CALL TTBLEX(TPAR(ISTA_2L,JSTA_2L,L),TTBL,ITB,JTB,KLRES             &
                    , PMID(ISTA_2L,JSTA_2L,L),PL,QQ,PP,RDP,THE0,STHE         &
                    , RDTHE,THESP,IPTB,ITHTB)
        ENDIF
!***
!***  COMPUTE PARCEL TEMPERATURE ALONG MOIST ADIABAT FOR PRESSURE>PLQ
!**
        IF(KNUMH > 0) THEN
          CALL TTBLEX(TPAR(ISTA_2L,JSTA_2L,L),TTBLQ,ITBQ,JTBQ,KHRES          &
                    , PMID(ISTA_2L,JSTA_2L,L),PLQ,QQ,PP,RDPQ                 &
                     ,THE0Q,STHEQ,RDTHEQ,THESP,IPTB,ITHTB)
        ENDIF

!------------SEARCH FOR EQ LEVEL----------------------------------------
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF(KHRES(I,J) > 0) THEN
              IF(TPAR(I,J,L) > T(I,J,L)) IEQL(I,J) = L
            ENDIF
          ENDDO
        ENDDO
!
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF(KLRES(I,J) > 0) THEN
              IF(TPAR(I,J,L) > T(I,J,L) .AND. &
               PMID(I,J,L)>100.) IEQL(I,J) = L
            ENDIF
          ENDDO
        ENDDO
!-----------------------------------------------------------------------
      ENDDO                  ! end of do l=lm,1,-1 loop
!------------COMPUTE CAPE AND CINS--------------------------------------
      LBEG = 1000
      LEND = 0
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          LBEG = MIN(IEQL(I,J),LBEG)
          LEND = MAX(LCL(I,J),LEND)
        ENDDO
      ENDDO
!
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IF(T(I,J,IEQL(I,J)) > 255.65) THEN
            THUNDER(I,J) = .FALSE.
          ENDIF
        ENDDO
      ENDDO
!
      DO L=LBEG,LEND

!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IDX(I,J) = 0
            IF(L >= IEQL(I,J).AND.L <= LCL(I,J)) THEN
              IDX(I,J) = 1
            ENDIF
          ENDDO
        ENDDO
!
!$omp  parallel do private(i,j,gdzkl,presk,thetaa,thetap,esatp,qsatp,tvp,tv)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF(IDX(I,J) > 0) THEN
              PRESK  = PMID(I,J,L)
              GDZKL   = (ZINT(I,J,L)-ZINT(I,J,L+1)) * G
              ESATP  = min(FPVSNEW(TPAR(I,J,L)),PRESK)
              QSATP  = EPS*ESATP/(PRESK-ESATP*ONEPS)
!              TVP    = TPAR(I,J,L)*(1+0.608*QSATP)
              TVP    = TVIRTUAL(TPAR(I,J,L),QSATP)
              THETAP = TVP*(H10E5/PRESK)**CAPA
!              TV     = T(I,J,L)*(1+0.608*Q(I,J,L)) 
              TV     = TVIRTUAL(T(I,J,L),Q(I,J,L))
              THETAA = TV*(H10E5/PRESK)**CAPA
              IF(THETAP < THETAA) THEN
                CINS(I,J) = CINS(I,J) + (LOG(THETAP)-LOG(THETAA))*GDZKL
              ELSEIF(THETAP > THETAA) THEN
                CAPE(I,J) = CAPE(I,J) + (LOG(THETAP)-LOG(THETAA))*GDZKL
                IF (THUNDER(I,J) .AND. T(I,J,L)  < 273.15                 &
                                 .AND. T(I,J,L)  > 253.15) THEN
                 CAPE20(I,J) = CAPE20(I,J) + (LOG(THETAP)-LOG(THETAA))*GDZKL
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!    
!     ENFORCE LOWER LIMIT OF 0.0 ON CAPE AND UPPER
!     LIMIT OF 0.0 ON CINS.
!
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          CAPE(I,J) = MAX(D00,CAPE(I,J))
          CINS(I,J) = MIN(CINS(I,J),D00)
! add equillibrium height
          ZEQL(I,J) = ZINT(I,J,IEQL(I,J))
          TEQL(I,J) = T(I,J,IEQL(I,J))
          IF (CAPE20(I,J) < 75.) THEN
            THUNDER(I,J) = .FALSE.
          ENDIF
          IF (THUNDER(I,J)) THEN
            THUND(I,J) = 1.0
          ELSE
            THUND(I,J) = 0.0
          ENDIF
        ENDDO
      ENDDO
!     
      DEALLOCATE(TPAR)
!     
      END SUBROUTINE CALCAPE
!
!-------------------------------------------------------------------------------------
!> calcape2() computes CAPE and CINS.
!>
!> This routine computes CAPE and CINS given temperature,
!> pressure, and specific humidty.  In "storm and cloud 
!> dynamics" (1989, academic press) cotton and anthes define
!> CAPE (equation 9.16, p501) as
!>
!> @code
!>                  EL
!>        CAPE =  SUM G * ln(THETAP/THETAA) DZ 
!>                 LCL
!>     
!>     Where,
!>      EL    = Equilibrium level,
!>     LCL    = Lifting condenstation level,
!>       G    = Gravitational acceleration,
!>     THETAP = Lifted parcel potential temperature,
!>     THETAA = Ambient potential temperature.
!> @endcode
!>     
!>     Note that the integrand ln(THETAP/THETAA) approximately
!>     equals (THETAP-THETAA)/THETAA.  This ratio is often used
!>     in the definition of CAPE/CINS.
!>     
!>     Two types of CAPE/CINS can be computed by this routine.  The
!>     summation process is the same For both cases.  What differs
!>     is the definition of the parcel to lift.  FOR ITYPE=1 the
!>     parcel with the warmest THETA-E in A DPBND pascal layer above
!>     the model surface is lifted.  the arrays P1D, T1D, and Q1D
!>     are not used.  For itype=2 the arrays P1D, T1D, and Q1D
!>     define the parcel to lift in each column.  Both types of
!>     CAPE/CINS may be computed in a single execution of the post
!>     processor.
!>     
!>     This algorithm proceeds as follows.
!>     For each column, 
!>        (1)  Initialize running CAPE and CINS SUM TO 0.0
!>        (2)  Compute temperature and pressure at the LCL using
!>             look up table (PTBL).  Use either parcel that gives
!>             max THETAE in lowest DPBND above ground (ITYPE=1)
!>             or given parcel from t1D,Q1D,...(ITYPE=2).
!>        (3)  Compute the temp of a parcel lifted from the LCL.
!>             We know that the parcel's
!>             equivalent potential temperature (THESP) remains
!>             constant through this process.  we can
!>             compute tpar using this knowledge using look
!>             up table (subroutine TTBLEX).
!>        (4)  Find the equilibrium level.  This is defined as the
!>             highest positively buoyant layer.
!>             (If there is no positively buoyant layer, CAPE/CINS
!>              will be zero)
!>        (5)  Compute CAPE/CINS.  
!>             (A) Compute THETAP.  We know TPAR and P.
!>             (B) Compute THETAA.  We know T and P.  
!>        (6)  Add G*(THETAP-THETAA)*DZ to the running CAPE or CINS sum.
!>             (A) If THETAP > THETAA, add to the CAPE sum.
!>             (B) If THETAP < THETAA, add to the CINS sum.
!>        (7)  Are we at equilibrium level? 
!>             (A) If yes, stop the summation.
!>             (b) if no, contiunue the summation.
!>        (8)  Enforce limits on CAPE and CINS (i.e. no negative CAPE)
!>
!> @param[in] ITYPE INTEGER Flag specifying how parcel to lift is identified.  See comments above.
!> @param[in] DPBND Depth over which one searches for most unstable parcel.
!> @param[in] P1D Array of pressure of parcels to lift.
!> @param[in] T1D Array of temperature of parcels to lift.
!> @param[in] Q1D Array of specific humidity of parcels to lift.
!> @param[in] L1D Array of model level of parcels to lift.
!> @param[out] CAPE Convective available potential energy (J/kg).
!> @param[out] CINS Convective inhibition (J/kg).
!> @param[out] LFC level of free convection (m).
!> @param[out] ESRHL Lower bound to account for effective helicity calculation.
!> @param[out] ESRHH Upper bound to account for effective helicity calculation.
!> @param[out] DCAPE downdraft CAPE (J/KG).
!> @param[out] DGLD Dendritic growth layer depth (m).
!> @param[out] ESP Enhanced stretching potential.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1993-02-10 | Russ Treadon  | Initial
!> 1993-06-19 | Russ Treadon  | Generalized routine to allow for type 2 CAPE/CINS calculations
!> 1994-09-23 | Mike Baldwin  | Modified to use look up tables instead of complicated equations
!> 1994-10-13 | Mike Baldwin  | Modified to continue CAPE/CINS calc up to at highest buoyant layer
!> 1998-06-12 | T Black       | Conversion from 1-D TO 2-D
!> 1998-08-18 | T Black       | Compute APE internally
!> 2000-01-04 | Jim Tuccillo  | MPI Version              
!> 2002-01-15 | Mike Baldwin  | WRF Version
!> 2003-08-24 | G Manikin     | Added level of parcel being lifted as output from the routine and added the depth over which one searches for the most unstable parcel as input
!> 2010-09-09 | G Manikin     | Changed computation to use virtual temp added eq lvl hght and thunder parameter    
!> 2015-??-?? | S Moorthi     | Optimization and threading
!> 2021-09-03 | J Meng        | Modified to add 0-3km CAPE/CINS, LFC, effective helicity, downdraft CAPE, dendritic growth layer depth, ESP
!> 2021-09-01 | E Colon       | Equivalent level height index for RTMA
!> 2022-08-27 | S Trahan      | Fixed bug in CALCAPE2 where extreme atmospheric conditions cause an out-of-bounds access
!> 2022-09-01 | S Trahan      | Fixed another bug in CALCAPE2 where extreme atmospheric conditions cause an out-of-bounds access
!>
!> @author Russ Treadon W/NP2 @date 1993-02-10
      SUBROUTINE CALCAPE2(ITYPE,DPBND,P1D,T1D,Q1D,L1D,    &  
                          CAPE,CINS,LFC,ESRHL,ESRHH,      &
                          DCAPE,DGLD,ESP)
      use vrbls3d,    only: pmid, t, q, zint
      use vrbls2d,    only: fis,ieql
      use gridspec_mod, only: gridtype
      use masks,      only: lmh
      use params_mod, only: d00, h1m12, h99999, h10e5, capa, elocp, eps,  &
                            oneps, g, tfrz
      use lookup_mod, only: thl, rdth, jtb, qs0, sqs, rdq, itb, ptbl,     &
                            plq, ttbl, pl, rdp, the0, sthe, rdthe, ttblq, &
                            itbq, jtbq, rdpq, the0q, stheq, rdtheq
      use ctlblk_mod, only: jsta_2l, jend_2u, lm, jsta, jend, im, jm, me, jsta_m, jend_m, spval,&
                            ista_2l, iend_2u,     ista, iend,             ista_m, iend_m
!     
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     INCLUDE/SET PARAMETERS.  CONSTANTS ARE FROM BOLTON (MWR, 1980).
      real,PARAMETER :: ISMTHP=2,ISMTHT=2,ISMTHQ=2
!     
!     DECLARE VARIABLES.
!
      integer,intent(in) :: ITYPE
      real,intent(in)    :: DPBND
      integer, dimension(ista:iend,Jsta:jend),intent(in)    :: L1D
      real,    dimension(ista:iend,Jsta:jend),intent(in)    :: P1D,T1D
!      real,    dimension(ista:iend,jsta:jend),intent(inout) :: Q1D,CAPE,CINS,PPARC,ZEQL
      real,    dimension(ista:iend,jsta:jend),intent(inout) :: Q1D,CAPE,CINS
      real,    dimension(ista:iend,jsta:jend)               :: PPARC,ZEQL
      real,    dimension(ista:iend,jsta:jend),intent(inout) :: LFC,ESRHL,ESRHH
      real,    dimension(ista:iend,jsta:jend),intent(inout) :: DCAPE,DGLD,ESP
      integer, dimension(ista:iend,jsta:jend) ::L12,L17,L3KM
!     
      integer, dimension(ista:iend,jsta:jend) :: IPTB, ITHTB, PARCEL, KLRES, KHRES, LCL, IDX
!     
      real,    dimension(ista:iend,jsta:jend) :: THESP, PSP, CAPE20, QQ, PP, THUND  
      integer, dimension(ista:iend,jsta:jend) :: PARCEL2 
      real,    dimension(ista:iend,jsta:jend) :: THESP2,PSP2
      real,    dimension(ista:iend,jsta:jend) :: CAPE4,CINS4 
      REAL, ALLOCATABLE :: TPAR(:,:,:)
      REAL, ALLOCATABLE :: TPAR2(:,:,:)

      LOGICAL THUNDER(ista:iend,jsta:jend), NEEDTHUN 
      real PSFCK,PKL,TBTK,QBTK,APEBTK,TTHBTK,TTHK,APESPK,TPSPK,        &
           BQS00K,SQS00K,BQS10K,SQS10K,BQK,SQK,TQK,PRESK,GDZKL,THETAP, &
           THETAA,P00K,P10K,P01K,P11K,TTHESK,ESATP,QSATP,TVP,TV
      real PRESK2, ESATP2, QSATP2, TVP2, THETAP2, TV2, THETAA2
!      real,external :: fpvsnew
      integer I,J,L,KNUML,KNUMH,LBEG,LEND,IQ, KB,ITTBK
      integer IE,IW,JN,JS,IVE(JM),IVW(JM),JVN,JVS
      integer ISTART,ISTOP,JSTART,JSTOP
      real,    dimension(ista:iend,jsta:jend) :: HTSFC

!     integer I,J,L,KNUML,KNUMH,LBEG,LEND,IQ,IT,LMHK, KB,ITTBK
!     
!**************************************************************
!     START CALCAPE HERE.
!     
      ALLOCATE(TPAR(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U,LM))
      ALLOCATE(TPAR2(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U,LM))
!
!     COMPUTE CAPE/CINS
!
!        WHICH IS: THE SUM FROM THE LCL TO THE EQ LEVEL OF
!             G * (LN(THETAP) - LN(THETAA)) * DZ
!
!             (POSITIVE AREA FOR CAPE, NEGATIVE FOR CINS)
!
!        WHERE:
!             THETAP IS THE PARCEL THETA
!             THETAA IS THE AMBIENT THETA
!             DZ IS THE THICKNESS OF THE LAYER
!
!         USING LCL AS LEVEL DIRECTLY BELOW SATURATION POINT
!         AND EQ LEVEL IS THE HIGHEST POSITIVELY BUOYANT LEVEL.
!  
!         IEQL = EQ LEVEL
!         P_thetaemax - real  pressure of theta-e max parcel (Pa)
!
!     INITIALIZE CAPE AND CINS ARRAYS
! 
!$omp  parallel do
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          CAPE(I,J)    = D00
          CAPE20(I,J)  = D00
          CAPE4(I,J)   = D00
          CINS(I,J)    = D00
          CINS4(I,J)   = D00
          LCL(I,J)     = 0
          THESP(I,J)   = D00
          IEQL(I,J)    = LM
          PARCEL(I,J)  = LM
          PSP(I,J)     = D00
          PPARC(I,J)   = D00
          THUNDER(I,J) = .TRUE.
          LFC(I,J)     = D00
          ESRHL(I,J)   = D00
          ESRHH(I,J)   = D00
          DCAPE(I,J)   = D00
          DGLD(I,J)    = D00
          ESP(I,J)     = D00
          THESP2(I,J)  = 1000.
          PSP2(I,J)    = D00
          PARCEL2(I,J) = LM
        ENDDO
      ENDDO
!
!$omp  parallel do
      DO L=1,LM
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            TPAR(I,J,L) = D00
            TPAR2(I,J,L) = D00
          ENDDO
        ENDDO
      ENDDO
!
!     FIND SURFACE HEIGHT
!
      IF(gridtype == 'E')THEN
        JVN =  1
        JVS = -1
        do J=JSTA,JEND
          IVE(J) = MOD(J,2)
          IVW(J) = IVE(J)-1
        enddo
        ISTART = ISTA_M
        ISTOP  = IEND_M
        JSTART = JSTA_M
        JSTOP  = JEND_M
      ELSE IF(gridtype == 'B')THEN
        JVN = 1
        JVS = 0
        do J=JSTA,JEND
          IVE(J)=1
          IVW(J)=0
        enddo
        ISTART = ISTA_M
        ISTOP  = IEND_M
        JSTART = JSTA_M
        JSTOP  = JEND_M
      ELSE
        JVN = 0
        JVS = 0
        do J=JSTA,JEND
          IVE(J) = 0
          IVW(J) = 0
        enddo
        ISTART = ISTA
        ISTOP  = IEND
        JSTART = JSTA
        JSTOP  = JEND
      END IF
!!$omp  parallel do private(htsfc,ie,iw)
      IF(gridtype /= 'A') CALL EXCH(FIS(ISTA:IEND,JSTA:JEND))
        DO J=JSTART,JSTOP
          DO I=ISTART,ISTOP
            IE = I+IVE(J)
            IW = I+IVW(J)
            JN = J+JVN
            JS = J+JVS
!mp          PDSLVK=(PD(IW,J)*RES(IW,J)+PD(IE,J)*RES(IE,J)+
!mp     1           PD(I,J+1)*RES(I,J+1)+PD(I,J-1)*RES(I,J-1))*0.25
!mp          PSFCK=AETA(LMV(I,J))*PDSLVK+PT
            IF (gridtype=='B')THEN
              HTSFC(I,J) = (0.25/g)*(FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(IE,JN))
            ELSE
              HTSFC(I,J) = (0.25/g)*(FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(I,JS))
            ENDIF
          ENDDO
        ENDDO
!
!     TYPE 2 CAPE/CINS:
!     NOTE THAT FOR TYPE 1 CAPE/CINS ARRAYS P1D, T1D, Q1D 
!     ARE DUMMY ARRAYS.
!     
      IF (ITYPE == 2) THEN
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            Q1D(I,J) = MIN(MAX(H1M12,Q1D(I,J)),H99999)
          ENDDO
        ENDDO
      ENDIF
!-------FOR ITYPE=1--FIND MAXIMUM THETA E LAYER IN LOWEST DPBND ABOVE GROUND-------
!-------FOR ITYPE=2--FIND THETA E LAYER OF GIVEN T1D, Q1D, P1D---------------------
!--------------TRIAL MAXIMUM BUOYANCY LEVEL VARIABLES-------------------

      DO KB=1,LM
!hc     IF (ITYPE==2.AND.KB>1) cycle
        IF (ITYPE == 1 .OR. (ITYPE == 2 .AND. KB == 1)) THEN

!$omp  parallel do private(i,j,apebtk,apespk,bqk,bqs00k,bqs10k,iq,ittbk,    &
!$omp &         p00k,p01k,p10k,p11k,pkl,psfck,qbtk,sqk,sqs00k,              &
!$omp &         sqs10k,tbtk,tpspk,tqk,tthbtk,tthesk,tthk)
          DO J=JSTA,JEND
            DO I=ISTA,IEND
              PSFCK  = PMID(I,J,NINT(LMH(I,J)))
              PKL    = PMID(I,J,KB)

!hc           IF (ITYPE==1.AND.(PKL<PSFCK-DPBND.OR.PKL>PSFCK)) cycle
              IF (ITYPE ==2 .OR.                                                &
                 (ITYPE == 1 .AND. (PKL >= PSFCK-DPBND .AND. PKL <= PSFCK)))THEN
                IF (ITYPE == 1) THEN
                  TBTK   = T(I,J,KB)
                  QBTK   = max(0.0, Q(I,J,KB))
                  APEBTK = (H10E5/PKL)**CAPA
                ELSE
                  PKL    = P1D(I,J)
                  TBTK   = T1D(I,J)
                  QBTK   = max(0.0, Q1D(I,J))
                  APEBTK = (H10E5/PKL)**CAPA
                ENDIF

!----------Breogan Gomez - 2009-02-06
! To prevent QBTK to be less than 0 which leads to a unrealistic value of PRESK
!  and a floating invalid.

!               if(QBTK < 0) QBTK = 0

!--------------SCALING POTENTIAL TEMPERATURE & TABLE INDEX--------------
                TTHBTK  =  TBTK*APEBTK
                TTHK    = (TTHBTK-THL)*RDTH
                QQ(I,J) = TTHK - AINT(TTHK)
                ITTBK   = INT(TTHK) + 1
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
                IF(ITTBK < 1)   THEN
                  ITTBK   = 1
                  QQ(I,J) = D00
                ENDIF
                IF(ITTBK >= JTB) THEN
                  ITTBK   = JTB-1
                  QQ(I,J) = D00
                ENDIF
!--------------BASE AND SCALING FACTOR FOR SPEC. HUMIDITY---------------
                BQS00K = QS0(ITTBK)
                SQS00K = SQS(ITTBK)
                BQS10K = QS0(ITTBK+1)
                SQS10K = SQS(ITTBK+1)
!--------------SCALING SPEC. HUMIDITY & TABLE INDEX---------------------
                BQK     = (BQS10K-BQS00K)*QQ(I,J) + BQS00K
                SQK     = (SQS10K-SQS00K)*QQ(I,J) + SQS00K
                TQK     = (QBTK-BQK)/SQK*RDQ
                PP(I,J) = TQK-AINT(TQK)
                IQ      = INT(TQK)+1
!--------------KEEPING INDICES WITHIN THE TABLE-------------------------
                IF(IQ < 1)    THEN
                  IQ      = 1
                  PP(I,J) = D00
                ENDIF
                IF(IQ >= ITB)  THEN
                  IQ      = ITB-1
                  PP(I,J) = D00
                ENDIF
!--------------SATURATION PRESSURE AT FOUR SURROUNDING TABLE PTS.-------
                P00K = PTBL(IQ  ,ITTBK  )
                P10K = PTBL(IQ+1,ITTBK  )
                P01K = PTBL(IQ  ,ITTBK+1)
                P11K = PTBL(IQ+1,ITTBK+1)
!--------------SATURATION POINT VARIABLES AT THE BOTTOM-----------------
                TPSPK = P00K + (P10K-P00K)*PP(I,J) + (P01K-P00K)*QQ(I,J)  &
                             + (P00K-P10K-P01K+P11K)*PP(I,J)*QQ(I,J)

!!from WPP::tgs          APESPK=(H10E5/TPSPK)**CAPA
                if (TPSPK > 1.0e-3) then
                  APESPK = (max(0.,H10E5/ TPSPK))**CAPA
                else
                  APESPK = 0.0
                endif

                TTHESK = TTHBTK * EXP(ELOCP*QBTK*APESPK/TTHBTK)
!--------------CHECK FOR MAXIMUM THETA E--------------------------------
                IF(TTHESK > THESP(I,J)) THEN
                  PSP  (I,J)  = TPSPK
                  THESP(I,J)  = TTHESK
                  PARCEL(I,J) = KB
                ENDIF
!--------------CHECK FOR MINIMUM THETA E--------------------------------
                IF(TTHESK < THESP2(I,J)) THEN
                  PSP2  (I,J)  = TPSPK
                  THESP2(I,J)  = TTHESK
                  PARCEL2(I,J) = KB
                ENDIF
              END IF 
            ENDDO  ! I  loop
          ENDDO    ! J  loop
        END IF
      ENDDO        ! KB loop

!----FIND THE PRESSURE OF THE PARCEL THAT WAS LIFTED
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            PPARC(I,J) = PMID(I,J,PARCEL(I,J))
          ENDDO
        ENDDO
!
!-----CHOOSE LAYER DIRECTLY BELOW PSP AS LCL AND------------------------
!-----ENSURE THAT THE LCL IS ABOVE GROUND.------------------------------
!-------(IN SOME RARE CASES FOR ITYPE=2, IT IS NOT)---------------------
      DO L=1,LM
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF (PMID(I,J,L) < PSP(I,J))    LCL(I,J) = L+1
          ENDDO
        ENDDO
      ENDDO
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IF (LCL(I,J) > NINT(LMH(I,J))) LCL(I,J) = NINT(LMH(I,J))
          IF (ITYPE  > 2) THEN
            IF (T(I,J,LCL(I,J)) < 263.15) THEN
              THUNDER(I,J) = .FALSE.
            ENDIF
          ENDIF

          ! Limit LCL to prevent out-of-bounds accesses later
          LCL(I,J) = max(min(LCL(I,J),LM-1),1)
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!---------FIND TEMP OF PARCEL LIFTED ALONG MOIST ADIABAT (TPAR)---------
!-----------------------------------------------------------------------
      DO L=LM,1,-1
!--------------SCALING PRESSURE & TT TABLE INDEX------------------------
        KNUML = 0
        KNUMH = 0
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            KLRES(I,J) = 0
            KHRES(I,J) = 0
            IF(L <= LCL(I,J)) THEN
              IF(PMID(I,J,L) < PLQ)THEN
                KNUML = KNUML + 1
                KLRES(I,J) = 1
              ELSE
                KNUMH = KNUMH + 1
                KHRES(I,J) = 1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
!***
!***  COMPUTE PARCEL TEMPERATURE ALONG MOIST ADIABAT FOR PRESSURE<PLQ
!**
        IF(KNUML > 0) THEN
          CALL TTBLEX(TPAR(ISTA_2L,JSTA_2L,L),TTBL,ITB,JTB,KLRES             &
                    , PMID(ISTA_2L,JSTA_2L,L),PL,QQ,PP,RDP,THE0,STHE         &
                    , RDTHE,THESP,IPTB,ITHTB)
        ENDIF
!***
!***  COMPUTE PARCEL TEMPERATURE ALONG MOIST ADIABAT FOR PRESSURE>PLQ
!**
        IF(KNUMH > 0) THEN
          CALL TTBLEX(TPAR(ISTA_2L,JSTA_2L,L),TTBLQ,ITBQ,JTBQ,KHRES          &
                    , PMID(ISTA_2L,JSTA_2L,L),PLQ,QQ,PP,RDPQ                 &
                     ,THE0Q,STHEQ,RDTHEQ,THESP,IPTB,ITHTB)
        ENDIF

!------------SEARCH FOR EQ LEVEL----------------------------------------
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF(KHRES(I,J) > 0) THEN
              IF(TPAR(I,J,L) > T(I,J,L)) IEQL(I,J) = L
            ENDIF
          ENDDO
        ENDDO
!
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF(KLRES(I,J) > 0) THEN
              IF(TPAR(I,J,L) > T(I,J,L)) IEQL(I,J) = L
            ENDIF
          ENDDO
        ENDDO
!-----------------------------------------------------------------------
      ENDDO                  ! end of do l=lm,1,-1 loop
!------------COMPUTE CAPE AND CINS--------------------------------------
      LBEG = 1000
      LEND = 0
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          LBEG = MIN(IEQL(I,J),LBEG)
          LEND = MAX(LCL(I,J),LEND)
        ENDDO
      ENDDO
!
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IF(T(I,J,IEQL(I,J)) > 255.65) THEN
            THUNDER(I,J) = .FALSE.
          ENDIF
        ENDDO
      ENDDO
!
!Ensure later calculations do not access LM+1
!
      LEND=MIN(LEND,LM-1)
!
!reverse L order from bottom up for ESRH calculation
!
      ESRHH = LCL
      ESRHL = LCL
!      DO L=LBEG,LEND
      DO L=LEND,LBEG,-1

!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IDX(I,J) = 0
            IF(L >= IEQL(I,J).AND.L <= LCL(I,J)) THEN
              IDX(I,J) = 1
            ENDIF
          ENDDO
        ENDDO
!
!$omp  parallel do private(i,j,gdzkl,presk,thetaa,thetap,esatp,qsatp,tvp,tv,&
!$omp &                    presk2,esatp2,qsatp2,tvp2,thetap2,tv2,thetaa2)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF(IDX(I,J) > 0) THEN
              PRESK  = PMID(I,J,L)
              GDZKL  = (ZINT(I,J,L)-ZINT(I,J,L+1)) * G
              ESATP  = min(FPVSNEW(TPAR(I,J,L)),PRESK)
              QSATP  = EPS*ESATP/(PRESK-ESATP*ONEPS)
!              TVP    = TPAR(I,J,L)*(1+0.608*QSATP)
              TVP    = TVIRTUAL(TPAR(I,J,L),QSATP)
              THETAP = TVP*(H10E5/PRESK)**CAPA
!              TV     = T(I,J,L)*(1+0.608*Q(I,J,L)) 
              TV     = TVIRTUAL(T(I,J,L),Q(I,J,L))
              THETAA = TV*(H10E5/PRESK)**CAPA
              IF(THETAP < THETAA) THEN
                CINS4(I,J) = CINS4(I,J) + (LOG(THETAP)-LOG(THETAA))*GDZKL
               IF(ZINT(I,J,L)-HTSFC(I,J) <= 3000.) THEN
                CINS(I,J) = CINS(I,J) + (LOG(THETAP)-LOG(THETAA))*GDZKL
               ENDIF
              ELSEIF(THETAP > THETAA) THEN
                CAPE4(I,J) = CAPE4(I,J) + (LOG(THETAP)-LOG(THETAA))*GDZKL
               IF(ZINT(I,J,L)-HTSFC(I,J) <= 3000.) THEN
                CAPE(I,J) = CAPE(I,J) + (LOG(THETAP)-LOG(THETAA))*GDZKL
               ENDIF
               IF (THUNDER(I,J) .AND. T(I,J,L)  < 273.15                 &
                                 .AND. T(I,J,L)  > 253.15) THEN
                 CAPE20(I,J) = CAPE20(I,J) + (LOG(THETAP)-LOG(THETAA))*GDZKL
               ENDIF
              ENDIF
              
! LFC
              IF (ITYPE /= 1) THEN
               PRESK2  = PMID(I,J,L+1)
               ESATP2  = min(FPVSNEW(TPAR(I,J,L+1)),PRESK2)
               QSATP2  = EPS*ESATP2/(PRESK2-ESATP2*ONEPS)
!               TVP2    = TPAR(I,J,L+1)*(1+0.608*QSATP2)
               TVP2    = TVIRTUAL(TPAR(I,J,L+1),QSATP2)
               THETAP2 = TVP2*(H10E5/PRESK2)**CAPA
!               TV2     = T(I,J,L+1)*(1+0.608*Q(I,J,L+1))
               TV2     = TVIRTUAL(T(I,J,L+1),Q(I,J,L+1))
               THETAA2 = TV2*(H10E5/PRESK2)**CAPA
               IF(THETAP >= THETAA .AND. THETAP2 <= THETAA2) THEN
                IF(LFC(I,J) == D00)THEN
                   LFC(I,J) = ZINT(I,J,L)
                ENDIF
               ENDIF
              ENDIF
!
! ESRH/CAPE threshold check
              IF(ZINT(I,J,L)-HTSFC(I,J) <= 3000.) THEN
                IF(CAPE4(I,J) >= 100. .AND. CINS4(I,J) >= -250.) THEN
                   IF(ESRHL(I,J) == LCL(I,J)) ESRHL(I,J)=L
                ENDIF
                ESRHH(I,J)=L
              ENDIF 

            ENDIF !(IDX(I,J) > 0)
          ENDDO
        ENDDO
      ENDDO

!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
              IF(ESRHH(I,J) > ESRHL(I,J)) ESRHH(I,J)=IEQL(I,J)
          ENDDO
        ENDDO
!    
!     ENFORCE LOWER LIMIT OF 0.0 ON CAPE AND UPPER
!     LIMIT OF 0.0 ON CINS.
!     ENFORCE LFC ABOVE LCL AND BELOW EL
!
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          CAPE(I,J) = MAX(D00,CAPE(I,J))
          CINS(I,J) = MIN(CINS(I,J),D00)
! equillibrium height
          ZEQL(I,J) = ZINT(I,J,IEQL(I,J))
          LFC(I,J)  = MIN(LFC(I,J),ZINT(I,J,IEQL(I,J)))
          LFC(I,J)  = MAX(ZINT(I,J, LCL(I,J)),LFC(I,J))
          IF (CAPE20(I,J) < 75.) THEN
            THUNDER(I,J) = .FALSE.
          ENDIF
          IF (THUNDER(I,J)) THEN
            THUND(I,J) = 1.0
          ELSE
            THUND(I,J) = 0.0
          ENDIF
        ENDDO
      ENDDO
!------------COMPUTE DCAPE--------------------------------------
!-----------------------------------------------------------------------
!---------FIND TEMP OF PARCEL DESCENDED ALONG MOIST ADIABAT (TPAR)---------
!-----------------------------------------------------------------------
      IF (ITYPE == 1) THEN

      DO L=LM,1,-1
!--------------SCALING PRESSURE & TT TABLE INDEX------------------------
        KNUML = 0
        KNUMH = 0
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            KLRES(I,J) = 0
            KHRES(I,J) = 0
            PSFCK  = PMID(I,J,NINT(LMH(I,J)))
            PKL    = PMID(I,J,L)
            IF(PKL >= PSFCK-DPBND) THEN
              IF(PMID(I,J,L) < PLQ)THEN
                KNUML = KNUML + 1
                KLRES(I,J) = 1
              ELSE
                KNUMH = KNUMH + 1
                KHRES(I,J) = 1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
!***
!***  COMPUTE PARCEL TEMPERATURE ALONG MOIST ADIABAT FOR PRESSURE<PLQ
!**
        IF(KNUML > 0) THEN
          CALL TTBLEX(TPAR2(ISTA_2L,JSTA_2L,L),TTBL,ITB,JTB,KLRES             &
                    , PMID(ISTA_2L,JSTA_2L,L),PL,QQ,PP,RDP,THE0,STHE         &
                    , RDTHE,THESP2,IPTB,ITHTB)
        ENDIF
!***
!***  COMPUTE PARCEL TEMPERATURE ALONG MOIST ADIABAT FOR PRESSURE>PLQ
!**
        IF(KNUMH > 0) THEN
          CALL TTBLEX(TPAR2(ISTA_2L,JSTA_2L,L),TTBLQ,ITBQ,JTBQ,KHRES          &
                    , PMID(ISTA_2L,JSTA_2L,L),PLQ,QQ,PP,RDPQ                 &
                    , THE0Q,STHEQ,RDTHEQ,THESP2,IPTB,ITHTB)
        ENDIF
      ENDDO                  ! end of do l=lm,1,-1 loop

      LBEG = 1
      LEND = LM

      DO L=LBEG,LEND
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IDX(I,J) = 0
            IF(L >= PARCEL2(I,J).AND.L < NINT(LMH(I,J))) THEN
              IDX(I,J) = 1
            ENDIF
          ENDDO
        ENDDO
!
!$omp  parallel do private(i,j,gdzkl,presk,thetaa,thetap,esatp,qsatp,tvp,tv)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF(IDX(I,J) > 0) THEN
              PRESK  = PMID(I,J,L)
              GDZKL  = (ZINT(I,J,L)-ZINT(I,J,L+1)) * G
              ESATP  = min(FPVSNEW(TPAR2(I,J,L)),PRESK)
              QSATP  = EPS*ESATP/(PRESK-ESATP*ONEPS)
!              TVP    = TPAR2(I,J,L)*(1+0.608*QSATP)
              TVP    = TVIRTUAL(TPAR2(I,J,L),QSATP)
              THETAP = TVP*(H10E5/PRESK)**CAPA
!              TV     = T(I,J,L)*(1+0.608*Q(I,J,L))
              TV     = TVIRTUAL(T(I,J,L),Q(I,J,L))
              THETAA = TV*(H10E5/PRESK)**CAPA
              !IF(THETAP < THETAA) THEN
                DCAPE(I,J) = DCAPE(I,J) + (LOG(THETAP)-LOG(THETAA))*GDZKL
              !ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          DCAPE(I,J) = MIN(D00,DCAPE(I,J))
        ENDDO
      ENDDO

      ENDIF !ITYPE=1 FOR DCAPE

!
! Dendritic Growth Layer depth
! the layer with temperatures from -12 to -17 C in meters
!
      L12=LM
      L17=LM
      DO L=LM,1,-1
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF(T(I,J,L) <= TFRZ-12. .AND. L12(I,J)==LM) L12(I,J)=L
            IF(T(I,J,L) <= TFRZ-17. .AND. L17(I,J)==LM) L17(I,J)=L
          ENDDO
        ENDDO
      ENDDO
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
           IF(L12(I,J)/=LM .AND. L17(I,J)/=LM) THEN
             DGLD(I,J)=ZINT(I,J,L17(I,J))-ZINT(I,J,L12(I,J))
             DGLD(I,J)=MAX(DGLD(I,J),0.)
           ENDIF
        ENDDO
      ENDDO
!
! Enhanced Stretching Potential
! ESP = (0-3 km MLCAPE / 50 J kg-1) * ((0-3 km lapse rate - 7.0) / 1.0 C (km-1)
! https://www.spc.noaa.gov/exper/soundings/help/params4.html
!
      L3KM=LM
      DO L=LM,1,-1
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            IF(ZINT(I,J,L)-HTSFC(I,J) <= 3000.) L3KM(I,J)=L
          ENDDO
        ENDDO
      ENDDO     
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
             ESP(I,J) = (CAPE(I,J) / 50.) * (T(I,J,LM) - T(I,J,L3KM(I,J)) - 7.0)
             IF((T(I,J,LM) - T(I,J,L3KM(I,J))) < 7.0) ESP(I,J) = 0.
!             IF(CAPE(I,J) < 250.) ESP(I,J) = 0.
        ENDDO
      ENDDO
!
      DEALLOCATE(TPAR)
      DEALLOCATE(TPAR2)
!     
      END SUBROUTINE CALCAPE2
!
!-------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------
!
      elemental function TVIRTUAL(T,Q)
!
!> TVIRTUAL() Computes virtual temperature
!>
!> @param[in] T real Temperature
!> @param[in] Q real Specific humidity
!> @return virtual temperature
!
      IMPLICIT NONE
      REAL TVIRTUAL
      REAL, INTENT(IN) :: T, Q

      TVIRTUAL = T*(1+0.608*Q)

      end function TVIRTUAL
!
!-------------------------------------------------------------------------------------
!
!> CALVOR() computes absolute vorticity. 
!>
!> @param[in] UWND U wind (m/s) mass-points.
!> @param[in] VWND V wind (m/s) mass-points.
!> @param[out] ABSV absolute vorticity (1/s) mass-points.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-22 | Russ Treadon | Initial
!> 1998-06-08 | T Black      | Convesion from 1-D to 2-D
!> 2000-01-04 | Jim Tuccillo | MPI Version            
!> 2002-01-15 | Mike Baldwin | WRF Version C-grid
!> 2005-03-01 | H Chuang     | Add NMM E grid
!> 2005-05-17 | H Chuang     | Add Potential vorticity calculation
!> 2005-07-07 | B Zhou       | Add RSM in computing DVDX, DUDY and UAVG
!> 2013-08-09 | S Moorthi    | Optimize the vorticity loop including threading
!> 2016-08-05 | S Moorthi    | add zonal filetering
!> 2019-10-17 | Y Mao        | Skip calculation when U/V is SPVAL
!> 2020-11-06 | J Meng       | Use UPP_MATH Module
!> 2022-05-26 | H Chuang     | Use GSL approach for FV3R
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22

      SUBROUTINE CALVOR(UWND,VWND,ABSV)
   
!
      use vrbls2d,      only: f
      use masks,        only: gdlat, gdlon, dx, dy
      use params_mod,   only: d00, dtr, small, erad
      use ctlblk_mod,   only: jsta_2l, jend_2u, spval, modelname, global, &
                              jsta, jend, im, jm, jsta_m, jend_m, gdsdegr,&
                              ista, iend, ista_m, iend_m, ista_2l, iend_2u, me, num_procs
      use gridspec_mod, only: gridtype, dyval
      use upp_math,     only: DVDXDUDY, DDVDX, DDUDY, UUAVG

      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in)    :: UWND, VWND
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(inout) :: ABSV
      REAL, dimension(IM,2) :: GLATPOLES, COSLPOLES, UPOLES, AVPOLES
      REAL, dimension(IM,JSTA:JEND) :: COSLTEMP, AVTEMP
!
      real,    allocatable ::  wrk1(:,:), wrk2(:,:), wrk3(:,:), cosl(:,:)
      INTEGER, allocatable ::  IHE(:),IHW(:), IE(:),IW(:)
!
      integer, parameter :: npass2=2, npass3=3
      integer I,J,ip1,im1,ii,iir,iil,jj,JMT2,imb2, npass, nn, jtem
      real    R2DX,R2DY,DVDX,DUDY,UAVG,TPH1,TPHI, tx1(im+2), tx2(im+2)
!     
!***************************************************************************
!     START CALVOR HERE.
!     
!     LOOP TO COMPUTE ABSOLUTE VORTICITY FROM WINDS.
!     
      IF(MODELNAME  == 'RAPR') then
!$omp  parallel do private(i,j)
        DO J=JSTA_2L,JEND_2U
          DO I=ISTA_2L,IEND_2U
            ABSV(I,J) = D00
          ENDDO
        ENDDO
      else
!$omp  parallel do private(i,j)
        DO J=JSTA_2L,JEND_2U
          DO I=ISTA_2L,IEND_2U
            ABSV(I,J) = SPVAL
          ENDDO
        ENDDO
      endif

      CALL EXCH(UWND)
      CALL EXCH(VWND)
!
      IF (MODELNAME == 'GFS' .or. global) THEN
        CALL EXCH(GDLAT(ISTA_2L,JSTA_2L))
        CALL EXCH(GDLON(ISTA_2L,JSTA_2L))

        allocate (wrk1(ista:iend,jsta:jend), wrk2(ista:iend,jsta:jend),          &
     &            wrk3(ista:iend,jsta:jend), cosl(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(iw(im),ie(im))

        imb2 = im/2
!$omp  parallel do private(i)
      do i=ista,iend
        ie(i) = i+1
        iw(i) = i-1
      enddo
!      iw(1)  = im
!      ie(im) = 1

!       if(1>=jsta .and. 1<=jend)then
!        if(cos(gdlat(1,1)*dtr)<small)poleflag=.T.
!       end if 	
!       call mpi_bcast(poleflag,1,MPI_LOGICAL,0,mpi_comm_comp,iret)

!$omp  parallel do private(i,j,ip1,im1)
        DO J=JSTA,JEND
          do i=ista,iend
            ip1 = ie(i)
            im1 = iw(i)
            cosl(i,j) = cos(gdlat(i,j)*dtr)
            IF(cosl(i,j) >= SMALL) then
              wrk1(i,j) = 1.0 / (ERAD*cosl(i,j))
            else
              wrk1(i,j) = 0.
            end if    
            if(i == im .or. i == 1) then
              wrk2(i,j) = 1.0 / ((360.+GDLON(ip1,J)-GDLON(im1,J))*DTR) !1/dlam
            else
              wrk2(i,j) = 1.0 / ((GDLON(ip1,J)-GDLON(im1,J))*DTR)      !1/dlam
            end if
          enddo
        enddo
!       CALL EXCH(cosl(1,JSTA_2L))
        CALL EXCH(cosl)

        call fullpole( cosl(ista_2l:iend_2u,jsta_2l:jend_2u),coslpoles)
        call fullpole(gdlat(ista_2l:iend_2u,jsta_2l:jend_2u),glatpoles)

!$omp  parallel do private(i,j,ii)
        DO J=JSTA,JEND
          if (j == 1) then
           if(gdlat(ista,j) > 0.) then ! count from north to south
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
          !     wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J+1)-GDLAT(II,J))*DTR) !1/dphi
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J+1)-GLATPOLES(ii,1))*DTR) !1/dphi
              enddo
            else ! count from south to north
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
          !     wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J+1)+GDLAT(II,J))*DTR) !1/dphi
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J+1)+GLATPOLES(ii,1))*DTR) !1/dphi
!
              enddo
            end if      
          elseif (j == JM) then
            if(gdlat(ista,j) < 0.) then ! count from north to south
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
          !      wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J-1)+GDLAT(II,J))*DTR)
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J-1)+GLATPOLES(ii,2))*DTR)
              enddo
            else ! count from south to north
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
          !     wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J-1)-GDLAT(II,J))*DTR)
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J-1)-GLATPOLES(ii,2))*DTR)
              enddo
            end if  
          else
            do i=ista,iend
              wrk3(i,j) = 1.0 / ((GDLAT(I,J-1)-GDLAT(I,J+1))*DTR) !1/dphi
            enddo
          endif
        enddo  

        npass = 0

        jtem = jm / 18 + 1
      
        call fullpole(uwnd(ista_2l:iend_2u,jsta_2l:jend_2u),upoles)

!$omp  parallel do private(i,j,ip1,im1,ii,jj,tx1,tx2)
        DO J=JSTA,JEND
!         npass = npass2
!         if (j > jm-jtem+1 .or. j < jtem) npass = npass3
          IF(J == 1) then                            ! Near North or South pole
            if(gdlat(ista,j) > 0.) then ! count from north to south
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VWND(ip1,J)==SPVAL .or. VWND(im1,J)==SPVAL .or. &
!                    UWND(II,J)==SPVAL .or. UWND(I,J+1)==SPVAL) cycle
                     UPOLES(II,1)==SPVAL .or. UWND(I,J+1)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)               &
!    &                      +  (UWND(II,J)*COSL(II,J)                            &
     &                      +  (upoles(II,1)*coslpoles(II,1)                            &
     &                      +   UWND(I,J+1)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)  &
     &                      + F(I,J)
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VWND(ip1,JJ)==SPVAL .or. VWND(im1,JJ)==SPVAL .or. &
                     UWND(I,J)==SPVAL .or. UWND(I,jj+1)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,JJ)-VWND(im1,JJ))*wrk2(i,jj)               &
     &                      -  (UWND(I,J)*COSL(I,J)                                 &
                            -   UWND(I,jj+1)*COSL(I,Jj+1))*wrk3(i,jj)) * wrk1(i,jj) &
     &                      + F(I,Jj)
                enddo
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VWND(ip1,J)==SPVAL .or. VWND(im1,J)==SPVAL .or. &
!                    UWND(II,J)==SPVAL .or. UWND(I,J+1)==SPVAL) cycle
                     UPOLES(II,1)==SPVAL .or. UWND(I,J+1)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)               &
!    &                      -  (UWND(II,J)*COSL(II,J)                            &
     &                      -  (upoles(II,1)*coslpoles(II,1)                            &
     &                      +   UWND(I,J+1)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)  &
     &                      + F(I,J)
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VWND(ip1,JJ)==SPVAL .or. VWND(im1,JJ)==SPVAL .or. &
                     UWND(I,J)==SPVAL .or. UWND(I,jj+1)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,JJ)-VWND(im1,JJ))*wrk2(i,jj)               &
     &                      +  (UWND(I,J)*COSL(I,J)                                 &
                            -   UWND(I,jj+1)*COSL(I,Jj+1))*wrk3(i,jj)) * wrk1(i,jj) &
     &                      + F(I,Jj)
                enddo
              ENDIF
            endif
          ELSE IF(J == JM) THEN                      ! Near North or South Pole
            if(gdlat(ista,j) < 0.) then ! count from north to south
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VWND(ip1,J)==SPVAL .or. VWND(im1,J)==SPVAL .or. &
!                    UWND(I,J-1)==SPVAL .or. UWND(II,J)==SPVAL) cycle
                     UWND(I,J-1)==SPVAL .or. UPOLES(II,2)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)              &
     &                      -  (UWND(I,J-1)*COSL(I,J-1)                         &
!    &                      +   UWND(II,J)*COSL(II,J))*wrk3(i,j)) * wrk1(i,j)   &
     &                      +   upoles(II,2)*coslpoles(II,2))*wrk3(i,j)) * wrk1(i,j)   &
     &                      + F(I,J)
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VWND(ip1,JJ)==SPVAL .or. VWND(im1,JJ)==SPVAL .or. &
                     UWND(I,jj-1)==SPVAL .or. UWND(I,J)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,JJ)-VWND(im1,JJ))*wrk2(i,jj)         &
     &                      -  (UWND(I,jj-1)*COSL(I,Jj-1)                     &
     &                      -   UWND(I,J)*COSL(I,J))*wrk3(i,jj)) * wrk1(i,jj) &
     &                      + F(I,Jj)
                enddo
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VWND(ip1,J)==SPVAL .or. VWND(im1,J)==SPVAL .or. &
!                    UWND(I,J-1)==SPVAL .or. UWND(II,J)==SPVAL) cycle
                     UWND(I,J-1)==SPVAL .or. UPOLES(II,2)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)              &
     &                      +  (UWND(I,J-1)*COSL(I,J-1)                         &
!    &                      +   UWND(II,J)*COSL(II,J))*wrk3(i,j)) * wrk1(i,j)   &
     &                      +   upoles(II,2)*coslpoles(II,2))*wrk3(i,j)) * wrk1(i,j)   &
     &                      + F(I,J)
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VWND(ip1,JJ)==SPVAL .or. VWND(im1,JJ)==SPVAL .or. &
                     UWND(I,jj-1)==SPVAL .or. UWND(I,J)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,JJ)-VWND(im1,JJ))*wrk2(i,jj)         &
     &                      +  (UWND(I,jj-1)*COSL(I,Jj-1)                     &
     &                      -   UWND(I,J)*COSL(I,J))*wrk3(i,jj)) * wrk1(i,jj) &
     &                      + F(I,Jj)
                enddo
              ENDIF
            endif
          ELSE
            DO I=ISTA,IEND
              ip1 = ie(i)
              im1 = iw(i)
              if(VWND(ip1,J)==SPVAL .or. VWND(im1,J)==SPVAL .or. &
                 UWND(I,J-1)==SPVAL .or. UWND(I,J+1)==SPVAL) cycle
              ABSV(I,J)   = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)               &
     &                    -  (UWND(I,J-1)*COSL(I,J-1)                          &
                          -   UWND(I,J+1)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)  &
                          + F(I,J)
            ENDDO
          END IF
          if (npass > 0) then
            do i=ista,iend
              tx1(i) = absv(i,j)
            enddo
            do nn=1,npass
              do i=ista,iend
                tx2(i+1) = tx1(i)
              enddo
              tx2(1)    = tx2(im+1)
              tx2(im+2) = tx2(2)
              do i=2,im+1
                tx1(i-1) = 0.25 * (tx2(i-1) + tx2(i+1)) + 0.5*tx2(i)
              enddo
            enddo
            do i=ista,iend
              absv(i,j) = tx1(i)
            enddo
          endif
        END DO                               ! end of J loop

!       deallocate (wrk1, wrk2, wrk3, cosl)
! GFS use lon avg as one scaler value for pole point

      ! call poleavg(IM,JM,JSTA,JEND,SMALL,COSL(1,jsta),SPVAL,ABSV(1,jsta))

        call exch(absv(ista_2l:iend_2u,jsta_2l:jend_2u))
        call fullpole(absv(ista_2l:iend_2u,jsta_2l:jend_2u),avpoles)     

        cosltemp=spval
        if(jsta== 1) cosltemp(1:im, 1)=coslpoles(1:im,1)
        if(jend==jm) cosltemp(1:im,jm)=coslpoles(1:im,2)
        avtemp=spval
        if(jsta== 1) avtemp(1:im, 1)=avpoles(1:im,1)
        if(jend==jm) avtemp(1:im,jm)=avpoles(1:im,2)
        
        call poleavg(IM,JM,JSTA,JEND,SMALL,cosltemp(1,jsta),SPVAL,avtemp(1,jsta))

        if(jsta== 1) absv(ista:iend, 1)=avtemp(ista:iend, 1)
        if(jend==jm) absv(ista:iend,jm)=avtemp(ista:iend,jm)
    
        deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)

      ELSE !(MODELNAME == 'GFS' .or. global)

      IF (GRIDTYPE == 'B')THEN
        CALL EXCH(VWND)
        CALL EXCH(UWND)
      ENDIF
     
      CALL DVDXDUDY(UWND,VWND)

      IF(GRIDTYPE == 'A')THEN
!$omp parallel do  private(i,j,jmt2,tphi,r2dx,r2dy,dvdx,dudy,uavg)
        DO J=JSTA_M,JEND_M
          JMT2 = JM/2+1
          TPHI = (J-JMT2)*(DYVAL/gdsdegr)*DTR
          DO I=ISTA_M,IEND_M
            IF(DDVDX(I,J)<SPVAL.AND.DDUDY(I,J)<SPVAL.AND. &
               UUAVG(I,J)<SPVAL.AND.UWND(I,J)<SPVAL.AND.  &
     &         UWND(I,J+1)<SPVAL.AND.UWND(I,J-1)<SPVAL) THEN
              DVDX   = DDVDX(I,J)
              DUDY   = DDUDY(I,J)
              UAVG   = UUAVG(I,J) 
!  is there a (f+tan(phi)/erad)*u term?
              IF(MODELNAME  == 'RAPR' .OR. MODELNAME  == 'FV3R') then
                 ABSV(I,J) = DVDX - DUDY + F(I,J)   ! for run RAP over north pole      
              else
                 ABSV(I,J) = DVDX - DUDY + F(I,J) + UAVG*TAN(GDLAT(I,J)*DTR)/ERAD  ! not sure about this???
              endif
            END IF
          END DO
        END DO

      ELSE IF (GRIDTYPE == 'E')THEN
       allocate(ihw(JSTA_2L:JEND_2U), IHE(JSTA_2L:JEND_2U))
!$omp  parallel do private(j)
        DO J=JSTA_2L,JEND_2U
          IHW(J) = -MOD(J,2)
          IHE(J) = IHW(J)+1
        ENDDO
!$omp parallel do  private(i,j,jmt2,tphi,r2dx,r2dy,dvdx,dudy,uavg)
        DO J=JSTA_M,JEND_M
          JMT2 = JM/2+1
          TPHI = (J-JMT2)*(DYVAL/1000.)*DTR
          TPHI = (J-JMT2)*(DYVAL/gdsdegr)*DTR
          DO I=ISTA_M,IEND_M
            IF(VWND(I+IHE(J),J) < SPVAL.AND.VWND(I+IHW(J),J) < SPVAL .AND.   &
     &         UWND(I,J+1) < SPVAL     .AND.UWND(I,J-1) < SPVAL) THEN
              DVDX   = DDVDX(I,J)
              DUDY   = DDUDY(I,J)
              UAVG   = UUAVG(I,J)
!  is there a (f+tan(phi)/erad)*u term?
              ABSV(I,J) = DVDX - DUDY + F(I,J) + UAVG*TAN(TPHI)/ERAD 
            END IF
          END DO
        END DO
       deallocate(ihw, IHE)
      ELSE IF (GRIDTYPE == 'B')THEN
!        CALL EXCH(VWND)      !done before dvdxdudy() Jesse 20200520
        DO J=JSTA_M,JEND_M
          JMT2 = JM/2+1
          TPHI = (J-JMT2)*(DYVAL/gdsdegr)*DTR
          DO I=ISTA_M,IEND_M         
            if(VWND(I,  J)==SPVAL .or. VWND(I,  J-1)==SPVAL .or. &
               VWND(I-1,J)==SPVAL .or. VWND(I-1,J-1)==SPVAL .or. &
               UWND(I,  J)==SPVAL .or. UWND(I-1,J)==SPVAL .or. &
               UWND(I,J-1)==SPVAL .or. UWND(I-1,J-1)==SPVAL) cycle
              DVDX   = DDVDX(I,J)
              DUDY   = DDUDY(I,J)
              UAVG   = UUAVG(I,J)
!  is there a (f+tan(phi)/erad)*u term?
           ABSV(I,J) = DVDX - DUDY + F(I,J) + UAVG*TAN(TPHI)/ERAD 
          END DO
        END DO 
      END IF 
      END IF
!     
!     END OF ROUTINE.
!     
      RETURN
      END

!> CALDIV computes divergence.
!>    
!> For GFS, this routine copmutes the horizontal divergence
!> using 2nd-order centered scheme on a lat-lon grid     
!>
!> @param[in] UWND U wind (m/s) mass-points.
!> @param[in] VWND V wind (m/s) mass-points.
!> @param[out] DIV divergence (1/s) mass-points.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2016-05-05 | Sajal Kar | Modified CALVORT to compute divergence from wind components
!> 2016-07-22 | S Moorthi | Modified polar divergence calculation
!>
!> @author Sajal Kar W/NP2 @date 2016-05-05

!-----------------------------------------------------------------------------
!> caldiv() computes divergence.
!
      SUBROUTINE CALDIV(UWND,VWND,DIV)
      use masks,        only: gdlat, gdlon
      use params_mod,   only: d00, dtr, small, erad
      use ctlblk_mod,   only: jsta_2l, jend_2u, spval, modelname, global, &
                              jsta, jend, im, jm, jsta_m, jend_m, lm,     &
                              ista, iend, ista_m, iend_m, ista_2l, iend_2u
      use gridspec_mod, only: gridtype

      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u,lm), intent(in)    :: UWND,VWND
      REAL, dimension(ista:iend,jsta:jend,lm),       intent(inout) :: DIV
      REAL, dimension(IM,2)         :: GLATPOLES, COSLPOLES, UPOLES, VPOLES, DIVPOLES
      REAL, dimension(IM,JSTA:JEND) :: COSLTEMP, DIVTEMP
!
      real,    allocatable ::  wrk1(:,:), wrk2(:,:), wrk3(:,:), cosl(:,:)
      INTEGER, allocatable ::  IHE(:),IHW(:), IE(:),IW(:)
!
      real                 :: dnpole, dspole, tem
      integer I,J,ip1,im1,ii,iir,iil,jj,imb2, l
!     
!***************************************************************************
!     START CALDIV HERE.
!     
!     LOOP TO COMPUTE DIVERGENCE FROM WINDS.
!     
      CALL EXCH(GDLAT(ISTA_2L,JSTA_2L))
      CALL EXCH(GDLON(ISTA_2L,JSTA_2L))

      allocate (wrk1(ista:iend,jsta:jend), wrk2(ista:iend,jsta:jend),          &
     &          wrk3(ista:iend,jsta:jend), cosl(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(iw(im),ie(im))

      imb2 = im/2
!$omp  parallel do private(i)
      do i=ista,iend
        ie(i) = i+1
        iw(i) = i-1
      enddo
!      iw(1)  = im
!      ie(im) = 1


!$omp  parallel do private(i,j,ip1,im1)
      DO J=JSTA,JEND
        do i=ista,iend
          ip1 = ie(i)
          im1 = iw(i)
          cosl(i,j) = cos(gdlat(i,j)*dtr)
          IF(cosl(i,j) >= SMALL) then
            wrk1(i,j) = 1.0 / (ERAD*cosl(i,j))
          else
            wrk1(i,j) = 0.
          end if    
          if(i == im .or. i == 1) then
            wrk2(i,j) = 1.0 / ((360.+GDLON(ip1,J)-GDLON(im1,J))*DTR) !1/dlam
          else
            wrk2(i,j) = 1.0 / ((GDLON(ip1,J)-GDLON(im1,J))*DTR)      !1/dlam
          end if
        enddo
      ENDDO

      CALL EXCH(cosl)
      CALL FULLPOLE(cosl,coslpoles)
      CALL FULLPOLE(gdlat(ista_2l:iend_2u,jsta_2l:jend_2u),glatpoles)
       
!$omp  parallel do private(i,j,ii)
      DO J=JSTA,JEND
        if (j == 1) then
          if(gdlat(ista,j) > 0.) then ! count from north to south
            do i=ista,iend
              ii = i + imb2
              if (ii > im) ii = ii - im
          !   wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J+1)-GDLAT(II,J))*DTR) !1/dphi
              wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J+1)-GLATPOLES(II,1))*DTR) !1/dphi
            enddo
          else ! count from south to north
            do i=ista,iend
              ii = i + imb2
              if (ii > im) ii = ii - im
          !   wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J+1)+GDLAT(II,J))*DTR) !1/dphi
              wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J+1)+GLATPOLES(II,1))*DTR) !1/dphi
            enddo
          end if      
        elseif (j == JM) then
          if(gdlat(ista,j) < 0.) then ! count from north to south
            do i=ista,iend
              ii = i + imb2
              if (ii > im) ii = ii - im
          !   wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J-1)+GDLAT(II,J))*DTR)
              wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J-1)+GLATPOLES(II,2))*DTR)
            enddo
          else ! count from south to north
            do i=ista,iend
              ii = i + imb2
              if (ii > im) ii = ii - im
          !   wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J-1)-GDLAT(II,J))*DTR)
              wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J-1)-GLATPOLES(II,2))*DTR)
            enddo
          end if  
        else
          do i=ista,iend
            wrk3(i,j) = 1.0 / ((GDLAT(I,J-1)-GDLAT(I,J+1))*DTR) !1/dphi
          enddo
        endif
      enddo  
      
      do l=1,lm
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            DIV(I,J,l) = SPVAL
          ENDDO
        ENDDO

        CALL EXCH(VWND(ista_2l,jsta_2l,l))
        CALL EXCH(UWND(ista_2l,jsta_2l,l))

        CALL FULLPOLE(VWND(ista_2l:iend_2u,jsta_2l:jend_2u,l),VPOLES)
        CALL FULLPOLE(UWND(ista_2l:iend_2u,jsta_2l:jend_2u,l),UPOLES)

!$omp  parallel do private(i,j,ip1,im1,ii,jj)
        DO J=JSTA,JEND
          IF(J == 1) then                          ! Near North pole
            if(gdlat(ista,j) > 0.) then ! count from north to south
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  DIV(I,J,l) = ((UWND(ip1,J,l)-UWND(im1,J,l))*wrk2(i,j)           &
     !&                    !  -  (VWND(II,J,l)*COSL(II,J)                          &
     &                       -  (VPOLES(II,1)*COSLPOLEs(II,1)                          &
     &                       +   VWND(I,J+1,l)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)
                enddo
!--
              ELSE                             !North pole point, compute at j=2
                jj = 2
                do i=ista,iend
                  ip1 = ie(i)
                  im1 = iw(i)
                  DIV(I,J,l) = ((UWND(ip1,jj,l)-UWND(im1,jj,l))*wrk2(i,jj)         &
     &                       +  (VWND(I,J,l)*COSL(I,J)                             &
                             -   VWND(I,jj+1,l)*COSL(I,jj+1))*wrk3(i,jj)) * wrk1(i,jj)
                enddo
!--
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  DIV(I,J,l) = ((UWND(ip1,J,l)-UWND(im1,J,l))*wrk2(i,j)           &
     !&                    !  +  (VWND(II,J,l)*COSL(II,J)                          &
     &                       +  (VPOLES(II,1)*COSLPOLES(II,1)                          &
     &                       +   VWND(I,J+1,l)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)
                enddo
!--
              ELSE                             !North pole point, compute at j=2
                jj = 2
                do i=ista,iend
                  ip1 = ie(i)
                  im1 = iw(i)
                  DIV(I,J,l) = ((UWND(ip1,jj,l)-UWND(im1,jj,l))*wrk2(i,jj)         &
     &                       -  (VWND(I,J,l)*COSL(I,J)                             &
                             -   VWND(I,jj+1,l)*COSL(I,jj+1))*wrk3(i,jj)) * wrk1(i,jj)
                enddo
              ENDIF
            endif
          ELSE IF(J == JM) THEN                    ! Near South pole
            if(gdlat(ista,j) < 0.) then ! count from north to south
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  DIV(I,J,l) = ((UWND(ip1,J,l)-UWND(im1,J,l))*wrk2(i,j)          &
     &                       +  (VWND(I,J-1,l)*COSL(I,J-1)                       &
     !&                    !  +   VWND(II,J,l)*COSL(II,J))*wrk3(i,j)) * wrk1(i,j)
     &                       +   VPOLES(II,2)*COSLPOLES(II,2))*wrk3(i,j)) * wrk1(i,j)
                enddo
!--
              ELSE                              !South pole point,compute at jm-1
                jj = jm-1
                do i=ista,iend
                  ip1 = ie(i)
                  im1 = iw(i)
                  DIV(I,J,l) = ((UWND(ip1,JJ,l)-UWND(im1,JJ,l))*wrk2(i,jj)       &
     &                       +  (VWND(I,jj-1,l)*COSL(I,Jj-1)                     &
     &                       -   VWND(I,J,l)*COSL(I,J))*wrk3(i,jj)) * wrk1(i,jj)

                enddo
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  DIV(I,J,l) = ((UWND(ip1,J,l)-UWND(im1,J,l))*wrk2(i,j)          &
     &                       -  (VWND(I,J-1,l)*COSL(I,J-1)                       &
     !&                    !  +   VWND(II,J,l)*COSL(II,J))*wrk3(i,j)) * wrk1(i,j)
     &                       +   VPOLES(II,2)*COSLPOLES(II,2))*wrk3(i,j)) * wrk1(i,j)
                enddo
!--
              ELSE                              !South pole point,compute at jm-1
                jj = jm-1
                do i=ista,iend
                  ip1 = ie(i)
                  im1 = iw(i)
                  DIV(I,J,l) = ((UWND(ip1,JJ,l)-UWND(im1,JJ,l))*wrk2(i,jj)       &
     &                       -  (VWND(I,jj-1,l)*COSL(I,Jj-1)                     &
     &                       -   VWND(I,J,l)*COSL(I,J))*wrk3(i,jj)) * wrk1(i,jj)

                enddo
              ENDIF
            endif
          ELSE
            DO I=ISTA,IEND
              ip1 = ie(i)
              im1 = iw(i)
              DIV(I,J,l) = ((UWND(ip1,J,l)-UWND(im1,J,l))*wrk2(i,j)           &
     &                   +  (VWND(I,J-1,l)*COSL(I,J-1)                        &
                         -   VWND(I,J+1,l)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)
            ENDDO
          ENDIF
        ENDDO                               ! end of J loop

! GFS use lon avg as one scaler value for pole point
!        call poleavg(IM,JM,JSTA,JEND,SMALL,COSL(1,jsta),SPVAL,DIV(1,jsta,l))

        call exch(div(ista_2l:iend_2u,jsta_2l:jend_2u,l))
        call fullpole(div(ista_2l:iend_2u,jsta_2l:jend_2u,l),divpoles)       

        COSLTEMP=SPVAL
        IF(JSTA== 1) COSLTEMP(1:IM, 1)=COSLPOLES(1:IM,1)
        IF(JEND==JM) COSLTEMP(1:IM,JM)=COSLPOLES(1:IM,2)
        DIVTEMP=SPVAL
        IF(JSTA== 1) DIVTEMP(1:IM, 1)=DIVPOLES(1:IM,1)
        IF(JEND==JM) DIVTEMP(1:IM,JM)=DIVPOLES(1:IM,2)
                
        call poleavg(IM,JM,JSTA,JEND,SMALL,COSLTEMP(1:IM,JSTA:JEND)  &
                    ,SPVAL,DIVTEMP(1:IM,JSTA:JEND))

        IF(JSTA== 1) DIV(ISTA:IEND, 1,L)=DIVTEMP(ISTA:IEND, 1)
        IF(JEND==JM) DIV(ISTA:IEND,JM,L)=DIVTEMP(ISTA:IEND,JM)

      enddo                        ! end of l looop
!--
      deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)
     

      END SUBROUTINE CALDIV

!------------------------------------------------------------------------
!> CALGRADPS computes gradients of a scalar field PS or LNPS.
!>
!> For GFS, this routine computes horizontal gradients of PS or LNPS.
!> Using 2nd-order centered scheme on a lat-lon grid.
!>
!> @param[in] PS Surface pressure (Pa) mass-points.
!> @param[out] PSX Zonal gradient of PS at mass-points.
!> @param[out] PSY Meridional gradient of PS at mass-points.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2016-05-05 | Sajal Kar | Reduced from CALVORT to zonal and meridional gradients of given surface pressure PS, or LNPS
!>
!> @author Sajal Kar W/NP2 @date 2016-05-05
!>
      SUBROUTINE CALGRADPS(PS,PSX,PSY)

      use masks,        only: gdlat, gdlon
      use params_mod,   only: dtr, d00, small, erad
      use ctlblk_mod,   only: jsta_2l, jend_2u, spval, modelname, global, &
                              jsta, jend, im, jm, jsta_m, jend_m,         &
                              ista, iend, ista_m, iend_m, ista_2l, iend_2u

      use gridspec_mod, only: gridtype

      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in)    :: PS
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(inout) :: PSX,PSY 
!
      real,    allocatable ::  wrk1(:,:), wrk2(:,:), wrk3(:,:), cosl(:,:)
      INTEGER, allocatable ::  IHE(:),IHW(:), IE(:),IW(:)
!
      integer I,J,ip1,im1,ii,iir,iil,jj,imb2
!     
!***************************************************************************
!     START CALGRADPS HERE.
!     
!     LOOP TO COMPUTE ZONAL AND MERIDIONAL GRADIENTS OF PS OR LNPS
!     
!sk06162016   DO J=JSTA_2L,JEND_2U
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          PSX(I,J) = SPVAL
          PSY(I,J) = SPVAL
!sk       PSX(I,J) = D00
!sk       PSY(I,J) = D00
        ENDDO
      ENDDO

      CALL EXCH(PS)

!     IF (MODELNAME == 'GFS' .or. global) THEN
        CALL EXCH(GDLAT(ISTA_2L,JSTA_2L))
        CALL EXCH(GDLON(ISTA_2L,JSTA_2L))

        allocate (wrk1(ista:iend,jsta:jend), wrk2(ista:iend,jsta:jend),          &
     &            wrk3(ista:iend,jsta:jend), cosl(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(iw(im),ie(im))

        imb2 = im/2
!$omp  parallel do private(i)
        do i=ista,iend
          ie(i) = i+1
          iw(i) = i-1
        enddo
!        iw(1)  = im
!        ie(im) = 1


!$omp  parallel do private(i,j,ip1,im1)
        DO J=JSTA,JEND
          do i=ista,iend
            ip1 = ie(i)
            im1 = iw(i)
            cosl(i,j) = cos(gdlat(i,j)*dtr)
            if(cosl(i,j) >= SMALL) then
              wrk1(i,j) = 1.0 / (ERAD*cosl(i,j))
            else
              wrk1(i,j) = 0.
            end if    
            if(i == im .or. i == 1) then
              wrk2(i,j) = 1.0 / ((360.+GDLON(ip1,J)-GDLON(im1,J))*DTR) !1/dlam
            else
              wrk2(i,j) = 1.0 / ((GDLON(ip1,J)-GDLON(im1,J))*DTR)      !1/dlam
            end if
          enddo
        ENDDO

        CALL EXCH(cosl)
       
!$omp  parallel do private(i,j,ii)
        DO J=JSTA,JEND
          if (j == 1) then
            if(gdlat(ista,j) > 0.) then ! count from north to south
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J+1)-GDLAT(II,J))*DTR) !1/dphi
              enddo
            else ! count from south to north
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J+1)+GDLAT(II,J))*DTR) !1/dphi
              enddo
            end if      
          elseif (j == JM) then
            if(gdlat(ista,j) < 0.) then ! count from north to south
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J-1)+GDLAT(II,J))*DTR)
              enddo
            else ! count from south to north
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J-1)-GDLAT(II,J))*DTR)
              enddo
            end if  
          else
            do i=ista,iend
              wrk3(i,j) = 1.0 / ((GDLAT(I,J-1)-GDLAT(I,J+1))*DTR) !1/dphi
            enddo
          endif
        ENDDO  

!$omp  parallel do private(i,j,ip1,im1,ii,jj)
        DO J=JSTA,JEND
          IF(J == 1) then                            ! Near North pole
            if(gdlat(ista,j) > 0.) then ! count from north to south
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSX(I,J) = (PS(ip1,J)-PS(im1,J))*wrk2(i,j)*wrk1(i,j)
                  PSY(I,J) = (PS(II,J)-PS(I,J+1))*wrk3(i,j)/ERAD 
                enddo
              ELSE                             !North pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSX(I,J) = (PS(ip1,jj)-PS(im1,jj))*wrk2(i,jj)*wrk1(i,jj)
                  PSY(I,J) = (PS(I,J)-PS(I,jj+1))*wrk3(i,jj)/ERAD
                enddo
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSX(I,J) = (PS(ip1,J)-PS(im1,J))*wrk2(i,j)*wrk1(i,j)
                  PSY(I,J) = - (PS(II,J)-PS(I,J+1))*wrk3(i,j)/ERAD
                enddo
              ELSE                             !North pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSX(I,J) = (PS(ip1,jj)-PS(im1,jj))*wrk2(i,jj)*wrk1(i,jj)
                  PSY(I,J) = - (PS(I,J)-PS(I,jj+1))*wrk3(i,jj)/ERAD
                enddo
              ENDIF
            endif
          ELSE IF(J == JM) THEN                      ! Near South pole
            if(gdlat(ista,j) < 0.) then ! count from north to south
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSX(I,J) = (PS(ip1,J)-PS(im1,J))*wrk2(i,j)*wrk1(i,j)
                  PSY(I,J) = (PS(I,J-1)-PS(II,J))*wrk3(i,j)/ERAD
                enddo
              ELSE                              !South pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSX(I,J) = (PS(ip1,JJ)-PS(im1,JJ))*wrk2(i,jj)*wrk1(i,jj)
                  PSY(I,J) = (PS(I,jj-1)-PS(I,J))*wrk3(i,jj)/ERAD
                enddo
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSX(I,J) = (PS(ip1,J)-PS(im1,J))*wrk2(i,j)*wrk1(i,j)
                  PSY(I,J) = - (PS(I,J-1)-PS(II,J))*wrk3(i,j)/ERAD
                enddo
              ELSE                              !South pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSX(I,J) = (PS(ip1,JJ)-PS(im1,JJ))*wrk2(i,jj)*wrk1(i,jj)
                  PSY(I,J) = - (PS(I,jj-1)-PS(I,J))*wrk3(i,jj)/ERAD
                enddo
              ENDIF
            endif
          ELSE
            DO I=ISTA,IEND
              ip1 = ie(i)
              im1 = iw(i)
              PSX(I,J)   = (PS(ip1,J)-PS(im1,J))*wrk2(i,j)*wrk1(i,j)
              PSY(I,J)   = (PS(I,J-1)-PS(I,J+1))*wrk3(i,j)/ERAD
            ENDDO
          END IF
!
        ENDDO                               ! end of J loop

        deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)
     
!     END IF 

      END SUBROUTINE CALGRADPS

!> calslr_roebber() computes snow solid-liquid-ratio slr using the Roebber algorithm.
!>
!> Obtained the code and data from WPC. WPC's SLR products include SLR computed from
!> GFS and NAM, SLR climotology, and averaged SLR. UPP computes SLR for GFS and RRFS. 
!> SLR climatology is not used in UPP calculation but the data is saved in fix directory 
!> for reference. Breadboard coefficients are included in this module to enhance the 
!> performance. Original Breadboard coefficients files are also saved in fix directory.
!> 
!> @param[in] tprs real Temperature on pressure levels.
!> @param[in] rhprs real Relative humidity on pressure levels.
!> @param[out] slr real Solid snow to liquid ratio.
!> 
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2022-07-11 | Jesse Meng | Initial
!> 2023-01-06 | Jesse Meng ! Import Breadboard coefficients into module
!>
!> @author Jesse Meng @date 2022-07-11

      SUBROUTINE CALSLR_ROEBBER(tprs,rhprs,slr)

      use masks,   only: lmh        
      use vrbls2d, only: slp, avgprec_cont, u10, v10, pshltr, tshltr, qshltr
      use vrbls3d, only: T, Q, PMID, PINT
      use ctlblk_mod, only: ista, iend, jsta, jend, &
                            ista_2l, iend_2u, jsta_2l, jend_2u, &
                            IM, JM, LM, LSM, SPL, MODELNAME, spval, me, idat
      use params_mod, only: CAPA, H1, H100
      use grib2_module, only: read_grib2_sngle

      implicit none

      real,dimension(ista_2l:iend_2u,jsta_2l:jend_2u,lsm),intent(in)  :: tprs
      real,dimension(ista_2l:iend_2u,jsta_2l:jend_2u,lsm),intent(in)  :: rhprs
      real,dimension(ista_2l:iend_2u,jsta_2l:jend_2u),    intent(out) :: slr !slr=snod/weasd=1000./sndens

! local variables
 
      character*256 :: climoFile
      logical file_exists
      integer :: ntot, height
      real,dimension(im,jm) :: CLIMO
      real,dimension(ista:iend,jsta:jend)    :: CLIMOSUB

      real,dimension(ista:iend,jsta:jend)    :: P1D,T1D,Q1D,RH1D
      real,dimension(ista:iend,jsta:jend)    :: T2M,RH2M

      type all_grids
           real :: grid
           real :: sigma
      end type all_grids

      real prob1, prob2, prob3
      real,dimension(0:14), parameter        :: sig = &
       (/0.0, 1.0, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, & 
                          0.8, 0.75, 0.7, 0.65, 0.6, 0.5, 0.4/)
      real,dimension(12), parameter          :: mf = &
       (/1.0, 0.67, 0.33, 0.0, -0.33, -0.67, -1.00, -0.67, -0.33, 0.0, 0.33, 0.67/)
      integer, dimension(0:37), parameter    :: levels = &
       (/2, 1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 725, 700, &
                  675, 650, 625, 600, 575, 550, 525, 500, 475, 450, 425, 400, &
                  375, 350, 325, 300, 275, 250, 225, 200, 175, 150, 125, 100/)

      real,dimension(0:14) :: tm, rhm

      real,dimension(0:30), parameter        :: co1 = &
       (/0.0, -.2926, .0070, -.0099, .0358, .0356, .0353, .0333, .0291, &
         .0235, .0169, .0060, -.0009, -.0052, -.0079, -.0093,&
        -.0116, -.0137, .0030, .0033, -.0005, -.0024, -.0023,&
        -.0021, -.0007, .0013, .0023, .0024, .0012, .0002, -.0010/)

      real,dimension(0:30), parameter        :: co2 = &
       (/0.0, -9.7961, .0099, -.0222, -.0036, -.0012, .0010, .0018, .0018,&
         .0011, -.0001, -.0016, -.0026, -.0021, -.0015, -.0010,&
         -.0008, -.0017, .0238, .0213, .0253, .0232, .0183, .0127,&
          .0041, -.0063, -.0088, -.0062, -.0029, .0002, .0019/)

      real,dimension(0:30), parameter        :: co3 = &
       (/0.0, 5.0037, -0.0097, -.0130, -.0170, -.0158, -.0141, -.0097,&
        -.0034, .0032, .0104, .0200, .0248, .0273, .0280, .0276,&
        .0285, .0308, -.0036, -.0042, -.0013, .0011, .0014, .0023,&
        .0011, -.0004, -.0022, -.0030, -.0033, -.0031, -.0019/)         

      real,dimension(0:30), parameter        :: co4 = &
       (/0.0, -5.0141, .0172, -.0267, .0015, .0026, .0033, .0015, -.0007,&
        -.0030, -.0063, -.0079, -.0074, -.0055, -.0035, -.0015,&
        -.0038, -.0093, .0052, .0059, .0019, -.0022, -.0077, -.0102,&
        -.0109, -.0077, .0014, .0160, .0217, .0219, .0190/)

      real,dimension(0:30), parameter        :: co5 = &
       (/0.0, -5.2807, -.0240, .0228, .0067, .0019, -.0010, -.0003, .0012,&
         .0027, .0056, .0067, .0067, .0034, .0005, -.0026, -.0039,&
         -.0033, -.0225, -.0152, -.0157, -.0094, .0049, .0138,&
         .0269, .0388, .0334, .0147, .0018, -.0066, -.0112/)

      real,dimension(0:30), parameter        :: co6 = &
       (/0.0, -2.2663, .0983, .3666, .0100, .0062, .0020, -.0008, -.0036,&
         -.0052, -.0074, -.0086, -.0072, -.0057, -.0040, -.0011,&
         .0006, .0014, .0012, -.0005, -.0019, .0003, -.0007, -.0008,&
         .0022, .0005, -.0016, -.0052, -.0024, .0008, .0037/)

      type(all_grids), dimension(ista:iend,jsta:jend,0:lsm) :: tmpk_grids, rh_grids
      integer,         dimension(ista:iend,jsta:jend,0:lsm) :: tmpk_levels, rh_levels

      real,dimension(ista:iend,jsta:jend)    :: hprob,mprob,lprob
      real,dimension(ista:iend,jsta:jend)    :: slrgrid, slrgrid2
      real,dimension(ista:iend,jsta:jend)    :: psfc,pres,qpf,swnd,prp

      character*20 nswFileName
      real :: psurf,p,sgw,sg1,sg2,dtds,rhds
      real :: f1,f2,f3,f4,f5,f6
      real :: p1,p2,p3
      real :: hprob_tot
      real :: mprob_tot
      real :: lprob_tot

      integer :: i,j,k,ks,L,LL,imo,iday
!
!***************************************************************************
!
! day and month of the year

      imo = idat(1)
      iday= idat(2)

! climatology
! currently not used, snoden climatology files saved in fix directory
!
!      climoFile='climo_snoden'
!      ntot=im*jm
!      CLIMO = spval
!      CLIMOSUB = spval
!      INQUIRE(FILE=climoFile, EXIST=file_exists)
!      if(file_exists) then
!         print*,trim(climoFile),' FOUND'
!         call read_grib2_sngle(climoFile,ntot,height,CLIMO)
!         do j=jsta,jend
!         do i=ista,iend
!            if(CLIMO(i,j).gt.0 .and. CLIMO(i,j).lt.1000) CLIMOSUB(i,j)=1000./CLIMO(i,j)
!         endif
!         end do
!         end do
!      else
!         print*,trim(climoFile),' NOT FOUND'
!      endif !if(file_exist)

! surface variables

!$omp parallel do private(i,j)
      DO J=JSTA,JEND
      DO I=ISTA,IEND
         PSFC(I,J)=PINT(I,J,NINT(LMH(I,J))+1)
         PRES(I,J)=SLP(I,J)
         QPF(I,J)=AVGPREC_CONT(I,J)*3600.*3.
         SWND(I,J)=SPVAL
         IF(U10(I,J)/=SPVAL .AND. V10(I,J)/=SPVAL) &
           SWND(I,J)=SQRT(U10(I,J)*U10(I,J)+V10(I,J)*V10(I,J))
      END DO
      END DO

! T2M and RH2M

!$omp parallel do private(i,j)
      DO J=JSTA,JEND
      DO I=ISTA,IEND
         IF(MODELNAME=='RAPR')THEN
            P1D(I,J) = PMID(I,J,NINT(LMH(I,J)))
            T1D(I,J) = T(I,J,NINT(LMH(I,J))) 
         ELSE
            P1D(I,J) = PINT(I,J,LM+1)*EXP(-0.068283/TSHLTR(I,J))
            T1D(I,J) = TSHLTR(I,J)*(PSHLTR(I,J)*1.E-5)**CAPA
         ENDIF
         Q1D(I,J) = QSHLTR(I,J)
         T2M(I,J) = T1D(I,J)
      ENDDO
      ENDDO
     
      CALL CALRH(P1D,T1D,Q1D,RH1D)

!$omp parallel do private(i,j)
      DO J=JSTA,JEND
      DO I=ISTA,IEND
         if(qshltr(i,j) /= spval)then
            RH2M(I,J) = min(H100,max(H1,RH1D(I,J)*100.))
         else
            RH2M(I,J) = spval 
         endif
      ENDDO
      ENDDO

!$omp parallel do private(i,j)
      do j=jsta,jend
      do i=ista,iend
         tmpk_grids(i,j,0)%grid=T2M(I,J)-273.15
         tmpk_levels(i,j,0)=pres(i,j)
         rh_grids(i,j,0)%grid=RH2M(I,J)
         rh_levels(i,j,0)=pres(i,j)
      end do
      end do

! T and RH all pressure levels

      DO L=1,LSM
        LL=LSM-L+1
!!!$omp parallel do private(i,j,ll)
      do j=jsta,jend
      do i=ista,iend
         tmpk_grids(i,j,LL)%grid=tprs(I,J,L)-273.15
         tmpk_levels(i,j,LL)=SPL(L)
         rh_grids(i,j,LL)%grid=rhprs(I,J,L)
         rh_levels(i,j,LL)=SPL(L)
      end do
      end do
      END DO

! convert to sigma

      tmpk_grids(:,:,0)%sigma = 1.0
      rh_grids(:,:,0)%sigma = 1.0

      DO L=1,LSM
        LL=LSM-L+1
!!!$omp parallel do private(i,j,ll)
        do j=jsta,jend
        do i=ista,iend
           if(pres(i,j) == spval) then
              tmpk_grids(i,j,LL)%sigma=spval
              rh_grids(i,j,LL)%sigma=spval
           else
              tmpk_grids(i,j,LL)%sigma=tmpk_levels(i,j,LL)/pres(i,j)     
              rh_grids(i,j,LL)%sigma=rh_levels(i,j,LL)/pres(i,j)
              prp(i,j)=pres(i,j)/psfc(i,j)
              prp(i,j)=prp(i,j)*100000./psfc(i,j)
           endif
        end do
        end do
      END DO

! main slr i/j loop starts

      do j=jsta,jend
      do i=ista,iend
         tm=spval
         rhm=spval
         slr(i,j)=spval
         slrgrid(i,j)=spval
         slrgrid2(i,j)=spval
         hprob(i,j)=spval
         mprob(i,j)=spval
         lprob(i,j)=spval

      if(pres(i,j)/=spval .and. qpf(i,j)/=spval .and. swnd(i,j)/=spval) then

! Interpolate T and RH to the 14 sigma levels      

      do ks=1,14
         psurf=pres(i,j)
         sgw=sig(ks)
         p=prp(i,j)
         do LL=0,LSM-1
           if(LL==0) then
              sg1 = psurf/psurf
           else
              sg1 = tmpk_levels(i,j,LL)/psurf
           endif  
           sg2 = tmpk_levels(i,j,LL+1)/psurf
           
           if(sg1 == sgw) then
              tm(ks) = tmpk_grids(i,j,LL)%grid
              rhm(ks)=   rh_grids(i,j,LL)%grid
           elseif (sg2 == sgw) then
              tm(ks) = tmpk_grids(i,j,LL+1)%grid
              rhm(ks)=   rh_grids(i,j,LL+1)%grid
           elseif ((sgw < sg1) .and. (sgw > sg2)) then
              dtds = (tmpk_grids(i,j,LL+1)%grid - tmpk_grids(i,j,LL)%grid)/(sg2-sg1)
              tm(ks) = ((sgw - sg1) * dtds) + tmpk_grids(i,j,LL)%grid
              rhds = (rh_grids(i,j,LL+1)%grid - rh_grids(i,j,LL)%grid)/(sg2-sg1)
              rhm(ks)= ((sgw - sg1) * rhds) + rh_grids(i,j,LL)%grid
           endif
         end do
      end do !loop ks

! Have surface wind, QPF, and temp/RH on the 14 sigma levels.
! Convert these data to the factors using regression equations

      f1 = co1(1)+co1(2)*qpf(i,j)+co1(3)*swnd(i,j)+co1(4)*tm(1)+co1(5)*tm(2)+co1(6)*tm(3)+ &
           co1(7)*tm(4)+co1(8)*tm(5)+co1(9)*tm(6)+co1(10)*tm(7)+co1(11)*tm(8)+ &
           co1(12)*tm(9)+co1(13)*tm(10)+co1(14)*tm(11)+co1(15)*tm(12)+co1(16)*tm(13)+ &
           co1(17)*tm(14)+co1(18)*rhm(1)+co1(19)*rhm(2)+co1(20)*rhm(3)+co1(21)*rhm(4)+ &
           co1(22)*rhm(5)+co1(23)*rhm(6)+co1(24)*rhm(7)+co1(25)*rhm(8)+co1(26)*rhm(9)+ &
           co1(27)*rhm(10)+co1(28)*rhm(11)+co1(29)*rhm(12)+co1(30)*rhm(13)

      f2 = co2(1)+co2(2)*qpf(i,j)+co2(3)*swnd(i,j)+co2(4)*tm(1)+co2(5)*tm(2)+co2(6)*tm(3)+ &
           co2(7)*tm(4)+co2(8)*tm(5)+co2(9)*tm(6)+co2(10)*tm(7)+co2(11)*tm(8)+ &
           co2(12)*tm(9)+co2(13)*tm(10)+co2(14)*tm(11)+co2(15)*tm(12)+co2(16)*tm(13)+ &
           co2(17)*tm(14)+co2(18)*rhm(1)+co2(19)*rhm(2)+co2(20)*rhm(3)+co2(21)*rhm(4)+ &
           co2(22)*rhm(5)+co2(23)*rhm(6)+co2(24)*rhm(7)+co2(25)*rhm(8)+co2(26)*rhm(9)+ &
           co2(27)*rhm(10)+co2(28)*rhm(11)+co2(29)*rhm(12)+co2(30)*rhm(13)

      f3 = co3(1)+co3(2)*qpf(i,j)+co3(3)*swnd(i,j)+co3(4)*tm(1)+co3(5)*tm(2)+co3(6)*tm(3)+ &
           co3(7)*tm(4)+co3(8)*tm(5)+co3(9)*tm(6)+co3(10)*tm(7)+co3(11)*tm(8)+ &
           co3(12)*tm(9)+co3(13)*tm(10)+co3(14)*tm(11)+co3(15)*tm(12)+co3(16)*tm(13)+ &
           co3(17)*tm(14)+co3(18)*rhm(1)+co3(19)*rhm(2)+co3(20)*rhm(3)+co3(21)*rhm(4)+ &
           co3(22)*rhm(5)+co3(23)*rhm(6)+co3(24)*rhm(7)+co3(25)*rhm(8)+co3(26)*rhm(9)+ &
           co3(27)*rhm(10)+co3(28)*rhm(11)+co3(29)*rhm(12)+co3(30)*rhm(13)

      f4 = co4(1)+co4(2)*qpf(i,j)+co4(3)*swnd(i,j)+co4(4)*tm(1)+co4(5)*tm(2)+co4(6)*tm(3)+ &
           co4(7)*tm(4)+co4(8)*tm(5)+co4(9)*tm(6)+co4(10)*tm(7)+co4(11)*tm(8)+ &
           co4(12)*tm(9)+co4(13)*tm(10)+co4(14)*tm(11)+co4(15)*tm(12)+co4(16)*tm(13)+ &
           co4(17)*tm(14)+co4(18)*rhm(1)+co4(19)*rhm(2)+co4(20)*rhm(3)+co4(21)*rhm(4)+ &
           co4(22)*rhm(5)+co4(23)*rhm(6)+co4(24)*rhm(7)+co4(25)*rhm(8)+co4(26)*rhm(9)+ &
           co4(27)*rhm(10)+co4(28)*rhm(11)+co4(29)*rhm(12)+co4(30)*rhm(13)

      f5 = co5(1)+co5(2)*qpf(i,j)+co5(3)*swnd(i,j)+co5(4)*tm(1)+co5(5)*tm(2)+co5(6)*tm(3)+ &
           co5(7)*tm(4)+co5(8)*tm(5)+co5(9)*tm(6)+co5(10)*tm(7)+co5(11)*tm(8)+ &
           co5(12)*tm(9)+co5(13)*tm(10)+co5(14)*tm(11)+co5(15)*tm(12)+co5(16)*tm(13)+ &
           co5(17)*tm(14)+co5(18)*rhm(1)+co5(19)*rhm(2)+co5(20)*rhm(3)+co5(21)*rhm(4)+ &
           co5(22)*rhm(5)+co5(23)*rhm(6)+co5(24)*rhm(7)+co5(25)*rhm(8)+co5(26)*rhm(9)+ &
           co5(27)*rhm(10)+co5(28)*rhm(11)+co5(29)*rhm(12)+co5(30)*rhm(13)
      
      f6 = co6(1)+co6(2)*qpf(i,j)+co6(3)*swnd(i,j)+co6(4)*tm(1)+co6(5)*tm(2)+co6(6)*tm(3)+ &
           co6(7)*tm(4)+co6(8)*tm(5)+co6(9)*tm(6)+co6(10)*tm(7)+co6(11)*tm(8)+ &
           co6(12)*tm(9)+co6(13)*tm(10)+co6(14)*tm(11)+co6(15)*tm(12)+co6(16)*tm(13)+ &
           co6(17)*tm(14)+co6(18)*rhm(1)+co6(19)*rhm(2)+co6(20)*rhm(3)+co6(21)*rhm(4)+ &
           co6(22)*rhm(5)+co6(23)*rhm(6)+co6(24)*rhm(7)+co6(25)*rhm(8)+co6(26)*rhm(9)+ &
           co6(27)*rhm(10)+co6(28)*rhm(11)+co6(29)*rhm(12)+co6(30)*rhm(13)

      hprob_tot = 0.
      mprob_tot = 0.
      lprob_tot = 0.
      do k=1,10
         if(k==1) then
            nswFileName='Breadboard1.nsw'
            call breadboard1_main(nswFileName,mf(imo),f1,f2,f3,f4,f5,f6,p1,p2,p3)
         elseif(k==2) then
            nswFileName='Breadboard2.nsw'
            call breadboard1_main(nswFileName,mf(imo),f1,f2,f3,f4,f5,f6,p1,p2,p3)
         elseif(k==3) then
            nswFileName='Breadboard3.nsw'
            call breadboard1_main(nswFileName,mf(imo),f1,f2,f3,f4,f5,f6,p1,p2,p3)
         elseif(k==4) then
            nswFileName='Breadboard4.nsw'
            call breadboard1_main(nswFileName,mf(imo),f1,f2,f3,f4,f5,f6,p1,p2,p3)
         elseif(k==5) then
            nswFileName='Breadboard5.nsw'
            call breadboard1_main(nswFileName,mf(imo),f1,f2,f3,f4,f5,f6,p1,p2,p3)
         elseif(k==6) then
            nswFileName='Breadboard6.nsw'
            call breadboard6_main(nswFileName,mf(imo),f1,f2,f3,f4,f5,f6,p1,p2,p3)
         elseif(k==7) then
            nswFileName='Breadboard7.nsw'
            call breadboard6_main(nswFileName,mf(imo),f1,f2,f3,f4,f5,f6,p1,p2,p3)
         elseif(k==8) then
            nswFileName='Breadboard8.nsw'
            call breadboard6_main(nswFileName,mf(imo),f1,f2,f3,f4,f5,f6,p1,p2,p3)
         elseif(k==9) then
            nswFileName='Breadboard9.nsw'
            call breadboard6_main(nswFileName,mf(imo),f1,f2,f3,f4,f5,f6,p1,p2,p3)
         elseif(k==10) then
            nswFileName='Breadboard10.nsw'
            call breadboard6_main(nswFileName,mf(imo),f1,f2,f3,f4,f5,f6,p1,p2,p3)
         endif
         hprob_tot = hprob_tot+p1
         mprob_tot = mprob_tot+p2
         lprob_tot = lprob_tot+p3
      enddo
      hprob(i,j) = hprob_tot/10.
      mprob(i,j) = mprob_tot/10.
      lprob(i,j) = lprob_tot/10.

      if(hprob(i,j) > mprob(i,j) .and. hprob(i,j) > lprob(i,j)) then
         slrgrid(i,j) = 8.0
      elseif(mprob(i,j) >= hprob(i,j) .and. mprob(i,j) >= lprob(i,j)) then
         slrgrid(i,j) = 13.0
      elseif(lprob(i,j) > hprob(i,j) .and. lprob(i,j) > mprob(i,j)) then   
          if(lprob(i,j) < .67) then
             slrgrid(i,j) = 18.0
          else 
             slrgrid(i,j) = 27.0
          endif
      endif

!     Weighted SLR

      if(lprob(i,j) < .67) then
         slrgrid2(i,j) = hprob(i,j)*8.0+mprob(i,j)*13.0+lprob(i,j)*18.0
         slrgrid2(i,j) = slrgrid2(i,j)*p/(hprob(i,j)+mprob(i,j)+lprob(i,j))
      else
         slrgrid2(i,j) = hprob(i,j)*8.0+mprob(i,j)*13.0+lprob(i,j)*27.0
         slrgrid2(i,j) = slrgrid2(i,j)*p/(hprob(i,j)+mprob(i,j)+lprob(i,j))
      endif
 
!      slr(i,j) = climosub(i,j)
!      slr(i,j) = slrgrid(i,j)
      slr(i,j) = slrgrid2(i,j)
      slr(i,j) = max(1.,min(25.,slr(i,j)))

      endif !if(pres(i,j), qpf(i,j), swnd(i,j) /= spval)
      enddo
      enddo

! main slr i/j loop ends

      END SUBROUTINE CALSLR_ROEBBER
!
!-------------------------------------------------------------------------------------
!> @brief breadboard1_main() _____ ???
      SUBROUTINE breadboard1_main(nswFileName,mf,f1,f2,f3,f4,f5,f6,p1,p2,p3)

      implicit none

      character*20 nswFileName
      real mf, f1, f2, f3, f4, f5, f6
      real p1, p2, p3

      real f(7)

      real inputFile(2,7)
      real inputAxon(7)
      real hidden1Axon(40)
      real outputAxon(3)
      real hidden1Synapse(7,40)
      real outputSynapse(40,3)
      real activeOutputProbe(2,3)

      real fgrid1(40), fgrid2(3), fgridsum

      integer i,j
!
      f(1) = mf
      f(2) = f1
      f(3) = f2
      f(4) = f3
      f(5) = f4
      f(6) = f5
      f(7) = f6

! Read nsw file and load weights

      inputFile(1,:)=1.
      inputFile(2,:)=0.
      inputAxon=0.
      hidden1Axon=0.
      outputAxon=0.
      hidden1Synapse=1.
      outputSynapse=1.
      activeOutputProbe(1,:)=1.
      activeOutputProbe(2,:)=0.

      if(trim(nswFileName)=='Breadboard1.nsw') then
        call Breadboard1(inputFile,hidden1Axon,hidden1Synapse,outputSynapse)
      elseif(trim(nswFileName)=='Breadboard2.nsw') then
        call Breadboard2(inputFile,hidden1Axon,hidden1Synapse,outputSynapse)
      elseif(trim(nswFileName)=='Breadboard3.nsw') then
        call Breadboard3(inputFile,hidden1Axon,hidden1Synapse,outputSynapse)
      elseif(trim(nswFileName)=='Breadboard4.nsw') then
        call Breadboard4(inputFile,hidden1Axon,hidden1Synapse,outputSynapse)
      elseif(trim(nswFileName)=='Breadboard5.nsw') then
        call Breadboard5(inputFile,hidden1Axon,hidden1Synapse,outputSynapse)
      endif

      if(activeOutputProbe(1,1)==1.) then
      do j=1,3
         activeOutputProbe(1,j)=8.999999761581421e-001
         activeOutputProbe(2,j)=5.000000074505806e-002
      enddo
      endif

! Run Network

      do j=1,7
         inputAxon(j) = inputFile(1,j) * f(j) + inputFile(2,j)
      enddo

      fgrid1=0.
!$omp parallel do private(i,j)
      do j=1,40
         do i=1,7
            fgrid1(j) = fgrid1(j) + hidden1Synapse(i,j) * inputAxon(i) 
         enddo
         fgrid1(j) = fgrid1(j) + hidden1Axon(j)
         fgrid1(j) = (exp(fgrid1(j))-exp(-fgrid1(j)))/(exp(fgrid1(j))+exp(-fgrid1(j)))
      enddo

      fgrid2=0.
      fgridsum=0.
      do j=1,3
         do i=1,40
            fgrid2(j) = fgrid2(j) + outputSynapse(i,j) * fgrid1(i)
         enddo
         fgrid2(j) = fgrid2(j) + outputAxon(j)
         fgrid2(j) = exp(fgrid2(j))
         fgridsum = fgridsum + fgrid2(j)
      enddo
      do j=1,3
         fgrid2(j) = fgrid2(j) / fgridsum
!         fgrid2(j) = activeOutputProbe(1,j) * fgrid2(j) + activeOutputProbe(2,j)  
      enddo

      p1 = fgrid2(1)
      p2 = fgrid2(2)
      p3 = fgrid2(3)

      END SUBROUTINE breadboard1_main
!
!-------------------------------------------------------------------------------------
!
      SUBROUTINE breadboard6_main(nswFileName,mf,f1,f2,f3,f4,f5,f6,p1,p2,p3)

      implicit none

      character*20 nswFileName
      real mf, f1, f2, f3, f4, f5, f6
      real p1, p2, p3

      real f(7)

      real inputFile(2,7)
      real inputAxon(7)
      real hidden1Axon(7)
      real hidden2Axon(4)
      real outputAxon(3)
      real hidden1Synapse(7,7)
      real hidden2Synapse(7,4)
      real outputSynapse(4,3)
      real activeOutputProbe(2,3)

      real fgrid1(7), fgrid2(4), fgrid3(3), fgridsum

      integer i,j
!
      f(1) = mf
      f(2) = f1
      f(3) = f2
      f(4) = f3
      f(5) = f4
      f(6) = f5
      f(7) = f6
!
      inputFile(1,:)=1.
      inputFile(2,:)=0.
      inputAxon=0.
      hidden1Axon=0.
      hidden2Axon=0.
      outputAxon=0.
      hidden1Synapse=1.
      hidden2Synapse=1.
      outputSynapse=1.
      activeOutputProbe(1,:)=1.
      activeOutputProbe(2,:)=0.

      if(trim(nswFileName)=='Breadboard6.nsw') then
        call Breadboard6(inputFile,hidden1Axon,hidden2Axon,&
             hidden1Synapse,hidden2Synapse,outputSynapse)
      elseif(trim(nswFileName)=='Breadboard7.nsw') then
        call Breadboard7(inputFile,hidden1Axon,hidden2Axon,&
             hidden1Synapse,hidden2Synapse,outputSynapse)
      elseif(trim(nswFileName)=='Breadboard8.nsw') then
        call Breadboard8(inputFile,hidden1Axon,hidden2Axon,&
             hidden1Synapse,hidden2Synapse,outputSynapse)
      elseif(trim(nswFileName)=='Breadboard9.nsw') then
        call Breadboard9(inputFile,hidden1Axon,hidden2Axon,&
             hidden1Synapse,hidden2Synapse,outputSynapse)
      elseif(trim(nswFileName)=='Breadboard10.nsw') then
        call Breadboard10(inputFile,hidden1Axon,hidden2Axon,&
             hidden1Synapse,hidden2Synapse,outputSynapse)
      endif

      if(activeOutputProbe(1,1)==1.) then
      do j=1,3
         activeOutputProbe(1,j)=8.999999761581421e-001 
         activeOutputProbe(2,j)=5.000000074505806e-002
      enddo
      endif

! Run Network

      do j=1,7
         inputAxon(j) = inputFile(1,j) * f(j) + inputFile(2,j)
      enddo

      fgrid1=0.
!$omp parallel do private(i,j)
      do j=1,7
         do i=1,7
            fgrid1(j) = fgrid1(j) + hidden1Synapse(i,j) * inputAxon(i) 
         enddo
         fgrid1(j) = fgrid1(j) + hidden1Axon(j)
         fgrid1(j) = (exp(fgrid1(j))-exp(-fgrid1(j)))/(exp(fgrid1(j))+exp(-fgrid1(j)))
      enddo

      fgrid2=0.
!$omp parallel do private(i,j)
      do j=1,4
         do i=1,7
            fgrid2(j) = fgrid2(j) + hidden2Synapse(i,j) * fgrid1(i)
         enddo
         fgrid2(j) = fgrid2(j) + hidden2Axon(j)
         fgrid2(j) = (exp(fgrid2(j))-exp(-fgrid2(j)))/(exp(fgrid2(j))+exp(-fgrid2(j)))
      enddo

      fgrid3=0.
      fgridsum=0.
      do j=1,3
         do i=1,4
            fgrid3(j) = fgrid3(j) + outputSynapse(i,j) * fgrid2(i)
         enddo
         fgrid3(j) = fgrid3(j) + outputAxon(j)
         fgrid3(j) = exp(fgrid3(j))
         fgridsum = fgridsum + fgrid3(j)
      enddo
      do j=1,3
         fgrid3(j) = fgrid3(j) / fgridsum
!         fgrid3(j) = activeOutputProbe(1,j) * fgrid3(j) + activeOutputProbe(2,j)
      enddo

      p1 = fgrid3(1)
      p2 = fgrid3(2)
      p3 = fgrid3(3)
      
      END SUBROUTINE breadboard6_main
!
!------------------------------------------------------------------------------
!
      SUBROUTINE Breadboard1(inputFile,hidden1Axon,hidden1Synapse,outputSynapse)
 
      implicit none
 
      real inputFile(2,7)
      real hidden1Axon(40)
      real hidden1Synapse(7,40)
      real outputSynapse(40,3)
 
      inputFile = reshape((/ &
         1.077844262123108E+00, -1.778443008661270E-01,&
         2.295625507831573E-01,  6.163756549358368E-02,&
         2.081887423992157E-01,  6.210270524024963E-01,&
         3.646677434444427E-01,  1.214343756437302E-01,&
         2.430133521556854E-01,  3.004860281944275E-01,&
         1.935067623853683E-01,  4.185551702976227E-01,&
         1.962280571460724E-01, -4.804643988609314E-01 &
        /), shape(inputFile))
 
      hidden1Axon = &
      (/-1.150484442710876E+00, -1.461968779563904E+00,  1.349107265472412E+00,  6.686212420463562E-01,&
        -8.486616015434265E-01, -1.908162593841553E+00, -1.514992356300354E+00, -1.632351636886597E+00,&
        -1.794843912124634E+00,  1.354879975318909E+00,  1.389558911323547E+00,  1.464104652404785E+00,&
         1.896052122116089E+00,  1.401677846908569E+00,  1.436681509017944E+00, -1.590880393981934E+00,&
        -1.070504426956177E+00,  2.047163248062134E+00,  1.564107656478882E+00,  1.298712372779846E+00,&
        -1.316817998886108E+00, -1.253177642822266E+00, -1.392926216125488E+00,  7.356406450271606E-01,&
         1.594561100006104E+00, -1.532955884933472E+00, -1.021214842796326E+00,  1.341110348701477E+00,&
         6.124811172485352E-01,  1.415654063224792E+00, -8.509962558746338E-01,  1.753035664558411E+00,&
         6.275475621223450E-01,  1.482843875885010E+00,  1.326028347015381E+00,  1.641556143760681E+00,&
         1.339018464088440E+00, -1.374068379402161E+00, -1.220067739486694E+00,  1.714797854423523E+00/)
 
      hidden1Synapse = reshape((/ &
        -4.612099826335907E-01, -3.177818655967712E-01, -2.800635099411011E-01, -6.984808295965195E-02,&
         6.583837419748306E-02, -5.769817233085632E-01,  3.955098092556000E-01, -1.624705344438553E-01,&
        -2.889076173305511E-01, -9.411631226539612E-01, -5.058886408805847E-01, -3.110982775688171E-01,&
        -3.723000884056091E-01,  8.419776558876038E-01,  2.598794996738434E-01, -1.364605724811554E-01,&
         9.416468143463135E-01, -4.025689139962196E-02,  4.176554381847382E-01,  1.196979433298111E-01,&
        -3.846398293972015E-01, -1.414917409420013E-01, -2.344214916229248E+00, -3.556166291236877E-01,&
        -7.762963771820068E-01, -1.243659138679504E+00,  4.907984733581543E-01, -1.891903519630432E+00,&
        -5.802390575408936E-01, -5.546363592147827E-01, -4.520095884799957E-01, -2.473797500133514E-01,&
        -7.757837772369385E-01, -5.350160598754883E-01,  1.817676275968552E-01, -1.932217180728912E-01,&
         5.944451093673706E-01, -6.568105518817902E-02, -1.562235504388809E-01,  4.926294833421707E-02,&
        -6.931540369987488E-01,  7.082754969596863E-01, -3.878217563033104E-02,  5.063381195068359E-01,&
        -7.642447352409363E-01, -2.539043128490448E-01, -4.328470230102539E-01, -4.773662984371185E-01,&
         6.699458956718445E-01, -1.670347154140472E-01,  6.986252665519714E-01, -6.806275844573975E-01,&
         1.059119179844856E-01,  5.320579931139946E-02, -4.806780517101288E-01,  7.601988911628723E-01,&
        -1.864496916532516E-01, -3.076690435409546E-01, -6.505665779113770E-01,  7.355872541666031E-02,&
        -4.033335149288177E-01, -2.168276757001877E-01,  5.354191064834595E-01,  2.991014420986176E-01,&
         4.275756180286407E-01,  6.465418934822083E-01, -1.401910781860352E-01,  5.381527543067932E-01,&
         9.247279167175293E-01, -3.687029778957367E-01,  6.354923844337463E-01, -1.423558890819550E-01,&
         9.430686831474304E-01,  1.187003701925278E-01,  5.426434278488159E-01,  7.573884129524231E-01,&
         3.361994773149490E-02,  3.300542756915092E-02, -4.439333379268646E-01,  5.953744649887085E-01,&
         3.412617444992065E-01,  1.421828866004944E-01,  5.224847793579102E-01, -8.267756700515747E-01,&
         5.009499788284302E-01,  2.736742198467255E-01,  8.603093624114990E-01,  9.373022615909576E-02,&
         1.714528501033783E-01,  9.114132076501846E-02, -1.638108491897583E-01,  5.879403948783875E-01,&
         5.585592240095139E-03,  8.149939179420471E-01, -1.340572237968445E-01,  3.880683779716492E-01,&
         3.857498764991760E-01, -8.105239868164062E-01,  5.239543914794922E-01,  7.420576363801956E-02,&
         7.694411277770996E-01, -3.954831138253212E-02,  5.615213513374329E-01,  4.560695886611938E-01,&
        -5.006425976753235E-01, -4.725854694843292E-01,  5.887325108051300E-02, -3.199687898159027E-01,&
        -5.229111015796661E-02, -6.034490466117859E-01, -8.414428234100342E-01,  1.826022863388062E-01,&
        -6.954011321067810E-01, -5.277091860771179E-01, -9.834931492805481E-01, -2.964940369129181E-01,&
         1.752081327140331E-02, -2.412298470735550E-01,  5.861807465553284E-01,  3.650662600994110E-01,&
        -1.846716850996017E-01,  3.277707397937775E-01,  1.213769540190697E-01,  1.398152709007263E-01,&
         1.624975651502609E-01, -7.172397375106812E-01, -4.065496101975441E-02, -1.131931394338608E-01,&
         7.050336003303528E-01,  3.453079611063004E-02,  5.642467141151428E-01,  7.171959280967712E-01,&
        -3.295499980449677E-01,  5.192958116531372E-01,  7.558688521385193E-01,  6.164067387580872E-01,&
        -1.597565859556198E-01,  1.512383669614792E-01,  5.231227278709412E-01, -2.199545800685883E-01,&
        -3.987313508987427E-01, -9.710572957992554E-01, -4.689137935638428E-01, -4.037811756134033E-01,&
        -4.528387784957886E-01, -4.784810543060303E-01,  1.759306043386459E-01,  7.449938654899597E-01,&
         1.120681285858154E+00, -5.609570741653442E-01,  1.393345594406128E+00,  1.374282408505678E-02,&
        -2.458193153142929E-01,  1.237058401107788E+00, -4.854794219136238E-02, -6.664386391639709E-01,&
        -8.786886334419250E-01, -3.208510577678680E-01, -4.315690398216248E-01, -5.186472535133362E-01,&
        -2.117208093404770E-01,  8.998587727546692E-02,  7.763032317161560E-01,  1.078992128372192E+00,&
         3.667660653591156E-01,  5.805531740188599E-01,  1.517073512077332E-01,  9.344519972801208E-01,&
         3.396262824535370E-01,  2.450248003005981E-01,  9.134629368782043E-01,  7.127542048692703E-02,&
        -1.287281513214111E-01,  3.953699469566345E-01, -4.097535610198975E-01, -5.983641743659973E-01,&
         4.500437378883362E-01, -8.147508651018143E-02, -7.916551083326340E-02, -1.505649089813232E-01,&
        -1.703914403915405E-01,  1.294612526893616E+00, -4.859757721424103E-01, -1.034098416566849E-01,&
        -6.859915256500244E-01,  4.521823674440384E-02,  3.100419938564301E-01, -9.373775720596313E-01,&
         5.841451883316040E-01,  7.020491957664490E-01, -1.681403964757919E-01,  6.397892832756042E-01,&
         1.168430075049400E-01,  4.124156236648560E-01,  5.404921174049377E-01, -3.311195969581604E-01,&
        -3.494578003883362E-01,  1.379718184471130E+00,  2.731607258319855E-01,  5.512273311614990E-01,&
         2.997024357318878E-01,  3.475511670112610E-01,  6.777516603469849E-01,  1.471205204725266E-01,&
         1.011002138257027E-01,  8.974244594573975E-01,  8.688372373580933E-02,  4.767233729362488E-01,&
         9.785303473472595E-01, -2.200428694486618E-01, -6.173372268676758E-01, -8.801369071006775E-01,&
        -1.111719012260437E+00, -3.223371803760529E-01, -6.491173505783081E-01, -3.894545435905457E-01,&
        -2.843862473964691E-01,  7.331426739692688E-01, -3.287445753812790E-02, -5.741032306104898E-03,&
         6.212961673736572E-01,  3.749484941363335E-02,  6.244438700377941E-03, -6.228777766227722E-01,&
        -4.667133837938309E-02,  2.016694307327271E+00,  2.834755480289459E-01,  6.229624748229980E-01,&
         6.552317738533020E-01, -9.771268069744110E-02,  7.506207823753357E-01,  6.942567825317383E-01,&
        -1.662521809339523E-01,  3.003259599208832E-01, -2.531996071338654E-01,  2.399661689996719E-01,&
         5.109554529190063E-01, -7.031706571578979E-01,  2.836774885654449E-01,  4.888223409652710E-01,&
         1.384589523077011E-01, -3.524579405784607E-01, -2.050135582685471E-01,  1.160808563232422E+00,&
        -4.008938968181610E-01,  1.656456440687180E-01, -5.116114616394043E-01,  8.800522685050964E-01,&
         6.836380064487457E-02, -5.902936309576035E-02,  5.672354102134705E-01, -7.219299674034119E-01,&
         3.463289514183998E-02, -1.044675827026367E+00, -8.341925591230392E-02, -3.036961853504181E-01,&
        -5.605638027191162E-01,  5.722484588623047E-01, -1.604338049888611E+00, -5.696258544921875E-01,&
        -2.531512081623077E-01, -4.675458073616028E-01, -6.486019492149353E-01, -2.437075823545456E-01,&
        -2.898264527320862E-01,  3.836293518543243E-01,  4.061043560504913E-01,  3.909072279930115E-01,&
        -8.113911151885986E-01,  1.260317683219910E+00, -3.924282491207123E-01,  3.586370870471001E-02,&
         7.703443765640259E-01,  6.714462637901306E-01, -4.909946396946907E-02,  3.536651730537415E-01,&
         1.900762617588043E-01,  3.638494014739990E-01,  2.248179465532303E-01, -6.255846619606018E-01 &
        /), shape(hidden1Synapse))
 
      outputSynapse = reshape((/ &
        -4.825605154037476E-01, -1.119017243385315E+00,  5.116804838180542E-01, -6.694142222404480E-01,&
        -5.718530416488647E-01, -7.233589291572571E-01, -8.200560212135315E-01, -6.121573448181152E-01,&
        -1.034205436706543E+00,  1.015549778938293E+00,  1.183975338935852E+00,  5.342597365379333E-01,&
         1.186208128929138E+00,  7.657266259193420E-01,  9.990772604942322E-01, -1.051267385482788E+00,&
        -7.288008332252502E-01,  9.447612762451172E-01,  6.943449974060059E-01,  5.248318314552307E-01,&
        -1.042970657348633E+00, -4.857340827584267E-04, -8.969252705574036E-01,  5.206210613250732E-01,&
         7.825390100479126E-01, -3.175100982189178E-01, -7.697273492813110E-01,  3.042222857475281E-01,&
         7.400255203247070E-01,  1.082547545433044E+00, -1.058874249458313E+00,  3.296852707862854E-01,&
         9.955985546112061E-01,  7.361931800842285E-01,  8.618848919868469E-01,  7.109408378601074E-01,&
         1.148022636771202E-01, -6.803723573684692E-01, -4.462003335356712E-02,  7.384030222892761E-01,&
        -2.215545326471329E-01, -8.702403903007507E-01,  8.234908580780029E-01,  6.819239258766174E-01,&
        -4.687527120113373E-01, -6.959788203239441E-01, -6.105158329010010E-01, -7.225347757339478E-01,&
        -7.860832810401917E-01,  5.608791112899780E-01,  9.937217235565186E-01,  6.797130703926086E-01,&
         8.231667280197144E-01,  1.115462303161621E+00,  5.290299654006958E-01, -4.602016210556030E-01,&
        -5.394889116287231E-01,  1.053055644035339E+00,  9.533493518829346E-01,  8.694807887077332E-01,&
        -4.802323281764984E-01, -1.070514082908630E+00, -8.236010670661926E-01,  7.932062149047852E-01,&
         1.111655592918396E+00, -1.025945305824280E+00, -2.268178462982178E-01,  6.432797908782959E-01,&
         2.442117929458618E-01,  7.986634969711304E-01, -3.561095297336578E-01,  1.058865070343018E+00,&
         6.459046602249146E-01,  4.042869210243225E-01,  2.976681292057037E-02,  1.033244490623474E+00,&
         9.110773205757141E-01, -6.528528332710266E-01, -8.971995115280151E-01,  1.046785235404968E+00,&
        -5.487565994262695E-01, -1.033755183219910E+00,  5.164890289306641E-01,  1.108534336090088E+00,&
        -2.507440149784088E-01, -1.150385260581970E+00, -1.040475010871887E+00, -1.114320755004883E+00,&
        -9.695596694946289E-01,  9.147439599037170E-01,  3.035557866096497E-01,  1.044997453689575E+00,&
         1.059857130050659E+00,  7.304399013519287E-01,  1.102171182632446E+00, -9.304327964782715E-01,&
        -5.997116565704346E-01,  1.120478868484497E+00,  6.444569826126099E-01,  2.137384265661240E-01,&
        -4.117920994758606E-01, -1.000458717346191E+00, -2.041520774364471E-01, -1.859422773122787E-01,&
         3.711319267749786E-01, -9.141649603843689E-01, -7.499164938926697E-01,  9.900025129318237E-01,&
        -2.189985066652298E-01,  8.942219614982605E-01, -3.195305764675140E-01,  6.445295810699463E-01,&
        -2.110123336315155E-01,  9.763143658638000E-01,  8.833498954772949E-01,  1.071311354637146E+00,&
         1.134591102600098E+00, -4.175429344177246E-01, -6.000540852546692E-01,  7.281569838523865E-01 &
        /), shape(outputSynapse))
 
      END SUBROUTINE Breadboard1         
!
!------------------------------------------------------------------------------
!
      SUBROUTINE Breadboard2(inputFile,hidden1Axon,hidden1Synapse,outputSynapse)
 
      implicit none
 
      real inputFile(2,7)
      real hidden1Axon(40)
      real hidden1Synapse(7,40)
      real outputSynapse(40,3)
 
      inputFile = reshape((/ &
         1.077844262123108E+00, -1.778443008661270E-01,&
         2.188449800014496E-01,  1.674167998135090E-02,&
         1.868382692337036E-01,  6.490761637687683E-01,&
         3.361344337463379E-01,  4.151264205574989E-02,&
         2.621995508670807E-01,  2.531536519527435E-01,&
         1.944894641637802E-01,  3.221717774868011E-01,&
         3.179650008678436E-01, -2.033386379480362E-01 &
        /), shape(inputFile))
 
      hidden1Axon = &
      (/-9.235364943742752E-02, -5.511198639869690E-01,  1.012191653251648E+00, -1.148184835910797E-01,&
        -8.399781584739685E-01, -4.726789295673370E-01,  7.570160627365112E-01, -3.985013365745544E-01,&
         1.164000511169434E+00,  2.212587594985962E-01,  9.570528268814087E-01, -1.504407286643982E+00,&
        -1.262813359498978E-01,  9.741528630256653E-01,  2.278975844383240E-01, -3.282702267169952E-01,&
         1.716251969337463E-01,  4.979004263877869E-01,  6.414948105812073E-01, -2.775986790657043E-01,&
        -6.721665859222412E-01,  7.226511836051941E-01, -1.020949006080627E+00, -9.638186097145081E-01,&
         4.050622135400772E-02, -8.287806510925293E-01, -2.900803685188293E-01,  1.004199028015137E+00,&
        -1.221053838729858E+00, -5.891714692115784E-01, -6.459002494812012E-01,  8.228222727775574E-01,&
         1.921370178461075E-01,  1.575044542551041E-01, -9.904603362083435E-01,  1.186665743589401E-01,&
         1.871918141841888E-01, -6.121324300765991E-01,  1.056765243411064E-01, -5.654883384704590E-01/)
 
      hidden1Synapse = reshape((/ &
        -5.215738341212273E-02,  6.958795785903931E-01, -3.700282871723175E-01,  4.440588057041168E-01,&
        -9.248711913824081E-02,  9.709199517965317E-02,  1.255098581314087E-01, -1.359838247299194E-01,&
         3.981630802154541E-01, -4.047442674636841E-01, -5.247595906257629E-01, -5.138890147209167E-01,&
         2.293408364057541E-01,  5.139534473419189E-01,  2.035804986953735E-01,  3.003124892711639E-01,&
        -2.340262830257416E-01,  3.037432730197906E-01,  4.666079878807068E-01,  3.753643631935120E-01,&
        -5.292671918869019E-02,  3.674933612346649E-01,  3.854512274265289E-01,  1.749511361122131E-01,&
         1.320011764764786E-01,  2.418431788682938E-01,  1.245125234127045E-01, -2.677426636219025E-01,&
         3.884479776024818E-02, -1.385747641324997E-01, -3.117613494396210E-01,  3.016934990882874E-01,&
        -2.856997251510620E-01, -4.838032424449921E-01,  4.488031566143036E-01, -3.862534165382385E-01,&
         2.520084977149963E-01, -6.066129356622696E-02, -2.037643343210220E-01, -9.749407321214676E-02,&
         1.909288167953491E-01, -2.689029574394226E-01,  8.022837042808533E-01,  4.543448388576508E-01,&
         1.268999278545380E-01,  2.794430553913116E-01,  4.331161379814148E-01, -1.717756092548370E-01,&
        -5.167780518531799E-01,  6.074145808815956E-02,  2.141399085521698E-01, -3.536535203456879E-01,&
        -2.548796236515045E-01, -4.349331259727478E-01,  3.771509276703000E-03,  1.351494044065475E-01,&
         8.080910146236420E-02, -2.638687789440155E-01,  1.792310923337936E-01, -5.317723155021667E-01,&
         6.300682574510574E-02,  1.391339004039764E-01, -6.581404209136963E-01,  1.574699729681015E-01,&
        -5.979638695716858E-01, -6.864693760871887E-01, -6.892689466476440E-01, -1.189238503575325E-01,&
        -1.904999166727066E-01, -4.838389158248901E-01,  4.585682973265648E-02,  3.201213181018829E-01,&
         5.204908251762390E-01, -3.531241044402122E-02,  4.392628967761993E-01,  4.307939708232880E-01,&
        -4.227218031883240E-02,  1.247199028730392E-01,  1.489800363779068E-01, -3.146159052848816E-01,&
         2.637389600276947E-01, -8.966535329818726E-02,  2.010040730237961E-01,  3.161593675613403E-01,&
        -8.221558481454849E-02, -4.601925909519196E-01, -3.832246661186218E-01,  2.877672016620636E-01,&
        -1.351716276258230E-02, -5.320604424923658E-03, -3.493662178516388E-02, -1.777663826942444E-01,&
        -1.865815520286560E-01,  6.387206912040710E-01, -4.405377805233002E-01,  4.452396631240845E-01,&
        -1.245370283722878E-01, -2.323225736618042E-01,  1.697962284088135E-01,  1.118463352322578E-01,&
        -2.475701570510864E-01, -3.791887685656548E-02,  5.509998202323914E-01,  1.247667223215103E-01,&
         3.189268708229065E-01, -3.584641516208649E-01,  8.915060758590698E-01,  9.720049053430557E-02,&
        -1.117252558469772E-01,  3.543806076049805E-01, -2.351483702659607E-01,  5.283502340316772E-01,&
         1.746209561824799E-01,  1.741478294134140E-01,  2.738423347473145E-01,  3.764865398406982E-01,&
         3.486587703227997E-01, -3.462808132171631E-01,  9.324266910552979E-01,  2.155355364084244E-01,&
        -5.171442404389381E-02,  6.311618685722351E-01, -1.088170856237411E-01,  4.840107262134552E-01,&
        -2.310744374990463E-01, -3.167505562305450E-01, -2.271509468555450E-01, -2.800688743591309E-01,&
         4.713648185133934E-02, -1.575807780027390E-01,  3.583298251032829E-02, -3.308865129947662E-01,&
        -2.662795484066010E-01,  1.894978582859039E-01,  7.474141567945480E-02, -1.493624746799469E-01,&
        -1.482628136873245E-01, -1.058527529239655E-01, -3.737696707248688E-01, -1.093639135360718E-01,&
        -4.270362555980682E-01,  1.249950975179672E-01, -1.971846818923950E-01,  3.135327398777008E-01,&
         4.604682624340057E-01, -4.614944458007812E-01,  4.820220768451691E-01,  3.806282877922058E-01,&
         3.629744052886963E-01,  3.986520171165466E-01, -2.283873707056046E-01,  1.246029064059258E-01,&
         3.940442204475403E-01,  2.390366494655609E-01,  8.402416110038757E-02,  3.498363792896271E-01,&
        -3.888027667999268E-01,  2.272991091012955E-01, -3.421411216259003E-01,  1.273499727249146E-01,&
         1.342627108097076E-01,  1.159043312072754E-01,  1.274240911006927E-01, -2.915177941322327E-01,&
         6.415430903434753E-01,  1.699399948120117E-01, -6.556300520896912E-01,  9.605846554040909E-02,&
         3.632318377494812E-01, -3.854629993438721E-01, -3.860571384429932E-01, -1.257066577672958E-01,&
        -1.186188161373138E-01, -1.368320286273956E-01, -2.300722897052765E-01, -4.762146174907684E-01,&
        -3.621844053268433E-01, -4.978014528751373E-02, -1.940275430679321E-01, -1.588442362844944E-02,&
        -1.519876420497894E-01,  1.312368810176849E-01,  1.862339228391647E-01,  6.462548375129700E-01,&
         5.544137358665466E-01, -3.416634351015091E-02,  9.995899349451065E-02, -6.969342380762100E-02,&
        -1.428494304418564E-01,  2.647481858730316E-01,  1.083492934703827E-01,  5.986538901925087E-02,&
        -1.576850377023220E-02,  1.962803453207016E-01,  6.334787011146545E-01, -1.408149152994156E-01,&
        -1.756295561790466E-01, -2.156554609537125E-01, -1.412229537963867E-01, -5.801249146461487E-01,&
        -5.700040608644485E-02, -3.019523918628693E-01, -1.161280944943428E-01, -3.032382726669312E-01,&
         1.140000447630882E-01, -2.648598253726959E-01, -2.016042023897171E-01, -3.181084990501404E-02,&
         7.931513339281082E-02,  5.399967432022095E-01, -4.595367014408112E-01,  9.602636098861694E-02,&
        -4.730868339538574E-01,  2.077568918466568E-01, -2.257115393877029E-01,  3.216529190540314E-01,&
         1.631081402301788E-01,  6.222640164196491E-03, -1.323710232973099E-01,  1.348871737718582E-01,&
         1.123578473925591E-01,  5.462109446525574E-01,  5.289056897163391E-01,  5.155519247055054E-01,&
         2.748569846153259E-01, -3.125837743282318E-01, -3.262098431587219E-01, -8.945185691118240E-03,&
        -4.980920553207397E-01,  5.064374208450317E-01, -1.056439951062202E-01, -3.115973472595215E-01,&
         3.343601152300835E-02, -7.157339155673981E-02,  5.459919571876526E-01,  2.175374031066895E-01,&
        -2.892075665295124E-02,  1.139620468020439E-01, -4.409461319446564E-01, -4.908669367432594E-02,&
        -2.098206430673599E-01,  3.024870157241821E-01, -3.447104394435883E-01, -2.666398882865906E-01,&
        -1.739841997623444E-01, -1.120999976992607E-01,  4.268572330474854E-01,  4.144327044487000E-01,&
         4.936498403549194E-01,  5.718982815742493E-01,  5.464938655495644E-02,  3.950506746768951E-01,&
        -1.432464718818665E-01, -8.016809076070786E-02,  5.947722792625427E-01, -1.419431418180466E-01,&
        -2.328271418809891E-01, -1.958254128694534E-01, -9.914696216583252E-03, -1.478249877691269E-01,&
         4.182004928588867E-01,  7.797469943761826E-02,  3.761124014854431E-01,  4.066407680511475E-01,&
         1.217691525816917E-01, -1.124059110879898E-01,  7.020493596792221E-02,  1.022125557065010E-01,&
        -5.025411844253540E-01, -2.482684552669525E-01, -5.819427594542503E-02, -1.587846502661705E-02,&
        -1.881837695837021E-01,  4.026338756084442E-01,  3.339109122753143E-01,  2.215891182422638E-01,&
         7.083265781402588E-01, -7.670203596353531E-02,  3.171359598636627E-01,  8.310161828994751E-01 &
        /), shape(hidden1Synapse))
 
      outputSynapse = reshape((/ &
         2.309078276157379E-01,  8.006124198436737E-02,  5.207773447036743E-01,  3.642434999346733E-02,&
        -5.444544181227684E-02, -2.300137132406235E-01,  4.965198636054993E-01, -3.590968847274780E-01,&
         1.392439752817154E-01, -2.941058278083801E-01,  6.655657291412354E-01, -4.931978881359100E-01,&
        -1.253394484519958E-01,  1.540697813034058E-01,  1.752252578735352E-01,  4.873855113983154E-01,&
         5.741749405860901E-01,  1.275441497564316E-01, -4.765471443533897E-02, -5.038099363446236E-02,&
        -8.334141224622726E-02,  5.842098593711853E-01, -4.490646719932556E-01, -5.416034907102585E-02,&
        -2.264686524868011E-01, -1.698177903890610E-01,  3.113179206848145E-01,  4.435532391071320E-01,&
        -5.240975022315979E-01,  1.108570247888565E-01,  2.321150526404381E-02,  2.374080866575241E-01,&
        -2.570592761039734E-01,  3.205819129943848E-01, -3.468126952648163E-01,  2.772298157215118E-01,&
         1.148034259676933E-01,  1.865169033408165E-03,  3.649827241897583E-01,  5.026416182518005E-01,&
        -2.502067089080811E-01, -6.028710007667542E-01, -6.978485733270645E-02,  8.656968921422958E-02,&
        -5.227651596069336E-01,  9.525942802429199E-02, -1.903700232505798E-01,  1.426358073949814E-01,&
         5.602359771728516E-01, -2.479453980922699E-01,  1.296138316392899E-01, -4.612154662609100E-01,&
        -4.198251068592072E-01,  6.053315401077271E-01, -1.160371229052544E-01, -4.044520258903503E-01,&
        -1.530461944639683E-02,  4.267008602619171E-01,  2.162231802940369E-01,  1.101492717862129E-01,&
        -9.195729345083237E-02, -3.771322593092918E-02,  3.320552408695221E-02, -4.979051947593689E-01,&
         1.581449210643768E-01, -5.021102428436279E-01,  1.184114068746567E-02,  4.836803376674652E-01,&
        -5.539562702178955E-01, -2.782657444477081E-01, -1.547775119543076E-01,  4.582551419734955E-01,&
         2.844007611274719E-01, -4.516306817531586E-01,  1.886052638292313E-02,  3.602048456668854E-01,&
         4.175081476569176E-02,  2.075715661048889E-01, -5.455711483955383E-01, -2.442489415407181E-01,&
        -2.680016458034515E-01,  2.636941149830818E-03,  4.164874255657196E-01,  8.120876550674438E-02,&
        -4.927250146865845E-01, -3.254565298557281E-01,  5.583248138427734E-01, -1.608870923519135E-01,&
         5.749610066413879E-01,  5.479150414466858E-01,  3.469662666320801E-01, -5.061987638473511E-01,&
         3.353976905345917E-01,  2.548734247684479E-01,  2.064624279737473E-01, -5.114225745201111E-01,&
        -4.629626572132111E-01, -1.936426460742950E-01,  2.327886223793030E-01, -4.583241790533066E-02,&
        -5.125665068626404E-01,  1.089363321661949E-01, -4.951449036598206E-01, -5.018569827079773E-01,&
         2.582837454974651E-02,  4.913705959916115E-02, -2.441505938768387E-01, -3.174663335084915E-02,&
        -1.644173413515091E-01, -2.947083115577698E-01, -5.097694396972656E-01,  7.136650383472443E-03,&
         1.942666023969650E-01,  1.587397605180740E-01, -4.691866040229797E-01, -4.862202703952789E-01,&
         1.432444006204605E-01, -4.405085742473602E-01,  3.072859644889832E-01, -4.172921180725098E-01 &
        /), shape(outputSynapse))
 
      END SUBROUTINE Breadboard2         
!
!------------------------------------------------------------------------------
!
      SUBROUTINE Breadboard3(inputFile,hidden1Axon,hidden1Synapse,outputSynapse)
 
      implicit none
 
      real inputFile(2,7)
      real hidden1Axon(40)
      real hidden1Synapse(7,40)
      real outputSynapse(40,3)
 
      inputFile = reshape((/ &
         1.077844262123108E+00, -1.778443008661270E-01,&
         2.442665100097656E-01,  3.212104737758636E-02,&
         2.107975035905838E-01,  6.168988943099976E-01,&
         3.646677434444427E-01,  1.214343756437302E-01,&
         2.485501170158386E-01,  2.868268489837646E-01,&
         1.976718604564667E-01,  4.469360709190369E-01,&
         3.208556175231934E-01, -2.509090602397919E-01 &
        /), shape(inputFile))
 
      hidden1Axon = &
       (/4.393131732940674E-01, -1.290386915206909E-01,  6.327351331710815E-01,  5.494017004966736E-01,&
         4.969031810760498E-01,  2.086368650197983E-01, -2.167895883321762E-01,  9.464725255966187E-01,&
         1.640024334192276E-01,  2.452306896448135E-01,  1.972979009151459E-01,  9.276027083396912E-01,&
         2.502645850181580E-01,  5.485208034515381E-01, -2.839279770851135E-01,  6.810981035232544E-01,&
        -2.170253098011017E-01, -3.821973502635956E-01,  8.861125111579895E-01, -6.720829606056213E-01,&
         2.960434183478355E-02, -3.987881243228912E-01, -1.057050973176956E-01,  6.963993310928345E-01,&
        -1.413413435220718E-01,  7.551014423370361E-01,  1.243001222610474E-02, -3.603826761245728E-01,&
         7.450697422027588E-01,  7.630060315132141E-01,  5.904716849327087E-01, -5.035977959632874E-01,&
         2.082890830934048E-03, -1.259811818599701E-01, -8.103467822074890E-01, -4.683765172958374E-01,&
        -3.666405081748962E-01, -5.880022794008255E-02, -5.269588828086853E-01, -1.594118028879166E-01/)
 
      hidden1Synapse = reshape((/ &
         2.258135080337524E-01, -8.417334407567978E-02, -6.296884268522263E-02, -1.971755474805832E-01,&
        -2.008096426725388E-01,  1.312222182750702E-01, -2.187249064445496E-01,  3.300825655460358E-01,&
        -1.458171010017395E-01, -2.447441816329956E-01,  2.373344898223877E-01, -3.369296491146088E-01,&
        -2.142974138259888E-01,  7.442125119268894E-03,  2.400149852037430E-01,  5.063241720199585E-01,&
         1.461273133754730E-01,  3.199279010295868E-01,  2.184794545173645E-01,  6.378577351570129E-01,&
         2.826454937458038E-01,  1.467282772064209E-01,  4.167218208312988E-01,  3.410821408033371E-02,&
        -1.507616639137268E-01,  1.607457697391510E-01,  1.063031926751137E-01,  4.860900044441223E-01,&
        -7.546984404325485E-02,  3.811344206333160E-01, -3.500247746706009E-02, -3.294828236103058E-01,&
        -2.355449087917805E-02,  3.319101631641388E-01,  1.341840159147978E-02, -2.975183129310608E-01,&
        -2.044427692890167E-01,  7.903610914945602E-02, -2.241216152906418E-01, -1.982768028974533E-01,&
         2.166045308113098E-01, -3.769811093807220E-01, -4.219292849302292E-02, -4.683617055416107E-01,&
         1.365721821784973E-01, -5.708352923393250E-01, -5.482509136199951E-01, -5.697317123413086E-01,&
         3.948671817779541E-01,  4.008982181549072E-01, -6.056785583496094E-01, -6.540334783494473E-03,&
        -4.144128859043121E-01, -9.239719808101654E-02,  1.977843493223190E-01, -2.407579571008682E-01,&
        -2.472878843545914E-01, -3.429937064647675E-01, -1.058190166950226E-01, -8.456809073686600E-02,&
         4.944565296173096E-01,  4.329789280891418E-01,  2.303941249847412E-01,  2.076211571693420E-01,&
         1.421037223190069E-02,  5.740813165903091E-02,  1.577541381120682E-01,  1.072699949145317E-01,&
         3.550452180206776E-03, -7.603026926517487E-02,  1.787180006504059E-01,  3.000865578651428E-01,&
        -4.790667295455933E-01, -1.263711899518967E-01, -1.886992603540421E-01, -1.971553862094879E-01,&
        -4.320513010025024E-01, -1.786982715129852E-01, -3.415124714374542E-01,  3.517304956912994E-01,&
         3.841716647148132E-01,  1.595797836780548E-01,  1.466515809297562E-01,  3.235963284969330E-01,&
         3.831133618950844E-02,  3.778985887765884E-02,  4.742037355899811E-01, -1.204959601163864E-01,&
        -6.766954064369202E-02,  4.763844013214111E-01,  2.847502529621124E-01, -2.614455521106720E-01,&
         4.211461246013641E-01,  2.459102123975754E-01, -3.291262984275818E-01,  4.159525930881500E-01,&
         1.433917880058289E-01,  5.506788492202759E-01, -4.396528601646423E-01,  3.432570993900299E-01,&
        -4.605481028556824E-01, -1.657515168190002E-01,  2.847986221313477E-01, -3.968485295772552E-01,&
         2.652311325073242E-01,  2.413431182503700E-03,  6.885899305343628E-01, -1.771224141120911E-01,&
        -2.605379931628704E-02,  1.681880354881287E-01,  4.201361536979675E-01, -2.905318737030029E-01,&
        -1.065197512507439E-01,  2.377779632806778E-01,  3.171224892139435E-01, -5.171843245625496E-02,&
         8.248845487833023E-02, -4.904226213693619E-02,  3.065647780895233E-01,  1.610077768564224E-01,&
         8.712385892868042E-01,  3.008154034614563E-01,  5.729283690452576E-01, -1.608658432960510E-01,&
        -3.810124993324280E-01,  6.462811827659607E-01, -2.662218213081360E-01, -5.297539830207825E-01,&
        -1.356185525655746E-01,  2.623566091060638E-01, -1.624718308448792E-01, -2.004417479038239E-01,&
        -3.377428650856018E-02,  3.970716595649719E-01, -1.560127288103104E-01,  4.747187346220016E-02,&
        -3.162815868854523E-01, -3.350041508674622E-01, -3.987393081188202E-01, -4.969080090522766E-01,&
        -1.142657846212387E-01, -7.119160890579224E-01,  1.153976768255234E-01, -6.001577973365784E-01,&
        -3.606468439102173E-01, -3.741255104541779E-01, -7.550917863845825E-01,  1.106901541352272E-01,&
        -1.475569456815720E-01, -2.016223073005676E-01, -2.226002812385559E-01,  2.520006597042084E-01,&
        -4.015582501888275E-01, -6.874573230743408E-01, -3.860632777214050E-01,  1.074488908052444E-01,&
        -3.594025373458862E-01, -2.556712925434113E-01,  2.491754293441772E-01, -1.749203801155090E-01,&
        -5.133146420121193E-03, -2.629097700119019E-01,  1.706630140542984E-01,  5.300921797752380E-01,&
         3.016012907028198E-01,  3.024738729000092E-01,  1.334729231894016E-02,  3.605858981609344E-01,&
        -3.797290921211243E-01,  2.125910073518753E-01, -3.324515819549561E-01, -2.657738924026489E-01,&
         8.549436926841736E-02,  2.843597829341888E-01, -1.628004312515259E-01,  4.068509638309479E-01,&
        -1.096388697624207E-01,  1.842555999755859E-01, -2.429902255535126E-01,  1.793259531259537E-01,&
         6.289024949073792E-01,  4.427114427089691E-01, -8.943214267492294E-02,  1.407862901687622E-01,&
        -4.747562706470490E-01,  1.607088744640350E-01,  2.691341638565063E-01, -1.326033025979996E-01,&
        -6.888723373413086E-02,  3.347525000572205E-01,  2.391179502010345E-01, -7.601787149906158E-02,&
         3.946174979209900E-01,  4.608300328254700E-01, -4.973608553409576E-01,  2.180006355047226E-02,&
        -2.155515551567078E-01,  4.018128812313080E-01,  5.872810482978821E-01, -2.970355451107025E-01,&
         6.164746284484863E-01, -2.832284271717072E-01, -7.214747369289398E-02,  3.505393862724304E-01,&
         3.504253327846527E-01, -3.037774860858917E-01, -3.341494500637054E-01, -2.143821418285370E-01,&
         3.230984508991241E-01, -6.691335439682007E-01, -1.196009963750839E-01,  2.609530091285706E-01,&
         6.332063078880310E-01, -2.495922595262527E-01, -1.421163380146027E-01,  4.370761811733246E-01,&
         2.344440817832947E-01, -4.770855009555817E-01, -1.213536486029625E-01, -4.947537779808044E-01,&
         2.018401175737381E-01, -3.219321966171265E-01, -1.836685538291931E-01,  6.838442683219910E-01,&
        -5.349717736244202E-01,  5.601373910903931E-01, -3.152181506156921E-01,  2.578000128269196E-01,&
         4.295753240585327E-01, -1.423847377300262E-01,  6.693964004516602E-01, -2.671292051672935E-02,&
        -2.906464338302612E-01, -6.406581997871399E-01, -5.139582753181458E-01,  2.622411847114563E-01,&
         2.534431815147400E-01, -1.518065035343170E-01, -4.292866215109825E-02,  4.628975689411163E-01,&
         1.969320774078369E-01,  4.264309704303741E-01, -4.475159347057343E-01, -5.727919340133667E-01,&
         5.388451814651489E-01, -2.982297539710999E-01, -3.593768924474716E-02, -1.298359930515289E-01,&
        -4.535509645938873E-01, -1.963836848735809E-01, -2.640297412872314E-01,  3.889253437519073E-01,&
        -2.371201291680336E-02,  5.441716909408569E-01, -3.557947278022766E-01, -1.912423074245453E-01,&
         3.168485462665558E-01, -3.096546828746796E-01,  2.481035888195038E-01,  2.293358147144318E-01,&
        -7.027690410614014E-01, -4.839945435523987E-01, -2.963027358055115E-01, -5.126427412033081E-01,&
         2.138081789016724E-01, -2.071801871061325E-01, -9.827529639005661E-02, -4.680003225803375E-01,&
        -3.230824470520020E-01, -2.535474896430969E-01,  2.779140770435333E-01, -5.119556188583374E-01,&
         1.893053054809570E-01, -5.211792513728142E-02,  4.212611019611359E-01, -5.767111182212830E-01,&
         3.436119556427002E-01,  1.560586243867874E-01, -1.338404417037964E-01,  2.465801686048508E-01 &
        /), shape(hidden1Synapse))
 
      outputSynapse = reshape((/ &
        -1.504478603601456E-01,  8.304652571678162E-02,  2.053809165954590E-01,  4.613898992538452E-01,&
         3.307471871376038E-01, -2.503668665885925E-01, -4.260648787021637E-01, -2.033478170633316E-01,&
         1.205723360180855E-01,  3.727485835552216E-01, -2.320208251476288E-01,  4.672348499298096E-01,&
        -1.567042618989944E-01,  4.181037843227386E-01, -2.018750756978989E-01,  2.649243474006653E-01,&
         2.292609065771103E-01,  2.745892405509949E-01,  2.554303109645844E-01, -3.891312777996063E-01,&
        -4.561745524406433E-01, -3.781261444091797E-01, -2.881123721599579E-01,  2.764029800891876E-01,&
         8.924255520105362E-02,  4.471623599529266E-01,  9.589984267950058E-02,  4.323486387729645E-01,&
         4.792469739913940E-01, -9.918873012065887E-02,  4.427296221256256E-01,  3.841804563999176E-01,&
         1.890532523393631E-01, -4.477364718914032E-01, -2.994475699961185E-02, -7.976207137107849E-02,&
         2.607934474945068E-01, -3.710708916187286E-01, -2.811897993087769E-01,  6.034602597355843E-02,&
         4.014556109905243E-01,  2.982565164566040E-01,  4.447779953479767E-01, -3.612459823489189E-02,&
        -2.895380258560181E-01,  2.155442684888840E-01, -3.415147066116333E-01,  4.278375506401062E-01,&
         1.896717213094234E-02, -9.841635823249817E-02,  1.671093255281448E-01,  3.151571452617645E-01,&
        -1.678100675344467E-01, -4.435905069112778E-02, -2.333792001008987E-01,  4.360995292663574E-01,&
         3.587894737720490E-01, -1.017290875315666E-01,  1.382773071527481E-01, -3.980610668659210E-01,&
        -2.268472909927368E-01, -2.996328286826611E-02,  2.546367645263672E-01,  1.532198935747147E-01,&
        -1.018586382269859E-02,  3.147244155406952E-01, -3.700032234191895E-01,  2.747226655483246E-01,&
         4.799823760986328E-01,  3.735623657703400E-01,  3.757937550544739E-01, -5.869687348604202E-02,&
         7.807171344757080E-02, -1.428240090608597E-01, -5.030028820037842E-01, -4.323083460330963E-01,&
        -2.643692195415497E-01, -4.277939200401306E-01,  3.172474205493927E-01, -4.587580561637878E-01,&
         4.488629996776581E-01, -1.273735053837299E-02,  2.275637537240982E-01,  2.276848852634430E-01,&
         1.995900124311447E-01, -1.224325075745583E-01, -1.321871429681778E-01,  4.938367307186127E-01,&
         3.713837862014771E-01,  4.943797290325165E-01, -8.973516523838043E-02,  3.630679845809937E-01,&
         3.118912279605865E-01,  3.763218820095062E-01, -2.658533453941345E-01,  5.210888572037220E-03,&
        -3.098636865615845E-01, -4.516429603099823E-01,  3.575363755226135E-01,  3.780608177185059E-01,&
         3.606519103050232E-01,  4.404914379119873E-01, -4.452764391899109E-01,  2.741447389125824E-01,&
         1.122588440775871E-01,  2.581178247928619E-01, -2.986721992492676E-01, -3.506239950656891E-01,&
        -4.466909915208817E-02,  1.343552619218826E-01, -2.677312493324280E-02, -5.070485472679138E-01,&
        -5.414816737174988E-01,  3.392856195569038E-02, -4.090670943260193E-01,  2.741051837801933E-02,&
         7.242175936698914E-02,  4.587205946445465E-01, -2.530987001955509E-02,  1.304957270622253E-02 &
        /), shape(outputSynapse))
 
      END SUBROUTINE Breadboard3         
!
!------------------------------------------------------------------------------
!
      SUBROUTINE Breadboard4(inputFile,hidden1Axon,hidden1Synapse,outputSynapse)
 
      implicit none
 
      real inputFile(2,7)
      real hidden1Axon(40)
      real hidden1Synapse(7,40)
      real outputSynapse(40,3)
 
      inputFile = reshape((/ &
         1.077844262123108E+00, -1.778443008661270E-01,&
         2.296211272478104E-01,  6.142363324761391E-02,&
         2.128665894269943E-01,  6.552034020423889E-01,&
         3.361344337463379E-01,  4.151264205574989E-02,&
         2.430133521556854E-01,  3.004860281944275E-01,&
         1.976718604564667E-01,  4.469360709190369E-01,&
         1.951007992029190E-01, -4.725341200828552E-01 &
        /), shape(inputFile))
 
      hidden1Axon = &
      (/-1.700838446617126E+00,  1.409139156341553E+00, -1.263895153999329E+00, -1.653346180915833E+00,&
        -1.753814935684204E+00,  1.510319232940674E+00, -1.652730584144592E+00,  1.968622922897339E+00,&
        -1.764715671539307E+00, -1.920537590980530E+00,  1.703584432601929E+00,  9.688673615455627E-01,&
         1.621924757957458E+00, -1.195185184478760E+00, -1.170735836029053E+00, -1.726262569427490E+00,&
         1.693020582199097E+00, -1.789734363555908E+00,  2.076834440231323E+00, -2.054785251617432E+00,&
         1.735462069511414E+00, -1.377997517585754E+00,  1.685962557792664E+00, -1.505226492881775E+00,&
         1.329061865806580E+00, -1.970339655876160E+00,  1.326048374176025E+00, -1.803932785987854E+00,&
        -1.356570959091187E+00, -7.451403737068176E-01,  1.977797389030457E+00,  1.962222456932068E+00,&
        -1.924186825752258E+00, -1.927103757858276E+00,  1.772511124610901E+00,  2.267752170562744E+00,&
         1.343345522880554E+00, -1.727791309356689E+00, -1.688525199890137E+00, -2.020093202590942E+00/)
 
      hidden1Synapse = reshape((/ &
        -3.217298686504364E-01, -1.535140275955200E-01, -9.374593496322632E-01, -3.773699328303337E-02,&
        -7.610699534416199E-01,  1.124547328799963E-03,  7.987623810768127E-01,  5.171887874603271E-01,&
         1.182283610105515E-01,  1.252476930618286E+00, -2.393243610858917E-01,  8.846385776996613E-02,&
         4.983871877193451E-01, -1.072657704353333E+00, -5.902777314186096E-01,  3.053096830844879E-01,&
        -1.245228290557861E+00, -9.408684819936752E-02, -1.261333227157593E+00,  7.626018673181534E-02,&
        -3.566111624240875E-01, -2.651087939739227E-01,  5.490935966372490E-02, -1.231116533279419E+00,&
        -3.552156984806061E-01, -4.995369017124176E-01, -1.970071047544479E-01,  6.921592950820923E-01,&
        -7.216929793357849E-01, -3.322352096438408E-02, -1.040984153747559E+00, -2.749272584915161E-01,&
        -3.936901688575745E-01, -5.485629439353943E-01,  2.315377295017242E-01,  3.925201594829559E-01,&
         2.289973348379135E-01,  9.091649055480957E-01, -2.400987595319748E-01,  2.274930775165558E-01,&
         7.657364010810852E-01, -4.531333744525909E-01, -3.045647442340851E-01, -1.612837314605713E-01,&
        -6.530205607414246E-01,  6.988145411014557E-02, -3.664937913417816E-01, -1.209497332572937E+00,&
         1.716423481702805E-01,  2.888691425323486E-01, -6.977611780166626E-01,  1.001697182655334E+00,&
        -3.773393929004669E-01, -3.817198425531387E-02,  3.071420192718506E-01, -1.018374800682068E+00,&
        -3.812201619148254E-01,  2.521711289882660E-01, -1.311386704444885E+00, -4.305998682975769E-01,&
        -2.096824795007706E-01, -6.536886692047119E-01,  9.946095943450928E-02, -8.006195425987244E-01,&
         6.314782798290253E-02, -9.162106513977051E-01,  1.249427199363708E-01, -1.967987567186356E-01,&
        -2.837883234024048E-01,  4.405716657638550E-01,  7.357195615768433E-01,  2.873047888278961E-01,&
         7.006355524063110E-01, -2.267676740884781E-01,  1.684177815914154E-01,  2.451081871986389E-01,&
        -6.897705197334290E-01, -1.359052062034607E-01, -1.217865824699402E+00,  6.268809437751770E-01,&
        -1.108817100524902E+00, -1.098538115620613E-01,  6.363938003778458E-02, -2.163156747817993E+00,&
         2.993230819702148E-01, -6.225543469190598E-02,  6.338689923286438E-01,  2.340336740016937E-01,&
         3.334980309009552E-01,  5.768545866012573E-01, -8.454492688179016E-01, -7.557854652404785E-01,&
        -6.227542161941528E-01, -1.105716824531555E+00,  2.116404175758362E-01, -2.117430865764618E-01,&
        -1.036560058593750E+00, -1.257222741842270E-01,  5.264365077018738E-01, -1.787502527236938E+00,&
        -6.102513074874878E-01, -1.036811590194702E+00, -1.041777491569519E+00,  6.762499362230301E-02,&
        -1.829331994056702E+00, -1.342972517013550E-01,  2.181535959243774E+00,  7.125011086463928E-01,&
         9.849542975425720E-01,  4.515964090824127E-01, -5.667360424995422E-01,  1.371907234191895E+00,&
         4.193291962146759E-01, -4.483173191547394E-01,  1.056447148323059E+00, -4.035096466541290E-01,&
         2.473213225603104E-01,  4.283659458160400E-01, -1.105738878250122E+00, -3.882422149181366E-01,&
         1.359030008316040E-01, -1.316889882087708E+00,  1.206199750304222E-01, -2.816296517848969E-01,&
        -3.856543898582458E-01, -1.341159194707870E-01,  2.931591272354126E-01, -8.115946650505066E-01,&
         1.549627929925919E-01, -3.494594991207123E-02,  1.392071247100830E-01,  8.500702381134033E-01,&
        -1.105314135551453E+00, -8.855208158493042E-01, -1.129539161920547E-01, -7.288187742233276E-01,&
         2.031663209199905E-01, -2.040854692459106E-01, -2.651244997978210E-01,  6.747405529022217E-01,&
         6.289814710617065E-01,  3.702930510044098E-01,  8.955963253974915E-01, -1.791490912437439E-01,&
         6.291658878326416E-01,  3.181912600994110E-01, -7.458741664886475E-01, -5.797970294952393E-01,&
         8.048549294471741E-03, -1.517996788024902E+00,  1.586797833442688E-02, -1.968807131052017E-01,&
        -6.696819067001343E-01,  2.561997175216675E-01,  1.585537791252136E-01, -3.939553797245026E-01,&
         1.001605153083801E+00, -3.178015723824501E-02,  2.169712930917740E-01,  7.597719430923462E-01,&
        -8.711787462234497E-01, -2.590858340263367E-01, -4.994206726551056E-01, -1.350332260131836E+00,&
        -1.754350513219833E-01, -5.298053622245789E-01, -1.044484019279480E+00, -5.103482306003571E-02,&
         8.845404386520386E-01,  4.584137201309204E-01,  1.076861619949341E+00,  1.874905377626419E-01,&
         2.787777185440063E-01,  8.369036912918091E-01, -8.217707276344299E-01, -2.826712131500244E-01,&
        -2.450734227895737E-01, -8.279343843460083E-01,  3.510917425155640E-01, -3.488889932632446E-01,&
        -7.627615332603455E-01,  3.606846034526825E-01,  5.258455872535706E-01, -5.099301040172577E-02,&
         6.352093815803528E-01, -1.835833787918091E-01,  1.247637987136841E+00,  5.917957425117493E-01,&
         1.019452288746834E-01, -5.673841834068298E-01,  1.377126276493073E-01, -1.055184245109558E+00,&
        -2.036373913288116E-01, -6.316062808036804E-01, -3.354403078556061E-01,  3.826665878295898E-01,&
        -6.721435189247131E-01, -6.410418748855591E-01, -1.417969822883606E+00, -8.955898880958557E-02,&
        -6.617363095283508E-01, -6.313887238502502E-01,  1.284139454364777E-01, -7.438000291585922E-02,&
         3.091568231582642E+00,  8.395515084266663E-01,  7.227233052253723E-01,  8.192335367202759E-01,&
        -2.106423974037170E-01,  2.122008800506592E+00,  7.060149908065796E-01,  3.394779860973358E-01,&
         6.117095947265625E-01, -3.271679580211639E-01,  1.616740077733994E-01,  1.569840312004089E-01,&
        -1.123665213584900E+00,  3.844760954380035E-01,  2.845884263515472E-01,  7.137780785560608E-01,&
         1.460106819868088E-01, -1.021391227841377E-01,  5.172263383865356E-01, -7.423986196517944E-01,&
        -2.789774909615517E-02, -1.258952766656876E-01, -1.325458526611328E+00, -5.270438194274902E-01,&
        -3.967397287487984E-02, -2.709308564662933E-01,  1.340401768684387E-01, -6.963784694671631E-01,&
        -3.221498429775238E-01, -8.531031608581543E-01,  3.377375304698944E-01,  1.652107536792755E-01,&
        -3.512997031211853E-01, -1.630981415510178E-01,  3.690161705017090E-01,  1.549807284027338E-02,&
         1.193455934524536E+00,  2.675475478172302E-01,  3.856497108936310E-01,  9.223973155021667E-01,&
        -8.005780726671219E-02,  7.949089407920837E-01,  1.678814589977264E-01,  5.589793920516968E-01,&
        -2.890521883964539E-01, -6.459630280733109E-02,  1.577395349740982E-01, -6.019581556320190E-01,&
         1.361452788114548E-01, -1.461234450340271E+00,  2.132855653762817E-01, -7.116237878799438E-01,&
        -1.837224513292313E-01,  6.981704831123352E-01, -1.456485867500305E+00, -8.896524459123611E-02,&
        -6.985316872596741E-01, -9.188821911811829E-01, -1.798982769250870E-01, -3.445543348789215E-01,&
        -9.767906665802002E-01,  6.575983762741089E-01, -5.698328614234924E-01,  2.794421613216400E-01,&
        -9.889149665832520E-01,  2.113757282495499E-01, -4.894487261772156E-01, -9.110729694366455E-01,&
         3.156659901142120E-01, -8.372070193290710E-01,  1.710339263081551E-02, -7.162731885910034E-01,&
        -9.848624467849731E-02, -2.407071143388748E-01, -4.630023241043091E-01,  5.028110146522522E-01 &
        /), shape(hidden1Synapse))
 
      outputSynapse = reshape((/ &
        -1.209702730178833E+00,  1.183213353157043E+00, -1.019356846809387E+00, -1.344744205474854E+00,&
        -1.445307731628418E+00,  1.024327754974365E+00, -1.584630727767944E+00,  1.083521246910095E+00,&
        -1.308865427970886E+00, -1.247952342033386E+00,  1.239847064018250E+00,  1.287056356668472E-01,&
         9.846584796905518E-01, -1.553632378578186E+00, -1.231866717338562E+00,  4.489912092685699E-02,&
         1.253254055976868E+00, -1.430614471435547E+00,  1.041161060333252E+00, -1.605084300041199E+00,&
         1.527578949928284E+00, -1.474965572357178E+00,  1.355290770530701E+00, -1.745877861976624E+00,&
         1.712602972984314E+00, -1.563431382179260E+00,  8.333104252815247E-01, -1.541154265403748E+00,&
        -1.556280970573425E+00,  7.898001670837402E-01,  1.451943874359131E+00,  1.376102089881897E+00,&
        -1.475358963012695E+00, -1.508958697319031E+00,  1.723131775856018E+00,  1.577485084533691E+00,&
         2.009120136499405E-01, -1.543342947959900E+00, -1.532042622566223E+00, -1.665173649787903E+00,&
        -1.577844977378845E+00,  1.509271860122681E+00, -1.648273229598999E+00, -1.399203181266785E+00,&
        -1.230364322662354E+00,  1.090018987655640E+00, -7.097014784812927E-01,  1.677408456802368E+00,&
        -1.743194699287415E+00, -1.423129081726074E+00,  7.856354713439941E-01,  1.262704372406006E+00,&
         1.029602646827698E+00, -8.157435655593872E-01, -1.168590903282166E+00, -1.007120013237000E+00,&
         1.498046159744263E+00, -1.094031929969788E+00,  1.288908720016479E+00, -1.570232629776001E+00,&
         1.331548571586609E+00, -1.591911792755127E+00,  1.173869848251343E+00, -1.569446206092834E+00,&
         1.071457147598267E+00, -1.386015534400940E+00,  1.319629669189453E+00, -1.251965403556824E+00,&
        -1.506981730461121E+00, -5.631150603294373E-01,  1.476744890213013E+00,  1.224819302558899E+00,&
        -1.190375804901123E+00, -4.876171946525574E-01,  1.674062848091125E+00,  1.343202710151672E+00,&
         8.375900387763977E-01, -1.624152183532715E+00, -1.477828741073608E+00, -1.320914030075073E+00,&
        -1.082759499549866E+00,  1.309733152389526E+00, -5.913071632385254E-01, -1.292264103889465E+00,&
        -1.440814852714539E+00,  1.020094513893127E+00, -1.208431601524353E+00,  1.691915869712830E+00,&
        -1.277797341346741E+00, -1.482174158096313E+00,  1.266713261604309E+00,  1.296367645263672E+00,&
         1.238657712936401E+00, -7.025628685951233E-01,  2.491326481103897E-01, -1.536825418472290E+00,&
         1.577931523323059E+00, -1.065637469291687E+00,  1.696800708770752E+00, -1.695444345474243E+00,&
         1.581656932830811E+00, -1.088520646095276E+00,  1.492973804473877E+00, -1.063908934593201E+00,&
         1.496415257453918E+00, -1.486176609992981E+00,  6.039925217628479E-01, -1.485497832298279E+00,&
        -1.147870540618896E+00, -1.266431331634521E+00,  1.607187867164612E+00,  1.494379520416260E+00,&
        -1.001191616058350E+00, -1.084854602813721E+00,  1.410489916801453E+00,  1.581320643424988E+00,&
         1.205576062202454E+00, -1.245357394218445E+00, -1.343545675277710E+00, -1.709581851959229E+00 &
        /), shape(outputSynapse))
 
      END SUBROUTINE Breadboard4         
!
!------------------------------------------------------------------------------
!
      SUBROUTINE Breadboard5(inputFile,hidden1Axon,hidden1Synapse,outputSynapse)
 
      implicit none
 
      real inputFile(2,7)
      real hidden1Axon(40)
      real hidden1Synapse(7,40)
      real outputSynapse(40,3)
 
      inputFile = reshape((/ &
         1.077844262123108E+00, -1.778443008661270E-01,&
         2.188449800014496E-01,  1.674167998135090E-02,&
         1.918158382177353E-01,  6.903452277183533E-01,&
         3.361344337463379E-01,  4.151264205574989E-02,&
         2.485501170158386E-01,  2.868268489837646E-01,&
         1.839550286531448E-01,  3.534696102142334E-01,&
         1.951007992029190E-01, -4.725341200828552E-01 &
        /), shape(inputFile))
 
      hidden1Axon = &
       (/3.177257776260376E-01, -3.444353640079498E-01,  5.270494818687439E-01, -5.221590399742126E-01,&
        -2.202716171741486E-01, -4.241476655006409E-01,  2.620704658329487E-02,  6.034846901893616E-01,&
        -3.619376122951508E-01, -3.380794525146484E-01,  4.901479184627533E-02,  4.951947927474976E-02,&
         1.800213754177094E-01, -2.407073378562927E-01, -3.286456167697906E-01, -6.795548200607300E-01,&
        -5.868792533874512E-01, -3.454326987266541E-01,  1.429300457239151E-01, -2.292728424072266E-01,&
         4.302643239498138E-01, -2.324737906455994E-01, -4.539224207401276E-01,  5.544423460960388E-01,&
        -4.054053127765656E-01, -1.476568281650543E-01, -2.141656428575516E-01,  1.077265888452530E-01,&
         5.846756696701050E-01,  3.272875547409058E-01,  1.847147941589355E-03, -4.990870654582977E-01,&
         1.531988829374313E-01,  1.791626960039139E-01, -6.736395359039307E-01, -5.093495845794678E-01,&
        -6.099227815866470E-02,  3.861090838909149E-01, -6.592265367507935E-01, -2.490588128566742E-01/)
 
      hidden1Synapse = reshape((/ &
         3.541271016001701E-02, -7.549672126770020E-01, -4.738137125968933E-01, -2.348672598600388E-03,&
        -2.733762562274933E-01, -8.357829414308071E-03, -8.771334886550903E-01, -2.402636408805847E-01,&
        -3.840126693248749E-01, -5.802615284919739E-01,  1.073393039405346E-03, -2.714654207229614E-01,&
        -1.682563573122025E-01,  2.412795424461365E-01,  6.722061038017273E-01, -2.907541096210480E-01,&
         1.961677670478821E-01, -3.303197622299194E-01,  1.424128562211990E-01,  5.971218943595886E-01,&
        -3.415485620498657E-01, -3.709296286106110E-01,  2.636498510837555E-01, -6.461778879165649E-01,&
        -4.282482266426086E-01, -1.192058548331261E-01, -7.758595943450928E-01, -4.671352729201317E-02,&
        -2.137460708618164E-01, -1.528403162956238E-02, -7.986806631088257E-01, -3.911508247256279E-02,&
        -5.328277871012688E-02, -6.519866585731506E-01,  3.402085006237030E-01,  1.100756451487541E-01,&
         6.820629835128784E-01,  7.288114726543427E-02,  2.484970390796661E-01, -1.383271068334579E-01,&
         1.246754452586174E-01,  6.508666276931763E-01,  3.158373534679413E-01, -5.986170172691345E-01,&
         6.103343367576599E-01, -6.012113094329834E-01, -1.359632611274719E-01, -2.586761862039566E-02,&
        -4.111338853836060E-01,  1.772232651710510E-01, -6.230232119560242E-01,  3.960133790969849E-01,&
        -6.472764015197754E-01, -3.764366805553436E-01, -9.892498701810837E-02, -9.984154999256134E-02,&
        -4.294761717319489E-01, -2.304461598396301E-01, -7.071238160133362E-01, -4.068204462528229E-01,&
        -4.626799225807190E-01, -3.020684123039246E-01,  6.521416902542114E-01,  1.521919965744019E-01,&
        -7.091572284698486E-01, -4.207086861133575E-01, -5.045717954635620E-01, -3.018378615379333E-01,&
        -4.485827982425690E-01, -5.111956596374512E-01, -8.567054569721222E-02,  4.856635630130768E-01,&
         2.459491789340973E-01, -1.496585756540298E-01, -1.183001995086670E-01,  4.713786244392395E-01,&
        -2.809847891330719E-01,  8.547450602054596E-02, -3.530589640140533E-01, -7.254429459571838E-01,&
        -1.860966980457306E-01, -6.639543771743774E-01,  4.769657552242279E-01, -7.412918210029602E-01,&
         3.024796843528748E-01, -6.272576451301575E-01, -5.452296733856201E-01, -2.242822349071503E-01,&
        -3.738160133361816E-01,  3.284691274166107E-01, -4.564896821975708E-01,  2.556349933147430E-01,&
         4.318492487072945E-02, -1.320876032114029E-01, -9.898099303245544E-02,  6.774403899908066E-02,&
         1.919083893299103E-01,  2.400640696287155E-01,  4.077304899692535E-01,  2.524036169052124E-01,&
         5.042297840118408E-01,  2.886471152305603E-01, -1.700776815414429E-01, -2.435589283704758E-01,&
        -2.057165205478668E-01,  1.996059715747833E-01,  2.711705565452576E-01,  3.861612975597382E-01,&
        -2.083975523710251E-01,  7.296724617481232E-02, -2.396509945392609E-01, -1.525006294250488E-01,&
        -4.502384066581726E-01, -5.351938009262085E-01, -3.890139460563660E-01,  1.700514107942581E-01,&
        -4.677065312862396E-01, -3.514041006565094E-01,  4.196007549762726E-01,  2.812465429306030E-01,&
        -2.938374876976013E-01, -3.160441517829895E-01, -4.980419874191284E-01,  3.127529323101044E-01,&
         2.271771281957626E-01, -1.466843336820602E-01, -6.397774219512939E-01,  4.446669816970825E-01,&
         8.942086249589920E-02,  9.681937843561172E-02, -5.533168092370033E-02, -4.528337121009827E-01,&
         6.882410049438477E-01, -3.133308887481689E-01, -2.058080136775970E-01, -2.226170003414154E-01,&
        -2.296325266361237E-01, -2.966837584972382E-01, -3.301460444927216E-01, -3.557955026626587E-01,&
         3.304032683372498E-01, -8.399857580661774E-02,  4.199078381061554E-01,  1.194518618285656E-02,&
         7.232509851455688E-01,  9.784302115440369E-02, -1.134829670190811E-01,  1.034526005387306E-01,&
        -8.523296117782593E-01,  5.190717577934265E-01,  5.323929339647293E-02,  1.697375029325485E-01,&
         5.581731796264648E-01, -9.171869754791260E-01, -1.815564483404160E-01,  3.742720186710358E-01,&
        -2.523972094058990E-01,  1.490504741668701E-01, -6.334505081176758E-01,  2.519290745258331E-01,&
         2.056387513875961E-01, -1.307390183210373E-01, -9.355121254920959E-01, -2.585434913635254E-01,&
        -4.636541008949280E-02, -1.257960349321365E-01,  1.712975054979324E-01, -7.756385207176208E-01,&
        -2.476336807012558E-01,  2.972539961338043E-01,  4.443784654140472E-01,  4.029458761215210E-02,&
        -2.695891633629799E-02, -1.858536303043365E-01, -1.682455986738205E-01, -1.443968862295151E-01,&
         3.042537868022919E-01, -4.171138703823090E-01, -1.896526068449020E-01,  1.934753060340881E-01,&
        -5.211362838745117E-01, -4.224704951047897E-02, -5.408123731613159E-01, -2.546814382076263E-01,&
        -3.727044463157654E-01, -4.361395835876465E-01,  1.507636755704880E-01,  8.203987777233124E-02,&
         1.366124451160431E-01,  5.710709095001221E-01,  3.028809726238251E-01,  9.636782407760620E-01,&
        -3.770071640610695E-02,  3.973050415515900E-01,  2.884645946323872E-03, -8.364310860633850E-01,&
         5.341901779174805E-01, -1.418879022821784E-03,  5.416565537452698E-01,  3.877540528774261E-01,&
        -1.585132908076048E-03,  1.770619601011276E-01,  4.701207578182220E-02,  4.187163114547729E-01,&
         9.934148788452148E-01,  2.260543704032898E-01,  7.113759517669678E-01,  4.728879332542419E-01,&
        -3.471966087818146E-01,  7.732371240854263E-02, -2.182047963142395E-01,  8.698941469192505E-01,&
         6.959328651428223E-01,  1.184082403779030E-01,  1.408587545156479E-01,  2.005882859230042E-01,&
         3.091167509555817E-01, -1.955157965421677E-01, -2.792426571249962E-02, -7.336559891700745E-02,&
         1.834385395050049E-01, -3.164150416851044E-01, -5.837532281875610E-01,  9.843266010284424E-01,&
        -5.053303837776184E-01,  9.432902336120605E-01,  2.762463316321373E-02,  3.678649663925171E-01,&
        -8.084134012460709E-02,  2.041484862565994E-01,  5.061163306236267E-01,  7.991071939468384E-01,&
         2.264233529567719E-01,  7.115226387977600E-01, -5.186138153076172E-01,  4.093891084194183E-01,&
        -1.001899018883705E-01, -1.933344826102257E-02,  1.815729439258575E-01, -1.810713559389114E-01,&
        -5.504883527755737E-01,  7.005249857902527E-01, -1.967341639101505E-02,  1.448700390756130E-02,&
         3.791421651840210E-01, -3.687309324741364E-01,  6.238684058189392E-01,  2.549594640731812E-02,&
         6.611171960830688E-01, -2.348230034112930E-01,  4.087108075618744E-01,  1.835047304630280E-01,&
         2.745413780212402E-01, -5.477424860000610E-01,  4.227129369974136E-02,  1.370747834444046E-01,&
        -1.771535575389862E-01,  2.915630638599396E-01,  8.117929100990295E-02, -5.147354602813721E-01,&
        -7.195407748222351E-01, -2.950702905654907E-01, -8.272841572761536E-01, -8.926602080464363E-03,&
         6.488984823226929E-01, -7.542604207992554E-01, -1.718278229236603E-01, -4.908424615859985E-02,&
        -3.619753718376160E-01, -9.747832268476486E-02, -9.625122696161270E-02, -1.545960754156113E-01,&
         4.842050671577454E-01, -9.618758410215378E-02,  1.017526090145111E-01, -1.527849882841110E-01,&
         5.150741338729858E-01, -2.614658325910568E-02, -4.681808650493622E-01,  6.698484718799591E-02 &
        /), shape(hidden1Synapse))
 
      outputSynapse = reshape((/ &
        -4.252142608165741E-01, -5.190939903259277E-01,  2.900628745555878E-01, -4.749988615512848E-01,&
        -2.432068884372711E-01,  2.475018054246902E-01,  1.508098654448986E-02, -1.032671928405762E-01,&
        -5.695398449897766E-01, -4.341589808464050E-01,  3.563072979450226E-01, -1.610363721847534E-01,&
        -1.529531776905060E-01,  3.572074323892593E-02, -1.639768481254578E-01, -2.103261351585388E-01,&
        -5.111085772514343E-01, -9.769214689731598E-02, -1.570120900869370E-01, -1.928524225950241E-01,&
         4.143640100955963E-01, -3.950143232941628E-02, -2.028328180313110E-01, -1.475265175104141E-01,&
        -2.296919003129005E-02, -3.979336936026812E-03, -3.908852040767670E-01,  4.192969501018524E-01,&
         2.397747188806534E-01,  4.962041378021240E-01,  4.480696618556976E-01, -2.336141020059586E-01,&
         3.938802778720856E-01,  2.352581322193146E-01,  1.772783696651459E-02, -5.289353057742119E-02,&
        -3.967223316431046E-02, -4.341553747653961E-01, -2.162312269210815E-01,  4.311326891183853E-02,&
         4.480128586292267E-01,  1.783114373683929E-01,  5.068565607070923E-01, -4.451447725296021E-01,&
        -5.096289515495300E-01, -4.807172119617462E-01,  1.144711822271347E-01,  3.887178003787994E-01,&
        -3.575057387351990E-01, -1.148879528045654E-01, -3.399987518787384E-02, -2.313354164361954E-01,&
        -7.217752188444138E-02,  3.657472431659698E-01,  3.738324940204620E-01,  4.177713990211487E-01,&
        -4.159389436244965E-01, -1.484509706497192E-01,  2.662932872772217E-01, -4.467738270759583E-01,&
         7.071519643068314E-02,  3.344006240367889E-01, -5.436876043677330E-02,  3.525221049785614E-01,&
        -2.395160868763924E-02, -3.141686320304871E-01,  3.852373957633972E-01,  4.932067096233368E-01,&
        -1.492380946874619E-01,  4.595996737480164E-01,  3.445216640830040E-02, -5.653984546661377E-01,&
        -4.437799155712128E-01,  1.460446715354919E-01, -4.742037057876587E-01,  1.456019878387451E-01,&
         3.867210447788239E-01,  4.871259629726410E-01, -4.954726397991180E-01,  1.770049333572388E-02,&
         2.028178423643112E-01, -3.220860958099365E-01,  2.971330881118774E-01, -1.783177554607391E-01,&
        -2.126741260290146E-01, -2.823735475540161E-01,  4.713099896907806E-01,  2.155631184577942E-01,&
        -3.713304102420807E-01,  2.199546098709106E-01,  2.943331003189087E-01,  4.534626007080078E-01,&
         3.414066731929779E-01, -1.535274535417557E-01, -1.036400645971298E-01, -4.483501911163330E-01,&
         8.723334968090057E-02, -1.368855964392424E-02, -5.010653138160706E-01,  4.472654759883881E-01,&
         1.106471717357635E-01,  5.139253139495850E-01, -2.296521663665771E-01,  4.545788764953613E-01,&
         1.664130948483944E-02,  2.438283525407314E-02, -1.943250745534897E-01,  4.952348470687866E-01,&
         3.839295804500580E-01, -3.456721901893616E-01, -1.650201976299286E-01, -3.892767727375031E-01,&
        -3.154349029064178E-01,  3.591218292713165E-01, -2.804268598556519E-01, -4.606449007987976E-01,&
         1.020256653428078E-01,  2.229744791984558E-01, -4.180959761142731E-01, -4.198006689548492E-01 &
        /), shape(outputSynapse))
 
      END SUBROUTINE Breadboard5         
!
!------------------------------------------------------------------------------
!
      SUBROUTINE Breadboard6(inputFile,hidden1Axon,hidden2Axon,&
                 hidden1Synapse,hidden2Synapse,outputSynapse)
 
      implicit none
 
      real inputFile(2,7)
      real hidden1Axon(7)
      real hidden2Axon(4)
      real hidden1Synapse(7,7)
      real hidden2Synapse(7,4)
      real outputSynapse(4,3)
 
      inputFile = reshape((/ &
         1.353383421897888E+00, -4.533834457397461E-01,&
         2.269289046525955E-01, -1.588500849902630E-02,&
         1.868382692337036E-01,  6.490761637687683E-01,&
         4.038590788841248E-01,  3.776083141565323E-02,&
         2.430133521556854E-01,  3.004860281944275E-01,&
         1.935067623853683E-01,  4.185551702976227E-01,&
         1.951007992029190E-01, -4.725341200828552E-01 &
        /), shape(inputFile))
 
     hidden1Axon = &
      (/ 7.384125608950853E-03, -2.202851057052612E+00,  2.003432661294937E-01, -2.467587143182755E-01,&
         5.973502993583679E-01,  3.834692537784576E-01,  2.687855064868927E-01/)
 
     hidden2Axon = &
      (/ 3.643750846385956E-01,  2.449363768100739E-01,  4.754272103309631E-01,  7.550075054168701E-01/)
 
     hidden1Synapse = reshape((/ &
         7.333400845527649E-01,  5.450296998023987E-01, -7.700046896934509E-01,  1.426693439483643E+00,&
        -1.024212338961661E-03, -6.459779292345047E-02,  1.028800487518311E+00, -2.116347402334213E-01,&
         3.591781139373779E+00,  2.435753583908081E+00, -6.687584519386292E-01,  1.201278567314148E+00,&
        -3.478864133358002E-01,  1.830960988998413E+00, -3.111673295497894E-01, -4.177703261375427E-01,&
        -3.920616805553436E-01, -5.040770769119263E-01, -5.354442000389099E-01, -1.534618530422449E-02,&
        -1.089364647865295E+00, -3.010036647319794E-01,  1.486289381980896E+00,  1.059559464454651E+00,&
         1.640596628189087E+00,  2.254628390073776E-01,  4.839954376220703E-01,  8.484285473823547E-01,&
        -6.926012784242630E-02,  4.926209524273872E-02,  2.834132313728333E-01,  3.028324842453003E-01,&
         2.161216735839844E-01,  7.251360416412354E-01,  2.851752638816833E-01, -5.653074979782104E-01,&
         3.640621304512024E-01,  1.341893225908279E-01,  7.511208057403564E-01, -1.088509336113930E-01,&
         1.044083759188652E-01,  6.529347300529480E-01, -6.885128021240234E-01, -1.003871187567711E-01,&
         9.337020665407181E-02, -4.425194561481476E-01, -3.668845295906067E-01, -2.661575675010681E-01,&
        -5.936880707740784E-01 &
        /), shape(hidden1Synapse))
 
     hidden2Synapse = reshape((/ &
        -5.461466908454895E-01, -1.490996479988098E+00,  7.721499800682068E-01, -3.842977285385132E-01,&
         1.134691461920738E-01, -7.171064615249634E-01,  4.990165829658508E-01, -4.233781099319458E-01,&
         5.502462983131409E-01, -1.000102013349533E-01,  1.481512188911438E+00,  1.637827455997467E-01,&
         5.879161506891251E-02, -3.256742060184479E-01,  4.237195849418640E-01,  1.471476674079895E+00,&
        -1.982609331607819E-01,  6.787789463996887E-01,  5.525223612785339E-01,  4.395257532596588E-01,&
         1.643348783254623E-01,  8.910947442054749E-01,  1.772162079811096E+00, -2.550726830959320E-01,&
         4.305597543716431E-01,  1.965346336364746E-01, -2.251276820898056E-01, -5.650298595428467E-01 &
        /), shape(hidden2Synapse))
 
     outputSynapse = reshape((/ &
         4.605286195874214E-02,  1.636024713516235E-01,  7.045555710792542E-01,  4.994805455207825E-01,&
         5.167593955993652E-01,  2.924540340900421E-01, -1.490857079625130E-02, -1.826021969318390E-01,&
         3.571106493473053E-01, -3.790216147899628E-01, -6.031348705291748E-01, -4.664786159992218E-01 &
        /), shape(outputSynapse))
 
      END SUBROUTINE Breadboard6         
!
!------------------------------------------------------------------------------
!
      SUBROUTINE Breadboard7(inputFile,hidden1Axon,hidden2Axon,&
                 hidden1Synapse,hidden2Synapse,outputSynapse)
 
      implicit none
 
      real inputFile(2,7)
      real hidden1Axon(7)
      real hidden2Axon(4)
      real hidden1Synapse(7,7)
      real hidden2Synapse(7,4)
      real outputSynapse(4,3)
 
      inputFile = reshape((/ &
         1.077844262123108E+00, -1.778443008661270E-01,&
         2.295625507831573E-01,  6.163756549358368E-02,&
         2.081165313720703E-01,  6.204994320869446E-01,&
         3.565062582492828E-01, -1.051693689078093E-02,&
         2.430133521556854E-01,  3.004860281944275E-01,&
         1.839550286531448E-01,  3.534696102142334E-01,&
         1.951007992029190E-01, -4.725341200828552E-01 &
        /), shape(inputFile))
 
     hidden1Axon = &
      (/-4.191969335079193E-01,  1.229978561401367E+00, -2.403785735368729E-01,  5.233071446418762E-01,&
         8.062141537666321E-01,  1.000604867935181E+00, -1.015548110008240E-01/)
 
     hidden2Axon = &
      (/-5.321261882781982E-01, -2.396449327468872E+00, -1.170158505439758E+00, -4.097367227077484E-01/)
 
     hidden1Synapse = reshape((/ &
         1.341468811035156E+00, -4.215665817260742E+00, -1.636691570281982E+00, -2.792109727859497E+00,&
        -1.489341259002686E+00,  4.075187742710114E-01, -2.091729402542114E+00, -5.029736161231995E-01,&
        -4.151493072509766E+00, -1.452428579330444E+00,  2.398953676223755E+00, -8.748555183410645E-01,&
         1.340690374374390E+00, -2.277854681015015E+00,  6.057588458061218E-01,  1.353034019470215E+00,&
        -1.214678883552551E+00, -3.864320814609528E-01,  1.148570895195007E+00,  5.792776346206665E-01,&
         1.344245020300150E-02, -8.885311484336853E-01, -1.594583272933960E+00,  4.960928857326508E-01,&
        -1.118881464004517E+00, -2.252289772033691E+00,  6.328870654106140E-01, -1.946701169013977E+00,&
        -2.910976111888885E-01,  2.447998225688934E-01,  2.001658976078033E-01, -1.229660585522652E-02,&
         6.969845890998840E-01, -5.897524300962687E-03, -5.688555836677551E-01,  2.619750201702118E-01,&
        -4.162483692169189E+00, -1.468571424484253E+00, -3.118389844894409E+00,  6.947994828224182E-01,&
        -2.687734663486481E-01, -2.110401153564453E+00,  3.224660456180573E-02,  8.378994464874268E-01,&
         9.896742701530457E-01, -7.354493737220764E-01,  6.684727072715759E-01,  1.465887904167175E+00,&
        -3.726872503757477E-01 &
        /), shape(hidden1Synapse))
 
     hidden2Synapse = reshape((/ &
        -3.395457863807678E-01, -5.815528631210327E-01,  2.929831743240356E-01, -5.629656314849854E-01,&
         4.701104387640953E-02, -9.300172328948975E-01, -1.461120098829269E-01, -3.458845615386963E-01,&
         1.266251802444458E-01,  6.342335790395737E-02,  1.869771480560303E-01, -1.476681977510452E-01,&
         5.144428834319115E-02, -3.145390946883708E-04,  8.697064518928528E-01,  1.057970225811005E-01,&
         2.603019773960114E-01,  4.393529295921326E-01, -2.832717299461365E-01,  5.771816968917847E-01,&
        -3.896601796150208E-01, -7.260112762451172E-01, -7.957320213317871E-01,  6.776907294988632E-02,&
        -3.073690235614777E-01, -1.540119051933289E-01, -6.733091473579407E-01,  2.009786069393158E-01 &
        /), shape(hidden2Synapse))
 
     outputSynapse = reshape((/ &
         3.156347572803497E-01, -8.236174583435059E-01, -9.946570396423340E-01,  4.212915897369385E-01,&
        -7.918102145195007E-01, -2.033229321241379E-01, -1.056663155555725E+00, -5.699685215950012E-01,&
        -9.666987657546997E-01, -5.505290031433105E-01,  8.724089711904526E-02, -9.536570906639099E-01 &
        /), shape(outputSynapse))
 
      END SUBROUTINE Breadboard7         
!
!------------------------------------------------------------------------------
!
      SUBROUTINE Breadboard8(inputFile,hidden1Axon,hidden2Axon,&
                 hidden1Synapse,hidden2Synapse,outputSynapse)
 
      implicit none
 
      real inputFile(2,7)
      real hidden1Axon(7)
      real hidden2Axon(4)
      real hidden1Synapse(7,7)
      real hidden2Synapse(7,4)
      real outputSynapse(4,3)
 
      inputFile = reshape((/ &
         1.353383421897888E+00, -4.533834457397461E-01,&
         2.188449800014496E-01,  1.674167998135090E-02,&
         1.906577646732330E-01,  6.807435750961304E-01,&
         3.361344337463379E-01,  4.151264205574989E-02,&
         2.491349428892136E-01,  3.307266235351562E-01,&
         1.839550286531448E-01,  3.534696102142334E-01,&
         1.951007992029190E-01, -4.725341200828552E-01 &
        /), shape(inputFile))
 
     hidden1Axon = &
      (/-3.274627029895782E-01,  2.668272238224745E-03, -3.019839525222778E-01, -4.557206928730011E-01,&
        -5.515558272600174E-02,  3.119016764685512E-04,  8.753398060798645E-02/)
 
     hidden2Axon = &
      (/ 2.733168303966522E-01, -3.423235416412354E-01,  8.666662573814392E-01, -6.124708056449890E-01/)
 
     hidden1Synapse = reshape((/ &
         2.732226848602295E-01,  1.847893238067627E+00, -1.084923520684242E-01,  1.385403037071228E+00,&
         2.885355055332184E-01, -3.135629594326019E-01,  1.057805895805359E+00, -5.868541821837425E-02,&
         3.278825521469116E+00,  4.641786217689514E-01,  4.461606740951538E-01, -1.952850073575974E-01,&
        -5.789646506309509E-01,  1.945697903633118E+00, -9.578172862529755E-02,  2.150904417037964E+00,&
         9.114052653312683E-01,  1.107189536094666E+00,  6.752110123634338E-01,  2.475811988115311E-01,&
         1.050705909729004E+00,  3.205673992633820E-01,  2.478840798139572E-01, -5.084273815155029E-01,&
        -2.407394796609879E-01, -1.702371835708618E-01,  1.456947028636932E-01,  3.221787512302399E-01,&
        -2.719256579875946E-01, -5.116361379623413E-01,  3.973563387989998E-02, -1.733802706003189E-01,&
        -1.649789661169052E-01, -4.471102654933929E-01, -4.071239829063416E-01, -1.492276042699814E-01,&
        -1.245773434638977E+00, -6.851593255996704E-01, -8.733592033386230E-01, -4.348643422126770E-01,&
        -3.520536422729492E-01, -9.930510520935059E-01,  1.956800930202007E-02, -9.781590104103088E-01,&
        -6.039583683013916E-01, -6.923800706863403E-01, -6.682770848274231E-01,  4.162513464689255E-02,&
        -1.004322052001953E+00 &
        /), shape(hidden1Synapse))
 
     hidden2Synapse = reshape((/ &
        -8.183520436286926E-01, -1.621446132659912E+00, -1.045793533325195E+00, -5.855653062462807E-02,&
         4.404523968696594E-01,  7.002395391464233E-01,  2.097517400979996E-01, -9.925779700279236E-02,&
        -8.263560533523560E-01, -1.043026208877563E+00,  4.524357020854950E-01,  2.231711596250534E-01,&
         8.736496567726135E-01,  8.797182440757751E-01,  6.963157653808594E-01,  2.816314399242401E-01,&
         1.525615751743317E-01,  1.936565339565277E-01,  1.900831162929535E-01,  1.180221140384674E-01,&
         1.027775928378105E-01,  9.149055480957031E-01,  1.129598617553711E+00,  6.131598353385925E-01,&
         2.547058761119843E-01,  2.556352131068707E-02, -3.627143800258636E-02, -6.722733378410339E-01 &
        /), shape(hidden2Synapse))
 
     outputSynapse = reshape((/ &
        -5.266965627670288E-01, -1.973343640565872E-01,  1.362649053335190E-01,  9.479679167270660E-02,&
         2.987665235996246E-01, -3.116582632064819E-01, -1.842434853315353E-01, -4.986568093299866E-01,&
         6.261917948722839E-01,  5.454919338226318E-01, -3.484728187322617E-02, -4.687039256095886E-01 &
        /), shape(outputSynapse))
 
      END SUBROUTINE Breadboard8         
!
!------------------------------------------------------------------------------
!
      SUBROUTINE Breadboard9(inputFile,hidden1Axon,hidden2Axon,&
                 hidden1Synapse,hidden2Synapse,outputSynapse)
 
      implicit none
 
      real inputFile(2,7)
      real hidden1Axon(7)
      real hidden2Axon(4)
      real hidden1Synapse(7,7)
      real hidden2Synapse(7,4)
      real outputSynapse(4,3)
 
      inputFile = reshape((/ &
         1.077844262123108E+00, -1.778443008661270E-01,&
         2.188449800014496E-01,  1.674167998135090E-02,&
         1.868382692337036E-01,  6.490761637687683E-01,&
         3.733665347099304E-01,  1.051026657223701E-01,&
         2.430133521556854E-01,  3.004860281944275E-01,&
         2.083092182874680E-01,  3.581876754760742E-01,&
         1.951007992029190E-01, -4.725341200828552E-01 &
        /), shape(inputFile))
 
     hidden1Axon = &
      (/ 1.012814998626709E+00, -3.782782554626465E-01, -2.220184087753296E+00, -3.424299955368042E-01,&
         1.449530482292175E+00, -2.592789530754089E-01, -4.670010507106781E-01/)
 
     hidden2Axon = &
      (/ 3.516010642051697E-01,  3.293374776840210E-01, -1.675553172826767E-01,  3.799068629741669E-01/)
 
     hidden1Synapse = reshape((/ &
         1.390573829412460E-01, -3.110583126544952E-01,  1.105552077293396E+00,  4.394045472145081E-01,&
         4.795211851596832E-01,  1.969023197889328E-01,  5.574952811002731E-02,  1.690310984849930E-01,&
         2.208244323730469E+00,  2.111947536468506E+00,  3.239532709121704E-01,  7.690296173095703E-01,&
         1.264077782630920E+00,  1.672740578651428E+00,  1.320844173431396E+00,  7.965675592422485E-01,&
        -7.341063618659973E-01,  3.702043294906616E+00,  1.716022133827209E+00, -6.642882823944092E-01,&
         1.686427950859070E+00, -4.863217473030090E-01,  1.285641908645630E+00,  1.281449794769287E+00,&
         2.356275558471680E+00, -1.406845331192017E+00,  6.027717590332031E-01,  6.652191877365112E-01,&
        -9.871492385864258E-01, -5.513690948486328E+00, -2.750334143638611E-01,  1.229651212692261E+00,&
        -2.504641294479370E+00, -3.219850361347198E-01, -2.744197607040405E+00, -4.023179113864899E-01,&
         9.932321496307850E-03, -6.916724443435669E-01, -2.260914087295532E+00,  1.261568814516068E-01,&
         3.248662948608398E-01,  6.963043808937073E-01,  1.830800414085388E+00, -2.054267644882202E+00,&
        -9.595731496810913E-01, -8.711494207382202E-01, -1.330682396888733E+00,  2.109736204147339E+00,&
        -6.145163774490356E-01 &
        /), shape(hidden1Synapse))
 
     hidden2Synapse = reshape((/ &
        -3.299105465412140E-01,  4.235435724258423E-01,  9.191738963127136E-01,  6.795659661293030E-01,&
        -1.440919041633606E+00,  4.634908214211464E-02, -1.265781879425049E+00,  2.394487708806992E-01,&
         1.205053567886353E+00,  5.790516138076782E-01,  1.087130665779114E+00, -6.723164916038513E-01,&
        -1.834900081157684E-01, -4.767680168151855E-01,  8.402896672487259E-02,  1.035530328750610E+00,&
         1.644443035125732E+00,  4.317290484905243E-01, -1.714672803878784E+00,  5.225644707679749E-01,&
        -5.602287650108337E-01,  1.068559288978577E+00, -2.211284125223756E-03, -2.943626642227173E-01,&
         1.341261714696884E-01,  4.324447214603424E-01, -5.482236146926880E-01, -4.985276758670807E-01 &
        /), shape(hidden2Synapse))
 
     outputSynapse = reshape((/ &
         3.726457059383392E-01,  7.749153375625610E-01,  4.159255921840668E-01,  5.234625935554504E-01,&
        -1.592817008495331E-01,  5.884559154510498E-01, -7.756121158599854E-01,  2.137655019760132E-01,&
        -6.172903776168823E-01, -4.417923986911774E-01, -4.576872885227203E-01,  4.440903961658478E-01 &
        /), shape(outputSynapse))
 
      END SUBROUTINE Breadboard9         
!
!------------------------------------------------------------------------------
!
      SUBROUTINE Breadboard10(inputFile,hidden1Axon,hidden2Axon,&
                 hidden1Synapse,hidden2Synapse,outputSynapse)
 
      implicit none
 
      real inputFile(2,7)
      real hidden1Axon(7)
      real hidden2Axon(4)
      real hidden1Synapse(7,7)
      real hidden2Synapse(7,4)
      real outputSynapse(4,3)
 
      inputFile = reshape((/ &
         1.077844262123108E+00, -1.778443008661270E-01,&
         2.269289046525955E-01, -1.588500849902630E-02,&
         1.906577646732330E-01,  6.807435750961304E-01,&
         3.703703582286835E-01, -4.592590779066086E-02,&
         2.611723542213440E-01,  3.901915252208710E-01,&
         1.911842674016953E-01,  4.027296602725983E-01,&
         1.951007992029190E-01, -4.725341200828552E-01 &
        /), shape(inputFile))
 
     hidden1Axon = &
      (/ 1.307985544204712E+00, -1.960705667734146E-01, -1.105142459273338E-01, -1.207442641258240E+00,&
        -1.665081620216370E+00,  1.251117825508118E+00, -7.307677268981934E-01/)
 
     hidden2Axon = &
      (/ 2.186001092195511E-02,  3.369570672512054E-01,  1.165086925029755E-01,  2.747000660747290E-03/)
 
     hidden1Synapse = reshape((/ &
        -3.375437259674072E-01, -3.020816326141357E+00, -1.435481071472168E+00,  1.473870635032654E+00,&
        -7.776365280151367E-01,  6.734371185302734E-01, -1.643768787384033E+00, -1.227448821067810E+00,&
        -7.365036606788635E-01, -4.473563134670258E-01, -5.696173906326294E-01, -2.562220990657806E-01,&
         8.557485342025757E-01, -8.057124614715576E-01,  4.266147911548615E-01,  2.171551227569580E+00,&
         3.776189982891083E-01,  5.574828386306763E-01,  3.814708292484283E-01,  2.591066062450409E-01,&
         1.959651827812195E+00,  1.003962755203247E-01, -1.228965446352959E-02, -3.882043361663818E-01,&
        -2.722288109362125E-02, -3.378733694553375E-01, -7.981095314025879E-01,  4.839731752872467E-01,&
         1.432798147201538E+00,  1.885666996240616E-01, -6.051751971244812E-01,  2.924412488937378E+00,&
         1.136252880096436E+00,  2.994727194309235E-01,  1.604383468627930E+00, -8.440219759941101E-01,&
         6.088087558746338E-01, -3.722844421863556E-01,  5.441566109657288E-01,  3.944540619850159E-01,&
         7.044004201889038E-01,  3.459328413009644E-01,  1.054268121719360E+00, -3.348083496093750E+00,&
        -7.199336886405945E-01, -1.489133596420288E+00, -4.090557992458344E-01,  8.203456401824951E-01,&
        -1.118073821067810E+00 &
        /), shape(hidden1Synapse))
 
     hidden2Synapse = reshape((/ &
        -6.871775984764099E-01, -1.148896694183350E+00, -2.102893590927124E-01, -5.890849828720093E-01,&
         5.899340510368347E-01,  7.098034024238586E-01, -1.422515869140625E+00, -1.206974506378174E+00,&
         4.104525446891785E-01,  3.567897081375122E-01,  2.746991515159607E-01,  1.193219542503357E+00,&
         3.167707324028015E-01, -1.222744822502136E+00, -9.918631613254547E-02,  4.355156719684601E-01,&
         2.938420772552490E-01, -1.012830615043640E+00, -1.290418803691864E-01,  7.479285597801208E-01,&
        -2.292920649051666E-01, -1.372484922409058E+00, -6.534293759614229E-03,  1.525195717811584E+00,&
         2.076585590839386E-01,  1.434590101242065E+00,  7.887706905603409E-02, -1.401232123374939E+00 &
        /), shape(hidden2Synapse))
 
     outputSynapse = reshape((/ &
         6.101396083831787E-01,  3.122945129871368E-01,  3.869898915290833E-01,  4.438063502311707E-01,&
         5.161536335945129E-01, -2.700618803501129E-01, -3.105166740715504E-02, -5.569267272949219E-01,&
        -5.549081563949585E-01, -3.867979049682617E-01,  1.623111665248871E-01, -6.052750945091248E-01 &
        /), shape(outputSynapse))
 
      END SUBROUTINE Breadboard10        
!
!-------------------------------------------------------------------------------------
!
!> calslr_uutah() computes snow solid-liquid-ratio slr using the Steenburgh algorithm.
!>
!> Obtained the code and data from U of Utah Jim Steenburgh and Peter Veals.
!> SLR = m1X1 + m2X2 + m3X3 + m4X4 + m5X5 + m6X6 + b.
!>
!>      X1 = wind speed at at 1km above ground level (AGL) in m/s
!>      m1 = -0.174848
!>
!>      X2 = temperature at 2km AGL in Kelvin
!>      m2 = -0.52644
!>
!>      X3 = wind speed at 2 km AGL in m/s
!>      m3 = 0.034911
!>
!>      X4 = wind speed at 500 m AGL in m/s
!>      m4 = -0.270473
!>
!>      X5 = temperature at 1 km AGL in Kelvin
!>      m5 = 0.028299
!>
!>      X6 = temperature at 500 m AGL in m/s
!>      m6 = 0.096273
!>
!>      b =  118.35844
!>
!> @param[out] SLR real Solid snow to liquid ratio
!> 
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2023-01-23 | Jesse Meng | Initial
!>
!> @author Jesse Meng @date 2023-01-03

      SUBROUTINE CALSLR_UUTAH(SLR)

      use vrbls3d,    only: ZINT,ZMID,PMID,T,Q,UH,VH
      use masks,      only: LMH,HTM
      use ctlblk_mod, only: ISTA,IEND,JSTA,JEND,ista_2l,iend_2u,jsta_2l,jend_2u,&
                            LM,SPVAL

      implicit none

      real,dimension(ista_2l:iend_2u,jsta_2l:jend_2u),intent(out) :: slr !slr=snod/weasd=1000./sndens

      integer, parameter :: NFL=3
      real,    parameter :: HTFL(NFL)=(/ 500., 1000., 2000. /)
      real,dimension(ISTA:IEND,JSTA:JEND,NFL) :: TFD,UFD,VFD

      real LHL(NFL),DZABH(NFL),SWND(NFL)
      real HTSFC,HTABH,DZ,RDZ,DELT,DELU,DELV

      real, parameter :: m1 = -0.174848
      real, parameter :: m2 = -0.52644
      real, parameter :: m3 =  0.034911
      real, parameter :: m4 = -0.270473
      real, parameter :: m5 =  0.028299
      real, parameter :: m6 =  0.096273
      real, parameter ::  b =118.35844

      integer,dimension(ISTA:IEND,JSTA:JEND) :: KARR
      integer,dimension(ISTA:IEND,JSTA:JEND) :: TWET05
      real,dimension(ISTA:IEND,JSTA:JEND)    :: ZWET

      REAL, ALLOCATABLE :: TWET(:,:,:)

      integer I,J,L,LLMH,LMHK,IFD
!
!***************************************************************************
!
      ALLOCATE(TWET(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U,LM))

      DO IFD = 1,NFL
!$omp parallel do private(i,j)      
        DO J=JSTA,JEND
          DO I=ISTA,IEND
             TFD(I,J,IFD)     = SPVAL
             UFD(I,J,IFD)     = SPVAL
             VFD(I,J,IFD)     = SPVAL
          ENDDO
        ENDDO
      ENDDO        

!        LOCATE VERTICAL INDICES OF T,U,V, LEVEL JUST
!        ABOVE EACH FD LEVEL.

      DO J=JSTA,JEND
      DO I=ISTA,IEND
      IF(ZINT(I,J,LM+1)<SPVAL) THEN
         HTSFC = ZINT(I,J,LM+1)
         LLMH  = NINT(LMH(I,J))
      IFD = 1
      DO L = LLMH,1,-1
         HTABH = ZMID(I,J,L)-HTSFC
         IF(HTABH>HTFL(IFD)) THEN
            LHL(IFD) = L
            DZABH(IFD) = HTABH-HTFL(IFD)
            IFD = IFD + 1
         ENDIF
         IF(IFD > NFL) exit
      ENDDO

!        COMPUTE T, U, V AT FD LEVELS.

      DO IFD = 1,NFL 
         L = LHL(IFD)
         IF (L<LM .AND. T(I,J,L)<SPVAL .AND. UH(I,J,L)<SPVAL .AND. VH(I,J,L)<SPVAL) THEN
           DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
           RDZ  = 1./DZ
           DELT = T(I,J,L)-T(I,J,L+1)
           TFD(I,J,IFD) = T(I,J,L) - DELT*RDZ*DZABH(IFD)
           DELU = UH(I,J,L)-UH(I,J,L+1)
           DELV = VH(I,J,L)-VH(I,J,L+1)
           UFD(I,J,IFD) = UH(I,J,L) - DELU*RDZ*DZABH(IFD)
           VFD(I,J,IFD) = VH(I,J,L) - DELV*RDZ*DZABH(IFD)
         ELSE
           TFD(I,J,IFD) = T(I,J,L)
           UFD(I,J,IFD) = UH(I,J,L)
           VFD(I,J,IFD) = VH(I,J,L)
         ENDIF
      ENDDO
      ENDIF !IF(ZINT(I,J,LM+1)<SPVAL)
      ENDDO !I loop
      ENDDO !J loop

!        COMPUTE SLR

      SLR = SPVAL

!$omp parallel do private(i,j)      
      DO J=JSTA,JEND
      DO I=ISTA,IEND
      IF(TFD(I,J,1)<SPVAL .AND. UFD(I,J,1)<SPVAL .AND. VFD(I,J,1)<SPVAL) THEN
         SWND(1)=sqrt(UFD(I,J,1)*UFD(I,J,1)+VFD(I,J,1)*VFD(I,J,1))
         SWND(2)=sqrt(UFD(I,J,2)*UFD(I,J,2)+VFD(I,J,2)*VFD(I,J,2))
         SWND(3)=sqrt(UFD(I,J,3)*UFD(I,J,3)+VFD(I,J,3)*VFD(I,J,3))
         SLR(I,J) = m1*SWND(2)+m2*TFD(I,J,3)+m3*SWND(3)+m4*SWND(1)+m5*TFD(I,J,2)+m6*TFD(I,J,1)+b
         SLR(I,J) = max(SLR(I,J),3.)
      ENDIF
      ENDDO
      ENDDO

!        COMPUTE WETBULB TEMPERATURE AND SEARCH FOR TWET > 0.5C

      KARR = 1
      CALL WETBULB(T,Q,PMID,HTM,KARR,TWET)

!$omp parallel do private(i,j)      
      DO J=JSTA,JEND
      DO I=ISTA,IEND
         ZWET(I,J)=ZMID(I,J,LM)
         TWET05(I,J)=-1
      ENDDO
      ENDDO

      DO L=LM,1,-1
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
      DO I=ISTA,IEND
         IF(TWET05(I,J) < 0) THEN
            IF(TWET(I,J,L) <= 273.15+0.5) THEN
               ZWET(I,J)=ZMID(I,J,L)
               TWET05(I,J)=1
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO

!$omp parallel do private(i,j,HTABH)      
      DO J=JSTA,JEND
      DO I=ISTA,IEND
         IF(TWET05(I,J) > 0 .AND. SLR(I,J)<SPVAL) THEN
            HTABH=ZWET(I,J)-ZINT(I,J,LM+1)
            IF(HTABH<0.) HTABH=0.
            SLR(I,J)=SLR(I,J)*(1.-HTABH/200.)
            IF(SLR(I,J)<0.) SLR(I,J)=0.
         ENDIF
      ENDDO
      ENDDO

      DEALLOCATE (TWET)

      END SUBROUTINE CALSLR_UUTAH
!
!-------------------------------------------------------------------------------------
!
  end module upp_physics

