!> @file
!>
!> @brief upp_physics is a collection of UPP subroutines for physics variables calculation.
!> @author Jesse Meng @date 2020-05-20

!> calcape() computes CAPE/CINS and other storm related variables.
!>
!> calcape2() computes additional storm related variables.
!>
!> calrh(), calrh_nam(), calrh_gfs(), calrh_gsd() compute RH using various algorithms.
!>
!> The NAM v4.1.18 algorithm (calrh_nam()) is selected as default for 
!> NMMB and FV3GFS, FV3GEFS, and FV3R for the UPP 2020 unification.
!>
!> calrh_pw() algorithm use at GSD for RUC and Rapid Refresh.
!>
!> fpvsnew() computes saturation vapor pressure.
!>
!> tvirtual() computes virtual temperature.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2020-05-20 | Jesse Meng | Initial
!>
!> @author Jesse Meng @date 2020-05-20
  module upp_physics

  implicit none

  private

  public :: CALCAPE, CALCAPE2
  public :: CALRH
  public :: CALRH_GFS, CALRH_GSD, CALRH_NAM
  public :: CALRH_PW
  public :: FPVSNEW
  public :: TVIRTUAL

  contains
!
!-------------------------------------------------------------------------------------
!
      SUBROUTINE CALRH(P1,T1,Q1,RH)

      use ctlblk_mod, only: im, jsta, jend, MODELNAME
      implicit none

      REAL,dimension(IM,jsta:jend),intent(in)    :: P1,T1
      REAL,dimension(IM,jsta:jend),intent(inout) :: Q1
      REAL,dimension(IM,jsta:jend),intent(out)   :: RH

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
!> @param[in] Q1 Specific humidity (kg/kg)
!> @param[out] RH Relative humidity  (decimal form)
!> @param[out] Q1 Specific humidity (kg/kg)
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
     use ctlblk_mod, only: jsta, jend, spval, im
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     SET PARAMETER.
!
!     DECLARE VARIABLES.
!     
      REAL,dimension(IM,jsta:jend),intent(in)    :: P1,T1
      REAL,dimension(IM,jsta:jend),intent(inout) :: Q1
      REAL,dimension(IM,jsta:jend),intent(out)   :: RH
      REAL QC
      integer I,J
!***************************************************************
!
!     START CALRH.
!
      DO J=JSTA,JEND
        DO I=1,IM
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
!> @param[in] Q1 Specific humidity (kg/kg)
!> @param[out] RH Relative humidity  (decimal form)
!> @param[out] Q1 Specific humidity (kg/kg)
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
      use ctlblk_mod, only: jsta, jend, spval, im
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
      REAL,dimension(IM,jsta:jend),intent(in)   :: P1,T1
      REAL,dimension(IM,jsta:jend),intent(inout):: Q1,RH
      REAL ES,QC
      integer :: I,J
!***************************************************************
!
!     START CALRH.
!
!$omp parallel do private(i,j,es,qc)
      DO J=JSTA,JEND
        DO I=1,IM
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
!      
      SUBROUTINE CALRH_GSD(P1,T1,Q1,RHB)
!
! Algorithm use at GSD for RUC and Rapid Refresh                           
!------------------------------------------------------------------
!

      use ctlblk_mod, only: jsta, jend, im, spval

      implicit none

      integer :: j, i
      real :: tx, pol, esx, es, e
      real, dimension(im,jsta:jend) :: P1, T1, Q1, RHB


      DO J=JSTA,JEND
        DO I=1,IM
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
      SUBROUTINE CALRH_PW(RHPW)
!
! Algorithm use at GSD for RUC and Rapid Refresh                           
!------------------------------------------------------------------
!

      use vrbls3d, only: q, pmid, t
      use params_mod, only: g
      use ctlblk_mod, only: lm, jsta, jend, lm, im, spval
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none

      real,PARAMETER :: svp1=6.1153,svp2=17.67,svp3=29.65

      REAL, dimension(im,jsta:jend):: PW, PW_SAT, RHPW
      REAL deltp,sh,qv,temp,es,qs,qv_sat
      integer i,j,l,k,ka,kb

      pw     = 0.
      pw_sat = 0.
      rhpw   = 0.

      DO L=1,LM
        k=lm-l+1
       DO J=JSTA,JEND
        DO I=1,IM
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

        if (i==120 .and. j==120 )                        &
          write (6,*)'pw-sat', temp, sh, qs, pmid(i,j,kb)    &
          ,pmid(i,j,ka),pw(i,j),pw_sat(i,j)

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
!> @param[out] fpvsnew Real(krealfp) Saturation vapor pressure in Pascals.
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
      use ctlblk_mod, only: jsta_2l, jend_2u, lm, jsta, jend, im, me, spval
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
      integer, dimension(IM,Jsta:jend),intent(in)    :: L1D
      real,    dimension(IM,Jsta:jend),intent(in)    :: P1D,T1D
      real,    dimension(IM,jsta:jend),intent(inout) :: Q1D,CAPE,CINS,PPARC,ZEQL
!     
      integer, dimension(im,jsta:jend) :: IPTB, ITHTB, PARCEL, KLRES, KHRES, LCL, IDX
!     
      real,    dimension(im,jsta:jend) :: THESP, PSP, CAPE20, QQ, PP, THUND  
      REAL, ALLOCATABLE :: TPAR(:,:,:)

      LOGICAL THUNDER(IM,jsta:jend), NEEDTHUN 
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
      ALLOCATE(TPAR(IM,JSTA_2L:JEND_2U,LM))
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
        DO I=1,IM
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
          DO I=1,IM
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
          DO I=1,IM
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
            DO I=1,IM
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
          DO I=1,IM
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
          DO I=1,IM
            IF (PMID(I,J,L) < PSP(I,J))    LCL(I,J) = L+1
          ENDDO
        ENDDO
      ENDDO
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=1,IM
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
          DO I=1,IM
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
          CALL TTBLEX(TPAR(1,JSTA_2L,L),TTBL,ITB,JTB,KLRES             &
                    , PMID(1,JSTA_2L,L),PL,QQ,PP,RDP,THE0,STHE         &
                    , RDTHE,THESP,IPTB,ITHTB)
        ENDIF
!***
!***  COMPUTE PARCEL TEMPERATURE ALONG MOIST ADIABAT FOR PRESSURE>PLQ
!**
        IF(KNUMH > 0) THEN
          CALL TTBLEX(TPAR(1,JSTA_2L,L),TTBLQ,ITBQ,JTBQ,KHRES          &
                    , PMID(1,JSTA_2L,L),PLQ,QQ,PP,RDPQ                 &
                     ,THE0Q,STHEQ,RDTHEQ,THESP,IPTB,ITHTB)
        ENDIF

!------------SEARCH FOR EQ LEVEL----------------------------------------
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            IF(KHRES(I,J) > 0) THEN
              IF(TPAR(I,J,L) > T(I,J,L)) IEQL(I,J) = L
            ENDIF
          ENDDO
        ENDDO
!
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
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
        DO I=1,IM
          LBEG = MIN(IEQL(I,J),LBEG)
          LEND = MAX(LCL(I,J),LEND)
        ENDDO
      ENDDO
!
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=1,IM
          IF(T(I,J,IEQL(I,J)) > 255.65) THEN
            THUNDER(I,J) = .FALSE.
          ENDIF
        ENDDO
      ENDDO
!
      DO L=LBEG,LEND

!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            IDX(I,J) = 0
            IF(L >= IEQL(I,J).AND.L <= LCL(I,J)) THEN
              IDX(I,J) = 1
            ENDIF
          ENDDO
        ENDDO
!
!$omp  parallel do private(i,j,gdzkl,presk,thetaa,thetap,esatp,qsatp,tvp,tv)
        DO J=JSTA,JEND
          DO I=1,IM
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
        DO I=1,IM
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
      use ctlblk_mod, only: jsta_2l, jend_2u, lm, jsta, jend, im, jm, me, jsta_m, jend_m, spval
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
      integer, dimension(IM,Jsta:jend),intent(in)    :: L1D
      real,    dimension(IM,Jsta:jend),intent(in)    :: P1D,T1D
!      real,    dimension(IM,jsta:jend),intent(inout) :: Q1D,CAPE,CINS,PPARC,ZEQL
      real,    dimension(IM,jsta:jend),intent(inout) :: Q1D,CAPE,CINS
      real,    dimension(IM,jsta:jend)               :: PPARC,ZEQL
      real,    dimension(IM,jsta:jend),intent(inout) :: LFC,ESRHL,ESRHH
      real,    dimension(IM,jsta:jend),intent(inout) :: DCAPE,DGLD,ESP
      integer, dimension(im,jsta:jend) ::L12,L17,L3KM
!     
      integer, dimension(im,jsta:jend) :: IPTB, ITHTB, PARCEL, KLRES, KHRES, LCL, IDX
!     
      real,    dimension(im,jsta:jend) :: THESP, PSP, CAPE20, QQ, PP, THUND  
      integer, dimension(im,jsta:jend) :: PARCEL2 
      real,    dimension(im,jsta:jend) :: THESP2,PSP2
      real,    dimension(im,jsta:jend) :: CAPE4,CINS4 
      REAL, ALLOCATABLE :: TPAR(:,:,:)
      REAL, ALLOCATABLE :: TPAR2(:,:,:)

      LOGICAL THUNDER(IM,jsta:jend), NEEDTHUN 
      real PSFCK,PKL,TBTK,QBTK,APEBTK,TTHBTK,TTHK,APESPK,TPSPK,        &
           BQS00K,SQS00K,BQS10K,SQS10K,BQK,SQK,TQK,PRESK,GDZKL,THETAP, &
           THETAA,P00K,P10K,P01K,P11K,TTHESK,ESATP,QSATP,TVP,TV
      real PRESK2, ESATP2, QSATP2, TVP2, THETAP2, TV2, THETAA2
!      real,external :: fpvsnew
      integer I,J,L,KNUML,KNUMH,LBEG,LEND,IQ, KB,ITTBK
      integer IE,IW,JN,JS,IVE(JM),IVW(JM),JVN,JVS
      integer ISTART,ISTOP,JSTART,JSTOP
      real,    dimension(IM,jsta:jend) :: HTSFC

!     integer I,J,L,KNUML,KNUMH,LBEG,LEND,IQ,IT,LMHK, KB,ITTBK
!     
!**************************************************************
!     START CALCAPE HERE.
!     
      ALLOCATE(TPAR(IM,JSTA_2L:JEND_2U,LM))
      ALLOCATE(TPAR2(IM,JSTA_2L:JEND_2U,LM))
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
        DO I=1,IM
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
          DO I=1,IM
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
        ISTART = 2
        ISTOP  = IM-1
        JSTART = JSTA_M
        JSTOP  = JEND_M
      ELSE IF(gridtype == 'B')THEN
        JVN = 1
        JVS = 0
        do J=JSTA,JEND
          IVE(J)=1
          IVW(J)=0
        enddo
        ISTART = 2
        ISTOP  = IM-1
        JSTART = JSTA_M
        JSTOP  = JEND_M
      ELSE
        JVN = 0
        JVS = 0
        do J=JSTA,JEND
          IVE(J) = 0
          IVW(J) = 0
        enddo
        ISTART = 1
        ISTOP  = IM
        JSTART = JSTA
        JSTOP  = JEND
      END IF
!!$omp  parallel do private(htsfc,ie,iw)
      IF(gridtype /= 'A') CALL EXCH(FIS(1:IM,JSTA:JEND))
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
          DO I=1,IM
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
            DO I=1,IM
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
          DO I=1,IM
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
          DO I=1,IM
            IF (PMID(I,J,L) < PSP(I,J))    LCL(I,J) = L+1
          ENDDO
        ENDDO
      ENDDO
!$omp  parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=1,IM
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
          DO I=1,IM
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
          CALL TTBLEX(TPAR(1,JSTA_2L,L),TTBL,ITB,JTB,KLRES             &
                    , PMID(1,JSTA_2L,L),PL,QQ,PP,RDP,THE0,STHE         &
                    , RDTHE,THESP,IPTB,ITHTB)
        ENDIF
!***
!***  COMPUTE PARCEL TEMPERATURE ALONG MOIST ADIABAT FOR PRESSURE>PLQ
!**
        IF(KNUMH > 0) THEN
          CALL TTBLEX(TPAR(1,JSTA_2L,L),TTBLQ,ITBQ,JTBQ,KHRES          &
                    , PMID(1,JSTA_2L,L),PLQ,QQ,PP,RDPQ                 &
                     ,THE0Q,STHEQ,RDTHEQ,THESP,IPTB,ITHTB)
        ENDIF

!------------SEARCH FOR EQ LEVEL----------------------------------------
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            IF(KHRES(I,J) > 0) THEN
              IF(TPAR(I,J,L) > T(I,J,L)) IEQL(I,J) = L
            ENDIF
          ENDDO
        ENDDO
!
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
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
        DO I=1,IM
          LBEG = MIN(IEQL(I,J),LBEG)
          LEND = MAX(LCL(I,J),LEND)
        ENDDO
      ENDDO
!
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=1,IM
          IF(T(I,J,IEQL(I,J)) > 255.65) THEN
            THUNDER(I,J) = .FALSE.
          ENDIF
        ENDDO
      ENDDO
!
!reverse L order from bottom up for ESRH calculation
!
      ESRHH = LCL
      ESRHL = LCL
!      DO L=LBEG,LEND
      DO L=LEND,LBEG,-1

!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
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
          DO I=1,IM
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
          DO I=1,IM
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
        DO I=1,IM
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
          DO I=1,IM
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
          CALL TTBLEX(TPAR2(1,JSTA_2L,L),TTBL,ITB,JTB,KLRES             &
                    , PMID(1,JSTA_2L,L),PL,QQ,PP,RDP,THE0,STHE         &
                    , RDTHE,THESP2,IPTB,ITHTB)
        ENDIF
!***
!***  COMPUTE PARCEL TEMPERATURE ALONG MOIST ADIABAT FOR PRESSURE>PLQ
!**
        IF(KNUMH > 0) THEN
          CALL TTBLEX(TPAR2(1,JSTA_2L,L),TTBLQ,ITBQ,JTBQ,KHRES          &
                    , PMID(1,JSTA_2L,L),PLQ,QQ,PP,RDPQ                 &
                    , THE0Q,STHEQ,RDTHEQ,THESP2,IPTB,ITHTB)
        ENDIF
      ENDDO                  ! end of do l=lm,1,-1 loop

      LBEG = 1
      LEND = LM

      DO L=LBEG,LEND
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            IDX(I,J) = 0
            IF(L >= PARCEL2(I,J).AND.L < NINT(LMH(I,J))) THEN
              IDX(I,J) = 1
            ENDIF
          ENDDO
        ENDDO
!
!$omp  parallel do private(i,j,gdzkl,presk,thetaa,thetap,esatp,qsatp,tvp,tv)
        DO J=JSTA,JEND
          DO I=1,IM
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
        DO I=1,IM
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
          DO I=1,IM
            IF(T(I,J,L) <= TFRZ-12. .AND. L12(I,J)==LM) L12(I,J)=L
            IF(T(I,J,L) <= TFRZ-17. .AND. L17(I,J)==LM) L17(I,J)=L
          ENDDO
        ENDDO
      ENDDO
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=1,IM
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
          DO I=1,IM
            IF(ZINT(I,J,L)-HTSFC(I,J) <= 3000.) L3KM(I,J)=L
          ENDDO
        ENDDO
      ENDDO     
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=1,IM
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
! COMPUTE VIRTUAL TEMPERATURE
!
      IMPLICIT NONE
      REAL TVIRTUAL
      REAL, INTENT(IN) :: T, Q

      TVIRTUAL = T*(1+0.608*Q)

      end function TVIRTUAL
!
!-------------------------------------------------------------------------------------
!
  end module upp_physics
