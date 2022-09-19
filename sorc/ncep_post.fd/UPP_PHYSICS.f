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
  public :: CALDIV
  public :: CALGRADPS
  public :: CALRH
  public :: CALRH_GFS, CALRH_GSD, CALRH_NAM
  public :: CALRH_PW
  public :: CALSLR_ROEBBER
  public :: CALVOR

  public :: FPVSNEW
  public :: TVIRTUAL

  contains
!
!-------------------------------------------------------------------------------------
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
!      
      SUBROUTINE CALRH_GSD(P1,T1,Q1,RHB)
!
! Algorithm use at GSD for RUC and Rapid Refresh                           
!------------------------------------------------------------------
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
      SUBROUTINE CALRH_PW(RHPW)
!
! Algorithm use at GSD for RUC and Rapid Refresh                           
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
!> @file
!> @brief Subroutine that computes absolute vorticity.
!>
!> This routine computes the absolute vorticity.
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

!      print*,'dyval in CALVOR= ',DYVAL 
  
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

        if(me==0          ) print*,'CALVOR ',me,glatpoles(ista,1),glatpoles(ista,2)
        if(me==num_procs-1) print*,'CALVOR ',me,glatpoles(ista,1),glatpoles(ista,2)

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
!          if(ABSV(I,J)>1.0)print*,'Debug CALVOR',i,j,VWND(ip1,J),VWND(im1,J), &
!          wrk2(i,j),UWND(I,J-1),COSL(I,J-1),UWND(I,J+1),COSL(I,J+1),wrk3(i,j),cosl(i,j),F(I,J),ABSV(I,J)
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
!sk06132016
              if(DIV(I,J,l)>1.0)print*,'Debug in CALDIV',i,j,UWND(ip1,J,l),UWND(im1,J,l), &
     &           wrk2(i,j),VWND(I,J-1,l),COSL(I,J-1),VWND(I,J+1,l),COSL(I,J+1),         &
     &           wrk3(i,j),wrk1(i,j),DIV(I,J,l)
!--
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

!sk06142016e
        if(DIV(ista,jsta,l)>1.0)print*,'Debug in CALDIV',jsta,DIV(ista,jsta,l)
!       print*,'Debug in CALDIV',' jsta= ',jsta,DIV(1,jsta,l)

      enddo                        ! end of l looop
!--
      deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)
     

      END SUBROUTINE CALDIV

      SUBROUTINE CALGRADPS(PS,PSX,PSY)
!> CALGRADPS computes gardients of a scalar field PS or LNPS.
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
!sk06142016A
              if(PSX(I,J)>100.0)print*,'Debug in CALGRADPS: PSX',i,j,PS(ip1,J),PS(im1,J), &
!             print*,'Debug in CALGRADPS',i,j,PS(ip1,J),PS(im1,J), &
     &           wrk2(i,j),wrk1(i,j),PSX(I,J)
              if(PSY(I,J)>100.0)print*,'Debug in CALGRADPS: PSY',i,j,PS(i,J-1),PS(i,J+1), &
!             print*,'Debug in CALGRADPS',i,j,PS(i,J-1),PS(i,J+1), &
     &           wrk3(i,j),ERAD,PSY(I,J)
!--
            ENDDO
          END IF
!
        ENDDO                               ! end of J loop

        deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)
     
!     END IF 

      END SUBROUTINE CALGRADPS

!> calslr_roebber() computes snow solid-liquid-ratio slr using the Roebber algorithm
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2022-07-11 | Jesse Meng | Initial
!>
!> @author Jesse Meng @date 2022-07-11

      SUBROUTINE CALSLR_ROEBBER(sno,si,slr)

      use vrbls3d, only: T, Q, PMID
      use vrbls2d, only: slp, prec, u10, v10
      use ctlblk_mod, only: ista, iend, jsta, jend, LM, spval

      implicit none

      real,dimension(ista:iend,jsta:jend),intent(in)    :: sno !weasd
      real,dimension(ista:iend,jsta:jend),intent(in)    :: si  !snod
      real,dimension(ista:iend,jsta:jend),intent(out)   :: slr !1/sndens=si/sno

! local variables

      real,dimension(ista:iend,jsta:jend)    :: P1D,T1D,Q1D,EGRID1
      real,dimension(ista:iend,jsta:jend,lm) :: RH3D

      type all_grids
           real :: grid
           real :: sigma
      end type all_grids

      real,dimension(0:14), parameter        :: sig = &
       (/0.0, 1.0, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85,&
         0.8, 0.75, 0.7, 0.65, 0.6, 0.5, 0.4/)

      real,dimension(12), parameter          :: mf = &
       (/1.0, 0.67, 0.33, 0.0, -0.33, -0.67, -1.00, -0.67, -0.33, 0.0, 0.33, 0.67/)

      real,dimension(0:14)                   :: tm, rhm

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

      type(all_grids), dimension(ista:iend,jsta:jend,lm) :: tmpk_grids, rh_grids
      integer,         dimension(ista:iend,jsta:jend,lm) :: tmpk_levels, rh_levels

      real,dimension(ista:iend,jsta:jend)    :: hprob,mprob,lprob
      real,dimension(ista:iend,jsta:jend)    :: slrgrid, slrgrid2
      real,dimension(ista:iend,jsta:jend)    :: pres,qpf,swnd

      character*20 nswFileName
      real :: psurf, sgw, sg1, sg2, dtds, rhds
      real :: f1, f2, f3, f4, f5, f6
      real :: p1, p2, p3
      real :: hprob_tot = 0.
      real :: mprob_tot = 0.
      real :: lprob_tot = 0.

      integer :: i, j, k, ks, L, LL, imo
!
!***************************************************************************
!
! calculate rh for all levels

      loop_lm: DO L=1,LM
        LL=LM-L+1
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
        DO I=ISTA,IEND
          P1D(I,J) = PMID(I,J,LL)
          T1D(I,J) = T(I,J,LL)
          Q1D(I,J) = Q(I,J,LL)
        ENDDO
        ENDDO
        CALL CALRH(P1D,T1D,Q1D,EGRID1)
        RH3D(:,:,LL)=EGRID1
      END DO loop_lm

! Load variables

      DO L=1,LM
        LL=LM-L+1
!$omp parallel do private(i,j)
      do j=jsta,jend
      do i=ista,iend
         tmpk_grids(i,j,LL)%grid=T(I,J,LL)-273.15
         tmpk_levels(i,j,LL)=PMID(I,J,LL)
         rh_grids(i,j,LL)%grid=RH3D(I,J,LL)
         rh_levels(i,j,LL)=PMID(I,J,LL)
      end do
      end do
      END DO

!$omp parallel do private(i,j)
      DO J=JSTA,JEND
      DO I=ISTA,IEND
         pres(i,j)=slp(i,j)
         qpf(i,j)=prec(i,j)
         swnd(i,j)=spval
         if(u10(i,j)/=spval .and. v10(i,j)/=spval) &
         swnd(i,j)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
      END DO
      END DO

! Convert to sigma

      tmpk_grids(:,:,1)%sigma = 1
      rh_grids(:,:,1)%sigma = 1

      DO L=1,LM-1
        LL=LM-L+1
!$omp parallel do private(i,j)
        do j=jsta,jend
        do i=ista,iend
           if(pres(i,j) == spval) then
              tmpk_grids(i,j,LL)%sigma=spval
              rh_grids(i,j,LL)%sigma=spval
           else
              tmpk_grids(i,j,LL)%sigma=tmpk_levels(i,j,LL) / pres(i,j)     
              rh_grids(i,j,LL)%sigma=rh_levels(i,j,LL) / pres(i,j)
           endif
        end do
        end do
      END DO

! main slr i/j loop start

!$omp parallel do private(i,j)
      loop_slr: do j=jsta,jend
      do i=ista,iend
         slr(i,j)=spval
        ! if(sno(i,j) /= spval .and. si(i,j) /= spval .and. si(i,j) > 0.) then
        !   slr(i,j) = si(i,j)/sno(i,j)
        ! endif
        ! slr(i,j) = RH3D(i,j,LM)

      if(pres(i,j)/=spval .and. qpf(i,j)/=spval .and. swnd(i,j)/=spval) then

! Interpolate to the 14 sigma levels      

      loop_ks15: do ks=1,14
         psurf = pres(i,j)
         sgw   = sig(ks)

         do L=LM,2,-1
           LL=LM-L+1
           if(LL==1) then
              sg1 = psurf/psurf
           else
              sg1 = tmpk_levels(i,j,LL)/psurf
           endif  
              sg2 = tmpk_levels(i,j,LL+1)/psurf
           
           if(sg1==sgw) then
              tm(ks) = tmpk_grids(i,j,LL)%grid
              rhm(ks)=   rh_grids(i,j,LL)%grid
           elseif (sg2==sgw) then
              tm(ks) = tmpk_grids(i,j,LL+1)%grid
              rhm(ks)=   rh_grids(i,j,LL+1)%grid
           elseif ((sgw < sg1) .and. (sgw > sg2)) then
              dtds = (tmpk_grids(i,j,LL+1)%grid - tmpk_grids(i,j,LL)%grid) / (sg2-sg1)
              tm(ks) = ((sgw - sg1) * dtds) + tmpk_grids(i,j,LL)%grid
              rhds = (rh_grids(i,j,LL+1)%grid - rh_grids(i,j,LL)%grid) / (sg2-sg1)
              rhm(ks)= ((sgw - sg1) * rhds) + rh_grids(i,j,LL)%grid
           endif    
         end do
      end do loop_ks15

! Have surface wind, QPF, and temp/RH on the 14 levels.
! Convert these data to the factors using regression equations

      f1 = co1(1)+co1(2)*qpf(i,j)+co1(3)*swnd(i,j)+co1(4)*tm(1)+co1(5)*tm(2)+co1(6)*tm(3)+ &
           co1(7)*tm(4)+co1(8)*tm(5)+co1(9)*tm(6)+co1(10)*tm(7)+co1(11)*tm(8)+ &
           co1(12)*tm(9)+co1(13)*tm(10)+co1(14)*tm(11)+co1(15)*tm(12)+co1(16)*tm(13)+ &
           co1(17)*tm(14)+co1(18)*rhm(1)+co1(19)*rhm(2)+co1(20)*rhm(3)+co1(21)*rhm(4)+ &
           co1(22)*rhm(5)+co1(23)*rhm(6)+co1(24)*rhm(7)+co1(25)*rhm(8)+co1(26)*rhm(9)+ &
           co1(27)*rhm(10)+co1(28)*rhm(11)+co1(29)*rhm(12)+co1(30)*rhm(13);

      f2 = co2(1)+co2(2)*qpf(i,j)+co2(3)*swnd(i,j)+co2(4)*tm(1)+co2(5)*tm(2)+co2(6)*tm(3)+ &
           co2(7)*tm(4)+co2(8)*tm(5)+co2(9)*tm(6)+co2(10)*tm(7)+co2(11)*tm(8)+ &
           co2(12)*tm(9)+co2(13)*tm(10)+co2(14)*tm(11)+co2(15)*tm(12)+co2(16)*tm(13)+ &
           co2(17)*tm(14)+co2(18)*rhm(1)+co2(19)*rhm(2)+co2(20)*rhm(3)+co2(21)*rhm(4)+ &
           co2(22)*rhm(5)+co2(23)*rhm(6)+co2(24)*rhm(7)+co2(25)*rhm(8)+co2(26)*rhm(9)+ &
           co2(27)*rhm(10)+co2(28)*rhm(11)+co2(29)*rhm(12)+co2(30)*rhm(13);

      f3 = co3(1)+co3(2)*qpf(i,j)+co3(3)*swnd(i,j)+co3(4)*tm(1)+co3(5)*tm(2)+co3(6)*tm(3)+ &
           co3(7)*tm(4)+co3(8)*tm(5)+co3(9)*tm(6)+co3(10)*tm(7)+co3(11)*tm(8)+ &
           co3(12)*tm(9)+co3(13)*tm(10)+co3(14)*tm(11)+co3(15)*tm(12)+co3(16)*tm(13)+ &
           co3(17)*tm(14)+co3(18)*rhm(1)+co3(19)*rhm(2)+co3(20)*rhm(3)+co3(21)*rhm(4)+ &
           co3(22)*rhm(5)+co3(23)*rhm(6)+co3(24)*rhm(7)+co3(25)*rhm(8)+co3(26)*rhm(9)+ &
           co3(27)*rhm(10)+co3(28)*rhm(11)+co3(29)*rhm(12)+co3(30)*rhm(13);

      f4 = co4(1)+co4(2)*qpf(i,j)+co4(3)*swnd(i,j)+co4(4)*tm(1)+co4(5)*tm(2)+co4(6)*tm(3)+ &
           co4(7)*tm(4)+co4(8)*tm(5)+co4(9)*tm(6)+co4(10)*tm(7)+co4(11)*tm(8)+ &
           co4(12)*tm(9)+co4(13)*tm(10)+co4(14)*tm(11)+co4(15)*tm(12)+co4(16)*tm(13)+ &
           co4(17)*tm(14)+co4(18)*rhm(1)+co4(19)*rhm(2)+co4(20)*rhm(3)+co4(21)*rhm(4)+ &
           co4(22)*rhm(5)+co4(23)*rhm(6)+co4(24)*rhm(7)+co4(25)*rhm(8)+co4(26)*rhm(9)+ &
           co4(27)*rhm(10)+co4(28)*rhm(11)+co4(29)*rhm(12)+co4(30)*rhm(13);

      f5 = co5(1)+co5(2)*qpf(i,j)+co5(3)*swnd(i,j)+co5(4)*tm(1)+co5(5)*tm(2)+co5(6)*tm(3)+ &
           co5(7)*tm(4)+co5(8)*tm(5)+co5(9)*tm(6)+co5(10)*tm(7)+co5(11)*tm(8)+ &
           co5(12)*tm(9)+co5(13)*tm(10)+co5(14)*tm(11)+co5(15)*tm(12)+co5(16)*tm(13)+ &
           co5(17)*tm(14)+co5(18)*rhm(1)+co5(19)*rhm(2)+co5(20)*rhm(3)+co5(21)*rhm(4)+ &
           co5(22)*rhm(5)+co5(23)*rhm(6)+co5(24)*rhm(7)+co5(25)*rhm(8)+co5(26)*rhm(9)+ &
           co5(27)*rhm(10)+co5(28)*rhm(11)+co5(29)*rhm(12)+co5(30)*rhm(13);
      
      f6 = co6(1)+co6(2)*qpf(i,j)+co6(3)*swnd(i,j)+co6(4)*tm(1)+co6(5)*tm(2)+co6(6)*tm(3)+ &
           co6(7)*tm(4)+co6(8)*tm(5)+co6(9)*tm(6)+co6(10)*tm(7)+co6(11)*tm(8)+ &
           co6(12)*tm(9)+co6(13)*tm(10)+co6(14)*tm(11)+co6(15)*tm(12)+co6(16)*tm(13)+ &
           co6(17)*tm(14)+co6(18)*rhm(1)+co6(19)*rhm(2)+co6(20)*rhm(3)+co6(21)*rhm(4)+ &
           co6(22)*rhm(5)+co6(23)*rhm(6)+co6(24)*rhm(7)+co6(25)*rhm(8)+co6(26)*rhm(9)+ &
           co6(27)*rhm(10)+co6(28)*rhm(11)+co6(29)*rhm(12)+co6(30)*rhm(13);

      do k=1,10
         p1 = 0
         p2 = 0
         p3 = 0
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
         hprob_tot = hprob_tot + p1
         mprob_tot = mprob_tot + p2
         lprob_tot = lprob_tot + p3
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
         slrgrid2(i,j) = hprob(i,j)*8.0 + mprob(i,j)*13.0 + lprob(i,j)*18.0
         slrgrid2(i,j) = slrgrid2(i,j)/(hprob(i,j)+mprob(i,j)+lprob(i,j))
      else
         slrgrid2(i,j) = hprob(i,j)*8.0 + mprob(i,j)*13.0 + lprob(i,j)*27.0
         slrgrid2(i,j) = slrgrid2(i,j)/(hprob(i,j)+mprob(i,j)+lprob(i,j))
      endif
               
      slr(i,j) = slrgrid2(i,j)
      slr(i,j) = RH3D(i,j,LM)

      endif !if(pres(i,j), qpf(i,j), swnd(i,j) /= spval)
      enddo
      enddo loop_slr

      END SUBROUTINE CALSLR_ROEBBER
!
!-------------------------------------------------------------------------------------
!
      SUBROUTINE breadboard1_main(nswFileName,mf,f1,f2,f3,f4,f5,f6,p1,p2,p3)

      implicit none

      character*20 nswFileName
      real mf, f1, f2, f3, f4, f5, f6
      real p1, p2, p3

      integer ieof
      character*100 bbstring
      integer datacount
      real inputFile(2,7)
      real inputAxon(7)
      real hidden1Axon(40)
      real outputAxon(3)
      real hidden1Synapse(7,40)
      real outputSynapse(40,3)
      real activeOutputProbe(2,3)

      integer i,j

      open(11,file=nswFileName,status='unknown')

      ieof = 0
      do while (ieof == 0)

      read(11,'(a)',iostat=ieof) bbstring

      if(trim(bbstring)=='#inputFile File') then
         print*,trim(bbstring)
         read(11,*) datacount
         do j=1,7
            read(11,*) inputFile(:,j)
         enddo
      endif

      if(trim(bbstring)=='#inputAxon Axon') then
         read(11,*) i,j
         read(11,*) j
      endif

      if(trim(bbstring)=='#hidden1Axon TanhAxon') then
         read(11,*) i,j
         read(11,*) j
         read(11,*) datacount, hidden1Axon
      endif

      if(trim(bbstring)=='#outputAxon SoftMaxAxon') then
         read(11,*) i,j
         read(11,*) j
         read(11,*) datacount, outputAxon
      endif

      if(trim(bbstring)=='#hidden1Synapse FullSynapse') then
         read(11,*) datacount, hidden1Synapse
      endif

      if(trim(bbstring)=='#outputSynapse FullSynapse') then
         read(11,*) datacount, outputSynapse
      endif

      if(trim(bbstring)=='#activeOutputProbe DataWriter') then
         read(11,*) datacount
         do j=1,3
            read(11,*) activeOutputProbe(:,j)
         enddo
      endif

      enddo

      close(11)

      END SUBROUTINE breadboard1_main
!
!-------------------------------------------------------------------------------------
!
      SUBROUTINE breadboard6_main(nswFileName,mf,f1,f2,f3,f4,f5,f6,p1,p2,p3)

      character*20 nswFileName
      real mf, f1, f2, f3, f4, f5, f6
      real p1, p2, p3

      integer ieof
      character*100 bbstring
      integer datacount
      real inputFile(2,7)
      real inputAxon(7)
      real hidden1Axon(40)
      real outputAxon(3)
      real hidden1Synapse(7,40)
      real outputSynapse(40,3)
      real activeOutputProbe(2,3)

      integer i,j

      open(11,file=nswFileName,status='unknown')

      ieof = 0
      do while (ieof == 0)

      read(11,'(a)',iostat=ieof) bbstring

      if(trim(bbstring)=='#inputFile File') then
         print*,trim(bbstring)
         read(11,*) datacount
         do j=1,7
            read(11,*) inputFile(:,j)
         enddo
      endif

      if(trim(bbstring)=='#inputAxon Axon') then
         read(11,*) i,j
         read(11,*) j
      endif

      if(trim(bbstring)=='#hidden1Axon TanhAxon') then
         read(11,*) i,j
         read(11,*) j
         read(11,*) datacount, hidden1Axon
      endif

      if(trim(bbstring)=='#outputAxon SoftMaxAxon') then
         read(11,*) i,j
         read(11,*) j
         read(11,*) datacount, outputAxon
      endif

      if(trim(bbstring)=='#hidden1Synapse FullSynapse') then
         read(11,*) datacount, hidden1Synapse
      endif

      if(trim(bbstring)=='#outputSynapse FullSynapse') then
         read(11,*) datacount, outputSynapse
      endif

      if(trim(bbstring)=='#activeOutputProbe DataWriter') then
         read(11,*) datacount
         do j=1,3
            read(11,*) activeOutputProbe(:,j)
         enddo
      endif

      enddo

      close(11)

      END SUBROUTINE breadboard6_main
!
!-------------------------------------------------------------------------------------
!
  end module upp_physics

