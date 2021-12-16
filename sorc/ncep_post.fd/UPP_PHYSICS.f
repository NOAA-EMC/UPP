!> @file
!
!> SUBPROGRAM: UPP_PHYSICS
!! @author JMENG @date 2020-05-20
!!
!! A collection of UPP subroutines for physics variables calculation.
!!
!! CALCAPE
!! Compute CAPE/CINS and other storm related variables.
!!
!! CALCAPE2
!! Compute additional storm related variables.
!!
!! CALRH
!! CALRH_NAM
!! CALRH_GFS
!! CALRH_GSD
!! Compute RH using various algorithms.
!! The NAM v4.1.18 ALGORITHM (CALRH_NAM) is selected as default for 
!! NMMB and FV3GFS, FV3GEFS, and FV3R for the UPP 2020 unification.
!!
!! CALRH_PW
!! Algorithm use at GSD for RUC and Rapid Refresh
!!
!! FPVSNEW
!! Compute saturation vapor pressure.
!!
!! TVIRTUAL
!! Compute virtual temperature.
!!
!! PROGRAM HISTORY LOG:
!!   MAY, 2020    Jesse Meng   Initial code
!!-------------------------------------------------------------------------------------
!!
  module upp_physics

  implicit none

  private

  public :: CALCAPE, CALCAPE2
  public :: CALDIV
  public :: CALGRADPS
  public :: CALRH
  public :: CALRH_GFS, CALRH_GSD, CALRH_NAM
  public :: CALRH_PW
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
      SUBROUTINE CALRH_NAM(P1,T1,Q1,RH)
!      SUBROUTINE CALRH(P1,T1,Q1,RH)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALRH       COMPUTES RELATIVE HUMIDITY
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-22       
!     
! ABSTRACT:  
!     THIS ROUTINE COMPUTES RELATIVE HUMIDITY GIVEN PRESSURE, 
!     TEMPERATURE, SPECIFIC HUMIDITY. AN UPPER AND LOWER BOUND
!     OF 100 AND 1 PERCENT RELATIVE HUMIDITY IS ENFORCED.  WHEN
!     THESE BOUNDS ARE APPLIED THE PASSED SPECIFIC HUMIDITY 
!     ARRAY IS ADJUSTED AS NECESSARY TO PRODUCE THE SET RELATIVE
!     HUMIDITY.
!   .     
!     
! PROGRAM HISTORY LOG:
!   ??-??-??  DENNIS DEAVEN
!   92-12-22  RUSS TREADON - MODIFIED AS DESCRIBED ABOVE.
!   98-06-08  T BLACK      - CONVERSION FROM 1-D TO 2-D
!   98-08-18  MIKE BALDWIN - MODIFY TO COMPUTE RH OVER ICE AS IN MODEL
!   98-12-16  GEOFF MANIKIN - UNDO RH COMPUTATION OVER ICE
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   02-06-11  MIKE BALDWIN - WRF VERSION
!   06-03-19  Wen Meng     - MODIFY TOP PRESSURE to 1 PA
!     
! USAGE:    CALL CALRH(P1,T1,Q1,RH)
!   INPUT ARGUMENT LIST:
!     P1     - PRESSURE (PA)
!     T1     - TEMPERATURE (K)
!     Q1     - SPECIFIC HUMIDITY (KG/KG)
!
!   OUTPUT ARGUMENT LIST: 
!     RH     - RELATIVE HUMIDITY  (DECIMAL FORM)
!     Q1     - ADJUSTED SPECIFIC HUMIDITY (KG/KG)
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!     LIBRARY:
!       NONE
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
!
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
      SUBROUTINE CALRH_GFS(P1,T1,Q1,RH)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALRH       COMPUTES RELATIVE HUMIDITY
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-22       
!     
! ABSTRACT:  
!     THIS ROUTINE COMPUTES RELATIVE HUMIDITY GIVEN PRESSURE, 
!     TEMPERATURE, SPECIFIC HUMIDITY. AN UPPER AND LOWER BOUND
!     OF 100 AND 1 PERCENT RELATIVE HUMIDITY IS ENFORCED.  WHEN
!     THESE BOUNDS ARE APPLIED THE PASSED SPECIFIC HUMIDITY 
!     ARRAY IS ADJUSTED AS NECESSARY TO PRODUCE THE SET RELATIVE
!     HUMIDITY.
!   .     
!     
! PROGRAM HISTORY LOG:
!   ??-??-??  DENNIS DEAVEN
!   92-12-22  RUSS TREADON - MODIFIED AS DESCRIBED ABOVE.
!   98-06-08  T BLACK      - CONVERSION FROM 1-D TO 2-D
!   98-08-18  MIKE BALDWIN - MODIFY TO COMPUTE RH OVER ICE AS IN MODEL
!   98-12-16  GEOFF MANIKIN - UNDO RH COMPUTATION OVER ICE
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   02-06-11  MIKE BALDWIN - WRF VERSION
!   13-08-13  S. Moorthi   - Threading
!   06-03-19  Wen Meng     - MODIFY TOP PRESSURE to 1 PA
!     
! USAGE:    CALL CALRH(P1,T1,Q1,RH)
!   INPUT ARGUMENT LIST:
!     P1     - PRESSURE (PA)
!     T1     - TEMPERATURE (K)
!     Q1     - SPECIFIC HUMIDITY (KG/KG)
!
!   OUTPUT ARGUMENT LIST: 
!     RH     - RELATIVE HUMIDITY  (DECIMAL FORM)
!     Q1     - ADJUSTED SPECIFIC HUMIDITY (KG/KG)
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!     LIBRARY:
!       NONE
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
!
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
!$$$     Subprogram Documentation Block
!
! Subprogram: fpvsnew         Compute saturation vapor pressure
!   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
!
! Abstract: Compute saturation vapor pressure from the temperature.
!   A linear interpolation is done between values in a lookup table
!   computed in gpvs. See documentation for fpvsx for details.
!   Input values outside table range are reset to table extrema.
!   The interpolation accuracy is almost 6 decimal places.
!   On the Cray, fpvs is about 4 times faster than exact calculation.
!   This function should be expanded inline in the calling routine.
!
! Program History Log:
!   91-05-07  Iredell             made into inlinable function
!   94-12-30  Iredell             expand table
! 1999-03-01  Iredell             f90 module
! 2001-02-26  Iredell             ice phase
!
! Usage:   pvs=fpvsnew(t)
!
!   Input argument list:
!     t          Real(krealfp) temperature in Kelvin
!
!   Output argument list:
!     fpvsnew       Real(krealfp) saturation vapor pressure in Pascals
!
! Attributes:
!   Language: Fortran 90.
!
!$$$
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
!

      SUBROUTINE CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,L1D,CAPE,    &  
                         CINS,PPARC,ZEQL,THUND)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALCAPE     COMPUTES CAPE AND CINS
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-02-10      
!     
! ABSTRACT:  
!     
!     THIS ROUTINE COMPUTES CAPE AND CINS GIVEN TEMPERATURE,
!     PRESSURE, AND SPECIFIC HUMIDTY.  IN "STORM AND CLOUD 
!     DYNAMICS" (1989, ACADEMIC PRESS) COTTON AND ANTHES DEFINE
!     CAPE (EQUATION 9.16, P501) AS
!
!                  EL
!         CAPE =  SUM G * LN(THETAP/THETAA) DZ 
!                 LCL
!     
!     WHERE,
!      EL    = EQUILIBRIUM LEVEL,
!     LCL    = LIFTING CONDENSTATION LEVEL,
!       G    = GRAVITATIONAL ACCELERATION,
!     THETAP = LIFTED PARCEL POTENTIAL TEMPERATURE,
!     THETAA = AMBIENT POTENTIAL TEMPERATURE.
!     
!     NOTE THAT THE INTEGRAND LN(THETAP/THETAA) APPROXIMATELY
!     EQUALS (THETAP-THETAA)/THETAA.  THIS RATIO IS OFTEN USED
!     IN THE DEFINITION OF CAPE/CINS.
!     
!     TWO TYPES OF CAPE/CINS CAN BE COMPUTED BY THIS ROUTINE.  THE
!     SUMMATION PROCESS IS THE SAME FOR BOTH CASES.  WHAT DIFFERS
!     IS THE DEFINITION OF THE PARCEL TO LIFT.  FOR ITYPE=1 THE
!     PARCEL WITH THE WARMEST THETA-E IN A DPBND PASCAL LAYER ABOVE
!     THE MODEL SURFACE IS LIFTED.  THE ARRAYS P1D, T1D, AND Q1D
!     ARE NOT USED.  FOR ITYPE=2 THE ARRAYS P1D, T1D, AND Q1D
!     DEFINE THE PARCEL TO LIFT IN EACH COLUMN.  BOTH TYPES OF
!     CAPE/CINS MAY BE COMPUTED IN A SINGLE EXECUTION OF THE POST
!     PROCESSOR.
!     
!     THIS ALGORITHM PROCEEDS AS FOLLOWS.
!     FOR EACH COLUMN, 
!        (1)  INITIALIZE RUNNING CAPE AND CINS SUM TO 0.0
!        (2)  COMPUTE TEMPERATURE AND PRESSURE AT THE LCL USING
!             LOOK UP TABLE (PTBL).  USE EITHER PARCEL THAT GIVES
!             MAX THETAE IN LOWEST DPBND ABOVE GROUND (ITYPE=1)
!             OR GIVEN PARCEL FROM T1D,Q1D,...(ITYPE=2).
!        (3)  COMPUTE THE TEMP OF A PARCEL LIFTED FROM THE LCL.
!             WE KNOW THAT THE PARCEL'S
!             EQUIVALENT POTENTIAL TEMPERATURE (THESP) REMAINS
!             CONSTANT THROUGH THIS PROCESS.  WE CAN
!             COMPUTE TPAR USING THIS KNOWLEDGE USING LOOK
!             UP TABLE (SUBROUTINE TTBLEX).
!        (4)  FIND THE EQUILIBRIUM LEVEL.  THIS IS DEFINED AS THE
!             HIGHEST POSITIVELY BUOYANT LAYER.
!             (IF THERE IS NO POSITIVELY BUOYANT LAYER, CAPE/CINS
!              WILL BE ZERO)
!        (5)  COMPUTE CAPE/CINS.  
!             (A) COMPUTE THETAP.  WE KNOW TPAR AND P.
!             (B) COMPUTE THETAA.  WE KNOW T AND P.  
!        (6)  ADD G*(THETAP-THETAA)*DZ TO THE RUNNING CAPE OR CINS SUM.
!             (A) IF THETAP > THETAA, ADD TO THE CAPE SUM.
!             (B) IF THETAP < THETAA, ADD TO THE CINS SUM.
!        (7)  ARE WE AT EQUILIBRIUM LEVEL? 
!             (A) IF YES, STOP THE SUMMATION.
!             (B) IF NO, CONTIUNUE THE SUMMATION.
!        (8)  ENFORCE LIMITS ON CAPE AND CINS (I.E. NO NEGATIVE CAPE)
!     
! PROGRAM HISTORY LOG:
!   93-02-10  RUSS TREADON
!   93-06-19  RUSS TREADON - GENERALIZED ROUTINE TO ALLOW FOR 
!                            TYPE 2 CAPE/CINS CALCULATIONS.     
!   94-09-23  MIKE BALDWIN - MODIFIED TO USE LOOK UP TABLES
!                            INSTEAD OF COMPLICATED EQUATIONS.
!   94-10-13  MIKE BALDWIN - MODIFIED TO CONTINUE CAPE/CINS CALC
!                            UP TO AT HIGHEST BUOYANT LAYER.
!   98-06-12  T BLACK      - CONVERSION FROM 1-D TO 2-D
!   98-08-18  T BLACK      - COMPUTE APE INTERNALLY
!   00-01-04  JIM TUCCILLO - MPI VERSION              
!   02-01-15  MIKE BALDWIN - WRF VERSION
!   03-08-24  G MANIKIN    - ADDED LEVEL OF PARCEL BEING LIFTED
!                            AS OUTPUT FROM THE ROUTINE AND ADDED
!                            THE DEPTH OVER WHICH ONE SEARCHES FOR
!                            THE MOST UNSTABLE PARCEL AS INPUT
!   10-09-09  G MANIKIN    - CHANGED COMPUTATION TO USE VIRTUAL TEMP
!                          - ADDED EQ LVL HGHT AND THUNDER PARAMETER    
!   15-xx-xx  S MOORTHI    - optimization and threading
!   21-07-28  W Meng       - Restrict computation from undefined grids.
!   21-09-01  E COLON      - equivalent level height index for RTMA
!
! USAGE:    CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,L1D,CAPE,
!                                CINS,PPARC)
!   INPUT ARGUMENT LIST:
!     ITYPE    - INTEGER FLAG SPECIFYING HOW PARCEL TO LIFT IS
!                IDENTIFIED.  SEE COMMENTS ABOVE.
!     DPBND    - DEPTH OVER WHICH ONE SEARCHES FOR MOST UNSTABLE PARCEL
!     P1D      - ARRAY OF PRESSURE OF PARCELS TO LIFT.
!     T1D      - ARRAY OF TEMPERATURE OF PARCELS TO LIFT.
!     Q1D      - ARRAY OF SPECIFIC HUMIDITY OF PARCELS TO LIFT.
!     L1D      - ARRAY OF MODEL LEVEL OF PARCELS TO LIFT.
!
!   OUTPUT ARGUMENT LIST: 
!     CAPE     - CONVECTIVE AVAILABLE POTENTIAL ENERGY (J/KG)
!     CINS     - CONVECTIVE INHIBITION (J/KG)
!     PPARC    - PRESSURE LEVEL OF PARCEL LIFTED WHEN ONE SEARCHES
!                  OVER A PARTICULAR DEPTH TO COMPUTE CAPE/CIN
!     
!   OUTPUT FILES:
!     STDOUT  - RUN TIME STANDARD OUT.
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       BOUND    - BOUND (CLIP) DATA BETWEEN UPPER AND LOWER LIMTS.
!       TTBLEX   - LOOKUP TABLE ROUTINE TO GET T FROM THETAE AND P
!
!     LIBRARY:
!       COMMON   - 
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN 90
!     MACHINE : CRAY C-90
!$$$  
!
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
!
      SUBROUTINE CALCAPE2(ITYPE,DPBND,P1D,T1D,Q1D,L1D,    &  
                          CAPE,CINS,LFC,ESRHL,ESRHH,      &
                          DCAPE,DGLD,ESP)
!      SUBROUTINE CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,L1D,CAPE,    &  
!                         CINS,PPARC,ZEQL,THUND)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALCAPE     COMPUTES CAPE AND CINS
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-02-10      
!     
! ABSTRACT:  
!     
!     THIS ROUTINE COMPUTES CAPE AND CINS GIVEN TEMPERATURE,
!     PRESSURE, AND SPECIFIC HUMIDTY.  IN "STORM AND CLOUD 
!     DYNAMICS" (1989, ACADEMIC PRESS) COTTON AND ANTHES DEFINE
!     CAPE (EQUATION 9.16, P501) AS
!
!                  EL
!         CAPE =  SUM G * LN(THETAP/THETAA) DZ 
!                 LCL
!     
!     WHERE,
!      EL    = EQUILIBRIUM LEVEL,
!     LCL    = LIFTING CONDENSTATION LEVEL,
!       G    = GRAVITATIONAL ACCELERATION,
!     THETAP = LIFTED PARCEL POTENTIAL TEMPERATURE,
!     THETAA = AMBIENT POTENTIAL TEMPERATURE.
!     
!     NOTE THAT THE INTEGRAND LN(THETAP/THETAA) APPROXIMATELY
!     EQUALS (THETAP-THETAA)/THETAA.  THIS RATIO IS OFTEN USED
!     IN THE DEFINITION OF CAPE/CINS.
!     
!     TWO TYPES OF CAPE/CINS CAN BE COMPUTED BY THIS ROUTINE.  THE
!     SUMMATION PROCESS IS THE SAME FOR BOTH CASES.  WHAT DIFFERS
!     IS THE DEFINITION OF THE PARCEL TO LIFT.  FOR ITYPE=1 THE
!     PARCEL WITH THE WARMEST THETA-E IN A DPBND PASCAL LAYER ABOVE
!     THE MODEL SURFACE IS LIFTED.  THE ARRAYS P1D, T1D, AND Q1D
!     ARE NOT USED.  FOR ITYPE=2 THE ARRAYS P1D, T1D, AND Q1D
!     DEFINE THE PARCEL TO LIFT IN EACH COLUMN.  BOTH TYPES OF
!     CAPE/CINS MAY BE COMPUTED IN A SINGLE EXECUTION OF THE POST
!     PROCESSOR.
!     
!     THIS ALGORITHM PROCEEDS AS FOLLOWS.
!     FOR EACH COLUMN, 
!        (1)  INITIALIZE RUNNING CAPE AND CINS SUM TO 0.0
!        (2)  COMPUTE TEMPERATURE AND PRESSURE AT THE LCL USING
!             LOOK UP TABLE (PTBL).  USE EITHER PARCEL THAT GIVES
!             MAX THETAE IN LOWEST DPBND ABOVE GROUND (ITYPE=1)
!             OR GIVEN PARCEL FROM T1D,Q1D,...(ITYPE=2).
!        (3)  COMPUTE THE TEMP OF A PARCEL LIFTED FROM THE LCL.
!             WE KNOW THAT THE PARCEL'S
!             EQUIVALENT POTENTIAL TEMPERATURE (THESP) REMAINS
!             CONSTANT THROUGH THIS PROCESS.  WE CAN
!             COMPUTE TPAR USING THIS KNOWLEDGE USING LOOK
!             UP TABLE (SUBROUTINE TTBLEX).
!        (4)  FIND THE EQUILIBRIUM LEVEL.  THIS IS DEFINED AS THE
!             HIGHEST POSITIVELY BUOYANT LAYER.
!             (IF THERE IS NO POSITIVELY BUOYANT LAYER, CAPE/CINS
!              WILL BE ZERO)
!        (5)  COMPUTE CAPE/CINS.  
!             (A) COMPUTE THETAP.  WE KNOW TPAR AND P.
!             (B) COMPUTE THETAA.  WE KNOW T AND P.  
!        (6)  ADD G*(THETAP-THETAA)*DZ TO THE RUNNING CAPE OR CINS SUM.
!             (A) IF THETAP > THETAA, ADD TO THE CAPE SUM.
!             (B) IF THETAP < THETAA, ADD TO THE CINS SUM.
!        (7)  ARE WE AT EQUILIBRIUM LEVEL? 
!             (A) IF YES, STOP THE SUMMATION.
!             (B) IF NO, CONTIUNUE THE SUMMATION.
!        (8)  ENFORCE LIMITS ON CAPE AND CINS (I.E. NO NEGATIVE CAPE)
!     
! PROGRAM HISTORY LOG:
!   93-02-10  RUSS TREADON
!   93-06-19  RUSS TREADON - GENERALIZED ROUTINE TO ALLOW FOR 
!                            TYPE 2 CAPE/CINS CALCULATIONS.     
!   94-09-23  MIKE BALDWIN - MODIFIED TO USE LOOK UP TABLES
!                            INSTEAD OF COMPLICATED EQUATIONS.
!   94-10-13  MIKE BALDWIN - MODIFIED TO CONTINUE CAPE/CINS CALC
!                            UP TO AT HIGHEST BUOYANT LAYER.
!   98-06-12  T BLACK      - CONVERSION FROM 1-D TO 2-D
!   98-08-18  T BLACK      - COMPUTE APE INTERNALLY
!   00-01-04  JIM TUCCILLO - MPI VERSION              
!   02-01-15  MIKE BALDWIN - WRF VERSION
!   03-08-24  G MANIKIN    - ADDED LEVEL OF PARCEL BEING LIFTED
!                            AS OUTPUT FROM THE ROUTINE AND ADDED
!                            THE DEPTH OVER WHICH ONE SEARCHES FOR
!                            THE MOST UNSTABLE PARCEL AS INPUT
!   10-09-09  G MANIKIN    - CHANGED COMPUTATION TO USE VIRTUAL TEMP
!                          - ADDED EQ LVL HGHT AND THUNDER PARAMETER    
!   15-xx-xx  S MOORTHI    - optimization and threading
!   19-09-03  J MENG       - MODIFIED TO ADD 0-3KM CAPE/CINS, LFC, 
!                            EFFECTIVE HELICITY, DOWNDRAFT CAPE,
!                            DENDRITIC GROWTH LAYER DEPTH, ESP
!   21-09-01  E COLON      - equivalent level height index for RTMA
!
! USAGE:    CALL CALCAPE2(ITYPE,DPBND,P1D,T1D,Q1D,L1D,    &
!                          CAPE,CINS,LFC,ESRHL,ESRHH,     &
!                          DCAPE,DGLD,ESP)
!
!   INPUT ARGUMENT LIST:
!     ITYPE    - INTEGER FLAG SPECIFYING HOW PARCEL TO LIFT IS
!                IDENTIFIED.  SEE COMMENTS ABOVE.
!     DPBND    - DEPTH OVER WHICH ONE SEARCHES FOR MOST UNSTABLE PARCEL
!     P1D      - ARRAY OF PRESSURE OF PARCELS TO LIFT.
!     T1D      - ARRAY OF TEMPERATURE OF PARCELS TO LIFT.
!     Q1D      - ARRAY OF SPECIFIC HUMIDITY OF PARCELS TO LIFT.
!     L1D      - ARRAY OF MODEL LEVEL OF PARCELS TO LIFT.
!
!   OUTPUT ARGUMENT LIST: 
!     CAPE     - CONVECTIVE AVAILABLE POTENTIAL ENERGY (J/KG)
!     CINS     - CONVECTIVE INHIBITION (J/KG)
!     LFC      - LEVEL OF FREE CONVECTION (M)
!     ESRHL    - LOWER BOUND TO ACCOUNT FOR EFFECTIVE HELICITY CALCULATION
!     ESRHH    - UPPER BOUND TO ACCOUNT FOR EFFECTIVE HELICITY CALCULATION
!     DCAPE    - DOWNDRAFT CAPE (J/KG)
!     DGLD     - DENDRITIC GROWTH LAYER DEPTH (M)
!     ESP      - ENHANCED STRETCHING POTENTIAL
!     
!   OUTPUT FILES:
!     STDOUT  - RUN TIME STANDARD OUT.
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       BOUND    - BOUND (CLIP) DATA BETWEEN UPPER AND LOWER LIMTS.
!       TTBLEX   - LOOKUP TABLE ROUTINE TO GET T FROM THETAE AND P
!
!     LIBRARY:
!       COMMON   - 
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN 90
!     MACHINE : CRAY C-90
!$$$  
!
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
!
!-------------------------------------------------------------------------------------
!

!> @file
!
!> SUBPROGRAM:    CALVOR      COMPUTES ABSOLUTE VORTICITY
!!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-22       
!!     
!! ABSTRACT:  
!!     THIS ROUTINE COMPUTES THE ABSOLUTE VORTICITY.
!!     
!! PROGRAM HISTORY LOG:
!!   92-12-22  RUSS TREADON
!!   98-06-08  T BLACK - CONVERSION FROM 1-D TO 2-D
!!   00-01-04  JIM TUCCILLO - MPI VERSION
!!   02-01-15  MIKE BALDWIN - WRF VERSION C-GRID
!!   05-03-01  H CHUANG - ADD NMM E GRID
!!   05-05-17  H CHUANG - ADD POTENTIAL VORTICITY CALCULATION
!!   05-07-07  B ZHOU   - ADD RSM IN COMPUTING DVDX, DUDY AND UAVG
!!   13-08-09  S MOORTHI - Optimize the vorticity loop including threading
!!   16-08-05  S Moorthi - add zonal filetering
!!   2019-10-17 Y Mao - Skip calculation when U/V is SPVAL
!!   2020-11-06 J Meng - USE UPP_MATH MODULE
!!   21-09-02  Bo Cui  - Decompose UPP in X direction, REPLACE EXCH_F to EXCH
!!   21-10-31  J MENG  - 2D DECOMPOSITION
!!     
!! USAGE:    CALL CALVOR(UWND,VWND,ABSV)
!!   INPUT ARGUMENT LIST:
!!     UWND     - U WIND (M/S) MASS-POINTS
!!     VWND     - V WIND (M/S) MASS-POINTS
!!
!!   OUTPUT ARGUMENT LIST: 
!!     ABSV     - ABSOLUTE VORTICITY (1/S) MASS-POINTS
!!     
!!   OUTPUT FILES:
!!     NONE
!!     
!!   SUBPROGRAMS CALLED:
!!     UTILITIES:
!!       NONE
!!     LIBRARY:
!!       COMMON   - CTLBLK
!!     
!!   ATTRIBUTES:
!!     LANGUAGE: FORTRAN
!!     MACHINE : WCOSS
!!
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
            IF(VWND(I+1,J)<SPVAL.AND.VWND(I-1,J)<SPVAL.AND.              &
     &         UWND(I,J+1)<SPVAL.AND.UWND(I,J-1)<SPVAL) THEN
              DVDX   = DDVDX(I,J)
              DUDY   = DDUDY(I,J)
              UAVG   = UUAVG(I,J) 
!  is there a (f+tan(phi)/erad)*u term?
              IF(MODELNAME  == 'RAPR') then
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

      SUBROUTINE CALDIV(UWND,VWND,DIV)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALDIV      COMPUTES DIVERGENCE
!   PRGRMMR: SAJAL KAR         ORG: W/NP2      DATE: 16-05-05
!     
! ABSTRACT:  
!     FOR GFS, THIS ROUTINE COMPUTES THE HORIZONTAL DIVERGENCE
!     USING 2ND-ORDER CENTERED SCHEME ON A LAT-LON GRID     
!
! PROGRAM HISTORY LOG:
!   16-05-05  SAJAL KAR MODIFIED CALVORT TO COMPUTE DIVERGENCE FROM
!             WIND COMPONENTS
!   16-07-22  S Moorthi modifying polar divergence calculation
!     
! USAGE:    CALL CALDIV(UWND,VWND,DIV)
!   INPUT ARGUMENT LIST:
!     UWND     - U WIND (M/S) MASS-POINTS
!     VWND     - V WIND (M/S) MASS-POINTS
!
!   OUTPUT ARGUMENT LIST: 
!     DIV     - DIVERGENCE (1/S) MASS-POINTS
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       NONE
!     LIBRARY:
!       COMMON   - CTLBLK
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : WCOSS
!$$$  
!     
!
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
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM: CALGRADPS COMPUTES GRADIENTS OF A SCALAR FIELD PS OR LNPS
!   PRGRMMR: SAJAL KAR         ORG: W/NP2      DATE: 16-05-05
!     
! ABSTRACT:  
!     FOR GFS, THIS ROUTINE COMPUTES  HRIZONTAL GRADIENTS OF PS OR LNPS
!     USING 2ND-ORDER CENTERED SCHEME ON A LAT-LON GRID
!     
! PROGRAM HISTORY LOG:
!   16-05-05  SAJAL KAR REDUCED FROM CALVORT TO ZONAL AND MERIDIONAL
!             GRADIENTS OF GIVEN SURFACE PRESSURE PS, OR LNPS
!     
! USAGE:    CALL CALGRADPS(PS,PSX,PSY)
!   INPUT ARGUMENT LIST:
!     PS       - SURFACE PRESSURE (PA) MASS-POINTS
!
!   OUTPUT ARGUMENT LIST: 
!     PSX     - ZONAL GRADIENT OF PS AT MASS-POINTS
!     PSY     - MERIDIONAL GRADIENT OF PS AT MASS-POINTS
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       NONE
!     LIBRARY:
!       COMMON   - CTLBLK
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : WCOSS
!$$$  
!     
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
!
!-------------------------------------------------------------------------------------
!
  end module upp_physics

