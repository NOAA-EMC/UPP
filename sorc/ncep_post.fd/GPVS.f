!> @file
!> @brief gpvs() computes saturation vapor pressure table.
!>
!> Compute saturation vapor pressure table as a function of
!> temperature for the table lookup function FPVS.
!> Exact saturation vapor pressures are calculated in subprogram FPVSX.
!> The current implementation computes a table with a length
!> of 7501 for temperatures ranging from 180.0 to 330.0 Kelvin.
!>
!> @param[out] pvu real (km) potential vorticity (10**-6*K*m**2/kg/s).
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1982-12-30 | N Phillips   | Initial
!> 1991-05-07 | Mark Iredell | Made into inlinable function
!> 1994-12-30 | Mark Iredell | Expand table
!> 1996-02-19 | Hong         | Ice effect
!>
!> @note Lookup tables for the saturation vapor pressure w/r/t water & ice.
!> @author N Phillips W/NP2 @date 1982-12-30
!-----------------------------------------------------------------------
!>
!> gpvs computes saturation vapor pressure table as a function of
!> temperature for the table lookup function FPVS.
!>
      SUBROUTINE GPVS
!     ******************************************************************

!----------------------------------------------------------------------
      use svptbl_mod, only: nx, c1xpvs, c2xpvs, c1xpvs0, c2xpvs0, tbpvs, tbpvs0
!- - - - - - - - - - -- - - -- - - -- - - -- - - - - -- - - -- - - -
      implicit none
!
      real xmin,xmax,xinc,x,t
      integer jx
      real,external :: fpvsx,fpvsx0
!----------------------------------------------------------------------
      XMIN=180.0
      XMAX=330.0
      XINC=(XMAX-XMIN)/(NX-1)
      C1XPVS=1.-XMIN/XINC
      C2XPVS=1./XINC
      C1XPVS0=1.-XMIN/XINC
      C2XPVS0=1./XINC
!
      DO JX=1,NX
        X=XMIN+(JX-1)*XINC
        T=X
        TBPVS(JX)=FPVSX(T)
        TBPVS0(JX)=FPVSX0(T)
      ENDDO
! 
      RETURN
      END
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!> @brief fpvs() computes saturation vapor pressure.
!> @note This function is mostly replaced by FPVSNEW in UPP_PHYSICS.f. 
!> The only routine that uses FPVS is CALMICT.f.
!>
!> Compute saturation vapor pressure from the temperature.
!> A linear interpolation is done between values in a lookup table
!> computed in GPVS. See documentation for FPVSX for details.
!> Input values outside table range are reset to table extrema.
!> The interpolation accuracy is almost 6 decimal places.
!> On the CRAY, FPVS is about 4 times faster than exact calculation.
!> This function should be expanded inline in the calling routine.
!>
!> @param[in] T real temperature in Kelvin.
!> @return FPVS real Saturation vapor pressure in kilopascals (CB).
!> 
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1982-12-30 | N Phillips   | Initial
!> 1991-05-07 | Mark Iredell | Made into inlinable function
!> 1994-12-30 | Mark Iredell | Expand table
!> 1996-02-19 | Hong         | Ice effect
!>
!> @author N Phillips W/NP2 @date 1982-12-30
!-----------------------------------------------------------------------
                           FUNCTION FPVS(T)
!-----------------------------------------------------------------------
      use svptbl_mod, only : NX,C1XPVS,C2XPVS,TBPVS
!
      implicit none
!
!      integer,parameter::NX=7501
!      real C1XPVS,C2XPVS,TBPVS(NX)

      real T
      real XJ
      integer JX
      real FPVS
!-----------------------------------------------------------------------
      XJ=MIN(MAX(C1XPVS+C2XPVS*T,1.),FLOAT(NX))
      JX=MIN(XJ,NX-1.)
      FPVS=TBPVS(JX)+(XJ-JX)*(TBPVS(JX+1)-TBPVS(JX))
!
      RETURN
      END
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> @brief FPVS0() computes saturation vapor pressure.
!> @note This function is no longer used. 
!>
!> @param T real Temperature (K).
!> @param NX integer Number of grid cells in the X direction _____.
!> @param C1XPVS0 real _____.
!> @param C2XPVS0 real _____.
!> @param TBPVS0 array _____.
!> @return FPVS0 real Saturation vapor pressure in kilopascals (CB).
!>

                       FUNCTION FPVS0(T,NX,C1XPVS0,C2XPVS0,TBPVS0)
!-----------------------------------------------------------------------
!      use svptbl_mod, only : NX,C1XPVS0,C2XPVS0,TBPVS0
      implicit none
!
      integer NX
      real C1XPVS0,C2XPVS0,TBPVS0(NX)
     
      real T
      real XJ1
      integer JX1
      real FPVS0
!-----------------------------------------------------------------------
      XJ1=MIN(MAX(C1XPVS0+C2XPVS0*T,1.),FLOAT(NX))
      JX1=MIN(XJ1,NX-1.)
      FPVS0=TBPVS0(JX1)+(XJ1-JX1)*(TBPVS0(JX1+1)-TBPVS0(JX1))
!
      RETURN
      END
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!> @brief fpvsx() computes saturation vapor pressure.
!>
!> Exactly compute saturation vapor pressure from temperature.
!> The water model assumes a perfect gas, constant specific heats
!> for gas and liquid, and neglects the volume of the liquid.
!> The model does account for the variation of the latent heat
!> of condensation with temperature. The ice option is not included.
!> The Clausius-Clapeyron equation is integrated from the triple point
!> To get the formula
!> @code
!>     PVS=PSATK*(TR**XA)*exp(XB*(1.-TR))
!> @endcode
!> where TR is TTP/T and other values are physical constants
!> This function should be expanded inline in the calling routine.
!>
!> @param[in] T real Temperature (K).
!> @return FPVSX real Saturation vapor pressure in kilopascals (CB).
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1982-12-30 | N Phillips   | Initial
!> 1991-05-07 | Mark Iredell | Made into inlinable function
!> 1994-12-30 | Mark Iredell | Exact computation
!> 1996-02-19 | Hong         | Ice effect
!>
!> @author N Phillips W/NP2 @date 1982-12-30
!-----------------------------------------------------------------------
                         FUNCTION FPVSX(T)
!-----------------------------------------------------------------------
    implicit none
!
    real,PARAMETER :: CP=1.0046E+3,RD=287.04,RV=4.6150E+2,              &
            TTP=2.7316E+2,HVAP=2.5000E+6,PSAT=6.1078E+2,                &
            CLIQ=4.1855E+3,CVAP= 1.8460E+3,CICE=2.1060E+3,HSUB=2.8340E+6
    real,PARAMETER :: PSATK=PSAT*1.E-3
    real,PARAMETER :: DLDT=CVAP-CLIQ,XA=-DLDT/RV,XB=XA+HVAP/(RV*TTP)
    real,PARAMETER :: DLDTI=CVAP-CICE,XAI=-DLDTI/RV,XBI=XAI+HSUB/(RV*TTP)
    real :: TR, T
    real :: FPVSX
!-----------------------------------------------------------------------
    TR=TTP/T
!
    IF(T>=TTP)THEN
      FPVSX=PSATK*(TR**XA)*EXP(XB*(1.-TR))
    ELSE
      FPVSX=PSATK*(TR**XAI)*EXP(XBI*(1.-TR))
    ENDIF
!
    RETURN
    END
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!> @brief fpvsx0() computes saturation vapor pressure.
!> 
!> @param T real Temperature (K).
!> @return FPVSX0 real Saturation vapor pressure in kilopascals (CB).
!
                        FUNCTION FPVSX0(T)
!-----------------------------------------------------------------------
    implicit none
!
    real,PARAMETER :: CP=1.0046E+3,RD=287.04,RV=4.6150E+2,            &
              TTP=2.7316E+2,HVAP=2.5000E+6,PSAT=6.1078E+2,            &
              CLIQ=4.1855E+3,CVAP=1.8460E+3,CICE=2.1060E+3,           &
              HSUB=2.8340E+6
    real,PARAMETER :: PSATK=PSAT*1.E-3
    real,PARAMETER :: DLDT=CVAP-CLIQ,XA=-DLDT/RV,XB=XA+HVAP/(RV*TTP)
    real,PARAMETER :: DLDTI=CVAP-CICE,XAI=-DLDT/RV,XBI=XA+HSUB/(RV*TTP)
    real TR
    real T,FPVSX0
!-----------------------------------------------------------------------
    TR=TTP/T
    FPVSX0=PSATK*(TR**XA)*EXP(XB*(1.-TR))
!
    RETURN
    END
