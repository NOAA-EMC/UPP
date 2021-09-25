!> @file
!> @brief Subroutines related to aviation.
!
!> Computes Low Level Wind Shear (0-2000feet) 
!!     
!! This program computes the low level wind shear(LLWS) over 0-2000 feet
!! (0-609.5m) layer. But because 10m wind represent sfc wind, 10-619.5 m
!! layer is used. (NOAA/NWS Instruction 10-813, 2004)
!!
!! Definition: LLWS(Z1,Z2) is vector difference of wind at z1 and z2
!! where:
!! - Z1 = 10m   + Surface height                
!! - Z2 = 619.5 + Surface height
!!
!! Algorithm: since Z2 is not defined in the model, so, first thing is
!! searching Z2 to see which layers it is located(ie between which two
!! pressure levels), then find the wind vector (U2,V2)at Z2 by
!! interpolating with the wind vectors of the at pressure levels above
!! and below then compute the vector difference between Z2 and Z1 (ie
!! U10,V10)
!!
!!<pre>                               
!!      ----------------------------------------- K2-1 ---------------------
!!                            ^
!!                            |
!!                            |
!!                            |            
!!                 ____       |  _____ Z2, U2=interpo[U(K2),U(K2-1)]
!!                  ^         |            V2=interpo[V(K2),V(K2-1)]
!!                  |         |                       
!!      ------------|---------|------------------ K2 ------------------------
!!                  |         |
!!                  |         |DH=SUM of all layers between K1-1 & K2-1 
!!                  |         |                                            .              
!!                  |609.5m   |                                            .
!!                  |(2000ft) |                                            .
!!                  |         v
!!      ------------|---------------------------------------------------- LSM-2
!!                  |               ^
!!                  |               |ZH1   
!!                  |               |
!!                 o-o 10m       ___v__ Z1,U10,V10                 
!!       FIS    ....|.....          ^
!!        ^   .            .        |
!!      --|-------------------------|------------ K1 -------------------- LSM-1
!!        | .                .      |
!!        |.                  .     |
!!       .|                    ...  |
!!      --|-------------------------|------------------------------------- LSM
!!      . |                         |
!!     ////////////////////////////////////////////////////////////////// Sea Level
!!</pre>                               
!!
!!
!! @param[in] U U wind profile (m/s) (at pressure level).
!! @param[in] V V wind (m/s) (at pressure level).
!! @param[in] H Height (m) (at pressure level).
!! @param[out] LLWS Low level wind shear (Knots/2000ft).
!!
!! Program History      
!! - 19-10-30  Bo CUI - REMOVE "GOTO" STATEMENT
!! - 21-04-01  Jesse Meng - computation on defined points only
!!     
!! @author Binbin Zhou NCEP/EMC  @date 2005-08-16       
      SUBROUTINE CALLLWS(U,V,H,LLWS)

!
      USE vrbls2d, only: fis, u10, v10
      use params_mod, only: gi
      use ctlblk_mod, only: jsta, jend, im, jm, lsm, spval
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL,DIMENSION(IM,JM,LSM),INTENT(IN)    :: U,V,H
      REAL,DIMENSION(IM,JM),INTENT(INOUT)     :: LLWS
      REAL    :: Z1,Z2,HZ1,DH,U2,V2,W2,RT 
      INTEGER :: K1,K2 
      integer I,J,LP

!***************************************************************
!
!

      DO 100 J=JSTA,JEND
        DO I=1,IM
 
          Z1 = 10.0 + FIS(I,J)*GI                              !Height of 10m levels geographic height (from sea level)
          
          IF(Z1<H(I,J,LSM)) THEN                            !First search location of 10m wind level
            K1 = LSM + 1                                       !to see it is in which pressure levels
          ELSE
            DO LP = LSM,2,-1                                   !If not found, keep searching upward                              
             IF(Z1>=H(I,J,LP).AND.Z1<H(I,J,LP-1)) THEN
               K1 = LP 
             END IF
            END DO
          END IF

          HZ1 = H(I,J,K1-1) - Z1                                !Distance between K1-1 and 10m level
 
          DH = 0.0

          IF((HZ1+10)>609.6) THEN                            !Then, search 2000ft(609.6m) location
            U2= U10(I,J) + (U(I,J,K1-1)-U10(I,J))*599.6/HZ1     !found it between K1-1 and K1, then linear
            V2= V10(I,J) + (V(I,J,K1-1)-V10(I,J))*599.6/HZ1     !interpolate to get wind at 2000ft U2,V2     
            Z2= FIS(I,J)*GI + 609.6
          ELSE                                                 !otherwise, keep on search upward
            DO LP = K1-1,2,-1
             DH=DH+(H(I,J,LP-1) - H(I,J,LP))
             IF((DH+HZ1+10)>609.6) THEN                      !found the 2000ft level 
               Z2=FIS(I,J)*GI+609.6   
               RT=(Z2-H(I,J,LP))/(H(I,J,LP-1)-H(I,J,LP))
               U2=U(I,J,LP)+RT*(U(I,J,LP-1)-U(I,J,LP))
               V2=V(I,J,LP)+RT*(V(I,J,LP-1)-V(I,J,LP))
               K2=LP
               exit
              END IF
             END DO
            END IF

!computer vector difference
         LLWS(I,J) = spval
         if(U10(I,J)<spval.and.V10(I,J)<spval)                   &
          LLWS(I,J)=SQRT((U2-U10(I,J))**2+(V2-V10(I,J))**2)/     &
                    609.6 * 1.943*609.6                         !unit: knot/2000ft
        ENDDO
 
100   CONTINUE     

      RETURN
      END

!> Computes In-Flight Icing.
!>     
!> This program computes the in-flight icing condition
!> with the T-RH-OMGA algorithm provided by S. Silberberg of
!> NCEP/AWC (improved new version).
!> 
!> According to S. Silberberg, Icing happens in following 
!> situation:
!> 1. -22C < T < 0C to      
!> 2.  RH > 70 %
!> 3. Ascent air, OMGA < 0 
!> 4. Equivalent Potential Vorticity (EPV) < 0
!> 5. Cloud water if SLD (supercooled large droplet)
!>
!> Current version dosn't consider SLD, so cloud water           
!> is not used. EPV computation is not available for current
!> NCEP/EMC models(NAM, WRF, RSM), so EPV is also not
!> used.
!>
!> @param[in] T1 TEMPERATURE (K)
!> @param[in] RH RELATIVE HUMIDITY  (DECIMAL FORM)
!> @param[in] OMGA Vertical velocity (Pa/sec)
!> @param[inout] ICING ICING CONDITION (1 or 0)
!>     
!> @author Binbin Zhou NCEP/EMC  @date 2005-08-16       
      SUBROUTINE CALICING (T1,RH,OMGA, ICING)
      use ctlblk_mod, only: jsta, jend, im, spval
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     DECLARE VARIABLES.
!     
      REAL, DIMENSION(IM,jsta:jend), INTENT(IN)    :: T1,RH,OMGA
      REAL, DIMENSION(IM,jsta:jend), INTENT(INOUT) :: ICING 
      integer I,J
!***************************************************************
!
!
      DO J=JSTA,JEND
        DO I=1,IM
        IF(OMGA(I,J)<SPVAL.AND.T1(I,J)<SPVAL.AND.RH(I,J)<SPVAL) THEN
         IF(OMGA(I,J) < 0.0 .AND.                       &
            (T1(I,J) <= 273.0 .AND. T1(I,J) >= 251.0)   &
              .AND. RH(I,J) >= 70.0) THEN

           ICING(I,J) = 1.0
         ELSE
           ICING(I,J) = 0.0
         END IF
        ELSE
           ICING(I,J) = SPVAL
        ENDIF
        ENDDO
      ENDDO

      RETURN
      END

!> Computes Clear Air Turbulence Index 
!>     
!> This program computes the Clear Air Turbulence condition which is
!> expressed as Index with Ellrod Algorithm (Gary P. Ellrod: Wea. and
!> Forecast,1992) and Ri number suggested by S. Silberberg of AWC. But Ri
!> number is still not classified into 3 level CAT, so current version
!> does not use Ri as suggested by S. Silberberg.
!>
!> PROGRAM HISTORY LOG:
!> - 05-09-19  H CHUANG - MODIFIED TO COMPUTE GRADIENTS FOR BOTH A AND E GRIDS
!>
!> According to Ellrod, the CAT is classied into 3 levels (index):
!> - Light:  CAT = 1     
!> - Middle: CAT = 2
!> - Severe: CAT = 3
!> - No CAT: CAT = 0 
!>
!> @param[in] U U wind profile (m/s) (at pressure level)
!> @param[in] V V wind (m/s) (at pressure level)
!> @param[in] H Height (m) (at pressure level)
!> @param[in] L # of pressure level
!> @param[inout] CAT CAT Index
!>     
!> @author Binbin Zhou NCEP/EMC @date 2005-08-16       
      SUBROUTINE CALCAT(U,V,H,U_OLD,V_OLD,H_OLD,CAT)
      use masks, only: dx, dy
      use ctlblk_mod, only: spval, jsta_2l, jend_2u, jsta_m, jend_m, &
              im, jm
      use gridspec_mod, only: gridtype
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none

!     
!     DECLARE VARIABLES.
!     
      REAL,DIMENSION(IM,jsta_2l:jend_2u),INTENT(IN)    :: U,V,H, &
                                                U_OLD,V_OLD,H_OLD
!      INTEGER,INTENT(IN)                      :: L
      REAL,DIMENSION(IM,jsta_2l:jend_2u),INTENT(INOUT) :: CAT

      REAL  DSH, DST, DEF, CVG, VWS, TRBINDX
      INTEGER  IHE(JM),IHW(JM)
      integer I,J
      integer ISTART,ISTOP,JSTART,JSTOP
      real VWS1,VWS2,VWS3,VWS4

!***************************************************************
!
!
      CAT=SPVAL 
      DO J=JSTA_2L,JEND_2U
       IF(GRIDTYPE == 'A')THEN
        IHW(J)=-1
        IHE(J)=1 
	ISTART=2
        ISTOP=IM-1
        JSTART=JSTA_M
        JSTOP=JEND_M
       ELSE IF(GRIDTYPE=='E')THEN
        IHW(J)=-MOD(J,2)
        IHE(J)=IHW(J)+1
	ISTART=2
        ISTOP=IM-1
        JSTART=JSTA_M
        JSTOP=JEND_M
       ELSE IF(GRIDTYPE=='B')THEN
        IHW(J)=-1
        IHE(J)=0 
	ISTART=2
        ISTOP=IM-1
        JSTART=JSTA_M
        JSTOP=JEND_M
       ELSE	
        print*,'no gridtype specified, exit calcat comp'
	return	
       END IF	
      ENDDO

      call exch_f(U)
      call exch_f(V)
      call exch_f(U_OLD)
      call exch_f(V_OLD)
      call exch_f(H)
      call exch_f(H_OLD)

      DO 100 J=JSTART,JSTOP
        DO I=ISTART,ISTOP
!
          IF(GRIDTYPE=='B')THEN
           IF(U(I,J)<spval.and.U(I,J-1)<spval.and.U(I-1,J)<spval.and.U(I-1,J-1)<spval.and.&
              V(I,J)<spval.and.V(I,J-1)<spval.and.V(I-1,J)<spval.and.V(I-1,J-1)<spval)THEN
!dsh=dv/dx+du/dy 
           DSH=(0.5*(V(I,J)+V(I,J-1))-0.5*(V(I-1,J)+V(I-1,J-1)))*10000./DX(I,J) &
	      +(0.5*(U(I,J)+U(I-1,J))-0.5*(U(I,J-1)+U(I-1,J-1)))*10000./DY(I,J)
!dst=du/dx-dv/dy
           DST =(0.5*(U(I,J)+U(I,J-1))-0.5*(U(I-1,J)+U(I-1,J-1)))*10000./DX(I,J) &
	      -(0.5*(V(I,J)+V(I-1,J))-0.5*(V(I,J-1)+V(I-1,J-1)))*10000./DY(I,J)
           DEF = SQRT (DSH*DSH + DST*DST)

!cvg=-(du/dx+dv/dy)
           CVG = -((0.5*(U(I,J)+U(I,J-1))-0.5*(U(I-1,J)+U(I-1,J-1)))*10000./DX(I,J) &
                +(0.5*(V(I,J)+V(I-1,J))-0.5*(V(I,J-1)+V(I-1,J-1)))*10000./DY(I,J))	   
           ELSE
            DEF = SPVAL
            CVG = SPVAL
           ENDIF
          ELSE
           IF(U(I,J+1)<spval.and.U(I,J-1)<spval.and.U(I+IHE(J),J)<spval.and.U(I+IHW(J),J)<spval.and.&
              V(I,J+1)<spval.and.V(I,J-1)<spval.and.V(I+IHE(J),J)<spval.and.V(I+IHW(J),J)<spval)THEN
!dsh=dv/dx+du/dy           
           DSH = (V(I+IHE(J),J) - V(I+IHW(J),J))*10000./(2*DX(I,J))   &  
              + (U(I,J+1) - U(I,J-1))*10000./(2*DY(I,J))

!dst=du/dx-dv/dy
           DST = (U(I+IHE(J),J) - U(I+IHW(J),J))*10000./(2*DX(I,J))   & 
              - (V(I,J+1) - V(I,J-1))*10000./(2*DY(I,J))

           DEF = SQRT (DSH*DSH + DST*DST)

!cvg=-(du/dx+dv/dy)
           CVG = -( (U(I+IHE(J),J) - U(I+IHW(J),J))*10000./(2*DX(I,J)) &
                  +(V(I,J+1) - V(I,J-1))*10000./(2*DY(I,J)) )
           ELSE
            DEF = SPVAL
            CVG = SPVAL
           ENDIF
          END IF
	  
          IF(GRIDTYPE == 'A')THEN
!vws=d|U|/dz
           IF(U_OLD(I,J)<spval.and.U(I,J)<spval.and.&
              V_OLD(I,J)<spval.and.V(I,J)<spval.and.&
              H_OLD(I,J)<spval.and.H(I,J)<spval)THEN
           VWS = ( SQRT(U_OLD(I,J)**2+V_OLD(I,J)**2 ) -               &
                  SQRT(U(I,J)**2+V(I,J)**2 )   ) *                    &
                  1000.0/(H_OLD(I,J) - H(I,J))
           ELSE
            VWS = SPVAL
           ENDIF
          else IF(GRIDTYPE == 'E')THEN
!vws=d|U|/dz
           IF(U_OLD(I+IHE(J),J)<spval.and.U(I+IHE(J),J)<spval.and.&
              V_OLD(I+IHE(J),J)<spval.and.V(I+IHE(J),J)<spval)THEN
            
	   VWS1 = ( SQRT(U_OLD(I+IHE(J),J)**2+V_OLD(I+IHE(J),J)**2 ) -&
                  SQRT(U(I+IHE(J),J)**2+V(I+IHE(J),J)**2 )   ) 
           ELSE
            VWS1 = SPVAL
           ENDIF
!vws=d|U|/dz
           IF(U_OLD(I+IHW(J),J)<spval.and.U(I+IHW(J),J)<spval.and.&
              V_OLD(I+IHW(J),J)<spval.and.V(I+IHW(J),J)<spval)THEN
           VWS2 = ( SQRT(U_OLD(I+IHW(J),J)**2+V_OLD(I+IHW(J),J)**2 ) -&   
                  SQRT(U(I+IHW(J),J)**2+V(I+IHW(J),J)**2 )   ) 
           ELSE
            VWS2 = SPVAL
           ENDIF
!vws=d|U|/dz
           IF(U_OLD(I,J-1)<spval.and.U(I,J-1)<spval.and.&
              V_OLD(I,J-1)<spval.and.V(I,J-1)<spval)THEN
           VWS3 = ( SQRT(U_OLD(I,J-1)**2+V_OLD(I,J-1)**2 ) -          & 
                  SQRT(U(I,J-1)**2+V(I,J-1)**2 )   ) 
           ELSE
            VWS3 = SPVAL
           ENDIF
!vws=d|U|/dz
           IF(U_OLD(I,J+1)<spval.and.U(I,J+1)<spval.and.&
              V_OLD(I,J+1)<spval.and.V(I,J+1)<spval)THEN
           VWS4 = ( SQRT(U_OLD(I,J+1)**2+V_OLD(I,J+1)**2 ) -          & 
                  SQRT(U(I,J+1)**2+V(I,J+1)**2 )   ) 
           ELSE
            VWS4 = SPVAL
           ENDIF

           IF(VWS1<spval.and.VWS2<spval.and.VWS3<spval.and.VWS4<spval.and.&
              H_OLD(I,J)<spval.and.H(I,J)<spval)THEN
           VWS=1000.0*(VWS1+VWS2+VWS3+VWS4)/4.0/(H_OLD(I,J) - H(I,J))
           ELSE
            VWS = SPVAL
           ENDIF
	  ELSE IF(GRIDTYPE == 'B')THEN
           IF(U_OLD(I+IHE(J),J)<spval.and.U(I+IHE(J),J)<spval.and.&
              V_OLD(I+IHE(J),J)<spval.and.V(I+IHE(J),J)<spval)THEN
	   VWS1 = ( SQRT(U_OLD(I+IHE(J),J)**2+V_OLD(I+IHE(J),J)**2 ) -&
                  SQRT(U(I+IHE(J),J)**2+V(I+IHE(J),J)**2 )   ) 
           ELSE
            VWS1 = SPVAL
           ENDIF
!vws=d|U|/dz
           IF(U_OLD(I+IHW(J),J)<spval.and.U(I+IHW(J),J)<spval.and.&
              V_OLD(I+IHW(J),J)<spval.and.V(I+IHW(J),J)<spval)THEN
           VWS2 = ( SQRT(U_OLD(I+IHW(J),J)**2+V_OLD(I+IHW(J),J)**2 ) -&   
                  SQRT(U(I+IHW(J),J)**2+V(I+IHW(J),J)**2 )   ) 
           ELSE
            VWS2 = SPVAL
           ENDIF
!vws=d|U|/dz
           IF(U_OLD(I,J-1)<spval.and.U(I,J-1)<spval.and.&
              V_OLD(I,J-1)<spval.and.V(I,J-1)<spval)THEN
           VWS3 = ( SQRT(U_OLD(I,J-1)**2+V_OLD(I,J-1)**2 ) -          & 
                  SQRT(U(I,J-1)**2+V(I,J-1)**2 )   ) 
           ELSE
            VWS3 = SPVAL
           ENDIF
!vws=d|U|/dz
           IF(U_OLD(I-1,J-1)<spval.and.U(I-1,J-1)<spval.and.&
              V_OLD(I-1,J-1)<spval.and.V(I-1,J-1)<spval)THEN
           VWS4 = ( SQRT(U_OLD(I-1,J-1)**2+V_OLD(I-1,J-1)**2 ) -          & 
                  SQRT(U(I-1,J-1)**2+V(I-1,J-1)**2 )   ) 
           ELSE
            VWS4 = SPVAL
           ENDIF

           IF(VWS1<spval.and.VWS2<spval.and.VWS3<spval.and.VWS4<spval.and.&
              H_OLD(I,J)<spval.and.H(I,J)<spval)THEN
           VWS=1000.0*(VWS1+VWS2+VWS3+VWS4)/4.0/(H_OLD(I,J) - H(I,J)) 
           ELSE
            VWS=SPVAL
           ENDIF
	  END IF  
          
         IF(VWS<spval.and.DEF<spval.and.CVG<spval)THEN 
          TRBINDX = ABS(VWS)*(DEF + ABS(CVG))
	  
          IF(TRBINDX<=4.) THEN
            CAT(I,J) = 0.0
          ELSE IF(TRBINDX<=8.) THEN
            CAT(I,J)=1.0
          ELSE IF(TRBINDX<=12.) THEN
            CAT(I,J)=2.0
          ELSE
            CAT(I,J)=3.0
          END IF        
         ELSE
          CAT(I,J)=SPVAL
         ENDIF
        ENDDO
 
100   CONTINUE     

      RETURN
      END

!> Computes ceiling.
!>     
!> This program computes the ceiling.  Definition: Ceiling is the cloud
!> base height for cloud fraction > 50% The cloud base is from sea level
!> in the model, while ceiling is from surface. If no ceiling, set
!> ceiling height = 20000 m
!>
!> @param[in] CLDZ CLOUD BASE HEIGHT from sea level(M)
!> @param[in] TCLD TOTAL CLOUD FRACTION (%)
!> @param[inout] CEILING CEILING HEIGHT from surface (m)
!>     
!> @author Binbin Zhou NCEP/EMC  @date 2005-08-18       
      SUBROUTINE CALCEILING (CLDZ,TCLD,CEILING)
      USE vrbls2d, only: fis
      use params_mod, only: small, gi
      use ctlblk_mod, only: jsta, jend, spval, im, modelname
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     DECLARE VARIABLES.
!     
      REAL, DIMENSION(IM,jsta:jend), INTENT(IN)    :: CLDZ, TCLD
      REAL, DIMENSION(IM,jsta:jend), INTENT(INOUT) :: CEILING
      integer I,J
!***************************************************************
!
!
      DO J=JSTA,JEND
        DO I=1,IM
          IF(ABS(TCLD(I,J)-SPVAL) <= SMALL) THEN
            CEILING(I,J)=SPVAL
          ELSE IF(TCLD(I,J) >= 50.) THEN
            if(MODELNAME == 'RAPR')then
              CEILING(I,J) = CLDZ(I,J) - FIS(I,J)*GI
            else
              CEILING(I,J) = CLDZ(I,J) ! for RAP/HRRR   - FIS(I,J)*GI
            endif
          ELSE
            CEILING(I,J) = 20000.0
          END IF

          IF(CEILING(I,J) < 0.0) CEILING(I,J)=20000.0

        ENDDO
      ENDDO

      RETURN
      END


!> Computes Ceiling.
!>     
!> This program computes the flight condition restriction 
!> which is defined as follow (NOAA/NWS/Instruction for TAF, 2004):
!>  
!> Ceiling(feet)| Visibility(miles) | FLTCND
!> -------------|-------------------|-------      
!> LIFR         | < 200 and/or < 1  |             1
!> IFR          | >= 500 to < 1000 and/or >=1 to <  3 |        2
!> MVFR         | >=1000 to <= 3000 and/or >=3 to <= 5|        3
!> VFR          | > 3000 > 5        |      5
!>
!>
!> @param[in] CEILING - CEILING HEIGHT from surface (m) NOTE: VIS -
!> Visibility is passed through COMMON /VISB/
!> @param[inout] FLTCND - FLIGHT CONDITION CATERGORY
!>      
!> @author Binbin Zhou NCEP/EMC @date 2005-08-18       
      SUBROUTINE CALFLTCND (CEILING,FLTCND)
      use vrbls2d, only: vis
      use ctlblk_mod, only: jsta, jend, im, spval
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     DECLARE VARIABLES.
!     
      REAL, DIMENSION(IM,jsta:jend), INTENT(IN)    :: CEILING
      REAL, DIMENSION(IM,jsta:jend), INTENT(INOUT) :: FLTCND
      REAL  CEIL,VISI
      integer I,J
!
!***************************************************************
!
!
      DO J=JSTA,JEND
        DO I=1,IM
 
        IF(CEILING(I,J)<spval.and.VIS(I,J)<spval)THEN
          CEIL = CEILING(I,J) * 3.2808               !from m -> feet
          VISI = VIS(I,J) / 1609.0                   !from m -> miles       

          IF(CEIL<500.0 .OR. VISI<1.0 ) THEN
             FLTCND(I,J) = 1.0

          ELSE IF( (CEIL>=500.AND.CEIL<1000.0) .OR.          &
                   (VISI>=1.0.AND.VISI<3.0) ) THEN
             FLTCND(I,J) = 2.0

          ELSE IF( (CEIL>=1000.AND.CEIL<=3000.0) .OR.         &
                   (VISI>=3.0.AND.VISI<=5.0) ) THEN
             FLTCND(I,J) = 3.0

          ELSE IF( CEIL>3000.0  .OR. VISI>5.0) THEN
             FLTCND(I,J) = 4.0

          END IF
        ELSE
          FLTCND(I,J) = SPVAL
        ENDIF
        ENDDO
      ENDDO

      RETURN
      END
