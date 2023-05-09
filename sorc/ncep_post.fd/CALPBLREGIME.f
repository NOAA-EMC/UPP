!> @file
!> @brief Subroutine that computes PBL height based on bulk RCH number.
!>
!> This routine computes the bulk Richardson number based on algorithms
!> from WRF surface layer and then derives PBL regime as follows:
!> 1. BR >= 0.2;
!> Represents nighttime stable conditions (Regime=1),
!>
!> 2. BR < 0.2 .AND. BR > 0.0;
!> Represents damped mechanical turbulent conditions
!> (Regime=2),
!>
!> 3. BR == 0.0
!> Represents forced convection conditions (Regime=3),
!>
!> 4. BR < 0.0
!> Represnets free convection conditions (Regime=4).    
!>    
!> @param[out] PBLREGIME real PBL Height above ground.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2007-04-27 | H Chuang | Initial
!> 2021-09-02 | Bo Cui   | Decompose UPP in X direction          
!>   
!> @author H Chuang @date 2007-04-27
!-----------------------------------------------------------------------
!> @brief Subroutine that computes PBL height based on WRF algorithm for 
!> bulk Richardson number and then derives PBL regime.
!>
!> @param[inout] PBLREGIME real PBL Height above ground.
!-----------------------------------------------------------------------
      SUBROUTINE CALPBLREGIME(PBLREGIME)

!
      use vrbls3d,      only: uh, vh, pmid, t, q, pint, zmid, zint
      use vrbls2d,      only: ths, qs, smstav, twbs, qwbs, pblh
      use masks,        only: dx
      use params_mod,   only: p1000, capa, d608, h1, g, rd, cp
      use ctlblk_mod,   only: jsta, jend, spval, lm, jsta_m, jend_m, im,    &
                              jsta_2l, jend_2u, ista, iend, ista_m, iend_m,ista_2l,iend_2u
      use gridspec_mod, only: gridtype
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!     
!     INCLUDE,DERIVE,SET PARAMETERS.
!     
      REAL    , PARAMETER ::  VCONVC=1.
!     
!     DECLARE VARIABLES.
!     
      REAL,dimension(ista_2l:iend_2u,jsta_2l:jend_2u),intent(inout) ::  PBLREGIME
!
      integer I,J,IE,IW,ii,jj
      real APE,THV,THVX,GOVRTH,UMASS,VMASS,WSPD,TSKV,DTHV,RHOX,fluxc,tsfc,  &
           VCONV,VSGD,BR,THX
     
!
!     
!*************************************************************************
!     
!     INITIALIZE ARRAYS.
!
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            PBLREGIME(I,J) = SPVAL
          ENDDO
        ENDDO
!
!     COMPUTE BULK RICHARDSON NUMBER AS CODED IN WRF module_sf_sfclay
!
!!$omp  parallel do
!!$omp& private(uhkl,ulkl,vhkl,vlkl,rib,ubot,utop,vbot,vtop,
!!$omp&         betta,ricr,ustarr,wmin,tvhtop,ztop,
!!$omp&         wndsl,wndslp,betta,ricr,ustarr,wmin 
!!$omp&       ,IFRSTLEV
!!$omp&       ,ICALPBL
!!$omp&       ,LVLP
!!$omp&       ,RIF
!!$omp&       ,RIBP
!!$omp&       ,UBOT1
!!$omp&       ,VBOT1
!!$omp&       ,ZBOT1
!!$omp&       ,THVBOT1)
!
      IF(GRIDTYPE /= 'A')THEN
         call exch(UH(1,jsta_2l,LM))
         call exch(VH(1,jsta_2l,LM))
      END IF
             
      DO J=JSTA_M,JEND_M
        DO I=ISTA_M,IEND_M
!
          IF(PMID(I,J,LM)<SPVAL .AND. QS(I,J)<SPVAL .AND. &
             SMSTAV(I,J)<SPVAL) THEN 
          APE  = (P1000/PMID(I,J,LM))**CAPA
          THX  = T(I,J,LM)*APE
          THVX = (Q(I,J,LM)*D608+H1)*THX
          GOVRTH = G/THX
          IF(GRIDTYPE == 'E')THEN
            IE=I+MOD(J+1,2) 
            IW=I+MOD(J+1,2)-1
            UMASS = (UH(I,J-1,LM)+UH(IW,J,LM)+UH(IE,J,LM)              &   
                  +  UH(I,J+1,LM))/4.0
            VMASS = (VH(I,J-1,LM)+VH(IW,J,LM)+VH(IE,J,LM)              &
                  +  VH(I,J+1,LM))/4.0
            WSPD= SQRT(UMASS*UMASS+VMASS*VMASS)
          ELSE IF(GRIDTYPE == 'B')THEN
            IE = I
            IW = I-1
            UMASS = (UH(IW,J-1,LM)+UH(IW,J,LM)+UH(IE,J-1,LM)              &   
                  +  UH(I,J,LM))/4.0
            VMASS = (VH(IW,J-1,LM)+VH(IW,J,LM)+VH(IE,J-1,LM)              &
                  +  VH(I,J,LM))/4.0
            WSPD= SQRT(UMASS*UMASS+VMASS*VMASS)  
          ELSE
            WSPD = SQRT(UH(I,J,LM)*UH(I,J,LM)+VH(I,J,LM)*VH(I,J,LM))
          END IF
                                                                                 
          TSKV = THS(I,J)*(1.+D608*QS(I,J)*SMSTAV(I,J))
          DTHV = (THVX-TSKV)
!  Convective velocity scale Vc and subgrid-scale velocity Vsg
!  following Beljaars (1995, QJRMS) and Mahrt and Sun (1995, MWR)
!                                ... HONG Aug. 2001
!
          rhox  = PINT(I,J,LM+1)/RD/(T(I,J,LM)*(Q(I,J,LM)*D608+H1)) !density
          fluxc = max(-twbs(i,j)/rhox/cp - d608*tskv*QWBS(i,j)/rhox,0.)
          tsfc  = THS(I,J)*(PINT(I,J,LM+1)/P1000)**CAPA
          VCONV = vconvc*(g/tsfc*pblh(i,j)*fluxc)**.33
! VCONV comes from Beljaars only
          VSGD  = 0.32 * (max(dx(i,j)/5000.-1.,0.))**.33
          WSPD  = SQRT(WSPD*WSPD+VCONV*VCONV+vsgd*vsgd)
          WSPD  = MAX(WSPD,0.1)
          BR    = GOVRTH*(ZMID(I,J,LM)-ZINT(I,J,LM+1))*DTHV/(WSPD*WSPD)
     
          IF(BR < 0.0) THEN
            PBLREGIME(I,J) = 4.0
          ELSE IF(BR == 0.0) THEN
            PBLREGIME(I,J) = 3.0
          ELSE IF(BR < 0.2) THEN
            PBLREGIME(I,J) = 2.0
          ELSE
            PBLREGIME(I,J) = 1.0 
          END IF

!         ii=im/2
!         jj=(jsta+jend)/2
!         if(i==ii.and.j==jj)print*,'Debug: CALPBLREGIME ',i,j,br,     &  
!         PBLREGIME(I,J)
          END IF !end IF PMID 
   
        ENDDO
      ENDDO
!      
!     END OF ROUTINE.
!     
      RETURN
      END

