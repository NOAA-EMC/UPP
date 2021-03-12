!> @file
!
!> SUBPROGRAM:    MISCLN      POSTS MISCELLANEOUS FIELDS
!!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-20
!!     
!! ABSTRACT:
!!     THIS ROUTINE HAS BECOME THE CATCH-ALL FOR MISCELLANEOUS
!!     OUTPUT FIELDS POSTED BY THE ETA POST PROCESSOR.  
!!     CURRENTLY THIS ROUTINE POSTS THE FOLLOWING FIELDS:
!!        (1) TROPOPAUSE LEVEL Z,P, T, U, V, AND VERTICAL WIND SHEAR,
!!        (2) MAX WIND LEVEL Z, P, U, AND V,
!!        (3) FD LEVEL T, Q, U, AND V,
!!        (4) FREEZING LEVEL Z AND RH,
!!        (5) CONSTANT MASS (BOUNDARY) FIELDS,
!!        (6) LFM LOOK-ALIKE FIELDS, AND
!!        (7) NGM LOOK-ALIKE FIELDS.
!!
!!     
!! PROGRAM HISTORY LOG:
!!   92-12-20  RUSS TREADON
!!   93-06-19  RUSS TREADON - ADDED TYPE 2 CAPE POSTING.
!!   94-11-07  MIKE BALDWIN - ADDED HELICITY POSTING.
!!   96-03-26  MIKE BALDWIN - CHANGE ETA BOUNDARY LAYER LABELS FOR GRIB
!!   96-11-19  MIKE BALDWIN - BACK OUT PREVIOUS CHANGE 
!!   97-04-25  MIKE BALDWIN - CHANGE ETA BOUNDARY LAYER LABELS FOR GRIB
!!   97-04-29  GEOFF MANIKIN - ADDED TROPOPAUSE HEIGHT AND
!!                             MAX WIND LEVEL FIELDS
!!   98-06-15  T BLACK       - CONVERSION FROM 1-D TO 2-D
!!   98-07-17  MIKE BALDWIN - REMOVED LABL84
!!   00-01-04  JIM TUCCILLO - MPI VERSION
!!   02-04-23  MIKE BALDWIN - WRF VERSION
!!   11-02-06  JUN WANG     - ADD GRIB2 OPTION
!!   11-10-16  SARAH LU     - ADD FD LEVEL DUST/ASH
!!   12-04-03  Jun Wang     - FIXED LVLSXML for fields at FD height (spec_hgt_lvl_above_grnd)
!!   13-05-3   Shrinivas Moorthi - Fix some bugs and make more efficient code
!!   14-02-21  Shrinivas Moorthi - Add more threading
!!   14-02-26  S Moorthi - threading datapd assignment and some cleanup &
!!                         bug fix
!!   15-11-18  S Moorthi - fixed some logical errors in the helicity and
!!   i                     storm motion part of the code
!!   17-06-01  Y Mao - ADD FD levels for GTG(EDPARM CATEDR MWTURB) and allow 
!!                     levels input from control file
!!   19-09-03  J Meng - ADD CAPE related variables for HRRR
!!   20-03-24  J Meng - remove grib1
!!   20-11-10  J Meng - USE UPP_PHYSICS MODULE
!!     
!! USAGE:    CALL MISCLN
!!   INPUT ARGUMENT LIST:
!!
!!   OUTPUT ARGUMENT LIST: 
!!     NONE
!!     
!!   SUBPROGRAMS CALLED:
!!     UTILITIES:
!!       TRPAUS  - COMPUTE TROPOPAUSE LEVEL FIELDS.
!!       CALMXW  - COMPUTE MAX WIND LEVEL FIELDS.
!!       SCLFLD  - SCALE ARRAY ELEMENTS BY CONSTANT.
!!       GRIBIT  - OUTPUT FIELD TO GRIB FILE.
!!       CALPOT  - CALCULATE POTENTIAL TEMPERATURE.
!!       FDLVL   - COMPUTE FD LEVEL DATA (AGL OR MSL).
!!       FRZLVL  - COMPUTE FREEZING LEVEL DATA.
!!       BOUND   - BOUND ARRAY ELEMENTS BETWEEN MINIMUM AND MAXIMUM VALUES.
!!       BNDLYR  - COMPUTE BOUNDARY LAYER FIELDS.
!!       CALDWP  - CALCULATE DEWPOINT TEMPERATURE.
!!       OTLFT   - COMPUTE LIFTED INDEX AT 500MB.
!!       CALLCL  - COMPUTE LCL DATA.
!!       LFMFLD  - COMPUTE LFM LOOK-ALIKE FIELDS.
!!       NGMFLD  - COMPUTE NGM LOOK-ALIKE FIELDS.
!!       CALTHTE - COMPUTE THETA-E.
!!       CALHEL  - COMPUTE HELICITY AND STORM MOTION.
!!
!!     LIBRARY:
!!       COMMON - RQSTFLD
!!                CTLBLK
!!     
!!   ATTRIBUTES:
!!     LANGUAGE: FORTRAN
!!     MACHINE : CRAY C-90
!!
      SUBROUTINE MISCLN

!
!
      use vrbls3d,    only: pmid, uh, vh, t, zmid, zint, pint, alpint, q, omga
      use vrbls3d,    only: catedr,mwt,gtg
      use vrbls2d,    only: pblh, cprate, fis
      use masks,      only: lmh
      use params_mod, only: d00, d50, h99999, h100, h1, h1m12, pq0, a2, a3, a4,    &
                            rhmin, rgamog, tfrz, small, g
      use ctlblk_mod, only: grib, cfld, fld_info, datapd, im, jsta, jend, jm, jsta_m, jend_m, &
                            nbnd, nbin_du, lm, htfd, spval, pthresh, nfd, petabnd, me,&
                            jsta_2l, jend_2u, MODELNAME, SUBMODELNAME
      use rqstfld_mod, only: iget, lvls, id, iavblfld, lvlsxml
      use grib2_module, only: pset
      use upp_physics, only: FPVSNEW, CALRH_PW, CALCAPE, CALCAPE2
      use gridspec_mod, only: gridtype
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
!     SET LOCAL PARAMETERS.  MAKE SURE NFD AND NBND AGREE
!     WITH THE VALUES SET IN SUBROUTINES FDLVL AND BNDLYR,
!     RESPECTIVELY.
      real,PARAMETER :: C2K=273.15
      real,parameter :: con_rd      =2.8705e+2 ! gas constant air    (J/kg/K)
      real,parameter :: con_rv      =4.6150e+2 ! gas constant H2O
      real,parameter :: con_eps     =con_rd/con_rv
      real,parameter :: con_epsm1   =con_rd/con_rv-1
      real,parameter :: cpthresh    =0.000004
!     
!     DECLARE VARIABLES.
!     
      LOGICAL NORTH, FIELD1,FIELD2
      LOGICAL, dimension(IM,JSTA:JEND) :: DONE, DONE1

      INTEGER, allocatable ::  LVLBND(:,:,:),LB2(:,:)
!     INTEGER LVLBND(IM,JM,NBND),LB2(IM,JM),LPBL(IM,JM)

      real,dimension(im,jm)        :: GRID1, GRID2
      real,dimension(im,jsta:jend) :: P1D, T1D, Q1D, U1D, V1D, SHR1D, Z1D,   &
                                      RH1D, EGRID1, EGRID2, EGRID3, EGRID4,  &
                                      EGRID5, EGRID6, EGRID7, EGRID8
      real, dimension(:,:,:),allocatable :: OMGBND, PWTBND, QCNVBND,   &
                                            PBND,   TBND,   QBND,      &
                                            UBND,   VBND,   RHBND,     &
                                            WBND,   T7D,    Q7D,       &
                                            U7D,    V6D,    P7D,       &
                                            ICINGFD,GTGFD,CATFD,MWTFD
      real, dimension(:,:,:,:),allocatable :: AERFD

      real, dimension(:,:),allocatable ::   QM8510, RH4710, RH8498,    &
                                            RH4796, RH1847, UST, VST,  &
                                            RH3310, RH6610, RH3366,    &
                                            PW3310, RH4410, RH7294,    &
                                            RH4472,                    &
                                            T78483, T89671, P78483, P89671

      REAL, dimension(:,:,:),allocatable :: HELI
      real, dimension(:,:),  allocatable :: USHR1, VSHR1, USHR6, VSHR6, &
                                            MAXWP, MAXWZ, MAXWU, MAXWV, &
                                            MAXWT
!                                           MAXWT, RHPW
      INTEGER,dimension(:,:),allocatable :: LLOW, LUPP, HTSFC
      REAL, dimension(:,:),allocatable   :: CANGLE,ESHR,UVECT,VVECT
!
      integer I,J,jj,L,ITYPE,ISVALUE,LBND,ILVL,IFD,ITYPEFDLVL(NFD),    &
              iget1, iget2, iget3
      real    DPBND,PKL1,PKU1,FAC1,FAC2,PL,TL,QL,QSAT,RHL,TVRL,TVRBLO, &
              ES1,ES2,QS1,QS2,RH1,RH2,ZSF,DEPTH(2), work1, work2, work3

!     Variables introduced to allow FD levels from control file - Y Mao
      integer :: N,NFDCTL
      REAL, allocatable :: HTFDCTL(:)
      integer, allocatable :: ITYPEFDLVLCTL(:)
      integer IE,IW,JN,JS,IVE(JM),IVW(JM),JVN,JVS
      integer ISTART,ISTOP,JSTART,JSTOP,MIDCAL

!     
!****************************************************************************
!     START MISCLN HERE.
!     
         allocate(USHR1(IM,jsta_2l:jend_2u),VSHR1(IM,jsta_2l:jend_2u), &
                  USHR6(IM,jsta_2l:jend_2u),VSHR6(IM,jsta_2l:jend_2u))
         allocate(UST(IM,jsta_2l:jend_2u),VST(IM,jsta_2l:jend_2u),     &
                  HELI(IM,jsta_2l:jend_2u,2))
!
!      HELICITY AND STORM MOTION.
       iget1 = IGET(162)
       iget2 = -1
       iget3 = -1
       if (iget1 > 0) then
         iget2 = LVLS(1,iget1)
         iget3 = LVLS(2,iget1)
       endif
       IF (iget1 > 0 .OR. IGET(163) > 0 .OR. IGET(164) > 0) THEN
         DEPTH(1) = 3000.0
         DEPTH(2) = 1000.0
         CALL CALHEL(DEPTH,UST,VST,HELI,USHR1,VSHR1,USHR6,VSHR6)
         IF (iget2 > 0) then 
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = HELI(I,J,1)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(iget1)
             fld_info(cfld)%lvl=LVLSXML(1,iget1)
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF

         IF (iget3 > 0) then
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = HELI(I,J,2)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(iget1)
             fld_info(cfld)%lvl=LVLSXML(2,iget1)
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF

         IF (IGET(163) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = UST(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(163))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
         IF (IGET(164) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = VST(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(164))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
       ENDIF

!     UPDRAFT HELICITY

       if (IGET(427) > 0) THEN
         CALL CALUPDHEL(GRID1(1,jsta_2l))
         if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(427))
!$omp parallel do private(i,j,jj)
           do j=1,jend-jsta+1
             jj = jsta+j-1
             do i=1,im
               datapd(i,j,cfld) = GRID1(i,jj)
             enddo
           enddo
         endif

       ENDIF

! CRA  0-1 KM AND 0-6 KM SHEAR

       IF(IGET(430) > 0 .OR. IGET(431) > 0 .OR. IGET(432) > 0      &
                                           .OR. IGET(433) > 0) THEN

         DEPTH = 6000.0
         CALL CALHEL(DEPTH,UST,VST,HELI,USHR1,VSHR1,USHR6,VSHR6)

         IF(IGET(430) > 0) THEN
!$omp parallel do private(i,j,jj)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = USHR1(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(430))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
         IF(IGET(431) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = VSHR1(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(431))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
         IF(IGET(432) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = USHR6(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(432))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
         IF(IGET(433) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
                GRID1(I,J) = VSHR6(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(433))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF
       ENDIF
!
       if (allocated(ushr1)) deallocate(ushr1)
       if (allocated(vshr1)) deallocate(vshr1)
       if (allocated(ushr6)) deallocate(ushr6)
       if (allocated(vshr6)) deallocate(vshr6)
       if (allocated(ust))   deallocate(ust)
       if (allocated(vst))   deallocate(vst)
       if (allocated(heli))  deallocate(heli)
! CRA
!     
!
!     ***BLOCK 1:  TROPOPAUSE P, Z, T, U, V, AND WIND SHEAR.
!    
      IF ((IGET(054)>0).OR.(IGET(055)>0).OR.       &
          (IGET(056)>0).OR.(IGET(057)>0).OR.       &
          (IGET(177)>0).OR.                           &
          (IGET(058)>0).OR.(IGET(108)>0) ) THEN
! Chuang: Use GFS algorithm per Iredell's and DiMego's decision on unification
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM

              if(PMID(I,J,1)<spval) then
! INPUT
              CALL TPAUSE(LM,PMID(I,J,1:LM),UH(I,J,1:LM)             & 
! INPUT
                         ,VH(I,J,1:LM),T(I,J,1:LM),ZMID(I,J,1:LM)    &
! OUTPUT
                         ,P1D(I,J),U1D(I,J),V1D(I,J),T1D(I,J)        &
! OUTPUT
                         ,Z1D(I,J),SHR1D(I,J))                       ! OUTPUT
              else
                P1D(I,J) = spval
                U1D(I,J) = spval
                V1D(I,J) = spval
                T1D(I,J) = spval
                Z1D(I,J) = spval
                SHR1D(I,J) = spval
              endif

            END DO
          END DO
!
!        TROPOPAUSE PRESSURE.
         IF (IGET(054) > 0) THEN
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = P1D(I,J)
               ENDDO
             ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(054))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF

!        ICAO HEIGHT OF TROPOPAUSE
         IF (IGET(399)>0) THEN
           CALL ICAOHEIGHT(P1D, GRID1(1,jsta))
!            print*,'sample TROPOPAUSE ICAO HEIGHTS',GRID1(im/2,(jsta+jend)/2)
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(399))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
          ENDIF

!        TROPOPAUSE HEIGHT.
         IF (IGET(177) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = Z1D(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(177))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
!
!        TROPOPAUSE TEMPERATURE.
         IF (IGET(055) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = T1D(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(055))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
!
!        TROPOPAUSE POTENTIAL TEMPERATURE.
         IF (IGET(108) > 0) THEN
           CALL CALPOT(P1D,T1D,GRID1(1,jsta))
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(108))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
!     
!        TROPOPAUSE U WIND AND/OR V WIND.
         IF ((IGET(056) > 0).OR.(IGET(057) > 0)) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=U1D(I,J)
               GRID2(I,J)=V1D(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            if(IGET(056)>0) then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(056))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
            if(IGET(057)>0) then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(057))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID2(i,jj)
                enddo
              enddo
            endif
           endif
         ENDIF
!
!        TROPOPAUSE WIND SHEAR.
         IF (IGET(058) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = SHR1D(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(058))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
      ENDIF
!
!
!     ***BLOCK 2:  MAX WIND LEVEL  P, Z, U, AND V
!
!        MAX WIND LEVEL CALCULATIONS
      IF ((IGET(173)>0) .OR. (IGET(174)>0) .OR.                &
          (IGET(175)>0) .OR. (IGET(176)>0)) THEN

          allocate(MAXWP(IM,jsta:jend), MAXWZ(IM,jsta:jend),         &
                   MAXWU(IM,jsta:jend), MAXWV(IM,jsta:jend),MAXWT(IM,jsta:jend))
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
            DO I=1,IM
             MAXWP(I,J)=SPVAL
             MAXWZ(I,J)=SPVAL
             MAXWU(I,J)=SPVAL
             MAXWV(I,J)=SPVAL
            ENDDO
           ENDDO

!            CALL CALMXW(MAXWP,MAXWZ,MAXWU,MAXWV,MAXWT)
! Chuang: Use GFS algorithm per Iredell's and DiMego's decision on unification
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
           loopI:DO I=1,IM
            DO L=1,LM
              IF (ABS(PMID(I,J,L)-SPVAL)<=SMALL .OR. &
                  ABS(UH(I,J,L)-SPVAL)<=SMALL .OR. &
                  ABS(UH(I,J,L)-SPVAL)<=SMALL .OR. &
                  ABS(VH(I,J,L)-SPVAL)<=SMALL .OR. &
                  ABS(T(I,J,L)-SPVAL)<=SMALL .OR. &
                  ABS(ZMID(I,J,L)-SPVAL)<=SMALL) cycle loopI
            ENDDO
! INPUT
            CALL MXWIND(LM,PMID(I,J,1:LM),UH(I,J,1:LM)               &
! INPUT
                       ,VH(I,J,1:LM),T(I,J,1:LM),ZMID(I,J,1:LM)      &
! OUTPUT
                       ,MAXWP(I,J),MAXWU(I,J),MAXWV(I,J)             &
! OUTPUT
                       ,MAXWT(I,J),MAXWZ(I,J))
           ENDDO loopI
          END DO 
!        PRESSURE OF MAX WIND LEVEL
         IF (IGET(173) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = MAXWP(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(173))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
!        ICAO HEIGHT OF MAX WIND LEVEL
         IF (IGET(398)>0) THEN
           CALL ICAOHEIGHT(MAXWP, GRID1(1,jsta))
!            print*,'sample MAX WIND ICAO HEIGHTS',GRID1(im/2,(jsta+jend)/2)
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(398))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF
!        HEIGHT OF MAX WIND LEVEL
         IF (IGET(174) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = MAXWZ(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(174))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         ENDIF

!        MAX WIND LEVEL U WIND AND/OR V WIND.
        IF ((IGET(175) > 0).OR.(IGET(176) > 0)) THEN
!$omp parallel do private(i,j)
          DO J=JSTA,JEND
            DO I=1,IM
              GRID1(I,J) = MAXWU(I,J)
              GRID2(I,J) = MAXWV(I,J)
            ENDDO
          ENDDO
          if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(175))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(176))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID2(i,jj)
              enddo
            enddo
          endif
        ENDIF
!        TEMPERATURE OF MAX WIND LEVEL
        IF (IGET(314) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=MAXWT(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(314))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
        ENDIF
          deallocate(MAXWP,MAXWZ,MAXWU,MAXWV,MAXWT)
      ENDIF

!
!     ***BLOCK 3-1:  FD LEVEL (selected) T, Q, U, AND V.
!     
      IF ( (IGET(059)>0.or.IGET(586)>0).OR.IGET(911)>0.OR.     &
           (IGET(060)>0.or.IGET(576)>0).OR.                     &
           (IGET(061)>0.or.IGET(577)>0).OR.                     &
           (IGET(601)>0.or.IGET(602)>0.or.IGET(603)>0).OR.      &
           (IGET(604)>0.or.IGET(605)>0).OR.                     &
           (IGET(451)>0.or.IGET(578)>0).OR.IGET(580)>0 ) THEN

         ALLOCATE(T7D(IM,JSTA:JEND,NFD), Q7D(IM,JSTA:JEND,NFD),    &
                  U7D(IM,JSTA:JEND,NFD), V6D(IM,JSTA:JEND,NFD),    &
                  P7D(IM,JSTA:JEND,NFD), ICINGFD(IM,JSTA:JEND,NFD) &
                 ,AERFD(IM,JSTA:JEND,NFD,NBIN_DU))

!
!     DETERMINE WHETHER TO DO MSL OR AGL FD LEVELS
!
         ITYPEFDLVL=1
         DO IFD = 1,NFD
           IF (IGET(059)>0) THEN
            IF (LVLS(IFD,IGET(059))>1) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(911)>0) THEN
            IF (LVLS(IFD,IGET(911))>1) ITYPEFDLVL(IFD)=2
           ENDIF
!for grib2, spec hgt only
           IF (IGET(586)>0) THEN
            IF(LVLS(IFD,IGET(586))>0) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(060)>0) THEN
            IF (LVLS(IFD,IGET(060))>1) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(576)>0) THEN
            IF(LVLS(IFD,IGET(576))>0) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(061)>0) THEN
            IF (LVLS(IFD,IGET(061))>1) ITYPEFDLVL(IFD)=2
           ENDIF
           IF (IGET(577)>0) then
            if(LVLS(IFD,IGET(577))>0) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(451)>0) THEN
	    IF (LVLS(IFD,IGET(451))>1) ITYPEFDLVL(IFD)=2
	   ENDIF
           IF (IGET(578)>0) then
            if(LVLS(IFD,IGET(578))>0) ITYPEFDLVL(IFD)=2
           ENDIF
	   
	   IF (IGET(580)>0) then
            if(LVLS(IFD,IGET(580))>1) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(587)>0) then
            if(LVLS(IFD,IGET(587))>0) ITYPEFDLVL(IFD)=2
           ENDIF

	   IF (IGET(601)>0) THEN
            IF (LVLS(IFD,IGET(601))>1) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(602)>0) THEN
            IF (LVLS(IFD,IGET(602))>1) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(603)>0) THEN
            IF (LVLS(IFD,IGET(603))>1) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(604)>0) THEN
            IF (LVLS(IFD,IGET(604))>1) ITYPEFDLVL(IFD)=2
           ENDIF
	   IF (IGET(605)>0) THEN
            IF (LVLS(IFD,IGET(605))>1) ITYPEFDLVL(IFD)=2
           ENDIF

         ENDDO
!         print *,'call FDLVL with ITYPEFDLVL: ', ITYPEFDLVL,'for tmp,lvls=',LVLS(1:15,iget(59)), &
!          'grib2tmp lvs=',LVLS(1:15,iget(586))

         CALL FDLVL(ITYPEFDLVL,T7D,Q7D,U7D,V6D,P7D,ICINGFD,AERFD)
!     
         DO 10 IFD = 1,NFD
!
!           FD LEVEL TEMPERATURE.
            iget1 = IGET(059)
            iget2 = IGET(586)
            if (iget1 > 0) then
              work1 = LVLS(IFD,iget1)
            else
              work1 = 0.0
            endif
            if (iget2 > 0) then
              work2 = LVLS(IFD,iget2)
            else
              work2 = 0.0
            endif
            IF (IGET1 > 0 .or. IGET2 > 0) THEN
              IF (work1 > 0 .or. work2 > 0) THEN
             
!$omp parallel do private(i,j)
                DO J=JSTA,JEND
                  DO I=1,IM
                    GRID1(I,J) = T7D(I,J,IFD)
                  ENDDO
                ENDDO
                IF(work1 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET1)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET1)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                ENDIF
                IF (work2 > 0) THEN
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET2)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET2)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                ENDIF
              ENDIF
            ENDIF

!           FD LEVEL VIRTUAL TEMPERATURE.
            IF (IGET(911)>0) THEN
              IF (LVLS(IFD,IGET(911))>0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 if ( T7D(I,J,IFD) > 600 ) then
                 GRID1(I,J)=SPVAL
                 else
                 GRID1(I,J)=T7D(I,J,IFD)*(1.+0.608*Q7D(I,J,IFD))
                 endif
                 !print *, "grid value ",T7D(I,J,IFD),Q7D(I,J,IFD),T7D(I,J,IFD)*(1.+0.608*Q7D(I,J,IFD)),GRID1(I,J)
               ENDDO
               ENDDO
               IF(LVLS(IFD,IGET(911))>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(911))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(911))
                   datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
                 endif
               ENDIF
              ENDIF
            ENDIF

!
!           FD LEVEL SPEC HUMIDITY.
            iget1 = IGET(451)
            iget2 = iget(578)
            if (iget1 > 0) then
              work1 = LVLS(IFD,iget1)
            else
              work1 = 0.0
            endif
            if (iget2 > 0) then
              work2 = LVLS(IFD,iget2)
            else
              work2 = 0.0
            endif
            IF (IGET1 > 0 .or. iget2 > 0) THEN
              IF (work1 > 0 .or. work2 > 0)THEN
!$omp parallel do private(i,j)
                DO J=JSTA,JEND
                  DO I=1,IM
                    GRID1(I,J) = Q7D(I,J,IFD)
                  ENDDO
                ENDDO
                if(work1 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET1)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET1)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                endif
                if(work2 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET2)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET2)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                endif
              ENDIF
            ENDIF
!
!           FD LEVEL PRESSURE
            iget1 = IGET(482)
            iget2 = iget(579)
            if (iget1 > 0) then
              work1 = LVLS(IFD,iget1)
            else
              work1 = 0.0
            endif
            if (iget2 > 0) then
              work2 = LVLS(IFD,iget2)
            else
              work2 = 0.0
            endif
            IF (IGET1 > 0 .or. IGET2 > 0) THEN
              IF (work1 > 0 .or. work2 > 0) THEN
!$omp parallel do private(i,j)
                DO J=JSTA,JEND
                  DO I=1,IM
                    GRID1(I,J) = P7D(I,J,IFD)
                  ENDDO
                ENDDO
                if(work1 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET1)
                   fld_info(cfld)%lvl   = LVLSXML(IFD,IGET1)
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                     jj = jsta+j-1
                     do i=1,im
                       datapd(i,j,cfld) = GRID1(i,jj)
                     enddo
                   enddo
                  endif
                endif
                if(work2 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET2)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET2)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                endif
              ENDIF
            ENDIF
!
!           FD LEVEL ICING
            iget1 = IGET(580)
            iget2 = iget(587)
            if (iget1 > 0) then
              work1 = LVLS(IFD,iget1)
            else
              work1 = 0.0
            endif
            if (iget2 > 0) then
              work2 = LVLS(IFD,iget2)
            else
              work2 = 0.0
            endif
            IF (IGET1 > 0 .or. IGET2 > 0) THEN
              IF (work1 > 0 .or. work2 > 0) THEN
!$omp parallel do private(i,j)
                DO J=JSTA,JEND
                  DO I=1,IM
                    GRID1(I,J) = ICINGFD(I,J,IFD)
                  ENDDO
                ENDDO
                if(iget1 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET1)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET1)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                endif
                if(work2 > 0) then
                  if(grib == 'grib2') then
                    cfld = cfld + 1
                    fld_info(cfld)%ifld = IAVBLFLD(IGET2)
                    fld_info(cfld)%lvl  = LVLSXML(IFD,IGET2)
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                  endif
                endif

              ENDIF
            ENDIF
!
!  ADD FD LEVEL DUST/ASH (GOCART)
            IF (IGET(601)>0) THEN                      ! DUST 1
	      IF (LVLS(IFD,IGET(601))>0) THEN
!$omp parallel do private(i,j)
	       DO J=JSTA,JEND
	       DO I=1,IM
	          GRID1(I,J)=AERFD(I,J,IFD,1) 
               ENDDO
               ENDDO
               if(iget(601)>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(601))
                  fld_info(cfld)%lvl=LVLSXML(IFD,IGET(601))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               endif
	      ENDIF
	    ENDIF

            IF (IGET(602)>0) THEN			! DUST 2
	      IF (LVLS(IFD,IGET(602))>0) THEN
!$omp parallel do private(i,j)
	       DO J=JSTA,JEND
	       DO I=1,IM
	          GRID1(I,J)=AERFD(I,J,IFD,2) 
               ENDDO
               ENDDO
               if(iget(602)>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(602))
                  fld_info(cfld)%lvl=LVLSXML(IFD,IGET(602))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               endif
	      ENDIF
	    ENDIF

            IF (IGET(603)>0) THEN			! DUST 3
	      IF (LVLS(IFD,IGET(603))>0) THEN
!$omp parallel do private(i,j)
	       DO J=JSTA,JEND
	       DO I=1,IM
	          GRID1(I,J)=AERFD(I,J,IFD,3) 
               ENDDO
               ENDDO
               if(iget(603)>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(603))
                  fld_info(cfld)%lvl=LVLSXML(IFD,IGET(603))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               endif
	      ENDIF
	    ENDIF

            IF (IGET(604)>0) THEN			! DUST 4
	      IF (LVLS(IFD,IGET(604))>0) THEN
!$omp parallel do private(i,j)
	       DO J=JSTA,JEND
	       DO I=1,IM
	          GRID1(I,J)=AERFD(I,J,IFD,4) 
               ENDDO
               ENDDO
               if(iget(604)>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(604))
                  fld_info(cfld)%lvl=LVLSXML(IFD,IGET(604))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               endif
	      ENDIF
	    ENDIF

            IF (IGET(605)>0) THEN			! DUST 5
	      IF (LVLS(IFD,IGET(605))>0) THEN
!$omp parallel do private(i,j)
	       DO J=JSTA,JEND
	       DO I=1,IM
	          GRID1(I,J)=AERFD(I,J,IFD,5) 
               ENDDO
               ENDDO
               if(iget(605)>0) then
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(605))
                  fld_info(cfld)%lvl=LVLSXML(IFD,IGET(605))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               endif
	      ENDIF
	    ENDIF

!
!
!           FD LEVEL U WIND AND/OR V WIND.
            IF ((IGET(060)>0).OR.(IGET(061)>0)) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=U7D(I,J,IFD)
                 GRID2(I,J)=V6D(I,J,IFD)
               ENDDO
               ENDDO
               IF (IGET(060)>0) THEN
                 IF (LVLS(IFD,IGET(060))>0) then
                  if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(060))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(060))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                  endif
                 ENDIF
               ENDIF
               IF (IGET(061)>0) THEN
                 IF (LVLS(IFD,IGET(061))>0) THEN
                  if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(061))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(061))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID2(i,jj)
                    enddo
                  enddo
                  endif
                 ENDIF
               ENDIF
            ENDIF
!
!           FD LEVEL U WIND AND/OR V WIND.
            IF ((IGET(576)>0).OR.(IGET(577)>0)) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = U7D(I,J,IFD)
                   GRID2(I,J) = V6D(I,J,IFD)
                 ENDDO
               ENDDO
               IF (IGET(576)>0) THEN
                 IF (LVLS(IFD,IGET(576))>0) then
                  if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(576))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(576))
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                     jj = jsta+j-1
                     do i=1,im
                       datapd(i,j,cfld) = GRID1(i,jj)
                     enddo
                   enddo
                  endif
                 ENDIF
               ENDIF
               IF (IGET(577)>0) THEN
                 IF (LVLS(IFD,IGET(577))>0) THEN
                  if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(577))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(577))
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                     jj = jsta+j-1
                     do i=1,im
                       datapd(i,j,cfld) = GRID2(i,jj)
                     enddo
                   enddo
                  endif
                 ENDIF
               ENDIF
            ENDIF

 10      CONTINUE
         DEALLOCATE(T7D,Q7D,U7D,V6D,P7D,ICINGFD,AERFD)
      ENDIF

!
!     ***BLOCK 3-2:  FD LEVEL (from control file) GTG
!     
      IF(IGET(467)>0.or.IGET(468)>0.or.IGET(469)>0) THEN
         if(IGET(467)>0) THEN          ! GTG
            N=IAVBLFLD(IGET(467))
            NFDCTL=size(pset%param(N)%level)
            if(allocated(ITYPEFDLVLCTL)) deallocate(ITYPEFDLVLCTL)
            allocate(ITYPEFDLVLCTL(NFDCTL))
            DO IFD = 1,NFDCTL
               ITYPEFDLVLCTL(IFD)=LVLS(IFD,IGET(467))
            enddo
            if(allocated(HTFDCTL)) deallocate(HTFDCTL)
            allocate(HTFDCTL(NFDCTL))
            HTFDCTL=pset%param(N)%level
!           print *, "GTG 467 levels=",pset%param(N)%level
            allocate(GTGFD(IM,JSTA:JEND,NFDCTL))
            call FDLVL_MASS(ITYPEFDLVLCTL,NFDCTL,HTFDCTL,GTG,GTGFD)
!           print *, "GTG 467 Done GTGFD=",me,GTGFD(IM/2,jend,1:NFDCTL)
            DO IFD = 1,NFDCTL
              IF (LVLS(IFD,IGET(467))>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                 DO I=1,IM
                    GRID1(I,J)=GTGFD(I,J,IFD)
                 ENDDO
                 ENDDO
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(467))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(467))
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                         datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                   enddo
                 endif
              ENDIF
            ENDDO
         endif

         if(IGET(468)>0) THEN          ! CAT
            N=IAVBLFLD(IGET(468))
            NFDCTL=size(pset%param(N)%level)
            if(allocated(ITYPEFDLVLCTL)) deallocate(ITYPEFDLVLCTL)
            allocate(ITYPEFDLVLCTL(NFDCTL))
            DO IFD = 1,NFDCTL
               ITYPEFDLVLCTL(IFD)=LVLS(IFD,IGET(468))
            enddo
            if(allocated(HTFDCTL)) deallocate(HTFDCTL)
            allocate(HTFDCTL(NFDCTL))
            HTFDCTL=pset%param(N)%level
            allocate(CATFD(IM,JSTA:JEND,NFDCTL))
            call FDLVL_MASS(ITYPEFDLVLCTL,NFDCTL,HTFDCTL,catedr,CATFD)
            DO IFD = 1,NFDCTL
              IF (LVLS(IFD,IGET(468))>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                 DO I=1,IM
                    GRID1(I,J)=CATFD(I,J,IFD)
                 ENDDO
                 ENDDO
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(468))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(468))
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                         datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                   enddo
                 endif
              ENDIF
            ENDDO
         endif

         if(IGET(469)>0) THEN          ! MWT
            N=IAVBLFLD(IGET(469))
            NFDCTL=size(pset%param(N)%level)
            if(allocated(ITYPEFDLVLCTL)) deallocate(ITYPEFDLVLCTL)
            allocate(ITYPEFDLVLCTL(NFDCTL))
            DO IFD = 1,NFDCTL
               ITYPEFDLVLCTL(IFD)=LVLS(IFD,IGET(469))
            enddo
            if(allocated(HTFDCTL)) deallocate(HTFDCTL)
            allocate(HTFDCTL(NFDCTL))
            HTFDCTL=pset%param(N)%level
            allocate(MWTFD(IM,JSTA:JEND,NFDCTL))
            call FDLVL_MASS(ITYPEFDLVLCTL,NFDCTL,HTFDCTL,MWT,MWTFD)
            DO IFD = 1,NFDCTL
              IF (LVLS(IFD,IGET(469))>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                 DO I=1,IM
                    GRID1(I,J)=MWTFD(I,J,IFD)
                 ENDDO
                 ENDDO
                 if(grib=='grib2') then
                   cfld=cfld+1
                   fld_info(cfld)%ifld=IAVBLFLD(IGET(469))
                   fld_info(cfld)%lvl=LVLSXML(IFD,IGET(469))
!$omp parallel do private(i,j,jj)
                   do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                         datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                   enddo
                 endif
              ENDIF
            ENDDO
         endif

         if(allocated(GTGFD)) deallocate(GTGFD)
         if(allocated(CATFD)) deallocate(CATFD)
         if(allocated(MWTFD)) deallocate(MWTFD)

         if(allocated(ITYPEFDLVLCTL)) deallocate(ITYPEFDLVLCTL)
         if(allocated(HTFDCTL)) deallocate(HTFDCTL)

      ENDIF
!     
!
!     ***BLOCK 4:  FREEZING LEVEL Z, RH AND P.
!     
      IF ( (IGET(062)>0).OR.(IGET(063)>0) ) THEN
         CALL FRZLVL(Z1D,RH1D,P1D)
!
!        FREEZING LEVEL HEIGHT.
         IF (IGET(062)>0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J)=Z1D(I,J)
              ENDDO
            ENDDO
            CALL BOUND (GRID1,D00,H99999)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(062))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
         ENDIF
!
!        FREEZING LEVEL RELATIVE HUMIDITY.
         IF (IGET(063)>0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J) = RH1D(I,J)
               ENDDO
            ENDDO
            CALL SCLFLD(GRID1,H100,IM,JM)
            CALL BOUND(GRID1,H1,H100)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(063))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
         ENDIF
!
!        FREEZING LEVEL PRESSURE
         IF (IGET(753)>0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J) = P1D(I,J)
              ENDDO
            ENDDO
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(753))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif

         ENDIF
      ENDIF

      IF (IGET(165)>0 .OR. IGET(350)>0.OR. IGET(756)>0) THEN
         CALL FRZLVL2(TFRZ,Z1D,RH1D,P1D)
!
!        HIGHEST FREEZING LEVEL HEIGHT.
          IF (IGET(165)>0)THEN  
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=Z1D(I,J)
               ENDDO
               ENDDO
            CALL BOUND (GRID1,D00,H99999)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(165))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST FREEZING LEVEL RELATIVE HUMIDITY
          IF (IGET(350)>0)THEN  
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH1D(I,J)*100.
               ENDDO
               ENDDO
            CALL BOUND (GRID1,H1,H100)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(350))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST FREEZING LEVEL PRESSURE
         IF (IGET(756)>0) THEN
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J) = P1D(I,J)
              ENDDO
            ENDDO
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(756))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
         ENDIF

      ENDIF
!
! HIGHEST -10 C ISOTHERM VALUES
!
      IF (IGET(776)>0 .OR. IGET(777)>0.OR. IGET(778)>0) THEN 
         CALL FRZLVL2(263.15,Z1D,RH1D,P1D)
!
!        HIGHEST -10C ISOTHERM LEVEL HEIGHT.
          IF (IGET(776)>0)THEN  
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=Z1D(I,J)
               ENDDO
               ENDDO
            CALL BOUND (GRID1,D00,H99999)
            if(grib=='grib2') then 
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(776))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST -10C ISOTHERM RELATIVE HUMIDITY
          IF (IGET(777)>0)THEN  
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH1D(I,J)*100.
               ENDDO
               ENDDO
            CALL BOUND (GRID1,H1,H100)
            if(grib=='grib2') then 
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(777))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST -10C ISOTHERM LEVEL PRESSURE
         IF (IGET(778)>0) THEN 
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=P1D(I,J)
               ENDDO
               ENDDO
            if(grib=='grib2') then 
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(778))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
         ENDIF

      ENDIF
!
! HIGHEST -20 C ISOTHERM VALUES
!
      IF (IGET(779)>0 .OR. IGET(780)>0.OR. IGET(781)>0) THEN
         CALL FRZLVL2(253.15,Z1D,RH1D,P1D)
!
!        HIGHEST -20C ISOTHERM LEVEL HEIGHT.
          IF (IGET(779)>0)THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=Z1D(I,J)
               ENDDO
               ENDDO
            CALL BOUND (GRID1,D00,H99999)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(779))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST -20C ISOTHERM RELATIVE HUMIDITY
          IF (IGET(780)>0)THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RH1D(I,J)*100.
               ENDDO
               ENDDO
            CALL BOUND (GRID1,H1,H100)
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(780))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
          END IF

!        HIGHEST -20C ISOTHERM LEVEL PRESSURE
         IF (IGET(781)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=P1D(I,J)
               ENDDO
               ENDDO
            if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(781))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
            endif
         ENDIF

      ENDIF
!     
      allocate(PBND(IM,jsta:jend,NBND), TBND(IM,jsta:jend,NBND),    &
               QBND(IM,jsta:jend,NBND), UBND(IM,jsta:jend,NBND),    &
               VBND(IM,jsta:jend,NBND), RHBND(IM,jsta:jend,NBND),   &
               WBND(IM,jsta:jend,NBND))

!
!     ***BLOCK 5:  BOUNDARY LAYER FIELDS.
!     
      IF ( (IGET(067)>0).OR.(IGET(068)>0).OR.       &
           (IGET(069)>0).OR.(IGET(070)>0).OR.       &
           (IGET(071)>0).OR.(IGET(072)>0).OR.       &
           (IGET(073)>0).OR.(IGET(074)>0).OR.       &
           (IGET(088)>0).OR.(IGET(089)>0).OR.       &
           (IGET(090)>0).OR.(IGET(075)>0).OR.       &
           (IGET(109)>0).OR.(IGET(110)>0).OR.       &
           (IGET(031)>0).OR.(IGET(032)>0).OR.       &
           (IGET(573)>0).OR.                           &
           (IGET(107)>0).OR.(IGET(091)>0).OR.       &
           (IGET(092)>0).OR.(IGET(093)>0).OR.       &
           (IGET(094)>0).OR.(IGET(095)>0).OR.       &
           (IGET(096)>0).OR.(IGET(097)>0).OR.       &
           (IGET(098)>0).OR.(IGET(221)>0) ) THEN
!
           allocate(OMGBND(IM,jsta:jend,NBND),PWTBND(IM,jsta:jend,NBND),  &
                    QCNVBND(IM,jsta:jend,NBND),LVLBND(IM,jsta:jend,NBND), &
                    LB2(IM,jsta:jend))

!        COMPUTE ETA BOUNDARY LAYER FIELDS.
         CALL BNDLYR(PBND,TBND,QBND,RHBND,UBND,VBND,      &
                     WBND,OMGBND,PWTBND,QCNVBND,LVLBND)

!$omp parallel do private(i,j)
         DO J=JSTA,JEND
           DO I=1,IM
             EGRID2(i,j) = SPVAL     
           ENDDO
         ENDDO

!     
!        LOOP OVER NBND BOUNDARY LAYERS.
         DO 20 LBND = 1,NBND
!     
!           BOUNDARY LAYER PRESSURE.
            IF (IGET(067)>0) THEN
              IF (LVLS(LBND,IGET(067))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = PBND(I,J,LBND)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(067))
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(067))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
                endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER TEMPERATURE.
            IF (IGET(068)>0) THEN
              IF (LVLS(LBND,IGET(068))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=TBND(I,J,LBND)
               ENDDO
               ENDDO
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(068))
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(068))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER POTENTIAL TEMPERATURE.
            IF (IGET(069)>0) THEN
              IF (LVLS(LBND,IGET(069))>0) THEN
               CALL CALPOT(PBND(1,jsta,LBND),TBND(1,jsta,LBND),GRID1(1,jsta))
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(069))
                 fld_info(cfld)%lvl=LVLSXML(IFD,IGET(069))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
                endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER RELATIVE HUMIDITY.
            IF (IGET(072)>0) THEN
              IF (LVLS(LBND,IGET(072))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=RHBND(I,J,LBND)
               ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(072))
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(072))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER DEWPOINT TEMPERATURE.
            IF (IGET(070)>0) THEN
              IF (LVLS(LBND,IGET(070))>0) THEN
               CALL CALDWP(PBND(1,jsta,LBND), QBND(1,jsta,LBND),     &
                           GRID1(1,jsta),     TBND(1,jsta,LBND))
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(070))
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(070))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER SPECIFIC HUMIDITY.
            IF (IGET(071)>0) THEN
              IF (LVLS(LBND,IGET(071))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J)=QBND(I,J,LBND)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,H1M12,H99999)
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(071))
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(071))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER MOISTURE CONVERGENCE.
            IF (IGET(088)>0) THEN
              IF (LVLS(LBND,IGET(088))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = QCNVBND(I,J,LBND)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(088))
                 fld_info(cfld)%lvl=LVLSXML(LBND,IGET(088))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER U WIND AND/OR V WIND.
!
            FIELD1=.FALSE.
            FIELD2=.FALSE.
!
            IF(IGET(073)>0)THEN
              IF(LVLS(LBND,IGET(073))>0)FIELD1=.TRUE.
            ENDIF
            IF(IGET(074)>0)THEN
              IF(LVLS(LBND,IGET(074))>0)FIELD2=.TRUE.
            ENDIF
!
            IF(FIELD1.OR.FIELD2)THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = UBND(I,J,LBND)
                   GRID2(I,J) = VBND(I,J,LBND)
                 ENDDO
               ENDDO
!
               IF (IGET(073)>0) THEN
                 IF (LVLS(LBND,IGET(073))>0) then
                   if(grib=='grib2') then
                    cfld=cfld+1
                    fld_info(cfld)%ifld=IAVBLFLD(IGET(073))
                    fld_info(cfld)%lvl=LVLSXML(LBND,IGET(073))
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID1(i,jj)
                      enddo
                    enddo
                   endif
                 ENDIF
               ENDIF
               IF (IGET(074)>0) THEN
                 IF (LVLS(LBND,IGET(074))>0) THEN
                   if(grib=='grib2') then
                    cfld=cfld+1
                    fld_info(cfld)%ifld=IAVBLFLD(IGET(074))
                    fld_info(cfld)%lvl=LVLSXML(LBND,IGET(074))
!$omp parallel do private(i,j,jj)
                    do j=1,jend-jsta+1
                      jj = jsta+j-1
                      do i=1,im
                        datapd(i,j,cfld) = GRID2(i,jj)
                      enddo
                    enddo
                   endif
                 ENDIF
               ENDIF
            ENDIF
!     
!           BOUNDARY LAYER OMEGA.
            IF (IGET(090)>0) THEN
              IF (LVLS(LBND,IGET(090))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = OMGBND(I,J,LBND)
                 ENDDO
               ENDDO
              if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(090))
               fld_info(cfld)%lvl=LVLSXML(LBND,IGET(090))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
              endiF
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER PRECIPITBLE WATER.
            IF (IGET(089)>0) THEN
              IF (LVLS(LBND,IGET(089))>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = PWTBND(I,J,LBND)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
              if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(089))
               fld_info(cfld)%lvl=LVLSXML(LBND,IGET(089))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
              endif
              ENDIF
            ENDIF
!     
!           BOUNDARY LAYER LIFTED INDEX.
            IF (IGET(075)>0 .OR. IGET(031)>0 .OR. IGET(573)>0) THEN
             CALL OTLFT(PBND(1,jsta,LBND),TBND(1,jsta,LBND),    &
                    QBND(1,jsta,LBND),GRID1(1,jsta))
             IF(IGET(075)>0)THEN
              IF (LVLS(LBND,IGET(075))>0) THEN
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(075))
                fld_info(cfld)%lvl=LVLSXML(LBND,IGET(075))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
              ENDIF
             END IF
             IF(IGET(031)>0 .or. IGET(573)>0)THEN
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                DO I=1,IM
                  EGRID2(I,J) = MIN(EGRID2(I,J),GRID1(I,J))
                END DO
              END DO
             END IF
            ENDIF
!
!        END OF ETA BOUNDARY LAYER LOOP.
 20      CONTINUE
         deallocate(OMGBND,PWTBND,QCNVBND)
!     
!        BEST LIFTED INDEX FROM BOUNDARY LAYER FIELDS.
!     
         IF (IGET(031)>0 .OR. IGET(573)>0 ) THEN
!           DO J=JSTA,JEND
!            DO I=1,IM
!              EGRID1(I,J) = H99999
!              EGRID2(I,J) = H99999
!            ENDDO
!            ENDDO
!
!            DO 50 LBND = 1,NBND
!               CALL OTLFT(PBND(1,1,LBND),TBND(1,1,LBND),      &
!                    QBND(1,1,LBND),EGRID2)
!               DO J=JSTA,JEND
!               DO I=1,IM
!                 EGRID1(I,J)=AMIN1(EGRID1(I,J),EGRID2(I,J))
!               ENDDO
!               ENDDO
! 50         CONTINUE
!$omp parallel do private(i,j)
            DO J=JSTA,JEND
              DO I=1,IM
                GRID1(I,J)=EGRID2(I,J)
              ENDDO
            ENDDO
!	    print*,'writting out best lifted index'

            if (IGET(031)>0) then
              if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(031))
               datapd(1:im,1:jend-jsta+1,cfld)=GRID1(1:im,jsta:jend)
              endif
            endif
   
            if(IGET(573)> 0 ) THEN
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(573))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            endif

	 END IF
!     
!        BEST BOUNDARY LAYER CAPE AND CINS.
!     
         FIELD1=.FALSE.
         FIELD2=.FALSE.
!
         IF(IGET(032)>0)THEN
           IF(LVLS(2,IGET(032))>0)FIELD1=.TRUE.
         ENDIF
         IF(IGET(107)>0)THEN
           IF(LVLS(2,IGET(107))>0)FIELD2=.TRUE.
         ENDIF
!
         IF(IGET(566)>0)THEN
           FIELD1=.TRUE.
         ENDIF
         IF(IGET(567)>0)THEN
           FIELD2=.TRUE.
         ENDIF
!
         !if(grib=="grib2") print *,'in MISCLN.f,iget(566)=',          &
         !  iget(566), 'iget(567)=',iget(567),'LVLSXML(1,IGET(566)=',  &
         !  LVLSXML(1,IGET(566)),'LVLSXML(1,IGET(567)=',               &
         !  LVLSXML(1,IGET(567)),'field1=',field1,'field2=',field2
!
         IF(FIELD1.OR.FIELD2)THEN
           ITYPE = 2
!
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               EGRID1(I,J) = -H99999
               EGRID2(I,J) = -H99999
             ENDDO
           ENDDO
!
           DO 80 LBND = 1,NBND
           CALL CALTHTE(PBND(1,jsta,LBND),TBND(1,jsta,LBND),        &
                        QBND(1,jsta,LBND),EGRID1)
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               IF (EGRID1(I,J) > EGRID2(I,J)) THEN
                 EGRID2(I,J) = EGRID1(I,J)
                 LB2(I,J)    = LVLBND(I,J,LBND)
                 P1D(I,J)    = PBND(I,J,LBND)
                 T1D(I,J)    = TBND(I,J,LBND)
                 Q1D(I,J)    = QBND(I,J,LBND)
               ENDIF
             ENDDO
           ENDDO
 80        CONTINUE
!
           DPBND = 0.
           CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,LB2,EGRID1,   &
                        EGRID2,EGRID3,EGRID4,EGRID5) 
!
           IF (IGET(566)>0) THEN
! dong add missing value for cape
              GRID1=spval
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                DO I=1,IM
                  IF(T1D(I,J) < spval) GRID1(I,J) = EGRID1(I,J)
                ENDDO
              ENDDO
             CALL BOUND(GRID1,D00,H99999)
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(566))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(566))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
           ENDIF
!
           IF (IGET(567) > 0) THEN
! dong add missing value for cape
              GRID1=spval
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(T1D(I,J) < spval) GRID1(I,J) = - EGRID2(I,J)
               ENDDO
             ENDDO
!
             CALL BOUND(GRID1,D00,H99999)
!
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(T1D(I,J) < spval) GRID1(I,J) = - GRID1(I,J)
               ENDDO
             ENDDO
!
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(567))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(567))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
           ENDIF
         ENDIF
!

!    PBL HEIGHT 
         IF(IGET(221) > 0) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = PBLH(I,J)
             ENDDO
           ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(221))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
         END IF
!        BOUNDARY LAYER LIFTING CONDENSATION PRESSURE AND HEIGHT.
!        EGRID1 IS LCL PRESSURE.  EGRID2 IS LCL HEIGHT.
!
         IF ( (IGET(109)>0).OR.(IGET(110)>0) ) THEN
            CALL CALLCL(PBND(1,jsta,1),TBND(1,jsta,1),          &
                        QBND(1,jsta,1),EGRID1,EGRID2)
            IF (IGET(109)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = EGRID2(I,J)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(109))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
            IF (IGET(110)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = EGRID1(I,J)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(110))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
         ENDIF
!     
!        NGM BOUNDARY LAYER FIELDS.
!     
         IF ( (IGET(091)>0).OR.(IGET(092)>0).OR.      &
              (IGET(093)>0).OR.(IGET(094)>0).OR.      &
              (IGET(095)>0).OR.(IGET(095)>0).OR.      &
              (IGET(096)>0).OR.(IGET(097)>0).OR.      &
              (IGET(098)>0) ) THEN

              allocate(T78483(im,jsta:jend), T89671(im,jsta:jend), &
                       P78483(im,jsta:jend), P89671(im,jsta:jend))
!
!  COMPUTE SIGMA 0.89671 AND 0.78483 TEMPERATURES
!    INTERPOLATE LINEAR IN LOG P
            IF (IGET(097)>0.OR.IGET(098)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   P78483(I,J) = LOG(PINT(I,J,NINT(LMH(I,J)))*0.78483)
                   P89671(I,J) = LOG(PINT(I,J,NINT(LMH(I,J)))*0.89671)
                 ENDDO
               ENDDO
               DONE =.FALSE.
               DONE1=.FALSE.
!!$omp  parallel do private(fac1,fac2,pkl1,pku1,t78483,t89671)
               DO L=2,LM
                 DO J=JSTA,JEND
                   DO I=1,IM                  
                     PKL1=0.5*(ALPINT(I,J,L)+ALPINT(I,J,L+1))
                     PKU1=0.5*(ALPINT(I,J,L)+ALPINT(I,J,L-1))
!                     IF(I==1 .AND. J==1)PRINT*,'L,P89671,PKL1,PKU1= ', &
!                                             L,P89671(I,J), PKL1, PKU1
                     IF(P78483(I,J) < PKL1.AND.P78483(I,J) > PKU1) THEN
                      FAC1 = (PKL1-P78483(I,J))/(PKL1-PKU1)
                      FAC2 = (P78483(I,J)-PKU1)/(PKL1-PKU1)
                      T78483(I,J) = T(I,J,L)*FAC2 + T(I,J,L-1)*FAC1
                      DONE1(I,J)=.TRUE.
                     ENDIF
                     IF(P89671(I,J) < PKL1.AND.P89671(I,J) > PKU1)THEN
                       FAC1 = (PKL1-P89671(I,J))/(PKL1-PKU1)
                       FAC2 = (P89671(I,J)-PKU1)/(PKL1-PKU1)
                       T89671(I,J) = T(I,J,L)*FAC2 + T(I,J,L-1)*FAC1
                       DONE(I,J)   = .TRUE.
!	  	       IF(I==1 .AND. J==1)PRINT*,'done(1,1)= ',done(1,1)
                     ENDIF
                   ENDDO
                 ENDDO
               ENDDO
!	       print*,'done(1,1)= ',done(1,1)
!$omp parallel do private(i,j,pl,tl,ql,qsat,rhl)
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(.NOT. DONE(I,J)) THEN
                   PL    = PINT(I,J,LM-1)
                   TL   = 0.5*(T(I,J,LM-2)+T(I,J,LM-1))
                   QL   = 0.5*(Q(I,J,LM-2)+Q(I,J,LM-1))
                   QSAT = PQ0/PL *EXP(A2*(TL-A3)/(TL-A4))
!
                   RHL = QL/QSAT
!
                   IF(RHL > 1.)THEN
                    RHL = 1.
                    QL  = RHL*QSAT
                   ENDIF
!
                   IF(RHL < RHmin)THEN
                    RHL = RHmin
                    QL  = RHL*QSAT
                   ENDIF
!
!                  TVRL   = TL*(1.+0.608*QL)
!                  TVRBLO = TVRL*(P89671(I,J)/PL)**RGAMOG
!                  T89671(I,J) = TVRBLO/(1.+0.608*QL)

                   T89671(I,J) = TL * (P89671(I,J)/PL)**RGAMOG
!     

!                   PKL1=0.5*(ALPINT(I,J,LM)+ALPINT(I,J,LM+1))
!                   PKU1=0.5*(ALPINT(I,J,LM-1)+ALPINT(I,J,LM))
!                   T89671(I,J)=T(I,J,LM)+(T(I,J,LM)-T(I,J,LM-1))*
!     +               (P89671(I,J)-PKL1)/(PKL1-PKU1)

!                   print*,'Debug T89671= ',i,j
!     +		     ,P89671(I,J), PKL1, PKU1  
!     +               ,T89671(I,J),T(I,J,LM-1),T(I,J,LM)
                 END IF
 
                 IF(.NOT. DONE1(I,J))THEN
                   PL   = PINT(I,J,LM-1)
                   TL   = 0.5*(T(I,J,LM-2)+T(I,J,LM-1))
                   QL   = 0.5*(Q(I,J,LM-2)+Q(I,J,LM-1))
                   QSAT = PQ0/PL *EXP(A2*(TL-A3)/(TL-A4))
!
                   RHL = QL/QSAT
!
                   IF(RHL>1.)THEN
                    RHL = 1.
                    QL  = RHL*QSAT
                   ENDIF
!
                   IF(RHL<RHmin)THEN
                     RHL = RHmin
                     QL  = RHL*QSAT
                   ENDIF
!
!                  TVRL   = TL*(1.+0.608*QL)
!                  TVRBLO = TVRL*(P78483(I,J)/PL)**RGAMOG
!                  T78483(I,J)  =TVRBLO/(1.+0.608*QL)

                   T78483(I,J)  = TL * (P78483(I,J)/PL)**RGAMOG
!     
                 END IF
 
               END DO
             END DO
!     
!           SIGMA 0.89671 TEMPERATURE
             IF (IGET(097) > 0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = T89671(I,J)
!                  IF(T89671(I,J)>350.)PRINT*,'LARGE T89671 ',   &
!                    I,J,T89671(I,J)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(097))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(097))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
             ENDIF
!     
!           SIGMA 0.78483 TEMPERATURE
             IF (IGET(098)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = T78483(I,J)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(098))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(098))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
             ENDIF
             deallocate(T78483, T89671, P78483, P89671)
            ENDIF
!     
!           NGM SIGMA LAYER 0.98230 FIELDS.  THESE FIELDS ARE 
!           THE FIRST ETA LAYER BOUNDARY LAYER FIELDS. 
!     
!     
            IF ( (IGET(091)>0).OR.(IGET(092)>0).OR.      &
                 (IGET(093)>0).OR.(IGET(094)>0).OR.      &
                 (IGET(095)>0).OR.(IGET(095)>0).OR.      &
                 (IGET(096)>0) ) THEN
!     
!     
!              PRESSURE.
               IF (IGET(091)>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                   DO I=1,IM
                     GRID1(I,J) = PBND(I,J,1)
                   ENDDO
                 ENDDO
                 if(grib=='grib2') then
                  cfld=cfld+1
                  fld_info(cfld)%ifld=IAVBLFLD(IGET(091))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               ENDIF
!     
!              TEMPERATURE.
               IF (IGET(092)>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                   DO I=1,IM
                     GRID1(I,J) = TBND(I,J,1)
                   ENDDO
                 ENDDO
                 if(grib=='grib2') then
                  cfld=cfld+1
                  fld_info(cfld)%ifld=IAVBLFLD(IGET(092))
                  fld_info(cfld)%lvl=LVLSXML(1,IGET(092))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               ENDIF
!     
!              SPECIFIC HUMIDITY.
               IF (IGET(093)>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                   DO I=1,IM
                     GRID1(I,J) = QBND(I,J,1)
                   ENDDO
                 ENDDO
                  CALL BOUND(GRID1,H1M12,H99999)
                 if(grib=='grib2') then
                  cfld=cfld+1
                  fld_info(cfld)%ifld=IAVBLFLD(IGET(093))
                  fld_info(cfld)%lvl=LVLSXML(1,IGET(093))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
                 endif
               ENDIF
!     
!              RELATIVE HUMIDITY.
               IF (IGET(094)>0) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                   DO I=1,IM
                     GRID1(I,J) = RHBND(I,J,1)
                   ENDDO
                 ENDDO
                  CALL SCLFLD(GRID1,H100,IM,JM)
                  CALL BOUND(GRID1,H1,H100)
                if(grib=='grib2') then
                  cfld=cfld+1
                  fld_info(cfld)%ifld=IAVBLFLD(IGET(094))
                  fld_info(cfld)%lvl=LVLSXML(1,IGET(094))
!$omp parallel do private(i,j,jj)
                  do j=1,jend-jsta+1
                    jj = jsta+j-1
                    do i=1,im
                      datapd(i,j,cfld) = GRID1(i,jj)
                    enddo
                  enddo
               endif
               ENDIF
!     
!              U AND/OR V WIND.
               IF ((IGET(095)>0).OR.(IGET(096)>0)) THEN
!$omp parallel do private(i,j)
                 DO J=JSTA,JEND
                   DO I=1,IM
                     GRID1(I,J) = UBND(I,J,1)
                     GRID2(I,J) = VBND(I,J,1)
                   ENDDO
                 ENDDO
                  IF (IGET(095)>0)  then
                    if(grib=='grib2') then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(IGET(095))
                     fld_info(cfld)%lvl=LVLSXML(1,IGET(095))
!$omp parallel do private(i,j,jj)
                     do j=1,jend-jsta+1
                       jj = jsta+j-1
                       do i=1,im
                         datapd(i,j,cfld) = GRID1(i,jj)
                       enddo
                     enddo
                    endif
                  ENDIF
                  IF (IGET(096)>0) then
                    if(grib=='grib2') then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(IGET(096))
                     fld_info(cfld)%lvl=LVLSXML(1,IGET(096))
!$omp parallel do private(i,j,jj)
                     do j=1,jend-jsta+1
                       jj = jsta+j-1
                       do i=1,im
                         datapd(i,j,cfld) = GRID2(i,jj)
                       enddo
                     enddo
                    endif
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!     
!     ENDIF FOR BOUNDARY LAYER BLOCK.
!
      ENDIF
!     
!
!     ***BLOCK 6:  MISCELLANEOUS LAYER MEAN LFM AND NGM FIELDS.
!     
      IF ( (IGET(066)>0).OR.(IGET(081)>0).OR.        &
           (IGET(082)>0).OR.(IGET(104)>0).OR.        &
           (IGET(099)>0).OR.(IGET(100)>0).OR.        &
           (IGET(101)>0).OR.(IGET(102)>0).OR.        &
           (IGET(103)>0) ) THEN
!     
!        LFM "MEAN" RELATIVE HUMIDITIES AND PRECIPITABLE WATER.
!     
         IF ( (IGET(066)>0).OR.(IGET(081)>0).OR.     &
              (IGET(082)>0).OR.(IGET(104)>0) ) THEN
            allocate(RH3310(IM,jsta:jend),RH6610(IM,jsta:jend),          &
                     RH3366(IM,jsta:jend),PW3310(IM,jsta:jend))
            CALL LFMFLD(RH3310,RH6610,RH3366,PW3310)
!     
!           SIGMA 0.33-1.00 MEAN RELATIVE HUMIIDITY.
            IF (IGET(066)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH3310(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
              if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(066))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(066))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
!               print *,'in miscln,RH0.33-1.0,cfld=',cfld,'fld=',  &
!                IAVBLFLD(IGET(066)),'lvl=',fld_info(cfld)%lvl
              endif
            ENDIF
!     
!           SIGMA 0.66-1.00 MEAN RELATIVE HUMIIDITY.
            IF (IGET(081)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH6610(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(081))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(081))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.33-0.66 MEAN RELATIVE HUMIIDITY.
            IF (IGET(082)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH3366(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(082))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(082))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.33-1.00 PRECIPITABLE WATER.
            IF (IGET(104)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = PW3310(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(104))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(104))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
            deallocate(RH3310,RH6610,RH3366,PW3310)
         ENDIF
!     
!        VARIOUS LAYER MEAN NGM SIGMA FIELDS.
!     
         IF ( (IGET(099)>0).OR.(IGET(100)>0).OR.    &
              (IGET(101)>0).OR.(IGET(102)>0).OR.    &
              (IGET(103)>0) ) THEN
            allocate(RH4710(IM,jsta_2l:jend_2u),RH4796(IM,jsta_2l:jend_2u), &
                     RH1847(IM,jsta_2l:jend_2u))
            allocate(RH8498(IM,jsta_2l:jend_2u),QM8510(IM,jsta_2l:jend_2u))

            CALL NGMFLD(RH4710,RH4796,RH1847,RH8498,QM8510)
!     
!           SIGMA 0.47191-1.00000 RELATIVE HUMIDITY.
            IF (IGET(099)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH4710(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(099))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(099))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.47191-0.96470 RELATIVE HUMIDITY.
            IF (IGET(100)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH4796(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(100))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(100))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.18019-0.47191 RELATIVE HUMIDITY.
            IF (IGET(101)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH1847(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(101))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(101))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.84368-0.98230 RELATIVE HUMIDITY.
            IF (IGET(102)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH8498(I,J)
                 ENDDO
               ENDDO
               CALL SCLFLD(GRID1,H100,IM,JM)
               CALL BOUND(GRID1,H1,H100)
              if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(102))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(102))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
              endif
            ENDIF
!     
!           SIGMA 0.85000-1.00000 MOISTURE CONVERGENCE.
            IF (IGET(103)>0) THEN
!           CONVERT TO DIVERGENCE FOR GRIB
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = -1.0*QM8510(I,J)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(103))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(103))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
            deallocate(RH4710,RH4796,RH1847)
            deallocate(RH8498,QM8510)
         ENDIF
      ENDIF

      IF ( (IGET(318)>0).OR.(IGET(319)>0).OR.     &
           (IGET(320)>0))THEN
       allocate(RH4410(IM,jsta:jend),RH7294(IM,jsta:jend),   &
                RH4472(IM,jsta:jend),RH3310(IM,jsta:jend))
       CALL LFMFLD_GFS(RH4410,RH7294,RH4472,RH3310)
!     
!           SIGMA 0.44-1.00 MEAN RELATIVE HUMIIDITY.
            IF (IGET(318)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH4410(I,J)*100.
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(318))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(318))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.72-0.94 MEAN RELATIVE HUMIIDITY.
            IF (IGET(319)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = RH7294(I,J)*100.
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(319))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(319))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
!     
!           SIGMA 0.44-0.72 MEAN RELATIVE HUMIIDITY.
            IF (IGET(320)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J)=RH4472(I,J)*100.
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H100)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(320))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(320))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
            ENDIF
            deallocate(RH4410,RH7294,RH4472,RH3310)
      END IF

! GFS computes sigma=0.9950 T, THETA, U, V from lowest two model level fields 
         IF ( (IGET(321)>0).OR.(IGET(322)>0).OR.     &
              (IGET(323)>0).OR.(IGET(324)>0).OR.     &
              (IGET(325)>0).OR.(IGET(326)>0)) THEN
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
	     DO I=1,IM
               EGRID2(I,J) = 0.995*PINT(I,J,LM+1)
               EGRID1(I,J) = LOG(PMID(I,J,LM)/EGRID2(I,J))   &
                           / LOG(PMID(I,J,LM)/PMID(I,J,LM-1))

        IF (MODELNAME == 'GFS' .OR. MODELNAME == 'FV3R') THEN
               EGRID1(I,J) = LOG(PMID(I,J,LM)/EGRID2(I,J))   &
                           / max(1.e-6,LOG(PMID(I,J,LM)/PMID(I,J,LM-1)))
               EGRID1(I,J) =max(-10.0,min(EGRID1(I,J), 10.0))
                   IF ( ABS(PMID(I,J,LM)-PMID(I,J,LM-1)) < 0.5 ) THEN
                     EGRID1(I,J) = -1.
                   ENDIF
        ENDIF

	     END DO
	   END DO
! Temperature	   
	   IF (IGET(321)>0) THEN
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = T(I,J,LM)+(T(I,J,LM-1)-T(I,J,LM)) &
                            * EGRID1(I,J)
               ENDDO
             ENDDO
           if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=IAVBLFLD(IGET(321))
            fld_info(cfld)%lvl=LVLSXML(1,IGET(321))
!$omp parallel do private(i,j,jj)
            do j=1,jend-jsta+1
              jj = jsta+j-1
              do i=1,im
                datapd(i,j,cfld) = GRID1(i,jj)
              enddo
            enddo
           endif
!            print *,'in miscln,iget(321,temp,sigmadata=',maxval(GRID1(1:im,jsta:jend)), &
!             minval(GRID1(1:im,jsta:jend)),'grib=',grib
            ENDIF
! Potential Temperature	    
            IF (IGET(322)>0) THEN
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID2(I,J) = T(I,J,LM)+(T(I,J,LM-1)-T(I,J,LM))     &
                            * EGRID1(I,J)
               ENDDO
             ENDDO
             CALL CALPOT(EGRID2,GRID2(1,jsta),GRID1(1,jsta))
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(322))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(322))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF
! RH	    
            IF (IGET(323)>0) THEN
!$omp parallel do private(i,j,es1,qs1,rh1,es2,qs2,rh2)
             DO J=JSTA,JEND
               DO I=1,IM
                 ES1 = min(PMID(I,J,LM),FPVSNEW(T(I,J,LM)))
                 QS1 = CON_EPS*ES1/(PMID(I,J,LM)+CON_EPSM1*ES1)
                 RH1 = Q(I,J,LM)/QS1
                 ES2 = min(PMID(I,J,LM-1),FPVSNEW(T(I,J,LM-1)))
                 QS2 = CON_EPS*ES2/(PMID(I,J,LM-1)+CON_EPSM1*ES2)
                 RH2 = Q(I,J,LM-1)/QS2
                 GRID1(I,J) = (RH1+(RH2-RH1)*EGRID1(I,J))*100.
               ENDDO
             ENDDO
             CALL BOUND(GRID1,D00,H100)
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(323))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(323))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF    
! U	   
            IF (IGET(324)>0) THEN
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = UH(I,J,LM)+(UH(I,J,LM-1)-UH(I,J,LM))    &
                            * EGRID1(I,J)
               ENDDO
             ENDDO
            if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(324))
             fld_info(cfld)%lvl=LVLSXML(1,IGET(324))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
            ENDIF
! V	   
            IF (IGET(325)>0) THEN
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = VH(I,J,LM)+(VH(I,J,LM-1)-VH(I,J,LM))    &
                            * EGRID1(I,J)
               ENDDO
             ENDDO
            if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(325))
             fld_info(cfld)%lvl=LVLSXML(1,IGET(325))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
              jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
            endif
           ENDIF
! OMGA	   
           IF (IGET(326)>0) THEN
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = OMGA(I,J,LM)+(OMGA(I,J,LM-1)-OMGA(I,J,LM))  &
                            * EGRID1(I,J)
               ENDDO
             ENDDO
            if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(IGET(326))
             fld_info(cfld)%lvl=LVLSXML(1,IGET(326))
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
              jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
            endif
           ENDIF
         END IF

!       MIXED LAYER CAPE AND CINS
!
         FIELD1=.FALSE.
         FIELD2=.FALSE.
!
         IF(IGET(032)>0)THEN
           IF(LVLS(3,IGET(032))>0)FIELD1=.TRUE.
         ENDIF
         IF(IGET(107)>0)THEN
           IF(LVLS(3,IGET(107))>0)FIELD2=.TRUE.
         ENDIF
!
         IF(IGET(582)>0)THEN
           FIELD1=.TRUE.
         ENDIF
         IF(IGET(583)>0)THEN
           FIELD2=.TRUE.
         ENDIF

         IF(FIELD1.OR.FIELD2)THEN
           ITYPE = 2
!
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               EGRID1(I,J) = -H99999
               EGRID2(I,J) = -H99999
!          ENDDO
!          ENDDO
!          DO J=JSTA,JEND
!          DO I=1,IM
               LB2(I,J)  = (LVLBND(I,J,1) + LVLBND(I,J,2) +           &
                            LVLBND(I,J,3))/3
               P1D(I,J)  = (PBND(I,J,1) + PBND(I,J,2) + PBND(I,J,3))/3
               T1D(I,J)  = (TBND(I,J,1) + TBND(I,J,2) + TBND(I,J,3))/3
               Q1D(I,J)  = (QBND(I,J,1) + QBND(I,J,2) + QBND(I,J,3))/3
             ENDDO
           ENDDO
!
           DPBND = 0.
           CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,LB2,EGRID1,           &
                        EGRID2,EGRID3,EGRID4,EGRID5)
 
           IF (IGET(582)>0) THEN
! dong add missing value for cape
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) GRID1(I,J) = EGRID1(I,J)
                 ENDDO
               ENDDO

               CALL BOUND(GRID1,D00,H99999)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(582))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(582))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                 jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif
           ENDIF
           IF (IGET(583)>0) THEN
! dong add missing value for cape
               GRID1=spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) GRID1(I,J) = - EGRID2(I,J)
                 ENDDO
               ENDDO
!
               CALL BOUND(GRID1,D00,H99999)
!
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) GRID1(I,J) = - GRID1(I,J)
                 ENDDO
               ENDDO
!
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(583))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(583))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif

           ENDIF
         ENDIF
              
!        MIXED LAYER LIFTING CONDENSATION PRESSURE AND HEIGHT.
!        EGRID1 IS LCL PRESSURE.  EGRID2 IS LCL HEIGHT.
!
!         IF ( (IGET(109)>0).OR.(IGET(110)>0) ) THEN
!            CALL CALLCL(P1D,T1D,Q1D,EGRID1,EGRID2)
!            IF (IGET(109)>0) THEN
!	       DO J=JSTA,JEND
!               DO I=1,IM
!                 GRID1(I,J)=EGRID2(I,J)
!               ENDDO
!               ENDDO
!           
!               ID(1:25) = 0
!	       
!	       CALL GRIBIT(IGET(109),1,
!     X              GRID1,IM,JM)
!            ENDIF
!	    
!            IF (IGET(110)>0) THEN
!	       DO J=JSTA,JEND
!               DO I=1,IM
!                 GRID1(I,J)=EGRID1(I,J)
!               ENDDO
!               ENDDO
!	       
!               ID(1:25) = 0
!	       
!	       CALL GRIBIT(IGET(110),1,
!     X              GRID1,IM,JM)
!            ENDIF
!         ENDIF
!
!       MOST UNSTABLE CAPE-LOWEST 300 MB
!
         FIELD1=.FALSE.
         FIELD2=.FALSE.
!
         IF(IGET(032)>0)THEN
           IF(LVLS(4,IGET(032))>0)FIELD1=.TRUE.
         ENDIF
              
         IF(IGET(107)>0)THEN
           IF(LVLS(4,IGET(107))>0)FIELD2=.TRUE.
         ENDIF
!
         IF(IGET(584)>0)THEN
           FIELD1=.TRUE.
         ENDIF
         IF(IGET(585)>0)THEN
           FIELD2=.TRUE.
         ENDIF
!
         IF(FIELD1.OR.FIELD2)THEN
           ITYPE = 1
!
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               EGRID1(I,J) = -H99999
               EGRID2(I,J) = -H99999
             ENDDO
           ENDDO
              
           DPBND = 300.E2
           CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,LB2,EGRID1,     &
                        EGRID2,EGRID3,EGRID4,EGRID5)
!
           IF (IGET(584)>0) THEN
! dong add missing value to cin
               GRID1 = spval
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) GRID1(I,J) = EGRID1(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(584))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(584))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif

           ENDIF
                
           IF (IGET(585)>0) THEN
! dong add missing value to cin
               GRID1 = spval
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) GRID1(I,J) = - EGRID2(I,J)
                 ENDDO
               ENDDO

               CALL BOUND(GRID1,D00,H99999)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF(T1D(I,J) < spval) GRID1(I,J) = - GRID1(I,J)
                 ENDDO
               ENDDO
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(585))
                 fld_info(cfld)%lvl=LVLSXML(1,IGET(585))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif

            ENDIF
              
!    EQUILLIBRIUM HEIGHT
           IF (IGET(443)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = EGRID4(I,J)
                 ENDDO
               ENDDO
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(443))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(443))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
           ENDIF

!      PRESSURE OF LEVEL FROM WHICH 300 MB MOST UNSTABLE CAPE
!      PARCEL WAS LIFTED (eq. PRESSURE OF LEVEL OF HIGHEST THETA-E)
           IF (IGET(246)>0) THEN
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = EGRID3(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
!               print *,'in miscln,PLPL=',maxval(grid1(1:im,jsta:jend)),  &
!                 minval(grid1(1:im,jsta:jend))
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(246))
                 fld_info(cfld)%lvl=LVLSXML(1,IGET(246))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
           ENDIF   ! 246

!    GENERAL THUNDER PARAMETER  ??? 458 ???
        IF (IGET(444)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   IF (CPRATE(I,J) > PTHRESH) THEN
                    GRID1(I,J) = EGRID5(I,J)
                   ELSE
                    GRID1(I,J) = 0
                   ENDIF
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               if(grib=='grib2') then
                 cfld=cfld+1
                 fld_info(cfld)%ifld=IAVBLFLD(IGET(444))
                 fld_info(cfld)%lvl=LVLSXML(1,IGET(444))
!$omp parallel do private(i,j,jj)
                 do j=1,jend-jsta+1
                   jj = jsta+j-1
                   do i=1,im
                     datapd(i,j,cfld) = GRID1(i,jj)
                   enddo
                 enddo
               endif
        ENDIF
      ENDIF
!
!       EXPAND HRRR CAPE/CIN RELATED VARIABLES
!       
!    CAPE AND CINS 0-3KM, FOLLOW ML PROCEDURE WITH HEIGHT 0-3KM

         FIELD1=.FALSE.
         FIELD2=.FALSE.
!
         IF(IGET(032)>0)THEN
           IF(LVLS(3,IGET(032))>0)FIELD1=.TRUE.
         ENDIF
         IF(IGET(107)>0)THEN
           IF(LVLS(3,IGET(107))>0)FIELD2=.TRUE.
         ENDIF
!
         IF(IGET(950)>0)THEN
           FIELD1=.TRUE.
         ENDIF
         IF(IGET(951)>0)THEN
           FIELD2=.TRUE.
         ENDIF
!
         IF(FIELD1.OR.FIELD2)THEN
           ITYPE = 2
!
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               EGRID1(I,J) = -H99999
               EGRID2(I,J) = -H99999
               EGRID3(I,J) = -H99999
               EGRID4(I,J) = -H99999
               EGRID5(I,J) = -H99999
               EGRID6(I,J) = -H99999
               EGRID7(I,J) = -H99999
               EGRID8(I,J) = -H99999
!          ENDDO
!          ENDDO
!          DO J=JSTA,JEND
!          DO I=1,IM
               LB2(I,J)  = (LVLBND(I,J,1) + LVLBND(I,J,2) +           &
                            LVLBND(I,J,3))/3
               P1D(I,J)  = (PBND(I,J,1) + PBND(I,J,2) + PBND(I,J,3))/3
               T1D(I,J)  = (TBND(I,J,1) + TBND(I,J,2) + TBND(I,J,3))/3
               Q1D(I,J)  = (QBND(I,J,1) + QBND(I,J,2) + QBND(I,J,3))/3
             ENDDO
           ENDDO
!
           DPBND = 0.
           CALL CALCAPE2(ITYPE,DPBND,P1D,T1D,Q1D,LB2,            &
                         EGRID1,EGRID2,EGRID3,EGRID4,EGRID5,     &
                         EGRID6,EGRID7,EGRID8)
!
!                        CAPE1, CINS2, LFC3,  ESRHL4,ESRHH5,
!                        DCAPE6,DGLD7, ESP8)
!
           IF (IGET(950)>0) THEN
! dong add missing value for cape
              GRID1=spval
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                DO I=1,IM
                  IF(T1D(I,J) < spval) GRID1(I,J) = EGRID1(I,J)
                ENDDO
              ENDDO
             CALL BOUND(GRID1,D00,H99999)
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(950))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(950))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
           ENDIF   !950
!
           IF (IGET(951)>0) THEN
! dong add missing value for cape
              GRID1=spval
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(T1D(I,J) < spval) GRID1(I,J) = - EGRID2(I,J)
               ENDDO
             ENDDO
!
             CALL BOUND(GRID1,D00,H99999)
!
!$omp parallel do private(i,j)
             DO J=JSTA,JEND
               DO I=1,IM
                 IF(T1D(I,J) < spval) GRID1(I,J) = - GRID1(I,J)
               ENDDO
             ENDDO
!
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(951))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(951))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF   !951

!    LFC HEIGHT

            IF (IGET(952)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = EGRID3(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(952))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(952))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
            ENDIF   !952

!    EFFECTIVE STORM RELATIVE HELICITY AND STORM MOTION.

         allocate(USHR1(IM,jsta_2l:jend_2u),VSHR1(IM,jsta_2l:jend_2u), &
                  USHR6(IM,jsta_2l:jend_2u),VSHR6(IM,jsta_2l:jend_2u))
         allocate(UST(IM,jsta_2l:jend_2u),VST(IM,jsta_2l:jend_2u),     &
                  HELI(IM,jsta_2l:jend_2u,2))
         allocate(LLOW(IM,jsta_2l:jend_2u),LUPP(IM,jsta_2l:jend_2u),   &
                  CANGLE(IM,jsta_2l:jend_2u))

       iget1 = IGET(953)
       iget2 = -1
       iget3 = -1
       if (iget1 > 0) then
         iget2 = LVLS(1,iget1)
         iget3 = LVLS(2,iget1)
       endif
      if(me==0) write(0,*) '953 ',iget1,iget2,iget3
       IF (iget1 > 0 .OR. IGET(162) > 0 .OR. IGET(953) > 0) THEN
         DEPTH(1) = 3000.0
         DEPTH(2) = 1000.0
!$omp parallel do private(i,j)
         DO J=JSTA,JEND
           DO I=1,IM
             LLOW(I,J) = INT(EGRID4(I,J))
             LUPP(I,J) = INT(EGRID5(I,J))
           ENDDO
         ENDDO

!         CALL CALHEL(DEPTH,UST,VST,HELI,USHR1,VSHR1,USHR6,VSHR6)
         CALL CALHEL2(LLOW,LUPP,DEPTH,UST,VST,HELI,CANGLE)
!                                                  
         IF (iget2 > 0) then
!$omp parallel do private(i,j)
           DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = HELI(I,J,1)
             !  GRID1(I,J) = HELI(I,J,2)
             ENDDO
           ENDDO
           if(grib=='grib2') then
             cfld=cfld+1
             fld_info(cfld)%ifld=IAVBLFLD(iget1)
             fld_info(cfld)%lvl=LVLSXML(1,iget1)
!$omp parallel do private(i,j,jj)
             do j=1,jend-jsta+1
               jj = jsta+j-1
               do i=1,im
                 datapd(i,j,cfld) = GRID1(i,jj)
               enddo
             enddo
           endif
         ENDIF

       ENDIF   !953

        IF (SUBMODELNAME == 'RTMA') THEN  !Start RTMA block

!EL field allocation

         allocate(ESHR(IM,jsta_2l:jend_2u),UVECT(IM,jsta_2l:jend_2u),&
                  VVECT(IM,jsta_2l:jend_2u),HTSFC(IM,jsta_2l:jend_2u))
!get surface height
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

      IF(gridtype /= 'A') CALL EXCH(FIS(1:IM,JSTA:JEND))
        DO J=JSTART,JSTOP
          DO I=ISTART,ISTOP
            IE = I+IVE(J)
            IW = I+IVW(J)
            JN = J+JVN
            JS = J+JVS
            IF (gridtype=='B')THEN
            HTSFC(I,J)=(0.25/g)*(FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(IE,JN))
            ELSE
            HTSFC(I,J)=(0.25/g)*(FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(I,JS))
            ENDIF
          ENDDO
        ENDDO

!Height of effbot
            IF (IGET(979)>0) THEN
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = ZINT(I,J,LLOW(I,J)) - HTSFC(I,J)
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(979))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(979))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF
!Height of effbot
            IF (IGET(980)>0) THEN
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = ZINT(I,J,LUPP(I,J)) - HTSFC(I,J)
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(980))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(980))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF
!Temperature of effbot
            IF (IGET(981)>0) THEN
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = T(I,J,LLOW(I,J))
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(981))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(981))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF
!Temperature of efftop
            IF (IGET(982)>0) THEN
             DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = T(I,J,LUPP(I,J))
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(982))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(982))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF

!U inflow based to 50% EL shear vector

            IF (IGET(983)>0) THEN
             DO J=JSTA,JEND
               DO I=1,IM
                       MIDCAL=INT(LLOW(I,J)+D50*(LUPP(I,J)-LLOW(I,J)))       
                                                            !mid-layer 
                                                            !vertical
                                                            !index
                       UVECT(I,J)=UH(I,J,MIDCAL)-UH(I,J,LLOW(I,J))
                       GRID1(I,J)=UVECT(I,J)
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(983))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(983))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF

!V inflow based to 50% EL shear vector
            IF (IGET(984)>0) THEN
             DO J=JSTA,JEND
               DO I=1,IM
                       MIDCAL=INT(LLOW(I,J)+D50*(LUPP(I,J)-LLOW(I,J)))
                                                            !mid-layer 
                                                            !vertical
                                                            !index
                       VVECT(I,J)=VH(I,J,MIDCAL)-VH(I,J,LLOW(I,J))
                       GRID1(I,J)=VVECT(I,J)
               ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(984))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(984))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF

!Inflow based (ESFC) to (50%) EL shear magnitude
            IF (IGET(985)>0) THEN
             DO J=JSTA,JEND
               DO I=1,IM
                       ESHR(I,J)=SQRT((UVECT(I,J)**2)+(VVECT(I,J))**2)
                                                               !effshear
                                                               !calc
                       GRID1(I,J)=ESHR(I,J) !Effective
             ENDDO
             ENDDO
             if(grib=='grib2') then
              cfld=cfld+1
              fld_info(cfld)%ifld=IAVBLFLD(IGET(985))
              fld_info(cfld)%lvl=LVLSXML(1,IGET(985))
!$omp parallel do private(i,j,jj)
              do j=1,jend-jsta+1
                jj = jsta+j-1
                do i=1,im
                  datapd(i,j,cfld) = GRID1(i,jj)
                enddo
              enddo
             endif
            ENDIF


        ENDIF   !END RTMA BLOCK


!    Critical Angle

            IF (IGET(957)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                    GRID1(I,J) = CANGLE(I,J)   
           !         IF(EGRID1(I,J)<100. .OR. EGRID2(I,J)>-250.) THEN
           !           GRID1(I,J) = 0.
           !         ENDIF
                 ENDDO
               ENDDO
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(957))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(957))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
            ENDIF   !957

!    Dendritic Layer Depth, -17C < T < -12C

            IF (IGET(955)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = EGRID7(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(955))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(955))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
            ENDIF   !955

!    Enhanced Stretching Potential

            IF (IGET(956)>0) THEN
!$omp parallel do private(i,j)
               DO J=JSTA,JEND
                 DO I=1,IM
                   GRID1(I,J) = EGRID8(I,J)
                 ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
             if(grib=='grib2') then
               cfld=cfld+1
               fld_info(cfld)%ifld=IAVBLFLD(IGET(956))
               fld_info(cfld)%lvl=LVLSXML(1,IGET(956))
!$omp parallel do private(i,j,jj)
               do j=1,jend-jsta+1
                 jj = jsta+j-1
                 do i=1,im
                   datapd(i,j,cfld) = GRID1(i,jj)
                 enddo
               enddo
             endif
            ENDIF   !956

!    Downdraft CAPE

           ITYPE = 1
          ! DO J=JSTA,JEND
          ! DO I=1,IM
          !     LB2(I,J)  = (LVLBND(I,J,1) + LVLBND(I,J,2) +           &
          !                  LVLBND(I,J,3))/3
          !     P1D(I,J)  = (PBND(I,J,1) + PBND(I,J,2) + PBND(I,J,3))/3
          !     T1D(I,J)  = (TBND(I,J,1) + TBND(I,J,2) + TBND(I,J,3))/3
          !     Q1D(I,J)  = (QBND(I,J,1) + QBND(I,J,2) + QBND(I,J,3))/3
          !   ENDDO
          ! ENDDO

           DPBND = 400.E2
          ! CALL CALCAPE2(ITYPE,DPBND,P1D,T1D,Q1D,LB2,            &
          !               EGRID1,EGRID2,EGRID3,EGRID4,EGRID5,     &
          !               EGRID6,EGRID7,EGRID8)

           IF (IGET(954)>0) THEN
               GRID1 = spval
!$omp parallel do private(i,j)
              DO J=JSTA,JEND
                 DO I=1,IM
                  IF(T1D(I,J) < spval) GRID1(I,J) = -EGRID6(I,J)
                 ENDDO
              ENDDO
               CALL BOUND(GRID1,D00,H99999)
               if(grib=='grib2') then
                cfld=cfld+1
                fld_info(cfld)%ifld=IAVBLFLD(IGET(954))
                fld_info(cfld)%lvl=LVLSXML(1,IGET(954))
!$omp parallel do private(i,j,jj)
                do j=1,jend-jsta+1
                  jj = jsta+j-1
                  do i=1,im
                    datapd(i,j,cfld) = GRID1(i,jj)
                  enddo
                enddo
               endif

           ENDIF   !954

       if (allocated(ushr1)) deallocate(ushr1)
       if (allocated(vshr1)) deallocate(vshr1)
       if (allocated(ushr6)) deallocate(ushr6)
       if (allocated(vshr6)) deallocate(vshr6)
       if (allocated(ust))   deallocate(ust)
       if (allocated(vst))   deallocate(vst)
       if (allocated(heli))  deallocate(heli)
       if (allocated(llow))  deallocate(llow)
       if (allocated(lupp))  deallocate(lupp)
       if (allocated(cangle))deallocate(cangle)

       ENDIF

      if (allocated(pbnd))   deallocate(pbnd)
      if (allocated(tbnd))   deallocate(tbnd)
      if (allocated(qbnd))   deallocate(qbnd)
      if (allocated(ubnd))   deallocate(ubnd)
      if (allocated(vbnd))   deallocate(vbnd)
      if (allocated(rhbnd))  deallocate(rhbnd)
      if (allocated(wbnd))   deallocate(wbnd)
      if (allocated(lvlbnd)) deallocate(lvlbnd)
      if (allocated(lb2))    deallocate(lb2)
!    
!
! RELATIVE HUMIDITY WITH RESPECT TO PRECIPITABLE WATER
       IF (IGET(749)>0) THEN
          CALL CALRH_PW(GRID1(1,jsta))
          if(grib=='grib2') then
           cfld=cfld+1
           fld_info(cfld)%ifld=IAVBLFLD(IGET(749))
!$omp parallel do private(i,j,jj)
           do j=1,jend-jsta+1
             jj = jsta+j-1
             do i=1,im
               datapd(i,j,cfld) = GRID1(i,jj)
             enddo
           enddo
          endif
       ENDIF       

 
!     END OF ROUTINE.
!     
      RETURN
      END
