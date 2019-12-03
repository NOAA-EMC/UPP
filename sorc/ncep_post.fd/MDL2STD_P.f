      SUBROUTINE MDL2STD_P()
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    MDL2STD_P       VERT INTRP OF MODEL LVLS TO STANDARD ATMOSPEHRIC PRESSURE
!   PRGRMMR: Y Mao           ORG: W/NP22     DATE: Sep 2019
!     
! ABSTRACT:
!     ORIGINATED FROM MISCLN.f. THIS ROUTINE INTERPOLATE TO STANDARD
!     ATMOSPHERIC PRESSURE, INSTEAD OF MODEL PRESSURE
!     
! PROGRAM HISTORY LOG:
!   19-09-24  Y Mao       - REWRITTEN FROM MISCLN.f
!
! USAGE:    CALL MDL2STD_P
!   INPUT ARGUMENT LIST:
!
!   OUTPUT ARGUMENT LIST: 
!     NONE       
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       FDLVL_UV   - COMPUTE FD LEVEL WIND (AGL OR MSL).
!       FDLVL_MASS - COMPUTE FD LEVEL MASS (AGL OR MSL).
!
!     LIBRARY:
!       COMMON   - CTLBLK
!                  RQSTFLD
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN 90
!     MACHINE : IBM SP
!$$$  
!
      use vrbls3d, only: pint, pmid, zmid
      use vrbls3d, only: t, q, uh, vh, omga, cwm, qqw, qqi, qqr, qqs, qqg
      
      use vrbls3d, only: ICING_GFIP, ICING_GFIS, catedr, mwt, gtg
      use ctlblk_mod, only: grib, cfld, fld_info, datapd, im, jsta, jend, jm, &
                            lm, htfd, spval, nfd, me,&
                            jsta_2l, jend_2u, MODELNAME
      use rqstfld_mod, only: iget, lvls, iavblfld, lvlsxml
      use grib2_module, only: pset

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
      implicit none

      real, external :: P2H

      real,dimension(im,jm)        :: GRID1
      real,dimension(im,jsta_2l:jend_2u) :: EGRID1,EGRID2,EGRID3,EGRID4

!
      integer I,J,jj,L,ITYPE,IFD,ITYPEFDLVL(NFD)

!     Variables introduced to allow FD levels from control file - Y Mao
      integer :: N,NFDCTL
      REAL, allocatable :: HTFDCTL(:)
      integer, allocatable :: ITYPEFDLVLCTL(:)
      real, allocatable :: QIN(:,:,:,:), QFD(:,:,:,:)
      real, allocatable :: VAR3D1(:,:,:), VAR3D2(:,:,:)

      integer, parameter :: NFDMAX=50 ! Max number of fields with the same HTFDCTL
      integer :: IDS(NFDMAX) ! All field IDs with the same HTFDCTL
      integer :: nFDS ! How many fields with the same HTFDCTL in the control file
      integer :: iID ! which field with HTFDCTL
      integer :: N1, N2
!     
!******************************************************************************
!
!     START MDL2STD_P. 
!

!     --------------WAFS block----------------------
!     450 ICIP
!     480 ICSEV
!     464 EDPARM
!     465 CAT
!     466 MWTURB
!     518 HGT
!     519 TMP
!     520 UGRD
!     521 VGRD
!     522 SPFH
!     523 RH
!     524 VVEL
!     525 ABSV
!     526 CLWMR
!     527 ICMR
!     528 RWMR
!     529 SNMR
!     530 GRLE
      IF(IGET(450)>0 .or. IGET(480)>0 .or. &
         IGET(464)>0 .or. IGET(465)>0 .or. IGET(466)>0 .or. &
         IGET(518)>0 .or. IGET(519)>0 .or. IGET(520)>0 .or. &
         IGET(521)>0 .or. IGET(522)>0 .or. IGET(523)>0 .or. &
         IGET(524)>0 .or. IGET(525)>0 .or. IGET(526)>0 .or. &
         IGET(527)>0 .or. IGET(528)>0 .or. IGET(529)>0 .or. &
         IGET(530)>0) then

!        STEP 1 -- U V (POSSIBLE FOR ABSV) INTERPLOCATION
         IF(IGET(520)>0 .or. IGET(521)>0 .or. IGET(525) > 0 ) THEN
!           U/V are always paired, use any for HTFDCTL          
            iID=520
            N = IAVBLFLD(IGET(iID))
            NFDCTL=size(pset%param(N)%level)
            if(allocated(ITYPEFDLVLCTL)) deallocate(ITYPEFDLVLCTL)
            allocate(ITYPEFDLVLCTL(NFDCTL))
            DO IFD = 1,NFDCTL
               ITYPEFDLVLCTL(IFD)=LVLS(IFD,IGET(iID))
            ENDDO
            if(allocated(HTFDCTL)) deallocate(HTFDCTL)
            allocate(HTFDCTL(NFDCTL))
            HTFDCTL=pset%param(N)%level
            DO i = 1, NFDCTL
               HTFDCTL(i)=P2H(HTFDCTL(i)/100.)
            ENDDO
            allocate(VAR3D1(IM,JSTA_2L:JEND_2U,NFDCTL))
            allocate(VAR3D2(IM,JSTA_2L:JEND_2U,NFDCTL))
            VAR3D1=SPVAL
            VAR3D2=SPVAL

            call FDLVL_UV(ITYPEFDLVLCTL,NFDCTL,HTFDCTL,VAR3D1,VAR3D2)

            DO IFD = 1,NFDCTL
               ! U
               IF (LVLS(IFD,IGET(520)) > 0) THEN
!$omp parallel do private(i,j)
                  DO J=JSTA,JEND
                  DO I=1,IM
                     GRID1(I,J)=VAR3D1(I,J,IFD)
                  ENDDO
                  ENDDO
                  if(grib=='grib2') then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(IGET(520))
                     fld_info(cfld)%lvl=LVLSXML(IFD,IGET(520))
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
               IF (LVLS(IFD,IGET(521)) > 0) THEN
!$omp parallel do private(i,j)
                  DO J=JSTA,JEND
                  DO I=1,IM
                     GRID1(I,J)=VAR3D2(I,J,IFD)
                  ENDDO
                  ENDDO
                  if(grib=='grib2') then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(IGET(521))
                     fld_info(cfld)%lvl=LVLSXML(IFD,IGET(521))
!$omp parallel do private(i,j,jj)
                     do j=1,jend-jsta+1
                        jj = jsta+j-1
                        do i=1,im
                           datapd(i,j,cfld) = GRID1(i,jj)
                        enddo
                     enddo
                  endif
               ENDIF
               ! ABSV
               IF (LVLS(IFD,IGET(525)) > 0) THEN
                  EGRID1=VAR3D1(1:IM,JSTA_2L:JEND_2U,IFD)
                  EGRID2=VAR3D2(1:IM,JSTA_2L:JEND_2U,IFD)
                  call CALVOR(EGRID1,EGRID2,EGRID3)
!$omp parallel do private(i,j)
                  DO J=JSTA,JEND
                  DO I=1,IM
                     GRID1(I,J)=EGRID3(I,J)
                  ENDDO
                  ENDDO
                  if(grib=='grib2') then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(IGET(525))
                     fld_info(cfld)%lvl=LVLSXML(IFD,IGET(525))
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

            deallocate(VAR3D1)
            deallocate(VAR3D2)

         ENDIF

!        STEP 2 -- MASS FIELDS INTERPOLATION EXCEPT:
!                  HGT(TO BE FIXED VALUES)
!                  RH ABSV (TO BE CACULATED)

         ALLOCATE(QIN(IM,JSTA:JEND,LM,NFDMAX))

!        INITIALIZE INPUTS
         nFDS = 0
         IF(IGET(450) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 450
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=icing_gfip(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(480) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 480
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=icing_gfis(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(464) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 464
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=gtg(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(465) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 465
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=catedr(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(466) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 466
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=mwt(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(519) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 519
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=T(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(522) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 522
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=Q(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(524) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 524
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=OMGA(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(526) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 526
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=QQW(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(527) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 527
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=QQR(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(528) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 528
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=QQS(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(529) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 529
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=QQG(1:IM,JSTA:JEND,1:LM)
         end if
         IF(IGET(530) > 0) THEN
            nFDS = nFDS + 1
            IDS(nFDS) = 530
            QIN(1:IM,JSTA:JEND,1:LM,nFDS)=QQI(1:IM,JSTA:JEND,1:LM)
         end if

!        FOR WAFS, ALL LEVLES OF DIFFERENT VARIABLES ARE THE SAME, USE ANY
         iID=IDS(1)
         N = IAVBLFLD(IGET(iID))
         NFDCTL=size(pset%param(N)%level)
         if(allocated(ITYPEFDLVLCTL)) deallocate(ITYPEFDLVLCTL)
         allocate(ITYPEFDLVLCTL(NFDCTL))
         DO IFD = 1,NFDCTL
            ITYPEFDLVLCTL(IFD)=LVLS(IFD,IGET(iID))
         ENDDO
         if(allocated(HTFDCTL)) deallocate(HTFDCTL)
         allocate(HTFDCTL(NFDCTL))
         HTFDCTL=pset%param(N)%level
         DO i = 1, NFDCTL
            HTFDCTL(i)=P2H(HTFDCTL(i)/100.)
         ENDDO

         ALLOCATE(QFD(IM,JSTA:JEND,NFDCTL,nFDS))
         QFD=SPVAL

         call FDLVL_MASS(ITYPEFDLVLCTL,NFDCTL,HTFDCTL,nFDS,QIN,QFD)

!        Adjust values before output
         N1 = -1
         DO N=1,nFDS
            iID=IDS(N)

!           Icing Potential
            if(iID==450) then
               N1=N
               DO IFD = 1,NFDCTL
                  DO J=JSTA,JEND
                  DO I=1,IM
                     if(QFD(I,J,IFD,N) < SPVAL) then
                        QFD(I,J,IFD,N)=max(0.0,QFD(I,J,IFD,N))
                        QFD(I,J,IFD,N)=min(1.0,QFD(I,J,IFD,N))
                     endif
                  ENDDO
                  ENDDO
               ENDDO
            endif

!           Icing severity categories
!              0 = none (0, 0.08)
!              1 = trace [0.08, 0.21]
!              2 = light (0.21, 0.37]
!              3 = moderate (0.37, 0.67]
!              4 = severe (0.67, 1]
!              http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-207.shtml
            if(iID==480) then
               DO IFD = 1,NFDCTL
                  DO J=JSTA,JEND
                  DO I=1,IM
                     if(N1 > 0) then
                        ! Icing severity is 0 when icing potential is too small
                        if(QFD(I,J,IFD,N1) < 0.001) QFD(I,J,IFD,N)=0.
                     endif
                     if(QFD(I,J,IFD,N) == SPVAL) cycle
                     if (QFD(I,J,IFD,N) < 0.08) then
                        QFD(I,J,IFD,N) = 0.0
                     elseif (QFD(I,J,IFD,N) <= 0.21) then
                        QFD(I,J,IFD,N) = 1.
                     else if(QFD(I,J,IFD,N) <= 0.37) then
                        QFD(I,J,IFD,N) = 2.0
                     else if(QFD(I,J,IFD,N) <= 0.67) then
                        QFD(I,J,IFD,N) = 3.0
                     else
                        QFD(I,J,IFD,N) = 4.0
                     endif
                  ENDDO
                  ENDDO
               ENDDO
            endif

!           GTG turbulence:  EDRPARM, CAT, MWTURB
            if(iID==464 .or. iID==465 .or. iID==466) then
               DO IFD = 1,NFDCTL
                  DO J=JSTA,JEND
                  DO I=1,IM
                     if(QFD(I,J,IFD,N) < SPVAL) then
                        QFD(I,J,IFD,N)=max(0.0,QFD(I,J,IFD,N))
                        QFD(I,J,IFD,N)=min(1.0,QFD(I,J,IFD,N))
                     endif
                  ENDDO
                  ENDDO
               ENDDO
            endif

!           SPFH
            if(iID==522) then
               DO IFD = 1,NFDCTL
                  DO J=JSTA,JEND
                  DO I=1,IM
                     if(QFD(I,J,IFD,N) < SPVAL) then
                        if(QFD(I,J,IFD,N)<1.0e-8) QFD(I,J,IFD,N)=0.
                     endif
                  ENDDO
                  ENDDO
               ENDDO
            endif

!           CLWMR ICMR RWMR SNMR GRLE
            if(iID==526 .or. iID==527 .or. iID==528 .or. &
               iID==529 .or. iID==530) then
               DO IFD = 1,NFDCTL
                  DO J=JSTA,JEND
                  DO I=1,IM
                     QFD(I,J,IFD,N)=max(0.0,QFD(I,J,IFD,N))
                  ENDDO
                  ENDDO
               ENDDO
            endif

         ENDDO

!        Output
         DO N=1,nFDS
            iID=IDS(N)
            DO IFD = 1,NFDCTL
               IF (LVLS(IFD,IGET(iID)) > 0) THEN
!$omp parallel do private(i,j)
                  DO J=JSTA,JEND
                  DO I=1,IM
                     GRID1(I,J)=QFD(I,J,IFD,N)
                  ENDDO
                  ENDDO
                  if(grib=='grib2') then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(IGET(iID))
                     fld_info(cfld)%lvl=LVLSXML(IFD,IGET(iID))
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
         ENDDO

         DEALLOCATE(QIN,QFD)

!        STEP 3 -  MASS FIELDS CALCULATION
!                  HGT(TO BE FIXED VALUES)
!                  RH ABSV (TO BE CACULATED)
         ! HGT
         IF(IGET(518) > 0) THEN
            iID=518
            N = IAVBLFLD(IGET(iID))
            NFDCTL=size(pset%param(N)%level)
            if(allocated(HTFDCTL)) deallocate(HTFDCTL)
            allocate(HTFDCTL(NFDCTL))
            HTFDCTL=pset%param(N)%level
            DO i = 1, NFDCTL
               HTFDCTL(i)=P2H(HTFDCTL(i)/100.)
            ENDDO

            DO IFD = 1,NFDCTL
               IF (LVLS(IFD,IGET(iID)) > 0) THEN
!$omp parallel do private(i,j)
                  DO J=JSTA,JEND
                  DO I=1,IM
                     GRID1(I,J)=HTFDCTL(IFD)
                  ENDDO
                  ENDDO
                  if(grib=='grib2') then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(IGET(iID))
                     fld_info(cfld)%lvl=LVLSXML(IFD,IGET(iID))
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
         ENDIF

         ! RH
         IF(IGET(523) > 0) THEN
            iID=523
            N = IAVBLFLD(IGET(iID))
            NFDCTL=size(pset%param(N)%level)
            if(allocated(ITYPEFDLVLCTL)) deallocate(ITYPEFDLVLCTL)
            allocate(ITYPEFDLVLCTL(NFDCTL))
            DO IFD = 1,NFDCTL
               ITYPEFDLVLCTL(IFD)=LVLS(IFD,IGET(iID))
            ENDDO
            if(allocated(HTFDCTL)) deallocate(HTFDCTL)
            allocate(HTFDCTL(NFDCTL))
            HTFDCTL=pset%param(N)%level
            DO i = 1, NFDCTL
               HTFDCTL(i)=P2H(HTFDCTL(i)/100.)
            ENDDO

            ALLOCATE(QIN(IM,JSTA:JEND,LM,2))
            QIN(1:IM,JSTA:JEND,1:LM,1)=T(1:IM,JSTA:JEND,1:LM)
            QIN(1:IM,JSTA:JEND,1:LM,2)=Q(1:IM,JSTA:JEND,1:LM)

            ALLOCATE(QFD(IM,JSTA:JEND,NFDCTL,2))
            QFD=SPVAL

            call FDLVL_MASS(ITYPEFDLVLCTL,NFDCTL,HTFDCTL,2,QIN,QFD)
            HTFDCTL=pset%param(N)%level ! Save back to pressure


            DO IFD = 1,NFDCTL
               IF (LVLS(IFD,IGET(iID)) > 0) THEN
!$omp parallel do private(i,j)
                  DO J=JSTA,JEND
                  DO I=1,IM
                     EGRID2(I,J) = HTFDCTL(IFD)                    ! P
                  ENDDO
                  ENDDO

                  EGRID3(1:IM,JSTA:JEND)=QFD(1:IM,JSTA:JEND,IFD,1) ! T
                  EGRID4(1:IM,JSTA:JEND)=QFD(1:IM,JSTA:JEND,IFD,2) ! Q
                  EGRID1 = SPVAL

                  IF(MODELNAME == 'GFS' .or. MODELNAME == 'FV3R')THEN
                     CALL CALRH_GFS(EGRID2(1,jsta),EGRID3(1,jsta),&
                       EGRID4(1,jsta), EGRID1(1,jsta))
                  ELSEIF (MODELNAME == 'RAPR')THEN 
                     CALL CALRH_GSD(EGRID2(1,jsta),EGRID3(1,jsta),&
                       EGRID4(1,jsta), EGRID1(1,jsta))
                  ELSE
                     CALL CALRH(EGRID2(1,jsta),EGRID3(1,jsta),&
                       EGRID4(1,jsta), EGRID1(1,jsta))
                  END IF

!$omp  parallel do private(i,j)
                  DO J=JSTA,JEND
                  DO I=1,IM
                     IF(EGRID1(I,J) < SPVAL) THEN
                        GRID1(I,J) = EGRID1(I,J)*100.
                     ELSE
                        GRID1(I,J)  = EGRID1(I,J)
                     ENDIF
                  ENDDO
                  ENDDO

                  if(grib=='grib2') then
                     cfld=cfld+1
                     fld_info(cfld)%ifld=IAVBLFLD(IGET(iID))
                     fld_info(cfld)%lvl=LVLSXML(IFD,IGET(iID))
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
            deallocate(QIN,QFD)

         ENDIF

      ENDIF
!
!     END OF ROUTINE.
!
      RETURN
      END

      FUNCTION P2H(p)
      implicit none
      real, intent(in) :: p
      real :: P2H
!     To convert pressure levels (hPa) to geopotantial heights
!     Uses ICAO standard atmosphere parameters as defined here:
!        https://www.nen.nl/pdfpreview/preview_29424.pdf
      real, parameter :: lapse = 0.0065
      real, parameter :: surf_temp = 288.15
      real, parameter :: gravity = 9.80665
      real, parameter :: moles_dry_air = 0.02896442
      real, parameter :: gas_const = 8.31432
      real, parameter :: surf_pres = 1013.25
      real, parameter :: power_const = (gravity * moles_dry_air) &
                                       / (gas_const * lapse)
      P2H = (surf_temp/lapse)*(1-(p/surf_pres)**(1/power_const))
      END
