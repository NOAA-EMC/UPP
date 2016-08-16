      SUBROUTINE CALVOR(UWND,VWND,ABSV)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALVOR      COMPUTES ABSOLUTE VORTICITY
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-22       
!     
! ABSTRACT:  
!     THIS ROUTINE COMPUTES THE ABSOLUTE VORTICITY.
!   .     
!     
! PROGRAM HISTORY LOG:
!   92-12-22  RUSS TREADON
!   98-06-08  T BLACK - CONVERSION FROM 1-D TO 2-D
!   00-01-04  JIM TUCCILLO - MPI VERSION
!   02-01-15  MIKE BALDWIN - WRF VERSION C-GRID
!   05-03-01  H CHUANG - ADD NMM E GRID
!   05-05-17  H CHUANG - ADD POTENTIAL VORTICITY CALCULATION
!   05-07-07  B ZHOU   - ADD RSM IN COMPUTING DVDX, DUDY AND UAVG
!   13-08-09  S MOORTHI - Optimize the vorticity loop including threading

!     
! USAGE:    CALL CALVOR(UWND,VWND,ABSV)
!   INPUT ARGUMENT LIST:
!     UWND     - U WIND (M/S) MASS-POINTS
!     VWND     - V WIND (M/S) MASS-POINTS
!
!   OUTPUT ARGUMENT LIST: 
!     ABSV     - ABSOLUTE VORTICITY (1/S) MASS-POINTS
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
      use vrbls2d,      only: f
      use masks,        only: gdlat, gdlon, dx, dy
      use params_mod,   only: d00, dtr, small, erad
      use ctlblk_mod,   only: jsta_2l, jend_2u, spval, modelname, global, &
                              jsta, jend, im, jm, jsta_m, jend_m, gdsdegr
      use gridspec_mod, only: gridtype, dyval

      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL, dimension(im,jsta_2l:jend_2u), intent(in)    :: UWND, VWND
      REAL, dimension(im,jsta_2l:jend_2u), intent(inout) :: ABSV
!
      real,    allocatable ::  wrk1(:,:), wrk2(:,:), wrk3(:,:), cosl(:,:)
      INTEGER, allocatable ::  IHE(:),IHW(:), IE(:),IW(:)
!
      integer I,J,ip1,im1,ii,iir,iil,jj,JMT2,imb2
      real    R2DX,R2DY,DVDX,DUDY,UAVG,TPH1,TPHI
!     
!***************************************************************************
!     START CALVOR HERE.
!     
!     LOOP TO COMPUTE ABSOLUTE VORTICITY FROM WINDS.
!     
      IF(MODELNAME  == 'RAPR') then
!$omp  parallel do private(i,j)
        DO J=JSTA_2L,JEND_2U
          DO I=1,IM
            ABSV(I,J) = D00
          ENDDO
        ENDDO
      else
!$omp  parallel do private(i,j)
        DO J=JSTA_2L,JEND_2U
          DO I=1,IM
            ABSV(I,J) = SPVAL
          ENDDO
        ENDDO
      endif

!      print*,'dyval in CALVOR= ',DYVAL 
  
      CALL EXCH_F(UWND)
!
      IF (MODELNAME == 'GFS' .or. global) THEN
        CALL EXCH(GDLAT(1,JSTA_2L))

        allocate (wrk1(im,jsta:jend), wrk2(im,jsta:jend),          &
     &            wrk3(im,jsta:jend), cosl(im,jsta_2l:jend_2u))
        allocate(iw(im),ie(im))

        imb2 = im/2
!$omp  parallel do private(i)
      do i=1,im
        ie(i) = i+1
        iw(i) = i-1
      enddo
      iw(1)  = im
      ie(im) = 1

!       if(1>=jsta .and. 1<=jend)then
!        if(cos(gdlat(1,1)*dtr)<small)poleflag=.T.
!       end if 	
!       call mpi_bcast(poleflag,1,MPI_LOGICAL,0,mpi_comm_comp,iret)

!$omp  parallel do private(i,j,ip1,im1)
        DO J=JSTA,JEND
          do i=1,im
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
       
!$omp  parallel do private(i,j,ii)
        DO J=JSTA,JEND
          if (j == 1) then
            if(gdlat(1,j) > 0.) then ! count from north to south
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J+1)-GDLAT(II,J))*DTR) !1/dphi
              enddo
            else ! count from south to north
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J+1)+GDLAT(II,J))*DTR) !1/dphi
!
              enddo
            end if      
          elseif (j == JM) then
            if(gdlat(1,j) < 0.) then ! count from north to south
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J-1)+GDLAT(II,J))*DTR)
              enddo
            else ! count from south to north
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J-1)-GDLAT(II,J))*DTR)
              enddo
            end if  
          else
            do i=1,im
              wrk3(i,j) = 1.0 / ((GDLAT(I,J-1)-GDLAT(I,J+1))*DTR) !1/dphi
            enddo
          endif
        enddo  

!$omp  parallel do private(i,j,ip1,im1,ii,jj)
        DO J=JSTA,JEND
          IF(J == 1) then                            ! Near North or South pole
            DO I=1,IM
              ip1 = ie(i)
              im1 = iw(i)
              IF(cosl(i,j) >= SMALL) THEN            !not a pole point
                ii = i + imb2
                if (ii > im) ii = ii - im
                ABSV(I,J) = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)               &
     &                    +  (UWND(II,J)*COSL(II,J)                            &
     &                    +   UWND(I,J+1)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)  &
     &                    + F(I,J)
              ELSE                                   !pole point, compute at j=2
                jj = 2
                ABSV(I,J) = ((VWND(ip1,JJ)-VWND(im1,JJ))*wrk2(i,jj)                   &
     &                    +  (UWND(I,J)*COSL(I,J)                                     &
                          +   UWND(I,jj+1)*COSL(I,Jj+1))*abs(wrk3(i,jj))) * wrk1(i,jj)&
     &                    + F(I,Jj)
              END IF
            ENDDO
          ELSE IF(J == JM) THEN                      ! Near North or South Pole
            DO I=1,IM
              ip1 = ie(i)
              im1 = iw(i)
              IF(cosl(i,j) >= SMALL) THEN            !not a pole point
                ii = i + imb2
                if (ii > im) ii = ii - im
                ABSV(I,J) = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)              &
     &                    +  (UWND(I,J-1)*COSL(I,J-1)                         &
     &                    +   UWND(II,J)*COSL(II,J))*wrk3(i,j)) * wrk1(i,j)   &
     &                    + F(I,J)
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                ABSV(I,J) = ((VWND(ip1,JJ)-VWND(im1,JJ))*wrk2(i,jj)             &
     &                    +  (UWND(I,jj-1)*COSL(I,Jj-1)                         &
     &                    +   UWND(I,J)*COSL(I,J))*abs(wrk3(i,jj))) * wrk1(i,jj)&
     &                    + F(I,Jj)
              END IF
            ENDDO
          ELSE
            DO I=1,IM
              ip1 = ie(i)
              im1 = iw(i)
              ABSV(I,J)   = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)               &
     &                    -  (UWND(I,J-1)*COSL(I,J-1)                          &
                          -   UWND(I,J+1)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)  &
                          + F(I,J)
            ENDDO
          END IF
!          if(ABSV(I,J)>1.0)print*,'Debug CALVOR',i,j,VWND(ip1,J),VWND(im1,J), &
!          wrk2(i,j),UWND(I,J-1),COSL(I,J-1),UWND(I,J+1),COSL(I,J+1),wrk3(i,j),cosl(i,j),F(I,J),ABSV(I,J)
        END DO                               ! end of J loop

!       deallocate (wrk1, wrk2, wrk3, cosl)
! GFS use lon avg as one scaler value for pole point

        call poleavg(IM,JM,JSTA,JEND,SMALL,COSL(1,jsta),SPVAL,ABSV(1,jsta))
        deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)
     
      ELSE IF(GRIDTYPE == 'A')THEN
!$omp parallel do  private(i,j,jmt2,tphi,r2dx,r2dy,dvdx,dudy,uavg)
        DO J=JSTA_M,JEND_M
          JMT2 = JM/2+1
          TPHI = (J-JMT2)*(DYVAL/gdsdegr)*DTR
          DO I=2,IM-1
            IF(VWND(I+1,J).LT.SPVAL.AND.VWND(I-1,J).LT.SPVAL.AND.              &
     &         UWND(I,J+1).LT.SPVAL.AND.UWND(I,J-1).LT.SPVAL) THEN
              R2DX   = 1./(2.*DX(I,J))
              R2DY   = 1./(2.*DY(I,J))
              DVDX   = (VWND(I+1,J)-VWND(I-1,J))*R2DX
              DUDY   = (UWND(I,J+1)-UWND(I,J-1))*R2DY
              UAVG   = 0.25*(UWND(I+1,J)+UWND(I-1,J)                           &
     &               +       UWND(I,J+1)+UWND(I,J-1))
!  is there a (f+tan(phi)/erad)*u term?
              ABSV(I,J) = DVDX - DUDY + F(I,J) + UAVG*TAN(GDLAT(I,J)*DTR)/ERAD  ! not sure about this???
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
          DO I=2,IM-1
            IF(VWND(I+IHE(J),J) < SPVAL.AND.VWND(I+IHW(J),J) < SPVAL .AND.   &
     &         UWND(I,J+1) < SPVAL     .AND.UWND(I,J-1) < SPVAL) THEN
              R2DX   = 1./(2.*DX(I,J))
              R2DY   = 1./(2.*DY(I,J))
              DVDX   = (VWND(I+IHE(J),J)-VWND(I+IHW(J),J))*R2DX
              DUDY   = (UWND(I,J+1)-UWND(I,J-1))*R2DY
              UAVG   = 0.25*(UWND(I+IHE(J),J)+UWND(I+IHW(J),J)                 &
     &               +       UWND(I,J+1)+UWND(I,J-1))
!  is there a (f+tan(phi)/erad)*u term?
              ABSV(I,J) = DVDX - DUDY + F(I,J) + UAVG*TAN(TPHI)/ERAD 
            END IF
          END DO
        END DO
       deallocate(ihw, IHE)
      ELSE IF (GRIDTYPE == 'B')THEN
        CALL EXCH_F(VWND)
        DO J=JSTA_M,JEND_M
          JMT2 = JM/2+1
          TPHI = (J-JMT2)*(DYVAL/gdsdegr)*DTR
          DO I=2,IM-1         
            R2DX = 1./DX(I,J)
            R2DY = 1./DY(I,J)
            DVDX = (0.5*(VWND(I,J)+VWND(I,J-1))-0.5*(VWND(I-1,J)               &
     &           +       VWND(I-1,J-1)))*R2DX
            DUDY = (0.5*(UWND(I,J)+UWND(I-1,J))-0.5*(UWND(I,J-1)               &
     &           +       UWND(I-1,J-1)))*R2DY
            UAVG = 0.25*(UWND(I-1,J-1)+UWND(I-1,J)                             &
     &           +       UWND(I,  J-1)+UWND(I,  J))
!  is there a (f+tan(phi)/erad)*u term?
           ABSV(I,J) = DVDX - DUDY + F(I,J) + UAVG*TAN(TPHI)/ERAD 
          END DO
        END DO 
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
! PROGRAM HISTORY LOG:
!   16-05-05  SAJAL KAR MODIFIED CALVORT TO COMPUTE DIVERGENCE FROM
!             WIND COMPONENTS
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
                              jsta, jend, im, jm, jsta_m, jend_m
      use gridspec_mod, only: gridtype

      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL, dimension(im,jsta_2l:jend_2u), intent(in)    :: UWND,VWND
      REAL, dimension(im,jsta_2l:jend_2u), intent(inout) :: DIV
!
      real,    allocatable ::  wrk1(:,:), wrk2(:,:), wrk3(:,:), cosl(:,:)
      INTEGER, allocatable ::  IHE(:),IHW(:), IE(:),IW(:)
!
      integer I,J,ip1,im1,ii,iir,iil,jj,imb2
!     
!***************************************************************************
!     START CALDIV HERE.
!     
!     LOOP TO COMPUTE DIVERGENCE FROM WINDS.
!     
!$omp  parallel do private(i,j)
!sk06162016    DO J=JSTA_2L,JEND_2U
      DO J=JSTA,JEND
        DO I=1,IM
          DIV(I,J) = SPVAL
!sk       DIV(I,J) = D00
        ENDDO
      ENDDO

      CALL EXCH_F(VWND)

!     IF (MODELNAME == 'GFS' .or. global) THEN
        CALL EXCH(GDLAT(1,JSTA_2L))

        allocate (wrk1(im,jsta:jend), wrk2(im,jsta:jend),          &
     &            wrk3(im,jsta:jend), cosl(im,jsta_2l:jend_2u))
        allocate(iw(im),ie(im))

        imb2 = im/2
!$omp  parallel do private(i)
        do i=1,im
          ie(i) = i+1
          iw(i) = i-1
        enddo
        iw(1)  = im
        ie(im) = 1


!$omp  parallel do private(i,j,ip1,im1)
        DO J=JSTA,JEND
          do i=1,im
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
       
!$omp  parallel do private(i,j,ii)
        DO J=JSTA,JEND
          if (j == 1) then
            if(gdlat(1,j) > 0.) then ! count from north to south
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J+1)-GDLAT(II,J))*DTR) !1/dphi
              enddo
            else ! count from south to north
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J+1)+GDLAT(II,J))*DTR) !1/dphi
              enddo
            end if      
          elseif (j == JM) then
            if(gdlat(1,j) < 0.) then ! count from north to south
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J-1)+GDLAT(II,J))*DTR)
              enddo
            else ! count from south to north
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J-1)-GDLAT(II,J))*DTR)
              enddo
            end if  
          else
            do i=1,im
              wrk3(i,j) = 1.0 / ((GDLAT(I,J-1)-GDLAT(I,J+1))*DTR) !1/dphi
            enddo
          endif
        enddo  

!$omp  parallel do private(i,j,ip1,im1,ii,jj)
        DO J=JSTA,JEND
          IF(J == 1) then                            ! Near North pole
            DO I=1,IM
              ip1 = ie(i)
              im1 = iw(i)
              IF(cosl(i,j) >= SMALL) THEN            !not a pole point
                ii = i + imb2
                if (ii > im) ii = ii - im
                DIV(I,J) = ((UWND(ip1,J)-UWND(im1,J))*wrk2(i,j)               &
     &                   -  (VWND(II,J)*COSL(II,J)                            &
     &                   +   VWND(I,J+1)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)
!sk06142016a
                if(DIV(I,J)>1.0)print*,'Debug in CALDIV',i,j,UWND(ip1,J),UWND(im1,J), &
     &             wrk2(i,j),VWND(II,J),COSL(II,J),VWND(I,J+1),COSL(I,J+1),        &
     &             wrk3(i,j),DIV(I,J)
!--
              ELSE                             !North pole point, compute at j=2
                jj = 2
                DIV(I,J) = ((UWND(ip1,jj)-UWND(im1,jj))*wrk2(i,jj)                   &
     &                   -  (VWND(I,J)*COSL(I,J)                                     &
                         +   VWND(I,jj+1)*COSL(I,jj+1))*abs(wrk3(i,jj))) * wrk1(i,jj)
!sk06142016b
                if(DIV(I,J)>1.0)print*,'Debug in CALDIV',i,j,UWND(ip1,jj),UWND(im1,jj), &
     &             wrk2(i,jj),VWND(I,J),COSL(I,J),VWND(I,jj+1),COSL(I,jj+1),         &
     &             wrk3(i,jj),wrk1(i,jj),DIV(I,J)
!--
              END IF
            ENDDO
          ELSE IF(J == JM) THEN                      ! Near South pole
            DO I=1,IM
              ip1 = ie(i)
              im1 = iw(i)
              IF(cosl(i,j) >= SMALL) THEN            !not a pole point
                ii = i + imb2
                if (ii > im) ii = ii - im
                DIV(I,J) = ((UWND(ip1,J)-UWND(im1,J))*wrk2(i,j)              &
     &                   -  (VWND(I,J-1)*COSL(I,J-1)                         &
     &                   +   VWND(II,J)*COSL(II,J))*wrk3(i,j)) * wrk1(i,j)
!sk06142016c
                if(DIV(I,J)>1.0)print*,'Debug in CALDIV',i,j,UWND(ip1,J),UWND(im1,J), &
     &             wrk2(i,j),VWND(I,J-1),COSL(I,J-1),VWND(II,J),COSL(II,J),        &
     &             wrk3(i,j),wrk1(i,j),DIV(I,J)
!--
              ELSE                              !South pole point,compute at jm-1
                jj = jm-1
                DIV(I,J) = ((UWND(ip1,JJ)-UWND(im1,JJ))*wrk2(i,jj)             &
     &                   -  (VWND(I,jj-1)*COSL(I,Jj-1)                         &
     &                   +   VWND(I,J)*COSL(I,J))*abs(wrk3(i,jj))) * wrk1(i,jj)
!sk06142016d
                if(DIV(I,J)>1.0)print*,'Debug in CALDIV',i,j,UWND(ip1,JJ),UWND(im1,JJ), &
     &             wrk2(i,jj),VWND(I,jj-1),COSL(I,jj-1),VWND(I,J),COSL(I,J),         &
     &             wrk3(i,jj),wrk1(i,jj),DIV(I,J)
!--
              END IF
            ENDDO
          ELSE
            DO I=1,IM
              ip1 = ie(i)
              im1 = iw(i)
              DIV(I,J)   = ((UWND(ip1,J)-UWND(im1,J))*wrk2(i,j)               &
     &                   +  (VWND(I,J-1)*COSL(I,J-1)                          &
                         -   VWND(I,J+1)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)
!sk06132016
              if(DIV(I,J)>1.0)print*,'Debug in CALDIV',i,j,UWND(ip1,J),UWND(im1,J), &
     &           wrk2(i,j),VWND(I,J-1),COSL(I,J-1),VWND(I,J+1),COSL(I,J+1),         &
     &           wrk3(i,j),wrk1(i,j),DIV(I,J)
!--
            ENDDO
          END IF
        ENDDO                               ! end of J loop

! GFS use lon avg as one scaler value for pole point
        call poleavg(IM,JM,JSTA,JEND,SMALL,COSL(1,jsta),SPVAL,DIV(1,jsta))
!sk06142016e
              if(DIV(1,jsta)>1.0)print*,'Debug in CALDIV',jsta,DIV(1,jsta)
!             print*,'Debug in CALDIV',' jsta= ',jsta,DIV(1,jsta)
!--
        deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)
     
!     END IF 

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
                              jsta, jend, im, jm, jsta_m, jend_m
      use gridspec_mod, only: gridtype

      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL, dimension(im,jsta_2l:jend_2u), intent(in)    :: PS
      REAL, dimension(im,jsta_2l:jend_2u), intent(inout) :: PSX,PSY 
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
!$omp  parallel do private(i,j)
!sk06162016   DO J=JSTA_2L,JEND_2U
      DO J=JSTA,JEND
        DO I=1,IM
          PSX(I,J) = SPVAL
          PSY(I,J) = SPVAL
!sk       PSX(I,J) = D00
!sk       PSY(I,J) = D00
        ENDDO
      ENDDO

      CALL EXCH_F(PS)

!     IF (MODELNAME == 'GFS' .or. global) THEN
        CALL EXCH(GDLAT(1,JSTA_2L))

        allocate (wrk1(im,jsta:jend), wrk2(im,jsta:jend),          &
     &            wrk3(im,jsta:jend), cosl(im,jsta_2l:jend_2u))
        allocate(iw(im),ie(im))

        imb2 = im/2
!$omp  parallel do private(i)
        do i=1,im
          ie(i) = i+1
          iw(i) = i-1
        enddo
        iw(1)  = im
        ie(im) = 1


!$omp  parallel do private(i,j,ip1,im1)
        DO J=JSTA,JEND
          do i=1,im
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
            if(gdlat(1,j) > 0.) then ! count from north to south
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J+1)-GDLAT(II,J))*DTR) !1/dphi
              enddo
            else ! count from south to north
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J+1)+GDLAT(II,J))*DTR) !1/dphi
              enddo
            end if      
          elseif (j == JM) then
            if(gdlat(1,j) < 0.) then ! count from north to south
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J-1)+GDLAT(II,J))*DTR)
              enddo
            else ! count from south to north
              do i=1,im
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J-1)-GDLAT(II,J))*DTR)
              enddo
            end if  
          else
            do i=1,im
              wrk3(i,j) = 1.0 / ((GDLAT(I,J-1)-GDLAT(I,J+1))*DTR) !1/dphi
            enddo
          endif
        ENDDO  

!$omp  parallel do private(i,j,ip1,im1,ii,jj)
        DO J=JSTA,JEND
          IF(J == 1) then                            ! Near North pole
            DO I=1,IM
              ip1 = ie(i)
              im1 = iw(i)
              IF(cosl(i,j) >= SMALL) THEN            !not a pole point
                ii = i + imb2
                if (ii > im) ii = ii - im
                PSX(I,J) = (PS(ip1,J)-PS(im1,J))*wrk2(i,j)*wrk1(i,j)
                PSY(I,J) = (PS(II,J)-PS(I,J+1))*wrk3(i,j)/ERAD 
              ELSE                             !North pole point, compute at j=2
                jj = 2
                PSX(I,J) = (PS(ip1,jj)-PS(im1,jj))*wrk2(i,jj)*wrk1(i,jj)
                PSY(I,J) = (PS(I,J)-PS(I,jj+1))*abs(wrk3(i,jj))/ERAD
              END IF
            ENDDO
          ELSE IF(J == JM) THEN                      ! Near South pole
            DO I=1,IM
              ip1 = ie(i)
              im1 = iw(i)
              IF(cosl(i,j) >= SMALL) THEN            !not a pole point
                ii = i + imb2
                if (ii > im) ii = ii - im
                PSX(I,J) = (PS(ip1,J)-PS(im1,J))*wrk2(i,j)*wrk1(i,j)
                PSY(I,J) = (PS(I,J-1)-PS(II,J))*wrk3(i,j)/ERAD
              ELSE                              !South pole point,compute at jm-1
                jj = jm-1
                PSX(I,J) = (PS(ip1,JJ)-PS(im1,JJ))*wrk2(i,jj)*wrk1(i,jj)
                PSY(I,J) = (PS(I,jj-1)-PS(I,J))*abs(wrk3(i,jj))/ERAD
              END IF
            ENDDO
          ELSE
            DO I=1,IM
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
        ENDDO                               ! end of J loop

        deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)
     
!     END IF 

      END SUBROUTINE CALGRADPS
