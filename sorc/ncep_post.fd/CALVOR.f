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
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE CALVOR(UWND,VWND,ABSV)

!     
!
      use vrbls2d,      only: f
      use masks,        only: gdlat, gdlon, dx, dy
      use params_mod,   only: d00, dtr, small, erad
      use ctlblk_mod,   only: jsta_2l, jend_2u, spval, modelname, global, &
                              jsta, jend, im, jm, jsta_m, jend_m, gdsdegr
      use gridspec_mod, only: gridtype, dyval
      use upp_math,     only: DVDXDUDY, DDVDX, DDUDY, UUAVG

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

        npass = 0

        jtem = jm / 18 + 1
!$omp  parallel do private(i,j,ip1,im1,ii,jj,tx1,tx2)
        DO J=JSTA,JEND
!         npass = npass2
!         if (j > jm-jtem+1 .or. j < jtem) npass = npass3
          IF(J == 1) then                            ! Near North or South pole
            if(gdlat(1,j) > 0.) then ! count from north to south
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VWND(ip1,J)==SPVAL .or. VWND(im1,J)==SPVAL .or. &
                     UWND(II,J)==SPVAL .or. UWND(I,J+1)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)               &
     &                      +  (UWND(II,J)*COSL(II,J)                            &
     &                      +   UWND(I,J+1)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)  &
     &                      + F(I,J)
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=1,IM
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
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VWND(ip1,J)==SPVAL .or. VWND(im1,J)==SPVAL .or. &
                     UWND(II,J)==SPVAL .or. UWND(I,J+1)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)               &
     &                      -  (UWND(II,J)*COSL(II,J)                            &
     &                      +   UWND(I,J+1)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)  &
     &                      + F(I,J)
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=1,IM
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
            if(gdlat(1,j) < 0.) then ! count from north to south
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VWND(ip1,J)==SPVAL .or. VWND(im1,J)==SPVAL .or. &
                     UWND(I,J-1)==SPVAL .or. UWND(II,J)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)              &
     &                      -  (UWND(I,J-1)*COSL(I,J-1)                         &
     &                      +   UWND(II,J)*COSL(II,J))*wrk3(i,j)) * wrk1(i,j)   &
     &                      + F(I,J)
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=1,IM
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
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VWND(ip1,J)==SPVAL .or. VWND(im1,J)==SPVAL .or. &
                     UWND(I,J-1)==SPVAL .or. UWND(II,J)==SPVAL) cycle
                  ABSV(I,J) = ((VWND(ip1,J)-VWND(im1,J))*wrk2(i,j)              &
     &                      +  (UWND(I,J-1)*COSL(I,J-1)                         &
     &                      +   UWND(II,J)*COSL(II,J))*wrk3(i,j)) * wrk1(i,j)   &
     &                      + F(I,J)
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=1,IM
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
            DO I=1,IM
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
            do i=1,im
              tx1(i) = absv(i,j)
            enddo
            do nn=1,npass
              do i=1,im
                tx2(i+1) = tx1(i)
              enddo
              tx2(1)    = tx2(im+1)
              tx2(im+2) = tx2(2)
              do i=2,im+1
                tx1(i-1) = 0.25 * (tx2(i-1) + tx2(i+1)) + 0.5*tx2(i)
              enddo
            enddo
            do i=1,im
              absv(i,j) = tx1(i)
            enddo
          endif
        END DO                               ! end of J loop

!       deallocate (wrk1, wrk2, wrk3, cosl)
! GFS use lon avg as one scaler value for pole point

        call poleavg(IM,JM,JSTA,JEND,SMALL,COSL(1,jsta),SPVAL,ABSV(1,jsta))
        deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)

      ELSE !(MODELNAME == 'GFS' .or. global)

      IF (GRIDTYPE == 'B')THEN
        CALL EXCH_F(VWND)
      ENDIF
     
      CALL DVDXDUDY(UWND,VWND)

      IF(GRIDTYPE == 'A')THEN
!$omp parallel do  private(i,j,jmt2,tphi,r2dx,r2dy,dvdx,dudy,uavg)
        DO J=JSTA_M,JEND_M
          JMT2 = JM/2+1
          TPHI = (J-JMT2)*(DYVAL/gdsdegr)*DTR
          DO I=2,IM-1
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
          DO I=2,IM-1
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
!        CALL EXCH_F(VWND)      !done before dvdxdudy() Jesse 20200520
        DO J=JSTA_M,JEND_M
          JMT2 = JM/2+1
          TPHI = (J-JMT2)*(DYVAL/gdsdegr)*DTR
          DO I=2,IM-1         
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
      use vrbls2d,      only: f
      use masks,        only: gdlat, gdlon
      use params_mod,   only: d00, dtr, small, erad
      use ctlblk_mod,   only: jsta_2l, jend_2u, spval, modelname, global, &
                              jsta, jend, im, jm, jsta_m, jend_m, lm
      use gridspec_mod, only: gridtype

      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL, dimension(im,jsta_2l:jend_2u,lm), intent(in)    :: UWND,VWND
      REAL, dimension(im,jsta:jend,lm),       intent(inout) :: DIV
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
      CALL EXCH(GDLAT(1,JSTA_2L))

      allocate (wrk1(im,jsta:jend), wrk2(im,jsta:jend),          &
     &          wrk3(im,jsta:jend), cosl(im,jsta_2l:jend_2u))
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
      
      do l=1,lm
!$omp  parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
            DIV(I,J,l) = SPVAL
          ENDDO
        ENDDO

        CALL EXCH_F(VWND(1,jsta_2l,l))

!$omp  parallel do private(i,j,ip1,im1,ii,jj)
        DO J=JSTA,JEND
          IF(J == 1) then                          ! Near North pole
            if(gdlat(1,j) > 0.) then ! count from north to south
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  DIV(I,J,l) = ((UWND(ip1,J,l)-UWND(im1,J,l))*wrk2(i,j)           &
     &                       -  (VWND(II,J,l)*COSL(II,J)                          &
     &                       +   VWND(I,J+1,l)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)
                enddo
!--
              ELSE                             !North pole point, compute at j=2
                jj = 2
                do i=1,im
                  ip1 = ie(i)
                  im1 = iw(i)
                  DIV(I,J,l) = ((UWND(ip1,jj,l)-UWND(im1,jj,l))*wrk2(i,jj)         &
     &                       +  (VWND(I,J,l)*COSL(I,J)                             &
                             -   VWND(I,jj+1,l)*COSL(I,jj+1))*wrk3(i,jj)) * wrk1(i,jj)
                enddo
!--
              ENDIF
            else
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  DIV(I,J,l) = ((UWND(ip1,J,l)-UWND(im1,J,l))*wrk2(i,j)           &
     &                       +  (VWND(II,J,l)*COSL(II,J)                          &
     &                       +   VWND(I,J+1,l)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)
                enddo
!--
              ELSE                             !North pole point, compute at j=2
                jj = 2
                do i=1,im
                  ip1 = ie(i)
                  im1 = iw(i)
                  DIV(I,J,l) = ((UWND(ip1,jj,l)-UWND(im1,jj,l))*wrk2(i,jj)         &
     &                       -  (VWND(I,J,l)*COSL(I,J)                             &
                             -   VWND(I,jj+1,l)*COSL(I,jj+1))*wrk3(i,jj)) * wrk1(i,jj)
                enddo
              ENDIF
            endif
          ELSE IF(J == JM) THEN                    ! Near South pole
            if(gdlat(1,j) < 0.) then ! count from north to south
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  DIV(I,J,l) = ((UWND(ip1,J,l)-UWND(im1,J,l))*wrk2(i,j)          &
     &                       +  (VWND(I,J-1,l)*COSL(I,J-1)                       &
     &                       +   VWND(II,J,l)*COSL(II,J))*wrk3(i,j)) * wrk1(i,j)
                enddo
!--
              ELSE                              !South pole point,compute at jm-1
                jj = jm-1
                do i=1,im
                  ip1 = ie(i)
                  im1 = iw(i)
                  DIV(I,J,l) = ((UWND(ip1,JJ,l)-UWND(im1,JJ,l))*wrk2(i,jj)       &
     &                       +  (VWND(I,jj-1,l)*COSL(I,Jj-1)                     &
     &                       -   VWND(I,J,l)*COSL(I,J))*wrk3(i,jj)) * wrk1(i,jj)

                enddo
              ENDIF
            else
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  DIV(I,J,l) = ((UWND(ip1,J,l)-UWND(im1,J,l))*wrk2(i,j)          &
     &                       -  (VWND(I,J-1,l)*COSL(I,J-1)                       &
     &                       +   VWND(II,J,l)*COSL(II,J))*wrk3(i,j)) * wrk1(i,j)
                enddo
!--
              ELSE                              !South pole point,compute at jm-1
                jj = jm-1
                do i=1,im
                  ip1 = ie(i)
                  im1 = iw(i)
                  DIV(I,J,l) = ((UWND(ip1,JJ,l)-UWND(im1,JJ,l))*wrk2(i,jj)       &
     &                       -  (VWND(I,jj-1,l)*COSL(I,Jj-1)                     &
     &                       -   VWND(I,J,l)*COSL(I,J))*wrk3(i,jj)) * wrk1(i,jj)

                enddo
              ENDIF
            endif
          ELSE
            DO I=1,IM
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
        call poleavg(IM,JM,JSTA,JEND,SMALL,COSL(1,jsta),SPVAL,DIV(1,jsta,l))
!sk06142016e
        if(DIV(1,jsta,l)>1.0)print*,'Debug in CALDIV',jsta,DIV(1,jsta,l)
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
!sk06162016   DO J=JSTA_2L,JEND_2U
!$omp  parallel do private(i,j)
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
            if(gdlat(1,j) > 0.) then ! count from north to south
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSX(I,J) = (PS(ip1,J)-PS(im1,J))*wrk2(i,j)*wrk1(i,j)
                  PSY(I,J) = (PS(II,J)-PS(I,J+1))*wrk3(i,j)/ERAD 
                enddo
              ELSE                             !North pole point, compute at j=2
                jj = 2
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSX(I,J) = (PS(ip1,jj)-PS(im1,jj))*wrk2(i,jj)*wrk1(i,jj)
                  PSY(I,J) = (PS(I,J)-PS(I,jj+1))*wrk3(i,jj)/ERAD
                enddo
              ENDIF
            else
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSX(I,J) = (PS(ip1,J)-PS(im1,J))*wrk2(i,j)*wrk1(i,j)
                  PSY(I,J) = - (PS(II,J)-PS(I,J+1))*wrk3(i,j)/ERAD
                enddo
              ELSE                             !North pole point, compute at j=2
                jj = 2
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSX(I,J) = (PS(ip1,jj)-PS(im1,jj))*wrk2(i,jj)*wrk1(i,jj)
                  PSY(I,J) = - (PS(I,J)-PS(I,jj+1))*wrk3(i,jj)/ERAD
                enddo
              ENDIF
            endif
          ELSE IF(J == JM) THEN                      ! Near South pole
            if(gdlat(1,j) < 0.) then ! count from north to south
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSX(I,J) = (PS(ip1,J)-PS(im1,J))*wrk2(i,j)*wrk1(i,j)
                  PSY(I,J) = (PS(I,J-1)-PS(II,J))*wrk3(i,j)/ERAD
                enddo
              ELSE                              !South pole point,compute at jm-1
                jj = jm-1
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSX(I,J) = (PS(ip1,JJ)-PS(im1,JJ))*wrk2(i,jj)*wrk1(i,jj)
                  PSY(I,J) = (PS(I,jj-1)-PS(I,J))*wrk3(i,jj)/ERAD
                enddo
              ENDIF
            else
              IF(cosl(1,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSX(I,J) = (PS(ip1,J)-PS(im1,J))*wrk2(i,j)*wrk1(i,j)
                  PSY(I,J) = - (PS(I,J-1)-PS(II,J))*wrk3(i,j)/ERAD
                enddo
              ELSE                              !South pole point,compute at jm-1
                jj = jm-1
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSX(I,J) = (PS(ip1,JJ)-PS(im1,JJ))*wrk2(i,jj)*wrk1(i,jj)
                  PSY(I,J) = - (PS(I,jj-1)-PS(I,J))*wrk3(i,jj)/ERAD
                enddo
              ENDIF
            endif
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
!
        ENDDO                               ! end of J loop

        deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)
     
!     END IF 

      END SUBROUTINE CALGRADPS
