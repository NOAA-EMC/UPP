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
!     MACHINE : CRAY C-90
!$$$  
!     
!
      use vrbls2d, only: f
      use masks, only: gdlat, gdlon, dx, dy
      use params_mod, only: d00, dtr, small, erad
      use ctlblk_mod, only: jsta_2l, jend_2u, spval, modelname, global, jsta, jend,&
              im, jm, jsta_m, jend_m
      use gridspec_mod, only: gridtype, dyval

      implicit none
!
!      PARAMETER (OMEGA=7.292E-5,TWOMG=2.*OMEGA)
!
!     DECLARE VARIABLES.
!     
      REAL, dimension(im,jm), intent(in) :: UWND, VWND
      REAL, dimension(im,jm), intent(inout) :: ABSV
!
      LOGICAL OLDRD,STRD
      real, allocatable ::  wrk1(:,:), wrk2(:,:), wrk3(:,:), cosl(:,:)
      real, parameter :: omega = 7.292E-5
      real, parameter :: twomg = 2.0 * omega
      INTEGER DUM 
      INTEGER, dimension(jm) :: IHE,IHW
!
      integer I,J,imb2,imb4,ip1,im1,ii,iir,iil,jj,JMT2
      real R2DX,R2DY,DVDX,DUDY,UAVG,TPH1,TPHI
!     
!***************************************************************************
!     START CALVOR HERE.
!     
!     LOOP TO COMPUTE ABSOLUTE VORTICITY FROM WINDS.
!     
!$omp  parallel do
      DO J=JSTA_2L,JEND_2U
      DO I=1,IM
        ABSV(I,J) = SPVAL
        IF(MODELNAME .EQ. 'RAPR') ABSV(I,J) = D00
      ENDDO
      ENDDO

      DO J=JSTA_2L,JEND_2U 
        IHW(J)=-MOD(J,2)
        IHE(J)=IHW(J)+1
      ENDDO
!      print*,'dyval in CALVOR= ',DYVAL 
  
      CALL EXCH_F(UWND)
!
!$omp  parallel do
!$omp& private(r2dx,r2dy,dvdx,dudy)
      IF (MODELNAME .EQ. 'GFS' .or. global)THEN
       CALL EXCH(GDLAT(1,JSTA_2L))

       allocate (wrk1(im,jsta:jend), wrk2(im,jsta:jend),          &
     &           wrk3(im,jsta:jend), cosl(im,jsta_2l:jend_2u))
!       if(1>=jsta .and. 1<=jend)then
!        if(cos(gdlat(1,1)*dtr)<small)poleflag=.T.
!       end if 	
!       call mpi_bcast(poleflag,1,MPI_LOGICAL,0,mpi_comm_comp,iret)
       imb2=im/2
       imb4=im/4
   
       DO J=JSTA,JEND
         do i=1,im
           ip1 = i + 1
           im1 = i - 1
           if (ip1 > im) ip1 = ip1 - im
           if (im1 < 1)  im1 = im1 + im
	   cosl(i,j) = cos(gdlat(i,j)*dtr)
           IF(cosl(i,j)>=SMALL)then
	     wrk1(i,j) = 1.0 / (ERAD*cosl(i,j))
	   else
	     wrk1(i,j) =0.
	   end if    
	   if(i == im .or. i == 1)then
             wrk2(i,j) = wrk1(i,j)                                   &
     &               / ((360.+GDLON(ip1,J)-GDLON(im1,J))*DTR) !dx
           else
	     wrk2(i,j) = wrk1(i,j)                                   &
     &               / ((GDLON(ip1,J)-GDLON(im1,J))*DTR) !dx
           end if
         enddo
       enddo
       CALL EXCH(cosl(1,JSTA_2L))
       
       DO J=JSTA,JEND
         if (j == 1) then
	   if(gdlat(1,j)>0.)then ! count from north to south
            do i=1,im
             ii = i + imb2
             if (ii > im) ii = ii - im
             wrk3(i,j) = 1.0                                        &
     &               / ((180.-GDLAT(i,J+1)-GDLAT(II,J))*DTR)/ERAD !dy
            enddo
	   else ! count from south to north
	    do i=1,im
             ii = i + imb2
             if (ii > im) ii = ii - im
             wrk3(i,j) = 1.0                                        &
     &               / ((180.+GDLAT(i,J+1)+GDLAT(II,J))*DTR)/ERAD 
            enddo
	   end if	      
         elseif (j == JM) then
	  if(gdlat(1,j)<0.)then ! count from north to south
           do i=1,im
             ii = i + imb2
             if (ii > im) ii = ii - im
             wrk3(i,j) = 1.0                                        &
     &                 / ((180.+GDLAT(i,J-1)+GDLAT(II,J))*DTR)/ERAD
           enddo
	  else ! count from south to north
	   do i=1,im
             ii = i + imb2
             if (ii > im) ii = ii - im
             wrk3(i,j) = 1.0                                        &
     &                 / ((180.-GDLAT(i,J-1)-GDLAT(II,J))*DTR)/ERAD
           enddo
	  end if  
         else
           do i=1,im
             wrk3(i,j) = 1.0 / ((GDLAT(I,J-1)-GDLAT(I,J+1))*DTR)/ERAD
           enddo
         endif
       enddo  

       DO J=JSTA,JEND
         DO I=1,IM
           ip1 = i + 1
           im1 = i - 1
           if (ip1 > im) ip1 = ip1 - im
           if (im1 < 1)  im1 = im1 + im
           ii = i + imb2
           if (ii > im) ii = ii - im
	         
	   IF(J == 1)then ! Near North or South pole
	    IF(cosl(i,j)>=SMALL)THEN !not a pole point
             ABSV(I,J) =(VWND(ip1,J)-VWND(im1,J)) * wrk2(i,j)          &
     &                 +(UWND(II,J)*COSL(II,J)                         &
     &                 + UWND(I,J+1)*COSL(I,J+1))*wrk3(i,j)/cosl(i,j)  &
     &                 + F(I,J)
            ELSE !pole point, compute at j=2
	     jj=2
             ABSV(I,J) =(VWND(ip1,JJ)-VWND(im1,JJ)) * wrk2(i,jj)  &
     &                 +(UWND(I,J)*COSL(I,J)+ UWND(I,Jj+1)  &
     &                 *COSL(I,Jj+1))*abs(wrk3(i,jj))/cosl(i,jj)  &
     &                 + F(I,Jj)
	    END IF 
           ELSE IF(J == JM)THEN ! Near North or South Pole
	    IF(cosl(i,j)>=SMALL)THEN !not a pole point
             ABSV(I,J) = (VWND(ip1,J)-VWND(im1,J))* wrk2(i,j)         &
     &                 +(UWND(I,J-1)*COSL(I,J-1)                      &
     &                 + UWND(II,J)*COSL(II,J))*wrk3(i,j)/cosl(i,j)   &
     &                 + F(I,J)
            ELSE !pole point, compute at jm-1
	     jj=jm-1
             ABSV(I,J) =(VWND(ip1,JJ)-VWND(im1,JJ)) * wrk2(i,jj) &
     &                 +(UWND(I,Jj-1)*COSL(I,Jj-1)                 &
     &                 + UWND(I,J)*COSL(I,J))*abs(wrk3(i,jj))/cosl(i,jj) &
     &                 + F(I,Jj)
	    END IF
           ELSE
	     ABSV(I,J) = (VWND(ip1,J)-VWND(im1,J))* wrk2(i,j)          &
     &	               -(UWND(I,J-1)*COSL(I,J-1)                       &
     &                 - UWND(I,J+1)*COSL(I,J+1))*wrk3(i,j)/cosl(i,j)  &
     &                 + F(I,J) 
           END IF
           if(ABSV(I,J)>1.0)print*,'Debug CALVOR',i,j,VWND(ip1,J),VWND(im1,J), &
           wrk2(i,j),UWND(I,J-1),COSL(I,J-1),UWND(I,J+1),COSL(I,J+1),wrk3(i,j),cosl(i,j),F(I,J),ABSV(I,J)
         END DO
       END DO  	  
!       deallocate (wrk1, wrk2, wrk3, cosl)
! GFS use lon avg as one scaler value for pole point
       call poleavg(IM,JM,JSTA,JEND,SMALL,COSL(1:IM,JSTA:JEND)         &
     &	      ,SPVAL,ABSV(1:IM,JSTA:JEND))
       deallocate (wrk1, wrk2, wrk3, cosl)
     
      ELSE IF(GRIDTYPE == 'A')THEN
       DO J=JSTA_M,JEND_M
        JMT2=JM/2+1
        TPHI=(J-JMT2)*(DYVAL/1000.)*DTR
        DO I=2,IM-1
         IF(VWND(I+1,J).LT.SPVAL.AND.VWND(I-1,J).LT.SPVAL.AND.         &
     &      UWND(I,J+1).LT.SPVAL.AND.UWND(I,J-1).LT.SPVAL) THEN          
          R2DX   = 1./(2.*DX(I,J))
          R2DY   = 1./(2.*DY(I,J))
          DVDX   = (VWND(I+1,J)-VWND(I-1,J))*R2DX
          DUDY   = (UWND(I,J+1)-UWND(I,J-1))*R2DY
	  UAVG=0.25*(UWND(I+1,J)+UWND(I-1,J)                           &
     &              +UWND(I,J+1)+UWND(I,J-1))
!  is there a (f+tan(phi)/erad)*u term?
          ABSV(I,J)=DVDX-DUDY+F(I,J)+UAVG*TAN(GDLAT(I,J)*DTR)/ERAD  ! not sure about this???
         END IF
	END DO
       END DO

      ELSE IF (GRIDTYPE == 'E')THEN
       DO J=JSTA_M,JEND_M
        JMT2=JM/2+1
        TPHI=(J-JMT2)*(DYVAL/1000.)*DTR
        DO I=2,IM-1
	 IF(VWND(I+IHE(J),J).LT.SPVAL.AND.VWND(I+IHW(J),J).LT.SPVAL    &
     &      .AND.UWND(I,J+1).LT.SPVAL.AND.UWND(I,J-1).LT.SPVAL) THEN          
          R2DX   = 1./(2.*DX(I,J))
          R2DY   = 1./(2.*DY(I,J))
	  DVDX   = (VWND(I+IHE(J),J)-VWND(I+IHW(J),J))*R2DX
          DUDY   = (UWND(I,J+1)-UWND(I,J-1))*R2DY
	  UAVG=0.25*(UWND(I+IHE(J),J)+UWND(I+IHW(J),J)                 &
     &              +UWND(I,J+1)+UWND(I,J-1))
!  is there a (f+tan(phi)/erad)*u term?
          ABSV(I,J)=DVDX-DUDY+F(I,J)+UAVG*TAN(TPHI)/ERAD 
       	 END IF
	END DO
       END DO
      ELSE IF (GRIDTYPE == 'B')THEN
       CALL EXCH_F(VWND)
       DO J=JSTA_M,JEND_M
        JMT2=JM/2+1
        TPHI=(J-JMT2)*(DYVAL/1000.)*DTR
        DO I=2,IM-1         
          R2DX   = 1./DX(I,J)
          R2DY   = 1./DY(I,J)
	  DVDX=(0.5*(VWND(I,J)+VWND(I,J-1))-0.5*(VWND(I-1,J) &
	       +VWND(I-1,J-1)))*R2DX
	  DUDY=(0.5*(UWND(I,J)+UWND(I-1,J))-0.5*(UWND(I,J-1) &
	       +UWND(I-1,J-1)))*R2DY
	  UAVG=0.25*(UWND(I-1,J-1)+UWND(I-1,J)                 &
     &              +UWND(I,J-1)+UWND(I,J))
!  is there a (f+tan(phi)/erad)*u term?
          ABSV(I,J)=DVDX-DUDY+F(I,J)+UAVG*TAN(TPHI)/ERAD 
	END DO
       END DO 
      END IF 
     
!     
!     END OF ROUTINE.
!     
      RETURN
      END
