  MODULE UPP_MATH
!------------------------------------------------------------------------
! A collection of UPP subroutines for numerical math functions calculation.
!
! DVDXDUDY
! computes dudy, dvdx, uwnd
!
! H2U, H2V, U2H, V2H
! interpolates variables between U, V, H, points
! adopted from UPP subroutine GRIDAVG.f
!
! program log:
!   MAY 20 2020    Jesse Meng   Initial code
!------------------------------------------------------------------------
  use masks,        only: dx, dy
  use ctlblk_mod,   only: im, jsta_2l, jend_2u, jsta_m, jend_m, spval
  use gridspec_mod, only: gridtype
!
  implicit none

  REAL, ALLOCATABLE :: DDVDX(:,:)
  REAL, ALLOCATABLE :: DDUDY(:,:)
  REAL, ALLOCATABLE :: UUAVG(:,:)
! 
!-------------------------------------------------------------------------------------
!
  contains
!
!-------------------------------------------------------------------------------------
  subroutine dvdxdudy(uwnd,vwnd)
!
      implicit none

      REAL, dimension(im,jsta_2l:jend_2u), intent(in)    :: UWND, VWND
!
!-- local variables
!--
      integer i, j
      real r2dx, r2dy
      INTEGER, allocatable ::  IHE(:),IHW(:)
!
      IF(GRIDTYPE == 'A')THEN
!$omp parallel do  private(i,j,r2dx,r2dy)
        DO J=JSTA_M,JEND_M
        DO I=2,IM-1
           IF(VWND(I+1,J).LT.SPVAL.AND.VWND(I-1,J).LT.SPVAL.AND.              &
     &        UWND(I,J+1).LT.SPVAL.AND.UWND(I,J-1).LT.SPVAL) THEN
              R2DX   = 1./(2.*DX(I,J))
              R2DY   = 1./(2.*DY(I,J))
              DDVDX(I,J)   = (VWND(I+1,J)-VWND(I-1,J))*R2DX
              DDUDY(I,J)   = (UWND(I,J+1)-UWND(I,J-1))*R2DY
              UUAVG(I,J)   = UWND(I,J)
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
!$omp parallel do  private(i,j,r2dx,r2dy)
        DO J=JSTA_M,JEND_M
          DO I=2,IM-1
            IF(VWND(I+IHE(J),J) < SPVAL.AND.VWND(I+IHW(J),J) < SPVAL .AND.   &
     &         UWND(I,J+1) < SPVAL     .AND.UWND(I,J-1) < SPVAL) THEN
               R2DX   = 1./(2.*DX(I,J))
               R2DY   = 1./(2.*DY(I,J))
               DDVDX(I,J)   = (VWND(I+IHE(J),J)-VWND(I+IHW(J),J))*R2DX
               DDUDY(I,J)   = (UWND(I,J+1)-UWND(I,J-1))*R2DY
               UUAVG(I,J)   = 0.25*(UWND(I+IHE(J),J)+UWND(I+IHW(J),J)               &
     &                      +       UWND(I,J+1)+UWND(I,J-1))
            END IF
          END DO
        END DO
        deallocate(ihw, IHE)
      ELSE IF (GRIDTYPE == 'B')THEN
        CALL EXCH_F(VWND)
!$omp parallel do  private(i,j,r2dx,r2dy)
        DO J=JSTA_M,JEND_M
          DO I=2,IM-1
            R2DX = 1./DX(I,J)
            R2DY = 1./DY(I,J)
            if(VWND(I,  J)==SPVAL .or. VWND(I,  J-1)==SPVAL .or. &
               VWND(I-1,J)==SPVAL .or. VWND(I-1,J-1)==SPVAL .or. &
               UWND(I,  J)==SPVAL .or. UWND(I-1,J  )==SPVAL .or. &
               UWND(I,J-1)==SPVAL .or. UWND(I-1,J-1)==SPVAL) cycle
            DDVDX(I,J) = (0.5*(VWND(I,J)+VWND(I,J-1))-0.5*(VWND(I-1,J) &
     &           +       VWND(I-1,J-1)))*R2DX
            DDUDY(I,J) = (0.5*(UWND(I,J)+UWND(I-1,J))-0.5*(UWND(I,J-1) &
     &           +       UWND(I-1,J-1)))*R2DY
            UUAVG(I,J) = 0.25*(UWND(I-1,J-1)+UWND(I-1,J)               &
     &           +       UWND(I,  J-1)+UWND(I,  J))
          END DO
        END DO
      END IF

  end subroutine dvdxdudy
!
!-------------------------------------------------------------------------------------
!
      subroutine H2U(ingrid,outgrid)
! This subroutine interpolates from H points onto U points
! Author: CHUANG, EMC, Dec. 2010

      use ctlblk_mod, only: spval, jsta, jend, jsta_m, jend, me, num_procs, jm,&
              im, jsta_2l, jend_2u , jend_m
      use gridspec_mod, only: gridtype

      implicit none

      INCLUDE "mpif.h"
      integer:: i,j,ie,iw
      real,dimension(IM,JSTA_2L:JEND_2U),intent(in)::ingrid
      real,dimension(IM,JSTA_2L:JEND_2U),intent(out)::outgrid
      outgrid=spval
      if(GRIDTYPE == 'A')THEN
       do j=jsta,jend
        do i=1,im
	 outgrid(i,j)=ingrid(i,j)
	end do
       end do
      else IF(GRIDTYPE == 'E')THEN
       call exch(ingrid(1,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=2,IM-1
	 IE=I+MOD(J,2)
         IW=IE-1
         outgrid(i,j)=(ingrid(IW,J)+ingrid(IE,J)+ingrid(I,J+1)+ingrid(I,J-1))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'B')THEN
       call exch(ingrid(1,jsta_2l))
       DO J=JSTA,JEND_M
        DO I=1,IM-1
	 outgrid(i,j)=(ingrid(i,j)+ingrid(i,j+1)+ingrid(i+1,j)+ingrid(i+1,j+1))/4.0
	end do
       end do
! Fill in boundary points because hysplit fails when 10 m wind has bitmaps
       do j=jsta,jend_m
        outgrid(im,j)=outgrid(im-1,j)
       end do
       IF(me == (num_procs-1) .and. jend_2u >= jm) then
        DO I=1,IM
         outgrid(i,jm) = outgrid(i,jm-1)
        END DO
       END IF      
      ELSE IF(GRIDTYPE == 'C')THEN
       DO J=JSTA,JEND
        DO I=1,IM-1
          outgrid(i,j)=(ingrid(i,j)+ingrid(i+1,j))/2.0
        end do
       end do
      end if 
      	 
!      return 
      end subroutine H2U
!
!-------------------------------------------------------------------------------------
!	
      subroutine H2V(ingrid,outgrid)
! This subroutine interpolates from H points onto V points
! Author: CHUANG, EMC, Dec. 2010
      use ctlblk_mod, only: spval, jsta, jend, jsta_m, jend_m, im, jsta_2l, jend_2u
      use gridspec_mod, only: gridtype
      implicit none
      INCLUDE "mpif.h"
      integer:: i,j,ie,iw
      real,dimension(IM,JSTA_2L:JEND_2U),intent(in)::ingrid
      real,dimension(IM,JSTA_2L:JEND_2U),intent(out)::outgrid
      outgrid=spval
      if(GRIDTYPE == 'A')THEN
       do j=jsta,jend
        do i=1,im
          outgrid(i,j)=ingrid(i,j)
        end do
       end do
      else IF(GRIDTYPE == 'E')THEN
       call exch(ingrid(1,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=2,IM-1
	 IE=I+MOD(J,2)
         IW=IE-1
         outgrid(i,j)=(ingrid(IW,J)+ingrid(IE,J)+ingrid(I,J+1)+ingrid(I,J-1))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'B')THEN
       call exch(ingrid(1,jsta_2l))
       DO J=JSTA,JEND_M
        DO I=1,IM-1
	 outgrid(i,j)=(ingrid(i,j)+ingrid(i,j+1)+ingrid(i+1,j)+ingrid(i+1,j+1))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'C')THEN
       call exch(ingrid(1,jsta_2l))
       DO J=JSTA,JEND_M
        DO I=1,IM
	 outgrid(i,j)=(ingrid(i,j)+ingrid(i,j+1))/2.0
	end do
       end do 
      end if 
      
!      return 
      end subroutine H2V
!
!-------------------------------------------------------------------------------------
!
      subroutine U2H(ingrid,outgrid)
! This subroutine interpolates from U points onto H points
! Author: CHUANG, EMC, Dec. 2010
      use ctlblk_mod, only: spval, jsta, jend, jsta_m, jend_m, im, jsta_2l, jend_2u
      use gridspec_mod, only: gridtype
      implicit none
      INCLUDE "mpif.h"
      integer:: i,j,ie,iw
      real,dimension(IM,JSTA_2L:JEND_2U),intent(in)::ingrid
      real,dimension(IM,JSTA_2L:JEND_2U),intent(out)::outgrid
      outgrid=spval
      if(GRIDTYPE == 'A')THEN
       do j=jsta,jend
        do i=1,im
	 outgrid(i,j)=ingrid(i,j)
	end do
       end do
      else IF(GRIDTYPE == 'E')THEN
       call exch(ingrid(1,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=2,IM-1
	 IE=I+MOD(J+1,2)
         IW=IE-1
         outgrid(i,j)=(ingrid(IW,J)+ingrid(IE,J)+ingrid(I,J+1)+ingrid(I,J-1))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'B')THEN
       call exch(ingrid(1,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=2,IM-1
	 outgrid(i,j)=(ingrid(i-1,j-1)+ingrid(i,j-1)+ingrid(i-1,j)+ingrid(i,j))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'C')THEN
       DO J=JSTA,JEND
        DO I=2,IM
	 outgrid(i,j)=(ingrid(i-1,j)+ingrid(i,j))/2.0
	end do
       end do       
      end if 
      	 
!      return 
      end subroutine U2H       
!
!-------------------------------------------------------------------------------------
!      	 
      subroutine V2H(ingrid,outgrid)
! This subroutine interpolates from V points onto H points
! Author: CHUANG, EMC, Dec. 2010
      use ctlblk_mod, only: spval, jsta, jend, jsta_m, jend_m, im, jsta_2l, jend_2u
      use gridspec_mod, only: gridtype
      implicit none
      INCLUDE "mpif.h"
      integer:: i,j,ie,iw
      real,dimension(IM,JSTA_2L:JEND_2U),intent(in)::ingrid
      real,dimension(IM,JSTA_2L:JEND_2U),intent(out)::outgrid
      outgrid=spval
      if(GRIDTYPE == 'A')THEN
       do j=jsta,jend
        do i=1,im
	 outgrid(i,j)=ingrid(i,j)
	end do
       end do
      else IF(GRIDTYPE == 'E')THEN
       call exch(ingrid(1,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=2,IM-1
	 IE=I+MOD(J,2)
         IW=IE-1
         outgrid(i,j)=(ingrid(IW,J)+ingrid(IE,J)+ingrid(I,J+1)+ingrid(I,J-1))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'B')THEN
       call exch(ingrid(1,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=2,IM-1
	 outgrid(i,j)=(ingrid(i-1,j-1)+ingrid(i,j-1)+ingrid(i-1,j)+ingrid(i,j))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'C')THEN
       call exch(ingrid(1,jsta_2l))
       DO J=JSTA_M,JEND
        DO I=1,IM
	 outgrid(i,j)=(ingrid(i,j-1)+ingrid(i,j))/2.0
	end do
       end do 
      end if 
      	 
!      return 
      end subroutine V2H	 	         
!
!-------------------------------------------------------------------------------------
!        
  end module UPP_MATH
