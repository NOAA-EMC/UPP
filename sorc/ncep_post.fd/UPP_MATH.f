!> @file
!>
!> @brief upp_math is a collection of UPP subroutines for numerical math functions calculation.
!> @author Jesse Meng @date 2020-05-20

!> dvdxdudy() computes dudy, dvdx, uwnd
!>
!> h2u(), h2v(), u2h(), v2h() interpolate variables between U, V, H, points
!> adopted from UPP subroutine GRIDAVG.f
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2020-05-20 | Jesse Meng | Initial
!>
!> @author Jesse Meng @date 2020-05-20
  module upp_math

  use masks,        only: dx, dy
  use ctlblk_mod,   only: im, jsta_2l, jend_2u, jsta_m, jend_m, spval,&
                              ista_2l, iend_2u, ista_m, iend_m
  use gridspec_mod, only: gridtype
!
  implicit none

  private

  public :: DDVDX, DDUDY, UUAVG
  public :: dvdxdudy
  public :: H2U, H2V, U2H, V2H

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

      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in)    :: UWND, VWND
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
        DO I=ISTA_M,IEND_M
           IF(VWND(I+1,J)<SPVAL.AND.VWND(I-1,J)<SPVAL.AND.              &
     &        UWND(I,J+1)<SPVAL.AND.UWND(I,J-1)<SPVAL) THEN
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
          DO I=ISTA_M,IEND_M
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
!$omp parallel do  private(i,j,r2dx,r2dy)
        DO J=JSTA_M,JEND_M
          DO I=ISTA_M,IEND_M
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

      use ctlblk_mod, only: spval, jsta, jend, jsta_m, jend_m, me, num_procs, jm,&
              im, jsta_2l, jend_2u, ista, iend, ista_m, iend_m, ista_2l, iend_2u 
      use gridspec_mod, only: gridtype

      implicit none

      INCLUDE "mpif.h"
      integer:: i,j,ie,iw
      real,dimension(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U),intent(in)::ingrid
      real,dimension(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U),intent(out)::outgrid
      outgrid=spval
      if(GRIDTYPE == 'A')THEN
       do j=jsta,jend
        do i=ista,iend
	 outgrid(i,j)=ingrid(i,j)
	end do
       end do
      else IF(GRIDTYPE == 'E')THEN
       call exch(ingrid(ista_2l,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=ISTA_M,IEND_M
	 IE=I+MOD(J,2)
         IW=IE-1
         outgrid(i,j)=(ingrid(IW,J)+ingrid(IE,J)+ingrid(I,J+1)+ingrid(I,J-1))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'B')THEN
       call exch(ingrid(ista_2l,jsta_2l))
       DO J=JSTA,JEND_M
        DO I=ISTA,IEND_M
	 outgrid(i,j)=(ingrid(i,j)+ingrid(i,j+1)+ingrid(i+1,j)+ingrid(i+1,j+1))/4.0
	end do
       end do
! Fill in boundary points because hysplit fails when 10 m wind has bitmaps
       do j=jsta,jend_m
        outgrid(iend,j)=outgrid(iend-1,j)
       end do
       IF(me == (num_procs-1) .and. jend_2u >= jm) then
        DO I=ISTA,IEND
         outgrid(i,jend) = outgrid(i,jend-1)
        END DO
       END IF      
      ELSE IF(GRIDTYPE == 'C')THEN
       DO J=JSTA,JEND
        DO I=ISTA,IEND_M
          outgrid(i,j)=(ingrid(i,j)+ingrid(i+1,j))/2.0
        end do
       end do
      end if 
      	 
      end subroutine H2U
!
!-------------------------------------------------------------------------------------
!	
      subroutine H2V(ingrid,outgrid)
! This subroutine interpolates from H points onto V points
! Author: CHUANG, EMC, Dec. 2010
      use ctlblk_mod, only: spval, jsta, jend, jsta_m, jend_m, im, jsta_2l, jend_2u,&
                                   ista, iend, ista_m, iend_m,     ista_2l, iend_2u
      use gridspec_mod, only: gridtype
      implicit none
      INCLUDE "mpif.h"
      integer:: i,j,ie,iw
      real,dimension(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U),intent(in)::ingrid
      real,dimension(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U),intent(out)::outgrid
      outgrid=spval
      if(GRIDTYPE == 'A')THEN
       do j=jsta,jend
        do i=ista,iend
          outgrid(i,j)=ingrid(i,j)
        end do
       end do
      else IF(GRIDTYPE == 'E')THEN
       call exch(ingrid(ista_2l,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=ISTA_M,IEND_M
	 IE=I+MOD(J,2)
         IW=IE-1
         outgrid(i,j)=(ingrid(IW,J)+ingrid(IE,J)+ingrid(I,J+1)+ingrid(I,J-1))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'B')THEN
       call exch(ingrid(ista_2l,jsta_2l))
       DO J=JSTA,JEND_M
        DO I=ISTA,IEND_M
	 outgrid(i,j)=(ingrid(i,j)+ingrid(i,j+1)+ingrid(i+1,j)+ingrid(i+1,j+1))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'C')THEN
       call exch(ingrid(ista_2l,jsta_2l))
       DO J=JSTA,JEND_M
        DO I=ISTA,IEND
	 outgrid(i,j)=(ingrid(i,j)+ingrid(i,j+1))/2.0
	end do
       end do 
      end if 
      
      end subroutine H2V
!
!-------------------------------------------------------------------------------------
!
      subroutine U2H(ingrid,outgrid)
! This subroutine interpolates from U points onto H points
! Author: CHUANG, EMC, Dec. 2010
      use ctlblk_mod, only: spval, jsta, jend, jsta_m, jend_m, im, jsta_2l, jend_2u,&
                                   ista, iend, ista_m, iend_m,     ista_2l, iend_2u
      use gridspec_mod, only: gridtype
      implicit none
      INCLUDE "mpif.h"
      integer:: i,j,ie,iw
      real,dimension(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U),intent(in)::ingrid
      real,dimension(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U),intent(out)::outgrid
      outgrid=spval
      if(GRIDTYPE == 'A')THEN
       do j=jsta,jend
        do i=ista,iend
	 outgrid(i,j)=ingrid(i,j)
	end do
       end do
      else IF(GRIDTYPE == 'E')THEN
       call exch(ingrid(ista_2l,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=ISTA_M,IEND_M
	 IE=I+MOD(J+1,2)
         IW=IE-1
         outgrid(i,j)=(ingrid(IW,J)+ingrid(IE,J)+ingrid(I,J+1)+ingrid(I,J-1))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'B')THEN
       call exch(ingrid(ista_2l,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=ISTA_M,IEND_M
	 outgrid(i,j)=(ingrid(i-1,j-1)+ingrid(i,j-1)+ingrid(i-1,j)+ingrid(i,j))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'C')THEN
       DO J=JSTA,JEND
        DO I=ISTA_M,IEND
	 outgrid(i,j)=(ingrid(i-1,j)+ingrid(i,j))/2.0
	end do
       end do       
      end if 
      	 
      end subroutine U2H       
!
!-------------------------------------------------------------------------------------
!      	 
      subroutine V2H(ingrid,outgrid)
! This subroutine interpolates from V points onto H points
! Author: CHUANG, EMC, Dec. 2010
      use ctlblk_mod, only: spval, jsta, jend, jsta_m, jend_m, im, jsta_2l, jend_2u,&
                                   ista, iend, ista_m, iend_m,     ista_2l, iend_2u
      use gridspec_mod, only: gridtype
      implicit none
      INCLUDE "mpif.h"
      integer:: i,j,ie,iw
      real,dimension(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U),intent(in)::ingrid
      real,dimension(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U),intent(out)::outgrid
      outgrid=spval
      if(GRIDTYPE == 'A')THEN
       do j=jsta,jend
        do i=ista,iend
	 outgrid(i,j)=ingrid(i,j)
	end do
       end do
      else IF(GRIDTYPE == 'E')THEN
       call exch(ingrid(ista_2l,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=ISTA_M,IEND_M
	 IE=I+MOD(J,2)
         IW=IE-1
         outgrid(i,j)=(ingrid(IW,J)+ingrid(IE,J)+ingrid(I,J+1)+ingrid(I,J-1))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'B')THEN
       call exch(ingrid(ista_2l,jsta_2l))
       DO J=JSTA_M,JEND_M
        DO I=ISTA_M,IEND_M
	 outgrid(i,j)=(ingrid(i-1,j-1)+ingrid(i,j-1)+ingrid(i-1,j)+ingrid(i,j))/4.0
	end do
       end do
      ELSE IF(GRIDTYPE == 'C')THEN
       call exch(ingrid(ista_2l,jsta_2l))
       DO J=JSTA_M,JEND
        DO I=ISTA,IEND
	 outgrid(i,j)=(ingrid(i,j-1)+ingrid(i,j))/2.0
	end do
       end do 
      end if 
      	 
      end subroutine V2H	 	         
!
!-------------------------------------------------------------------------------------
!        
  end module upp_math
