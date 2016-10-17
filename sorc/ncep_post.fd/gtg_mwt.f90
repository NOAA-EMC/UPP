module gtg_mountainwave

  use ctlblk_mod, only: jsta,jend,jsta_2l, jend_2u, jsta_m, jend_m, &
       jsta_m2, jend_m2,im,jm,lm, modelname,global
  use ctlblk_mod, only: SPVAL
  use params_mod, only: D608,H1,SMALL
   use physcons, only : RD=>con_rd

  use gtg_config, only : SMALL1,SMALL2,DRADDEG
  use gtg_config, only : icoord,isentropic_coord,sigma_coord,p_coord,z_coord
  use gtg_filter, only : fillybdys2d,filt2d

  implicit none

contains

!-----------------------------------------------------------------------
  subroutine mwt_init(zm,ugm,vgm,wm,Tm,pm,qvm,Rim,Nsqm, &
       truelat1,truelat2,stand_lon,latg,long,hmean,msfx,msfy,dx,dy,mwfilt,mws)
!     --- Computes low-level (lowest 1500 m) parameters used in MWT algorithms.

    implicit none

    real,dimension(im,jsta_2l:jend_2u,LM),intent(in) :: zm,ugm,vgm,wm,Tm,pm,qvm,Rim,Nsqm
    real,intent(in) :: truelat1,truelat2,stand_lon
    real,dimension(im,jsta_2l:jend_2u),intent(in) :: latg,long,hmean
    real,dimension(im,jsta_2l:jend_2u),intent(in) :: msfx,msfy,dx,dy
    real,dimension(im,jsta_2l:jend_2u),intent(inout) :: mwfilt
    real,dimension(im,jsta_2l:jend_2u),intent(inout) :: mws

!   --- local arrays
    real, dimension(im,jsta_2l:jend_2u,LM) :: um,vm
    integer, parameter :: Nmwtd=20
    real, dimension(im,jsta_2l:jend_2u,Nmwtd) :: mwtd
!-----------------------------------------------------------------------
!     --- Computes low-level (lowest 1500 m) parameters used in MWT
!     --- algorithms.  On output these are stored in the mwtd(nx,ny,17)
!     --- defined as follows:
!     ---   mwtd(i,j,1)=hmaxt (m)
!     ---   mwtd(i,j,2)=Umaxt (m/s)
!     ---   mwtd(i,j,3)=dragxt (Nt/m^2)
!     ---   mwtd(i,j,4)=speedmaxt (m/s)
!     ---   mwtd(i,j,5)=gradht (m/km)
!     ---   mwtd(i,j,6)=Nmaxt (1/s)
!     ---   mwtd(i,j,7)=Rimaxt (ND)
!     ---   mwtd(i,j,8)=wmaxt (m/s)
!     ---   mwtd(i,j,9)=avg wind direction (deg)
!     ---   mwtd(i,j,10)=unsmoothed umaxt*hmaxt (m^2/s)
!     ---   mwtd(i,j,11)=unsmoothed wmaxt*hmaxt (m^2/s)
!     ---   mwtd(i,j,12)=mwfilt (ND - 0 or 1)
!     ---   mwtd(i,j,13)=avg u (m/s)
!     ---   mwtd(i,j,14)=avg v (m/s)
!     ---   mwtd(i,j,15)=tausx=rho*ux*Nsqm)avg  (Nt/m^3)
!     ---   mwtd(i,j,16)=tausy=rho*vy*Nsqm)avg  (Nt/m^3)
!     ---   mwtd(i,j,17)=tauss=rho*speed*Nsqm)avg  (Nt/m^3)
!     --- On output mwfilt(i,j)=1 if in designated mtn wave region, 0 otherwise
!     --- On output mws(i,j) is the mtn wave multiplier (m^2/s)

    integer :: i,j,k,ii,iii,ip1,im1,jp1,jm1,kp1,km1,idx
    integer :: ii1,ii2,jj1,jj2
    real :: umax,smax,unmax,stmax,Rimax,wmax,htmax
    real :: dxm,dym
    real :: dragx,dhdx,dhdy,gradh
    real :: ux,vy,uavg,vavg,beta0
    real :: N,rhouNavg,rhovNavg,rhosNavg
    real :: dqv,Tvk,pijk,rho
    real :: ht,speed
    !
    integer :: ijk,jj,kk
    integer :: idel,jdel
    integer :: im3,ip3
    real :: aream
    real :: cm,cms,cmw,mwf,sms,hms,UNms,ums,wms
    integer :: idir
    integer :: Filttype, nsmooth

!   --- define min topographic height and search depth for max mtn top winds
    real,parameter :: hmin=500.    ! m
    real,parameter :: gradhmin=2.5 ! m/km

!   --- define min topographic height and search depth for max mtn top winds
    real,parameter :: bldepth=1500. ! m  RDS 04-08-2014
    integer,parameter ::  mwsflag = 1 ! 1=speed, 2=w

!     --- If computing MWT indices get surface parameters
    write(*,*) 'enter mwt_init'

!   --- Initialization
    mwtd = SPVAL

!   --- Get geographic um,vm from input grid relative ugm, vgm
    idir=+1	! geographic winds from grid-relative winds
    call rotu(truelat1,truelat2,stand_lon,latg,long,idir,ugm,vgm,um,vm)

!   --- Get mountain top pbl parameters
    do j=jsta,jend
       jp1=j-1
       jm1=j+1
       if(jp1<1) jp1=1
       if(jm1>jm) jm1=jm
    do i=1,IM
       ip1=i+1
       im1=i-1
       if(im1<1) then
          if(modelname == 'GFS' .or. global) then
             im1=im1+IM
          else
             im1=1
          end if
       end if
       if(ip1>IM) then
          if(modelname == 'GFS' .or. global) then
             ip1=ip1-IM
          else
             ip1=im
          end if
       endif

       htmax=0.
       umax=0.
       smax=0.
       unmax=0.
       stmax=0.
       Rimax=0.
       wmax=0.
       ijk=0
       uavg=0.
       vavg=0.
       rhouNavg=0.
       rhovNavg=0.
       rhosNavg=0.

!      --- define search perimeter for maximum mtn top winds (idel,jdel)
       idel=1
       jdel=1
       ii1=i-idel
       ii2=i+idel
       jj1=MAX(j-jdel,1+1)
       jj2=MIN(j+jdel,JM-1)
       do jj=jj2,jj1,-1 ! post is north-south, original GTG is south-north
          do iii=ii1,ii2
             ii = iii
             if(ii < 1) then
                if(modelname == 'GFS' .or. global) then
                   ii = ii + IM
                else
                   ii = 1
                endif
             elseif(ii > IM) then
                if(modelname == 'GFS' .or. global) then
                   ii = ii - IM
                else
                   ii = IM
                endif
             end if
             ht = hmean(ii,jj)
             htmax = MAX(htmax,ht)
             dxm=dx(ii,jj)/msfx(ii,jj)
             dym=dy(ii,jj)/msfy(ii,jj)
             do k=LM,1,-1  ! GFS is top-bottom, original GTG is bottom-top
                kp1=k-1
                km1=k+1
                if(k==LM) km1=LM
                if(k==1) kp1=1

!               --- Record the maxiumum values of U,N,U*N,Ri in the lowest
!               --- bldepth meters above the terrain
                if(zm(ii,jj,k)>ht+bldepth) exit
                ijk=ijk+1
                ux = um(ii,jj,k)
                vy = vm(ii,jj,k)
                speed = SQRT(ux**2 + vy**2)
                smax=MAX(smax,speed)
                umax=MAX(umax,ABS(ux))
                wmax=MAX(wmax,ABS(wm(ii,jj,k)))
                uavg=uavg+ux
                vavg=vavg+vy
!               --- derive density from p and Tv
                dqv = MAX(qvm(ii,jj,k),0.)
                Tvk = Tm(ii,jj,k)*(H1+D608*dqv) !Tv from specific humidity
                pijk = pm(ii,jj,k)
                rho = pijk/(Rd*Tvk)
                if(Nsqm(ii,jj,k)>0.) then
                   N=SQRT(Nsqm(ii,jj,k))
                   stmax=MAX(stmax,N)
                   rhouNavg=rhouNavg+rho*ux*N
                   rhovNavg=rhovNavg+rho*vy*N
                   rhosNavg=rhosNavg+rho*speed*N
                endif
                Rimax=MAX(Rimax,Rim(ii,jj,k))
             enddo  ! k loop
             UNmax=umax*stmax
          enddo ! jj loop
       enddo ! ii loop
!      --- Save mountain top parameters in the (idel,jdel) box
!      --- surrounding (i,j) in the mwtd array
       ux = uavg/MAX(FLOAT(ijk),1.)
       vy = vavg/MAX(FLOAT(ijk),1.)
!       --- get average wind direction in the box surrounding i,j
       if((ABS(ux)<SMALL).and.(ABS(vy)<SMALL)) then
          beta0=0.             ! wind dir indeterminate
       else
          beta0=ATAN2(-ux,-vy)  ! wind dir (radians) 
       endif
       beta0 = beta0/DRADDEG  ! wind dir (deg)
       if(beta0<0.)   beta0=beta0+360.
       if(beta0>=360.) beta0=beta0-360.
!      --- Compute |grad(ht)|
       dxm=dx(i,j)/msfx(i,j)
       dym=dy(i,j)/msfy(i,j)
       dhdx=dreg(hmean(im1,j),hmean(i,j),hmean(ip1,j),dxm)
       dhdy=dreg(hmean(i,jm1),hmean(i,j),hmean(i,jp1),dym)
       gradh=SQRT(dhdx**2+dhdy**2)
!
       rhouNavg=rhouNavg/MAX(FLOAT(ijk),1.)
       rhovNavg=rhovNavg/MAX(FLOAT(ijk),1.)
       rhosNavg=rhosNavg/MAX(FLOAT(ijk),1.)
!      --- Save the maximum U,speed,w,N,UN,Ri in the box surrounding i,j
       mwtd(i,j,1)=htmax     ! m
       mwtd(i,j,2)=Umax      ! E-W wind m/s
       mwtd(i,j,4)=smax      ! speed m/s
!      mwtd(i,j,5)=UNmax     ! U*N m/s^2
       mwtd(i,j,5)=1000.*gradh  ! m/km
       mwtd(i,j,6)=stmax     ! N s^-1
       mwtd(i,j,7)=Rimax     ! Ri ND
       mwtd(i,j,8)=wmax      ! w m/s
       mwtd(i,j,9)=beta0     ! deg
!      mwtd(i,j,12)=mwfilt
       mwtd(i,j,13)=ux       ! avg u (m/s)
       mwtd(i,j,14)=vy       ! avg v (m/s)
       mwtd(i,j,15)=rhouNavg ! tausx=rho*ux*Nsqm)avg  (Nt/m^3)
       mwtd(i,j,16)=rhovNavg ! tausy=rho*vy*Nsqm)avg  (Nt/m^3)
       mwtd(i,j,17)=rhosNavg ! tauss=rho*speed*Nsqm)avg  (Nt/m^3)
    enddo  ! i loop
    enddo  ! j loop

    call fillybdys2d(mwtd(1:IM,jsta_2l:jend_2u,5))
!    --- Use 1-2-1 smoother nsmooth times on gradht.  This is designed to
!    --- fill in regions where terrain differences between grid points are small.
    Filttype=1  ! 1-2-1 smoother
    nsmooth=5
    call filt2d(nsmooth,Filttype,mwtd(1:IM,jsta_2l:jend_2u,5))

!   --- Compute wave drag in x =integral p*dhdx and store in mwtd(i,j,3)
    do j=jsta_m2,jend_m2
    do i=1,IM
       dragx=0.
!       --- At each (i,j) point within the MWT region compute the
!       --- local dragx based on the integral of pdh/dx from 3 points
!       --- to the left to 3 points to the right of the (i,j) pt.  Also
!       --- average over 3 y points.
       im3=i-3
       ip3=i+3
       jp1=MAX(j-1,1+1)  ! post is north-south, original GTG is south-north
       jm1=MIN(j+1,JM-1) ! post is north-south, original GTG is south-north
       aream=0.
       do jj=jm1,jp1,-1  ! post is north-south, original GTG is south-north
          dym=dy(i,jj)/msfy(i,jj)
          aream=aream+dym
          do iii=im3,ip3
             ii = iii
             if(ii < 1) then
                if(modelname == 'GFS' .or. global) then
                   ii = ii + IM
                else
                   ii = 1
                endif
             elseif(ii > IM) then
                if(modelname == 'GFS' .or. global) then
                   ii = ii - IM
                else
                   ii = IM
                endif
             end if

             dxm=dx(ii,jj)/msfx(ii,jj)
             ip1=ii+1
             im1=ii-1
             if(im1<1) then
                if(modelname == 'GFS' .or. global) then
                   im1=im1+IM
                else
                   im1=1
                end if
             end if
             if(ip1>IM) then
                if(modelname == 'GFS' .or. global) then
                   ip1=ip1-IM
                else
                   ip1=im
                end if
             endif
             if(ABS(pm(ii,jj,1)-SPVAL)<SMALL1 .or. &
                ABS(hmean(ip1,jj)-SPVAL)<SMALL1 .or. &
                ABS(hmean(im1,jj)-SPVAL)<SMALL1) cycle
             dhdx = (hmean(ip1,jj)-hmean(im1,jj))/(2.*dxm)
             dragx = dragx + (pm(ii,jj,1)*dhdx)*dxm
          enddo  ! ii loop
       enddo  ! jj loop
       mwtd(i,j,3)=dragx/MAX(aream,1.)
    enddo  ! j loop
    enddo  ! i loop

    call fillybdys2d( mwtd(1:IM,jsta_2l:jend_2u,3))

!   --- Get MWT flag defined as an (i,j) where h>hmin and gradh>gradhmin 
    do j=jsta,jend
    do i=1,IM
       mwfilt(i,j)=0
       hms=mwtd(i,j,1)    ! hmaxt
       hms=MIN(hms,3000.) ! ~10,000 ft
!      hms=MIN(hms,2750.) ! 9,000 ft
       gradh=mwtd(i,j,5)
       if(hms>=hmin .and. gradh>=gradhmin) then
          mwfilt(i,j)=1
       endif
    enddo
    enddo
!   --- Smooth mwfilt
    nsmooth=2
    Filttype=1  ! 1-2-1 smoother
    call filt2d(nsmooth,Filttype,mwfilt)
!
    do j=jsta,jend
    do i=1,IM
       mwtd(i,j,12)=mwfilt(i,j)
    end do
    end do

!   --- Get MWT diagnostic multiplier
    do j=jsta,jend
    do i=1,IM
       cm =0.
       cms=0.
       cmw=0.
       mws(i,j)=0.
       mwf=mwfilt(i,j)
       hms=mwtd(i,j,1)  ! hmaxt
       if(mwf >= SMALL1 .or. ABS(hms-SPVAL) > SMALL1) then
          sms=mwtd(i,j,4)  ! speedmaxt
          wms=mwtd(i,j,8)  ! wmaxt
          if(ABS(sms-SPVAL) > SMALL1) cms=MAX(sms*hms,0.)
          if(ABS(wms-SPVAL) > SMALL1) cmw=MAX(wms*hms,0.)
       endif
       mwtd(i,j,10)=cms
       mwtd(i,j,11)=cmw

       if(mwsflag==1) then ! 1=speed, 2=w
          cm=cms
       else
          cm=cmw
       end if
       mws(i,j)=cm   ! Multiplier for MWT diagnostics

    end do
    end do

!   --- Smooth mws
    nsmooth=1
    Filttype=1  ! 1-2-1 smoother
    call filt2d(nsmooth,Filttype,mws)

    write(*,*) "mwtd 1=",mwtd(IM/2,jsta,1)
    write(*,*) "mwtd 2=",mwtd(IM/2,jsta,2)
    write(*,*) "mwtd 3=",mwtd(IM/2,jsta,3)
    write(*,*) "mwtd 4=",mwtd(IM/2,jsta,4)
    write(*,*) "mwtd 5=",mwtd(IM/2,jsta,5)
    write(*,*) "mwtd 6=",mwtd(IM/2,jsta,6)
    write(*,*) "mwtd 7=",mwtd(IM/2,jsta,7)
    write(*,*) "mwtd 8=",mwtd(IM/2,jsta,8)
    write(*,*) "mwtd 9=",mwtd(IM/2,jsta,9)
    write(*,*) "mwtd 10=",mwtd(IM/2,jsta,10)
    write(*,*) "mwtd 11=",mwtd(IM/2,jsta,11)
    write(*,*) "mwtd 12=",mwtd(IM/2,jsta,12)
    write(*,*) "mwtd 13=",mwtd(IM/2,jsta,13)
    write(*,*) "mwtd 14=",mwtd(IM/2,jsta,14)
    write(*,*) "mwtd 15=",mwtd(IM/2,jsta,15)
    write(*,*) "mwtd 16=",mwtd(IM/2,jsta,16)
    write(*,*) "mwtd 17=",mwtd(IM/2,jsta,17)

    return
  end subroutine mwt_init

!-----------------------------------------------------------------------
  subroutine ucritl(kmin,kmax,hgt,mwfilt,zm,um,vm,ucrit)
!     --- Computes proximity to critical level where speed~0.

    implicit none

    integer,intent(in) :: kmin,kmax
    real,dimension(1:IM,JSTA:JEND),intent(in) :: hgt
    real,dimension(im,jsta_2l:jend_2u),intent(in) :: mwfilt
    real,dimension(1:IM,JSTA:JEND,1:LM),intent(in) :: zm,um,vm
    real,dimension(1:IM,JSTA:JEND,1:LM),intent(inout) :: ucrit

    real :: spdmaxt(im,jsta_2l:jend_2u) ! work array

    integer :: i,j,k,k1,k2,kcrit,kk
    real :: zk,zkk,ucl,ucritmax,speedk,speedkk
    real,parameter :: bldepth=1500.    ! m
    logical :: break_k

    write(*,*) 'in ucritl'

    ucritmax=0.
    k1=MIN(kmax,LM-1) ! GFS is top-bottom, original GTG is bottom-top
    k2=MAX(kmin,2)
    do j=jsta,jend
    do i=1,IM
       if(mwfilt(i,j)<=0.) cycle
!      --- Look for critical levels only above the PBL (>~1500 m)
       k1=LM-1
       do k=LM-1,1,-1 ! GFS is top-bottom, original GTG is bottom-top
          speedk=SQRT(um(i,j,k)**2+vm(i,j,k)**2)
          spdmaxt(i,j)=MAX(spdmaxt(i,j),speedk)
          if(zm(i,j,k)>=hgt(i,j)+bldepth) then
             k1=k
             exit
          endif
       enddo
       k1=MIN(k1,LM-1)  ! GFS is top-bottom, original GTG is bottom-top
       k1=MAX(k1,k2)
       do k=LM,k1+1,-1  ! GFS is top-bottom, original GTG is bottom-top
          ucrit(i,j,k)=0.
       enddo
       break_k=.false.
       do k=k1,k2,-1  ! GFS is top-bottom, original GTG is bottom-top
          zk = zm(i,j,k)
          ucrit(i,j,k)=0.
!         --- Don't include uncomputed (i,j,k) or pts below terrain
          if(ABS(um(i,j,k)-SPVAL)<=SMALL1 .or. &
             ABS(vm(i,j,k)-SPVAL)<=SMALL1 .or. &
             ABS(zm(i,j,k)-SPVAL)<=SMALL1) cycle
!         --- Look for critical levels only above the PBL (>~1500 m)
          if(zm(i,j,k)<hgt(i,j)+bldepth) cycle
          speedk=SQRT(um(i,j,k)**2+vm(i,j,k)**2)
          kcrit=-1
          if(speedk<=0.5) then
!            --- speed < 0.5 m/s
             kcrit=k
             if(kcrit<=1 .or. kcrit>LM) exit
!            --- Extend influence of CL downward by cldepth
             zk=zm(i,j,kcrit)
             do kk=kcrit,k1 ! GFS is top-bottom, original GTG is bottom-top
                zkk=zm(i,j,kk)
                speedkk=SQRT(um(i,j,kk)**2+vm(i,j,kk)**2)
                if(zk-zkk<bldepth) then
                   ucl=spdmaxt(i,j)**2/(MAX(0.1,speedkk**2))
                   ucl=MAX(ucl,SMALL)
                   ucl=SQRT(ucl)
                   ucrit(i,j,kk)=ucl
                   if(ucrit(i,j,kk)>ucritmax) ucritmax=ucrit(i,j,k)
                else
                   break_k=.true.
                   exit
                endif
             enddo
!            --- Use only the lowest cl above the BL
             if(break_k) exit
          endif
       enddo ! k loop

    enddo  ! i loop
    enddo  ! j loop

    return
  end subroutine ucritl

!-----------------------------------------------------------------------
! Not used any more
!
!  subroutine tke_gwbM(pm,zm,Tm,qvm,qcm,ugm,vgm,thetav,hgt,
!     1  hvar,hpbl,mwfilt,nx,ny,nz,imin,imax,jmin,jmax,GWB,
!     2  printflag,ic,jc,iprt)
!     --- Computes gravity wave drag according to the original
!     --- Palmer et al formulation (QJRMS,112,1001-1039,1986).
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Not used any more
!
!      subroutine tke_gwbMz(printflag,iprt)
!     --- Computes gravity wave drag according to the original
!     --- Palmer et al formulation (QJRMS,112,1001-1039,1986).
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Not used any more
!      function Rifc(Ri2,Ric)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 'lwfinder' is not used since its ouputs 'ckmax,lwflag' are not used
! Then the related 'funclw' and 'MULLER' are deleted.
!      subroutine lwfinder(lsqi,zi,ni,dzi,ckmax,lwflag,printflag,iprt)
!      subroutine funclw(k,cerr,ierr)
!      SUBROUTINE MULLER(KN,N,RTS,MAXIT,EP1,EP2,FUNC,FNREAL,KOUNT)
!----------------------------------------------------------------------

  SUBROUTINE CQAUPIX(NN,X,Y,X0,Y0,INSIDE)
!***********************************************************************
!* PURPOSE:        DETERMINES IF A POINT (X0,Y0) IS INSIDE OR OUTSIDE
!*                 OF A POLYGON WITH VERTICES (X(J),Y(J),J=1,NN).
!*
!* DESCRIPTION:    THIS ROUTINE USES CAUCHY'S RESIDUE THEREOM.
!*                 IT INTEGRATES ALONG THE (ASSUMED) STRAIGHT LINE
!*                 SEGMENTS OF THE POLYGON.  IF THE EVALUATED LINE
!*                 INTEGRAL IS AN INTEGRAL MULTIPLE OF 2*PI, THE POINT
!*                 IN QUESTION IS INTERIOR, OTHERWISE IT IS EXTERIOR.
!*
!* CALLING INTERFACE:
!*
!*   INPUT:
!*       X,Y   - ARRAY OF POLYGON VERTICES
!*       NN    - NUMBER OF VERTICES IN THE POLYGON
!*       X0,Y0 - POINT IN QUESTION
!*
!*   OUTPUT:
!*       INSIDE - TRUE IF POINT (X0,Y0) IS INSIDE, FALSE IF OUTSIDE
!***********************************************************************
    IMPLICIT NONE

    INTEGER,intent(in) :: NN
    REAL(kind=4),intent(in) ::  X(NN), Y(NN), X0, Y0
    LOGICAL,intent(out) :: INSIDE

    INTEGER :: J, ISUM
    REAL(kind=4) ::  XX(NN),YY(NN), YTOP, XBTM, YSUM, ANGLE
!************************** BEGIN LOGIC ********************************
!
!     --- INITIALIZATIONS
    INSIDE = .FALSE.
    if(NN > 4000) then
       write(*, *) 'error in CQAUPIX: NN exceeds dimension'
       return
    endif

    ISUM = 0
!   --- TRANSLATE THE VERTICES OF THE POLYGON WITH RESPECT TO
!   --- (X0,Y0) AS THE ORIGIN
    DO J=1,NN
       XX(J)=X(J)-X0
       YY(J)=Y(J)-Y0
    end DO

!   --- LOOP THROUGH THE VERTICES OF THE POLYGON, SUMMING THE ANGLES
!   --- AT EACH VERTEX.  NOTE USING CAUCHY'S LEMMA REQUIRES COMPLEX
!   --- NUMBER NOTATION, SO THAT Z(J) = X(J) + I*Y(J), I=SQRT(-1).
    YSUM = 0.
    DO  J=1,NN-1
       YTOP = XX(J)*YY(J+1) - XX(J+1)*YY(J)
       XBTM = XX(J+1)*XX(J) + YY(J+1)*YY(J)
!      --- REV 001. PREVENT 0/0 ARGUMENT TO ATAN2 FUNCTION.
       IF((ABS(YTOP)<=SMALL) .AND. (ABS(XBTM)<=SMALL)) THEN
          ANGLE = 0.
       ELSE
          ANGLE = ATAN2(YTOP,XBTM)
       ENDIF
       YSUM = YSUM + ANGLE
    end DO
    ISUM = NINT(YSUM)

    IF(ISUM==0) THEN
       INSIDE = .FALSE.
    ELSE
       INSIDE = .TRUE.
    ENDIF

    RETURN
  END SUBROUTINE CQAUPIX

!-----------------------------------------------------------------------
  subroutine rotu(truelat1,truelat2,stand_lon,gdlat,gdlon,idir,ui,vi,uo,vo)
!     --- If idir>0, input is grid-relative winds, rotate to 
!     --- Earth-relative winds
!     --- If idir=<0, input is Earth-relative winds, rotate to 
!     --- grid-relative winds
!     --- Note: Only rotate winds for Lambert conformal, polar stereographic

    implicit none

    real,intent(in) :: truelat1,truelat2,stand_lon
    real,dimension(im,jsta_2l:jend_2u),intent(in) :: gdlat,gdlon
    integer,intent(in) :: idir
    real,dimension(im,jsta_2l:jend_2u,LM),intent(in) :: ui,vi
    real,dimension(im,jsta_2l:jend_2u,LM),intent(inout) :: uo,vo

    integer :: i,j,k
    real :: ugrid,vgrid,umet,vmet
    real :: cone,diff,alpha,cosalpha,sinalpha,hemi

    write(*,*) 'enter rotu'

!   --- Only rotate winds for Lambert conformal, polar stereographic
!   GFS is on Gaussian latlon, no rotation
    if(modelname == 'GFS') then
       uo = ui
       vo = vi
       write(*,*) "No rotation for GFS"
       return
    end if

!   --- Compute the cone factor
    cone=1.
    ! Lambert conformal projection
    if(modelname=='RAP' .or. modelname=='NAM') then
       call lc_cone(truelat1, truelat2, cone)
    endif

    if_dir: if(idir<0) then
!      --- Convert true east-west, north south velocity components to
!      --- grid relative components

       do j=jsta,jend
       do i=1,IM
          diff = stand_lon-gdlon(i,j)
          if(diff>180.) diff=diff-360.
          if(diff<-180.) diff=diff+360.
!         --- Calculate the rotation angle, alpha, in radians.  If PS cone=1
          if(gdlat(i,j) < 0.) then
             hemi =-1.0
          ELSE
             hemi = 1.0
          ENDIF
          alpha = hemi*cone*diff*DRADDEG
          cosalpha=cos(alpha)
          sinalpha=sin(alpha)
          do k=1,LM
             uo(i,j,k) = SPVAL
             vo(i,j,k) = SPVAL
             umet=ui(i,j,k)  ! true east-west
             vmet=vi(i,j,k)  ! true north-south
             if(ABS(umet-SPVAL) < SMALL1 .or. &
                ABS(vmet-SPVAL) < SMALL1) cycle
             uo(i,j,k)= umet*cosalpha + vmet*sinalpha
             vo(i,j,k)=-umet*sinalpha + vmet*cosalpha
          enddo ! k loop
       enddo ! i loop
       enddo ! j loop
    else  ! idir>=0
!       --- Compute true east-west, north-south (met) velocity components
!       --- from grid relative velocity components

       do j=jsta,jend
       do i=1,IM
!         --- Compute the grid relative velocities from east-west,
!         --- north=south (met) velocities
          diff = stand_lon-gdlon(i,j)
          if(diff>180.) diff=diff-360.
          if(diff<-180.) diff=diff+360.
!         --- Calculate the rotation angle, alpha, in radians.  If PS cone=1
          if(gdlat(i,j) < 0.) then
             hemi =-1.0
          ELSE
             hemi = 1.0
          ENDIF
          alpha = hemi*cone*diff*DRADDEG
          cosalpha=cos(alpha)
          sinalpha=sin(alpha)
          do k=1,LM
             uo(i,j,k)= SPVAL
             vo(i,j,k)= SPVAL
             ugrid=ui(i,j,k)
             vgrid=vi(i,j,k)
             if(ABS(ugrid-SPVAL) < SMALL1 .or. &
                ABS(vgrid-SPVAL) < SMALL1) cycle
             uo(i,j,k)= ugrid*cosalpha - vgrid*sinalpha
             vo(i,j,k)= ugrid*sinalpha + vgrid*cosalpha
          enddo ! k loop
       enddo ! i loop
       enddo ! j loop

    endif if_dir

    return
  end subroutine rotu

!-----------------------------------------------------------------------
  subroutine lc_cone(truelat1, truelat2, cone)
!     --- Subroutine to compute the cone factor of a Lambert Conformal projection
!     --- F77 translation from wrfv3.01 F90 routine of the same name

    IMPLICIT NONE

    real :: truelat1, truelat2
    real :: cone
    real :: rad_per_deg

!     --- Initializations
    cone=1.
    rad_per_deg=DRADDEG

!     --- First, see if this is a secant or tangent projection.  For tangent
!     --- projections, truelat1 = truelat2 and the cone is tangent to the 
!     --- Earth's surface at this latitude.  For secant projections, the cone
!     --- intersects the Earth's surface at each of the distinctly different
!     --- latitudes.  Ref: Haltiner and Martin, p. 14.
    IF (ABS(truelat1-truelat2) > 0.1) THEN
       cone = ALOG(COS(truelat1*rad_per_deg)) - &
              ALOG(COS(truelat2*rad_per_deg))
       cone = cone/(ALOG(TAN((45.0-ABS(truelat1)/2.0)*rad_per_deg)) - &
                    ALOG(TAN((45.0-ABS(truelat2)/2.0)*rad_per_deg)))
    ELSE
       cone = SIN(ABS(truelat1)*rad_per_deg )
    ENDIF

    RETURN
  END subroutine lc_cone

!-----------------------------------------------------------------------
  function dreg(f1,f2,f3,dx)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK 
!     --- Computes df/dx on a regular grid.
!     --- Uses centered differences unless one of f1,f2,f3 is
!     --- missing, in which case one-sided differences are attempted.
!     --- If that is unsuccessful a missing value is returned.
!     --- f1,f2,f3 = f(i-1), f(i), f(i+1)
!$$$
    implicit none
    real, intent(in) :: f1,f2,f3,dx
    real :: dreg

    dreg=SPVAL

!   --- Don't include uncomputed (i,j,k) or pts below terrain 
    if(ABS(dx)<SMALL) return

!   --- Try two-sided difference
    if( (ABS(f1-SPVAL)>SMALL1) .and. &
        (ABS(f3-SPVAL)>SMALL1) ) then
       dreg = (f3-f1)/(2.*dx)
!   --- try 1-sided differences if near terrain
    else
       if(ABS(f2-SPVAL) < SMALL1) return
       if(ABS(f3-SPVAL)>SMALL1) then
          dreg = (f3-f2)/dx
       elseif(ABS(f1-SPVAL)>SMALL1) then
          dreg = (f2-f1)/dx
       endif
    end if
    return
  end function dreg

!-----------------------------------------------------------------------
  function dirreg(f1,f2,f3,x1,x2,x3)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK 
!     --- Computes df/dx on an irregular grid of 1d values.
!     --- Checks for missing values  
!     --- f1,f2,f3 = f(i-1), f(i), f(i+1)
!     --- x1,x2,x3 = x(i-1), x(i), x(i+1)
!$$$
    implicit none

    real, intent(in) :: f1,f2,f3,x1,x2,x3
    real :: dirreg

    real :: dx1,dx2

    dirreg=SPVAL

!   --- Don't include uncomputed (i,j,k) or pts below terrain 
    if(ABS(x1-SPVAL) < SMALL1 .or. &
       ABS(x2-SPVAL) < SMALL1 .or. &
       ABS(x3-SPVAL) < SMALL1) return

    if(ABS(f2-SPVAL) < SMALL1) return

!   --- Try two-sided difference
    dx1 = x2-x1
    dx2 = x3-x2
    if( (ABS(f1-SPVAL)>SMALL1) .and. &
        (ABS(f3-SPVAL)>SMALL1) .and. &
        (ABS(dx1)>SMALL2) .and. &
        (ABS(dx2)>SMALL2) .and. &
        (ABS(dx1+dx2)>SMALL2) ) then
       dirreg = (f3-f2)*(dx1/(dx2*(dx1+dx2))) + &
                (f2-f1)*(dx2/(dx1*(dx1+dx2)))
!   --- try 1-sided differences if near terrain
    elseif((ABS(f3-SPVAL)>SMALL1).and.(ABS(dx2)>SMALL2)) then
       dirreg = (f3-f2)/dx2
    elseif((ABS(f1-SPVAL)>SMALL1).and.(ABS(dx1)>SMALL2)) then
       dirreg = (f2-f1)/dx1
    endif

    return
  end function dirreg

end module gtg_mountainwave
