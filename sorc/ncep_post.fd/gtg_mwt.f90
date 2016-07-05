module gtg_mountainwave

  use ctlblk_mod, only: jsta,jend,jsta_2l, jend_2u, jsta_m, jend_m, &
       jsta_m2, jend_m2,im,jm,lm, modelname,global
  use ctlblk_mod, only: SPVAL
  use params_mod, only: D608,H1,SMALL
   use physcons, only : RD=>con_rd

  use gtg_config, only : SMALL1,SMALL2,DRADDEG
  use gtg_config, only : icoord,isentropic_coord,sigma_coord,p_coord,z_coord

  implicit none

contains

!-----------------------------------------------------------------------
  subroutine get_mw_regions(latg,long,hmean,msfx,msfy,dx,dy, &
    npolygons,polypts,maxpolygons,maxpolypts,Polygonlatlon,&
    use_input_polygons,mwfilt)
!     --- Defines mountain wave regions either by using input polygons
!     --- or terrain characteristics (terrain height and gradient)
!     --- Output is in mwfilt(i,j) =1 if in mw region, 0 if not

    use gtg_filter, only : filt2d

    implicit none

    real,dimension(im,jsta_2l:jend_2u),intent(in) :: latg,long,hmean
    real,dimension(im,jsta_2l:jend_2u),intent(in) :: msfx,msfy,dx,dy
!   --- mtw area polygon declarations
    integer,intent(in) :: npolygons
    integer,intent(in) ::  polypts(maxpolygons)
    integer,intent(in) :: maxpolygons,maxpolypts
    real,intent(in) :: Polygonlatlon(maxpolygons,maxpolypts,2)
    logical,intent(in) :: use_input_polygons
!   --- output
    real,dimension(im,jsta_2l:jend_2u),intent(inout)  :: mwfilt

    real,dimension(im,jsta_2l:jend_2u) :: gradh ! grad(h) (m/km)
    integer :: i,j,k,ip1,im1,jp1,jm1
    integer :: nb
    real :: dxm,dym
    real :: dhdx,dhdy
    integer :: npts,nmpts
    logical :: inmtnbox
    integer :: ip,ib,np
    integer :: Filttype,nsmooth
    integer,parameter :: maxb=15
    real :: blat(maxb),blon(maxb)
!   --- define min topographic height and search depth for max mtn top winds
    real,parameter :: hmin=500.    ! m
    real,parameter :: gradhmin=2.5 ! m/km

    write(*,*) 'enter get_mw_regions'

!   --- Initialize the output array
    mwfilt=-1
    gradh=SPVAL

!   --- Compute terrain ht gradient (m/km)
    do j=jend_m2,jsta_m2,-1 ! post is north-south, original GTG is south-north
       jp1=j-1
       jm1=j+1
       if(jp1<1) jp1=1
       if(jm1>JM) jm1=JM
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
                ip1=IM
             end if
          endif

          dxm=dx(i,j)/msfx(i,j)
          dym=dy(i,j)/msfy(i,j)
          dhdx=dreg(hmean(im1,j),hmean(i,j),hmean(ip1,j),dxm)
          dhdy=dreg(hmean(i,jm1),hmean(i,j),hmean(i,jp1),dym)
          gradh(i,j)=1000.*SQRT(dhdx**2+dhdy**2)   ! m/km
       enddo
    enddo

    if_use: if(use_input_polygons) then
!      --- For conus models use predefined mwt polygon regions to assign
!      --- mwt zones. Get lat,lon boundaries of mtn wave region
       np=nPolygons
       nb=0
       do ip=1,np
          nb=nb+polypts(ip)
       enddo
       if(np<=0 .or. nb<=0) then
          write(*,*) 'No input MWT polygons specified: using terrain'
          ! use_input_polygons=.FALSE.
       else
          do j=jsta,jend
          do i=1,IM
             mwfilt(i,j)=0
             np=nPolygons
!            -- Polygon loop
             do ip=1,np
                nb=polypts(ip)
                do ib=1,nb
                   blat(ib)=Polygonlatlon(ip,ib,1)
                   blon(ib)=Polygonlatlon(ip,ib,2)
                enddo
!               --- Determine if input lat,lon is within the defined polygon
                inmtnbox=.FALSE.
                call cqaupix(nb,blat,blon,latg(i,j),long(i,j),inmtnbox)
                if(inmtnbox) then
                   mwfilt(i,j)=1
                   exit
                endif
             enddo  ! polygon loop
          enddo  ! j loop
          enddo  ! i loop
       endif
    else
!      --- If input polygons are unavailable, assume mwt regions are defined
!      --- where the mean terrain ht over a grid cell is >= hmin (nominally
!      --- 500 m and grad(h) >= 3.0E-3).
!      --- Use 1-2-1 smoother nsmooth times.  This is designed to fill in regions
!      --- where terrain differences between grid points are small.
       Filttype=1  ! 1-2-1 smoother
       nsmooth=6
       call filt2d(nsmooth,Filttype,gradh)
!      --- Determine which grid points satify the criteria for a mw region
       do j=jsta,jend
       do i=1,IM
          mwfilt(i,j)=0
          if(hmean(i,j)>=hmin .and. gradh(i,j)>=gradhmin) then
             mwfilt(i,j)=1
          endif
       enddo
       enddo
    endif if_use

    return
  end subroutine get_mw_regions

!-----------------------------------------------------------------------
  subroutine mwt_init(zm,ugm,vgm,wm,Tm,pm,qvm,Rim,Nsqm, &
       truelat1,truelat2,stand_lon,latg,long,hmean,msfx,msfy,dx,dy,mwfilt,mws)
!     --- Computes low-level (lowest 1500 m) parameters used in MWT algorithms.

    implicit none

    real,dimension(im,jsta_2l:jend_2u,LM),intent(in) :: zm,ugm,vgm,wm,Tm,pm,qvm,Rim,Nsqm
    real,intent(in) :: truelat1,truelat2,stand_lon
    real,dimension(im,jsta_2l:jend_2u),intent(in) :: latg,long,hmean
    real,dimension(im,jsta_2l:jend_2u),intent(in) :: msfx,msfy,dx,dy
    real,dimension(im,jsta_2l:jend_2u),intent(in) :: mwfilt
    real,dimension(im,jsta_2l:jend_2u),intent(inout) :: mws

!   --- local arrays
    real, dimension(im,jsta_2l:jend_2u,LM) :: mwtd,um,vm
!     --- On output these are stored in the mwtd(nx,ny,17) defined as follows:
!     ---   mwtd(i,j,1)=hmaxt (m)
!     ---   mwtd(i,j,2)=Umaxt (m/s)
!     ---   mwtd(i,j,3)=dragxt (Nt/m^2)
!     ---   mwtd(i,j,4)=speedmaxt (m/s)
!     ---   mwtd(i,j,5)=UNmaxt (m/s^2)
!     ---   mwtd(i,j,6)=Nmaxt (1/s)
!     ---   mwtd(i,j,7)=Rimaxt (ND)
!     ---   mwtd(i,j,8)=wmaxt (m/s)
!     ---   mwtd(i,j,9)=avg wind direction (deg)
!     ---   mwtd(i,j,10)=Txmaxt (K/m)
!     ---   mwtd(i,j,11)=gradTmaxt (K/m)
!     ---   mwtd(i,j,12)=mwfilt (ND - 0 or 1)
!     ---   mwtd(i,j,13)=avg u (m/s)
!     ---   mwtd(i,j,14)=avg v (m/s)
!     ---   mwtd(i,j,15)=tausx=rho*ux*Nsqm)avg  (Nt/m^3)
!     ---   mwtd(i,j,16)=tausy=rho*vy*Nsqm)avg  (Nt/m^3)
!     ---   mwtd(i,j,17)=tauss=rho*speed*Nsqm)avg  (Nt/m^3)

    integer :: i,j,k,ii,iii,ip1,im1,jp1,jm1,kp1,km1,idx
    integer :: ii1,ii2,jj1,jj2
    real :: umax,smax,unmax,stmax,Rimax,wmax,htmax
    real :: dTdx,dTdy,dTdz,dzdx,dzdy,Tx,Ty,gradT,Txmax,gradTmax
    real :: dxm,dym
    real :: dragx,dhdx
    real :: ux,vy,uavg,vavg,beta0
    real :: N,rhouNavg,rhovNavg,rhosNavg
    real :: dqv,Tvk,pijk,rho
    real :: ylast,ysave
    real :: ht,speed
    !
    integer :: immin,immax,jmmin,jmmax
    integer :: it1,it2,jt1,jt2,iquadrant
    integer :: ijk,jj,kk
    integer :: idel,jdel
    integer :: im3,ip3
    real :: aream
    real :: cm,mwf,sms,hms,UNms,ums,wms
    integer :: idir

!   --- define min topographic height and search depth for max mtn top winds
    real,parameter :: bldepth=1500. ! m  RDS 04-08-2014
!   --- define number of upstream points to backup to look for max terrain and w.
    integer,parameter :: nbackup=2  ! RDS 04-08-2014
    integer,parameter ::  mwtmultflag = 1 ! 1=speed, 2=w

!     --- If computing MWT indices get surface parameters
    write(*,*) 'enter mwt_init'

!   --- Initialization
    mwtd = SPVAL

!   --- Get geographic um,vm from input grid relative ugm, vgm
    idir=+1	! geographic winds from grid-relative winds
    call rotu(truelat1,truelat2,stand_lon,latg,long,idir,ugm,vgm,um,vm)

!   --- Find grid boundaries
!    immin=IM
!    immax=1
!    jmmin=jend
!    jmmax=jsta
!    do j=jsta,jend
!    do i=1,IM
!       if(mwfilt(i,j)>0) then
!          immin=MIN(immin,i)
!          immax=MAX(immax,i)
!          jmmin=MIN(jmmin,j)
!          jmmax=MAX(jmmax,j)
!       endif
!    enddo  ! j loop
!    enddo  ! i loop

!   --- Get mountain top pbl parameters
!     do i=immin,immax
!     do j=jmmin,jmmax
    do j=jsta_m2,jend_m2
    do i=1,IM
       htmax=0.
       umax=0.
       smax=0.
       unmax=0.
       stmax=0.
       Rimax=0.
       wmax=0.
       Txmax=-1.0E10
       gradTmax=0.
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
          jp1=jj-1
          jm1=jj+1
          if(jp1<1) jp1=1
          if(jm1>jm) jm1=jm
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
             ip1=ii+1
             im1=ii-1
             if(im1<1) then
                im1=1
                if(modelname == 'GFS' .or. global) im1=im1+IM
             end if
             if(ip1>IM) then
                ip1=im
                if(modelname == 'GFS' .or. global) ip1=ip1-IM
             endif

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
                Tvk = Tm(i,j,k)*(H1+D608*dqv) !Tv from specific humidity
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
!               --- Compute low-level temperature gradients
!               --- dT/dx, dT/dy on native grid
                dTdx=dreg(Tm(im1,jj,k),Tm(ii,jj,k),Tm(ip1,jj,k),dxm)
                dTdy=dreg(Tm(ii,jm1,k),Tm(ii,jj,k),Tm(ii,jp1,k),dym)
!               --- Don't include uncomputed (i,j,k) or pts below terrain v17
                if(ABS(dTdx-SPVAL)<SMALL1 .or. &
                   ABS(dTdy-SPVAL)<SMALL1)cycle
                Tx=dTdx
                Ty=dTdy
!               --- If native eta grid is not a constant z coordinate, 
!               --- transform T gradients to constant z surface by using
!                   dT/dx)z = dT/dx)eta - (dT/dz)*dz/dx)eta 
!                   dT/dy)z = dT/dy)eta - (dT/dz)*dz/dy)eta
!               --- see e.g. Haltiner and Williams p. 15.
                dTdz=SPVAL
                dzdx=SPVAL
                dzdy=SPVAL
                ! models NOT on const z coordinate
                if(icoord /= z_coord) then
                   dTdz = dirreg(Tm(ii,jj,km1),Tm(ii,jj,k),Tm(ii,jj,kp1),&
                                 zm(ii,jj,km1),zm(ii,jj,k),zm(ii,jj,kp1) )
                   if(ABS(dTdz-SPVAL)<SMALL1) cycle
                   dzdx=dreg(zm(im1,jj,k),zm(ii,jj,k),zm(ip1,jj,k),dxm)
                   dzdy=dreg(zm(ii,jm1,k),zm(ii,jj,k),zm(ii,jp1,k),dym)
!                  --- Don't include uncomputed (i,j,k) or pts below terrain v17
                   if(ABS(dzdx-SPVAL)<SMALL1 .or. &
                      ABS(dzdy-SPVAL)<SMALL1) cycle
                   Tx = dTdx - dTdz*dzdx
                   Ty = dTdy - dTdz*dzdy
                endif
!               --- Compute |delT|
                gradT=SQRT(Tx**2 + Ty**2)
                Txmax = MAX(Txmax,Tx)
                gradTmax = MAX(gradTmax,gradT)
             enddo  ! k loop
             UNmax=umax*stmax
          enddo ! jj loop
       enddo ! ii loop
!      --- Save the maximum U,speed,w,N,UN,Ri in the box surrounding i,j 
       mwtd(i,j,2)=Umax
       mwtd(i,j,4)=smax
       mwtd(i,j,18)=wmax  ! temporary storage for wmax
       mwtd(i,j,19)=htmax  ! temporary storage for hmax
       mwtd(i,j,6)=stmax
       mwtd(i,j,5)=UNmax
       mwtd(i,j,7)=Rimax
       mwtd(i,j,10)=Txmax
       mwtd(i,j,11)=gradTmax
       ux = uavg/MAX(FLOAT(ijk),1.)
       vy = vavg/MAX(FLOAT(ijk),1.)
       rhouNavg=rhouNavg/MAX(FLOAT(ijk),1.)
       rhovNavg=rhovNavg/MAX(FLOAT(ijk),1.)
       rhosNavg=rhosNavg/MAX(FLOAT(ijk),1.)
       mwtd(i,j,13)=ux
       mwtd(i,j,14)=vy
       mwtd(i,j,15)=rhouNavg
       mwtd(i,j,16)=rhovNavg
       mwtd(i,j,17)=rhosNavg
!       --- get average wind direction in the box surrounding i,j
       if((ABS(ux)<SMALL).and.(ABS(vy)<SMALL)) then
          beta0=0.             ! wind dir indeterminate
       else
          beta0=ATAN2(-ux,-vy)  ! wind dir (radians) 
       endif
       beta0 = beta0/DRADDEG  ! wind dir (deg)
       if(beta0<0.)   beta0=beta0+360.
       if(beta0>=360.) beta0=beta0-360.
       mwtd(i,j,9)=beta0
    enddo  ! i loop
    enddo  ! j loop

!   --- Look for max terrain/w in the downstream quadrant
!   --- First determine downstream quadrant 
!     do i=immin,immax
!     do j=jmmin,jmmax
    do j=jsta,jend
    do i=1,IM
       iquadrant=0
       beta0=mwtd(i,j,9)
       if((beta0>=337.7 .and. beta0<=360.0) .or. &
          (beta0>=0. .and. beta0<22.5) ) then
!         --- from north
          iquadrant=1
          it1=i
          it2=i
          jt1=j         ! post is north-south, original GTG is south-north
          jt2=j+nbackup ! post is north-south, original GTG is south-north
       elseif(beta0>=22.5 .and. beta0<112.5)then
!         --- from east
          iquadrant=2
          it1=i-nbackup
          it2=i
          jt1=j
          jt2=j
       elseif(beta0>=112.5 .and.beta0<202.5)then
!         --- from south
          iquadrant=3
          it1=i
          it2=i
          jt1=j-nbackup ! post is north-south, original GTG is south-north
          jt2=j         ! post is north-south, original GTG is south-north
        elseif(beta0>=202.5 .and.beta0<337.7)then
!         --- from west
          iquadrant=4
          it1=i
          it2=i+nbackup
          jt1=j
          jt2=j
       endif
       wmax=SPVAL
       if(iquadrant>0) then
          jt1=MAX(jt1,1)
          jt2=MIN(jt2,JM)
!         --- Get max terrain and w using temporary values above, i.e.,
!         --- max of (i,j),(i+1,j),(i,j+1),(i+1,j+1)
          htmax=0.
          wmax=0.
          do jj=jt1,jt2
          do iii=it1,it2
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
             htmax=MAX(mwtd(ii,jj,19),htmax)
             wmax=MAX(mwtd(ii,jj,18),wmax)
          enddo
          enddo
       endif
       mwtd(i,j,1)=htmax
       mwtd(i,j,8)=wmax
    enddo
    enddo

!   --- GTG doens't use wave drag
!   --- wave drag=integral p*dhdx
!     do i=immin,immax
!     do j=jmmin,jmmax
    do j=jsta_m2,jend_m2
    do i=1,IM
       dragx=0.
!       --- At each (i,j) point within the MWT region compute the
!       --- local dragx based on the integral of pdh/dx from 3 points
!       --- to the left to 3 points to the right of the (i,j) pt.  Also
!       --- average over 3 y points.
       im3=i-3
       ip3=i+3
       jp1=MAX(j-1,1+1)    ! post is north-south, original GTG is south-north
       jm1=MIN(j+1,JM-1) ! post is north-south, original GTG is south-north
       aream=0.
       do jj=jp1,jm1,-1  ! post is north-south, original GTG is south-north
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

             dxm=dx(ii,j)/msfx(ii,j)
             ip1=ii+1
             im1=ii-1
             if(im1<1) then
                im1=1
                if(modelname == 'GFS' .or. global) im1=im1+IM
             end if
             if(ip1>IM) then
                ip1=im
                if(modelname == 'GFS' .or. global) ip1=ip1-IM
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

!   --- Get MWT diagnostic multiplier
    do j=jsta,jend
    do i=1,IM
       mws(i,j)=SPVAL
       mwf=mwfilt(i,j)
       if(ABS(mwf-SPVAL)<=SMALL1) cycle
       sms=mwtd(i,j,4)  ! speedmaxt
       hms=mwtd(i,j,1)  ! hmaxt
       UNms=mwtd(i,j,5) ! UNmaxt
       hms=MIN(hms,2750.)  ! 9,000 ft
       ums=mwtd(i,j,2)  ! umaxt
       wms=mwtd(i,j,8)  ! wmaxt
       if(ABS(hms-SPVAL)<=SMALL1) cycle
       if(mwtmultflag==1) then ! 1=speed, 2=w
          if(ABS(sms-SPVAL)<=SMALL1) cycle
          cm=MAX(mwf*sms*hms,0.)
       else
          if(ABS(wms-SPVAL)<=SMALL1) cycle
          cm=MAX(mwf*wms*hms,0.)
       endif
       mws(i,j)=cm   ! Multiplier for MWT diagnostics
    enddo
    enddo

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
