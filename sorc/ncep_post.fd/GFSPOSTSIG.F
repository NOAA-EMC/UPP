! Add Iredells subroutine to read sigma files
!-------------------------------------------------------------------------------
subroutine rtsig(lusig,head,k1,k2,kgds,ijo,levs,ntrac,jcap,lnt2,me,     &
                 h,p,px,py,t,u,v,d,trc,iret)
!$$$  Subprogram documentation block
!
! Subprogram: rtsig      Read and transform sigma file
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram reads a sigma file and transforms
!   the fields to a designated global grid.
!
! Program history log:
!   1999-10-18  Mark Iredell
!   2013-04-19  Jun  Wang:  add option to get tmp and ps(in pascal) 
!                           from enthalpy and ps(cb) option
!   2013-05-06 Shrinivas Moorthi: Initialize midea to 0
!   2013-05-07 Shrinivas Moorthi: Remove mo3, mct, midea and define io3, ict etc
!                                 correctly and get correct cloud condensate.
!   2013-08-02 Shrinivas Moorthi: Rewrote the whole routine to read the sigma
!                                 file differently and to read all tracers
!                                 Addedd sptezj for two 2d fields
!   2014-02-20 Shrinivas Moorthi: Modified conversion from spectral to grid
!                                 taking advantage of threding in SP library.
!                                 This really speeds up the code
!                                 Also threaded loop for Temperature from Tv

!
! Usage:  call rtsig(lusig,head,k1,k2,kgds,ijo,nct,                     &
!                    h,p,px,py,t,tx,ty,u,v,d,z,sh,o3,ct,iret,o,o2)
!   Input argument list:
!     lusig    integer(sigio_intkind) sigma file unit number
!     head     type(sigio_head) sigma file header
!     k1       integer first model level to return
!     k2       integer last model level to return
!     kgds     integer (200) GDS to which to transform
!     ijo      integer dimension of output fields
!     levs     integer number of total vertical levels
!     ntrac    integer number of output tracers
!     jcap     integer number of waves
!     lnt2     integer (jcap+1)*(jcap+2)
!   Output argument list:
!     h        real (ijo) surface orography (m)
!     p        real (ijo) surface pressure (Pa)
!     px       real (ijo) log surface pressure x-gradient (1/m)
!     py       real (ijo) log surface pressure y-gradient (1/m)
!     t        real (ijo,k1:k2) temperature (K)
!     tx       real (ijo,k1:k2) virtual temperature x-gradient (K/m)
!     ty       real (ijo,k1:k2) virtual temperature y-gradient (K/m)
!     u        real (ijo,k1:k2) x-component wind (m/s)
!     v        real (ijo,k1:k2) y-component wind (m/s)
!     d        real (ijo,k1:k2) wind divergence (1/s)
!     trc      real (ijo,k1:k2,ntrac) tracers
!                                1 = specific humidity (kg/kg)
!                                2 = Ozone mixing ratio (kg/kg)
!                                3 = cloud condensate mixing ratio (kg/kg)
!                                .
!                                .
!                                    atomic oxyge, oxygen etc
!
!     iret     integer return code
!
! Modules used:
!   sigio_r_module sigma file I/O
!
! Subprograms called:
!   sigio_rrdati   read sigma single data field
!   sptez          scalar spectral transform
!   sptezd         gradient spectral transform
!   sptezm         multiple scalar spectral transform
!   sptezmv        multiple vector spectral transform
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use sigio_module,   only : sigio_intkind, sigio_head
  use sigio_r_module, only : sigio_dati, sigio_rrdati
  use physcons_post,       only : con_omega, con_fvirt
  use omp_lib
  implicit none
  integer(sigio_intkind),intent(in)    :: lusig
  type(sigio_head),      intent(in)    :: head
  integer,intent(in)                   :: k1,k2,kgds(200),ijo,levs,ntrac,jcap,lnt2,me
  real,dimension(ijo),   intent(out)   :: h,p,px,py
  real,dimension(ijo,k1:k2),intent(out):: t,u,v,d
  real,dimension(ijo,k1:k2,ntrac),intent(out),target :: trc
  integer,intent(out) :: iret
!
  integer idrt,io,jo,iidea
! integer idrt,io,jo,mo3,mct,iidea,midea
  integer(sigio_intkind):: irets
! type(sigio_datm):: datm
  type(sigio_dati) dati
! type griddata
! real,dimension(:,:),pointer :: datm
! endtype griddata
! type(griddata),dimension(:),pointer :: datatrc
  real, target ::  trisca(lnt2,k1:k2+1), triscb(lnt2,k1:k2)
  real,dimension(:),  allocatable :: cpi
  real,dimension(:,:),allocatable :: wrk
  integer io3,ict,jct,n,i,k,jc,nt
  integer idvm, klen
  real    pmean,sumq,xcp
! integer, parameter :: latch=20
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Determine output grid
  idrt = kgds(1)
  if(kgds(1) == 0 .and. kgds(4) < 90000) idrt = 256
  io = kgds(2)
  jo = kgds(3)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Read and transform surface fields
  iret = 1
  if (me == 0) then
    print*,'Debug rtsig= ',lusig,k1,k2,ijo,kgds(1:20)
  endif

  idvm  = head%idvm
! jc = omp_get_num_threads()
! write(0,*)' in RTSIG lnt2=',lnt2,' threads=',jc,' latch=',latch,   &
!          ' jcap=',jcap,' io=',io,' jo=',jo,' ijo=',ijo
!
  if (k2 < k1) return

  dati%i = 1                                           ! hs
  dati%f => trisca(:,k1)
  call sigio_rrdati(lusig,head,dati,irets)
  if(irets /= 0) return

! call sptez(0,jcap,idrt,io,jo,trisca(1,k1),h,1)
! call sptez(0,jcap,idrt,io,jo,dats%hs,h,1)
! call sptez(0,jcap,idrt,io,jo,dats%ps,p,1)
! call sptezj(jcap,lnt2,1,idrt,io,jo,jc,trisca,h,latch,1)
!
  dati%i = 2                               ! Surface pressure
  dati%f => trisca(:,k1+1)
  call sigio_rrdati(lusig,head,dati,irets)
  if(irets /= 0) return
!
! call sptez(0,jcap,idrt,io,jo,trisca(1,k1),p,1)
! call sptezj(jcap,lnt2,1,idrt,io,jo,jc,trisca,p,latch,1)
!--
  allocate(wrk(ijo,2))
  call sptezm(0,jcap,idrt,io,jo,2,trisca(1,k1),wrk,1)
  if( mod(idvm,10) < 2) then
!$omp parallel do private(i)
    do i=1,ijo
      h(i) = wrk(i,1)
      p(i) = 1.e3*exp(wrk(i,2))
!     p(i) = 1.e3*exp(p(i))
    enddo
  elseif(mod(idvm,10) == 2) then
!$omp parallel do private(i)
    do i=1,ijo
      h(i) = wrk(i,1)
      p(i) = 1000.*wrk(i,2)
!     p(i) = 1000.*p(i)
    enddo
  endif
  if (allocated(wrk)) deallocate(wrk)

  call sptezd(0,jcap,idrt,io,jo,trisca(1,k1+1),pmean,px,py,1)
  iret = 0

! if (k2 < k1) return

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Read and transform fields on levels k1 through k2
  iret = 2
  if (k2 >= k1) then
    klen = k2-k1+1
    do k=k1,k2
      write(0,*)' retriving T for k=',k,' k1=',k1,' k2=',k2
      dati%i = k + 2                        ! Virtual Temperature or CpT
      dati%f => trisca(:,k)

      call sigio_rrdati(lusig,head,dati,iret)
    enddo
    call sptezm(0,jcap,idrt,io,jo,klen,trisca(1,k1),t(1,k1),1)
!   call sptezm(0,jcap,idrt,io,jo,klen,trisca,t,1)
    do k=k1,k2
      dati%i = levs + 2 + (k-1) * 2 + 1     ! Divergence
      dati%f => trisca(:,k)
      call sigio_rrdati(lusig,head,dati,irets)
      if(irets /= 0) return
      dati%i = levs + 2 + (k-1) * 2 + 2     ! Vorticity
      dati%f => triscb(:,k)
      call sigio_rrdati(lusig,head,dati,irets)
      if(irets /= 0) return
    enddo
    call sptezmv(0,jcap,idrt,io,jo,klen,trisca(1,k1),triscb(1,k1),  &
                 u(1,k1),v(1,k1),1)
    call sptezm(0,jcap,idrt,io,jo,klen,trisca(1,k1),d(1,k1),1)

!   call sptezm(0,jcap,idrt,io,jo,1,triscb,z(1,k),1)
    write(0,*)' retriving d/z for k=',k,' k1=',k1,' k2=',k2
!   datm%z(3,:) = datm%z(3,:)+2*con_omega/sqrt(1.5)
!   call sptezm(0,jcap,idrt,io,jo,klen,datm%z,z,1)
   write(0,*)' start get tracer'
    do nt=1,ntrac
      do k=k1,k2
        dati%i = levs * (2+nt) + 2 + k      ! Tracers starting with q
        dati%f => trisca(:,k)
        call sigio_rrdati(lusig,head,dati,irets)
      enddo
      call sptezm(0,jcap,idrt,io,jo,klen,trisca(1,k1),trc(1,k1,nt),1)
      write(0,*)' retriving d/z for nt=',nt,'ntrac=',ntrac,'k=',k,' k1=',k1,' k2=',k2
    enddo
    !t=t/(1+con_fvirt*sh)
   write(0,*)' end get tracer,idvm=',idvm,'ijo=',ijo,'ntrac=',ntrac
!
!-- get temp 
    if (mod(idvm/10,10) == 3) then ! Enthalpy case
      allocate(cpi(0:ntrac))
!     write(0,*)'aft read sig, cpi=',head%cpi
      cpi(0:ntrac) = head%cpi(1:ntrac+1)
!     write(0,*)'cpi=',cpi(0:ntrac)
!$omp parallel do private(k,i,xcp,sumq,n)
      do k=k1,k2
        do i=1,ijo
          xcp  = 0.0
          sumq = 0.0
          do n=1,ntrac
            if( cpi(n) /= 0.0 ) then
              xcp  = xcp  + cpi(n)*trc(i,k,n)
              sumq = sumq + trc(i,k,n)
            endif
          enddo
          xcp    = (1.-sumq)*cpi(0) + xcp
          t(i,k) = t(i,k) / xcp   ! Now g1 contains T
        enddo
      enddo
      if (allocated(cpi)) deallocate(cpi)
    else
!$omp parallel do private(i,k)
      do k=k1,k2
        do i=1,ijo
          t(i,k) = t(i,k) / (1+con_fvirt*trc(i,k,1)) !get temp from virtual temp
        enddo
      enddo
    endif
  endif
! write(0,*)'end comput t'
  iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine  

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine modstuff(km,idvc,idsl,nvcoord,vcoord,ps,psx,psy,d,u,v,&
                      pi,pm,om)
!$$$  Subprogram documentation block
!
! Subprogram: modstuff   Compute model coordinate dependent functions
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes fields which depend on the model coordinate
!           such as pressure thickness and vertical velocity.
!
! Program history log:
!   1999-10-18  Mark Iredell
!   2013-04-19  Jun  Wang:  add option to get pi by using 8byte real computation
!
! Usage:  call modstuff(km,idvc,idsl,nvcoord,vcoord,ps,psx,psy,d,u,v,&
!                       pd,pi,pm,os,om,px,py)
!   Input argument list:
!     km       integer number of levels
!     idvc     integer vertical coordinate id (1 for sigma and 2 for hybrid)
!     idsl     integer type of sigma structure (1 for phillips or 2 for mean)
!     nvcoord  integer number of vertical coordinates
!     vcoord   real (km+1,nvcoord) vertical coordinates
!     ps       real surface pressure (Pa)
!     psx      real log surface pressure x-gradient (1/m)
!     psy      real log surface pressure y-gradient (1/m)
!     d        real (km) wind divergence (1/s)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!   Output argument list:
!     pi       real (km+1) interface pressure (Pa)
!     pm       real (km) mid-layer pressure (Pa)
!     om       real (km) vertical velocity (Pa/s)
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use sigio_module, only: sigio_modprd
    implicit none
    integer,intent(in):: km,idvc,idsl,nvcoord
    real,intent(in):: vcoord(km+1,nvcoord)
    real,intent(in):: ps,psx,psy
    real,intent(in):: u(km),v(km),d(km)
!   real,intent(out):: pi(km+1),pm(km)
    real*8, intent(out):: pi(km+1),pm(km)
    real,intent(out):: om(km)
    real*8 ps8,pm8(km),pd8(km),vcoord8(km+1,nvcoord)
    real*8 dpmdps(km),dpddps(km),dpidps(km+1),pi8(km+1)
    real vgradp,pd(km),px(km),py(km),os
    integer k,iret,logk
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ps8=ps
    vcoord8=vcoord
    call sigio_modprd(1,1,km,nvcoord,idvc,idsl,vcoord8,iret,&
                     ps=(/ps8/),&
                     pm=pm8,pd=pd8,dpmdps=dpmdps,dpddps=dpddps)
!
!jw: has to be 8 real for wam     
    pi8(1)=ps
    pm=pm8
!   pd=pd8
    dpidps(1)=1.
    do k=1,km
      pi8(k+1)=pi8(k)-pd8(k)
      dpidps(k+1)=dpidps(k)-dpddps(k)
!     if(pi(8)<0.) then
!        print *,'in modstuff,pi8=',pi8(k)
!     endif
    enddo
    pi=pi8
!
    os=0
    do k=km,1,-1
      vgradp=u(k)*psx+v(k)*psy
      os=os-vgradp*ps*(dpmdps(k)-dpidps(k+1))-d(k)*(pm(k)-pi(k+1))
      om(k)=vgradp*ps*dpmdps(k)+os
      os=os-vgradp*ps*(dpidps(k)-dpmdps(k))-d(k)*(pi(k)-pm(k))
    enddo
    px=ps*dpmdps*psx
    py=ps*dpmdps*psy
  end subroutine

!-------------------------------------------------------------------------------
  subroutine modstuff2(im,ix,km,idvc,idsl,nvcoord,vcoord,ps,psx,psy,d,u,v,&
                       pi,pm,om,me)
!$$$  Subprogram documentation block
!
! Subprogram: modstuff   Compute model coordinate dependent functions
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes fields which depend on the model coordinate
!           such as pressure thickness and vertical velocity.
!
! Program history log:
!   1999-10-18  Mark Iredell
!   2013-04-19  Jun  Wang:  add option to get pi by using 8byte real computation
!   2013-08-13  Shrinivas Moorthi - Modified to include im points and thread
!
! Usage:  call modstuff(km,idvc,idsl,nvcoord,vcoord,ps,psx,psy,d,u,v,&
!                       pd,pi,pm,os,om,px,py)
!   Input argument list:
!     im       integer - inner computational domain
!     ix       integer - maximum inner dimension
!     km       integer number of levels
!     idvc     integer vertical coordinate id (1 for sigma and 2 for hybrid)
!     idsl     integer type of sigma structure (1 for phillips or 2 for mean)
!     nvcoord  integer number of vertical coordinates
!     vcoord   real (km+1,nvcoord) vertical coordinates
!     ps       real surface pressure (Pa)
!     psx      real log surface pressure x-gradient (1/m)
!     psy      real log surface pressure y-gradient (1/m)
!     d        real (km) wind divergence (1/s)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!   Output argument list:
!     pi       real (km+1) interface pressure (Pa)
!     pm       real (km) mid-layer pressure (Pa)
!     om       real (km) vertical velocity (Pa/s)
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use sigio_module, only : sigio_modprd
    implicit none
    integer,                    intent(in)  :: im,ix,km,idvc,idsl,nvcoord,me
    real,                       intent(in)  :: vcoord(km+1,nvcoord)
    real,   dimension(ix),      intent(in)  :: ps,psx,psy
    real,   dimension(ix,km),   intent(in)  :: u,v,d
    real*8, dimension(ix,km+1), intent(out) :: pi
    real*8, dimension(ix,km),   intent(out) :: pm
    real,   dimension(ix,km),   intent(out) :: om
!   real*8, allocatable :: ps8(:), pm8(:,:), pd8(:,:),dpmdps(:,:),dpddps(:,:), &
!                          dpidps(:,:),pi8(:,:),vcoord8(:,:)
!   real,   allocatable :: os(:)
!   real,   allocatable :: pd(:,:),px(:,:), py(:,:), os(:)

!   real vgradpps

    real*8 ps8(ix),pm8(ix,km),pd8(ix,km),vcoord8(km+1,nvcoord)
    real*8 dpmdps(ix,km),dpddps(ix,km),dpidps(ix,km+1),pi8(ix,km+1)
    real vgradpps,pd(im,km),os(im)
!   real vgradpps,pd(im,km),px(im,km),py(im,km),os(im),tem
    integer i,k,iret,logk
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ps8     = ps
    vcoord8 = vcoord
    call sigio_modprd(im,ix,km,nvcoord,idvc,idsl,vcoord8,iret,                 &
                     ps=ps8,pd=pd8,dpddps=dpddps,pm=pm8,dpmdps=dpmdps)

!
!   if (me == 0) then
!     write(0,*)' pd8=',pd8(1,60:64)
!     write(0,*)' pm8=',pm8(1,60:64)
!    endif
!jw: has to be 8 real for wam     

!$omp parallel do private(i)
    do i=1,im
      pi8(i,1)    = ps(i)
      dpidps(i,1) = 1.
      os(i)       = 0
      pi(i,1)     = pi8(i,1)
    enddo
    do k=1,km
!$omp parallel do private(i)
      do i=1,im
        pi8(i,k+1)    = pi8(i,k)    - pd8(i,k)
        dpidps(i,k+1) = dpidps(i,k) - dpddps(i,k)
!       if(pi(i,8)<0.) then
!          print *,'in modstuff,pi8=',pi8(i,k),' i=',i,' k=',k,' me=',me
!       endif
        pi(i,k+1) = pi8(i,k+1)
        pm(i,k)   = pm8(i,k)
      enddo
    enddo
!
    do k=km,1,-1
!$omp parallel do private(i,vgradpps)
      do i=1,im
        vgradpps = (u(i,k)*psx(i) + v(i,k)*psy(i)) * ps(i)

        os(i)    = os(i) - vgradpps*(dpmdps(i,k)-dpidps(i,k+1))             &
                         - d(i,k)*(pm(i,k)-pi(i,k+1))

        om(i,k)  = os(i) + vgradpps*dpmdps(i,k)

        os(i)    = os(i) - vgradpps*(dpidps(i,k)-dpmdps(i,k))               &
                         - d(i,k)*(pi(i,k)-pm(i,k))
      enddo
    enddo
  end subroutine

!-----------------------------------------------------------------------
      subroutine trssc(jcap,nc,km,ntrac,idvc,idvm,idsl,nvcoord,vcoord,     &
                       cpi,idrt,lonb,latb,ijl,ijn,j1,j2,jc,chgq0,          &
                       szs,sps,st,sd,sz,sq,gfszs,gfsps,gfsp,gfsdp,         &
                       gfst,gfsu,gfsv,gfsq,gfsw)
!$$$  subprogram documentation block
!
! subprogram:    trssc       transform sigma spectral fields to grid
!   prgmmr: iredell          org: w/nmc23     date: 92-10-31
!
! abstract: transforms sigma spectral fields to grid and converts
!   log surface pressure to surface pressure and virtual temperature
!   to temperature.
!
! program history log:
!   91-10-31  mark iredell
!
! usage:    call trssc(jcap,nc,km,ntrac,idvm,
!    &                 idrt,lonb,latb,ijl,j1,j2,jc,
!    &                 szs,sps,st,sd,sz,sq,zs,ps,t,u,v,q)
!   input argument list:
!     jcap         integer spectral truncation
!     nc           integer first dimension (nc>=(jcap+1)*(jcap+2))
!     km           integer number of levels
!     ntrac        integer number of tracers
!     idvm         integer mass variable id
!     idrt         integer data representation type
!     lonb         integer number of longitudes
!     latb         integer number of latitudes
!     ijl          integer horizontal dimension
!     j1           integer first latitude
!     j2           integer last latitude
!     jc           integer number of cpus
!     szs          real (nc) orography
!     sps          real (nc) log surface pressure
!     st           real (nc,levs) virtual temperature
!     sd           real (nc,levs) divergence
!     sz           real (nc,levs) vorticity
!     sq           real (nc,levs*ntrac) tracers
!   output argument list:
!     zs           real (ijl) orography
!     ps           real (ijl) surface pressure
!     t            real (ijl,km) temperature
!     u            real (ijl,km) zonal wind
!     v            real (ijl,km) meridional wind
!     q            real (ijl,km*ntrac) tracers
!
! subprograms called:
!   sptran       perform a scalar spherical transform
!
! attributes:
!   language: fortran
!
!c$$$
      use gfsio_module
!      use gfsio_rst
      implicit none
      integer,intent(in)::jcap,nc,km,ntrac,idvc,idvm,idsl,nvcoord,idrt,lonb,latb
      integer,intent(in)::ijl,ijn,j1,j2,jc,chgq0
      real,intent(in):: szs(nc),sps(nc),st(nc,km),sd(nc,km),sz(nc,km),sq(nc,km*ntrac)
      real,intent(in):: cpi(0:ntrac)
      real*8,intent(in):: vcoord(km+1,nvcoord)
      real,dimension(ijn),intent(inout):: gfszs,gfsps
      real,dimension(ijn,km),intent(inout):: gfsp,gfsdp,gfst,gfsu,gfsv,gfsw
      real,dimension(ijn,km*ntrac),intent(inout):: gfsq
      real zs(ijl),ps(ijl),t(ijl,km),u(ijl,km),v(ijl,km),q(ijl,km*ntrac)
      real wi(ijl,km),pi(ijl,km),dpo(ijl,km)
      real tvcon,xcp,sumq
      integer thermodyn_id,jn,js,is,in
      integer jj,jjm,k,n,j,i,ij,lonb2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  spectral transforms
      if(j1==732)print*,'sample input to trssc= ',jcap,nc,km,ntrac, &
                          idvc,idvm,idsl,nvcoord,    &
                       idrt,lonb,latb,ijl,ijn,j1,j2,jc,chgq0
      lonb2=lonb*2
      ij=lonb2*(j2-j1+1)
      in=1
      is=1+lonb
      call sptran(0,jcap,idrt,lonb,latb,1,1,1,lonb2,lonb2,nc,ijl,           &
                  j1,j2,jc,szs,zs(in),zs(is),1)
      call sptran(0,jcap,idrt,lonb,latb,1,1,1,lonb2,lonb2,nc,ijl,           &
                  j1,j2,jc,sps,ps(in),ps(is),1)
      call sptran(0,jcap,idrt,lonb,latb,km,1,1,lonb2,lonb2,nc,ijl,          &
                  j1,j2,jc,st,t(in,1),t(is,1),1)
      call sptranv(0,jcap,idrt,lonb,latb,km,1,1,lonb2,lonb2,nc,ijl,         &
                   j1,j2,jc,sd,sz,u(in,1),u(is,1),v(in,1),v(is,1),1)
      call sptran(0,jcap,idrt,lonb,latb,km*ntrac,1,1,lonb2,lonb2,nc,ijl,    &
                  j1,j2,jc,sq,q(in,1),q(is,1),1)
      if(j1==732)then
       do k=1,km
        do i=1,ijl
	 if(t(i,k)>400. .or. t(i,k)<100.)print*,'bad T from sptran',i,k,t(i,k)
	 if(q(i,k)>1.)print*,'bad Q  from sptran',i,k,q(i,k)
	 if(q(i,2*k)>1.)print*,'bad Q  from sptran',i,k,q(i,2*k)
	 if(q(i,3*k)>1.)print*,'bad Q  from sptran',i,k,q(i,3*k)
	end do
       end do	 
      end if  
      select case(mod(idvm,10))
      case(0,1)
        do i=1,ij
          ps(i)=1.e3*exp(ps(i))
        enddo
      case(2)
        do i=1,ij
          ps(i)=1.e3*ps(i)
        enddo
      case default
        do i=1,ij
          ps(i)=1.e3*exp(ps(i))
        enddo
      end select
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      thermodyn_id=mod(idvm/10,10)
      if (thermodyn_id == 3) then
        do k=1,km
         do i=1,ij
            t(i,k) = t(i,k)/cpi(0)   ! enthalpy (cpt/cpd)
         end do
        end do
!
      endif

      call getomega(jcap,nc,km,idvc,idvm,idrt,idsl,                    &
        nvcoord,vcoord,lonb,latb,ijl,j1,j2,1,sd,sps,                   &
        ps,t,u,v,wi,pi,dpo)
      if(j1==732)then
       do k=1,km
        do i=1,ijl
	 if(t(i,k)>400. .or. t(i,k)<100.)print*,'bad T after getomega',i,k,t(i,k)
	 if(q(i,k)>1. )print*,'bad Q  after getomega',i,k,q(i,k)
	 if(q(i,2*k)>1. )print*,'bad Q  after getomega',i,2*k,q(i,2*k)
	end do
       end do	 
      end if  
      if(thermodyn_id /= 2)then
!  convert to surface pressure and temperature
         if (thermodyn_id == 3) then
           do k=1,km
            do i=1,ij
               xcp  = 0.0
               sumq = 0.0
               do n=1,ntrac
                 if( cpi(n) .ne. 0.0 ) then
                   xcp  = xcp  + cpi(n)*q(i,k+(n-1)*km)
                   sumq = sumq + q(i,k+(n-1)*km)
                 endif
               enddo
               t(i,k)  = t(i,k)/((1.-sumq)*cpi(0)+xcp)
            end do
           end do

          else 
           tvcon=(461.50/287.05-1.)
           t(:,:) = t(:,:)/(1.+tvcon*q(:,1:km))
          endif
      end if
      if(j1==732)then
       do k=1,km
        do i=1,ijl
	 if(t(i,k)>400. .or. t(i,k)<100.)print*,'bad T after Tv to T',i,k,t(i,k)
	 if(q(i,k)>1.)print*,'bad Q  after Tv to T',i,k,q(i,k)
	 if(q(i,2*k)>1. )print*,'bad Q  after Tv to T',i,k,q(i,2*k)
	end do
       end do	 
      end if  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!----force tracers to be positive
      if (chgq0 == 1) q = max(q, 0.0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  pass data to gfsdatao
      do j=1,j2-j1+1
        jn=j+j1-1
        js=latb+1-jn
        jn=(jn-1)*lonb
        js=(js-1)*lonb
        jj=j*lonb
        jjm=(j-1)*lonb
        do i=1,lonb
          gfszs(i+jn) = zs(i+jjm)
          gfsps(i+jn) = ps(i+jjm)
        enddo
        do i=1,lonb
          gfszs(i+js) = zs(i+jj)
          gfsps(i+js) = ps(i+jj)
        enddo
        do k=1,km
         do i=1,lonb
          gfsdp(i+jn,k) = dpo(i+jjm,k)
          gfsp(i+jn,k)  = pi(i+jjm,k)
          gfst(i+jn,k)  = t(i+jjm,k)
          gfsu(i+jn,k)  = u(i+jjm,k)
          gfsv(i+jn,k)  = v(i+jjm,k)
          gfsw(i+jn,k)  = wi(i+jjm,k)
         enddo
         do i=1,lonb
          gfsdp(i+js,k) = dpo(i+jj,k)
          gfsp(i+js,k)  = pi(i+jj,k)
          gfst(i+js,k)  = t(i+jj,k)
          gfsu(i+js,k)  = u(i+jj,k)
          gfsv(i+js,k)  = v(i+jj,k)
          gfsw(i+js,k)  = wi(i+jj,k)
         enddo
        enddo
        do k=1,km*ntrac
          do i=1,lonb
            gfsq(i+jn,k) = q(i+jjm,k)
          enddo
          do i=1,lonb
            gfsq(i+js,k) = q(i+jj,k)
          enddo
        enddo
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine getomega(jcap,nc,km,idvc,idvm,idrt,idsl,nvcoord,vcoord,  &
            lonb,latb,ijn,j1,j2,jc,sd,sps,psi,ti,ui,vi,wi,pm,pd)
!!!!!
      use sigio_module, only : sigio_modprd
      implicit none

      integer,intent(in):: jcap,nc,km,idvc,idvm,idrt,idsl,nvcoord
      integer,intent(in):: lonb,latb,j1,j2,jc,ijn
      real*8,intent(in):: vcoord(km+1,nvcoord)
      real,intent(in):: sd(nc,km),sps(nc)
      real,intent(in):: psi(ijn),ti(ijn,km),ui(ijn,km),vi(ijn,km)
      real,intent(out):: wi(ijn,km),pm(ijn,km),pd(ijn,km)
      real :: pi(ijn,km+1)
      real :: os
      real*8 psi8(ijn),ti8(ijn,km),pm8(ijn,km),pd8(ijn,km)
      real*8 dpmdps(ijn,km),dpddps(ijn,km),dpidps(ijn,km+1),vgradp,psmean
      real di(ijn,km),psx(ijn),psy(ijn)
      integer k,i,ij,lonb2,iret,is,in
!----1. spectral transform
      lonb2=lonb*2
      ij=lonb2*(j2-j1+1)
      in=1
      is=1+lonb
      call sptrand(0,jcap,idrt,lonb,latb,1,1,1,lonb2,lonb2,nc,ijn,    &
           j1,j2,jc,sps,psmean,psx(in),psx(is),psy(in),psy(is),1)

      call sptran(0,jcap,idrt,lonb,latb,km,1,1,lonb2,lonb2,nc,ijn,    &
                  j1,j2,jc,sd,di(in,1),di(is,1),1)
      psi8=psi	
      ti8=ti  

      call sigio_modprd(ijn,ijn,km,nvcoord,idvc,idsl,vcoord,iret,     &
                   ps=psi8,t=ti8,pm=pm8,pd=pd8,dpmdps=dpmdps,dpddps=dpddps)
      pm=pm8
      pd=pd8		   

      select case(mod(idvm,10))
      case(0,1)
          continue
      case(2)
          do i=1,ijn
           psx(i)=psx(i)/(psi(i)*1.0e-3)
           psy(i)=psy(i)/(psi(i)*1.0e-3)
          enddo
      case default
          do i=1,ijn
           psx(i)=psx(i)/psi(i)
           psy(i)=psy(i)/psi(i)
          enddo
      end select

!----3.omeda from modstuff
      do i=1,ijn
       pi(i,1)=psi(i)
       dpidps(i,1)=1.
      enddo
      do k=1,km
       do i=1,ijn
         pi(i,k+1)=pi(i,k)-pd(i,k)
         dpidps(i,k+1)=dpidps(i,k)-dpddps(i,k)
       enddo
      enddo
      do i=1,ijn
       os=0.
       do k=km,1,-1
        vgradp=ui(i,k)*psx(i)+vi(i,k)*psy(i)
        os=os-vgradp*psi(i)*(dpmdps(i,k)-dpidps(i,k+1))-                 &
           di(i,k)*(pm(i,k)-pi(i,k+1))
        wi(i,k)=vgradp*psi(i)*dpmdps(i,k)+os
        os=os-vgradp*psi(i)*(dpidps(i,k)-dpmdps(i,k))-                   &
           di(i,k)*(pi(i,k)-pm(i,k))
       enddo
!
      enddo
!---
       return
       end subroutine 
