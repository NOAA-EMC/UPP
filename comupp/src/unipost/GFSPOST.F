  subroutine pvetc(km,p,px,py,t,tx,ty,h,u,v,av,hm,s,bvf2,pvn,theta,sigma,pvu)
!$$$  Subprogram documentation block
!
! Subprogram: pvetc      Compute potential vorticity, etc
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes
!             the Montgomery streamfunction
!               hm=cp*t+g*z
!             the specific entropy
!               s=cp*log(t/t0)-r*log(p/p0)
!             the Brunt-Vaisala frequency squared
!               bvf2=g/cp*ds/dz
!             the potential vorticity defined as
!               pvn=(av*ds/dz-dv/dz*ds/dx+du/dz*ds/dy)/rho/cp
!             the potential temperature
!               theta=t0*exp(s/cp)
!             the static stability
!               sigma=t/g*bvf2
!             and the potential vorticity in PV units
!               pvu=10**-6*theta*pvn
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call pvetc(km,p,px,py,t,tx,ty,h,u,v,av,s,bvf2,pvn,theta,sigma,pvu)
!   Input argument list:
!     km       integer number of levels
!     p        real (km) pressure (Pa)
!     px       real (km) pressure x-gradient (Pa/m)
!     py       real (km) pressure y-gradient (Pa/m)
!     t        real (km) (virtual) temperature (K)
!     tx       real (km) (virtual) temperature x-gradient (K/m)
!     ty       real (km) (virtual) temperature y-gradient (K/m)
!     h        real (km) height (m)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     av       real (km) absolute vorticity (1/s)
!   Output argument list:
!     hm       real (km) Montgomery streamfunction (m**2/s**2)
!     s        real (km) specific entropy (J/K/kg)
!     bvf2     real (km) Brunt-Vaisala frequency squared (1/s**2)
!     pvn      real (km) potential vorticity (m**2/kg/s)
!     theta    real (km) (virtual) potential temperature (K)
!     sigma    real (km) static stability (K/m)
!     pvu      real (km) potential vorticity (10**-6*K*m**2/kg/s)
!
! Modules used:
!   physcons       Physical constants
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons
    implicit none
    integer,intent(in):: km
    real,intent(in),dimension(km):: p,px,py,t,tx,ty,h,u,v,av
    real,intent(out),dimension(km):: hm,s,bvf2,pvn,theta,sigma,pvu
!   real,parameter:: hhmin=500.,t0=2.e2,p0=1.e5
    real,parameter:: hhmin=5.,t0=2.e2,p0=1.e5
    integer k,kd,ku,k2(2)
    real cprho,sx,sy,sz,uz,vz
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    do k=1,km
      hm(k)=con_cp*t(k)+con_g*h(k)
      s(k)=con_cp*log(t(k)/t0)-con_rd*log(p(k)/p0)
    enddo
    do k=1,km

#ifdef ibm_sp
      call rsearch1(km,h,2,(/h(k)-hhmin,h(k)+hhmin/),k2)
#else
      call rsearch(1,km,1,1,h,2,1,1,(/h(k)-hhmin,h(k)+hhmin/),1,1,k2)
#endif
!      kd=max(k2(1),1)
!      ku=min(k2(2)+1,km)
!      kd=min(k2(1),km) ! Chuang: post counts from top down, redefine lower bound
      kd=min(k2(1)+1,km) ! Chuang: post counts from top down,
!      ku=max(k2(2)-1,1)
      ku=max(k2(2),1)
      if(ku==1)kd=2 ! Chuang: make sure ku ne kd at model top
      cprho=p(k)/(con_rocp*t(k))
      sx=con_cp*tx(k)/t(k)-con_rd*px(k)/p(k)
      sy=con_cp*ty(k)/t(k)-con_rd*py(k)/p(k)
      sz=(s(ku)-s(kd))/(h(ku)-h(kd))
      uz=(u(ku)-u(kd))/(h(ku)-h(kd))
      vz=(v(ku)-v(kd))/(h(ku)-h(kd))
      bvf2(k)=con_g/con_cp*sz
      pvn(k)=(av(k)*sz-vz*sx+uz*sy)/cprho
      theta(k)=t0*exp(s(k)/con_cp)
      sigma(k)=t(k)/con_g*bvf2(k)
      pvu(k)=1.e6*theta(k)*pvn(k)
    enddo
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine p2th(km,theta,u,v,h,t,pvu,sigma,rh,omga,kth,th &
,lth,uth,vth,hth,tth,zth,sigmath,rhth,oth)
!$$$  Subprogram documentation block
!
! Subprogram: p2th       Interpolate to isentropic level
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram interpolates fields to given isentropic levels.
!   The interpolation is linear in entropy.
!   Outside the domain the bitmap is set to false.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call p2th(km,theta,u,v,h,t,puv,kth,th,uth,vth,tth)
!   Input argument list:
!     km       integer number of levels
!     theta    real (km) potential temperature (K)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     h        real (km) height (m)
!     t        real (km) temperature (K)
!     pvu      real (km) potential vorticity in PV units (10**-6*K*m**2/kg/s)
!     kth      integer number of isentropic levels
!     th       real (kth) isentropic levels (K)
!   Output argument list:
!     lpv      logical*1 (kth) bitmap
!     uth      real (kth) x-component wind (m/s)
!     vth      real (kth) y-component wind (m/s)
!     hth      real (kth) height (m)
!     tth      real (kth) temperature (K)
!     zth      real (kth) potential vorticity in PV units (10**-6*K*m**2/kg/s)
!
! Subprograms called:
!   rsearch       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    implicit none
    integer,intent(in):: km,kth
    real,intent(in),dimension(km):: theta,u,v,h,t,pvu,sigma,rh,omga
    real,intent(in):: th(kth)
    logical*1,intent(out),dimension(kth):: lth
    real,intent(out),dimension(kth):: uth,vth,hth,tth,zth &
    ,sigmath,rhth,oth
    real w
    integer loc(kth),l
    integer k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef ibm_sp
    call rsearch1(km,theta(1),kth,th(1),loc(1))
#else
    call rsearch(1,km,1,1,theta(1),kth,1,1,th(1),1,1,loc(1))
#endif
    do k=1,kth
      l=loc(k)
      lth(k)=l.gt.0.and.l.lt.km
      if(lth(k)) then
        w=log(th(k)/theta(l))/log(theta(l+1)/theta(l))
        uth(k)=u(l)+w*(u(l+1)-u(l))
        vth(k)=v(l)+w*(v(l+1)-v(l))
        hth(k)=h(l)+w*(h(l+1)-h(l))
        tth(k)=t(l)+w*(t(l+1)-t(l))
        zth(k)=pvu(l)+w*(pvu(l+1)-pvu(l))
	sigmath(k)=sigma(l)+w*(sigma(l+1)-sigma(l))
	rhth(k)=rh(l)+w*(rh(l+1)-rh(l))
!	pth(k)=p(l)+w*(p(l+1)-p(l))
	oth(k)=omga(l)+w*(omga(l+1)-omga(l))
      endif
    enddo
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine p2pv(km,pvu,h,t,p,u,v,kpv,pv,pvpt,pvpb,&
                  lpv,upv,vpv,hpv,tpv,ppv,spv)
!$$$  Subprogram documentation block
!
! Subprogram: p2pv       Interpolate to potential vorticity level
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram interpolates fields to given potential vorticity
!   levels within given pressure limits.
!   The output level is the first  encountered from the top pressure limit.
!   If the given potential vorticity level is not found, the outputs are zero
!   and the bitmap is false. The interpolation is linear in potential vorticity.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call p2pv(km,pvu,h,t,p,u,v,kpv,pv,pvpt,pvpb,&
!                   lpv,upv,vpv,hpv,tpv,ppv,spv)
!   Input argument list:
!     km       integer number of levels
!     pvu      real (km) potential vorticity in PV units (10**-6*K*m**2/kg/s)
!     h        real (km) height (m)
!     t        real (km) temperature (K)
!     p        real (km) pressure (Pa)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     kpv      integer number of potential vorticity levels
!     pv       real (kpv) potential vorticity levels (10**-6*K*m**2/kg/s)
!     pvpt     real (kpv) top pressures for PV search (Pa)
!     pvpb     real (kpv) bottom pressures for PV search (Pa)
!   Output argument list:
!     lpv      logical*1 (kpv) bitmap
!     upv      real (kpv) x-component wind (m/s)
!     vpv      real (kpv) y-component wind (m/s)
!     hpv      real (kpv) temperature (K)
!     tpv      real (kpv) temperature (K)
!     ppv      real (kpv) pressure (Pa)
!     spv      real (kpv) wind speed shear (1/s)
!
! Subprograms called:
!   rsearch       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons
    implicit none
    integer,intent(in):: km,kpv
    real,intent(in),dimension(km):: pvu,h,t,p,u,v
    real,intent(in):: pv(kpv),pvpt(kpv),pvpb(kpv)
    logical*1,intent(out),dimension(kpv):: lpv
    real,intent(out),dimension(kpv):: upv,vpv,hpv,tpv,ppv,spv
    real,parameter:: pd=2500.
    real w,spdu,spdd
    integer k,l1,l2,lu,ld,l
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    do k=1,kpv
#ifdef ibm_sp
      call rsearch1(km,p,1,pvpb(k),l1)
      call rsearch1(km,p,1,pvpt(k),l2)
#else
      call rsearch(1,km,1,1,p,1,1,1,pvpb(k),1,1,l1)
      call rsearch(1,km,1,1,p,1,1,1,pvpt(k),1,1,l2)
#endif
!      l1=l1+1
      l=0
      if(pv(k).ge.0.) then
!        do lu=l2-1,l1,-1
!        do lu=l2,l1-1 ! Chuang: post counts top down	
        do lu=l2+2,l1 ! Chuang: post counts top down
!          if(pv(k).lt.pvu(lu+1).and.pv(k).ge.pvu(lu)) then
          if(pv(k).ge.pvu(lu+1).and.pv(k).lt.pvu(lu)) then	  
#ifdef ibm_sp
            call rsearch1(km,p,1,p(lu)+pd,ld)
#else
            call rsearch(1,km,1,1,p,1,1,1,p(lu)+pd,1,1,ld)
#endif
!            if(all(pv(k).ge.pvu(ld:lu-1))) then
	    if(all(pv(k).ge.pvu(lu+1:ld))) then
              l=lu
              exit
            endif
          endif
        enddo
      else
!        do lu=l2-1,l1,-1
!        do lu=l2,l1-1 ! Chuang: post counts top down	
        do lu=l2+2,l1 ! Chuang: post counts top down
!          if(pv(k).gt.pvu(lu+1).and.pv(k).le.pvu(lu)) then
	  if(pv(k).le.pvu(lu+1).and.pv(k).gt.pvu(lu)) then
#ifdef ibm_sp
            call rsearch1(km,p,1,p(lu)+pd,ld)
#else
            call rsearch(1,km,1,1,p,1,1,1,p(lu)+pd,1,1,ld)
#endif
!            if(all(pv(k).le.pvu(ld:lu-1))) then
	    if(all(pv(k).le.pvu(lu+1:ld))) then
              l=lu
              exit
            endif
          endif
        enddo
      endif
      lpv(k)=l.gt.0
      if(lpv(k)) then
        w=(pv(k)-pvu(l))/(pvu(l+1)-pvu(l))
        upv(k)=u(l)+w*(u(l+1)-u(l))
        vpv(k)=v(l)+w*(v(l+1)-v(l))
        hpv(k)=h(l)+w*(h(l+1)-h(l))
        tpv(k)=t(l)+w*(t(l+1)-t(l))
        ppv(k)=p(l)*exp((h(l)-hpv(k))*(1-0.5*(tpv(k)/t(l)-1))/(con_rog*t(l)))
        spdu=sqrt(u(l+1)**2+v(l+1)**2)
        spdd=sqrt(u(l)**2+v(l)**2)
        spv(k)=(spdu-spdd)/(h(l+1)-h(l))
      endif
    enddo
  end subroutine
!-------------------------------------------------------------------------------
#ifdef ibm_sp
!-------------------------------------------------------------------------------
subroutine rsearch1(km1,z1,km2,z2,l2)
!$$$  subprogram documentation block
!
! subprogram:    rsearch1    search for a surrounding real interval
!   prgmmr: iredell    org: w/nmc23     date: 98-05-01
!
! abstract: this subprogram searches a monotonic sequences of real numbers
!   for intervals that surround a given search set of real numbers.
!   the sequences may be monotonic in either direction; the real numbers
!   may be single or double precision.
!
! program history log:
! 1999-01-05  mark iredell
!
! usage:    call rsearch1(km1,z1,km2,z2,l2)
!   input argument list:
!     km1    integer number of points in the sequence
!     z1     real (km1) sequence values to search
!            (z1 must be monotonic in either direction)
!     km2    integer number of points to search for
!     z2     real (km2) set of values to search for
!            (z2 need not be monotonic)
!
!   output argument list:
!     l2     integer (km2) interval locations from 0 to km1
!            (z2 will be between z1(l2) and z1(l2+1))
!
! subprograms called:
!   sbsrch essl binary search
!   dbsrch essl binary search
!
! remarks:
!   returned values of 0 or km1 indicate that the given search value
!   is outside the range of the sequence.
!
!   if a search value is identical to one of the sequence values
!   then the location returned points to the identical value.
!   if the sequence is not strictly monotonic and a search value is
!   identical to more than one of the sequence values, then the
!   location returned may point to any of the identical values.
!
!   if l2(k)=0, then z2(k) is less than the start point z1(1)
!   for ascending sequences (or greater than for descending sequences).
!   if l2(k)=km1, then z2(k) is greater than or equal to the end point
!   z1(km1) for ascending sequences (or less than or equal to for
!   descending sequences).  otherwise z2(k) is between the values
!   z1(l2(k)) and z1(l2(k+1)) and may equal the former.
!
! attributes:
!   language: fortran
!
!$$$
  use machine,only:kint_mpi
  implicit none
  integer,intent(in):: km1,km2
  real,intent(in):: z1(km1),z2(km2)
  integer,intent(out):: l2(km2)
  integer(kint_mpi) incx,n,incy,m,indx(km2),rc(km2),iopt
  integer k2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find the surrounding input interval for each output point.
  if(z1(1).le.z1(km1)) then
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  input coordinate is monotonically ascending.
    incx=1
    n=km2
    incy=1
    m=km1
    iopt=1
    if(digits(1.).lt.digits(1._8)) then
      call sbsrch(z2,incx,n,z1,incy,m,indx,rc,iopt)
    else
      call dbsrch(z2,incx,n,z1,incy,m,indx,rc,iopt)
    endif
    do k2=1,km2
      l2(k2)=indx(k2)-rc(k2)
    enddo
  else
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  input coordinate is monotonically descending.
    incx=1
    n=km2
    incy=-1
    m=km1
    iopt=0
    if(digits(1.).lt.digits(1._8)) then
      call sbsrch(z2,incx,n,z1,incy,m,indx,rc,iopt)
    else
      call dbsrch(z2,incx,n,z1,incy,m,indx,rc,iopt)
    endif
    do k2=1,km2
      l2(k2)=indx(k2)-rc(k2)
    enddo
  else
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  input coordinate is monotonically descending.
    incx=1
    n=km2
    incy=-1
    m=km1
    iopt=0
    if(digits(1.).lt.digits(1._8)) then
      call sbsrch(z2,incx,n,z1,incy,m,indx,rc,iopt)
    else
      call dbsrch(z2,incx,n,z1,incy,m,indx,rc,iopt)
    endif
    do k2=1,km2
      l2(k2)=km1+1-indx(k2)
    enddo
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
#endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine tpause(km,p,u,v,t,h,ptp,utp,vtp,ttp,htp,shrtp)
!$$$  Subprogram documentation block
!
! Subprogram: tpause     Compute tropopause level fields
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram finds the tropopause level and computes fields 
!   at the tropopause level.  The tropopause is defined as the lowest level
!   above 500 mb which has a temperature lapse rate of less than 2 K/km.
!   The lapse rate must average less than 2 K/km over a 2 km depth.
!   If no such level is found below 50 mb, the tropopause is set to 50 mb.
!   The tropopause fields are interpolated linearly in lapse rate.
!   The tropopause pressure is found hydrostatically.
!   The tropopause wind shear is computed as the partial derivative
!   of wind speed with respect to height at the tropopause level.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call tpause(km,p,u,v,t,h,ptp,utp,vtp,ttp,htp,shrtp)
!   Input argument list:
!     km       integer number of levels
!     p        real (km) pressure (Pa)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     t        real (km) temperature (K)
!     h        real (km) height (m)
!   Output argument list:
!     ptp      real tropopause pressure (Pa)
!     utp      real tropopause x-component wind (m/s)
!     vtp      real tropopause y-component wind (m/s)
!     ttp      real tropopause temperature (K)
!     htp      real tropopause height (m)
!     shrtp    real tropopause wind shear (1/s)
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   rsearch       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons
    implicit none
    integer,intent(in):: km
    real,intent(in),dimension(km):: p,u,v,t,h
    real,intent(out):: ptp,utp,vtp,ttp,htp,shrtp
    real,parameter:: ptplim(2)=(/500.e+2,50.e+2/),gamtp=2.e-3,hd=2.e+3
    real gamu,gamd,td,gami,wtp,spdu,spdd
    integer klim(2),k,kd,ktp
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find tropopause level
#ifdef ibm_sp
    call rsearch1(km-2,p(2),2,ptplim(1),klim(1))
#else
    call rsearch(1,km-2,1,1,p(2),2,1,1,ptplim(1),1,1,klim(1))
#endif
!    klim(1)=klim(1)+2
    klim(1)=klim(1)+1
!    klim(2)=klim(2)+1
    klim(2)=klim(2)+2
    gamd=1.e+9
    ktp=klim(2)
    wtp=0
!    do k=klim(1),klim(2)
    do k=klim(1),klim(2),-1
!      gamu=(t(k-1)-t(k+1))/(h(k+1)-h(k-1))
      gamu=(t(k+1)-t(k-1))/(h(k-1)-h(k+1))
      if(gamu.le.gamtp) then
!        call rsearch1(km-k-1,h(k+1),1,h(k)+hd,kd)
#ifdef ibm_sp
	call rsearch1(k-2,h(2),1,h(k)+hd,kd)
#else
	call rsearch(1,k-2,1,1,h(2),1,1,1,h(k)+hd,1,1,kd)
#endif
!        td=t(k+kd)+(h(k)+hd-h(k+kd))/(h(k+kd+1)-h(k+kd))*(t(k+kd+1)-t(k+kd))
	td=t(kd+2)+(h(k)+hd-h(2+kd))/(h(kd+1)-h(2+kd))*(t(kd+1)-t(2+kd))
        gami=(t(k)-td)/hd
        if(gami.le.gamtp) then
          ktp=k
          wtp=(gamtp-gamu)/(max(gamd,gamtp+0.1e-3)-gamu)
          exit
        endif
      endif
      gamd=gamu
    enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  compute tropopause level fields
    utp=u(ktp)-wtp*(u(ktp)-u(ktp-1))
    vtp=v(ktp)-wtp*(v(ktp)-v(ktp-1))
    ttp=t(ktp)-wtp*(t(ktp)-t(ktp-1))
    htp=h(ktp)-wtp*(h(ktp)-h(ktp-1))
    ptp=p(ktp)*exp((h(ktp)-htp)*(1-0.5*(ttp/t(ktp)-1))/(con_rog*t(ktp)))
    spdu=sqrt(u(ktp)**2+v(ktp)**2)
    spdd=sqrt(u(ktp-1)**2+v(ktp-1)**2)
    shrtp=(spdu-spdd)/(h(ktp)-h(ktp-1))
    
    utp=u(ktp)-wtp*(u(ktp)-u(ktp+1))
    vtp=v(ktp)-wtp*(v(ktp)-v(ktp+1))
    ttp=t(ktp)-wtp*(t(ktp)-t(ktp+1))
    htp=h(ktp)-wtp*(h(ktp)-h(ktp+1))
    ptp=p(ktp)*exp((h(ktp)-htp)*(1-0.5*(ttp/t(ktp)-1))/(con_rog*t(ktp)))
    spdu=sqrt(u(ktp)**2+v(ktp)**2)
    spdd=sqrt(u(ktp+1)**2+v(ktp+1)**2)
    shrtp=(spdu-spdd)/(h(ktp)-h(ktp+1))
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine mxwind(km,p,u,v,t,h,pmw,umw,vmw,tmw,hmw)
!$$$  Subprogram documentation block
!
! Subprogram: mxwind     Compute maximum wind level fields
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram finds the maximum wind level and computes fields 
!   at the maximum wind level.  The maximum wind level is searched for
!   between 500 mb and 100 mb.  The height and wind speed at the maximum wind
!   speed level is calculated by assuming the wind speed varies quadratically
!   in height in the neighborhood of the maximum wind level.  The other fields
!   are interpolated linearly in height to the maximum wind level.
!   The maximum wind level pressure is found hydrostatically.
!
! Program history log:
!   1999-10-18  Mark Iredell
!   2005-02-02  Mark Iredell  changed upper limit to 100 mb
!
! Usage:  call mxwind(km,p,u,v,t,h,pmw,umw,vmw,tmw,hmw)
!   Input argument list:
!     km       integer number of levels
!     p        real (km) pressure (Pa)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     t        real (km) temperature (K)
!     h        real (km) height (m)
!   Output argument list:
!     pmw      real maximum wind level pressure (Pa)
!     umw      real maximum wind level x-component wind (m/s)
!     vmw      real maximum wind level y-component wind (m/s)
!     tmw      real maximum wind level temperature (K)
!     hmw      real maximum wind level height (m)
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   rsearch       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons
    implicit none
    integer,intent(in):: km
    real,intent(in),dimension(km):: p,u,v,t,h
    real,intent(out):: pmw,umw,vmw,tmw,hmw
    real,parameter:: pmwlim(2)=(/500.e+2,100.e+2/)
    integer klim(2),k,kmw
    real spd(km),spdmw,wmw,dhd,dhu,shrd,shru,dhmw,ub,vb,spdb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find maximum wind level
#ifdef ibm_sp
    call rsearch1(km,p(1),2,pmwlim(1),klim(1))
#else
    call rsearch(1,km,1,1,p(1),2,1,1,pmwlim(1),1,1,klim(1))
#endif
!    klim(1)=klim(1)+1
    klim(2)=klim(2)+1
!    spd(klim(1):klim(2))=sqrt(u(klim(1):klim(2))**2+v(klim(1):klim(2))**2)
    spd(klim(2):klim(1))=sqrt(u(klim(2):klim(1))**2+v(klim(2):klim(1))**2)
    spdmw=spd(klim(1))
    kmw=klim(1)
!    do k=klim(1)+1,klim(2)
    do k=klim(1)-1,klim(2),-1
      if(spd(k).gt.spdmw) then
        spdmw=spd(k)
        kmw=k
      endif
    enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find speed and height at the maximum wind level
    if(kmw.eq.klim(1).or.kmw.eq.klim(2)) then
      hmw=h(kmw)
      spdmw=spd(kmw)
      wmw=0.
    else
!      dhd=h(kmw)-h(kmw-1)
      dhd=h(kmw)-h(kmw+1) !post counts top down
!      dhu=h(kmw+1)-h(kmw)
      dhu=h(kmw-1)-h(kmw)
!      shrd=(spd(kmw)-spd(kmw-1))/(h(kmw)-h(kmw-1))
      shrd=(spd(kmw)-spd(kmw+1))/(h(kmw)-h(kmw+1))
!      shru=(spd(kmw)-spd(kmw+1))/(h(kmw+1)-h(kmw))
      shru=(spd(kmw)-spd(kmw-1))/(h(kmw-1)-h(kmw))
      dhmw=(shrd*dhu-shru*dhd)/(2*(shrd+shru))
      hmw=h(kmw)+dhmw
      spdmw=spd(kmw)+dhmw**2*(shrd+shru)/(dhd+dhu)
!      if(dhmw.gt.0) kmw=kmw+1
      if(dhmw.gt.0) kmw=kmw-1
!      wmw=(h(kmw)-hmw)/(h(kmw)-h(kmw-1))
      wmw=(h(kmw)-hmw)/(h(kmw)-h(kmw+1))
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  compute maximum wind level fields
!    ub=u(kmw)-wmw*(u(kmw)-u(kmw-1))
    ub=u(kmw)-wmw*(u(kmw)-u(kmw+1))
!    vb=v(kmw)-wmw*(v(kmw)-v(kmw-1))
    vb=v(kmw)-wmw*(v(kmw)-v(kmw+1))
    spdb=max(sqrt(ub**2+vb**2),1.e-6)
    umw=ub*spdmw/spdb
    vmw=vb*spdmw/spdb
!    tmw=t(kmw)-wmw*(t(kmw)-t(kmw-1))
    tmw=t(kmw)-wmw*(t(kmw)-t(kmw+1))
    pmw=p(kmw)*exp((h(kmw)-hmw)*(1-0.5*(tmw/t(kmw)-1))/(con_rog*t(kmw)))
  end subroutine 
