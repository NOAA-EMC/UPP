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
    use physcons_post, only: con_cp, con_g, con_rd, con_rocp
!
    implicit none
    integer,intent(in):: km
    real,intent(in), dimension(km):: p,px,py,t,tx,ty,h,u,v,av
    real,intent(out),dimension(km):: hm,s,bvf2,pvn,theta,sigma,pvu
!   real,parameter:: hhmin=500.,t0=2.e2,p0=1.e5
    real,parameter:: hhmin=5.,t0=2.e2,p0=1.e5
    integer k,kd,ku,k2(2)
    real cprho,sx,sy,sz,uz,vz
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    do k=1,km
      hm(k) = con_cp*t(k) + con_g*h(k)
      s(k)  = con_cp*log(t(k)/t0) - con_rd*log(p(k)/p0)
    enddo
    do k=1,km
      call rsearch1(km,h,2,(/h(k)-hhmin,h(k)+hhmin/),k2)
!     kd = max(k2(1),1)
!     ku = min(k2(2)+1,km)
!     kd = min(k2(1),km)   ! Chuang: post counts from top down, redefine lower bound
      kd = min(k2(1)+1,km) ! Chuang: post counts from top down,
!     ku = max(k2(2)-1,1)
      ku = max(k2(2),1)
      if(ku==1) kd=2       ! Chuang: make sure ku ne kd at model top
      cprho    = p(k)/(con_rocp*t(k))
      sx       = con_cp*tx(k)  / t(k)-con_rd*px(k)/p(k)
      sy       = con_cp*ty(k)  / t(k)-con_rd*py(k)/p(k)
      sz       = (s(ku)-s(kd)) / (h(ku)-h(kd))
      uz       = (u(ku)-u(kd)) / (h(ku)-h(kd))
      vz       = (v(ku)-v(kd)) / (h(ku)-h(kd))
      bvf2(k)  = con_g/con_cp*sz
      pvn(k)   = (av(k)*sz - vz*sx + uz*sy) / cprho
      theta(k) = t0*exp(s(k)/con_cp)
      sigma(k) = t(k)/con_g*bvf2(k)
      pvu(k)   = 1.e6*theta(k)*pvn(k)
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
!   rsearch1       search for a surrounding real interval
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
    call rsearch1(km,theta(1),kth,th(1),loc(1))
    do k=1,kth
      l = loc(k)
      lth(k) = l > 0 .and.l < km
      if(lth(k)) then
        w = log(th(k)/theta(l)) / log(theta(l+1)/theta(l))
        uth(k)     = u(l)     + w*(u(l+1)-u(l))
        vth(k)     = v(l)     + w*(v(l+1)-v(l))
        hth(k)     = h(l)     + w*(h(l+1)-h(l))
        tth(k)     = t(l)     + w*(t(l+1)-t(l))
        zth(k)     = pvu(l)   + w*(pvu(l+1)-pvu(l))
        sigmath(k) = sigma(l) + w*(sigma(l+1)-sigma(l))
        rhth(k)    = rh(l)    + w*(rh(l+1)-rh(l))
!	pth(k)     = p(l)     + w*(p(l+1)-p(l))
        oth(k)     = omga(l)  + w*(omga(l+1)-omga(l))
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
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons_post, only: con_rog
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
      call rsearch1(km,p,1,pvpb(k),l1)
      call rsearch1(km,p,1,pvpt(k),l2)
!     l1=l1+1
      l = 0
      if(pv(k) >= 0.) then
!       do lu=l2-1,l1,-1
!       do lu=l2,l1-1 ! Chuang: post counts top down	
        do lu=l2+2,l1 ! Chuang: post counts top down
!         if(pv(k).lt.pvu(lu+1).and.pv(k).ge.pvu(lu)) then
          if(pv(k) >= pvu(lu+1).and.pv(k) < pvu(lu)) then
            call rsearch1(km,p,1,p(lu)+pd,ld)
!           if(all(pv(k).ge.pvu(ld:lu-1))) then
            if(all(pv(k) >= pvu(lu+1:ld))) then
              l = lu
              exit
            endif
          endif
        enddo
      else
!       do lu=l2-1,l1,-1
!       do lu=l2,l1-1 ! Chuang: post counts top down	
        do lu=l2+2,l1 ! Chuang: post counts top down
!         if(pv(k).gt.pvu(lu+1).and.pv(k).le.pvu(lu)) then
          if(pv(k) <= pvu(lu+1).and.pv(k) > pvu(lu)) then
            call rsearch1(km,p,1,p(lu)+pd,ld)
!           if(all(pv(k).le.pvu(ld:lu-1))) then
            if(all(pv(k) <= pvu(lu+1:ld))) then
              l = lu
              exit
            endif
          endif
        enddo
      endif
      lpv(k) = l > 0
      if(lpv(k)) then
        w = (pv(k)-pvu(l))/(pvu(l+1)-pvu(l))
        upv(k) = u(l) + w*(u(l+1)-u(l))
        vpv(k) = v(l) + w*(v(l+1)-v(l))
        hpv(k) = h(l) + w*(h(l+1)-h(l))
        tpv(k) = t(l) + w*(t(l+1)-t(l))
        ppv(k) = p(l)*exp((h(l)-hpv(k))*(1-0.5*(tpv(k)/t(l)-1))/(con_rog*t(l)))

        spdu   = sqrt(u(l+1)*u(l+1) + v(l+1)*v(l+1))
        spdd   = sqrt(u(l)*u(l)     + v(l)*v(l))
        spv(k) = (spdu-spdd) / (h(l+1)-h(l))
      endif
    enddo
  end subroutine
!-------------------------------------------------------------------------------

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
  implicit none
  integer,intent(in):: km1,km2
  real,intent(in):: z1(km1),z2(km2)
  integer,intent(out):: l2(km2)
  integer k1,k2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find the surrounding input interval for each output point.
  if(z1(1) <= z1(km1)) then
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  input coordinate is monotonically ascending.
    do k2=1,km2
      if (z1(1) >= z2(k2)) then
        l2(k2) = 1
      else
      l2(k2)=km1
      do k1=1,km1-1
        if(z1(k1) <= z2(k2) .and. z1(k1+1) > z2(k2)) then
          l2(k2) = k1
          exit
        endif
      enddo
     endif
    enddo
  else
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  input coordinate is monotonically descending.
    do k2=1,km2
     if (z1(1) <= z2(k2)) then
        l2(k2) = 1
      else
      l2(k2)=km1
      do k1=km1,2,-1
        if(z2(k2) >= z1(k1) .and. z2(k2) < z1(k1-1)) then
          l2(k2) = k1-1
          exit
        endif
      enddo
     endif
    enddo
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
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
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons_post, only: con_rog
    implicit none
    integer,intent(in):: km
    real,intent(in),dimension(km):: p,u,v,t,h
    real,intent(out):: ptp,utp,vtp,ttp,htp,shrtp
    real,parameter:: ptplim(2)=(/500.e+2,50.e+2/),gamtp=2.e-3,hd=2.e+3
    real gamu,gamd,td,gami,wtp,spdu,spdd
    integer klim(2),k,kd,ktp
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find tropopause level
    call rsearch1(km-2,p(2),2,ptplim,klim)
    klim(1)=klim(1)+1
    klim(2)=klim(2)+2
    ! klim(1) > klim(2) or loops does not run ; klim(2) has a
    ! minimum value of 3 to insure k-2 != 0 in called subprogram
    gamd=1.e+9
    ktp=klim(2)
    wtp=0
!    do k=klim(1),klim(2)
    do k=klim(1),klim(2),-1
!      gamu=(t(k-1)-t(k+1))/(h(k+1)-h(k-1))
      gamu=(t(k+1)-t(k-1))/(h(k-1)-h(k+1))
      if(gamu.le.gamtp) then
!        call rsearch1(km-k-1,h(k+1),1,h(k)+hd,kd)
	call rsearch1(k-2,h(2),1,h(k)+hd,kd)
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
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons_post, only: con_rog
    implicit none
    integer,intent(in):: km
    real,intent(in),dimension(km):: p,u,v,t,h
    real,intent(out):: pmw,umw,vmw,tmw,hmw
    real,parameter:: pmwlim(2)=(/500.e+2,100.e+2/)
    integer klim(2),k,kmw
    real spd(km),spdmw,wmw,dhd,dhu,shrd,shru,dhmw,ub,vb,spdb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find maximum wind level
    call rsearch1(km,p(1),2,pmwlim,klim)
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
  
subroutine mptgen(mpirank,mpisize,nd,jt1,jt2,j1,j2,jx,jm,jn)
!$$$  Subprogram documentation block
!
! Subprogram:    mptgen      Generate grid decomposition dimensions
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram decomposes total dimensions of a problem
!   into smaller domains to be managed on a distributed memory system.
!   The last dimension given is decomposed first.  If more decompositions
!   are possible, the next to last dimension is decomposed next, and so on.
!   The transpositions between decompositions should be done by mptran*.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:    call mptgen(mpirank,mpisize,nd,jt1,jt2,j1,j2,jx,jm,jn)
!   Input argument list:
!     mpirank  integer(kint_mpi) rank of the process (from mpi_comm_rank)
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     nd       integer(kint_mpi) number of dimensions to decompose
!     jt1      integer(kint_mpi) (nd) lower bounds of total dimensions
!     jt2      integer(kint_mpi) (nd) upper bounds of total dimensions
!   Output argument list:
!     j1       integer(kint_mpi) (nd) lower bounds of local decompositions
!     j2       integer(kint_mpi) (nd) upper bounds of local decompositions
!     jx       integer(kint_mpi) (nd) local size of decompositions
!     jm       integer(kint_mpi) (nd) maximum size of decompositions
!     jn       integer(kint_mpi) (nd) number of decompositions
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use machine_post,only:kint_mpi
  implicit none
  integer(kint_mpi),intent(in):: mpirank,mpisize,nd,jt1(nd),jt2(nd)
  integer(kint_mpi),intent(out):: j1(nd),j2(nd),jx(nd),jm(nd),jn(nd)
  integer msize,mrank,msn,mrn,n
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  msize=mpisize
  mrank=mpirank
  do n=nd,1,-1
    if(jt2(n).ge.jt1(n)) then
      jm(n)=(jt2(n)-jt1(n))/msize+1
      msn=max(msize/(jt2(n)-jt1(n)+1),1)
      if(n.eq.1) msn=1
      jn(n)=msize/msn
      mrn=mrank/msn
      j1(n)=min(jt1(n)+jm(n)*mrn,jt2(n)+1)
      j2(n)=min(jt1(n)+jm(n)*mrn+jm(n)-1,jt2(n))
      jx(n)=j2(n)-j1(n)+1
      msize=msn
      mrank=mod(mrank,msn)
      write(0,*)' mrank=',mrank,' j1=',j1(n),' j2=',j2(n),' jx=',jx(n),' jm=',jm(n)
    else
      jm(n)=0
      jn(n)=1
      j1(n)=jt1(n)
      j2(n)=jt2(n)
      jx(n)=0
    endif
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mptranr4(mpicomm,mpisize,im,ida,idb,&
                    jm,jma,jmb,jda,km,kma,kmb,kdb,a,b,ta,tb)
!$$$  Subprogram documentation block
!
! Subprogram:    mptranr4    Transpose grid decompositions
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram transposes an array of data from one
!   grid decomposition to another by using message passing.
!   The grid decompositions should be generated by mptgen.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:  call mptranr4(mpicomm,mpisize,im,ida,idb,&
!                       jm,jma,jmb,jda,km,kma,kmb,kdb,a,b)
!   Input argument list:
!     mpicomm  integer(kint_mpi) mpi communicator
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     im       integer(kint_mpi) undecomposed range
!     ida      integer(kint_mpi) undecomposed input dimension
!     idb      integer(kint_mpi) undecomposed output dimension
!     jm       integer(kint_mpi) output grid decomposition size
!     jma      integer(kint_mpi) input grid undecomposed range
!     jmb      integer(kint_mpi) output grid decomposed range
!     jda      integer(kint_mpi) input grid undecomposed dimension
!     km       integer(kint_mpi) input grid decomposition size
!     kma      integer(kint_mpi) input grid decomposed range
!     kmb      integer(kint_mpi) output grid undecomposed range
!     kdb      integer(kint_mpi) output grid undecomposed dimension
!     a        real(4) (ida,jda,kma) input array
!   Output argument list:
!     b        real(4) (idb,kdb,jmb) output array
!     ta,tb    real(4) (im,jm,km,mpisize)  work arrays
!
! Subprograms called:
!   mpi_alltoall  mpi exchange of data between every process pair
!
! Remarks:
!   While this routine serves a wide variety of scalable transpose functions
!   for multidimensional grids,
!     (a) it does not work with nonrectanguloid grids;
!     (b) it does not do any load balancing;
!     (c) it does not do any communication hiding.
!
!   This subprogram must be used rather than mpi_alltoall
!   in any of the following cases:
!     (a) The undecomposed range is less than the respective dimension
!         (either im.lt.ida or im.lt.idb)
!     (b) The decomposition size is greater than one
!         (either km.gt.1 or jm.gt.1)
!     (c) The decomposed range is ever zero
!         (either kma.eq.0 or jmb.eq.0 for any process)
!     (d) The output grid range is not the full extent
!         (either kmb.lt.mpisize or kmb.lt.kda or jma.lt.mpisize or jma.lt.jda)
!   If none of these conditions apply, mpi_alltoall could be used directly
!   rather than this subprogram and would be more efficient.
!
!   Example 1.  Transpose a 1000 x 10000 matrix.
!
!!!   include 'mpif.h'                                     ! use mpi
!!!   parameter(jt=1000,kt=10000)                          ! set problem size
!!!   real,allocatable:: a(:,:),b(:,:)                     ! declare arrays
!!!   call mpi_init(ierr)                                  ! initialize mpi
!!!   call mpi_comm_rank(MPI_COMM_WORLD,mpirank,ierr)      ! get mpi rank
!!!   call mpi_comm_size(MPI_COMM_WORLD,mpisize,ierr)      ! get mpi size
!!!   call mptgen(mpirank,mpisize,1,1,jt,j1,j2,jx,jm,jn)   ! decompose output
!!!   call mptgen(mpirank,mpisize,1,1,kt,k1,k2,kx,km,kn)   ! decompose input
!!!   allocate(a(jt,k1:k2),b(kt,j1:j2))                    ! allocate arrays
!!!   a=reshape((/((j+k,j=1,jt),k=k1,k2)/),(/jt,k2-k1+1/)) ! initialize input
!!!   call mptranr4(MPI_COMM_WORLD,mpisize,1,1,1,          ! transpose arrays
!!!  &              jm,jt,j2-j1+1,jt,km,k2-k1+1,kt,kt,a,b)
!!!   print '(2i8,f16.1)',((k,j,b(k,j),k=2000,kt,2000),    ! print some values
!!!  &                    j=((j1-1)/200+1)*200,j2,200)
!!!   call mpi_finalize(ierr)                              ! finalize mpi
!!!   end
!
!   This transpose took 0.6 seconds on 4 2-way winterhawk nodes.
!   A 20000x10000 transpose took 3.4 seconds on 16 2-way winterhawk nodes.
!   Thus a transpose may take about 1 second for every 16 Mb per node.
!
! Attributes:
!   language: fortran 90
!
!$$$
  use machine_post,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: mpicomm,mpisize
  integer(kint_mpi),intent(in):: im,ida,idb
  integer(kint_mpi),intent(in):: jm,jma,jmb,jda
  integer(kint_mpi),intent(in):: km,kma,kmb,kdb
  real(4),dimension(ida,jda,kma),intent(in):: a
  real(4),dimension(idb,kdb,jmb),intent(out):: b
  real(4),dimension(im,jm,km,mpisize),intent(inout):: ta,tb
  integer(4) jmb1(1),jmbf(mpisize),kma1(1),kmaf(mpisize)
  integer(kint_mpi)::i,j,k,l,ierr,ja,kb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  jmb1(1) = jmb
  call mpi_allgather(jmb1,1,MPI_INTEGER4,jmbf,1,MPI_INTEGER4,mpicomm,ierr)
  kma1(1) = kma
  call mpi_allgather(kma1,1,MPI_INTEGER4,kmaf,1,MPI_INTEGER4,mpicomm,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose input array
!$omp parallel do private(i,j,k,l,ja)
  do l=1,mpisize
    do k=1,kma
      do j=1,jm
        ja = j + sum(jmbf(1:l-1))
        if(ja <= jma) then
          do i=1,im
            ta(i,j,k,l) = a(i,ja,k)
          enddo
        endif
      enddo
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  externally transpose data
  call mpi_alltoall(ta,im*jm*km,MPI_REAL4,tb,im*jm*km,MPI_REAL4,mpicomm,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose output array
!$omp parallel do private(i,j,k,l,kb)
  do l=1,mpisize
    do k=1,km
      kb = k + sum(kmaf(1:l-1))
      if(kb <= kmb) then
        do j=1,jmb
          do i=1,im
            b(i,kb,j) = tb(i,j,k,l)
          enddo
        enddo
      endif
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-----------------------------------------------------------------------
