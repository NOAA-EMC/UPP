!---------------------------------------------------------------------
! The GFIP algorithm computes the probability of icing within a model
!    column given the follow input data
!
! 2-D data: convective precip (xacp),
!           non-convective precip (xncp)
!           the topography height (xalt)
!           the latitude and longitude (xlat and xlon)
!           d  (xlat and xlon)
!           the number of vertical levels (num_z) = 47 in my GFS file
! 3-D data
!            pressure(num_z) in PA
!            temperature(num_z) in K
!            rh(num_z) in %
!            hgt(num_z) in GPM)
!            cwat(num_z) in kg/x kg
!            vv(num_z) in Pa/sec
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------+
! This subroutine calculates the icing potential for a GFS column
! First the topography is mapped to the model's vertical coordinate
!      in (topoK) 
! Then cloud layers are defined in (cloud_layer) 
!      up to 12 layers are possible 
! The icing is computed in (icing_pot). The icing  
! output should range from 0 - 1. 
!
!-----------------------------------------------------------------------+
!234567890


      subroutine icing_algo(i,j,pres,temp,rh,hgt,cwat,vv,num_z,xlat,xlon, &
                      xalt,xncp,xacp,ice_pot)
      implicit none

      integer i,j, l      
      integer num_z,surface,region,num_lyr
      integer,dimension(12) :: cld_top, cld_base

      real xlat, xlon, xalt, xncp, xacp, topok
      real,dimension(12) :: ctt, cbt, thick
      real,dimension(num_z) :: pres, hgt, rh, temp, cwat, vv, ice_pot


      if(i==50 .and. j==50)then
       print*,'sample input to FIP ',i,j,num_z,xlat,xlon,xalt,xncp, &
       xacp
       do l=1,num_z
        print*,'l,P,T,RH,H,CWM,VV',l,pres(l),temp(l),rh(l),hgt(l),cwat(l),vv(l)
       end do
      end if
       	
       surface = topoK(hgt,xalt,num_z)
       call cloud_layer(xlat,xlon,hgt,rh,temp,surface,num_z,cwat, &
                     region,num_lyr,cld_top,ctt,cld_base,cbt,thick)
       call icing_pot(hgt,rh,temp,surface,num_z,cwat,vv, &
          region,num_lyr,cld_top,ctt,cld_base,ice_pot)


      return


      end


!-----------------------------------------------------------------------+
       real function dew_pt(t,rh)
      real tc, eps, vapr, top, bottom, td  

      tc = t - 273.15
      if (rh.lt.0.0001) then 
       rh = 0.0001
      endif 

      eps = 6.112
      vapr = rh * (vapor_press(tc)/100.0)
      top = 243.5 * (log(eps) - log(vapr))
      bottom = log(vapr) - log(eps) - 17.67
      dew_pt = top/bottom

      return 
      end

!-----------------------------------------------------------------------+
!  press in mb, T and Td in degrees C
      real function thetae(press, t, td)
      real rmix, e, thtam
      real mix_ratio 
        press = press/100.0
        t = t - 273.15
        td = td - 273.15
        rmix = mix_ratio(td,press)
        e = (2.0/7.0) * ( 1 - (0.001 * 0.28 * rmix))
        thtam = (t + 273.15) * ((1000.0/press)**e) 
        thetae = thtam * exp( (3.376/tLCL(t,td) - 0.00254) *  &
               (rmix * (1 + (81 * 0.001 * rmix)) ))

      return
      end


!-----------------------------------------------------------------------+

      real function vapor_press(t);
         vapor_press = 6.112 * exp( (17.67*t)/(t+243.5));

      return
      end


!-----------------------------------------------------------------------+

      real function  tLCL(t, td);
      real tk, dk 

         tk = t + 273.15;
         dk = td + 273.15;
         tLCL = ( 1 / (1/(dk - 56) + log(tk/dk)/800.0)) + 56.0;

      return
      end

!-----------------------------------------------------------------------+

       real function  mix_ratio(td, press);
         real corr, e

         corr = 1.001 + ( (press - 100) /900) * 0.0034;
         e = vapor_press(td) * corr;
         mix_ratio = 0.62197 * (e/(press-e)) * 1000;


      return
      end



!-----------------------------------------------------------------------+

      real function topoK(hgt,alt,num_z) 

      real hgt(num_z)
      integer k

      topoK = num_z 
      do 110 k=2,num_z
       if ((hgt(k-1).gt.alt).and.(hgt(k).le.alt)) then
        topoK = k-1
       endif
  110 continue

      return
      end  

!-----------------------------------------------------------------------+

      subroutine cloud_layer(xlat,xlon,hgt,rh,temp,surface,num_z,cwat, &
                       region,num_lyr,cld_top,ctt,cld_base,cbt,thick)

      real hgt(num_z),rh(num_z),temp(num_z),cwat(num_z)
      real ctt(12),cbt(12),t_rh,thick(12) 
      integer cloud(num_z),in_cld,cur_base,surface 
      integer region,cld_top(12),cld_base(12),num_lyr


! get the global region and set the RH thresholds           

      if (abs(xlat).lt.23.5) then
        t_rh = 80.0
        region = 1
      elseif ( (abs(xlat).ge.23.5).and.(abs(xlat).lt.66)) then
        t_rh = 75.0
        region =2 
      else
        t_rh = 70.0
        region =2 
      endif 

! zero the clouid_on flag 
       do 115 k=1,num_z
        cloud(k) = 0
 115  continue

! loop from the top (start at 2 so we can have n+1) )
! bottom is num+z and top is 1

      num_lyr = 0
      in_cld = 0
      cur_base = 1 

      do 120 k=2,surface 
       if (cur_base.le.k) then
        if ((rh(k-1).lt.t_rh).and.(in_cld.eq.0)) then
         if ((rh(k).ge.t_rh).and.(rh(k+1).ge.t_rh)) then
         num_lyr = num_lyr + 1 
         in_cld = 1
         cld_top(num_lyr) = k 
         ctt(num_lyr) = temp(k)

! find the cloud base
       do 125 kk=cld_top(num_lyr),surface-1
        if  ((rh(kk-1).ge.t_rh).and.(rh(kk).lt.t_rh)) then  
         if ((rh(kk+1).lt.t_rh).and.(in_cld.eq.1)) then 
         cld_base(num_lyr) = kk
         cur_base = kk
         cbt(num_lyr) = temp(kk)
         in_cld = 0
        endif
        endif
 125    continue
        endif
        endif
        endif

 120   continue



! if the loop exits still in cloud make the cloud base the surface
       if (in_cld.eq.1) then 
        cld_base(num_lyr) = surface
        cbt(num_lyr) = temp(surface) 
       endif

!  get the cloud thickness
       do 130 n=1,num_lyr
         thick(n) = hgt(cld_top(n)) - hgt(cld_base(n))
 130   continue



      return
      end

!-----------------------------------------------------------------------+

      subroutine icing_pot(hgt,rh,temp,surface,num_z,cwat,vv, &
          region,num_lyr,cld_top,ctt,cld_base,ice_pot)

      real hgt(num_z),rh(num_z),temp(num_z),cwat(num_z),vv(num_z) 
      real ctt(num_lyr),cbt(num_lyr),t_rh,thick(num_lyr)
      real ice_pot(num_z),ice_ctt(num_z), slw(num_z), slw_fac(num_z) 
      real vv_fac(num_z)
      integer surface
      integer region,cld_top(num_lyr),cld_base(num_lyr),num_lyr

      do 200 k=1,num_z
       ice_pot(k) = 0.0
       ice_ctt(k) = 0.0
       vv_fac(k) = 0.0
       slw_fac(k) = 0.0
 200  continue

! apply the cloud top temperature to the cloud layer

      do 205 n=1,num_lyr
       do 210 k=cld_top(n),cld_base(n)  
        ice_ctt(k) = ctt(n)
 210   continue
 205  continue
      
! convert the cwater to slw if the CTT > -14C

      do 212 k=1,num_z
        if((ice_ctt(k).ge.259.15).and.(ice_ctt(k).le.273.15))  then
         slw(k) = cwat(k) * 1000.0
        else 
         slw(k) = 0.0
        endif
 212   continue 



! run the icing algorithm


      do 215 n=1,num_lyr
       do 220 k=cld_top(n),cld_base(n)  
        ice_pot(k) = tmap(temp(k))*rh_map(rh(k))*ctt_map(ice_ctt(k))

! add the VV map
        if (vv_map(vv(k)).ge.0.0) then
          vv_fac(k) = (1.0-ice_pot(k))*(0.25*vv_map(vv(k))) 
        else
          vv_fac(k) = ice_pot(k) * (0.25 * vv_map(vv(k)))  
        endif

! add the slw
       slw_fac(k) = (1.0 - ice_pot(k))*(0.4*slw_map(slw(k)))  

! calculate the final icing potential
      if (ice_pot(k).gt.0.001) then 
       ice_pot(k) = ice_pot(k) + vv_fac(k) + slw_fac(k)
      endif 

! make sure the values don't exceed 1.0
      if (ice_pot(k).gt.1.0) then
       ice_pot(k) = 1.0
      endif 

!        print *,ice_pot(k),temp(k),tmap(temp(k)),rh(k),rh_map(rh(k))
!     *          , ice_ctt(k),ctt_map(ice_ctt(k)),vv(k),vv_map(vv(k))
!     *          , cwat(k), slw(k), slw_map(slw(k))  

 220   continue
 215  continue

      return
      end

!-----------------------------------------------------------------------+
       real function tmap(temp)
       real temp

        if((temp.ge.248.15).and.(temp.le.263.15)) then 
          tmap=((temp-248.15)/15.0)**1.5
        elseif((temp.gt.263.15).and.(temp.lt.268.15)) then 
          tmap=1.0
        elseif((temp.gt.268.15).and.(temp.lt.273.15)) then 
         tmap=((273.15 - temp)/5.0)**0.75
        else
         tmap=0.0;
       endif 

       return
       end
 

!-----------------------------------------------------------------------+
       real function rh_map(rh)
       real rhmap

      
       if (rh.gt.95.0) then
          rh_map=1.0
       elseif ( (rh.le.95.0).and.(rh.ge.50.0) ) then 
          rh_map=((20/9) * ((rh/100.0)-0.5))**3.0
       else
          rh_map=0.0
       endif

       return
       end


!-----------------------------------------------------------------------+

       real function ctt_map(cldtops)

        if((cldtops.ge.261.0).and.(cldtops.lt.280.)) then 
         ctt_map = 1.0
       elseif((cldtops.gt.223.0).and.(cldtops.lt.261.0 )) then 
         ctt_map = 0.2 + 0.8*(((cldtops - 223.0)/38.0)**2.0)
       elseif(cldtops.lt.223.0) then 
         ctt_map = 0.2
       else
         ctt_map = 0.0
       endif 


       return
       end


!-----------------------------------------------------------------------+

      real function vv_map(vv)

        if (vv.gt.0.0) then 
          vv_map = 0.0
        elseif (vv.lt.-0.5) then 
          vv_map = 1.0
        else
          vv_map = -1.0 * (vv/0.5);
        endif 

       return
       end


!-----------------------------------------------------------------------+

       real function slw_map(slw)

       if(slw.gt.0.2) then 
        slw_map = 1.0
       elseif (slw.le.0.001) then 
        slw_map = 0.0
       else
        slw_map = slw/0.2;
       endif 

      return
      end

