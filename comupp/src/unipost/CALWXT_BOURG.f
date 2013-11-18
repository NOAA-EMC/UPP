!$$$  Subprogram documentation block
!
! Subprogram: calwxt_bourg    Calculate precipitation type (Bourgouin)
!   Prgmmr: Baldwin      Org: np22        Date: 1999-07-06
!
! Abstract: This routine computes precipitation type
!    using a decision tree approach that uses the so-called
!    "energy method" of Bourgouin of AES (Canada) 1992
!
! Program history log:
!   1999-07-06  M Baldwin
!   1999-09-20  M Baldwin  make more consistent with bourgouin (1992)
!   2005-08-24  G Manikin  added to wrf post
!   2007-06-19  M Iredell  mersenne twister, best practices
!
! Usage:    call calwxt_bourg(im,jm,jsta_2l,jend_2u,jsta,jend,lm,lp1,   &
!    &                        iseed,g,pthresh,                          &
!    &                        t,q,pmid,pint,lmh,prec,zint,ptype)
!   Input argument list:
!     im       integer i dimension
!     jm       integer j dimension
!     jsta_2l  integer j dimension start point (including haloes)
!     jend_2u  integer j dimension end point (including haloes)
!     jsta     integer j dimension start point (excluding haloes)
!     jend     integer j dimension end point (excluding haloes)
!     lm       integer k dimension
!     lp1      integer k dimension plus 1
!     iseed    integer random number seed
!     g        real gravity (m/s**2)
!     pthresh  real precipitation threshold (m)
!     t        real(im,jsta_2l:jend_2u,lm) mid layer temp (K)
!     q        real(im,jsta_2l:jend_2u,lm) specific humidity (kg/kg)
!     pmid     real(im,jsta_2l:jend_2u,lm) mid layer pressure (Pa)
!     pint     real(im,jsta_2l:jend_2u,lp1) interface pressure (Pa)
!     lmh      real(im,jsta_2l:jend_2u) max number of layers
!     prec     real(im,jsta_2l:jend_2u) precipitation (m)
!     zint     real(im,jsta_2l:jend_2u,lp1) interface height (m)
!   Output argument list:
!     ptype    real(im,jm) instantaneous weather type ()
!              acts like a 4 bit binary
!                1111 = rain/freezing rain/ice pellets/snow
!                where the one's digit is for snow
!                      the two's digit is for ice pellets
!                      the four's digit is for freezing rain
!                  and the eight's digit is for rain
!              in other words...
!                ptype=1 snow
!                ptype=2 ice pellets/mix with ice pellets
!                ptype=4 freezing rain/mix with freezing rain
!                ptype=8 rain
!
! Modules used:
!   mersenne_twister pseudo-random number generator
!
! Subprograms called:
!   random_number    pseudo-random number generator
!
! Attributes:
!   Language: Fortran 90
!
! Remarks: vertical order of arrays must be layer   1 = top
!                                       and layer lmh = bottom
!
!$$$
      subroutine calwxt_bourg_post(im,jm,jsta_2l,jend_2u,jsta,jend,lm,lp1,   &
     &                        iseed,g,pthresh,                          &
     &                        t,q,pmid,pint,lmh,prec,zint,ptype)
      use mersenne_twister, only: random_number
      implicit none
!
!    input:
      integer,intent(in):: im,jm,jsta_2l,jend_2u,jsta,jend,lm,lp1
      integer,intent(in):: iseed
      real,intent(in):: g,pthresh
      real,intent(in):: t(im,jsta_2l:jend_2u,lm)
      real,intent(in):: q(im,jsta_2l:jend_2u,lm)
      real,intent(in):: pmid(im,jsta_2l:jend_2u,lm)
      real,intent(in):: pint(im,jsta_2l:jend_2u,lp1)
      real,intent(in):: lmh(im,jsta_2l:jend_2u)
      real,intent(in):: prec(im,jsta_2l:jend_2u)
      real,intent(in):: zint(im,jsta_2l:jend_2u,lp1)
!
!    output:
      real,intent(out):: ptype(im,jm)
!
      integer i,j,ifrzl,iwrml,l,lhiwrm,lmhk
      real pintk1,areane,tlmhk,areape,pintk2,surfw,area1,dzkl,psfck
      real rn(im*jm*2),r1,r2
!
!     initialize weather type array to zero (ie, off).
!     we do this since we want ptype to represent the
!     instantaneous weather type on return.
!     
!$omp  parallel do
      do j=jsta,jend
      do i=1,im
        ptype(i,j) = 0
      enddo
      enddo
!
      call random_number(rn,iseed)
!
!$omp  parallel do
!$omp& private(a,lmhk,tlmhk,iwrml,psfck,lhiwrm,pintk1,pintk2,area1,
!$omp&         areape,dzkl,surfw,r1,r2)
      do j=jsta,jend
      do i=1,im
      lmhk=nint(lmh(i,j))
      psfck=pint(i,j,lmhk+1)
!
!   skip this point if no precip this time step 
!
      if (prec(i,j).le.pthresh) cycle
!     find the depth of the warm layer based at the surface
!     this will be the cut off point between computing
!     the surface based warm air and the warm air aloft
!
!
!     lowest layer t
!
      tlmhk = t(i,j,lmhk)
      iwrml = lmhk + 1
      if (tlmhk.ge.273.15) then
        do l = lmhk, 2, -1
         if (t(i,j,l).ge.273.15.and.t(i,j,l-1).lt.273.15.and.           &
     &            iwrml.eq.lmhk+1) iwrml = l
          end do
      end if
!
!     now find the highest above freezing level
!
! gsm  added 250 mb check to prevent stratospheric warming situations
!       from counting as warm layers aloft
      lhiwrm = lmhk + 1
      do l = lmhk, 1, -1
          if (t(i,j,l).ge.273.15 .and. pmid(i,j,l).gt.25000.) lhiwrm = l
      end do

!     energy variables
!     surfw is the positive energy between the ground
!     and the first sub-freezing layer above ground
!     areane is the negative energy between the ground
!     and the highest layer above ground
!     that is above freezing
!     areape is the positive energy "aloft"
!     which is the warm energy not based at the ground
!     (the total warm energy = surfw + areape)
!
!     pintk1 is the pressure at the bottom of the layer
!     pintk2 is the pressure at the top of the layer
!     dzkl is the thickness of the layer
!     ifrzl is a flag that tells us if we have hit
!     a below freezing layer
!
      pintk1 = psfck
      ifrzl = 0
      areane = 0.0
      areape = 0.0
      surfw = 0.0                                         

      do l = lmhk, 1, -1
          if (ifrzl.eq.0.and.t(i,j,l).le.273.15) ifrzl = 1
          pintk2=pint(i,j,l)
          dzkl=zint(i,j,l)-zint(i,j,l+1)
          area1 = alog(t(i,j,l)/273.15) * g * dzkl
          if (t(i,j,l).ge.273.15.and. pmid(i,j,l).gt.25000.) then
              if (l.lt.iwrml) areape = areape + area1
              if (l.ge.iwrml) surfw = surfw + area1
          else
              if (l.gt.lhiwrm) areane = areane + abs(area1)
          end if
          pintk1 = pintk2
      end do
      
!
!     decision tree time
!
      if (areape.lt.2.0) then
!         very little or no positive energy aloft, check for
!         positive energy just above the surface to determine rain vs. snow
          if (surfw.lt.5.6) then
!             not enough positive energy just above the surface
!             snow = 1
              ptype(i,j) = 1
          else if (surfw.gt.13.2) then
!             enough positive energy just above the surface
!             rain = 8
              ptype(i,j) = 8
          else
!             transition zone, assume equally likely rain/snow
!             picking a random number, if <=0.5 snow
              r1 = rn(i+im*(j-1))
              if (r1.le.0.5) then
!                 snow = 1
                  ptype(i,j) = 1
              else
!                 rain = 8
                  ptype(i,j) = 8
              end if
          end if
!
      else
!         some positive energy aloft, check for enough negative energy
!         to freeze and make ice pellets to determine ip vs. zr
          if (areane.gt.66.0+0.66*areape) then
!             enough negative area to make ip,
!             now need to check if there is enough positive energy
!             just above the surface to melt ip to make rain
              if (surfw.lt.5.6) then
!                 not enough energy at the surface to melt ip
!                 ice pellets = 2
                  ptype(i,j) = 2
              else if (surfw.gt.13.2) then
!                 enough energy at the surface to melt ip
!                 rain = 8
                  ptype(i,j) = 8
              else
!                 transition zone, assume equally likely ip/rain
!                 picking a random number, if <=0.5 ip
                  r1 = rn(i+im*(j-1))
                  if (r1.le.0.5) then
!                     ice pellets = 2
                      ptype(i,j) = 2
                  else
!                     rain = 8
                      ptype(i,j) = 8
                  end if
              end if
          else if (areane.lt.46.0+0.66*areape) then
!             not enough negative energy to refreeze, check surface temp
!             to determine rain vs. zr
              if (tlmhk.lt.273.15) then
!                 freezing rain = 4
                  ptype(i,j) = 4
              else
!                 rain = 8
                  ptype(i,j) = 8
              end if
          else
!             transition zone, assume equally likely ip/zr
!             picking a random number, if <=0.5 ip
              r1 = rn(i+im*(j-1))
              if (r1.le.0.5) then
!                 still need to check positive energy
!                 just above the surface to melt ip vs. rain
                  if (surfw.lt.5.6) then
!                     ice pellets = 2
                      ptype(i,j) = 2
                  else if (surfw.gt.13.2) then
!                     rain = 8
                      ptype(i,j) = 8
                  else
!                     transition zone, assume equally likely ip/rain
!                     picking a random number, if <=0.5 ip
                      r2 = rn(i+im*(j-1)+im*jm)
                      if (r2.le.0.5) then
!                         ice pellets = 2
                          ptype(i,j) = 2
                      else
!                         rain = 8
                          ptype(i,j) = 8
                      end if
                  end if
              else
!                 not enough negative energy to refreeze, check surface temp
!                 to determine rain vs. zr
                  if (tlmhk.lt.273.15) then
!                     freezing rain = 4
                      ptype(i,j) = 4
                  else
!                     rain = 8
                      ptype(i,j) = 8
                  end if
              end if
          end if
      end if
      end do
      end do
      end
