module gtg_itfa
  use gtg_config
  implicit none
contains

!-----------------------------------------------------------------------
  subroutine ITFA_MWT(nids,indxpicked,kregions,cat,qitfam)
!qix,qit,qitfam,iditfam,nx,ny,nz,imin,imax,
!     1  jmin,jmax,kimin,kimax,mask,printflag,ic,jc,iprt,ioutputflag,
!     2  Qdir)

!     --- Computes an ITFA combination for MWT and outputs as qitfa
!     --- If ioutputflag > 0 qitfa is also stored on disk in directory Qdir as idQ.Q. 
!     --- qix and qit are work arrays.
!     --- The sum is only computed between i=imin,imax, j=jmin,jmax,
!     --- k=kmin,kmax.
!     --- The individual indices are read in from the interpolated 4XX.Q
!     --- files on disk in the Qdir.

    implicit none

    integer,intent(in) :: nids
    integer,intent(in) :: indxpicked(nids)
    integer,intent(in) :: kregions(IM,jsta:jend,MAXREGIONS,2)
    real,intent(in) :: cat(IM,jsta:jend,LM,nids)
    real,intent(inout) :: qitfam(IM,jsta:jend,LM)

    integer :: kpickitfa(nids)
    integer :: nidsmwt


    real :: qs(IM,jsta:jend,LM) ! work array

    integer :: idx
    integer :: n,i,j,k,idx,iregion,kk
    real :: wtindx,wqs,qijk

      real    qix(nx,ny,nz),qit(nx,ny,nz),qitfam(nx,ny,nz)
      integer iditfam
      integer imin,imax,jmin,jmax,kimin,kimax
      integer ioutputflag
      integer ic,jc,printflag,iprt
      integer mask(nx,ny)
      character*200 Qdir
!-----------------------------------------------------------------------
      integer i,j,ki,idx,idQ,ierr
      integer iregion
      real    qmax,qmin
      character*24 cnamei
      integer nimax
      parameter (nimax=101)
      integer kpickitfa(nimax)
!-----------------------------------------------------------------------

    write(iprt,*) 'enter ITFA_MWT'

    qitfam = 0.

    nidsmwt = 0
    kpickitfa = 0
    do i = 1, nids
       if(indxpicked(i) >= 476 .and. indxpicked(i) <= 490) then ! MWT
          nidsmwt = nidsmwt + 1
          kpickitfa(nidsmwt) = indxpicked(i)
       end if
    end do
    if (nidsmwt <= 0) then
       write(*,*) "There is no MWT indices picked"
       return
    end if

!   --- Perform the combination using altitude-dependent static weights
    call itfacompQ(nidsmwt,kpickitfa(1:nidsmwt),kregions,static_wts,&
         cat,clampitfaL,clampitfaH,qitfam)

!   --- Merge the regions
    call MergeItfaRegions(qitfam,nx,ny,nz,imin,imax,jmin,jmax,
     1  kregion,minregion,maxregion,mask,printflag,ic,jc,iprt)
      if(printflag>=2) then
        write(iprt,*) 'after merge'
        do ki=1,nz
          write(iprt,*) 'i,j,k,itfamwt=',ic,jc,ki,qitfam(ic,jc,ki)
        enddo
!         ki=31
!         do i=1,nx
!           write(iprt,*) 'i,j,k,itfamwt=',i,jc,ki,qitfam(i,jc,ki)
!         enddo
      endif
!     --- Write out ITFAMWT to .Q file
      if(ioutputflag>0) then
        idQ=iditfam+399
        cnamei=cname(iditfam)
        call putqix(qitfam,nx,ny,nz,idQ,cnamei,printflag,iprt,Qdir,ierr)
      endif
      call TIstats(qitfam,nx,ny,nz,1,nx,1,ny,kimin,kimax,qmax,qmin,idQ,
     1  iprt)
!
      return
      end
!
      subroutine ITFA_static(qix,qit,qitfad,iditfad,nx,ny,nz,imin,imax,
     1  jmin,jmax,kimin,kimax,mask,printflag,ic,jc,iprt,ioutputflag,
     2  Qdir)
!     --- Computes an ITFA combination using static weights and outputs as qitfa.
!     --- If ioutputflag > 0 qitfa is also stored on disk in directory Qdir as idQ.Q. 
!     --- qix and qit are work arrays.
!     --- The sum is only computed between i=imin,imax, j=jmin,jmax,
!     --- k=kmin,kmax.
!     --- The individual indices are read in from the interpolated 4XX.Q
!     --- files on disk in the Qdir.
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nz
      real    qix(nx,ny,nz),qit(nx,ny,nz),qitfad(nx,ny,nz)
      integer iditfad
      integer nitfa
      integer imin,imax,jmin,jmax,kimin,kimax
      integer ioutputflag
      integer ic,jc,printflag,iprt
      integer mask(nx,ny)
      character*200 Qdir
!-----------------------------------------------------------------------
      integer ki,idx,ierr,i
      integer idQ,ni
      integer iregion
      real    qmin,qmax
      character*24 cnamei
      integer nimax
      parameter (nimax=101)
      integer kpickitfa(nimax)
      include 'scoreparams.inc'
      include 'scorei.inc'
      include 'pireps.inc'   ! npireps
      include 'consts.inc'   ! SPVAL
!-----------------------------------------------------------------------
!
!     --- Initializations
      if(printflag>=1) then
        write(iprt,*)'enter ITFA_static: idQ,remap_option=',iditfad+399,
     1    remap_option
        write(iprt,*) 'minregion,maxregion=',minregion,maxregion
        do iregion=minregion,maxregion
          do idx=1,idmax
            if(ipickitfa(iregion,idx)>0) then
              write(iprt,*)'iregion,idQ,ipickitfa=',iregion,idx+399,
     1          ipickitfa(iregion,idx)
            endif
          enddo
        enddo
      endif
      call TIinit(qitfad,nx,ny,nz)
      call TIinit(qix   ,nx,ny,nz)
      call TIinit(qit   ,nx,ny,nz)
!
!     --- Perform the combination using altitude-dependent static weights
      do iregion=minregion,maxregion
        write(iprt,*) 'computing ITFADEF for iregion=',iregion
!       --- Compute the ITFADEF combination for this region
        do idx=1,idmax
          kpickitfa(idx)=0
          if(ifcsttype(idx)==CAT)kpickitfa(idx)=ipickitfa(iregion,idx)
        enddo
        if(printflag>=1) then
          ni=0
          do idx=1,idmax
            if(kpickitfa(idx)>0) then
              ni=ni+1
              write(iprt,*) 'iregion,idQ,kpickitfa,wt=',iregion,idx+399,
     1         kpickitfa(idx),static_wts(iregion,idx)
            endif
          enddo
          write(iprt,*) 'nindices to use=',ni
        endif
        nitfa=0
!       --- qix and qit are work arrays.  Output is in qitfad
        call itfacompQ(qix,qit,qitfad,static_wts,kpickitfa,
     1    remap_option,iregion,kregion,minregion,maxregion,
     2    clampitfaL,clampitfaH,nx,ny,nz,idmax,nitfa,imin,imax,
     3    jmin,jmax,kimin,kimax,mask,iditfad,printflag,ic,jc,iprt,Qdir)
        write(iprt,*) 'after ITFACOMP: iregion,nitfa=',iregion,nitfa
        if(printflag>=2) then
          write(iprt,*) 'after itfacomp'
          do ki=1,nz
            write(iprt,*) 'i,j,k,itfadef=',ic,jc,ki,qitfad(ic,jc,ki)
          enddo
          ki=31
          do i=1,nx
            write(iprt,*) 'i,j,k,itfadef=',i,jc,ki,qitfad(i,jc,ki)
          enddo
        endif
      enddo  ! iregion loop
!
!     --- Merge the regions
      call MergeItfaRegions(qitfad,nx,ny,nz,imin,imax,jmin,jmax,
     1  kregion,minregion,maxregion,mask,printflag,ic,jc,iprt)
      if(printflag>=2) then
        write(iprt,*) 'after merge'
        do ki=1,nz
          write(iprt,*) 'i,j,k,itfadef=',ic,jc,ki,qitfad(ic,jc,ki)
        enddo
      endif
!
!     --- Write out ITFADEF to .Q file, and optionally netcdf file
      if(ioutputflag>0) then
        idQ=iditfad+399
        cnamei=cname(iditfad)
        call putqix(qitfad,nx,ny,nz,idQ,cnamei,printflag,iprt,Qdir,ierr)
      endif
      if(printflag>=1) then
        do ki=1,nz
          write(iprt,*) 'i,j,k,itfadef=',ic,jc,ki,qitfad(ic,jc,ki)
        enddo
        idQ=iditfad+399
        call TIstats(qitfad,nx,ny,nz,1,nx,1,ny,kimin,kimax,qmax,qmin,
     1    idQ,iprt)
        write(iprt,*)'exit  ITFA_static: nitfa=',nitfa
      endif
!
      return
      end

!-----------------------------------------------------------------------
  subroutine itfacompQ(nids,indxpicked,kregions,wts,cat,calmpL,clampH,qitfa)
!     --- Given a set of weights wts, computes the itfa combination
!     --- stored in qitfa(nx,ny,nz).  The individual indices are read 
!     --- in from the interpolated 4XX.Q files on disk in the Qdir.
!     --- qix and qs are work arrays.
!     --- The sum is only computed between i=imin,imax, j=jmin,jmax,
!     --- k=kmin,kmax.
!     --- nids is the output number of indices used.

    implicit none

    integer,intent(in) :: nids
    integer,intent(in) :: indxpicked(nids)
    integer,intent(in) :: kregions(IM,jsta:jend,MAXREGIONS,2)
    real,intent(in) :: wts(MAXREGIONS,IDMAX)
    real,intent(in) :: cat(IM,jsta:jend,LM,nids)
    real,intent(in) :: calmpL,clampH
    real,intent(inout) :: qitfa(IM,jsta:jend,LM)

    real :: qs(IM,jsta:jend,LM) ! work array

    integer :: idx
    integer :: n,i,j,k,idx,iregion,kk
    real :: wtindx,wqs,qijk

    write(*,*) 'enter itfacompQ'

!   --- Loop over all indices in the sum
    qitfa = SPVAL
    loop_n_idx: do n=1,nids
       idx = indxpicked(n)-399
       call remapi(idx,kregions,cat(1:IM,jsta:jend,1:LM,n),qs)

       ! --- Compute the weighted sum and store in qitfa
       do k=LM,1,-1
       do j=jsta,jend
       do i=1,IM

          ! --- Determine which k region
          iregion = -1
          do kk = 1, MAXREGIONS
             if (k <= kregions(i,j,kk,1) .and. k > kregions(i,j,kk,2)) then
                iregion = kk
                exit
             endif
          end do
          if (iregion < 0) cycle

          ! --- Get the right weight
          wtindx=wts(iregion,idx)

          if(ABS(qs(i,j,k)-SPVAL)<1.0E-3) cycle
          wqs=wtindx*MAX(qs(i,j,k),0.)
          qitfa(i,j,k)=qitfa(i,j,k)+wqs
       enddo
       enddo
       enddo
    end do loop_n_idx

!   --- Clamp the resultant sum between clampL and clampH
    do k=LM,1,-1
    do j=jsta,jend
    do i=1,IM
        qijk=qitfa(i,j,k)
        if(ABS(qijk-SPVAL)<1.0E-3) cycle
        qijk=MAX(qijk,clampL)
        qijk=MIN(qijk,clampH)
        qitfa(i,j,k)=qijk
     enddo
     enddo
     enddo

     return
   end subroutine itfacompQ

!-----------------------------------------------------------------------
  subroutine MergeItfaRegions(kregions,qitfa)
!     --- Merge the boundaries of the regions of multiple regions

    implicit none
     
    integer,intent(in) :: kregions(IM,jsta:jend,MAXREGIONS,2)
    real,intent(inout) :: qitfa(IM,jsta:jend,LM) 

    integer :: i,j,k,iregion,nregions
    real :: qijk,qsum,qk(LM)
    integer :: ksum
    integer :: kbdy,kbdy_m,kbdy_p
    write(*,*) 'enter MergeItfaRegions'


    do j=jsta,jend
    do i=1,IM

       nregions = 0
       do kk = 1, MAXREGIONS
          if(kregions(i,j,kk,1) > 0) nregions = nregions + 1
       end do
       if(nregions <= 1) return 

       do kk = 1, MAXREGIONS

          kbdy=kregions(i,j,kk,2)
          if(kbdy < 0) cycle

          do k=1,LM
             qk(k)=qitfa(i,j,k) ! save qitfa to qk because qitfa will be changed later
          enddo
          qsum=0.
          ksum=0
          kbdy_m=max(kbdy-1,1)
          kbdy_p=min(kbdy+1,LM)
          do k=kbdy_m,kbdy_p
             qijk=qk(k)
             if(ABS(qijk-SPVAL)>1.0E-3) then
                qsum=qsum+qijk
                ksum=ksum+1
             endif
          enddo
          if(ksum>0) then
             qitfa(i,j,kbdy)=qsum/FLOAT(ksum)
!            --- Merge in points above kbdy
             qsum=0.
             ksum=0
             kbdy_m=kbdy
             kbdy_p=min(kbdy+2,LM)
             do k=kbdy_m,kbdy_p
                qijk=qitfa(i,j,k)
                if(ABS(qijk-SPVAL)>1.0E-3) then
                   qsum=qsum+qijk
                   ksum=ksum+1
                endif
             enddo
             kbdy_p=min(kbdy+1,LM)
             if(ksum>0) then
                qitfa(i,j,kbdy_p)=qsum/FLOAT(ksum)
             endif
!            --- Merge in points below kbdy
             qsum=0.
             ksum=0
             kbdy_m=max(kbdy-2,1)
             kbdy_p=kbdy
             do k=kbdy_m,kbdy_p
                qijk=qitfa(i,j,k)
                if(ABS(qijk-SPVAL)>1.0E-3) then
                   qsum=qsum+qijk
                   ksum=ksum+1
                endif
             enddo
             if(ksum>0) then
                kbdy_m=max(kbdy-1,1)
                qitfa(i,j,kbdy_m)=qsum/FLOAT(ksum)
             endif
          endif
       enddo
    enddo
    enddo

    return
  end subroutine MergeItfaRegions

      subroutine itfamaxQ(qitfa,qitfam,qitfax,qit,iditfa,iditfam,
     1  iditfax,imin,imax,jmin,jmax,kimin,kimax,mbc,mask,
     2  nx,ny,nz,printflag,ic,jc,iprt,ioutputflag,Qdir)
!     --- Reads in ITFA (CAT) and ITFA (MWT) and takes max of two, 
!     --- grid point by grid point, and outputs as qitfax.
!     --- qit is a work array.
!     --- If ioutputflag > 0 qitfa is also stored on disk in directory Qdir as idQ.Q. 
!     --- On input ITFADEF or ITFADYN is in qitfa and ITFAMWT is in qitfam.
!     --- The max is only computed between i=imin,imax, j=jmin,jmax,
!     --- k=kmin,kmax and where mask(i,j)>0.
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nz
      real    qitfa(nx,ny,nz),qitfam(nx,ny,nz),qitfax(nx,ny,nz)
      real    qit(nx,ny,nz)
      integer iditfa,iditfam,iditfax
      integer imin,imax,jmin,jmax,kimin,kimax,mbc
      integer ic,jc,printflag,iprt
      real    qi,qm
      integer ioutputflag
      integer mask(nx,ny)
      character*200 Qdir
!-----------------------------------------------------------------------
      integer i,j,k,ki,idQ,jdQ,kdQ
      real    qijk,qmin,qmax
      integer Filttype,nftxy,nftz
      integer ierr
      character*24 cnamei
      include 'scorei.inc'
      include 'consts.inc'   ! SPVAL
!-----------------------------------------------------------------------
!
!     --- Initializations
      if(printflag>=1) then
        write(iprt,*) 'enter itfamaxQ: iditfa,iditfam,iditfax=',
     1    iditfa,iditfam,iditfax
        write(iprt,*) 'kimin,kimax=',kimin,kimax
        if(printflag>=2) then
          do k=1,nzi
            write(iprt,*) 'i,j,k,qitfa,qitfam=',ic,jc,k,qitfa(ic,jc,k),
     1        qitfam(ic,jc,k)
          enddo
        endif
      endif
      idQ=iditfa+399  ! input ITFA
      jdQ=iditfam+399 ! input ITFAMWT
      kdQ=iditfax+399 ! output ITFAMAX
      call TIinit(qitfax,nx,ny,nz)
!
!     --- Now obtain ITFAMAX=MAX(ITFA,ITFAMWT)
      qmin=+1.0E30
      qmax=-1.0E30
      do i=imin,imax
      do j=jmin,jmax
        do ki=kimin,kimax
          qitfax(i,j,ki)=SPVAL
          if(mask(i,j)<=0) go to 44
          qi=SPVAL
          qm=SPVAL
          if(ABS(qitfa(i,j,ki)-SPVAL)<1.0E-3) go to 44
          qi=qitfa(i,j,ki)  ! CAT index
          qitfa(i,j,ki)=qi
          if(ABS(qitfam(i,j,ki)-SPVAL)<1.0E-3) go to 44
          qm=qitfam(i,j,ki)  ! MWT index
          qitfax(i,j,ki)=MAX(qi,qm)
   44     continue
          if(printflag>=1 .and. i==ic .and. j==jc) then
            write(iprt,*) 'i,j,ki,z,ITFA,ITFAMWT,MAX=',i,j,ki,
     1       zi(ki),qi,qm,qitfax(i,j,ki)
          endif
        enddo
      enddo
      enddo
!
!     --- Clamp the resultant sum between clampL and clampH
      qmin=+1.0E30
      qmax=-1.0E30
      do ki=kimin,kimax
      do j=jmin,jmax
      do i=imin,imax
        qijk=qitfax(i,j,ki)
        if(ABS(qijk-SPVAL)<1.0E-3) go to 45
        qijk=MAX(qijk,clampitfaL)
        qijk=MIN(qijk,clampitfaH)
        qitfax(i,j,ki)=qijk
        qmin=MIN(qmin,qijk)
        qmax=MAX(qmax,qijk)
   45   continue
      enddo
      enddo
      enddo
!
!     --- Merge the regions
      call MergeItfaRegions(qitfax,nx,ny,nz,imin,imax,jmin,jmax,
     1  kregion,minregion,maxregion,mask,printflag,ic,jc,iprt)
      if(printflag>=2) then
        write(iprt,*) 'after merge'
        do ki=1,nzi
          write(iprt,*) 'i,j,k,itfamax=',ic,jc,ki,qitfax(ic,jc,ki)
        enddo
      endif
!
!     --- Perform one smoothing for the blend, and for consistency
!     --- also smooth ITFA
      Filttype=1
      nftxy=1
      nftz=1
      call filt3d(qitfax,qit,nx,ny,nzi,imin,imax,jmin,jmax,
     1  kimin,kimax,mbc,nftxy,nftz,Filttype)
!     call filt3d(qitfa ,qit,nx,ny,nzi,imin,imax,jmin,jmax,
!    1  kimin,kimax,mbc,nftxy,nftz,Filttype)
!     call filt3d(qitfam,qit,nx,ny,nzi,imin,imax,jmin,jmax,
!    1  kimin,kimax,mbc,nftxy,nftz,Filttype)
      if(printflag>=1) then
        do ki=1,nzi
          write(iprt,*)'i,j,k,itfamax smooth=',ic,jc,ki,qitfax(ic,jc,ki)
        enddo
        do ki=1,nzi
          write(iprt,*)'i,j,k,itfamwt smooth=',ic,jc,ki,qitfam(ic,jc,ki)
        enddo
      endif
      if(ioutputflag>0) then
!       --- Write out ITFAMAX to .Q file
        idQ=iditfax+399
        cnamei=cname(iditfax)
        call putqix(qitfax,nx,ny,nz,idQ,cnamei,printflag,iprt,Qdir,ierr)
!       idQ=iditfa+399
!       cnamei=cname(iditfa)
!       call putqix(qitfa ,nx,ny,nz,idQ,cnamei,printflag,iprt,Qdir,ierr)
!       idQ=iditfam+399
!       cnamei=cname(iditfam)
!       call putqix(qitfam,nx,ny,nz,idQ,cnamei,printflag,iprt,Qdir,ierr)
      endif
!
      if(printflag>=1) then
        do ki=1,nzi
          write(iprt,*) 'i,j,k,z,itfamaxQ=',ic,jc,ki,zi(ki),
     1     qitfax(ic,jc,ki)
        enddo
        write(iprt,*) 'exit  itfamaxQ: qmin,qmax=',qmin,qmax
      endif
      idQ=iditfax+399
      call TIstats(qitfax,nx,ny,nzi,1,nx,1,ny,kimin,kimax,qmax,qmin,idQ,
     1  iprt)
!
      return
      end

!-----------------------------------------------------------------------
  subroutine remapi(idx,kregions,q,qs)
!     --- Performs mapping of input variable q into output 
!     --- variable qs (0-1).  Input index specific thresholds
!     --- (null,light,moderate,severe) are in timap.  Corresponding
!     --- (0-1) thresholds are in tis. remap_option is in scorei.inc
!     --- If remap_option=1 use piecewise linear mapping
!     --- If remap_option=2 use fit to log-normal PDF

    implicit none

    integer,intent(in) :: idx
    integer,intent(in) :: kregions(IM,jsta:jend,MAXREGIONS,2)
    real,intent(in) :: q(IM,jsta:jend,LM)
    real,intent(inout) :: qs(IM,jsta:jend,LM)

    integer :: i,j,k,idQ,iregion,kk
    real :: qijk,qsijk
    real :: LS,ceps,fact,e32,epsilon

    idQ=idx+399

!   --- Do the remap according to remap_option

    do k=LM,1,-1
    do j=jsta,jend
    do i=1,IM

!      --- Determine which k region
       iregion = -1
       do kk = 1, MAXREGIONS
          if (k <= kregions(i,j,kk,1) .and. k > kregions(i,j,kk,2)) then
             iregion = kk
             exit
          endif
       end do
       if (iregion < 0) cycle

       qijk=q(i,j,k)
       if(ABS(qijk-SPVAL)<=1.0E-3) cycle
       qs(i,j,k)=0.
       if(idQ==465) then ! SGSTKE
!         --- Map according to epsilon=(C/L)e^3/2 with L~dx
!         LS=dxm  ! assume length scale ~ horizontal grid spacing
          LS=100. ! assumed vertical length scale
!         LS= 20. ! assumed vertical length scale
          ceps = 0.845
!         fact=5. ! multiply by a factor to get better agreement with observations
!         fact=10.! multiply by a factor to get better agreement with observations
!         fact=100.! multiply by a factor to get better agreement with observations
          fact=3.  ! multiply by a factor to get better agreement with observations
          ceps=fact*ceps
          e32 = qijk**1.5
          epsilon=(ceps/LS)*e32
          qsijk=epsilon**(1./3.)
       else
          call remapq(iregion,idx,qijk,qsijk)
       endif
       qs(i,j,k)=qsijk
    enddo
    enddo
    enddo

    return
  end subroutine remapi

!-----------------------------------------------------------------------
  subroutine remapq(iregion,idx,q,qs)
!     --- Performs mapping of input variable q into output 
!     --- variable qs (0-1).  Input index specific thresholds
!     --- (null,light,moderate,severe) are in timap.  Corresponding
!     --- (0-1) thresholds are in tis.
!     --- If remap_option=1 use piecewise linear mapping
!     --- If remap_option=2 use fit to log-normal PDF

    implicit none

    integer,intent(in) :: iregion,idx
    real,intent(in) :: q
    real,intent(out) :: qs

    real :: tii(NTI), A(2)

    qs=SPVAL
    if(ABS(q-SPVAL)<=1.0E-3) return

    if(remap_option==2) then
!   --- Use fit of indices to log-normal PDF
!   --- Log(epsilon^(1/3)) = a + b Log(I)
!   --- lnedrfits is saved in timap(1:2)
       A(1)=timap(iregion,idx,1)
       A(2)=timap(iregion,idx,2)
       call remap2(A,q,qs)
    else
!   --- Default: use piecewise linear remap of entire grid
       tii(1:NTI)=timap(iregion,idx,1:NTI)
       call remap1(NTI,tii,tis,q,qs)
    endif

!   --- clamp the remapped index between clampL and clampH
    qs=MAX(qs,clampidxL)
    qs=MIN(qs,clampidxH)

    return
  end subroutine remapq

!-----------------------------------------------------------------------
  subroutine remap1(n,ti,tis,q,qs)
!     --- Performs linear remap of index value contained in q into
!     --- scaled value in the range (0-1)

    implicit none

    integer, intent(in) :: n
    real, intent(in) :: ti(n),tis(n)
    real, intent(in) :: q
    real, intent(out) :: qs
    real :: slope

    qs=SPVAL
    if(q<ti(2)) then
       slope=(tis(2)-tis(1))/(ti(2)-ti(1))
       qs=tis(1)+slope*(q-ti(1))
    elseif(q<ti(3)) then
       slope=(tis(3)-tis(2))/(ti(3)-ti(2))
       qs=tis(2)+slope*(q-ti(2))
    elseif(q<ti(4)) then
       slope=(tis(4)-tis(3))/(ti(4)-ti(3))
       qs=tis(3)+slope*(q-ti(3))
    else
       slope=(tis(5)-tis(4))/(ti(5)-ti(4))
       qs=tis(4)+slope*(q-ti(4))
    endif
    return
  end subroutine remap1

!-----------------------------------------------------------------------
  subroutine remap2(A,q,qs)
!     --- Fit q to log-normal PDF, output is epsilon^1/3=qs
!     --- Log(epsilon^(1/3)) = a + b Log(I), where a=A(1), b=A(2)

    implicit none

    real,intent(in) :: q
    real,intent(in) :: A(2)
    real,intent(out) :: qs

    real :: ai,bi,qq,logqs

!   --- use fit of indices to log-normal PDF
!   --- Log(epsilon^(1/3)) = a + b Log(I)

    qs=0.
    if(q<1.0E-20) return
    ai=A(1)
    bi=A(2)
    qq=MAX(q,1.0E-20)  ! protect against 0 or neg values
    logqs = ai + bi*ALOG(qq)
    qs = EXP(logqs)

    return
  end subroutine remap2
