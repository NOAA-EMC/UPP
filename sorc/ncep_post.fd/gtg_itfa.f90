!
      subroutine ITFA_MWT(qix,qit,qitfam,iditfam,nx,ny,nz,imin,imax,
     1  jmin,jmax,kimin,kimax,mask,printflag,ic,jc,iprt,ioutputflag,
     2  Qdir)
!     --- Computes an ITFA combination for MWT and outputs as qitfa
!     --- If ioutputflag > 0 qitfa is also stored on disk in directory Qdir as idQ.Q. 
!     --- qix and qit are work arrays.
!     --- The sum is only computed between i=imin,imax, j=jmin,jmax,
!     --- k=kmin,kmax.
!     --- The individual indices are read in from the interpolated 4XX.Q
!     --- files on disk in the Qdir.
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nz
      real    qix(nx,ny,nz),qit(nx,ny,nz),qitfam(nx,ny,nz)
      integer iditfam
      integer nwts
      integer imin,imax,jmin,jmax,kimin,kimax
      integer ioutputflag
      integer ic,jc,printflag,iprt
      integer mask(nx,ny)
      character*200 Qdir
!-----------------------------------------------------------------------
      integer i,j,ki,idx,idQ,ierr
      integer iregion
      real    qmax,qmin
      real    wsum
      character*24 cnamei
      integer nimax
      parameter (nimax=101)
      integer kpickitfa(nimax)
      include 'scoreparams.inc'
      include 'scorei.inc'
      include 'pireps.inc'   ! npireps
      include 'consts.inc'   ! RMISSD
!-----------------------------------------------------------------------
!
!     --- Initializations
      if(printflag.ge.1) then
        write(iprt,*) 'enter ITFA_MWT: iditfa,remap_option=',
     1    iditfam,remap_option
        write(iprt,*) 'minregion,maxregion=',minregion,maxregion
        do iregion=minregion,maxregion
          do idx=1,idmax
            if(ifcsttype(idx).eq.MWT .and. ipickitfa(iregion,idx).gt.0)
     1      then
              write(iprt,*)'iregion,idQ,ipickitfa=',iregion,idx+399,
     1          ipickitfa(iregion,idx)
            endif
          enddo
        enddo
      endif
      call TIinit(qix,nx,ny,nz)
      call TIinit(qit,nx,ny,nz)
!     --- Initialize ITFAMWT to 0.
      do ki=1,nz
      do j=1,ny
      do i=1,nx
        qitfam(i,j,ki)=0.
      enddo
      enddo
      enddo
!
!     --- Perform the combination using altitude-dependent static weights
      do iregion=minregion,maxregion
!       --- Compute the ITFAMWT combination for this region
        write(iprt,*) 'computing ITFAMWT for iregion=',iregion
        nwts=0
        wsum=0.
        do idx=1,idmax
          kpickitfa(idx)=0
          if(ifcsttype(idx).eq.MWT .and. ipickitfa(iregion,idx).gt.0)
     1    then
            kpickitfa(idx)=ipickitfa(iregion,idx)
            nwts=nwts+1
            wsum=wsum+static_wts(iregion,idx)
            if(printflag.ge.2) then
              write(iprt,*) 'iregion,idQ,kpickitfa=',iregion,idx+399,
     1         kpickitfa(idx)
            endif
          endif
        enddo
!       --- Check the weights
        write(iprt,*) 'nwts,wsum=',nwts,wsum
        if(nwts.le.0) then
          write(iprt,*) 'No weights specificied for MWT region=',iregion
          go to 153
        endif
        nwts=0
!       --- qix and qit are work arrays.  Output is in qitfam
        call itfacompQ(qix,qit,qitfam,static_wts,kpickitfa,
     1    remap_option,iregion,kregion,minregion,maxregion,
     2    clampitfaL,clampitfaH,nx,ny,nz,idmax,nwts,imin,imax,
     3    jmin,jmax,kimin,kimax,mask,iditfam,printflag,ic,jc,iprt,Qdir)
        write(iprt,*) 'after ITFACOMP: iregion,nitfa=',iregion,nwts
        if(printflag.ge.2) then
          write(iprt,*) 'after itfacomp'
          do ki=1,nz
            write(iprt,*) 'i,j,k,itfamwt=',ic,jc,ki,qitfam(ic,jc,ki)
          enddo
!         ki=31
!         do i=1,nx
!           write(iprt,*) 'i,j,k,itfamwt=',i,jc,ki,qitfam(i,jc,ki)
!         enddo
        endif
  153   continue
      enddo  ! iregion loop
!     --- Merge the regions
      call MergeItfaRegions(qitfam,nx,ny,nz,imin,imax,jmin,jmax,
     1  kregion,minregion,maxregion,mask,printflag,ic,jc,iprt)
      if(printflag.ge.2) then
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
      if(ioutputflag.gt.0) then
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
      integer nwts
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
      include 'consts.inc'   ! RMISSD
!-----------------------------------------------------------------------
!
!     --- Initializations
      if(printflag.ge.1) then
        write(iprt,*)'enter ITFA_static: idQ,remap_option=',iditfad+399,
     1    remap_option
        write(iprt,*) 'minregion,maxregion=',minregion,maxregion
        do iregion=minregion,maxregion
          do idx=1,idmax
            if(ipickitfa(iregion,idx).gt.0) then
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
          if(ifcsttype(idx).eq.CAT)kpickitfa(idx)=ipickitfa(iregion,idx)
        enddo
        if(printflag.ge.1) then
          ni=0
          do idx=1,idmax
            if(kpickitfa(idx).gt.0) then
              ni=ni+1
              write(iprt,*) 'iregion,idQ,kpickitfa,wt=',iregion,idx+399,
     1         kpickitfa(idx),static_wts(iregion,idx)
            endif
          enddo
          write(iprt,*) 'nindices to use=',ni
        endif
        nwts=0
!       --- qix and qit are work arrays.  Output is in qitfad
        call itfacompQ(qix,qit,qitfad,static_wts,kpickitfa,
     1    remap_option,iregion,kregion,minregion,maxregion,
     2    clampitfaL,clampitfaH,nx,ny,nz,idmax,nwts,imin,imax,
     3    jmin,jmax,kimin,kimax,mask,iditfad,printflag,ic,jc,iprt,Qdir)
        write(iprt,*) 'after ITFACOMP: iregion,nitfa=',iregion,nwts
        if(printflag.ge.2) then
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
      if(printflag.ge.2) then
        write(iprt,*) 'after merge'
        do ki=1,nz
          write(iprt,*) 'i,j,k,itfadef=',ic,jc,ki,qitfad(ic,jc,ki)
        enddo
      endif
!
!     --- Write out ITFADEF to .Q file, and optionally netcdf file
      if(ioutputflag.gt.0) then
        idQ=iditfad+399
        cnamei=cname(iditfad)
        call putqix(qitfad,nx,ny,nz,idQ,cnamei,printflag,iprt,Qdir,ierr)
      endif
      if(printflag.ge.1) then
        do ki=1,nz
          write(iprt,*) 'i,j,k,itfadef=',ic,jc,ki,qitfad(ic,jc,ki)
        enddo
        idQ=iditfad+399
        call TIstats(qitfad,nx,ny,nz,1,nx,1,ny,kimin,kimax,qmax,qmin,
     1    idQ,iprt)
        write(iprt,*)'exit  ITFA_static: nitfa=',nwts
      endif
!
      return
      end
!
      subroutine itfacompQ(qix,qs,qitfa,wts,ipickitfa,
     1  remap_option,iregion,kregion,minregion,maxregion,
     2  clampitfaL,clampitfaH,nx,ny,nz,idmax,nitfa,imin,imax,
     3  jmin,jmax,kmin,kmax,mask,iditfa,printflag,ic,jc,iprt,Qdir)
!     --- Given a set of weights wts, computes the itfa combination
!     --- stored in qitfa(nx,ny,nz).  The individual indices are read 
!     --- in from the interpolated 4XX.Q files on disk in the Qdir.
!     --- qix and qs are work arrays.
!     --- The sum is only computed between i=imin,imax, j=jmin,jmax,
!     --- k=kmin,kmax.
!     --- nitfa is the output number of indices used.
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nz,idmax
      real    qix(nx,ny,nz),qs(nx,ny,nz),qitfa(nx,ny,nz)
      real    wts(3,idmax)    ! 3=low,mid,high
      integer kregion(3,2)
      integer iregion,minregion,maxregion
      integer ipickitfa(idmax)
      integer remap_option
      real    clampitfaL,clampitfaH
      integer iditfa
      integer nitfa
      integer imin,imax,jmin,jmax,kmin,kmax,ic,jc,printflag,iprt
      character*200 Qdir
      integer mask(nx,ny)
!-----------------------------------------------------------------------
      integer i,j,k,idx,idQ,ierr,ni
      real    qijk,qmin,qmax,indxwt,wqs,qitfalast
      integer kmaxi,kmini
      integer krmin,krmax
      integer nregions
      include 'consts.inc'   ! RMISSD
!-----------------------------------------------------------------------
!
!     --- Initializations
      if(printflag.ge.1) then
        idQ=iditfa+399
        write(iprt,*) 'enter itfacompQ: idQ,remap_option=',
     1    idQ,remap_option
        write(iprt,*) 'iregion,minregion,maxregion=',iregion,minregion,
     1    maxregion
        if(printflag.ge.2) then
          ni=0
          do idx=1,idmax
            indxwt=wts(iregion,idx)
            if(indxwt.gt.0.) then
              ni=ni+1
              write(iprt,*) 'iregion,idx,ipickitfa,wts=',iregion,idx,
     1        ipickitfa(idx),indxwt
            endif
          enddo
          write(iprt,*) 'iregion,nitfa=',iregion,ni
        endif
      endif
!
!     --- Initializations
      nregions=maxregion-minregion+1
      krmin=kregion(iregion,1)
      krmax=kregion(iregion,2)
      if(nregions.eq.1) then
        kmini=kmin
        kmaxi=kmax
      else
        kmini=krmin
        kmaxi=krmax
      endif
      if(printflag.ge.2) then
        write(iprt,*) 'iregion,kmini,kmaxi=',iregion,kmini,kmaxi
      endif
!
!     --- Loop over all indices in the sum
      nitfa=0
      do idx=1,idmax
        if(ipickitfa(idx).le.0) go to 33
        idQ=idx+399
!       --- Retrieve the index 
        call getqix(qix,nx,ny,nz,idQ,printflag,iprt,Qdir,ierr)
        if(ierr.ne.0) then
          write(iprt,*) 'error in itfacompQ reading index idQ=',idQ
          go to 33
        endif
        if(printflag.ge.2) then
          write(iprt,*) 'retrieving qix iregion,idQ=',iregion,idQ
          do k=1,nz
            write(iprt,*) 'i,j,k,qix=',ic,jc,k,qix(ic,jc,k)
          enddo
        endif
!       --- Remap index qix to qs
        call TIinit(qs,nx,ny,nz)
        call remapi(qix,qs,nx,ny,nz,iregion,idx,
     1    imin,imax,jmin,jmax,kmini,kmaxi,mask,printflag,ic,jc,iprt)
!       --- Compute the weighted sum and store in qitfa for the
!       --- current region
        indxwt=wts(iregion,idx)
        if(printflag.ge.2) then
          write(iprt,*) 'iregion,idQ,wts=',iregion,idQ,indxwt
        endif
        qmin=1.0E30
        qmax=-1.0E30
        do i=imin,imax
        do j=jmin,jmax
          if(mask(i,j).le.0) go to 32
        do k=kmini,kmaxi
           if(nitfa.le.0) then
            qitfalast=0.
          else
            qitfalast=qitfa(i,j,k)
          endif
          wqs=RMISSD
          if(ABS(qitfalast-RMISSD).LT.1.0E-3) go to 30
          if(ABS(qs(i,j,k)-RMISSD).LT.1.0E-3) go to 30
          wqs=indxwt*MAX(qs(i,j,k),0.)
          qitfa(i,j,k)=qitfalast+wqs
          qmin=MIN(qmin,qitfa(i,j,k))
          qmax=MAX(qmax,qitfa(i,j,k))
   30     continue
          if(printflag.ge.2 .and. i.eq.ic .and. j.eq.jc) then
            write(iprt,*) 'idQ,i,j,k,qiftalast,qs,wqs,qitfa=',
     1       idQ,ic,jc,k,qitfalast,qs(ic,jc,k),wqs,qitfa(ic,jc,k)
          endif
        enddo
   32   continue   
        enddo
        enddo
!       --- Increment the number of indices used
        nitfa=nitfa+1
   33   continue
      enddo  ! idx loop
!
!     --- Clamp the resultant sum between clampL and clampH
      qmin=+1.0E30
      qmax=-1.0E30
      do k=kmini,kmaxi
      do j=jmin,jmax
      do i=imin,imax
        qijk=qitfa(i,j,k)
        if(ABS(qijk-RMISSD).LT.1.0E-3) go to 35
        qijk=MAX(qijk,clampitfaL)
        qijk=MIN(qijk,clampitfaH)
        qitfa(i,j,k)=qijk
        qmin=MIN(qmin,qijk)
        qmax=MAX(qmax,qijk)
   35   continue
      enddo
      enddo
      enddo
!
      if(printflag.ge.1) then
        write(iprt,*)  'exit itfacompQ: iregion=',iregion
        if(printflag.ge.2) then
          do k=1,nz
            write(iprt,*)  'i,j,k,qitfa=',ic,jc,k,qitfa(ic,jc,k)
          enddo
        endif
        write(iprt,*) 'exit itfacompQ: nitfa=',nitfa
        write(iprt,*) 'exit itfacompQ: clamped qitfa min,max=',qmin,qmax
      endif
!
      return
      end
!
      subroutine MergeItfaRegions(qitfa,nx,ny,nzi,imin,imax,jmin,jmax,
     1  kregion,minregion,maxregion,mask,printflag,ic,jc,iprt)
!     --- Merge the boundaries of the regions of multiple regions
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nzi
      real    qitfa(nx,ny,nzi)
      integer kregion(3,2)
      integer minregion,maxregion
      integer ic,jc,printflag,iprt
      integer imin,imax,jmin,jmax
      integer mask(nx,ny)
      integer i,j,k,iregion,nregions
      real    qijk,qsum
      integer ksum
      integer nzmax
      parameter (nzmax=201)
      real    qk(nzmax)
      integer kbdy
      include 'consts.inc'   ! RMISSD
!-----------------------------------------------------------------------
!
      nregions=maxregion-minregion+1
      if(printflag.ge.1) then
        write(iprt,*) 'enter MergeItfaRegions: nregions=',nregions
        if(printflag.ge.2) then
          do k=1,nzi
            write(iprt,*) 'i,j,k,inputqitfa=',ic,jc,k,qitfa(ic,jc,k)
          enddo
        endif
      endif
!
      if(nregions.le.1) return
        do iregion=minregion,maxregion-1
          kbdy=kregion(iregion,2)
          if(printflag.ge.2) then
            write(iprt,*) 'iregion,kbdy=',iregion,kbdy
          endif
          do j=jmin,jmax
          do i=imin,imax
            if(mask(i,j).le.0) go to 10
            do k=1,nzi
              qk(k)=qitfa(i,j,k)
            enddo
            qsum=0.
            ksum=0
            do k=kbdy-1,kbdy+1
              qijk=qk(k)
              if(ABS(qijk-RMISSD).GT.1.0E-3) then
                qsum=qsum+qijk
                ksum=ksum+1
              endif
              if(printflag.ge.2 .and. i.eq.ic .and. j.eq.jc) then
                write(iprt,*) 'iregion,i,j,k,qitfa,qsum,ksum=',
     1            iregion,ic,jc,k,qijk,qsum,ksum
              endif
            enddo
            if(ksum.gt.0) then
              qitfa(i,j,kbdy)=qsum/FLOAT(ksum)
!             --- Merge in points above kbdy
              qsum=0.
              ksum=0
              do k=kbdy,kbdy+2
                qijk=qitfa(i,j,k)
                if(ABS(qijk-RMISSD).GT.1.0E-3) then
                  qsum=qsum+qijk
                  ksum=ksum+1
                endif
              enddo
              if(ksum.gt.0) then
                qitfa(i,j,kbdy+1)=qsum/FLOAT(ksum)
              endif
!             --- Merge in points below kbdy
              qsum=0.
              ksum=0
              do k=kbdy,kbdy-2
                qijk=qitfa(i,j,k)
                if(ABS(qijk-RMISSD).GT.1.0E-3) then
                  qsum=qsum+qijk
                  ksum=ksum+1
                endif
              enddo
              if(ksum.gt.0) then
                qitfa(i,j,kbdy-1)=qsum/FLOAT(ksum)
              endif
              if(printflag.ge.3 .and. i.eq.ic .and. j.eq.jc) then
                do k=1,nzi
                  write(iprt,*) 'iregion,i,j,k,qk,qitfa=',
     1              iregion,ic,jc,k,qk(k),qitfa(i,j,k)
                enddo
              endif
          endif
   10     continue
          enddo
          enddo
        enddo
!
      if(printflag.ge.1) then
        write(iprt,*) 'exit MergeItfaRegions'
        if(printflag.ge.2) then
          do k=1,nzi
            write(iprt,*) 'i,j,k,mergedqitfa=',ic,jc,k,qitfa(ic,jc,k)
          enddo
        endif
      endif
!
      return
      end
!
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
      include 'consts.inc'   ! RMISSD
!-----------------------------------------------------------------------
!
!     --- Initializations
      if(printflag.ge.1) then
        write(iprt,*) 'enter itfamaxQ: iditfa,iditfam,iditfax=',
     1    iditfa,iditfam,iditfax
        write(iprt,*) 'kimin,kimax=',kimin,kimax
        if(printflag.ge.2) then
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
          qitfax(i,j,ki)=RMISSD
          if(mask(i,j).le.0) go to 44
          qi=RMISSD
          qm=RMISSD
          if(ABS(qitfa(i,j,ki)-RMISSD).lt.1.0E-3) go to 44
          qi=qitfa(i,j,ki)  ! CAT index
          qitfa(i,j,ki)=qi
          if(ABS(qitfam(i,j,ki)-RMISSD).lt.1.0E-3) go to 44
          qm=qitfam(i,j,ki)  ! MWT index
          qitfax(i,j,ki)=MAX(qi,qm)
   44     continue
          if(printflag.ge.1 .and. i.eq.ic .and. j.eq.jc) then
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
        if(ABS(qijk-RMISSD).LT.1.0E-3) go to 45
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
      if(printflag.ge.2) then
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
      if(printflag.ge.1) then
        do ki=1,nzi
          write(iprt,*)'i,j,k,itfamax smooth=',ic,jc,ki,qitfax(ic,jc,ki)
        enddo
        do ki=1,nzi
          write(iprt,*)'i,j,k,itfamwt smooth=',ic,jc,ki,qitfam(ic,jc,ki)
        enddo
      endif
      if(ioutputflag.gt.0) then
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
      if(printflag.ge.1) then
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
