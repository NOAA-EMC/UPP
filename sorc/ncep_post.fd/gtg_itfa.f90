module gtg_itfa
  use ctlblk_mod, only: jsta,jend, IM,JM,LM, SPVAL
  use gtg_config, only : static_wts,remap_option
  use gtg_config, only : clampidxL,clampidxH,clampitfaL,clampitfaH
  use gtg_filter

  implicit none
contains

  subroutine ITFAcompF(ncat,ipickitfa,kregions,cat,comp_ITFAMWT,comp_ITFADYN,qitfax)
! Computes ITFA combinations as the weighted sum of turbulence indices.
! static_wgt(MAXREGIONS, IDMAX) is from gtg_config

! idQ=491 for itfa CAT combination using dynamic weights
! idQ=492 for itfa CAT combination using static weights
! idQ=493 for itfa MWT combination using static weights
! idQ=494 for MAX(CAT, MWT) itfa combinations using static weights
! idQ=495 for MAX(CAT, MWT) itfa combinations using dynamic weights

    implicit none

    integer,intent(in) :: ncat
    integer,intent(in) :: ipickitfa(MAXREGIONS,ncat)
    integer,intent(in) :: kregions(MAXREGIONS,2)
    real,intent(in) :: cat(IM,jsta:jend,LM,ncat)
    logical, intent(in) :: comp_ITFAMWT,comp_ITFADYN
    real,intent(inout) :: qitfax(IM,jsta:jend,LM)

    ! for one kregion, save for MWT/CAT
    integer :: kpickitfa(ncat)
    real :: wts(ncat)

    integer :: iregion,i

!   qitfa=Output itfa CAT combination using static/dynamic weights
!   qitfam=Output itfa MWT combination using static weights
!   qitfax=Output MAX(CAT, MWT) itfa combinations
    real :: qitfa(IM,jsta:jend,LM), qitfam(IM,jsta:jend,LM)

    write(*,*) 'enter ITFAcompF'

    qitfa = 0.
    qitfam = 0.
    qitfax = 0.

    loop_iregion: do iregion=1,MAXREGIONS

       if(comp_ITFADYN) then 
          ! Comute an ITFA based on dynamic weights
          ! not available for current version
          qitfa = SPVAL
       else
          kpickitfa = 0
          do i = 1, ncat
             if(ipickitfa(iregion,i)<=475 .and. ipickitfa(iregion,i)>0) then ! CAT
                kpickitfa(i) = ipickitfa(iregion,i)
                wts(i)=static_wts(iregion,i)
             end if
          end do
          if (sum(kpickitfa) == 0) then
             write(*,*) "There is no CAT static indices picked"
             return
          end if

          ! Compute an ITFA combination using a set of default weight
          call itfasum()
       end if

       if(comp_ITFAMWT) then
          kpickitfa = 0
          do i = 1, ncat
             if(ipickitfa(iregion,i) >= 476 .and. ipickitfa(iregion,i) <= 490) then ! MWT
                kpickitfa(i) = ipickitfa(iregion,i)
                wts(i)=static_wts(iregion,i)
             end if
          end do
          if (sum(kpickitfa) == 0) then
             write(*,*) "There is no MWT indices picked"
             return
          end if

          ! Compute an ITFA combination using a set of default wei

          call itfasum()
       else
          qitfam = 0.
       endif

       ! Now obtain ITFAMAX=MAX(ITFA,ITFAMWT)
       call itfamax()

       write(iprt,*)'exit  ITFAcompF'

    end do loop_idx
    end do loop_iregion
    return
  end subroutine ITFAcompF


!-----------------------------------------------------------------------
  subroutine ITFA_wts(ncat,kpickitfa,kregions,wts,cat,miss,qitfa)
    implicit none

    integer,intent(in) :: ncat
    integer,intent(in) :: kpickitfa(MAXREGIONS,ncat)
    integer,intent(in) :: kregions(MAXREGIONS,2)
    real, intent(in) :: wts(MAXREGIONS, IDMAX)
    real,intent(in) :: cat(IM,jsta:jend,LM,ncat)
    real, intent(in) :: miss
    real,intent(inout) :: qitfaa(IM,jsta:jend,LM)

    integer :: i, iregion

    write(*,*) 'enter ITFA_wts'

    qitfa = miss


    loop_iregion: do iregion = 1, MAXREGIONS
      ! --- Compute the ITFAMWT combination for this region 
       nitfa = 0
       kpickitfa = 0
       do i = 1, ncat
          if(ipickitfa(iregion,i) >= 476 .and. ipickitfa(iregion,i) <= 490) then ! MWT
             nitfa = nitfa + 1
             kpickitfa(nitfa) = ipickitfa(iregion,i)
             nids_itfa(nitfa) = i ! new array for index of cat(:,:,:,n) 
          end if
       end do
       if (nitfa <= 0) then
          write(*,*) "There is no MWT indices picked"
          return
       end if

!  Perform the combination using altitude-dependent static weights
       call itfacompQ(nitfa,kpickitfa(1:nitfa),nids_itfa(1:nitfa),kregions,static_wgt,&
                      cat,clampitfaL,clampitfaH,qitfam)

       write(*,*) "ITFA_wts: iregion,n,kpickitfa(:)=",iregion,nitfa,kpickitfa(1:nitfa),qitfam(IM/2,jsta,1:LM)

    end do loop_iregion

!   --- Merge the regions
    call MergeRegions(kregions,qitfam)

    return
  end subroutine ITFA_wts

!-----------------------------------------------------------------------
  subroutine ITFA_MWT(ncat,ipickitfa,kregions,cat,qitfam)
! Computes an ITFA combination for MWT and outputs as qitfa
! If ioutputflag > 0 qitfa is also stored on disk in directory Qdir as idQ.Q. 
! qix and qit are work arrays.
! The sum is only computed between i=imin,imax, j=jmin,jmax,
! k=kmin,kmax.
! The individual indices are read in from the interpolated 4XX.Q
! files on disk in the Qdir.

    implicit none

    integer,intent(in) :: ncat
    integer,intent(in) :: ipickitfa(MAXREGIONS,ncat)
    integer,intent(in) :: kregions(MAXREGIONS,2)
    real,intent(in) :: cat(IM,jsta:jend,LM,ncat)
    real,intent(inout) :: qitfam(IM,jsta:jend,LM)

    integer :: kpickitfa(ncat),nids_itfa(ncat)
    integer :: nitfa,i
    integer :: iregion

    write(*,*) 'enter ITFA_MWT'

    qitfam = 0.


    loop_iregion: do iregion = 1, MAXREGIONS
      ! --- Compute the ITFAMWT combination for this region 
       nitfa = 0
       kpickitfa = 0
       do i = 1, ncat
          if(ipickitfa(iregion,i) >= 476 .and. ipickitfa(iregion,i) <= 490) then ! MWT
             nitfa = nitfa + 1
             kpickitfa(nitfa) = ipickitfa(iregion,i)
             nids_itfa(nitfa) = i ! new array for index of cat(:,:,:,n) 
          end if
       end do
       if (nitfa <= 0) then
          write(*,*) "There is no MWT indices picked"
          return
       end if

!  Perform the combination using altitude-dependent static weights
       call itfacompQ(nitfa,kpickitfa(1:nitfa),nids_itfa(1:nitfa),kregions,static_wgt,&
                      cat,clampitfaL,clampitfaH,qitfam)

       write(*,*) "ITFA_MWT: iregion,n,kpickitfa(:)=",iregion,nitfa,kpickitfa(1:nitfa),qitfam(IM/2,jsta,1:LM)

    end do loop_iregion

!   --- Merge the regions
    call MergeRegions(kregions,qitfam)

    return
  end subroutine ITFA_MWT


!-----------------------------------------------------------------------
  subroutine ITFA_static(ncat,ipickitfa,kregions,cat,qitfad)
! Computes an ITFA combination using static weights and outputs as qitfa.

    implicit none

    integer,intent(in) :: ncat
    integer,intent(in) :: ipickitfa(MAXREGIONS,ncat)
    integer,intent(in) :: kregions(MAXREGIONS,2)
    real,intent(in) :: cat(IM,jsta:jend,LM,ncat)
    real,intent(inout) :: qitfad(IM,jsta:jend,LM)

    integer :: kpickitfa(ncat),nids_itfa(ncat)
    integer :: nitfa,i
    integer :: iregion

    write(*,*) 'enter ITFA_static'

    qitfad = 0.

    loop_iregion: do iregion = 1, MAXREGIONS
       nitfa = 0
       kpickitfa = 0
       do i = 1, ncat
          if(ipickitfa(iregion,i) <= 475 .and. ipickitfa(iregion,i) > 0) then ! CAT
             nitfa = nitfa + 1
             kpickitfa(nitfa) = ipickitfa(iregion,i)
             nids_itfa(nitfa) = i ! new array for index of cat(:,:,:,n) 
          end if
       end do
       if (nitfa <= 0) then
          write(*,*) "There is no ITFA_static indices picked"
          return
       end if

!  Perform the combination using altitude-dependent static weights
       call itfacompQ(nitfa,kpickitfa(1:nitfa),nids_itfa(1:nitfa),kregions,static_wgt,&
                      cat,clampitfaL,clampitfaH,qitfad)
    end do loop_iregion

!   --- Merge the regions
    call MergeRegions(kregions,qitfad)

    return
  end subroutine ITFA_static


!-----------------------------------------------------------------------
  subroutine itfasum(iregion,kmin,kmax,ncat,kpickitfa,wts,cat,qitfa)
! Given a set of weights wts, computes the itfa combination stored in cat

    implicit none

    integer,intent(in) :: iregion
    integer,intent(in) :: kmin,kmax
    integer,intent(in) :: ncat
    integer,intent(in) :: kpickitfa(ncat)
    real,intent(in)    :: wts(ncat)
    real,intent(in)    :: cat(IM,jsta:jend,LM,ncat)
    real,intent(inout) :: qitfa(IM,jsta:jend,LM)

    integer :: i,j,k,idx
    real :: weight
    real :: qitfalast,qijk,qs,wqs
    ! nitfa is the output number of indices used.
    integer :: nitfa

    write(*,*) 'enter itfasum'

    qitfa = SPVAL

!   --- Loop over all 'picked' indices in the sum
    nitfa = 0
    loop_n_idx: do idx=1,ncat
       if(kpickitfa(idx) <= 0) cycle

       weight=wts(idx)

       ! --- Compute the weighted sum and store in qitfa for the current region
       do k=kmin,kmax
       do j=jsta,jend
       do i=1,IM
          if(nitfa == 0) then
             qitfalast = 0.
          else
             qitfalast = qitfa(i,j,k)
             if(ABS(qitfalast-SPVAL) < SMALL1) qitfalast = 0.
          endif

!         remap the raw index value to edr
          qijk = cat(i,j,k,idx)
          call remapq(iregion,idx,qijk,qs)
          if(ABS(qs-SPVAL)<SMALL1) cycle
          wqs=weight*MAX(qs,0.)
          qitfa(i,j,k)=qitfa(i,j,k)+wqs
       enddo
       enddo
       enddo
       nitfa = nitfa+1
    end do loop_n_idx

!   --- Clamp the resultant sum between clampL and clampH
    do k=LM,1,-1
    do j=jsta,jend
    do i=1,IM
       qijk=qitfa(i,j,k)
       if(ABS(qijk-SPVAL)<SMALL1) cycle
       qijk=MAX(qijk,clampL)
       qijk=MIN(qijk,clampH)
       qitfa(i,j,k)=qijk
    enddo
    enddo
    enddo

    return
  end subroutine itfasum

!-----------------------------------------------------------------------
  subroutine MergeRegions(nregions,iregion,kmax,qitfa)
! Merge the boundaries of the regions of multiple regions

    implicit none
     
    integer,intent(in) :: iregion,nregions
    integer,intent(in) :: kmax
    real,intent(inout) :: qitfa(IM,jsta:jend,LM) 

    integer :: i,j,k
    real :: qijk,qsum,qk(LM)
    integer :: ksum
    integer :: kbdy
    write(*,*) 'enter MergeRegions'

    if(nregions <= 1 .or. iregion == 1) return

    kbdy = kmax ! GFS is top-bottom
    kbdy = MAX(kbdy,3)
    kbdy = MIN(kbdy,LM-2)

    do j=jsta,jend
    do i=1,IM

       do k=1,LM
          qk(k)=qitfa(i,j,k) ! save qitfa to qk because qitfa will be changed later
       enddo

       qsum=0.
       ksum=0
       do k=kbdy-1,kbdy+1
          qijk=qk(k)
          if(ABS(qijk-SPVAL)>SMALL1) then
             qsum=qsum+qijk
             ksum=ksum+1
          endif
       enddo
       if(ksum>0) then
          qitfa(i,j,kbdy)=qsum/FLOAT(ksum)
       end if

!      Merge in points below kbdy
       qsum=0.
       ksum=0
       do k=kbdy,kbdy+2
          qijk=qk(k)
          if(ABS(qijk-SPVAL)>SMALL1) then
             qsum=qsum+qijk
             ksum=ksum+1
          endif
       enddo
       if(ksum>0) then
          qitfa(i,j,kbdy+1)=qsum/FLOAT(ksum)
       endif

!      Merge in points above kbdy
       qsum=0.
       ksum=0
       do k=kbdy-2,kbdy
          qijk=qk(k)
          if(ABS(qijk-SPVAL)>SMALL1) then
             qsum=qsum+qijk
             ksum=ksum+1
          endif
       enddo
       if(ksum>0) then
          qitfa(i,j,kbdy-1)=qsum/FLOAT(ksum)
       endif
    enddo
    enddo

    return
  end subroutine MergeRegions


!-----------------------------------------------------------------------
  subroutine itfamaxQ(kmin,kmax,kregions,qitfa,qitfam,qitfax)
! Reads in ITFA (CAT) and ITFA (MWT) and takes max of two, 
! grid point by grid point, and outputs as qitfax.
! qit is a work array.
! If ioutputflag > 0 qitfa is also stored on disk in directory Qdir as idQ.Q. 
! On input ITFADEF or ITFADYN is in qitfa and ITFAMWT is in qitfam.
! The max is only computed between i=imin,imax, j=jmin,jmax,
! k=kmin,kmax and where mask(i,j)>0.

    implicit none

    integer,intent(in) :: kmin,kmax
    integer,intent(in) :: kregions(MAXREGIONS,2)
    real,intent(in) :: qitfa(IM,jsta:jend,LM) 
    real,intent(in) :: qitfam(IM,jsta:jend,LM) 
    real,intent(inout) :: qitfax(IM,jsta:jend,LM) 

    integer :: i,j,k
    real :: qi,qm,qijk

    integer :: Filttype,nftxy,nftz

    qitfax = SPVAL

!   --- Now obtain ITFAMAX=MAX(ITFA,ITFAMWT)
    do k=1,LM
       do j=JSTA,JEND
       do i=1,IM
          qitfax(i,j,k)=SPVAL
          qi=SPVAL
          qm=SPVAL
          if(ABS(qitfa(i,j,k)-SPVAL)<SMALL1) cycle
          qi=qitfa(i,j,k)  ! CAT index
          if(ABS(qitfam(i,j,k)-SPVAL)<SMALL1) cycle
          qm=qitfam(i,j,k)  ! MWT index
          qitfax(i,j,k)=MAX(qi,qm)
       enddo
       enddo
    enddo

!   --- Clamp the resultant sum between clampL and clampH
    do k=1,LM
       do j=JSTA,JEND
       do i=1,IM
          qijk=qitfax(i,j,k)
          if(ABS(qijk-SPVAL)<SMALL1) cycle
          qijk=MAX(qijk,clampitfaL)
          qijk=MIN(qijk,clampitfaH)
          qitfax(i,j,k)=qijk
       enddo
       enddo
    enddo

!   --- Merge the regions
    call MergeRegions(kregions,qitfax)
!
!   --- Perform one smoothing for the blend, and for consistency
!   --- also smooth ITFA
    Filttype=1
    nftxy=1
    nftz=1
    call filt3d(kmin,kmax,nftxy,nftz,Filttype,qitfax)

    return
  end subroutine itfamaxQ


!-----------------------------------------------------------------------
  subroutine remapq(iregion,idx,q,qs)
! Performs mapping of input variable q into output 
! variable qs (0-1).  Input index specific thresholds
! (null,light,moderate,severe) are in timap.  Corresponding
! (0-1) thresholds are in tis.
! If remap_option=1 use piecewise linear mapping
! If remap_option=2 use fit to log-normal PDF

    implicit none

    integer,intent(in) :: iregion,idx
    real,intent(in) :: q
    real,intent(out) :: qs

    real, parameter :: clampidxL=0.0,clampidxH=1.5

    real :: tii(NTI), A(2)

    qs=SPVAL
    if(ABS(q-SPVAL)<=SMALL1) return

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
! Performs linear remap of index value contained in q into
! scaled value in the range (0-1)

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
! Fit q to log-normal PDF, output is epsilon^1/3=qs
! Log(epsilon^(1/3)) = a + b Log(I), where a=A(1), b=A(2)

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

end module gtg_itfa

      subroutine itfasum(qix,qitfa,iditfa,wts,kpickitfa,
     1  nregions,remap_option,iregion,
     2  nx,ny,nz,idmax,nitfa,imin,imax,jmin,jmax,kmin,kmax,
     3  mask,printflag,ic,jc,iprt,Fdir,Qdir,iFQflag)
! Given a set of weights wts, computes the itfa combination
! and return in qitfa(nx,ny,nz).
! The individual indices are read in from idQ.F if iFQflag=1
! or idQ.Q if iFQflag=2
! qix and qs are work arrays.
! The sum is only computed between i=imin,imax, j=jmin,jmax,
! k=kmin,kmax.
! nitfa is the output number of indices used.
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nz,idmax
      real    qix(nx,ny,nz),qitfa(nx,ny,nz)
      real    wts(idmax)
      integer kpickitfa(idmax)
      integer iregion,nregions
      integer remap_option
      integer iFQflag
      integer iditfa
      integer nitfa
      integer imin,imax,jmin,jmax,kmin,kmax
      integer ic,jc,printflag,iprt
      character*200 Fdir,Qdir
      integer mask(nx,ny)
!-----------------------------------------------------------------------
      integer i,j,k,idx,idQ,ierr,ni
      real    qijk,qs,qmin,qmax,indxwt,wqs,qitfalast
      integer nzmax
      parameter (nzmax=201)
      real    zk(nzmax)
      include 'consts.inc'   ! RMISSD
      real    clampitfaL,clampitfaH
      parameter (clampitfaL=0.0,clampitfaH=1.0)
!-----------------------------------------------------------------------
!
! Initializations
      if(printflag.ge.2) then
        idQ=iditfa+399
        write(iprt,*) 'enter itfasum: iregion,iditfa,remap_option=',
     1    iregion,idQ,remap_option
        write(iprt,*) 'iFQflag=',iFQflag
        write(iprt,*) 'imin,imax,jmin,jmax,kmin,kmax=',
     1   imin,imax,jmin,jmax,kmin,kmax
        write(iprt,*) 'nx,ny,nz=',nx,ny,nz
        write(iprt,*) 'nregions=',nregions
        write(iprt,*) 'clampitfaL,clampitfaH=',clampitfaL,clampitfaH
        if(printflag.ge.3) then
          ni=0
          do idx=1,idmax
            indxwt=wts(idx)
            if(kpickitfa(idx).gt.0.) then
              ni=ni+1
              write(iprt,*) 'iregion,idx,ipickitfa,wts=',iregion,idx,
     1        kpickitfa(idx),indxwt
            endif
          enddo
          write(iprt,*) 'iregion,nitfa=',iregion,ni
        endif
      endif
!
! Initializations
      do k=1,nz
        zk(k)=RMISSD
      enddo
      nitfa=0
!
! Loop over all indices in the sum
      do idx=1,idmax
        if(kpickitfa(idx).le.0) go to 33
        idQ=idx+399
!   Retrieve the index 
        if(iFQflag.eq.1) then ! .F
          call readF(qix,nx,ny,nz,idQ,printflag,iprt,FDir,ierr)
        else
          call getqix(qix,nx,ny,nz,idQ,printflag,iprt,QDir,ierr)
        endif
        if(ierr.ne.0) then
          write(iprt,*) 'error in itfasum reading index idQ=',idQ
          go to 33
        endif
        if(printflag.ge.2) then
          write(iprt,*) 'retrieved qix for idQ=',idQ
          if(printflag.ge.3) then
            do k=1,nz
              write(iprt,*) 'i,j,k,qix=',ic,jc,k,qix(ic,jc,k)
            enddo
          endif
        endif
!   Compute the weighted sum and store in qitfa for the
!   current region
        indxwt=wts(idx)
        if(printflag.ge.2) then
          write(iprt,*) 'iregion,idQ,wts=',iregion,idQ,indxwt
        endif
        do i=imin,imax
        do j=jmin,jmax
          if(mask(i,j).le.0) go to 32
          do k=kmin,kmax
            if(nitfa.le.0) then
              qitfalast=0.
            else
              qitfalast=qitfa(i,j,k)
              if(ABS(qitfalast-RMISSD).LT.SMALL1) qitfalast=0.
            endif
!       remap the raw index qix value to edr
            qijk=qix(i,j,k)
            call remapq(qijk,qs,iregion,idx,printflag,iprt)
            wqs=RMISSD
            if(ABS(qitfalast-RMISSD).LT.SMALL1) go to 30
            if(ABS(qs-RMISSD).LT.SMALL1) go to 30
            wqs=indxwt*MAX(qs,0.)
            qitfa(i,j,k)=qitfalast+wqs
   30       continue
            if(printflag.ge.2 .and. i.eq.ic .and. j.eq.jc) then
              write(iprt,*) 'idQ,i,j,k,qiftalast,qs,wqs,qitfa=',
     1         idQ,ic,jc,k,qitfalast,qs,wqs,qitfa(ic,jc,k)
            endif
          enddo
   32     continue   
        enddo
        enddo
!   Increment the number of indices used
        nitfa=nitfa+1
   33   continue
      enddo  ! idx loop
!
! Clamp the resultant sum between clampL and clampH
      qmin=+1.0E30
      qmax=-1.0E30
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
        qijk=qitfa(i,j,k)
        if(ABS(qijk-RMISSD).LT.SMALL1) go to 35
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
! Merge the regions
      call MergeRegions(qitfa,nx,ny,nz,imin,imax,jmin,jmax,kmin,
     1  iregion,nregions,mask,printflag,ic,jc,iprt)
!
      if(printflag.ge.1) then 
        if(printflag.ge.2) then
          do k=1,nz
            write(iprt,*)  'i,j,k,qitfa=',ic,jc,k,qitfa(ic,jc,k)
          enddo
        endif
        write(iprt,*) 'exit  itfasum: nitfa=',nitfa
        write(iprt,*) 'exit  itfasum: clamped qitfa min,max=',qmin,qmax
      endif
!
      return
      end
!
      subroutine itfamax(qitfa,qitfam,qitfax,qit,iditfa,iditfam,
     1  iditfax,minrgn,maxrgn,maxregions,imin,imax,
     2  jmin,jmax,kminrgn,kmaxrgn,mbc,mask,nx,ny,nz,
     3  printflag,ic,jc,iprt,
     4  ioutputflag,cnamei,Fdir,Qdir,iFQflag)
! Reads in ITFA (CAT) and ITFA (MWT) and takes max of two, 
! grid point by grid point, and outputs as qitfax.
! qit is a work array.
! If ioutputflag > 0 qitfax is also stored on disk in directory
! Fdir as idQ.F (iFQflag=1) or Qdir as idQ.Q (iFQflag=2), where
! idQ=iditfax+399. 
! On input ITFADEF or ITFADYN is in qitfa and ITFAMWT is in qitfam.
! The max is only computed between i=imin,imax, j=jmin,jmax,
! k=kminrgn(iregion),kmaxrgn(iregion), iregion=minrgn,maxrgn,
! and where mask(i,j)>0.
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nz,maxregions
      real    qitfa(nx,ny,nz),qitfam(nx,ny,nz),qitfax(nx,ny,nz)
      real    qit(nx,ny,nz)
      integer kminrgn(maxregions),kmaxrgn(maxregions)
      integer iditfa,iditfam,iditfax
      integer imin,imax,jmin,jmax,mbc
      integer minrgn,maxrgn
      integer ic,jc,printflag,iprt
      character*24 cnamei
      real    qi,qm
      integer ioutputflag
      integer iFQflag
      integer mask(nx,ny)
      character*200 Fdir,Qdir
!-----------------------------------------------------------------------
      integer i,j,k,idQ,jdQ,kdQ
      real    qijk,qmin,qmax
      integer kminr,kmaxr
      integer iregion,nregions
      integer Filttype,nftxy,nftz
      integer ierr
      include 'consts.inc'   ! RMISSD
      real    clampitfaL,clampitfaH
      parameter (clampitfaL=0.0,clampitfaH=1.0)
!-----------------------------------------------------------------------
!
! Initializations
      if(printflag.ge.1) then
        write(iprt,*) 'enter itfamax: iditfa,iditfam,iditfax=',
     1    iditfa,iditfam,iditfax
        write(iprt,*) 'nx,ny,nz=',nx,ny,nz
        if(printflag.ge.2) then
          do k=1,nz
            write(iprt,*) 'i,j,k,qitfa,qitfam=',ic,jc,k,qitfa(ic,jc,k),
     1        qitfam(ic,jc,k)
          enddo
        endif
      endif
      idQ=iditfa+399  ! input ITFA
      jdQ=iditfam+399 ! input ITFAMWT
      kdQ=iditfax+399 ! output ITFAMAX
      call TIinit(qitfax,nx,ny,nz)
      nregions=maxrgn-minrgn+1
!
! Now obtain ITFAMAX=MAX(ITFA,ITFAMWT) region by region
      qmin=+1.0E30
      qmax=-1.0E30
      do iregion=minrgn,maxrgn
        kminr=kminrgn(iregion)
        kmaxr=kmaxrgn(iregion)
        qmin=+1.0E30
        qmax=-1.0E30
        do i=imin,imax
        do j=jmin,jmax
        do k=kminr,kmaxr
          if(mask(i,j).le.0) go to 44
          qi=RMISSD
          qm=RMISSD
          if(ABS(qitfa(i,j,k)-RMISSD).lt.SMALL1) go to 44
          qi=qitfa(i,j,k)  ! CAT index
          qitfa(i,j,k)=qi
          if(ABS(qitfam(i,j,k)-RMISSD).lt.SMALL1) go to 44
          qm=qitfam(i,j,k)  ! MWT index
          qijk=MAX(qi,qm)
          qijk=MAX(qijk,clampitfaL)
          qijk=MIN(qijk,clampitfaH)
          qitfax(i,j,k)=qijk
          qmin=MIN(qmin,qijk)
          qmax=MAX(qmax,qijk)
   44     continue
          if(printflag.ge.2 .and. i.eq.ic .and. j.eq.jc) then
            write(iprt,*) 'i,j,k,ITFA,ITFAMWT,MAX=',i,j,k,
     1        qi,qm,qitfax(i,j,k)
          endif
        enddo
        enddo
        enddo
!
!   Merge the regions
        call MergeRegions(qitfax,nx,ny,nz,imin,imax,jmin,jmax,kminr,
     1    iregion,nregions,mask,printflag,ic,jc,iprt)
!
      enddo  ! region loop
!
! Perform one smoothing for the blend, and for consistency
! also smooth ITFA
      Filttype=1
      nftxy=1
      nftz=1
      call filt3d(qitfax,qit,nx,ny,nz,imin,imax,jmin,jmax,
     1  kminr,kmaxr,mbc,nftxy,nftz,Filttype)
!     call filt3d(qitfa ,qit,nx,ny,nz,imin,imax,jmin,jmax,
!    1  kminr,kmaxr,mbc,nftxy,nftz,Filttype)
      call filt3d(qitfam,qit,nx,ny,nz,imin,imax,jmin,jmax,
     1  kminr,kmaxr,mbc,nftxy,nftz,Filttype)
      if(printflag.ge.1) then
        do k=1,nz
          write(iprt,*)'i,j,k,itfamax smooth=',ic,jc,k,qitfax(ic,jc,k)
        enddo
        do k=1,nz
          write(iprt,*)'i,j,k,itfamwt smooth=',ic,jc,k,qitfam(ic,jc,k)
        enddo
      endif
!
      if(ioutputflag.gt.0) then
!   Write out ITFAMAX to output binary file
        idQ=iditfax+399
        if(iFQflag.eq.1) then ! .F
          call writeF(qitfax,nx,ny,nz,idQ,cnamei,Fdir,printflag,
     1      iprt,ierr)
        else ! .Q
          call putqix(qitfax,nx,ny,nz,idQ,cnamei,printflag,iprt,
     1      Qdir,ierr)
        endif
      endif
!
      if(printflag.ge.1) then
        if(printflag.ge.2) then
          do k=1,nz
            write(iprt,*) 'i,j,k,ITFAMAX=',ic,jc,k,qitfax(ic,jc,k)
          enddo
        endif
        idQ=iditfax+399
        call TIstats(qitfax,nx,ny,nz,1,nx,1,ny,kminr,kmaxr,qmax,qmin,
     1    idQ,iprt)
        write(iprt,*) 'exit  itfamax: qmin,qmax=',qmin,qmax
      endif
!
      return
      end
!
