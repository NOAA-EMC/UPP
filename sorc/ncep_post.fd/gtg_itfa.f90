module gtg_itfa
  use ctlblk_mod, only: jsta,jend, IM,JM,LM, SPVAL
  use gtg_config
  use gtg_filter

  implicit none
contains
!-----------------------------------------------------------------------
  subroutine ITFA_MWT(nids,indxpicked,kregions,cat,qitfam)
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
    integer :: nids_itfa,i

    write(*,*) 'enter ITFA_MWT'

    qitfam = 0.

    nids_itfa = 0
    kpickitfa = 0
    do i = 1, nids
       if(indxpicked(i) >= 476 .and. indxpicked(i) <= 490) then ! MWT
          nids_itfa = nids_itfa + 1
          kpickitfa(nids_itfa) = indxpicked(i)
       end if
    end do
    if (nids_itfa <= 0) then
       write(*,*) "There is no MWT indices picked"
       return
    end if

!   --- Perform the combination using altitude-dependent static weights
    call itfacompQ(nids_itfa,kpickitfa(1:nids_itfa),kregions,static_wgt,&
         cat,clampitfaL,clampitfaH,qitfam)

!   --- Merge the regions
    call MergeItfaRegions(kregions,qitfam)

    return
  end subroutine ITFA_MWT


!-----------------------------------------------------------------------
  subroutine ITFA_static(nids,indxpicked,kregions,cat,qitfad)
!     --- Computes an ITFA combination using static weights and outputs as qitfa.
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
    real,intent(inout) :: qitfad(IM,jsta:jend,LM)

    integer :: kpickitfa(nids)
    integer :: nids_itfa,i

    write(*,*) 'enter ITFA_static'

    qitfad = 0.

    nids_itfa = 0
    kpickitfa = 0
    do i = 1, nids
       if(indxpicked(i) <= 475) then ! CAT
          nids_itfa = nids_itfa + 1
          kpickitfa(nids_itfa) = indxpicked(i)
       end if
    end do
    if (nids_itfa <= 0) then
       write(*,*) "There is no ITFA_static indices picked"
       return
    end if

!   --- Perform the combination using altitude-dependent static weights
    call itfacompQ(nids_itfa,kpickitfa(1:nids_itfa),kregions,static_wgt,&
         cat,clampitfaL,clampitfaH,qitfad)

!   --- Merge the regions
    call MergeItfaRegions(kregions,qitfad)

    return
  end subroutine ITFA_static


!-----------------------------------------------------------------------
  subroutine itfacompQ(nids,indxpicked,kregions,wts,cat,clampL,clampH,qitfa)
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
    real,intent(in) :: clampL,clampH
    real,intent(inout) :: qitfa(IM,jsta:jend,LM)

    real :: qs(IM,jsta:jend,LM) ! work array

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

    integer :: i,j,k,kk,iregion,nregions
    real :: qijk,qsum,qk(LM)
    integer :: ksum
    integer :: kbdy,kbdy_m,kbdy_p
    write(*,*) 'enter MergeItfaRegions'


    do j=jsta,jend
    do i=1,IM

       nregions = 0
       do kk = 1, MAXREGIONS-1
          if(kregions(i,j,kk,1) > 0) nregions = nregions + 1
       end do
       if(nregions <= 1) cycle

       do kk = 1, MAXREGIONS-1

          kbdy=kregions(i,j,kk,2)
          if(kbdy <= 0) cycle

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
                qijk=qk(k)
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
                qijk=qk(k)
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


!-----------------------------------------------------------------------
  subroutine itfamaxQ(kmin,kmax,kregions,qitfa,qitfam,qitfax)
!     --- Reads in ITFA (CAT) and ITFA (MWT) and takes max of two, 
!     --- grid point by grid point, and outputs as qitfax.
!     --- qit is a work array.
!     --- If ioutputflag > 0 qitfa is also stored on disk in directory Qdir as idQ.Q. 
!     --- On input ITFADEF or ITFADYN is in qitfa and ITFAMWT is in qitfam.
!     --- The max is only computed between i=imin,imax, j=jmin,jmax,
!     --- k=kmin,kmax and where mask(i,j)>0.

    implicit none

    integer,intent(in) :: kmin,kmax
    integer,intent(in) :: kregions(IM,jsta:jend,MAXREGIONS,2)
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
          if(ABS(qitfa(i,j,k)-SPVAL)<1.0E-3) cycle
          qi=qitfa(i,j,k)  ! CAT index
          if(ABS(qitfam(i,j,k)-SPVAL)<1.0E-3) cycle
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
          if(ABS(qijk-SPVAL)<1.0E-3) cycle
          qijk=MAX(qijk,clampitfaL)
          qijk=MIN(qijk,clampitfaH)
          qitfax(i,j,k)=qijk
       enddo
       enddo
    enddo

!   --- Merge the regions
    call MergeItfaRegions(kregions,qitfax)
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

end module gtg_itfa
