module gtg_itfa
  use ctlblk_mod, only: jsta,jend, IM,JM,LM, SPVAL
  use gtg_config
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
          call itfa_comp()
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

          call itfa_comp()
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
    call MergeItfaRegions(kregions,qitfam)

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
    call MergeItfaRegions(kregions,qitfam)

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
    call MergeItfaRegions(kregions,qitfad)

    return
  end subroutine ITFA_static


!-----------------------------------------------------------------------
  subroutine itfa_comp(iregion,kmin,kmax,ncat,kpickitfa,wts,cat,qitfa)
! Given a set of weights wts, computes the itfa combination stored in cat

    implicit none

    integer,intent(in) :: iregion
    integer,intent(in) :: kmin,kmax
    integer,intent(in) :: ncat
    integer,intent(in) :: kpickitfa(ncat)
    real,intent(in)    :: wts(ncat)
    real,intent(in)    :: cat(IM,jsta:jend,LM,ncat)
    real,intent(inout) :: qitfa(IM,jsta:jend,LM)

    real,parameter :: clampL=0.0,clampH=1.0

    integer :: i,j,k,idx
    real :: weight
    real :: qitfalast,qijk,qs,wqs
    ! nitfa is the output number of indices used.
    integer :: nitfa

    write(*,*) 'enter itfa_comp'

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
          call remapq(qijk,qs,iregion,idx)
       write(*,*) "itfa_comp cat nn, idx=",nn,idx
       call remapi(idx,kregions,cat(1:IM,jsta:jend,1:LM,nn),qs)
          if(ABS(qs(i,j,k)-SPVAL)<1.0E-3) cycle
          wqs=weight*MAX(qs(i,j,k),0.)
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
        if(ABS(qijk-SPVAL)<1.0E-3) cycle
        qijk=MAX(qijk,clampL)
        qijk=MIN(qijk,clampH)
        qitfa(i,j,k)=qijk
     enddo
     enddo
     enddo

     return
   end subroutine itfa_comp

!-----------------------------------------------------------------------
  subroutine MergeItfaRegions(kregions,qitfa)
! Merge the boundaries of the regions of multiple regions

    implicit none
     
    integer,intent(in) :: kregions(MAXREGIONS,2)
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
          if(kregions(kk,1) > 0) nregions = nregions + 1
       end do
       if(nregions <= 1) cycle

       do kk = 1, MAXREGIONS-1

          kbdy=kregions(kk,2)
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
!        Merge in points above kbdy
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
!        Merge in points below kbdy
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
! Performs mapping of input variable q into output 
! variable qs (0-1).  Input index specific thresholds
! (null,light,moderate,severe) are in timap.  Corresponding
! (0-1) thresholds are in tis. remap_option is in scorei.inc
! If remap_option=1 use piecewise linear mapping
! If remap_option=2 use fit to log-normal PDF

    implicit none

    integer,intent(in) :: idx
    integer,intent(in) :: kregions(MAXREGIONS,2)
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

!  Determine which k region
       iregion = -1
       do kk = 1, MAXREGIONS
          if (k <= kregions(kk,1) .and. k > kregions(kk,2)) then
             iregion = kk
             exit
          endif
       end do
       if (iregion < 0) cycle

       qijk=q(i,j,k)
       if(ABS(qijk-SPVAL)<=1.0E-3) cycle
       qs(i,j,k)=0.
       if(idQ==465) then ! SGSTKE
!     Map according to epsilon=(C/L)e^3/2 with L~dx
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

      subroutine ITFA_MWT(qix,qit,qitfam,iditfam,nx,ny,nz,
     1 minrgn,maxrgn,maxregions,ipickitfa,idmax,
     2 static_wts,iFcstType,remap_option,imin,imax,jmin,jmax,
     3 kminrgn,kmaxrgn,mask,printflag,ic,jc,iprt,
     4 ioutputflag,cnamei,Fdir,Qdir,iFQflag)
! Computes an ITFA MWT combination based on a set of dynamic weights, and
! outputs as qitfam.
! qix and qit are work arrays.
! The sum is only computed between i=imin,imax, j=jmin,jmax,
! k=kminrgn(iregion),kmaxrgn(iregion).
! The individual indices are read in from idQ.F if iFQflag=1
! or idQ.Q if iFQflag=2
! Output values are clamped between clampitfaL,clampitfaH
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nz,idmax,maxregions
      real    qix(nx,ny,nz),qit(nx,ny,nz),qitfam(nx,ny,nz)
      integer kminrgn(maxregions),kmaxrgn(maxregions)
      integer minrgn,maxrgn
      integer ipickitfa(maxregions,idmax)
      real    static_wts(maxregions,idmax)
      integer iFcstType(idmax)
      integer iditfam
      integer remap_option
      integer nwts
      integer imin,imax,jmin,jmax
      integer ioutputflag,iFQflag
      integer ic,jc,printflag,iprt
      integer mask(nx,ny)
      character*24 cnamei
      character*200 Fdir,Qdir
!-----------------------------------------------------------------------
      integer i,j,k,idx,idQ
      integer iregion,ni,nregions
      real    wsum
      real    qmin,qmax
      integer kminr,kmaxr
      integer ierr
      integer nimax
      parameter (nimax=101)
      integer kpickitfa(nimax)
      real    wts(nimax)
      include 'consts.inc'   ! RMISSD
      integer BLT,CAT,CIT,MWT,WAKE
      parameter (BLT=1, CAT=2, CIT=3, MWT=4, WAKE=5)
!-----------------------------------------------------------------------
!
! Initializations
      idQ=iditfam+399
      if(printflag.ge.1) then
        write(iprt,*) 'enter ITFA_MWT: idQ,iFQflag,remap_option=',
     1    idQ,iFQflag,remap_option
        write(iprt,*) 'nx,ny,nz=',nz,ny,nz
        write(iprt,*) 'minregion,maxregion=',minrgn,maxrgn
        do iregion=minrgn,maxrgn
          do idx=1,idmax
            if(iFcstType(idx).eq.MWT .and. ipickitfa(iregion,idx).gt.0)
     1      then
              write(iprt,*)'iregion,idQ,ipickitfa,wts=',iregion,idx+399,
     1          ipickitfa(iregion,idx),static_wts(iregion,idx)
            endif
          enddo
        enddo
      endif
      nregions=maxrgn-minrgn+1
      call TIinit(qix,nx,ny,nz)
      call TIinit(qit,nx,ny,nz)
! Initialize ITFAMWT to 0.
      do k=1,nz
      do j=1,ny
      do i=1,nx
        qitfam(i,j,k)=0.
      enddo
      enddo
      enddo
!
! Perform the combination using altitude-dependent static weights
      do iregion=minrgn,maxrgn
        kminr=kminrgn(iregion)
        kmaxr=kmaxrgn(iregion)
!   Get the static MWT weights
        if(printflag.ge.2) then 
          write(iprt,*) 'computing ITFAMWT for iregion=',iregion
          write(iprt,*) 'kminr,kmaxr=',kminr,kmaxr
        endif
        nwts=0
        wsum=0.
        do idx=1,idmax
          kpickitfa(idx)=0
          if(iFcstType(idx).eq.MWT .and. ipickitfa(iregion,idx).gt.0)
     1    then
            kpickitfa(idx)=ipickitfa(iregion,idx)
            wts(idx)=static_wts(iregion,idx)
            nwts=nwts+1
            wsum=wsum+wts(idx)
            if(printflag.ge.2) then
              write(iprt,*) 'iregion,idQ,kpickitfa=',iregion,idx+399,
     1         kpickitfa(idx)
            endif
          endif
        enddo
!   Check the weights
        if(printflag.ge.2) write(iprt,*) 'nwts,wsum=',nwts,wsum
        if(nwts.le.0) then
          write(iprt,*) 'No weights specificied for MWT region=',iregion
          go to 153
        endif
        nwts=0
        if(printflag.ge.2) then
          ni=0
          do idx=1,idmax
            if(kpickitfa(idx).gt.0) then
              ni=ni+1
              write(iprt,*) 'iregion,idQ,kpickitfa,wt=',iregion,idx+399,
     1         kpickitfa(idx),wts(idx)
            endif
          enddo
          write(iprt,*) 'nindices to use=',ni
        endif
        nwts=0
!   Compute the ITFAMWT combination for this region
!   qix and qit are work arrays.  Output is in qitfam
        call itfasum(qix,qitfam,iditfam,wts,kpickitfa,
     1    nregions,remap_option,iregion,
     2    nx,ny,nz,idmax,nwts,imin,imax,jmin,jmax,kminr,kmaxr,
     3    mask,printflag,ic,jc,iprt,Fdir,Qdir,iFQflag)
        if(printflag.ge.2) then
          write(iprt,*) 'after itfasum: iregion,nitfa=',iregion,nwts
          do k=1,nz
            write(iprt,*) 'i,j,k,itfamwt=',ic,jc,k,qitfam(ic,jc,k)
          enddo
        endif
  153   continue
      enddo  ! iregion loop
!
      if(ioutputflag.gt.0) then
!   Write out ITFAMWT to output binary file
        idQ=iditfam+399
        if(iFQflag.eq.1) then ! .F
          call writeF(qitfam,nx,ny,nz,idQ,cnamei,Fdir,printflag,
     1      iprt,ierr)
        else ! .Q
          call putqix(qitfam,nx,ny,nz,idQ,cnamei,printflag,iprt,
     1      Qdir,ierr)
        endif
      endif
!
      if(printflag.ge.1) then
        if(printflag.ge.2) then
          do k=1,nz
            write(iprt,*) 'i,j,k,itfamwt=',ic,jc,k,qitfam(ic,jc,k)
          enddo
        endif
        call TIstats(qitfam,nx,ny,nz,1,nx,1,ny,kminr,kmaxr,qmax,qmin,
     1    idQ,iprt)
        write(iprt,*)'exit  ITFA_MWT: idQ,nitfa,qmin,qmax=',
     1    iditfam+399,nwts,qmin,qmax
      endif
!
      return
      end
!
      subroutine ITFA_static(qix,qit,qitfad,iditfad,nx,ny,nz,
     1  minrgn,maxrgn,maxregions,ipickitfa,idmax,
     2  static_wts,iFcstType,remap_option,imin,imax,jmin,jmax,
     3  kminrgn,kmaxrgn,mask,printflag,ic,jc,iprt,
     4  ioutputflag,cnamei,Fdir,Qdir,iFQflag)
! Computes an ITFA combination using static weights and outputs as qitfad.
! qix and qit are work arrays.
! The sum is only computed between i=imin,imax, j=jmin,jmax,
! k=kminrgn(iregion),kmaxrgn(iregion).
! The individual indices are read in from idQ.F if iFQflag=1
! or idQ.Q if iFQflag=2
! Output values are clamped between clampitfaL,clampitfaH
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nz,idmax,maxregions
      real    qix(nx,ny,nz),qit(nx,ny,nz),qitfad(nx,ny,nz)
      integer kminrgn(maxregions),kmaxrgn(maxregions)
      integer minrgn,maxrgn
      integer ipickitfa(maxregions,idmax)
      real    static_wts(maxregions,idmax)
      integer iFcstType(idmax)
      integer remap_option
      integer iditfad
      integer nwts
      integer imin,imax,jmin,jmax
      integer ioutputflag,iFQflag
      integer ic,jc,printflag,iprt
      integer mask(nx,ny)
      character*24 cnamei
      character*200 Fdir,Qdir
!-----------------------------------------------------------------------
      integer k,idx,idQ
      integer ni
      integer iregion,nregions
      real    wsum
      real    qmin,qmax
      integer kminr,kmaxr
      integer ierr
      integer nimax
      parameter (nimax=101)
      integer kpickitfa(nimax)
      real    wts(nimax)
      include 'consts.inc'   ! RMISSD
      integer BLT,CAT,CIT,MWT,WAKE
      parameter (BLT=1, CAT=2, CIT=3, MWT=4, WAKE=5)
!-----------------------------------------------------------------------
!
! Initializations
      idQ=iditfad+399
      if(printflag.ge.1) then
        write(iprt,*)'enter ITFA_static: idQ,iFQflag,remap_option=',
     1    idQ,iFQflag,remap_option
        write(iprt,*) 'minregion,maxregion=',minrgn,maxrgn
        do iregion=minrgn,maxrgn
          do idx=1,idmax
            if(ipickitfa(iregion,idx).gt.0) then
              write(iprt,*)'iregion,idQ,ipickitfa,wts=',iregion,idx+399,
     1          ipickitfa(iregion,idx),static_wts(iregion,idx)
            endif
          enddo
        enddo
      endif
      nregions=maxrgn-minrgn+1
!
      call TIinit(qitfad,nx,ny,nz)
      call TIinit(qix   ,nx,ny,nz)
      call TIinit(qit   ,nx,ny,nz)
!
! Perform the combination using altitude-dependent static weights
      do iregion=minrgn,maxrgn
        kminr=kminrgn(iregion)
        kmaxr=kmaxrgn(iregion)
        if(printflag.ge.2) then 
          write(iprt,*) 'computing ITFASTATIC for iregion=',iregion
          write(iprt,*) 'kminr,kmaxr=',kminr,kmaxr
        endif
!   Get the static CAT weights
        nwts=0
        wsum=0.
        do idx=1,idmax
          kpickitfa(idx)=0
          if(iFcstType(idx).eq.CAT .and. ipickitfa(iregion,idx).gt.0)
     1    then
            kpickitfa(idx)=ipickitfa(iregion,idx)
            wts(idx)=static_wts(iregion,idx)
            nwts=nwts+1
            wsum=wsum+wts(idx)
            if(printflag.ge.2) then
              write(iprt,*) 'iregion,idQ,kpickitfa=',iregion,idx+399,
     1         kpickitfa(idx)
            endif
          endif
        enddo
!   Check the weights
        if(printflag.ge.2) write(iprt,*) 'nwts,wsum=',nwts,wsum
        if(nwts.le.0) then
          write(iprt,*) 'No weights specificied for CAT region=',iregion
          go to 153
        endif
        nwts=0
        if(printflag.ge.2) then
          ni=0
          do idx=1,idmax
            if(kpickitfa(idx).gt.0) then
              ni=ni+1
              write(iprt,*) 'iregion,idQ,kpickitfa,wt=',iregion,idx+399,
     1         kpickitfa(idx),wts(idx)
            endif
          enddo
          write(iprt,*) 'nindices to use=',ni
        endif
        nwts=0
!   Compute the ITFADEF combination for this region
!   qix and qit are work arrays.  Output is in qitfad
        call itfasum(qix,qitfad,iditfad,wts,kpickitfa,
     1    nregions,remap_option,iregion,
     2    nx,ny,nz,idmax,nwts,imin,imax,jmin,jmax,kminr,kmaxr,
     3    mask,printflag,ic,jc,iprt,Fdir,Qdir,iFQflag)
        if(printflag.ge.2) then
          write(iprt,*) 'after itfasum: iregion,nitfa=',iregion,nwts
          do k=1,nz
            write(iprt,*) 'i,j,k,itfadef=',ic,jc,k,qitfad(ic,jc,k)
          enddo
        endif
  153   continue
      enddo  ! iregion loop
!
      if(ioutputflag.gt.0) then
!   Write out ITFADEF to output binary file
        idQ=iditfad+399
        if(iFQflag.eq.1) then ! .F
          call writeF(qitfad,nx,ny,nz,idQ,cnamei,Fdir,printflag,
     1      iprt,ierr)
        else ! .Q
          call putqix(qitfad,nx,ny,nz,idQ,cnamei,printflag,iprt,
     1      Qdir,ierr)
        endif
      endif
!
      if(printflag.ge.1) then
        if(printflag.ge.2) then
          do k=1,nz
            write(iprt,*) 'i,j,k,itfadef=',ic,jc,k,qitfad(ic,jc,k)
          enddo
        endif
        call TIstats(qitfad,nx,ny,nz,1,nx,1,ny,kminr,kmaxr,qmax,qmin,
     1    idQ,iprt)
        write(iprt,*)'exit  ITFA_static: idQ,nitfa,qmin,qmax=',
     1    iditfad+399,nwts,qmin,qmax
      endif
!
      return
      end
!
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
              if(ABS(qitfalast-RMISSD).LT.1.0E-3) qitfalast=0.
            endif
!       remap the raw index qix value to edr
            qijk=qix(i,j,k)
            call remapq(qijk,qs,iregion,idx,printflag,iprt)
            wqs=RMISSD
            if(ABS(qitfalast-RMISSD).LT.1.0E-3) go to 30
            if(ABS(qs-RMISSD).LT.1.0E-3) go to 30
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
          if(ABS(qitfa(i,j,k)-RMISSD).lt.1.0E-3) go to 44
          qi=qitfa(i,j,k)  ! CAT index
          qitfa(i,j,k)=qi
          if(ABS(qitfam(i,j,k)-RMISSD).lt.1.0E-3) go to 44
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
      subroutine GetkBdys(zm,zi,nx,ny,nz,nzi,imin,imax,jmin,jmax,
     1  kminrgn,kmaxrgn,zregion,minrgn,maxrgn,maxregions,printflag,
     2  iprt,iFQflag)
! Merge the boundaries of the regions of multiple regions
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nz,nzi,maxregions
      real    zm(nx,ny,nz)
      real    zi(nzi)
      real    zregion(maxregions,2)
      integer kminrgn(maxregions),kmaxrgn(maxregions)
      integer minrgn,maxrgn
      integer imin,imax,jmin,jmax
      integer iprt,printflag
      integer iFQflag
      integer i,j,k,iregion
      integer kminr,kmaxr
      real    zmin
      integer nregions
      integer nzk,izfmin,jzfmin
      integer nzmax
      parameter (nzmax=201)
      real    zk(nzmax)
      include 'consts.inc'   ! RMISSD
!-----------------------------------------------------------------------
!
! Initializations
      nregions=maxrgn-minrgn+1
      if(printflag.ge.2) then
        write(iprt,*) 'enter GetkBdys: iFQflag=',iFQflag
        write(iprt,*) 'nx,ny,nz,nzi=',nx,ny,nz,nzi
        write(iprt,*) 'nregions=',nregions
      endif
      do k=1,nz
        zk(k)=RMISSD
      enddo
!
! Restrict vertical indices in the sum to only be those in 
! input region iregion
      if(iFQflag.eq.1 ) then !.F
!   If computing on native grid (iFQflag=1) search left boundary
!   points for the minimum altitude, i.e., lowest terrain.  Use
!   this location to define the vertical limits for the input region.
        zmin=1.0E20
        izfmin=1
        jzfmin=1
        do i=imin,imin+1
        do j=jmin,jmax
          if(zm(i,j,1).lt.zmin) then
            zmin=zm(i,j,1)
            izfmin=i
            jzfmin=j
          endif
        enddo
        enddo
        do k=1,nz
          zk(k)=zm(izfmin,jzfmin,k)*3.28  ! m -> ft
          zk(k)=MAX(zk(k),0.)
        enddo
        nzk=nz
      else
!   This is interpolated grid (.Q, iFQflag=2).  Find the vertical
!   limits for the input region
        do k=1,nzi
          zk(k)=zi(k) ! ft
        enddo
        nzk=nzi
      endif
      if(printflag.ge.2) then
        if(iFQflag.eq.1) then
          write(iprt,*) 'iFQflag,z0=',iFQflag,zm(izfmin,jzfmin,1)*3.28
        endif
        do k=1,nzk
          write(iprt,*) 'k,zk(ft)=',k,zk(k)
        enddo
      endif
!
! Find the vertical indices k corresponding to the minimum
! and maximum extents of each region
      do iregion=1,nregions
      kminr=1
      kmaxr=nzk
      if(iregion.eq.1) go to 12
      do k=1,nzk
        kminr=k
        if(zk(k).ge.zregion(iregion,1)) go to 12
      enddo
   12 continue
      do k=kminr+1,nzk
        kmaxr=k
        if(zk(k).ge.zregion(iregion,2)) go to 14
      enddo
   14 continue
      if(printflag.ge.2) then
        write(iprt,*) 'zregion12,kminr,kmaxr=',zregion(iregion,1),
     1    zregion(iregion,2),kminr,kmaxr
      endif
!
      kminrgn(iregion)=kminr
      kmaxrgn(iregion)=kmaxr
        
      enddo ! iregion loop
!
      if(printflag.ge.2) then
        do iregion=1,nregions
          write(iprt,*) 'iregion,kmin,kmax=',iregion,kminrgn(iregion),
     1    kmaxrgn(iregion)
        enddo
        write(iprt,*) 'exit  GetkBdys'
      endif
!
      return
      end
!
      subroutine MergeRegions(qitfa,nx,ny,nz,imin,imax,jmin,jmax,kminr,
     1  iregion,nregions,mask,printflag,ic,jc,iprt)
! Merge the boundaries of the regions of multiple regions
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nz
      real    qitfa(nx,ny,nz)
      integer mask(nx,ny)
      integer imin,imax,jmin,jmax,kminr
      integer printflag,iprt
      integer ic,jc
      integer i,j,k,iregion,nregions
      real    qijk,qsum
      integer ksum,kbdy
      integer nzmax
      parameter (nzmax=201)
      real    qk(nzmax)
      include 'consts.inc'   ! RMISSD
!-----------------------------------------------------------------------
!
      if(printflag.ge.2) then
        write(iprt,*) 'enter MergeRegions: iregion,nregions=',
     1    iregion,nregions
        if(printflag.ge.2) then
          do k=1,nz
            write(iprt,*) 'i,j,k,inputqitfa=',ic,jc,k,qitfa(ic,jc,k)
          enddo
        endif
      endif
!
      if(nregions.le.1 .or. iregion.eq.1) return
      kbdy=kminr
      kbdy=MAX(kbdy,3)
      kbdy=MIN(kbdy,nz-2)
      if(printflag.ge.2) then
        write(iprt,*) 'iregion,kbdy=',iregion,kbdy
      endif
      do j=jmin,jmax
      do i=imin,imax
        if(mask(i,j).le.0) go to 40
        do k=1,nz
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
          if(printflag.ge.3 .and. i.eq.ic .and. j.eq.jc) then
            write(iprt,*) 'iregion,i,j,k,qitfa,qsum,ksum=',
     1        iregion,ic,jc,k,qijk,qsum,ksum
          endif
        enddo
        if(ksum.gt.0) then
          qitfa(i,j,kbdy)=qsum/FLOAT(ksum)
        endif
!   Merge in points above kbdy
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
!   Merge in points below kbdy
        qsum=0.
        ksum=0
        do k=kbdy-2,kbdy
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
          do k=1,nz
            write(iprt,*) 'iregion,i,j,k,qk,qitfa=',
     1        iregion,ic,jc,k,qk(k),qitfa(i,j,k)
          enddo
        endif
   40   continue
      enddo
      enddo
!
      if(printflag.ge.2) then
        if(printflag.ge.2) then
          do k=1,nz
            write(iprt,*) 'i,j,k,mergedqitfa=',ic,jc,k,qitfa(ic,jc,k)
          enddo
        endif
        write(iprt,*) 'exit  MergeRegions'
      endif
!
      return
      end
!
      subroutine get_dynamic_wts(qix,qs,nx,ny,nz,imin,imax,jmin,jmax,
     1  minrgn,maxrgn,mask,dxm,latg,long,msfx,msfy,usedynamicwts,nwts,
     2  printflag,ic,jc,iprt,Qdir)
! Computes the dynamic weight for index idx proportional to a 
! score based on current agreement with observations.  The 
! observation list is in pireps.inc. npregion is the output 
! number of observations within the input region iregion.
! The index to be scored and assigned a weight is in the
! interpolated qix.  The remapped interpolated index is
! output in qs.  The indices to be scored are in kpickitfa(idx).
! The output weight is in dynamic_wts(iregion,idx) in scorei.inc.
! The number of wts computed for the input region is output in nwts.
! usedynamicwts is output and is set to TRUE when the number
! of observations available for scoring > minObsCalcWts.
! When the number of observations < minObsCalcWts the dynamic_wts
! are set to static_wts. When this routine is called at an 
! The individual indices are read in from the interpolated 4XX.Q
! files on disk in the Qdir.
! qix and qs are work arrays
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nz,maxrgn,minrgn
      real    qix(nx,ny,nz),qs(nx,ny,nz)
      real    dxm
      real    latg(nx,ny),long(nx,ny),msfx(nx,ny),msfy(nx,ny)
      integer iregion
      integer nwts
      logical usedynamicwts
      integer imin,imax,jmin,jmax
      integer ic,jc,printflag,iprt
      character*200 Qdir
      integer mask(nx,ny)
!-----------------------------------------------------------------------
      integer k,idx,idQ,ierr
      integer kmini,kmaxi
      integer npregion
      real    zpmin,zpmax,dzp
      real    TSSmin,TSSmax
      real    wsum
      real    optwts(101)
      include 'scoreparams.inc'
      include 'scorei.inc'
      include 'pireps.inc'   ! npireps
      include 'consts.inc'   ! RMISSD
!-----------------------------------------------------------------------
!
! Initializations
      if(printflag.ge.1) then
        write(iprt,*) 'enter get_dynamic_wts: remap_option=',
     1    remap_option
        write(iprt,*) 'npireps,minObsCalcWts=',npireps,minObsCalcWts
        do iregion=minrgn,maxrgn
        do idx=1,idmax
          if(ipickitfa(iregion,idx).gt.0) then
            write(iprt,*)'iregion,idQ,ipickitfa=',iregion,idx+399,
     1        ipickitfa(iregion,idx)
          endif
        enddo
        enddo
      endif
! Initialize dynamic weights to static weights
      do iregion=minrgn,maxrgn
      do idx=1,idmax
        dynamic_wts(iregion,idx)=RMISSD
        if(ipickitfa(iregion,idx).gt.0) then
          dynamic_wts(iregion,idx)=static_wts(iregion,idx)
        endif
      enddo
      enddo
      npregion=0
      nwts=0
!
! If no pireps over all regions skip
      if(npireps.le.0) then
        if(printflag.ge.1) then
          write(iprt,*) '## No obs for this case: skipping'
        endif
        go to 50
      endif
!
! Loop over all regions
      do iregion=minrgn,maxrgn

!   Set boundries of each index interpolated grid to be scored (.Q)
        kmini=kregion(iregion,1)
        kmaxi=kregion(iregion,2)
        if(printflag.ge.1) then
          write(iprt,*) 'iregion,kmini,kmaxi=',iregion,kmini,kmaxi
        endif
!
!   Filter observations according to the altitude within
!   the input region.  Use observations within dzp ft of
!   the lower and upper boundaries of the region.
        dzp=1000.
        zpmin=MAX(zi(kmini)-dzp,0.)
        zpmax=zi(kmaxi)+dzp
!
        TSSmax=-2.
        TSSmin=+2.
        do idx=1,idmax
          optwts(idx)=RMISSD
          idQ=idx+399
          NNi(idx)=-1
          YNi(idx)=-1
          NYi(idx)=-1
          YYi(idx)=-1
          PODYi(idx)=-1.
          PODNi(idx)=-1.
          TSSi(idx)=-1.
          funci(idx)=+9999.
          if(ipickitfa(iregion,idx).le.0) go to 25
!     Retrieve the index in the appropriate .Q file
          call getqix(qix,nx,ny,nz,idQ,printflag,iprt,Qdir,ierr)
          if(ierr.ne.0) then
            write(iprt,*) 'error reading qix index from disk'
            go to 25
          endif
          if(printflag.ge.3) then
            write(iprt,*) 'after getqix for idQ=',idQ
            do k=1,nz
              write(iprt,*) 'i,j,k,qix=',ic,jc,k,qix(ic,jc,k)
            enddo
          endif
!     Remap the interpploated index to edr and output in the qs array
          iremapd=0  ! remapped flag
          call remapi(qix,qs,nx,ny,nz,iregion,idx,
     1      imin,imax,jmin,jmax,kmini,kmaxi,mask,printflag,ic,jc,iprt)
          iremapd=1  ! indicate already remapped to 0-1
          npregion=0
!     Score the index
          call scoreindx(qs,nx,ny,nz,iregion,idx,zpmin,zpmax,
     1      indx_scoring_option,timog,dxm,latg,long,msfx,msfy,mask,
     2      npregion,printflag,iprt)
          TSSmax=MAX(TSSmax,TSSi(idx))
          TSSmin=MIN(TSSmin,TSSi(idx))
!         if(printflag.ge.2) then
          if(printflag.ge.1) then
            write(iprt,*) 'after scoreindx for iregion,idQ,npregion=',
     1       iregion,idQ,npregion
            if(printflag.ge.3) then
            write(iprt,*) 'idQ tcount,YY,YN,NY,NN,PODY,PODN,TSS,FUNC=',
     1        idQ,tcounti(idx),YYi(idx),YNi(idx),NYi(idx),NNi(idx),
     2        PODYi(idx),PODNi(idx),TSSi(idx),funci(idx)
!           write(iprt,*) 'idQ tcountMS,YYMS,YNMS,NYMS,NNMS=',
!    1        idQ,tcountMSi(idx),YYMSi(idx),YNMSi(idx),NYMSi(idx),NNMSi(idx)
            endif
          endif
          if(funci(idx).gt.0.) optwts(idx)=funci(idx)
   25     continue
!
        enddo  ! index loop
!
!   Use dynamic weights only if the number of observations used
!   in calculating them is > minObsCalcWts.  Otherwise use the static weights.
        if(npregion.lt.minObsCalcWts) then
          if(printflag.ge.1) then
            write(iprt,*)'iregion,npregion,minObsCalcWts=',
     1       iregion,npregion,minObsCalcWts
            write(iprt,*)
     1       '## npregion<minObsCalcWts: using static weights'
          endif
          wsum=0.
          nwts=0
          usedynamicwts=.FALSE.
          do idx=1,idmax
            if(ipickitfa(iregion,idx).gt.0) then
              dynamic_wts(iregion,idx)=static_wts(iregion,idx)
              wsum=wsum+static_wts(iregion,idx)
              nwts=nwts+1
            endif
          enddo  ! index loop
        else
          usedynamicwts=.TRUE.
          wsum=0.
          nwts=0
          do idx=1,idmax
            if(ipickitfa(iregion,idx).gt.0) then
              dynamic_wts(iregion,idx)=optwts(idx)
              wsum=wsum+optwts(idx)
              nwts=nwts+1
            endif
          enddo  ! index loop
!     Normalize the weights
          do idx=1,idmax
            if(dynamic_wts(iregion,idx).ge.0.) 
     1        dynamic_wts(iregion,idx)=dynamic_wts(iregion,idx)/wsum
            if(printflag.ge.3)write(iprt,*)'iregion,idx,wts,normwts=',
     1        idx,optwts(idx),dynamic_wts(iregion,idx)
          enddo
!     Check sums
          wsum=0.
          do idx=1,idmax
            if(dynamic_wts(iregion,idx).ge.0.) 
     1      wsum=wsum+dynamic_wts(iregion,idx)
          enddo
          if(printflag.ge.3) write(iprt,*) 'wsum (should=1)=',wsum
        endif
!
      enddo  ! region loop
!
   50 continue
      if(printflag.ge.1) then
        do iregion=minrgn,maxrgn
        do idx=1,idmax
          idQ=idx+399
          write(iprt,*) 'iregion,idQ,static,dynamic_wts=',iregion,idQ,
     1      static_wts(iregion,idx),dynamic_wts(iregion,idx)
        enddo
        enddo
        write(iprt,*) 'exit  get_dynamic_wts'
      endif
!
      return
      end
!
      subroutine InterpITFAFtoQ(qix,qitfa,zm,pm,Tm,Tv,qvm,hgt,k3d,kta,
     1  nx,ny,nz,nziq,imin,imax,jmin,jmax,kimin,kimax,mask,
     2  comp_ITFAMWT,comp_ITFADYN,printflag,ic,jc,iprt,ioutputflag,
     3  Fdir,Qdir)
! Interpolates each of the ITFAs from native coordinates (stored as
! idQ.F in Fdir to flight levels and stores as idQ.Q in the Qdir.
! idQ=491 for itfa CAT combination using dynamic weights
! idQ=492 for itfa CAT combination using static weights
! idQ=493 for itfa MWT combination using static weights
! idQ=494 for MAX(CAT, MWT) itfa combinations using static weights
! idQ=495 for MAX(CAT, MWT) itfa combinations using dynamic weights
! qitfa, qix, and qit are work arrays.
      implicit none
!-----------------------------------------------------------------------
      include 'scoreparams.inc'
      include 'scorei.inc'
      integer nx,ny,nz,nziq
      real    zm(nx,ny,nz),pm(nx,ny,nz),Tm(nx,ny,nz),qvm(nx,ny,nz),
     1        Tv(nx,ny,nz)
      real    hgt(nx,ny)
      real    qix(nx,ny,nziq)
      integer k3d(nx,ny,nziq)        ! interpolation array
      real    qitfa(nx,ny,nz)
      integer kta
      integer mask(nx,ny)
      logical comp_ITFAMWT,comp_ITFADYN
      integer imin,imax,jmin,jmax,kimin,kimax
      integer ioutputflag
      integer ic,jc,printflag,iprt
      character*200 Fdir,Qdir
      include 'consts.inc'   ! RMISSD
!-----------------------------------------------------------------------
      integer idx,idQ,ierr
      integer nFtoQ
      character*24 cnamei
!-----------------------------------------------------------------------
!
! Initializations
      if(printflag.ge.1) then
        write(iprt,*) 
     1   'enter InterpITFAFtoQ: comp_ITFAMWT,comp_ITFADYN,=',
     2   comp_ITFAMWT,comp_ITFADYN
        write(iprt,*) 'imin,imax,jmin,jmax,kimin,kimax=',
     1   imin,imax,jmin,jmax,kimin,kimax
        write(iprt,*) 'nx,ny,nz,nzi,nziq=',nx,ny,nz,nzi,nziq
      endif
      nFtoQ=0
!
! ITFADYN (idQ=491)
      idx=iditfa
      idQ=399+idx
      cnamei=cname(iditfa)
      if(comp_ITFADYN) then 
!   Read ITFADYN from disk into qitfa
        if(printflag.ge.1) then
          write(iprt,*) 'interpolating ITFADYN: idQ=',idQ
        endif
        call TIinit(qitfa,nx,ny,nz)
        call readF(qitfa,nx,ny,nz,idQ,printflag,iprt,Fdir,ierr)
!   Interpolate native grid ITFADYN in qitfa to flight
!   levels (zi) and output the interpolated grid as qix and
!   also write as idQ.Q file (where idQ=iditfa+399).
        call FtoQ(zm,pm,Tm,qvm,Tv,zi,hgt,pstdi,idQ,imin,imax,
     1    jmin,jmax,kimin,kimax,kta,zsmin,nx,ny,nz,nzi,cnamei,
     2    qitfa,k3d,qix,mask,ic,jc,printflag,iprt,ierr,nFtoQ,Qdir)
!   At low levels the interpolation/extrapolation procedure
!   may provide unreasonable values.  Ensure the values fall
!   within clampitfaL-clampitfaH.  Overwrite the idQ.Q file with the
!   clamped values.  qix is a work array.
        call checkITFAQ(qix,idQ,nx,ny,nzi,imin,imax,jmin,jmax,
     1    kimin,kimax,mask,hgt,zi,printflag,ic,jc,iprt,
     2    ioutputflag,cnamei,Qdir)
      else
!   write as missing
        call TIinit(qix,nx,ny,nzi)
        call putqix(qix,nx,ny,nzi,idQ,cnamei,printflag,iprt,
     1    Qdir,ierr)
      endif
!
! ITFADEF (idQ=492)
      idx=iditfad
      idQ=399+idx
      cnamei=cname(iditfad)
      if(printflag.ge.1) then
        write(iprt,*)'interpolating ITFADEF: idQ=',idQ
      endif
! Read ITFADEF from disk into qitfa
      call TIinit(qitfa,nx,ny,nz)
      call readF(qitfa,nx,ny,nz,idQ,printflag,iprt,Fdir,ierr)
        write(iprt,*) 'after readF for itfadef ic,jc,printflag=',
     1   ic,jc,printflag
! Interpolate native grid ITFADEF in qitfa to flight
! levels (zi) and output the interpolated grid as qix and
! also write as idQ.Q file (where idQ=iditfa+399).
      call FtoQ(zm,pm,Tm,qvm,Tv,zi,hgt,pstdi,idQ,imin,imax,
     1  jmin,jmax,kimin,kimax,kta,zsmin,nx,ny,nz,nzi,cnamei,
     2  qitfa,k3d,qix,mask,ic,jc,printflag,iprt,ierr,nFtoQ,Qdir)
        write(iprt,*) 'after FtoQ for itfadef ic,jc,printflag=',
     1   ic,jc,printflag
! At low levels the interpolation/extrapolation procedure
! may provide unreasonable values.  Ensure the values fall
! within clampitfaL-clampitfaH.  Overwrite the idQ.Q file with the
! clamped values.  qix is a work array.
      call checkITFAQ(qix,idQ,nx,ny,nzi,imin,imax,jmin,jmax,
     1  kimin,kimax,mask,hgt,zi,printflag,ic,jc,iprt,
     2  ioutputflag,cnamei,Qdir)
        write(iprt,*) 'after checkITFAQ for itfadef ic,jc,printflag=',
     1   ic,jc,printflag
!
! ITFAMWT (idQ=493)
      idx=iditfam
      idQ=399+idx
      cnamei=cname(iditfam)
      if(comp_ITFAMWT) then
!   Read ITFAMWT from disk into qitfa
        if(printflag.ge.1) then
          write(iprt,*)'interpolating ITFAMWT: idQ=',idQ
        endif
        call TIinit(qitfa,nx,ny,nz)
        call readF(qitfa,nx,ny,nz,idQ,printflag,iprt,Fdir,ierr)
        write(iprt,*) 'after readF for itfamwt ic,jc,printflag=',
     1   ic,jc,printflag
!   Interpolate native grid ITFAMWT in qitfa to flight
!   levels (zi) and output the interpolated grid as qix and
!   also write as idQ.Q file (where idQ=iditfa+399).
        call FtoQ(zm,pm,Tm,qvm,Tv,zi,hgt,pstdi,idQ,imin,imax,
     1    jmin,jmax,kimin,kimax,kta,zsmin,nx,ny,nz,nzi,cnamei,
     2    qitfa,k3d,qix,mask,ic,jc,printflag,iprt,ierr,nFtoQ,Qdir)
        write(iprt,*) 'after FtoQ for itfamwt ic,jc,printflag=',
     1   ic,jc,printflag
!   At low levels the interpolation/extrapolation procedure
!   may provide unreasonable values.  Ensure the values fall
!   within clampitfaL-clampitfaH.  Overwrite the idQ.Q file with the
!   clamped values.  qix is a work array.
        call checkITFAQ(qix,idQ,nx,ny,nzi,imin,imax,jmin,jmax,
     1    kimin,kimax,mask,hgt,zi,printflag,ic,jc,iprt,
     2    ioutputflag,cnamei,Qdir)
        write(iprt,*) 'after checkITFAQ for itfamwt ic,jc,printflag=',
     1   ic,jc,printflag
      else
!   write as missing
        call TIinit(qix,nx,ny,nzi)
        call putqix(qix,nx,ny,nzi,idQ,cnamei,printflag,iprt,
     1    Qdir,ierr)
      endif
!
! ITFAMAX (idQ=494)
      idx=iditfax
      idQ=399+idx
      cnamei=cname(iditfax)
      if(comp_ITFAMWT) then
!   Read ITFAMAX from disk into qitfa
        if(printflag.ge.1) then
          write(iprt,*)'interpolating ITFAMAX: idQ=',idQ
        endif
        call TIinit(qitfa,nx,ny,nz)
        call readF(qitfa,nx,ny,nz,idQ,printflag,iprt,Fdir,ierr)
        write(iprt,*) 'after readF for itfamax ic,jc,printflag=',
     1   ic,jc,printflag
!   Interpolate native grid ITFAMAX in qitfa to flight
!   levels (zi) and output the interpolated grid as qix and
!   also write as idQ.Q file (where idQ=iditfa+399).
        call FtoQ(zm,pm,Tm,qvm,Tv,zi,hgt,pstdi,idQ,imin,imax,
     1    jmin,jmax,kimin,kimax,kta,zsmin,nx,ny,nz,nzi,cnamei,
     2    qitfa,k3d,qix,mask,ic,jc,printflag,iprt,ierr,nFtoQ,Qdir)
        write(iprt,*) 'after FtoQ for itfamax ic,jc,printflag=',
     1   ic,jc,printflag
!   At low levels the interpolation/extrapolation procedure
!   may provide unreasonable values.  Ensure the values fall
!   within clampitfaL-clampitfaH.  Overwrite the idQ.Q file with the
!   clamped values.  qix is a work array.
        call checkITFAQ(qix,idQ,nx,ny,nzi,imin,imax,jmin,jmax,
     1    kimin,kimax,mask,hgt,zi,printflag,ic,jc,iprt,
     2    ioutputflag,cnamei,Qdir)
        write(iprt,*) 'after checkITFAQ for itfamax ic,jc,printflag=',
     1   ic,jc,printflag
      else
!   write out interpolated ITFADEF
        call getqix(qix,nx,ny,nzi,399+iditfad,printflag,iprt,Qdir,ierr)
        call putqix(qix,nx,ny,nzi,idQ,cnamei,printflag,iprt,
     1    Qdir,ierr)
      endif
!
! ITFAMAXD (idQ=495)
      idx=iditfaxd
      idQ=399+idx
      cnamei=cname(iditfaxd)
      if(.NOT.comp_ITFADYN) then
!   Read ITFAMAXD from disk into qitfa
        if(printflag.ge.1) then
          write(iprt,*)'interpolating ITFAMAXD: idQ=',idQ
        endif
        call TIinit(qitfa,nx,ny,nz)
        call readF(qitfa,nx,ny,nz,idQ,printflag,iprt,Fdir,ierr)
!   Interpolate native grid ITFAMAXD in qitfa to flight
!   levels (zi) and output the interpolated grid as qix and
!   also write as idQ.Q file (where idQ=iditfa+399).
        call FtoQ(zm,pm,Tm,qvm,Tv,zi,hgt,pstdi,idQ,imin,imax,
     1    jmin,jmax,kimin,kimax,kta,zsmin,nx,ny,nz,nzi,cnamei,
     2    qitfa,k3d,qix,mask,ic,jc,printflag,iprt,ierr,nFtoQ,Qdir)
!   At low levels the interpolation/extrapolation procedure
!   may provide unreasonable values.  Ensure the values fall
!   within clampitfaL-clampitfaH.  Overwrite the idQ.Q file with the
!   clamped values.  qix is a work array.
        call checkITFAQ(qix,idQ,nx,ny,nzi,imin,imax,jmin,jmax,
     1    kimin,kimax,mask,hgt,zi,printflag,ic,jc,iprt,
     2    ioutputflag,cnamei,Qdir)
      else
!   write as missing
        call TIinit(qix,nx,ny,nzi)
        call putqix(qix,nx,ny,nzi,idQ,cnamei,printflag,iprt,
     1    Qdir,ierr)
      endif
!
      if(printflag.ge.1) then
        write(iprt,*)'exit  InterpITFAFtoQ'
      endif
!
      return
      end
!
      subroutine checkITFAQ(qitfa,idQ,nx,ny,nzi,imin,imax,
     1 jmin,jmax,kmin,kmax,mask,hgt,zi,printflag,ic,jc,iprt,
     2 ioutputflag,cnamei,Qdir)
! Checks ITFA forecasted values to make sure they are in the
! expected range (clampitfaL,clampitfaH).  Also checks to make
! values below the terrain are set to missing.
      implicit none
!-----------------------------------------------------------------------
      integer nx,ny,nzi
      real    qitfa(nx,ny,nzi)
      real    zi(nzi)
      real    hgt(nx,ny)
      integer mask(nx,ny)
      integer idQ
      integer imin,imax,jmin,jmax,kmin,kmax
      integer ioutputflag
      integer ic,jc,printflag,iprt
      character*24 cnamei
      character*200 Qdir
!-----------------------------------------------------------------------
      integer i,j,k,ierr
      real    qijk,qmin,qmax
      real    hft
      integer nclamped
      include 'consts.inc'   ! RMISSD
      real    clampitfaL,clampitfaH
      parameter (clampitfaL=0.0,clampitfaH=1.0)
!-----------------------------------------------------------------------
!
! Initializations
      nclamped=0
      if(printflag.ge.2) then
        write(iprt,*) 'enter checkITFAQ: idQ,ioutputflag=',
     1    idQ,ioutputflag
      endif
!
      if(ioutputflag.le.0) then
        write(iprt,*) 'ioutputflag.le.0 - no files written - return'
        return
      endif
!
! Retrieve the itfa idQ.Q file 
      call TIinit(qitfa,nx,ny,nzi)
      call getqix(qitfa,nx,ny,nzi,idQ,printflag,iprt,QDir,ierr)
      if(ierr.ne.0) then
        write(iprt,*) 'error in checkITFAQ reading idQ=',idQ
        go to 50
      endif
!
! Check the readin values
      qmin=+1.0E20
      qmax=-1.0E20
      do j=jmin,jmax
      do i=imin,imax
        if(mask(i,j).le.0) go to 40
        do k=kmin,kmax
          qijk=qitfa(i,j,k)
          if(ABS(qijk-RMISSD).LT.1.0E-3) go to 35
          if(qijk.lt.clampitfaL) then
            qijk=clampitfaL
            nclamped=nclamped+1
          endif
          if(qijk.gt.clampitfaH) then
            qijk=clampitfaH
            nclamped=nclamped+1
          endif
          qitfa(i,j,k)=qijk
          qmin=MIN(qmin,qijk)
          qmax=MAX(qmax,qijk)
   35     continue
        enddo
!   Ensure values below the terrain are set to missing
        do k=1,kmax
!     zi is in ft, terrain hgt(i,j) is in m, so convert
          hft=hgt(i,j)*3.28
          if(zi(k).lt.hft-100.) then
            qijk=qitfa(i,j,k)
            if(ABS(qijk-RMISSD).GT.1.0E-3) then
              qitfa(i,j,k)=RMISSD
              nclamped=nclamped+1
            endif
          else
            go to 40
          endif
        enddo
   40 continue
      enddo
      enddo
!
! Write back if adjustments were made
      if(nclamped.gt.0) then
        call putqix(qitfa,nx,ny,nzi,idQ,cnamei,printflag,iprt,
     1    Qdir,ierr)
      endif
!
   50 continue
      if(printflag.ge.2) then
        write(iprt,*) 'exit  checkITFAQ: idQ,nclamped=',idQ,nclamped
        write(iprt,*) 'qmin,qmax=',qmin,qmax
      endif
!
      return
      end
