module gtg_filter
  use ctlblk_mod, only: SPVAL
  use ctlblk_mod, only: jsta,jend,jsta_2l, jend_2u, jsta_m, jend_m, &
       jsta_m2, jend_m2,im,jm,lm, modelname,global
  use gtg_config, only : SMALL1

  implicit none

  integer, parameter :: grid_extent=1
  logical, parameter :: overwrite=.true.
  integer, parameter :: nmax = 100

! --- Specify x boundary conditions options (mbc)
! --- 1=No smoothing of end points,
! --- 2=Assume 0 gradients at end points,
! --- 3=cyclic in x ((nx+1,j,k)=u(1,j,k), u(0,j,k)=u(nx,j,k))
  integer, parameter :: DIRICHLET=1, NEUMANN=2, CYCLIC=3

contains

!-----------------------------------------------------------------------
  subroutine filt3d(kmin,kmax,nftxy,nftz,filttype,u)

!$$$ SUBPROGRAM DOCUMENTATION BLOCK 
! ABSTRACT: Applies filter to input array u and overwrites u with
!    filtered um on output.  um is a work array.
!    filttype=1 specifies 1-2-1 smoothing in each direction (default)
!    filttype=2 specifies median filter
!$$$

    implicit none

    integer, intent(in) :: kmin,kmax
    integer, intent(in) :: nftxy,nftz
    integer, intent(in) :: filttype
    real, intent(inout) :: u(im,jsta_2l:jend_2u,lm)

    real :: um(im,jsta_2l:jend_2u,lm) ! a work array
    
!   --- filttype=1 specifies 1-2-1 smoothing in each direction (default)
!   --- filttype=2 specifies median filter
    if(filttype == 2) then
       call medianFilter3D(kmin,kmax,u,um)
    else
       call meanFilter3D(kmin,kmax,nftxy,nftz,u,um)
    endif

!   --- overwrite input u with um for output.  This is for compability 
!   --- with previous versions
    if(overwrite) then
       u = um
    endif

    return
  end subroutine filt3d

!-----------------------------------------------------------------------
  subroutine meanFilter3D(kmin,kmax,nftxy,nftz,u,um)

!$$$  SUBPROGRAM DOCUMENTATION BLOCK
! ABSTRACT: Applies 2d 1-2-1 smoother nftxy times in x,y at each vertical
!   level and nftz times in z at each (i,j) column.  Equivalent 
!   to a simple 27-point filter if nftxy=nftz=1.
!   E.g., Haltiner and Williams (1980) eqs. (11-107,11-108).
!$$$

    implicit none
    integer, intent(in) :: kmin,kmax
    integer, intent(in) :: nftxy,nftz
    real, dimension(IM,jsta_2l:jend_2u,LM), intent(in) :: u
    real, dimension(IM,jsta_2l:jend_2u,LM), intent(inout) :: um

    integer :: i,j,k,ift
    integer kmin1,kmax1,imin1,imax1,jmin1,jmax1
    integer :: im1,ip1
    real :: umm1,ump1
    real :: usave,ulast

!   --- Copy u into work array um
    um=u

!   --- Ensure smoothing function extends only to within grid extent pts
!   --- from the grid boundaries.
    imin1 = 1  + grid_extent
    imax1 = IM - grid_extent
    if(grid_extent == 2) then
       jmin1=jsta_m2
       jmax1=jend_m2
    else
       jmin1=jsta_m
       jmax1=jend_m
    end if
    kmin1=MAX(kmin,1+grid_extent)
    kmax1=MIN(kmax,LM-grid_extent)

!   --- xy filter loop
    do ift = 1,nftxy
       do k=kmin,kmax

!         --- Smooth in y direction. Note no smoothing on first and last points.
          do i=1,IM
             ulast=um(i,jmin1-1,k)
             do j=jmin1,jmax1
                usave=um(i,j,k)
!               --- Don't include uncomputed u(i,j)
                if(.not.(ABS(um(i,j+1,k)-SPVAL)<SMALL1 .or. &
                         ABS(usave-SPVAL)<SMALL1 .or. &
                         ABS(ulast-SPVAL)<SMALL1)) then
                   um(i,j,k)=.25*(um(i,j+1,k) +2.*usave + ulast)
                end if
                ulast=usave
             enddo
!            --- Smooth j=1 point assuming u(0)=u(1)
             if(jsta==1) then
                if(.not. (ABS(um(i,1,k)-SPVAL)<SMALL1 .or. &
                          ABS(um(i,2,k)-SPVAL)<SMALL1)) then
                   umm1=um(i,1,k)
                   um(i,1,k)=.25*(umm1 + 2.*um(i,1,k) + um(i,2,k))
                endif
             end if
!            --- Smooth j=jm point assuming u(jm+1)=u(jm)
             if(jend==jm) then
                if(.not.(ABS(um(i,jm,  k)-SPVAL)<SMALL1 .or. &
                         ABS(um(i,jm-1,k)-SPVAL)<SMALL1)) then
                   ump1=um(i,jm,k)
                   um(i,jm,k)=.25*(um(i,jm-1,k) + 2.*um(i,jm,k) + ump1)
                endif
             end if
          enddo  ! i loop

!         --- Smooth in x direction. Smooth first and last points only if cyclic.
          do j=jsta, jend
             im1=imin1-1 ! will be >= 1 since imin1 = 1  + grid_extent
             ulast=um(im1,j,k)
             do i=imin1,imax1
                usave=um(i,j,k)
!               --- u(i,j) = .25*(u(i+1,j) +2.*u(i,j) + u(i-1,j))
!               --- Don't include uncomputed u(i,j)
                if(.not. (ABS(um(i+1,j,k)-SPVAL)<SMALL1 .or. &
                          ABS(usave      -SPVAL)<SMALL1 .or. &
                          ABS(ulast      -SPVAL)<SMALL1)) then
                   um(i,j,k) = .25*(um(i+1,j,k) +2.*usave + ulast)
                endif
                ulast=usave
             enddo  ! i loop
!            --- Smooth i=1 pt with BC depending on the value of mbc
             if(modelname == 'GFS' .or. global) then  ! assumes cyclic u(0)=u(im)
                umm1=um(im,j,k)
             else               ! assumes u(0)=u(1)
                umm1=um(1,j,k)
             endif
             if(.not. (ABS(um(1,j,k)-SPVAL)<SMALL1 .or. &
                       ABS(um(2,j,k)-SPVAL)<SMALL1 .or. &
                       ABS(umm1-SPVAL)<SMALL1)) then
                um(1,j,k)=.25*(umm1 + 2.*um(1,j,k) + um(2,j,k))
             end if
!           --- Smooth i=im pt with BC depending on the value of mbc
             if(modelname == 'GFS' .or. global) then  ! assumes cyclic u(im+1)=u(1)
                ump1=um(1,j,k)
             else               ! assumes u(im+1)=u(im)
                ump1=um(im,j,k)
             endif
             if(.not. (ABS(um(im,  j,k)-SPVAL)<SMALL1 .or. &
                       ABS(um(im-1,j,k)-SPVAL)<SMALL1 .or. &
                       ABS(ump1-SPVAL)<SMALL1)) then
                um(im,j,k)=.25*(um(im-1,j,k) + 2.*um(im,j,k) + ump1)
             endif
          enddo  ! j loop

       enddo  ! k loop
    enddo  ! nftxy loop

!   --- z filter loop
    do ift = 1,nftz
       do j=jsta,jend
       do i=1,im

!         --- Smooth in z direction. Note 2pt smoothing on first and last points.
          ulast=um(i,j,kmin1-1)
          do k=kmin1,kmax1
             usave=um(i,j,k)
!            --- Don't include uncomputed u(i,j)
             if(.not. (ABS(ulast      -SPVAL)<SMALL1 .or. &
                       ABS(usave      -SPVAL)<SMALL1 .or. &
                       ABS(um(i,j,k+1)-SPVAL)<SMALL1)) then
!               --- u(i,j) = .25*(u(i,j,k+1) +2.*u(i,j,k) + u(i,j,k-1))
                um(i,j,k)=.25*(um(i,j,k+1) + 2.*usave + ulast)
             endif
             ulast=usave
          enddo  ! k loop

       enddo  ! i loop
       enddo  ! j loop
    enddo  ! nftz loop

    return
  end subroutine meanFilter3D

!-----------------------------------------------------------------------
  subroutine filt2d(nftxy,filttype,u)

!$$$  SUBPROGRAM DOCUMENTATION BLOCK 
! ABSTRACT: Applies filter to input 2d array u and overwrites u with
!   filtered um on output.  um is a work array.
!   filttype=1 specifies 1-2-1 smoothing in each direction (default)
!   filttype=2 specifies median filter
!$$$

    implicit none
    integer, intent(in) :: nftxy
    integer, intent(in) :: filttype
    real, intent(inout) :: u(im,jsta_2l:jend_2u)

    real :: um(im,jsta_2l:jend_2u) ! a work array

!   --- filttype=1 specifies 1-2-1 smoothing in each direction (default)
!   --- filttype=2 specifies median filter
    if(filttype==2) then
       call medianFilter2D(u,um)
    else
       call meanFilter2D(nftxy,u,um)
    endif

!   --- overwrite input u with um for output.  This is for compability 
!   --- with previous versions
    if(overwrite) then
       u = um
    endif
!
    return
  end subroutine filt2d

!-----------------------------------------------------------------------
  subroutine meanFilter2D(nftxy,u,um)

!$$$ SUBPROGRAM DOCUMENTATION BLOCK 
! ABSTRACT: uses a simple 3 point smoothing function on on 2-D arrays
!     --- Note smoothing in x followed by a smoothing in y is equivalent
!     --- to a 9-pt filter.
!     --- E.g., Haltiner and Williams (1980) eqs. (11-107,11-108).
!$$$

    implicit none

    integer, intent(in) :: nftxy
    real, intent(inout) :: u(im,jsta_2l:jend_2u)
    real, intent(inout) :: um(im,jsta_2l:jend_2u) ! a work array

    integer :: i,j,ift
    integer :: imin1,imax1,jmin1,jmax1
    integer :: im1,ip1
    real :: umm1,ump1
    real :: usave,ulast

!   --- Initializations
    imin1 = 1  + grid_extent
    imax1 = IM - grid_extent
    if(grid_extent == 2) then
       jmin1=jsta_m2
       jmax1=jend_m2
    else
       jmin1=jsta_m
       jmax1=jend_m
    end if

!   --- Copy u into work array um
    um = u
  
    do ift = 1,nftxy    ! xy filter loop

!      --- Smooth in y direction.
       do i=1,im
          ulast=um(i,jmin1-1)
          do j=jmin1,jmax1
             usave=um(i,j)
!            --- Don't include uncomputed u(i,j)
             if(.not. (ABS(um(i,j+1)-SPVAL)<SMALL1 .or. &
                       ABS(um(i,j  )-SPVAL)<SMALL1 .or. &
                       ABS(ulast    -SPVAL)<SMALL1)) then
                um(i,j)=.25*(um(i,j+1) +2.*um(i,j) + ulast)
             endif
             ulast=usave
          enddo
!         --- Smooth j=1 point assuming u(0)=u(1)
          if( jsta==1) then
             if(.not. (ABS(um(i,1)-SPVAL)<SMALL1 .or. &
                       ABS(um(i,2)-SPVAL)<SMALL1)) then
                umm1=um(i,1)
                um(i,1)=.25*(umm1 + 2.*um(i,1) + um(i,2))
             endif
          endif
!         --- Smooth j=jm point assuming u(jm+1)=u(jm)
          if(jend==jm) then
             if(.not. (ABS(um(i,jm)-SPVAL)<SMALL1 .or. &
                       ABS(um(i,jm-1)-SPVAL)<SMALL1)) then
                ump1=um(i,jm)
                um(i,jm)=.25*(um(i,jm-1) + 2.*um(i,jm) + ump1)
             endif
          endif
       enddo  ! i loop

!      --- Smooth in x direction. Account for possible cyclic BC.
       do j=jsta, jend
          im1=imin1-1 ! will be >= 1 since imin1 = 1  + grid_extent
          ulast=um(im1,j)
          do i=imin1,imax1
             usave=um(i,j)
!            --- u(i,j) = .25*(u(i+1,j) +2.*u(i,j) + u(i-1,j))
!            --- Don't include uncomputed u(i,j)
             if(.not. (ABS(um(i+1,j)-SPVAL)<SMALL1 .or. &
                       ABS(um(i,j  )-SPVAL)<SMALL1 .or. &
                       ABS(ulast    -SPVAL)<SMALL1)) then
                um(i,j) = .25*(um(i+1,j) +2.*um(i,j) + ulast)
             endif
             ulast=usave
          enddo
!         --- Smooth i=1 pt with BC depending on the value of mbc
          if(modelname == 'GFS' .or. global) then  ! assumes cyclic u(0)=u(im)
             umm1=um(im,j)
          else               ! assumes u(0)=u(1)
             umm1=um(1,j)
          endif
          if(.not. (ABS(um(1,j)-SPVAL)<SMALL1 .or. &
                    ABS(um(2,j)-SPVAL)<SMALL1 .or. &
                    ABS(umm1-SPVAL)<SMALL1)) then
             um(1,j)=.25*(umm1 + 2.*um(1,j) + um(2,j))
          end if
!         --- Smooth i=im pt with BC depending on the value of mbc
          if(modelname == 'GFS' .or. global) then  ! assumes cyclic u(im+1)=u(1)
             ump1=um(1,j)
          else               ! assumes u(im+1)=u(im)
             ump1=um(im,j)
          endif
          if(.not. (ABS(um(im,j  )-SPVAL)<SMALL1 .or. &
                    ABS(um(im-1,j)-SPVAL)<SMALL1 .or. &
                    ABS(ump1-SPVAL)<SMALL1)) then
             um(im,j)=.25*(um(im-1,j) + 2.*um(im,j) + ump1)
          endif
       enddo  ! j loop

    enddo  ! nftxy loop

    return
  end subroutine meanFilter2D

!-----------------------------------------------------------------------
  subroutine medianFilter3D(kmin,kmax,u,um)
! --- Applies median filter to input array u and outputs in array um

    use Quicksort

    implicit none

    integer, intent(in) :: kmin,kmax
    real, dimension(IM,jsta_2l:jend_2u,LM), intent(in) :: u
    real, dimension(IM,jsta_2l:jend_2u,LM), intent(inout) :: um
    
    integer :: i,j,k,ii,iii,jj,kk
    integer :: imin1,imax1,jmin1,jmax1,kmin1,kmax1
    real :: ua(nmax)
    integer :: n,nh
    real :: umedian

!     --- Initializations
    imin1=1
    imax1=IM
    jmin1=jsta
    jmax1=jend
    kmin1=kmin
    kmax1=kmax
!   --- Copy u into work array um
    um = u
    do n=1,nmax
       ua(n)=SPVAL
    enddo

    do k=kmin1,kmax1
    do j=jmin1,jmax1
    do i=imin1,imax1
       n=0
       do kk=k-grid_extent,k+grid_extent
       do jj=j-grid_extent,j+grid_extent
       do iii=i-grid_extent,i+grid_extent

          if ( jj < jsta_2l .or. jj > jend_2u .or. &
               kk < 1 .or. kk > LM ) cycle

          ii=iii
          if(iii<=0) then
!            --- cyclic boundary conditions u(0)=u(im), etc.
             if(modelname == 'GFS' .or. global) then
                ii=iii+IM
             else
                cycle
             endif
          endif
          if(iii>IM) then
!            --- cyclic boundary conditions u(im+1)=u(1), etc.
             if(modelname == 'GFS' .or. global) then 
                ii=iii-IM
             else
                cycle
             endif
          endif
          if(.not. (ABS(u(ii,jj,kk)-SPVAL)<SMALL1)) then
             n=n+1
             ua(n)=u(ii,jj,kk)
          end if
       enddo
       enddo
       enddo
       if(n<=1) then
          um(i,j,k)=u(i,j,k)
       else
!         --- Sort ua into ascending order
          call quick_sort(ua,n)
!         --- Capture median
          nh=n/2
          if(2*nh==n) then
             umedian=0.5*(ua(nh)+ua(nh+1))
          else
             umedian=ua(nh+1)
          endif
          um(i,j,k)=umedian
       endif
    enddo
    enddo
    enddo

    return
  end subroutine medianFilter3D

!-----------------------------------------------------------------------
  subroutine medianFilter2D(u,um)
!   --- Applies median filter to input array u and outputs in array um

    use Quicksort

    implicit none

    real, dimension(IM,jsta_2l:jend_2u), intent(in) :: u
    real, dimension(IM,jsta_2l:jend_2u), intent(inout) :: um

    integer :: i,j,ii,iii,jj
    integer :: imin1,imax1,jmin1,jmax1
    real :: ua(nmax)
    integer :: n,nh
    real  ::  umedian

!   --- Initializations
    imin1=1
    imax1=IM
    jmin1=jsta
    jmax1=jend
!   --- Copy u into work array um
    um =u
    do n=1,nmax
       ua(n)=SPVAL
    enddo

    do j=jmin1,jmax1
    do i=imin1,imax1
       n=0
       do jj=j-grid_extent,j+grid_extent
       do iii=i-grid_extent,i+grid_extent

          if ( jj < jsta_2l .or. jj > jend_2u) cycle

          ii=iii
          if(iii<=0) then
!            --- cyclic boundary conditions u(0)=u(im), etc.
             if(modelname == 'GFS' .or. global) then
                ii=iii+IM
             else
                cycle
             endif
          endif
          if(iii>IM) then
!           --- cyclic boundary conditions u(im+1)=u(1), etc.
             if(modelname == 'GFS' .or. global) then
                ii=iii-IM
             else
                cycle
            endif
          endif
          if(.not. (u(ii,jj)-SPVAL)<SMALL1) then
             n=n+1
             ua(n)=u(ii,jj)
          end if
       enddo
       enddo
       if(n<=1) then
          um(i,j)=u(i,j)
       else
!         --- Sort ua into ascending order
          call quick_sort(ua,n)
!         --- Capture median
          nh=n/2
          if(2*nh==n) then
             umedian=0.5*(ua(nh)+ua(nh+1))
          else
             umedian=ua(nh+1)
          endif
          um(i,j)=umedian
       endif
    enddo
    enddo

    return
  end subroutine medianFilter2D

end module gtg_filter
