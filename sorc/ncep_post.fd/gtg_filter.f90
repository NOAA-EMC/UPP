module gtg_filter
  use ctlblk_mod, only: SPVAL
  use ctlblk_mod, only: jsta,jend,jsta_2l, jend_2u, jsta_m, jend_m, &
       jsta_m2, jend_m2,im,jm,lm, modelname,global
  use gtg_config, only : SMALL1

  implicit none

  integer, parameter :: nmax = 100

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
    real, dimension(IM,jsta_2l:jend_2u,LM), intent(inout) :: u,um

    integer :: i,j,k,ift
    integer :: im1,ip1,jm1,jp1

!   --- xy filter loop
    do ift = 1,nftxy
       do k=kmin,kmax

!         --- Initialize dummy um
          do j = jsta, jend
          do i = 1,IM
             um(i,j,k)=u(i,j,k)
          end do
          end do

!         --- 1-2-1 smoothing in y direction first 
          do j = jsta_m, jend_m ! No smoothing on j=1 or j=JM point
             jm1=j-1
             jp1=j+1
             do i=1,IM
!               --- Don't include uncomputed u(i,j)
                if(.not.(ABS(u(i,jp1,k)-SPVAL)<SMALL1 .or. &
                         ABS(u(i,j  ,k)-SPVAL)<SMALL1 .or. &
                         ABS(u(i,jm1,k)-SPVAL)<SMALL1)) then
                   um(i,j,k)=0.25*(u(i,jp1,k)+2.*u(i,j,k)+u(i,jm1,k))
                end if
             enddo
          enddo  ! j loop

!         --- Smooth in x direction. Smooth first and last points only if cyclic.
          do i=1,IM

             im1=i-1
             ip1=i+1
             if(im1<1) then
                if(modelname == 'GFS' .or. global) then
                   im1=im1+IM
                else
                   im1=1
                end if
             end if
             if(ip1>IM) then
                if(modelname == 'GFS' .or. global) then
                   ip1=ip1-IM
                else
                   ip1=IM
                end if
             endif

             do j=jsta, jend
!               --- Don't include uncomputed u(i,j)
                if(.not. (ABS(um(im1,j,k)-SPVAL)<SMALL1 .or. &
                          ABS(um(i,  j,k)-SPVAL)<SMALL1 .or. &
                          ABS(um(ip1,j,k)-SPVAL)<SMALL1)) then
                   u(i,j,k) = 0.25*(um(ip1,j,k) +2.*um(i,j,k) + um(im1,j,k))
                endif
             enddo  ! j loop
          end do ! i loop

       enddo  ! k loop
    enddo  ! nftxy loop

!   --- z filter loop
    do ift = 1,nftz

       do j=jsta,jend
       do i=1,im

!         --- Initialize dummy um
          do k=kmin,kmax
             um(i,j,k)=u(i,j,k)
          enddo

!         --- Smooth in z direction.
          do k=kmin+1,kmax-1 ! Don't include boundary points
!            --- Don't include uncomputed u(i,j)
             if(.not. (ABS(um(i,j,k-1)-SPVAL)<SMALL1 .or. &
                       ABS(um(i,j,k)  -SPVAL)<SMALL1 .or. &
                       ABS(um(i,j,k+1)-SPVAL)<SMALL1)) then
                u(i,j,k)=0.25*(um(i,j,k+1) + 2.*um(i,j,k) + um(i,j,k-1))
             endif
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
    real, dimension(IM,jsta_2l:jend_2u), intent(inout) :: u, um

    integer :: i,j,ift
    integer :: im1,ip1,jm1,jp1

    do ift = 1,nftxy    ! xy filter loop

!      --- Initialize dummy um
       do j = jsta, jend
       do i = 1,IM
          um(i,j)=u(i,j)
       end do
       end do

!      --- 1-2-1 smoothing in y direction first
       do j = jsta_m, jend_m ! No smoothing on j=1 or j=JM point
          jm1=j-1
          jp1=j+1
          do i=1,IM
!            --- Don't include uncomputed u(i,j)
             if(.not. (ABS(u(i,jp1)-SPVAL)<SMALL1 .or. &
                       ABS(u(i,j  )-SPVAL)<SMALL1 .or. &
                       ABS(u(i,jm1)-SPVAL)<SMALL1)) then
                um(i,j)=0.25*(u(i,jp1)+2.*u(i,j)+u(i,jm1))
             endif
          enddo
       enddo  ! j loop

!      --- Smooth in x direction. Smooth first and last points only if cyclic.
       do i=1,IM

          im1=i-1
          ip1=i+1
          if(im1<1) then
             if(modelname == 'GFS' .or. global) then
                im1=im1+IM
             else
                im1=1
             end if
          end if
          if(ip1>IM) then
             if(modelname == 'GFS' .or. global) then
                ip1=ip1-IM
             else
                ip1=IM
             end if
          endif

          do j=jsta, jend
!            --- Don't include uncomputed u(i,j)
             if(.not. (ABS(um(im1,j)-SPVAL)<SMALL1 .or. &
                       ABS(um(i,  j)-SPVAL)<SMALL1 .or. &
                       ABS(um(ip1,j)-SPVAL)<SMALL1)) then
                u(i,j) = 0.25*(um(ip1,j) +2.*um(i,j) + um(im1,j))
             endif
          enddo ! j loop
       enddo    ! i loop

    enddo  ! nftxy loop

    return
  end subroutine meanFilter2D

!-----------------------------------------------------------------------
  subroutine medianFilter3D(kmin,kmax,u,um)
! --- Applies median filter to input array u and outputs in array um

    use Quicksort

    implicit none

    integer, intent(in) :: kmin,kmax
    real, dimension(IM,jsta_2l:jend_2u,LM), intent(inout) :: u,um

    integer, parameter :: grid_extent=1
    
    integer :: i,j,k,ii,iii,jj,kk
    real :: ua(nmax)
    integer :: n,nh
    real :: umedian

!   --- Copy u into work array um
    um = u

    do k=kmin,kmax
    do j=jsta,jend
    do i=1,IM

       do n=1,nmax
          ua(n)=SPVAL
       enddo

       n=0
       do kk=k-grid_extent,k+grid_extent
       do jj=j-grid_extent,j+grid_extent
       do iii=i-grid_extent,i+grid_extent

          if ( jj < 1 .or. jj > JM .or. &
               kk < 1 .or. kk > LM ) cycle

          ii=iii
          if(iii<1) then
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
          u(i,j,k)=um(i,j,k)
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
          u(i,j,k)=umedian
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

    real, dimension(IM,jsta_2l:jend_2u), intent(inout) :: u, um

    integer, parameter :: grid_extent=1

    integer :: i,j,ii,iii,jj
    real :: ua(nmax)
    integer :: n,nh
    real  ::  umedian

!   --- Copy u into work array um
    um =u

    do j=jsta,jend
    do i=1,IM

       do n=1,nmax
          ua(n)=SPVAL
       enddo

       n=0
       do jj=j-grid_extent,j+grid_extent
       do iii=i-grid_extent,i+grid_extent

          if (jj < 1 .or. jj > JM) cycle

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
          u(i,j)=um(i,j)
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
          u(i,j)=umedian
       endif
    enddo
    enddo

    return
  end subroutine medianFilter2D

end module gtg_filter
