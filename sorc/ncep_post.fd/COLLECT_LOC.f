!> @file
!> @brief Subroutine that collect gathers from all MPI tasks.
!>
!> @param[in] A Array being gathered.
!> @param[out] A gathered array - only valid on task 0.
!>
!> Gather "A" from all MPI tasks onto task 0.
!>
!> ### Program history log:
!> Date       | Programmer          | Comments
!> -----------|---------------------|----------
!> 2000-01-06 | Jim Tuccillo        | Initial
!> 2021-06-01 | George Vandenberghe | 2D Decomposition             
!>
!> @author Jim Tuccillo IBM @date 2000-01-06
      SUBROUTINE COLLECT_LOC ( A, B ) 


      use CTLBLK_mod, only: num_procs, jsta, icnt, idsp, mpi_comm_comp, im,&
                            jsta_2l, jend_2u, jm, me,    &
                            buff,ista_2l,iend_2u,jexa,iexa,jsxa,isxa,ista,iend
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      include 'mpif.h'
      integer ii,jj,isum
      real, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in) :: a
      real, dimension(im,jm), intent(out) :: b
      integer ierr,n
      real, allocatable :: rbufs(:)
      allocate(buff(im*jm))
      jj=( jexa(me)-jsxa(me)+1) * (iexa(me)-isxa(me)+1)  
      allocate( rbufs(( jexa(me)-jsxa(me)+1) * (iexa(me)-isxa(me)+1)) )
!
      if ( num_procs <= 1 ) then
         b = a
      else
          
!GWV   reshape the receive subdomain

        isum=1
        do jj=jsxa(me),jexa(me)
          do ii=isxa(me),iexa(me)
            if(isum .gt. im*jm .or. ii .gt. im .or. ii .lt. 1 .or. jj .gt. jm .or. jj.lt. 1) &
              write(*,901)' BOUNDS2 FAIL in reshape ',isum,ii,jj,im*jm,im,im*jm
              rbufs(isum)=a(ii,jj)
              isum=isum+1
          end do
        end do

!GWV  end reshape

        call mpi_gatherv(rbufs,icnt(me),MPI_REAL,  buff,icnt,idsp,MPI_REAL,0,MPI_COMM_WORLD, ierr ) 

!GWV   reshape the gathered array 

        if(me .eq. 0) then
          isum=1 
          do n=0,num_procs-1
            do jj=jsxa(n),jexa(n)
              do ii=isxa(n),iexa(n)
                if(isum .gt. im*jm .or. ii .gt. im .or. ii .lt. 1 .or. jj .gt. jm .or. jj .lt. 1) &
                 write(*,901)' BOUNDS FAIL in reshape ',isum,ii,jj,im*jm,im,im*jm
                 b(ii,jj)=buff(isum)
                 isum=isum+1
              end do
            end do
          end do
        end if

      endif   !  num_procs <= 1

 901  format(a30,10i10)

      deallocate(buff)
      deallocate(rbufs)

      end               
!
!-----------------------------------------------------------------------
!
      SUBROUTINE COLLECT_ALL ( A, B )

      use CTLBLK_mod, only: num_procs, jsta, icnt, idsp, mpi_comm_comp, im,&
              jsta_2l, jend_2u, jm, me,    &
       buff,ista_2l,iend_2u,jexa,iexa,jsxa,isxa,ista,iend,jend
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      include 'mpif.h'
      integer ii,jj,isum
      real, dimension(ista:iend,jsta:jend), intent(in) :: a
      real, dimension(im,jm), intent(out) :: b
      integer ierr,n
      real, allocatable :: rbufs(:)
      allocate(buff(im*jm))
      jj=( jexa(me)-jsxa(me)+1) * (iexa(me)-isxa(me)+1)
      allocate( rbufs(( jexa(me)-jsxa(me)+1) * (iexa(me)-isxa(me)+1)) )
!
      if ( num_procs <= 1 ) then
         b = a
      else

!GWV   reshape the receive subdomain
        isum=1
        do jj=jsxa(me),jexa(me)
          do ii=isxa(me),iexa(me)
            if(isum .gt. im*jm .or. ii .gt. im .or. ii .lt. 1 .or. jj .gt. jm .or. jj.lt. 1) &
               write(*,901)' BOUNDS2 FAIL in reshape',isum,ii,jj,im*jm,im,im*jm
            rbufs(isum)=a(ii,jj)
            isum=isum+1
          end do
        end do
!GWV  end reshape

        call mpi_allgatherv(rbufs,icnt(me),MPI_REAL,buff,icnt,idsp,MPI_REAL, mpi_comm_comp, ierr ) 
        call mpi_barrier(mpi_comm_comp,ierr)

!GWV   reshape the gathered array and collect in all procs
        isum=1
        do n=0,num_procs-1
          do jj=jsxa(n),jexa(n)
            do ii=isxa(n),iexa(n)
              if(isum .gt. im*jm .or. ii .gt. im .or. ii .lt. 1 .or. jj .gt. jm .or. jj .lt. 1) &
              write(*,901)' BOUNDS FAIL in reshape',isum,ii,jj,im*jm,im,im*jm
              b(ii,jj)=buff(isum)
              isum=isum+1
            end do
          end do
        end do

      endif   ! num_procs <= 1

 901  format(a30,10i10)

      deallocate(buff)
      deallocate(rbufs)

      end

