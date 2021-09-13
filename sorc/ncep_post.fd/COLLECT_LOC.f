!> @file
!
!> SUBPROGRAM:    COLLECT     GATHERS FROM ALL MPI TASKS
!!   PRGRMMR: TUCCILLO        ORG: IBM
!!
!! ABSTRACT:
!!     GATHER "A" FROM ALL MPI TASKS ONTO TASK 0
!!
!! PROGRAM HISTORY LOG:
!!   00-01-06  TUCCILLO - ORIGINAL
!!
!! USAGE:    CALL COLLECT(A)
!!   INPUT ARGUMENT LIST:
!!     A        - ARRAY BEING GATHERED
!!
!!   OUTPUT ARGUMENT LIST:
!!     A        - GATHERED ARRAY - ONLY VALID ON TASK 0
!!
!!   OUTPUT FILES:
!!     STDOUT  - RUN TIME STANDARD OUT.
!!
!!   SUBPROGRAMS CALLED:
!!       MPI_GATHERV
!!     UTILITIES:
!!       NONE
!!     LIBRARY:
!!       COMMON   - CTLBLK.comm
!!
!!   ATTRIBUTES:
!!     LANGUAGE: FORTRAN
!!     MACHINE : IBM RS/6000 SP
!!
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
      write(0,*) ' GWVX COLL CALL'
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
               write(0,901)' GWVX BOUNDS2 FAIL in reshape ',isum,ii,jj,im*jm,im,im*jm
            rbufs(isum)=a(ii,jj)
            isum=isum+1
           end do
           end do
!GWV  end reshape

!UNCOMMENT    POST TEST      call mpi_gatherv(rbufs,icnt(me),MPI_REAL, buff,icnt,idsp,MPI_REAL,0,MPI_COMM_COMP, ierr )
        call mpi_gatherv(rbufs,icnt(me),MPI_REAL,  buff,icnt,idsp,MPI_REAL,0,MPI_COMM_WORLD, ierr )  !GWVX COMMENT

!GWV   reshape the gathered array 
             if(me .eq. 0) then
             isum=1 
             do n=0,num_procs-1
            do jj=jsxa(n),jexa(n)
            do ii=isxa(n),iexa(n)
            if(isum .gt. im*jm .or. ii .gt. im .or. ii .lt. 1 .or. jj .gt. jm .or. jj .lt. 1) &
               write(0,901)' GWVX BOUNDS FAIL in reshape ',isum,ii,jj,im*jm,im,im*jm
 901    format(a30,10i10)
            b(ii,jj)=buff(isum)
            isum=isum+1
            end do
            end do
            end do
            

      end if
           endif
             deallocate(buff)

      end               
