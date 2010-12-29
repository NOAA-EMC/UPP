c Stub versions of MPI F77 routines (single processor) - most do nothing
c 20101228 slovacek add mpi_file_open and mpi_file_read_at declarations
c          these files will do nothing, but report that serial bbuilds
c          don't support the functionalty and return a failure code

c timer routine
c can replace etime with standard UNIX call on a particular system

      double precision function mpi_wtime()
      real array(2)
      
      mpi_wtime = etime(array)
      
      return
      end


      subroutine mpi_init(ierror)

      return
      end

      subroutine mpi_initialized(mpi_inited, ierror)
      logical mpi_inited

      return
      end



      subroutine mpi_finalize(ierror)

      return
      end


      subroutine mpi_abort(mpi_comm,ierror)

      stop
      end

c return me = 0

      subroutine mpi_comm_rank(mpi_comm,me,ierror)

      me = 0

      return
      end

c return nprocs = 1

      subroutine mpi_comm_size(mpi_comm,nprocs,ierror)

      nprocs = 1

      return
      end


      subroutine mpi_barrier(mpi_comm,ierror)

      return
      end

c warn against sending message to self, since no data copy is done

      subroutine mpi_send(data,n,mpi_datatype,iproc,itag,
     $     mpi_comm,ierror)
      
      write (6,*)
     $     'MPI Stub WARNING: Should not send message to self'

      return
      end

        subroutine mpi_comm_dup(oldcomm, newcomm,ierror)
        integer oldcomm, newcomm,ierror
        newcomm = oldcomm
        return
        end

        subroutine mpi_isend(buf,count,datatype,source,
     &  tag,comm,request,ierror)
        integer buf(*), count,datatype,source,tag,comm,
     &  request,ierror
        call mpi_error()
        return
        end

        subroutine mpi_error()
        print *, 'mpi_error called'
        stop
        end

        subroutine mpi_comm_group (com,group,ierr )
        integer com,group,ierr
        write(6,*) 'dont want to see this'
        group=com
        return
        end

        subroutine mpi_group_excl(intin,n,ranks,ngroup,ierr)
        integer intin,n,ranks(n),ngroup,ierr
        ngroup=intin
        return
        end

        subroutine mpi_group_free (intin,ierr )
        integer intin,ierr
        return
        end

        subroutine mpi_intercomm_create (lcom,llead,ipeer,irem,itag,
     &                  newcom,ierr )
        integer lcom,llead,ipeer,irem,itag,newcom,ierr
        newcom=lcom
        write(6,*) 'shouldnt be calling this!'
        return
        end

        subroutine mpi_comm_create (com,group,ncom,ierr )
        integer com,group,ncom,ierr
        ncom=com
        return
        end

        subroutine mpi_sendrecv(sendbuf,sendcount,sendtype,
     +  dest,sendtag,recvbuf,recvcount,recvtype,source,recvtag,
     +  comm,status,ierr)

        integer sendcount,sendtype
        integer dest,sendtag,recvcount
        integer recvtype,source,recvtag,comm,ierr,status(*)

        write(6,*) 'should never execute this when running with one CPU'

        return
        end

      subroutine mpi_gatherv(sendbuf,sendcount,sendtype,
     +  recvbuf,recvcounts,displs,recvtype,root,comm,ierr)

        integer sendcount,sendtype,recvcounts(*),displs(*)
        integer recvtype,root,comm,ierr
        write(6,*) 'should never execute this when running with one CPU'

      return
      end


c warn against sending message to self, since no data copy is done

      subroutine mpi_rsend(data,n,mpi_datatype,iproc,itag,
     $     mpi_comm,ierror)

      write (6,*)
     $     'MPI Stub WARNING: Should not send message to self'

      return
      end

c warn against receiving message from self, since no data copy is done

      subroutine mpi_probe(iproc,itag,mpi_comm,istatus,ierror)

      return
      end


      subroutine mpi_recv(data,n,mpi_datatype,iproc,itag,
     $     mpi_comm,istatus,ierror)

      write (6,*)
     $     'MPI Stub WARNING: Should not recv message from self'

      return
      end

c warn against querying message from self, since no data copy is done

      subroutine mpi_get_count(istatus,mpi_datatype,icount,ierror)

      write (6,*)
     $     'MPI Stub WARNING: Should not query message from self'

      return
      end


c warn against receiving message from self, since no data copy is done

      subroutine mpi_irecv(data,n,mpi_datatype,iproc,itag,
     $     mpi_comm,irequest,ierror)

      write (6,*)
     $     'MPI Stub WARNING: Should not recv message from self'

      return
      end

c warn against waiting on message from self, since no data copy is done

      subroutine mpi_wait(irequest,istatus,ierror)

      write (6,*)
     $     'MPI Stub WARNING: Should not wait on message from self'

      return
      end

c warn against waiting on message from self, since no data copy is done

      subroutine mpi_waitall(icount,irequest,istatus,ierror)

      write (6,*)
     $     'MPI Stub WARNING: Should not wait on message from self'

      return
      end

c warn against waiting on message from self, since no data copy is done

      subroutine mpi_waitany(icount,array_of_requests,
     $     index,istatus,ierror)

      write (6,*)
     $     'MPI Stub WARNING: Should not wait on message from self'

      return
      end


      subroutine mpi_bcast(data,n,mpi_datatype,node,mpi_comm,ierror)

      return
      end

      subroutine mpi_type_size(mpi_datatype,isize,ierror)
     
      return
      end

c copy values from data1 to data2

      subroutine mpi_reduce(data1,data2,n,mpi_datatype,
     $     mpi_operation,io_task,mpi_comm,ierror)
      include "mpif.h"

      if (mpi_datatype.eq.mpi_integer) then
        call mpi_copy_integer(data1,data2,n)
      else if (mpi_datatype.eq.mpi_real) then
        call mpi_copy_real(data1,data2,n)
      else if (mpi_datatype.eq.mpi_double_precision) then
        call mpi_copy_double_precision(data1,data2,n)
      else if (mpi_datatype.eq.mpi_real8) then
        call mpi_copy_double_precision(data1,data2,n)
      endif

      return
      end


      subroutine mpi_allreduce(data1,data2,n,mpi_datatype,
     $     mpi_operation,mpi_comm,ierror)
      include "mpif.h"

      print *,'allreduce ',mpi_datatype,mpi_integer,mpi_real
      if (mpi_datatype.eq.mpi_integer) then
        call mpi_copy_integer(data1,data2,n)
      else if (mpi_datatype.eq.mpi_real) then
        call mpi_copy_real(data1,data2,n)
      else if (mpi_datatype.eq.mpi_double_precision) then
        call mpi_copy_double_precision(data1,data2,n)
      else if (mpi_datatype.eq.mpi_real8) then
        call mpi_copy_double_precision(data1,data2,n)
      endif

      return
      end

      subroutine mpi_gather(data1,nsend,mpi_sendtype,data2,
     $     nrecv,mpi_recvtype,io_task,mpi_comm,ierror)
      include "mpif.h"

      if (mpi_sendtype.eq.mpi_integer) then
        call mpi_copy_integer(data1,data2,nsend)
      else if (mpi_sendtype.eq.mpi_real) then
        call mpi_copy_real(data1,data2,nsend)
      else if (mpi_sendtype.eq.mpi_double_precision) then
        call mpi_copy_double_precision(data1,data2,nsend)
      else if (mpi_sendtype.eq.mpi_real8) then
        call mpi_copy_double_precision(data1,data2,nsend)
      endif

      return
      end


c copy values from data1 to data2

      subroutine mpi_allgather(data1,nsend,mpi_sendtype,
     $     data2,nrecv,mpi_recvtype,mpi_comm,ierror)
      include "mpif.h"

      if (mpi_sendtype.eq.mpi_integer) then
        call mpi_copy_integer(data1,data2,nsend)
      else if (mpi_sendtype.eq.mpi_real) then
        call mpi_copy_real(data1,data2,nsend)
      else if (mpi_sendtype.eq.mpi_double_precision) then
        call mpi_copy_double_precision(data1,data2,nsend)
      else if (mpi_sendtype.eq.mpi_real8) then
        call mpi_copy_double_precision(data1,data2,nsend)
      endif

      return
      end


c copy values from data1 to data2

      subroutine mpi_allgatherv(data1,nsend,mpi_sendtype,
     $     data2,nrecv,ndispls,mpi_recvtype,mpi_comm,ierror)
      include "mpif.h"

      if (mpi_sendtype.eq.mpi_integer) then
        call mpi_copy_integer(data1,data2,nsend)
      else if (mpi_sendtype.eq.mpi_real) then
        call mpi_copy_real(data1,data2,nsend)
      else if (mpi_sendtype.eq.mpi_double_precision) then
        call mpi_copy_double_precision(data1,data2,nsend)
      else if (mpi_sendtype.eq.mpi_real8) then
        call mpi_copy_double_precision(data1,data2,nsend)
      endif

      return
      end


c copy values from data1 to data2

      subroutine mpi_reduce_scatter(data1,data2,n,mpi_datatype,
     $     mpi_operation,mpi_comm,ierror)
      include "mpif.h"

      if (mpi_datatype.eq.mpi_integer) then
        call mpi_copy_integer(data1,data2,n)
      else if (mpi_datatype.eq.mpi_real) then
        call mpi_copy_real(data1,data2,n)
      else if (mpi_datatype.eq.mpi_double_precision) then
        call mpi_copy_double_precision(data1,data2,n)
      else if (mpi_datatype.eq.mpi_real8) then
        call mpi_copy_double_precision(data1,data2,n)
      endif

      return
      end


      subroutine mpi_cart_create(mpi_comm,ndims,dims,
     $     periods,reorder,mpi_comm_cart,ierror)
      logical periods(*),reorder
      integer dims(*)

      return
      end


c set all coords = 0

      subroutine mpi_cart_get(mpi_comm,ndims,dims,periods,
     $     coords,ierror)
      logical periods(*)
      integer dims(*),coords(*)

      do i = 1,ndims
        coords(i) = 0
      enddo

      return
      end


c set isource = idest = self = 0

      subroutine mpi_cart_shift(mpi_comm,idir,idisp,
     $     isource,idest,ierror)

      isource = 0
      idest = 0

      return
      end


      subroutine mpi_comm_split(mpi_comm,icolor,ikey,new_comm,ierror)

      return
      end


      subroutine mpi_comm_free(mpi_comm,ierror)

      return
      end


c -------------------
c Added routines for data copying
c -------------------

      subroutine mpi_copy_integer(data1,data2,n)
      integer data1(*),data2(*)

      do i = 1,n
        data2(i) = data1(i)
      enddo

      return
      end


      subroutine mpi_copy_real(data1,data2,n)
      real data1(*),data2(*)

      do i = 1,n
        data2(i) = data1(i)
      enddo

      return
      end


      subroutine mpi_copy_double_precision(data1,data2,n)
      double precision data1(*),data2(*)

      do i = 1,n
        data2(i) = data1(i)
      enddo

      return
      end

      subroutine mpi_cart_coords(mpi_comm, rank, maxdims, coords, ierr)
      integer coords(*)

      return
      end

      subroutine mpi_scatterv(globbuf, counts, displs, gettype, 
     $     field, my_count, sendtype, root, comm, ierr)
      return
      end

      subroutine mpi_group_incl(group,n,ranks,newgroup,ierror)
      return
      end

c-----------------------------
c Add file IO routine which are currently unsuported
c-----------------------------
      subroutine mpi_file_open(comm, filename, amode, info, fh, ierr)
        integer ierr

        print *, "MPI STUB - mpi_file_open not supported"
        ierr = 255

        return
      end 
  
      subroutine mpi_file_read_at(fh, offset, buf, count, dtype, status,
     $      ierr)
        integer ierr

        print *, "MPI STUB - mpi_file_read_at not supported"
        ierr = 255

        return
      end

      subroutine mpi_file_close(fh, ierr)
      return
      end
