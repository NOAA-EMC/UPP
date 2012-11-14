/* Stub versions of MPI C routines (single processor) - most do nothing */

#include "mpi.h"
#include "stdio.h"

/* function prototypes */

void mpi_copy_integer(int *, int *, int);
void mpi_copy_float(float *, float *, int);
void mpi_copy_double(double *, double *, int);
void mpi_copy_byte(char *, char *, int);

/* MPI Functions */

void MPI_Comm_rank(MPI_Comm comm, int *me)
{
  *me = 0;
}

void MPI_Comm_size(MPI_Comm comm, int *nprocs)
{
  *nprocs = 1;
}

void MPI_Send(void *buf, int count, MPI_Datatype datatype,
	      int dest, int tag, MPI_Comm comm)
{
  printf("MPI Stub WARNING: Should not send message to self\n");
}

void MPI_Recv(void *buf, int count, MPI_Datatype datatype,
	      int source, int tag, MPI_Comm comm, MPI_Status *status)
{
  printf("MPI Stub WARNING: Should not recv message from self\n");
}

void MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
	       int source, int tag, MPI_Comm comm, MPI_Request *request)
{
  printf("MPI Stub WARNING: Should not recv message from self\n");
}

void MPI_Isend(void *buf, int count, MPI_Datatype datatype,
               int source, int tag, MPI_Comm comm, MPI_Request *request)
{
  printf("MPI Stub WARNING: Should not send message from self\n");
}


void MPI_Wait(MPI_Request *request, MPI_Status *status)
{
  printf("MPI Stub WARNING: Should not wait on message from self\n");
}

void MPI_Abort(MPI_Comm comm, MPI_Status *status)
{
  printf("MPI Stub : from Abort\n");
}

void  MPI_Finalize( void )
{
 printf("MPI Stub WARNING: Should not MPI_Finalize\n");
}

void  MPI_Initialized( int * tag )
{
 printf("MPI Stub WARNING: Should not MPI_Initialize\n");
}


void MPI_Waitany(int count, MPI_Request *request, int *index, 
		 MPI_Status *status)
{
  printf("MPI Stub WARNING: Should not wait on message from self\n");
}

void MPI_Comm_dup(MPI_Comm comm, MPI_Comm *comm_out) { }

void MPI_Comm_free(MPI_Comm *comm) { }

int MPI_Bcast(void *sendbuf, int sendcount, MPI_Datatype datatype,
              int task, MPI_Comm comm)
{
  return 0;
}

/* copy values from data1 to data2 */

void MPI_Allreduce(void *sendbuf, void *recvbuf, int sendcount,
		   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  if (datatype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,sendcount);
  else if (datatype == MPI_FLOAT)
    mpi_copy_float(sendbuf,recvbuf, sendcount);
  else if (datatype == MPI_DOUBLE)
    mpi_copy_double(sendbuf,recvbuf, sendcount);
  else if (datatype == MPI_BYTE)
    mpi_copy_byte(sendbuf,recvbuf, sendcount);
}

/* copy values from data1 to data2 */

void MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, int recvcount, MPI_Datatype recvtype,
                int task, MPI_Comm comm)
{
  if (sendtype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_FLOAT)
    mpi_copy_float(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_DOUBLE)
    mpi_copy_double(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_BYTE)
    mpi_copy_byte(sendbuf,recvbuf,sendcount);
}

/* copy values from data1 to data2 */

void MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                 void *recvbuf, int *recvcounts, MPI_Datatype recvtype,
                 int *task, MPI_Comm comm)
{
  if (sendtype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_FLOAT)
    mpi_copy_float(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_DOUBLE)
    mpi_copy_double(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_BYTE)
    mpi_copy_byte(sendbuf,recvbuf,sendcount);
}

/* copy values from data1 to data2 */

void MPI_Scatterv(void *sendbuf, int sendcount, int *displs, 
                  MPI_Datatype sendtype, void *recvbuf, 
                  int *recvcounts, MPI_Datatype recvtype,
                  int *task, MPI_Comm comm)
{
  if (sendtype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_FLOAT)
    mpi_copy_float(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_DOUBLE)
    mpi_copy_double(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_BYTE)
    mpi_copy_byte(sendbuf,recvbuf,sendcount);
}

/* copy values from data1 to data2 */

void MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		   void *recvbuf, int recvcount, MPI_Datatype recvtype,
		   MPI_Comm comm)
{
  if (sendtype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_FLOAT)
    mpi_copy_float(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_DOUBLE)
    mpi_copy_double(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_BYTE)
    mpi_copy_byte(sendbuf,recvbuf,sendcount);
}

/* copy values from data1 to data2 */

void MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		    void *recvbuf, int *recvcounts, int *displs,
		    MPI_Datatype recvtype, MPI_Comm comm)
{
  if (sendtype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_FLOAT)
    mpi_copy_float(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_DOUBLE)
    mpi_copy_double(sendbuf,recvbuf,sendcount);
  else if (sendtype == MPI_BYTE)
    mpi_copy_byte(sendbuf,recvbuf,sendcount);
}

/* copy values from data1 to data2 */

void MPI_Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts,
			MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  if (datatype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,*recvcounts);
  else if (datatype == MPI_FLOAT)
    mpi_copy_float(sendbuf,recvbuf,*recvcounts);
  else if (datatype == MPI_DOUBLE)
    mpi_copy_double(sendbuf,recvbuf,*recvcounts);
  else if (datatype == MPI_BYTE)
    mpi_copy_byte(sendbuf,recvbuf,*recvcounts);
}


/*
-------------------
Added routines for data copying
-------------------
*/

void mpi_copy_integer(int *data1, int *data2, int n)
{
  int i;
  for (i = 0; i < n; i++) data2[i] = data1[i];
}

void mpi_copy_float(float *data1, float *data2, int n)
{
  int i;
  for (i = 0; i < n; i++) data2[i] = data1[i];
}

void mpi_copy_double(double *data1, double *data2, int n)
{
  int i;
  for (i = 0; i < n; i++) data2[i] = data1[i];
}

void mpi_copy_byte(char *data1, char *data2, int n)
{
  int i;
  for (i = 0; i < n; i++) data2[i] = data1[i];
}

void MPI_Cart_shift(MPI_Comm comm,int direction,int disp,
               int *rank_source,int *rank_dest)
{

}


int MPI_Cart_rank ( MPI_Comm comm, int *coords, int *rank )
{
  return 0;
}

int MPI_Alltoall ( void *sendbuf, int *sendcnts, MPI_Datatype sendtype,
        void *recvbuf, int *recvcnts, MPI_Datatype recvtype, MPI_Comm comm ) 
{
  if (sendtype == MPI_INT)
    mpi_copy_integer(sendbuf,recvbuf,sendcnts[0]);
  else if (sendtype == MPI_FLOAT)
    mpi_copy_float((float *)sendbuf,(float *)recvbuf, sendcnts[0]);
  else if (sendtype == MPI_DOUBLE)
    mpi_copy_double((double *)sendbuf,(double *)recvbuf, sendcnts[0]);
  else if (sendtype == MPI_BYTE)
    mpi_copy_byte((char *)sendbuf,(char *)recvbuf, sendcnts[0]);
  return 0;
}

int MPI_Alltoallv ( void *sendbuf, int *sendcnts, int *sdispls, MPI_Datatype sendtype,
        void *recvbuf, int *recvcnts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm ) 
{
  int *recv_offset;

  recv_offset = recvbuf + *sdispls;

  if (sendtype == MPI_INT)
    mpi_copy_integer((int *)sendbuf,(int *)recv_offset,sendcnts[0]);
  else if (sendtype == MPI_FLOAT)
    mpi_copy_float((float *)sendbuf,(float *)recv_offset, sendcnts[0]);
  else if (sendtype == MPI_DOUBLE)
    mpi_copy_double((double *)sendbuf,(double *)recv_offset, sendcnts[0]);
  else if (sendtype == MPI_BYTE)
    mpi_copy_byte((char *)sendbuf,(char *)recv_offset, sendcnts[0]);
  return 0;
}

int MPI_Cart_coords ( MPI_Comm comm, int rank, int maxdims, int *coords )
{
  return 0;
}
