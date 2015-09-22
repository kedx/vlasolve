#ifndef __MPI_COM_HXX__
#define __MPI_COM_HXX__

#if HAVE_MPI
#include<mpi.h>
#endif

class mpiCommunication {
private:
  int myRank;
  int nProcs;
public:

#if HAVE_MPI

  mpiCommunication(int argc, char **argv)
  {
    MPI_Init(&argc,&argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    //printf("Hello from %d / %d.\n",myRank,nProcs);
  }

  ~mpiCommunication()
  {
    MPI_Finalize();
  }

#else

  mpiCommunication(int argc, char **argv)
  {
    myRank=0;
    nProcs=1;
  }

  ~mpiCommunication()
  {

  }

#endif

  int rank() {return myRank;}
  int size() {return nProcs;}
};


#endif
