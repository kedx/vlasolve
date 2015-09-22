#ifndef __MPI_COM_HXX__
#define __MPI_COM_HXX__

#include "myMPI.hxx"
#include <sys/types.h>
#include <unistd.h>

class mpiCommunication {
private:
  int myRank;
  int nProcs;
  MPI_Comm com;
  bool finalize;
  bool deleteCom;

#if HAVE_MPI

private:
  void init(MPI_Comm myCom)
  {
    com=myCom;
    MPI_Comm_rank(com, &myRank);
    MPI_Comm_size(com, &nProcs);
  }
  
  mpiCommunication(MPI_Comm myCom):deleteCom(true),finalize(false)
  {
    init(myCom);
  }

public:

  mpiCommunication(int argc, char **argv):deleteCom(false),finalize(true)
  {
    MPI_Init(&argc,&argv);
    init(MPI_COMM_WORLD);
    //printf("Hello from %d / %d.\n",myRank,nProcs);
  }

  ~mpiCommunication()
  {
    if (deleteCom) MPI_Comm_free(&com);
    if (finalize) MPI_Finalize();   
  }

  mpiCommunication split(int color, int key)
  {
    MPI_Comm myCom;
    MPI_Comm_split(com,color,key,&myCom);
    return mpiCommunication(myCom);
  }

  void debug()
  {
    int wait = 1;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("Process %d on %s has PID %d.\n",rank(),hostname,getpid());
    fflush(stdout);
    while (wait)
      {
	sleep(10);
	MPI_Allreduce(MPI_IN_PLACE,&wait,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
      }
  }
  
  template <typename T>
  T max(const T &val)
  {
    T tmp=val;
    MPI_Allreduce(MPI_IN_PLACE, &tmp, 1, MPI_Type<T>::get(), MPI_MAX, com);
    return tmp;
  }

  void barrier()
  {
    MPI_Barrier(com);
  }

  template <class container>
  int Bcast(container &buffer, int root,long count=-1)
  {
    if (count<0) count=buffer.size();
    return MPI_Bcast(&buffer[0],count,MPI_Type<typename container::value_type>::get(),root,com);
  }

  template <typename T>
  int Allgather_inplace(T* buffer, long count)
  {
    return MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
			 buffer,count,MPI_Type<T>::get(),com);
		  
  }

  template <typename T>
  int Allreduce_inplace(T* buffer, long count, MPI_Op op)
  {
    return MPI_Allreduce(MPI_IN_PLACE,buffer,count,MPI_Type<T>::get(),op,com);
		  
  }

  template <class container>
  int Allgather_inplace(container &buffer)
  {
    return MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,
			 &buffer[0],buffer.size()/size(),
			 MPI_Type<typename container::value_type>::get(),com);
		  
  }

  template <class container>
  int Allreduce_inplace(container &buffer, MPI_Op op)
  {
    return MPI_Allreduce(MPI_IN_PLACE,&buffer[0],buffer.size(),
			 MPI_Type<typename container::value_type>::get(), op, com);		  
  }

  

#else

public:
  mpiCommunication(int argc, char **argv)
  {
    myRank=0;
    nProcs=1;
  }

  ~mpiCommunication()
  {

  }

  mpiCommunication split(int color, int key)
  {
    return mpiCommunication(0,NULL);
  }

  void debug() {}

  double max(double &val) {return val;}

  void barrier() {}
 
#endif

  int rank() {return myRank;}
  int size() {return nProcs;}
};


#endif
