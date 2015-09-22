#ifndef __MY_MPI_HXX__
#define __MY_MPI_HXX__

#if HAVE_MPI
#include<mpi.h>


template <typename T>
struct MPI_Type
{
  static inline MPI_Datatype get(){return MPI_DATATYPE_NULL;}
};

#define ASSOCIATE_MPI_TYPE(T,M)			\
  template <>					\
  struct MPI_Type<T> {				\
    static inline MPI_Datatype get(){return M;}	\
  };
  
  ASSOCIATE_MPI_TYPE(char, MPI_CHAR)
  ASSOCIATE_MPI_TYPE(unsigned char,MPI_UNSIGNED_CHAR)
  ASSOCIATE_MPI_TYPE(short,MPI_SHORT)
  ASSOCIATE_MPI_TYPE(unsigned short,MPI_UNSIGNED_SHORT)
  ASSOCIATE_MPI_TYPE(int,MPI_INT)
  ASSOCIATE_MPI_TYPE(unsigned int,MPI_UNSIGNED)
  ASSOCIATE_MPI_TYPE(long,MPI_LONG)
  ASSOCIATE_MPI_TYPE(unsigned long,MPI_UNSIGNED_LONG)
  ASSOCIATE_MPI_TYPE(float,MPI_FLOAT)
  ASSOCIATE_MPI_TYPE(double,MPI_DOUBLE)
  ASSOCIATE_MPI_TYPE(long double,MPI_LONG_DOUBLE)
  
#else

  //#define 

  template <typename T>
  struct MPI_Type<T>
  {
    static inline MPI_Datatype get(){return 0;}
  }

#endif
#endif

