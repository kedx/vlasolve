#ifndef __MY_MPI_HXX__
#define __MY_MPI_HXX__

#if HAVE_MPI
#include<mpi.h>

#ifndef MPID_TAG_UB
#define MPID_TAG_UB (1<<((sizeof(int)*8) -2)) // maximum tag value
#endif

template <typename T>
struct MPI_Type
{
  static inline MPI_Datatype get(){return MPI_DATATYPE_NULL;}
};

#define ASSOCIATE_MPI_TYPE(T,M)			\
  template <>					\
  struct MPI_Type<T> {				\
    static inline MPI_Datatype get(){return M;}	\
  };						\

 
ASSOCIATE_MPI_TYPE(char              ,MPI_CHAR)
ASSOCIATE_MPI_TYPE(unsigned char     ,MPI_UNSIGNED_CHAR)
ASSOCIATE_MPI_TYPE(short             ,MPI_SHORT)
ASSOCIATE_MPI_TYPE(unsigned short    ,MPI_UNSIGNED_SHORT)
ASSOCIATE_MPI_TYPE(int               ,MPI_INT)
ASSOCIATE_MPI_TYPE(unsigned int      ,MPI_UNSIGNED)
ASSOCIATE_MPI_TYPE(long              ,MPI_LONG)
ASSOCIATE_MPI_TYPE(unsigned long     ,MPI_UNSIGNED_LONG)
ASSOCIATE_MPI_TYPE(long long         ,MPI_LONG_LONG)
ASSOCIATE_MPI_TYPE(unsigned long long,MPI_UNSIGNED_LONG_LONG)
ASSOCIATE_MPI_TYPE(float             ,MPI_FLOAT)
ASSOCIATE_MPI_TYPE(double            ,MPI_DOUBLE)
ASSOCIATE_MPI_TYPE(long double       ,MPI_LONG_DOUBLE)

/*
ASSOCIATE_MPI_TYPE(int8_t  ,MPI_INT8_T)
ASSOCIATE_MPI_TYPE(int16_t ,MPI_INT16_T)
ASSOCIATE_MPI_TYPE(int32_t ,MPI_INT32_T)
ASSOCIATE_MPI_TYPE(int64_t ,MPI_INT64_T)

ASSOCIATE_MPI_TYPE(uint8_t  ,MPI_UINT8_T)
ASSOCIATE_MPI_TYPE(uint16_t ,MPI_UINT16_T)
ASSOCIATE_MPI_TYPE(uint32_t ,MPI_UINT32_T)
ASSOCIATE_MPI_TYPE(uint64_t ,MPI_UINT64_T)
*/

#else // we don't have MPI

#include <time.h>

#ifndef MPID_TAG_UB
#define MPID_TAG_UB (1<<((sizeof(int)*8) -2))
#endif
#define MPI_ANY_TAG 0

#ifndef MPI_ORDER_FORTRAN
#define MPI_ORDER_FORTRAN 1
#endif
#ifndef MPI_ORDER_C
#define MPI_ORDER_C 2
#endif

typedef void* MPI_Comm;
typedef int MPI_Datatype;
typedef int MPI_Request;
typedef void* MPI_Op;
typedef int MPI_Aint;
struct MPI_Status {static const int MPI_SOURCE=0;};

MPI_Op MPI_MAX = NULL;
MPI_Op MPI_MIN = NULL;
MPI_Op MPI_SUM = NULL;
MPI_Datatype MPI_BYTE = 1;
MPI_Datatype MPI_DATATYPE_NULL = 0;


template <typename T>
struct MPI_Type
{
  static inline MPI_Datatype get(){return MPI_DATATYPE_NULL;}
};

  
int MPI_Type_free(MPI_Datatype *datatype) {return 0;}
void MPI_Type_create_struct(long,int*,MPI_Aint*,MPI_Datatype*,MPI_Datatype*) {}
void MPI_Type_commit(MPI_Datatype*) {}
int MPI_Type_create_subarray(int ndims,
			     int array_of_sizes[],
			     int array_of_subsizes[],
			     int array_of_starts[],
			     int order,
			     MPI_Datatype oldtype,
			     MPI_Datatype *newtype) {return 0;}
int MPI_Type_contiguous(int count,
			MPI_Datatype old_type,
			MPI_Datatype *new_type_p) {return 0;}


#endif

#endif
