#ifndef DUMMY_CELL_21_HXX
#define DUMMY_CELL_21_HXX

#include "myIO.hxx"

// for delayed boundary conditions
template <class dataType>
struct dummyCell21
{
  typedef dataType dataT;

  double R;
  double U;
  dataT *val;
  int nVal;

  dummyCell21(const dummyCell21 &other):
    R(other.R),U(other.U),nVal(other.nVal)
  {
    if (nVal>0)
      {
	val=(dataT*)malloc(sizeof(dataT)*nVal);
	std::copy(other.val,other.val+nVal,val);
      }
    else
      {
	val=NULL;
	nVal=-1;
      }
  }

  dummyCell21& operator=(const dummyCell21& other)
  {
    R=other.R;
    U=other.U;

    if (other.nVal!=nVal)
      {
	nVal=other.nVal;
	if (nVal>0) val=(dataT*)realloc(val,sizeof(dataT)*nVal);
	else 
	  {
	    free(val); 
	    val=NULL;
	  }
      }
     
    if (nVal>0) std::copy(other.val,other.val+nVal,val);
  }

  dummyCell21()
  {
    nVal=-1;
    val=NULL;
  }
 
  template <class inputIterator>
  dummyCell21(double R_, double U_,const inputIterator &start,const inputIterator &stop):
    nVal(std::distance(start,stop)),R(R_),U(U_)
  {
    val=(dataT*)malloc(sizeof(dataT)*nVal);
    std::copy(start,stop,val);
  }

  dummyCell21(double R_, double U_, int ct):
    nVal(ct),R(R_),U(U_)
  {
    val=(dataT*)calloc(nVal,sizeof(dataT));
  }
      
  ~dummyCell21() 
  {
    free(val);
  }

  void read(FILE *f, bool swap)
  {
    int nv=nVal;
    myIO::fread(&R,sizeof(double),1,f,swap);
    myIO::fread(&U,sizeof(double),1,f,swap);
    myIO::fread(&nVal,sizeof(int),1,f,swap);
    if (nv!=nVal) val=(dataT*)realloc(val,sizeof(dataT)*nVal);
    myIO::fread(val,sizeof(dataT),nVal,f,swap);
  }

  void write(FILE *f)
  {
    fwrite(&R,sizeof(double),1,f);
    fwrite(&U,sizeof(double),1,f);
    fwrite(&nVal,sizeof(int),1,f);
    fwrite(val,sizeof(dataT),nVal,f);
  }

  void print()
  {
    printf("(R=%e U=%e V=%e)",R,U,val[0]);
  }
};

#endif
