#ifndef DUMMY_CELL_21_HXX
#define DUMMY_CELL_21_HXX

#include <vector>
#include "myIO.hxx"

template <class dataType>
struct dummyCell21
{
  typedef dataType dataT;

  double R;
  double U;
  std::vector<dataT> val;
  
  dummyCell21()
  {
   
  }
 
  template <class inputIterator>
  dummyCell21(double R_, double U_,const inputIterator &start,const inputIterator &stop):
    R(R_),U(U_)
  {
    val.assign(start,stop);
  }

  dummyCell21(double R_, double U_, int ct):
    R(R_),U(U_)
  {
    val.assign(ct,0);
  }
      
  ~dummyCell21() 
  {
   
  }

  void read(FILE *f, bool swap)
  {
    int nv;
    myIO::fread(&R,sizeof(double),1,f,swap);
    myIO::fread(&U,sizeof(double),1,f,swap);
    myIO::fread(&nv,sizeof(int),1,f,swap);  
    val.resize(nv);
    myIO::fread(&val[0],sizeof(dataT),nv,f,swap);  
  }

  void write(FILE *f)
  {
    fwrite(&R,sizeof(double),1,f);
    fwrite(&U,sizeof(double),1,f);
    int nv=val.size();
    fwrite(&nv,sizeof(int),1,f);
    fwrite(&val[0],sizeof(dataT),nv,f);
  }

  void print()
  {
    printf("(R=%e U=%e V=%e)",R,U,val[0]);
  }

  dataT operator[](int id) const
  {
    return val[id];
  }
};

#endif
