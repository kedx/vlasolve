#ifndef __GRID_NN_HXX__
#define __GRID_NN_HXX__

#include <utility>
#include <string>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "sliceNNIterator.hxx"
#include "grid_base.hxx"
#include "defines.h"

template <
  int PU_DIMS,
  typename dataType,
  class valLocationTraits
  >
class gridNN : public gridBase<dimTraits<PU_DIMS,PU_DIMS>,dataType,valLocationTraits> {
public:

  static const std::string getTag() {char tmp[255];sprintf(tmp,"%s%d%s","GRID_",PU_DIMS," v0.10");return std::string(tmp);}

  typedef dataType value_type;
  typedef gridNN<PU_DIMS,value_type,valLocationTraits> myT;  
  typedef gridBase<dimTraits<PU_DIMS,PU_DIMS>,value_type,valLocationTraits> bT;

  typedef typename bT::iterator iterator;
  friend class sliceNNIterator<myT>;
  typedef sliceNNIterator<myT> sliceIterator;
  
  typedef typename bT::scaleT scaleT;
  typedef typename bT::scaleTypeT scaleTypeT;
  typedef typename bT::scaleTypeV scaleTypeV; 
  typedef typename bT::valLocationT valLocationT;
  typedef typename bT::valLocationV valLocationV;
  typedef typename bT::valLocationTr valLocationTr;
  typedef typename bT::paramT paramT;
  
  typedef typename bT::dirT dirT;

  static const int P_DIMS = bT::P_DIMS;
  static const int U_DIMS = bT::U_DIMS;
  static const int J_DIMS = bT::J_DIMS;
  static const int DIMS = bT::DIMS;

protected:
  long sliceCount[2];
  long sliceCountFull[2];

  long sliceStride[2][P_DIMS+1];
  long sliceStrideFull[2][U_DIMS+1];

public:
  gridNN():bT()
  {}

   ~gridNN()
  {
    
  }
  
  virtual void write(FILE *f)
  {
    bT::write(f);
    myIO::writeTag(f,getTag());
    std::vector<char> spare(256*8,0);
    fwrite(&spare[0],sizeof(char),spare.size(),f);
  }

  virtual void read(FILE *f, bool swap)
  {
    bT::read(f,swap);
    myIO::checkTag(f,getTag()); 
    std::vector<char> spare(256*8,0);
    myIO::fread(&spare[0],sizeof(char),spare.size(),f,swap);
    init();
  }

  virtual void init(const paramT &p_)
  {
    bT::init(p_);
    init();
  }

private:
  void init()
  {
    int i,t;
    sliceStride[0][0]=sliceStride[1][0]=1;
    sliceStrideFull[0][0]=sliceStrideFull[1][0]=1;
   
    for (i=0;i<P_DIMS;i++) 
      sliceStride[1][i+1]=sliceStride[1][i]*
	(bT::getValCoord_P(i).size()-bT::p.lowPmargin[i]-bT::p.highPmargin[i]);
    for (i=0;i<U_DIMS;i++) 
      sliceStride[0][i+1]=sliceStride[0][i]*
	(bT::getValCoord_U(i).size()-bT::p.lowUmargin[i]-bT::p.highUmargin[i]);

    for (i=0;i<P_DIMS;i++) sliceStrideFull[1][i+1]=sliceStrideFull[1][i]*bT::getValCoord_P(i).size();
    for (i=0;i<U_DIMS;i++) sliceStrideFull[0][i+1]=sliceStrideFull[0][i]*bT::getValCoord_U(i).size();
    
    sliceCount[1]=sliceStride[1][P_DIMS];
    sliceCount[0]=sliceStride[0][U_DIMS];
    sliceCountFull[1]=sliceStrideFull[1][P_DIMS];
    sliceCountFull[0]=sliceStrideFull[0][U_DIMS];
    //printf("[%ld][%ld][%ld][%ld]\n",sliceCount[0],sliceCount[1],sliceCountFull[0],sliceCountFull[1]);
    //printf("lpomm %ld %ld -> %ld %ld\n",sliceCount[0],sliceCount[1],sliceCountFull[0],sliceCountFull[1]);
  }

public:

  long sliceIndex(const int coord[DIMS], int type=0, bool full=false)
  {
    long res=0;
    long *stride=(full)?(sliceStrideFull[type]):(sliceStride[type]);
    const int *co=(type)?coord:(&coord[P_DIMS]);

    for (int i=0;i<P_DIMS;i++) res+=co[i]*stride[i];
  
    return res;
  }

  long nSlices(int type=0, bool full=false)
  {
    if ((type>1)||(type<0)) return 0;

    if (full)
      return sliceCountFull[type];
    else
      return sliceCount[type];
  }

  sliceIterator slice_begin(long j, int type=0, bool full=false)
  {
    std::vector<int> coord(DIMS,-1);
    dirT dir=(full)?gridNav::undefined():gridNav::dir();
    long *str=(full)?sliceStrideFull[type]:sliceStride[type];
    int *co=&coord[(1-type)*P_DIMS];
  
    co[P_DIMS-1]=j/str[P_DIMS-1];
    for (int i=P_DIMS-2;i>=0;i--) 
      {
	j-=co[i+1]*str[i+1];
	co[i]=j/str[i];
      }
    
    return sliceIterator(this,&coord[0],dir);
  }

  sliceIterator slice_end(long j, int type=0, bool full=false)
  {
    return sliceIterator(this,NULL,(full)?gridNav::undefined():gridNav::dir(),true);
  }

  
  value_type *fullSlicePtr(long j, int type)
  {
    if (type!=0) return NULL;
    return &bT::arr[j*bT::stride_Val[P_DIMS+1]];
  }

};

#endif
