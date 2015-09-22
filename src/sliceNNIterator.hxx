#ifndef __SLICENN_ITERATOR_HXX__
#define __SLICENN_ITERATOR_HXX__

#include <iterator>
#include "subgridIterator.hxx"
#include "fieldIterator.hxx"
#include "helpers.hxx"

template <class gridType>
class sliceNNIterator : public subgridIterator<gridType>
{
public:
  template <class gridType2> friend class subgridIterator;
  

  typedef subgridIterator<gridType> bT;

  typedef sliceNNIterator<gridType> self_type;
  typedef typename bT::container_type container_type;
  typedef typename bT::iterator_category iterator_category;
  typedef typename bT::value_type value_type;
  typedef typename bT::pointer pointer;
  typedef typename bT::reference reference;
  typedef typename bT::const_pointer const_pointer;
  typedef typename bT::const_reference const_reference;

  typedef typename bT::valLocationT valLocationT;
  typedef typename bT::valLocationV valLocationV;
  typedef typename gridNav::dirT dirT;

  static const int P_DIMS = bT::P_DIMS;
  static const int U_DIMS = bT::U_DIMS;
  static const int DIMS = bT::DIMS;

  static const valLocationT valLocation = bT::valLocation;
  

protected:
  int dirID;
  int dim[DIMS];
  int nd;
  long stride0;

  int *pmin[DIMS];
  int *pmax[DIMS];
  
  void construct(const int *coord)
  {
    if (bT::data==NULL) return;
    if (coord==NULL) return;

    long i,j,k;
    int ndf=0;
    int dimF[DIMS];
    int cF[DIMS];
    long stride[DIMS+1]; 
    long strideF[DIMS+1];    
    dirID=0;
    for (i=0,j=0,k=0;i<DIMS;i++)
      {
	if (coord[i]>=0) 
	  {
	    cF[i]=coord[i];
	    strideF[j]=bT::g->getValStride(i)/bT::nfields;
	    strideF[j+1]=bT::g->getValStride(i+1)/bT::nfields;
	    dimF[j++]=i;
	    dirID|=(1<<i);	    
	  }
	else 
	  {
	    stride[k]=bT::g->getValStride(i)/bT::nfields;
	    stride[k+1]=bT::g->getValStride(i+1)/bT::nfields;
	    dim[k++]=i;	    
	  }
      }
    nd=k;ndf=j;

    for (i=nd;i<DIMS;i++) 
      {
	dim[i]=DIMS;
	stride[i+1]=0;
      }
   
    for (i=ndf;i<DIMS;i++) 
      {
	dimF[i]=DIMS;
	strideF[i+1]=0;
      }
    
    for (i=0;i<ndf;i++)
      {
	bT::i+=(coord[dimF[i]]-bT::w[dimF[i]])*strideF[i];
	bT::data+=(coord[dimF[i]]-bT::w[dimF[i]])*strideF[i]*bT::nfields;
	bT::w[dimF[i]]=coord[dimF[i]];
	bT::min[dimF[i]]=bT::w[dimF[i]];
	bT::max[dimF[i]]=bT::w[dimF[i]]+1;
	//bT::delta[i]=stride[i+1]-(bT::max[i]-bT::min[i])*stride[i];
      }

    for (i=0;i<nd;i++)
      {	
	//bT::min[i]=bT::min[dim[i]];
	//bT::max[i]=bT::max[dim[i]];
	pmin[i]=&bT::min[dim[i]];
	pmax[i]=&bT::max[dim[i]];
	bT::delta[i]=stride[i+1]-(*pmax[i]-*pmin[i])*stride[i];
      }
    for (i=nd;i<DIMS;i++) 
      {
	//bT::min[i]=cF[i-nd];
	//bT::max[i]=cF[i-nd]+1;
	pmin[i]=pmax[i]=NULL;
      }

    stride0=stride[0];
  }

public:

  sliceNNIterator():bT()
  {
  }

  sliceNNIterator(const bT &it, const int *coord):bT(it) // fixer dim et bT::w !!!!
  {
    construct(coord);
  }

  sliceNNIterator(container_type *g, const int *coord, dirT region=gridNav::undefined(), bool end=false,const int *lowPmargin=NULL, const int *highPmargin=NULL,const int *lowUmargin=NULL, const int *highUmargin=NULL):
    bT(g,region,end,lowPmargin,highPmargin,lowUmargin,highUmargin)
  {
    construct(coord);
  }

  self_type &operator++()
  {     
    if (bT::data==NULL) return *this;
    bT::w[dim[0]]++;bT::i+=stride0;bT::data+=stride0*bT::nfields;

    if (bT::w[dim[0]]>=*pmax[0])
      {
	bT::w[dim[0]]=*pmin[0];
	if (dim[1]>=bT::DIMS) {bT::i=-1;bT::data=NULL;}
	else 
	  {
	    bT::w[dim[1]]++;
	    bT::i+=bT::delta[0];bT::data+=bT::delta[0]*bT::nfields;
	    if (bT::w[dim[1]]>=*pmax[1])
	      {
		bT::w[dim[1]]=*pmin[1];
		int ct=2;
		while (dim[ct]<bT::DIMS)
		  {
		    bT::w[dim[ct]]++;bT::i+=bT::delta[ct-1];bT::data+=bT::delta[ct-1]*bT::nfields;
		    if (bT::w[dim[ct]]>=*pmax[ct])
		      {
			bT::w[dim[ct]]=*pmin[ct];
			ct++;
		      }
		    else break;
		  };
		if (dim[ct]>=bT::DIMS) {bT::i=-1;bT::data=NULL;}
	      }
	  }
      }
     
    return *this;   
  }

  void print() const {
    
    double tmp[4];
    this->get_C(tmp);
    bT::print();
    printf("data=%ld,i=%ld: P=[%g %g] U=[%g %g]\n",(long)bT::data,bT::i,tmp[0],tmp[1],tmp[2],tmp[3]);
  }

  self_type operator++(int)
  {
    self_type it(*this);
    ++(*this);
    return it;
  }

  int ndims()const
  {
    return nd;
  }

  double slice_cellVol(bool ortho=false) const
  {
    double result=1.;

    if (ortho)
      {
	int i=0;
	for (int id=0;id<DIMS;id++)
	  {
	    if (dirID&(1<<id)) 
	      result*=bT::g->getCellDelta(id)[bT::w[id]];
	  }
      }
    else
      {
	for (int id=0;id<nd;id++)
	  result*=bT::g->getCellDelta(dim[id])[bT::w[dim[id]]];
	
      }

    return result;
  }

  double slice_cellAvg(int s=0, bool ortho=false) const
  {    
    const std::vector<long> &itPt=bT::g->getIntegrationPointsPtr((!ortho)?dirID:((~dirID)&((1<<DIMS)-1)));    
    double res=0;

    for (int ii=0;ii<itPt.size();ii+=2)
      res+=bT::data[itPt[ii]+s]+bT::data[itPt[ii+1]+s];
  
    return res/itPt.size();
  }

  /*
  double get_sliceAvg(const int s=0, bool ortho=false) const
  {    
   
    const std::vector<long> &itPt=bT::g->getIntegrationPointsPtr((!ortho)?dirID:((~dirID)&((1<<DIMS)-1)));
    double res=0;
   
    for (int ii=0;ii<itPt.size();ii+=2)
      res+=bT::data[itPt[ii]+s]+bT::data[itPt[ii+1]+s];//itPt[i][delta]+itPt[i+1][delta];
    return res/itPt.size();
  }
  */

};

template <>
template <class gridT>
class field_iterator< sliceNNIterator<gridT> >
  :public sliceNNIterator<gridT>
{
private:
  typedef sliceNNIterator<gridT> bT;
  typedef typename bT::bT bbT;
  int fct;
public:
  typedef field_iterator<bT> self_type;

  field_iterator(const bT &it):bT(it),fct(0)
  {
    
  }

  self_type &operator++()
  {     
    if (bbT::data==NULL) return *this;
    ++bbT::data;
    if ((++fct)==bbT::nfields) {fct=0;bbT::w[bT::dim[0]]++;bbT::i+=bT::stride0;bbT::data+=(bT::stride0-1)*bbT::nfields;} 
    else return *this;

    if (bbT::w[bT::dim[0]]>=*bT::pmax[0])
      {
	bbT::w[bT::dim[0]]=*bT::pmin[0];
	
	if (bT::dim[1]>=bbT::DIMS) {bbT::i=-1;bbT::data=NULL;}
	else 
	  {
	    bbT::w[bT::dim[1]]++;
	    bbT::i+=bbT::delta[0];bbT::data+=bbT::delta[0]*bbT::nfields;
	    if (bbT::w[bT::dim[1]]>=*bT::pmax[1])
	      {
		bbT::w[bT::dim[1]]=*bT::pmin[1];
		int ct=2;
		while (bT::dim[ct]<bbT::DIMS)
		  {
		    bbT::w[bT::dim[ct]]++;bbT::i+=bbT::delta[ct-1];bT::data+=bbT::delta[ct-1]*bbT::nfields;
		    if (bbT::w[bT::dim[ct]]>=*bT::pmax[ct])
		      {
			bbT::w[bT::dim[ct]]=*bT::pmin[ct];
			ct++;
		      }
		    else break;
		  };
		if (bT::dim[ct]>=bbT::DIMS) {bbT::i=-1;bbT::data=NULL;}
	      }
	  }
      }
     
    return *this;   
  }

  self_type operator++(int)
  {
    self_type it(*this);
    ++(*this);
    return it;
  }
  
};

template <>
template <class gridT>
class oneField_iterator< sliceNNIterator<gridT> > 
  :public sliceNNIterator<gridT>
{
private:
  typedef sliceNNIterator<gridT> bT;
  typedef typename bT::bT bbT;
 
public:
  typedef oneField_iterator<bT> self_type;

  oneField_iterator(const bT &it, int id):bT(it)
  {
    if (bbT::data!=NULL) bbT::data+=id;
  }
}; 


#endif
