#ifndef __SLICE21_ITERATOR_HXX__
#define __SLICE21_ITERATOR_HXX__

#include <iterator>
#include "fieldIterator.hxx"

template <class gridType>
class slice21Iterator
{

public:

  typedef slice21Iterator<gridType> self_type;
  typedef gridType container_type;

  typedef std::forward_iterator_tag iterator_category;
  typedef typename container_type::value_type value_type;
  typedef value_type* pointer;
  typedef value_type& reference;
  typedef const value_type* const_pointer;
  typedef const value_type& const_reference;

  typedef typename container_type::valLocationT valLocationT;
  typedef typename container_type::valLocationV valLocationV;

  static const int P_DIMS = container_type::P_DIMS;
  static const int U_DIMS = container_type::U_DIMS;
  static const int DIMS = container_type::DIMS;
  //static const valLocationT valLocation = container_type::valLocation;

protected: 
  long j;  
  long r,u;
  long rmax,umax;
  long nfields; 
  long rval_delta;
  long i;
  int lrm,hrm,lum,hum;

  value_type *curD;
  container_type *g; 

public:
    
  slice21Iterator():curD(NULL),i(-1),g(NULL){}

  slice21Iterator(container_type *g_, long j_, bool full=false, bool end=false):
    g(g_),curD(g_->getDataPtr()),nfields(g_->getNFields()),j(j_)
  {  
    rval_delta=g->R_Val.size()*nfields;
    if (end) {i=-1;curD=NULL;}
    else {
      lrm=(full)?0:g->lowRmargin;
      hrm=(full)?0:g->highRmargin;
      lum=(full)?0:g->lowUmargin;
      hum=(full)?0:g->highUmargin;

      r=lrm;
      u=lum;
      
      rmax=g->R_Val.size();
      umax=g->U_Val.size();
      i=(j*g->getValStride(2)+u*g->getValStride(1)+r);
      curD+=i*nfields;
    }
  }
  /*
  slice21Iterator(const self_type &it):
    g(it.g),curD(it.curD),j(it.j),i(it.i),r(it.r),u(it.u),rmax(it.rmax),umax(it.umax)
  {}
  
  self_type& operator=(const self_type& it)
  {g=it.g;curD=it.curD;j=it.j;i=it.i;r=it.r;u=it.u;rmax=it.rmax;umax=it.umax;}
  */
  ~slice21Iterator(){}

  void print() const {printf("i=%ld: R[%ld]=%e, U[%ld]=%e, J[%ld]=%e.\n",get_i(),get_r(),get_R(),get_u(),get_U(),get_j(),get_J());}
  long get_r() const {return r/*-lrm*/;}
  long get_u() const {return u/*-lum*/;}
  long get_j() const {return j;}
  double get_R() const {return g->R_Val[r];}
  double get_U() const {return g->U_Val[u];}
  double get_J() const {return g->J_Val[j];}
  double get_dR() const {return (r>=g->R_V.size()-1)?0:(g->R_V[r+1]-g->R_V[r]);}
  double get_dU() const {return (u>=g->U_V.size()-1)?0:(g->U_V[u+1]-g->U_V[u]);}
  double get_dJ() const {return (j>=g->U_V.size()-1)?0:(g->J_V[j+1]-g->J_V[j]);}
 
  long get_i() const {return i;}

  double get_avg() const {
    assert(r<g->R_Val.size()-hrm-1);
    assert(u<g->U_Val.size()-hum-1);
    return ((*curD)+curD[nfields]+curD[rval_delta]+curD[rval_delta+nfields])*0.25;
  }

  double get_Ravg() const {
    assert(r<g->R_Val.size()-hrm-1);
    return ((*curD)+curD[nfields])*0.5;
  }

  double get_Uavg() const {
    assert(u<g->U_Val.size()-hum-1);
    return ((*curD)+curD[rval_delta])*0.5;
  }

  double get_avg(int id) const {
    assert(r<g->R_Val.size()-hrm-1);
    assert(u<g->U_Val.size()-hum-1);
    return (curD[id]+curD[nfields+id]+curD[rval_delta+id]+curD[rval_delta+nfields+id])*0.25;
  }

  double get_Ravg(int id) const {
    assert(r<g->R_Val.size()-hrm-1);
    return (curD[id]+curD[nfields+id])*0.5;
  }

  double get_Uavg(int id) const {
    assert(u<g->U_Val.size()-hum-1);
    return (curD[id]+curD[rval_delta+id])*0.5;
  }


  self_type &operator++()
  { 
    if (curD==NULL) return *this;
    r++;curD+=nfields;i++;
    if (r>=rmax-hrm) 
      {
	r=lrm;
	u++;
	if (u>=umax-hum) {i=-1;curD=NULL;}
	else 
	  {
	    i+=(hrm+lrm);
	    curD+=nfields*(hrm+lrm);
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

  reference operator*() const
  {
    return *curD;
  }

  reference operator[](int id)
  {
    return curD[id];
  }

  bool operator==(const self_type& r) const
  {return (curD==r.curD);}

  bool operator!=(const self_type& r) const
  {return (curD!=r.curD);}

  bool operator<(const self_type& r) const
  {return (curD<r.curD);}

  bool operator>(const self_type& r) const
  {return (curD>r.curD);}	  

  bool operator<=(const self_type& r) const
  {return (curD<=r.curD);}	  	   

  bool operator>=(const self_type& r) const
  {return (curD>=r.curD);}
};

template <>
template <class gridType>
class field_iterator< slice21Iterator<gridType> >
  :public slice21Iterator<gridType>
{
private:
  typedef slice21Iterator<gridType> bT;
  int fct;
public:
  typedef field_iterator<bT> self_type;

  field_iterator(const bT &it):bT(it),fct(0)
  {
    
  }

  self_type &operator++()
  {     
    if (bT::curD==NULL) return *this;
    
    ++bT::curD;
    if ((++fct)==bT::nfields) {fct=0;bT::r++;bT::i++;}
    else return *this;

    if (bT::r>=bT::rmax-bT::hrm) 
      {
	bT::r=bT::lrm;
	bT::u++;
	if (bT::u>=bT::umax-bT::hum) {bT::i=-1;bT::curD=NULL;}
	else 
	  {
	    bT::i+=(bT::hrm+bT::lrm);
	    bT::curD+=bT::nfields*(bT::hrm+bT::lrm);
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
template <class gridType>
class oneField_iterator< slice21Iterator<gridType> > 
  :public slice21Iterator<gridType>
{
private:
  typedef slice21Iterator<gridType> bT;
 
public:
  typedef oneField_iterator<bT> self_type;

  oneField_iterator(const bT &it, int id):bT(it)
  {
    if (bT::curD!=NULL) bT::curD+=id;
  }
}; 


#endif
