#ifndef __GRID_21_HXX__
#define __GRID_21_HXX__

#include <iterator>
#include <utility>
#include <string>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "grid_base.hxx"
#include "types.h"

template <
  typename dataType,
  valLocationType valLocationTemplate
  >
class grid21 : public gridBase<dimTraits<1,2>,dataType,valLocationTemplate> {
public:
  static std::string getTag() {return "GRID_21 v0.1";}

  typedef dataType dataT;

  typedef grid21<dataT,valLocationTemplate> myT;
  typedef gridBase<dimTraits<1,2>,dataT,valLocationTemplate> bT;
  
  typedef typename bT::scaleT scaleT;
  typedef typename bT::scaleTypeT scaleTypeT;
  typedef typename bT::scaleTypeV scaleTypeV; 
  typedef typename bT::valLocationT valLocationT;
  typedef typename bT::valLocationV valLocationV;
  typedef typename bT::paramT paramT;

  static const int P_DIMS = bT::P_DIMS;
  static const int U_DIMS = bT::U_DIMS;
  static const int J_DIMS = bT::J_DIMS;
  static const valLocationT valLocation = bT::valLocation;

protected:
  std::vector<double> &R_C;
  std::vector<double> &U_C;
  std::vector<double> &J_C;
  std::vector<double> &R_V;
  std::vector<double> &U_V;
  std::vector<double> &J_V;  
  std::vector<double> &R_Val;
  std::vector<double> &U_Val;
  std::vector<double> &J_Val;  

  long highRmargin;
  long lowRmargin;
  long highUmargin;
  long lowUmargin;
  long highJmargin;
  long lowJmargin;

public:
  grid21():bT(),
	   R_C(bT::P_Cell[0]),
	   U_C(bT::U_Cell[0]),
	   J_C(bT::U_Cell[1]),
	   R_V(bT::P_Vert[0]),
	   U_V(bT::U_Vert[0]),
	   J_V(bT::U_Vert[1]),
	   R_Val(bT::P_Val[0]),
	   U_Val(bT::U_Val[0]),
	   J_Val(bT::U_Val[1])
  {
    
  }

  ~grid21()
  {
    
  }

  virtual void write(FILE *f)
  {
    bT::write(f);
    myIO::writeTag(f,getTag());  
  }

  virtual void read(FILE *f, bool swap)
  {
    bT::read(f,swap);
    myIO::checkTag(f,getTag()); 
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
    highRmargin=bT::p.highPmargin[0];
    lowRmargin=bT::p.lowPmargin[0];
    highUmargin=bT::p.highUmargin[0];
    lowUmargin=bT::p.lowUmargin[0];
    highJmargin=bT::p.highUmargin[1];
    lowJmargin=bT::p.lowUmargin[1];
  }

public:
  friend class sliceIterator;
  class sliceIterator
  {
  protected: 
    long j,r,u;
    long rmax,umax;
    long i;
    myT *g; 
    dataT *data;
  public:

    typedef std::forward_iterator_tag iterator_category;
    typedef dataT value_type;
    typedef dataT* pointer;
    typedef dataT& reference;
    
    sliceIterator():i(-1),g(NULL){}

    sliceIterator(myT *g_, long j_, valLocationT vlt=valLocation, bool end=false):
      g(g_),data(g_->getDataPtr()),j(j_)
    {
      if (end) i=-1;
      else {
	r=g->lowRmargin;
	u=g->lowUmargin;
	if (valLocation==valLocationV::VERTEX)
	  {
	    rmax=g->R_V.size();
	    umax=g->U_V.size();
	    i=j*g->getVertStride(2)+u*g->getVertStride(1)+r;
	  }
	else
	  {
	    rmax=g->R_C.size();
	    umax=g->U_C.size();
	    i=j*g->getCellStride(2)+u*g->getCellStride(1)+r;
	  }
	//printf("new(%ld) : i=%ld\n",j,i);
      }
    }

    sliceIterator(const sliceIterator &it):
      g(it.g),data(it.data),j(it.j),i(it.i),r(it.r),u(it.u),rmax(it.rmax),umax(it.umax)
    {}

    sliceIterator& operator=(const sliceIterator& it)
    {g=it.g;data=it.data;j=it.j;i=it.i;r=it.r;u=it.u;rmax=it.rmax;umax=it.umax;}
	
    ~sliceIterator(){}

    void print() const {printf("i=%ld: R[%ld]=%e, U[%ld]=%e, J[%ld]=%e.\n",i,get_r(),get_R(),get_u(),get_U(),get_j(),get_J());}
    long get_r() const {return r/*-g->lowRmargin*/;}
    long get_u() const {return u/*-g->lowUmargin*/;}
    long get_j() const {return j;}
    double get_R() const {return g->R_Val[r];}
    double get_U() const {return g->U_Val[u];}
    double get_J() const {return g->J_Val[j];}
    double get_dR() const {return (r>=g->R_V.size()-g->highRmargin-1)?0:(g->R_V[r+1]-g->R_V[r]);}
    double get_dU() const {return (u>=g->U_V.size()-g->highUmargin-1)?0:(g->U_V[u+1]-g->U_V[u]);}
    double get_dJ() const {return g->J_V[j+1]-g->J_V[j];}
  
    long get_i() const {return i;}

    double get_avg() const {
      assert(r<g->R_Val.size()-g->highRmargin-1);
      assert(u<g->U_Val.size()-g->highUmargin-1);
      return (data[i]+data[i+1]+data[i+g->R_Val.size()]+data[i+g->R_Val.size()+1])*0.25;
    }


    sliceIterator operator++()
    { 
      if (i<0) return *this;
      i++;r++;
      if (r>=rmax-g->highRmargin) 
	{
	  r=g->lowRmargin;
	  u++;
	  if (u>=umax-g->highUmargin) i=-1;
	  else i+=g->highRmargin+g->lowRmargin;
	}
      return *this;
    }

    sliceIterator operator++(int)
    {
      sliceIterator it(*this);
      ++(*this);
      return it;
    }

    reference operator*() const
    { return data[i]; }

    bool operator==(const sliceIterator& r) const
    {return (i==r.i);}

    bool operator!=(const sliceIterator& r) const
    {return (i!=r.i);}

    bool operator<(const sliceIterator& r) const
    {return (i<r.i);}

    bool operator>(const sliceIterator& r) const
    {return (i>r.i);}	  

    bool operator<=(const sliceIterator& r) const
    {return (i<=r.i);}	  	   

    bool operator>=(const sliceIterator& r) const
    {return (i>=r.i);}
  };

  sliceIterator slice_begin(long j, valLocationT vlt=valLocation)
  {
    return sliceIterator(this,j,vlt);
  }

  sliceIterator slice_end(long j, valLocationT vlt=valLocation)
  {
    return sliceIterator(this,j,vlt,true);
  }

  dataT *fullSlicePtr(long j)
  {
    return &bT::arr[j*bT::stride_Val[2]];
  }

  long nSlices()
  {
    return J_Val.size();
  }

public:

  const std::vector<double>& getValCoord21_R() const {return R_Val;}
  const std::vector<double>& getValCoord21_U() const {return U_Val;}
  const std::vector<double>& getValCoord21_J() const {return J_Val;}
  const std::vector<double>& getCellCoord21_R() const {return R_C;}
  const std::vector<double>& getCellCoord21_U() const {return U_C;}
  const std::vector<double>& getCellCoord21_J() const {return J_C;}
  const std::vector<double>& getVertCoord21_R() const {return R_V;}
  const std::vector<double>& getVertCoord21_U() const {return U_V;}
  const std::vector<double>& getVertCoord21_J() const {return J_V;}
  
  std::vector<double>& computeM_R(std::vector<double> &Mr,double massBelowR0=0)
  {
    Mr.resize(R_Val.size());
    computeM_R(&Mr[0],massBelowR0);
    return Mr;
  }

  double *computeM_R(double *Mr=NULL,double massBelowR0=0)
  { 
  
    long j,u,r,i;
    double R,U,J;
    double dJ,dU,dR;
    
    dataT *data=bT::arr;
    const long Msize=bT::P_Vert[0].size();

    if (Mr==NULL) Mr=new double[Msize];
    for (r=0;r<Msize;r++) Mr[r]=0;
    
#pragma omp parallel for shared(data,Mr) private(J,U,R,j,u,r,i,dJ,dR,dU)
    for (j=0;j<J_C.size();j++)
      {
	i=j*bT::stride_Cell[P_DIMS+U_DIMS-1];
	J=J_C[j];dJ=(J_V[j+1]-J_V[j]);	
	double Ct1=8*M_PI*M_PI*J*dJ;
	for (u=0;u<U_C.size();u++)
	  {
	    dU=U_V[u+1]-U_V[u];
	    double Ct=Ct1*dU;
	    for (r=0;r<R_C.size();r++,i++) 
	      {
		dR=R_V[r+1]-R_V[r];
		
		double tmp=data[i]*dR*Ct;
		//if (tmp<0) tmp=0;
		#pragma omp atomic
		Mr[r+1]+= tmp;
	      }
	  }
      }
    Mr[0]=0*massBelowR0;
    for (r=0;r<Msize-1;r++) Mr[r+1]+=Mr[r];
    
    return Mr;
  }

  std::vector<double>& computeRho_R(std::vector<double> &Rhor,double massBelowR0=0)
  {
    Rhor.resize(bT::P_Cell[0].size());
    computeRho_R(&Rhor[0],massBelowR0);
    return Rhor;
  }

  double *computeRho_R(double *Rhor=NULL,double massBelowR0=0)
  {     
    long j,u,r,i;
    double R,U,J;
    double dJ,dU,dR;
    
    
    dataT *data=bT::arr;
    const long RhoSize = bT::P_Cell[0].size();
    
    if (Rhor==NULL) Rhor=new double[RhoSize];
    for (r=0;r<RhoSize;r++) Rhor[r]=0;
     
    #pragma omp parallel for shared(data,Rhor) private(J,U,R,j,u,r,i,dJ,dR,dU)
    for (j=0;j<J_C.size();j++)
      {
	i=j*bT::stride_Cell[P_DIMS+U_DIMS-1];
	J=J_C[j];dJ=(J_V[j+1]-J_V[j]);	
	double Ct1=2*M_PI*J*dJ;
	for (u=0;u<U_C.size();u++)
	  {
	    dU=U_V[u+1]-U_V[u];
	    double Ct=Ct1*dU;
	    for (r=0;r<R_C.size();r++,i++) 
	      {
		/*
		double tmp = 
		  data[i]+
		  data[i+stride[1]]+
		  data[i+stride[2]]+
		  data[i+stride[2]+stride[1]];
		tmp = tmp * 0.25*Ct;
		*/
		double tmp=data[i]*Ct;
		//if (tmp<0) tmp=0;
		#pragma omp atomic
		Rhor[r]+=tmp;
	      }
	  }
	
      }
    /*
    if (R_V[0]>0) {
      R=R_V[0];
      Rhor[0]=massBelowR0/(4./3.*M_PI*R*R*R);
      }*/
    for (r=0;r<RhoSize;r++) {
      //R=0.5*(R_V[r+1]+R_V[r]);//printf("%f %f %f\n",R_V[r-1],R,R_V[r]);
      R=R_C[r];
      Rhor[r]/=(R*R);
    } 
    return Rhor;
  }
  
};

#endif
