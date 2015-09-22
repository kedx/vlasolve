#ifndef __GRID_21_HXX__
#define __GRID_21_HXX__

#include <utility>
#include <string>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "slice21Iterator.hxx"
#include "grid_base.hxx"
#include "defines.h"

template <
  typename dataType,
  class valLocationTraits
  >
class grid21 : public gridBase<dimTraits<1,2>,dataType,valLocationTraits> {
public:

  static std::string getTag() {return "GRID_21 v0.11";}

  typedef dataType value_type;
  typedef grid21<value_type,valLocationTraits> myT;  
  typedef gridBase<dimTraits<1,2>,value_type,valLocationTraits> bT;

  typedef typename bT::iterator iterator;
  friend class slice21Iterator<myT>;
  typedef slice21Iterator<myT> sliceIterator;
  
  typedef typename bT::scaleT scaleT;
  typedef typename bT::scaleTypeT scaleTypeT;
  typedef typename bT::scaleTypeV scaleTypeV; 
  typedef typename bT::valLocationT valLocationT;
  typedef typename bT::valLocationV valLocationV;
  typedef typename bT::paramT paramT;

  static const int P_DIMS = bT::P_DIMS;
  static const int U_DIMS = bT::U_DIMS;
  static const int J_DIMS = bT::J_DIMS;
  static const int DIMS = bT::DIMS;

  //static const valLocationT valLocation_P = bT::valLocation_P;
  //static const valLocationT valLocation_U = bT::valLocation_U;
  //static const valLocationT valLocation_J = bT::valLocation_J;

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
    highRmargin=bT::p.highPmargin[0];
    lowRmargin=bT::p.lowPmargin[0];
    highUmargin=bT::p.highUmargin[0];
    lowUmargin=bT::p.lowUmargin[0];
    highJmargin=bT::p.highUmargin[1];
    lowJmargin=bT::p.lowUmargin[1];
  }

public:

  long nSlices(int type=0, bool full=false)
  {
    if (type!=0) return 0;
    return J_Val.size();
  }

  sliceIterator slice_begin(long j,int type=0, bool full=false)
  {
    if (type!=0) return slice_end(j,0,full);

    if ((j<0)||(j>=J_Val.size()-highJmargin-lowJmargin))
      return slice_end(j,full);
      
    return sliceIterator(this,j+lowJmargin,full);
  }
  
  sliceIterator slice_end(long j, int type=0, bool full=false)
  {
    return sliceIterator(this,j+lowJmargin,full,true);
  }
  
  value_type *fullSlicePtr(long j, int type=0)
  {
    if (type!=0) return NULL;
    return &bT::arr[j*bT::stride_Val[2]];
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

  /*
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
    
    value_type *data=bT::arr;
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
    
    
    value_type *data=bT::arr;
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
		
		// double tmp = 
		//   data[i]+
		//   data[i+stride[1]]+
		//   data[i+stride[2]]+
		//   data[i+stride[2]+stride[1]];
		// tmp = tmp * 0.25*Ct;
		
		double tmp=data[i]*Ct;
		//if (tmp<0) tmp=0;
		#pragma omp atomic
		Rhor[r]+=tmp;
	      }
	  }
	
      }
    
    // if (R_V[0]>0) {
    //   R=R_V[0];
    //   Rhor[0]=massBelowR0/(4./3.*M_PI*R*R*R);
    //   }
    for (r=0;r<RhoSize;r++) {
      //R=0.5*(R_V[r+1]+R_V[r]);//printf("%f %f %f\n",R_V[r-1],R,R_V[r]);
      R=R_C[r];
      Rhor[r]/=(R*R);
    } 
    return Rhor;
  }
  */
};

#endif
