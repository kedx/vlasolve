#ifndef __GRID_HXX__
#define __GRID_HXX__

#include <map>
#include <string>
#include "scale.hxx"

template <class dimensionTraits,
	  typename dataType=double, 
	  class scaleClass=scale<dataType>
	  >
class regularGrid
{ 
public:
  typedef dimensionTraits dimTr;

  static const int DIMS = dimensionTraits::P_DIMS;
  static const int P_DIMS= dimensionTraits::P_DIMS;
  static const int V_DIMS= dimensionTraits::V_DIMS;
  static const int J_DIMS= V_DIMS-P_DIMS;
  
  typedef dataType dataT;
  typedef scaleClass scaleT;  
  typedef typename scaleT::scaleType scaleType;
  typedef typename scaleT::scaleTypeVal scaleTypeVal;

  struct gridParams {
    double Pmin[P_DIMS];
    double Pmax[P_DIMS];
    double Vmin[V_DIMS];
    double Vmax[V_DIMS];
    int Pres[P_DIMS];
    int Vres[V_DIMS];
    scaleType Pscale[P_DIMS];
    scaleType Vscale[V_DIMS];
  };

  typedef gridParams paramT;
  
private:  
  typedef std::map<std::string,double> valueMapT;
  typedef typename valueMapT::iterator valueMapItT;
    
  long n_Vert;
  long n_Cell;
  
  std::vector<double> P_Vert[P_DIMS];
  std::vector<double> U_Vert[V_DIMS];
  std::vector<double> P_Cell[P_DIMS];
  std::vector<double> U_Cell[V_DIMS];

  long stride_Vert[P_DIMS+V_DIMS+1];
  long stride_Cell[P_DIMS+V_DIMS+1];
  valueMapT valueMap;

  paramT p;
  dataT *arr;

  void init()
  {    
    int i,j;

    n_Vert=1; 
    n_Cell=1; 
    stride_Vert[0]=1;
    stride_Cell[0]=1;

    for (i=0;i<P_DIMS;i++ ) {
      
      P_Vert[i]=scaleT::genScale(p.Pmin[i],p.Pmax[i],p.Pres[i]+1,p.Pscale[i]);
      P_Cell[i]=scaleT::genScale(p.Pmin[i],p.Pmax[i],p.Pres[i]+1,p.Pscale[i],
				 scaleT::valLocationTypeVal::PIXEL);
      n_Vert*=(p.Pres[i]+1);
      n_Cell*=p.Pres[i];   
      //printf("%ld %ld \n",P_Vert[i].size(),P_Cell[i].size());
      stride_Cell[i+1]=stride_Cell[i]*p.Pres[i];
      stride_Vert[i+1]=stride_Vert[i]*(p.Pres[i]+1);
      //P_Cell[i].resize(P_Vert[i].size()-1);
      //for (j=0;j<P_Cell[i].size();j++)  P_Cell[i][j]=(P_Vert[i][j]+P_Vert[i][j+1])*0.5;
    }
    for (i=0;i<V_DIMS;i++ ) {
      
      U_Vert[i]=scaleT::genScale(p.Vmin[i],p.Vmax[i],p.Vres[i]+1,p.Vscale[i]);
      U_Cell[i]=scaleT::genScale(p.Vmin[i],p.Vmax[i],p.Vres[i]+1,p.Vscale[i],
				 scaleT::valLocationTypeVal::PIXEL);
      //printf("%ld %ld \n",U_Vert[i].size(),U_Cell[i].size());
      n_Vert*=(p.Vres[i]+1);
      n_Cell*=p.Vres[i];  
      stride_Cell[i+P_DIMS+1]=stride_Cell[i+P_DIMS]*p.Vres[i];
      stride_Vert[i+P_DIMS+1]=stride_Vert[i+P_DIMS]*(p.Vres[i]+1);
      //U_Cell[i].resize(U_Vert[i].size()-1);
      //for (j=0;j<U_Cell[i].size();j++)  U_Cell[i][j]=(U_Vert[i][j]+U_Vert[i][j+1])*0.5;
    }
    
    arr = new dataT[n_Cell];
    //printf("allocated %ld nodes (%.3fMo)\n",nNodes,(double)(8*nNodes)/(1024*1024));
  }

public:
  bool registerValue(const char *name, double value, bool replace=true)
  {
    std::string pname(name);
    std::pair<valueMapItT,bool> res=
      valueMap.insert(std::make_pair(pname,value));

    if ((!res.second)&&(replace)) {
      valueMap.erase(res.first);
      res=valueMap.insert(std::make_pair(pname,value));
      return false;
    }
    
    return true;
  }
  
  bool getRegisteredValue(const char *name, double &value)
  {
    std::string pname(name);
    valueMapItT res=valueMap.find(pname);
    if (res==valueMap.end()) return false;
    value=res->second;
    return true;
  }

  void set(paramT &p_){    
    p=p_;
    if (arr!=NULL) delete arr;
    arr=NULL;
    init();
  }

  //std::vector<double> get_dP() const {return std::vector<double>(dP,dP+P_DIMS);}
  //std::vector<double> get_dV() const {return std::vector<double>(dV,dV+V_DIMS);}

  long getCellStride(int i) const {return stride_Cell[i];}
  long getVertStride(int i) const {return stride_Vert[i];}

  long getNVertex() const {return stride_Vert[P_DIMS+V_DIMS];}
  long getNCell() const {return stride_Cell[P_DIMS+V_DIMS];}

  std::vector< std::vector<double> > getCellCoord_P() const 
  {
    return std::vector< std::vector<double> > (P_Cell,P_Cell+P_DIMS);
  }
  std::vector< std::vector<double> > getCellCoord_U() const 
  {
    return std::vector< std::vector<double> > (U_Cell,U_Cell+V_DIMS);
  }
  std::vector< std::vector<double> > getVertCoord_P() const 
  {
    return std::vector< std::vector<double> > (P_Vert,P_Vert+P_DIMS);
  }
  std::vector< std::vector<double> > getVertCoord_U() const 
  {
    return std::vector< std::vector<double> > (U_Vert,U_Vert+V_DIMS);
  }

  paramT get_Params() const {return p;}

  //double get_nNodes() const {return nNodes;} 
  dataT *getDataPtr() {return arr;}

public:
  
  regularGrid() {
    arr=NULL;
  }

  regularGrid(paramT &p_) {
    arr=NULL;
    p=p_;
    init();
  }

  ~regularGrid() {
    delete arr;
  }

public:
  std::vector<double>& computeM_R(std::vector<double> &Mr,double massBelowR0=0)
  {
    Mr.resize(P_Vert[0].size());
    computeM_R(&Mr[0],massBelowR0);
    return Mr;
  }

  double *computeM_R(double *Mr=NULL,double massBelowR0=0)
  {     
    long j,u,r,i;
    double R,U,J;
    double dJ,dU,dR;
    std::vector<double> &R_C=P_Cell[0];
    std::vector<double> &U_C=U_Cell[0];
    std::vector<double> &J_C=U_Cell[1];
    std::vector<double> &R_V=P_Vert[0];
    std::vector<double> &U_V=U_Vert[0];
    std::vector<double> &J_V=U_Vert[1];
    dataT *data=arr;
    const long Msize = P_Vert[0].size();

    if (Mr==NULL) Mr=new double[Msize];
    for (r=0;r<Msize;r++) Mr[r]=0;
    
#pragma omp parallel for shared(R_C,U_C,J_C,R_V,U_V,J_V,data,Mr) private(J,U,R,j,u,r,i,dJ,dR,dU)
    for (j=0;j<J_C.size();j++)
      {
	i=j*stride_Cell[P_DIMS+V_DIMS-1];
	J=J_C[j];dJ=(J_V[j+1]-J_V[j]);	
	double Ct1=8*M_PI*M_PI*J*dJ;
	for (u=0;u<U_C.size();u++)
	  {
	    dU=U_V[u+1]-U_V[u];
	    double Ct=Ct1*dU;
	    for (r=0;r<R_C.size();r++,i++) 
	      {
		dR=R_V[r+1]-R_V[r];
		/*
		double tmp = 
		  data[i]+
		  data[i+stride[0]]+
		  data[i+stride[1]]+
		  data[i+stride[1]+stride[0]]+
		  data[i+stride[2]]+
		  data[i+stride[2]+stride[0]]+
		  data[i+stride[2]+stride[1]]+
		  data[i+stride[2]+stride[1]+stride[0]];
		tmp = tmp * 0.125*dR*Ct;
		*/
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
    
    Rhor.resize(P_Cell[0].size());
    computeRho_R(&Rhor[0],massBelowR0);
    return Rhor;
  }

  double *computeRho_R(double *Rhor=NULL,double massBelowR0=0)
  {     
    long j,u,r,i;
    double R,U,J;
    double dJ,dU,dR;
    std::vector<double> &R_C=P_Cell[0];
    std::vector<double> &U_C=U_Cell[0];
    std::vector<double> &J_C=U_Cell[1];
    std::vector<double> &R_V=P_Vert[0];
    std::vector<double> &U_V=U_Vert[0];
    std::vector<double> &J_V=U_Vert[1];
    dataT *data=arr;
    const long RhoSize = P_Cell[0].size();
    
    if (Rhor==NULL) Rhor=new double[RhoSize];
    for (r=0;r<RhoSize;r++) Rhor[r]=0;
     
    #pragma omp parallel for shared(R_C,U_C,J_C,R_V,U_V,J_V,data,Rhor) private(J,U,R,j,u,r,i,dJ,dR,dU)
    for (j=0;j<J_C.size();j++)
      {
	i=j*stride_Cell[P_DIMS+V_DIMS-1];
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
