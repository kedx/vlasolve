#ifndef __INIT_GEN_HERNQUIST_HXX__
#define __INIT_GEN_HERNQUIST_HXX__

#include <vector>
#include <map> 
#include <string>
#include <math.h>

#include "global.h"
#include "initGen.hxx"

template <class gridType>
class Hernquist : public initGen<gridType> {

public:
  typedef initGen<gridType> base; 
  typedef typename base::dataT dataT;
  typedef typename base::gridT gridT;
  typedef typename base::dimTr dimTr;
  typedef typename base::scaleT scaleT;
  typedef typename base::scaleType scaleType;
  typedef typename base::scaleTypeVal scaleTypeVal;  
  typedef typename base::gridParamT gridParamT;
  typedef typename base::solverParamT solverParamT;

  struct Hernquist_params {
    double G;
    double M;
    double R0;
    double Ra;
  };

  typedef Hernquist_params paramT; 

private:
void registerParams() {
    p.G=1;
    p.R0=1;
    p.Ra=-1;// anisotropy radius, set to -1 for no anisotropy
    p.M=1;
    
    pMapD.insert(std::make_pair(std::string("G"),&p.G));
    pMapD.insert(std::make_pair(std::string("R0"),&p.R0));
    pMapD.insert(std::make_pair(std::string("M"),&p.M));    
  }

public:

  std::string getTypeStr() {
    return std::string("Hernquist");
  }

  gridParamT defaultGridParams() {
    gridParamT gp;

    gp.Pmin[0]=0.01;
    gp.Pmax[0]=25*10;

    gp.Vmin[0]=-1.45;
    gp.Vmax[0]=1.45;
    gp.Vmin[1]=1.E-4;
    gp.Vmax[1]=3.24*4;
    
    gp.Pres[0]=45*2;
    gp.Vres[0]=40*2;
    gp.Vres[1]=20*4;//20*2*5; 

    gp.Pscale[0]=scaleTypeVal::LOGARITHMIC;
    gp.Vscale[0]=scaleTypeVal::REGULAR;
    gp.Vscale[1]=scaleTypeVal::LOGARITHMIC;
    return gp;
  }

  solverParamT defaultSolverParams() {
    solverParamT sp;   
    sp.T0=0;
    sp.Tmax=62; // no Tmax
    sp.dt=0.025;
    return sp;
  }

private:
  paramT p;
  std::map<std::string, double*> pMapD;
  std::map<std::string, char *> pMapC;
  double Vg;
  
private:

  void ComputeConstants() {
    Vg=sqrt((p.G*p.M)/(p.R0));
  }

  double computePhi(double R) {
    return -(p.G*p.M)/(p.R0+R);
  }

  double computeRho(double R) {
    return (p.M*p.R0)/(2.*M_PI*R)*pow(p.R0+R,-3);
  }

  std::vector<double> computePhi(std::vector<double> &R_) 
  {
    std::vector<double> Phi;
    
    Phi.resize(R_.size());
    long i;
    for (i=0;i<R_.size();i++)
      Phi[i]=computePhi(R_[i]);

    return Phi;
  }

  double computeF(double phi, double r, double u, double j)
  {
    double eps=phi+(u*u+(j*j)/(r*r))*0.5;
    
    if (p.Ra>0) {
      eps = eps+(j*j)/(2*p.Ra*p.Ra);
    }
    if (eps>=0) return 0;

    double q=sqrt(-(p.R0)/(p.G*p.M) * eps); 
    double q2=q*q;
    double A1=(p.M)/(8.*sqrt(2.)*pow(M_PI*p.R0*Vg,3.))*pow(1.-q2,-2.5);
    double A2=3*asin(q)+q*sqrt(1.-q2)*(1.-2*q2)*(8*q2*q2-8*q2-3);
    double result = A1*A2;

    if (p.Ra>0) {
      A1=(p.M)/(sqrt(2)*pow(M_PI*p.R0*Vg,3))*pow(p.R0/p.Ra,2);
      A2=q*(1.-2*q2);
      result+=A1*A2;
    }
    
    return result;
  } 
  
  template <class DT>
  dataT valueAt(std::vector<double>& pos, dimTraits<1,2>)
  {
    double phi=computePhi(pos[0]);
    return computeF(phi,pos[0],pos[1],pos[2]);
  }
  

  template <class DT>
  gridT *generate(gridParamT &gp,solverParamT &sp, dimTraits<1,2>)
  {
    printf("Generating initial conditions (%s) ... ",getTypeStr().c_str());fflush(0);
    gridT *g=new gridT(gp);
    dataT *d=g->getDataPtr();    
    
    long r,u,j,i;
    double R,U,J;    
    std::vector<std::vector<double> > P_Cell = g->getCellCoord_P();
    std::vector<std::vector<double> > P_Vert = g->getVertCoord_P();
    std::vector<std::vector<double> > U_Cell = g->getCellCoord_U();
    std::vector<std::vector<double> > U_Vert = g->getVertCoord_U();
    std::vector<double> &R_C=P_Cell[0];
    std::vector<double> &U_C=U_Cell[0];
    std::vector<double> &J_C=U_Cell[1];
    std::vector<double> &R_V=P_Vert[0];
    std::vector<double> &U_V=U_Vert[0];
    std::vector<double> &J_V=U_Vert[1];
    long delta = g->getCellStride(2);
    
    std::vector<double> Phi=computePhi(R_C);

    #pragma omp parallel for shared(R_C,U_C,J_C,R_V,U_V,J_V,delta,d,Phi) private(J,U,R,j,u,r,i)
    for (j=0;j<J_C.size();j++)
      {
	i=j*delta;
	J=J_C[j];
	for (u=0;u<U_C.size();u++)
	  {
	    U=U_C[u];
	    for (r=0;r<R_C.size();r++,i++) 
	      {
		R=R_C[r];
    		d[i]=computeF(Phi[r],R,U,J);		
	      }
	  }
      }
    
    double massBelowR0=0;
    double dr=R_V[0]/10000.;
    for (R=0;R<R_V[0];R+=dr)
      {
	massBelowR0 += 4./3.*M_PI*((R+dr)*(R+dr)*(R+dr) - R*R*R)*computeRho(R+0.5*dr);
      }
    
    g->registerValue("massBelowR0",massBelowR0);
    g->registerValue("R0",p.R0);
    
    std::vector<double> Mr;
    g->computeM_R(Mr);
    for (i=0;i<R_C.size()*U_C.size()*J_C.size();i++) d[i] *= p.M/Mr.back();
    g->registerValue("INIGEN_NormFactor",p.M/Mr.back());

    printf("done.\n");
    return g;
  }



  /* NO NEED TO MODIFY BELOW HERE */

public:

  dataT valueAt(std::vector<double> &pos) {
    return valueAt<dimTr>(pos);
  }

  gridT *generate(gridParamT &gp, solverParamT &sp) {
    return generate<dimTr>(gp,sp);
  }

  gridT *generate(solverParamT &sp) {
    gridParamT gp=defaultGridParams();
    return generate(gp,sp);
  }

  gridT *generate(gridParamT &gp) {
    solverParamT sp=defaultSolverParams();
    return generate(gp,sp);
  }

  gridT *generate() {
    gridParamT gp=defaultGridParams();
    solverParamT sp=defaultSolverParams();
    return generate(gp,sp);
  }

  void setParams(std::string &str) {
    
  }

  void setParams(paramT &p_) {
    p=p_;
    ComputeConstants();
  }

  Hernquist():base() {
    registerParams();
    /*
    p.G=1;
    p.R0=1;
    p.M=1;
    */
    ComputeConstants();
  }
  /*
  Hernquist(paramT &p_) : base(){
    registerParams();
    p=p_;
    ComputeConstants();
  }
  */
  ~Hernquist() {
  }



private:
  template <class DT> 
  dataT valueAt(std::vector<double>& pos){
    return valueAt<DT>(pos,DT());
  }

  template <class DT>
  dataT valueAt(std::vector<double>& pos, DT){
    fprintf(stderr,"ERROR generating initial conditions :\n");
    fprintf(stderr,"  %s: not implemented in (%d+%d)D.\n",getTypeStr().c_str(),DT::P_DIMS,DT::V_DIMS);
    exit(0);
  }

   template <class DT> 
   gridT *generate(gridParamT &gp,solverParamT &sp){
     return generate<DT>(gp,sp,DT());
  }
  
  template <class DT>
  gridT *generate(gridParamT &gp,solverParamT &sp, DT){
    fprintf(stderr,"ERROR generating initial conditions :\n");
    fprintf(stderr,"  %s: not implemented in (%d+%d)D.\n",getTypeStr().c_str(),DT::P_DIMS,DT::V_DIMS);
    exit(0);
  }
};


#endif
  


