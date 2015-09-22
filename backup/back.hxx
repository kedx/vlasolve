#ifndef __INIT_GEN_UNIFORM_DENSITY_SPHERE_HXX__
#define __INIT_GEN_UNIFORM_DENSITY_SPHERE_HXX__

#include <vector>
#include <map>
#include <string>

#include "global.h"
#include "initGen.hxx"

template <class gridType>
class UniformDensitySphere : public initGen<gridType> {

public:
  typedef initGen<gridType> base; 
  typedef typename base::dataT dataT;
  typedef typename base::gridT gridT;
  typedef typename base::dimTr dimTr;
  typedef typename base::scaleT scaleT;
  typedef typename base::scaleTypeT scaleTypeT;
  typedef typename base::scaleTypeValT scaleTypeValT;  
  typedef typename base::gridParamT gridParamT;
  typedef typename base::solverParamT solverParamT;
  
  struct UDS_params {
    double G;
    double M;
    double R0;
    double alpha;
  };

  typedef UDS_params paramT; 

private:
void registerParams() {
    p.G=1;
    p.R0=2;
    p.M=1;
    p.alpha=0.5; // alpha=1 -> equilibrium (virial theorem satisfied)

    pMapD.insert(std::make_pair(std::string("G"),&p.G));
    pMapD.insert(std::make_pair(std::string("R0"),&p.R0));
    pMapD.insert(std::make_pair(std::string("M"),&p.M));
    pMapD.insert(std::make_pair(std::string("alpha"),&p.alpha));
  }

public:

  std::string getTypeStr() {
    return std::string("UDF");
  }

  gridParamT defaultGridParams() {
    gridParamT gp;

    gp.Pmin[0]=0.01;//*0+0.0937081;//0.2/10;
    gp.Pmax[0]=25;//15.9*0+16.9676;//20;

    gp.Umin[0]=-2.0;
    gp.Umax[0]=2.0;
    gp.Umin[1]=0;
    gp.Umax[1]=1.6;
    
    gp.Pres[0]=45*4;//256;
    gp.Ures[0]=40*4;//256;
    gp.Ures[1]=40*4;//128;
 
    gp.Pscale[0]=scaleTypeValT::LOGARITHMIC;//scaleTypeVal::LOGARITHMIC;
    gp.Uscale[0]=scaleTypeValT::REGULAR;
    gp.Uscale[1]=scaleTypeValT::QUADRATIC;
    return gp;
  }

  solverParamT defaultSolverParams() {
    solverParamT sp;
    sp.T0=0;    
    sp.dt=0.1;//0.05/10;    
    sp.Tmax=Tc*10;
    return sp;
  }

private:
  paramT p;
  std::map<std::string, double*> pMapD;
  std::map<std::string, char *> pMapC;

  double Tc;
  double rho0;
  double C1;
  double C2;
  double alpha;
  double sig0;

private:

  double computeRho(double R) {
    if (R>p.R0) return 0;
    return (3.*p.M)/(4.*M_PI*pow(p.R0,3));
  }

  void computeConstants() {
    rho0=(3*p.M)/(4*M_PI*pow(p.R0,3));
    Tc=M_PI*sqrt(pow(p.R0,3)/(2*p.G*p.M));

    // alpha=-2T/W
    // cf henon(64) / fujiwara (83) 
    // -->  1/(2*sig0*sig0) =5R/(2alpha)
    // for a=0.5 :
    sig0 = sqrt((p.alpha)/(5*p.R0));
    C1=rho0*pow(2*M_PI*sig0*sig0,-1.5);
    C2=2*sig0*sig0;
  }

  double computeF(double r, double u, double j)
  {
    if (r>p.R0) return 0;
    //return C1*exp(-(u*u+(j*j)/(r*r))/C2);
    return (exp(-pow(r,15)/pow(p.R0,15)) - exp(-1))/(1.-exp(-1))*C1*exp(-(u*u+(j*j)/(r*r))/C2);
  } 

  template <class DT>
  dataT valueAt(std::vector<double>& pos, dimTraits<1,2>)
  {
    return computeF(pos[0],pos[1],pos[2]);
  }
  
  template <class DT>
  gridT *generate(gridParamT &gp,solverParamT &sp, dimTraits<1,2>)
  {
    printf("Generating initial conditions (%s) ... ",getTypeStr().c_str());fflush(0);
    gridT *g=new gridT(gp);
    dataT *d=g->getDataPtr();    
    
    long r,u,j,i;
    double R,U,J;   
    double dJ,dU,dR;
        
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
   
#pragma omp parallel for shared(R_C,U_C,J_C,R_V,U_V,J_V,delta,d) private(J,U,R,j,u,r,i,dJ,dR,dU)
    //for (j=0;j<gp.Vres[1];j++)
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
		
		d[i]=computeF(R,U,J);
		/*
		d[i]=1./8. * (
		  computeF(R_V[r],U_V[u],J_V[j])+
		  computeF(R_V[r+1],U_V[u],J_V[j])+
		  computeF(R_V[r],U_V[u+1],J_V[j])+
		  computeF(R_V[r+1],U_V[u+1],J_V[j])+
		  computeF(R_V[r],U_V[u],J_V[j+1])+
		  computeF(R_V[r+1],U_V[u],J_V[j+1])+
		  computeF(R_V[r],U_V[u+1],J_V[j+1])+
		  computeF(R_V[r+1],U_V[u+1],J_V[j+1]));
		*/
	      }
		
	  }
      }

    double massBelowR0=computeRho(0)*4./3.*M_PI*pow(R_V[0],3);
    
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

   UniformDensitySphere():base() {
    registerParams();
    /*
    p.G=1;
    p.R0=2;
    p.M=1;
    p.alpha=0.5; // alpha=1 -> equilibrium (virial theorem satisfied)
    */
    computeConstants();
  }

  UniformDensitySphere(paramT &p_) : base(){
    registerParams();
    p=p_;
    computeConstants();
  }

  ~UniformDensitySphere() {
  }



private:
  template <class DT> 
  dataT valueAt(std::vector<double>& pos){
    return valueAt<DT>(pos,DT());
  }

  template <class DT>
  dataT valueAt(std::vector<double>& pos, DT){
    fprintf(stderr,"ERROR generating initial conditions :\n");
    fprintf(stderr,"  %s: not implemented in (%d+%d)D.\n",getTypeStr().c_str(),DT::P_DIMS,DT::U_DIMS);
    exit(0);
  }

  template <class DT> 
  gridT *generate(gridParamT &gp,solverParamT &sp){
    return generate<DT>(gp,sp,DT());
  }
  
  template <class DT>
  gridT *generate(gridParamT &gp,solverParamT &sp, DT){
    fprintf(stderr,"ERROR generating initial conditions :\n");
    fprintf(stderr,"  %s: not implemented in (%d+%d)D.\n",getTypeStr().c_str(),DT::P_DIMS,DT::U_DIMS);
    exit(0);
  }
};


#endif
  


