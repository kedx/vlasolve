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
  typedef initGen<gridType> baseT; 
  typedef typename baseT::dataT dataT;
  typedef typename baseT::gridT gridT;
  typedef typename baseT::dimTr dimTr;
  typedef typename baseT::scaleT scaleT;
  typedef typename baseT::scaleTypeT scaleTypeT;
  typedef typename baseT::scaleTypeV scaleTypeV;  
  typedef typename baseT::gridParamT gridParamT;
  typedef typename baseT::solverParamT solverParamT;
  typedef typename baseT::gridSliceItT gridSliceItT;

  typedef typename baseT::queryT queryT;

  struct Hernquist_params {
    double M;
    double R0;
    double Ra;
  };

  typedef Hernquist_params paramT; 


void setupParameters(const paramsParser &parser,const solverParamT &sp_) {
  sp=sp_;
   
  p.R0=1;
  p.Ra=-1;// anisotropy radius, set to -1 for no anisotropy
  p.M=1;

  p.R0 = parser.get<>("R0",baseT::parserCategory(),p.R0);
  p.M = parser.get<>("M",baseT::parserCategory(),p.M);
  p.Ra = parser.get<>("Ra",baseT::parserCategory(),p.Ra);

  Vg=sqrt((sp.G*p.M)/(p.R0));
}

public:

  std::string getTypeStr() {
    return std::string("Hernquist");
  }

  gridParamT defaultGridParams() {
    gridParamT gp;

    gp.Pmin[0]=0.01;
    gp.Pmax[0]=25*10;

    gp.Umin[0]=-1.45;
    gp.Umax[0]=1.45;
    gp.Umin[1]=1.E-4;
    gp.Umax[1]=3.24*4;
    
    gp.Pres[0]=45*2;
    gp.Ures[0]=40*2;
    gp.Ures[1]=20*4;//20*2*5; 

    gp.Pscale[0]=scaleTypeV::LOGARITHMIC;
    gp.Uscale[0]=scaleTypeV::LINEAR;
    gp.Uscale[1]=scaleTypeV::LOGARITHMIC;
    return gp;
  }

  solverParamT defaultSolverParams() {
    solverParamT sp;   
    sp.G=1;
    sp.T0=0;
    sp.Tmax=62; // no Tmax
    sp.dt=0.025;
    return sp;
  }

private:
  paramT p;
  solverParamT sp;
  double Vg;
  
private:

  double computePhi(double R) {
    return -(sp.G*p.M)/(p.R0+R);
  }

  double computeRho(double R) {
    return (p.M*p.R0)/(2.*M_PI*R)*pow(p.R0+R,-3);
  }

   double kernelMass(double r0)
  {
    double result=0;
    long N=10000;
    double dr=r0/N;
    long i;

    for (i=0;i<N;i++)
      {
	double R=i*dr;
	result += 4./3.*M_PI*((R+dr)*(R+dr)*(R+dr) - R*R*R)*computeRho(R+0.5*dr);
      }
		  
    return result;
  }

  std::vector<double> computePhi(const std::vector<double> &R_) 
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

    double q=sqrt(-(p.R0)/(sp.G*p.M) * eps); 
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
  dataT valueAt(const std::vector<double>& pos, dimTraits<1,2>)
  {
    double phi=computePhi(pos[0]);
    return computeF(phi,pos[0],pos[1],pos[2]);
  }
  
  template <class DT>
  gridT *generate(const gridParamT &gp, dimTraits<1,2>)
  {
    printf("Generating initial conditions (%s) ... ",getTypeStr().c_str());fflush(0); 
    gridT *g=new gridT();
    g->init(gp);
    long j;
    std::vector<double> Phi=computePhi(g->getValCoord21_R());

    const long nSlices=g->nSlices();
    #pragma omp parallel for 
    for (j=0;j<nSlices;j++)
      {
	const gridSliceItT sit_end=g->slice_end(j);
	gridSliceItT sit=g->slice_begin(j);
	double J=sit.get_J();
	for (;sit!=sit_end;++sit) 
	  {
	    (*sit)=computeF(Phi[sit.get_r()],sit.get_R(),sit.get_U(),J);		
	  }
      }
    
    /*
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
    */
    printf("done.\n");
    return g;
  }

public:

  virtual double query(const queryT &what, double param=0)
  {
    switch (what)
      {
      case baseT::queryV::RENORMALIZE: 
	return p.M;
      case baseT::queryV::KERNELMASS:
	return kernelMass(param);
      default:
	printf("WARNING in initGen_UDF::query : invalid query (%d)!\n",(int)what);
	return -1;
      }
  }

  Hernquist():baseT() {
  }

  ~Hernquist() {
  }

  /* NO NEED TO MODIFY BELOW HERE */
  /* COPY/PASTE to create a new class */

  dataT valueAt(const std::vector<double> &pos) {
    assert(baseT::initialized);
    return valueAt<dimTr>(pos,dimTr());
  }

  gridT *generate(const gridParamT &gp) {
    assert(baseT::initialized);
    return generate<dimTr>(gp,dimTr());
  }
 
private:
  template <class DT>
  dataT valueAt(const std::vector<double>& pos, DT){
    fprintf(stderr,"ERROR generating initial conditions :\n");
    fprintf(stderr,"  %s: not implemented in (%d+%d)D.\n",getTypeStr().c_str(),DT::P_DIMS,DT::U_DIMS);
    exit(0);
  }
  
  template <class DT>
  gridT *generate(const gridParamT &gp, DT){
    fprintf(stderr,"ERROR generating initial conditions :\n");
    fprintf(stderr,"  %s: not implemented in (%d+%d)D.\n",getTypeStr().c_str(),DT::P_DIMS,DT::U_DIMS);
    exit(0);
  }
};


#endif
  


