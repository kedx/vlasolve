#ifndef __INIT_GEN_UNIFORM_DENSITY_SPHERE_HXX__
#define __INIT_GEN_UNIFORM_DENSITY_SPHERE_HXX__

#include <vector>
#include <map>
#include <string>
#include <assert.h>

//#include "quadrature21.hxx"
#include "global.h"
#include "initGen.hxx"


template <class gridType>
class UniformDensitySphere : public initGen<gridType> {

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
  
  
  struct UDS_params {
    double M;
    double R0;
    double alpha;
    double smoothingRadius;
  };

  typedef UDS_params paramT; 

public:

  void setupParameters(const paramsParser &parser,const solverParamT &sp_)
  {
    //default parameters
    sp=sp_;

    p.R0=2;
    p.M=1;
    p.alpha=0.5; // alpha=1 -> equilibrium (virial theorem satisfied)
    p.smoothingRadius=0.5; // in pixels

    p.R0 = parser.get<>("R0",baseT::parserCategory(),p.R0);
    p.M = parser.get<>("M",baseT::parserCategory(),p.M);
    p.alpha = parser.get<>("alpha",baseT::parserCategory(),p.alpha);
    p.smoothingRadius=parser.get<>("smoothingRadius",baseT::parserCategory(),p.smoothingRadius);

    rho0=(3*p.M)/(4*M_PI*pow(p.R0,3));
    Tc=M_PI*sqrt(pow(p.R0,3)/(2*sp.G*p.M));

    // alpha=-2T/W
    // cf henon(64) / fujiwara (83) 
    // -->  1/(2*sig0*sig0) =5R/(2alpha)
    // for a=0.5 :
    sig0 = sqrt((p.alpha)/(5*p.R0));
    C1=rho0*pow(2*M_PI*sig0*sig0,-1.5);
    C2=2*sig0*sig0;
  }

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
    
    gp.Pres[0]=200;//256;
    gp.Ures[0]=200;//256;
    gp.Ures[1]=200;//128;
 
    gp.Pscale[0]=scaleTypeV::LOGARITHMIC;//scaleTypeVal::LOGARITHMIC;
    gp.Uscale[0]=scaleTypeV::LINEAR;
    gp.Uscale[1]=scaleTypeV::QUADRATIC;
    return gp;
  }

  solverParamT defaultSolverParams() {
    solverParamT sp_;
    sp_.G=1;
    sp_.T0=0;    
    sp_.dt=0.005;
    sp_.Tmax=Tc*10;
    return sp_;
  }

private:
  paramT p;
  solverParamT sp;

  double Tc;
  double rho0;
  double C1;
  double C2;
  double alpha;
  double sig0;

private:

  double computeRho(double r) {
    double sm=0.5*(1.+erf((p.R0-r)/p.smoothingRadius));
    return sm*(3.*p.M)/(4.*M_PI*pow(p.R0,3));
    //if (R>p.R0) return 0;
    //return (3.*p.M)/(4.*M_PI*pow(p.R0,3));
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

  /*  
  double smoothingF(double r, double u, double j)
  {
    double h=0.5*C1*exp(-(u*u+(j*j)/(p.R0*p.R0))/C2);
    double sgn=(r<=p.R0)?-1.0:1.0;
    return h*sgn*exp(-fabs(r-p.R0)*1./p.smoothingRadius);
  }
  */
  double computeF(const double r,const double u,const double j)
  {
    double sm=0.5*(1.+erf((p.R0-r)/p.smoothingRadius));
    double res=sm*C1*exp(-(u*u+(j*j)/(r*r))/C2);
    if (res<1.E-200) res=1.E-200;
    return res;
  } 

  template <class DT>
  dataT valueAt(const std::vector<double>& pos, dimTraits<1,2>)
  {
    return computeF(pos[0],pos[1],pos[2]);
  }
  
  template <class DT>
  gridT *generate(const gridParamT &gp, dimTraits<1,2>)
  {
    printf("Generating initial conditions (%s) ... ",getTypeStr().c_str());fflush(0);
    
    gridT *g=new gridT();
    g->init(gp);
    long j;
    
    const long nSlices=g->nSlices();
    #pragma omp parallel for 
    for (j=0;j<nSlices;j++)
      {
	const gridSliceItT sit_end=g->slice_end(j);
	gridSliceItT sit=g->slice_begin(j);
	double J=sit.get_J();
	for (;sit!=sit_end;++sit) 
	  (*sit)=computeF(sit.get_R(),sit.get_U(),J);
      }
  
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

  UniformDensitySphere():baseT() {
  }

  ~UniformDensitySphere() {
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
  


