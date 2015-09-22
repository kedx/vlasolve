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
  static const int MAX_SPECIES=baseT::MAX_SPECIES; 
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
    double M[MAX_SPECIES];
    double R0[MAX_SPECIES];
    double alpha[MAX_SPECIES];
    double smoothingRadius[MAX_SPECIES];
    double gap[2][MAX_SPECIES];
  };

  typedef UDS_params paramT; 

public:

  void setupParameters(const paramsParser &parser,const solverParamT &sp_)
  {
    int i;
    //default parameters
    sp=sp_;

    for (i=0;i<sp.nSpecies;i++)
      {
	p.R0[i]=2;
	p.M[i]=1;
	p.alpha[i]=0.5; // alpha=1 -> equilibrium (virial theorem satisfied)
	p.smoothingRadius[i]=0.5; // in pixels
	p.gap[0][i]=p.gap[1][i]=0;

	p.R0[i] = parser.get<>("R0",baseT::parserCategory(),p.R0[i],i);
	p.M[i] = parser.get<>("M",baseT::parserCategory(),p.M[i],i);
	p.alpha[i] = parser.get<>("alpha",baseT::parserCategory(),p.alpha[i],i);
	p.smoothingRadius[i]=parser.get<>("smoothingRadius",baseT::parserCategory(),p.smoothingRadius[i],i);
	p.gap[0][i]=parser.get<>("gap",baseT::parserCategory(),p.gap[0][i],2*i);
	p.gap[1][i]=parser.get<>("gap",baseT::parserCategory(),p.gap[1][i],2*i+1);


	rho0[i]=(3*p.M[i])/(4*M_PI*pow(p.R0[i],3));
	Tc[i]=M_PI*sqrt(pow(p.R0[i],3)/(2*sp.G*p.M[i]));
   
	// alpha=-2T/W
	// cf henon(64) / fujiwara (83) 
	// -->  1/(2*sig0*sig0) =5R/(2alpha)
	// for a=0.5 :
	sig0[i] = sqrt((p.alpha[i])/(5*p.R0[i]));
	C1[i]=rho0[i]*pow(2*M_PI*sig0[i]*sig0[i],-1.5);
	C2[i]=2*sig0[i]*sig0[i];
      }
  }

  std::string getTypeStr() {
    return std::string("UDF");
  }

  gridParamT defaultGridParams() {
    gridParamT gp = {};

    gp.Pmin[0]=0.01;//*0+0.0937081;//0.2/10;
    gp.Pmax[0]=25;//15.9*0+16.9676;//20;

    gp.Umin[0]=-2.0;
    gp.Umax[0]=2.0;
    gp.Umin[1]=0;
    gp.Umax[1]=1.6;
    
    gp.Pres[0]=200;//256;
    gp.Ures[0]=200;//256;
    gp.Ures[1]=200;//128;
 
    gp.Pscale[0]=scaleTypeV::LOGARITHMIC;
    gp.Uscale[0]=scaleTypeV::LINEAR;
    gp.Uscale[1]=scaleTypeV::QUADRATIC;

    return gp;
  }

  solverParamT defaultSolverParams() {
    solverParamT sp_;
    sp_.G=1;
    sp_.T0=0;    
    sp_.dt=0.005;
    sp_.Tmax=60.;
    sp_.nSpecies=1;
    return sp_;
  }

private:
  paramT p;
  solverParamT sp;

  double Tc[MAX_SPECIES];
  double rho0[MAX_SPECIES];
  double C1[MAX_SPECIES];
  double C2[MAX_SPECIES];
  double alpha[MAX_SPECIES];
  double sig0[MAX_SPECIES];

private:

  double computeRho(double r, int s=0) {
    double sm=0.5*(1.+erf((p.R0[s]-r)/p.smoothingRadius[s]));
    return sm*(3.*p.M[s])/(4.*M_PI*pow(p.R0[s],3));
    //if (R>p.R0) return 0;
    //return (3.*p.M)/(4.*M_PI*pow(p.R0,3));
  }

  double kernelMass(double r0, int s=0)
  {
    double result=0;
    long N=100;
    double dr=r0/N;
    long i;

    for (i=0;i<N;i++)
      {
	double R=i*dr;
	result += 4./3.*M_PI*((R+dr)*(R+dr)*(R+dr) - R*R*R)*computeRho(R+0.5*dr,s);
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
  double computeF(const double r,const double u,const double j, int s=0)
  {
    double sm=0.5*(1.+erf((p.R0[s]-r)/p.smoothingRadius[s]));
    double res=sm*C1[s]*exp(-(u*u+(j*j)/(r*r))/C2[s]);
    
    if (p.gap[1][s]>0)
      {	
	double sm0=0.5*(1.+erf((p.gap[0][s]-r)/0.15));
	double sm1=0.5*(1.+erf(-(p.gap[1][s]-r)/0.15));
	//printf("sm0+sm1 = %lg %lg<[%lg,%lg]\n",sm0+sm1,r,p.gap[0][s],p.gap[1][s]);
	res *= (sm0+sm1);
      }

    if (res<1.E-200) res=1.E-200;
    return res;
  } 

  template <class DT>
  dataT valueAt(const std::vector<double>& pos, int s, dimTraits<1,2>)
  {
    return computeF(pos[0],pos[1],pos[2],s);
  }
  
  template <class DT>
  gridT *generate(const gridParamT &gp, dimTraits<1,2>)
  {
    printf("Generating initial conditions (%s[%d]) ... ",getTypeStr().c_str(),sp.nSpecies);fflush(0);
    
    gridT *g=new gridT();
    g->init(gp);
    long s,j;
    
    const long nSlices=g->nSlices();

    for (s=0;s<sp.nSpecies;s++)
      {
#pragma omp parallel for 
	for (j=0;j<nSlices;j++)
	  {	   
	    const oneField_iterator<gridSliceItT> sit_end(g->slice_end(j),s);
	    oneField_iterator<gridSliceItT> sit(g->slice_begin(j),s);
	  
	    double J=sit.get_J();
	    for (;sit!=sit_end;++sit) 
	      (*sit)=computeF(sit.get_R(),sit.get_U(),J,s);
	  }
      }
  
    printf("done.\n");
    return g;
  }

public:

  virtual double getMass(int s)
  {
    if (s<0)
      {
	double mtot=0;
	for (int i=0;i<sp.nSpecies;i++) 
	  mtot+=p.M[i];
	return mtot;
      }
    else return p.M[s];    
  }
  
  virtual double query(const queryT &what, double param=0)
  {
    switch (what)
      {
      case baseT::queryV::RENORMALIZE:
	{
	  double mtot=0;
	  for (int i=0;i<sp.nSpecies;i++) 
	    mtot+=p.M[i];
	  return mtot;
	}
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
  
  dataT valueAt(const std::vector<double> &pos, int s) {
    assert(baseT::initialized);
    return valueAt<dimTr>(pos,s,dimTr());
  }

   gridT *generate(const gridParamT &gp) {
    assert(baseT::initialized);
    return generate<dimTr>(gp,dimTr());
  }
 
private:

  template <class DT>
  dataT valueAt(const std::vector<double>& pos, int s, DT){
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
  


