#ifndef __INIT_GEN_TEST_NN_HXX__
#define __INIT_GEN_TEST_NN_HXX__

#include <vector>
#include <map>
#include <string>
#include <assert.h>

//#include "quadrature21.hxx"
#include "global.h"
#include "initGen.hxx"

template <class gridType>
class testNN : public initGen<gridType> {

public:
  typedef initGen<gridType> bT; 

  static const int P_DIMS= bT::P_DIMS;
  static const int U_DIMS= bT::U_DIMS;
  static const int J_DIMS= bT::J_DIMS;
  static const int DIMS = bT::DIMS;

  static const int MAX_SPECIES=bT::MAX_SPECIES; 
  typedef typename bT::dataT dataT;
  typedef typename bT::gridT gridT;
  typedef typename bT::dimTr dimTr;
  typedef typename bT::scaleT scaleT;
  typedef typename bT::scaleTypeT scaleTypeT;
  typedef typename bT::scaleTypeV scaleTypeV;  
  typedef typename bT::gridParamT gridParamT;
  typedef typename bT::solverParamT solverParamT;
  typedef typename bT::gridSliceItT gridSliceItT;

  typedef typename bT::queryT queryT;
  
  struct TNN_params {
    double M[MAX_SPECIES];
    double R0[MAX_SPECIES];
    double alpha[MAX_SPECIES];
    double smoothingRadius[MAX_SPECIES];
  };

  typedef TNN_params paramT; 

public:

  void setupParameters(const paramsParser &parser,const solverParamT &sp_)
  {
    int i;
    //default parameters
    sp=sp_;

    for (i=0;i<sp.nSpecies;i++)
      {
	p.R0[i]=1;
	p.M[i]=1;
	p.alpha[i]=0.5;
	p.smoothingRadius[i]=0.5; 

	p.R0[i] = parser.get<>("R0",bT::parserCategory(),p.R0[i],i);
	p.M[i] = parser.get<>("M",bT::parserCategory(),p.M[i],i);
	p.alpha[i] = parser.get<>("alpha",bT::parserCategory(),p.alpha[i],i);
	p.smoothingRadius[i]=parser.get<>("smoothingRadius",bT::parserCategory(),p.smoothingRadius[i],i);
      }
  }

  std::string getTypeStr() {
    return std::string("TNN");
  }

  gridParamT defaultGridParams() {
    //gridParamT gp = bT::defaultGridParams();
    gridParamT gp = {};
    int i;
    for (i=0;i<P_DIMS;i++)
      {
	gp.Pmin[i]=-2;
	gp.Pmax[i]=2;
	gp.Pres[i]=128/2;
	gp.Pscale[i]=scaleTypeV::LINEAR;
      }

     for (i=0;i<U_DIMS;i++)
      {
	gp.Umin[i]=-2;
	gp.Umax[i]=2;
	gp.Ures[i]=32*2;
	gp.Uscale[i]=scaleTypeV::LINEAR;
      }

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

private:

  double computeF(const double *r,const double *u,int s=0)
  {
    //double res=cos(9*r[0]*sin(8*u[0]))*100+80*sin(15*(u[1]*r[0]-2*u[0]*r[1]));
    double rr=0;
    double uu=0;
    for (int i=0;i<P_DIMS;i++) rr+=r[i]*r[i];
    for (int i=0;i<U_DIMS;i++) uu+=u[i]*u[i];
    
    //double res=erfc(1000.*(sqrt(rr)-p.R0[s]))/2;
    //double res=(rr==0)?1:0;
    double res=exp(-rr/p.R0[s])*exp(-uu/(0.1))*erfc(5*(sqrt(rr)-p.R0[s]));

    return res;
  } 

  template <class DT>
  dataT valueAt(const std::vector<double>& pos, int s, dimTraits<P_DIMS,P_DIMS>)
  {
    return 0;
  }
  
  template <class DT>
  gridT *generate(const gridParamT &gp, dimTraits<P_DIMS,P_DIMS>)
  {
    printf("Generating initial conditions (%s[%d]) ... ",getTypeStr().c_str(),sp.nSpecies);fflush(0);
    
    gridT *g=new gridT();
    g->init(gp);
    long s,j;
    
    const long nSlices=g->nSlices();
    //printf("nSlice = %ld\n",nSlices);
    for (s=0;s<sp.nSpecies;s++)
      {
#pragma omp parallel for 
	for (j=0;j<nSlices;j++)
	  {
	    double pos[DIMS];
	    const oneField_iterator<gridSliceItT> sit_end(g->slice_end(j),s);
	    oneField_iterator<gridSliceItT> sit(g->slice_begin(j),s);
	    
	    //sit.print();
	    //sit.get_C(pos);
	    //sit.get_C(pos);
	    
	    for (;sit!=sit_end;++sit) 
	      {		
		//sit.print();
		sit.C(pos);
		(*sit)=computeF(pos,&pos[P_DIMS],s);
	      }
	    

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
      default:
	printf("WARNING in initGen_TNN::query : invalid query (%d)!\n",(int)what);
	return -1;
      }
  }

  testNN():bT() {
  }

  ~testNN() {
  }
 
  /* NO NEED TO MODIFY BELOW HERE */
  /* COPY/PASTE to create a new class */
  
  dataT valueAt(const std::vector<double> &pos, int s) {
    assert(bT::initialized);
    return valueAt<dimTr>(pos,s,dimTr());
  }

   gridT *generate(const gridParamT &gp) {
    assert(bT::initialized);
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
  


