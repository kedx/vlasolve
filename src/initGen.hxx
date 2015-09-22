#ifndef __INITGEN_BASE__HXX__
#define __INITGEN_BASE__HXX__

#include <string>
#include "solverParam.hxx"

//#include "solverParamsBase.hxx"
#include "dimTraits.hxx"
#include "paramsParser.hxx"

template <class gridType>
class initGen {

public:  

  static const int P_DIMS= gridType::P_DIMS;
  static const int U_DIMS= gridType::U_DIMS;
  static const int J_DIMS= U_DIMS-P_DIMS;
  static const int DIMS = U_DIMS+P_DIMS;

  static const int MAX_SPECIES=5; 

  struct queryV {
    enum type {RENORMALIZE=0, KERNELMASS=1};
  };

  typedef typename queryV::type queryT;

  static const std::string parserCategory() {return "init";}
  typedef solverParam solverParamT;
  
  typedef gridType gridT;
  typedef typename gridT::dimTr dimTr;
  typedef typename gridT::value_type dataT;  
  typedef typename gridT::scaleT scaleT;
  typedef typename gridT::scaleTypeT scaleTypeT;
  typedef typename gridT::scaleTypeV scaleTypeV;
  typedef typename gridT::paramT gridParamT; 
  typedef typename gridT::sliceIterator gridSliceItT;
  
  virtual std::string getTypeStr()=0;
  virtual void setupParameters(const paramsParser &parser, const solverParamT &sp_)=0;
  virtual gridParamT defaultGridParams()=0;
  /*
  {
    int i;
    for (i=0;i<P_DIMS;i++)
      {
	gp.Pmin[i]=0;
	gp.Pmax[i]=1;
	gp.Pres[i]=128;
	gp.Pscale[i]=scaleTypeV::LINEAR;
      }

    for (i=0;i<U_DIMS;i++)
      {
	gp.Umin[i]=-1;
	gp.Umax[i]=1;
	gp.Ures[i]=128;
	gp.Uscale[i]=scaleTypeV::LINEAR;
      }
    gp.haveParentGrid=0;
    gp.parent_P_DIMS=P_DIMS;
    gp.parent_U_DIMS=U_DIMS;  
    for (i=0;i<DIMS;i++) gp.parent_which[DIMS]=i;

  }
  */
  virtual solverParamT defaultSolverParams()=0;

  virtual dataT valueAt(const std::vector<double> &pos, int s=0)=0;
  virtual gridT *generate(const gridParamT &gp)=0;

  //virtual double renormalize() {return -1;}
  //virtual double getValue(const std::string &what, double param)
  virtual double getMass(int s=-1)=0;

  virtual double query(const queryT &what, double param=0)
  {
    printf("WARNING in initGen_UDF::query : invalid query (=%d)!\n",(int)what);
    return -1;
  }

  initGen():initialized(false),nSpecies(0)
  {}

  virtual ~initGen()
  {}

  gridT *generate() {
    gridParamT gp=defaultGridParams();
    return generate(gp);
  }

  void init(const paramsParser &parser, const solverParamT &sp_) {
    nSpecies = sp_.nSpecies;//parser.get<>("nSpecies",parserCategory(),1);
    setupParameters(parser,sp_);
    initialized=true;
  }

  int getNSpecies() {return nSpecies;}

protected : 
  int nSpecies;
  bool initialized;
};


#endif
 
