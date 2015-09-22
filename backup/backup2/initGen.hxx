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

  struct queryV {
    enum type {RENORMALIZE=0, KERNELMASS=1};
  };

  typedef typename queryV::type queryT;

  static const std::string parserCategory() {return "init";}
  typedef solverParam solverParamT;
  
  typedef gridType gridT;
  typedef typename gridT::dimTr dimTr;
  typedef typename gridT::dataT dataT;  
  typedef typename gridT::scaleT scaleT;
  typedef typename gridT::scaleTypeT scaleTypeT;
  typedef typename gridT::scaleTypeV scaleTypeV;
  typedef typename gridT::paramT gridParamT; 
  typedef typename gridT::sliceIterator gridSliceItT;
  
  virtual std::string getTypeStr()=0;
  virtual void setupParameters(const paramsParser &parser, const solverParamT &sp_)=0;
  virtual gridParamT defaultGridParams()=0;
  virtual solverParamT defaultSolverParams()=0;

  virtual dataT valueAt(const std::vector<double> &pos)=0;
  virtual gridT *generate(const gridParamT &gp)=0;

  //virtual double renormalize() {return -1;}
  //virtual double getValue(const std::string &what, double param)
  virtual double query(const queryT &what, double param=0)
  {
    printf("WARNING in initGen_UDF::query : invalid query (=%d)!\n",(int)what);
    return -1;
  }

  initGen():initialized(false)
  {}

  virtual ~initGen()
  {}

  gridT *generate() {
    gridParamT gp=defaultGridParams();
    return generate(gp);
  }

  void init(const paramsParser &parser, const solverParamT &sp_) {
    setupParameters(parser,sp_);
    initialized=true;
  }

protected : 

  bool initialized;
};


#endif
 
