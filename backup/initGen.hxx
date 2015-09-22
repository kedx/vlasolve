#ifndef __INITGEN_BASE__HXX__
#define __INITGEN_BASE__HXX__

#include <string>

#include "solverParams.hxx"
#include "dimTraits.hxx"

template <class gridType>
class initGen {

public: 
 
  typedef solverParams solverParamT;  
  typedef gridType gridT;
  typedef typename gridT::dimTr dimTr;
  typedef typename gridT::dataT dataT;  
  typedef typename gridT::scaleT scaleT;
  typedef typename gridT::scaleType scaleType;
  typedef typename gridT::scaleTypeVal scaleTypeVal;
  typedef typename gridT::paramT gridParamT;  
  
  virtual std::string getTypeStr()=0;
  virtual void setParams(std::string &str)=0;
  virtual gridParamT defaultGridParams()=0;
  virtual solverParamT defaultSolverParams()=0;

  virtual dataT valueAt(std::vector<double> &pos)=0;
  virtual gridT *generate(solverParamT &sp)=0;
  virtual gridT *generate(gridParamT &gp, solverParamT &sp)=0;
  virtual gridT *generate(gridParamT &gp)=0;
  virtual gridT *generate()=0;
  
  initGen()
  {}

  virtual ~initGen()
  {}
};


#endif
 
