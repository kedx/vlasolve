#ifndef __INIT_GEN_ALL_HXX__
#define __INIT_GEN_ALL_HXX__

//#include "setup.hxx"
#include <vector>
#include "initGen.hxx"
#include "initGen_UniformDensitySphere.hxx"
#include "initGen_Plummer.hxx"
#include "initGen_Hernquist.hxx"

template <class gridT>
class initGenAll {
public:
typedef initGen<gridT> initGenT;

private:
  typedef std::vector<initGenT *> listT;
  typedef typename listT::iterator listItT;

  static void registerTypes(std::vector<initGenT *> &initAll)
  {
    initAll.push_back((initGenT *)new UniformDensitySphere<gridT>());
    initAll.push_back((initGenT *)new Plummer<gridT>());
    initAll.push_back((initGenT *)new Hernquist<gridT>());
  }

public:  
  
  static initGenT *get(std::string &str)
  {
    std::vector<initGenT *> initAll;
    initGenT *result=NULL;
    registerTypes(initAll);    
    for (listItT it=initAll.begin();it!=initAll.end();it++)
      {
	if ((*it)->getTypeStr() == str) result=*it;
	else delete *it;
      } 
    
    if (result==NULL) {
      fprintf(stderr,"ERROR: unknown initial conditions type %s\n.",str.c_str());
      exit(0);
    }
   
    return result;
  }

  static initGenT *get(const char *strc)
  {
    std::string str(strc);
    return get(str);
  }
};

#endif
