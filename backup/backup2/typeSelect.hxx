#ifndef __TYPE_SELECT_HXX__
#define __TYPE_SELECT_HXX__

#include <map>
#include <string>

template <class S>
struct typeSelect {
private:
  typedef typename S:: type T;
  typedef typename std::map<std::string,T> mapT;
  typedef typename mapT::iterator mapItT;
  typedef typename std::map<T,std::string> mapInvT;
  typedef typename mapInvT::iterator mapInvItT;
  
  mapT m_;
  mapInvT minv_;

public:
  typeSelect() {}
  virtual ~typeSelect() {}

  void insert(const std::string &str,const T &val)
  {
    m_.insert(std::make_pair(str,val));
    minv_.insert(std::make_pair(val,str));
  }

  T getVal(const std::string &str, bool errorIfUndefined=true) 
  { 
    mapItT it = m_.find(str);
    if (it==m_.end()) 
      {
	if (errorIfUndefined)
	  {
	    fprintf(stderr,"ERROR: unknown %s value: '%s'\n",name().c_str(),str.c_str());
	    exit(-1);
	  }
	return S::UNDEFINED;
      }
    else return it->second;
  }

  std::string getString(const T &val, bool errorIfUndefined=true) 
  { 
    mapInvItT it = minv_.find(val);
    if (it==minv_.end())
      {
	if (errorIfUndefined)
	  {
	    fprintf(stderr,"ERROR in typeselect::getString : no string associated to type '%s'\n",name().c_str());
	    exit(-1);
	  }
	return "undefined";
      }
    else return it->second;
  }

  virtual std::string name()=0;

};

#endif
