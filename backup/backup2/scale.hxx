#ifndef __SCALE_HXX__
#define __SCALE_HXX__

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"

#include <algorithm>
#include <string>

#include "valLocationType.hxx"

template <typename T>
struct scale {

  struct scaleTypeV {
    enum type {LINEAR=0,LOGARITHMIC=1,QUADRATIC=2,CUBIC=3,QUARTIC=4,UNDEFINED=-1};
  };
  typedef typename scaleTypeV::type scaleTypeT;

  struct scaleTypeSelect : public typeSelect<scaleTypeV> {
    scaleTypeSelect()
    {
      insert("linear",scaleTypeV::LINEAR);
      insert("logarithmic",scaleTypeV::LOGARITHMIC);
      insert("quadratic",scaleTypeV::QUADRATIC);
      insert("cubic",scaleTypeV::CUBIC);
      insert("quartic",scaleTypeV::QUARTIC);
    }
    std::string name() {return "scale type";}
  };

  typedef valLocationVal valLocationV;
  typedef valLocationType valLocationT;  
  typedef T dataT;

  /*
  static scaleTypeT str2Type(std::string str)
  {
    
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);

    if (str=="LINEAR") return scaleTypeV::LINEAR;
    if (str=="LOGARITHMIC") return scaleTypeV::LOGARITHMIC;
    if (str=="QUADRATIC") return scaleTypeV::QUADRATIC;
    if (str=="CUBIC") return scaleTypeV::CUBIC;
    if (str=="QUARTIC") return scaleTypeV::QUARTIC;
    
    fprintf(stderr,"ERROR in scale::str2Type : unknow type '%s'.\n",str.c_str());
    exit(-1);
  }

  static std::string type2Str(scaleTypeT type)
  {
    if (type==scaleTypeV::LINEAR) return "LINEAR";
    if (type==scaleTypeV::LOGARITHMIC) return "LOGARITHMIC";
    if (type==scaleTypeV::QUADRATIC) return "QUADRATIC";
    if (type==scaleTypeV::CUBIC) return "CUBIC";
    if (type==scaleTypeV::QUARTIC) return "QUARTIC";
    

    fprintf(stderr,"ERROR in scale::str2Type : invalid type .\n");
    exit(-1);
  }
  
  static std::string valLocation2Str(valLocationT vl)
  {
    if (vl==valLocationV::CELL) return "C";
    if (vl==valLocationV::VERTEX) return "V";
    
    fprintf(stderr,"ERROR in valLocation2Str : invalid type .\n");
    exit(-1);
  }
  */

  static T valueAt(double start, double stop, long N, long index, scaleTypeT sct=scaleTypeV::LINEAR, valLocationT ltype=valLocationV::VERTEX)
  {
    double val;

    if (sct==scaleTypeV::LINEAR) { 
      double dx=(stop-start)/(N);
      val =start+dx*index;
    }
    else if (sct==scaleTypeV::LOGARITHMIC) {
      if (start*stop<=0) {
	fprintf(stderr,"error: cannot generate logscale from a 0-crossing range.\n");
	exit(0);
      }
      double s=(start<0)?-1:1;
      val=valueAt(log(fabs(start)),log(fabs(stop)),N,index,scaleTypeV::LINEAR,valLocationV::VERTEX);
      val = s*exp(val);
    }
    else if ((int)sct>=2) {
      int pp=(int)sct;
      val=valueAt(pow(start,1./pp),pow(stop,1./pp),N,index,scaleTypeV::LINEAR,valLocationV::VERTEX);
      val=pow(val,pp);
    }

    if (ltype==valLocationV::CELL) {
      double val2 = valueAt(start,stop,N,index+1,sct,valLocationV::VERTEX);
      val = (val+val2)*0.5;
    }
    else
      {
	if (index==0) return start;
	if (index==N) return stop;
      }

    return val;
  }

  static std::vector<T> genScale(double start, double stop, long N, scaleTypeT sct=scaleTypeV::LINEAR, valLocationT ltype=valLocationV::VERTEX)
  {
    std::vector<double> scale;
    long i,j;
     
    if (sct==scaleTypeV::LINEAR) { 
      double dx=(stop-start)/(N);
      scale.resize(N+1);
      for (i=0;i<N+1;i++) scale[i]=start+dx*i;
    } 
    else if (sct==scaleTypeV::LOGARITHMIC) {
      if (start*stop<=0) {
	fprintf(stderr,"error: cannot generate logscale from a 0-crossing range.\n");
	exit(0);
      }
      double s=(start<0)?-1:1;
      scale=genScale(log(fabs(start)),log(fabs(stop)),N,scaleTypeV::LINEAR,valLocationV::VERTEX);
      for (i=0;i<N+1;i++) scale[i]=s*exp(scale[i]);
    }
    else if ((int)sct>=2) {
      int pp=(int)sct;
      scale=genScale(pow(start,1./pp),pow(stop,1./pp),N,scaleTypeV::LINEAR,valLocationV::VERTEX);
      for (i=0;i<N+1;i++) scale[i]=pow(scale[i],pp);
    }
  
    if (ltype==valLocationV::CELL) {
      int j;
      for (j=0;j<scale.size()-1;j++)  scale[j]=(scale[j]+scale[j+1])*0.5;
      scale.resize(scale.size()-1);
    }
    else
      {
	scale.front()=start;
	scale.back()=stop;	
      }

    return scale;
  }

  template <typename dT>
  static bool rescale(dT &x,scaleTypeT sct)
  {
    if (sct==scaleTypeV::LOGARITHMIC) x=log((x));
    else if (sct>=2) x=pow(x,1./sct);
    else return false;
    
    return true;
  }

  template <typename dT>
  static void rescale(dT &x,dT &y,scaleTypeT xsct,scaleTypeT ysct)
  {
    rescale(x,xsct);
    rescale(y,ysct);
  }

  template <typename dT>
  static void rescale(dT &x,dT &y,dT &z,scaleTypeT xsct,scaleTypeT ysct,scaleTypeT zsct)
  {
    rescale(x,xsct);
    rescale(y,ysct);
    rescale(z,zsct);
  }

  template <typename dT>
  static void unrescale(dT &x,scaleTypeT sct)
  {
    if (sct==scaleTypeV::LOGARITHMIC) x=exp(x);
    else if (sct>=2) x=pow(x,sct);
  }

  template <typename dT>
  static void unrescale(dT &x,dT &y,scaleTypeT xsct,scaleTypeT ysct)
  {
    unrescale(x,xsct);
    unrescale(y,ysct);
  }

  template <typename dT>
  static void unrescale(dT &x,dT &y,dT &z,scaleTypeT xsct,scaleTypeT ysct,scaleTypeT zsct)
  {
    unrescale(x,xsct);
    unrescale(y,ysct);
    unrescale(z,zsct);
  }
};

#endif
