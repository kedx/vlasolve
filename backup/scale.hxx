#ifndef __SCALE_HXX__
#define __SCALE_HXX__

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"

template <typename T>
struct scale {
  struct scaleTypeVal {
    enum type {REGULAR=0,LOGARITHMIC=1,QUADRATIC=2,CUBIC=3,QUARTIC=4};
  };
  struct valLocationTypeVal {
    enum type {PIXEL, VERTEX};
  };

  typedef typename scaleTypeVal::type scaleType;
  typedef typename valLocationTypeVal::type valLocationType;
  typedef T dataT;

  static std::vector<T> genScale(double start, double stop, long N, scaleType sct=scaleTypeVal::REGULAR, valLocationType ltype=valLocationTypeVal::VERTEX)
  {
    std::vector<double> scale;
    long i,j;
    
   
    if (sct==scaleTypeVal::REGULAR) { 
      double dx=(stop-start)/(N-1);
      scale.resize(N);
      for (i=0;i<N;i++) scale[i]=start+dx*i;
    } 
    else if (sct==scaleTypeVal::LOGARITHMIC) {
      if (start*stop<0) {
	fprintf(stderr,"error: cannot generate logscale from a 0-crossing range.\n");
	exit(0);
      }
      double s=(start<0)?-1:1;
      scale=genScale(log(fabs(start)),log(fabs(stop)),N);
      for (i=0;i<N;i++) scale[i]=s*exp(scale[i]);
    }
    else if ((int)sct>=2) {
      int pp=(int)sct;
      //scale.resize(N);
      scale=genScale(pow(start,1./pp),pow(stop,1./pp),N);
      for (i=0;i<N;i++) scale[i]=pow(scale[i],pp);
      /*
      for (i=0;i<N;i++) printf("%g -",scale[i]-(start + (double)pow(i,pp)/(double)pow((N-1),pp) * (stop-start)));printf("\n\n");fflush(0);
      */
    }
    /*
    else if (sct==scaleTypeVal::QUADRATIC) {
      scale.resize(N);
      for (i=0;i<N;i++) scale[i]=start + (double)(i*i)/(double)((N-1)*(N-1)) * (stop-start);
    }
    else if (sct==scaleTypeVal::CUBIC) {
      scale.resize(N);
      for (i=0;i<N;i++) scale[i]=start + (double)(i*i*i)/(double)((N-1)*(N-1)*(N-1)) * (stop-start);
    }
    else if (sct==scaleTypeVal::QUARTIC) {
      scale.resize(N);
      for (i=0;i<N;i++) scale[i]=start + (double)(i*i*i*i)/(double)((N-1)*(N-1)*(N-1)*(N-1)) * (stop-start);
    }
    */
    if (ltype==valLocationTypeVal::PIXEL) {
      int j;
      for (j=0;j<scale.size()-1;j++)  scale[j]=(scale[j]+scale[j+1])*0.5;
      scale.resize(scale.size()-1);
    }

    return scale;
  }

  template <typename dT>
  static bool rescale(dT &x,scaleType sct)
  {
    if (sct==scaleTypeVal::LOGARITHMIC) x=log((x));
    else if (sct>=2) x=pow(x,1./sct);
    else return false;
    
    return true;
  }

  template <typename dT>
  static void rescale(dT &x,dT &y,scaleType xsct,scaleType ysct)
  {
    rescale(x,xsct);
    rescale(y,ysct);
  }

  template <typename dT>
  static void rescale(dT &x,dT &y,dT &z,scaleType xsct,scaleType ysct,scaleType zsct)
  {
    rescale(x,xsct);
    rescale(y,ysct);
    rescale(z,zsct);
  }

  template <typename dT>
  static void unrescale(dT &x,scaleType sct)
  {
    if (sct==scaleTypeVal::LOGARITHMIC) x=exp(x);
    else if (sct>=2) x=pow(x,sct);
  }

  template <typename dT>
  static void unrescale(dT &x,dT &y,scaleType xsct,scaleType ysct)
  {
    unrescale(x,xsct);
    unrescale(y,ysct);
  }

  template <typename dT>
  static void unrescale(dT &x,dT &y,dT &z,scaleType xsct,scaleType ysct,scaleType zsct)
  {
    unrescale(x,xsct);
    unrescale(y,ysct);
    unrescale(z,zsct);
  }
};

#endif
