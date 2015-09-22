#ifndef __INTERPOLATION_INTERFACE_HXX__
#define __INTERPOLATION_INTERFACE_HXX__

//#include "scale.hxx"

struct interpType {
  enum type {DUMMY,SPLINE,PFC2};
};

template <interpType::type interpTraits, typename dataT> struct interpol;

#include "interpol_DUMMY.hxx"
#include "interpol_spline.hxx"
#include "interpol_PFC2.hxx"



/*
template <interpType::type interpTraits, typename dataT> struct interpol{
  enum bTypeVal {BT_PERIODICC, BT_DERIV1, 
		 BT_DERIV2, BT_FLAT, 
		 BT_NATURALL, BT_ANTIPERIODIC};
  typedef int initT;
  typedef void interpT;
  typedef scale<double> scaleT;
  typedef scaleT::scaleType scaleType; 
  typedef scaleT::scaleTypeVal scaleTypeVal;
  static initT createInit(double xstart,double xstop,int Nx,
			  double ystart,double ystop,int Ny,
			  bTypeVal xbctL, bTypeVal xbctR,
			  bTypeVal ybctL, bTypeVal ybctR,	
			  scaleType xscaleType=scaleTypeVal::REGULAR,
			  scaleType yscaleType=scaleTypeVal::REGULAR,
			  double xbcvL=0,double xbcvR=0,
			  double ybcvL=0,double ybcvR=0)
  {};

  static interpT *create(initT &ini, dataT *data){};
  static void evaluate(interpT *interp, double x, double y, dataT *v){};
  static void evaluate(interpT *interp, double x, double y, dataT *v, dataT *g){};
  static void evaluate(interpT *interp, double x, double y, dataT *v, dataT *g, dataT *h){};
  static void destroy(interpT *interp){};
};
*/

#endif
