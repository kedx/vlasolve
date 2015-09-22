#ifndef __INTERPOL_DUMMY2_INTERFACE_HXX__
#define __INTERPOL_DUMMY2_INTERFACE_HXX__

#include "interpol.hxx"
#include "scale.hxx"

template <typename dataT> 
struct interpol<interpType::DUMMY,dataT>
{     
  
  typedef scale<dataT> scaleT;
  typedef typename scaleT::scaleType scaleType; 
  typedef typename scaleT::scaleTypeVal scaleTypeVal;
  typedef typename scaleT::valLocationType valLocationType;
  typedef typename scaleT::valLocationTypeVal valLocationTypeVal;
  
  struct boundaryTypeVal {
    enum type {BT_PERIODIC, BT_DERIV1, 
	       BT_DERIV2, BT_FLAT, 
	       BT_NATURAL, BT_ANTIPERIODIC};
  };
  typedef typename boundaryTypeVal::type boundaryType;
  
  struct ugridT{
    double start;
    double end;
    int num;
  };

  struct BCType {
    boundaryType lCode;
    boundaryType rCode;
    double lVal;
    double rVal;
  };

  struct initT {
    ugridT xgrid; ugridT ygrid; 
    BCType xbc; BCType ybc;     
    scaleType xscaleType; scaleType yscaleType;
    bool rescaleNeeded;
  };

  struct interpT {
    scaleType xscaleType;
    scaleType yscaleType;
    bool rescaleNeeded;
  };

  
  static initT createInit(double xstart,double xstop,int Nx,
			  double ystart,double ystop,int Ny,
			  boundaryType xbctL, boundaryType xbctR,
			  boundaryType ybctL, boundaryType ybctR,			  
			  scaleType xscaleType=scaleTypeVal::REGULAR,
			  scaleType yscaleType=scaleTypeVal::REGULAR,
			  valLocationType vltype= valLocationTypeVal::PIXEL,
			  double xbcvL=0,double xbcvR=0,
			  double ybcvL=0,double ybcvR=0)
  {
    initT g;
    //ugridT xgrid,ygrid;
    g.xgrid.start=xstart;g.ygrid.start=ystart;
    g.xgrid.end=xstop;g.ygrid.end=ystop;
    g.xgrid.num=Nx;g.ygrid.num=Ny;
    
    //bcT xbc, ybc;
    g.xbc.lCode = xbctL;g.xbc.rCode = xbctR;
    g.ybc.lCode = ybctL;g.ybc.rCode = ybctR;
    g.xbc.lVal = xbcvL;g.xbc.rVal = xbcvR;
    g.ybc.lVal = ybcvL;g.ybc.rVal = ybcvR;

    g.xscaleType=xscaleType;
    g.yscaleType=yscaleType;

    g.rescaleNeeded=false;

    g.rescaleNeeded |= scaleT::rescale(g.xgrid.start,xscaleType);
    scaleT::rescale(g.xgrid.end,xscaleType);
    g.rescaleNeeded |= scaleT::rescale(g.ygrid.start,yscaleType);
    scaleT::rescale(g.ygrid.end,yscaleType);
    
    return g;
  }

  static interpT *create(initT &ini, dataT *data){
    interpT *dummy=new interpT;
    dummy->xscaleType = ini.xscaleType;
    dummy->yscaleType = ini.yscaleType;
    dummy->rescaleNeeded = ini.rescaleNeeded;
    return dummy;
  }

  static void evaluate(interpT *dummy, double x, double y, dataT *v){
    if (dummy->rescaleNeeded) scaleT::rescale(x,y,dummy->xscaleType,dummy->yscaleType);        
  }

  static void evaluate(interpT *dummy, double x, double y, dataT *v, dataT *g){
  }

  static void evaluate(interpT *dummy, double x, double y, dataT *v, dataT *g, dataT *h){
  }

  static void destroy(interpT *dummy){
    delete dummy;
    dummy=NULL;
  }
};

#endif
