#ifndef __INTERPOL_SPLINE_INTERFACE_HXX__
#define __INTERPOL_SPLINE_INTERFACE_HXX__

#include "interpol.hxx"

#include "einspline/bspline.h"
#include "interpol.hxx"
#include "scale.hxx"

template <> struct interpol<interpType::SPLINE,double>
{   
  typedef Ugrid ugridT;
  typedef UBspline_2d_d UBsplineT;
  
  typedef scale<double> scaleT;
  typedef scaleT::scaleType scaleType; 
  typedef scaleT::scaleTypeVal scaleTypeVal;
  typedef double dataT;
  typedef BCtype_d bcT;

  enum bTypeVal {BT_PERIODIC=PERIODIC, BT_DERIV1= DERIV1, 
		 BT_DERIV2=DERIV2, BT_FLAT=FLAT, 
		 BT_NATURAL=NATURAL, BT_ANTIPERIODIC=ANTIPERIODIC};

  struct initT {
    ugridT xgrid; ugridT ygrid; 
    bcT xbc; bcT ybc; 
    scaleType xscaleType; scaleType yscaleType;
    bool rescaleNeeded;
  };

  struct interpT {
    UBsplineT *spline;
    scaleType xscaleType;
    scaleType yscaleType;
    bool rescaleNeeded;
  };

  
  static initT createInit(double xstart,double xstop,int Nx,
			  double ystart,double ystop,int Ny,
			  bTypeVal xbctL, bTypeVal xbctR,
			  bTypeVal ybctL, bTypeVal ybctR,	
			  scaleType xscaleType=scaleTypeVal::REGULAR,
			  scaleType yscaleType=scaleTypeVal::REGULAR,
			  double xbcvL=0,double xbcvR=0,
			  double ybcvL=0,double ybcvR=0)
  {
    initT g;
    //ugridT xgrid,ygrid;
    g.xgrid.start=xstart;g.ygrid.start=ystart;
    g.xgrid.end=xstop;g.ygrid.end=ystop;
    g.xgrid.num=Nx;g.ygrid.num=Ny;
    
    //bcT xbc, ybc;
    g.xbc.lCode = (bc_code)xbctL;g.xbc.rCode = (bc_code)xbctR;
    g.ybc.lCode = (bc_code)ybctL;g.ybc.rCode = (bc_code)ybctR;
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
    interpT *spline=new interpT;
    spline->spline=create_UBspline_2d_d (ini.ygrid,ini.xgrid,ini.ybc,ini.xbc,data);
    spline->xscaleType = ini.xscaleType;
    spline->yscaleType = ini.yscaleType;
    spline->rescaleNeeded = ini.rescaleNeeded;
    return spline;
  }

  /*
  static interpT *create(ugridT &x, ugridT &y, bcT &xbc, bcT &ybc, dataT *data,
			 scaleType xscaleType=scaleTypeVal::REGULAR,
			 scaleType yscaleType=scaleTypeVal::REGULAR){
    initT g;
    g.xgrid=x;
    g.ygrid=y;
    g.xbc=xbc;
    g.ybc=ybc;
    g.xscaleType=xscaleType;
    g.yscaleType=yscaleType;
    g.rescaleNeeded=false;

    g.rescaleNeeded=false;

    g.rescaleNeeded |= scaleT::rescale(g.xgrid.start,xscaleType);
    scaleT::rescale(g.xgrid.end,xscaleType);
    g.rescaleNeeded |= scaleT::rescale(g.ygrid.start,yscaleType);
    scaleT::rescale(g.ygrid.end,yscaleType);
    
    return create(g,data);
  }
  */
  static void evaluate(interpT *spline, double x, double y, dataT *v){
    if (!spline->rescaleNeeded) return eval_UBspline_2d_d(spline->spline,y,x,v);    
    scaleT::rescale(x,y,spline->xscaleType,spline->yscaleType);        
    if (x>spline->spline->y_grid.end) x=spline->spline->y_grid.end;
    if (x<spline->spline->y_grid.start) x=spline->spline->y_grid.start;
    if (y<spline->spline->x_grid.start) y=spline->spline->x_grid.start;
    if (y>spline->spline->x_grid.end) y=spline->spline->x_grid.end;
    
    eval_UBspline_2d_d(spline->spline,y,x,v);
    //scaleT::unrescale(x,y,spline->xscaleType,spline->yscaleType);
  }

  static void evaluate(interpT *spline, double x, double y, dataT *v, dataT *g){
    if (!spline->rescaleNeeded) return eval_UBspline_2d_d_vg(spline->spline,y,x,v,g);
    scaleT::rescale(x,y,spline->xscaleType,spline->yscaleType);
    if (x<spline->spline->y_grid.start) x=spline->spline->y_grid.start;
    if (x>spline->spline->y_grid.end) x=spline->spline->y_grid.end;
    if (y<spline->spline->x_grid.start) y=spline->spline->x_grid.start;
    if (y>spline->spline->x_grid.end) y=spline->spline->x_grid.end;
    
    eval_UBspline_2d_d_vg(spline->spline,y,x,v,g);
    //scaleT::unrescale(x,y,spline->xscaleType,spline->yscaleType);
  }

  static void evaluate(interpT *spline, double x, double y, dataT *v, dataT *g, dataT *h){
    if (!spline->rescaleNeeded) return eval_UBspline_2d_d_vgh(spline->spline,y,x,v,g,h);
    scaleT::rescale(x,y,spline->xscaleType,spline->yscaleType);
    if (x<spline->spline->y_grid.start) x=spline->spline->y_grid.start;
    if (x>spline->spline->y_grid.end) x=spline->spline->y_grid.end;
    if (y<spline->spline->x_grid.start) y=spline->spline->x_grid.start;
    if (y>spline->spline->x_grid.end) y=spline->spline->x_grid.end;
    
    eval_UBspline_2d_d_vgh(spline->spline,y,x,v,g,h);
    //scaleT::unrescale(x,y,spline->xscaleType,spline->yscaleType);
  }

  static void destroy(interpT *spline){
    destroy_Bspline(spline->spline);
    delete spline;
    spline=NULL;
  }
};

#endif
