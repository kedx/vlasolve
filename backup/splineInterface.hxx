#ifndef __SPLINE_INTERFACE_HXX__
#define __SPLINE_INTERFACE_HXX__

#include "einspline/bspline.h"
#include "scale.hxx"

template <int D, typename T> struct einspline;

template <> struct einspline<2,double>
{    
  typedef Ugrid ugridT;
  typedef UBspline_2d_d UBsplineT;
  
  typedef scale<double> scaleT;
  typedef scaleT::scaleType scaleType; 
  typedef scaleT::scaleTypeVal scaleTypeVal;
  typedef double dataT;
  typedef BCtype_d bcT;

  struct initT {
    ugridT xgrid; ugridT ygrid; 
    bcT xbc; bcT ybc; 
    scaleType xscaleType; scaleType yscaleType;
    bool rescaleNeeded;
  };

  struct splineT {
    UBsplineT *spline;
    scaleType xscaleType;
    scaleType yscaleType;
    bool rescaleNeeded;
  };
  
  static initT createInit(double xstart,double xstop,int Nx,
			  double ystart,double ystop,int Ny,			  
			  bc_code xbctL, bc_code xbctR,
			  bc_code ybctL, bc_code ybctR,	
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

  static splineT *create(initT &ini, dataT *data){
    splineT *spline=new splineT;
    spline->spline=create_UBspline_2d_d (ini.ygrid,ini.xgrid,ini.ybc,ini.xbc,data);
    spline->xscaleType = ini.xscaleType;
    spline->yscaleType = ini.yscaleType;
    spline->rescaleNeeded = ini.rescaleNeeded;
    return spline;
  }

  
  static splineT *create(ugridT &x, ugridT &y, bcT &xbc, bcT &ybc, dataT *data,
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

  static void evaluate(splineT *spline, double x, double y, dataT *v){
    if (!spline->rescaleNeeded) return eval_UBspline_2d_d(spline->spline,y,x,v);    
    scaleT::rescale(x,y,spline->xscaleType,spline->yscaleType);        
    if (x>spline->spline->y_grid.end) x=spline->spline->y_grid.end;
    if (x<spline->spline->y_grid.start) x=spline->spline->y_grid.start;
    if (y<spline->spline->x_grid.start) y=spline->spline->x_grid.start;
    if (y>spline->spline->x_grid.end) y=spline->spline->x_grid.end;
    
    eval_UBspline_2d_d(spline->spline,y,x,v);
    //scaleT::unrescale(x,y,spline->xscaleType,spline->yscaleType);
  }

  static void evaluate(splineT *spline, double x, double y, dataT *v, dataT *g){
    if (!spline->rescaleNeeded) return eval_UBspline_2d_d_vg(spline->spline,y,x,v,g);
    scaleT::rescale(x,y,spline->xscaleType,spline->yscaleType);
    if (x<spline->spline->y_grid.start) x=spline->spline->y_grid.start;
    if (x>spline->spline->y_grid.end) x=spline->spline->y_grid.end;
    if (y<spline->spline->x_grid.start) y=spline->spline->x_grid.start;
    if (y>spline->spline->x_grid.end) y=spline->spline->x_grid.end;
    
    eval_UBspline_2d_d_vg(spline->spline,y,x,v,g);
    //scaleT::unrescale(x,y,spline->xscaleType,spline->yscaleType);
  }

  static void evaluate(splineT *spline, double x, double y, dataT *v, dataT *g, dataT *h){
    if (!spline->rescaleNeeded) return eval_UBspline_2d_d_vgh(spline->spline,y,x,v,g,h);
    scaleT::rescale(x,y,spline->xscaleType,spline->yscaleType);
    if (x<spline->spline->y_grid.start) x=spline->spline->y_grid.start;
    if (x>spline->spline->y_grid.end) x=spline->spline->y_grid.end;
    if (y<spline->spline->x_grid.start) y=spline->spline->x_grid.start;
    if (y>spline->spline->x_grid.end) y=spline->spline->x_grid.end;
    
    eval_UBspline_2d_d_vgh(spline->spline,y,x,v,g,h);
    //scaleT::unrescale(x,y,spline->xscaleType,spline->yscaleType);
  }

  static void destroy(splineT *spline){
    destroy_Bspline(spline->spline);
    delete spline;
    spline=NULL;
  }
};

#endif
