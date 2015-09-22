#ifndef __INTERPOL_SPLINE_HXX__
#define __INTERPOL_SPLINE_HXX__

#include "einspline/bspline.h"
#include "scale.hxx"

template <
  typename dataType,
  class scaleType=scale<dataType>
  >
struct interpolSpline
{   
  typedef Ugrid ugridT;
  typedef UBspline_2d_d UBsplineT;

  typedef dataType dataT;
  typedef scaleType scaleT;

  typedef typename scaleT::scaleTypeT scaleTypeT; 
  typedef typename scaleT::scaleTypeV scaleTypeV;
  typedef typename scaleT::valLocationT valLocationT;
  typedef typename scaleT::valLocationV valLocationV;  

  struct boundaryTypeV {
    enum type {BT_PERIODIC=PERIODIC, BT_DERIV1= DERIV1, 
	       BT_DERIV2=DERIV2, BT_FLAT=FLAT, 
	       BT_NATURAL=NATURAL, BT_ANTIPERIODIC=ANTIPERIODIC};
  };
  typedef typename boundaryTypeV::type boundaryTypeT;

  typedef BCtype_d bcT;

  struct initT {
    ugridT xgrid; ugridT ygrid; 
    bcT xbc; bcT ybc; 
    scaleTypeT xscaleType; scaleTypeT yscaleType;
    bool rescaleNeeded;
    bool useLog;
  };

  struct interpT {
    UBsplineT *spline;
    scaleTypeT xscaleType;
    scaleTypeT yscaleType;
    bool rescaleNeeded;
    bool useLog;
  };

  
  static initT createInit(double xstart,double xstop,int Nx,
			  double ystart,double ystop,int Ny,
			  boundaryTypeT xbctL, boundaryTypeT xbctR,
			  boundaryTypeT ybctL, boundaryTypeT ybctR,
			  scaleTypeT xscaleType,
			  scaleTypeT yscaleType,
			  valLocationT vltype,
			  bool useLog=false,
			  double xbcvL=0,double xbcvR=0,
			  double ybcvL=0,double ybcvR=0)
  {
    initT g;
    std::vector<double> tmp;

    tmp=scaleT::genScale(xstart,xstop,Nx,xscaleType,vltype);
    tmp.back()+=(tmp[tmp.size()-1]-tmp[tmp.size()-2])*1.E-10;
    g.xgrid.start=tmp.front();
    g.xgrid.end=tmp.back();
    g.xgrid.num=tmp.size();
        
    tmp=scaleT::genScale(ystart,ystop,Ny,yscaleType,vltype);
    tmp.back()+=(tmp[tmp.size()-1]-tmp[tmp.size()-2])*1.E-10;
    g.ygrid.start=tmp.front();
    g.ygrid.end=tmp.back();
    g.ygrid.num=tmp.size();    
        
    g.useLog = useLog;

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
    dataT *slice=data;
    std::vector<dataT> tmp;
    if (ini.useLog)
      {
	tmp.assign(data,data+(ini.ygrid.num*ini.xgrid.num));
	long i;
	for (i=0;i<tmp.size();i++) tmp[i]=log(tmp[i]);
	slice =&tmp[0];
      }
    spline->spline=create_UBspline_2d_d (ini.ygrid,ini.xgrid,ini.ybc,ini.xbc,slice);
    spline->xscaleType = ini.xscaleType;
    spline->yscaleType = ini.yscaleType;
    spline->rescaleNeeded = ini.rescaleNeeded;
    spline->useLog = ini.useLog;
    return spline;
  }

  static void evaluate(interpT *spline, double x, double y, dataT *v){
    if (!spline->rescaleNeeded) eval_UBspline_2d_d(spline->spline,y,x,v);    
    else {
      scaleT::rescale(x,y,spline->xscaleType,spline->yscaleType);        
      if (x>spline->spline->y_grid.end) x=spline->spline->y_grid.end;
      if (x<spline->spline->y_grid.start) x=spline->spline->y_grid.start;
      if (y<spline->spline->x_grid.start) y=spline->spline->x_grid.start;
      if (y>spline->spline->x_grid.end) y=spline->spline->x_grid.end;
    
      eval_UBspline_2d_d(spline->spline,y,x,v);
    }

    if (spline->useLog) 
      {
	//printf("%e ->",*v);
	if (*v<-700) *v=-700;
	*v=exp(*v);
	//printf("%e.\n",*v);
      }
    //scaleT::unrescale(x,y,spline->xscaleType,spline->yscaleType);
  }

  static void evaluate(interpT *spline, double x, double y, dataT *v, dataT *g){
    if (!spline->rescaleNeeded) eval_UBspline_2d_d_vg(spline->spline,y,x,v,g);
    else {
      scaleT::rescale(x,y,spline->xscaleType,spline->yscaleType);
      if (x<spline->spline->y_grid.start) x=spline->spline->y_grid.start;
      if (x>spline->spline->y_grid.end) x=spline->spline->y_grid.end;
      if (y<spline->spline->x_grid.start) y=spline->spline->x_grid.start;
      if (y>spline->spline->x_grid.end) y=spline->spline->x_grid.end;
    
      eval_UBspline_2d_d_vg(spline->spline,y,x,v,g);
    }

    if (spline->useLog) 
      {
	if (*v<-700) *v=-700;
	*v=exp(*v);
	
	g[0] = -g[0]*v[0];
	g[1] = -g[1]*v[0];
      }
    //scaleT::unrescale(x,y,spline->xscaleType,spline->yscaleType);
  }

  static void evaluate(interpT *spline, double x, double y, dataT *v, dataT *g, dataT *h){
    if (!spline->rescaleNeeded) eval_UBspline_2d_d_vgh(spline->spline,y,x,v,g,h);
    else {
      scaleT::rescale(x,y,spline->xscaleType,spline->yscaleType);
      if (x<spline->spline->y_grid.start) x=spline->spline->y_grid.start;
      if (x>spline->spline->y_grid.end) x=spline->spline->y_grid.end;
      if (y<spline->spline->x_grid.start) y=spline->spline->x_grid.start;
      if (y>spline->spline->x_grid.end) y=spline->spline->x_grid.end;
      
      eval_UBspline_2d_d_vgh(spline->spline,y,x,v,g,h);
    }

    if (spline->useLog) 
    {
      if (*v<-700) *v=-700;
      *v=exp(*v);

      g[0] = -g[0]*v[0];
      g[1] = -g[1]*v[0];
      fprintf(stderr,"ERROR: second derivative evaluation for spline with log valued function not implemented.\n");
      exit(-1);
    }
    //scaleT::unrescale(x,y,spline->xscaleType,spline->yscaleType);
  }

  static void destroy(interpT *spline){
    destroy_Bspline(spline->spline);
    delete spline;
    spline=NULL;
  }
};

#endif
