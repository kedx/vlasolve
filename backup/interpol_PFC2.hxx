#ifndef __INTERPOL_PFC2_INTERFACE_HXX__
#define __INTERPOL_PFC2_INTERFACE_HXX__

#include "interpol.hxx"
#include "scale.hxx"

template <typename dataT> 
struct interpol<interpType::PFC2,dataT>
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

    double dxh;
    double dyh;

    std::vector<double> xh;
    std::vector<double> yh;  

    std::vector<double> x;
    std::vector<double> y;  
  };

  struct interpT {
    initT ini;
    std::vector<dataT> data;
    /*
    scaleType xscaleType;
    scaleType yscaleType;
    bool rescaleNeeded;

    double dxh;
    double dyh;

    std::vector<dataT> data;
    std::vector<double> xh;
    std::vector<double> yh;
    */
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
    long i;

    if (vltype!=valLocationTypeVal::PIXEL)
      {
	fprintf(stderr,"ERROR in interpol_pfc2.hxx: Value location type can only be PIXELS for this type of interpolator.\n");
	exit(0);
      }

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

    
    g.x = scaleT::genScale(xstart,xstop,Nx,xscaleType,valLocationTypeVal::VERTEX);
    g.xh = scaleT::genScale(xstart,xstop,Nx,xscaleType,valLocationTypeVal::PIXEL);
    g.y = scaleT::genScale(ystart,ystop,Ny,yscaleType,valLocationTypeVal::VERTEX);
    g.yh = scaleT::genScale(ystart,ystop,Ny,yscaleType,valLocationTypeVal::PIXEL);

    return g;
  }

  static interpT *create(initT &ini, dataT *data){
    interpT *pfc=new interpT;
    long Nx=ini.xh.size();
    long Ny=ini.yh.size();    
    long N=Nx*Ny;

    pfc->ini=ini;
    pfc->data.resize(Nx*Ny);

    std::vector<dataT> &D = pfc->data;
    long xi,yi,i;
    double deltaX,deltaY;

    
    for (yi=0,i=0;yi<Ny;yi++)
      {
	for (xi=0;xi<Nx;xi++,i++)
	  {
	    D[i]=data[i];
	    
	  }
      }

    

    return pfc;
  }

  static void evaluate(interpT *pfc, double x, double y, dataT *v){
    initT &ini = pfc->ini;
    

    /*
    initT &ini = pfc->ini;
    if (ini.rescaleNeeded) scaleT::rescale(x,y,ini.xscaleType,ini.yscaleType);  
    
    
    long Nhx=ini.xh.size();
    long Nhy=ini.yh.size();
    long Nht=Nx*Ny;    

    long xi=(long)((x-ini.xh[0])/(ini.dxh));
    long yi=(long)((y-ini.yh[0])/(ini.dyh)) * Nhx;

    double xxi = ini.xh[0] + xi*ini.dxh;
    double yyi = ini.yh[0] + yi*ini.dyh;
    
    double xxd[5]={xi-3*ini.dxh,xi-ini.dxh,xi,xi+ini.dxh,xi+3*ini.dxh};
    double yyd[5]={yi-3*ini.dyh,yi-ini.dyh,yi,yi+ini.dyh,yi+3*ini.dyh};

    double f=pfc->data[(xi>>1)+(yi>>1)*ini.xgrid.num];

    double fm1;
    if (xi==0) {
    }
    else xm1= xi-1;
    if (xi==0) {
    }
    else xp1= xi+1;

    if (yi==0) {
    }
    else ym1= yi-1;
    if (yi==0) {
    }
    else yp1= yi+1;
    */
    
  }

  static void evaluate(interpT *pfc, double x, double y, dataT *v, dataT *g){
  }

  static void evaluate(interpT *pfc, double x, double y, dataT *v, dataT *g, dataT *h){
  }

  static void destroy(interpT *pfc){
    delete pfc;
    pfc=NULL;
  }
};

#endif
