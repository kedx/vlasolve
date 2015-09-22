#ifndef __SLICER_HXX__
#define __SLICER_HXX__

#include "gridNav.hxx"
#include "valLocationType.hxx"

template <class gridType>
class slicer {
public:
  typedef slicer<gridType> myT;
  typedef gridType gridT;
  typedef typename gridT::paramT paramT;
  typedef typename gridT::scaleT scaleT;
  typedef typename gridT::valLocationT valLocationT;
  typedef typename gridT::valLocationV valLocationV;
  typedef typename scaleT::scaleTypeV scaleTypeV;
  typedef typename scaleT::scaleTypeT scaleTypeT;
  typedef typename gridNav::dirT dirT;
  

  static const int P_DIMS = gridT::P_DIMS;
  static const int U_DIMS = gridT::U_DIMS;
  static const int J_DIMS = gridT::J_DIMS;
  static const int DIMS = gridT::DIMS;

  //static const valLocationT valLocation_P = gridT::valLocation_P;
  //static const valLocationT valLocation_U = gridT::valLocation_U;
  //static const valLocationT valLocation_J = gridT::valLocation_J;

protected:

  static paramT divide(const paramT &gp, long index, long N, int lowMargin,int highMargin, bool periodic, int dir, int &fullgridpos)
  {
    paramT out=gp;  
    long deltaI;
    
    long min,max;  
    const bool useJ=(dir==(DIMS-J_DIMS));
    const bool useU=(dir>=P_DIMS);
    valLocationT vl;

    if (useJ) vl=gridT::valLocation_J;
    else if (useU) vl=gridT::valLocation_U;
    else vl=gridT::valLocation_P;
    

    int &out_lowMargin=(useU)?out.lowUmargin[dir-P_DIMS]:out.lowPmargin[dir];
    int &out_highMargin=(useU)?out.highUmargin[dir-P_DIMS]:out.highPmargin[dir];
    const int &gp_res=(useU)?gp.Ures[dir-P_DIMS]:gp.Pres[dir];
    int &out_res=(useU)?out.Ures[dir-P_DIMS]:out.Pres[dir];

    if ((index!=0)||(periodic)) out_lowMargin=lowMargin;
    if ((index!=N-1)||(periodic)) out_highMargin=highMargin;
   
    if (vl == gridT::valLocationV::VERTEX) 
      deltaI=(long)((gp_res+1)/N);
    else
      deltaI=(long)((gp_res)/N);

    min=(index*deltaI) - out_lowMargin;

    if (index==N-1)
      {
	out_res = (gp_res - deltaI*index) + out_lowMargin + out_highMargin;
	max= gp_res + out_highMargin;
      }
    else 
      {
	out_res = deltaI + out_lowMargin + out_highMargin;
	max=((index+1)*deltaI) + out_highMargin;

	if (vl == gridT::valLocationV::VERTEX)
	  {
	    out_res--;
	    max--;
	  }
      }

    const double &gp_min=(useU)?gp.Umin[dir-P_DIMS]:gp.Pmin[dir];
    double &out_min=(useU)?out.Umin[dir-P_DIMS]:out.Pmin[dir];
    const double &gp_max=(useU)?gp.Umax[dir-P_DIMS]:gp.Pmax[dir];
    double &out_max=(useU)?out.Umax[dir-P_DIMS]:out.Pmax[dir];
    const scaleTypeT &gp_scale=(useU)?gp.Uscale[dir-P_DIMS]:gp.Pscale[dir];

    out_min=scaleT::valueAt(gp_min,gp_max,gp_res,min,gp_scale,gridT::valLocationV::VERTEX);
    out_max=scaleT::valueAt(gp_min,gp_max,gp_res,max,gp_scale,gridT::valLocationV::VERTEX);
    
    if (min==0) out_min=gp_min;
    if (max==gp_res) out_max=gp_max;
    
    fullgridpos = min;
   
    printf("index(%ld) :(%e,%e,%d)==> [%ld,%ld]->[%g,%g] res=%d\n",index,gp_min,gp_max,gp_res,min,max,out_min,out_max,out_res);

    return out;
  }
  
  
};


#endif
