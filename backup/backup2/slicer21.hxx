#ifndef __SLICER_21_TRAITS_HXX__
#define __SLICER_21_TRAITS_HXX__

template <class gridType>
struct slicer21 {
  typedef gridType gridT;
  typedef typename gridType::paramT paramT;
  typedef typename gridType::scaleT scaleT;
  typedef typename scaleT::scaleTypeV scaleTypeV;

  static long defaultPmargin(int i) {return (i==0)?3:0;};
  static long defaultUmargin(int i) {return 0;};

  static void setFullGridMargin(paramT &gp)
  {
    long i;

    for (i=0;i<gridT::P_DIMS;i++) gp.lowPmargin[i]=gp.highPmargin[i]=0;
    for (i=0;i<gridT::U_DIMS;i++) gp.lowUmargin[i]=gp.highUmargin[i]=0;
  }

  static paramT getSubGridParams(const paramT &gp, long index,long nGrids)
  {
    paramT out=gp;
    long deltaI=(long)((gp.Pres[0])/nGrids);
    double resD=gp.Pres[0];
    long xmin,xmax;
    long defaultMarginSize = defaultPmargin(0);

    long i;
    for (i=0;i<gridT::P_DIMS;i++) 
      {
	if (index!=0) out.lowPmargin[i]=defaultPmargin(i);
	else out.lowPmargin[i]=0;
	if (index!=nGrids-1) out.highPmargin[i]=defaultPmargin(i);
	else out.highPmargin[i]=0;
      }
    for (i=0;i<gridT::U_DIMS;i++) 
      {
	out.lowUmargin[i]=defaultUmargin(i);
	out.highUmargin[i]=defaultUmargin(i);
      }
    
    xmin=(index*deltaI) - out.lowPmargin[0];
   
    if (index==nGrids-1) 
      {
	out.Pres[0] = (gp.Pres[0] - deltaI*index) + out.lowPmargin[i] + out.highPmargin[i];
	xmax=gp.Pres[0] + out.highPmargin[i];
      }
    else 
      {
	out.Pres[0] = deltaI + out.lowPmargin[i] + out.highPmargin[i];
	xmax=((index+1)*deltaI) + out.highPmargin[i];
      }
    
    out.Pmin[0]=scaleT::valueAt(gp.Pmin[0],gp.Pmax[0],gp.Pres[0],xmin,gp.Pscale[0],scaleT::valLocationV::VERTEX);
    out.Pmax[0]=scaleT::valueAt(gp.Pmin[0],gp.Pmax[0],gp.Pres[0],xmax,gp.Pscale[0],scaleT::valLocationV::VERTEX);
    
    if (xmin==0) out.Pmin[0]=gp.Pmin[0];
    if (xmax==gp.Pres[0]) out.Pmax[0]=gp.Pmax[0];

    //printf("index(%ld) : [%ld,%ld]->[%g,%g] res=%d\n",index,xmin,xmax,out.Pmin[0],out.Pmax[0],out.Pres[0]);

    return out;
  }

};

#endif
