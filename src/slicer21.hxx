#ifndef __SLICER_21_HXX__
#define __SLICER_21_HXX__

#include "slicer.hxx"

template <class gridType>
class slicer21 : public slicer<gridType> {
public:
  typedef slicer<gridType> bT;
  typedef gridType gridT;

  typedef typename bT::paramT paramT;
  typedef typename bT::scaleT scaleT;
  typedef typename bT::scaleTypeV scaleTypeV;
  typedef typename bT::dirT dirT;

  static const int P_DIMS = bT::P_DIMS;
  static const int U_DIMS = bT::U_DIMS;
  static const int DIMS = bT::DIMS;

  static long defaultHPmargin(int i) {return 1;}
  static long defaultHUmargin(int i) {return 0;}
  static long defaultLPmargin(int i) {return 0;}
  static long defaultLUmargin(int i) {return 0;}

private:
  paramT fullGridParams;
  paramT gridParams;
  std::vector<double> fullGrid[3];
    
  int slicerSize[2];  // number of slices in each dim
  int myID;     
  int myCoords[2]; // coords within the slices 
  int myFullGridPos[2]; // coords within the fullgrid

public:

  /*
  static void setDefaultFullGridMargin(paramT &gp)
  {
    long i;

    for (i=0;i<P_DIMS;i++) gp.lowPmargin[i]=gp.highPmargin[i]=0;
    for (i=0;i<U_DIMS;i++) gp.lowUmargin[i]=gp.highUmargin[i]=0;
  }
  */

  void init(const paramT &gp, long index,long nGrids,long nthreads, bool periodic=false)
  {
    int nSlicePerGrid=(int)(gp.Ures[1]/nthreads);

    if (nSlicePerGrid<1) nSlicePerGrid=1;

    int nNodes=nGrids;
    int nSlices=nSlicePerGrid;

    if (nSlices>=nNodes)
      {
	nSlices=nNodes;
	nNodes=1;
      }
    else
      {
	while(nSlices>0)
	  {
	    int ip=(int)(nNodes/nSlices);
	    if (nNodes==ip*nSlices) 
	      {
		nNodes=ip;
		break;
	      }
	  };
      }
    printf("SLICER : grid will be sliced into %dx%d chuncks.\n",nSlices,nNodes);
  
    // COMMENTER ICI
    // -> sinon, ca force la paralellisation en P uniquement !!!
    /**
    nNodes=nGrids;
    nSlices=1;
    **/

    fullGridParams=gp;
    fullGrid[0]=scaleT::genScale(gp.Pmin[0],gp.Pmax[0],gp.Pres[0],gp.Pscale[0],gridT::valLocation_P);
    fullGrid[1]=scaleT::genScale(gp.Umin[0],gp.Umax[0],gp.Ures[0],gp.Uscale[0],gridT::valLocation_U);
    fullGrid[2]=scaleT::genScale(gp.Umin[1],gp.Umax[1],gp.Ures[1],gp.Uscale[1],gridT::valLocation_J);

    myID=index;   
    slicerSize[0]=nNodes;
    slicerSize[1]=nSlices;    
    myCoords[0]=myID%slicerSize[0];
    myCoords[1]=myID/slicerSize[0];   

    paramT gp2=this->divide(gp,myCoords[1],slicerSize[1],0,0,periodic,2,myFullGridPos[1]);
    gridParams=this->divide(gp2,myCoords[0],slicerSize[0],defaultLPmargin(0),defaultHPmargin(0),periodic,0,myFullGridPos[0]);    
  }

  slicer21()
  {
    
  }

  slicer21(const paramT &gp, long index,long nGrids,long nthreads)
  {
    init(gp,index,nGrids,nthreads);
  }

  ~slicer21()
  {
    
  }

  paramT getGridParams()
  {
    return gridParams;
  }

  long getMyID()
  {
    return myID;
  }

  long getNeiID(dirT dir)
  {
    if (slicerSize[0]<=1) return -1;
		      
    if (dir==gridNav::dir(0,1))
      {
	if (myCoords[0]+1 >= slicerSize[0]) return -1;
	else return myID+1;
      }
    else if (dir==gridNav::dir(0,-1))
      {
	if (myCoords[0]<=0) return -1;
	return myID-1;
      }
    return -1;
  }

  int getMyCoord(int dir)
  {
    if (dir==0) return myCoords[0];
    else if (dir==2) return myCoords[1];
    else return 0;
  }

  int getSlicerSize(int dir)
  {
    if (dir==0) return slicerSize[0];
    else if (dir==2) return slicerSize[1];
    else return 0;
  }

  int getMyFullGridPos(int dir)
  {
    if (dir==0) return myFullGridPos[0];
    else if (dir==2) return myFullGridPos[1];
    else return 0;
  }

  int getFullGridSize(int dir)
  {
    return fullGrid[dir].size();
  }

  const std::vector<double> &getFullGrid(int dir)
  {
    return fullGrid[dir];
  }

  /*
  static paramT getSubGridParams(const paramT &gp, long index,long nGrids,long nthreads)
  {
    int nSlicePerGrid=(int)(gp.Ures[1]/nthreads);

    if (nSlicePerGrid<1) nSlicePerGrid=1;

    int nNodes=nGrids;
    int nSlices=nSlicePerGrid;

    if (nSlices>=nNodes)
      {
	nSlices=nNodes;
	nNodes=1;
      }
    else
      {
	while(nSlices>0)
	  {
	    int ip=(int)(nNodes/nSlices);
	    if (nNodes==ip*nSlices) 
	      {
		nNodes=ip;
		break;
	      }
	  };
      }
    printf("SLICER : grid will be sliced into %dx%d chuncks.\n",nSlices,nNodes);
    
    paramT gp2=divide(gp,index%nSlices,nSlices,0,0,false,2);
    return divide(gp2,index/nSlices,nNodes,defaultLPmargin(0),defaultHPmargin(0),false,0);    
  }
  */
};

#endif
