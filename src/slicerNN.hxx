#ifndef __SLICER_NN_HXX__
#define __SLICER_NN_HXX__

#include "slicer.hxx"

template <class gridType>
class slicerNN : public slicer<gridType> {
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
  std::vector<double> fullGrid[DIMS];
    
  int slicerSize[DIMS];  // number of slices in each dim
  int myID;     
  int myCoords[DIMS]; // coords within the slices 
  int myFullGridPos[DIMS]; // coords within the fullgrid

public:

  void init(const paramT &gp, long index,long nGrids,long nthreads, bool periodic=false)
  {
 
   
  
    int i;
    fullGridParams=gp;
    for (i=0;i<P_DIMS;i++)
      fullGrid[i]=scaleT::genScale(gp.Pmin[i],gp.Pmax[i],gp.Pres[i],gp.Pscale[i],gridT::valLocation_P);
    for (i=0;i<U_DIMS;i++)
      fullGrid[P_DIMS+i]=scaleT::genScale(gp.Umin[i],gp.Umax[i],gp.Ures[i],gp.Uscale[i],gridT::valLocation_U);
    

    myID=index;   

    slicerSize[0]=nGrids;
    for (int i=1;i<DIMS;i++) slicerSize[i]=1;

    myCoords[0]=myID%slicerSize[0];
    myCoords[1]=myID/slicerSize[0];   
    for (int i=1;i<DIMS;i++) myCoords[i]=0;

    for (i=0;i<DIMS;i++) myFullGridPos[i]=0;
    //paramT gp2=this->divide(gp,myCoords[1],slicerSize[1],0,0,periodic,1,myFullGridPos[1]);
    gridParams=this->divide(gp,myCoords[0],slicerSize[0],defaultLPmargin(0),defaultHPmargin(0),periodic,0,myFullGridPos[0]);    

    printf("SLICER : grid %ld will be sliced into %ldx%ld chuncks.\n",index,nGrids,1l);
    /*
    printf(" margin: [%d %d][%d %d]--[%d %d][%d %d].\n",
	   gridParams.lowPmargin[0],gridParams.highPmargin[0],
	   gridParams.lowPmargin[1],gridParams.highPmargin[1],
	   gridParams.lowUmargin[0],gridParams.highUmargin[0],
	   gridParams.lowUmargin[1],gridParams.highUmargin[1]);
    */
  }

  slicerNN()
  {
    
  }

  slicerNN(const paramT &gp, long index,long nGrids,long nthreads)
  {
    init(gp,index,nGrids,nthreads);
  }

  ~slicerNN()
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
    return (dir>=DIMS)?0:myCoords[0];
  }

  int getSlicerSize(int dir)
  {
    return (dir>=DIMS)?0:slicerSize[dir];
  }

  int getMyFullGridPos(int dir)
  {
    return (dir>=DIMS)?0:myFullGridPos[dir];
  }

  int getFullGridSize(int dir)
  {
    return (dir>=DIMS)?0:fullGrid[dir].size();
  }

  const std::vector<double> &getFullGrid(int dir)
  {
    return fullGrid[dir];
  }
};

#endif
