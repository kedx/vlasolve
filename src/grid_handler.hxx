#ifndef __GRID_HANDLER_HXX__
#define __GRID_HANDLER_HXX__

#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>

#include "mpiCommunication.hxx"
#include "paramsParser.hxx"
#include "myIO.hxx"
#include "initGen.hxx"

//template <class gridHandlerType> class snapshot;

template <
  class gridType,
  class slicerType
  >
class gridHandler {
public:

  static std::string getTag() {return "GRID_HANDLER v0.11";}
  template<class gh> friend class snapshot;

  typedef gridType gridT;
  typedef mpiCommunication mpiComT;
  typedef slicerType slicerT;
  //typedef paramsParserType paramsParserT;
 
  static const int P_DIMS= gridType::P_DIMS;
  static const int U_DIMS= gridType::U_DIMS;
  static const int J_DIMS= gridType::J_DIMS;
  static const int DIMS= gridType::DIMS;
  
  typedef typename gridT::value_type dataT;
  typedef typename gridT::scaleT scaleT;
  typedef typename gridT::scaleTypeT scaleTypeT;
  typedef typename gridT::scaleTypeV scaleTypeV;
  typedef typename gridT::paramT paramT;
  typedef typename gridT::dirT dirT;

  typedef typename gridT::iterator gridItT;
  typedef typename gridT::fieldIterator gridFieldItT;

  typedef initGen<gridT> initGenT;

protected:
  std::vector< std::pair<long,dirT> > sendDir;
  //std::vector< std::pair<gridItT,gridItT> > sendIterators;
  std::vector< std::pair<gridFieldItT,gridFieldItT> > sendIterators;
  std::vector< std::vector<dataT> > sendBuffer;

  std::vector< std::pair<long,dirT> > receiveDir;
  //std::vector< std::pair<gridItT,gridItT> > receiveIterators;
  std::vector< std::pair<gridFieldItT,gridFieldItT> > receiveIterators;
  std::vector< std::vector<dataT> > receiveBuffer;
  
  void initSync()
  {
    long i,j,k;

    std::vector <std::vector<dirT> > marginDir(DIMS);
    std::vector<long> buffer((2+DIMS)*hlp::IntPower<3,DIMS>::value);
 
    for (i=0;i<DIMS;i++) 
      {
	marginDir[i].push_back(gridNav::dir(i,0));
	if (grid->getLowMargin(i)>0) marginDir[i].push_back(gridNav::dir(i,-1));
	if (grid->getHighMargin(i)>0) marginDir[i].push_back(gridNav::dir(i,+1));
      }
    
    std::vector<int> w(DIMS,0);
    while (true)
      {
	dirT direction=0;
	int delta=1;
	for (i=0;i<DIMS;i++) 
	  {
	    direction |= marginDir[i][w[i]];
	    w[i]+=delta;
	    if (w[i]>=marginDir[i].size()) {w[i]=0;delta=1;}
	    else delta=0;
	  }

	if (direction) 
	  {	   
	    long nei = slicer.getNeiID(direction);//slicerT::getNeiID(mpiCom->rank(), mpiCom->size(), direction);
	    receiveDir.push_back(std::make_pair(nei,direction));
	   
	    gridItT bit=grid->margin_begin(receiveDir.back().second);
	    gridItT eit=grid->margin_end(receiveDir.back().second);
	    
	    //receiveIterators.push_back(std::make_pair(bit,eit));
	    receiveIterators.push_back(std::make_pair(gridFieldItT(bit),gridFieldItT(eit)));
	  }

	if (delta) break;
      }
   
    receiveBuffer.resize(receiveDir.size());
    for (i=0;i<receiveDir.size();i++)
      receiveBuffer[i].resize(grid->getNFields()*grid->margin_size(receiveDir[i].second));      
      
    for (i=0;i<mpiCom->size();i++)
      {
	if (i==mpiCom->rank())
	  {
	    buffer[0]=receiveDir.size();
	    for (j=0;j<receiveDir.size();j++)
	      {
		buffer[1+(2+DIMS)*j]=(long)receiveDir[j].first;
		buffer[1+(2+DIMS)*j+1]=(long)gridNav::reverse(receiveDir[j].second);
		gridItT it=grid->margin_begin(receiveDir[j].second);
		for (k=0;k<DIMS;k++) buffer[1+(2+DIMS)*j+2+k]=gridItT::dim(it,k);
	      }
	  }

	mpiCom->Bcast(buffer,i);
	//MPI_Bcast(&buffer[0],buffer.size(),MPI_LONG,i,MPI_COMM_WORLD);

	int sendInfoP[P_DIMS];
	int sendInfoU[U_DIMS];
	for (j=0;j<buffer[0];j++)
	  {
	    if (buffer[1+(2+DIMS)*j] == mpiCom->rank())
	      {
		sendDir.push_back(std::make_pair(i,buffer[1+(2+DIMS)*j+1]));
		for (k=0;k<P_DIMS;k++) sendInfoP[k]=buffer[1+(2+DIMS)*j+2+k];
		for (k=0;k<U_DIMS;k++) sendInfoU[k]=buffer[1+(2+DIMS)*j+2+k+P_DIMS];
		//printf("send : %d %d %d %ld\n",sendInfoP[0],sendInfoU[0],sendInfoU[1],sendDir.back().second);	
		gridItT bit=grid->innerMargin_begin(sendInfoP,sendInfoU,sendDir.back().second);
		gridItT eit=grid->innerMargin_end(sendInfoP,sendInfoU,sendDir.back().second);
		//sendIterators.push_back(std::make_pair(bit,eit));
		sendIterators.push_back(std::make_pair(gridFieldItT(bit),gridFieldItT(eit)));
	      }
	  }
      }
   
    sendBuffer.resize(sendDir.size());
    for (i=0;i<sendDir.size();i++)
      sendBuffer[i].resize(grid->getNFields()*gridItT::boxSize(sendIterators[i].first));      
    
    if (debug)
      {
	char tmpS[10000];
	char tmpR[10000];

	strcpy(tmpS,"");
	strcpy(tmpR,"");

	for (i=0;i<sendDir.size();i++) 
	  sprintf(tmpS,"%s (%ld,%ld:%ld)",tmpS,sendDir[i].first,sendDir[i].second,gridItT::boxSize(sendIterators[i].first));
	for (i=0;i<receiveDir.size();i++) 
	  sprintf(tmpR,"%s (%ld,%ld:%ld)",tmpR,receiveDir[i].first,receiveDir[i].second,gridItT::boxSize(receiveIterators[i].first));

	printf("  (%d:) I send data to : %s.\n  (%d:) I receive data from: %s\n",mpiCom->rank(),tmpS,mpiCom->rank(),tmpR);
      }
   
  }

public:

  long getNeighbor(int dim, int dir)
  {
    return slicer.getNeiID(gridNav::dir(dim,dir));
  }

  long getNeighbor(dirT dir)
  {
    return slicer.getNeiID(dir);
  }
  
  void synchronize()
  {
    MPI_Request requestS[sendDir.size()];
    MPI_Request requestR[receiveDir.size()];
    long i;
    
#pragma omp parallel for 
    for (i=0;i<sendDir.size();i++)
      std::copy(sendIterators[i].first,sendIterators[i].second,sendBuffer[i].begin());
     
    for (i=0;i<receiveDir.size();i++)
      mpiCom->Irecv(&receiveBuffer[i][0],receiveBuffer[i].size(),receiveDir[i].first,&requestR[i],0);
    //MPI_Irecv(&receiveBuffer[i][0],receiveBuffer[i].size(),MPI_DOUBLE,receiveDir[i].first,0,MPI_COMM_WORLD,&requestR[i]);
  
    for (i=0;i<sendDir.size();i++)
      mpiCom->Isend(&sendBuffer[i][0],sendBuffer[i].size(),sendDir[i].first,&requestS[i],0);
    //MPI_Isend(&sendBuffer[i][0],sendBuffer[i].size(),MPI_DOUBLE,sendDir[i].first,0,MPI_COMM_WORLD,&requestS[i]);
    
    mpiCom->Waitall(sendDir.size(),requestS);
    mpiCom->Waitall(receiveDir.size(),requestR);
    /*
    MPI_Status status;
    for (i=0;i<sendDir.size();i++)
      MPI_Wait(&requestS[i], &status);
      
    for (i=0;i<receiveDir.size();i++)
      MPI_Wait(&requestR[i], &status);
    */

    //MPI_Barrier(MPI_COMM_WORLD);
#pragma omp parallel for   
    for (i=0;i<receiveDir.size();i++)
      std::copy(receiveBuffer[i].begin(),receiveBuffer[i].end(),receiveIterators[i].first);
 
 
  }

  gridHandler(const paramsParser &params,const paramT &gp, mpiComT &mpiCom_, initGenT *init):
    mpiCom(&mpiCom_),grid(NULL)
  {
    long i;

    nGrids = mpiCom->size();
    fullGridP = gp;

    typename scaleT::scaleTypeSelect scaleSelect;

    for (i=0;i<P_DIMS;i++)
      {
	fullGridP.Pmin[i]=params.template get<>("Pmin",gridT::parserCategory(),fullGridP.Pmin[i],i);
	fullGridP.Pmax[i]=params.template get<>("Pmax",gridT::parserCategory(),fullGridP.Pmax[i],i);
	fullGridP.Pres[i]=params.template get<>("Pres",gridT::parserCategory(),fullGridP.Pres[i],i);
	fullGridP.Pscale[i]=scaleSelect.getVal(params.template get<>("Pscale",gridT::parserCategory(),scaleSelect.getString(fullGridP.Pscale[i]),i));
      }
    for (i=0;i<U_DIMS;i++)
      {
	fullGridP.Umin[i]=params.template get<>("Umin",gridT::parserCategory(),fullGridP.Umin[i],i);
	fullGridP.Umax[i]=params.template get<>("Umax",gridT::parserCategory(),fullGridP.Umax[i],i);
	fullGridP.Ures[i]=params.template get<>("Ures",gridT::parserCategory(),fullGridP.Ures[i],i);
	fullGridP.Uscale[i]=scaleSelect.getVal(params.template get<>("Uscale",gridT::parserCategory(),scaleSelect.getString(fullGridP.Uscale[i]),i));
      }

    fullGridP.haveParentGrid=0;
    fullGridP.parent_P_DIMS=P_DIMS;
    fullGridP.parent_U_DIMS=U_DIMS;  
    for (i=0;i<DIMS;i++) fullGridP.parent_which[i]=i;
   
    //slicerT::setDefaultFullGridMargin(fullGridP);
    slicer.init(fullGridP, mpiCom->rank(), nGrids, num_omp_threads);
    //slicerT::getSubGridParams(fullGridP, mpiCom->rank(), nGrids, num_omp_threads)
    //printf("xres = %d %d\n",slicer.getGridParams().Pres[0],fullGridP.Pres[0]);   
    grid=init->generate(slicer.getGridParams());
    initSync();
    //printf("Fullgrid position (%d): %d %d %d %d\n",mpiCom->rank(),
    //getMyFullGridPos(0),getMyFullGridPos(1),
    //getMyFullGridPos(2),getMyFullGridPos(3));
    
    //for (i=0;i<DIMS;i++) 
    //fullGridCoord[i]=slicer.getMyFullGridCoord(i);
    //mpiCom->barrier();
    //exit(0);
  }

  gridHandler(FILE *f, mpiComT &mpiCom_):
    mpiCom(&mpiCom_),grid(NULL)
  {
    read(f);    
    slicer.init(fullGridP, mpiCom->rank(), nGrids, num_omp_threads);
    //for (i=0;i<DIMS;i++) 
    //fullGridCoord[i]=slicer.getMyFullGridCoord(i);
    initSync();
  }

  ~gridHandler()
  {
    if (grid!=NULL) 
      delete grid;
  }
  
  //gridT *getGrid() {return grid;}
  //void setGrid(gridT *g) {grid=g;}
  gridT *getGrid() {return grid;}

  //const paramT &getGridParams() const {return grid->getParams();}
  const paramT &getFullGridParams() const {return fullGridP;}
  const std::vector<double> &getFullGridCoords(int dir) {return slicer.getFullGrid(dir);}
  int getMyFullGridPos(int dir) {return slicer.getMyFullGridPos(dir);}
  
  bool isFullGridBoundary_Plow(int i) {return (grid->getParams().Pmin[i]<=fullGridP.Pmin[i]);}
  bool isFullGridBoundary_Phigh(int i) {return (grid->getParams().Pmax[i]>=fullGridP.Pmax[i]);}
  bool isFullGridBoundary_Ulow(int i) {return (grid->getParams().Umin[i]<=fullGridP.Umin[i]);}
  bool isFullGridBoundary_Uhigh(int i) {return (grid->getParams().Umax[i]>=fullGridP.Umax[i]);}
  bool isFullGridBoundary_low(int i) {return (i<P_DIMS)?isFullGridBoundary_Plow(i):isFullGridBoundary_Ulow(i-P_DIMS);}
  bool isFullGridBoundary_high(int i) {return (i<P_DIMS)?isFullGridBoundary_Phigh(i):isFullGridBoundary_Uhigh(i-P_DIMS);} 

  template <class subGridT>
  void getSubGrid(subGridT &result,int which[subGridT::DIMS], int nFields, bool global=true, long minStorage=0)
  {
    typename subGridT::paramT subP;

    if (global) 
      getFullGridParams().getSubGridParams(subP,which,minStorage);
    else
      grid->getParams().getSubGridParams(subP,which,minStorage);
    if (nFields>=1) subP.nFields=nFields;
    result.init(subP);
  }


  /*
  void write(const std::string &fname)
  {
    char filename[1024];
    if (mpiCom->size()>1)
      sprintf(filename,"%s_%6.6d.gh",fname.c_str(),mpiCom->rank());
    else
      sprintf(filename,"%s.gh",fname.c_str());
    
    FILE *f=fopen(fname.c_str(),"w");

    if (!(f = fopen(filename,"r")))
    {
      fprintf(stderr,"ERROR: opening file '%s' for writing.\n",fname.c_str());
      exit(-1);
    }

    write(f);  
    fclose(f);
  }
  */

  void write(FILE *f)
  {
    long i;
    int dummy;

    myIO::writeTag(f,getTag());
     
    dummy=1;fwrite(&dummy,sizeof(int),1,f);// used to test for endianness
    int pp[3];
    pp[0]=P_DIMS;pp[1]=U_DIMS;pp[2]=J_DIMS;
    fwrite(pp,sizeof(int),3,f);
    
    i=sizeof(dataT);
    fwrite(&i,sizeof(long),1,f);
    fwrite(&nGrids,sizeof(long),1,f);
    i=mpiCom->rank();
    fwrite(&i,sizeof(long),1,f);
    fullGridP.write(f);
    grid->write(f);

    std::vector<char> spare(256*8,0);
    fwrite(&spare[0],sizeof(char),spare.size(),f);
  }

  template <class iterator>
  std::vector<double> gather_all(const iterator &tab_start, long dir, long addSpare=0)
  {
    int fg_size=slicer.getFullGridSize(dir);
    std::vector<double> result(fg_size+addSpare,0);
    long i;   
    
    std::vector<double> Total_in(slicer.getSlicerSize(dir)+addSpare,0);
    std::vector<double> Total_out(slicer.getSlicerSize(dir)+addSpare,0);
    
    const long r=slicer.getMyCoord(dir);

    const long imin=grid->getLowMargin(dir);
    const long imax=grid->getValCoord(dir).size()-grid->getHighMargin(dir);       

    
    int pos=getMyFullGridPos(dir);
  
    std::copy(tab_start+imin,tab_start+imax,&result[addSpare+pos+imin]);
    mpiCom->Allreduce_inplace(result,MPI_SUM);
    //MPI_Allreduce( MPI_IN_PLACE, &result[0],result.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    
    return result;
  }

  template <class iterator>
  void gather(const iterator &tab_start, long dir, long addSpare=0)
  {
    std::vector<double> tmp=gather_all(tab_start,dir,addSpare);

    int pos=getMyFullGridPos(dir);
    long size = grid->getValCoord(dir).size();
    
    std::copy(&tmp[pos],&tmp[pos+size],tab_start);    
  }


  template <class iterator>
  std::vector<double> accumulate_all(const iterator &tab_start, long dir, long addSpare=0)
  {
    std::vector<double> result=gather_all(tab_start,dir,addSpare);
    long i;
   
    for (i=1;i<result.size();i++) result[i]+=result[i-1];
      
    return result;
  }

  template <class iterator>
  double accumulate(const iterator &tab_start, long dir, long addSpare=0)
  {
    std::vector<double> tmp=accumulate_all(tab_start,dir,addSpare);

    int pos=slicer.getMyFullGridPos(dir);
    long size = grid->getValCoord(dir).size();

    std::copy(&tmp[pos],&tmp[pos+size],tab_start);
    return tmp.back();
  }

  template < class destGridT >
  void gatherGrids(destGridT &dst)
  {
    if (mpiCom->size()==1) return;
   
    const int rank=mpiCom->rank();
    const int size=mpiCom->size();

    int gridPos[size][destGridT::DIMS];
    int gridHM[size][destGridT::DIMS];
    int gridDim[size][destGridT::DIMS];

    if (!dst.findGridPos(*grid,gridPos[rank],gridDim[rank],true))
      {
	fprintf(stderr,"ERROR in gatherGrids: incompatible grids!\n");
	exit(-1);
      }
    //printf("(%d:) size= %d %d\n",rank,gridDim[rank][0],gridDim[rank][1]);
    
    mpiCom->Allgather_inplace(&gridPos[0][0],destGridT::DIMS);
    mpiCom->Allgather_inplace(&gridDim[0][0],destGridT::DIMS);
    
    /*
    MPI_Allgather(&gridPos[rank],destGridT::DIMS,MPI_INT,
		  gridPos,size*destGridT::DIMS,MPI_INT,
		  MPI_COMM_WORLD);
    MPI_Allgather(&gridDim[rank],destGridT::DIMS,MPI_INT,
		  gridDim,size*destGridT::DIMS,MPI_INT,
		  MPI_COMM_WORLD);
    */
    //for (int i=0;i<size;i++) printf("(%d:) size[%d]= %d %d\n",rank,i,gridDim[i][0],gridDim[i][1]);
    long gridSize[size];
    long maxSize=0;
    for (int j=0;j<size;j++)
      {
	gridSize[j]=1;
	for (int i=0;i<destGridT::DIMS;i++) 
	  {
	    gridSize[j]*=gridDim[j][i];
	    gridHM[j][i]=dst.getValCoord(i).size()-(gridPos[j][i]+gridDim[j][i]);
	  }
	if (gridSize[j]>maxSize) maxSize=gridSize[j];
      }
   
    std::vector<double> buffer(maxSize);
   
    for (int i=0;i<size;i++)
      {
	if (rank==i)
	  {
	    /* printf("process %d writing [%d %d][%d %d][%d %d]\n",rank,
		   gridPos[i][0],gridPos[i][1],
		   gridDim[i][0],gridDim[i][1],
		   gridHM[i][0],gridHM[i][1]);*/
	    typename destGridT::iterator it=dst.subbox_begin(gridPos[i],gridHM[i]);	    
	    const typename destGridT::iterator it_end=dst.subbox_end(gridPos[i],gridHM[i]);
	    
	    long j=0;
	    for (;it!=it_end;it++)
	      buffer[j++]=(*it);
	    //printf("process %d writing... done : %ld/%ld sent\n",rank,j,maxSize);
	  }
	
	//if (rank==i) printf("process %d broadcasting ...\n",rank);
	mpiCom->Bcast(buffer,i);
	//printf("process %d done broadcasting ... done.\n",rank);

	
	if (rank!=i)
	  {
	    //printf("process %d reading [%d %d][%d %d]\n",rank,gridPos[i][0],gridPos[i][1],gridHM[i][0],gridHM[i][1]);
	    typename destGridT::iterator it=dst.subbox_begin(gridPos[i],gridHM[i]);
	    //const typename destGridT::iterator it_end=dst.subbox_end(gridPos[i],gridHM[i]);
	    const long N=(long)gridSize[i];
	    long j=0;
	    for (;j<N;j++,it++)
	      {
		(*it)=buffer[j];
	      }
	    //printf("process %d reading... done: %ld/%ld read\n",rank,j,maxSize);
	  }
	
	//mpiCom->barrier();
      }
    
  }

  /*
  double accumulate(std::vector<double> &tab, long dir)
  {
    long i;
    assert(tab.size() == grid->getValCoord(dir).size());
    // std::vector<double> Total_in(mpiCom->size()+1,0);
    // std::vector<double> Total_out(mpiCom->size()+1,0);
    // const long r=mpiCom->rank();
    std::vector<double> Total_in(slicer.getSlicerSize(dir)+1,0);
    std::vector<double> Total_out(slicer.getSlicerSize(dir)+1,0);
    
    const long r=slicer.getMyCoord(dir);

    const long imin=grid->getLowMargin(dir);
    const long imax=grid->getValCoord(dir).size()-grid->getHighMargin(dir);       

    int size=slicer.getFullGridSize(dir);
    int pos=slicer.getFullGridCoords(dir);
    //printf("size = %d, pos=%d+%ld\n",size,pos,tab.size());
    std::vector<double> tmp(size,0);
    std::copy(&tab[imin],&tab[imax],&tmp[pos+imin]);
    MPI_Allreduce( MPI_IN_PLACE, &tmp[0],tmp.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    //for (i=std::max(pos,1);i<pos+tab.size();i++) tmp[i]+=tmp[i-1];
    for (i=1;i<tmp.size();i++) tmp[i]+=tmp[i-1];
    tab.assign(&tmp[pos],&tmp[pos+tab.size()]);

    return tmp.back();
    
    
    // pour du MPI en J seulement
    // faire des groupes de com et sommer dessus pour P et J
    // std::vector<double> tmp(tab.size(),0);
    // MPI_Allreduce( &tab[0], &tmp[0],tab.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    // tab=tmp;
    // for (i=1;i<tab.size();i++) tab[i]+=tab[i-1];
    

    
    // pour du MPI en P seulement 
    
    // for (i=imin+1;i<imax;i++) tab[i]+=tab[i-1];
    // Total_in[r+1]+=tab[imax-1];
    
    // MPI_Allreduce( &Total_in[1], &Total_out[1], slicer.getSize(dir), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
       
    // for (i=1;i<=r;i++) Total_out[i]+=Total_out[i-1];   
    // for (i=imin;i<imax;i++) tab[i]+=Total_out[r];
    
    //if (r>0) for (i=0;i<imin;i++) tab[i]+=Total_out[r-1];
    //for (i=imax;i<tab.size();i++) tab[i]+=tab[i-1];
    
  }
*/
  /*
  void accumulate(std::vector<double> &tab, long dir)
  {
    long i;
    assert(tab.size() == grid->getValCoord(dir).size());
    std::vector<double> Total_in(mpiCom->size()+1,0);
    std::vector<double> Total_out(mpiCom->size()+1,0);
    const long r=mpiCom->rank();

    const long imin=grid->getLowMargin(dir);
    const long imax=grid->getValCoord(dir).size()-grid->getHighMargin(dir);

    for (i=imin+1;i<imax;i++) tab[i]+=tab[i-1];
    
    Total_in[r+1]+=tab[imax-1];

    MPI_Allreduce( &Total_in[1], &Total_out[1], mpiCom->size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    for (i=1;i<=r;i++) Total_out[i]+=Total_out[i-1];
    for (i=imin;i<imax;i++) tab[i]+=Total_out[r];
    }
  */
  //double sum(std::vector<double> &tab, long dir)

  template <class iterator>
  double sum(const iterator &tab_start, long dir)
  {
    long i;
    //assert(tab.size() == grid->getValCoord(dir).size());
    double sum=0;
    const long imin=grid->getLowMargin(dir);
    const long imax=grid->getValCoord(dir).size()-grid->getHighMargin(dir);
    iterator it=tab_start+imin;
    for (i=imin;i<imax;i++,++it) sum+=(*it);
    //for (i=imin;i<imax;i++) sum+=tab[i];
    
    mpiCom->Allreduce_inplace(&sum,1,MPI_SUM);
    //MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return sum;
  }

  double sum(double val)
  {
    double res=val;
    mpiCom->Allreduce_inplace(&res,(long)1,MPI_SUM);
    //MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return res;
  }

public:
  //mpiComT *mpiCom;
private:
  paramT fullGridP;
  
  mpiComT *mpiCom;
  long nGrids;
  gridT *grid;
  slicerT slicer;

  void read(FILE *f)
  {
    long i;
    int dummy[3];
    bool swap=false;
    /*
    FILE *f=fopen(fname.c_str(),"r");

    if (!f)
    {
      fprintf(stderr,"ERROR: File %s does not exist.\n",fname.c_str());
      exit(-1);
    }
    */
    myIO::checkTag(f,getTag());
   
    int ret=fread(dummy,sizeof(int),1,f);
    if ((*dummy)!=1) swap=true;
    myIO::fread(dummy,sizeof(int),3,f,swap);
    if ((dummy[0]!=P_DIMS)||(dummy[1]!=U_DIMS)||(dummy[2]!=J_DIMS))
      {
	fprintf(stderr,"ERROR: container and file have different dimensions.\n");
	exit(-1);
      }
    myIO::fread(&i,sizeof(long),1,f,swap);
    if (i!=sizeof(dataT))
      {
	fprintf(stderr,"ERROR: container and file have different data types.\n");
	exit(-1);
      }
    myIO::fread(&nGrids,sizeof(long),1,f,swap);
    if (nGrids != mpiCom->size())
      {
	fprintf(stderr,"ERROR: this file can only be loaded on a %ld nodes MPI run.\n",nGrids);
	exit(-1);
      }
    myIO::fread(&i,sizeof(long),1,f,swap);
    if (i != mpiCom->rank())
      {
	fprintf(stderr,"ERROR: wrong rank (%ld, should be %d)!!! .\n",i,mpiCom->rank());
	exit(-1);
      }
    
    fullGridP.read(f,swap);
    if (grid!=NULL) delete grid;
    grid=new gridT();
    grid->read(f,swap);
    
    std::vector<char> spare(256*8,0);
    myIO::fread(&spare[0],sizeof(char),spare.size(),f,swap);   
  }
};


#endif
