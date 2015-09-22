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

  static std::string getTag() {return "GRID_HANDLER v0.1";}

  template<class gh> friend class snapshot;

  typedef gridType gridT;
  typedef mpiCommunication mpiComT;
  typedef slicerType slicerT;
  //typedef paramsParserType paramsParserT;
 
  static const int P_DIMS= gridType::P_DIMS;
  static const int U_DIMS= gridType::U_DIMS;
  static const int J_DIMS= gridType::J_DIMS;
  static const int DIMS= gridType::DIMS;
  
  typedef typename gridT::dataT dataT;
  typedef typename gridT::scaleT scaleT;
  typedef typename gridT::scaleTypeT scaleTypeT;
  typedef typename gridT::scaleTypeV scaleTypeV;
  typedef typename gridT::paramT paramT;

  typedef initGen<gridT> initGenT;

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
							     

	  //fullGridP.Pscale[i]=scaleT::str2Type(params.template get<>("Pscale",gridT::parserCategory(),scaleT::type2Str(fullGridP.Pscale[i]),i));
	//fullGridP.lowPmargin[i]=slicerTr::defaultPMargin(i);
      }
    for (i=0;i<U_DIMS;i++)
      {
	fullGridP.Umin[i]=params.template get<>("Umin",gridT::parserCategory(),fullGridP.Umin[i],i);
	fullGridP.Umax[i]=params.template get<>("Umax",gridT::parserCategory(),fullGridP.Umax[i],i);
	fullGridP.Ures[i]=params.template get<>("Ures",gridT::parserCategory(),fullGridP.Ures[i],i);
	fullGridP.Uscale[i]=scaleSelect.getVal(params.template get<>("Uscale",gridT::parserCategory(),scaleSelect.getString(fullGridP.Uscale[i]),i));
	  //fullGridP.Uscale[i]=scaleT::str2Type(params.template get<>("Uscale",gridT::parserCategory(),scaleT::type2Str(fullGridP.Uscale[i]),i));
	//fullGridP.Umargin[i]=slicerTr::defaultUMargin(i);
      }

    slicerT::setFullGridMargin(fullGridP);
    grid=init->generate(slicerT::getSubGridParams(fullGridP, mpiCom->rank(), nGrids));
    //gridP = slicerT::getSubGridParams(fullGridP, mpiCom->rank(), nGrids); 
    //grid=new gridT(gridP);   
    //grid=NULL;
  }

  gridHandler(FILE *f, mpiComT &mpiCom_):
    mpiCom(&mpiCom_),grid(NULL)
  {
    read(f);
  }

  ~gridHandler()
  {
    if (grid!=NULL) 
      delete grid;
  }
  
  //gridT *getGrid() {return grid;}
  //void setGrid(gridT *g) {grid=g;}
  gridT *getGrid() {return grid;}

  const paramT &getGridParams() const {return grid->getParams();}
  const paramT &getFullGridParams() const {return fullGridP;}

  bool isFullGridBoundary_Plow(int i) {return (grid->getParams().Pmin[i]<=fullGridP.Pmin[i]);}
  bool isFullGridBoundary_Phigh(int i) {return (grid->getParams().Pmax[i]>=fullGridP.Pmax[i]);}
  bool isFullGridBoundary_Ulow(int i) {return (grid->getParams().Umin[i]<=fullGridP.Umin[i]);}
  bool isFullGridBoundary_Uhigh(int i) {return (grid->getParams().Umax[i]>=fullGridP.Umax[i]);}
  bool isFullGridBoundary_low(int i) {return (i<P_DIMS)?isFullGridBoundary_Plow(i):isFullGridBoundary_Ulow(i-P_DIMS);}
  bool isFullGridBoundary_high(int i) {return (i<P_DIMS)?isFullGridBoundary_Phigh(i):isFullGridBoundary_Uhigh(i-P_DIMS);} 
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
  }

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
    for (i=1;i<r;i++) Total_out[r]+=Total_out[r-1];
    for (i=imin;i<imax;i++) tab[i]+=Total_out[r];
  }

private:
  paramT fullGridP;
  //paramT gridP;
  mpiComT *mpiCom;
  long nGrids;
  gridT *grid;

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
	fprintf(stderr,"ERROR: wrong rank !!! .\n",nGrids);
	exit(-1);
      }

    fullGridP.read(f,swap);
    if (grid!=NULL) delete grid;
    grid=new gridT();
    grid->read(f,swap);
    //gridP = grid->getParams();
    //fclose(f);
  }
};


#endif
