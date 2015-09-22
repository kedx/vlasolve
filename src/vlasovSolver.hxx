#ifndef __VLASOV_SOLVER_HXX__
#define __VLASOV_SOLVER_HXX__

#include <sys/times.h>

#include "solverParam.hxx"
#include "paramsParser.hxx"
#include "initGen_all.hxx"

#include "myIO.hxx"

template <class gridHandlerType>
class vlasovSolver {
public:

  static const std::string parserCategory() {return "solver";}
  static const std::string getTag() {return "VLASOV_SOLVER v0.13";}

  typedef vlasovSolver<gridHandlerType> myT;

  typedef solverParam solverParamT;
  typedef gridHandlerType gridHandlerT;

  typedef typename gridHandlerT::gridT gridT;
  typedef typename gridHandlerT::mpiComT mpiComT;
  typedef typename gridHandlerT::slicerT slicerT;

  //typedef snapshot<gridHandlerT> snapshotT;
  typedef initGenAll<gridT> initGenAllT;

  //typedef typename gridHandlerT::paramsParserT paramsParserT;
  
  typedef typename initGenAllT::initGenT initGenT;

  static const int P_DIMS= gridT::P_DIMS;
  static const int U_DIMS= gridT::U_DIMS;
  static const int J_DIMS= gridT::J_DIMS;
  static const int PU_DIMS = P_DIMS+U_DIMS-J_DIMS;
  static const int DIMS= gridT::DIMS;

  typedef typename gridT::value_type dataT;
  typedef typename gridT::scaleT scaleT;
  typedef typename gridT::scaleTypeT scaleTypeT;
  typedef typename gridT::scaleTypeV scaleTypeV;
  typedef typename gridT::paramT gridParamT;
  typedef typename gridT::iterator gridItT;
  typedef typename gridT::dirT dirT;
    
  vlasovSolver(mpiComT &mpiCom_):
    mpiCom(&mpiCom_),initialized(false),gh(NULL),time(0),deltaTime(0),stepCount(0)
  {
   
  }

  virtual ~vlasovSolver()
  {
    delete gh;
    //delete snap;
  }

  void init(const paramsParser &params)
  {
    std::string initTypeStr;
    initTypeStr=params.template get<std::string>("type",initGenT::parserCategory(),"UDF");
    initGenT *init = initGenAllT::get(initTypeStr); 

    sp=init->defaultSolverParams();
    sp.G = params.template get<>("G",parserCategory(),sp.G);
    sp.T0 = params.template get<>("T0",parserCategory(),sp.T0);
    sp.Tmax = params.template get<>("Tmax",parserCategory(),sp.Tmax);
    sp.dt = params.template get<>("dt",parserCategory(),sp.dt);
    sp.nSpecies = params.template get<>("nSpecies",parserCategory(),sp.nSpecies);

    sp.statsEvery = params.template get<>("statsEvery",parserCategory(),25);
    sp.snapshotEvery = params.template get<>("snapshotEvery",parserCategory(),500);
    sp.restartEvery = params.template get<>("restartEvery",parserCategory(),5000);
    
    noLinking=params.template get<>("noLinking",parserCategory(),false);

    time = sp.T0;
    deltaTime = sp.dt;

    init->init(params,sp);

    gridParamT gp=init->defaultGridParams();
    gp.nFields = init->getNSpecies() * nSolverFields();
    
    gh = new gridHandlerT(params,gp,*mpiCom,init);
    //grid=gh->getGrid(); 
    
    setup(init,params);  
    gh->synchronize();
    
    for (int i=0;i<sp.nSpecies;i++)
      renormalize(init->getMass(i),i);
    
    delete init;
    
    //snap = new snapshotT(gh);
    initialized=true;
    
  }

  void init(const paramsParser &params, const std::string &fname)
  {
    if (fname=="") 
      {
	init(params);
	return;
      }

    std::string filename;
    filename = myIO::genFname(fname,"vls",mpiCom->rank(), mpiCom->size());

    FILE *f=fopen(filename.c_str(),"r");
    if (!f)
      {
	fprintf(stderr,"ERROR: File %s does not exist.\n",filename.c_str());
	exit(-1);
      }

    read(f);
    fclose(f);
    
    sp.Tmax = params.template get<>("Tmax",parserCategory(),sp.Tmax);
    sp.dt = params.template get<>("dt",parserCategory(),sp.dt);
    sp.statsEvery = params.template get<>("statsEvery",parserCategory(),sp.statsEvery);
    sp.snapshotEvery = params.template get<>("snapshotEvery",parserCategory(),sp.snapshotEvery);
    sp.restartEvery = params.template get<>("restartEvery",parserCategory(),sp.restartEvery);
    //snap = new snapshotT(gh);
    initialized=true;  
  }

  void run()
  {
    struct tms tmsT;
    long deltaT;

    if (!initialized) 
      {
	fprintf(stderr,"ERROR in valsovSolver: init() must be called before runnning the solver.\n");
	exit(-1);
      }
    
    if (!checkTime())
      {
	printf("WARNING: nothing to do, T0 = %.3f >= Tmax = %.3f.\n",time,sp.Tmax);
	return;
      }

    beforeRunning();

    if (mpiCom->rank()==0) printf("Starting at T=%.5f, dt=%g, tmax = %g\n",time,deltaTime,sp.Tmax);fflush(0);  
    
    deltaT=times(&tmsT);
    long stepCount0=stepCount;

    do {
      if (checkTimer(sp.restartEvery,stepCount))
	{
	  char fname[255];
	  sprintf(fname,"restart_%.4f",time);
	  write(fname);
	}
      if (checkTimer(sp.snapshotEvery,stepCount)) 
	{
	  gh->synchronize();
	  snapshot();
	}
      if (checkTimer(sp.statsEvery,stepCount)) 
	{
	  gh->synchronize();
	  statistics(stepCount==0);
	}
     
      if (mpiCom->rank()==0) {printf("\rT=%.4f -> %.4f   ",time,time+deltaTime);fflush(0);}
      makeOneStep(time);
      time+=deltaTime; 
      stepCount++;
    } while (checkTime());

    if (mpiCom->rank()==0) printf("\n");   

    deltaT-=times(&tmsT);
    double T=-(double)deltaT/sysconf(_SC_CLK_TCK);
    if (mpiCom->rank()==0) printf("Computed %ld timesteps in %gs. (%3.3gs / step)",stepCount,T,T/(stepCount-stepCount0));

    afterRunning();
   
  }

  virtual void makeOneStep(double t)=0;
  virtual void renormalize(double norm, int s=-1)=0;
  
protected:
  gridHandlerT *gh;
  //gridT *grid;
  solverParamT sp;
  mpiComT *mpiCom;
  //snapshotT *snap;

  double time;
  double deltaTime;
  long stepCount;
  bool noLinking;
  bool initialized;
  int nSpecies;

  virtual void snapshot()=0;
  virtual void statistics(bool reset=false)=0;
  
  bool checkTimer(long every, long n=-1) 
  {
    if (n<0) n=stepCount;
    return (every>=0)&&((n%every)==0);
  }

  bool checkTime()
  {
    return ((sp.Tmax<=0)||((time<sp.Tmax)&&(fabs(time-sp.Tmax)>deltaTime/100)));
  }

  virtual int nSolverFields()=0;
  
  virtual void setup(initGenT *init,const paramsParser &params)
  { 
  
  }

  virtual void beforeRunning() 
  {
    gh->synchronize();
  }

  virtual void afterRunning() 
  {
    if (checkTimer(sp.snapshotEvery,stepCount)) 
      {
	gh->synchronize();
	snapshot();
      }

    char fname[255];
    sprintf(fname,"final_%.4f",time);
    write(fname);   
  }

  void write(const std::string &fname)
  {
    char filename[1024];
    
    if (mpiCom->size()>1)
      sprintf(filename,"%s_%6.6d.vls",fname.c_str(),mpiCom->rank());
    else
      sprintf(filename,"%s.vls",fname.c_str());
    
    FILE *f=fopen(filename,"w");

    if (!f)
      {
	fprintf(stderr,"ERROR: opening file '%s' for writing.\n",filename);
	exit(-1);
      }
   
    if (verbose) {printf("\nDumping solver state (%s) ... ",filename);fflush(0);}    
    write(f); 
    fclose(f); 
    if (verbose) printf("done.\n");
  }
 

  virtual void write(FILE *f)
  {
    int dummy=1;fwrite(&dummy,sizeof(int),1,f);// used to test for endianness

    myIO::writeTag(f,getTag());         
    sp.write(f);
    fwrite(&time,sizeof(double),1,f);
    fwrite(&deltaTime,sizeof(double),1,f);
    fwrite(&stepCount,sizeof(long),1,f);  
    dummy=(int)noLinking;
    fwrite(&dummy,sizeof(int),1,f);  
    gh->write(f);  

    std::vector<char> spare(256*8,0);
    fwrite(&spare[0],sizeof(char),spare.size(),f);  
  }

  void read(FILE *f)
  {
    bool swap=false;
    int dummy;
    int ret=fread(&dummy,sizeof(int),1,f);
    if ((dummy)!=1) swap=true;
    
    read(f,swap);
  }

  virtual void read(FILE *f,bool swap)
  {
    myIO::checkTag(f,getTag());
    sp.read(f,swap);
    myIO::fread(&time,sizeof(double),1,f,swap);
    myIO::fread(&deltaTime,sizeof(double),1,f,swap);
    myIO::fread(&stepCount,sizeof(long),1,f,swap);
    int dummy;
    myIO::fread(&dummy,sizeof(int),1,f,swap);
    noLinking=(bool)dummy;
    gh = new gridHandlerT(f, *mpiCom);
    //grid=gh->getGrid();

    std::vector<char> spare(256*8,0);
    myIO::fread(&spare[0],sizeof(char),spare.size(),f,swap);
  }


};

#endif
