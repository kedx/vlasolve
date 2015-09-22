#ifndef __VLASOV_SOLVER_NN_HXX__
#define __VLASOV_SOLVER_NN_HXX__

#ifdef HAVE_FFTW3

#include <algorithm>
#include "global.h"
#include "vlasovSolver.hxx"
#include "dimTraits.hxx"
#include "NDfield.hxx"
#include "typeSelect.hxx"
#include "helpers.hxx"
#include "grid_base.hxx"

#include "quadratureNN.hxx"
#include "mpiRequests.hxx"
#include "poissonNN.hxx"

using namespace NDF;

template <class gridHandlerType>
class vlasovSolverNN : public vlasovSolver<gridHandlerType> {
  public:

  typedef vlasovSolver<gridHandlerType> bT;
  typedef vlasovSolverNN<gridHandlerType> myT;

  typedef typename bT::gridHandlerT gridHandlerT;
  typedef typename bT::gridT gridT;
  typedef typename bT::mpiComT mpiComT;
  typedef typename bT::slicerT slicerT;
  typedef typename bT::initGenAllT initGenAllT;

  static const int P_DIMS= bT::P_DIMS;
  static const int U_DIMS= bT::U_DIMS;
  static const int J_DIMS= bT::J_DIMS;
  static const int PU_DIMS = bT::PU_DIMS;
  static const int DIMS= bT::DIMS;

  typedef typename bT::initGenT initGenT;

  typedef typename bT::dataT dataT;
  typedef typename bT::scaleT scaleT;
  typedef typename bT::scaleTypeT scaleTypeT;
  typedef typename bT::scaleTypeV scaleTypeV;
  typedef typename bT::gridParamT gridParamT;
  typedef typename bT::dirT dirT;
  typedef typename scaleT::valLocationV valLocationV;
  typedef typename gridT::valLocationTr valLocationTr;

  typedef gridBase< dimTraits<P_DIMS,0>, dataT, valLocationTr > projectedPSlicesT;
  typedef typename projectedPSlicesT::paramT projectedPSlicesParamT;
  //typedef gridBase< dimTraits<0,U_DIMS>, dataT, gridT::valLocationTr > porjectedUSlicesT;

  typedef typename gridT::sliceIterator gridSliceItT;
  typedef mpiValAtCoordsRequestT<typename gridSliceItT::pointer, mpiComT, DIMS> requestT;

  typedef quadratureNN<gridHandlerT> quadratureT;  
  
  typedef poissonFFTSolver<dataT,P_DIMS> poissonFFTSolverT;

  //static const std::string getTag() {std::string out("VLASOV_SOLVER"); return out<<P_DIMS<<" v0.10";}
  static const std::string getTag() {char tmp[255];sprintf(tmp,"%s%d%s","VLASOV_SOLVER_",PU_DIMS," v0.10");return std::string(tmp);}

  vlasovSolverNN(mpiComT &mpiCom_) : 
    bT(mpiCom_),poissonSolver(),quickInit(false)
  {
    //static_checkCompatibility<P_DIMS,U_DIMS>();
  }

  virtual ~vlasovSolverNN()
  {
    
  } 

  // WARNING! buffer must be allocated for DIMS values (not P_DIMS !!)
  bool getDrift(const gridSliceItT &it, double *buffer)
  {
    bool out=false;
    it.C(buffer);
        
    for (int i=0;i<P_DIMS;i++)
      {
	buffer[i]-=buffer[i+P_DIMS]*drift_dt;
	if ((buffer[i]<Pmin[i])||(buffer[i]>Pmax[i])) out=true;
      }
    
    return out;
  }

  bool updateDrift(const gridSliceItT &it, double *buffer)
  {
    bool out=false;
    it.CP(buffer);
    //bool pr=false;
    //if ((buffer[0]<-1.7)||(buffer[1]<-1.7)) pr=true;
    //if (pr) printf("[%g %g %g %g] ->",buffer[0],buffer[1],buffer[2],buffer[3]);
    for (int i=0;i<P_DIMS;i++)
      {
	buffer[i]-=buffer[i+P_DIMS]*drift_dt;
	if ((buffer[i]<Pmin[i])||(buffer[i]>Pmax[i])) out=true;
      }
    //if (pr) printf("[%g %g %g %g] (%d)\n",buffer[0],buffer[1],buffer[2],buffer[3],(int)(out));
    return out;
  }

 protected:
  projectedPSlicesT density;
  projectedPSlicesT forcePoisson[P_DIMS];
  projectedPSlicesT force;

  //for communication with surrounding grids
  requestT myRequests;
  requestT neiRequests;
  double drift_dt;
  poissonFFTSolverT poissonSolver;

  void setupRequests()
  {   
    gridT *grid=bT::gh->getGrid();
    
    myRequests.clear();
    neiRequests.clear();

    if (bT::mpiCom->size()<2) return;
    
    long i,j;
    const long nSlices=grid->nSlices();

    // no openmp here
    for (j=0;j<nSlices;j++)
      {
	const gridSliceItT it_end=grid->slice_end(j);
	gridSliceItT it=grid->slice_begin(j);

	double disp_drift[DIMS];
	getDrift(it,disp_drift);

	for (;it!=it_end;++it)
	  {
	    dirT dir=0;
	    updateDrift(it,disp_drift);

	    for (int i=0;i<P_DIMS;i++)
	      {
		if (disp_drift[i]>P_Val[i].back()) 
		  dir+=gridNav::dir(i,1);
		else if (disp_drift[i]<P_Val[i].front())
		  dir+=gridNav::dir(i,-1);
	      }	  
	    //printf("HOP\n");
	    int neiID=bT::gh->getNeighbor(dir);
	    //printf("BOP\n");
	    if (neiID<0) continue;
	    else 
	      {
		myRequests.procID.push_back(neiID);
		myRequests.ptr.push_back(&(*it));
		myRequests.coords.insert(myRequests.coords.end(),disp_drift,disp_drift+DIMS);
		myRequests.coords[myRequests.coords.size()-1]=j;
		//myRequests.coords.push_back(disp_driftR[i]);
	        //myRequests.coords.push_back(disp_driftU[i]);
		//myRequests.coords.push_back((double)j);
	      }
	  }
      }
    
    myRequests.setupLocal(grid->getNFields());  
    neiRequests.setupRemote(myRequests,bT::mpiCom);
    
    //  char txt[255];
    //  strcpy(txt,"");
    //  for (i=0;i<neiRequests.nProcs;i++) sprintf(txt,"%s%d ",txt,neiRequests.procID[i]);  
    //  printf("OK! rank %d sends to %d procs: %s\n",bT::mpiCom->rank(),neiRequests.nProcs,txt);
    //  strcpy(txt,"");
    
    //  for (i=0;i<myRequests.nProcs;i++) sprintf(txt,"%s%d ",txt,myRequests.procID[i]);  
    // printf("OK! rank %d receive from %d procs: %s\n",bT::mpiCom->rank(),myRequests.nProcs,txt);
    
  }


  
 private:
  void setup()
  { 
    gridT *grid=bT::gh->getGrid();
    const gridParamT &gp=bT::gh->getGrid()->getParams();
    const gridParamT &fgp=bT::gh->getFullGridParams();
    
    for (int i=0;i<P_DIMS;i++)
      {
	P_C[i]=grid->getCellCoord_P(i);
	U_C[i]=grid->getCellCoord_U(i);
	P_V[i]=grid->getVertCoord_P(i);
	U_V[i]=grid->getVertCoord_U(i);
	P_Val[i]=grid->getValCoord_P(i);
	U_Val[i]=grid->getValCoord_U(i);

	Pscale[i] = gp.Pscale[i];
	Uscale[i] = gp.Uscale[i];
	Pmin[i]=P_Val[i].front();
	Pmax[i]=P_Val[i].back();
	Umin[i]=U_Val[i].front();
	Umax[i]=U_Val[i].back();
      }
    drift_dt=(bT::noLinking)?bT::deltaTime*0.5:bT::deltaTime;

    int which[P_DIMS];
    for (int i=0;i<P_DIMS;i++) which[i]=i;
       
    bT::gh->template getSubGrid<projectedPSlicesT>(density,which,1,true);
    //bT::gh->template getSubGrid<projectedPSlicesT>(force,which,P_DIMS);    

    density.setMinStorage(poissonFFTSolverT::inplaceArraySize(density.getArrDims()));    

    poissonSolver.setG(bT::sp.G);
    poissonSolver.setQuickInit(quickInit);
    poissonSolver.setup(density.getArrDims(),density.getBoxSize(),
			density.getDataPtr(),density.getDataPtr(),
			0*poissonFFTSolverT::PERIODIC);

    typename projectedPSlicesT::paramT p=density.getParams(); 
    for (int i=0;i<P_DIMS;i++)
      forcePoisson[i].init(p,poissonSolver.getForcePtr(i));    
    
    /*
      typename projectedPSlicesT::paramT p=density.getParams(); 
      p.minStorageSize=poissonFFTSolverT::inplaceArraySize(density.getArrDims());
      forcePoisson[0].init(p);    
      poissonSolver.setup(density.getArrDims(),density.getBoxSize(),
      density.getDataPtr(),forcePoisson[0].getDataPtr(),
      0*poissonFFTSolverT::PERIODIC);
    
      p.minStorageSize=0;
      for (int i=1;i<P_DIMS;i++)
      forcePoisson[i].init(p,poissonSolver.getForcePtr(i));    
    */
    setupRequests();
    
  }

protected:
  bool quickInit;

  std::vector<double> P_C[P_DIMS];
  std::vector<double> U_C[U_DIMS];
  std::vector<double> P_V[P_DIMS];
  std::vector<double> U_V[U_DIMS];
  std::vector<double> P_Val[P_DIMS];
  std::vector<double> U_Val[U_DIMS];

  scaleTypeT Pscale[P_DIMS];
  scaleTypeT Uscale[P_DIMS];

  double Pmin[P_DIMS];
  double Pmax[P_DIMS];
  double Umin[P_DIMS];
  double Umax[P_DIMS];

  virtual void renormalize(double norm, int which)
  {
    if (norm<=0) return;
    gridT *grid=bT::gh->getGrid();
    double mass;
    
    if (which<0)
      mass=quadratureT::integrate(bT::gh,intNN::f<P_DIMS>(0,bT::sp.nSpecies),true);  
    else
      mass=quadratureT::integrate(bT::gh,intNN::f<P_DIMS>(which,1),true);  
    //printf("mass = %lg\n",mass);
    double fac=norm/mass;

    if (!bT::mpiCom->rank()) 
      {
	if (which<0) 
	  printf("Renormalizing mass : M=%g -> %g.\n",mass,norm);
	else
	  printf("Renormalizing mass : M[%d]=%g -> %g.\n",which,mass,norm);
      }

    const long nSlices=grid->nSlices();
    long j;
    long k=0;

    if (which<0)
      {
#pragma omp parallel for
	for (j=0;j<nSlices;j++)
	  {
	    const field_iterator<gridSliceItT> it_end(grid->slice_end(j));
	    field_iterator<gridSliceItT> it(grid->slice_begin(j));
	    for (;it!=it_end;++it)
	      (*it)*=fac;
	  }
      }
    else
      {
#pragma omp parallel for
	for (j=0;j<nSlices;j++)
	  {
	    const oneField_iterator<gridSliceItT> it_end(grid->slice_end(j),which);
	    oneField_iterator<gridSliceItT> it(grid->slice_begin(j),which);
	    for (;it!=it_end;++it)
	      (*it)*=fac;
	  }
      }
    
  }

  virtual void beforeRunning()
  {
    bT::beforeRunning();
    
  }
  
  virtual void setup(initGenT *init,const paramsParser &params)
  {
    bT::setup(init,params);
    quickInit = params.template get<>("quickInit",bT::parserCategory(),false);
    setup();
  }

  virtual void read(FILE *f, bool swap)
  {
    bT::read(f,swap);
    myIO::checkTag(f,getTag());
    
    std::vector<char> spare(256*8,0);
    myIO::fread(&spare[0],sizeof(char),spare.size(),f,swap);
 
    setup();
  }

  virtual void write(FILE *f)
  {
    bT::write(f);

    myIO::writeTag(f,getTag());
    std::vector<char> spare(256*8,0);
    fwrite(&spare[0],sizeof(char),spare.size(),f);
  }

  virtual void snapshot()
  {
    static long nsnap=0;
    printf("\n");
    /*
    char fname[1024];
    sprintf(fname,"snapshot_%4.4f",bT::time);
 
    if (bT::mpiCom->size()>1)
      sprintf(fname,"%s_%6.6d.ND",fname,bT::mpiCom->rank());
    else
      sprintf(fname,"%s.ND",fname);
    
    NDfield *f=bT::gh->getGrid()->toNDfield();
    Save_NDfield(f,fname);
    free(f);
    
    sprintf(fname,"density_%4.4f",bT::time);
 
    if (bT::mpiCom->size()>1)
      sprintf(fname,"%s_%6.6d.ND",fname,bT::mpiCom->rank());
    else
      sprintf(fname,"%s.ND",fname);
    
    f=density.toNDfield();
    Save_NDfield(f,fname);
    free(f);
    
    sprintf(fname,"force_%4.4f",bT::time);
 
    if (bT::mpiCom->size()>1)
      sprintf(fname,"%s_%6.6d.ND",fname,bT::mpiCom->rank());
    else
      sprintf(fname,"%s.ND",fname);
    
    f=force.toNDfield();
    Save_NDfield(f,fname);
    free(f);
    */
    char dir[][2]={"X","Y","Z"};
    for (int i=0;i<=P_DIMS;i++)
      {
	char fname[1024];
	NDfield *f=NULL;
	if (i==P_DIMS)
	  {
	    sprintf(fname,"snapshot_%4.4f",bT::time);
 
	    if (bT::mpiCom->size()>1)
	      sprintf(fname,"%s_%6.6d.ND",fname,bT::mpiCom->rank());
	    else
	      sprintf(fname,"%s.ND",fname);
    
	    NDfield *f=bT::gh->getGrid()->toNDfield();
	    Save_NDfield(f,fname);
	    free(f);
	  }
	else if (nsnap)
	  {
	    sprintf(fname,"force%s_%4.4f",dir[i],bT::time);
 
	    if (bT::mpiCom->size()>1)
	      sprintf(fname,"%s_%6.6d.ND",fname,bT::mpiCom->rank());
	    else
	      sprintf(fname,"%s.ND",fname);
    
	    f=forcePoisson[i].toNDfield();
	    Save_NDfield(f,fname);
	    free(f);
	  }
      }
    nsnap++;
  }

  double potentialEnergy(double &Ep)
  {
    return 0;
    printf("WARNING potentialEnergy: NOT IMPLEMENTED\n");
  }

  
  virtual void statistics(bool reset)
  {
    gridT *grid=bT::gh->getGrid();
    long i;
    char fname[1024];
    int nsp=bT::sp.nSpecies;
    
    double Ep_tot=0;
    double M_tot=0;
    double Ek_tot=0;
    double E_tot=0;
    double L1_tot=0;

    for (int s=0;s<nsp;s++)
      {
	if (nsp>1)
	  sprintf(fname,"statistics_S%d.txt",s);
	else 
	  sprintf(fname,"statistics.txt");
     
	double Ep=0;
	double M=potentialEnergy(Ep);
	double Ek=quadratureT::integrate(bT::gh,intNN::kineticEnergy<P_DIMS>(),true);
	double E=Ek+Ep;
	double S=quadratureT::integrate(bT::gh,intNN::entropy<P_DIMS>(),true);  
	double L1=quadratureT::integrate(bT::gh,intNN::L1Norm<P_DIMS>(),true);  
	double L2=quadratureT::integrate(bT::gh,intNN::L2Norm<P_DIMS>(),true);

	Ep_tot+=Ep;
	M_tot+=M;
	Ek_tot+=Ek;
	E_tot+=E;
	L1_tot+=L1;

	if (bT::mpiCom->rank()!=0) return;
	FILE *f;
	if (reset) 
	  {
	    f=fopen(fname,"w");
	    fprintf(f,"# Time Step Mass E Ek Ep S L1 L2\n");
	  }
	else f=fopen(fname,"a");

	fprintf(f,"%.10lg %ld %10.10lg %10.10lg %10.10lg %10.10lg %10.10lg %10.10lg %10.10lg\n",
		bT::time,bT::stepCount,M,E,Ek,Ep,S,L1,L2);
    
	fclose(f);
      }

    if (nsp>1)
      {
	sprintf(fname,"statistics.txt");
	double S_tot=quadratureT::integrate(bT::gh,intNN::entropy<P_DIMS>(0,nsp),true);
	double L2_tot=quadratureT::integrate(bT::gh,intNN::L2Norm<P_DIMS>(0,nsp),true);
      
	if (bT::mpiCom->rank()!=0) return;
	FILE *f;
	if (reset) 
	  {
	    f=fopen(fname,"w");
	    fprintf(f,"# Time Step Mass E Ek Ep S L1 L2\n");
	  }
	else f=fopen(fname,"a");

	fprintf(f,"%.10g %ld %10.10g %10.10g %10.10g %10.10g %10.10g %10.10g %10.10g\n",bT::time,bT::stepCount,M_tot,E_tot,Ek_tot,Ep_tot,S_tot,L1_tot,L2_tot);
    
	fclose(f);
      }
  }
  
private:  
  /*
  NDfield *NDfieldFromData(int sliceAtJ=-1)
  {
    long i;
    int useSpecies=(bT::sp.nSpecies>1)?1:0;
    int dims[P_DIMS+U_DIMS+useSpecies];
    double x0[P_DIMS+U_DIMS+useSpecies];
    double delta[P_DIMS+U_DIMS+useSpecies];
    int ndims;
    long stride[P_DIMS+U_DIMS+1+useSpecies];
    int index;
    dataT *d;
    char comment[80];
    gridT *grid=bT::gh->getGrid();
    const gridParamT &gp=bT::gh->getGrid()->getParams();

    stride[0]=1;    
    if (useSpecies)
      {
	dims[0]=bT::sp.nSpecies;
	x0[0]=0;
	delta[0]=bT::sp.nSpecies;
	stride[1]=stride[0]*dims[0];
      }

    long dp=useSpecies;
    for (i=0;i<DIMS;i++) {
      dims[i+dp]=grid->getValCoord(i).size();
      x0[i+dp]=grid->getValCoord(i).front();
      delta[i+dp]=grid->getValCoord(i).back()-grid->getValCoord(i).front();
      stride[i+1+dp]=grid->getValStride(i+1);
    }
    
    strcpy(comment,"");
    for (i=0;i<P_DIMS;i++) 
      {
	char tmp[256];
	sprintf(tmp,"%d %d ",Pscale[i],(int)gridT::valLocation_P);
	strcat(comment,tmp);
      }
    for (i=0;i<U_DIMS;i++) 
      {
	char tmp[256];
	sprintf(tmp,"%d %d ",Uscale[i],(int)gridT::valLocation_U);
	strcat(comment,tmp);
      }
    
    if (sliceAtJ<0) {
      d = grid->getDataPtr();
      ndims=P_DIMS+U_DIMS+useSpecies;
    }
    else {
      printf("NDfieldfromData: NOT IMPLEMENTED for slice>=0\n");exit(0);
      //ndims=P_DIMS+U_DIMS-1+useSpecies;
      //d=grid->fullSlicePtr(sliceAtJ);
    }
    
    return Create_NDfield(dims,ndims,0,ND_DOUBLE,x0,delta,d,comment);
  }
  */
  /*
  template <int NP, int NU>
  void static_checkCompatibility()
  {    
    gridDimensionsAreCompatible< hlp::isEqual<NP,NU>::result >(dimTraits<NP,NU>());
  }
  
  //template <int NP, int NU,bool dimEq> 
  //void gridDimensionsAreCompatible(dimTraits<NP,NU>);
  //template <int NP, int NU> 
  //void gridDimensionsAreCompatible<NP,NU,true>(dimTraits<NP,NU>){};
  
  template <bool dims_are_compatible> 
  void gridDimensionsAreCompatible();
  template <> 
  void gridDimensionsAreCompatible<true>(){};
  */
};

#endif

#endif
