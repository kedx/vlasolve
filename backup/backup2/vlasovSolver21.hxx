#ifndef __VLASOV_SOLVER_21_HXX__
#define __VLASOV_SOLVER_21_HXX__

#include "vlasovSolver.hxx"
#include "cellChain.hxx"
#include "dimTraits.hxx"
#include "find_unordered_maps.hxx"
#include "NDfield.hxx"
#include "typeSelect.hxx"

#include "quadrature21.hxx"

using namespace NDF;

template <class gridHandlerType>
class vlasovSolver21 : public vlasovSolver<gridHandlerType> {
public:

  static const std::string getTag() {return "VLASOV_SOLVER21 v0.11";}

  typedef vlasovSolver<gridHandlerType> bT;

  typedef typename bT::gridHandlerT gridHandlerT;
  typedef typename bT::gridT gridT;
  typedef typename bT::mpiComT mpiComT;
  typedef typename bT::slicerT slicerT;
  //typedef typename bT::paramsParserT paramsParserT;
  typedef typename bT::initGenAllT initGenAllT;

  static const int P_DIMS= bT::P_DIMS;
  static const int U_DIMS= bT::U_DIMS;
  static const int J_DIMS= bT::J_DIMS;
  static const int DIMS= bT::DIMS;

  typedef typename bT::initGenT initGenT;

  typedef typename bT::dataT dataT;
  typedef typename bT::scaleT scaleT;
  typedef typename bT::scaleTypeT scaleTypeT;
  typedef typename bT::scaleTypeV scaleTypeV;
  typedef typename bT::gridParamT gridParamT;
  typedef typename bT::gridSliceItT gridSliceItT;

  typedef quadrature21<gridHandlerT> quadratureT;

  vlasovSolver21(mpiComT &mpiCom_) : 
    bT(mpiCom_)
  {
    static_checkCompatibility<P_DIMS,U_DIMS>();
  }

  virtual ~vlasovSolver21()
  {
    
  }  

protected:

  struct kernelTypeV {
    enum type {DELAYED=0, REFLECTIVE=1, UNDEFINED=-1};
  };

  typedef typename kernelTypeV::type kernelTypeT;

  struct kernelTypeSelect : public typeSelect<kernelTypeV> {
    kernelTypeSelect()
    {
      insert("delayed",kernelTypeV::DELAYED);
      insert("reflective",kernelTypeV::REFLECTIVE);
    }
    std::string name() {return "kernel type";}
  };
  kernelTypeT kernelType;
  
  struct dummyCellT
  {
    typedef dataT valT;

    double R;
    double U;
    valT val;

    dummyCellT(double R_=0, double U_=0,valT val_=0):
      R(R_),U(U_),val(val_)
    {}
   
    ~dummyCellT() {}

    void read(FILE *f, bool swap)
    {
      myIO::fread(&R,sizeof(double),1,f,swap);
      myIO::fread(&U,sizeof(double),1,f,swap);
      myIO::fread(&val,sizeof(valT),1,f,swap);
    }

    void write(FILE *f)
    {
      fwrite(&R,sizeof(double),1,f);
      fwrite(&U,sizeof(double),1,f);
      fwrite(&val,sizeof(valT),1,f);
    }

    void print()
    {
      printf("(R=%e U=%e V=%e)",R,U,val);
    }
  };

  typedef cellChain<dummyCellT> cellChainT;
  typedef std::list<cellChainT> cellChainListT;
  typedef typename cellChainListT::iterator cellChainListItT;
  typedef typename my_unordered_map<long,cellChainListItT>::type kernelCellsT;
  typedef typename kernelCellsT::iterator kernelCellsItT;

  cellChainListT allChains;  
  std::vector<kernelCellsT> kernelBoundary;
  double kernelMass;
  
  cellChainListItT newKernelChain(long j, long index)
  {
    cellChainListItT it=allChains.insert(allChains.end(),cellChainT());
   
    if (!kernelBoundary[j].insert(std::make_pair(index,it)).second)
      {

	allChains.erase(it);
	return allChains.end();
      }

    return it;
  }

  std::vector<double> disp_driftR;
  std::vector<double> disp_driftU;  

  void generateDisplacement(double dt, bool reflectiveKernel=true)
  {
    long j;
    const long nSlices=bT::grid->nSlices();
    const double dth=dt/2;

    double R0=(reflectiveKernel)?Rmin:-1;

    //FILE *f2=fopen("testUV.i","wb");
    //fprintf(f2,"uv=array(double,[3,4,%ld,%ld]);",nSlices,bT::grid->getNVal()/nSlices);

    //if (!bT::grid->getRegisteredValue("INIGEN_NormFactor",normFactor)) normFactor=1;
    disp_driftR.resize(bT::grid->getNVal());
    disp_driftU.resize(bT::grid->getNVal());
    
    #pragma omp parallel for
    for (j=0;j<nSlices;j++)
      {
	std::pair<double,double> newVal;
	const gridSliceItT it_end=bT::grid->slice_end(j);
	gridSliceItT it=bT::grid->slice_begin(j);
	double J=it.get_J();
	//fprintf(f2,"id=1;\n");
	for (;it!=it_end;++it)
	  {
	    newVal=advanceParticleDrift(it.get_R(),it.get_U(),J,dth,R0);
	    disp_driftR[it.get_i()]=newVal.first;
	    disp_driftU[it.get_i()]=newVal.second;
	    
	    //fprintf(f2,"uv(,%ld,id++)=[%g,%g,%g,%g];\n",j+1,it.get_R(),it.get_U(),disp_driftR[it.get_i()],disp_driftU[it.get_i()]);
	  }
      }
    //fclose(f2);
  }

  void setupBoundaryConditions(double dt, initGenT *init)
  {
    //double normFactor=1;
    bool isLowB=bT::gh->isFullGridBoundary_Plow(0);
    long j;
    const double dth=dt/2;
    //if (!bT::grid->getRegisteredValue("INIGEN_NormFactor",normFactor)) normFactor=1;
    allChains.clear();
    kernelBoundary.clear();
    kernelBoundary.resize(J_Val.size());

    if (kernelType!=kernelTypeV::DELAYED) return;

    const long nSlices=bT::grid->nSlices();
    // no openMP cause cellchains are not thread safe ...
    for (j=0;j<nSlices;j++)
      {
	std::pair<double,double> newVal;

	const gridSliceItT sit_end=bT::grid->slice_end(j);
	gridSliceItT sit=bT::grid->slice_begin(j);
	double J=sit.get_J();
	for (;sit!=sit_end;++sit)
	  {
	    const long i=sit.get_i();

	    if (disp_driftR[i]<Rmin)
	      {
		if (!isLowB)
		  {
		    fprintf(stderr,"ERROR: time step is too high for this grid.\n");
		    fprintf(stderr,"    R=%g => R_drift=%g<%g.\n",sit.get_R(),disp_driftR[i],R_Val[0]);
		    exit(-1);
		  }

		double deltaT=dth;
		cellChainListItT cit = newKernelChain(j,i);		    
		std::pair<double,double> newVal=std::make_pair(disp_driftR[i],disp_driftU[i]);
		std::vector<double> pos(3);
		pos[2]=J;
		while (newVal.first<R_Val[0])
		  {
		    pos[0]=newVal.first;pos[1]=newVal.second;
		    dataT val=init->valueAt(pos);
			
		    cit->push(dummyCellT(newVal.first,newVal.second,val));
		    deltaT+=dth;
		    newVal=advanceParticleDrift(sit.get_R(),sit.get_U(),J,deltaT);
		  }	
		cit->push(dummyCellT(newVal.first,newVal.second,0));
	      }
	  }
      }
   
    printf("Nchains : %ld\n",allChains.size());
  }

  std::pair<double,double> advanceParticleDrift(double R, double U, double J,double dt, double R0=-1)
  {
    double H0; 
    double newU;
    double newR;
    double tlim,tR0;
    double sgn;
    double deltaT;
    
    H0= 0.5* (U*U + (J*J)/(R*R));
    sgn= (U<=0)?-1:1;
    deltaT = dt;		
    newR=R;
    newU=U;

    if (sgn<0) {
      tlim=10*dt;
      tR0=10*dt;
    }
    else {
      tlim=sqrt(2*R*R*H0-J*J)/(2*H0);
      if (R0>0) {
	double A,B,C,D;
	A=2*H0;
	B=-2*sgn*sqrt(2*R*R*H0-J*J);
	C=(R*R-R0*R0);
	D=B*B-4*A*C;
	if (D<0) tR0=10*dt;
	else tR0=(-B-sqrt(D))/(2*A);
	if (tR0<0) {
	  tR0=(-B+sqrt(D))/(2*A);
	  if (tR0<0) tR0=10*dt;
	}		  
      } else tR0=10*dt;
    }	

    
    if ((tR0<deltaT)&&(tR0<tlim))
      {
	newR=sqrt(pow(fabs(newR*newU)-2*sgn*H0*tR0,2)+J*J)*pow(2*H0,-0.5);
	newU=(2*H0-(J*J)/(newR*newR));
	newU=-sqrt(newU);
	
	sgn=-1;
	H0 = 0.5* (newU*newU + (J*J)/(newR*newR));
	deltaT -= tR0;
	tlim=10*dt;
	tR0=10*dt;
      }
    
    if (tlim<deltaT) {
      newR = sqrt(pow(fabs(newR*newU)-2*sgn*H0*tlim,2)+J*J)*pow(2*H0,-0.5);
      newU = 0;
      
      H0 = 0.5* (J*J)/(newR*newR);
      sgn=-1;
      
      deltaT -= tlim;
      tR0 -= tlim;
      tlim=10*dt;
      
      if (tR0<deltaT) {
	printf("H0 R=%f U=%f (Tr0=%f)",R,U,tR0);exit(0);
	// This will probably never happen
	newR = sqrt(pow(2*H0*tR0,2)+J*J)*pow(2*H0,-0.5);
	newU=(2*H0-(J*J)/(newR*newR));
	newU=sgn*sqrt(newU);
	H0 = 0.5* (newU*newU + (J*J)/(newR*newR));
	sgn=-1;	    
	deltaT-=tR0;
	tR0=10*dt;
      }
    }		
    
    newR = sqrt(pow(fabs(newR*newU)-2*sgn*H0*deltaT,2)+J*J)*pow(2*H0,-0.5);
    newU = (2*H0-(J*J)/(newR*newR));
    if (newU<=0) newU=0;
    else newU = sgn*sqrt(newU);
    
    return std::make_pair(newR,newU);
  }

private:
  void setup()
  { 
    R_C=bT::grid->getCellCoord21_R();
    U_C=bT::grid->getCellCoord21_U();
    J_C=bT::grid->getCellCoord21_J();

    R_V=bT::grid->getVertCoord21_R();
    U_V=bT::grid->getVertCoord21_U();
    J_V=bT::grid->getVertCoord21_J();

    R_Val=bT::grid->getValCoord21_R();
    U_Val=bT::grid->getValCoord21_U();
    J_Val=bT::grid->getValCoord21_J();
    
    const gridParamT &gp=bT::gh->getGridParams();
    const gridParamT &fgp=bT::gh->getFullGridParams();

    Rscale = gp.Pscale[0];
    Uscale = gp.Uscale[0];
    Jscale = gp.Uscale[1];

    Rmin=R_Val[0];

    generateDisplacement(bT::deltaTime,(kernelType==kernelTypeV::REFLECTIVE));
    /*
    std::vector<double> Mr = quadrature21<gridT>::M_of_R(bT::grid);
    //for (long j=1;j<Mr.size();j++) Mr[j]+=Mr[j-1];
    bT::gh->accumulate(Mr,0);
    //for (long j=0;j<Mr.size();j++) printf("(%e,%e) ",R_V[j],Mr[j]);
    if (bT::mpiCom->size() == bT::mpiCom->rank()+1) printf("M = %e\n",Mr.back());
    printf("\n");
    */
    
  }


protected:

  std::vector<double> R_C;
  std::vector<double> U_C;
  std::vector<double> J_C;
  std::vector<double> R_V;
  std::vector<double> U_V;
  std::vector<double> J_V;
  std::vector<double> R_Val;
  std::vector<double> U_Val;
  std::vector<double> J_Val;

  scaleTypeT Rscale;
  scaleTypeT Uscale;
  scaleTypeT Jscale;

  double Rmin;

  virtual void renormalize(double norm)
  {
    if (norm<=0) return;
    std::vector<double> Mr = quadratureT::M_of_less_than_R(bT::gh);
    
    double fac = norm/(Mr.back()+kernelMass);
    
    const long nSlices=bT::grid->nSlices();
    long j;

#pragma omp parallel for
    for (j=0;j<nSlices;j++)
      {
	const gridSliceItT it_end=bT::grid->slice_end(j);
	for (gridSliceItT it=bT::grid->slice_begin(j);it!=it_end;++it)
	  (*it)*=fac;
      }

    if (kernelType==kernelTypeV::DELAYED)
      {
	for (cellChainListItT it=allChains.begin();it!=allChains.end();it++)
	  it->multiply(fac);
      }
  }
  

  virtual void beforeRunning()
  {
    printf("Using %s kernel with M0=%e.\n",kernelTypeSelect().getString(kernelType).c_str(),kernelMass);
  }
  
  virtual void setup(initGenT *init,const paramsParser &params)
  {
    bT::setup(init,params);

    const std::string kernelTypeStr = params.template get<std::string>("kernelType",bT::parserCategory(),"delayed");
    kernelType = kernelTypeSelect().getVal(kernelTypeStr);

    setup();

    if (kernelType!=kernelTypeV::DELAYED) kernelMass=0;
    else kernelMass=init->query(initGenT::queryV::KERNELMASS, Rmin);
    //printf("kernel mass : %e\n",kernelMass);
    setupBoundaryConditions(bT::deltaTime,init);
  }

  virtual void read(FILE *f, bool swap)
  {
    bT::read(f,swap);
    myIO::checkTag(f,getTag());
    
    kernelType = myIO::readEnum<kernelTypeV>(f,swap);
    myIO::fread(&kernelMass,sizeof(double),1,f,swap);
    //printf("kernel mass (%d) : %e\n",(int)kernelType,kernelMass);
    setup();

    // load boundary conditions
    if (kernelType == kernelTypeV::DELAYED)
      {
	long size,size2,i,j,key;
	myIO::fread(&size,sizeof(long),1,f,swap);
    
	kernelBoundary.resize(size);
	allChains.clear();
    
	for (i=0;i<size;i++)
	  {
	    myIO::fread(&size2,sizeof(long),1,f,swap);
	    //printf("size2 = %ld\n",size2);
	    for (j=0;j<size2;j++)
	      {
		myIO::fread(&key,sizeof(long),1,f,swap);
		kernelBoundary[i][key] = allChains.insert(allChains.end(),cellChainT(f,swap));	      
	      }
	  }   
      }
  }

  virtual void write(FILE *f)
  {
    bT::write(f);

    long i;

    myIO::writeTag(f,getTag());
    myIO::writeEnum<kernelTypeV>(f,kernelType);  
    fwrite(&kernelMass,sizeof(double),1,f);

    // save boundary conditions
    if (kernelType == kernelTypeV::DELAYED)
      {
	long size=kernelBoundary.size();
	fwrite(&size,sizeof(long),1,f);
	for (i=0;i<size;i++)
	  {
	    long size2=kernelBoundary[i].size();
	    fwrite(&size2,sizeof(long),1,f);
	    for (kernelCellsItT it = kernelBoundary[i].begin(); it != kernelBoundary[i].end(); it++)
	      {
		fwrite(&it->first,sizeof(long),1,f);
		it->second->write(f);
	      }
	  }
      }
  }

  virtual void snapshot()
  {
    char fname[255];
    sprintf(fname,"snapshot_%4.4f.ND",bT::time);
    NDfield *f=NDfieldFromData();
    Save_NDfield(f,fname);
    free(f);
  }

private:

  typedef typename scaleT::valLocationV valLocationV;

  NDfield *NDfieldFromData(int sliceAtJ=-1)
  {
    long i;
    int dims[P_DIMS+U_DIMS];
    double x0[P_DIMS+U_DIMS];
    double delta[P_DIMS+U_DIMS];
    int ndims;
    long stride[P_DIMS+U_DIMS+1];
    int index;
    dataT *d;
    char comment[80];

    const gridParamT &gp=bT::gh->getGridParams();

    stride[0]=1;
    for (i=0;i<P_DIMS;i++) {
      dims[i]=gp.Pres[i]+1;
      x0[i]=gp.Pmin[i];
      delta[i]=gp.Pmax[i]-gp.Pmin[i];
      stride[i+1]=stride[i]*dims[i];
    }
    long dp=P_DIMS;
    for (i=0;i<U_DIMS;i++) {
      if (i<U_DIMS-1) dims[i+dp]=gp.Ures[i]+1;
      else dims[i+dp]=gp.Ures[i];
      x0[i+dp]=gp.Umin[i];  
      delta[i+dp]=gp.Umax[i]-gp.Umin[i];
      stride[i+dp+1]=stride[i+dp]*dims[i+dp];
    }
    /*
    sprintf(comment,"%s %s %s %s %s %s",
	    scaleT::type2Str(Rscale).c_str(),scaleT::valLocation2Str(gridT::valLocation).c_str(),
	    scaleT::type2Str(Uscale).c_str(),scaleT::valLocation2Str(gridT::valLocation).c_str(),
	    scaleT::type2Str(Jscale).c_str(),scaleT::valLocation2Str(valLocationV::CELL).c_str());
    */
    sprintf(comment,"%d %d %d %d %d %d",
	    (int)Rscale,(int)gridT::valLocation,
	    (int)Uscale,(int)gridT::valLocation,
	    (int)Jscale,(int)valLocationV::CELL);
   
    
    if (sliceAtJ<0) {
      d = bT::grid->getDataPtr();
      ndims=P_DIMS+U_DIMS;
    }
    else {
      ndims=P_DIMS+U_DIMS-1;
      d=bT::grid->fullSlicePtr(sliceAtJ);
    }
    
    return Create_NDfield(dims,ndims,0,ND_DOUBLE,x0,delta,d,comment);
  }

  
  template <int NP, int NU>
  void static_checkCompatibility()
  {
    gridDimensionsAreCompatible<NP,NU>(dimTraits<NP,NU>());
  }

  template <int NP, int NU> void gridDimensionsAreCompatible(dimTraits<1,2>){}
};

#endif
