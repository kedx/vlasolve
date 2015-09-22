#ifndef __VLASOV_SOLVER_21_HXX__
#define __VLASOV_SOLVER_21_HXX__

#include <algorithm>

#include "vlasovSolver.hxx"
#include "dimTraits.hxx"
#include "NDfield.hxx"
#include "typeSelect.hxx"

#include "quadrature21.hxx"
#include "dummyCell21_vect.hxx"
#include "sphericalKernel.hxx"
#include "mpiRequests.hxx"

using namespace NDF;

template <class gridHandlerType>
class vlasovSolver21 : public vlasovSolver<gridHandlerType> {
public:

  static const std::string getTag() {return "VLASOV_SOLVER21 v0.12";}

  typedef vlasovSolver<gridHandlerType> bT;
  typedef vlasovSolver21<gridHandlerType> myT;

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

  typedef typename gridT::sliceIterator gridSliceItT;
  typedef mpiValAtCoordsRequestT<typename gridSliceItT::pointer, mpiComT, DIMS> requestT;

  typedef quadrature21<gridHandlerT> quadratureT;

  vlasovSolver21(mpiComT &mpiCom_) : 
    bT(mpiCom_)
  {
    static_checkCompatibility<P_DIMS,U_DIMS>();
  }

  virtual ~vlasovSolver21()
  {
    
  }  

protected: // Kernel related definitions 
  
  struct kernelTypeV {
    enum type {DELAYED=0, REFLECTIVE=1, UNDEFINED=-1};
  };

  typedef typename kernelTypeV::type kernelTypeT;

  struct kernelTypeSelect : public typeSelect<kernelTypeV> {
    kernelTypeSelect()
    {
      this->insert("delayed",kernelTypeV::DELAYED);
      this->insert("reflective",kernelTypeV::REFLECTIVE);
    }
    std::string name() {return "kernel type";}
  };

  typedef dummyCell21<dataT> dummyCellT;
  typedef sphericalKernel<dummyCellT> sphericalKernelT;
  typedef typename sphericalKernelT::cellChainListItT cellChainListItT;
  typedef typename sphericalKernelT::kernelCellsItT kernelCellsItT;

  kernelTypeT kernelType;
  double kernelMass;
  sphericalKernelT kernel; //delayed kernel

protected:

  std::vector<double> disp_driftR;
  std::vector<double> disp_driftU;  
  std::vector<bool> disp_driftOut; 
  std::vector<int> CFL_R;
  std::vector<int> CFL_U;

  //for communication with surrounding grids
  requestT myRequests;
  requestT neiRequests;
  
  void generateDisplacement(double dt, bool reflectiveKernel=true)
  {
    long j;
    gridT *grid=bT::gh->getGrid();
    const long nSlices=grid->nSlices();
    const double dth=0.5*dt;
    int maxdU=0, maxdR=0;
    double R0=(reflectiveKernel)?kernelRadius:-1;
    
    gridParamT fgp=bT::gh->getFullGridParams();   
    std::vector<double> fgR_Val=scaleT::genScale(fgp.Pmin[0],fgp.Pmax[0],fgp.Pres[0],fgp.Pscale[0],gridT::valLocation_P);
    std::vector<double> fgU_Val=scaleT::genScale(fgp.Umin[0],fgp.Umax[0],fgp.Ures[0],fgp.Uscale[0],gridT::valLocation_U);   
    typedef std::vector<double>::iterator itT;
    
    //FILE *f2=fopen("testUV.i","wb");
    //fprintf(f2,"uv=array(double,[3,4,%ld,%ld]);",nSlices,grid->getNVal()/nSlices);

    //if (!grid->getRegisteredValue("INIGEN_NormFactor",normFactor)) normFactor=1;
    disp_driftR.resize(grid->getNVal());
    disp_driftU.resize(grid->getNVal());
    disp_driftOut.assign(grid->getNVal(),false);
    //CFL_R.resize(grid->getNVal());
    //CFL_U.resize(grid->getNVal());
   
#pragma omp parallel for
    for (j=0;j<nSlices;j++)
      {
	std::pair<double,double> newVal;
	const gridSliceItT it_end=grid->slice_end(j);
	gridSliceItT it=grid->slice_begin(j);
	double J=it.get_J();
	//fprintf(f2,"id=1;\n");
	for (;it!=it_end;++it)
	  {
	    long i=it.get_i();	    
	    newVal=advanceParticleDrift(it.get_R(),it.get_U(),J,dth,R0);
	    disp_driftR[i]=newVal.first;
	    disp_driftU[i]=newVal.second;
	    if ((disp_driftR[i]>=R_Val.back())||
		(disp_driftR[i]<=R_Val.front())||
		(disp_driftU[i]>=U_Val.back())||
		(disp_driftU[i]<=U_Val.front())) disp_driftOut[i]=true;
	    // to check for CFL condition 
	    
	    if ((disp_driftR[i]>=fgR_Val.back())||
		(disp_driftR[i]<=fgR_Val.front())||
		(disp_driftU[i]>=fgU_Val.back())||
		(disp_driftU[i]<=fgU_Val.front())) {/*CFL_R[i]=CFL_U[i]=0;*/}
	    else
	      {
	        itT itR=std::upper_bound(fgR_Val.begin(),fgR_Val.end(),disp_driftR[i]);	
    
		if (itR!=fgR_Val.begin()) itR--;
		int deltaR = (int)(itR-fgR_Val.begin())-it.get_r();
		if (deltaR<0) deltaR=-deltaR; else deltaR++;
		
		itT itU;	
		if (disp_driftU[i] * it.get_U() >= 0)
		  itU=std::upper_bound(fgU_Val.begin(),fgU_Val.end(),disp_driftU[i]);	
		else
		  itU=std::upper_bound(fgU_Val.begin(),fgU_Val.end(),-disp_driftU[i]);	

		if (itU!=fgU_Val.begin()) itU--;
		int deltaU = (int)(itU-fgU_Val.begin())-it.get_u();
		if (deltaU<0) deltaU=-deltaU; else deltaU++;
	    
		//CFL_R[i]=deltaR;
		//CFL_U[i]=deltaU;

		/*
		if ((CFL_R[i]>maxdR)||(CFL_U[i]>maxdU))
		  {
		    printf("delta[j(%ld)=%g] = [dR<%d,DU<%d]\n",it.get_j(),it.get_J(),deltaR,deltaU);
		    printf("[%g,%g] -> [%g,%g]\n",it.get_R(),it.get_U(),disp_driftR[i],disp_driftU[i]);
		    if (deltaR>maxdR) maxdR=deltaR;
		    if (deltaU>maxdU) maxdU=deltaU;	    
		  }
		*/
	      }
	    
	  }
      }

    //printf("Maximum displacement for CFL condition: [dR<%d,DU<%d]\n",maxdR,maxdU);
    //fclose(f2);
    /*
    FILE *f=fopen("testDisp.dat","wb");
    for (j=0;j<nSlices;j++)
      {
	std::pair<double,double> newVal;
	const gridSliceItT it_end=grid->slice_end(j);
	gridSliceItT it=grid->slice_begin(j);
	double J=it.get_J();
	//fprintf(f2,"id=1;\n");
	for (;it!=it_end;++it)
	  {
	    long i=it.get_i();	    
	    fprintf(f,"%ld %g %g %g %g\n",j,it.get_R(),it.get_U(),disp_driftR[i],disp_driftU[i]);
	  }
      }
    fclose(f);
    */
  }

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
	for (;it!=it_end;++it)
	  {
	    long i=it.get_i();
	    dirT dir=0;
	    
	    if (disp_driftR[i]>R_Val.back()) 
	      dir+=gridNav::dir(0,1);
	    else if (disp_driftR[i]<R_Val.front())
	      dir+=gridNav::dir(0,-1);
	    
	    if (disp_driftU[i]>U_Val.back()) 
	      dir+=gridNav::dir(1,1);
	    else if (disp_driftU[i]<U_Val.front()) 
	      dir+=gridNav::dir(1,-1);
	    
	    int neiID=bT::gh->getNeighbor(dir);
	    if (neiID<0) continue;
	    else 
	      {
		myRequests.procID.push_back(neiID);
		myRequests.ptr.push_back(&(*it));
		myRequests.coords.push_back(disp_driftR[i]);
		myRequests.coords.push_back(disp_driftU[i]);
		myRequests.coords.push_back((double)j);
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

  void setupBoundaryConditions(double dt, initGenT *init)
  {
    if (kernelType!=kernelTypeV::DELAYED) return;

    gridT *grid=bT::gh->getGrid();
    const long nSlices=grid->nSlices();  
    bool isLowB=bT::gh->isFullGridBoundary_Plow(0); 
    long j;
    double dth=0.5*dt;    
    std::vector<double> &driftR=disp_driftR;
    std::vector<double> &driftU=disp_driftU;

    kernel.reset();
    kernel.reserve(nSlices);    

    // !! no openMP cause cellchains are not thread safe ... !!
    for (j=0;j<nSlices;j++)
      {
	std::pair<double,double> newVal;

	const gridSliceItT sit_end=grid->slice_end(j);
	gridSliceItT sit=grid->slice_begin(j);
	double J=sit.get_J();
	for (;sit!=sit_end;++sit)
	  {
	    const long i=sit.get_i();

	    if (driftR[i]<kernelRadius)
	      {
		if (!isLowB)
		  {
		    fprintf(stderr,"ERROR: time step is too high for this grid.\n");
		    fprintf(stderr,"    R=%g => R_drift=%g<%g (%g).\n",sit.get_R(),driftR[i],R_Val[0],kernelRadius);
		    exit(-1);
		  }

		double deltaT=dth;
		cellChainListItT cit = kernel.newChain(j,i);		    
		std::pair<double,double> newVal=advanceParticleDrift(sit.get_R(),sit.get_U(),J,deltaT);
		std::vector<double> pos(3);
		pos[2]=J;
		//int parity=0;
		while (newVal.first<kernelRadius)
		  {
		    pos[0]=newVal.first;pos[1]=newVal.second;
		    dataT val[bT::sp.nSpecies];
		    for (int s=0;s<bT::sp.nSpecies;s++) 
		      val[s]=init->valueAt(pos,s);
		    //dataT val=init->valueAt(pos);
		    dummyCellT cell(newVal.first,newVal.second,&val[0],&val[0]+bT::sp.nSpecies);
		    cit->push(cell);
		    deltaT+=dth;
		    newVal=advanceParticleDrift(sit.get_R(),sit.get_U(),J,deltaT);
		    //parity++;
		  }	
		cit->push(dummyCellT(newVal.first,newVal.second,1));
	      }
	  }
      }
   
    //printf("Nchains : %ld\n",allChains.size());
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
    gridT *grid=bT::gh->getGrid();

    R_C=grid->getCellCoord21_R();
    U_C=grid->getCellCoord21_U();
    J_C=grid->getCellCoord21_J();

    R_V=grid->getVertCoord21_R();
    U_V=grid->getVertCoord21_U();
    J_V=grid->getVertCoord21_J();

    R_Val=grid->getValCoord21_R();
    U_Val=grid->getValCoord21_U();
    J_Val=grid->getValCoord21_J();
    
    const gridParamT &gp=bT::gh->getGrid()->getParams();
    const gridParamT &fgp=bT::gh->getFullGridParams();

    Rscale = gp.Pscale[0];
    Uscale = gp.Uscale[0];
    Jscale = gp.Uscale[1];
    
    bool refKernel=(kernelType==kernelTypeV::REFLECTIVE);
    
    if (bT::gh->isFullGridBoundary_Plow(0)) 
    	kernelRadius = R_Val[kernelBufferSize];
    else kernelRadius=0;

    kernelRadius=bT::mpiCom->max(kernelRadius);

    printf("kernel radius: %e\n",kernelRadius);
    if (bT::noLinking) generateDisplacement(bT::deltaTime,refKernel);
    else generateDisplacement(bT::deltaTime*2,refKernel);
    setupRequests();
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

  double kernelRadius;
  int kernelBufferSize;

  virtual void renormalize(double norm, int which)
  {
    if (norm<=0) return;
    gridT *grid=bT::gh->getGrid();
    double mass;
    
    if (which<0)
      mass=quadratureT::integrate(bT::gh,int21::f(0,bT::sp.nSpecies),true)+kernelMass;  
    else
      mass=quadratureT::integrate(bT::gh,int21::f(which,1),true)+kernelMass;  
    
    double fac=norm/mass;
    //double fac=norm/(mass+kernelMass);
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
	    //printf("j=%ld\n",j);MPI_Barrier(MPI_COMM_WORLD);fflush(0);MPI_Barrier(MPI_COMM_WORLD);
	    const field_iterator<gridSliceItT> it_end(grid->slice_end(j));
	    field_iterator<gridSliceItT> it(grid->slice_begin(j));
	    for (;it!=it_end;++it)
	      (*it)*=fac;
	  }

	//MPI_Barrier(MPI_COMM_WORLD);exit(0);
	if (kernelType==kernelTypeV::DELAYED) 
	  for (j=0;j<bT::sp.nSpecies;j++) kernel.multiply(fac,j);
      }
    else
      {
	#pragma omp parallel for
	for (j=0;j<nSlices;j++)
	  {
	    //printf("j=%ld\n",j);MPI_Barrier(MPI_COMM_WORLD);fflush(0);MPI_Barrier(MPI_COMM_WORLD);
	    const oneField_iterator<gridSliceItT> it_end(grid->slice_end(j),which);
	    oneField_iterator<gridSliceItT> it(grid->slice_begin(j),which);
	    for (;it!=it_end;++it)
	      (*it)*=fac;
	  }
	
	//MPI_Barrier(MPI_COMM_WORLD);exit(0);
	if (kernelType==kernelTypeV::DELAYED) 
	  kernel.multiply(fac,which);
      }
  }
  

  virtual void beforeRunning()
  {
    bT::beforeRunning();
    if (bT::gh->isFullGridBoundary_Plow(0)) printf("Using %s kernel with M0=%e.\n",kernelTypeSelect().getString(kernelType).c_str(),kernelMass);
  }
  
  virtual void setup(initGenT *init,const paramsParser &params)
  {
    bT::setup(init,params);

    const std::string kernelTypeStr = params.template get<std::string>("kernelType",bT::parserCategory(),"delayed");
    kernelType = kernelTypeSelect().getVal(kernelTypeStr);
    kernelBufferSize = params.template get<>("kernelBufferSize",bT::parserCategory(),0);
    setup();

    if (kernelType!=kernelTypeV::DELAYED) kernelMass=0;
    else kernelMass=0;//init->query(initGenT::queryV::KERNELMASS, kernelRadius);
    printf("kernel mass : %e\n",kernelMass);
    
    if (bT::noLinking) setupBoundaryConditions(bT::deltaTime,init);
    else setupBoundaryConditions(bT::deltaTime*2,init);
  }

  virtual void read(FILE *f, bool swap)
  {
    bT::read(f,swap);
    myIO::checkTag(f,getTag());
    
    kernelType = myIO::readEnum<kernelTypeV>(f,swap);
    myIO::fread(&kernelMass,sizeof(double),1,f,swap);
    myIO::fread(&kernelBufferSize,sizeof(int),1,f,swap);
    //printf("kernel mass (%d) : %e\n",(int)kernelType,kernelMass);

    std::vector<char> spare(256*8,0);
    myIO::fread(&spare[0],sizeof(char),spare.size(),f,swap);
 
    setup();

    if (kernelType == kernelTypeV::DELAYED) kernel.read(f,swap);    
  }

  virtual void write(FILE *f)
  {
    bT::write(f);

    long i;

    myIO::writeTag(f,getTag());
    myIO::writeEnum<kernelTypeV>(f,kernelType);  
    fwrite(&kernelMass,sizeof(double),1,f);
    fwrite(&kernelBufferSize,sizeof(int),1,f);

    std::vector<char> spare(256*8,0);
    fwrite(&spare[0],sizeof(char),spare.size(),f);

    if (kernelType == kernelTypeV::DELAYED) kernel.write(f);
  }

  virtual void snapshot()
  {
    char fname[1024];
    sprintf(fname,"snapshot_%4.4f",bT::time);
 
    if (bT::mpiCom->size()>1)
      sprintf(fname,"%s_%6.6d.ND",fname,bT::mpiCom->rank());
    else
      sprintf(fname,"%s.ND",fname);
    
    NDfield *f=NDfieldFromData();
    Save_NDfield(f,fname);
    free(f);
  }
  /*
  template <class iterator>
  static double kineticEnergy(double r, double u, double j, const iterator &it)
  {
    return 0.5*it.get_avg()*(u*u+j*j/(r*r));
  }
  
  template <class iterator>
  static double entropy(double r, double u, double j, const iterator &it)
  {
    double val=it.get_avg();
    return (val<=0)?0:-val*log(val);
  }

  template <class iterator>
  static double L1Norm(double r, double u, double j, const iterator &it)
  {
    return fabs(it.get_avg());
  }

  template <class iterator>
  static double L2Norm(double r, double u, double j, const iterator &it)
  {
    double val=it.get_avg();
    return val*val;
  }
  */

  double potentialEnergy(double &Ep)
  {
    //std::vector<double> Mr = quadratureT::M_of_R(bT::gh,true);
    std::vector<double> Mr = quadratureT::integrate1D(bT::gh,0,int21::f(0,bT::sp.nSpecies),true);
    const std::vector<double> &R=bT::gh->getFullGridCoords(0);
    double M=Mr[0]+kernelMass;
    Ep=0;

    for (int i=1;i<Mr.size();i++) 
      {
	double r=0.5*(R[i-1]+R[i]);	
	Ep-=(Mr[i]*bT::sp.G*M)/r;
	M+=Mr[i];
      }
    return M;
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
     
	double Ep;
	double M=potentialEnergy(Ep);
	double Ek=quadratureT::integrate(bT::gh,int21::kineticEnergy(),true);
	double E=Ek+Ep;
	double S=quadratureT::integrate(bT::gh,int21::entropy(),true);  
	double L1=quadratureT::integrate(bT::gh,int21::L1Norm(),true);  
	double L2=quadratureT::integrate(bT::gh,int21::L2Norm(),true);

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

	fprintf(f,"%.10g %ld %10.10g %10.10g %10.10g %10.10g %10.10g %10.10g %10.10g\n",bT::time,bT::stepCount,M,E,Ek,Ep,S,L1,L2);
    
	fclose(f);
      }

    if (nsp>1)
      {
	sprintf(fname,"statistics.txt");
	double S_tot=quadratureT::integrate(bT::gh,int21::entropy(0,nsp),true);
	double L2_tot=quadratureT::integrate(bT::gh,int21::L2Norm(0,nsp),true);
      
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

  typedef typename scaleT::valLocationV valLocationV;

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
    for (i=0;i<P_DIMS;i++) {
      dims[i+dp]=gp.Pres[i]+1;
      x0[i+dp]=gp.Pmin[i];
      delta[i+dp]=gp.Pmax[i]-gp.Pmin[i];
      stride[i+1+dp]=stride[i+dp]*dims[i+dp];
    }
    dp+=P_DIMS;
    for (i=0;i<U_DIMS;i++) {
      if (i<U_DIMS-1) dims[i+dp]=gp.Ures[i]+1;
      else dims[i+dp]=gp.Ures[i];
      x0[i+dp]=gp.Umin[i];  
      delta[i+dp]=gp.Umax[i]-gp.Umin[i];
      stride[i+dp+1]=stride[i+dp]*dims[i+dp];
    }

    sprintf(comment,"%d %d %d %d %d %d",
	    (int)Rscale,(int)gridT::valLocation_P,
	    (int)Uscale,(int)gridT::valLocation_U,
	    (int)Jscale,(int)gridT::valLocation_J);
   
    
    if (sliceAtJ<0) {
      d = grid->getDataPtr();
      ndims=P_DIMS+U_DIMS+useSpecies;
    }
    else {
      ndims=P_DIMS+U_DIMS-1+useSpecies;
      d=grid->fullSlicePtr(sliceAtJ);
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
