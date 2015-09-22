#ifndef __VLASOV_SOLVERNN_SPLINE_HXX__
#define __VLASOV_SOLVERNN_SPLINE_HXX__

#ifdef HAVE_FFTW3

#include <algorithm>
#include "vlasovSolverNN.hxx"
#include "quadratureNN.hxx"
#include "interpol_spline.hxx"
#include "dimTraits.hxx"

#include "localSpline.hxx"
#include "openMP_interface.hxx"

template <class gridHandlerType>
class vlasovSolverNN_spline : public vlasovSolverNN<gridHandlerType> {
public:
   
  typedef vlasovSolverNN<gridHandlerType> bNNT;
  typedef typename bNNT::bT bT;
  typedef vlasovSolverNN_spline<gridHandlerType> myT;
  
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
  typedef typename bT::gridItT gridItT;
  typedef typename bT::dirT dirT;
  typedef typename bNNT::gridSliceItT gridSliceItT;
  typedef typename bNNT::requestT requestT;
  //typedef typename bNNT::projectedPSlicesT densityGridT;
  
  typedef typename bNNT::quadratureT quadratureT;

  static const std::string getTag() {char tmp[255];sprintf(tmp,"%s%d%s","VLASOV_SOLVER_",PU_DIMS,"_SPLINE v0.10");return std::string(tmp);}

  vlasovSolverNN_spline(mpiComT &mpiCom_) : 
    bNNT(mpiCom_)
  {
    lbCount[0]=lbCount[1]=0;

    for (int i=0;i<2;i++)
      for (int j=0;j<P_DIMS;j++)
	for (int k=0;k<2;k++)
	  localBoundarySize[i][j][k]=0;
  }

  virtual ~vlasovSolverNN_spline()
  {
    
  }  
  
private:

  // lookup tables for requests an drift
  std::vector<int> driftR_I;
  std::vector<double> driftR_F; 
  std::vector<int> neiReqCoord_I; //interger coordinate for requests

  void decomposeDisplacement()
  {
    long i,j,k;
    interpolT &spl=spline[0][0][0];
    gridT *grid=bT::gh->getGrid();
    long nVal=grid->getNVal();
    
    // This uses a lot of memory, should have an option for skipping it ?
    //printf("Precomputing lookup tables ... ");fflush(0);
    driftR_I.resize(nVal*P_DIMS);
    driftR_F.resize(nVal*P_DIMS);
	
    const long nSlices=grid->nSlices();
    
#pragma omp parallel for 
    for (j=0;j<nSlices;j++)
      {
	const gridSliceItT it_end=grid->slice_end(j);
	gridSliceItT it=grid->slice_begin(j);
	double disp_drift[DIMS];
	getDrift(it,disp_drift);

	for (;it!=it_end;++it)
	  {
	    i=it.get_i()*P_DIMS;
	    updateDrift(it,disp_drift);
	    for (k=0;k<P_DIMS;k++)
	      {
		std::pair<int,double> r=spl.pos2Coord(disp_drift[k],k);
		driftR_I[i+k]=r.first;
		driftR_F[i+k]=r.second;
	      }
	  }
      }

    requestT &req = bNNT::neiRequests;
    if (!req.nProcs) {return;}
    // NOTE : SUPPOSE NO MPI MARGIN IN V !!
    neiReqCoord_I.clear();
    neiReqCoord_I.reserve((req.coords.size()*(P_DIMS))/DIMS);

    for (i=0;i<req.coords.size();i+=DIMS)
      {
	//int w[DIMS];
	//grid->pos2Coord(&req.coords[i],w);
	
	for (j=0;j<P_DIMS;j++)
	  {	    
	    std::pair<int,double> r=spl.pos2Coord(req.coords[i+j],j);
	    neiReqCoord_I.push_back(r.first);	    
	    req.coords[i+j]=r.second;
	  }
	//printf("%d == %d\n",(int)req.coords[i+DIMS-1],(int)grid->sliceIndex(w));
	  //neiReqCoord_I.push_back(grid->sliceIndex(w));
      }
    //printf("done.\n");
    
  }

  typedef localSpline<P_DIMS,dataT> interpolT;
  typedef typename interpolT::initT interpInitT;
  typedef typename interpolT::boundaryTypeT boundaryTypeT;
  typedef typename interpolT::boundaryTypeV boundaryTypeV; 
  typedef typename interpolT::localBoundaryT localBoundaryT;
  std::vector< std::vector<interpolT> > spline[2];  

  int lbCount[2];
  long localBoundarySize[2][P_DIMS][2];
  localBoundaryT localBoundaryBufferI;
  localBoundaryT localBoundaryBufferO;
  std::vector<localBoundaryT>  localBoundaryPerThread;

  void setupSplines(int type)
  {
    std::vector< std::vector<interpolT> > &spline = myT::spline[type];
    gridT *grid=bT::gh->getGrid();
    const int nSpecies=bT::sp.nSpecies;     
    const long nSlices=grid->nSlices(type,true);  
    long j,s;  
    
    //printf("Setting up splines[%d] : %ld slices.\n",type,nSlices);
    
    for (s=0;s<nSpecies;s++)
      {
	if (!lbCount[type])
	  {
	    #pragma omp parallel for 
	    for (j=0;j<nSlices;j++)
	      {
		const oneField_iterator<gridSliceItT> start(grid->slice_begin(j,type,true),s);
		const oneField_iterator<gridSliceItT> stop(grid->slice_end(j,type,true),s);
		spline[s][j].assign(start,stop);
		spline[s][j].build();
	      }
	    continue;
	  }
	
#pragma omp parallel for 
	for (j=0;j<nSlices;j++)
	  {
	    interpolT &spl=spline[s][j];
	    localBoundaryT &lb=localBoundaryPerThread[omp_get_thread_num()];
	    const oneField_iterator<gridSliceItT> start(grid->slice_begin(j,type,true),s);
	    const oneField_iterator<gridSliceItT> stop(grid->slice_end(j,type,true),s);
	    spl.assign(start,stop);
	    spl.getLocalBoundaries(lb);
	    //if (j==0) printf("(%d:) lb = [[%ld,%ld][%ld,%ld]]\n",bT::mpiCom->rank(),
	    //	   lb.B[0][0].size(),lb.B[0][1].size(),lb.B[1][0].size(),lb.B[1][1].size());
	    for (int i=0;i<P_DIMS;i++)
	      {
		//printf("N=%ld*%ld : %ld %ld\n",spline[j].boundarySize(i),j,lb.B[i][0].size(),lb.B[i][1].size());
		if (lb.B[i][0].size()) 
		  std::copy(lb.B[i][0].begin(),lb.B[i][0].end(),&localBoundaryBufferO.B[i][0][spl.boundarySize(i)*j]);
		if (lb.B[i][1].size())
		  std::copy(lb.B[i][1].begin(),lb.B[i][1].end(),&localBoundaryBufferO.B[i][1][spl.boundarySize(i)*j]);
	      }
	  }
	
	int Nreq=0;
	MPI_Request request[4*2];
      
	for (int i=0;i<P_DIMS;i++)
	  {
	    if (localBoundarySize[type][i][0]) 
	      {
		int procID=bT::gh->getNeighbor(i,-1);
		//printf("(%d:) [%d %d %d]: send:%ld==%ld rec:%ld==%ld\n",bT::mpiCom->rank(),type,i,0,
		//       localBoundaryBufferO.B[i][0].size(),localBoundarySize[type][i][0],
		//       localBoundaryBufferI.B[i][0].size(),localBoundarySize[type][i][0]);
		// MPI_Isend(&localBoundaryBufferO.B[i][0][0],localBoundarySize[type][i][0],
		// 	  MPI_DOUBLE,procID,0,MPI_COMM_WORLD,&request[Nreq++]);
		// MPI_Irecv(&localBoundaryBufferI.B[i][0][0],localBoundarySize[type][i][0],
		// 	  MPI_DOUBLE,procID,0,MPI_COMM_WORLD,&request[Nreq++]);
		bT::mpiCom->Isend(&localBoundaryBufferO.B[i][0][0],localBoundaryBufferO.B[i][0].size(),
				  procID,&request[Nreq++],0);
		bT::mpiCom->Irecv(&localBoundaryBufferI.B[i][0][0],localBoundaryBufferI.B[i][0].size(),
				     procID,&request[Nreq++],0);
		
	      }
	    if (localBoundarySize[type][i][1])
	      {
		int procID=bT::gh->getNeighbor(i,+1);
		//printf("(%d:) [%d %d %d]: send:%ld==%ld rec:%ld==%ld\n",bT::mpiCom->rank(),type,i,1,
		//       localBoundaryBufferO.B[i][1].size(),localBoundarySize[type][i][1],
		//       localBoundaryBufferI.B[i][1].size(),localBoundarySize[type][i][1]);

		// MPI_Isend(&localBoundaryBufferO.B[i][1][0],localBoundarySize[type][i][1],
		// 	  MPI_DOUBLE,procID,0,MPI_COMM_WORLD,&request[Nreq++]);
		// MPI_Irecv(&localBoundaryBufferI.B[i][1][0],localBoundarySize[type][i][1],
		// 	  MPI_DOUBLE,procID,0,MPI_COMM_WORLD,&request[Nreq++]);
		
		bT::mpiCom->Isend(&localBoundaryBufferO.B[i][1][0],localBoundaryBufferO.B[i][1].size(),
				  procID,&request[Nreq++],0);
		bT::mpiCom->Irecv(&localBoundaryBufferI.B[i][1][0],localBoundaryBufferI.B[i][1].size(),
				     procID,&request[Nreq++],0);
	      }
	  }
	// MPI_Status status;
	// for (int i=0;i<Nreq;i++) MPI_Wait(&request[i], &status);
	bT::mpiCom->Waitall(Nreq,&request[0]);

	bT::mpiCom->barrier();
#pragma omp parallel for 
	for (j=0;j<nSlices;j++)
	  {
	    interpolT &spl=spline[s][j];
	    localBoundaryT &lb=localBoundaryPerThread[omp_get_thread_num()];
	
	    for (int i=0;i<P_DIMS;i++)
	      {
	    
		if (localBoundarySize[type][i][0])
		  {
		    dataT *it = &localBoundaryBufferI.B[i][0][spl.boundarySize(i)*j];
		    std::copy(it,it+spl.boundarySize(i),lb.B[i][0].begin());
		  }
		if (localBoundarySize[type][i][1])
		  {
		    dataT *it = &localBoundaryBufferI.B[i][1][spl.boundarySize(i)*j];
		    std::copy(it,it+spl.boundarySize(i),lb.B[i][1].begin());
		  }
	      }

	    spl.setLocalBoundaries(lb);
	    spl.build();
	  }	

      }     
      
  }

  void updateRequests(int type=0)
  {
    if (type!=0) return;

    gridT *grid=bT::gh->getGrid();
    long nfields=grid->getNFields();
    requestT &neiReq=bNNT::neiRequests;
    requestT &myReq=bNNT::myRequests;  
    const long N=neiReq.val.size()/nfields;    
    const int nSpecies=bT::sp.nSpecies;
    
    long i,ii,iii,s;
    std::vector<double> &val=neiReq.val;
    std::vector<double> &coords=neiReq.coords;
    //std::vector<double> &myVal=myReq.val;

    std::vector< std::vector<interpolT> > &spline = myT::spline[type];
    
    const long nSlices=grid->nSlices(type,0);

    for (s=0;s<nfields;s++)
      {
#pragma omp parallel for
	for (i=0;i<N;i++)
	  {
	    //int j=grid->sliceIndex(&neiReqCoord_I[P_DIMS*i]);if (j>=nSlices) printf("j=%d > %ld\n",j,nSlices);
	    int *reqC=&neiReqCoord_I[(P_DIMS)*i];
	    //int j=reqC[P_DIMS];
	    int j=(int)coords[DIMS*(i+1)-1];//printf("j=%d > %ld == %ld\n",j,nSlices,spline[s].size());
	    val[i*nfields+s]=spline[s][j].evaluateC(reqC,&coords[DIMS*i]);
	  }
      }
    //printf("HELLO %d\n",bT::mpiCom->rank());bT::mpiCom->barrier();exit(0);
    myReq.synchronize(neiReq);
  }
  
  void setup()
  {
    const int nSpecies=bT::sp.nSpecies;
    interpInitT ini;
    gridT *grid=bT::gh->getGrid();
    //long nval;
    long i;
    printf("Initializing splines ... ");fflush(0);
    for (int type=0;type<2;type++)
      {
	const long nSlices=grid->nSlices(type);
	std::vector< std::vector<interpolT> > &spl = myT::spline[type];
	//nval=1;
	lbCount[type]=0;
	if (type==0)
	  {
	    for (i=0;i<P_DIMS;i++)
	      {
		ini.start[i]=bNNT::P_V[i].front();
		ini.stop[i]=bNNT::P_V[i].back();
		ini.N[i]=bNNT::P_C[i].size();
		ini.scaleType[i]=bNNT::Pscale[i];
		ini.vlType[i]=grid->valLocation_P;

		if (bT::gh->isFullGridBoundary_Plow(i))
		  ini.bCond[i][0] = boundaryTypeV::BT_NATURAL;//boundaryTypeV::BT_NOT_A_KNOT;//boundaryTypeV::BT_NATURAL;
		else
		  {
		    //printf("BT:: Low local!!\n");
		    ini.bCond[i][0] = boundaryTypeV::BT_LOCAL;
		    lbCount[type]++;
		  }

		if (bT::gh->isFullGridBoundary_Phigh(i))
		  ini.bCond[i][1] = boundaryTypeV::BT_NATURAL;//boundaryTypeV::BT_NOT_A_KNOT;//boundaryTypeV::BT_NATURAL;
		else
		  {	
		    //printf("BT:: high local!!\n");
		    ini.bCond[i][1] = boundaryTypeV::BT_LOCAL;
		    lbCount[type]++;
		  }
	      }
	    //for (i=0;i<P_DIMS;i++) nval*=bNNT::U_Val[i].size();
	  }
	else
	  {
	    for (i=0;i<U_DIMS;i++)
	      {
		ini.start[i]=bNNT::U_V[i].front();
		ini.stop[i]=bNNT::U_V[i].back();
		ini.N[i]=bNNT::U_C[i].size();
		ini.scaleType[i]=bNNT::Uscale[i];
		ini.vlType[i]=grid->valLocation_U;
	      
		if (bT::gh->isFullGridBoundary_Ulow(i))
		  ini.bCond[i][0] = boundaryTypeV::BT_NATURAL;
		else
		  {
		    ini.bCond[i][0] = boundaryTypeV::BT_LOCAL;
		    lbCount[type]++;
		  }

		if (bT::gh->isFullGridBoundary_Uhigh(i))
		  ini.bCond[i][1] = boundaryTypeV::BT_NATURAL;
		else
		  {
		    ini.bCond[i][1] = boundaryTypeV::BT_LOCAL;
		    lbCount[type]++;
		  }
		//for (i=0;i<U_DIMS;i++) nval*=bNNT::P_Val[i].size();
	      }
	  }

	spl.resize(nSpecies);
	for (long s=0;s<nSpecies;s++)
	  spl[s].assign(nSlices,interpolT(ini));
	//printf("lbcount=%d, nv=%ld\n",lbCount[type],nSlices);
	if (lbCount[type]) 
	  {	    
	    localBoundaryPerThread.resize(num_omp_threads);	    

	    for (int i=0;i<P_DIMS;i++)
	      {
		long sz=spl[0][0].boundarySize(i)*nSlices;
		
		if (ini.bCond[i][0]==boundaryTypeV::BT_LOCAL) 
		  {
		    //printf("bsize(%d) : %d %d %ld -> %ld %ld\n",type, i,0,spl[0][0].boundarySize(i),nSlices,sz);
		    localBoundaryBufferI.B[i][0].resize(std::max(sz,(long)localBoundaryBufferI.B[i][0].size()));
		    localBoundaryBufferO.B[i][0].resize(std::max(sz,(long)localBoundaryBufferO.B[i][0].size()));
		    localBoundarySize[type][i][0]=sz;
		  } else localBoundarySize[type][i][0]=0;

		if (ini.bCond[i][1]==boundaryTypeV::BT_LOCAL)
		  {
		    //printf("bsize(%d) : %d %d %ld -> %ld %ld\n",type,i,1,spl[0][0].boundarySize(i),nSlices,sz);
		    localBoundaryBufferI.B[i][1].resize(std::max(sz,(long)localBoundaryBufferI.B[i][1].size()));
		    localBoundaryBufferO.B[i][1].resize(std::max(sz,(long)localBoundaryBufferO.B[i][1].size()));
		    localBoundarySize[type][i][1]=sz;
		  } else localBoundarySize[type][i][1]=0;
		
	      }
	  }
      }
    
    decomposeDisplacement();
    printf("done.\n");
  }

protected:

  virtual int nSolverFields() {return 1;}

  virtual void setup(initGenT *init,const paramsParser &params)
  {
    bNNT::setup(init,params);
    setup();
  }

  virtual void read(FILE *f, bool swap)
  {
    bNNT::read(f,swap);
    myIO::checkTag(f,getTag());

    std::vector<char> spare(256*8,0);
    myIO::fread(&spare[0],sizeof(char),spare.size(),f,swap);
    setup();
  }

  virtual void write(FILE *f)
  {
    bNNT::write(f);
    myIO::writeTag(f,getTag());

    std::vector<char> spare(256*8,0);
    fwrite(&spare[0],sizeof(char),spare.size(),f);
  }

   void makeOneStep(double t)
  {   
    gridT *grid=bT::gh->getGrid();
    int nFieldsTotal=grid->getNFields();
    long j,s;
    double dth=0.5*bT::deltaTime;
    long delta = grid->getCellStride(2);
    //printf("NOT INMPLEMENT makeonestep\n");exit(0);
    const long nSlices=grid->nSlices(); 
    const long nSlices1=grid->nSlices(1);
    const int nSpecies=bT::sp.nSpecies;

    if (bT::noLinking)
      {
	bT::gh->synchronize();
	setupSplines(0);
	
	for (s=0;s<nSpecies;s++)
	  {
#pragma omp parallel for
	    for (j=0;j<nSlices;j++)
	      { 
		
		interpolT &spl=spline[0][s][j];
		double drift[P_DIMS];
		oneField_iterator<gridSliceItT> it(grid->slice_begin(j),s);
		const oneField_iterator<gridSliceItT> it_end(grid->slice_end(j),s);
		
		getDrift(it,drift);

		for (;it!=it_end;++it)
		  {
		    if (updateDrift(it,drift)) (*it=0);
		    else (*it)=spl.evaluate(drift);
		  }
		
	      }
	  }
	updateRequests();
      }

    quadratureT::project(bT::gh,intNN::f<P_DIMS>(0,nSpecies),bNNT::density);  
    //bNNT::density.writeToNDfield(std::string("density"));
    bNNT::poissonSolver.solve();
    //bNNT::force.interleave(bNNT::forcePoisson,bNNT::forcePoisson+P_DIMS,bT::deltaTime);    
    //bNNT::force.writeToNDfield(std::string("force"));exit(0);
    bT::gh->synchronize(); //bNNT::snapshot();exit(0);
    setupSplines(1);
    
    for (s=0;s<nSpecies;s++)
      {	
#pragma omp parallel for
	for (j=0;j<nSlices1;j++)
	  { 
	    interpolT &spl=spline[1][s][j];
	    oneField_iterator<gridSliceItT> it(grid->slice_begin(j,1),s);
	    const oneField_iterator<gridSliceItT> it_end(grid->slice_end(j,1),s);
	    int i;
	    double newU[P_DIMS];
	    
	    double dsp[P_DIMS];
	    for (i=0;i<P_DIMS;i++) 
	      dsp[i]=bT::deltaTime*(*bNNT::forcePoisson[i].getDataPtr(it.get_w()));
	    
	    //dataT *disp=bNNT::force.getDataPtr(it.get_w());
	    
	    for (;it!=it_end;++it)
	      {	
		it.CU(newU);
		//printf("[%g,%g] -> ",newU[0],newU[1]);
		for (i=0;i<P_DIMS;i++) 
		  {
		    newU[i]-=dsp[i];
		    if ((newU[i]<bNNT::Umin[i])||(newU[i]>bNNT::Umax[i]))
		      {
			(*it=0);
			break;
		      }
		  }
		if (i==P_DIMS) (*it)=spl.evaluate(newU);
		//printf("[%g,%g] ([%g,%g]<..<[%g,%g])\n",newU[0],newU[1],
		//bNNT::Umin[0],bNNT::Umin[1],
		//bNNT::Umax[0],bNNT::Umax[1]);
	      }
	    
	  }
      }    
    //exit(0);
    
    bT::gh->synchronize();  
    setupSplines(0);

    for (s=0;s<nSpecies;s++)
      {
#pragma omp parallel for
	for (j=0;j<nSlices;j++)
	  { 
	    interpolT &spl=spline[0][s][j];
	    double drift[DIMS];
	    //double coords[DIMS];
	    oneField_iterator<gridSliceItT> it(grid->slice_begin(j),s);
	    const oneField_iterator<gridSliceItT> it_end(grid->slice_end(j),s);
	    
	    getDrift(it,drift);
	    
	    for (;it!=it_end;++it)
	      {	
		//it.get_C(coords);
		//fprintf(f1,"drift(,i)=[%lg,%lg,",coords[0],coords[1]);
		//printf("[%lg %lg][%lg %lg]",coords[0],coords[1],coords[2],coords[3]);

		if (updateDrift(it,drift)) (*it=0);
		//if (getDrift(it,drift)) (*it=0);
		else (*it)=spl.evaluate(drift);
		
		//fprintf(f1,"%lg,%lg,%lg,%lg];i++;\n",drift[0],drift[1],coords[2],coords[3]);
		  //printf(" ==> [%lg %lg]\n",drift[0],drift[1]);
	      }	
	  }
      }
    
    // fclose(f1);exit(0);
      updateRequests();
      //quadratureT::project(bT::gh,intNN::f<P_DIMS>(0,nSpecies),bNNT::density);
  }
};

#endif 
#endif
