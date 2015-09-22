#ifndef __VLASOV_SOLVER_HXX__
#define __VLASOV_SOLVER_HXX__

#include <math.h>

#include "global.h"
#include "dimTraits.hxx"
#include "snapshots.hxx"

#include <sys/time.h>
#include <vector>

//#include "splineInterface.hxx"
#include "initGen.hxx"
#include "solverParams.hxx"
#include "NDfield.h"
#include "string.h"

#ifdef USE_OPENMP
#include <parallel/algorithm>
#endif

#ifdef USE_THREADS
#include <errno.h>
#include <pthread.h>
#endif

template <class gridType,class interpolT> class vlSolver;

template <class gridType,class interpolT>
class vlSolver {
public:
  
  friend class snapshots< vlSolver<gridType,interpolT> >;
  typedef snapshots< vlSolver<gridType,interpolT> > snapshotT;

  typedef gridType gridT;
  typedef typename gridT::dimTr dimTr;
  typedef typename gridT::dataT dataT;

private:
  snapshotT *snap;

  struct dummyCell
  {
    double R;
    double U;
    dataT val;

    dummyCell(dataT _val=0, double _R=0, double _U=0)
    {
      R=_R;
      U=_U;
      val=_val;
    }
  };
  
  typedef std::vector<dummyCell> dCellChain; 
 
public:

  static const int DIMS = gridT::DIMS;
  static const int P_DIMS = gridT::P_DIMS;
  static const int V_DIMS = gridT::V_DIMS;
  static const int J_DIMS = gridT::J_DIMS;
  
  typedef solverParams paramT;
  typedef initGen<gridT> initT;
  typedef typename gridT::paramT gridParamT;
  
  
  
private:  

  paramT p;
  gridT *grid;
  gridParamT gp; 
  dataT *data;
  typename interpolT::initT interpP; 
  std::vector<bool> mask;
  std::vector<double> timeStep;
 
  std::vector< std::map<long,dCellChain> > coreCells;
  std::vector<double> disp_kick;
  std::vector<double> disp_driftR;
  std::vector<double> disp_driftU;  
 

  std::vector<std::vector<double> > P_Cell;
  std::vector<std::vector<double> > P_Vert;
  std::vector<std::vector<double> > U_Cell;
  std::vector<std::vector<double> > U_Vert;

  std::pair<double,double> advanceParticleDrift(double R, double U, double J,double dt, double R0=-1)
  {
    double H0; 
    double newU;
    double newR;
    double tlim,tR0;
    double sgn;
    double deltaT;
    
    H0= 0.5* (U*U + (J*J)/(R*R));
    sgn= (U<0)?-1:1;
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
	//newR = sqrt(pow(fabs(newR*newU)-2*sgn*H0*(dth-tR0),2)+J*J)*pow(2*H0,-0.5);	       
	
	//printf("J=%f HI R=%f U=%f (Tr0=%f/tlim=%f)",J,R,U,tR0,tlim);
	//printf(" --> R=%f U=%f\n",newR,sgn*sqrt((2*H0-(J*J)/(newR*newR))));
	
	deltaT -= tR0;
	//tlim -= tR0;
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
    newU = sgn*sqrt(newU);

    return std::make_pair(newR,newU);
  }

  std::vector<double> EvaluateDt(double dt)
  {
    std::vector<double> &R_C=P_Cell[0];
    std::vector<double> &U_C=U_Cell[0];
    std::vector<double> &J_C=U_Cell[1];
    std::vector<double> &R_V=P_Vert[0];
    std::vector<double> &U_V=U_Vert[0];
    std::vector<double> &J_V=U_Vert[1];
    double R,U,J;
    long r,u,j,i;
    double R0=R_C[0];
    long delta = grid->getCellStride(2);
    double dth=dt/2;

    std::vector<double> deltaT(J_C.size(),dth);
    
    typename interpolT::interpT *interp;
    dataT *data = grid->getDataPtr();

    mask.assign(grid->getNCell(),false);

    printf("Evaluating Dt ...");fflush(0);

#pragma omp parallel for shared(R_C,U_C,J_C,R_V,U_V,J_V,deltaT) private(J,U,R,j,u,r,i,interp)  
    for (j=0;j<J_C.size();j++)
      {
	std::pair<double,double> newVal;
	i=j*delta;
	interp = interpolT::create(interpP,&data[i]);
	J=J_C[j];
	for (u=0;u<U_C.size();u++)
	  {
	    U=U_C[u];
	    for (r=0;r<R_C.size();r++,i++) 
	      {
		double curdt = dth;
		long nsteps=1;
		//if ((mask.size())&&(mask[i])) continue;

		R=R_C[r];
		newVal=advanceParticleDrift(R,U,J,curdt,R0);
		double Dr = fabs( (newVal.first-R_C[r])/(R_V[r+1]-R_V[r]) );
		double Du = fabs( (fabs(newVal.second)-fabs(U_C[u]))/(U_V[u+1]-U_V[u]) );
		
		if (Dr>1)
		  {
		    do {
		      curdt*=0.5;nsteps++;
		      newVal=advanceParticleDrift(R,U,J,curdt,R0);
		      if ((newVal.first>=R_C.back())||
			    (newVal.second<=U_C[0])||
			    (newVal.second>=U_C.back())) break;
		      Dr = fabs( (newVal.first-R)/(R_V[r+1]-R_V[r]) );
		    } while (Dr>1);
		  }
		
		if (Du>1)
		  {
		    do {
		      curdt*=0.5;nsteps++;
		      newVal=advanceParticleDrift(R,U,J,curdt,R0);
		      if ((newVal.first>=R_C.back())||
			    (newVal.second<=U_C[0])||
			    (newVal.second>=U_C.back())) break;
		      Du = fabs( (fabs(newVal.second)-fabs(U))/(U_V[u+1]-U_V[u]) );
		    } while (Du>1);
		  }

		if ((curdt<dth)&&(data[i]<1.E-10))
		  {
		    double v,t;
		    int step;

		    mask[i]=true;
		    for (step=1;step<=nsteps;step++)
		      {
			v=data[i];
			t=step*curdt;
			newVal=advanceParticleDrift(R,U,J,t,R0);
			if ((newVal.first>=R_C.back())||
			    (newVal.second<=U_C[0])||
			    (newVal.second>=U_C.back())) break;
			else interpolT::evaluate (interp,newVal.first,newVal.second,&v);

			if (v>1.E-5) {
			  mask[i]=false;	       
			  break;
			}

			/*
			if ((data[i]==0)&&(v==0)) continue; 
			if (fabs(2*(data[i]-v)/(data[i]+v)) > 1.E-20)  
			  {
			    mask[i]=false;
			    break;
			  }
			*/
		      }		     
		  }
		
		if (!mask[i]) 
		  if (deltaT[j]>curdt)
		    {
		      //newVal=advanceParticleDrift(R,U,J,dth/pow(2,10),R0);
		      //Dr = fabs(newVal.first-R_C[r])/(R_V[r+1]-R_V[r]);
		      //Du = fabs(fabs(newVal.second)-fabs(U_C[u]))/(U_V[u+1]-U_V[u]);
		      //printf("%e (%ld steps): RU=(%g,%g) NRU=(%g,%g) DRU=(%g,%g)\n",curdt,nsteps,R,U,newVal.first,newVal.second,Dr,Du);
		      deltaT[j]=curdt;
		    }
	      }
	  }

	interpolT::destroy(interp);
      }

    for (j=0;j<J_C.size();j++) deltaT[j]*=2;
    
    printf("done.\n");

    return deltaT;
  }

  
  void generateDisplacement(double dt, initT *ini=NULL, bool reflectiveKernel=true)
  {
    std::vector<double> &R_C=P_Cell[0];
    std::vector<double> &U_C=U_Cell[0];
    std::vector<double> &J_C=U_Cell[1];
    std::vector<double> &R_V=P_Vert[0];
    std::vector<double> &U_V=U_Vert[0];
    std::vector<double> &J_V=U_Vert[1];
    long delta = grid->getCellStride(2);

    long r,u,j,i;
    double R,U,J;
    double newR;
    double newU;
    double dth=dt/2;
    double R0=(reflectiveKernel)?R_C[0]:-1;
    
    double normFactor;
    if (!grid->getRegisteredValue("INIGEN_NormFactor",normFactor)) normFactor=1;
    

    disp_driftR.resize(grid->getNCell());
    disp_driftU.resize(grid->getNCell());
    
    coreCells.clear();
    //if (!reflectiveKernel) 
    coreCells.resize(J_C.size());
    
    FILE *f=fopen("testXY.i","wb");
    FILE *f2=fopen("testUV.i","wb");
    FILE *f3=fopen("maxRU.dat","wb");
    fprintf(f,"xy=[];\n");
    fprintf(f2,"uv=array(double,[2,4,%ld]);id=1;\n",R_C.size()*U_C.size());

    
    std::vector<double> maxDr(J_C.size());
    std::vector<double> maxDu(J_C.size());

//#pragma omp parallel for shared(R_C,U_C,J_C,R_V,U_V,J_V) private(J,U,R,j,u,r,i)  
    for (j=0;j<J_C.size();j++)
      {
	std::pair<double,double> newVal;
	i=j*delta;
	J=J_C[j];	
	
	for (u=0;u<U_C.size();u++)
	  {
	    U=U_C[u];
	    
	    for (r=0;r<R_C.size();r++,i++) 
	      {
		R=R_C[r];
		double deltaT=dth;
		
		newVal=advanceParticleDrift(R,U,J,deltaT,R0);
		disp_driftR[i]=newVal.first;
		disp_driftU[i]=newVal.second;	
		
		double Dr = (disp_driftR[i]-R_C[r])/(R_V[r+1]-R_V[r]);
		double Du = (disp_driftU[i]-U_C[u])/(U_V[u+1]-U_V[u]);
		if (fabs(Dr)>fabs(maxDr[j])) maxDr[j]=Dr;
		if (fabs(Du)>fabs(maxDu[j])) maxDu[j]=Du;

		//fprintf(f3,"%ld %g %g %g %g %g\n",j,R_C[r],U_C[u],J_C[j],Dr,Du);

		if (j==0) fprintf(f2,"uv(,id++)=[%g,%g,%g,%g];\n",R,U,disp_driftR[i],disp_driftU[i]);
		if ((!reflectiveKernel)&&(newVal.first<R_C[0]))
		  {
		    // generate dummy cells
		    double val;		    
		    dCellChain &curChain = 
		    coreCells[j].insert(std::make_pair(i,dCellChain())).first->second;
		    //dCellChain curChain(0);

		    //if (j==0) printf("%ld(%g/%g) \n",i,disp_driftR[i],newVal.first);
		    std::vector<double> pos(3);
		    pos[2]=J;	
		    if (j==0) fprintf(f,"grow,xy,&[");
		    while (newVal.first<R_C[0])
		      {
			pos[0]=newVal.first;pos[1]=newVal.second;
			val=normFactor*ini->valueAt(pos);
			if (j==0) fprintf(f,"[%g,%g,%g,%g],",pos[0],pos[1],pos[2],val);
			curChain.push_back(dummyCell(val,newVal.first,newVal.second));
			deltaT+=dth;
			newVal=advanceParticleDrift(R,U,J,deltaT,R0);
		      }	
		    // For that cell,  value will be computed when solving --> 0
		    curChain.push_back(dummyCell(0,newVal.first,newVal.second));
		    if (j==0) fprintf(f,"[%g,%g,%g,%g]]\n",newVal.first,newVal.second,J,0.);
		    //if (j==0) fprintf(f,"plg,(*xy(0))(2,),(*xy(0))(1,);\n");
		    //coreCells[j].insert(std::make_pair(i,curChain));
		    //if (j==0) printf("  SIZE = %ld\n",curChain.size());
		  }
	      }
	  }
	//if (j==0) printf("\n\n");
      }
    //exit(0);
    for (j=0;j<J_C.size();j++)
      fprintf(f3,"%g %g %g\n",J_C[j],maxDr[j],maxDu[j]);

    fclose(f);fclose(f2);fclose(f3);
  }
  
  void init(initT *ini) {
    grid=ini->generate(p);
    data=grid->getDataPtr();
    gp=grid->get_Params();

    P_Cell = grid->getCellCoord_P();
    P_Vert = grid->getVertCoord_P();
    U_Cell = grid->getCellCoord_U();
    U_Vert = grid->getVertCoord_U();  

    std::vector<double> &R_C=P_Cell[0];
    std::vector<double> &U_C=U_Cell[0];
    std::vector<double> &J_C=U_Cell[1];
    std::vector<double> &R_V=P_Vert[0];
    std::vector<double> &U_V=U_Vert[0];
    std::vector<double> &J_V=U_Vert[1];

    interpP=interpolT::createInit(R_V.front(),R_V.back(),R_V.size(),
				  U_V.front(),U_V.back(),U_V.size(),
				  interpolT::boundaryTypeVal::BT_NATURAL,
				  interpolT::boundaryTypeVal::BT_FLAT,
				  interpolT::boundaryTypeVal::BT_FLAT,
				  interpolT::boundaryTypeVal::BT_FLAT,
				  gp.Pscale[0],gp.Vscale[0],
				  interpolT::valLocationTypeVal::PIXEL);
  
    
    //generateDisplacement(p.dt,ini,false); // no reflective kernel  

    timeStep = EvaluateDt(p.dt);
    FILE *f=fopen("deltaT.dat","wb");
    for (long j=0;j<J_C.size();j++)
      {
	fprintf(f,"%10.10g %10.10g\n",J_C[j],timeStep[j]);
	printf("%10.10g %10.10g\n",J_C[j],timeStep[j]);
      }
    fclose(f);
    generateDisplacement(p.dt,ini);  

    snap=new snapshotT(this);   
    snap->addSliceTrigger(p.T0,p.T0-1,0.25,10);    
    snap->addTrigger(p.T0,p.T0-1,0.5);    
  
  }

  double massBelowR0;
  double R0;
  
public:
   
  vlSolver(initT *ini, paramT *params=NULL )
  {   
    if (params==NULL) p=ini->defaultSolverParams();
    else p=*params;
        
    init(ini);
  }

  ~vlSolver() {
    delete grid;
    delete snap;
  }

  void solve()
  {    
    double t=p.T0;
    double dt=p.dt;

    initSolver<dimTr>();
    
    printf("Starting at T=%.5f\n",t);fflush(0);
    snap->update(t);
    do {
      makeOneStep<dimTr>(t,dt);
      t+=dt;           
      printf("\rT=%.3f   ",t);fflush(0);
      snap->update(t);
    } while ((t<p.Tmax)||(p.Tmax==0));
  }
  /*
  template <class splineT>
  void testSpline(typename splineT::splineT *spline, long j, long Nu=10, long Nr=10)
  {
    
    std::vector<double> &R_=Pcoord[0];
    std::vector<double> &U_=Vcoord[0];
    std::vector<double> &J_=Vcoord[1];
    
    int ndims=2;
    int dims[ndims];
    double x0[ndims];
    double delta[ndims];
    dims[0]=(Nr*(gp.Pres[0]-1)+1);
    dims[1]=(Nu*(gp.Vres[0]-1)+1);
    x0[0]=gp.Pmin[0];delta[0]=gp.Pmax[0]-gp.Pmin[0];  
    x0[1]=gp.Vmin[0];delta[1]=gp.Vmax[0]-gp.Vmin[0];    
    dataT *test=new dataT[dims[0]*dims[1]];
    NDfield *f=Create_NDfield(dims,ndims,0,ND_DOUBLE,x0,delta,test,"radial density");

    long u,uu,r,rr,i;
    double U,R;
    i=0;
    for (u=0;u<gp.Vres[0]-1;u++)
      {
	if (u==gp.Vres[0]-2) Nu++;
	for (uu=0;uu<Nu;uu++)
	  {
	    U=U_[u] + uu * (U_[u+1]-U_[u])/(Nu);
	    for (r=0;r<gp.Pres[0]-1;r++)
	      {
		if (r==gp.Pres[0]-2) Nr++;
		for (rr=0;rr<Nr;rr++,i++)
		  {
		    R=R_[r] + rr * (R_[r+1]-R_[r])/(Nr);
		    splineT::evaluate(spline,R,U,&test[i]);
		  }
		if (r==gp.Pres[0]-2) Nr--;
	      }
	  }
	if (u==gp.Vres[0]-2) Nu--;
      }

    char fname[255];
    sprintf(fname,"Test-interpolation.ND");
    Save_NDfield(f,fname);
    Free_NDfield(&f);
    
    f=NDfieldFromData(j);
    sprintf(fname,"Test-original.ND");
    Save_NDfield(f,fname);
    Free_NDfield(&f);
    
  }
    */

private:
  template <class DTr> 
  void initSolver(){
    initSolver<DTr>(DTr());
  };

  template <class DT>
  void initSolver(dimTraits<1,2>){
    massBelowR0=0;
    grid->getRegisteredValue("R0",R0);
    grid->getRegisteredValue("massBelowR0",massBelowR0);
  };

  template <class DTr> 
  void makeOneStep(double t, double dt){
    makeOneStep<DTr>(t,dt,DTr());
  };
  
  template <class DT>
  void makeOneStep(double t, double dt,dimTraits<1,2>)
  {
    //Ugrid p_grid, v_grid;
    //BCtype_d p_bc, v_bc;
    dataT *data = grid->getDataPtr();
    /*******/
    double G=1; // A CHANGER !!!!
    /*******/    
    
    long r,u,j,i;
    double R,U,J;
    double newR;
    double newU;
    double dth=0.5*dt;

    
    std::vector<double> &R_C=P_Cell[0];
    std::vector<double> &U_C=U_Cell[0];
    std::vector<double> &J_C=U_Cell[1];
    std::vector<double> &R_V=P_Vert[0];
    std::vector<double> &U_V=U_Vert[0];
    std::vector<double> &J_V=U_Vert[1];
    long delta = grid->getCellStride(2);

    //typedef einspline<2,dataT> splineT;
    typename interpolT::interpT *interp;
    /*
    typename interpolT::initT interpP =
      interpolT::createInit(R_V.front(),R_V.back(),R_V.size(),
			    U_V.front(),U_V.back(),U_V.size(),
			    interpolT::boundaryTypeVal::BT_NATURAL,
			    interpolT::boundaryTypeVal::BT_FLAT,
			    interpolT::boundaryTypeVal::BT_FLAT,
			    interpolT::boundaryTypeVal::BT_FLAT,
			    gp.Pscale[0],gp.Vscale[0],
			    interpolT::valLocationTypeVal::PIXEL);
    */

    /*
    j=0;
    i=j*gp.Vres[0]*gp.Pres[0];
    spline = splineT::create(splineP,&data[i]);
    testSpline<splineT>(spline,j);
    exit(0);
    */    

#pragma omp parallel for shared(R_C,U_C,J_C,R_V,U_V,J_V,data) private(J,U,R,j,u,r,i,interp,newR,newU)
    for (j=0;j<J_C.size();j++)
      {
	i=j*delta;
	const long imax=i+U_C.size()*R_C.size();
	interp = interpolT::create(interpP,&data[i]);

	if (coreCells[j].size())
	  {
	    typename std::map<long,dCellChain>::iterator chain_it;
	    typename dCellChain::iterator it;
	    for (chain_it=coreCells[j].begin();chain_it!=coreCells[j].end();chain_it++)
	      {
		dummyCell &c = chain_it->second.back();
		if ((c.U<=U_C[0])||(c.U>=U_C.back())) c.val=0;
		else interpolT::evaluate (interp,c.R,c.U,&c.val);
	      }
	  }

	for (;i<imax;i++)
	  {	
	    if (disp_driftR[i]<R_C[0]) data[i]=coreCells[j].find(i)->second[0].val;
	    else if ((disp_driftR[i]>=R_C.back())||
		     (disp_driftU[i]<=U_C[0])||
		     (disp_driftU[i]>=U_C.back())) data[i]=0;
	    else interpolT::evaluate (interp,disp_driftR[i],disp_driftU[i],&data[i]);
	  }

	interpolT::destroy(interp);

	if (coreCells[j].size())
	  {
	    typename std::map<long,dCellChain>::iterator chain_it;
	    for (chain_it=coreCells[j].begin();chain_it!=coreCells[j].end();chain_it++)
	      {
		for (u=0;u<chain_it->second.size()-1;u++)
		  chain_it->second[u]=chain_it->second[u+1];
	      }
	  }
      }

    std::vector<double> Mr;
    grid->computeM_R(Mr,massBelowR0);
    for (r=0;r<R_C.size();r++)  
      Mr[r] = (0.5*(Mr[r]+Mr[r+1])*G/(R_C[r]*R_C[r]))*dt;
    
#pragma omp parallel for shared(R_C,U_C,J_C,R_V,U_V,J_V,data) private(J,U,R,j,u,r,i,interp,newR,newU)
    for (j=0;j<J_C.size();j++)
      {
	i=j*delta;
	J=J_C[j];
	interp = interpolT::create(interpP,&data[i]);
	for (u=0;u<U_C.size();u++)
	  {
	    U=U_C[u];
	    for (r=0;r<R_C.size();r++,i++) 
	      {
		//R=R_C[r];
		//newU = U+(0.5*(Mr[r]+Mr[r+1])*G/(R*R))*dt;
		newU = U+Mr[r];
		newR = R_C[r];
		if ((newU<=U_C[0])||(newU>=U_C.back())) data[i]=0;
		else interpolT::evaluate (interp,newR,newU,&data[i]);
		//if (data[i]<0) data[i]=0;
	      }
	  }
	interpolT::destroy(interp);
      }

#pragma omp parallel for shared(R_C,U_C,J_C,R_V,U_V,J_V,data) private(J,U,R,j,u,r,i,interp,newR,newU)
    for (j=0;j<J_C.size();j++)
      {	
	i=j*delta;
	const long imax=i+U_C.size()*R_C.size();
	interp = interpolT::create(interpP,&data[i]);

	if (coreCells[j].size())
	  {
	    typename std::map<long,dCellChain>::iterator chain_it;
	    typename dCellChain::iterator it;
	    for (chain_it=coreCells[j].begin();chain_it!=coreCells[j].end();chain_it++)
	      {
		dummyCell &c = chain_it->second.back();
		if ((c.U<=U_C[0])||(c.U>=U_C.back())) c.val=0;
		else interpolT::evaluate (interp,c.R,c.U,&c.val);
	      }
	  }

	for (;i<imax;i++)
	  {
	    if (disp_driftR[i]<R_C[0]) data[i]=coreCells[j].find(i)->second[0].val;
	    else if ((disp_driftR[i]>=R_C.back())||
		(disp_driftU[i]<=U_C[0])||
		(disp_driftU[i]>=U_C.back())) data[i]=0;
	    else interpolT::evaluate (interp,disp_driftR[i],disp_driftU[i],&data[i]);
	  }

	interpolT::destroy(interp);

	if (coreCells[j].size())
	  {
	    typename std::map<long,dCellChain>::iterator chain_it;
	    for (chain_it=coreCells[j].begin();chain_it!=coreCells[j].end();chain_it++)
	      {
		for (u=0;u<chain_it->second.size()-1;u++)
		  chain_it->second[u]=chain_it->second[u+1];
	      }
	  }

      }
  }

};



#endif
