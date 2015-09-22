#ifndef __SNAPSHOTS_HXX__
#define __SNAPSHOTS_HXX__

#include <iostream>
#include <queue>
#include <functional>
#include <algorithm>
#include <set>
#include <map>

#include "NDfield.h"
#include "scale.hxx"

template <class vlSolver>
class snapshots {
public:
  typedef typename vlSolver::gridT gridT;
  typedef typename vlSolver::dataT dataT;
  typedef scale<double> scaleT;
  
  /*
  struct dumpTypeVal {
    enum type {FULL=0, SLICE=1};
  };

  typedef typename dumpTypeVal::type dumpType;
*/
  struct sliceType {
    std::set<int> J;
  };

  struct regDumpType {
    double t0;
    double t1;
    double dt;
  };

private:
  vlSolver *solver;

  std::set<double> fullDump;  
  std::map<double, sliceType> sliceDump;

  std::vector<regDumpType> fullDump_reg;
  std::vector<std::pair<regDumpType,sliceType> > sliceDump_reg;

  double lastUpdate;

  long sliceCount;
  long fullCount;
  
  char fullFName[255];
  char densityFName[255];
  char sliceFName[255];

public:  
  snapshots(vlSolver *sl)
  {
    solver=sl;
    sliceCount=0;
    fullCount=0;
    sprintf(fullFName,"fullDump");
    sprintf(sliceFName,"sliceDump");
    sprintf(densityFName,"density");

    lastUpdate=-1.E100;
  }

  template <class inputIteratorT>
  void addTrigger(inputIteratorT Tstart, inputIteratorT Tstop)
  {
    for (inputIteratorT it = Tstart; it != Tstop; it++)
      {
	fullDump.insert((double) *it);
      }
  }
  
  void addTrigger(double t0, double t1, double dt)
  {
    regDumpType t={t0,t1,dt};
    fullDump_reg.push_back(t);
  }

  template <class inputIteratorT,class inputIteratorS>
  void addSliceTrigger(inputIteratorT Tstart, inputIteratorT Tstop,inputIteratorS Sstart, inputIteratorS Sstop)
  {
    sliceType t;
    
    for (inputIteratorS itS = Sstart; itS != Sstop; itS++)
      t.insert(*itS);

    for (inputIteratorT it = Tstart; it != Tstop; it++)
      {
	sliceDump.push(std::make_pair((double)*it,t));
      }
  }

  template <class inputIteratorS>
  void addSliceTrigger(double t0, double t1, double dt,inputIteratorS Sstart, inputIteratorS Sstop)
  {
    sliceType s;    
    regDumpType t={t0,t1,dt};

    for (inputIteratorS itS = Sstart; itS != Sstop; itS++)
      s.J.insert(*itS);
    
    sliceDump_reg.push_back(std::make_pair(t,s));
  }  

  void addSliceTrigger(double t0, double t1, double dt, int nSlices)
  {
    sliceType s;    
    long ntot = solver->U_Cell[solver->V_DIMS-1].size();
    regDumpType t={t0,t1,dt};

    double d,j;
    long i;

    for (i=0;i<nSlices;i++) {
      modf (i * (double)(ntot-1)/(nSlices-1) , &j);
      s.J.insert((int)j);
    }

    sliceDump_reg.push_back(std::make_pair(t,s));
  }  

  void update(double t)
  {
    bool newFull=false;
    bool newSlice=false;
        
    if (fullDump.begin()!=fullDump.end())
    while (*fullDump.rbegin() <= t)
      {
	newFull=true;
	doFullDump(t);
	fullDump.erase(fullDump.rbegin().base());	
      };    

    if (sliceDump.begin()!=sliceDump.end())
    while (sliceDump.rbegin()->first <= t)
      {
	newSlice=true;
	doSliceDump(t,sliceDump.rbegin()->second);
	sliceDump.erase(sliceDump.rbegin().base());
      };    

    for (typename std::vector<regDumpType>::iterator it
	   =fullDump_reg.begin();it != fullDump_reg.end(); it++)
      {
	if (((t>=it->t0)&&(t<=it->t1))
	    ||(it->t0>it->t1)) { 
	  double intpart1;
	  double intpart2;
	  modf ((t-it->t0)/(it->dt) , &intpart1);
	  modf ((lastUpdate-it->t0)/(it->dt) , &intpart2);
	  if (intpart1!=intpart2)
	    {
	      newFull=true;
	      doFullDump(t);
	    }
	}
      }
	
    for (typename std::vector< std::pair<regDumpType,sliceType> >::iterator it
	   =sliceDump_reg.begin();it != sliceDump_reg.end(); it++)
      {
	if (((t>=it->first.t0)&&(t<=it->first.t1))
	    ||(it->first.t0>it->first.t1)) { 
	  double intpart1;
	  double intpart2;
	  modf ((t-it->first.t0)/(it->first.dt) , &intpart1);
	  modf ((lastUpdate-it->first.t0)/(it->first.dt) , &intpart2);
	  if (intpart1!=intpart2)
	    {
	      newSlice=true;
	      doSliceDump(t,it->second);
	    }
	}
      }
	
    if (newFull) fullCount++;
    if (newSlice) sliceCount++;

    lastUpdate = t;
  }

private:

  void doFullDump(double t)
  {
    char fname[255];
    NDfield *f=NDfieldFromData();
    sprintf(fname,"%s_%5.5ld_T%4.4f.ND",fullFName,fullCount,t);
    Save_NDfield(f,fname);
    free(f);

    f=NDfieldFromDensity();
    sprintf(fname,"%s_%5.5ld_T%4.4f.ND",densityFName,fullCount,t);
    Save_NDfield(f,fname);
    Free_NDfield(&f);
  }

  void doSliceDump(double t, sliceType &sl)
  {
    char fname[255];
    std::set<int> &Ji = sl.J;
    for (std::set<int>::iterator it=Ji.begin(); it!=Ji.end();it++) 
      {
	NDfield *f=NDfieldFromData(*it);
	double J = solver->U_Cell[solver->V_DIMS-1][*it];
	
	sprintf(fname,"%s_%5.5ld_J%5.5f_T%4.4f.ND",
		sliceFName,sliceCount,J,t);
	Save_NDfield(f,fname);
	free(f);
      }
  }

  NDfield *NDfieldFromDensity()
  {
    std::vector<double> rhoV;
    solver->grid->computeRho_R(rhoV,solver->massBelowR0);
    double *rho = (double *)malloc(sizeof(double)*rhoV.size());
    std::copy(rhoV.begin(),rhoV.end(),rho);
    //memcpy(rho,&rhoV[0],sizeof(double)*rhoV.size());

    int dims[1];
    double x0[1];
    double delta[1];
    int ndims;
    
    dims[0]=rhoV.size();    
    x0[0]=solver->gp.Pmin[0];
    delta[0]=solver->gp.Pmax[0]-solver->gp.Pmin[0];  
    ndims=1;

    return Create_NDfield(dims,ndims,0,ND_DOUBLE,x0,delta,rho,"radial density");
  }


  NDfield *NDfieldFromData(int sliceAtJ=-1)
  {
    long i;
    int dims[solver->P_DIMS+solver->V_DIMS];
    double x0[solver->P_DIMS+solver->V_DIMS];
    double delta[solver->P_DIMS+solver->V_DIMS];
    int ndims;
    long stride[solver->P_DIMS+solver->V_DIMS+1];
    int index;
    dataT *d;

    stride[0]=1;
    for (i=0;i<solver->P_DIMS;i++) {
      dims[i]=solver->gp.Pres[i];
      x0[i]=solver->gp.Pmin[i];
      delta[i]=solver->gp.Pmax[i]-solver->gp.Pmin[i];
      stride[i+1]=stride[i]*dims[i];
    }
    long dp=solver->P_DIMS;
    for (i=0;i<solver->V_DIMS;i++) {
      dims[i+dp]=solver->gp.Vres[i];
      x0[i+dp]=solver->gp.Vmin[i];  
      delta[i+dp]=solver->gp.Vmax[i]-solver->gp.Vmin[i];
      stride[i+dp+1]=stride[i+dp]*dims[i+dp];
    }
    
    if (sliceAtJ<0) {
      d=solver->data;
      ndims=solver->P_DIMS+solver->V_DIMS;
    }
    else {
      ndims=solver->P_DIMS+solver->V_DIMS-1;
      d=&solver->data[sliceAtJ*stride[ndims]];
    }
    
    return Create_NDfield(dims,ndims,0,ND_DOUBLE,x0,delta,d,"no comment");
  }

};

#endif
