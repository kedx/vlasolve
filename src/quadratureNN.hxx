#ifndef __QUADRATURE_NN_HXX__
#define __QUADRATURE_NN_HXX__

#include "integrandNN.hxx"

template <class gridHandlerType>
class quadratureNN
{
private:
  typedef quadratureNN<gridHandlerType> myT;
public:
  typedef gridHandlerType gridHandlerT;
  typedef typename gridHandlerT::gridT gridT;

  typedef typename gridT::sliceIterator gridSliceItT;
  typedef typename gridT::value_type dataT;
  typedef typename gridT::valLocationT valLocationT;
  typedef typename gridT::valLocationV valLocationV;

  static const int P_DIMS= gridT::P_DIMS;
  static const int U_DIMS= gridT::U_DIMS;
  static const int J_DIMS= gridT::J_DIMS;
  static const int DIMS= gridT::DIMS;
  
  static const valLocationT valLocation = gridT::valLocation;

  //typedef gridBase< dimTraits<P_DIMS,0>, dataT, valLocation > projectedSliceT;

  typedef intNN::f<P_DIMS> defaultIntT;

  template <class funcT>
  static double integrate(gridHandlerT *gh,const funcT &myFunc, bool global=false)
  {
    double result=integrate<funcT>(gh->getGrid(), myFunc);
    if (global) return gh->sum(result);
    else return result;	
  }
  
  // projectedSliceT should be a gridBase<> type
  template <class funcT, class subGridT> 
  static void project(gridHandlerT *gh, const funcT &myFunc, subGridT &result, bool full=false, bool global=true, int *which=NULL)
  {
    if (which!=NULL)
      gh->template getSubGrid<subGridT>(result,which,0,global);	
    else if (!result.isInitialized())
      {
	fprintf(stderr,"ERROR in quadratureNN::project: subGrid must be initialized or a 'which' parameter must be given.\n");
	exit(0);
      }
    
    if (global) full=false;
    project<funcT,subGridT>(gh->getGrid(),myFunc,result,full);
   
    if (global) gh->gatherGrids(result);
    
  }
  
  template <class funcT>
  static std::vector<double> integrate1D(gridHandlerT *gh,const int dir, const funcT &myFunc, bool global=false)
  {
    std::vector<double> result=myT::integrate1D<funcT>(gh->getGrid(),dir,myFunc);
    if (global) result=gh->gather_all(result.begin(),dir,1);
    else gh->gather(result.begin(),dir,1);	
    
    return result;
  }

  template <class funcT>
  static std::vector<double> integrate1D_acc(gridHandlerT *gh,const int dir, double &tot, const funcT &myFunc, bool global=false)
  {
    std::vector<double> result=myT::integrate1D<funcT>(gh->getGrid(),dir,myFunc);
    if (global)
      tot=gh->accumulate_all(result.begin(),dir,1).back();
    else
      tot=gh->accumulate(result.begin(),dir,1);
      
    return result;
  }

private:
  
  template <class funcT>
  static double integrate(gridT *g,const funcT &myFunc=funcT())
  {     
    const long nSlices=g->nSlices();
    std::vector<double> result(num_omp_threads,0);
    long j;
    
#pragma omp parallel for 
    for (j=0;j<nSlices;j++)
      {
	gridSliceItT it=g->slice_begin(j);
	const gridSliceItT it_end=g->slice_end(j);
	double p[DIMS];
	double *res=&result[omp_get_thread_num()];

	it.cellC(p);	 	
	for (;it!=it_end;++it)
	  {	
	    it.cellCP(p);
	    double dv=it.cellVol();
	    if (dv==0) continue;
	    (*res)+=myFunc(p,it)*dv;
	  }
      }
    
    for (int i=1;i<result.size();i++) result[0]+=result[i];
     
    return result[0];
  }

  template <class funcT>
  static std::vector<double> integrate1D(gridT *g, const int dir, const funcT &myFunc=funcT())
  {      
    std::vector<double> Mrv;    
    const long nSlices=g->nSlices();
    long j;
    
    const long tabSize=g->getValCoord(dir).size();
    //const std::vector<double>& Rtab= g->getValCoord21_R();
    std::vector<double> tmp[num_omp_threads];
    printf("IntegrandNN::integrate1D not implemented yet!\n");exit(-1);
     /*
    Mrv.assign(tabSize,0);

#pragma omp parallel for 
    for (j=0;j<num_omp_threads;j++) 
      tmp[j].assign(tabSize,0);
    
#pragma omp parallel for 
    for (j=0;j<nSlices;j++)
      {
	double *Mr=&tmp[omp_get_thread_num()][0];
	gridSliceItT it=g->slice_begin(j);
	const gridSliceItT it_end=g->slice_end(j);

	double J=it.get_J();
	double dJ=it.get_dJ();
	double Ct1=8*M_PI*M_PI*J*dJ;
	
	for (;it!=it_end;++it)
	  {
	    double dU=it.get_dU();if (dU<=0) continue;
	    double dR=it.get_dR();if (dR<=0) continue;
	    double tmp=Ct1*dU*dR*myFunc(it.get_R()+0.5*dR,it.get_U()+0.5*dU,it.get_J(),it); 

	    switch (dir)
	      {
	      case 0:
		Mr[it.get_r()]+=tmp;
		break;
	      case 1:
		Mr[it.get_u()]+=tmp;
		break;
	      case 2:
		Mr[it.get_j()]+=tmp;
		break;
	      }
	  }
      }

     #pragma omp parallel for 
    for (j=0;j<num_omp_threads;j++)
      for (long i=0;i<tabSize;i++) Mrv[i]+=tmp[j][i]; 
     */
    return Mrv; 

  }
  
  template <class funcT, class subGridT> 
  static void project(gridT *g, const funcT &myFunc, subGridT &result, bool full)
  {
    int type;
    
    std::pair<int,int> parentD=result.getParentDims();
    if ((subGridT::P_DIMS == P_DIMS)&&(subGridT::U_DIMS == 0)) type=1;
    else if ((subGridT::P_DIMS == 0)&&(subGridT::U_DIMS == U_DIMS)) type=0;
    else
      {
	fprintf(stderr,"ERROR in project: cannot project a [%d,%d] grid on a [%d,%d] subgrid.\n",
		P_DIMS,U_DIMS,
		subGridT::P_DIMS,subGridT::U_DIMS);
	fprintf(stderr,"      NOT IMPLEMENTED YET.\n");
	exit(0);
      }

    const long nSlices=g->nSlices(type,full);    
   
    FILE *f;
  
#pragma omp parallel for 
    for (long j=0;j<nSlices;j++)
      {
	gridSliceItT it=g->slice_begin(j,type,full);
	const gridSliceItT it_end=g->slice_end(j,type,full);
	double p[DIMS];

	it.cellC(p);
		
	typename subGridT::iterator res=result.importIterator(it);
	*res=0;

	for (;it!=it_end;++it)
	  {	
	    double dv=it.slice_cellVol(false);	    
	    if (dv==0) continue;
	    (*res)+=myFunc(p,it,false)*dv;	    
	  }
      }
    
  }
  
};

#endif
