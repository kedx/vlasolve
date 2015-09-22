#ifndef __QUADRATURE_21_HXX__
#define __QUADRATURE_21_HXX__

#include "integrand21.hxx"

template <class gridHandlerType>
class quadrature21
{
private:
  typedef quadrature21<gridHandlerType> myT;
public:
  typedef gridHandlerType gridHandlerT;
  typedef typename gridHandlerT::gridT gridT;

  typedef typename gridT::sliceIterator gridSliceItT;
  typedef typename gridT::value_type dataT;
  typedef typename gridT::valLocationT valLocationT;
  typedef typename gridT::valLocationV valLocationV;
  
  //static const valLocationT valLocation = gridT::valLocation;

  typedef int21::f defaultIntT;
  
  template <class funcT>
  static double integrate(gridHandlerT *gh,const funcT &myFunc, bool global=false)
  {
    double result=integrate<funcT>(gh->getGrid(), myFunc);
    if (global) return gh->sum(result);
    else return result;	
  }

  template <class funcT>
  static std::vector<double> integrate1D(gridHandlerT *gh,const int dir, const funcT &myFunc, bool global=false)
  {
    std::vector<double> result=myT::template integrate1D<funcT>(gh->getGrid(),dir,myFunc);
    if (global) result=gh->gather_all(result.begin(),dir,1);
    else gh->gather(result.begin(),dir,1);	
    
    return result;
  }

  template <class funcT>
  static std::vector<double> integrate1D_acc(gridHandlerT *gh, const int dir, double &tot, const funcT &myFunc, bool global=false)
  {
    std::vector<double> result=myT::template integrate1D<funcT>(gh->getGrid(),dir,myFunc);
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
	double *res=&result[omp_get_thread_num()];

	double J=it.get_J();
	double dJ=it.get_dJ();
	double Ct1=8*M_PI*M_PI*J*dJ;
	
	for (;it!=it_end;++it)
	  {
	    double dU=it.get_dU();if (dU<=0) continue;
	    double dR=it.get_dR();if (dR<=0) continue;

	    (*res)+=Ct1*dU*dR*myFunc(it.get_R()+0.5*dR,it.get_U()+0.5*dU,it.get_J(),it);
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
    //long j;
    
    const long tabSize=g->getValCoord(dir).size();
    //const std::vector<double>& Rtab= g->getValCoord21_R();
    std::vector<double> tmp[num_omp_threads];

    Mrv.assign(tabSize,0);

#pragma omp parallel for 
    for (long j=0;j<num_omp_threads;j++) 
      tmp[j].assign(tabSize,0);
    
#pragma omp parallel for 
    for (long j=0;j<nSlices;j++)
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
    for (long i=0;i<tabSize;i++) 
      for (long j=0;j<num_omp_threads;j++)
	Mrv[i]+=tmp[j][i];
    
    return Mrv; 
  }
  /*
  static std::vector<double> M_of_R(gridT *g)
  {
    return myT::integrate< 0,defaultIntT >(g);
  }
  */
};

#endif
