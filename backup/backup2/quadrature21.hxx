#ifndef __QUADRATURE_21_HXX__
#define __QUADRATURE_21_HXX__

template <class gridHandlerType>
class quadrature21
{
public:
  typedef gridHandlerType gridHandlerT;
  typedef typename gridHandlerT::gridT gridT;

  typedef typename gridT::sliceIterator gridSliceItT;
  typedef typename gridT::dataT dataT;
  typedef typename gridT::valLocationT valLocationT;
  typedef typename gridT::valLocationV valLocationV;
  
  static const valLocationT valLocation = gridT::valLocation;

  static std::vector<double> M_of_less_than_R(gridHandlerT *gh)
  {
    std::vector<double> Mr = M_of_R(gh->getGrid());
    gh->accumulate(Mr,0);
    return Mr;
  }

private:
// REMOVE ATOMIC FOR SPEEDUP !!!
  static std::vector<double> M_of_R(gridT *g)
  {
    std::vector<double> Mr;
    const std::vector<double>& Rtab= g->getValCoord21_R();
    const long nSlices=g->nSlices();
    long j;

    Mr.assign(Rtab.size(),0);
    
#pragma omp parallel for 
    for (j=0;j<nSlices;j++)
      {
	const gridSliceItT it_begin=g->slice_begin(j);
	double J=it_begin.get_J();
	double dJ=it_begin.get_dJ();
	double Ct1=8*M_PI*M_PI*J*dJ;
	const gridSliceItT it_end=g->slice_end(j);
	for (gridSliceItT it=it_begin;it!=it_end;++it)
	  {
	    double dU=it.get_dU();if (dU<=0) continue;
	    double dR=it.get_dR();if (dR<=0) continue;
	    double tmp=Ct1*dU*dR*it.get_avg();
	    #pragma omp atomic
	    Mr[it.get_r()+1]+=tmp;	    
	  }
      }
    
    return Mr;
  }

  
};

#endif
