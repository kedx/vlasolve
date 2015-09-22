#ifndef __QUADRATURE_21_HXX__
#define __QUADRATURE_21_HXX__

template <class gridHandlerType>
class quadrature21
{
public:
  typedef gridHandlerType gridHandlerT;
  typedef typename gridHandlerT::gridT gridT;

  typedef typename gridT::sliceIterator gridSliceItT;
  typedef typename gridT::value_type dataT;
  typedef typename gridT::valLocationT valLocationT;
  typedef typename gridT::valLocationV valLocationV;
  
  static const valLocationT valLocation = gridT::valLocation;

  static double M(gridHandlerT *gh)
  {
    std::vector<double> Mr = M_of_R(gh->getGrid());
    return gh->sum(Mr,0);
  }

  static std::vector<double> M_of_less_than_R(gridHandlerT *gh, double &mTot)
  {
    std::vector<double> Mr = M_of_R(gh->getGrid());
    mTot=gh->accumulate(Mr,0);
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
	gridSliceItT it=g->slice_begin(j);
	const gridSliceItT it_end=g->slice_end(j);

	double J=it.get_J();
	double dJ=it.get_dJ();
	double Ct1=8*M_PI*M_PI*J*dJ;
	
	for (;it!=it_end;++it)
	  {
	    double dU=it.get_dU();if (dU<=0) continue;
	    double dR=it.get_dR();if (dR<=0) continue;
	    double tmp=Ct1*dU*dR*it.get_avg();
	    #pragma omp atomic
	    Mr[it.get_r()+1]+=tmp;	    
	  }
      }
    for (j=0;j<Mr.size();j++) printf("(%ld,%g,%g)",j,Rtab[j],Mr[j]);
    printf("\n");
    return Mr;
  }

  
};

#endif
