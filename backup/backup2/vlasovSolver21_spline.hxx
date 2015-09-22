#ifndef __VLASOV_SOLVER21_SPLINE_HXX__
#define __VLASOV_SOLVER21_SPLINE_HXX__

#include "vlasovSolver21.hxx"
#include "quadrature21.hxx"
#include "interpol_spline.hxx"
#include "dimTraits.hxx"

template <class gridHandlerType>
class vlasovSolver21_spline : public vlasovSolver21<gridHandlerType> {
public:

  static const std::string getTag() {return "VLASOV_SOLVER21_SPLINE v0.1";}
  
  typedef vlasovSolver21<gridHandlerType> b21T;

  typedef typename b21T::bT bT;
  typedef typename b21T::kernelTypeT kernelTypeT;

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

  vlasovSolver21_spline(mpiComT &mpiCom_) : 
    b21T(mpiCom_)
  {
    
  }

  virtual ~vlasovSolver21_spline()
  {
    
  }  

private:
  typedef typename b21T::cellChainListItT cellChainListItT;
  typedef typename b21T::kernelCellsItT kernelCellsItT;
  typedef typename b21T::dummyCellT dummyCellT;

  void setup()
  {
    boundaryTypeT bType_Rhigh,bType_Rlow;
    boundaryTypeT bType_Uhigh,bType_Ulow;

    bType_Rlow = boundaryTypeV::BT_NATURAL;

    if (bT::gh->isFullGridBoundary_Phigh(0))
      bType_Rhigh = boundaryTypeV::BT_FLAT;
    else
      bType_Rhigh = boundaryTypeV::BT_NATURAL;

    if (bT::gh->isFullGridBoundary_Ulow(0))
      bType_Ulow = boundaryTypeV::BT_FLAT;
    else
      bType_Ulow = boundaryTypeV::BT_NATURAL;

    if (bT::gh->isFullGridBoundary_Uhigh(0))
      bType_Uhigh = boundaryTypeV::BT_FLAT;
    else
      bType_Uhigh = boundaryTypeV::BT_NATURAL;
    
    interpInit = interpolT::createInit(b21T::R_V.front(),b21T::R_V.back(),b21T::R_C.size(),
				       b21T::U_V.front(),b21T::U_V.back(),b21T::U_C.size(),
				       bType_Rlow,bType_Rhigh,
				       bType_Ulow,bType_Uhigh,
				       b21T::Rscale,b21T::Uscale,
				       bT::grid->valLocation);  
  }

  
protected:
  typedef interpolSpline<dataT,scaleT> interpolT;
  typedef typename interpolT::interpT interpT;
  typedef typename interpolT::initT interpInitT;
  typedef typename interpolT::boundaryTypeT boundaryTypeT;
  typedef typename interpolT::boundaryTypeV boundaryTypeV;
  //typedef typename interpolT::valLocationV valLocationV;

  interpInitT interpInit;

  virtual void setup(initGenT *init,const paramsParser &params)
  {
    b21T::setup(init,params);
    setup();
  }

  virtual void read(FILE *f, bool swap)
  {
    b21T::read(f,swap);
    myIO::checkTag(f,getTag());
    setup();
  }

  virtual void write(FILE *f)
  {
    b21T::write(f);
    myIO::writeTag(f,getTag());
  }


  void makeOneStep(double t)
  {   
    long r,u,j,i;
    //double newR;
    //double newU;
    double dth=0.5*bT::deltaTime;
    long delta = bT::grid->getCellStride(2);
    const long nSlices=bT::grid->nSlices();

#pragma omp parallel for 
    for (j=0;j<nSlices;j++)
      { 		
	interpT *interp = interpolT::create(interpInit,bT::grid->fullSlicePtr(j));
	
	if (b21T::kernelType==b21T::kernelTypeV::DELAYED)
	  {
	    //long tmp=0;
	    const kernelCellsItT kit_end = b21T::kernelBoundary[j].end();
	    for (kernelCellsItT kit=b21T::kernelBoundary[j].begin();kit!=kit_end;++kit)
	      {
		//if ((j==0)&&(tmp==0)) kit->second->print("\nbefore:","\n");
		dummyCellT &c = kit->second->back();
		if ((c.U<=b21T::U_Val.front())||(c.U>=b21T::U_Val.back())) c.val=0;
		else interpolT::evaluate(interp,c.R,c.U,&c.val);
		//if ((j==0)&&((tmp++)==0)) kit->second->print("\n\nbefore1:","\n\n");
	      }
	  }

	const gridSliceItT it_end=bT::grid->slice_end(j);
	for (gridSliceItT it=bT::grid->slice_begin(j);it!=it_end;++it)
	  {
	    long i=it.get_i();
	    if (b21T::disp_driftR[i]<b21T::R_Val[0]) *it=b21T::kernelBoundary[j].find(i)->second->readRoll().val;
	    else if ((b21T::disp_driftR[i]>=b21T::R_Val.back())||
		     (b21T::disp_driftU[i]<=b21T::U_Val.front())||
		     (b21T::disp_driftU[i]>=b21T::U_Val.back())) (*it)=0;
	    else interpolT::evaluate(interp,b21T::disp_driftR[i],b21T::disp_driftU[i],&(*it));
	  }
	/*
	if (b21T::kernelType==b21T::kernelTypeV::DELAYED)
	  {
	    long tmp=0;
	    const kernelCellsItT kit_end = b21T::kernelBoundary[j].end();
	    for (kernelCellsItT kit=b21T::kernelBoundary[j].begin();kit!=kit_end;++kit,++tmp)
	      {
		if ((j==0)&&((tmp++)==0)) kit->second->print("after1:","\n");
	      }
	  }
	*/
	interpolT::destroy(interp);
      }
    
    //printf("data = %e\n",data[180*80+10]);
    
    //std::vector<double> Mr = quadratureT::M_of_R(bT::grid);
    //bT::gh->accumulate(Mr,0);
    std::vector<double> Mr = quadratureT::M_of_less_than_R(bT::gh);
    for (j=0;j<Mr.size();j++) Mr[j]+=b21T::kernelMass;
    const std::vector<double> &R = bT::grid->getValCoord21_R();
    
    if (bT::mpiCom->size() == bT::mpiCom->rank()+1) printf("M = %e\n",Mr.back());
    for (i=0;i<Mr.size();i++) Mr[i]*=bT::sp.G/(R[i]*R[i])*bT::deltaTime;
    
    
#pragma omp parallel for
    for (j=0;j<nSlices;j++)
      {
	const gridSliceItT it_begin=bT::grid->slice_begin(j);
	const gridSliceItT it_end=bT::grid->slice_end(j);

	interpT *interp = interpolT::create(interpInit,bT::grid->fullSlicePtr(j));
	
	for (gridSliceItT it=it_begin;it!=it_end;++it)
	  {
	    double newU=it.get_U()+Mr[it.get_r()];
	    double newR=it.get_R();
	    if ((newU<=b21T::U_Val.front())||(newU>=b21T::U_Val.back())) (*it)=0;
	    else interpolT::evaluate(interp,newR,newU,&(*it));
	    if (debug)
	      {
		if (it.get_R()==b21T::R_Val.front()) 
		  {
		    long i=it.get_i();
		    if ((b21T::disp_driftR[i]<b21T::R_Val[0])&&(fabs(b21T::disp_driftU[i]/(newU-it.get_U()))<100))
		      printf("%3.3ld: r=%e U=%e dU=%e drift: %e (ratio:%.3f)\n",j,it.get_R(),it.get_U(),newU-it.get_U(),b21T::disp_driftU[i],b21T::disp_driftU[i]/(newU-it.get_U()));
		  }
	      }
	  }
	interpolT::destroy(interp);
      }
    
    //printf("data = %e %e\n",data[180*80+10],Mr[10]);
    
#pragma omp parallel for
    for (j=0;j<nSlices;j++)
      { 
	interpT *interp = interpolT::create(interpInit,bT::grid->fullSlicePtr(j));

	if (b21T::kernelType==b21T::kernelTypeV::DELAYED)
	  {
	    //long tmp=0;
	    const kernelCellsItT kit_end = b21T::kernelBoundary[j].end();
	    for (kernelCellsItT kit=b21T::kernelBoundary[j].begin();kit!=kit_end;++kit)
	      {
		dummyCellT &c = kit->second->back();
		if ((c.U<=b21T::U_Val.front())||(c.U>=b21T::U_Val.back())) c.val=0;
		else interpolT::evaluate(interp,c.R,c.U,&c.val);
		//if ((j==0)&&((tmp++)==0)) kit->second->print("\n\nbefore2:","\n\n");
	      }
	  }

	const gridSliceItT it_end=bT::grid->slice_end(j);
	for (gridSliceItT it=bT::grid->slice_begin(j);it!=it_end;++it)
	  {
	    long i=it.get_i();	 
	    if (b21T::disp_driftR[i]<b21T::R_Val.front()) *it=b21T::kernelBoundary[j].find(i)->second->readRoll().val;
	    else if ((b21T::disp_driftR[i]>=b21T::R_Val.back())||
		     (b21T::disp_driftU[i]<=b21T::U_Val.front())||
		     (b21T::disp_driftU[i]>=b21T::U_Val.back())) (*it)=0;
	    else interpolT::evaluate (interp,b21T::disp_driftR[i],b21T::disp_driftU[i],&(*it));
	  }
	/*
	if (b21T::kernelType==b21T::kernelTypeV::DELAYED)
	  {
	    long tmp=0;
	    const kernelCellsItT kit_end = b21T::kernelBoundary[j].end();
	    for (kernelCellsItT kit=b21T::kernelBoundary[j].begin();kit!=kit_end;++kit,++tmp)
	      {
		if ((j==0)&&((tmp++)==0)) kit->second->print("after2:","\n");
	      }
	  }
	*/
	interpolT::destroy(interp);
      }
    //printf("data = %e\n",data[180*80+10]);
    
  }

};

#endif
