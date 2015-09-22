#ifndef __VLASOV_SOLVER21_SPLINE_HXX__
#define __VLASOV_SOLVER21_SPLINE_HXX__

#include "vlasovSolver21.hxx"
#include "quadrature21.hxx"
#include "interpol_spline.hxx"
#include "dimTraits.hxx"

#include "localSpline.hxx"


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
  typedef typename bT::gridItT gridItT;
  typedef typename bT::dirT dirT;
  typedef typename b21T::gridSliceItT gridSliceItT;

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

  std::vector<int> driftR_I;
  std::vector<int> driftU_I;  
  std::vector<double> driftR_F;
  std::vector<double> driftU_F;  

  void decomposeDisplacement()
  {
    std::vector<double> &driftR=b21T::disp_driftR;
    std::vector<double> &driftU=b21T::disp_driftU;
    long i;

    driftR_I.resize(driftR.size());
    driftR_F.resize(driftR.size());
    driftU_I.resize(driftU.size());
    driftU_F.resize(driftU.size());
    
    for (i=0;i<driftR.size();i++)
      {
	std::pair<int,double> r=spline[0].pos2Coord(driftR[i],0);
	driftR_I[i]=r.first;
	driftR_F[i]=r.second;
      }

    for (i=0;i<driftU.size();i++)
      {
	std::pair<int,double> r=spline[0].pos2Coord(driftU[i],1);
	driftU_I[i]=r.first;
	driftU_F[i]=r.second;
      }
  }
  
  struct localSplineDataT
  {
    gridItT it;
    int dim;
    int dir;
    
    localSplineDataT(gridT *g, int dim_, int dir_)
    {
      int dp[P_DIMS];
      int du[U_DIMS];
     
      memset(dp,0,sizeof(int)*P_DIMS);
      memset(du,0,sizeof(int)*U_DIMS);

      dir=dir_;
      dim=dim_;
      int delta=(dir>0)?10:11;
      if (dim<P_DIMS) 
	dp[dim]=delta; 
      else 
	du[dim-P_DIMS]=delta;

      it = g->innerMargin_begin(dp,du,dir);
    }
  };

  std::vector<localSplineDataT> localSplineData;

  void setupLocalSplines()
  {
    long i;

    for (i=0;i<P_DIMS;i++)
      {
	if (!bT::gh->isFullGridBoundary_Plow(i)) 
	  localSplineData.push_back(localSplineDataT(bT::grid,i,-1));
	if (!bT::gh->isFullGridBoundary_Phigh(i))
	  localSplineData.push_back(localSplineDataT(bT::grid,i,+1));
      }
    for (i=0;i<U_DIMS;i++)
      {
	if (!bT::gh->isFullGridBoundary_Ulow(i)) 
	  localSplineData.push_back(localSplineDataT(bT::grid,i+P_DIMS,-1));
	if (!bT::gh->isFullGridBoundary_Uhigh(i))
	  localSplineData.push_back(localSplineDataT(bT::grid,i+P_DIMS,+1));
      }
  }
  /**
  //EINSPLINE
  typedef interpolSpline<dataT> interpolT;
  typedef typename interpolT::interpT interpT;
  typedef typename interpolT::initT interpInitT;
  typedef typename interpolT::boundaryTypeT boundaryTypeT;
  typedef typename interpolT::boundaryTypeV boundaryTypeV; 
  interpInitT interpInit;
  
  
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
    setupLocalSplines();
  }

  /**/
  /**/
  
  typedef localSpline<2,dataT> interpolT;
  typedef typename interpolT::initT interpInitT;
  typedef typename interpolT::boundaryTypeT boundaryTypeT;
  typedef typename interpolT::boundaryTypeV boundaryTypeV; 
  std::vector<interpolT> spline;
  
  void setup()
  {
    interpInitT ini;

    ini.start[0]=b21T::R_V.front();
    ini.start[1]=b21T::U_V.front();
    ini.stop[0]=b21T::R_V.back();
    ini.stop[1]=b21T::U_V.back();
    ini.N[0]=b21T::R_C.size();
    ini.N[1]=b21T::U_C.size();
    ini.scaleType[0]=b21T::Rscale;
    ini.scaleType[1]=b21T::Uscale;
    ini.vlType[0]=bT::grid->valLocation;
    ini.vlType[1]=bT::grid->valLocation;

    ini.bCond[0][0] = boundaryTypeV::BT_NOT_A_KNOT;
    if (bT::gh->isFullGridBoundary_Phigh(0))
      ini.bCond[0][1] = boundaryTypeV::BT_NOT_A_KNOT;
    else
      ini.bCond[0][1] = boundaryTypeV::BT_LOCAL;

    if (bT::gh->isFullGridBoundary_Ulow(0))
      ini.bCond[1][0] = boundaryTypeV::BT_NOT_A_KNOT;
    else
      ini.bCond[1][0] = boundaryTypeV::BT_LOCAL;

    if (bT::gh->isFullGridBoundary_Uhigh(0))
      ini.bCond[1][1] = boundaryTypeV::BT_NOT_A_KNOT;
    else
      ini.bCond[1][1] = boundaryTypeV::BT_LOCAL;
   
    spline.assign(b21T::J_Val.size(),interpolT(ini));
    decomposeDisplacement();
    setupLocalSplines(); 
  }
  
  /**/
protected:

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

  void testSpline1D()
  {
      FILE *ff1=fopen("test_spline1D_X1.txt","w");
      FILE *ff2=fopen("test_spline1D_X2.txt","w");
      FILE *ff3=fopen("test_spline1DALL_X1.txt","w");
      FILE *ff4=fopen("test_spline1DALL_X2.txt","w");
      typedef localSpline<1,double> splineT;
      splineT::initT ini;
      
      long i=0;
      long j,k;
      const long nSlices=bT::grid->nSlices();
      for (i=0;i<1;i++)
	{
	  ini.start[i]=1;
	  ini.stop[i]=5;
	  ini.N[i]=50;
	  ini.scaleType[i]=splineT::scaleTypeV::LINEAR;
	  ini.vlType[i]=splineT::valLocationV::VERTEX;
	  ini.bCond[i][0]=splineT::boundaryTypeV::BT_NOT_A_KNOT;
	  if (i==0) ini.bCond[i][1]=splineT::boundaryTypeV::BT_LOCAL;
	  else ini.bCond[i][1]=splineT::boundaryTypeV::BT_NOT_A_KNOT;
	  //ini.bCond[i][1]=splineT::boundaryTypeV::BT_NOT_A_KNOT;
	}

      splineT::initT ini2(ini);
      ini2.start[0]=ini.stop[0];
      ini2.stop[0]=ini.stop[0] + (ini.stop[0]-ini.start[0]);
      ini2.bCond[0][0]=ini.bCond[0][1];
      ini2.bCond[0][1]=ini.bCond[0][0];

    splineT::initT ini3(ini);
    for (i=0;i<1;i++)
      {
	ini3.bCond[i][0]=ini3.bCond[i][1]=splineT::boundaryTypeV::BT_NOT_A_KNOT;
	ini3.start[i]=ini.start[i];
	ini3.stop[i]=ini2.stop[i];
	if (i==0) ini3.N[i]=ini.N[i]+ini2.N[i];
      } 


      std::vector<double> tabX=scaleT::genScale(ini.start[0],ini.stop[0],ini.N[0],ini.scaleType[0],ini.vlType[0]);
      std::vector<double> tabY=tabX;
      std::vector<double> tabX2=scaleT::genScale(ini.start[0],ini.stop[0],ini.N[0]*88,ini.scaleType[0],ini.vlType[0]);
      std::vector<double> tabY2=tabX2;

      std::vector<double> batX=scaleT::genScale(ini2.start[0],ini2.stop[0],ini2.N[0],ini2.scaleType[0],ini2.vlType[0]);
      std::vector<double> batY=batX;
      std::vector<double> batX2=scaleT::genScale(ini2.start[0],ini2.stop[0],ini2.N[0]*88,ini2.scaleType[0],ini2.vlType[0]);
      std::vector<double> batY2=batX2;

      std::vector<double> allX=scaleT::genScale(ini3.start[0],ini3.stop[0],ini3.N[0],ini3.scaleType[0],ini3.vlType[0]);
      std::vector<double> allY=allX;
      std::vector<double> allX2=scaleT::genScale(ini3.start[0],ini3.stop[0],ini3.N[0]*88,ini3.scaleType[0],ini3.vlType[0]);
      std::vector<double> allY2=allX2;

      for (i=0;i<tabX.size();i++) tabY[i]=cos(tabX[i]*6.59)+sin(tabX[i]*4.84+0.46)*cos(tabX[i]*3.56-0.289);
      for (i=0;i<batX.size();i++) batY[i]=cos(batX[i]*6.59)+sin(batX[i]*4.84+0.46)*cos(batX[i]*3.56-0.289);
      for (i=0;i<allX.size();i++) allY[i]=cos(allX[i]*6.59)+sin(allX[i]*4.84+0.46)*cos(allX[i]*3.56-0.289);
      /*
      for (i=0;i<tabX.size();i++) tabY[i]=exp(-tabX[i]);
      for (i=0;i<batX.size();i++) batY[i]=exp(-batX[i]);
      for (i=0;i<allX.size();i++) allY[i]=exp(-allX[i]);
      */
      splineT spline(ini,tabY.begin(),tabY.end());    
      splineT spline2(ini2,batY.begin(),batY.end());    
      splineT spline3(ini3,allY.begin(),allY.end());
      
      typename splineT::localBoundaryT lb=spline.getLocalBoundaries();
      typename splineT::localBoundaryT lb2=spline2.getLocalBoundaries();

      spline.setLocalBoundaries(lb2);
      spline2.setLocalBoundaries(lb);
      
      spline.build();
      spline2.build();
      spline3.build(); 

      for (i=0;i<tabX.size();i++) fprintf(ff1,"%ld %e %e %e %e\n",i,tabX[i],tabY[i],batX[i],batY[i]);
      for (i=0;i<tabX2.size();i++) 
	{
	  tabY2[i]=spline.evaluate(&tabX2[i]);
	  batY2[i]=spline2.evaluate(&batX2[i]);
	  fprintf(ff2,"%ld %e %e %e %e\n",i,tabX2[i],tabY2[i],batX2[i],batY2[i]);
	}
      for (i=0;i<allX.size();i++) fprintf(ff3,"%ld %e %e \n",i,allX[i],allY[i]);
      for (i=0;i<allX2.size();i++) 
	{
	  allY2[i]=spline3.evaluate(&allX2[i]);
	  fprintf(ff4,"%ld %e %e\n",i,allX2[i],allY2[i]);
	}
      
      fclose(ff1);
      fclose(ff2);     
      fclose(ff3);     
      fclose(ff4);     
    }

  void testSpline2D()
  {
    FILE *ff1=fopen("test_spline2D_X1.txt","w");
    FILE *ff2=fopen("test_spline2D_X2.txt","w");
    FILE *ff3=fopen("test_spline2DALL_X1.txt","w");
    FILE *ff4=fopen("test_spline2DALL_X2.txt","w");

    typedef localSpline<2,double> splineT;
    splineT::initT ini;
    
    long i=0;
    long j,k;
    const long nSlices=bT::grid->nSlices();
    for (i=0;i<2;i++)
      {
	ini.start[i]=0;
	ini.stop[i]=3.14159*2;
	ini.N[i]=20;
	ini.scaleType[i]=splineT::scaleTypeV::LINEAR;
	ini.vlType[i]=splineT::valLocationV::VERTEX;
	ini.bCond[i][0]=splineT::boundaryTypeV::BT_NOT_A_KNOT;	
	if (i==0) ini.bCond[i][1]=splineT::boundaryTypeV::BT_LOCAL;
	else ini.bCond[i][1]=splineT::boundaryTypeV::BT_NOT_A_KNOT;
	//ini.bCond[i][1]=splineT::boundaryTypeV::BT_NOT_A_KNOT;
      }

    splineT::initT ini2(ini);
    ini2.start[0]=ini.stop[0];
    ini2.stop[0]=ini.stop[0] + (ini.stop[0]-ini.start[0]);
    ini2.bCond[0][0]=ini.bCond[0][1];
    ini2.bCond[0][1]=ini.bCond[0][0];

    splineT::initT ini3(ini);
    for (i=0;i<2;i++)
      {
	ini3.bCond[i][0]=ini3.bCond[i][1]=splineT::boundaryTypeV::BT_NOT_A_KNOT;
	ini3.start[i]=ini.start[i];
	ini3.stop[i]=ini2.stop[i];
	if (i==0) ini3.N[i]=ini.N[i]+ini2.N[i];
      } 

    std::vector<double> tabX=scaleT::genScale(ini.start[0],ini.stop[0],ini.N[0],ini.scaleType[0],ini.vlType[0]);
    std::vector<double> tabY=scaleT::genScale(ini.start[1],ini.stop[1],ini.N[1],ini.scaleType[1],ini.vlType[1]);
    std::vector<double> tabZ(tabX.size()*tabY.size());
    std::vector<double> tabX2=scaleT::genScale(ini.start[0],ini.stop[0],ini.N[0]*88,ini.scaleType[0],ini.vlType[0]);
    std::vector<double> tabY2=scaleT::genScale(ini.start[1],ini.stop[1],ini.N[1]*88,ini.scaleType[1],ini.vlType[0]);
    std::vector<double> tabZ2(tabX2.size()*tabY2.size());

    std::vector<double> batX=scaleT::genScale(ini2.start[0],ini2.stop[0],ini2.N[0],ini2.scaleType[0],ini2.vlType[0]);
    std::vector<double> batY=scaleT::genScale(ini2.start[1],ini2.stop[1],ini2.N[1],ini2.scaleType[1],ini2.vlType[1]);
    std::vector<double> batZ(batX.size()*batY.size());
    std::vector<double> batX2=scaleT::genScale(ini2.start[0],ini2.stop[0],ini2.N[0]*88,ini2.scaleType[0],ini2.vlType[0]);
    std::vector<double> batY2=scaleT::genScale(ini2.start[1],ini2.stop[1],ini2.N[1]*88,ini2.scaleType[1],ini2.vlType[0]);
    std::vector<double> batZ2(batX2.size()*batY2.size());

    std::vector<double> allX=scaleT::genScale(ini3.start[0],ini3.stop[0],ini3.N[0],ini3.scaleType[0],ini3.vlType[0]);
    std::vector<double> allY=scaleT::genScale(ini3.start[1],ini3.stop[1],ini3.N[1],ini3.scaleType[1],ini3.vlType[1]);
    std::vector<double> allZ(allX.size()*allY.size());
    std::vector<double> allX2=scaleT::genScale(ini3.start[0],ini3.stop[0],ini3.N[0]*88,ini3.scaleType[0],ini3.vlType[0]);
    std::vector<double> allY2=scaleT::genScale(ini3.start[1],ini3.stop[1],ini3.N[1]*88,ini3.scaleType[1],ini3.vlType[0]);
    std::vector<double> allZ2(allX2.size()*allY2.size());

    k=0;
    for (i=0;i<tabY.size();i++)
      for (j=0;j<tabX.size();j++,k++)
	{
	  tabZ[k]=cos(tabX[j])+sin(tabX[j]*0.84+tabY[i]*0.46-0.13)*cos(tabX[j]*0.56-tabY[i]);
	  batZ[k]=cos(batX[j])+sin(batX[j]*0.84+batY[i]*0.46-0.13)*cos(batX[j]*0.56-batY[i]);	  

	  //tabZ[k]=cos(tabX[j])+cos(tabX[j]-tabY[i]);	  
	  //batZ[k]=cos(batX[j])+cos(batX[j]-batY[i]);	  

	}
    
    k=0;
    for (i=0;i<allY.size();i++)
      for (j=0;j<allX.size();j++,k++)
	{
	  allZ[k]=cos(allX[j])+sin(allX[j]*0.84+allY[i]*0.46-0.13)*cos(allX[j]*0.56-allY[i]);
	  //allZ[k]=cos(allX[j])+cos(allX[j]-allY[i]);
	}
    

    splineT spline(ini,tabZ.begin(),tabZ.end());    
    splineT spline2(ini2,batZ.begin(),batZ.end());    
    splineT spline3(ini3,allZ.begin(),allZ.end());
  
    typename splineT::localBoundaryT lb=spline.getLocalBoundaries();
    typename splineT::localBoundaryT lb2=spline2.getLocalBoundaries();

    spline.setLocalBoundaries(lb2);
    spline2.setLocalBoundaries(lb);
    
    spline.build();
    spline2.build();
    spline3.build();
    
    k=0;
    for (i=0;i<tabY.size();i++)
      for (j=0;j<tabX.size();j++,k++)
	{
	  fprintf(ff1,"%ld %e %e %e %e %e %e\n",k,tabX[j],tabY[i],tabZ[k],batX[j],batY[i],batZ[k]);
	}
    k=0;
    for (i=0;i<tabY2.size();i++)
      for (j=0;j<tabX2.size();j++,k++)
	{
	  tabZ2[k]=spline.evaluate(tabX2[j],tabY2[i]);
	  batZ2[k]=spline2.evaluate(batX2[j],batY2[i]);
	  fprintf(ff2,"%ld %e %e %e %e %e %e\n",k,tabX2[j],tabY2[i],tabZ2[k],batX2[j],batY2[i],batZ2[k]);
	}
    
    k=0;
    for (i=0;i<allY.size();i++)
      for (j=0;j<allX.size();j++,k++)
	{
	  fprintf(ff3,"%ld %e %e %e\n",k,allX[j],allY[i],allZ[k]);
	}
    k=0;
    for (i=0;i<allY2.size();i++)
      for (j=0;j<allX2.size();j++,k++)
	{
	  allZ2[k]=spline3.evaluate(allX2[j],allY2[i]);
	  fprintf(ff4,"%ld %e %e %e\n",k,allX2[j],allY2[i],allZ2[k]);
	}
    
    fclose(ff1);
    fclose(ff2);
    fclose(ff3);
    fclose(ff4);
  }
  
  /**/
  void makeOneStep(double t)
  {  
    
    // testSpline1D();
    // testSpline2D();
    // exit(-1);

    long r,u,j,i;
    double dth=0.5*bT::deltaTime;
    long delta = bT::grid->getCellStride(2);
    const long nSlices=bT::grid->nSlices();
    double *driftR=&b21T::disp_driftR[0];
    double *driftU=&b21T::disp_driftU[0];

    if (bT::noLinking)
      {
	bT::gh->synchronize();
#pragma omp parallel for 
	for (j=0;j<nSlices;j++)
	  { 	
	    spline[j].assign(bT::grid->fullSlicePtr(j));
	    spline[j].build();
	    if (b21T::kernelType==b21T::kernelTypeV::DELAYED)
	      {
		const kernelCellsItT kit_end = b21T::kernel.end(j);
		for (kernelCellsItT kit=b21T::kernel.begin(j);kit!=kit_end;++kit)
		  {
		    dummyCellT &c = kit->second->back();
		    if ((c.U<=b21T::U_Val.front())||(c.U>=b21T::U_Val.back())) c.val=0;
		    else c.val=spline[j].evaluate(c.R,c.U);
		  }
	      }
	  
	    const gridSliceItT it_end=bT::grid->slice_end(j);
	    for (gridSliceItT it=bT::grid->slice_begin(j);it!=it_end;++it)
	      {
		long i=it.get_i();
		if (driftR[i]<b21T::kernelRadius) *it=b21T::kernel.readRoll(j,i).val;
		else if ((driftR[i]>=b21T::R_Val.back())||
			 (driftU[i]<=b21T::U_Val.front())||
			 (driftU[i]>=b21T::U_Val.back())) (*it)=0;
		else (*it)=spline[j].evaluate(driftR[i],driftU[i]);
	      }
	  }	
      }

    bT::gh->synchronize();
    
    std::vector<double> Mr = quadratureT::M_of_less_than_R(bT::gh);
    for (j=0;j<Mr.size();j++) Mr[j]+=b21T::kernelMass;
    const std::vector<double> &R = bT::grid->getValCoord21_R();

    //printf("(%d) M = %e\n",bT::mpiCom->rank(),Mr.back());
    for (i=0;i<Mr.size();i++) Mr[i]*=bT::sp.G/(R[i]*R[i])*bT::deltaTime;
        
#pragma omp parallel for
    for (j=0;j<nSlices;j++)
      {
	const gridSliceItT it_begin=bT::grid->slice_begin(j);
	const gridSliceItT it_end=bT::grid->slice_end(j);

	spline[j].assign(bT::grid->fullSlicePtr(j));
	spline[j].build();
	for (gridSliceItT it=it_begin;it!=it_end;++it)
	  {
	    double newU=it.get_U()+Mr[it.get_r()];
	    if ((newU<=b21T::U_Val.front())||(newU>=b21T::U_Val.back())) (*it)=0;
	    else (*it)=spline[j].evaluateNxPy(it.get_r(),newU);
	  }
      }

    bT::gh->synchronize();
   
#pragma omp parallel for
    for (j=0;j<nSlices;j++)
      { 
	spline[j].assign(bT::grid->fullSlicePtr(j));
	spline[j].build();
	if (b21T::kernelType==b21T::kernelTypeV::DELAYED)
	  {
	    const kernelCellsItT kit_end = b21T::kernel.end(j);
	    for (kernelCellsItT kit=b21T::kernel.begin(j);kit!=kit_end;++kit)
	      {
		dummyCellT &c = kit->second->back();
		if ((c.U<=b21T::U_Val.front())||(c.U>=b21T::U_Val.back())) c.val=0;
		else c.val=spline[j].evaluate(c.R,c.U);
	      }
	  }

	const gridSliceItT it_end=bT::grid->slice_end(j);
	for (gridSliceItT it=bT::grid->slice_begin(j);it!=it_end;++it)
	  {
	    long i=it.get_i();	 
	    if (driftR[i]<b21T::kernelRadius) *it=b21T::kernel.readRoll(j,i).val;
	    else if ((driftR[i]>=b21T::R_Val.back())||
		     (driftU[i]<=b21T::U_Val.front())||
		     (driftU[i]>=b21T::U_Val.back())) (*it)=0;
	    else (*it)=spline[j].evaluateC(driftR_I[i],driftR_F[i],driftU_I[i],driftU_F[i]);
	  }
      }
       
  }
  /**/
  /**
  // EINSPLINE
  void makeOneStep(double t)
  {  
   

    long r,u,j,i;
    double dth=0.5*bT::deltaTime;
    long delta = bT::grid->getCellStride(2);
    const long nSlices=bT::grid->nSlices();
    double *driftR=&b21T::disp_driftR[0];
    double *driftU=&b21T::disp_driftU[0];

    // for (j=0;j<nSlices;j++)
    //   { 
    // 	interpT *interp = interpolT::create(interpInit,bT::grid->fullSlicePtr(j));
    // 	const gridSliceItT it_end=bT::grid->slice_end(j);
    // 	for (gridSliceItT it=bT::grid->slice_begin(j);it!=it_end;++it)
    // 	  {
    // 	    long i=it.get_i();
    // 	    double R=it.get_R();
    // 	    double U=it.get_U();
    // 	    if (U<=b21T::U_Val.front()) continue;
    // 	    else if (U>=b21T::U_Val.back()) continue;
    // 	    //else if (R<=b21T::R_Val.front()) continue;
    // 	    else if (R>=b21T::R_Val.back()) continue;
    // 	    else interpolT::evaluate(interp,it.get_R(),it.get_U(),&(*it));
	   
    // 	  }
    // 	  interpolT::destroy(interp);
    //   }
    //   return;
    
    if (bT::noLinking)
      {
	bT::gh->synchronize();
#pragma omp parallel for 
	for (j=0;j<nSlices;j++)
	  { 		
	    interpT *interp = interpolT::create(interpInit,bT::grid->fullSlicePtr(j));
	    //printf("HI\n");
	    if (b21T::kernelType==b21T::kernelTypeV::DELAYED)
	      {
		const kernelCellsItT kit_end = b21T::kernel.end(j);
		for (kernelCellsItT kit=b21T::kernel.begin(j);kit!=kit_end;++kit)
		  {
		    dummyCellT &c = kit->second->back();
		    if ((c.U<=b21T::U_Val.front())||(c.U>=b21T::U_Val.back())) c.val=0;
		    else interpolT::evaluate(interp,c.R,c.U,&c.val);
		  }
	      }
	    //printf("HI2\n");
	    const gridSliceItT it_end=bT::grid->slice_end(j);
	    for (gridSliceItT it=bT::grid->slice_begin(j);it!=it_end;++it)
	      {
		long i=it.get_i();
		if (driftR[i]<b21T::kernelRadius) *it=b21T::kernel.readRoll(j,i).val;
		else if ((driftR[i]>=b21T::R_Val.back())||
			 (driftU[i]<=b21T::U_Val.front())||
			 (driftU[i]>=b21T::U_Val.back())) (*it)=0;
		else interpolT::evaluate(interp,driftR[i],driftU[i],&(*it));
	      }
	    //printf("BI\n");
	    interpolT::destroy(interp);
	  }	
      }

    bT::gh->synchronize();

    std::vector<double> Mr = quadratureT::M_of_less_than_R(bT::gh);
    for (j=0;j<Mr.size();j++) Mr[j]+=b21T::kernelMass;
    const std::vector<double> &R = bT::grid->getValCoord21_R();

    //printf("(%d) M = %e\n",bT::mpiCom->rank(),Mr.back());
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
	    
	    // if (debug)
	    //   {
	    // 	if (it.get_R()==b21T::R_Val.front()) 
	    // 	  {
	    // 	    long i=it.get_i();
	    // 	    if ((driftR[i]<b21T::kernelRadius)&&(fabs(driftU[i]/(newU-it.get_U()))<100))
	    // 	      printf("%3.3ld: r=%e U=%e dU=%e drift: %e (ratio:%.3f)\n",j,it.get_R(),it.get_U(),newU-it.get_U(),driftU[i],driftU[i]/(newU-it.get_U()));
	    // 	  }
	    //   }

	  }
	interpolT::destroy(interp);
      }

    
    bT::gh->synchronize();
    
#pragma omp parallel for
    for (j=0;j<nSlices;j++)
      { 
	interpT *interp = interpolT::create(interpInit,bT::grid->fullSlicePtr(j));

	if (b21T::kernelType==b21T::kernelTypeV::DELAYED)
	  {
	    const kernelCellsItT kit_end = b21T::kernel.end(j);
	    for (kernelCellsItT kit=b21T::kernel.begin(j);kit!=kit_end;++kit)
	      {
		dummyCellT &c = kit->second->back();
		if ((c.U<=b21T::U_Val.front())||(c.U>=b21T::U_Val.back())) c.val=0;
		else interpolT::evaluate(interp,c.R,c.U,&c.val);
	      }
	  }

	const gridSliceItT it_end=bT::grid->slice_end(j);
	for (gridSliceItT it=bT::grid->slice_begin(j);it!=it_end;++it)
	  {
	    long i=it.get_i();	 
	    if (driftR[i]<b21T::kernelRadius) *it=b21T::kernel.readRoll(j,i).val;
	    else if ((driftR[i]>=b21T::R_Val.back())||
		     (driftU[i]<=b21T::U_Val.front())||
		     (driftU[i]>=b21T::U_Val.back())) (*it)=0;
	    else interpolT::evaluate(interp,driftR[i],driftU[i],&(*it));
	  }

	interpolT::destroy(interp);
      }
       
  }
  /***/
};

#endif
