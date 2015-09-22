#ifndef __LOCAL_SPLINE_HXX__
#define __LOCAL_SPLINE_HXX__

#include <algorithm>
#include "helpers.hxx"
#include "scale.hxx"

template <
  int NDIMS,
  typename dataType = double
  >
class localSpline
{
public:
  typedef localSpline<NDIMS,dataType> myT;
  typedef dataType dataT;
  typedef scale<dataT> scaleT;
  
  typedef typename scaleT::scaleTypeT scaleTypeT; 
  typedef typename scaleT::scaleTypeV scaleTypeV;
  typedef typename scaleT::valLocationT valLocationT;
  typedef typename scaleT::valLocationV valLocationV;  

  // Boundary conditions:
  // PERIODIC: 
  // FLAT: f'(x0)=0
  // NATURAL: f''(x0)=0
  // SUPRA_NATURAL: f'''(x0)=0
  // NOT_A_KNOT: f''''(x0)=0
  // LOCAL: f'(x0) is computed from the spline interpolations over two adjacent domains  
  //        This is used for MPI parallelisation

  struct boundaryTypeV {
    enum type {BT_PERIODIC=0, BT_FLAT=1, BT_NATURAL=2, BT_SUPRA_NATURAL=3, BT_NOT_A_KNOT=4, BT_LOCAL=5};
  };
  typedef typename boundaryTypeV::type boundaryTypeT;

  struct initT {
    double start[NDIMS];
    double stop[NDIMS];
    int N[NDIMS];
    scaleTypeT scaleType[NDIMS];
    valLocationT vlType[NDIMS];

    boundaryTypeT bCond[NDIMS][2];
  };

  struct localBoundaryT {
    std::vector<dataT> B[NDIMS][2];
    void clear() 
    {
      for (int i=0;i<NDIMS;i++)
	{
	  B[i][0].clear();
	  B[i][1].clear();
	}
    }
  };
  
  void init(const initT &ini)
  {
    long i,i0;
    long maxS=0;
    
    purge();
    nVal=1;
    stride[0]=1;
    for (i=0;i<NDIMS;i++) 
      {
	sclType[i]=ini.scaleType[i];
	bCond[i][0]=ini.bCond[i][0];
	bCond[i][1]=ini.bCond[i][1];
	std::vector<double> scl=scaleT::genScale(ini.start[i],ini.stop[i],ini.N[i],sclType[i],ini.vlType[i]);	
	dim[i]=scl.size();
	stride[i+1]=stride[i]*(dim[i]+2);
	nVal*=dim[i];

	start[i]=scl.front();
	stop[i]=scl.back();

	scaleT::rescale(start[i],sclType[i]);
	scaleT::rescale(stop[i],sclType[i]);
	step[i]=(stop[i]-start[i])/(dim[i]-1);
	step_inv[i]=1./step[i];
	/*
	m3h6[i]=-1./(step[i]*2.);
	p3h6[i]=1./(step[i]*2.);
	m3h[i]=-3./(step[i]);
	p3h[i]=3./(step[i]);
	mh3[i]=-step[i]/3.;
	*/
	if (dim[i]>maxS) maxS=dim[i];
      } 
    //printf("start= (%g %g) = (%g %g %ld %ld %ld)\n",start[0],start[1],ini.start[0],ini.start[1],(long)ini.N[0],(long)sclType[0],(long)ini.vlType[0]);
    for (i=0;i<NDIMS;i++)
      {
	bStride[i][0]=1;
	
	long k=0;
	for (long j=0;j<NDIMS;j++)
	  {
	    if (i!=j) 
	      {
		bStride[i][k+1]=bStride[i][k]*(dim[j]+2);
		bDim[i][k]=dim[j];
		k++;
	      }
	    else 
	      {
		bStride[i][k]*=(dim[j]+2);
	      }
	  }
      }
    purge();
    computeLocalCoefficients();
    /*
    D.clear();D.reserve(maxS);
    L.clear();L.reserve(maxS-1);
    
    D.push_back(4);
    D.push_back(3.5);
    L.push_back(0.25);
    
    // D.push_back(6);
    // D.push_back(4);
    // L.push_back(1.0/6.0);
    
    i0=D.size();

    for (i=i0;i<maxS;i++)
      {
	L.push_back(1./D.back());
	D.push_back(4.-L.back());
      }
    
    for (i=0;i<NDIMS;i++) 
      {
	L_back[i]=1./(D[D.size()-1]*D[D.size()-2]);
	D_back[i]=1.-L_back[i];
      }
      
    // for (i=0;i<NDIMS;i++) 
    //   {
    // 	L_back[i]=-3.0/(D[D.size()-1]*D[D.size()-2]);
    // 	D_back[i]=1.0-L_back[i];
    //   }
      
    for (i=0;i<D.size();i++) D[i]=1./D[i];
*/
  }

  void init(const initT &ini, const dataT *dataP)
  {
    init(ini);
    assign(dataP);
  }

  template <class inputIterator>
  void init(const initT &ini, const inputIterator &start,const inputIterator &stop)
  {
    init(ini);
    assign(start,stop);
  }
  
  localSpline(const initT &ini):data(NULL)
  {
    init(ini);
  }

  localSpline(const initT &ini, const dataT *dataP):data(NULL)
  {
    init(ini,dataP);
  }

  template <class inputIterator>
  localSpline(const initT &ini, const inputIterator &start,const inputIterator &stop):data(NULL)
  {
    init(ini,start,stop);
  }
  /*
  localSpline(const myT& b)
  {
    (*this)=b;
    }

  myT & operator=(const myT& b)
  {
    memcpy(dim,b.dim,sizeof(int)*NDIMS);
    memcpy(stride,b.stride,sizeof(long)*(NDIMS+1));
    memcpy(bStride,b.bStride,sizeof(long)*NDIMS*NDIMS);
    memcpy(bDim,b.bDim,sizeof(int)*NDIMS*NDIMS);
    memcpy(start,b.start,sizeof(double)*NDIMS);
    memcpy(stop,b.stop,sizeof(double)*NDIMS);
    memcpy(step,b.step,sizeof(double)*NDIMS);
    memcpy(sclType,b.sclType,sizeof(scaleTypeT)*NDIMS);
    memcpy(bCond,b.bCond,sizeof(boundaryTypeT)*NDIMS*2);
    data=b.data;
  }
  */
  ~localSpline()
  {
    if (data!=NULL) free(data);
  }

  void assign(const dataT *dataP)
  {
    assign(dataP,dataP+nVal);
  }

  template <class inputIterator>
  void assign(const inputIterator &start,const inputIterator &stop)
  {
    int w[NDIMS];
    long i,k;
    long strd[NDIMS+1];

    for (i=0;i<=NDIMS;i++) strd[i]=2*stride[i];
    
    //data.assign(stride[NDIMS],0);
    if (data==NULL)
      {
	int res=posix_memalign((void**)&data,32,stride[NDIMS]*sizeof(dataT));
	memset(data,0,stride[NDIMS]*sizeof(dataT));
	//data=(dataT*)calloc(stride[NDIMS],sizeof(dataT));
	//printf("(long)data = %p(%ld) -> %ld\n",data,(long)data,((long)data)%32);
      }
    for (i=0,k=0;i<NDIMS;i++) 
      {
	w[i]=0;
	k+=stride[i];
      }
    
    for (inputIterator it=start;it!=stop;++it)
      {
	data[k]=*it;
	
	++w[0];++k;
	if (w[0]==dim[0]) 
	  {	   
	    k+=strd[0];
	    hlp::getNext<NDIMS>(w,k,dim,strd);
	  }
      }
   
    computeLocalBoundaries();
  }
  
  void build()
  {
    // check if boundary was assigned !!!!

    switch (NDIMS)
      {
      case 1: create_1D(); break;
      case 2: create_2D(); break;
      case 3: create_3D(); break;
      default: create_ND();
      }
  }

  long boundarySize(int dm)
  {
    return stride[NDIMS]/(dim[dm]+2);
  }
  
  void getLocalBoundaries(localBoundaryT &lb)
  {
    int dm;
    long i;

    for (int dm=0;dm<NDIMS;dm++)
      {
	for (int dr=0;dr<=1;dr++)
	  {
	    if (bCond[dm][dr]==boundaryTypeV::BT_LOCAL)
	      {
		std::vector<dataT> &dest=lb.B[dm][dr];
	
		int w[NDIMS]= {0};
		const long bN=stride[NDIMS]/(dim[dm]+2);
		long k=(dr>0)?((dim[dm]+1)*stride[dm]):0;

		dest.resize(bN);
		for (i=0;i<bN;i++)
		  {
		    dest[i]=data[k];
		    w[0]++;k+=bStride[dm][0];
		    if (w[0]==bDim[dm][0]) hlp::getNext<NDIMS-1>(w,k,bDim[dm],bStride[dm]);	
		  }
	      }
	    else lb.B[dm][dr].clear();
	  }
      }
  }
  
  localBoundaryT getLocalBoundaries()
  {
    localBoundaryT lb;
    getLocalBoundaries(lb);
    return lb;
  }

  bool setLocalBoundaries(const localBoundaryT &lb,const bool reverse=false)
  {
    long i;

    for (int dm=0;dm<NDIMS;dm++)
      {
	for (int dr=0;dr<=1;dr++)
	  {
	    if (bCond[dm][dr]==boundaryTypeV::BT_LOCAL)
	      {
		const std::vector<dataT> &src=lb.B[dm][(reverse)?(1-dr):dr];
		
		int w[NDIMS]= {0};
		const long bN=stride[NDIMS]/(dim[dm]+2);
		long k=(dr>0)?((dim[dm]+1)*stride[dm]):0;
		//assert(src.size()==bN);

		if (src.size()!=bN) 
		  {
		    fprintf(stderr,"ERROR in setLocalBoundaries: boundary[%d][%d] has invalid size (%ld!=%ld)!\n",dm,dr,src.size(),bN);
		    exit(-1);
		  }
												
		for (i=0;i<bN;i++)
		  {
		    data[k]+=src[i];
		    w[0]++;k+=bStride[dm][0];
		    if (w[0]==bDim[dm][0]) hlp::getNext<NDIMS-1>(w,k,bDim[dm],bStride[dm]);	
		  }
	      }
	  }
      }
    return true;
  }
  
  dataT evaluate(const double pos[NDIMS])
  {
    switch (NDIMS)
      {
      case 1: return evaluate(pos[0]);
      case 2: return evaluate(pos[0],pos[1]);
      case 3: return evaluate(pos[0],pos[1],pos[2]);
      default: 
	fprintf(stderr,"ERROR: spline::evaluate not implemented in %dD!\n",NDIMS);
	exit(-1);
      }
    return 0; //dummy
  }

  std::pair<int,double> pos2Coord(double p, int dir)
  {
    scaleT::rescale(p,sclType[dir]);
    double x=(p-start[dir])*step_inv[dir];
    x=modf(x,&p);
    return std::make_pair((int)p,x);
  }

  dataT evaluateC(int *N, double *x)
  {
    switch (NDIMS)
      {
      case 1: return evaluateC(*N,*x);
      case 2: return evaluateC(N[0],x[0],N[1],x[1]);
      case 3: return evaluateC(N[0],x[0],N[1],x[1],N[2],x[2]);
      
      default: 
	fprintf(stderr,"ERROR: spline::evaluateC not implemented in %dD!\n",NDIMS);
	exit(-1);
      }
    return 0; //dummy
  }

  dataT evaluateC(int N, double x)
  {
    return (1./6.)*(data[N]*BV2(1.-x)+data[N+1]*BV1(1.-x)+data[N+2]*BV1(x)+data[N+3]*BV2(x));
  }

  // x/y are fractional parts, Nx and Ny coords
  dataT evaluateC(int Nx, double x, int Ny, double y)
  {
    Nx++;Ny++;
    const long ref[4]={Nx+(Ny-1)*stride[1],Nx+Ny*stride[1],Nx+(Ny+1)*stride[1],Nx+(Ny+2)*stride[1]};    
    const double b2xI=BV2(1.-x);
    const double b1xI=BV1(1.-x);
    const double b2x=BV2(x);
    const double b1x=BV1(x);
    double result=1.0/36.0;
    
    result *= 
      BV2(1.-y)*(data[ref[0]-1]*b2xI+data[ref[0]]*b1xI+data[ref[0]+1]*b1x+data[ref[0]+2]*b2x) +
      BV1(1.-y)*(data[ref[1]-1]*b2xI+data[ref[1]]*b1xI+data[ref[1]+1]*b1x+data[ref[1]+2]*b2x) +
      BV1(y)*(data[ref[2]-1]*b2xI+data[ref[2]]*b1xI+data[ref[2]+1]*b1x+data[ref[2]+2]*b2x) +
      BV2(y)*(data[ref[3]-1]*b2xI+data[ref[3]]*b1xI+data[ref[3]+1]*b1x+data[ref[3]+2]*b2x);

    return result;
  }

  // x/y are fractional parts, Nx and Ny coords
  dataT evaluateC(int Nx, double x, int Ny, double y, int Nz, double z)
  {
    printf("spline::evaluateC not implemented in 3D!\n");exit(0);
    /*
    Nx++;Ny++;Nz++;
    const long ref[4]={Nx+(Ny-1)*stride[1],Nx+Ny*stride[1],Nx+(Ny+1)*stride[1],Nx+(Ny+2)*stride[1]};    
    const double b2xI=BV2(1.-x);
    const double b1xI=BV1(1.-x);
    const double b2x=BV2(x);
    const double b1x=BV1(x);
    double result=1.0/36.0;
    
    result *= 
      BV2(1.-y)*(data[ref[0]-1]*b2xI+data[ref[0]]*b1xI+data[ref[0]+1]*b1x+data[ref[0]+2]*b2x) +
      BV1(1.-y)*(data[ref[1]-1]*b2xI+data[ref[1]]*b1xI+data[ref[1]+1]*b1x+data[ref[1]+2]*b2x) +
      BV1(y)*(data[ref[2]-1]*b2xI+data[ref[2]]*b1xI+data[ref[2]+1]*b1x+data[ref[2]+2]*b2x) +
      BV2(y)*(data[ref[3]-1]*b2xI+data[ref[3]]*b1xI+data[ref[3]+1]*b1x+data[ref[3]+2]*b2x);

    return result;
    */
  }

  dataT evaluateNxPy(int Nx,double yy)
  {
    Nx++;
    scaleT::rescale(yy,sclType[1]);
    double y=(yy-start[1])*step_inv[1];
    y=modf(y,&yy);
    int Ny=1+(int)yy;
    
    const long ref[4]={Nx+(Ny-1)*stride[1],Nx+Ny*stride[1],Nx+(Ny+1)*stride[1],Nx+(Ny+2)*stride[1]};
    double result=1.0/36.0;

    result *= 
      BV2(1.-y)*(data[ref[0]-1]+data[ref[0]]*4+data[ref[0]+1]) +
      BV1(1.-y)*(data[ref[1]-1]+data[ref[1]]*4+data[ref[1]+1]) +
      BV1(y)*(data[ref[2]-1]+data[ref[2]]*4+data[ref[2]+1]) +
      BV2(y)*(data[ref[3]-1]+data[ref[3]]*4+data[ref[3]+1]);

    return result;
  }

  dataT evaluatePxNy(double xx,int Ny)
  {
    Ny++;
    scaleT::rescale(xx,sclType[0]);
    double x=(xx-start[0])*step_inv[0];
    x=modf(x,&xx);
    int Nx=1+(int)xx;
    
    const long ref[4]={Nx+(Ny-1)*stride[1],Nx+Ny*stride[1],Nx+(Ny+1)*stride[1],Nx+(Ny+2)*stride[1]};
    const double b2xI=BV2(1.-x);
    const double b1xI=BV1(1.-x);
    const double b2x=BV2(x);
    const double b1x=BV1(x);
    double result=1.0/36.0;

    result *= 
      (data[ref[0]-1]*b2xI+data[ref[0]]*b1xI+data[ref[0]+1]*b1x+data[ref[0]+2]*b2x) +
      4*(data[ref[1]-1]*b2xI+data[ref[1]]*b1xI+data[ref[1]+1]*b1x+data[ref[1]+2]*b2x) +
      (data[ref[2]-1]*b2xI+data[ref[2]]*b1xI+data[ref[2]+1]*b1x+data[ref[2]+2]*b2x);
      
    return result;
  }

  dataT evaluate(double xx)
  {
    scaleT::rescale(xx,sclType[0]);

    double x=(xx-start[0])*step_inv[0];
    double result=1./6.;

    x=modf(x,&xx);
    return evaluateC((int)xx,x);
  }

  dataT evaluate(double xx,double yy)
  {    
    scaleT::rescale(xx,sclType[0]);
    scaleT::rescale(yy,sclType[1]);
    
    double x=(xx-start[0])*step_inv[0];
    double y=(yy-start[1])*step_inv[1];
    x=modf(x,&xx);
    y=modf(y,&yy);
    
    return evaluateC((int)xx,x,(int)yy,y);
  }
 
  dataT evaluate(double xx,double yy,double zz)
  {
    scaleT::rescale(xx,sclType[0]);
    scaleT::rescale(yy,sclType[1]);
    scaleT::rescale(zz,sclType[2]);
    
    double x=(xx-start[0])*step_inv[0];
    double y=(yy-start[1])*step_inv[1];
    double z=(zz-start[2])*step_inv[2];
    x=modf(x,&xx);
    y=modf(y,&yy);
    z=modf(z,&zz);
    return evaluateC((int)xx,x,(int)yy,y,(int)zz,z);
  }

  void purge()
  {
    //data.clear();
    if (data!=NULL) free(data);
    int k;
    for (k=0;k<NDIMS;k++)
      {
	lookup[k].clear();
      }
  }

private:
  // double m3h6[NDIMS],p3h6[NDIMS];
  // double m3h[NDIMS],p3h[NDIMS];
  // double mh3[NDIMS];
  // static const double A[3];

  //std::vector<double> D;
  //std::vector<double> L;
  //double D_back[NDIMS];
  //double L_back[NDIMS];

  static const int WS3;
  static const double WM3[];

  long nVal;
  
  int dim[NDIMS];
  long stride[NDIMS+1];
  long bStride[NDIMS][NDIMS];  
  int bDim[NDIMS][NDIMS];  
  double start[NDIMS];
  double stop[NDIMS];
  double step[NDIMS];
  double step_inv[NDIMS];
  scaleTypeT sclType[NDIMS];
  boundaryTypeT bCond[NDIMS][2];
  std::vector<double> lookup[NDIMS];
  //std::vector<dataT> data;
  dataT *data;
  //std::vector<double> tmpArr;
  dataT *tmpArr;
  std::vector<double> localCoef;

  double BV2(double x) {return x*x*x;}
  double BV1(double x) {return 1.0+3.0*x+3.0*x*x-3.0*x*x*x;}

 
  void computeLocalCoefficients(int order=3)
  {   
    if (order<1) order=1;
    const int Nmax=(1<<order)+2;
    
    if (Nmax==localCoef.size()) return;
    /*
    else if (Nmax==WS3)
      {
	localCoef.assign(WM3,WM3+WS3);
	return;
      }
    */
    localCoef.assign(Nmax,0);
    localCoef[0]=6./7.;
    localCoef[1]=-3./14.;    
    std::vector<double> tmp(localCoef);
    localCoef[2]=1./14.;

    int Ncur=2;
    double alpha=1.0;
    long i,j;
    for (i=1;i<order;i++)
      {   
	for (j=0;j<Ncur-1;j++) tmp[j]-=localCoef[Ncur]*localCoef[Ncur-j-2]/alpha;
	for (j=0;j<Ncur;j++) tmp[j+Ncur]+=localCoef[Ncur]*localCoef[j]/alpha;
	double td=localCoef[Ncur]*localCoef[Ncur]/alpha;
	Ncur*=2;	
	std::copy(tmp.begin(),tmp.begin()+Ncur,localCoef.begin());
	localCoef[Ncur]=td;
	//localCoef[Ncur]=1./(alpha*196.0);
	
	alpha -= (2.0/(alpha*196.0));
	/*
	printf("order %ld: ",i);
	for (j=0;j<=Ncur;j++) printf(" %g",localCoef[j]/alpha);
	printf(" ; alpha=%g\n",alpha);	
	*/
      }
    for (i=0;i<Nmax;i++) 
      localCoef[i]/=alpha;
    alpha=localCoef[Nmax-2];
    localCoef[Nmax-2]=0;
    
    localCoef[Nmax-5]+=alpha/12.;
    localCoef[Nmax-4]-=2.*alpha/3.;
    localCoef[Nmax-2]+=2.*alpha/3.;
    localCoef[Nmax-1]-=alpha/12.;
    for (j=0;j<Nmax;j++) localCoef[j]=-localCoef[j];
    std::reverse(localCoef.begin(),localCoef.end());
    /*
    printf("Final result for W-: ");
    for (j=0;j<Nmax;j++) printf(" %g",localCoef[j]);
    printf("\n");
    */
    
    //exit(0);
  }
  /*
  void create_record_1D(dataT *tab,dataT bVal0,dataT bVal1, const int dir, const long delta)
  {
    std::vector<double> &lk = lookup[dir];
    lk.clear();
    const long M=(dim[dir]+2);
    const double hi=step_inv[dir];
 
    //double *A=&tmpArr[0];
    double A[4*M];
    double *a=A;
    dataT *t=&tab[delta];
    double e[5];
    double f[2];

    bool thirdOrder[2]={false,false};
    int row;
    int stopRow=0;
    
    switch (bCond[dir][0])
      {
      case boundaryTypeV::BT_PERIODIC:
    	a[0]=1.0/6.0;a[1]=4.0/6.0;a[2]=1.0/6.0;a[3]=tab[dim[dir]*delta];
	break;
      case boundaryTypeV::BT_FLAT:
    	a[0]=-0.5*hi;a[1]=0;a[2]=0.5*hi;a[3]=0;
	break;
      case boundaryTypeV::BT_NATURAL:
	{
	  double hi2=hi*hi;
	  a[0]=hi2;a[1]=-2.0*hi2;a[2]=hi2;a[3]=0;
	}
	break;
      case boundaryTypeV::BT_SUPRA_NATURAL:
	{
	  double hi3=hi*hi*hi;
	  a[0]=hi3;a[1]=-3.0*hi3;a[2]=3.0*hi3;e[0]=-hi3;e[1]=0;a[3]=0;
	}
	thirdOrder[0]=true;
	break;
      case boundaryTypeV::BT_NOT_A_KNOT:
	a[0]=-1.0;a[1]=4.0;a[2]=-6.0;e[0]=a[1];e[1]=a[0];a[3]=0;
	thirdOrder[0]=true;
	break;
      case boundaryTypeV::BT_LOCAL:
    	a[0]=-0.5*hi;a[1]=0;a[2]=0.5*hi;a[3]=bVal0;
	break;
      }

    a+=4;
    for (int i=1; i<M-1; i++, t+=delta, a+=4) {
      a[0] = 1.0/6.0;
      a[1] = 4.0/6.0;
      a[2] = 1.0/6.0;
      a[3] = *t;
    }

    a=&A[(M-1)*4];
    switch (bCond[dir][1])
      {
      case boundaryTypeV::BT_PERIODIC:
    	a[0]=1.0/6.0;a[1]=4.0/6.0;a[2]=1.0/6.0;a[3]=tab[delta];
	break;
      case boundaryTypeV::BT_FLAT:
    	a[0]=-0.5*hi;a[1]=0;a[2]=0.5*hi;a[3]=0;
	break;
      case boundaryTypeV::BT_NATURAL:
	{
	  double hi2=hi*hi;
	  a[0]=hi2;a[1]=-2.0*hi2;a[2]=hi2;a[3]=0;
	}
	break;
      case boundaryTypeV::BT_SUPRA_NATURAL:
	{
	  double hi3=hi*hi*hi;
	  a[0]=hi3;a[1]=-3.0*hi3;a[2]=3.0*hi3;e[0]=-hi3;e[1]=0;a[3]=0;
	}
	thirdOrder[1]=true;
	break;
      case boundaryTypeV::BT_NOT_A_KNOT:
	a[0]=-1.0;a[1]=4.0;a[2]=-6.0;f[0]=a[1];f[1]=a[0];a[3]=0;
	thirdOrder[1]=true;
	break;
      case boundaryTypeV::BT_LOCAL:
    	a[0]=-0.5*hi;a[1]=0;a[2]=0.5*hi;a[3]=bVal1;
	break; 
      }

    if (!thirdOrder[0])
      {
	A[4*(0)+1] /= A[4*(0)+0];
	A[4*(0)+2] /= A[4*(0)+0];
	A[4*(0)+3] /= A[4*(0)+0];lk.push_back(1./A[4*(0)+0]);
	A[4*(0)+0] = 1.0;
	A[4*(1)+1] -= A[4*(1)+0]*A[4*(0)+1];
	A[4*(1)+2] -= A[4*(1)+0]*A[4*(0)+2];
	A[4*(1)+3] -= A[4*(1)+0]*A[4*(0)+3];lk.push_back(-A[4*(1)+0]);
	A[4*(1)+0] = 0.0;
	A[4*(1)+2] /= A[4*(1)+1];
	A[4*(1)+3] /= A[4*(1)+1];lk.push_back(1./A[4*(1)+1]);
	A[4*(1)+1] = 1.0;
	row=2;
      }
    else
      {
	A[4*(0)+1] /= A[4*(0)+0];
	A[4*(0)+2] /= A[4*(0)+0];
	A[4*(0)+3] /= A[4*(0)+0];lk.push_back(1./A[4*(0)+0]);
	e[0]/=A[4*(0)+0];
	e[1]/=A[4*(0)+0];
	A[4*(0)+0] = 1.0;

	A[4*(1)+1] -= A[4*(1)+0]*A[4*(0)+1];
	A[4*(1)+2] -= A[4*(1)+0]*A[4*(0)+2];
	A[4*(1)+3] -= A[4*(1)+0]*A[4*(0)+3];lk.push_back(-A[4*(1)+0]);
	e[2]= -e[0]*A[4*(1)+0];
	e[3]= -e[1]*A[4*(1)+0];
	A[4*(1)+0] = 0.0;
	A[4*(1)+2] /= A[4*(1)+1];
	A[4*(1)+3] /= A[4*(1)+1];lk.push_back(1./A[4*(1)+1]);
	e[2]/=A[4*(1)+1];
	e[3]/=A[4*(1)+1];
	A[4*(1)+1] = 1.0;

	A[4*(2)+1] -= A[4*(2)+0]*A[4*(2-1)+2];
	A[4*(2)+2] -= A[4*(2)+0]*e[2];
	A[4*(2)+3] -= A[4*(2)+0]*A[4*(2-1)+3];lk.push_back(-A[4*(2)+0]);
	e[4] = -e[3]*A[4*(2)+0];
	A[4*(2)+2] /= A[4*(2)+1];
	A[4*(2)+3] /= A[4*(2)+1];lk.push_back(1.0/A[4*(2)+1]);
	e[4]/=A[4*(2)+1];
	A[4*(2)+0] = 0.0;
	A[4*(2)+1] = 1.0;

	A[4*(3)+1] -= A[4*(3)+0]*A[4*(3-1)+2];
	A[4*(3)+2] -= A[4*(3)+0]*e[4];
	A[4*(3)+3] -= A[4*(3)+0]*A[4*(3-1)+3];lk.push_back(-A[4*(3)+0]);
	A[4*(3)+2] /= A[4*(3)+1];
	A[4*(3)+3] /= A[4*(3)+1];lk.push_back(1./A[4*(3)+1]);
	A[4*(3)+0] = 0.0;
	A[4*(3)+1] = 1.0;

	row=4;stopRow=2;
      }
    
    // The central part does not depend on boundary conditions
    a=&A[4*(row-1)+0];
    
    for (; row <M-1; row++,a+=4) {
      a[4+1] -= a[4+0]*a[2];
      a[4+3] -= a[4+0]*a[3];lk.push_back(-a[4+0]);
      a[4+2] /= a[4+1];
      a[4+3] /= a[4+1];lk.push_back(1./a[4+1]);
      a[4+0] = 0.0;
      a[4+1] = 1.0;
    }
    
    // but the last row does 
    if (thirdOrder[1])
      {
	f[1] -= f[0] * A[4*(M-5)+2];
	A[4*(M-1)+3] -= f[0]*A[4*(M-5)+3];lk.push_back(-f[0]);

	A[4*(M-1)+0] -= f[1]*A[4*(M-4)+2];
	A[4*(M-1)+3] -= f[1]*A[4*(M-4)+3];lk.push_back(-f[1]);

	A[4*(M-1)+1] -= A[4*(M-1)+0]*A[4*(M-3)+2];
	A[4*(M-1)+3] -= A[4*(M-1)+0]*A[4*(M-3)+3];lk.push_back(-A[4*(M-1)+0]);

	A[4*(M-1)+2] -= A[4*(M-1)+1]*A[4*(M-2)+2];
	A[4*(M-1)+3] -= A[4*(M-1)+1]*A[4*(M-2)+3];lk.push_back(-A[4*(M-1)+1]);

	A[4*(M-1)+3] /= A[4*(M-1)+2];lk.push_back(1.0/A[4*(M-1)+2]);
	A[4*(M-1)+2] = 1.0;
      }
    else
      {
	A[4*(M-1)+1] -= A[4*(M-1)+0]*A[4*(M-3)+2];
	A[4*(M-1)+3] -= A[4*(M-1)+0]*A[4*(M-3)+3];lk.push_back(-A[4*(M-1)+0]);

	A[4*(M-1)+2] -= A[4*(M-1)+1]*A[4*(M-2)+2];
	A[4*(M-1)+3] -= A[4*(M-1)+1]*A[4*(M-2)+3];lk.push_back(-A[4*(M-1)+1]);

	A[4*(M-1)+3] /= A[4*(M-1)+2];lk.push_back(1.0/A[4*(M-1)+2]);
	A[4*(M-1)+2] = 1.0;
      }

    t=&tab[delta*(M-1)];
    *t = A[4*(M-1)+3];
    t-=delta;
    
    for (row=M-2; row>stopRow; row--,t-=delta)
      {
	(*t) = A[4*(row)+3] - A[4*(row)+2]*t[delta];lk.push_back(- A[4*(row)+2]);
      }

    if (!thirdOrder[0])
      {
	(*t) = A[4*(0)+3] - A[4*(0)+1]*t[delta] - A[4*(0)+2]*t[2*delta];
	lk.push_back(-A[4*(0)+1]);lk.push_back(-A[4*(0)+2]);
      }
    else
      {
	a=&A[4*stopRow];

	(*t) = a[4*(0)+3] - a[4*(0)+2]*t[delta]  - e[4] * t[2*delta];
	lk.push_back(-a[4*(0)+2]);lk.push_back(-e[4]);
	t-=delta;a-=4;
	(*t) = a[4*(0)+3] - a[4*(0)+2]*t[delta]  - e[2] * t[2*delta]- e[3] * t[3*delta];
	lk.push_back(-a[4*(0)+2]);lk.push_back(-e[2]);lk.push_back(-e[3]);
	t-=delta;a-=4;
	(*t) = a[4*(0)+3] - a[4*(0)+1]*t[delta] - a[4*(0)+2]*t[2*delta] - e[0] * t[3*delta]- e[1] * t[4*delta];	
	lk.push_back(-a[4*(0)+1]);lk.push_back(-a[4*(0)+2]);
	lk.push_back(-e[0]);lk.push_back(-e[1]);
      }
  
  }
  
  void create_multi_1D(dataT *tab,dataT bVal0,dataT bVal1, const int dir, const long delta)
  {
    const std::vector<double> &lk = lookup[dir];

    if (lk.size()==0)
      {
	fprintf(stderr,"in localspline: create_record_1D not implemented correctly yet ! \n");
	exit(-1);
	create_record_1D(tab,bVal0,bVal1,dir,delta);
	return;
      }

    bool thirdOrder[2]={false,false};
    
    switch (bCond[dir][0])
      {
      case boundaryTypeV::BT_PERIODIC:
	bVal0=tab[dim[dir]*delta];
	break;
      case boundaryTypeV::BT_FLAT:
    	bVal0=0;
	break;
      case boundaryTypeV::BT_NATURAL:
	bVal0=0;
	break;
      case boundaryTypeV::BT_SUPRA_NATURAL:
	bVal0=0;
	thirdOrder[0]=true;
	break;
      case boundaryTypeV::BT_NOT_A_KNOT:
	bVal0=0;
	thirdOrder[0]=true;
	break;
      }
    switch (bCond[dir][1])
      {
      case boundaryTypeV::BT_PERIODIC:
    	bVal1=tab[delta];
	break;
      case boundaryTypeV::BT_FLAT:
    	bVal1=0;
	break;
      case boundaryTypeV::BT_NATURAL:
	bVal1=0;
	break;
      case boundaryTypeV::BT_SUPRA_NATURAL:
	bVal1=0;
	thirdOrder[1]=true;
	break;
      case boundaryTypeV::BT_NOT_A_KNOT:
	bVal1=0;
	thirdOrder[1]=true;
	break;
      }

    

    int row;
    int stopRow=0;
    const long M=(dim[dir]+2);
    double *t=tab;
    const double *l=&lk[0];

    t[delta*(M-1)]=bVal1;

    if (!thirdOrder[0])
      {
	t[0]=bVal0*l[0];
	t[delta]+=l[1]*t[0];
	t[delta]*=l[2];
	l+=3;
	row=2;
      }
    else
      {
	t[0]=bVal0*l[0];
	t[delta]+=l[1]*t[0];
	t[delta]*=l[2];
	t[2*delta]+=t[delta]*l[3];
	t[2*delta]*=l[4];
	t[3*delta]+=t[2*delta]*l[5];
	t[3*delta]*=l[6];
	l+=7;
	row=4;stopRow=2;
      }

    t=&tab[delta*(row-1)];    
    
    for (; row <M-1; row++,t+=delta,l+=2)     
      t[delta]=l[1]*(t[delta]+t[0]*l[0]);
    
    if (thirdOrder[1])
      {
	t=&tab[delta*(M-5)];
	t[4]=l[4]*(t[4]+l[0]*t[0]+l[1]*t[1]+l[2]*t[2]+l[3]*t[3]);
	l+=5;
      }
    else
      {
	t=&tab[delta*(M-3)];
	t[2]=l[2]*(t[2]+l[0]*t[0]+l[1]*t[1]);
	l+=3;
      }
    
    t=&tab[delta*(M-2)];
    
    for (row=M-2; row>stopRow; row--,t-=delta,l++)
      (*t) += t[delta]*l[0];
   
    if (!thirdOrder[0])
      {
	(*t) += l[0]*t[delta] + l[1]*t[2*delta];
      }
    else
      {
	*t += l[0]*t[delta]+l[1]*t[2*delta];
	t-=delta;l+=2;
	*t += l[0]*t[delta]+l[1]*t[2*delta]+l[2]*t[3*delta];
	t-=delta;l+=3;
	*t += l[0]*t[delta]+l[1]*t[2*delta]+l[2]*t[3*delta]+l[3]*t[4*delta];
      }
  
  }
  */
  void create_1D(dataT *tab,dataT bVal0,dataT bVal1, const int dir, const long delta)
  {
    const long M=(dim[dir]+2);
    const double hi=step_inv[dir];
 
    double *A=&tmpArr[0];
    double *a=A;
    dataT *t=&tab[delta];
    double e[5];
    double f[2];

    bool thirdOrder[2]={false,false};
    int row;
    int stopRow=0;
    
    switch (bCond[dir][0])
      {
      case boundaryTypeV::BT_PERIODIC:
    	a[0]=1.0/6.0;a[1]=4.0/6.0;a[2]=1.0/6.0;a[3]=tab[dim[dir]*delta];
	break;
      case boundaryTypeV::BT_FLAT:
    	a[0]=-0.5*hi;a[1]=0;a[2]=0.5*hi;a[3]=0;
	break;
      case boundaryTypeV::BT_NATURAL:
	{
	  double hi2=hi*hi;
	  a[0]=hi2;a[1]=-2.0*hi2;a[2]=hi2;a[3]=0;
	}
	break;
      case boundaryTypeV::BT_SUPRA_NATURAL:
	{
	  fprintf(stderr,"spline: BT_SUPRA_NATURAL not implemented correctly yet, sorry!\n");
	  exit(-1);
	  double hi3=hi*hi*hi;
	  a[0]=hi3;a[1]=-3.0*hi3;a[2]=3.0*hi3;e[0]=-hi3;e[1]=0;a[3]=0;
	}
	thirdOrder[0]=true;
	break;
      case boundaryTypeV::BT_NOT_A_KNOT:
	fprintf(stderr,"spline: BT_NOT_A_KNOT not implemented correctly yet, sorry!\n");
	exit(-1);
	a[0]=-1.0;a[1]=4.0;a[2]=-6.0;e[0]=a[1];e[1]=a[0];a[3]=0;
	thirdOrder[0]=true;
	break;
      case boundaryTypeV::BT_LOCAL:
    	a[0]=-0.5*hi;a[1]=0;a[2]=0.5*hi;a[3]=bVal0;
	break;
      }

    a+=4;
    for (int i=1; i<M-1; i++, t+=delta, a+=4) {
      a[0] = 1.0/6.0;
      a[1] = 4.0/6.0;
      a[2] = 1.0/6.0;
      a[3] = *t;
    }

    a=&A[(M-1)*4];
    switch (bCond[dir][1])
      {
      case boundaryTypeV::BT_PERIODIC:
    	a[0]=1.0/6.0;a[1]=4.0/6.0;a[2]=1.0/6.0;a[3]=tab[delta];
	break;
      case boundaryTypeV::BT_FLAT:
    	a[0]=-0.5*hi;a[1]=0;a[2]=0.5*hi;a[3]=0;
	break;
      case boundaryTypeV::BT_NATURAL:
	{
	  double hi2=hi*hi;
	  a[0]=hi2;a[1]=-2.0*hi2;a[2]=hi2;a[3]=0;
	}
	break;
      case boundaryTypeV::BT_SUPRA_NATURAL:
	{
	  fprintf(stderr,"spline: BT_SUPRA_NATURAL not implemented correctly yet, sorry!\n");
	  exit(-1);
	  double hi3=hi*hi*hi;
	  a[0]=hi3;a[1]=-3.0*hi3;a[2]=3.0*hi3;e[0]=-hi3;e[1]=0;a[3]=0;
	}
	thirdOrder[1]=true;
	break;
      case boundaryTypeV::BT_NOT_A_KNOT:
	fprintf(stderr,"spline: BT_NOT_A_KNOT not implemented correctly yet, sorry!\n");
	exit(-1);
	a[0]=-1.0;a[1]=4.0;a[2]=-6.0;f[0]=a[1];f[1]=a[0];a[3]=0;
	thirdOrder[1]=true;
	break;
      case boundaryTypeV::BT_LOCAL:
    	a[0]=-0.5*hi;a[1]=0;a[2]=0.5*hi;a[3]=bVal1;
	break; 
      }

    if (!thirdOrder[0])
      {
	A[4*(0)+1] /= A[4*(0)+0];
	A[4*(0)+2] /= A[4*(0)+0];
	A[4*(0)+3] /= A[4*(0)+0];
	A[4*(0)+0] = 1.0;
	A[4*(1)+1] -= A[4*(1)+0]*A[4*(0)+1];
	A[4*(1)+2] -= A[4*(1)+0]*A[4*(0)+2];
	A[4*(1)+3] -= A[4*(1)+0]*A[4*(0)+3];
	A[4*(1)+0] = 0.0;
	A[4*(1)+2] /= A[4*(1)+1];
	A[4*(1)+3] /= A[4*(1)+1];
	A[4*(1)+1] = 1.0;
	row=2;
      }
    else
      {
	A[4*(0)+1] /= A[4*(0)+0];
	A[4*(0)+2] /= A[4*(0)+0];
	A[4*(0)+3] /= A[4*(0)+0];
	e[0]/=A[4*(0)+0];
	e[1]/=A[4*(0)+0];
	A[4*(0)+0] = 1.0;

	A[4*(1)+1] -= A[4*(1)+0]*A[4*(0)+1];
	A[4*(1)+2] -= A[4*(1)+0]*A[4*(0)+2];
	A[4*(1)+3] -= A[4*(1)+0]*A[4*(0)+3];
	e[2]= -e[0]*A[4*(1)+0];
	e[3]= -e[1]*A[4*(1)+0];
	A[4*(1)+0] = 0.0;
	A[4*(1)+2] /= A[4*(1)+1];
	A[4*(1)+3] /= A[4*(1)+1];
	e[2]/=A[4*(1)+1];
	e[3]/=A[4*(1)+1];
	A[4*(1)+1] = 1.0;

	A[4*(2)+1] -= A[4*(2)+0]*A[4*(2-1)+2];
	A[4*(2)+2] -= A[4*(2)+0]*e[2];
	A[4*(2)+3] -= A[4*(2)+0]*A[4*(2-1)+3];
	e[4] = -e[3]*A[4*(2)+0];
	A[4*(2)+2] /= A[4*(2)+1];
	A[4*(2)+3] /= A[4*(2)+1];
	e[4]/=A[4*(2)+1];
	A[4*(2)+0] = 0.0;
	A[4*(2)+1] = 1.0;

	A[4*(3)+1] -= A[4*(3)+0]*A[4*(3-1)+2];
	A[4*(3)+2] -= A[4*(3)+0]*e[4];
	A[4*(3)+3] -= A[4*(3)+0]*A[4*(3-1)+3];
	A[4*(3)+2] /= A[4*(3)+1];
	A[4*(3)+3] /= A[4*(3)+1];
	A[4*(3)+0] = 0.0;
	A[4*(3)+1] = 1.0;

	row=4;stopRow=2;
      }
   
    // The central part does not depend on boundary conditions
    
    a=&A[4*(row-1)+0];
    for (; row <M-1; row++,a+=4) {
      a[4+1] -= a[4+0]*a[2];
      a[4+3] -= a[4+0]*a[3];
      a[4+2] /= a[4+1];
      a[4+3] /= a[4+1];
      a[4+0] = 0.0;
      a[4+1] = 1.0;
    }
    
    /*
    for (; row <M-1; row++) {
      A[4*(row)+1] -= A[4*(row)+0]*A[4*(row-1)+2];
      A[4*(row)+3] -= A[4*(row)+0]*A[4*(row-1)+3];
      A[4*(row)+2] /= A[4*(row)+1];
      A[4*(row)+3] /= A[4*(row)+1];
      A[4*(row)+0] = 0.0;
      A[4*(row)+1] = 1.0;
    }
    */

    // but the last row does 
    if (thirdOrder[1])
      {
	f[1] -= f[0] * A[4*(M-5)+2];
	A[4*(M-1)+3] -= f[0]*A[4*(M-5)+3];

	A[4*(M-1)+0] -= f[1]*A[4*(M-4)+2];
	A[4*(M-1)+3] -= f[1]*A[4*(M-4)+3];

	A[4*(M-1)+1] -= A[4*(M-1)+0]*A[4*(M-3)+2];
	A[4*(M-1)+3] -= A[4*(M-1)+0]*A[4*(M-3)+3];

	A[4*(M-1)+2] -= A[4*(M-1)+1]*A[4*(M-2)+2];
	A[4*(M-1)+3] -= A[4*(M-1)+1]*A[4*(M-2)+3];

	A[4*(M-1)+3] /= A[4*(M-1)+2];
	A[4*(M-1)+2] = 1.0;
      }
    else
      {
	A[4*(M-1)+1] -= A[4*(M-1)+0]*A[4*(M-3)+2];
	A[4*(M-1)+3] -= A[4*(M-1)+0]*A[4*(M-3)+3];

	A[4*(M-1)+2] -= A[4*(M-1)+1]*A[4*(M-2)+2];
	A[4*(M-1)+3] -= A[4*(M-1)+1]*A[4*(M-2)+3];

	A[4*(M-1)+3] /= A[4*(M-1)+2];
	A[4*(M-1)+2] = 1.0;
      }

    t=&tab[delta*(M-1)];
    *t = A[4*(M-1)+3];
    t-=delta;

    // now we have an U type matrix -> substitute back
    /*
    a=&A[4*(M-2)+0];
    for (row=M-2; row>stopRow; row--,t-=delta, a-=4)
      (*t) = a[3] - a[2]*t[delta];
*/
    
    for (row=M-2; row>stopRow; row--,t-=delta)
      (*t) = A[4*(row)+3] - A[4*(row)+2]*t[delta];
    
    if (!thirdOrder[0])
      {
	(*t) = A[4*(0)+3] - A[4*(0)+1]*t[delta] - A[4*(0)+2]*t[2*delta];
      }
    else
      {
	a=&A[4*stopRow];
	(*t) = a[4*(0)+3] - a[4*(0)+2]*t[delta]  - e[4] * t[2*delta];
	t-=delta;a-=4;
	(*t) = a[4*(0)+3] - a[4*(0)+2]*t[delta]  - e[2] * t[2*delta]- e[3] * t[3*delta];
	t-=delta;a-=4;
	(*t) = a[4*(0)+3] - a[4*(0)+1]*t[delta] - a[4*(0)+2]*t[2*delta] - e[0] * t[3*delta]- e[1] * t[4*delta];	
      }
  }

  /*
  void create_1D_LU(dataT *tab, dataT *b0, dataT *b1, const int dir, const long delta)
  {
    long i,j;
    const double h=step[dir];
    const long imax=dim[dir];
    const double MH3=-h/3.;
    const double P3H=3./h;
    const double M3H=-P3H;
    const double P3H6=P3H/6.;
    const double M3H6=M3H/6.;
    
    dataT *t=tab;
       
    // L.Phi = F  
    (*t) -= MH3*(*b0);
    for (i=1;i<imax;i++,t+=delta)
      t[delta]-=L[i-1]*(*t);
    
    (*b1)-= (*t)*P3H*L_back[dir] + (*(t-delta))*M3H*L[imax-2];

    // U.n = F
    (*b1) /= (P3H6*D_back[dir]);
    (*t)=(6*(*t)-(*b1))*D[imax-1];
    t-=delta;
    for (i=imax-2;i>0;i--,t-=delta)
      (*t)=(6*(*t)-t[delta])*D[i];

    (*t)=(6*(*t)-2*t[delta])*D[0];
    (*b0)=((*b0)-P3H6*t[delta])/M3H6;
  }
  */

#if defined(__ICC) || defined(__INTEL_COMPILER)
// 175 - disable spurious subscript out of range warning 
#pragma warning push
#pragma warning disable 175
#endif

  void create_1D()
  {
    long size=4*(dim[0]+2);

    //tmpArr.resize(size);
    int res=posix_memalign((void**)&tmpArr,32,size*sizeof(dataT));
    //tmpArr=(dataT*)malloc(size*sizeof(dataT));
    
    create_1D(&data[0],data[0],data[dim[0]+1],0,stride[0]);

    //tmpArr.clear();
    free(tmpArr);
  }


  void create_2D()
  {
    long i;
    dataT *tab;
    long size=4l*(2l+(long)std::max(dim[0],dim[1]));

    //tmpArr.resize(size);
    int res=posix_memalign((void**)&tmpArr,32,size*sizeof(dataT));
    //tmpArr=(dataT*)malloc(size*sizeof(dataT));

    tab=&data[0+stride[1]];
    for (i=-1+1;i<=dim[1]-1;i++,tab+=stride[1])
      create_1D(tab,tab[0],tab[(dim[0]+1)*stride[0]],0,stride[0]);
      
    tab=&data[0];
    for (i=-1;i<=dim[0];i++,tab+=stride[0])
      create_1D(tab,tab[0],tab[(dim[1]+1)*stride[1]],1,stride[1]);
    
    //tmpArr.clear();
    free(tmpArr);
  }

  void create_3D()
  {
    long i,j;
    dataT *tab;
    long size=4l*(2l+(long)std::max(std::max(dim[0],dim[1]),dim[2]));

    //tmpArr.resize(size);
    int res=posix_memalign((void**)&tmpArr,32,size*sizeof(dataT));
    //tmpArr=(dataT*)malloc(size*sizeof(dataT));

    tab=&data[0+stride[2]+stride[1]];
    for (i=-1+1;i<=dim[2]-1;i++,tab+=2*stride[1])
      for (j=-1+1;j<=dim[1]-1;j++,tab+=stride[1])
	create_1D(tab,tab[0],tab[(dim[0]+1)*stride[0]],0,stride[0]);
   
    tab=&data[0+stride[2]];
    for (i=-1+1;i<=dim[2]-1;i++,tab+=stride[2]-stride[1])
      for (j=-1;j<=dim[0];j++,tab+=stride[0])
      	create_1D(tab,tab[0],tab[(dim[1]+1)*stride[1]],1,stride[1]);

    tab=&data[0];
    for (i=-1;i<=dim[1];i++)
      for (j=-1;j<=dim[0];j++,tab+=stride[0])
      	create_1D(tab,tab[0],tab[(dim[2]+1)*stride[2]],2,stride[2]);
   
    //tmpArr.clear();
    free(tmpArr);
  }

#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma warning pop
#endif

  void create_ND()
  {
    fprintf(stderr,"localspline with NDIMS=%d are not implemented yet!\n",NDIMS);
    exit(-1);
  }

  void computeLocalBoundary(const int dm,const int dr)
  {
    long i,j;
    dataT *dataPtr=&data[0];
    const long bN=stride[NDIMS]/(dim[dm]+2);
    const long delta=stride[dm];
    const double fac=(dr>0)?(step_inv[dm]):(-step_inv[dm]);

    int w[NDIMS]= {0};
    long k=(dr>0)?((dim[dm]+1)*delta):0;
    dataT *ref=dataPtr+k;

    const int WS=localCoef.size();
    const double *WM=&localCoef[0];
    
    if (dr>0)
      {
	for (i=0;i<bN;i++)
	  {	
	    //printf("k=%ld\n",k);
	    dataT *cur=ref-delta*(WS+1+1);
	
	    (*ref)=
	      WM[0]*cur[1*delta]+WM[1]*cur[2*delta]+WM[2]*cur[3*delta]+WM[3]*cur[4*delta]+WM[4]*cur[5*delta]+
	      WM[5]*cur[6*delta]+WM[6]*cur[7*delta]+WM[7]*cur[8*delta]+WM[8]*cur[9*delta]+WM[9]*cur[10*delta];
	    
	    (*ref)*=fac;
	    //printf("ref=%e\n",*ref);
	    //(*ref)=-5.7;
	
	    w[0]++;k+=bStride[dm][0];
	    if (w[0]==bDim[dm][0]) hlp::getNext<NDIMS-1>(w,k,bDim[dm],bStride[dm]);	
	    ref=dataPtr+k;
	  }
      }
    else
      {
	for (i=0;i<bN;i++)
	  {	
	    //printf("k=%ld\n",k);
	    dataT *cur=ref+1*delta;

	    (*ref)=
	      WM[9]*cur[1*delta]+WM[8]*cur[2*delta]+WM[7]*cur[3*delta]+WM[6]*cur[4*delta]+WM[5]*cur[5*delta]+
	      WM[4]*cur[6*delta]+WM[3]*cur[7*delta]+WM[2]*cur[8*delta]+WM[1]*cur[9*delta]+WM[0]*cur[10*delta];
	    
	    (*ref)*=fac;
	    //printf("ref=%e\n",*ref);
	    //(*ref)=-5.7;
	
	    w[0]++;k+=bStride[dm][0];
	    if (w[0]==bDim[dm][0]) hlp::getNext<NDIMS-1>(w,k,bDim[dm],bStride[dm]);	
	    ref=dataPtr+k;
	  }
      }
  }

  void computeLocalBoundaries()
  {
    for (int dm=0;dm<NDIMS;dm++)
      {
	if (bCond[dm][0] == boundaryTypeV::BT_LOCAL)
	  computeLocalBoundary(dm,0);
	if (bCond[dm][1] == boundaryTypeV::BT_LOCAL)
	  computeLocalBoundary(dm,1);
      }
  }

};

template <int NDIMS,typename dataType>
const int localSpline<NDIMS,dataType>::WS3 = 10;

// template <int NDIMS,typename dataType>
// const double localSpline<NDIMS,dataType>::WM3[localSpline<NDIMS,dataType>::WS3]=
//   {0.2214309755E-5, -1.771447804E-5, 7.971515119E-5, -3.011461267E-4, 1.113797807E-3,
//    -4.145187862E-3, 0.01546473933  , -0.05771376946, 0.2153903385   , -0.8038475846 };

template <int NDIMS,typename dataType>
const double localSpline<NDIMS,dataType>::WM3[10]=
  {0.2214309755E-5, -1.771447804E-5, 7.971515119E-5, -3.011461267E-4, 1.113797807E-3,
   -4.145187862E-3, 0.01546473933  , -0.05771376946, 0.2153903385   , -0.8038475846 };


#endif
