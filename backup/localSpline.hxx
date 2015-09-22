#ifndef __LOCAL_SPLINE_HXX__
#define __LOCAL_SPLINE_HXX__

#include <algorithm>

#include "scale.hxx"

template <
  int NDIMS,
  typename dataType = double
  >
class localSpline
{
public:
  typedef dataType dataT;
  typedef scale<dataT> scaleT;
  
  typedef typename scaleT::scaleTypeT scaleTypeT; 
  typedef typename scaleT::scaleTypeV scaleTypeV;
  typedef typename scaleT::valLocationT valLocationT;
  typedef typename scaleT::valLocationV valLocationV;  

  struct boundaryTypeV {
    enum type {BT_PERIODIC, BT_FLAT, BT_NATURAL, BT_LOCAL};
  };
  typedef typename boundaryTypeV::type boundaryTypeT;

  struct initT {
    double start[NDIMS];
    double stop[NDIMS];
    int N[NDIMS];
    scaleTypeT scaleType[NDIMS];
    valLocationT vlType[NDIMS];

    boundaryTypeT bCondL[NDIMS];
    boundaryTypeT bCondR[NDIMS];
  };

  void init(const initT &ini)
  {
    long i,i0;
    long maxS=0;

    purge();
    stride[0]=1;
    for (i=0;i<NDIMS;i++) 
      {
	sclType[i]=ini.scaleType[i];
	bCondL[i]=ini.bCondL[i];
	bCondR[i]=ini.bCondR[i];
	scl[i]=scaleT::genScale(ini.start[i],ini.stop[i],ini.N[i],sclType[i],ini.vlType[i]);	
	dim[i]=scl[i].size();
	stride[i+1]=stride[i]*dim[i];

	start[i]=scl[i].front();
	stop[i]=scl[i].back();
	rescale[i]=scaleT::rescale(start[i],sclType[i]);
	scaleT::rescale(stop[i],sclType[i]);
	step[i]=(stop[i]-start[i])/(dim[i]-1);

	m3h6[i]=-1./(step[i]*2.);
	p3h6[i]=1./(step[i]*2.);
	m3h[i]=-3./(step[i]);
	p3h[i]=3./(step[i]);
	mh3[i]=-step[i]/3.;

	if (dim[i]>maxS) maxS=dim[i];
      } 

    for (i=0;i<NDIMS;i++) 
      {
	bStride[i][0]=1;
	long k=0;
	for (long j=0;j<NDIMS;j++)
	  {
	    if (i!=j) 
	      {
		if (k<i)
		  bStride[i][k+1]=bStride[i][k]*(dim[j]+2);
		else
		  bStride[i][k+1]=bStride[i][k]*dim[j];

		k++;
	      }
	  }
      }

    if (D.size()==0) 
      {
	D.push_back(4);
	D.push_back(3.5);
	L.push_back(0.25);
      }
    else
      {
	D.resize(std::min((long)D.size(),maxS));
	L.resize(std::min((long)L.size(),maxS-1));
      }

    D.reserve(maxS);
    L.reserve(maxS-1);

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

    
    
    //printf("\n");
  }
  
  localSpline(const initT &ini)
  {
    init(ini);
  }
  
  ~localSpline()
  {

  }

  void create(dataT *dataP)
  {
    create(dataP,dataP+stride[NDIMS]);
  }

  template <class inputIterator>
  void create(const inputIterator &start,const inputIterator &stop)
  {
    data.assign(start,stop);
    //printf("stride[NDIMS]=%ld == %ld\n",stride[NDIMS],std::distance(start,stop));
    //data.resize(stride[NDIMS]);
    for (long i=0;i<NDIMS;i++)
      {
	b0[i].resize(bStride[i][NDIMS-1]);
	b1[i].resize(bStride[i][NDIMS-1]);

	b0[i].assign(bStride[i][NDIMS-1],0);
	b1[i].assign(bStride[i][NDIMS-1],0);
      }

    //std::copy(start,stop,data.begin());

    switch (NDIMS)
      {
      case 1: create_1D(); break;
      case 2: create_2D(); break;
      case 3: create_3D(); break;
      default: create_ND();
      }
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

  dataT evaluate(double x)
  {
    double h=step[0];    
    int N=(x-start[0])/h;    
    double result;
  
    if (x>=stop[0]) 
      {
	return (1./6.)*(data[dim[0]-2]+4*data[dim[0]-1]+b1[0][0]);
      }
    if (x<start[0]) return 0;

    x=x-(start[0]+N*h);
    
    if (N==0) result=b0[0][0]*BV2(h,h-x)+data[0]*BV1(h,h-x)+data[1]*BV1(h,x)+data[2]*BV2(h,x);
    else if (N>=dim[0]-2) result=data[N-1]*BV2(h,h-x)+data[N]*BV1(h,h-x)+data[N+1]*BV1(h,x)+b1[0][0]*BV2(h,x);
    else result=data[N-1]*BV2(h,h-x)+data[N]*BV1(h,h-x)+data[N+1]*BV1(h,x)+data[N+2]*BV2(h,x);

    return (1./(6.*h*h*h))*result;   
  }

  dataT evaluate(double x,double y)
  {
    double hx=step[0];    
    double hy=step[1];    
    int Nx=(x-start[0])/hx;    
    int Ny=(y-start[1])/hy;    
    double result;

    if (x>=stop[0]) 
      {
	if (y>=stop[1]) return 0;
	else if (y<start[1]) return 0;
	return 0;
      }

    if (x<start[0]) 
      {
	if (y>=stop[1]) return 0;
	else if (y<start[1]) return 0;
	return 0;
      }

    x=x-(start[0]+Nx*hx);
    y=y-(start[1]+Ny*hy);

    const long ref[4]={Nx+(Ny-1)*stride[1],Nx+Ny*stride[1],Nx+(Ny+1)*stride[1],Nx+(Ny+2)*stride[1]};
    const double d0[4]={data[ref[0]-1],data[ref[1]-1],data[ref[2]-1],data[ref[3]-1]};
    const double d1[4]={data[ref[0]+2],data[ref[1]+2],data[ref[2]+2],data[ref[3]+2]};

    if (Nx==0) return 0;
    if (Ny==0) return 0;
    if (Nx>=dim[0]-2) return 0;
    if (Ny>=dim[1]-2) return 0;
    
    result = 
      BV2(hy,hy-y)*(data[ref[0]-1]*BV2(hx,hx-x)+data[ref[0]]*BV1(hx,hx-x)+data[ref[0]+1]*BV1(hx,x)+data[ref[0]+2]*BV2(hx,x)) +
      BV1(hy,hy-y)*(data[ref[1]-1]*BV2(hx,hx-x)+data[ref[1]]*BV1(hx,hx-x)+data[ref[1]+1]*BV1(hx,x)+data[ref[1]+2]*BV2(hx,x)) +
      BV1(hy,y)*(data[ref[2]-1]*BV2(hx,hx-x)+data[ref[2]]*BV1(hx,hx-x)+data[ref[2]+1]*BV1(hx,x)+data[ref[2]+2]*BV2(hx,x)) +
      BV2(hy,y)*(data[ref[3]-1]*BV2(hx,hx-x)+data[ref[3]]*BV1(hx,hx-x)+data[ref[3]+1]*BV1(hx,x)+data[ref[3]+2]*BV2(hx,x));
	
    return 1./(36.*hx*hx*hx*hy*hy*hy)*result;
    
  }

  dataT evaluate(double x,double y,double z)
  {
    return 0;
  }

  void purge()
  {
    data.clear();
    for (long i=0;i<NDIMS;i++)
      {
	b0[i].clear();
	b1[i].clear();
      }
  }

private:
  double m3h6[NDIMS],p3h6[NDIMS];
  double m3h[NDIMS],p3h[NDIMS];
  double mh3[NDIMS];
  //static const double A[3];

  std::vector<double> D;
  std::vector<double> L;
  double D_back[NDIMS];
  double L_back[NDIMS];

private:
  int dim[NDIMS];
  long stride[NDIMS+1];
  long bStride[NDIMS][NDIMS];  
  bool rescale[NDIMS];
  double start[NDIMS];
  double stop[NDIMS];
  double step[NDIMS];
  scaleTypeT sclType[NDIMS];
  boundaryTypeT bCondL[NDIMS];
  boundaryTypeT bCondR[NDIMS];

  std::vector<double> scl[NDIMS];
  std::vector<dataT> data;
  std::vector<dataT> b0[NDIMS];
  std::vector<dataT> b1[NDIMS];

  double BV2(double h, double x) {return x*x*x;}
  double BV1(double h, double x) {return h*h*h+3*h*h*x+3*h*x*x-3*x*x*x;}

  void create_1D()
  {
    create_1D(&data[0],&b0[0][0],&b1[0][0],0,stride[0]);
  }

  //void create_1D(dataT *tab, dataT *b0, dataT *b1, int dir)
  void create_1D(dataT *tab, dataT *b0, dataT *b1, const int dir, const long delta)
  {
    long i,j;
    const double h=step[dir];
    const long imax=dim[dir];
    const double MH3=-h/3.;
    const double P3H=3./h;
    const double M3H=-P3H;
    const double P3H6=P3H/6.;
    const double M3H6=M3H/6.;
    /*
    double MH3=mh3[dir];
    double P3H=p3h[dir];
    double M3H=m3h[dir];
    double P3H6=p3h6[dir];
    double M3H6=m3h6[dir];
    const long imax=dim[dir];
    const long delta=stride[dir];
    */

    //double old=tab[1];

    //b0[0]=b1[0]=0;

    // L.Phi = F  
    //printf("L(mh3=%g):\n (%ld,%g)",MH3,(long)0,b0[0]);
    //for (i=0;i<dim[dir];i++) printf(" (%ld,%g,%g)",i+1,tab[i],(i>0)?L[i-1]:0);
    //printf(" (%ld,%g)\n",i+1,b1[0]);
 
    tab[0] -= MH3*b0[0];
    for (i=1,j=delta;i<imax;i++,j+=delta)
      tab[j]-=L[i-1]*tab[j-delta];

    //printf("b1[0]-= tab[j-delta]*P3H*L_back[dir] + tab[j-2*delta]*M3H*L[imax-2]\n");
    //printf("%g -= %g*%g*%g+%g*%g*%g\n",b1[0],tab[j-delta],P3H,L_back[dir],tab[j-2*delta],M3H,L[imax-2]);

    b1[0]-= tab[j-delta]*P3H*L_back[dir] + tab[j-2*delta]*M3H*L[imax-2];
    

    //printf("L(mh3=%g):\n (%ld,%g)",MH3,(long)0,b0[0]);
    //for (i=0;i<dim[dir];i++) printf(" (%ld,%g)",i+1,tab[i]);
    //printf(" (%ld,%g)\n",i+1,b1[0]);

    // U.n = F
    b1[0]/=(P3H6*D_back[dir]);
    j=(imax-1)*delta;
    //printf("tab=(%g-%g)/%g\n",tab[j],b1[0],D[imax-1]);
    tab[j]=(6*tab[j]-b1[0])/D[imax-1];
    
    for (i=imax-2,j=(imax-2)*delta;i>0;i--,j-=delta)
      tab[j]=(6*tab[j]-tab[j+delta])/D[i];

    tab[0]=(6*tab[0]-2*tab[delta])/D[0];
    b0[0]=(b0[0]-P3H6*tab[delta])/M3H6;

    //printf("U: (%ld,%g)",(long)0,b0[0]);
    //for (i=0;i<dim[dir];i++) printf(" (%ld,%g)",i+1,tab[i]);
    //printf(" (%ld,%g)\n",i+1,b1[0]);

    //printf("old=%g == %g/6+2.%g/6+%g/6 = %g\n",old,data[0],data[1],data[2],data[0]/6+2*data[1]/3+data[2]/6);
  }

  void create_2D(dataT *tab, int dir)
  {
    
  }

  void create_2D()
  {
    long i;
    dataT *tab=&data[0];
    
    for (i=0;i<dim[1];i++,tab+=stride[0])
      create_1D(tab,&b0[0][i],&b1[0][i],0,stride[1]);
    
    /*
    // SOLVE 2 equations for second derivatives on boundary here !
    */
    tab=&data[0];
    for (i=0;i<dim[0];i++,tab+=stride[1])
      create_1D(tab,&b0[1][i+1],&b1[1][i+1],1,stride[0]);
     
    create_1D(&b0[1][1],&b0[1][0],&b0[1].back(),1,1);
    create_1D(&b1[1][1],&b1[1][0],&b1[1].back(),1,1);
    
  }

  void create_3D()
  {
    
  }

  void create_ND()
  {
    fprintf(stderr,"localspline with NDIMS=%d are not implemented yet!\n",NDIMS);
    exit(-1);
  }
};

//template <int NDIMS,typename dataType>
//const double localSpline<NDIMS,dataType>::A[3]={1./6.,4./6.,1./6.};

#endif
