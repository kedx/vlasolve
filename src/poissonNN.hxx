#ifndef __POISSON_NN_HXX__
#define __POISSON_NN_HXX__

#ifdef HAVE_FFTW3

#include "FFTW_tools.hxx"
#include "helpers.hxx"
#include "global.h"
#include <string.h>


template <typename dataType,int NDIMS>
class poissonFFTSolver
{
public:
  typedef poissonFFTSolver<dataType,NDIMS> myT;
  typedef FFTW_tools<dataType,NDIMS> FFTWT;
  typedef dataType dataT;
  typedef typename FFTWT::cplxT cplxDataT;
  typedef typename FFTWT::planT planT;

  static const int PERIODIC = ((1<<NDIMS)-1);
  
private:
  bool initialized;
  
  bool threadsInitialized;
  int nthreads;
  int periodicity;

  struct solverDataT 
  {
    planT pf;
    planT pb;
    planT pb3;
    bool havePlans;
    bool havePB3;
    dataT *kernelR;
    fftw_complex *kernelF;
  };

  solverDataT data;
  dataType *arr;
  dataType *arr_in;
  fftw_complex *arrF;
  int initType;

  double boxSize[NDIMS];
  double boxSizeP[NDIMS];//padded boxsize
  
  long dim[NDIMS]; //dims of real space, no padding at all
  int dim4[NDIMS]; //dims of real space (not padding at all, int version for fftw)

  long dimP[NDIMS]; //dims of real space (padded, no extra space for fftw_complex)
  int dimP4[NDIMS]; //dims of real space (int version for fftw)

  long stride[NDIMS+1]; //stride of real space (no padding at all)
  long strideP[NDIMS+1]; //stride of real space (padding, not extra space for fftw_complex)  
  long stridePE[NDIMS+1]; //stride of real space (padding + extra space for fftw_complex)  
  long deltaP[NDIMS+1]; //stride difference due to extra fftw_complex space
  long deltaPF[NDIMS+1]; //stride difference due to extra fftw_complex space + symetry 

  long dimF[NDIMS]; //dims of fourier space (padded + symetric matrix)
  long strideF[NDIMS+1]; //stride of fourier space (padded + symetric matrix)
  std::vector<dataT*> forceDataPtr;

  long arraySize;

  double G;

public:
    
  poissonFFTSolver(double G_=1.0, int periodicity_=0, int nthreads_=num_omp_threads):
    threadsInitialized(false),nthreads(nthreads_),
    initialized(false),periodicity(periodicity_),
    initType(FFTW_MEASURE),G(G_)
  {
    memset(&data,0,sizeof(solverDataT));
    if (nthreads>1) 
      {
	fftw_init_threads();
	threadsInitialized=true;	
      }

  }

  ~poissonFFTSolver()
  {
    destroyPlans();
    if (threadsInitialized) fftw_cleanup_threads();
    if (data.kernelR!=NULL) free(data.kernelR);
    //if (data.kernelF!=NULL) free(data.kernelF);
  }

  void setG(double G_)
  {
    G=G_;
  }
  
  void setQuickInit(bool fast=false)
  {
    if (fast) initType=FFTW_ESTIMATE;
    else initType=FFTW_MEASURE;
  }

  template <typename T>
  static long inplaceArraySize(const T dim[NDIMS], int periodicity=0, bool potentialOnly=false)
  {
    int pad=~periodicity;
    long result=FFTWT::actualSize(dim,pad,true);
    return (potentialOnly)?result:NDIMS*result;//*2;//factor of 2 because fiueld is not hermitian after force comp.
  }

  template <typename T>
  long getArraySize(const T dim[NDIMS], int periodicity_=(1<<31), bool potentialOnly=false)
  {
    int pad=(periodicity_!=(1<<31))?~periodicity_:~periodicity;
    long result=FFTWT::actualSize(dim,pad,true);
    return (potentialOnly)?result:NDIMS*result;//*2;
  }

  template <typename dimT, typename dimS>
  void setup(dimT dim_[NDIMS],dimS boxSize_[NDIMS], dataT *source, dataT *dest, int periodicity_=(1<<31), bool potentialOnly=false, int nthreads_=0)
  {   

    printf("Setting-up FFT poisson solver (%s) ... ",
	   (initType==FFTW_ESTIMATE)?"fast":"slow");fflush(0);
    periodicity=(periodicity_!=(1<<31))?periodicity_:periodicity;
    nthreads=(nthreads_>0)?nthreads_:nthreads;
    int pad=~periodicity;

    arr_in=source;
    arr=dest;
    arrF=(fftw_complex *)arr;

    stride[0]=1;
    for (int i=0;i<NDIMS;i++) 
      {
	dim4[i]=(int)dim_[i];
	dim[i]=dim_[i];	
	boxSize[i]=boxSize_[i];	
	boxSizeP[i]=(periodicity&(1<<i))?boxSize_[i]:2*boxSize_[i];	
	stride[i+1]=stride[i]*dim[i];
      }

    arraySize=FFTWT::actualSize(dim,dimP,deltaP,stridePE,pad,true);
    FFTWT::fourierSize(dim,dimF,strideF,pad,true);

    strideP[0]=1;   
    for (int i=0;i<NDIMS;i++) 
      {
	dimP4[i]=(int)dimP[i];  
	strideP[i+1]=strideP[i]*dimP[i];	
      }

    //for (int i=0;i<NDIMS+1;i++)
    //printf("stride(,P,F)[%d]=[%ld %ld %ld], delta[%d]=%ld\n",i,stride[i],stridePE[i],strideF[i],i,deltaP[i]);
    
    int nth=(nthreads_<1)?nthreads:nthreads_;
    destroyPlans();
    fftw_plan_with_nthreads(nth);
    data.pf=fftw_plan_dft_r2c(NDIMS,dimP4,arr,(cplxDataT*)arr,initType);//,FFTW_MEASURE);
    data.pb=fftw_plan_dft_c2r(NDIMS,dimP4,(cplxDataT*)arr,arr,initType);//,FFTW_MEASURE);
    if (!potentialOnly)
      {
	
	data.pb3=fftw_plan_many_dft_c2r(NDIMS,dimP4,NDIMS,
					(cplxDataT*)arr,NULL,1,arraySize/2,
					arr,NULL,1,arraySize,
					initType);//,FFTW_MEASURE);
	  /*
	data.pb3=fftw_plan_many_dft(NDIMS,dimP4,NDIMS,
				    (cplxDataT*)arr,NULL,1,strideP[NDIMS], 
				    arr,NULL,1,strideP[NDIMS],
				    FFTW_BACKWARD,FFTW_MEASURE);
	  */
	data.havePB3=true;
	forceDataPtr.clear();
	for (int i=0;i<NDIMS;i++) forceDataPtr.push_back(&arr[i*strideF[NDIMS]*2]);
      }
    else forceDataPtr.clear();
    data.havePlans=true;   
    
    createKernels();
    initialized=true;
    printf("done.\n");
  }

  dataT *solve(bool potential=false)
  { 
    if (!data.havePlans)
      {
	fprintf(stderr,"ERROR in poissonFFTSolver::solve(): you must call setup(...) before solve().\n");
	exit(-1);
      }
    if ((!potential)&&(!data.havePB3))
      {
	fprintf(stderr,"ERROR in poissonFFTSolver::solve():\n");
	fprintf(stderr,"     solver is configured for potential only, cannot compute force.\n");
	exit(-1);
      }
  
    FFTWT::rearrangeToFFTW(arr_in,arr,dim,~periodicity,true);
    fftw_execute(data.pf);

    
    // multiply by green's function for potential
    if ((periodicity&PERIODIC)!=PERIODIC)
      {
#pragma omp parallel for
	for (long i=0;i<strideF[NDIMS];++i)
	  arrF[i]*=data.kernelF[i];
      }
    else
      {
#pragma omp parallel for
	for (long i=0;i<strideF[NDIMS];++i)
	  arrF[i]*=data.kernelR[i];
      }
    
    if (potential)
      {
	fftw_execute(data.pb);
	FFTWT::rearrangeFromFFTW(arr,arr,dim,~periodicity,true);
	return arr;
      }
  
#pragma omp parallel for
    for (long i=1;i<NDIMS;i++)
      {
	std::copy(forceDataPtr[0],forceDataPtr[0]+arraySize,forceDataPtr[i]);
      }
    
#pragma omp parallel for
    for (long i=0;i<NDIMS;i++)
      {
	long j;
	cplxDataT k[dimF[i]];
	long p[NDIMS];
	long &ki=p[i];
	cplxDataT *a=(cplxDataT *)forceDataPtr[i];
	double fac=-TWOPI/boxSizeP[i]; 
	
	for (j=0;j<NDIMS;j++) p[j]=0; 
	for (j=0;j<dimF[i];j++) k[j]=fac*FFTWT::indGen(j,dimP[i])*I;//F=-grad(phi)


	for (j=0;j<strideF[NDIMS];)
	  {
	    /*
	    double norm=0;
	    for (int n=0;n<NDIMS;n++)
	      {
		norm+=FFTWT::indGen(p[n],dimP[n])*FFTWT::indGen(p[n],dimP[n]);
	      }
	    norm = 2*PI*sqrt(norm)/dimP[0];
	    norm = (norm>PI)?0:(0.5*(1+cos(norm)));
	    */
	    a[j]*=k[ki];//*norm;
	    ++j;++p[0];
	    if (p[0]==dimF[0]) hlp::getNext<NDIMS>(p,dimF);
	  }
      }
    
    //fftw_execute(data.pb);
    fftw_execute(data.pb3);
    
#pragma omp parallel for
    for (long i=0;i<NDIMS;i++)
      {
	FFTWT::rearrangeFromFFTW(forceDataPtr[i],forceDataPtr[i],dim,~periodicity,true);
      }    
    
    return arr;
  }

  dataT *solvePotential()
  {
    return solve(true);
  }

  dataT *getForcePtr(long i)
  {
    if ((i<0)||(i>=forceDataPtr.size()))
      return NULL;

    return forceDataPtr[i];
  }

  const std::vector<dataT*> &getForcePtr()
  {
    return forceDataPtr;
  }

private:

  void destroyPlans()
  {
    if (data.havePlans)
      {
	fftw_destroy_plan(data.pf);
	fftw_destroy_plan(data.pb);
	if (data.havePB3) fftw_destroy_plan(data.pb3);
      }
  }

  void createKernels()
  {
    long p[NDIMS];    
    long i,j;
    double r,u;   
    double fac=1;
  
    if(data.kernelR!=NULL) free(data.kernelR);      

    if ((periodicity&PERIODIC)!=PERIODIC)
      {
	for (i=0;i<NDIMS;i++) 
	  {
	    //fac *= 1.0d/dimP[i];
	    fac *= boxSize[i]/(((double)dim[i]-1)*dimP[i]);
	    p[i]=0;  
	  }
	fac = -G*fac;
	data.kernelF=(cplxDataT*)calloc(strideF[NDIMS],sizeof(cplxDataT));
	data.kernelR=(dataT*)data.kernelF;

	//green
	for (i=0;i<stridePE[NDIMS];)
	  {
	    double d;
	    d=0;
	    for (j=0;j<NDIMS;j++)
	      {
		double q=FFTWT::norm(p[j],dimP[j],boxSizeP[j]);
		d+=q*q;
	      }
	    if (d==0) 
	      {
		data.kernelR[i]= fac/sqrt(boxSize[0]/(dim[0]-1));
		//data.kernelR[i]= fac;
	      }
	    else 
	      {
		data.kernelR[i]= fac/sqrt(d);
	      }
	
	    ++i;++p[0];
	    if (p[0]==dimP[0]) 
	      {
		i+=deltaP[1];
		hlp::getNext<NDIMS>(p,dimP);
	      }
	  }

	planT kp=fftw_plan_dft_r2c(NDIMS,dimP4,data.kernelR,data.kernelF,FFTW_ESTIMATE);	
	fftw_execute(kp);
	fftw_destroy_plan(kp);	

	//hanning
	for (i=0;i<NDIMS;i++) p[i]=0;
	for (i=0;i<strideF[NDIMS];)
	  {
	    double norm = 2.0F*PI*FFTWT::norm(p,dimP);
	    norm = (norm>PI)?0:(0.5*(1.0F+cos(norm)));

	    data.kernelF[i]*=norm;
	    ++i;++p[0];
	    if (p[0]==dimF[0]) hlp::getNext<NDIMS>(p,dimF);
	  }

      }
    else
      {
	for (i=0;i<NDIMS;i++) 
	  {
	    fac*=1.0/(dimP[i]);
	    //fac *= 1.0/(dimP[i]); // correct for volume element + FFT(FFT_inv()))
	    p[i]=0;  
	  }
	fac = -pow(fac,1.5)*4*PI*G;
	data.kernelR=(dataT*)calloc(strideF[NDIMS],sizeof(dataT)); //kernel is real here ...
	data.kernelF=(cplxDataT*)data.kernelR;
   
	for (i=0;i<strideF[NDIMS];)
	  {
	    double k2=-FFTWT::norm2(p,dimP)*PISQ*4;
	    //double k2=0;
	    //for (j=0;j<NDIMS;j++)
	    //k2+=pow(FFTWT::indgen(p[j],dimP[j]),2)*PISQ;
	    

	    //k2=-FFTWT::norm2(p,dimP,boxSizeP)*PISQ;//fac compensates for the factor in ffts
	    //k2=-FFTWT::norm2(p,dimP,boxSizeP)*PISQ;//fac compensates for the factor in ffts
	
	    if (k2!=0)
	      {
		data.kernelR[i]=fac/k2;	   
		//data.kernelF[i]=4*PI*G/k2;	    
	      }
	    else 
	      {
		data.kernelR[i]=fac;
		//data.kernelF[i]=0;	   
	      }

	    ++i;++p[0];
	    if (p[0]==dimF[0]) hlp::getNext<NDIMS>(p,dimF);
	  }
      }
  }
};

#endif
#endif
