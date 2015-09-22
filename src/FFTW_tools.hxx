#ifndef __FFTW_TOOLS_HXX__
#define __FFTW_TOOLS_HXX__

#ifdef HAVE_FFTW3

#include <complex.h>
#include <fftw3.h>
#include <string.h>

#include "helpers.hxx"

template <
  typename dataType,
  int NDIMS>
class FFTW_tools
{
public:
  typedef dataType dataT;
  typedef fftw_complex cplxT;
  typedef fftw_plan planT;

  template <typename T>
  static long actualSize(const T dim[NDIMS], T dimP[NDIMS], T deltaP[NDIMS+1],
			 T strideP[NDIMS+1], int padding, bool inplace)
  {
    long i;
    long delta=1;

    for (i=0;i<NDIMS+1;i++) deltaP[i]=0;    
    if (padding&(1<<0)) delta=2;
    else if (!(dim[0]&1)) delta=2;
    deltaP[1]=delta;
    /*
    if ((padding&(1<<0))&(dim[0]&1))
      delta=2;
    else 
      delta=1;
    */
    

    strideP[0]=1;
    if (padding&(1<<0)) 
      {
	dimP[0]=2*(dim[0]-1);
	strideP[1]=strideP[0]*(delta+dimP[0]);
      }
    else
      {
	dimP[0]=dim[0];
	strideP[1]=strideP[0]*(delta+dimP[0]);
      }
    
    for (i=2;i<=NDIMS;i++) 
      if (padding&(1<<(i-1))) 
	{
	  dimP[i-1]=2*(dim[i-1]-1);
	  strideP[i]=strideP[i-1]*dimP[i-1];
	}
      else
	{
	  dimP[i-1]=dim[i-1];
	  strideP[i]=strideP[i-1]*dimP[i-1];
	}
    
    return strideP[NDIMS];
  }

  template <typename T>
  static long actualSize(const T dim[NDIMS], int padding, bool inplace)
  {
    T strideP[NDIMS+1]; 
    T dimP[NDIMS+1]; 
    T deltaP[NDIMS+1]; 
    return actualSize(dim,dimP,deltaP,strideP,padding,inplace);   
  }

  template <typename T>
  static long fourierSize(const T dim[NDIMS], T dimF[NDIMS], T strideF[NDIMS+1], int padding, bool inplace)
  {
    int i;
    int fac[NDIMS];

    for (i=0;i<NDIMS;i++)
      {
	if (padding&(1<<i)) fac[i]=2;
	else fac[i]=1;
      }

    if (inplace) 
      {
	if (padding&(1<<i))
	  dimF[0]=(2*(dim[0]-1))/2+1;
	else
	  dimF[0]=fac[0]*dim[0]/2+1;
      }
    else
      {
	if (padding&(1<<i))
	  dimF[0]=2*(dim[0]-1);
	else
	  dimF[0]=dim[0];
      }
      
    strideF[0]=1;strideF[1]=dimF[0];
    for (i=1;i<NDIMS;i++) 
      {
	if (padding&(1<<i))
	  {
	    dimF[i]=2*(dim[i]-1);
	    strideF[i+1]=strideF[i]*dimF[i];
	  }
	else
	  {
	    dimF[i]=dim[i];
	    strideF[i+1]=strideF[i]*dimF[i];
	  }
      }

    return strideF[NDIMS];
  }

  template <typename T>
  static dataT *rearrange(dataT *data, dataT *dest,const T dim[NDIMS], int padding, bool inplace, bool toFFTW=true)
  {
    dataT *newData;
    T strideP[NDIMS+1];
    T dimP[NDIMS];
    T deltaP[NDIMS];
    T x[NDIMS];
    long nval=1;
    long i,j,k,m;

    for (i=0;i<NDIMS;i++) nval*=dim[i];   
    
    long size = actualSize(dim,dimP,deltaP,strideP,padding,inplace);

    if (dest==NULL)
      {
	if (toFFTW) 
	  dest=(dataT *) malloc (size*sizeof(dataT));	
	else
	  dest=(dataT *) malloc (nval*sizeof(dataT));	
      }
    newData=dest;

    if (toFFTW)
      {
	
	/*
	for (i=0;i<NDIMS;i++) x[i]=dim[i]-1;
	// reverse order so that data can be rearranged inplace
	for (i=nval-1;i>=0;i--)
	  {	
	    for (j=0,m=0;j<NDIMS;j++) m+=x[j]*strideP[j];

	    newData[m]=(double)data[i];
	    if (i!=m) data[i]=0;
		
	    x[0]--;
	    if (x[0]<0)
	      {
		j=1;

		do
		  {
		    x[j-1]=dim[j-1]-1;
		    x[j++]--; 
		  } while ((j<NDIMS)&&(x[j-1]<0));
	      }	
	  }
	*/
	long sz[NDIMS];
	sz[0]=(deltaP[1]+(dimP[0]-dim[0])*strideP[0]);
	for (i=1;i<NDIMS;i++) sz[i]=(dimP[i]-dim[i])*strideP[i];
	for (i=0;i<NDIMS;i++) x[i]=0;

	j=nval-1;
	for (i=strideP[NDIMS]-1;i>=0;)
	  {
	    //i-=(deltaP[1]+(dimP[0]-dim[0])*strideP[0]);
	    //for (int l=0;l<(deltaP[1]+(dimP[0]-dim[0])*strideP[0]);l++) newData[i--]=0;	
	    //i-=sz[0];
	    //memset(&newData[i],0,sizeof(dataT)*sz[0]);
	    
	    for (k=0;k<NDIMS;k++)
	      if (x[k]==0) 
		{
		  i-=sz[k];
		  memset(&newData[i],0,sizeof(dataT)*sz[k]);
		  //i-=(dimP[k]-dim[k])*strideP[k];
		  //for (int l=0;l<(dimP[k]-dim[k])*strideP[k];l++) newData[i--]=0;	
		}
	      else break;
	    	  
	    for (k=0;k<dim[0];k++) newData[i--]=data[j--];
	    hlp::getNext<NDIMS>(x,dim);
	  }
	//printf("i=%ld j=%ld\n",i,j);
	
      }
    else
      {
	const long Ncpy=dim[0];
	const size_t mmvSz=Ncpy*sizeof(dataT);
	const long Ntot=nval;
	T deltaP[NDIMS+1];
	
	for (i=0;i<NDIMS;i++) x[i]=0;
	 
	deltaP[0]=strideP[1];
	for (i=1;i<=NDIMS;i++) deltaP[i]=strideP[i]*(dimP[i]-dim[i]);

	for (i=0,j=0;i<Ntot;i+=Ncpy)
	  {	   
	    memmove(&newData[i],&data[j],mmvSz);

	    j+=deltaP[0];
	    hlp::getNext<NDIMS>(x,j,dim,deltaP);
	    }
      }

    return newData;
  }

  template <typename T>
  static dataT *rearrangeToFFTW(dataT *src,dataT *dst, const T dim[NDIMS], int padding, bool inplace)
  {
    return rearrange(src,dst,dim,padding,inplace,true);
  }

  template <typename T>
  static dataT *rearrangeFromFFTW(dataT *src,dataT *dst, const T dim[NDIMS], int padding, bool inplace)
  {
    return rearrange(src,dst,dim,padding,inplace,false);
  }
  /*
  template <typename T>
  static dataT *unHermitian(cplxT *data, const T dim[NDIMS], const T dimF[NDIMS])
  {
    long nval=1;
    long nvalF=1;
    long i,j,iF;
    T x[NDIMS];    
        
    for (i=0;i<NDIMS;i++) 
      {
	nvalF*=dimF[i];
	nval*=dim[i];
	x[i]=0;	
      }

    i=0;iF=nvalF;
    for (j=nval-1;j>=0;)
      {
	if (x[0]<dimF[0])
	  data[j]=conj(data[i]);
	else
	  data[j]=data[iF];

	if (x[0]>=(dim[0]-dimF[0])) iF--;

	x[0]++;i++;j--;
	if (x[0]==dim[0]) hlp::getNext<NDIMS>(x,dim);
      }
      
  }
  */
  
  template <typename T>
  static T indGen(T k,T dim)
  {
    return (k<=(dim*0.5))?k:(k-dim);  
  }

  template <typename T>
  static T indGen(T k[NDIMS],T dim[NDIMS], T out[NDIMS])
  {
    for (int i=0;i<NDIMS;i++)
      out[i]=(k[i]<=(dim[i]*0.5))?k[i]:(k[i]-dim[i]);  
  }

  template <typename T>
  static T indGen(T *k,T *dim, T *out, int ndims)
  {
    for (int i=0;i<ndims;i++)
      out[i]=(k[i]<=(dim[i]*0.5))?k[i]:(k[i]-dim[i]);  
  }

  template <typename T,typename T2,typename T3>
  static double norm(const T k, const T2 delta, const T3 size)
  {
    double a=(k<=(delta*0.5))?k:(k-delta);
    return (a*size)/delta;
  }
 

  template <typename T,typename T2>
  static double norm(const T k[NDIMS], const T2 delta[NDIMS])
  {
    double a[NDIMS];
    double result=0;

    for (long i=0;i<NDIMS;i++)
      {
	a[i]=(k[i]<=(delta[i]*0.5))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i])/(delta[i]*delta[i]);
      }

    return sqrt(result);
  }

  template <typename T,typename T2>
  static double norm(const T *k, const T2 *delta, int ndims)
  {
    double a[ndims];
    double result=0;

    for (long i=0;i<ndims;i++)
      {
	a[i]=(k[i]<=(delta[i]*0.5))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i])/(delta[i]*delta[i]);
      }

    return sqrt(result);
  }

  template <typename T,typename T2,typename T3>
  static double norm(const T k[NDIMS], const T2 delta[NDIMS], const T3 size[NDIMS])
  {
    double a[NDIMS];
    double result=0;

    for (long i=0;i<NDIMS;i++)
      {
	a[i]=(k[i]<=(delta[i]*0.5))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i]*size[i]*size[i])/(delta[i]*delta[i]);
      }

    return sqrt(result);
  }

  template <typename T,typename T2,typename T3>
  static double norm(const T *k, const T2 *delta, const T3 *size,int ndims)
  {
    double a[ndims];
    double result=0;

    for (long i=0;i<ndims;i++)
      {
	a[i]=(k[i]<=(delta[i]*0.5))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i]*size[i]*size[i])/(delta[i]*delta[i]);
      }

    return sqrt(result);
  }

  template <typename T,typename T2>
  static double norm2(const T k[NDIMS], const T2 delta[NDIMS])
  {
    double a[NDIMS];
    double result=0;

    for (long i=0;i<NDIMS;i++)
      {
	a[i]=(k[i]<=(delta[i]*0.5))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i])/(delta[i]*delta[i]);
      }

    return result;
  }

  template <typename T,typename T2>
  static double norm2(const T *k, const T2 *delta, int ndims)
  {
    double a[ndims];
    double result=0;

    for (long i=0;i<ndims;i++)
      {
	a[i]=(k[i]<=(delta[i]*0.5))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i])/(delta[i]*delta[i]);
      }

    return result;
  }

  template <typename T,typename T2,typename T3>
  static double norm2(const T k[NDIMS], const T2 delta[NDIMS], const T3 size[NDIMS])
  {
    double a[NDIMS];
    double result=0;

    for (long i=0;i<NDIMS;i++)
      {
	a[i]=(k[i]<=(delta[i]*0.5))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i]*size[i]*size[i])/(delta[i]*delta[i]);
      }

    return result;
  }

  template <typename T,typename T2,typename T3>
  static double norm2(const T *k, const T2 *delta, const T3 *size, int ndims)
  {
    double a[ndims];
    double result=0;

    for (long i=0;i<ndims;i++)
      {
	a[i]=(k[i]<=(delta[i]*0.5))?k[i]:(k[i]-delta[i]);
	result+=(a[i]*a[i]*size[i]*size[i])/(delta[i]*delta[i]);
      }

    return result;
  }
};

#endif
#endif
