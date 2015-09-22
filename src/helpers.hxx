#ifndef __HELPERS_HXX__
#define __HELPERS_HXX__

#include <algorithm>

namespace hlp {

  template<int ND,typename T1,typename T2>
  inline void getNext(T1 w[], T2 &k, const T1 dm[], const T2 delta[])
  {
    w[0]=0;++w[1];
    if (w[1]==dm[1])
      {
	k+=delta[1];	
	if (ND==2) return;
	w[1]=0;++w[2];
	if (w[2]==dm[2])
	  {
	    k+=delta[2];
	    if (ND==3) return;
	    w[2]=0;++w[3];
	    if (w[3]==dm[3])
	      {
		k+=delta[3];
		if (ND==4) return;
		w[3]=0;++w[4];
		if (w[4]==dm[4])
		  {
		    k+=delta[4];
		    if (ND==5) return;
		    w[4]=0;++w[5];
		    if (w[5]==dm[5])
		      {
			k+=delta[5];
			if (ND==6) return;
			w[5]=0;++w[6];
		      
			int ct=6;
			for (;;)
			  {
			    if (w[ct]==dm[ct])
			      {
				k+=delta[ct];
				if (ND==ct+1) return;
				w[ct]=0;++w[++ct];
			      } 
			    else return;
			  } 
		      
		      }
		  }
	      }
	  }
      }
  }

  template<int ND,typename T>
  inline void getNext(T w[], const T dm[])
  {   
    w[0]=0;++w[1];
    if (w[1]==dm[1])
      {
	if (ND==2) return;
	w[1]=0;++w[2];
	if (w[2]==dm[2])
	  {
	    if (ND==3) return;
	    w[2]=0;++w[3];
	    if (w[3]==dm[3])
	      {
		if (ND==4) return;
		w[3]=0;++w[4];
		if (w[4]==dm[4])
		  {
		    if (ND==5) return;
		    w[4]=0;++w[5];
		    if (w[5]==dm[5])
		      {
			if (ND==6) return;
			w[5]=0;++w[6];
		      
			int ct=6;
			for (;;)
			  {
			    if (w[ct]==dm[ct])
			      {
				if (ND==ct+1) return;
				w[ct]=0;++w[++ct];
			      } 
			    else return;
			  } 
		      
		      }
		  }
	      }
	  }
      }
  }

  template<>
  inline void getNext<1>(int w[], long &k, const int dm[], const long strd[])
  {
    return;
  }

  template <int A, int B> 
  struct isEqual
  {
  public:
    static const bool result=(A==B);
  };

  template <class container>
  typename container::const_iterator findValue(const container &C,const typename container::value_type &v, double tolerance=1.E-10)
  {
    typedef typename container::const_iterator const_iteratorT;
    const_iteratorT pos_it=std::lower_bound(C.begin(),C.end(),v);

    if ((*pos_it)!=v)
      {
	if (fabs(((*pos_it)-v)/((*pos_it)+v))>tolerance)
	  {
	    pos_it--;
	    if ((*pos_it)!=v)
	      {
		if (fabs(((*pos_it)-v)/((*pos_it)+v))>tolerance) return C.end();		
	      }
	  }
      }
    return pos_it;
  }

   //**** CONSTANT_VALUE 

  template<long N>
  struct ConstantValue
  {
    //static const long value = N; 
    enum { value = N };
  };

  template<typename T>
  struct ConstantType
  {
    typedef T Type;
  };
  

  // INTEGER power
  template <long A, long B, long value>
  struct IntPowerHelper
  {
    typedef typename IntPowerHelper<A,B-1,value*A>::Result Result;
  };

  template <long A, long value>
  struct IntPowerHelper<A,0,value>
  {
    typedef ConstantValue<value> Result;
  };

  template <long A, long B>
  struct IntPower
  {
    typedef typename IntPowerHelper<A,B,1>::Result Result;
    enum { value = Result::value };
  };

  // Get the sign of any type
  template <typename T> T sgn(T val) {
    return (T(0) < val) - (val < T(0));
  }


  /*
  template <int ND, class T, class T2> 
  void get_CT(T *res,const T *t1,const  T* t2,const  T2* w)
  {
    int i;
    for (i=0;i<ND;i++) res[i]=t1[w[i]];
    for (i=ND;i<ND+ND;i++) res[i]=t2[w[i]];
  }

  template <>
  void get_CT<1,double,int>(double *res,const double *t1,const double* t2,const int* w)
  {
    res[0]=t1[w[0]];
    res[1]=t2[w[1]];
  }

  template <>
  void get_CT<2,double,int>(double *res,const double *t1,const double* t2,const int* w)
  {
    res[0]=t1[w[0]];res[1]=t1[w[1]];
    res[2]=t2[w[2]];res[3]=t2[w[3]];
  }

  template <>
  void get_CT<3,double,int>(double *res,const double *t1,const double* t2,const int* w)
  {
    res[0]=t1[w[0]];res[1]=t1[w[1]];res[2]=t1[w[2]];
    res[3]=t2[w[3]];res[4]=t2[w[4]];res[5]=t2[w[5]];
  }
  */

};

#endif
