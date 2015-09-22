#ifndef __INTEGRAND_21_HXX__
#define __INTEGRAND_21_HXX__

namespace int21 {
  //template <class iterator>
  class integrand
  {
  protected:
    typedef integrand myT;
    const int nval;
    const int val;
    const std::vector<int> which;

    static std::vector<int> range(const int start,const int count)
    {
      std::vector<int> result;
      for (int i=start;i<start+count;i++) result.push_back(i);
      return result;
    }

    template <class iterator>
    double it_avg(const iterator &it) const
    {
     
      switch (nval)
	{
	case 0:
	  return it.get_avg();
	case 1:
	  return it.get_avg(val);
	}
      
      double result=0;
      for (int i=0;i<nval;i++) {
	result+=it.get_avg(which[i]);
	//printf("+ avg[%d]=%g\n",which[i],it.get_avg(which[i]));
      }
      return result;
    }
      
  public:
    integrand()
      :nval(0),
       val(0)
    {
    }

    integrand(int start, int count)
      :nval(((count==1)&&(start==0))?0:count),
       val(start),
       which(myT::range(start,count))
    {
    }

    template <class container>
    integrand(const container &c)
      :nval(c.size()),
       val(*c.begin()),
       which(c.begin(),c.end())
    {
    }
    
  };

  //template <class iterator>
  struct f : public integrand
  {
    typedef integrand bT;
    f():bT(){}
    f(int start,int count):bT(start,count){}
    template <class container>
    f(const container &c):bT(c){}

    template <class iterator>
    double operator()(double, double, double, const iterator &it) const
    {
      return it_avg(it);
    }
  };

  //template <class iterator>
  struct kineticEnergy : public integrand
  {
    typedef integrand bT;
    kineticEnergy():bT(){}
    kineticEnergy(int start,int count):bT(start,count){}
    template <class container>
    kineticEnergy(const container &c):bT(c){}

    template <class iterator>
    double operator()(double r, double u, double j, const iterator &it) const
    {
      return 0.5*it_avg(it)*(u*u+j*j/(r*r));
    }
  };
  
  //template <class iterator>
  struct entropy : public integrand
  {
    typedef integrand bT;
    entropy():bT(){}
    entropy(int start,int count):bT(start,count){}
    template <class container>
    entropy(const container &c):bT(c){}

    template <class iterator>
    double operator()(double, double, double, const iterator &it) const
    {
      double val=it_avg(it);
      return (val<=0)?0:-val*log(val);
    }
  };

  //template <class iterator>
  struct L1Norm : public integrand
  {
    typedef integrand bT;
    L1Norm():bT(){}
    L1Norm(int start,int count):bT(start,count){}
    template <class container>
    L1Norm(const container &c):bT(c){}

    template <class iterator>
    double operator()(double, double, double, const iterator &it) const
    {
      return fabs(it_avg(it));
    }
  };

  //template <class iterator>
  struct L2Norm : public integrand
  {
    typedef integrand bT;
    L2Norm():bT(){}
    L2Norm(int start,int count):bT(start,count){}
    template <class container>
    L2Norm(const container &c):bT(c){}
    template <class iterator>

    double operator()(double, double, double, const iterator &it) const
    {
      double val=it_avg(it);
      return val*val;
    }
  };
};


#endif
