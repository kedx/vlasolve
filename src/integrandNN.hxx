#ifndef __INTEGRAND_NN_HXX__
#define __INTEGRAND_NN_HXX__

namespace intNN {
  
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
	  return it.cellAvg();
	case 1:
	  return it.cellAvg(val);
	}
      
      double result=0;
      for (int i=0;i<nval;i++) {
	result+=it.cellAvg(which[i]);
      }
      return result;
    }

    template <class iterator>
    double it_sliceAvg(const iterator &it, bool ortho) const
    {     
      switch (nval)
	{
	case 0:
	  return it.slice_cellAvg(0,ortho);
	case 1:
	  return it.slice_cellAvg(val,ortho);
	}
      
      double result=0;
      for (int i=0;i<nval;i++) {
	result+=it.slice_cellAvg(which[i],ortho);
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

  template <int ND>
  struct f : public integrand
  {
    typedef integrand bT;
    f():bT(){}
    f(int start,int count):bT(start,count){}
    template <class container>
    f(const container &c):bT(c){}

    template <class iterator>
    double operator()(double p[2*ND], const iterator &it) const
    {
      return it_avg(it);
    }

    template <class iterator>
    double operator()(double p[2*ND], const iterator &it, bool ortho) const
    {
      return it_sliceAvg(it,ortho);
    }
  };

  template <int ND>
  struct kineticEnergy : public integrand
  {
    typedef integrand bT;
    kineticEnergy():bT(){}
    kineticEnergy(int start,int count):bT(start,count){}
    template <class container>
    kineticEnergy(const container &c):bT(c){}

    template <class iterator>
    double operator()(double p[2*ND], const iterator &it) const
    {
      double v=p[ND]*p[ND];
      for (int i=ND+1;i<2*ND;i++) v+=p[i]*p[i];
      return 0.5*it_avg(it)*v;
      //printf("IntegrandNN::kineticEnergy not implemented!\n");exit(-1);
    }
  };
  
  template <int ND>
  struct entropy : public integrand
  {
    typedef integrand bT;
    entropy():bT(){}
    entropy(int start,int count):bT(start,count){}
    template <class container>
    entropy(const container &c):bT(c){}

    template <class iterator>
    double operator()(double p[2*ND], const iterator &it) const
    {
      double val=it_avg(it);
      return (val<=0)?0:-val*log(val);
    }
  };

  template <int ND>
  struct L1Norm : public integrand
  {
    typedef integrand bT;
    L1Norm():bT(){}
    L1Norm(int start,int count):bT(start,count){}
    template <class container>
    L1Norm(const container &c):bT(c){}

    template <class iterator>
    double operator()(double p[2*ND], const iterator &it) const
    {
      return fabs(it_avg(it));
    }
  };

  template <int ND>
  struct L2Norm : public integrand
  {
    typedef integrand bT;
    L2Norm():bT(){}
    L2Norm(int start,int count):bT(start,count){}
    template <class container>
    L2Norm(const container &c):bT(c){}
    template <class iterator>

    double operator()(double p[2*ND], const iterator &it) const
    {
      double val=it_avg(it);
      return val*val;
    }
  };
};


#endif
