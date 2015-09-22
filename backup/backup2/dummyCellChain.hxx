#ifndef __CELL_CHAIN_HXX__
#define __CELL_CHAIN_HXX__

#include <vector>
#include <list>

struct cellType
{
  const double R;
  const double U;
  dataT val;

  cell(dataT val_=0,double R_=0, double U_=0):
    R(R_),U(U_),val(val_)
  {
     
  }

  ~cell() {}
};

template <typename dataT, struct cellType>
class CellChains {
public:
  typedef cellType cellT;
  
  template <class cellType>
  class chainType
  {
  public:
    typedef cellType cellT;
    chainT()
    {

    }
    ~chainT() {}

    void insert_back(const cellT &c)
    {
      chain.push_back(c);
    }
    
    void push_back(const cellT &c)
    {
      chain.erase(chain.begin());
      chain.push_back(c);
    }
    
    const cellT &front() const
    {
      return chain.front();
    }

    const cellT &back() const
    {
      return chain.front();
    }

  private:
    typedef std::list<cellT> chainListT;
    typedef typename chainListT::iterator chainListItT;
    chainListT chain;
  };


  typedef chainType<cellT> chainT; 

  //typedef std::vector<cellT> chainT; 
  typedef std::vector<chainT> chainArrayT; 
  //typedef typename chainT::iterator chainItT; 
  typedef typename chainArrayT::iterator chainArrayItT; 

  dummyCellChains()
  {
    
  }

  ~dummyCellChains()
  {
    
  }

  std::pair<chainItT,long> addChain()
  {
    chainItT it = chains.insert(chains.end(),chainT());
    return std::makePair(it,chains.size());
  }

  push(long chainId, )

  chainT &chainAt(long index)
  {
    return chains[index];
  }
  
private:

  chainArrayT chains;

}

#endif

