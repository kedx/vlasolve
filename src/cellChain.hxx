#ifndef __CELL_CHAIN_HXX__
#define __CELL_CHAIN_HXX__

#include <vector>
#include <list>
#include <map>

template <class cellType>
class cellChain {
public:
  typedef cellType cellT;
  typedef typename cellT::dataT cellDataT;

  cellChain()
  {
    
  }

  ~cellChain() {}

  cellChain(FILE *f, bool swap)
  {
    read(f,swap);
  }
  
  void push(const cellT &c)
  {
    chain.push_back(c);
  }

  void print(const std::string &before="", const std::string &after="")
  {
    printf("%s",before.c_str());
    for (chainListItT it=chain.begin();it!=chain.end();it++) 
      it->print();
    printf("%s",after.c_str());
  }

  void multiply(double fac, int index=0)
  {
    
    for (chainListItT it=chain.begin();it!=chain.end();it++)
      {
	//printf("%g->%g\n",it->val[index],it->val[index]*fac);
	it->val[index]*=fac;
      }
  }

  void roll()
  {
    chainListItT it=chain.begin();
    chainListItT cur=it++;
    
    for (;it!=chain.end();it++,cur++) 
      cur->val.swap(it->val);
  }
  
  const cellDataT readRoll(int s=0)
  {
    cellDataT result=chain.front()[s];
    roll();
    return result;
  }
  
  const cellDataT readFront(int s=0)
  {
    return chain.front()[s];    
  }
   
  cellT &back()
  {
    return chain.back();
  }

  cellT &front()
  {
    return chain.back();
  }
  
  void read(FILE *f, bool swap)
  {
    long size;
    myIO::fread(&size,sizeof(long),1,f,swap);
    chain.resize(size);
   
    for (chainListItT it=chain.begin();it!=chain.end();it++)
      it->read(f,swap);
   
  }

  void write(FILE *f)
  {
    long size=chain.size();
    fwrite(&size,sizeof(long),1,f);
    for (chainListItT it=chain.begin();it!=chain.end();it++)
      it->write(f);
  }

private:
  typedef std::list<cellT> chainListT;
  typedef typename chainListT::iterator chainListItT;
  chainListT chain;
};


#endif

