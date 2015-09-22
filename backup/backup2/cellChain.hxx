#ifndef __CELL_CHAIN_HXX__
#define __CELL_CHAIN_HXX__

#include <vector>
#include <list>
#include <map>

template <class cellType>
class cellChain {
public:
  typedef cellType cellT;
  typedef typename cellT::valT cellValT;

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

  void multiply(double fac)
  {
    for (chainListItT it=chain.begin();it!=chain.end();it++) 
      it->val *= fac;
  }

  const cellT readRoll()
  {
    cellT result=chain.front();
    chainListItT it=chain.begin();
    chainListItT cur=it++;
    
    for (;it!=chain.end();it++,cur++) cur->val=it->val;
    
    //chain.splice(chain.end(),chain,chain.begin());
    return result;
    //const cellT &result=chain.front();
    //chainListItT it=chain.end();  
    //chain.splice(chain.begin(),chain,--it);
    //return result;
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

