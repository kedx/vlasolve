/*
#ifndef __SPHERICAL_KERNEL_HXX__
#define __SPHERICAL_KERNEL_HXX__

#include "cellChain.hxx"
#include "find_unordered_maps.hxx"
#include "typeSelect.hxx"

namespace sphericalKernel {

  struct typeV {
    enum type {DELAYED=0, REFLECTIVE=1, UNDEFINED=-1};
  };
  
  typedef typename typeV::type typeT;
  
  struct select : public typeSelect<typeV> {
    select()
    {
      insert("delayed",typeV::DELAYED);
      insert("reflective",typeV::REFLECTIVE);
    }
    std::string name() {return "kernel type";}
  };

  kernelTypeT kernelType;
  kernel::typeT kernelType;

  template <class cellType>
  class cells {
  public:
    static const std::string getTag() {return "KERNEL_CELLS v0.10";}

    typedef cellType cellT;
    typedef cellChain<cellT> cellChainT;
    typedef std::list<cellChainT> cellChainListT;
    typedef typename cellChainListT::iterator cellChainListItT;
    typedef typename my_unordered_map<long,cellChainListItT>::type kernelCellsT;
    typedef typename kernelCellsT::iterator kernelCellsItT;
    
  private:
    cellChainListT allChains;  
    std::vector<kernelCellsT> kernelBoundary;
    bool initialized;

  public:
    cells():initialized(false) {}
    ~cells() {}

    kernelCellsItT begin(long j) {return kernelBoundary[j].begin();}
    kernelCellsItT end(long j) {return kernelBoundary[j].end();}  
    const kernelCellsItT begin(long j) const {return kernelBoundary[j].begin();}
    const kernelCellsItT end(long j) const {return kernelBoundary[j].end();}
    
    const cellT readRoll(long j, long index)
    {
      return kernelBoundary[j].find(i)->second->readRoll();
    }

    kernelCellsItT find(long j, long index) {return kernelBoundary[j].find(index);}
    const kernelCellsItT find(long j, long index) const {return kernelBoundary[j].find(index);}
    

    cellChainListItT newChain(long j, long index)
    {
      cellChainListItT it=allChains.insert(allChains.end(),cellChainT());

      if (j>=kernelBoundary.size()) kernelBoundary.resize(j+1);

      if (!kernelBoundary[j].insert(std::make_pair(index,it)).second)
	{
	  allChains.erase(it);
	  return allChains.end();
	}

      return it;
    }

    void write(FILE *f)
    {
      long size=kernelBoundary.size();

      myIO::writeTag(f,getTag());

      fwrite(&size,sizeof(long),1,f);
      for (i=0;i<size;i++)
	{
	  long size2=kernelBoundary[i].size();
	  fwrite(&size2,sizeof(long),1,f);
	  for (kernelCellsItT it = kernelBoundary[i].begin(); it != kernelBoundary[i].end(); it++)
	    {
	      fwrite(&it->first,sizeof(long),1,f);
	      it->second->write(f);
	    }
	}
    }

    void read(FILE *f, bool swap)
    {
      kernelBoundary.clear();
      allChains.clear();

      myIO::checkTag(f,getTag());
      long size,size2,i,j,key;
      myIO::fread(&size,sizeof(long),1,f,swap);
    
      kernelBoundary.resize(size);
      allChains.clear();
    
      for (i=0;i<size;i++)
	{
	  myIO::fread(&size2,sizeof(long),1,f,swap);
	  //printf("size2 = %ld\n",size2);
	  for (j=0;j<size2;j++)
	    {
	      myIO::fread(&key,sizeof(long),1,f,swap);
	      kernelBoundary[i][key] = allChains.insert(allChains.end(),cellChainT(f,swap));	      
	    }
	}   
    }
  };

}

#endif
*/
