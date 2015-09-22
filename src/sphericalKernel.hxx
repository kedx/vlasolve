#ifndef __SPHERICAL_KERNEL_HXX__
#define __SPHERICAL_KERNEL_HXX__

#include "cellChain.hxx"
#include "find_unordered_maps.hxx"
#include "typeSelect.hxx"

template <class cellType>
class sphericalKernel {
public:
  static const std::string getTag() {return "KERNEL_CELLS v0.10";}

  typedef cellType cellT;
  typedef typename cellT::dataT dataT;
  typedef cellChain<cellT> cellChainT;
  typedef std::list<cellChainT> cellChainListT;
  typedef typename cellChainListT::iterator cellChainListItT;
  typedef typename my_unordered_map<long,cellChainListItT>::type kernelCellsT;
  typedef typename kernelCellsT::iterator kernelCellsItT;
    
private:
  cellChainListT allChains;  
  std::vector<kernelCellsT> kernelBoundary;
  kernelCellsT dummyK;
    
public:
  sphericalKernel() {}
  ~sphericalKernel() {}

  kernelCellsItT begin(long j) {return ((j>=kernelBoundary.size())||(j<0))?dummyK.end():kernelBoundary[j].begin();}
  kernelCellsItT end(long j) {return ((j>=kernelBoundary.size())||(j<0))?dummyK.end():kernelBoundary[j].end();}  
  const kernelCellsItT begin(long j) const {return ((j>=kernelBoundary.size())||(j<0))?dummyK.end():kernelBoundary[j].begin();}
  const kernelCellsItT end(long j) const {return ((j>=kernelBoundary.size())||(j<0))?dummyK.end():kernelBoundary[j].end();}
  /*
  const cellT readRoll(long j, long index, bool noRoll=false)
  {
    if (noRoll)
      return kernelBoundary[j].find(index)->second->readFront();
    else
      return kernelBoundary[j].find(index)->second->readRoll();
  }
  */
  void roll(long j, long index)
  {
    kernelBoundary[j].find(index)->second->roll();
  }

  dataT readAt(long j, long index, int s=0, bool roll=false)
  {
    if (roll) 
      return kernelBoundary[j].find(index)->second->readRoll(s);
    else 
      return kernelBoundary[j].find(index)->second->readFront(s);
  }

  void multiply(double fac, int index=0)
  {
    for (cellChainListItT it=allChains.begin();it!=allChains.end();it++)
      it->multiply(fac,index);
  }

  void reserve(long j)
  {
    kernelBoundary.resize(j);
  }

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

  void reset()
  {
    kernelBoundary.clear();
    allChains.clear();
  }

  void write(FILE *f)
  {
    long size=kernelBoundary.size();
    long i;

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
    reset();
     
    myIO::checkTag(f,getTag());
    long size,size2,i,j,key;
    myIO::fread(&size,sizeof(long),1,f,swap);
    
    kernelBoundary.resize(size);
    allChains.clear();
    
    for (i=0;i<size;i++)
      {
	myIO::fread(&size2,sizeof(long),1,f,swap);
	for (j=0;j<size2;j++)
	  {
	    myIO::fread(&key,sizeof(long),1,f,swap);
	    
	    kernelBoundary[i][key] = allChains.insert(allChains.end(),cellChainT(f,swap));
	    
	  }
      }   
  }
};

#endif
