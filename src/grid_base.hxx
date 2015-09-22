#ifndef __GRID_BASE_HXX__
#define __GRID_BASE_HXX__

#include <map>
#include <vector>
#include <stdio.h>
#include <malloc.h>

#include "NDfield.hxx"
#include "subgridIterator.hxx"
#include "gridNav.hxx"
#include "scale.hxx"
#include "dimTraits.hxx"
#include "myIO.hxx"
#include "valLocation_traits.hxx"
#include "helpers.hxx"

template <
  class dimensionTraits,
  typename dataType,
  class valLocationTraits
  >
class gridBase {
public:
  //template template <class dimensionTraits2,typename dataType2,class valLocationTraits2> 
  //friend class gridBase;
  
  typedef dataType value_type;
  typedef gridBase<dimensionTraits,dataType,valLocationTraits> myT;  

  friend class subgridIterator<myT>;
  typedef subgridIterator<myT> iterator;
  typedef field_iterator<iterator> fieldIterator;
  typedef oneField_iterator<iterator> oneFieldIterator;

  static std::string getTag() {return "GRID_BASE v0.12";}
  
  typedef scale<value_type> scaleT;
  typedef dimensionTraits dimTr;

  typedef typename gridNav::dirT dirT;

  typedef typename scaleT::scaleTypeT scaleTypeT;
  typedef typename scaleT::scaleTypeV scaleTypeV;
  typedef typename scaleT::valLocationT valLocationT;
  typedef typename scaleT::valLocationV valLocationV;

  static const int P_DIMS= dimensionTraits::P_DIMS;
  static const int U_DIMS= dimensionTraits::U_DIMS;
  static const int J_DIMS= U_DIMS-P_DIMS;
  static const int DIMS = U_DIMS+P_DIMS;

  typedef valLocationTraits valLocationTr;
  static const valLocationT valLocation_P = valLocationTraits::VP;   
  static const valLocationT valLocation_U = valLocationTraits::VU;
  static const valLocationT valLocation_J = valLocationTraits::VJ;
  

  struct paramT {
    static const int P_DIMS= myT::P_DIMS;
    static const int U_DIMS= myT::U_DIMS;
    static const int DIMS= myT::DIMS;

    double Pmin[P_DIMS];
    double Pmax[P_DIMS];
    double Umin[U_DIMS];
    double Umax[U_DIMS];
    int Pres[P_DIMS];
    int Ures[U_DIMS];
    scaleTypeT Pscale[P_DIMS];
    scaleTypeT Uscale[U_DIMS];

    int lowPmargin[P_DIMS];
    int highPmargin[P_DIMS];
    int lowUmargin[U_DIMS];
    int highUmargin[U_DIMS];
    
    int nFields;

    int haveParentGrid;
    int parent_P_DIMS;
    int parent_U_DIMS;  
    int parent_which[DIMS];
    long minStorageSize;
  
    template <class subParamT>
    void getSubGridParams(subParamT &sp, int which[subParamT::DIMS], long minStorage=0) const
    {
      int i;
      sp.nFields=nFields;      
      for (i=0;i<subParamT::P_DIMS;i++)
	{
	  bool useU=(which[i])>=P_DIMS;
	  int id=(useU)?(which[i]-P_DIMS):which[i];
	  
	  sp.Pmin[i]=(useU)?Umin[id]:Pmin[id];
	  sp.Pmax[i]=(useU)?Umax[id]:Pmax[id];
	  sp.Pres[i]=(useU)?Ures[id]:Pres[id];
	  sp.Pscale[i]=(useU)?Uscale[id]:Pscale[id];
	  sp.lowPmargin[i]=(useU)?lowUmargin[id]:lowPmargin[id];
	  sp.highPmargin[i]=(useU)?highUmargin[id]:highPmargin[id];
	}
      for (i=0;i<subParamT::U_DIMS;i++)
	{
	  bool useU=(which[i+subParamT::P_DIMS])>=P_DIMS;
	  int id=(useU)?(which[i+subParamT::P_DIMS]-P_DIMS):which[i+subParamT::P_DIMS];
	  
	  sp.Umin[i]=(useU)?Umin[id]:Pmin[id];
	  sp.Umax[i]=(useU)?Umax[id]:Pmax[id];
	  sp.Ures[i]=(useU)?Ures[id]:Pres[id];
	  sp.Uscale[i]=(useU)?Uscale[id]:Pscale[id];
	  sp.lowUmargin[i]=(useU)?lowUmargin[id]:lowPmargin[id];
	  sp.highUmargin[i]=(useU)?highUmargin[id]:highPmargin[id];
	}

      sp.haveParentGrid=1;
      sp.parent_P_DIMS=P_DIMS;
      sp.parent_U_DIMS=U_DIMS;
      sp.minStorageSize=minStorage;
      for (i=0;i<DIMS;i++) if (i<subParamT::DIMS) sp.parent_which[i]=which[i];
    }
    
    void write(FILE *f)
    {
      long i;
      int t;

      fwrite(Pmin,sizeof(double),P_DIMS,f);
      fwrite(Pmax,sizeof(double),U_DIMS,f);
      fwrite(Umin,sizeof(double),P_DIMS,f);
      fwrite(Umax,sizeof(double),U_DIMS,f);
      fwrite(Pres,sizeof(int),P_DIMS,f);
      fwrite(Ures,sizeof(int),U_DIMS,f);

      for (i=0;i<P_DIMS;i++) {t=(int)Pscale[i];fwrite(&t,sizeof(int),1,f);}
      for (i=0;i<U_DIMS;i++) {t=(int)Uscale[i];fwrite(&t,sizeof(int),1,f);}

      fwrite(lowPmargin,sizeof(int),P_DIMS,f);
      fwrite(highPmargin,sizeof(int),P_DIMS,f);
      fwrite(lowUmargin,sizeof(int),U_DIMS,f);
      fwrite(highUmargin,sizeof(int),U_DIMS,f);

      fwrite(&nFields,sizeof(int),1,f);

      fwrite(&haveParentGrid,sizeof(int),1,f);
      fwrite(&parent_P_DIMS,sizeof(int),1,f);
      fwrite(&parent_U_DIMS,sizeof(int),1,f);
      fwrite(parent_which,sizeof(int),DIMS,f);
      fwrite(&minStorageSize,sizeof(long),1,f);

      std::vector<char> spare(256*8 - sizeof(int)*(3+DIMS)-sizeof(long),0);
      fwrite(&spare[0],sizeof(char),spare.size(),f);
    }

    void read(FILE *f,bool swap)
    {
      long i;
      myIO::fread(Pmin,sizeof(double),P_DIMS,f,swap);
      myIO::fread(Pmax,sizeof(double),U_DIMS,f,swap);
      myIO::fread(Umin,sizeof(double),P_DIMS,f,swap);
      myIO::fread(Umax,sizeof(double),U_DIMS,f,swap);
      myIO::fread(Pres,sizeof(int),P_DIMS,f,swap);
      myIO::fread(Ures,sizeof(int),U_DIMS,f,swap);

      int tmp[P_DIMS+U_DIMS];
      myIO::fread(tmp,sizeof(int),P_DIMS+U_DIMS,f,swap);
      for (i=0;i<P_DIMS;i++) Pscale[i]=(scaleTypeT)tmp[i];
      for (i=P_DIMS;i<P_DIMS+U_DIMS;i++) Uscale[i-P_DIMS]=(scaleTypeT)tmp[i];

      myIO::fread(lowPmargin,sizeof(int),P_DIMS,f,swap);
      myIO::fread(highPmargin,sizeof(int),P_DIMS,f,swap);
      myIO::fread(lowUmargin,sizeof(int),U_DIMS,f,swap);
      myIO::fread(highUmargin,sizeof(int),U_DIMS,f,swap);   

      myIO::fread(&nFields,sizeof(int),1,f,swap);

      myIO::fread(&haveParentGrid,sizeof(int),1,f,swap);
      myIO::fread(&parent_P_DIMS,sizeof(int),1,f,swap);
      myIO::fread(&parent_U_DIMS,sizeof(int),1,f,swap);
      myIO::fread(parent_which,sizeof(int),DIMS,f,swap);
      myIO::fread(&minStorageSize,sizeof(long),1,f,swap);

      std::vector<char> spare(256*8- sizeof(int)*(3+DIMS)-sizeof(long),0);
      myIO::fread(&spare[0],sizeof(char),spare.size(),f,swap);
    }
  };

  static const std::string parserCategory() {return "grid";}

  
protected:
  paramT p; 
  long n_Vert;
  long n_Cell;

  long n_Val;
  long n_Fields;

  long arr_dims[DIMS];
  double boxSize[DIMS];

  std::vector<double> P_Vert[P_DIMS];
  std::vector<double> U_Vert[U_DIMS];

  std::vector<double> P_Cell[P_DIMS];
  std::vector<double> U_Cell[U_DIMS];

  std::vector<double> P_Val[P_DIMS];
  std::vector<double> U_Val[U_DIMS];

  std::vector<double> P_CellDelta[P_DIMS];
  std::vector<double> U_CellDelta[U_DIMS];

  std::vector<double> P_VertDelta[P_DIMS];
  std::vector<double> U_VertDelta[U_DIMS];

  long stride_Vert[P_DIMS+U_DIMS+1];
  long stride_Cell[P_DIMS+U_DIMS+1];
  long stride_Val[P_DIMS+U_DIMS+1];

  std::vector<long> integrationPoints[1<<DIMS];
  
  value_type *arr;
  bool ownArr;
  
private:
  bool initialized;
  typedef std::map<std::string,double> valueMapT;
  typedef typename valueMapT::iterator valueMapItT;
  valueMapT valueMap;

public:
  bool registerValue(const char *name, double value, bool replace=true)
  {
    std::string pname(name);
    std::pair<valueMapItT,bool> res=
      valueMap.insert(std::make_pair(pname,value));

    if ((!res.second)&&(replace)) {
      valueMap.erase(res.first);
      res=valueMap.insert(std::make_pair(pname,value));
      return false;
    }
    
    return true;
  }

  bool getRegisteredValue(const char *name, double &value)
  {
    std::string pname(name);
    valueMapItT res=valueMap.find(pname);
    if (res==valueMap.end()) return false;
    value=res->second;
    return true;
  }
  
private:
  void init(value_type *myArr=NULL)
  {    
    int i,j;
    if (initialized) fprintf(stderr,"WARNING : gridBase should be initialized once only!\n");
    n_Vert=1; 
    n_Cell=1; 
    n_Val=1;
    n_Fields=p.nFields;
    //printf("n_F=%ld\n",n_Fields);
    stride_Vert[0]=1;//1*n_Fields;
    stride_Cell[0]=1;//1*n_Fields;
    stride_Val[0]=1;//1*n_Fields;

    for (i=0;i<P_DIMS;i++) {
      
      P_Vert[i]=scaleT::genScale(p.Pmin[i],p.Pmax[i],p.Pres[i],p.Pscale[i],valLocationV::VERTEX);
      P_Cell[i]=scaleT::genScale(p.Pmin[i],p.Pmax[i],p.Pres[i],p.Pscale[i],valLocationV::CELL);
      P_Val[i]=scaleT::genScale(p.Pmin[i],p.Pmax[i],p.Pres[i],p.Pscale[i],valLocation_P);
      
      P_VertDelta[i]=scaleT::genScaleDelta(p.Pmin[i],p.Pmax[i],p.Pres[i],p.Pscale[i],valLocationV::VERTEX);
      P_CellDelta[i]=scaleT::genScaleDelta(p.Pmin[i],p.Pmax[i],p.Pres[i],p.Pscale[i],valLocationV::CELL);
      P_CellDelta[i].push_back(0);

      n_Vert*=P_Vert[i].size();
      n_Cell*=P_Cell[i].size();
      n_Val*=P_Val[i].size();

      stride_Cell[i+1]=stride_Cell[i]*P_Cell[i].size();
      stride_Vert[i+1]=stride_Vert[i]*P_Vert[i].size();
      stride_Val[i+1]=stride_Val[i]*P_Val[i].size();

    }
    for (i=0;i<U_DIMS;i++) {
      
      U_Vert[i]=scaleT::genScale(p.Umin[i],p.Umax[i],p.Ures[i],p.Uscale[i],valLocationV::VERTEX);
      U_Cell[i]=scaleT::genScale(p.Umin[i],p.Umax[i],p.Ures[i],p.Uscale[i],valLocationV::CELL);

      U_VertDelta[i]=scaleT::genScaleDelta(p.Umin[i],p.Umax[i],p.Ures[i],p.Uscale[i],valLocationV::VERTEX);
      U_CellDelta[i]=scaleT::genScaleDelta(p.Umin[i],p.Umax[i],p.Ures[i],p.Uscale[i],valLocationV::CELL);
      U_CellDelta[i].push_back(0);
      
      if (i==(U_DIMS-J_DIMS)) 
	U_Val[i]=scaleT::genScale(p.Umin[i],p.Umax[i],p.Ures[i],p.Uscale[i],valLocation_J);
      else 
	U_Val[i]=scaleT::genScale(p.Umin[i],p.Umax[i],p.Ures[i],p.Uscale[i],valLocation_U);

      n_Vert*=U_Vert[i].size();
      n_Cell*=U_Cell[i].size();
      n_Val*=U_Val[i].size();

      stride_Cell[i+P_DIMS+1]=stride_Cell[i+P_DIMS]*U_Cell[i].size();
      stride_Vert[i+P_DIMS+1]=stride_Vert[i+P_DIMS]*U_Vert[i].size();
      stride_Val[i+P_DIMS+1]=stride_Val[i+P_DIMS]*U_Val[i].size();
    }
    
    long storageSize = std::max(p.minStorageSize,n_Val*n_Fields);
   
    

    if ((arr!=NULL)&&(ownArr)) free(arr);
    
   
    if (myArr!=NULL)
      {
	ownArr=false;
	arr=myArr;
      }
    else
      {
	ownArr=true;
	int res=posix_memalign((void**)&arr,32,storageSize*sizeof(value_type));
	memset(arr,0,storageSize*sizeof(value_type));
	//arr=(value_type*)calloc(storageSize,sizeof(value_type));
      }

    /*
    if (arr!=NULL) delete[] arr;
    arr = new value_type[storageSize];
    */

    for (int k=0;k<(1<<DIMS);k++) integrationPoints[k].clear();
    for (int k=0;k<(1<<DIMS);k++)
      {	
	for (i=0;i<(1<<DIMS);i++)
	  {	
	    long ptr=0;	    
	    for (j=0;j<DIMS;j++)
	      {
		if (i&(1<<j))
		  {
		    if (k&(1<<j)) break;
		    ptr+=stride_Val[j];
		  }
	      }
	    if (j==DIMS) integrationPoints[k].push_back(ptr*n_Fields);
	  }
      }
    
    for (i=0;i<DIMS;i++)
      {
	arr_dims[i]=getValCoord(i).size();
	boxSize[i]=getVertCoord(i).back()-getVertCoord(i).front();	
      }

    /*
      for (i=0;i<(1<<DIMS);i++)
      {	
      integrationPoints[i]=arr;
      for (j=0;j<DIMS;j++)
      {
      if (i&(1<<j))
      integrationPoints[i]+=stride_Val[j];
      }
      }
    */

    //printf("res : %d %d %d",p.Pres[0],p.Ures[0],p.Ures[1]);
    //printf("nval = %ld (%ld %ld)\n",n_Val, n_Cell, n_Vert);

    //if (arr!=NULL) free(arr);
    // arr = (value_type *) malloc((long)sizeof(value_type)*n_Val*n_Fields);
    initialized=true;
  }

public:
  value_type *getDataPtr() {return arr;}
  template <typename T>
  value_type *getDataPtr(T i) {return &arr[i*p.nFields];}
  template <typename T>
  value_type *getDataPtr(T* w) 
  {
    long index=w[0]*stride_Val[0];
    for (int i=1;i<DIMS;i++) 
      index+=stride_Val[i]*w[i];
    return getDataPtr(index);
  }

  long getCellStride(int i) const {return stride_Cell[i];}
  long getVertStride(int i) const {return stride_Vert[i];}
  long getValStride(int i) const {return stride_Val[i];}
  const long *getValStridePtr() const {return stride_Val;}
  const std::vector<long> &getIntegrationPointsPtr(int id=0) const {return integrationPoints[id];}

  virtual long getLowMargin(int i) const {return (i<P_DIMS)?p.lowPmargin[i]:p.lowUmargin[i-P_DIMS];}
  virtual long getHighMargin(int i) const {return (i<P_DIMS)?p.highPmargin[i]:p.highUmargin[i-P_DIMS];}

  long getNVertex() const {return n_Vert;}
  long getNCell() const {return n_Cell;}
  long getNVal() const {return n_Val;}
  long getNFields() const {return n_Fields;}
  long getArrDim(int i) const {return arr_dims[i];}
  double getBoxSize(int i) const {return boxSize[i];}

  const std::vector<double>& getCellCoord(long i) const {return (i<P_DIMS)?P_Cell[i]:U_Cell[i-P_DIMS];}
  const std::vector<double>& getVertCoord(long i) const {return (i<P_DIMS)?P_Vert[i]:U_Vert[i-P_DIMS];}
  const std::vector<double>& getValCoord(long i) const {return (i<P_DIMS)?P_Val[i]:U_Val[i-P_DIMS];}

  const std::vector<double>& getCellCoord_P(long i) const {return P_Cell[i];}
  const std::vector<double>& getVertCoord_P(long i) const {return P_Vert[i];}
  const std::vector<double>& getValCoord_P(long i) const {return P_Val[i];}

  const std::vector<double>& getCellCoord_U(long i) const {return U_Cell[i];}
  const std::vector<double>& getVertCoord_U(long i) const {return U_Vert[i];}
  const std::vector<double>& getValCoord_U(long i) const {return U_Val[i];}

  const std::vector<double>& getCellDelta_P(long i) const {return P_CellDelta[i];}
  const std::vector<double>& getVertDelta_P(long i) const {return P_VertDelta[i];}
  const std::vector<double>& getCellDelta_U(long i) const {return U_CellDelta[i];}
  const std::vector<double>& getVertDelta_U(long i) const {return U_VertDelta[i];}
  const std::vector<double>& getCellDelta(long i) const {return (i<P_DIMS)?P_CellDelta[i]:U_CellDelta[i-P_DIMS];}
  const std::vector<double>& getVertDelta(long i) const {return (i<P_DIMS)?P_VertDelta[i]:U_VertDelta[i-P_DIMS];}

  const std::vector<double>& getCellCoord_J() const {return U_Cell[U_DIMS-J_DIMS];}
  const std::vector<double>& getVertCoord_J() const {return U_Vert[U_DIMS-J_DIMS];}
  const std::vector<double>& getValCoord_J() const {return U_Val[U_DIMS-J_DIMS];}

  const long *getArrDims() const {return arr_dims;} 
  const double *getBoxSize() const {return boxSize;} 

  const paramT &getParams() const {return p;}
  /*
  template <class subParamT>
  void getSubParams(subParamT &sp, int whichP[subParamT::P_DIMS], int whichU[subParamT::U_DIMS]=NULL) const 
  {
    p.getSubGridParams(sp,whichP,whichU);
  }

  template <class subGridT>
  subGridT getSubGrid(int whichP[subGridT::P_DIMS], int whichU[subGridT::U_DIMS]=NULL)
  {
    subGridT result;

    typename subGridT::paramT subP;
    getSubParams(subP,whichP,whichU);
    result.init(subP);

    return result;
  }
  */
  bool isInitialized() const {return initialized;}
  
public:
  
  gridBase():initialized(false),arr(NULL),ownArr(false)
  {
    
  }

  virtual ~gridBase() {
    //if (arr!=NULL) delete[] arr;
    if ((arr!=NULL)&&(ownArr)) free(arr);
  }

private:
  // disable copy-construction
  gridBase( const gridBase& other );
  gridBase& operator=(const gridBase&);

public:

  virtual void init(const paramT &p_, value_type* myArr=NULL)
  {
    p=p_;
    /*
    parent_P_DIMS=parent_P_DIMS_;
    parent_U_DIMS=parent_U_DIMS_;
    if (parent_which_==NULL) 
      for (i=0;i<DIMS;i++) parent_which[i]=-1;
    else
      for (i=0;i<DIMS;i++) parent_which[i]=parent_which_[i];
    */
    init(myArr);
  }

  virtual void write(FILE *f)
  {
    myIO::writeTag(f,getTag());
    p.write(f);   
    fwrite(arr,sizeof(value_type),n_Val*n_Fields,f);
  }
  
  virtual void read(FILE *f, bool swap)
  {
    myIO::checkTag(f,getTag());
    p.read(f,swap);
    init();       
    myIO::fread(arr,sizeof(value_type),n_Val*n_Fields,f,swap);
  }

  iterator begin()
  {
    return iterator(this);
  }

  iterator end()
  {
    return iterator(this,gridNav::undefined(),true);
  }

  iterator margin_begin(dirT dir=gridNav::dir())
  {
    return iterator(this,dir);
  }

  iterator margin_end(dirT dir=gridNav::dir())
  {
    return iterator(this,dir,true);
  }

  long margin_size(dirT dir=gridNav::dir())
  {
    return iterator::boxSize(margin_begin(dir));
  }

  iterator subbox_begin(const int *lp,const int *hp,const int *lu,const int *hu,dirT dir=gridNav::dir())
  {
    return iterator(this,dir,false,lp,hp,lu,hu);
  }

  iterator subbox_end(const int *lp,const int *hp,const int *lu,const int *hu,dirT dir=gridNav::dir())
  {
    return iterator(this,dir,true,lp,hp,lu,hu);
  }

  iterator subbox_begin(const int *l,const int *h,dirT dir=gridNav::dir())
  {
    return iterator(this,dir,false,l,h,(U_DIMS)?(&l[P_DIMS]):l,(U_DIMS)?(&h[P_DIMS]):h);
  }

  iterator subbox_end(const int *l,const int *h,dirT dir=gridNav::dir())
  {
    return iterator(this,dir,true,l,h,&l[P_DIMS],&h[P_DIMS]);
  }

  long subbox_size(const int *lp=NULL,const int *hp=NULL,const int *lu=NULL,const int *hu=NULL,dirT dir=gridNav::dir())
  {
    return iterator::boxSize(subbox_begin(lp,hp,lu,hu,dir));
  }

  iterator innerMargin_begin(const int *dp_, const int *du_,dirT dir=gridNav::dir())
  {
    int lpm[P_DIMS];
    int hpm[P_DIMS];
    int lum[U_DIMS];
    int hum[U_DIMS];
    int i;
  
    for (i=0;i<P_DIMS;i++)
      {
	int dp=dp_[i];
	if (dp==0) dp=P_Val[i].size()-p.highPmargin[i]-p.lowPmargin[i];
	else if (dp<0) dp=P_Val[i].size();

	if (dir&gridNav::dir(i,1)) 
	  {
	    lpm[i]=P_Val[i].size()-dp-p.highPmargin[i];
	    hpm[i]=p.highPmargin[i];
	  }
	else
	  {
	    lpm[i]=p.lowPmargin[i];
	    hpm[i]=P_Val[i].size()-dp-p.lowPmargin[i];
	  }
      }

    for (i=0;i<U_DIMS;i++)
      {
	int du=du_[i];
	if (du<=0) du=U_Val[i].size()-p.highUmargin[i]-p.lowUmargin[i];
	else if (du<0) du=U_Val[i].size();

	if (dir&gridNav::dir(i+P_DIMS,1)) 
	  {
	    lum[i]=U_Val[i].size()-du-p.highUmargin[i];
	    hum[i]=p.highUmargin[i];
	  }
	else
	  {
	    lum[i]=p.lowUmargin[i];
	    hum[i]=U_Val[i].size()-du-p.lowUmargin[i];
	  }
      }
   
    return subbox_begin(lpm,hpm,lum,hum);
  }

  iterator innerMargin_end(const int *dp, const int *du,dirT dir=gridNav::dir())
  {
    return iterator(this,dir,true);
  }

  long innerMargin_size(const int *dp, const int *du,dirT dir=gridNav::dir())
  {
    return iterator::boxSize(innerMargin_begin(dp,du,dir));
  }
  /*
  void pos2Coord(const double p[DIMS], int c[DIMS])
  {
    for (int i=0;i<P_DIMS;i++)
      {
	c[i]=std::distance(P_Val[i].begin(),std::lower_bound(P_Val[i].begin(),P_Val[i].end(),p[i]));	
	//if (p[i]==P_Val[i].back()) c[i]=P_Val[i].size();
      }
    for (int i=0;i<U_DIMS;i++)
      {
    	c[i+P_DIMS]=std::distance(U_Val[i].begin(),std::lower_bound(U_Val[i].begin(),U_Val[i].end(),p[i+P_DIMS]));
	//if (p[i+P_DIMS]==U_Val[i].back()) c[i+P_DIMS]=U_Val[i].size();
      }
  }
  */
  template <class otherItT>
  iterator importIterator(const otherItT &it)
  {
   
    if (p.haveParentGrid)
      return iterator(this,it,p.parent_which);
    else
      return iterator(this,it);
  }

  std::pair<int,int> getParentDims() {return std::make_pair(p.parent_P_DIMS,p.parent_U_DIMS);}

  template <class subGridT>
  bool findGridPos(subGridT &sub, int pos[DIMS], int size[DIMS], bool inner=false)
  {
    typedef typename std::vector<double>::const_iterator const_iteratorT;    
    int which[DIMS];
    //typename subGridT::paramsT &sP=sub.getParams();

    if ((P_DIMS!=subGridT::P_DIMS)||(U_DIMS!=subGridT::U_DIMS))
      {
	if ((p.haveParentGrid)&&(p.parent_P_DIMS==subGridT::P_DIMS)&&(p.parent_U_DIMS==subGridT::U_DIMS))
	  {
	    for (int i=0;i<DIMS;i++) which[i]=p.parent_which[i];
	  }
	else if ((p.parent_P_DIMS<=subGridT::P_DIMS)&&(p.parent_U_DIMS<=subGridT::U_DIMS))
	  {
	    for (int i=0;i<DIMS;i++) 
	      which[i]=(i<P_DIMS)?i:i-P_DIMS+subGridT::P_DIMS;
	  }
	else return false;
      }
    else for (int i=0;i<DIMS;i++) which[i]=i;

    double C[subGridT::DIMS];
    sub.margin_begin().C(C);
      
    for (int ct=0;ct<DIMS;ct++)
      {
	int oct=which[ct];
	const_iteratorT pos_it=hlp::findValue(getValCoord(ct),C[oct]);
	if (pos_it==getValCoord(ct).end()) return false;
	pos[ct]=std::distance(getValCoord(ct).begin(),pos_it);
	size[ct]=sub.getValCoord(oct).size();
	if (inner) 
	  {
	    pos[ct]+=sub.getLowMargin(oct);
	    size[ct]-=(sub.getLowMargin(oct)+sub.getHighMargin(oct));
	  }
	
      }
   
    return true;
  }
  /*
  template <class containerT>
  void interleave(const containerT &other, double fac=1.)
  {
    

    if (other.size()!=p.nFields)
      {
	fprintf(stderr,"ERROR in gridBase::interleave : nFields(%d) and number of conatiners(%ld) do not match.\n",nFields,other.size);
	exit(-1)
      }
    for (int i=0;i<;i++)
    interleave()
  }
  */
  template <class itT>
  void interleave(const itT start,const itT stop, double fac=1.)
  {    
    if (std::distance(start,stop)!=p.nFields)
      {
	fprintf(stderr,"ERROR in gridBase::interleave : nFields(%d) and number of grids(%ld) do not match.\n",p.nFields,std::distance(start,stop));
	exit(-1);
      }

    long i=0;
    value_type *other[p.nFields];

    for (itT it=start;it!=stop;it++)
      other[i++]=&(*it);

    return interleave(other,&other[p.nFields],fac);
  }

  void interleave(myT* start, myT* stop, double fac=1.)
  {
    if (std::distance(start,stop)!=p.nFields)
      {
	fprintf(stderr,"ERROR in gridBase::interleave : nFields(%d) and number of grids(%ld) do not match.\n",p.nFields,std::distance(start,stop));
	exit(-1);
      }
    
    iterator it=begin();
    #pragma omp parallel for
    for (long s=0; s<p.nFields;s++)
      {
	oneField_iterator<iterator> it(begin(),s);
	const oneField_iterator<iterator> it_end(end(),s);
	iterator oit=start[s].importIterator(it);
	
	for (;it!=it_end;it++,oit++)
	  (*it)=fac*(*oit);
      }
  }
  
  void multiply(double factor)
  {
    #pragma omp parallel for 
    for (long i=0;i<n_Val;i++)
      {
	arr[i]*=factor;
      }
  }
  
  void writeToNDfield(const std::string &fname_)
  {
    char fname[255];
    sprintf(fname,"%s.ND",fname_.c_str());
    /*
      if (bT::mpiCom->size()>1)
      sprintf(fname,"%s_%6.6d.ND",fname_.c_str(),bT::mpiCom->rank());
      else
      sprintf(fname,"%s.ND",fname_);
    */
    NDF::NDfield *f=toNDfield();
    NDF::Save_NDfield(f,fname);
    free(f);
  }

  void copyRaw(value_type *d, long N)
  {
    if (p.minStorageSize<N)
      {
	fprintf(stderr,"ERROR in gridBase::copyDataNoCheck: insufficient space to copy, please call set setMinStorage ...\n");
	exit(-1);
      }
    std::copy(d,d+N,arr);
  }

  template <class otherGridT>
  void copyRaw(otherGridT &g)
  {
    copyRaw(g.getDataPtr(),g.getNVal());
  }

  void setMinStorage(long N)
  {
    long after=std::max(N,n_Val*n_Fields);
    long before=std::max(p.minStorageSize,n_Val*n_Fields);
    
    if (before==after) return;
    //printf("Storage before: (%ld,%ld)\n",p.minStorageSize,before);
    if (after<=n_Val*n_Fields)
      p.minStorageSize=0;
    else 
      p.minStorageSize=after;
    //printf("Storage : (%ld,%ld)\n",p.minStorageSize,after);
    arr=(value_type*)realloc(arr,sizeof(value_type)*after);
  }
  
  NDF::NDfield *toNDfield()
  {
    long i;
    const paramT &gp=getParams();
    int useSpecies=(gp.nFields>1)?1:0;
    int dims[P_DIMS+U_DIMS+useSpecies];
    double x0[P_DIMS+U_DIMS+useSpecies];
    double delta[P_DIMS+U_DIMS+useSpecies];
    int ndims;
    long stride[P_DIMS+U_DIMS+1+useSpecies];
    int index;
    value_type *d;
    char comment[80];
    //gridT *grid=bT::gh->getGrid();
    

    stride[0]=1;    
    if (useSpecies)
      {
	dims[0]=gp.nFields;
	x0[0]=0;
	delta[0]=gp.nFields;
	stride[1]=stride[0]*dims[0];
      }

    long dp=useSpecies;
    for (i=0;i<DIMS;i++) {
      dims[i+dp]=getValCoord(i).size();
      x0[i+dp]=getValCoord(i).front();
      delta[i+dp]=getValCoord(i).back()-getValCoord(i).front();
      stride[i+1+dp]=getValStride(i+1);
    }
    
    strcpy(comment,"");
    for (i=0;i<P_DIMS;i++) 
      {
	char tmp[256];
	sprintf(tmp,"%d %d ",gp.Pscale[i],(int)valLocation_P);
	strcat(comment,tmp);
      }
    for (i=0;i<U_DIMS;i++) 
      {
	char tmp[256];
	sprintf(tmp,"%d %d ",gp.Uscale[i],(int)valLocation_U);
	strcat(comment,tmp);
      }
    
   
    d = arr;
    ndims=P_DIMS+U_DIMS+useSpecies;
 
    return NDF::Create_NDfield(dims,ndims,0,ND_DOUBLE,x0,delta,d,comment);
  }

};

#endif
