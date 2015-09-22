#ifndef __GRID_BASE_HXX__
#define __GRID_BASE_HXX__

#include <map>
#include <vector>
#include <stdio.h>

#include "scale.hxx"
#include "dimTraits.hxx"
#include "myIO.hxx"

template <
  class dimensionTraits,
  typename dataType,
  valLocationType valLocationTemplate
  >
class gridBase {
public:

  typedef gridBase<dimensionTraits,dataType,valLocationTemplate> myT;

  static std::string getTag() {return "GRID_BASE v0.1";}

  typedef dataType dataT;
  typedef scale<dataType> scaleT;
  typedef dimensionTraits dimTr;

  typedef typename scaleT::scaleTypeT scaleTypeT;
  typedef typename scaleT::scaleTypeV scaleTypeV;
  typedef typename scaleT::valLocationT valLocationT;
  typedef typename scaleT::valLocationV valLocationV;

  static const int P_DIMS= dimensionTraits::P_DIMS;
  static const int U_DIMS= dimensionTraits::U_DIMS;
  static const int J_DIMS= U_DIMS-P_DIMS;
  static const int DIMS = U_DIMS+P_DIMS;
  static const valLocationT valLocation = valLocationTemplate;

  struct paramT {
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
    }
  };

  static const std::string parserCategory() {return "grid";}

protected:
  long n_Vert;
  long n_Cell;
  long n_Val;
  
  std::vector<double> P_Vert[P_DIMS];
  std::vector<double> U_Vert[U_DIMS];
  std::vector<double> P_Cell[P_DIMS];
  std::vector<double> U_Cell[U_DIMS];

  std::vector<double> P_Val[P_DIMS];
  std::vector<double> U_Val[U_DIMS];

  long stride_Vert[P_DIMS+U_DIMS+1];
  long stride_Cell[P_DIMS+U_DIMS+1];
  long stride_Val[P_DIMS+U_DIMS+1];

  paramT p;
  dataT *arr;
  
private:
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
  void init()
  {    
    int i,j;
    
    n_Vert=1; 
    n_Cell=1; 
    n_Val=1;
    stride_Vert[0]=1;
    stride_Cell[0]=1;
    stride_Val[0]=1;

    for (i=0;i<P_DIMS;i++) {
      
      P_Vert[i]=scaleT::genScale(p.Pmin[i],p.Pmax[i],p.Pres[i],p.Pscale[i],valLocationV::VERTEX);
      P_Cell[i]=scaleT::genScale(p.Pmin[i],p.Pmax[i],p.Pres[i],p.Pscale[i],valLocationV::CELL);
      P_Val[i]=scaleT::genScale(p.Pmin[i],p.Pmax[i],p.Pres[i],p.Pscale[i],valLocation);
      
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

      if (i==(U_DIMS-J_DIMS)) 
	U_Val[i]=scaleT::genScale(p.Umin[i],p.Umax[i],p.Ures[i],p.Uscale[i],valLocationV::CELL);
      else 
	U_Val[i]=scaleT::genScale(p.Umin[i],p.Umax[i],p.Ures[i],p.Uscale[i],valLocation);

      n_Vert*=U_Vert[i].size();
      n_Cell*=U_Cell[i].size();
      n_Val*=U_Val[i].size();

      stride_Cell[i+P_DIMS+1]=stride_Cell[i+P_DIMS]*U_Cell[i].size();
      stride_Vert[i+P_DIMS+1]=stride_Vert[i+P_DIMS]*U_Vert[i].size();
      stride_Val[i+P_DIMS+1]=stride_Val[i+P_DIMS]*U_Val[i].size();
    }
    //printf("res : %d %d %d",p.Pres[0],p.Ures[0],p.Ures[1]);
    //printf("nval = %ld (%ld %ld)\n",n_Val, n_Cell, n_Vert);
    
    arr = new dataT[n_Val];   
  }

public:
  dataT *getDataPtr() {return arr;}

  long getCellStride(int i) const {return stride_Cell[i];}
  long getVertStride(int i) const {return stride_Vert[i];}
  long getValStride(int i) const {return stride_Val[i];}

  virtual long getLowMargin(int i) const {return (i<P_DIMS)?p.lowPmargin[i]:p.lowUmargin[i-P_DIMS];}
  virtual long getHighMargin(int i) const {return (i<P_DIMS)?p.highPmargin[i]:p.highUmargin[i-P_DIMS];}

  long getNVertex() const {return n_Vert;}
  long getNCell() const {return n_Cell;}
  long getNVal() const {return n_Val;}

  const std::vector<double>& getCellCoord(long i) const {return (i<P_DIMS)?P_Cell[i]:U_Cell[i-P_DIMS];}
  const std::vector<double>& getVertCoord(long i) const {return (i<P_DIMS)?P_Vert[i]:U_Vert[i-P_DIMS];}
  const std::vector<double>& getValCoord(long i) const {return (i<P_DIMS)?P_Val[i]:U_Val[i-P_DIMS];}

  const std::vector<double>& getCellCoord_P(long i) const {return P_Cell[i];}
  const std::vector<double>& getVertCoord_P(long i) const {return P_Vert[i];}
  const std::vector<double>& getValCoord_P(long i) const {return P_Val[i];}

  const std::vector<double>& getCellCoord_U(long i) const {return U_Cell[i];}
  const std::vector<double>& getVertCoord_U(long i) const {return U_Vert[i];}
  const std::vector<double>& getValCoord_U(long i) const {return U_Val[i];}

  const std::vector<double>& getCellCoord_J() const {return U_Cell[U_DIMS-J_DIMS];}
  const std::vector<double>& getVertCoord_J() const {return U_Vert[U_DIMS-J_DIMS];}
  const std::vector<double>& getValCoord_J() const {return U_Val[U_DIMS-J_DIMS];}

  const paramT &getParams() const {return p;}
  
public:
  
  gridBase()
  {
    arr=NULL;
  }

  virtual ~gridBase() {
    delete arr;
  }

  virtual void init(const paramT &p_)
  {
    p=p_;
    init();
  }

  virtual void write(FILE *f)
  {
    myIO::writeTag(f,getTag());
    p.write(f);   
    fwrite(arr,sizeof(dataT),n_Val,f);
    
  }
  
  virtual void read(FILE *f, bool swap)
  {
    myIO::checkTag(f,getTag());
    p.read(f,swap);
    init();       
    myIO::fread(arr,sizeof(dataT),n_Val,f,swap);
    
  }
  
};

#endif
