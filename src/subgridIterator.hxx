#ifndef __SUBGRID_ITERATOR_HXX_
#define __SUBGRID_ITERATOR_HXX_

#include <iterator>
#include "gridNav.hxx"
#include "fieldIterator.hxx"
#include "helpers.hxx"

//template <class gridType> class subgridIterator;

template <class gridType>
class subgridIterator
{
public:
  template <class gridType2> friend class subgridIterator;  

  typedef std::forward_iterator_tag iterator_category;
  typedef subgridIterator<gridType> self_type;

  typedef gridType container_type;
  typedef typename container_type::value_type value_type;
  typedef value_type* pointer;
  typedef value_type& reference;
  typedef const value_type* const_pointer;
  typedef const value_type& const_reference;
  typedef long  difference_type;

  typedef typename container_type::valLocationT valLocationT;
  typedef typename container_type::valLocationV valLocationV;

  static const int P_DIMS = container_type::P_DIMS;
  static const int U_DIMS = container_type::U_DIMS;
  static const int DIMS = container_type::DIMS;

  //static const valLocationT valLocation = container_type::valLocation;

  typedef typename gridNav::dirT dirT;
  
  static const long N_INT_P = (1<<DIMS);
  static const double N_INT_P_INV = 1.0F/(1<<DIMS);

protected:
  long delta[DIMS-1];
  long nfields;
  long i; 
  value_type *data;
  container_type *g;
  int w[DIMS];
  int min[DIMS];
  int max[DIMS];

public:
    
  subgridIterator():
    g(NULL),
    i(-1),data(NULL)
  {
    
  }

  subgridIterator(container_type *g_, dirT region=gridNav::undefined(), bool end=false,const int *lowPmargin=NULL, const int *highPmargin=NULL,const int *lowUmargin=NULL, const int *highUmargin=NULL):
    g(g_),
    data(g_->getDataPtr())
  {
    long ct;
    const int *lpm=(lowPmargin==NULL)?g->getParams().lowPmargin:lowPmargin;
    const int *hpm=(highPmargin==NULL)?g->getParams().highPmargin:highPmargin;
    const int *lum=(lowUmargin==NULL)?g->getParams().lowUmargin:lowUmargin;
    const int *hum=(highUmargin==NULL)?g->getParams().highUmargin:highUmargin;
    nfields = g->getNFields();
    if (end) {i=-1;data=NULL;}
    else 
      {
	i=0;
	for (ct=0;ct<DIMS;ct++)
	  {
	    if (region==gridNav::undefined())
	      {
		min[ct]=0;
		max[ct]=g->getArrDim(ct);//(ct<P_DIMS)?g->getValCoord_P(ct).size():g->getValCoord_U(ct-P_DIMS).size();
	      }
	    else if (region&gridNav::dir(ct,-1))
	      {
		min[ct]=0;
		max[ct]=(ct<P_DIMS)?lpm[ct]:lum[ct-P_DIMS];
	      }
	    else if (region&gridNav::dir(ct,+1))
	      {
		min[ct]=max[ct]=g->getArrDim(ct);//(ct<P_DIMS)?g->getValCoord_P(ct).size():g->getValCoord_U(ct-P_DIMS).size();
		min[ct]-=(ct<P_DIMS)?hpm[ct]:hum[ct-P_DIMS];
	      }
	    else
	      {
		min[ct]=(ct<P_DIMS)?lpm[ct]:lum[ct-P_DIMS];
		max[ct]=g->getArrDim(ct);//(ct<P_DIMS)?g->getValCoord_P(ct).size():g->getValCoord_U(ct-P_DIMS).size();
		max[ct]-=(ct<P_DIMS)?hpm[ct]:hum[ct-P_DIMS];
	      }

	    if (max[ct]<=min[ct])
	      {
		i=-1;
		break;
	      }

	    w[ct] = min[ct];
	    i+=min[ct]*g->getValStride(ct);
	  }
	for (ct=0;ct<DIMS-1;ct++) delta[ct]=(g->getValStride(ct+1)-(max[ct]-min[ct])*g->getValStride(ct))/nfields;
	//printf("margin: (%d %d)(%d %d)(%d %d) min/max: (%d %d)(%d %d)(%d %d)\n",lpm[0],hpm[0],lum[0],hum[0],lum[1],hum[1],min[0],max[0],min[1],max[1],min[2],max[2]);
	data+=i;
	i/=nfields;
      }
  }

  
  template <class otherItT>
  subgridIterator(container_type *g_, const otherItT &it,const int *which=NULL):
    g(g_),
    data(g_->getDataPtr())
  {
    typedef typename std::vector<double>::const_iterator const_iteratorT;
    double C[otherItT::DIMS];
    
    it.C(C);
    nfields = g->getNFields();

    if (it.data==NULL) {i=-1;data=NULL;}
    else
      {	
	i=0;
	for (int ct=0;ct<DIMS;ct++)
	  {
	    int oct=(which!=NULL)?(which[ct]):((ct>P_DIMS)?ct+(otherItT::P_DIMS-P_DIMS):ct);
	    const_iteratorT pos_it=hlp::findValue(g->getValCoord(ct),C[oct]);
	    
	    if (pos_it==g->getValCoord(ct).end())
	      {
		    
		fprintf(stderr,"ERROR in subgridIterator constructor : trying to build an iterator from incompatible grid type.\n");
		fprintf(stderr," this->coord[%d]=%g != other->Coord[%d]=%g.(delta=%g)\n",ct,*pos_it,oct,C[oct],C[oct]-(*pos_it));
		exit(-1);
	      }
	    /*
	    const_iteratorT pos_it=std::lower_bound(g->getValCoord(ct).begin(),g->getValCoord(ct).end(),C[oct]);
	    if ((*pos_it)!=C[oct])
	      {
		if (fabs(((*pos_it)-C[oct])/((*pos_it)+C[oct]))>1.E-10)
		  {
		    pos_it--;
		    if ((*pos_it)!=C[oct])
		      {
			if (fabs(((*pos_it)-C[oct])/((*pos_it)+C[oct]))>1.E-10)
			  {
		    
			    fprintf(stderr,"ERROR in subgridIterator constructor : trying to build an iterator from incompatible grid type.\n");
			    fprintf(stderr," this->coord[%d]=%g != other->Coord[%d]=%g.(delta=%g)\n",ct,*pos_it,oct,C[oct],C[oct]-(*pos_it));
			    exit(-1);
			  }
		      }
		  }
	      }
	    */
	    w[ct]=min[ct]=std::distance(g->getValCoord(ct).begin(),pos_it);	    
	    max[ct]=min[ct]+it.max[oct]-it.min[oct];
	    i+=min[ct]*g->getValStride(ct);
	  }
	for (int ct=0;ct<DIMS-1;ct++) delta[ct]=(g->getValStride(ct+1)-(max[ct]-min[ct])*g->getValStride(ct))/nfields;
	data+=i;
	i/=nfields;
      }
    //print();
    //it.print();
  }

  virtual ~subgridIterator(){}

  void print() const
  {
    int ii;

    char str[2000];
    sprintf(str,"delta = [");
    for (ii=0;ii<DIMS-1;ii++) sprintf(str,"%s%ld ",str,delta[ii]);
    sprintf(str,"%s];nfields=%ld, i=%ld, data=%ld(+%ld); ",str,nfields,i,(long)data,(long)(data-g->getDataPtr()));
    sprintf(str,"%sw = [",str);
    for (ii=0;ii<DIMS;ii++) sprintf(str,"%s%d ",str,w[ii]);sprintf(str,"%s]; ",str);
    sprintf(str,"%smin = [",str);
    for (ii=0;ii<DIMS;ii++) sprintf(str,"%s%d ",str,min[ii]);sprintf(str,"%s]; ",str);
    sprintf(str,"%smax = [",str);
    for (ii=0;ii<DIMS;ii++) sprintf(str,"%s%d ",str,max[ii]);sprintf(str,"%s];",str);
    printf("%s\n",str);
  }
    
  self_type &operator++()
  {     
    if (data==NULL) return *this;
    w[0]++;data+=nfields;i++;

    if (w[0]>=max[0])
      {
	if (DIMS==1) {i=-1;data=NULL;}
	else
	  {
	    w[0]=min[0];
	    w[1]++;
	    i+=delta[0];data+=nfields*delta[0];
	    if (w[1]>=max[1])
	      {
		w[1]=min[1];
		int ct=2;
		while (ct<DIMS)
		  {
		    w[ct]++;i+=delta[ct-1];data+=nfields*delta[ct-1];
		    if (w[ct]>=max[ct])
		      {
			w[ct]=min[ct];
			ct++;
		      }
		    else break;
		  };
		if (ct==DIMS) {i=-1;data=NULL;}
	      }
	  }
      }
     
    return *this;
  }

  self_type operator++(int)
  {
    self_type it(*this);
    ++(*this);
    return it;
  }

  reference operator*() const
  {    
    return *data;
    //return (g->getDataPtr())[i];
  }

  long get_i() const {return i;}
  const int *get_w() const {return w;}

  bool operator==(const self_type& r) const
  {return (data==r.data);}

  bool operator!=(const self_type& r) const
  {return (data!=r.data);}

  bool operator<(const self_type& r) const
  {return (data<r.data);}

  bool operator>(const self_type& r) const
  {return (data>r.data);}	  

  bool operator<=(const self_type& r) const
  {return (data<=r.data);}	  	   

  bool operator>=(const self_type& r) const
  {return (data>=r.data);}

  static long boxSize(const self_type& it)
  {
    long result=1;
    int s;
    if (it.data==NULL) return -1;
    for (s=0;s<DIMS;s++) result *=it.max[s]-it.min[s];
    return result;
  }

  static long dim(const self_type& it, int d)
  {
    return (d<DIMS)?(it.max[d]-it.min[d]):0;
  }

  
  void cellC(double *res) const 
  {
    for (int id=0;id<P_DIMS;id++) res[id]=g->getCellCoord_P(id)[w[id]];
    for (int id=P_DIMS;id<DIMS;id++) res[id]=g->getCellCoord_U(id-P_DIMS)[w[id]];
  }

  void cellCP(double *res) const 
  {
    for (int id=0;id<P_DIMS;id++) res[id]=g->getCellCoord_P(id)[w[id]];
  }

  void cellCU(double *res) const 
  {
    for (int id=0;id<U_DIMS;id++) res[id]=g->getCellCoord_U(id)[w[id+P_DIMS]];
  }

  bool cellD(double *res) const
  {
    bool nonNull=true;
    for (int id=0;id<P_DIMS;id++) 
      if ((res[id]=g->getCellDelta_P(id)[w[id]])==0) nonNull=false;
    for (int id=P_DIMS;id<DIMS;id++)
      if ((res[id]=g->getCellDelta_U(id-P_DIMS)[w[id]])==0) nonNull=false; 

    return nonNull;
  }

  double cellVol() const
  {
    double result=1.;

    for (int id=0;id<P_DIMS;id++) 
      result*=g->getCellDelta_P(id)[w[id]];
    for (int id=P_DIMS;id<DIMS;id++)
      result*=g->getCellDelta_U(id-P_DIMS)[w[id]];

    return result;
  }

  /*
  bool get_DC(double *res) const 
  {
    bool nonNull=true;
    for (int id=0;id<P_DIMS;id++) 
      {
	const int ww=w[id];	
	const std::vector<double>&p=g->getValCoord_P(id);
	if (ww>=p.size()-1)
	   {
	    res[id]=0;
	    nonNull=false;
	  }
	else res[id]=(p[ww+1]-p[ww]);
	//res[id]=(ww>=p.size()-1)?(nonNull=false):(p[ww+1]-p[ww]);
      }
    for (int id=P_DIMS;id<DIMS;id++)
      {
	const int ww=w[id];	
	const std::vector<double>&p=g->getValCoord_U(id-P_DIMS);
	if (ww>=p.size()-1)
	   {
	    res[id]=0;
	    nonNull=false;
	  }
	else res[id]=(p[ww+1]-p[ww]);
	//res[id]=(ww>=p.size()-1)?(nonNull=false):(p[ww+1]-p[ww]);
      }
    return nonNull;
  }
  */

  bool cellDP(double *res) const 
  {
    bool nonNull=true;
    for (int id=0;id<P_DIMS;id++) 
      if ((res[id]=g->getCellDelta_P(id)[w[id]])==0) nonNull=false;
    return nonNull;
  }

  /*
  bool get_DCP(double *res) const 
  {
    bool nonNull=true;
    for (int id=0;id<P_DIMS;id++) 
      {
	const int ww=w[id];	
	const std::vector<double>&p=g->getValCoord_P(id);
	if (ww>=p.size()-1)
	   {
	    res[id]=0;
	    nonNull=false;
	  }
	else res[id]=(p[ww+1]-p[ww]);
      }
    
    return nonNull;
  }
  */

  bool cellDU(double *res) const 
  {
    bool nonNull=true;
    for (int id=0;id<U_DIMS;id++) 
      if ((res[id]=g->getCellDelta_U(id)[w[id+P_DIMS]])==0) nonNull=false;
    return nonNull;
  }

  /*
  bool get_DCU(double *res) const 
  {
    bool nonNull=true;   
    for (int id=0;id<U_DIMS;id++)
      {
	const int ww=w[id+P_DIMS];	
	const std::vector<double>&p=g->getValCoord_U(id);
	if (ww>=p.size()-1)
	  {
	    res[id]=0;
	    nonNull=false;
	  }
	else res[id]=(p[ww+1]-p[ww]);
      }
    return nonNull;
  }
  */

  double cellAvg(int s=0) const
  {    
    const std::vector<long> &itPt=g->getIntegrationPointsPtr();
    //long delta=(long)(data-itPt[0])+s;
    double res=0;

    for (int i=0;i<N_INT_P;i+=4)
      res+=data[itPt[i]]+data[itPt[i+1]]+data[itPt[i+2]]+data[itPt[i+3]];
    return res*N_INT_P_INV;
  }

  void C(double *res) const 
  {
    for (int id=0;id<P_DIMS;id++) res[id]=g->getValCoord_P(id)[w[id]];
    for (int id=P_DIMS;id<DIMS;id++) res[id]=g->getValCoord_U(id-P_DIMS)[w[id]];
  }

  void CP(double *res) const 
  {
    for (int id=0;id<P_DIMS;id++) res[id]=g->getValCoord_P(id)[w[id]];
  }

  void CU(double *res) const 
  {
    for (int id=0;id<U_DIMS;id++) res[id]=g->getValCoord_U(id)[w[id+P_DIMS]];
  }

  
  double vol(int s=0) const
  {
    double result=1;
    for (int id=0;id<P_DIMS;id++) 
      result*=g->getVertDelta_P(id)[w[id]];
    for (int id=P_DIMS;id<DIMS;id++)
      result*=g->getVertDelta_U(id-P_DIMS)[w[id]];
    return result;    
  }

};

template <>
template <class gridType>
class field_iterator< subgridIterator<gridType> >
  :public subgridIterator<gridType>
{
private:
  typedef subgridIterator<gridType> bT;
  int fct;
public:
  typedef field_iterator<bT> self_type;

  field_iterator(const bT &it):bT(it),fct(0)
  {
    
  }

  self_type &operator++()
  {     
    if (bT::data==NULL) return *this;
    ++bT::data;
    if ((++fct)==bT::nfields) {fct=0;bT::w[0]++;bT::i++;}
    else return *this;

    if (bT::w[0]>=bT::max[0])
      {
	if (bT::DIMS==1) {bT::i=-1;bT::data=NULL;}
	else
	  {
	    bT::w[0]=bT::min[0];
	    bT::w[1]++;
	    bT::i+=bT::delta[0];bT::data+=bT::nfields*bT::delta[0];
	    if (bT::w[1]>=bT::max[1])
	      {
		bT::w[1]=bT::min[1];
		int ct=2;
		while (ct<bT::DIMS)
		  {
		    bT::w[ct]++;bT::i+=bT::delta[ct-1];bT::data+=bT::nfields*bT::delta[ct-1];
		    if (bT::w[ct]>=bT::max[ct])
		      {
			bT::w[ct]=bT::min[ct];
			ct++;
		      }
		    else break;
		  };
		if (ct==bT::DIMS) {bT::i=-1;bT::data=NULL;}
	      }
	  }
      }
     
    return *this;
  }

  self_type operator++(int)
  {
    self_type it(*this);
    ++(*this);
    return it;
  }
  
};

template <>
template <class gridType>
class oneField_iterator< subgridIterator<gridType> > 
  :public subgridIterator<gridType>
{
private:
  typedef subgridIterator<gridType> bT;
 
public:
  typedef oneField_iterator<bT> self_type;

  oneField_iterator(const bT &it, int id):bT(it)
  {
    if (bT::data!=NULL) bT::data+=id;
  }
}; 

/*

template <class T>
class field_iterator : public T
{
private:
  typedef T bT;
  int fct;
public:
  typedef field_iterator<T> self_type;

  field_iterator(T &it):T(it),fct(0)
  {
    
  }

  self_type operator++()
  {     
    if (bT::i<0) return *this;
    ++bT::i;
    if ((++fct)==bT::nfields) {fct=0;bT::w[0]++;}
    else return *this;

    if (bT::w[0]>=bT::max[0])
      {
	bT::w[0]=bT::min[0];
	bT::w[1]++;
	bT::i+=bT::delta[0];
	if (bT::w[1]>=bT::max[1])
	  {
	    bT::w[1]=bT::min[1];
	    int ct=2;
	    while (ct<bT::DIMS)
	      {
		bT::w[ct]++;bT::i+=bT::delta[ct-1];
		if (bT::w[ct]>=bT::max[ct])
		  {
		    bT::w[ct]=bT::min[ct];
		    ct++;
		  }
		else break;
	      };
	    if (ct==bT::DIMS) bT::i=-1;
	  }
      }
     
    return *this;
  }

  self_type operator++(int)
  {
    self_type it(*this);
    ++(*this);
    return it;
  }
  
};

template <class T>
class oneField_iterator : public T
{
private:
  typedef T bT;
 
public:
  typedef oneField_iterator<T> self_type;

  oneField_iterator(T &it, int id):T(it)
  {
    bT::i+=id;
  }
}; 
*/

#endif
