#ifndef _GRIDNAV__HXX_
#define _GRIDNAV__HXX_

#include <vector>

struct gridNav {
  typedef long dirT;
  static const dirT UNDEFINED=(1l <<63);
  
  static dirT dir(const char *dim,const char *dir)
  {
    int i=0;
    dirT result=0;

    do {
      int d=2*(dim[i]-'0');
      if (dir[i]=='+') result|= 1<<(d);
      else if (dir[i]=='-') result|= 1<<(d+1);
      else break;
      i++;
    } while(true);

    return result;
  }

  static dirT undefined()
  {
    return UNDEFINED;
  }

  static dirT dir()
  {
    return 0;
  }
  
  static dirT dir(const int dim,const int dir)
  {
    if (!dir) return 0;
    return (dir>0)?(1<<(dim*2)):(1<<(dim*2+1));
  }
  
  static void setDim(dirT &which, const int dim, const int dir)
  {
    which &= (~((1<<(2*dim))|(1<<(2*dim+1))));
    if (!dir) return;
    if (dir>0) which |=  (1<<(dim*2));
    else which |=  (1<<(dim*2+1));     
  }

  static dirT dir(const std::vector<int> &dim,const std::vector<int> &dir)
  {  
    dirT tmp=0;
    long i;
    for (i=0;i<dir.size();i++)
      {
	if (!dir[i]) continue;
	if (dir[i]>0) tmp |=  (1<<(dim[i]*2));
	else tmp |=  (1<<(dim[i]*2+1));   
      }
   
    return tmp;
  }

  static dirT reverse(dirT dir)
  {
    if (dir==UNDEFINED) return dir;

    dirT result=0;
    dirT test=3;
    int i;
    long bef=(long)dir;
    const long m[4]={0,2,1,0};
    for (i=0;i<sizeof(dirT)*8;i+=2)
      result|=(m[(dir&(test<<i))>>i]<<i);
    //printf("%ld -> %ld\n",bef,result);
    return result;
  }

};

#endif
