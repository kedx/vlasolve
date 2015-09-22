#ifndef __SNAPSHOT_HXX__
#define __SNAPSHOT_HXX__

#include <iostream>
#include <queue>
#include <functional>
#include <algorithm>
#include <set>
#include <map>

template <class gridHandlerType>
class snapshot {
public:
  typedef gridHandlerType gridHandlerT;
  
  typedef typename gridHandlerT::dataT dataT;
  typedef typename gridHandlerT::gridT gridT;
  
  snapshot(const gridHandlerT *gh_):
  gh(gh_)
  {
    
  }
  
  ~snapshot()
  {
    
  }

  
private:
  const gridHandlerT *gh;

};

#endif
