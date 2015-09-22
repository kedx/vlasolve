#ifndef __VAL_LOCATION_TRAITS_HXX__
#define __VAL_LOCATION_TRAITS_HXX__

#include "valLocationType.hxx"

template <valLocationType P, valLocationType U, valLocationType J=valLocationVal::UNDEFINED>
struct valLocation_traits
{
  static const valLocationType VP = P;
  static const valLocationType VU = U;
  static const valLocationType VJ = J;
};


#endif
