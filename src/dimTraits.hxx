#ifndef __DIMS_TRAITS__HXX__
#define __DIMS_TRAITS__HXX__

template <int PD, int UD>
struct dimTraits {
  static const int P_DIMS=PD;
  static const int U_DIMS=UD;
};

#endif
