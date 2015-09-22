#ifndef __DIMS_TRAITS__HXX__
#define __DIMS_TRAITS__HXX__

template <int PD, int VD>
struct dimTraits {
  static const int P_DIMS=PD;
  static const int V_DIMS=VD;
};

#endif
