#ifndef __OPENMP_INTERFACE_HXX__
#define __OPENMP_INTERFACE_HXX__

#ifdef USE_OPENMP

#include "openmp.h"

#else

int omp_get_thread_num() {return 0;}
int omp_get_num_threads() {return 1;}
void omp_set_num_threads(int n) {return;}

#endif //USE_OPENMP

#endif
