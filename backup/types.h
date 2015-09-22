#ifndef __MYTYPES_H__
#define __MYTYPES_H__

#define PI ((double)3.141592653589793116)
#define RAD2DEG (((double)180.)/((double)3.141592653589793116))
#define DEG2RAD (((double)3.141592653589793116)/((double)180.))

#define TWOPI 6.283185307179586232
#define PI_INV 0.31830988618379069122
#define TWOPI_INV 0.15915494309189534561

#ifdef SIMPLEPRECISION
#define FLOAT float
#define FFTW_COMPLEX fftwf_complex
#define FFTW_PLAN fftwf_plan
#define FFTWNAME(x) fftwf_##x
#else
#define FLOAT double
#define FFTW_COMPLEX fftw_complex
#define FFTW_PLAN fftw_plan
#define FFTWNAME(x) fftw_##x
#endif

#ifdef USELONGINT
#define INT long int
#else
#define INT int
#endif

typedef unsigned int uint;
typedef unsigned long ulong;

#endif
