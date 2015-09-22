#ifndef __MYTYPES_H__
#define __MYTYPES_H__

#define PI ((double)3.141592653589793116)
#define RAD2DEG (((double)180.)/((double)3.141592653589793116))
#define DEG2RAD (((double)3.141592653589793116)/((double)180.))

#define TWOPI 6.283185307179586232
#define PI_INV 0.31830988618379069122
#define TWOPI_INV 0.15915494309189534561

#ifdef USELONGINT
typedef long int INT;
typedef unsigned long int UINT;
#else
typedef int INT;
typedef unsigned int UINT;
#endif

#ifdef USESIMPLEPRECISION
typedef float FLOAT;
#else
typedef double FLOAT;
#endif

typedef unsigned int uint;
typedef unsigned long ulong;

#endif
