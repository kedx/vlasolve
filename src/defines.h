#ifndef __MYTYPES_H__
#define __MYTYPES_H__

#define PI (3.141592653589793116F)
#define RAD2DEG ((180.0d)/(3.141592653589793116F))
#define DEG2RAD ((3.141592653589793116F)/(180.0F))

#define PISQ 9.869604401089357992F
#define TWOPI 6.283185307179586232F
#define PI_INV 0.31830988618379069122F
#define TWOPI_INV 0.15915494309189534561F

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
