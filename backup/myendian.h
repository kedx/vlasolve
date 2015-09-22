#ifndef __ENDIAN_SWAP
#define __ENDIAN_SWAP

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif


int fread_sw(void *,size_t,size_t,FILE *,int);

int swapI(int);
float swapF(float);
double swapD(double);
void Dswap2B(void*);
void Dswap4B(void*);
void Dswap8B(void*);
void Dswap2BArr(void*,size_t);
void Dswap4BArr(void*,size_t);
void Dswap8BArr(void*,size_t);

#ifdef __cplusplus
}
#endif


#endif
