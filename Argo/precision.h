#ifndef _PRECISION_H_
#define _PRECISION_H_
#ifndef USE_DOUBLE
#define USE_DOUBLE 1
#endif

#if USE_DOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double scalar_t;
#else
typedef float scalar_t;
#endif

#endif _PRECISION_H_