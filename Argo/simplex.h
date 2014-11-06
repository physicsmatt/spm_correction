#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#define S_MODE_SEQUENTIAL 0x0000
#define S_MODE_OPENCL 0x0001
#define S_MODE_TBB 0x0002


#define S_MODE_TYPE S_MODE_TBB

#if ( S_MODE_TYPE == S_MODE_OPENCL )
	#undef S_MODE_SEQUENTIAL
	#undef S_MODE_TBB

	#include <CL\cl.hpp>
	#include "Kernels.h"
#else
	#if ( S_MODE_TYPE == S_MODE_TBB )
		#include <tbb\task_scheduler_init.h>
		#include <tbb\parallel_for.h>
		#include <tbb\critical_section.h>
		#undef S_MODE_SEQUENTIAL
		#undef S_MODE_OPENCL


	#elif ( S_MODE_TYPE == S_MODE_SEQUENTIAL )
		#undef S_MODE_OPENCL
		#undef S_MODE_TBB
	#endif
#endif

#include "FImage.h"


#define S_WRITE_INTERP F_BSPLINE
#define S_INTERP_CONST_A -0.5
#if ( S_WRITE_INTERP == F_BSPLINE )	// Cubic B-spline
	#define S_INTERP_CONST_B 1
	#define S_INTERP_CONST_C 0
#elif ( S_WRITE_INTERP == F_CATMULL_BSPLINE )	// Catmull-Rom Spline
	#define S_INTERP_CONST_B 0
	#define S_INTERP_CONST_C 0.5
#else // Mitchell-Netravali Cubic Filter ( Not sure if defines can be fractions. Test this later. )
	#define	S_INTERP_CONST_B 1.0 / 3.0
	#define	S_INTERP_CONST_C 1.0 / 3.0
#endif


int simplex ( FImage *baseImage, FImage *sliverImage, double x[], double z[], bool _fastZ, double precision[], double alpha, double beta, double gamma,
		int maxi );


#endif _SIMPLEX_H_
