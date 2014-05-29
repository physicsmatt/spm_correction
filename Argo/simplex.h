#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include "FImage.h"
int simplex ( FImage *baseImage, FImage *sliverImage, double x[], double z[], int n, double precision[], double alpha, double beta, double gamma,
		double errhalt, int maxi );
///new things to add for regions
//,int nr, int enreg, int ROH, int RH, double InitC3

#endif _SIMPLEX_H_
