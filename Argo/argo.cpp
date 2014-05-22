/*********************************************************************
 *Image Master procedural image registration and stabilization routine*
 *Coded by: Brian Salmons & Dr. Matthew Trawick                       *
 **********************************************************************
 *Dr. Matthew Trawick                                                 *
 *University of Richmond                                              *
 **********************************************************************
 */

/** Version 1.1
 *
 *Known Issues:
 *-Sliverdrift is SIGNIFICANTLY slower.  Need to try to speed this up.
 *-Images still have to be flipped.  It seems... random.  I need to find the cause of this.
 *
 *-Significant speed improvements for non-sliverdrift version.
 *-Fixed annoying division bug resulting from changes to the way areas are calculated
 *-Program now outputs the version built on (up to the coder to change it)
 *-Initial simplex values are now no longer output... they aren't needed.
 *-Significant changes to the way that sliver drift computations occur.
 *-Printouts from simplex supressed.
 */

//#define WINVER 0x0501 
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <math.h>
#include <time.h>
#include "Eigen\Eigen"
#include <stdio.h>
#include <valarray>
#include <algorithm>
#include "image_basic.h"
#include "simplex.h"
#include <tbb\parallel_for.h>
#include <tbb\critical_section.h>

using namespace tbb;
//#define SANITYCHECKS
#define VERSION "1.1"
#define CODEFOLD

bool LOG_OUTPUT = false;

struct times	{
	double combos_time = 0, diffs_time = 0,
			grid_time = 0, simplex_time = 0,
			image_read_time = 0, image_write_time = 0;
	double total_time = 0;
};

/**
 * Struct for holding the 6 Beta values used to evaluate the C coefficients.
 */
struct beta_values {
	double B0, B1, B2, B3, B4, B5;
};

/**
 * Struct for holding the 3 gamma values representing the 3 points on our parabola used in the creation of diffs.
 */
struct gamma_values {
	double G0, G1, G2;
};

/**
 * Parameter Combo Struct.
 */
struct param_combo {
	double P1, P2, P3;
};

struct input_params	{

	// File names for sliver and base images.
	std::string rundata1, rundata2;

	// Size of a block in pixels.
	int blocksize = 0;

	// Maximum possible constant shift in pixels. Defaults to 5.
	int MaxA0 = 5, MaxB0 = 5;

	/* Maximum linear shift, as fraction of image size
	 * (if a1_Max = 0.02, and image height is 500 pixels,
	 * this term causes max shift at top of image of 10 pixels)
	 */
	double a1_Max = 0.1, b1_Max = 0.1;

	/* Fractional shift, relative to a1_max.
	 * If a1_max implies max shift of 10 pixels at top, a2_Multiplier of 0.5 causes deviation
	 * from linearity of 5 pixels at midpoint.
	 */
	double a2_Multiplier = ( 1.0 / 9.0 ), b2_Multiplier = ( 1.0 / 9.0 );

	// Fractional shift for cubic term a3, also relative to shift caused by a1_max.
	double a3_Multiplier = .1, b3_Multiplier = .1;  //(0.0/9.0);

	/**
	 * For initial grid search, image is rescaled (downsampled) by this factor.
	 * This parameter should be set to the size, in pixels, of the smallest real feature on the image.
	 */
	int precision = 1;

	/**
	 * These parameters control the simplex routine.
	 * These parameters were found, using trial and error on not a whole lot of different test cases.
	 * To avoid being trapped in some rut, it seemse to be VERY important to have the reflection
	 * parameter slightly less than 1.0, and the growth parameter less than 2.0.
	 */
	double growthParam = 1.5, contractParam = 0.5, reflectParam = 0.9, haltParam = 1e-11;
	int maxRefineIterationsParam = 5000;
};

struct calculated_params	{
	int SWIDTH, SHEIGHT, BWIDTH, BHEIGHT, blocksize, numblocks, MaxA0, MaxB0;
	double a1_Max, a2_Max, a3_Max, b1_Max, b2_Max, b3_Max;
	double a1_step, a2_step, a3_step, b1_step, b2_step, b3_step;
	long totaldriftA, totaldriftB;

	int Asliverdrift, Bsliverdrift;
	int Adiffletpoints, Bdiffletpoints;
	int Apoints, Bpoints;
	int num_sliver_blocks;
	int sliver_block_width;
	unsigned long int dyn_diff_length;
	long int number_of_A_combos, number_of_B_combos;
};

struct images_store	{
	image_basic orig_sliver, orig_base, resamp_sliver, resamp_base;
};

struct results_params	{
	double bestA0, bestA1, bestA2, bestA3;
	double bestB0, bestB1, bestB2, bestB3;
	double bestC0, bestC1, bestC2, bestC3;
	
	// While performing diff summations, if diff goes past bestdiff, break and continue.
	double bestdiff = 1e37;
	
	unsigned long long count = 0;
	unsigned long long iterationsIgnored = 0;
};

struct beta_gamma_store	{
	gamma_values * difflets;
	beta_values * dynamic_diffs;
};

struct combos_store	{
	param_combo * A_combos;
	param_combo * B_combos;
};

/**
 * Basically just return floor() for negatives and ceil() for positives. Gosh why such a crappy implementation?
 * @param a
 * @return
 */
int rnd(double a) {
	if ( a >= 0 )
		return (int( a + 0.5 ));
	else
		return (int( a - 0.5 ));
}

/**
 * Returns < for param_combo based upon comparisons of P1, P2, and P3 with priorities in that order.
 *
 * @param a
 * @param b
 * @return The smaller param_combo.
 */
bool operator<(const param_combo& a, const param_combo& b) {
	if ( a.P1 != b.P1 )
		return (a.P1 < b.P1);
	else if ( a.P2 != b.P2 )
		return (a.P2 < b.P2);
	else
		return (a.P3 < b.P3);
}

/**
 * Returns > for param_combo based upon comparisons of P1, P2, and P3 with priorities in that order.
 *
 * @param a
 * @param b
 * @return The larger param_combo
 */
bool operator>(const param_combo& a, const param_combo& b) {
	if ( a.P1 != b.P1 )
		return (a.P1 > b.P1);
	else if ( a.P2 != b.P2 )
		return (a.P2 > b.P2);
	else
		return (a.P3 > b.P3);
}

/**
 * Returns <= for param_combo based upon comparisons of P1, P2, and P3 with priorities in that order.
 *
 * @param a
 * @param b
 * @return The smaller param_combo if different.
 */
bool operator<=(const param_combo& a, const param_combo& b) {
	if ( a.P1 != b.P1 )
		return (a.P1 < b.P1);
	else if ( a.P2 != b.P2 )
		return (a.P2 < b.P2);
	else
		return (a.P3 <= b.P3);
}

/**
 * Returns >= for param_combo based upon comparisons of P1, P2, and P3 with priorities in that order.
 *
 * @param a
 * @param b
 * @return The larger param_combo if different.
 */
bool operator>=(const param_combo& a, const param_combo& b) {
	if ( a.P1 != b.P1 )
		return (a.P1 > b.P1);
	else if ( a.P2 != b.P2 )
		return (a.P2 > b.P2);
	else
		return (a.P3 >= b.P3);
}

/**
 * Returns pixel value from float x,y by bilerping (could use trilinear algorithm...)
 * @param baseImage pointer to image
 * @param x coordinate
 * @param y coordinate
 * @return pixel value as float
 */
float interp_pixel_float(image_basic *baseImage, float x, float y) {

	float s, t;
	int image_width = baseImage->width;
	int image_height = baseImage->height;
	int left_val, right_val;
	float top_val, bottom_val;

	// Get integer coordinates for top left corner
	int left_index = (int) x, top_index = (int) y;
	int right_index;
	int bottom_index;

	// Get integer coordinates of bottom and right edges
	if ( left_index > image_width )
		right_index = left_index;             // Don't fall off the edge!
	else
		right_index = (int) ceil( x );             //left_index + 1;
	if ( top_index > image_height )
		bottom_index = top_index;             // Don't fall off the edge!
	else
		bottom_index = (int) ceil( y );             //top_index + 1;

	// Interpolate across the top edge
	s = x - (float) left_index;
	t = y - (float) top_index;
	left_val = baseImage->get( left_index, top_index );
	right_val = baseImage->get( right_index, top_index );
	// Linear interpolation recoded to only one multiply
	top_val = s * (float) right_val + ( 1 - s ) * (float) left_val; //right_val + s * (left_val-right_val);

	// Interpolate across the bottom edge
	left_val = baseImage->get( left_index, bottom_index );
	right_val = baseImage->get( right_index, bottom_index );
	//	bottom_val = right_val + s * (left_val-right_val);
	//Matt changed above to below
	bottom_val = s * (float) right_val + ( 1 - s ) * (float) left_val;

	// Interpolate between top and bottom
	return ( t * bottom_val + ( 1 - t ) * top_val );	//(bottom_val + t * (top_val-bottom_val));

}

/**
 * Returns pixel value from double x,y by bilerping. Has more precision and similar runtime as interp_pixel_float.
 * @param baseImage pointer to image
 * @param x coordinate
 * @param y coordinate
 * @return pixel value as double
 */
double interp_pixel(image_basic *baseImage, double x, double y) {

	double s, t;
	int image_width = baseImage->width;
	int image_height = baseImage->height;
	int left_val, right_val;
	double top_val, bottom_val;

	// Get integer coordinates for top left corner
	int left_index = (int) x, top_index = (int) y;
	int right_index;
	int bottom_index;

	// Get integer coordinates of bottom and right edges
	if ( left_index > image_width )
		right_index = left_index;             // Don't fall off the edge!
	else
		right_index = (int) ceil( x );             //left_index + 1;
	if ( top_index > image_height )
		bottom_index = top_index;             // Don't fall off the edge!
	else
		bottom_index = (int) ceil( y );             //top_index + 1;

	// Interpolate across the top edge
	s = x - (double) left_index;
	t = y - (double) top_index;
	left_val = baseImage->get( left_index, top_index );
	right_val = baseImage->get( right_index, top_index );
	// Linear interpolation recoded to only one multiply
	top_val = s * right_val + ( 1 - s ) * ( left_val ); //right_val + s * (left_val-right_val);

	// Interpolate across the bottom edge
	left_val = baseImage->get( left_index, bottom_index );
	right_val = baseImage->get( right_index, bottom_index );
	//	bottom_val = right_val + s * (left_val-right_val);
	//Matt changed above to below
	bottom_val = s * right_val + ( 1 - s ) * ( left_val );

	// Interpolate between top and bottom
	return ( t * bottom_val + ( 1 - t ) * top_val );	//(bottom_val + t * (top_val-bottom_val));

}

/**
 * Based upon given parameter values, the base image is warped.
 *
 * @param baseImage pointer to the image to warp
 * @param sizex width of base image
 * @param sizey height of base image
 * @param aterms A parameter values
 * @param bterms B parameter values
 * @param cterms C parameter values
 * @param dominantAxis Determines which axis (x or y) to apply the warping to.
 * @param warpedImage pointer to the image to write to
 */
void warp_image(image_basic *baseImage, int sizex, int sizey, double aterms[], double bterms[], double cterms[], int dominantAxis, image_basic *warpedImage) {
	//we have to retain the integrity of the old image, so we have to make a new, unwarped version.
	//image_basic warpedImage(sizex,sizey,1);
	int pixelx = 0;
	int pixely = 0;

	double newx = 0;
	double newy = 0;
	double newz = 0;
	int pixelvalue = 0;

	if ( dominantAxis == 1 ) {
		for ( pixely = 0; pixely < sizey; pixely++ ) {
			int pixely2 = pixely * pixely;
			int pixely3 = pixely2 * pixely;
			for ( pixelx = 0; pixelx < sizex; pixelx++ ) {
				//math.pow() is slightly slower than explicit multiplication.  Since this will be happening thousands of times, 
				//this IS slightly faster.  Please note that since
				//Nathan's edit (saving the value is even faster...moved calculation up to reduce redundancy by a factor of 1000.)

				newx = pixelx + aterms[ 0 ] + aterms[ 1 ] * pixely + aterms[ 2 ] * pixely2 + aterms[ 3 ] * pixely3;
				newy = bterms[ 0 ] + bterms[ 1 ] * pixely + bterms[ 2 ] * pixely2 + bterms[ 3 ] * pixely3;
				newz = cterms[ 0 ] + cterms[ 1 ] * pixely + cterms[ 2 ] * pixely2 + cterms[ 3 ] * pixely3;
				if ( newx < 0 || newx >= baseImage->width || newy < 0 || newy >= baseImage->height ) {
					warpedImage->set( pixelx, pixely, -1 );
				}
				else {
					pixelvalue = rnd( interp_pixel( baseImage, (float) newx, (float) newy ) + newz );     //where the new value is actually computed.
					//note,above, the cheap-ass implementation of rounding
					//ceil(newz);
					//pixelvalue+= newz;
					warpedImage->set( pixelx, pixely, pixelvalue );
				}
			}
		}
	}
	else {
		for ( pixelx = 0; pixelx < sizex; pixelx++ ) {
			int pixelx2 = pixelx * pixelx;
			int pixelx3 = pixelx2 * pixelx;
			for ( pixely = 0; pixely < sizey; pixely++ ) {
				// math.pow() is slightly slower than explicit multiplication.  Since this will be happening thousands of times
				// this IS slightly faster.  Please note that since
				// Nathan's edit. Again. see above. Reversed order of for loops 5/23/2011
				newx = aterms[ 0 ] + ( aterms[ 1 ] + 1 ) * pixelx + aterms[ 2 ] * pixelx2 + aterms[ 3 ] * pixelx3;
				newy = pixely + bterms[ 0 ] + ( bterms[ 1 ] - 1 ) * pixelx + bterms[ 2 ] * pixelx2 + bterms[ 3 ] * pixelx3;
				if ( newx < 0 || newx >= baseImage->width || newy < 0 || newy >= baseImage->height ) {
					warpedImage->set( pixelx, pixely, -1 );
				}
				else {
					pixelvalue = rnd( interp_pixel( baseImage, (float) newx, (float) newy ) );     //where the new value is actually computed.
					warpedImage->set( pixelx, pixely, pixelvalue );
				}
			}
		}
	}
}

/**
 * This method makes an array of pixelvalues, rather than the image_basic inhouse valarray usage, which, while fast, is also highly confusing.
 * @param baseImage pointer to the image
 * @param image pointer to two dimensional array
 */
void val2Array(image_basic *baseImage, int **image) {
	int sizex = baseImage->width;
	int sizey = baseImage->height;
	int pixelx = 0;
	int pixely = 0;

	for ( pixely = 0; pixely < sizey; pixely++ ) {
		for ( pixelx = 0; pixelx < sizex; pixelx++ ) {
			image[ pixelx ][ pixely ] = (int) interp_pixel( baseImage, (float) pixelx, (float) pixely );
			//note: since pixelx and pixely are integers, I have no idea why Brian is using interp_pixel for this.
		}
	}
}

/**
 * Compute the array index in an x-contiguous array.
 * @param sizex total size of the first dimension
 * @param x the index within the first dimension
 * @param y the index within the second dimension
 * @return the index within a one dimensional array representing two dimensions
 */
int compute1D(int sizex, int x, int y) {
	return ( x + sizex * y );
}

/**
 * This method makes an array of pixel values, rather than the image_basic in-house valarray usage, which, while fast, is also highly confusing.
 * @param baseImage
 * @param image
 */
void val21DArray(image_basic *baseImage, int image[]) {
	int sizex = baseImage->width;
	int sizey = baseImage->height;
	int pixelx = 0;
	int pixely = 0;

	for ( pixely = 0; pixely < sizey; pixely++ ) {
		for ( pixelx = 0; pixelx < sizex; pixelx++ ) {
			image[ compute1D( sizex, pixelx, pixely ) ] = baseImage->get( pixelx, pixely );
		}
	}
}

/**
 * Returns the nearest neighbor values
 * @param x current x
 * @param y current y
 * @param newx pointer to nearest neighbor value for x
 * @param newy pointer to nearest neighbor value for y
 */
void nearestneighbor(double x, double y, double* newx, double* newy) {
	*newx = rnd( x );
	*newy = rnd( y );
}

/**
 * This function down-samples an image by an integer factor.
 * Assumes values of images are integers
 * This is NOT a particularly fast function.
 * @param resamp
 * @param base
 * @param factor
 */
void resample(image_basic *resamp, image_basic base, int factor) {
	unsigned int rx_size = base.width / factor;
	unsigned int ry_size = base.height / factor;
	//resamp->set_color_mode(base.get_color_mode());
	resamp->initialize( rx_size, ry_size, base.get_color_mode() );
	//int cm = base.get_color_mode();
	//image_basic *resamp(rx_size,ry_size,base.get_color_mode()); 

	double inv_area = 1.0 / ( factor * factor );
	for ( unsigned int ry = 0; ry < ry_size; ++ry ) {
		for ( unsigned int rx = 0; rx < rx_size; ++rx ) {
			int pixel_value = 0;
			for ( unsigned int by = ry * factor; by < ( ry + 1 ) * factor; ++by )
				for ( unsigned int bx = rx * factor; bx < ( rx + 1 ) * factor; ++bx )
					pixel_value += base.fast_get( bx, by );
			resamp->set( rx, ry, rnd( pixel_value * inv_area ) );
		}
	}
}

#ifdef CODEFOLD

void logInputParams(input_params* p)	{
	printf("******** Input Parameters ********");
	printf("MaxA0 %d\n", p->MaxA0);
	printf("MaxB0 %d\n", p->MaxB0);
	printf("Maxa1 %f\n", p->a1_Max);
	printf("Maxb1 %f\n", p->b1_Max);
	printf("a2 %f\n", p->a2_Multiplier);
	printf("b2 %f\n", p->b2_Multiplier);
	printf("a3 %f\n", p->a3_Multiplier);
	printf("b3 %f\n", p->b3_Multiplier);
	printf("Precision %d\n", p->precision);
	printf("Growth %f\n", p->growthParam);
	printf("Contraction %f\n", p->contractParam);
	printf("Reflection %f\n", p->reflectParam);
	printf("Error Halt %f\n", p->haltParam);
	printf("Maximum Iterations %d\n", p->maxRefineIterationsParam);
	printf("\n");
}

void logCalculatedParams(calculated_params* c)	{
	printf("******** Calculated Parameters ********");
	printf("a1_Max = %f, a1_step = %f\n", c->a1_Max, c->a1_step);
	printf("a2_Max = %f, a2_step = %f\n", c->a2_Max, c->a2_step);
	printf("a3_Max = %f, a3_step = %f\n", c->a3_Max, c->a3_step);

	printf("b1_Max = %f, b1_step = %f\n", c->b1_Max, c->b1_step);
	printf("b2_Max = %f, b2_step = %f\n", c->b2_Max, c->b2_step);
	printf("b3_Max = %f, b3_step = %f\n", c->b3_Max, c->b3_step);

	printf("a1_num_steps = %f, a2_num_steps = %f, a3_num_steps = %f \n", c->a1_Max / c->a1_step, c->a2_Max / c->a2_step, c->a3_Max / c->a3_step);
	printf("b1_num_steps = %f, b2_num_steps = %f, b3_num_steps = %f \n", c->b1_Max / c->b1_step, c->b2_Max / c->b2_step, c->b3_Max / c->b3_step);
	printf("drifts: %d, %d\n", c->MaxA0, c->MaxB0);
	printf("total drifts: %ld, %ld\n", c->totaldriftA, c->totaldriftB);
	printf("height removed %ld\n", (c->totaldriftB + 2 * c->Bsliverdrift + c->MaxB0));
	printf("\n");
}

void logBetaGammaInfo(calculated_params* c)	{
	printf(
		"blocklet width: %d blockHeight (blocksize): %d number of blocklets (num_sliver_blocks):"
		" %d number of blocks (numblocks): %d size of beta array: %ld\n",
		c->SWIDTH / c->num_sliver_blocks, c->blocksize, c->num_sliver_blocks, c->numblocks, c->dyn_diff_length);
	printf("\n");
}

void logCurrentBest(results_params* r)	{
	printf("\n");
	printf("BEST: %e   COUNT: %I64u\n", r->bestdiff, r->count);
	printf("           %f, %f, %e, %e\n", r->bestA0, r->bestA1, r->bestA2, r->bestA3);
	printf("           %f, %f, %e, %e\n", r->bestB0, r->bestB1, r->bestB2, r->bestB3);
	printf("           %f, %f, %e, %e\n", r->bestC0, r->bestC1, r->bestC2, r->bestC3);
	printf("\n");
}

void logGridSearchInfo(times* t, results_params* r)	{
	printf("******** Finished Grid Search ********\n");
	printf("Time for Diffs : \t %f \n", t->diffs_time);
	printf("Time for Combos : \t %f \n", t->combos_time);
	printf("Time for Grid Search: \t %lf \n", t->grid_time);
	printf("Total Count : %I64u\n", r->count);
	printf("Iterations Ignored : %I64u\n", r->iterationsIgnored);
	printf("Total Difference : %f\n", r->bestdiff);
	printf("\n");
}

void logSimplexRoutineInfo(times* t, images_store* images, double* z)	{
	printf("******** Finished Simplex Routine ********\n");
	printf("Drift Coefficients :\n");
	printf("A0: %f\n", z[0]);
	printf("A1: %f\n", z[2]);
	printf("A2: %e\n", z[4]);
	printf("A3: %e\n", z[6]);
	printf("B0: %f\n", z[1]);
	printf("B1: %f\n", z[3] - 1);
	printf("B2: %e\n", z[5]);
	printf("B3: %e\n", z[7]);
	printf("C0: %e\nC1: %e\nC2: %e\nC3: %e\n", z[8], z[9], z[10], z[11]);
	printf("Time for Simplex Routine : \t %f seconds\n", t->simplex_time);
	printf("\n");
}

void logProgramInformation(double* z, double* x, times* t)	{
	FILE *paramfile = fopen("Output_Parameters.txt", "w");
	FILE *timefile = fopen("Runtimes.txt", "w");

	if (paramfile == NULL || timefile == NULL) {
		printf("File Not opened!");
	}
	else {
		fprintf(paramfile, "Simplex Routine Parameters : \n");
		fprintf(paramfile, "A Params\n");
		fprintf(paramfile, "%e %e %e %e\n", z[0], z[2], z[4], z[6]);
		fprintf(paramfile, "B Params\n");
		fprintf(paramfile, "%e %e %e %e\n", z[1], z[3] - 1.0, z[5], z[7]);
		fprintf(paramfile, "C Params\n");
		fprintf(paramfile, "%e %e %e %e\n", z[8], z[9], z[10], z[11]);
		fprintf(paramfile, "\n");
		fprintf(paramfile, "Grid Search Parameters : \n");
		fprintf(paramfile, "A Params\n");
		fprintf(paramfile, "%e %e %e %e\n", x[0], x[2], x[4], x[6]);
		fprintf(paramfile, "B Params\n");
		fprintf(paramfile, "%e %e %e %e\n", x[1], x[3] - 1.0, x[5], x[7]);
		fprintf(paramfile, "C Params\n");
		fprintf(paramfile, "%e %e %e %e\n", x[8], x[9], x[10], x[11]);
		fclose(paramfile);

		fprintf(timefile, "Image Read Time : \t%f\n", t->image_read_time);
		fprintf(timefile, "Combos Time : \t%f\n", t->combos_time);
		fprintf(timefile, "Diffs Time : \t%f\n", t->diffs_time);
		fprintf(timefile, "Grid Search Time : \t%f\n", t->grid_time);
		fprintf(timefile, "Simplex Routine Time : \t%f\n", t->simplex_time);
		fprintf(timefile, "Image Write Time\t%f\n", t->image_write_time);
		fprintf(timefile, "Total Program Time : \t%f\n", t->total_time);
		fclose(timefile);
		printf("********End of Program********\n");
	}
	getchar();  //puts a getchar in only if not run from commandline.
}

#endif

void readInputParams(int argc, char *argv[], input_params* p)	{

	if ( ( argc < 33 ) || ( argc > 33 ) ) {
		//perror("Error: Usage is ImageMaster -i: ImagePath -s: SliverPath -A0: value -A1: value -A2: value -B0: value -B1: value -B2: value -p: precision -b: blocksize -g: growth -c: contract -r: reflection\n");
		//ACTUALLY, there should be at least 31 arguments now, since I added cubic terms
		//Also, usage is -0 A0value, -1 B0value, -2 a1value -3 b1value -4 a2multiplier -5 b2multiplier, etc.
		printf( "Amount of arguments provided : %d\n", argc );
	}

	char command_arg;
	for ( int i = 1; i < argc; ++i ) {
		if ( argv[ i ][ 0 ] == '-' ) {
			command_arg = argv[ i ][ 1 ];

			// Look at argument after option to get corresponding value
			switch ( command_arg ) {
				case 'i':
					p->rundata1 = argv[ i + 1 ];
					for ( unsigned int x = 0; x < p->rundata1.length(); x++ ) {
						if ( p->rundata1[ x ] == '|' )
							p->rundata1[ x ] = ' ';
					}
					break;
				case 's':
					p->rundata2 = argv[ i + 1 ];
					for ( unsigned int x = 0; x < p->rundata2.length(); x++ ) {
						if ( p->rundata2[ x ] == '|' )
							p->rundata2[ x ] = ' ';
					}
					break;
				case '0':
					//A0percent = atof(argv[i+1]);
					p->MaxA0 = atoi( argv[ i + 1 ] );
					break;
				case '1':
					//A1percent = atof(argv[i+1]);
					p->MaxB0 = atoi( argv[ i + 1 ] );
					break;
				case '2':
					//A2percent = atof(argv[i+1]);
					p->a1_Max = atof( argv[ i + 1 ] );
					break;
				case '3':
					//B0percent = atof(argv[i+1]);
					p->b1_Max = atof( argv[ i + 1 ] );
					break;
				case '4':
					//B1percent = atof(argv[i+1]);
					p->a2_Multiplier = atof( argv[ i + 1 ] );
					break;
				case '5':
					//B2percent = atof(argv[i+1]);
					p->b2_Multiplier = atof( argv[ i + 1 ] );
					break;
				case '6':
					p->a3_Multiplier = atof( argv[ i + 1 ] );
					break;
				case '7':
					p->b3_Multiplier = atof( argv[ i + 1 ] );
					break;
				case 'p':
					p->precision = atoi( argv[ i + 1 ] );
					break;
				case 'b':
					p->blocksize = atoi( argv[ i + 1 ] );
					break;
				case 'g':
					p->growthParam = atof( argv[ i + 1 ] );
					break;
				case 'c':
					p->contractParam = atof( argv[ i + 1 ] );
					break;
				case 'r':
					p->reflectParam = atof( argv[ i + 1 ] );
					break;
				case 'h':
					p->haltParam = atof( argv[ i + 1 ] );
					break;
				case 't':
					p->maxRefineIterationsParam = atoi( argv[ i + 1 ] );
					break;
				default:
					break;
			}
		}
	}
	// If precision has not been specified, use maximum precision.
	if ( p->precision == 0 )	{
		p->precision = 1;
	}
}

void initCalculatedParams(calculated_params* c, input_params* p)	{
	// Set precision adjusted MaxA0 and MaxB0.
	c->MaxA0 = ( int ) ceil( ( float ) p->MaxA0 / ( float ) p->precision );
	c->MaxB0 = ( int ) ceil( ( float ) p->MaxB0 / ( float ) p->precision );

	// Set a1, a2, a3, b1, b2, b3 step sizes.
	c->a1_step = 1.0 / c->BHEIGHT; //that's 1 pixel
	double a2_Max = p->a2_Multiplier * ( p->a1_Max * 4 / c->BHEIGHT );  //one third of max deviation from linear part.  This is totally arbitrary.
	c->a2_step = 4.0 / ( c->BHEIGHT * c->BHEIGHT );
	double a3_Max = p->a3_Multiplier * 20.78 * p->a1_Max / ( c->BHEIGHT * c->BHEIGHT );
	c->a3_step = 20.78 / ( c->BHEIGHT * c->BHEIGHT * c->BHEIGHT ); //changed constant from 10.4 on 1/15/2007 <--- Nathan would like to know where magic 20.78 comes from.

	c->b1_step = 1.0 / c->BHEIGHT;
	double b2_Max = p->b2_Multiplier * ( p->b1_Max * 4 / c->BHEIGHT );
	c->b2_step = 4.0 / ( c->BHEIGHT * c->BHEIGHT );
	double b3_Max = p->b3_Multiplier * 20.78 * p->b1_Max / ( c->BHEIGHT * c->BHEIGHT );
	c->b3_step = 20.78 / ( c->BHEIGHT * c->BHEIGHT * c->BHEIGHT );

	// Make a1, a2, a3, b1, b2, b3 maximums exact ratios (rather than relative) of their respective steps
	c->a1_Max = ceil( p->a1_Max / c->a1_step ) * c->a1_step;
	c->a2_Max = ceil( a2_Max / c->a2_step ) * c->a2_step;
	c->a3_Max = ceil( a3_Max / c->a3_step ) * c->a3_step;

	c->b1_Max = ceil( p->b1_Max / c->b1_step ) * c->b1_step;
	c->b2_Max = ceil( b2_Max / c->b2_step ) * c->b2_step;
	c->b3_Max = ceil( b3_Max / c->b3_step ) * c->b3_step;

	// This set of variables record how many pixels of drift in either axis that will be encountered by the program.
	// Note that the factor of 1.38 is hard-coded for a2_multiplier = 1/3, a3_multiplier = 1/9
	c->totaldriftA = c->MaxA0 + ( int ) ceil( 1.38 * c->a1_Max * c->BHEIGHT );
	c->totaldriftB = c->MaxB0 + ( int ) ceil( 1.38 * c->b1_Max * c->BHEIGHT );

	// Compute maximum slope.
	double maxslope, maxslopeA, maxslopeB = 0;
	maxslopeA = c->a1_Max +
			c->a2_Max * c->BHEIGHT +
			0.5 * c->BHEIGHT * c->BHEIGHT * c->a3_Max;
	maxslopeB = c->b1_Max +
			c->b2_Max * c->BHEIGHT +
			0.5 * c->BHEIGHT * c->BHEIGHT * c->b3_Max;
	maxslope = maxslopeB > maxslopeA ? maxslopeB : maxslopeA;

	// Recalculate block size based upon calculated maximum slope if user does not provide a block size.
	if (p->blocksize <= 0)
	{
		c->blocksize = (int)(1 / maxslope);
	}
	else	{
		c->blocksize = p->blocksize;
	}

	// This adds one if there is a remainder;
	c->numblocks = 0; 

	c->Asliverdrift = ( int ) ( maxslopeA * c->SWIDTH );
	c->Bsliverdrift = ( int ) ( maxslopeB * c->SWIDTH );
	c->Adiffletpoints = ( c->totaldriftA + c->Asliverdrift ) * 2 + 1;  //get factor of two correct here.
	c->Bdiffletpoints = ( c->totaldriftB + c->Bsliverdrift ) * 2 + 1;
	
	c->numblocks = ( c->BHEIGHT - c->totaldriftB - 2 * c->Bsliverdrift - c->MaxB0 ) / c->blocksize;

	c->Apoints = c->totaldriftA * 2 + 1;
	c->Bpoints = c->totaldriftB * 2 + 1;
	c->num_sliver_blocks = c->SWIDTH / c->blocksize;
	c->sliver_block_width = c->SWIDTH / c->num_sliver_blocks;

	c->dyn_diff_length = c->numblocks * c->Apoints * c->Bpoints;
}

void readImages(images_store* images, input_params* p, calculated_params* c)	{

	// Read image files and resample based upon provided precision.
	images->orig_base.file_read_tiff( p->rundata1 );
	images->orig_sliver.file_read_tiff( p->rundata2 );
	resample( &images->resamp_base, images->orig_base, p->precision );
	resample( &images->resamp_sliver, images->orig_sliver, p->precision );

	// Set dimension parameters.
	c->BWIDTH = images->resamp_base.width;
	c->BHEIGHT = images->resamp_base.height;
	c->SWIDTH = images->resamp_sliver.width;
	c->SHEIGHT = images->resamp_sliver.height;
}

void initBetaGamma(calculated_params* c, images_store* images, beta_gamma_store* bg)	{
	int smin = c->BWIDTH / 2 - c->SWIDTH / 2;
	long long partial_sum1;
	long long partial_sum2;
	long long partial_sum3;
	long long diff;

	// Figure out maximum sliver drift, make difflets larger by that much.  (Apoints+Asliverpoints)
	bg->difflets = new gamma_values[c->numblocks * c->Adiffletpoints *
									c->Bdiffletpoints * c->num_sliver_blocks ];
	gamma_values difflet;
	for ( long int sb = 0; sb < c->num_sliver_blocks; ++sb ) {
		long int sliver_block_x = sb * c->sliver_block_width;
		for ( long int j = 0; j <= 2 * ( c->totaldriftB + c->Bsliverdrift ); ++j ) {
			for ( long int i = 0; i <= 2 * ( c->totaldriftA + c->Asliverdrift ); ++i ) {
				long int x_initial = smin - c->totaldriftA - c->Asliverdrift + i;
				x_initial += sliver_block_x;
				for ( long int b = 0; b < c->numblocks; b++ ) {
					partial_sum1 = 0;
					partial_sum2 = 0;
					partial_sum3 = 0;
					for ( long int y = c->MaxB0 - c->totaldriftB + j + b * c->blocksize,
							ys = c->MaxB0 + c->Bsliverdrift + b * c->blocksize;
							ys < c->MaxB0 + ( b + 1 ) * c->blocksize + c->Bsliverdrift;
							++y, ++ys ) {
						// Skip all negative elements.
						if ( y < 0 ) {
							continue;
						}
						for ( long int x = x_initial, xs = sliver_block_x;
								xs < sliver_block_x + c->sliver_block_width;
								++x, ++xs ) {
							diff = -images->resamp_sliver.fast_get( xs, ys ) + images->resamp_base.fast_get( x, y ) - 1;
							partial_sum1 += diff * diff;
							++diff;
							partial_sum2 += diff * diff;
							++diff;
							partial_sum3 += diff * diff;
						}
					}

					double y1 = ( double ) partial_sum1;
					double y2 = ( double ) partial_sum2;
					double y3 = ( double ) partial_sum3;

					double g2 = 0.5 * ( y1 + y3 ) - y2;
					double g1 = 0.5 * ( y3 - y1 );
					double g0 = y2;

					difflet.G0 = g0;
					difflet.G1 = g1;
					difflet.G2 = g2;

					bg->difflets[b + i * c->numblocks + j * c->numblocks * c->Adiffletpoints +
					          sb * c->numblocks * c->Adiffletpoints * c->Bdiffletpoints ] = difflet;
				} //end for b
			} //end for i
		} //end for j
	} //end for sb

	// Initialize array for dynamic_diffs once.
	bg->dynamic_diffs = new beta_values[c->dyn_diff_length];
}

/**
  * We use Matt's trick here of getting a list of indexes here, even though we have to use a for loop.  Pre-calculating the positions IS faster,
  * and Matt's trick precomputes as much as possible.
  */
void initCombos(calculated_params* c, combos_store* ab)	{
	// Make lists of all vaues of A1,A2,A3, B1,B2,B3 to try
	c->number_of_A_combos = ((int)(c->a1_Max / c->a1_step) * 2 + 1) *
		((int)(c->a2_Max / c->a2_step) * 2 + 1) *
		((int)(c->a3_Max / c->a3_step) * 2 + 1);
	ab->A_combos = new param_combo[c->number_of_A_combos];
	long int A_combos_counter = 0;
	for (double a3 = -c->a3_Max; a3 <= c->a3_Max + c->a3_step / 2; a3 += c->a3_step) {
		double p2_term1 = -1.5 * c->BHEIGHT * a3;
		double p1_term1 = 0.5 * c->BHEIGHT * c->BHEIGHT * a3;
		for (double a2 = -c->a2_Max; a2 <= c->a2_Max + c->a2_step / 2; a2 += c->a2_step) {
			double p1_term2 = -a2 * c->BHEIGHT;
			for (double a1 = -c->a1_Max; a1 <= c->a1_Max + c->a1_step / 2; a1 += c->a1_step) {
				ab->A_combos[A_combos_counter].P1 = a1 + p1_term2 + p1_term1;
				ab->A_combos[A_combos_counter].P2 = a2 + p2_term1;
				ab->A_combos[A_combos_counter].P3 = a3;
				++A_combos_counter;
			}
		}
	}

	c->number_of_B_combos = ((int)(c->b1_Max / c->b1_step) * 2 + 1) * 
								((int)(c->b2_Max / c->b2_step) * 2 + 1) * 
								((int)(c->b3_Max / c->b3_step) * 2 + 1);
	ab->B_combos = new param_combo[c->number_of_B_combos];
	long int B_combos_counter = 0;
	for ( double b3 = -c->b3_Max; b3 <= c->b3_Max + c->b3_step / 2; b3 += c->b3_step ) {
		double p2_term1 = - 1.5 * c->BHEIGHT * b3;
		double p1_term1 = 0.5 * c->BHEIGHT * c->BHEIGHT * b3;
		for ( double b2 = -c->b2_Max; b2 <= c->b2_Max + c->b2_step / 2; b2 += c->b2_step ) {
			double p1_term2 = - b2 * c->BHEIGHT;
			for ( double b1 = -c->b1_Max; b1 <= c->b1_Max + c->b1_step / 2; b1 += c->b1_step ) {
				ab->B_combos[B_combos_counter].P1 = b1 + p1_term2 + p1_term1;
				ab->B_combos[B_combos_counter].P2 = b2 + p2_term1;
				ab->B_combos[B_combos_counter].P3 = b3;
				++B_combos_counter;
			}
		}
	}
	std::sort(ab->A_combos, (ab->A_combos + (c->number_of_A_combos)));
	std::sort(ab->B_combos, (ab->B_combos + (c->number_of_B_combos)));
}

void performGridSearch(calculated_params* c, results_params* r, combos_store* ab, beta_gamma_store* bg)	{

	double * yb_arr = new double[c->numblocks];
	double * yb2_arr = new double[c->numblocks];
	double * yb3_arr = new double[c->numblocks];
	double * yb4_arr = new double[c->numblocks];
	double * yb5_arr = new double[c->numblocks];
	double * yb6_arr = new double[c->numblocks];
	int mult1, mult2, mult3;
	mult3 = c->numblocks * c->Adiffletpoints * c->Bdiffletpoints;
	mult2 = c->numblocks * c->Adiffletpoints;
	mult1 = c->numblocks * c->Apoints;

	for (int i = 0; i < c->numblocks; ++i) {
		yb_arr[i] = i * c->blocksize + c->MaxB0; //+ Bsliverdrift;
		yb2_arr[i] = yb_arr[i] * yb_arr[i];
		yb3_arr[i] = yb2_arr[i] * yb_arr[i];
		yb4_arr[i] = yb3_arr[i] * yb_arr[i];
		yb5_arr[i] = yb4_arr[i] * yb_arr[i];
		yb6_arr[i] = yb5_arr[i] * yb_arr[i];
	}

	//TO DO for tomorrow: think about how to do initial and final B1 sensibly and symmetrically
	//Incorporate actual sliver_width into calculation, no hard-wiring to 32.
	//Refill DIFFS
	long int Aindex_f = 0, Bindex_f = 0;
	for (long int Bindex_i = 0; Bindex_i < c->number_of_B_combos - 1; Bindex_i = Bindex_f) {
		double B1_i = ab->B_combos[Bindex_i].P1;
		double B1_f = B1_i + .01;	//(1.0/32);
		//search ahead for index Bi_f
		for (Bindex_f = Bindex_i + 1; (ab->B_combos[Bindex_f].P1 < B1_f) && (Bindex_f < c->number_of_B_combos); ++Bindex_f);

		for (long int Aindex_i = 0; Aindex_i < c->number_of_A_combos - 1; Aindex_i = Aindex_f) {
			double A1_i = ab->A_combos[Aindex_i].P1;
			double A1_f = A1_i + .01;		//(1.0/32);
			//search ahead for index Ai_f
			for (Aindex_f = Aindex_i + 1; (ab->A_combos[Aindex_f].P1 < A1_f) && (Aindex_f < c->number_of_A_combos); ++Aindex_f);
			
			// Empty out array of beta values.
			memset(bg->dynamic_diffs, 0, c->dyn_diff_length * sizeof(beta_values));

			// Populate array with new beta values.
			long int diffs_address;
			long int difflets_address;
			int Xc = -c->sliver_block_width;
			int totX;
			int totY;
			int blocktimesApoints = c->Apoints * c->numblocks;
			long A0_B0_adj_inc = c->numblocks;  //A0increase is now assumed to be 1.
			long MaxA0_times_numblocks = c->MaxA0 * c->numblocks;
			long B0_adj_min = -c->MaxB0 * blocktimesApoints;
			long B0_adj_max = c->MaxB0 * blocktimesApoints;
			//long BO_adj_inc = B0increase*numblocks*Apoints;
			long BO_adj_inc = blocktimesApoints;  //B0increase is now assumed to be 1.
			int finalAddress = -mult3;

			for (int sb = 0; sb < c->num_sliver_blocks; ++sb) {
				//Now calculate offsets sx and symeh
				finalAddress += mult3;
				Xc += c->sliver_block_width;
				int Xc2 = Xc * Xc;
				int doubleXc = 2 * Xc;
				int sx = c->Asliverdrift - (int)(Xc * A1_i);  //which A1 do I use?  first?  last?  average?
				int sy = c->Bsliverdrift - (int)(Xc * B1_i);
				//printf("%f, %d, %d, %d\n",A1_i,sb,sx,sy);
				int YPart = -mult1;
				for (int dy = 0; dy < c->Bpoints; ++dy) {
					totY = ( dy + sy ) * mult2;
					YPart += mult1;
					int XPart = -c->numblocks;
					for (int dx = 0; dx < c->Apoints; ++dx) {
						totX = (dx + sx) * c->numblocks;
						XPart += (int)c->numblocks;
						for (int b = 0; b < c->numblocks; ++b) {  //partially unroll this loop?
							//diffs[dy,dx,b] += difflets[sb,dy+sy,dx+sx,b]

							//Xc=0;
							//printf("X position: %f \n",Xc);

							//difflets_address =(b + (dx+sx)*numblocks + (dy+sy)*numblocks*Adiffletpoints + sb*numblocks*Adiffletpoints*Bdiffletpoints);
							difflets_address = ( b + totX + totY + finalAddress );
							//diffs_address = b+dx*numblocks+dy*numblocks*Apoints;
							diffs_address = b + XPart + YPart;
							bg->dynamic_diffs[diffs_address].B0 += bg->difflets[difflets_address].G0;
							bg->dynamic_diffs[diffs_address].B1 += bg->difflets[difflets_address].G1;
							bg->dynamic_diffs[diffs_address].B2 += bg->difflets[difflets_address].G2;
							bg->dynamic_diffs[diffs_address].B3 += (bg->difflets[difflets_address].G1 * Xc);	// Adjust for new Beta3 which is sum of negatives
							bg->dynamic_diffs[diffs_address].B4 += bg->difflets[difflets_address].G2 * Xc2;
							bg->dynamic_diffs[diffs_address].B5 += (bg->difflets[difflets_address].G2 * doubleXc);// Adjust for new Beta5 which is sum of negatives
							//}
						}
					}
				}
			}

			critical_section cs;
			parallel_for((long int)Bindex_i, Bindex_f, (long int)1, [&](long int Bindex)
			{
				for (long int Aindex = Aindex_i; Aindex < Aindex_f; ++Aindex) {
					double B1 = ab->B_combos[Bindex].P1;
					double B2 = ab->B_combos[Bindex].P2;
					double B3 = ab->B_combos[Bindex].P3;
					double A1 = ab->A_combos[Aindex].P1;
					double A2 = ab->A_combos[Aindex].P2;
					double A3 = ab->A_combos[Aindex].P3;

					for (long B0_adj = B0_adj_min; B0_adj <= B0_adj_max; B0_adj += BO_adj_inc) {
						long A0_B0_adj_max = MaxA0_times_numblocks + B0_adj;

						for (long A0_B0_adj = -MaxA0_times_numblocks + B0_adj; A0_B0_adj <= A0_B0_adj_max; A0_B0_adj += A0_B0_adj_inc) {

							Eigen::Matrix4d A;
							Eigen::MatrixXd B(4, 1);
							Eigen::MatrixXd C(4, 1);
							Eigen::Matrix4d D;
							Eigen::MatrixXd E(4, 1);

							A << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
							D << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
							B << 0, 0, 0, 0;
							C << 0, 0, 0, 0;
							E << 0, 0, 0, 0;

							for (int i = 0; i<c->numblocks; i++) {
								double beta1, double_beta2, beta3, beta4, beta5;

								long xpos = (long)(A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + c->totaldriftA);
								long ypos = (long)(B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + c->totaldriftB);

								int index = i + xpos*c->numblocks + ypos*blocktimesApoints + A0_B0_adj;

								double yb = yb_arr[i];
								double yb2 = yb2_arr[i];
								double yb3 = yb3_arr[i];
								double yb4 = yb4_arr[i];
								double yb5 = yb5_arr[i];
								double yb6 = yb6_arr[i];

								beta1 = bg->dynamic_diffs[index].B1;
								double_beta2 = bg->dynamic_diffs[index].B2 * 2;
								beta3 = bg->dynamic_diffs[index].B3;
								beta4 = bg->dynamic_diffs[index].B4;
								beta5 = bg->dynamic_diffs[index].B5;

								D << double_beta2, double_beta2*yb + (beta5), double_beta2*yb2, double_beta2*yb3,//1st row
										double_beta2*yb + (beta5), double_beta2*yb2 + (2 * beta4 + 2 * beta5*yb), double_beta2*yb3 + (beta5*yb2), double_beta2*yb4 + (beta5*yb3),//2nd Row
										double_beta2*yb2, double_beta2*yb3 + (beta5*yb2), double_beta2*yb4, double_beta2*yb5,//3rd Row
										double_beta2*yb3, double_beta2*yb4 + (beta5*yb3), double_beta2*yb5, double_beta2*yb6;//4th row
								E << -beta1,
										-beta1 * yb - beta3,
										-beta1 * yb2,
										-beta1 * yb3;
								A += D;
								B += E;
							}

							C = A.inverse()*B;
							double C0 = C.data()[0];
							double C1 = C.data()[1];
							double C2 = C.data()[2];
							double C3 = C.data()[3];
							

							// Critical section lock to prevent data overwriting.
							cs.lock();
							double function = 0;
							for (int i = 0; i<c->numblocks; i++){
								long xpos = (long)(A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + c->totaldriftA);
								long ypos = (long)(B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + c->totaldriftB);
								//int index = diffs_indices[i] + A0_B0_adj;//*numBytes
								int index = i + xpos*c->numblocks + ypos*blocktimesApoints + A0_B0_adj;
								//int index = diffs_indices[i] + A0_B0_adj;
								double yb = yb_arr[i];
								double yb2 = yb2_arr[i];
								double yb3 = yb3_arr[i];

								double Beta0 = bg->dynamic_diffs[index].B0;
								double Beta1 = bg->dynamic_diffs[index].B1;
								double Beta2 = bg->dynamic_diffs[index].B2;
								double Beta3 = bg->dynamic_diffs[index].B3;
								double Beta4 = bg->dynamic_diffs[index].B4;
								double Beta5 = bg->dynamic_diffs[index].B5;
								double delz = C0 + C1*yb + C2*yb2 + C3*yb3;

								double sum = Beta0 + Beta1*delz + Beta2*delz*delz
									+ Beta3*C1 + Beta4*C1*C1 + Beta5*C1*delz;

								if (sum < 0)	{
									printf("Negative Sum!\n");
								}

								function += sum;
								if (function > r->bestdiff)	{
									r->iterationsIgnored += c->numblocks - (i + 1);
									break;
								}
							}
							// Note: this "if" block doesn't seem to slow down the program at all, compared to simpler statements.
							if (function <= r->bestdiff){
								r->bestdiff = function;
								r->bestA0 = (A0_B0_adj - B0_adj) / c->numblocks;
								r->bestB0 = B0_adj / blocktimesApoints;
								r->bestA1 = A1;
								r->bestA2 = A2;
								r->bestA3 = A3;
								r->bestB1 = B1;
								r->bestB2 = B2;
								r->bestB3 = B3;
								r->bestC0 = C0;
								r->bestC1 = C1;
								r->bestC2 = C2;
								r->bestC3 = C3;

								if (LOG_OUTPUT)	{
									logCurrentBest(r);
								}
							}
							r->count++;
							// Unlock here.
							cs.unlock();
						}
					}
				}
			} );
		}
	}
}

void performSimplexRoutine(input_params* pparams, results_params* rparams, images_store* images, double* x, double* z, double* precisionArr)	{
	rparams->bestC1 /= pparams->precision;
	rparams->bestC2 /= pparams->precision * pparams->precision;
	rparams->bestC3 /= pparams->precision * pparams->precision * pparams->precision;

	// These are how you should fill the x array. Use this initialization 
	// if you haven't skipped the grid search
	x[0] = rparams->bestA0 * pparams->precision;
	x[1] = rparams->bestB0 * pparams->precision;
	x[2] = rparams->bestA1;
	x[3] = rparams->bestB1 + 1;
	x[4] = rparams->bestA2 / pparams->precision;
	x[5] = rparams->bestB2 / pparams->precision;
	x[6] = rparams->bestA3 / (pparams->precision * pparams->precision);
	x[7] = rparams->bestB3 / (pparams->precision * pparams->precision);
	x[8] = rparams->bestC0;
	x[9] = rparams->bestC1;
	x[10] = rparams->bestC2;
	x[11] = rparams->bestC3;

	// Setup precision for simplex.  
	// These values are all based on the changes in the parameters that will cause at most a shift of 
	// 0.1 pixel.
	for (int i = 0; i < 12; ++i)	{
		z[i] = x[i];
	}
	
	precisionArr[0] = 0.1;
	precisionArr[1] = 0.1;
	precisionArr[2] = 0.1 / images->orig_base.height;
	precisionArr[3] = 0.1 / images->orig_base.height;
	precisionArr[4] = 0.1 / pow(images->orig_base.height, 2.0);
	precisionArr[5] = 0.1 / pow(images->orig_base.height, 2.0);
	precisionArr[6] = 0.1 / pow(images->orig_base.height, 3.0);
	precisionArr[7] = 0.1 / pow(images->orig_base.height, 3.0);
	precisionArr[8] = 0.1;
	precisionArr[9] = 0.1 / images->orig_base.height;
	precisionArr[10] = 0.1 / pow(images->orig_base.height, 2.0);
	precisionArr[11] = 0.1 / pow(images->orig_base.height, 3.0);

	simplex(&images->orig_base, &images->orig_sliver, x, z, 12, precisionArr, pparams->reflectParam,
		pparams->contractParam, pparams->growthParam, pparams->haltParam,
		pparams->maxRefineIterationsParam);
}

void performImageCorrection(double* z, images_store* images, input_params* pparams)	{
	// Perform final warp, and write the output tiff
	double aterms[4] = { z[0], z[2], z[4], z[6] };
	double bterms[4] = { z[1], z[3], z[5], z[7] };
	double cterms[4] = { z[8], z[9], z[10], z[11] };

	image_basic final(images->orig_base.width, images->orig_base.height, images->orig_base.get_color_mode());
	warp_image(&images->orig_base, images->orig_base.width, images->orig_base.height, aterms, bterms, cterms, 1, &final);

	std::string finalstring = pparams->rundata1;
	std::basic_string < char > finalmarkerstring("_corrected");
	finalstring.insert(finalstring.rfind("."), finalmarkerstring);
	final.file_write_tiff(finalstring);
}

/**
 * Main method for our program.
 */
int main(int argc, char *argv[]) {
	clock_t begintime = clock();

	bool debugmode = false;
	LOG_OUTPUT = true;

	// Various clock variables for timing.
	clock_t time0, time1;

	input_params* pparams = new input_params;
	calculated_params* cparams = new calculated_params;
	images_store* images = new images_store;
	results_params* rparams = new results_params;
	beta_gamma_store* bgparams = new beta_gamma_store;
	combos_store* abparams = new combos_store;
	times* t = new times;


	time0 = clock();
	readInputParams(argc, argv, pparams);
	if ( debugmode ) {
		pparams->precision = 1;
		pparams->blocksize = 8;
	}
	readImages(images, pparams, cparams);
	initCalculatedParams(cparams, pparams);
	time1 = clock();
	t->image_read_time = ((double)time1 - (double)time0) / 1000;

	if (LOG_OUTPUT)	{
		logInputParams(pparams);
		logCalculatedParams(cparams);
	}

	
	time0 = clock();
	initBetaGamma(cparams, images, bgparams);
	time1 = clock();
	t->diffs_time = ((double)time1 - (double)time0) / 1000;

	if (LOG_OUTPUT)	{
		logBetaGammaInfo(cparams);
	}

	time0 = clock();
	initCombos(cparams, abparams);
	time1 = clock();
	t->combos_time = ((double)time1 - (double)time0) / 1000;

	time0 = clock();
	performGridSearch(cparams, rparams, abparams, bgparams);
	time1 = clock();
	t->grid_time = ((double)time1 - (double)time0) / 1000;

	if (LOG_OUTPUT)	{
		logGridSearchInfo(t, rparams);
	}

	// Number of parameters to be optimized 
	int n = 12;
	double x[12];
	double z[12];
	double precisionArr[12];

	time0 = clock();
	performSimplexRoutine(pparams, rparams, images, x, z, precisionArr);
	time1 = clock();
	t->simplex_time = ((double)time1 - (double)time0) / 1000;

	if (LOG_OUTPUT)	{
		logSimplexRoutineInfo(t, images, z);
	}

	time0 = clock();
	performImageCorrection(z, images, pparams);
	time1 = clock();
	t->image_write_time = ((double)time1 - (double)time0) / 1000;

	int endtime = clock();
	t->total_time = (double)(endtime - begintime) * .001;

	if (LOG_OUTPUT)	{
		logProgramInformation(z, x, t);
	}
}
