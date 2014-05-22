/*********************************************************************
 *Image Master procedural image registration and stabilization routine*
 *Coded by: Brian Salmons & Dr. Matthew Trawick                       *
 **********************************************************************
 *Dr. Matthew Trawick                                                 *
 *University of Richmond                                              *
 **********************************************************************
 */

/** Version 1.2
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
#include "argo.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <valarray>
#include <algorithm>
#include <tbb\parallel_for.h>
#include <tbb\critical_section.h>
#include "Eigen\Eigen"
#include "simplex.h"

//#define SANITYCHECKS
#define VERSION "1.2"
#define CODEFOLD

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

argo::argo()	{
	blocksize = 0;
	MaxA0 = 5, MaxB0 = 5;
	a1_Max = 0.1, b1_Max = 0.1;
	a2_Multiplier = (1.0 / 9.0), b2_Multiplier = (1.0 / 9.0);
	a3_Multiplier = .1, b3_Multiplier = .1;
	precision = 1;

	/**
	* These parameters were found, using trial and error on not a whole lot of different test cases.
	* To avoid being trapped in some rut, it seems to be VERY important to have the reflection
	* parameter slightly less than 1.0, and the growth parameter less than 2.0.
	*/
	growthParam = 1.5, contractParam = 0.5, reflectParam = 0.9, haltParam = 1e-11;
	maxRefineIterationsParam = 5000;

	rsliver_width, rsliver_height, rbase_width, rbase_height, numblocks;
	a2_Max, a3_Max, b2_Max, b3_Max;
	a1_step, a2_step, a3_step, b1_step, b2_step, b3_step;
	totaldriftA, totaldriftB;

	Asliverdrift, Bsliverdrift;
	Adiffletpoints, Bdiffletpoints;
	Apoints, Bpoints;
	num_sliver_blocks;
	sliver_block_width;
	dyn_diff_length;
	number_of_A_combos, number_of_B_combos;

	results.bestdiff = 1e37;
	results.count = 0;
	results.iterationsIgnored = 0;
}

argo::~argo()	{
	if (beta_gamma_store.difflets)	{
		delete[] beta_gamma_store.difflets;
	}
	if (beta_gamma_store.dynamic_diffs)	{
		delete[] beta_gamma_store.dynamic_diffs;
	}
	if (combos_store.A_combos)	{
		delete[] combos_store.A_combos;
	}
	if (combos_store.B_combos)	{
		delete[] combos_store.B_combos;
	}
}

#ifdef CODEFOLD

void argo::logInputParams()	{
	printf("******** Input Parameters ********\n");
	printf("Input Max A0 : \t %d\n", MaxA0);
	printf("Input Max B0 : \t %d\n", MaxB0);
	printf("Input Max a1 : \t %f\n", a1_Max);
	printf("Input Max b1 : \t %f\n", b1_Max);
	printf("Input a2 Multiplier : \t %f\n", a2_Multiplier);
	printf("Input b2 Multiplier : \t %f\n", b2_Multiplier);
	printf("Input a3 Multiplier : \t %f\n", a3_Multiplier);
	printf("Input b3 Multiplier : \t %f\n", b3_Multiplier);

	printf("******** Simplex Routine Information ********\n");
	printf("Input Precision : \t %d\n", precision);
	printf("Input Growth : \t %f\n", growthParam);
	printf("Input Contraction : \t %f\n", contractParam);
	printf("Input Reflection : \t %f\n", reflectParam);
	printf("Input Error Halt : \t %f\n", haltParam);
	printf("Input Maximum Iterations : \t %d\n", maxRefineIterationsParam);
	printf("\n");
}

void argo::logCalculatedParams()	{
	printf("******** Calculated Parameters ********\n");
	printf("a1 Max : %f, a1 Step Size : %f, a1 Steps : %f\n", a1_Max, a1_step, a1_Max / a1_step);
	printf("a2 Max : %f, a2 Step Size : %f, a2 Steps : %f\n", a2_Max, a2_step, a2_Max / a2_step);
	printf("a3 Max : %f, a3 Step Size : %f, a3 Steps : %f\n", a3_Max, a3_step, a3_Max / a3_step);

	printf("b1 Max : %f, b1 Step Size : %f, b1 Steps : %f\n", b1_Max, b1_step, b1_Max / b1_step);
	printf("b2 Max : %f, b2 Step Size : %f, b2 Steps : %f\n", b2_Max, b2_step, b2_Max / b2_step);
	printf("b3 Max : %f, b3 Step Size : %f, b3 Steps : %f\n", b3_Max, b3_step, b3_Max / b3_step);

	printf("A0 Max : %d, B0 Max %d\n", MaxA0, MaxB0);
	printf("total drifts : \t %ld, %ld\n", totaldriftA, totaldriftB);
	printf("height removed : \t %ld\n", (totaldriftB + 2 * Bsliverdrift + MaxB0));
	printf("\n");
}

void argo::logComboInfo()	{

}

void argo::logBetaGammaInfo()	{
	printf("Number of Blocks : \t %d\n", numblocks);
	printf("Number of Blocklets : \t %d\n", num_sliver_blocks);
	printf("Blocklet Pixel Width : \t %d\n", rsliver_width / num_sliver_blocks);
	printf("Block and Blocklet Pixel Height : \t %d\n", blocksize);
	printf("Size of Beta Array : \t %d\n", dyn_diff_length);
	printf("\n");
}

void argo::logCurrentBest()	{
	printf("\n");
	printf("BEST: %e   COUNT: %I64u\n", results.bestdiff, results.count);
	printf("           %f, %f, %e, %e\n", results.bestA0, results.bestA1, results.bestA2, results.bestA3);
	printf("           %f, %f, %e, %e\n", results.bestB0, results.bestB1, results.bestB2, results.bestB3);
	printf("           %f, %f, %e, %e\n", results.bestC0, results.bestC1, results.bestC2, results.bestC3);
	printf("\n");
}

void argo::logGridSearchInfo()	{
	printf("******** Finished Grid Search ********\n");
	printf("Time for Diffs : \t %f \n", times_store.diffs_time);
	printf("Time for Combos : \t %f \n", times_store.combos_time);
	printf("Time for Grid Search: \t %lf \n", times_store.grid_time);
	printf("Total Count : \t %I64u\n", results.count);
	printf("Iterations Ignored : \t %I64u\n", results.iterationsIgnored);
	printf("Total Difference : \t %f\n", results.bestdiff);
	printf("\n");
}

void argo::logSimplexRoutineInfo()	{
	printf("******** Finished Simplex Routine ********\n");
	printf("Drift Coefficients :\n");
	printf("Best A0 : %f\n", z[0]);
	printf("Best A1 : %f\n", z[2]);
	printf("Best A2 : %e\n", z[4]);
	printf("Best A3 : %e\n", z[6]);
	printf("Best B0 : %f\n", z[1]);
	printf("Best B1 : %f\n", z[3] - 1);
	printf("Best B2 : %e\n", z[5]);
	printf("Best B3 : %e\n", z[7]);
	printf("Best C0 : %e\nBest C1: %e\nBest C2: %e\nBest C3: %e\n", z[8], z[9], z[10], z[11]);
	printf("Time for Simplex Routine : \t %f seconds\n", times_store.simplex_time);
	printf("\n");
}

void argo::logProgramInformation()	{
	FILE *paramfile = fopen("Output_Parameters.txt", "w");

	if (paramfile == NULL) {
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

		fprintf(paramfile, "Image Read Time : \t%f\n", times_store.image_read_time);
		fprintf(paramfile, "Combos Time : \t%f\n", times_store.combos_time);
		fprintf(paramfile, "Diffs Time : \t%f\n", times_store.diffs_time);
		fprintf(paramfile, "Grid Search Time : \t%f\n", times_store.grid_time);
		fprintf(paramfile, "Simplex Routine Time : \t%f\n", times_store.simplex_time);
		fprintf(paramfile, "Image Write Time\t%f\n", times_store.image_write_time);
		fprintf(paramfile, "Total Program Time : \t%f\n", times_store.total_time);

		fclose(paramfile);
		printf("********End of Program********\n");
	}
}

#endif

void argo::readImages()	{

	// Read image files and resample based upon provided precision.
	images_store.orig_base.file_read_tiff(rundata1);
	images_store.orig_sliver.file_read_tiff(rundata2);
	resample(&images_store.resamp_base, images_store.orig_base, precision);
	resample(&images_store.resamp_sliver, images_store.orig_sliver, precision);

	// Set dimension parameters.
	rbase_width = images_store.resamp_base.width;
	rbase_height = images_store.resamp_base.height;
	rsliver_width = images_store.resamp_sliver.width;
	rsliver_height = images_store.resamp_sliver.height;
}

void argo::readInputParams(int argc, char *argv[])	{

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
					rundata1 = argv[ i + 1 ];
					for ( unsigned int x = 0; x < rundata1.length(); x++ ) {
						if ( rundata1[ x ] == '|' )
							rundata1[ x ] = ' ';
					}
					break;
				case 's':
					rundata2 = argv[ i + 1 ];
					for ( unsigned int x = 0; x < rundata2.length(); x++ ) {
						if ( rundata2[ x ] == '|' )
							rundata2[ x ] = ' ';
					}
					break;
				case '0':
					MaxA0 = atoi( argv[ i + 1 ] );
					break;
				case '1':
					MaxB0 = atoi( argv[ i + 1 ] );
					break;
				case '2':
					a1_Max = atof( argv[ i + 1 ] );
					break;
				case '3':
					b1_Max = atof( argv[ i + 1 ] );
					break;
				case '4':
					a2_Multiplier = atof( argv[ i + 1 ] );
					break;
				case '5':
					b2_Multiplier = atof( argv[ i + 1 ] );
					break;
				case '6':
					a3_Multiplier = atof( argv[ i + 1 ] );
					break;
				case '7':
					b3_Multiplier = atof( argv[ i + 1 ] );
					break;
				case 'p':
					precision = atoi( argv[ i + 1 ] );
					break;
				case 'b':
					blocksize = atoi( argv[ i + 1 ] );
					break;
				case 'g':
					growthParam = atof( argv[ i + 1 ] );
					break;
				case 'c':
					contractParam = atof( argv[ i + 1 ] );
					break;
				case 'r':
					reflectParam = atof( argv[ i + 1 ] );
					break;
				case 'h':
					haltParam = atof( argv[ i + 1 ] );
					break;
				case 't':
					maxRefineIterationsParam = atoi( argv[ i + 1 ] );
					break;
				default:
					break;
			}
		}
	}
	// If precision has not been specified, use maximum precision.
	if ( precision == 0 )	{
		precision = 1;
	}
}

void argo::initCalculatedParams()	{
	// Set precision adjusted MaxA0 and MaxB0.
	MaxA0 = ( int ) ceil( ( float ) MaxA0 / ( float ) precision );
	MaxB0 = ( int ) ceil( ( float ) MaxB0 / ( float ) precision );

	// Set a1, a2, a3, b1, b2, b3 step sizes.
	a1_step = 1.0 / rbase_height; //that's 1 pixel
	a2_Max = a2_Multiplier * ( a1_Max * 4 / rbase_height );  //one third of max deviation from linear part.  This is totally arbitrary.
	a2_step = 4.0 / ( rbase_height * rbase_height );
	a3_Max = a3_Multiplier * 20.78 * a1_Max / ( rbase_height * rbase_height );
	a3_step = 20.78 / ( rbase_height * rbase_height * rbase_height ); //changed constant from 10.4 on 1/15/2007 <--- Nathan would like to know where magic 20.78 comes from.

	b1_step = 1.0 / rbase_height;
	b2_Max = b2_Multiplier * ( b1_Max * 4 / rbase_height );
	b2_step = 4.0 / ( rbase_height * rbase_height );
	b3_Max = b3_Multiplier * 20.78 * b1_Max / ( rbase_height * rbase_height );
	b3_step = 20.78 / ( rbase_height * rbase_height * rbase_height );

	// Make a1, a2, a3, b1, b2, b3 maximums exact ratios (rather than relative) of their respective steps
	a1_Max = ceil( a1_Max / a1_step ) * a1_step;
	a2_Max = ceil( a2_Max / a2_step ) * a2_step;
	a3_Max = ceil( a3_Max / a3_step ) * a3_step;

	b1_Max = ceil( b1_Max / b1_step ) * b1_step;
	b2_Max = ceil( b2_Max / b2_step ) * b2_step;
	b3_Max = ceil( b3_Max / b3_step ) * b3_step;

	// This set of variables record how many pixels of drift in either axis that will be encountered by the program.
	// Note that the factor of 1.38 is hard-coded for a2_multiplier = 1/3, a3_multiplier = 1/9
	totaldriftA = MaxA0 + ( int ) ceil( 1.38 * a1_Max * rbase_height );
	totaldriftB = MaxB0 + ( int ) ceil( 1.38 * b1_Max * rbase_height );

	// Compute maximum slope.
	double maxslope, maxslopeA, maxslopeB = 0;
	maxslopeA = a1_Max +
			a2_Max * rbase_height +
			0.5 * rbase_height * rbase_height * a3_Max;
	maxslopeB = b1_Max +
			b2_Max * rbase_height +
			0.5 * rbase_height * rbase_height * b3_Max;
	maxslope = maxslopeB > maxslopeA ? maxslopeB : maxslopeA;

	// Recalculate block size based upon calculated maximum slope if user does not provide a block size.
	if (blocksize <= 0)
	{
		blocksize = (int)(1 / maxslope);
	}

	// This adds one if there is a remainder;
	numblocks = 0; 

	Asliverdrift = ( int ) ( maxslopeA * rsliver_width );
	Bsliverdrift = ( int ) ( maxslopeB * rsliver_width );

	Adiffletpoints = ( totaldriftA + Asliverdrift ) * 2 + 1;  //get factor of two correct here.
	Bdiffletpoints = ( totaldriftB + Bsliverdrift ) * 2 + 1;
	
	numblocks = ( rbase_height - totaldriftB - 2 * Bsliverdrift - MaxB0 ) / blocksize;

	Apoints = totaldriftA * 2 + 1;
	Bpoints = totaldriftB * 2 + 1;
	num_sliver_blocks = rsliver_width / blocksize;
	sliver_block_width = rsliver_width / num_sliver_blocks;
	dyn_diff_length = numblocks * Apoints * Bpoints;
}

/**
* We use Matt's trick here of getting a list of indexes here, even though we have to use a for loop.  Pre-calculating the positions IS faster,
* and Matt's trick precomputes as much as possible.
*/
void argo::initCombos()	{
	// Make lists of all vaues of A1,A2,A3, B1,B2,B3 to try
	number_of_A_combos = ((int)(a1_Max / a1_step) * 2 + 1) *
		((int)(a2_Max / a2_step) * 2 + 1) *
		((int)(a3_Max / a3_step) * 2 + 1);
	combos_store.A_combos = new param_combo[number_of_A_combos];
	long int A_combos_counter = 0;
	for (double a3 = -a3_Max; a3 <= a3_Max + a3_step / 2; a3 += a3_step) {
		double p2_term1 = -1.5 * rbase_height * a3;
		double p1_term1 = 0.5 * rbase_height * rbase_height * a3;
		for (double a2 = -a2_Max; a2 <= a2_Max + a2_step / 2; a2 += a2_step) {
			double p1_term2 = -a2 * rbase_height;
			for (double a1 = -a1_Max; a1 <= a1_Max + a1_step / 2; a1 += a1_step) {
				combos_store.A_combos[A_combos_counter].P1 = a1 + p1_term2 + p1_term1;
				combos_store.A_combos[A_combos_counter].P2 = a2 + p2_term1;
				combos_store.A_combos[A_combos_counter].P3 = a3;
				++A_combos_counter;
			}
		}
	}

	number_of_B_combos = ((int)(b1_Max / b1_step) * 2 + 1) *
		((int)(b2_Max / b2_step) * 2 + 1) *
		((int)(b3_Max / b3_step) * 2 + 1);
	combos_store.B_combos = new param_combo[number_of_B_combos];
	long int B_combos_counter = 0;
	for (double b3 = -b3_Max; b3 <= b3_Max + b3_step / 2; b3 += b3_step) {
		double p2_term1 = -1.5 * rbase_height * b3;
		double p1_term1 = 0.5 * rbase_height * rbase_height * b3;
		for (double b2 = -b2_Max; b2 <= b2_Max + b2_step / 2; b2 += b2_step) {
			double p1_term2 = -b2 * rbase_height;
			for (double b1 = -b1_Max; b1 <= b1_Max + b1_step / 2; b1 += b1_step) {
				combos_store.B_combos[B_combos_counter].P1 = b1 + p1_term2 + p1_term1;
				combos_store.B_combos[B_combos_counter].P2 = b2 + p2_term1;
				combos_store.B_combos[B_combos_counter].P3 = b3;
				++B_combos_counter;
			}
		}
	}
	std::sort(combos_store.A_combos, (combos_store.A_combos + (number_of_A_combos)));
	std::sort(combos_store.B_combos, (combos_store.B_combos + (number_of_B_combos)));
}

void argo::initBetaGamma()	{
	// Initialize array for dynamic_diffs once.
	beta_gamma_store.dynamic_diffs = new beta_values[dyn_diff_length];


	int smin = rbase_width / 2 - rsliver_width / 2;
	long long partial_sum1;
	long long partial_sum2;
	long long partial_sum3;
	long long diff;

	// Figure out maximum sliver drift, make difflets larger by that much.  (Apoints+Asliverpoints)
	beta_gamma_store.difflets = new gamma_values[numblocks * Adiffletpoints *
									Bdiffletpoints * num_sliver_blocks ];

	gamma_values difflet;
	for ( long int sb = 0; sb < num_sliver_blocks; ++sb ) {
		long int sliver_block_x = sb * sliver_block_width;
		for ( long int j = 0; j <= 2 * ( totaldriftB + Bsliverdrift ); ++j ) {
			for ( long int i = 0; i <= 2 * ( totaldriftA + Asliverdrift ); ++i ) {
				long int x_initial = smin - totaldriftA - Asliverdrift + i;
				x_initial += sliver_block_x;
				for ( long int b = 0; b < numblocks; b++ ) {
					partial_sum1 = 0;
					partial_sum2 = 0;
					partial_sum3 = 0;
					for ( long int y = MaxB0 - totaldriftB + j + b * blocksize,
							ys = MaxB0 + Bsliverdrift + b * blocksize;
							ys < MaxB0 + ( b + 1 ) * blocksize + Bsliverdrift;
							++y, ++ys ) {
						// Skip all negative elements.
						if ( y < 0 ) {
							continue;
						}
						for ( long int x = x_initial, xs = sliver_block_x;
								xs < sliver_block_x + sliver_block_width;
								++x, ++xs ) {
							diff = -images_store.resamp_sliver.fast_get( xs, ys ) + images_store.resamp_base.fast_get( x, y ) - 1;
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

					beta_gamma_store.difflets[b + i * numblocks + j * numblocks * Adiffletpoints +
					          sb * numblocks * Adiffletpoints * Bdiffletpoints ] = difflet;
				} //end for b
			} //end for i
		} //end for j
	} //end for sb
}

void argo::performGridSearch(bool verbose)	{

	double * yb_arr = new double[numblocks];
	double * yb2_arr = new double[numblocks];
	double * yb3_arr = new double[numblocks];
	double * yb4_arr = new double[numblocks];
	double * yb5_arr = new double[numblocks];
	double * yb6_arr = new double[numblocks];
	int mult1, mult2, mult3;
	mult3 = numblocks * Adiffletpoints * Bdiffletpoints;
	mult2 = numblocks * Adiffletpoints;
	mult1 = numblocks * Apoints;

	for (int i = 0; i < numblocks; ++i) {
		yb_arr[i] = i * blocksize + MaxB0; //+ Bsliverdrift;
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
	for (long int Bindex_i = 0; Bindex_i < number_of_B_combos - 1; Bindex_i = Bindex_f) {
		double B1_i = combos_store.B_combos[Bindex_i].P1;
		double B1_f = B1_i + .01;	//(1.0/32);
		//search ahead for index Bi_f
		for (Bindex_f = Bindex_i + 1; (combos_store.B_combos[Bindex_f].P1 < B1_f) && (Bindex_f < number_of_B_combos); ++Bindex_f);

		for (long int Aindex_i = 0; Aindex_i < number_of_A_combos - 1; Aindex_i = Aindex_f) {
			double A1_i = combos_store.A_combos[Aindex_i].P1;
			double A1_f = A1_i + .01;		//(1.0/32);
			//search ahead for index Ai_f
			for (Aindex_f = Aindex_i + 1; (combos_store.A_combos[Aindex_f].P1 < A1_f) && (Aindex_f < number_of_A_combos); ++Aindex_f);
			
			// Empty out array of beta values.
			if (beta_gamma_store.dynamic_diffs != 0)	{
				delete[] beta_gamma_store.dynamic_diffs;
			}
			beta_gamma_store.dynamic_diffs = new beta_values[dyn_diff_length];

			// Populate array with new beta values.
			long int diffs_address;
			long int difflets_address;
			int Xc = -sliver_block_width;
			int totX;
			int totY;
			int blocktimesApoints = Apoints * numblocks;
			long A0_B0_adj_inc = numblocks;  //A0increase is now assumed to be 1.
			long MaxA0_times_numblocks = MaxA0 * numblocks;
			long B0_adj_min = -MaxB0 * blocktimesApoints;
			long B0_adj_max = MaxB0 * blocktimesApoints;
			//long BO_adj_inc = B0increase*numblocks*Apoints;
			long BO_adj_inc = blocktimesApoints;  //B0increase is now assumed to be 1.
			int finalAddress = -mult3;

			for (int sb = 0; sb < num_sliver_blocks; ++sb) {
				//Now calculate offsets sx and symeh
				finalAddress += mult3;
				Xc += sliver_block_width;
				int Xc2 = Xc * Xc;
				int doubleXc = 2 * Xc;
				int sx = Asliverdrift - (int)(Xc * A1_i);  //which A1 do I use?  first?  last?  average?
				int sy = Bsliverdrift - (int)(Xc * B1_i);
				//printf("%f, %d, %d, %d\n",A1_i,sb,sx,sy);
				int YPart = -mult1;
				for (int dy = 0; dy < Bpoints; ++dy) {
					totY = ( dy + sy ) * mult2;
					YPart += mult1;
					int XPart = -numblocks;
					for (int dx = 0; dx < Apoints; ++dx) {
						totX = (dx + sx) * numblocks;
						XPart += (int)numblocks;
						for (int b = 0; b < numblocks; ++b) {  //partially unroll this loop?
							//diffs[dy,dx,b] += difflets[sb,dy+sy,dx+sx,b]

							//Xc=0;
							//printf("X position: %f \n",Xc);

							//difflets_address =(b + (dx+sx)*numblocks + (dy+sy)*numblocks*Adiffletpoints + sb*numblocks*Adiffletpoints*Bdiffletpoints);
							difflets_address = ( b + totX + totY + finalAddress );
							//diffs_address = b+dx*numblocks+dy*numblocks*Apoints;
							diffs_address = b + XPart + YPart;
							beta_gamma_store.dynamic_diffs[diffs_address].B0 += beta_gamma_store.difflets[difflets_address].G0;
							beta_gamma_store.dynamic_diffs[diffs_address].B1 += beta_gamma_store.difflets[difflets_address].G1;
							beta_gamma_store.dynamic_diffs[diffs_address].B2 += beta_gamma_store.difflets[difflets_address].G2;
							beta_gamma_store.dynamic_diffs[diffs_address].B3 += (beta_gamma_store.difflets[difflets_address].G1 * Xc);	// Adjust for new Beta3 which is sum of negatives
							beta_gamma_store.dynamic_diffs[diffs_address].B4 += beta_gamma_store.difflets[difflets_address].G2 * Xc2;
							beta_gamma_store.dynamic_diffs[diffs_address].B5 += (beta_gamma_store.difflets[difflets_address].G2 * doubleXc);// Adjust for new Beta5 which is sum of negatives
							//}
						}
					}
				}
			}

			tbb::critical_section cs;
			tbb::parallel_for((long int)Bindex_i, Bindex_f, (long int)1, [&](long int Bindex)
			{
				for (long int Aindex = Aindex_i; Aindex < Aindex_f; ++Aindex) {
					double B1 = combos_store.B_combos[Bindex].P1;
					double B2 = combos_store.B_combos[Bindex].P2;
					double B3 = combos_store.B_combos[Bindex].P3;
					double A1 = combos_store.A_combos[Aindex].P1;
					double A2 = combos_store.A_combos[Aindex].P2;
					double A3 = combos_store.A_combos[Aindex].P3;

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

							for (int i = 0; i<numblocks; i++) {
								double beta1, double_beta2, beta3, beta4, beta5;

								long xpos = (long)(A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + totaldriftA);
								long ypos = (long)(B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + totaldriftB);

								int index = i + xpos*numblocks + ypos*blocktimesApoints + A0_B0_adj;

								double yb = yb_arr[i];
								double yb2 = yb2_arr[i];
								double yb3 = yb3_arr[i];
								double yb4 = yb4_arr[i];
								double yb5 = yb5_arr[i];
								double yb6 = yb6_arr[i];

								beta1 = beta_gamma_store.dynamic_diffs[index].B1;
								double_beta2 = beta_gamma_store.dynamic_diffs[index].B2 * 2;
								beta3 = beta_gamma_store.dynamic_diffs[index].B3;
								beta4 = beta_gamma_store.dynamic_diffs[index].B4;
								beta5 = beta_gamma_store.dynamic_diffs[index].B5;

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
							for (int i = 0; i<numblocks; i++){
								long xpos = (long)(A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + totaldriftA);
								long ypos = (long)(B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + totaldriftB);
								//int index = diffs_indices[i] + A0_B0_adj;//*numBytes
								int index = i + xpos*numblocks + ypos*blocktimesApoints + A0_B0_adj;
								//int index = diffs_indices[i] + A0_B0_adj;
								double yb = yb_arr[i];
								double yb2 = yb2_arr[i];
								double yb3 = yb3_arr[i];

								double Beta0 = beta_gamma_store.dynamic_diffs[index].B0;
								double Beta1 = beta_gamma_store.dynamic_diffs[index].B1;
								double Beta2 = beta_gamma_store.dynamic_diffs[index].B2;
								double Beta3 = beta_gamma_store.dynamic_diffs[index].B3;
								double Beta4 = beta_gamma_store.dynamic_diffs[index].B4;
								double Beta5 = beta_gamma_store.dynamic_diffs[index].B5;
								double delz = C0 + C1*yb + C2*yb2 + C3*yb3;

								double sum = Beta0 + Beta1*delz + Beta2*delz*delz
									+ Beta3*C1 + Beta4*C1*C1 + Beta5*C1*delz;

								if (sum < 0)	{
									printf("Negative Sum!\n");
								}

								function += sum;
								if (function > results.bestdiff)	{
									results.iterationsIgnored += numblocks - (i + 1);
									break;
								}
							}
							// Note: this "if" block doesn't seem to slow down the program at all, compared to simpler statements.
							if (function <= results.bestdiff){
								results.bestdiff = function;
								results.bestA0 = (A0_B0_adj - B0_adj) / numblocks;
								results.bestB0 = B0_adj / blocktimesApoints;
								results.bestA1 = A1;
								results.bestA2 = A2;
								results.bestA3 = A3;
								results.bestB1 = B1;
								results.bestB2 = B2;
								results.bestB3 = B3;
								results.bestC0 = C0;
								results.bestC1 = C1;
								results.bestC2 = C2;
								results.bestC3 = C3;

								if (verbose)	{
									logCurrentBest();
								}
							}
							results.count++;
							// Unlock here.
							cs.unlock();
						}
					}
				}
			} );
		}
	}
}

void argo::performSimplexRoutine()	{
	results.bestC1 /= precision;
	results.bestC2 /= precision * precision;
	results.bestC3 /= precision * precision * precision;

	// These are how you should fill the x array. Use this initialization 
	// if you haven't skipped the grid search
	x[0] = results.bestA0 * precision;
	x[1] = results.bestB0 * precision;
	x[2] = results.bestA1;
	x[3] = results.bestB1 + 1;
	x[4] = results.bestA2 / precision;
	x[5] = results.bestB2 / precision;
	x[6] = results.bestA3 / (precision * precision);
	x[7] = results.bestB3 / (precision * precision);
	x[8] = results.bestC0;
	x[9] = results.bestC1;
	x[10] = results.bestC2;
	x[11] = results.bestC3;

	// Setup precision for simplex.  
	// These values are all based on the changes in the parameters that will cause at most a shift of 
	// 0.1 pixel.
	for (int i = 0; i < 12; ++i)	{
		z[i] = x[i];
	}
	
	precisionArr[0] = 0.1;
	precisionArr[1] = 0.1;
	precisionArr[2] = 0.1 / images_store.orig_base.height;
	precisionArr[3] = 0.1 / images_store.orig_base.height;
	precisionArr[4] = 0.1 / pow(images_store.orig_base.height, 2.0);
	precisionArr[5] = 0.1 / pow(images_store.orig_base.height, 2.0);
	precisionArr[6] = 0.1 / pow(images_store.orig_base.height, 3.0);
	precisionArr[7] = 0.1 / pow(images_store.orig_base.height, 3.0);
	precisionArr[8] = 0.1;
	precisionArr[9] = 0.1 / images_store.orig_base.height;
	precisionArr[10] = 0.1 / pow(images_store.orig_base.height, 2.0);
	precisionArr[11] = 0.1 / pow(images_store.orig_base.height, 3.0);

	simplex(&images_store.orig_base, &images_store.orig_sliver, x, z, 12, precisionArr, reflectParam,
		contractParam, growthParam, haltParam,
		maxRefineIterationsParam);
}

void argo::performImageCorrection()	{
	// Perform final warp, and write the output tiff
	double aterms[4] = { z[0], z[2], z[4], z[6] };
	double bterms[4] = { z[1], z[3], z[5], z[7] };
	double cterms[4] = { z[8], z[9], z[10], z[11] };

	image_basic final(images_store.orig_base.width, images_store.orig_base.height, images_store.orig_base.get_color_mode());
	warp_image(&images_store.orig_base, images_store.orig_base.width, images_store.orig_base.height, aterms, bterms, cterms, 1, &final);

	std::string finalstring = rundata1;
	std::basic_string < char > finalmarkerstring("_corrected");
	finalstring.insert(finalstring.rfind("."), finalmarkerstring);
	final.file_write_tiff(finalstring);
}

void argo::correctImages()	{

}

void argo::correctImages(bool verbose)	{

}

void argo::correctImages(int argc, char* argv[])	{
	clock_t begintime = clock();
	bool debugmode = false;
	clock_t time0, time1;

	time0 = clock();
	readInputParams(argc, argv);
	if (debugmode) {
		precision = 1;
		blocksize = 8;
	}
	readImages();
	initCalculatedParams();
	time1 = clock();
	times_store.image_read_time = ((double)time1 - (double)time0) / 1000;

	time0 = clock();
	initCombos();
	time1 = clock();
	times_store.combos_time = ((double)time1 - (double)time0) / 1000;

	time0 = clock();
	initBetaGamma();
	time1 = clock();
	times_store.diffs_time = ((double)time1 - (double)time0) / 1000;

	time0 = clock();
	performGridSearch(false);
	time1 = clock();
	times_store.grid_time = ((double)time1 - (double)time0) / 1000;

	time0 = clock();
	performSimplexRoutine();
	time1 = clock();
	times_store.simplex_time = ((double)time1 - (double)time0) / 1000;

	time0 = clock();
	performImageCorrection();
	time1 = clock();
	times_store.image_write_time = ((double)time1 - (double)time0) / 1000;

	int endtime = clock();
	times_store.total_time = (double)(endtime - begintime) * .001;
}

void argo::correctImages(int argc, char* argv[], bool verbose)	{
	clock_t begintime = clock();
	bool debugmode = false;
	clock_t time0, time1;

	time0 = clock();
	readInputParams(argc, argv);
	if (debugmode) {
		precision = 1;
		blocksize = 8;
	}
	if (verbose)	{
		logInputParams();
	}

	readImages();
	initCalculatedParams();
	time1 = clock();
	times_store.image_read_time = ((double)time1 - (double)time0) / 1000;

	time0 = clock();
	initCombos();
	time1 = clock();
	times_store.combos_time = ((double)time1 - (double)time0) / 1000;

	time0 = clock();
	initBetaGamma();
	time1 = clock();
	times_store.diffs_time = ((double)time1 - (double)time0) / 1000;

	if (verbose)	{
		logCalculatedParams();
		logComboInfo();
		logBetaGammaInfo();
	}

	time0 = clock();
	performGridSearch(verbose);
	time1 = clock();
	times_store.grid_time = ((double)time1 - (double)time0) / 1000;

	if (verbose)	{
		logGridSearchInfo();
	}

	time0 = clock();
	performSimplexRoutine();
	time1 = clock();
	times_store.simplex_time = ((double)time1 - (double)time0) / 1000;

	if (verbose)	{
		logSimplexRoutineInfo();
	}

	time0 = clock();
	performImageCorrection();
	time1 = clock();
	times_store.image_write_time = ((double)time1 - (double)time0) / 1000;

	int endtime = clock();
	times_store.total_time = (double)(endtime - begintime) * .001;

	if (verbose)	{
		logProgramInformation();
	}
}

/**
 * Main method for our program.
 */
int main(int argc, char *argv[]) {

	argo* program = new argo();
	program->correctImages(argc, argv, true);
	delete program;

	getchar();

	//////////////////// OLD CODE ////////////////////
	/*

	// Size of a block in pixels.
	int blocksize = 0;

	// Maximum possible constant shift in pixels. Defaults to 5.
	int MaxA0 = 5;
	int MaxB0 = 5;

	// Maximum linear shift, as fraction of image size
	// (if a1_Max = 0.02, and image height is 500 pixels,
	// this term causes max shift at top of image of 10 pixels)

	double a1_Max = 0.1;
	double b1_Max = 0.1;

	// Fractional shift, relative to a1_max.
	// If a1_max implies max shift of 10 pixels at top, a2_Multiplier of 0.5 causes deviation
	// from linearity of 5 pixels at midpoint.
	//
	double a2_Multiplier = ( 1.0 / 9.0 );
	double b2_Multiplier = ( 1.0 / 9.0 );

	// Fractional shift for cubic term a3, also relative to shift caused by a1_max.
	double a3_Multiplier = .1;  //(0.0/9.0);
	double b3_Multiplier = .1;  //(0.0/9.0);

	//
	// For initial grid search, image is rescaled (downsampled) by this factor.
	// This parameter should be set to the size, in pixels, of the smallest real feature on the image.
	//
	int precision = 1;

	//
	// These parameters control the simplex routine.
	// These parameters were found, using trial and error on not a whole lot of different test cases.
	// To avoid being trapped in some rut, it seemse to be VERY important to have the reflection
	// parameter slightly less than 1.0, and the growth parameter less than 2.0.
	//
	double growthParam = 1.5;
	double contractParam = 0.5;
	double reflectParam = 0.9;
	double haltParam = 1e-11;
	int maxRefineIterationsParam = 5000;



	printf( "Running on version %s\n", VERSION );

	if ( ( argc < 33 ) || ( argc > 33 ) ) {
		//perror("Error: Usage is ImageMaster -i: ImagePath -s: SliverPath -A0: value -A1: value -A2: value -B0: value -B1: value -B2: value -p: precision -b: blocksize -g: growth -c: contract -r: reflection\n");
		//ACTUALLY, there should be at least 31 arguments now, since I added cubic terms
		//Also, usage is -0 A0value, -1 B0value, -2 a1value -3 b1value -4 a2multiplier -5 b2multiplier, etc.
		printf( "Amount of arguments provided : %d\n", argc );
		//getchar();
		//return 0;
	}

	std::string rundata1, rundata2;
	char command_arg;
	for ( int i = 1; i < argc; ++i ) {
		if ( argv[ i ][ 0 ] == '-' ) {
			command_arg = argv[ i ][ 1 ];

			// found an option, so look at next
			// argument to get the value of
			// the option

			switch ( command_arg ) {
				case 'i':
					rundata1 = argv[ i + 1 ];
					for ( unsigned int x = 0; x < rundata1.length(); x++ ) {
						if ( rundata1[ x ] == '|' )
							rundata1[ x ] = ' ';
					}
					break;
				case 's':
					rundata2 = argv[ i + 1 ];
					for ( unsigned int x = 0; x < rundata2.length(); x++ ) {
						if ( rundata2[ x ] == '|' )
							rundata2[ x ] = ' ';
					}
					break;
				case '0':
					//A0percent = atof(argv[i+1]);
					MaxA0 = atoi( argv[ i + 1 ] );
					break;
				case '1':
					//A1percent = atof(argv[i+1]);
					MaxB0 = atoi( argv[ i + 1 ] );
					break;
				case '2':
					//A2percent = atof(argv[i+1]);
					a1_Max = atof( argv[ i + 1 ] );
					break;
				case '3':
					//B0percent = atof(argv[i+1]);
					b1_Max = atof( argv[ i + 1 ] );
					break;
				case '4':
					//B1percent = atof(argv[i+1]);
					a2_Multiplier = atof( argv[ i + 1 ] );
					break;
				case '5':
					//B2percent = atof(argv[i+1]);
					b2_Multiplier = atof( argv[ i + 1 ] );
					break;
				case '6':
					a3_Multiplier = atof( argv[ i + 1 ] );
					break;
				case '7':
					b3_Multiplier = atof( argv[ i + 1 ] );
					break;
				case 'p':
					precision = atoi( argv[ i + 1 ] );
					break;
				case 'b':
					blocksize = atoi( argv[ i + 1 ] );
					break;
				case 'g':
					growthParam = atof( argv[ i + 1 ] );
					break;
				case 'c':
					contractParam = atof( argv[ i + 1 ] );
					break;
				case 'r':
					reflectParam = atof( argv[ i + 1 ] );
					break;
				case 'h':
					haltParam = atof( argv[ i + 1 ] );
					break;
				case 't':
					maxRefineIterationsParam = atoi( argv[ i + 1 ] );
					break;

				default:
					break;
			}
		}
	}


	if ( debugmode ) {
		precision = 1;
		blocksize = 8;
	}

	printf( "MaxA0 %d\n", MaxA0 );
	printf( "MaxB0 %d\n", MaxB0 );
	printf( "Maxa1 %f\n", a1_Max );
	printf( "Maxb1 %f\n", b1_Max );
	printf( "a2 %f\n", a2_Multiplier );
	printf( "b2 %f\n", b2_Multiplier );
	printf( "a3 %f\n", a3_Multiplier );
	printf( "b3 %f\n", b3_Multiplier );
	printf( "precision %d\n", precision );
	printf( "growth %f\n", growthParam );
	printf( "contraction %f\n", contractParam );
	printf( "reflection %f\n", reflectParam );
	printf( "errhalt %f\n", haltParam );
	printf( "maxi %d\n", maxRefineIterationsParam );

	image_basic base_Image, sliver_Image, resamp_base_Image, resamp_sliver_Image;

	// Read files to our image_basic objects
	base_Image.file_read_tiff( rundata1 );
	sliver_Image.file_read_tiff( rundata2 );

	// Resample images based upon precision provided
	resample( &resamp_base_Image, base_Image, precision );
	resample( &resamp_sliver_Image, sliver_Image, precision );

	// HEIGHT will only be used within diffs creation and the main grid search.
	int HEIGHT = resamp_base_Image.height;


	MaxA0 = ( int ) ceil( ( float ) MaxA0 / ( float ) precision );
	MaxB0 = ( int ) ceil( ( float ) MaxB0 / ( float ) precision );

	double bestA0, bestA1, bestA2, bestA3;
	double bestB0, bestB1, bestB2, bestB3;

	// Keep track of best difference during calculations.
	double bestdiff = 1e37;

	// START COMMENTING TO SKIP GRID SEARCH HERE

	//This is the start of code block for dealing with new, smarter way of parametrizing the variables.
	//double a1_Max = (.05 * HEIGHT) / HEIGHT ; //hard coded for now.  This is as a percentage of the image height.
	double a1_step = 1.0 / HEIGHT; //that's 1 pixel
	double a2_Max = a2_Multiplier * ( a1_Max * 4 / HEIGHT );  //one third of max deviation from linear part.  This is totally arbitrary.
	double a2_step = 4.0 / ( HEIGHT * HEIGHT );
	double a3_Max = a3_Multiplier * 20.78 * a1_Max / ( HEIGHT * HEIGHT );
	double a3_step = 20.78 / ( HEIGHT * HEIGHT * HEIGHT ); //changed constant from 10.4 on 1/15/2007 <--- Nathan would like to know where magic 20.78 comes from.


	// kludge below:
	// a3_step *= 2.0;
	// a2_step *= 1.3;
	// a1_step *= 1.0;
	// force "0.0" to be one of the actual values tested for all parameters

	// Make A1, A2, A3 exact ratios (rather than relative) of their respective steps
	a1_Max = ceil( a1_Max / a1_step ) * a1_step;
	a2_Max = ceil( a2_Max / a2_step ) * a2_step;
	a3_Max = ceil( a3_Max / a3_step ) * a3_step;

	//double b1_Max = (.05 * HEIGHT) / HEIGHT ; //hard coded for now.  This is as a percentage of the image height.
	double b1_step = 1.0 / HEIGHT; //that's 1 pixel
	double b2_Max = b2_Multiplier * ( b1_Max * 4 / HEIGHT );
	double b2_step = 4.0 / ( HEIGHT * HEIGHT );
	double b3_Max = b3_Multiplier * 20.78 * b1_Max / ( HEIGHT * HEIGHT );
	double b3_step = 20.78 / ( HEIGHT * HEIGHT * HEIGHT ); //changed constant from 10.4 on 1/15/2007

	// kludge below:
	// b3_step *= 2.0;
	// b2_step *= 1.3;
	// b1_step *= 1.0;
	// force "0.0" to be one of the actual values tested for all parameters

	// Make B1, B2, B3 exact ratios (rather than relative) of their respective steps
	b1_Max = ceil( b1_Max / b1_step ) * b1_step;
	b2_Max = ceil( b2_Max / b2_step ) * b2_step;
	b3_Max = ceil( b3_Max / b3_step ) * b3_step;

	printf("a1_Max = %f, a1_step = %f\n", a1_Max, a1_step);
	printf("a2_Max = %f, a2_step = %f\n", a2_Max, a2_step);
	printf("a3_Max = %f, a3_step = %f\n", a3_Max, a3_step);

	printf("b1_Max = %f, b1_step = %f\n", b1_Max, b1_step);
	printf("b2_Max = %f, b2_step = %f\n", b2_Max, b2_step);
	printf("b3_Max = %f, b3_step = %f\n", b3_Max, b3_step);

	printf( "%f %f %f \n", a1_Max / a1_step, a2_Max / a2_step, a3_Max / a3_step );
	printf( "%f %f %f \n", b1_Max / b1_step, b2_Max / b2_step, b3_Max / b3_step );

	//This set of variables record how many pixels of drift in either axis that will be encountered by the program.
	long totaldriftA = MaxA0 + ( int ) ceil( 1.38 * a1_Max * HEIGHT );
	long totaldriftB = MaxB0 + ( int ) ceil( 1.38 * b1_Max * HEIGHT );
	//note that the factor of 1.38 is hard-coded for a2_multiplier = 1/3, a3_multiplier = 1/9

	printf( "drifts: %d, %d\n", MaxA0, MaxB0 );
	printf( "total drifts: %ld, %ld\n", totaldriftA, totaldriftB );

	//******************
	//* Block Calculations
	//*******************
	// Block size is determined by the maximum drift that may occur at a maximum.
	// This may be changed to become more dynamic in later versions, however, as of
	// now it is separate from the loop, and so worst case must be assumed.
	double maxslope, maxslopeA, maxslopeB = 0;
	maxslopeA = a1_Max + a2_Max * HEIGHT + 0.5 * HEIGHT * HEIGHT * a3_Max;
	maxslopeB = b1_Max + b2_Max * HEIGHT + 0.5 * HEIGHT * HEIGHT * b3_Max;

	// Set maxslope to larger value (inline if)
	maxslope = maxslopeB > maxslopeA ? maxslopeB : maxslopeA;

	if (blocksize <= 0)
	{
		//double maxslope=max((MaxA1+2*MaxA2*(subsize[1]-1)),MaxB1+2*MaxB2*(subsize[1]-1)); //An approximation
		//maxslope = a1_Max + a2_Max*HEIGHT + 0.5*HEIGHT*HEIGHT*a3_Max;
		//This is the correct calculation, though it leads to an inconveniently small blocksize.
		//Better to hardcode at 8 or something?
		//Note: need to adjust this in case b1's are bigger.
		blocksize = (int)(1 / maxslope);
	}
	//int remainderexists=((HEIGHT % blocksize) > 0); //we do this as an int so that it may be included in calculations
	int numblocks = 0;//= HEIGHT/blocksize + (remainderexists);  //this adds one if there is a remainder;
	//int remainderarea = HEIGHT % blocksize;

	int Asliverdrift = ( int ) ( maxslopeA * resamp_sliver_Image.width );
	int Bsliverdrift = ( int ) ( maxslopeB * resamp_sliver_Image.width );
	int Adiffletpoints = ( totaldriftA + Asliverdrift ) * 2 + 1;  //get factor of two correct here.
	int Bdiffletpoints = ( totaldriftB + Bsliverdrift ) * 2 + 1;
	printf( "height removed %ld\n", ( totaldriftB + 2 * Bsliverdrift + MaxB0 ) );
	numblocks = ( HEIGHT - totaldriftB - 2 * Bsliverdrift - MaxB0 ) / blocksize;

	//Hard code the block information, for now.
	//numblocks = subsize[1] / blocksize * 0.9;
	//remainderexists = 0;
	//remainderarea = 0;
	int Apoints = totaldriftA * 2 + 1;
	int Bpoints = totaldriftB * 2 + 1;
	int num_sliver_blocks = resamp_sliver_Image.width/blocksize; // Was hardcoded to 8
	int sliver_block_width = resamp_sliver_Image.width / num_sliver_blocks;
	time0 = clock();

	//************************
	//*Difflets Array Creation
	//************************

	//int smin=imagesize[0]/2-sliver_Image.width/2;
	int smin = resamp_base_Image.width / 2 - resamp_sliver_Image.width / 2;
	long long partial_sum1;
	long long partial_sum2;
	long long partial_sum3;
	long long diff;

	//double * infoC = new double[4];// Never used
	double bestC0, bestC1, bestC2, bestC3;
	//double changeInZForBlock;	Value is never used so commenting out for now.
	//int smin=resamp_base_Image.width /2- resamp_sliver_Image.width/2;
	//unsigned long partial_sum;

	// Compute dyn_diff_length here to avoid computing it again inside loop
	unsigned long int dyn_diff_length = numblocks * Apoints * Bpoints;
	beta_values * dynamic_diffs = new beta_values[ dyn_diff_length ];

	//figure out maximum sliver drift, make difflets larger by that much.  (Apoints+Asliverpoints)
	gamma_values * difflets = new gamma_values[ numblocks * Adiffletpoints * Bdiffletpoints * num_sliver_blocks ];
	gamma_values difflet;
	printf(
			"blocklet width: %d blockHeight (blocksize): %d number of blocklets (num_sliver_blocks): %d number of blocks (numblocks): %d size of beta array: %ld\n",
			resamp_sliver_Image.width / num_sliver_blocks, blocksize, num_sliver_blocks, numblocks, dyn_diff_length );
	
	for ( long int sb = 0; sb < num_sliver_blocks; ++sb ) {
		long int sliver_block_x = sb * sliver_block_width;
		for ( long int j = 0; j <= 2 * ( totaldriftB + Bsliverdrift ); ++j ) {
			for ( long int i = 0; i <= 2 * ( totaldriftA + Asliverdrift ); ++i ) {

				long int x_initial = smin - totaldriftA - Asliverdrift + i;

				x_initial += sliver_block_x;
				for ( long int b = 0; b < numblocks; b++ ) {
					partial_sum1 = 0;
					partial_sum2 = 0;
					partial_sum3 = 0;
					//partial_sum4=0;
					//partial_sum5=0;
					//partial_sum6=0;
					//for(long int y=MaxB0 - totaldriftB - Bsliverdrift + j + b*blocksize, ys=MaxB0+b*blocksize; ys<MaxB0+(b+1)*blocksize; y++, ys++) {
					for ( long int y = MaxB0 - totaldriftB + j + b * blocksize, ys = MaxB0 + Bsliverdrift + b * blocksize;
							ys < MaxB0 + ( b + 1 ) * blocksize + Bsliverdrift; ++y, ++ys ) {
						if ( y < 0 ) {
							continue;
						}					// just skip all negative elements.
											//for(unsigned int x= x_initial, xs=0; xs < resamp_sliver_Image.width; x++, xs++) {
						for ( long int x = x_initial, xs = sliver_block_x; xs < sliver_block_x + sliver_block_width; ++x, ++xs ) {
							//capture 3 points with the differents in z;
							//this gives us three points on our perfect parabola
							//***************************
							///**Changed from main - sliver
							///********june 22 2010********
							//diff = resamp_sliver_Image.fast_get(xs,ys)-resamp_base_Image.fast_get(x,y)-1; //start at -1 delZ
							diff = -resamp_sliver_Image.fast_get( xs, ys ) + resamp_base_Image.fast_get( x, y ) - 1;
							partial_sum1 += diff * diff;
							//printf("X %d Y %d XS %d YS %d ps %lld diff %lld\n",x,y,xs,ys,partial_sum1,diff);
							//getchar();
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

					//diffs[b+i*numblocks+j*numblocks*Apoints] = partial_sum;
					difflets[ b + i * numblocks + j * numblocks * Adiffletpoints + sb * numblocks * Adiffletpoints * Bdiffletpoints ] = difflet;

					//printf("diffs: %d, %d, %d, %lu\n",b,i,j,diffs[b+i*numblocks+j*numblocks*Apoints]);
				} //end for b
			} //end for i
		} //end for j
	} //end for sb


	//unsigned long long sum0, sum1,term0,term1;
	//double function;

		//We use Matt's trick here of getting a list of indexes here, even though we have to use a for loop.  Pre-calculating the positions IS faster,
	//and Matt's trick precomputes as much as possible.

	unsigned long long count = 0;
	unsigned long long iterationsIgnored = 0;
	//float delz;

	//long * diffs_indices1=(long *)&(diffs_indices[1]);  //interestingly, using this never bought me much!
	double * yb_arr = new double[ numblocks ];  //for a 2nd order polynomial, this precalculation of yb and yb^2 made it slower.
	double * yb2_arr = new double[ numblocks ];
	double * yb3_arr = new double[ numblocks ];
	double * yb4_arr = new double[ numblocks ];
	double * yb5_arr = new double[ numblocks ];
	double * yb6_arr = new double[ numblocks ];
	int mult1, mult2, mult3;
	mult3 = numblocks * Adiffletpoints * Bdiffletpoints;
	mult2 = numblocks * Adiffletpoints;
	mult1 = numblocks * Apoints;

	time1 = clock();
	for ( int i = 0; i < numblocks; ++i ) {
		yb_arr[ i ] = i * blocksize + MaxB0; //+ Bsliverdrift;
		yb2_arr[ i ] = yb_arr[ i ] * yb_arr[ i ];
		yb3_arr[ i ] = yb2_arr[ i ] * yb_arr[ i ];
		yb4_arr[ i ] = yb3_arr[ i ] * yb_arr[ i ];
		yb5_arr[ i ] = yb4_arr[ i ] * yb_arr[ i ];
		yb6_arr[ i ] = yb5_arr[ i ] * yb_arr[ i ];
	}

	//*****************************************************
	//* Make lists of all vaues of A1,A2,A3, B1,B2,B3 to try
	//*****************************************************

	long int number_of_B_combos = ((int)(b1_Max / b1_step) * 2 + 1)*((int)(b2_Max / b2_step) * 2 + 1)*((int)(b3_Max / b3_step) * 2 + 1);
	param_combo * B_combos = new param_combo[ number_of_B_combos ];
	long int B_combos_counter = 0;

	for ( double b3 = -b3_Max; b3 <= b3_Max + b3_step / 2; b3 += b3_step ) {
		double p2_term1 = - 1.5 * HEIGHT * b3;
		double p1_term1 = 0.5 * HEIGHT * HEIGHT * b3;
		for ( double b2 = -b2_Max; b2 <= b2_Max + b2_step / 2; b2 += b2_step ) {
			double p1_term2 = - b2 * HEIGHT;
			for ( double b1 = -b1_Max; b1 <= b1_Max + b1_step / 2; b1 += b1_step ) {
				B_combos[ B_combos_counter ].P1 = b1 + p1_term2 + p1_term1;
				B_combos[ B_combos_counter ].P2 = b2 + p2_term1;
				B_combos[ B_combos_counter ].P3 = b3;
				++B_combos_counter;
			}
		}
	}

	//long int number_of_B_combos = ((int)(b1_Max / b1_step) * 2 + 1)*((int)(b2_Max / b2_step) * 2 + 1)*((int)(b3_Max / b3_step) * 2 + 1);
	//param_combo * B_combos = new param_combo[number_of_B_combos];
	//long int B_combos_counter = 0;
	//for (double b1 = -b1_Max; b1 <= b1_Max + b1_step / 2; b1 += b1_step){
	//	for (double b2 = -b2_Max; b2 <= b2_Max + b2_step / 2; b2 += b2_step){
	//		for (double b3 = -b3_Max; b3 <= b3_Max + b3_step / 2; b3 += b3_step){
	//			B_combos[B_combos_counter].P1 = b1 - b2*HEIGHT + 0.5*HEIGHT*HEIGHT*b3;
	//			B_combos[B_combos_counter].P2 = b2 - 1.5 * HEIGHT * b3;
	//			B_combos[B_combos_counter].P3 = b3;
	//			++B_combos_counter;
	//		}
	//	}
	//}

	long int number_of_A_combos = ((int)(a1_Max / a1_step) * 2 + 1)*((int)(a2_Max / a2_step) * 2 + 1)*((int)(a3_Max / a3_step) * 2 + 1);
	param_combo * A_combos = new param_combo[ number_of_A_combos ];
	long int A_combos_counter = 0;
	for ( double a3 = -a3_Max; a3 <= a3_Max + a3_step / 2; a3 += a3_step ) {
		double p2_term1 = -1.5 * HEIGHT * a3;
		double p1_term1 = 0.5 * HEIGHT * HEIGHT * a3;
		for ( double a2 = -a2_Max; a2 <= a2_Max + a2_step / 2; a2 += a2_step ) {
			double p1_term2 = - a2 * HEIGHT;
			for ( double a1 = -a1_Max; a1 <= a1_Max + a1_step / 2; a1 += a1_step ) {
				A_combos[ A_combos_counter ].P1 = a1 + p1_term2 + p1_term1;
				A_combos[ A_combos_counter ].P2 = a2 + p2_term1;
				A_combos[ A_combos_counter ].P3 = a3;
				++A_combos_counter;
			}
		}
	}

	std::sort( B_combos, ( B_combos + ( number_of_B_combos ) ) );
	std::sort( A_combos, ( A_combos + ( number_of_A_combos ) ) );

	time2 = clock();

	FILE *data = fopen("data.txt", "w");

	//TO DO for tomorrow: think about how to do initial and final B1 sensibly and symmetrically
	//Incorporate actual sliver_width into calculation, no hard-wiring to 32.
	//Refill DIFFS
	long int Aindex_f = 0, Bindex_f = 0;

	for ( long int Bindex_i = 0; Bindex_i < number_of_B_combos - 1; Bindex_i = Bindex_f ) {
		double B1_i = B_combos[ Bindex_i ].P1;
		double B1_f = B1_i + .01;	//(1.0/32);
		//search ahead for index Bi_f
		for ( Bindex_f = Bindex_i + 1; ( B_combos[ Bindex_f ].P1 < B1_f ) && ( Bindex_f < number_of_B_combos ); ++Bindex_f );

		for ( long int Aindex_i = 0; Aindex_i < number_of_A_combos - 1; Aindex_i = Aindex_f ) {
			double A1_i = A_combos[ Aindex_i ].P1;
			double A1_f = A1_i + .01;		//(1.0/32);
			//search ahead for index Ai_f
			for ( Aindex_f = Aindex_i + 1; ( A_combos[ Aindex_f ].P1 < A1_f ) && ( Aindex_f < number_of_A_combos ); ++Aindex_f );

			// CALCULATE NEW DIFFS HERE
			// Empty out beta array
			memset( dynamic_diffs, 0, dyn_diff_length * sizeof(beta_values) );

			long int diffs_address;
			long int difflets_address;
			int Xc = -sliver_block_width;
			int totX;
			int totY;
			int blocktimesApoints = Apoints * numblocks;
			long A0_B0_adj_inc = numblocks;  //A0increase is now assumed to be 1.
			long MaxA0_times_numblocks = MaxA0 * numblocks;
			long B0_adj_min = -MaxB0 * blocktimesApoints;
			long B0_adj_max = MaxB0 * blocktimesApoints;
			//long BO_adj_inc = B0increase*numblocks*Apoints;
			long BO_adj_inc = blocktimesApoints;  //B0increase is now assumed to be 1.
			int finalAddress = -mult3;

			for ( int sb = 0; sb < num_sliver_blocks; ++sb ) {
				//Now calculate offsets sx and symeh
				finalAddress += mult3;
				Xc += sliver_block_width;
				int Xc2 = Xc * Xc;
				int doubleXc = 2 * Xc;
				int sx = Asliverdrift - ( int ) ( Xc * A1_i );  //which A1 do I use?  first?  last?  average?
				int sy = Bsliverdrift - ( int ) ( Xc * B1_i );
				//printf("%f, %d, %d, %d\n",A1_i,sb,sx,sy);
				int YPart = -mult1;
				for ( int dy = 0; dy < Bpoints; ++dy ) {
					totY = ( dy + sy ) * mult2;
					YPart += mult1;
					int XPart = -numblocks;
					for ( int dx = 0; dx < Apoints; ++dx ) {
						totX = ( dx + sx ) * numblocks;
						XPart += ( int ) numblocks;
						for ( int b = 0; b < numblocks; ++b ) {  //partially unroll this loop?
							//diffs[dy,dx,b] += difflets[sb,dy+sy,dx+sx,b]

							//Xc=0;
							//printf("X position: %f \n",Xc);

							//difflets_address =(b + (dx+sx)*numblocks + (dy+sy)*numblocks*Adiffletpoints + sb*numblocks*Adiffletpoints*Bdiffletpoints);
							difflets_address = ( b + totX + totY + finalAddress );
							//diffs_address = b+dx*numblocks+dy*numblocks*Apoints;
							diffs_address = b + XPart + YPart;
							dynamic_diffs[ diffs_address ].B0 += difflets[ difflets_address ].G0;
							dynamic_diffs[ diffs_address ].B1 += difflets[ difflets_address ].G1;
							dynamic_diffs[ diffs_address ].B2 += difflets[ difflets_address ].G2;
							dynamic_diffs[ diffs_address ].B3 += ( difflets[ difflets_address ].G1 * Xc );	// Adjust for new Beta3 which is sum of negatives
							dynamic_diffs[ diffs_address ].B4 += difflets[ difflets_address ].G2 * Xc2;
							dynamic_diffs[ diffs_address ].B5 += ( difflets[ difflets_address ].G2 * doubleXc );// Adjust for new Beta5 which is sum of negatives
							//}
						}
					}
				}
			}

			//printf("Betas for block 5\n B0: %g B1: %f B2 %f B3 %f B4 %f B5 %f",
			//dynamic_diffs[5].B0, dynamic_diffs[5].B1,dynamic_diffs[5].B2
			//				,dynamic_diffs[5].B3,dynamic_diffs[5].B4,dynamic_diffs[5].B5);
			//getchar();

			//#define numBytes 8 //change to 4 if you make 
			//difflets, base_diffs_adress, diffs indices, and dynamic diffs 32 bit instead of 64

			critical_section cs;

			

			parallel_for((long int)Bindex_i, Bindex_f, (long int)1, [&](long int Bindex)
			// for ( long int Bindex = Bindex_i; Bindex < Bindex_f; ++Bindex )
			{
				//long * diffs_indices = new  long[numblocks]; Never used
				for (long int Aindex = Aindex_i; Aindex < Aindex_f; ++Aindex) {
					//					printf("combo indices: %d %d \n",Bindex,Aindex);
					double B1 = B_combos[Bindex].P1;
					double B2 = B_combos[Bindex].P2;
					double B3 = B_combos[Bindex].P3;
					double A1 = A_combos[Aindex].P1;
					double A2 = A_combos[Aindex].P2;
					double A3 = A_combos[Aindex].P3;
					//for (long int i=0;i<numblocks;i++){
					//long yb=i*blocksize + MaxB0drift;// + (blocksize/2);
					//long xpos = (long) (     A1*yb + A2*yb*yb + A3*yb*yb*yb + totaldriftA);
					//long ypos = (long) (     B1*yb + B2*yb*yb + B3*yb*yb*yb + totaldriftB);
					//int xpos = (int) (      totaldriftA);
					//int ypos = (int) (   totaldriftB);

					//long xpos = (long) (     A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + totaldriftA);
					//long ypos = (long) (     B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + totaldriftB);

					//long yb=i*blocksize + MaxB0drift;// + (blocksize/2);
					//int xpos = (int) (     A1*yb_arr[i] + A2*yb2_arr[i] + totaldriftA);
					//int ypos = (int) ( (B1-1)*yb_arr[i] + B2*yb2_arr[i] + totaldriftB);
					//diffs_indices[i]=xpos*numblocks+ypos*numblocks*Apoints+i;

					//diffs_indices[i]=(int) ( (      A1*yb + A2*yb*yb + totaldriftA  )
					//						+( (  (B1-1)*yb + B2*yb*yb + totaldriftB)  * Apoints) * numblocks);

					//}
					//long A0_B0_adj_inc = A0increase*numblocks;

					//for (long int i=0;i<numblocks;i++){
					//
					//	long xpos = (long) (     A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + totaldriftA);
					//	long ypos = (long) (     B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + totaldriftB);
					//
					//
					//	diffs_indices[i]=(i + xpos*numblocks + ypos*blocktimesApoints);
					//
					// }


					for (long B0_adj = B0_adj_min; B0_adj <= B0_adj_max; B0_adj += BO_adj_inc) {
						long A0_B0_adj_max = MaxA0_times_numblocks + B0_adj;
						//for (double A0=-MaxA0;A0<=MaxA0;A0+=A0increase){

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

							for (int i = 0; i<numblocks; i++) {

								double beta1, double_beta2, beta3, beta4, beta5;

								//term0 = * (unsigned long long *)((char*)base_diffs_address+(* ((long long *) (((char*)diffs_indices)+i))));
								//term1 = * (unsigned long long *)((char*)base_diffs_address+(* ((long long *) (((char*)diffs_indices)+i+numBytes))));
								//sum0 += term0;
								//sum1 += term1;
								//j + A0_B0_adj;
								//if(i<30 && i %10 ==0){

								//int index = diffs_indices[i] + A0_B0_adj;//*numBytes

								long xpos = (long)(A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + totaldriftA);
								long ypos = (long)(B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + totaldriftB);

								int index = i + xpos*numblocks + ypos*blocktimesApoints + A0_B0_adj;

								double yb = yb_arr[i];
								double yb2 = yb2_arr[i];
								double yb3 = yb3_arr[i];
								double yb4 = yb4_arr[i];
								double yb5 = yb5_arr[i];
								double yb6 = yb6_arr[i];

								beta1 = dynamic_diffs[index].B1;
								double_beta2 = dynamic_diffs[index].B2 * 2;
								beta3 = dynamic_diffs[index].B3;
								beta4 = dynamic_diffs[index].B4;
								beta5 = dynamic_diffs[index].B5;

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
								//	}

							}

							C = A.inverse()*B;
							double C0 = C.data()[0];
							double C1 = C.data()[1];
							double C2 = C.data()[2];
							double C3 = C.data()[3];
							

							
							// Critical section lock to prevent data overwriting.
							cs.lock();
							double function = 0;
							for (int i = 0; i<numblocks; i++){
								long xpos = (long)(A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + totaldriftA);
								long ypos = (long)(B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + totaldriftB);
								//int index = diffs_indices[i] + A0_B0_adj;//*numBytes
								int index = i + xpos*numblocks + ypos*blocktimesApoints + A0_B0_adj;
								//int index = diffs_indices[i] + A0_B0_adj;
								double yb = yb_arr[i];
								double yb2 = yb2_arr[i];
								double yb3 = yb3_arr[i];

								double Beta0 = dynamic_diffs[index].B0;
								double Beta1 = dynamic_diffs[index].B1;
								double Beta2 = dynamic_diffs[index].B2;
								double Beta3 = dynamic_diffs[index].B3;
								double Beta4 = dynamic_diffs[index].B4;
								double Beta5 = dynamic_diffs[index].B5;
								double delz = C0 + C1*yb + C2*yb2 + C3*yb3;

								double sum = Beta0 + Beta1*delz + Beta2*delz*delz
									+ Beta3*C1 + Beta4*C1*C1 + Beta5*C1*delz;

								if (sum < 0)	{
									std::cout << "NEGATIVE SUM!!!" << std::endl;
								}

								function += sum;
								if (function > bestdiff)	{
									iterationsIgnored += numblocks - (i + 1);
									break;
								}
							}
							// Note: this "if" block doesn't seem to slow down the program at all, compared to simpler statements.
							if (function <= bestdiff){
								bestdiff = function;
								bestA0 = (A0_B0_adj - B0_adj) / numblocks;
								bestB0 = B0_adj / blocktimesApoints;
								bestA1 = A1;
								bestA2 = A2;
								bestA3 = A3;
								bestB1 = B1;
								bestB2 = B2;
								bestB3 = B3;
								bestC0 = C0;
								bestC1 = C1;
								bestC2 = C2;
								bestC3 = C3;

								printf("\n");
								printf("BEST: %e   COUNT: %I64u\n", bestdiff, count);
								printf("           %f, %f, %e, %e\n", bestA0, A1, A2, A3);
								printf("           %f, %f, %e, %e\n", bestB0, B1, B2, B3);
								printf("           %f, %f, %e, %e\n", bestC0, bestC1, bestC2, bestC3);

								fprintf(data, "\n");
								fprintf(data, "BEST: %e   COUNT: %I64u\n", bestdiff, count);
								fprintf(data, "           %f, %f, %e, %e\n", bestA0, A1, A2, A3);
								fprintf(data, "           %f, %f, %e, %e\n", bestB0, B1, B2, B3);
								fprintf(data, "           %f, %f, %e, %e\n", bestC0, bestC1, bestC2, bestC3);
							}
							count++;
							// Unlock here.
							cs.unlock();
						}
					}
				}
			} );
		}
	}

	fclose(data);

	time3 = clock();

		printf( "Time 0: \t %I64u \n", ( long long ) time0 );
	printf( "Time 1: \t %I64u \n", ( long long ) time1 );
	printf( "Time 2: \t %I64u \n", ( long long ) time2 );
	printf( "Time 3: \t %I64u \n", ( long long ) time3 );
	printf( "Time for creating diffs array: \t %f \n", ( time1 - time0 ) * 0.001 );
	printf( "Time for creating and sorting A and B param arrays: \t %f \n", ( time2 - time1 ) * 0.001 );
	double grid_time = ( ( double ) time3 - ( double ) time2 ) / 1000;
	printf( "Time for grid search: \t %lf \n", grid_time );
	printf( "count: %I64u\n", count );
	printf( "Grid search calculations took %lf clockcycles per set of parameters.\n", ( grid_time * 3460000000.0 ) / ( double ) count ); //!~! clock line

	bestC1 /= precision;
	bestC2 /= precision * precision;
	bestC3 /= precision * precision * precision;
	//bestC0=bestC1=bestC2=bestC3=0;

	// END COMMENTING TO SKIP GRID SEARCH HERE

	printf( "Done grid search\n" );

	int n = 12;         // number of parameters to be optimized 

	//now calculate the starting parameters for the simplex routine
	double x[ 12 ];


	// These are how you should fill the x array. Use this initialization 
	// if you haven't skipped the grid search
	x[ 0 ] = bestA0 * precision;
	x[ 1 ] = bestB0 * precision;
	x[ 2 ] = bestA1;
	x[ 3 ] = bestB1 + 1;
	x[ 4 ] = bestA2 / precision;
	x[ 5 ] = bestB2 / precision;
	x[ 6 ] = bestA3 / ( precision * precision );
	x[ 7 ] = bestB3 / ( precision * precision );
	x[ 8 ] = bestC0;
	x[ 9 ] = bestC1;
	x[ 10 ] = bestC2;
	x[ 11 ] = bestC3;

	// x[0] = bestA0 ;          
	// x[1] = bestB0;
	// x[2] = bestA1;          
	// x[3] = bestB1+1;
	// x[4] = bestA2;
	// x[5] = bestB2;
	// x[6] = bestA3;
	// x[7] = bestB3;
	// x[8] = bestC0;
	// x[9] = bestC1;
	// x[10] = bestC2;
	// x[11] = bestC3;

	//intentionally screw up the best guess, as a torture test for the simplex refinement routine
	//x[2] *=1.2
	//x[4] *=1.2
	//x[1] +=1

	// Setup precision for simplex.  
	// These values are all based on the changes in the parameters that will cause at most a shift of 
	// 0.1 pixel.
	double z[ 12 ];
	for ( int i = 0; i < n; ++i )
		z[ i ] = x[ i ];
	
	// printf("Coefficients of drift were\n");
	// printf("Totaldiff: %d\n",bestdiff);
	// printf("A0   %f\n",z[0]);
	// printf("A1   %f   %f\n",z[2],z[2]*base_Image.height);
	// printf("A2   %e   %f\n",z[4],z[4]*base_Image.height*base_Image.height);
	// printf("A3   %e   %f\n",z[6],z[6]*base_Image.height*base_Image.height*base_Image.height);
	// printf("B0   %f\n",z[1]);
	// printf("B1   %f   %f\n",z[3]-1,(z[3]-1)*base_Image.height);
	// printf("B2   %e   %f\n",z[5],z[5]*base_Image.height*base_Image.height);
	// printf("B3   %e   %f\n",z[7],z[7]*base_Image.height*base_Image.height*base_Image.height);
	// printf("C0: %e  %f\nC1: %e  %f\nC2: %e %f\nC3: %e % f\n", z[8],z[8],z[9],z[9]*base_Image.height,z[10],z[10]*base_Image.height*base_Image.height,z[11],z[11]*base_Image.height*base_Image.height*base_Image.height);
	 
	//printf("time: %f\n",time2-time3);
	printf( "Beginning simplex\n" );
	//double ycorrect; Never gets used so commenting out for now.
	double precisionArr[ 12 ];
	precisionArr[ 0 ] = 0.1;
	precisionArr[ 1 ] = 0.1;
	precisionArr[ 2 ] = 0.1 / base_Image.height;
	precisionArr[ 3 ] = 0.1 / base_Image.height;
	precisionArr[ 4 ] = 0.1 / pow( base_Image.height, 2.0 );
	precisionArr[ 5 ] = 0.1 / pow( base_Image.height, 2.0 );
	precisionArr[ 6 ] = 0.1 / pow( base_Image.height, 3.0 );
	precisionArr[ 7 ] = 0.1 / pow( base_Image.height, 3.0 );
	precisionArr[ 8 ] = 0.1;	///base_Image.height;
	precisionArr[ 9 ] = 0.1 / base_Image.height;
	precisionArr[ 10 ] = 0.1 / pow( base_Image.height, 2.0 );
	precisionArr[ 11 ] = 0.1 / pow( base_Image.height, 3.0 );

	simplex( &base_Image, &sliver_Image, x, z, n, precisionArr, reflectParam,
			contractParam, growthParam, haltParam,
			maxRefineIterationsParam );

	printf( "Coefficients of drift were\n" );
	printf( "Totaldiff: %f\n", bestdiff );
	printf( "A0   %f\n", z[ 0 ] );
	printf( "A1   %f   %f\n", z[ 2 ], z[ 2 ] * base_Image.height );
	printf( "A2   %e   %f\n", z[ 4 ], z[ 4 ] * base_Image.height * base_Image.height );
	printf( "A3   %e   %f\n", z[ 6 ], z[ 6 ] * base_Image.height * base_Image.height * base_Image.height );
	printf( "B0   %f\n", z[ 1 ] );
	printf( "B1   %f   %f\n", z[ 3 ] - 1, ( z[ 3 ] - 1 ) * base_Image.height );
	printf( "B2   %e   %f\n", z[ 5 ], z[ 5 ] * base_Image.height * base_Image.height );
	printf( "B3   %e   %f\n", z[ 7 ], z[ 7 ] * base_Image.height * base_Image.height * base_Image.height );
	printf( "C0: %e  %f\nC1: %e  %f\nC2: %e %f\nC3: %e % f\n", z[ 8 ], z[ 8 ], z[ 9 ], z[ 9 ] * base_Image.height, z[ 10 ],
			z[ 10 ] * base_Image.height * base_Image.height, z[ 11 ], z[ 11 ] * base_Image.height * base_Image.height * base_Image.height );

	time4 = clock();
	printf( "Time for simplex refinement: \t %f seconds\n", ( time4 - time3 ) * 0.001 );

	//now do a final warp, and write the output tiff
	double aterms[ 4 ] = { z[ 0 ], z[ 2 ], z[ 4 ], z[ 6 ] };  //-1 is temporary, needed to fix an offset error.
	double bterms[ 4 ] = { z[ 1 ], z[ 3 ], z[ 5 ], z[ 7 ] };
	double cterms[ 4 ] = { z[ 8 ], z[ 9 ], z[ 10 ], z[ 11 ] };
	image_basic final( base_Image.width, base_Image.height, base_Image.get_color_mode() );
	warp_image( &base_Image, base_Image.width, base_Image.height, aterms, bterms, cterms, 1, &final );
	int endtime = clock();

	// std::string finalstring = rundata1;

	std::string finalstring = rundata1;
	std::basic_string < char > finalmarkerstring( "_corrected" );
	finalstring.insert( finalstring.rfind( "." ), finalmarkerstring );

	final.file_write_tiff( finalstring );
	double totalelapsed = ( double ) ( endtime - begintime ) * .001;

	HEIGHT = final.height;
	printf( "Trying to output log data\n\n" );

	FILE *simplexfile = fopen("simplex_parameters.txt", "w");
	FILE *gridfile = fopen("grid_search_parameters.txt", "w");
	FILE *timefile = fopen( "runtimes.txt", "a" );

	if (simplexfile == NULL || timefile == NULL) {
		printf( "file Not opened" );
	}
	else {
		fprintf(simplexfile, "A Params\n");
		fprintf(simplexfile, "%e %e %e %e\n", z[0], z[2], z[4], z[6]);
		fprintf(simplexfile, "B Params\n");
		fprintf(simplexfile, "%e %e %e %e\n", z[1], z[3] - 1.0, z[5], z[7]);
		fprintf(simplexfile, "C Params\n");
		fprintf(simplexfile, "%e %e %e %e\n", z[8], z[9], z[10], z[11]);

		fprintf(gridfile, "A Params\n");
		fprintf(gridfile, "%e %e %e %e\n", x[0], x[2], x[4], x[6]);
		fprintf(gridfile, "B Params\n");
		fprintf(gridfile, "%e %e %e %e\n", x[1], x[3] - 1.0, x[5], x[7]);
		fprintf(gridfile, "C Params\n");
		fprintf(gridfile, "%e %e %e %e\n", x[8], x[9], x[10], x[11]);

		fprintf( timefile, "Total Program Time =\t%f\n", totalelapsed );
		fprintf( timefile, "Total number of iterations ignored = %I64u\n", iterationsIgnored );
		printf( "%e %e %e %e %e %e %e %e\n", z[ 0 ], z[ 2 ], z[ 4 ], z[ 6 ], z[ 1 ], z[ 3 ] - 1.0, z[ 5 ], z[ 7 ] );
		fclose(simplexfile);
		fclose( timefile );
		fclose(gridfile);
		printf( "------End of Program------\n" );
	}
	//getchar();  //puts a getchar in only if not run from commandline.
	
	*/
	////////////////// OLD CODE END //////////////////
}
