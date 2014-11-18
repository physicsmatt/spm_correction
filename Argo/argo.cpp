#define _CRT_SECURE_NO_WARNINGS

#include <algorithm>
#include <time.h>
#include <iostream>
#include <tbb\parallel_for.h>
#include <tbb\critical_section.h>
#include "Eigen\Eigen"
#include "simplex.h"
#include "argo.h"

//#define SANITYCHECKS


/**
 *	Returns < for param_combo based upon comparisons of P1, P2, and P3 with priorities in that order.
 *
 *	@param a
 *	@param b
 *	@return The smaller param_combo.
 */
bool operator< ( const param_combo& a, const param_combo& b ) {
	if ( a.P1 != b.P1 )
		return ( a.P1 < b.P1 );
	else if ( a.P2 != b.P2 )
		return ( a.P2 < b.P2 );
	else
		return ( a.P3 < b.P3 );
}

/**
 *	Returns > for param_combo based upon comparisons of P1, P2, and P3 with priorities in that order.
 *
 *	@param a
 *	@param b
 *	@return The larger param_combo
 */
bool operator> ( const param_combo& a, const param_combo& b ) {
	if ( a.P1 != b.P1 )
		return ( a.P1 > b.P1 );
	else if ( a.P2 != b.P2 )
		return ( a.P2 > b.P2 );
	else
		return ( a.P3 > b.P3 );
}

/**
 *	Returns <= for param_combo based upon comparisons of P1, P2, and P3 with priorities in that order.
 *
 *	@param a
 *	@param b
 *	@return The smaller param_combo if different.
 */
bool operator<= ( const param_combo& a, const param_combo& b ) {
	if ( a.P1 != b.P1 )
		return ( a.P1 < b.P1 );
	else if ( a.P2 != b.P2 )
		return ( a.P2 < b.P2 );
	else
		return ( a.P3 <= b.P3 );
}

/**
 *	Returns >= for param_combo based upon comparisons of P1, P2, and P3 with priorities in that order.
 *
 *	@param a
 *	@param b
 *	@return The larger param_combo if different.
 */
bool operator>= ( const param_combo& a, const param_combo& b ) {
	if ( a.P1 != b.P1 )
		return ( a.P1 > b.P1 );
	else if ( a.P2 != b.P2 )
		return ( a.P2 > b.P2 );
	else
		return ( a.P3 >= b.P3 );
}

/**
 *	Constructor. Initializes default program values.
 */
argo::argo () {
	blocksize = 0;
	A0_max = 5, B0_max = 5;
	a1_max = 0.1, b1_max = 0.1;
	a2_mult = ( 1.0 / 9.0 ), b2_mult = ( 1.0 / 9.0 );
	a3_mult = .1, b3_mult = .1;
	precision = 1;

	/**
	 *	These parameters were found, using trial and error on not a whole lot of different test cases.
	 *	To avoid being trapped in some rut, it seems to be VERY important to have the reflection
	 *	parameter slightly less than 1.0, and the growth parameter less than 2.0.
	 */
	simplex_growth = 1.5, simplex_contract = 0.5, simplex_reflect = 0.9;
	simplex_iterations = 5000;

	results.bestdiff = INFINITY;
	results.count = 0;
	results.iterations_ignored = 0;

	images_store.flipped = 0;
}

/**
 *	Destructor. Deletes any allocations that had taken place.
 */
argo::~argo () {
	if ( beta_gamma_store.difflets ) {
		delete[] beta_gamma_store.difflets;
	}
	if ( beta_gamma_store.dynamic_diffs ) {
		delete[] beta_gamma_store.dynamic_diffs;
	}
	if ( combos.A_combos ) {
		delete[] combos.A_combos;
	}
	if ( combos.B_combos ) {
		delete[] combos.B_combos;
	}
	if ( images_store.orig_base ) {
		delete images_store.orig_base;
	}
	if ( images_store.resamp_base ) {
		delete images_store.resamp_base;
	}
	if ( images_store.orig_sliver ) {
		delete images_store.orig_sliver;
	}
	if ( images_store.resamp_sliver ) {
		delete images_store.resamp_sliver;
	}
}

/**
 *	Log to screen the input parameters.
 */
void argo::logInputParams () {
	printf( "******** Input Parameters ********\n" );
	printf( "Input Max A0 : \t %d\n", A0_max );
	printf( "Input Max B0 : \t %d\n", B0_max );
	printf( "Input Max a1 : \t %f\n", a1_max );
	printf( "Input Max b1 : \t %f\n", b1_max );
	printf( "Input a2 Multiplier : \t %f\n", a2_mult );
	printf( "Input b2 Multiplier : \t %f\n", b2_mult );
	printf( "Input a3 Multiplier : \t %f\n", a3_mult );
	printf( "Input b3 Multiplier : \t %f\n", b3_mult );

	printf( "******** Simplex Routine Information ********\n" );
	printf( "Input Precision : \t %d\n", precision );
	printf( "Input Growth : \t %f\n", simplex_growth );
	printf( "Input Contraction : \t %f\n", simplex_contract );
	printf( "Input Reflection : \t %f\n", simplex_reflect );
	printf( "Input Maximum Iterations : \t %d\n", simplex_iterations );
	printf( "\n" );
}

/**
 *	Log to screen the calculated parameters.
 */
void argo::logCalculatedParams () {
	printf( "******** Calculated Parameters ********\n" );
	printf( "a1 Max : %f, a1 Step Size : %f, a1 Steps : %f\n", a1_max, a1_step, a1_max / a1_step );
	printf( "a2 Max : %f, a2 Step Size : %f, a2 Steps : %f\n", a2_max, a2_step, a2_max / a2_step );
	printf( "a3 Max : %f, a3 Step Size : %f, a3 Steps : %f\n", a3_max, a3_step, a3_max / a3_step );

	printf( "b1 Max : %f, b1 Step Size : %f, b1 Steps : %f\n", b1_max, b1_step, b1_max / b1_step );
	printf( "b2 Max : %f, b2 Step Size : %f, b2 Steps : %f\n", b2_max, b2_step, b2_max / b2_step );
	printf( "b3 Max : %f, b3 Step Size : %f, b3 Steps : %f\n", b3_max, b3_step, b3_max / b3_step );

	printf( "A0 Max : %d, B0 Max %d\n", A0_max, B0_max );
	printf( "total drifts : \t %ld, %ld\n", A_total_drift, B_total_drift );
	printf( "height removed : \t %ld\n", ( B_total_drift + 2 * B_sliver_drift + B0_max ) );
	printf( "\n" );
}

/**
 *	Log to screen the beta gamma diffs information.
 */
void argo::logBetaGammaInfo () {
	printf( "Number of Blocks : \t %d\n", numblocks );
	printf( "Number of Blocklets : \t %d\n", num_sliver_blocks );
	printf( "Blocklet Pixel Width : \t %d\n", rsliver_width / num_sliver_blocks );
	printf( "Block and Blocklet Pixel Height : \t %d\n", blocksize );
	printf( "Size of Beta Array : \t %d\n", dyn_diffs_size );
	printf( "\n" );
}

/**
 *	Log to screen information about the current best value grid-search has computed.
 */
void argo::logCurrentBest () {
	printf( "\n" );
	printf( "BEST: %e   COUNT: %I64u\n", results.bestdiff, results.count );
	printf( "           %e, %e, %e, %e\n", results.bestA0, results.bestA1, results.bestA2, results.bestA3 );
	printf( "           %e, %e, %e, %e\n", results.bestB0, results.bestB1, results.bestB2, results.bestB3 );
	printf( "           %e, %e, %e, %e\n", results.bestC0, results.bestC1, results.bestC2, results.bestC3 );
	printf( "\n" );
}

/**
 *	Log to screen the final information from grid-search.
 */
void argo::logGridSearchInfo () {
	printf( "******** Finished Grid Search ********\n" );
	printf( "Time for Diffs : \t %f \n", times.diffs_time );
	printf( "Time for Combos : \t %f \n", times.combos_time );
	printf( "Time for Grid Search: \t %lf \n", times.grid_time );
	printf( "Total Count : \t %I64u\n", results.count );
	printf( "Iterations Ignored : \t %I64u\n", results.iterations_ignored );
	printf( "Total Difference : \t %f\n", results.bestdiff );
	printf( "\n" );
}

/**
 *	Log to screen the final information for Simplex.
 */
void argo::logSimplexRoutineInfo () {
	printf( "******** Finished Simplex Routine ********\n" );
	printf( "Drift Coefficients :\n" );
	printf( "Best A0 : %f\n", simplex_best[ 0 ] );
	printf( "Best A1 : %f\n", simplex_best[ 2 ] );
	printf( "Best A2 : %e\n", simplex_best[ 4 ] );
	printf( "Best A3 : %e\n", simplex_best[ 6 ] );
	printf( "Best B0 : %f\n", simplex_best[ 1 ] );
	printf( "Best B1 : %f\n", simplex_best[ 3 ] - 1 );
	printf( "Best B2 : %e\n", simplex_best[ 5 ] );
	printf( "Best B3 : %e\n", simplex_best[ 7 ] );
	if ( !simplex_mode ) {
		printf( "Best C0 : %e\nBest C1 : %e\nBest C2 : %e\nBest C3 : %e\n", simplex_best[ 8 ], simplex_best[ 9 ], simplex_best[ 10 ], simplex_best[ 11 ] );
	}
	else {
		printf( "Best Sliver A1 : %e\nBest Sliver B1 : %e\n", simplex_best[ 8 ], simplex_best[ 9 ] );
	}
	printf( "Time for Simplex Routine : \t %f seconds\n", times.simplex_time );
	printf( "\n" );
}

/**
 *	Log to file both grid-search and simplex information.
 */
void argo::logProgramInformation () {
	FILE *paramfile = fopen( "output.txt", "w" );

	if ( paramfile == NULL ) {
		printf( "File Not opened!" );
	}
	else {
		fprintf( paramfile, "Simplex Routine Parameters : \n" );
		fprintf( paramfile, "A0 = %e;\nA1 = %e;\nA2 = %e;\nA3 = %e;\n", simplex_best[ 0 ], simplex_best[ 2 ], simplex_best[ 4 ], simplex_best[ 6 ] );
		fprintf( paramfile, "B0 = %e;\nB1 = %e;\nB2 = %e;\nB3 = %e;\n", simplex_best[ 1 ], simplex_best[ 3 ] - 1.0, simplex_best[ 5 ], simplex_best[ 7 ] );
		if ( !simplex_mode ) {
			fprintf( paramfile, "C0 = %e;\nC1 = %e;\nC2 = %e;\nC3 = %e;\n", simplex_best[ 8 ], simplex_best[ 9 ], simplex_best[ 10 ], simplex_best[ 11 ] );
			fprintf( paramfile, "Sliver A1 = %e;\nSliver B1 = %e;\nSliver C1 = %e'\n", simplex_best[ 12 ], simplex_best[ 13 ] - 1.0, simplex_best[ 14 ] );
		}
		else {
			fprintf( paramfile, "Sliver A1 = %e;\nSliver B1 = %e;\n", simplex_best[ 8 ], simplex_best[ 9 ] - 1.0 );
		}
		
		fprintf( paramfile, "\n" );
		fprintf( paramfile, "Grid Search Parameters : \n" );
		fprintf( paramfile, "A0 = %e;\nA1 = %e;\nA2 = %e;\nA3 = %e;\n", grid_best[ 0 ], grid_best[ 2 ], grid_best[ 4 ], grid_best[ 6 ] );
		fprintf( paramfile, "B0 = %e;\nB1 = %e;\nB2 = %e;\nB3 = %e;\n", grid_best[ 1 ], grid_best[ 3 ] - 1.0, grid_best[ 5 ], grid_best[ 7 ] );
		fprintf( paramfile, "C0 = %e;\nC1 = %e;\nC2 = %e;\nC3 = %e;\n", grid_best[ 8 ], grid_best[ 9 ], grid_best[ 10 ], grid_best[ 11 ] );
		fprintf( paramfile, "\n" );
		fprintf( paramfile, "Time Information : \n" );
		fprintf( paramfile, "Image Read Time : \t%f\n", times.image_read_time );
		fprintf( paramfile, "Combos Time : \t%f\n", times.combos_time );
		fprintf( paramfile, "Diffs Time : \t%f\n", times.diffs_time );
		fprintf( paramfile, "Grid Search Time : \t%f\n", times.grid_time );
		fprintf( paramfile, "Simplex Routine Time : \t%f\n", times.simplex_time );
		fprintf( paramfile, "Image Write Time\t%f\n", times.image_write_time );
		fprintf( paramfile, "Total Program Time : \t%f\n", times.total_time );
		fprintf( paramfile, "\n");
		fprintf( paramfile, "Program Information : \n" );
		fprintf( paramfile, "Total Parameters Generated : %I64u\n", results.count );
		fprintf( paramfile, "Best Difference Calculated by Grid Search : %e\n", results.bestdiff );
		fprintf( paramfile, "Total Iterations Ignored : %I64u\n", results.iterations_ignored );

		fclose( paramfile );
		printf( "********End of Program********\n" );
	}
}

/**
*	Read in image files.
*	
*	@param	verbose	Boolean flag to print debug info.
*/
void argo::readImages ( bool verbose ) {
	// Read image files and resample based upon provided precision.
	if ( verbose ) {
		printf( "Reading in images : '%s' and '%s'.\n", base_name.c_str(), sliver_name.c_str() );
	}
	images_store.orig_base = new FImage( base_name, images_store.flipped );
	images_store.orig_sliver = new FImage( sliver_name, images_store.flipped );

	if ( verbose ) {
		printf( "Images read succesfully!\n" );
		printf( "Image Type : %d.\n", images_store.orig_base->metadata.type );
		printf( "Image Format : %d.\n", images_store.orig_base->metadata.format );
	}



	images_store.resamp_base = new FImage( images_store.orig_base->metadata );
	images_store.resamp_sliver = new FImage( images_store.orig_base->metadata );

	images_store.orig_base->resample( images_store.resamp_base, precision );
	images_store.orig_sliver ->resample( images_store.resamp_sliver, precision );

	// Set dimension parameters.
	rbase_width = images_store.resamp_base->width;
	rbase_height = images_store.resamp_base->height;
	rsliver_width = images_store.resamp_sliver->width;
	rsliver_height = images_store.resamp_sliver->height;
}

/**
 *	Reads and parses a string of values. This is for command-line argo invocation.
 *	Basic argo call : ./argo.exe -i "base.tif" -s "sliver.tif" -f 0 -0 3 -1 3 -2 .1 -3 .1 -4 .04 -5 .04 -6 .04 -7 .04 -p 2 -x 2 -b 0 -g 1.8 -c 0.7 -r 1.3 -t 15000
 *
 *	-i		The base image.
 *	-s		The sliver image.
 *	-f		Boolean flag indicating if the image needs to be flipped horizontally.
 *	-0		The maximum A0 value.
 *	-1		The maximum B0 value.
 *	-2		The maximum a1 (not A1) value.
 *	-3		The maximum b1 (not B1) value.
 *	-4		The multiplier for a2.
 *	-5		The multiplier for b2.
 *	-6		The multiplier for a3.
 *	-7		The multiplier for b3.
 *	-p		The precision to down-scale the image (i.e -p 2 means down-scale the image by 50%)
 *	-b		The block size.
 *	-g		The simplex growth parameter.
 *	-c		The simplex contraction parameter.
 *	-r		The simplex reflection parameter.
 *	-t		The maximum amount of iterations the simplex routine should run regardless of convergence.
 *	-z		The z-correction mode, 0 for slow-z. Everything else results in fast-z.
 *
 *	@param	argc	The number of arguments.
 *	@param	argv	The arguments.
 */
void argo::readInputParams ( int argc, char *argv[] ) {

	char command_arg;
	for ( int i = 1; i < argc; ++i ) {
		if ( argv[ i ][ 0 ] == '-' ) {
			command_arg = argv[ i ][ 1 ];

			// Look at argument after option to get corresponding value
			switch ( command_arg ) {
				case 'i':
					base_name = argv[ i + 1 ];
					for ( unsigned int x = 0; x < base_name.length(); x++ ) {
						if ( base_name[ x ] == '|' )
							base_name[ x ] = ' ';
					}
					break;
				case 's':
					sliver_name = argv[ i + 1 ];
					for ( unsigned int x = 0; x < sliver_name.length(); x++ ) {
						if ( sliver_name[ x ] == '|' )
							sliver_name[ x ] = ' ';
					}
					break;
				case 'f':
					images_store.flipped = atoi( argv[ i + 1 ] ) != 0 ? true : false;
					break;
				case '0':
					A0_max = atoi( argv[ i + 1 ] );
					break;
				case '1':
					B0_max = atoi( argv[ i + 1 ] );
					break;
				case '2':
					a1_max = atof( argv[ i + 1 ] );
					break;
				case '3':
					b1_max = atof( argv[ i + 1 ] );
					break;
				case '4':
					a2_mult = atof( argv[ i + 1 ] );
					break;
				case '5':
					b2_mult = atof( argv[ i + 1 ] );
					break;
				case '6':
					a3_mult = atof( argv[ i + 1 ] );
					break;
				case '7':
					b3_mult = atof( argv[ i + 1 ] );
					break;
				case 'p':
					precision = atoi( argv[ i + 1 ] ) != 0 ? atoi( argv[ i + 1 ] ) : 1;
					break;
				case 'b':
					blocksize = atoi( argv[ i + 1 ] );
					break;
				case 'g':
					simplex_growth = atof( argv[ i + 1 ] );
					break;
				case 'c':
					simplex_contract = atof( argv[ i + 1 ] );
					break;
				case 'r':
					simplex_reflect = atof( argv[ i + 1 ] );
					break;
				case 't':
					simplex_iterations = atoi( argv[ i + 1 ] );
					break;
				case 'z':
					simplex_mode = atoi( argv[ i + 1 ] ) != 0 ? true : false;
				default:
					break;
			}
		}
	}
}

/**
 *	Initialize parameters that will be necessary for further grid-search functioning.
 */
void argo::initCalculatedParams () {
	// Set precision adjusted MaxA0 and MaxB0.
	A0_max = ( int ) ceil( ( float ) A0_max / ( float ) precision );
	B0_max = ( int ) ceil( ( float ) B0_max / ( float ) precision );

	// Set a1, a2, a3, b1, b2, b3 step sizes.
	a1_step = 1.0 / rbase_height; //that's 1 pixel
	a2_max = a2_mult * ( a1_max * 4 / rbase_height );  //one third of max deviation from linear part.  This is totally arbitrary.
	a2_step = 4.0 / ( rbase_height * rbase_height );
	a3_max = a3_mult * 20.78 * a1_max / ( rbase_height * rbase_height );
	a3_step = 20.78 / ( rbase_height * rbase_height * rbase_height ); //changed constant from 10.4 on 1/15/2007 <--- Nathan would like to know where magic 20.78 comes from.

	b1_step = 1.0 / rbase_height;
	b2_max = b2_mult * ( b1_max * 4 / rbase_height );
	b2_step = 4.0 / ( rbase_height * rbase_height );
	b3_max = b3_mult * 20.78 * b1_max / ( rbase_height * rbase_height );
	b3_step = 20.78 / ( rbase_height * rbase_height * rbase_height );

	// Make a1, a2, a3, b1, b2, b3 maximums exact ratios (rather than relative) of their respective steps
	a1_max = ceil( a1_max / a1_step ) * a1_step;
	a2_max = ceil( a2_max / a2_step ) * a2_step;
	a3_max = ceil( a3_max / a3_step ) * a3_step;

	b1_max = ceil( b1_max / b1_step ) * b1_step;
	b2_max = ceil( b2_max / b2_step ) * b2_step;
	b3_max = ceil( b3_max / b3_step ) * b3_step;

	// This set of variables record how many pixels of drift in either axis that will be encountered by the program.
	// Note that the factor of 1.38 is hard-coded for a2_multiplier = 1/3, a3_multiplier = 1/9
	A_total_drift = A0_max + ( int ) ceil( 1.38 * a1_max * rbase_height );
	B_total_drift = B0_max + ( int ) ceil( 1.38 * b1_max * rbase_height );

	// Compute maximum slope.
	double max_slope, maxslopeA, maxslopeB;
	maxslopeA = a1_max + a2_max * rbase_height + 0.5 * rbase_height * rbase_height * a3_max;
	maxslopeB = b1_max + b2_max * rbase_height + 0.5 * rbase_height * rbase_height * b3_max;
	max_slope = maxslopeB > maxslopeA ? maxslopeB : maxslopeA;

	// Recalculate block size based upon calculated maximum slope if user does not provide a block size.
	if ( blocksize <= 0 ) {
		blocksize = ( int ) ( 1 / max_slope );
	}

	numblocks = 0;

	A_sliver_drift = ( int ) ( maxslopeA * rsliver_width );
	B_sliver_drift = ( int ) ( maxslopeB * rsliver_width );

	A_difflet_points = ( A_total_drift + A_sliver_drift ) * 2 + 1;  //get factor of two correct here.
	B_difflet_points = ( B_total_drift + B_sliver_drift ) * 2 + 1;

	numblocks = ( rbase_height - B_total_drift - 2 * B_sliver_drift - B0_max ) / blocksize;

	A_points = A_total_drift * 2 + 1;
	B_points = B_total_drift * 2 + 1;
	num_sliver_blocks = rsliver_width / blocksize;
	sliver_block_width = rsliver_width / num_sliver_blocks;
	dyn_diffs_size = numblocks * A_points * B_points;
}

/**
 *	We use Matt's trick here of getting a list of indexes here, even though we have to use a for loop.  Pre-calculating the positions IS faster,
 *	and Matt's trick precomputes as much as possible.
 */
void argo::initCombos () {
	// Make lists of all vaues of A1,A2,A3, B1,B2,B3 to try
	A_combos_size = ( ( int ) ( a1_max / a1_step ) * 2 + 1 ) * ( ( int ) ( a2_max / a2_step ) * 2 + 1 ) * ( ( int ) ( a3_max / a3_step ) * 2 + 1 );
	combos.A_combos = new param_combo[ A_combos_size ];
	long int count = 0;
	for ( double a3 = -a3_max; a3 <= a3_max + a3_step / 2; a3 += a3_step ) {
		double p2_term1 = -1.5 * rbase_height * a3;
		double p1_term1 = 0.5 * rbase_height * rbase_height * a3;
		for ( double a2 = -a2_max; a2 <= a2_max + a2_step / 2; a2 += a2_step ) {
			double p1_term2 = -a2 * rbase_height;
			for ( double a1 = -a1_max; a1 <= a1_max + a1_step / 2; a1 += a1_step ) {
				combos.A_combos[ count ].P1 = a1 + p1_term1 + p1_term2;
				combos.A_combos[ count ].P2 = a2 + p2_term1;
				combos.A_combos[ count ].P3 = a3;
				++count;
			}
		}
	}

	B_combos_size = ( ( int ) ( b1_max / b1_step ) * 2 + 1 ) * ( ( int ) ( b2_max / b2_step ) * 2 + 1 ) * ( ( int ) ( b3_max / b3_step ) * 2 + 1 );
	combos.B_combos = new param_combo[ B_combos_size ];
	count = 0;
	for ( double b3 = -b3_max; b3 <= b3_max + b3_step / 2; b3 += b3_step ) {
		double p2_term1 = -1.5 * rbase_height * b3;
		double p1_term1 = 0.5 * rbase_height * rbase_height * b3;
		for ( double b2 = -b2_max; b2 <= b2_max + b2_step / 2; b2 += b2_step ) {
			double p1_term2 = -b2 * rbase_height;
			for ( double b1 = -b1_max; b1 <= b1_max + b1_step / 2; b1 += b1_step ) {
				combos.B_combos[ count ].P1 = b1 + p1_term1 + p1_term2;
				combos.B_combos[ count ].P2 = b2 + p2_term1;
				combos.B_combos[ count ].P3 = b3;
				++count;
			}
		}
	}

	// Sort generated parameters
	std::sort( combos.A_combos, combos.A_combos + A_combos_size );
	std::sort( combos.B_combos, combos.B_combos + B_combos_size );
}

/**
 *	Initialize the difflets (gamma values) which will later be used to generate beta values.
 */
void argo::initBetaGamma()	{
	// Initialize array for dynamic_diffs once.
	beta_gamma_store.dynamic_diffs = new beta_values[dyn_diffs_size];

	int smin = rbase_width / 2 - rsliver_width / 2;
	double partial_sum1;
	double partial_sum2;
	double partial_sum3;
	double diff;

	data_range = images_store.resamp_base->getRange();

	// Hardcoded for stepping over image. For a byte scale image thus, we will have that z-step which is approximately 1.
	double zstep = data_range / 256;	
	double inv2zstep = 1 / ( 2 * zstep );
	double invzstep2 = 1 / ( zstep * zstep );

	// Figure out maximum sliver drift, make difflets larger by that much.  (Apoints+Asliverpoints)
	beta_gamma_store.difflets = new gamma_values[numblocks * A_difflet_points *
									B_difflet_points * num_sliver_blocks ];

	gamma_values difflet;
	for ( long int sb = 0; sb < num_sliver_blocks; ++sb ) {
		long int sliver_block_x = sb * sliver_block_width;
		for ( long int j = 0; j <= 2 * ( B_total_drift + B_sliver_drift ); ++j ) {
			for ( long int i = 0; i <= 2 * ( A_total_drift + A_sliver_drift ); ++i ) {
				long int x_initial = smin - A_total_drift - A_sliver_drift + i;
				x_initial += sliver_block_x;
				for ( long int b = 0; b < numblocks; b++ ) {
					partial_sum1 = 0;
					partial_sum2 = 0;
					partial_sum3 = 0;
					for ( long int y = B0_max - B_total_drift + j + b * blocksize,
							ys = B0_max + B_sliver_drift + b * blocksize;
							ys < B0_max + ( b + 1 ) * blocksize + B_sliver_drift;
							++y, ++ys ) {
						// Skip all negative elements.
						if ( y < 0 ) {
							continue;
						}
						for ( long int x = x_initial, xs = sliver_block_x;
								xs < sliver_block_x + sliver_block_width;
								++x, ++xs ) {
							diff = -images_store.resamp_sliver->fastGet( xs, ys ) + images_store.resamp_base->fastGet( x, y ) - zstep;
							partial_sum1 += diff * diff;
							diff += zstep;
							partial_sum2 += diff * diff;
							diff += zstep;
							partial_sum3 += diff * diff;
						}
					}
					
					double g2 = ( 0.5 * ( partial_sum1 + partial_sum3 ) - partial_sum2 ) * invzstep2;
					double g1 = ( partial_sum3 - partial_sum1 ) * inv2zstep;
					double g0 = partial_sum2;

					difflet.G0 = g0;
					difflet.G1 = g1;
					difflet.G2 = g2;

					beta_gamma_store.difflets[b + i * numblocks + j * numblocks * A_difflet_points +
					          sb * numblocks * A_difflet_points * B_difflet_points ] = difflet;
				} //end for b
			} //end for i
		} //end for j
	} //end for sb
}

/**
 *	The bulk of grid-search calculations. Given index values, we generate appropriate beta
 *	values using the pre-computed gammas. We then go through our combos range and perform
 *	a matrix multiply to generate C parameters and calculate the final difference.
 *
 *	@param	verbose	Boolean flag indicating whether to display debugging information.
 */
void argo::performGridSearch(bool verbose)	{

	double * yb_arr = new double[numblocks];
	double * yb2_arr = new double[numblocks];
	double * yb3_arr = new double[numblocks];
	double * yb4_arr = new double[numblocks];
	double * yb5_arr = new double[numblocks];
	double * yb6_arr = new double[numblocks];
	int mult1, mult2, mult3;
	mult3 = numblocks * A_difflet_points * B_difflet_points;
	mult2 = numblocks * A_difflet_points;
	mult1 = numblocks * A_points;

	for (int i = 0; i < numblocks; ++i) {
		yb_arr[i] = i * blocksize + B0_max; //+ Bsliverdrift;
		yb2_arr[i] = yb_arr[i] * yb_arr[i];
		yb3_arr[i] = yb2_arr[i] * yb_arr[i];
		yb4_arr[i] = yb3_arr[i] * yb_arr[i];
		yb5_arr[i] = yb4_arr[i] * yb_arr[i];
		yb6_arr[i] = yb5_arr[i] * yb_arr[i];
	}

	tbb::critical_section cs;
	//TO DO for tomorrow: think about how to do initial and final B1 sensibly and symmetrically
	//Refill DIFFS
	long int Aindex_f = 0, Bindex_f = 0;
	for (long int Bindex_i = 0; Bindex_i < B_combos_size - 1; Bindex_i = Bindex_f) {
		double B1_i = combos.B_combos[Bindex_i].P1;
		double B1_f = B1_i + .01;
		//search ahead for index Bi_f
		for (Bindex_f = Bindex_i + 1; (combos.B_combos[Bindex_f].P1 < B1_f) && (Bindex_f < B_combos_size); ++Bindex_f);

		for (long int Aindex_i = 0; Aindex_i < A_combos_size - 1; Aindex_i = Aindex_f) {
			double A1_i = combos.A_combos[Aindex_i].P1;
			double A1_f = A1_i + .01;
			//search ahead for index Ai_f
			for (Aindex_f = Aindex_i + 1; (combos.A_combos[Aindex_f].P1 < A1_f) && (Aindex_f < A_combos_size); ++Aindex_f);
			
			// Empty out array of dynamic diffs.
			memset(beta_gamma_store.dynamic_diffs, 0, dyn_diffs_size * sizeof(beta_values));

			// Populate array with new beta values.
			// **** NOTE: Parallelize in the future. The loops should be decoupled enough to allow for this. Determine which of the loops
			// would be best to parallelize. ****
			// long int diffs_address;
			// long int difflets_address;
			// int Xc = -sliver_block_width;
			// int totX;
			// int totY;
			// int finalAddress = -mult3;
			for (int sb = 0; sb < num_sliver_blocks; ++sb) {
				//Now calculate offsets sx and symeh
				//finalAddress += mult3;
				//Xc += sliver_block_width;
				int finalAddress = mult3 * sb;
				int Xc = sliver_block_width * sb;
				int Xc2 = Xc * Xc;
				int doubleXc = 2 * Xc;
				int sx = A_sliver_drift - (int)(Xc * A1_i);  //which A1 do I use?  first?  last?  average?
				int sy = B_sliver_drift - (int)(Xc * B1_i);
				//printf("%f, %d, %d, %d\n",A1_i,sb,sx,sy);
				//int YPart = -mult1;
				for (int dy = 0; dy < B_points; ++dy) {
					int totY = ( dy + sy ) * mult2;
					//YPart += mult1;
					int YPart = mult1 * dy;
					//int XPart = -numblocks;
					for (int dx = 0; dx < A_points; ++dx) {
						int totX = (dx + sx) * numblocks;
						int XPart = ( ( int ) numblocks ) * dx;
						//XPart += (int)numblocks;
						for (int b = 0; b < numblocks; ++b) {  //partially unroll this loop?
							//diffs[dy,dx,b] += difflets[sb,dy+sy,dx+sx,b]

							//Xc=0;
							//printf("X position: %f \n",Xc);

							//difflets_address =(b + (dx+sx)*numblocks + (dy+sy)*numblocks*Adiffletpoints + sb*numblocks*Adiffletpoints*Bdiffletpoints);
							long int difflets_address = ( b + totX + totY + finalAddress );
							//diffs_address = b+dx*numblocks+dy*numblocks*Apoints;
							long int diffs_address = b + XPart + YPart;
							beta_gamma_store.dynamic_diffs[diffs_address].B0 += beta_gamma_store.difflets[difflets_address].G0;
							beta_gamma_store.dynamic_diffs[diffs_address].B1 += beta_gamma_store.difflets[difflets_address].G1;
							beta_gamma_store.dynamic_diffs[diffs_address].B2 += beta_gamma_store.difflets[difflets_address].G2;
							beta_gamma_store.dynamic_diffs[diffs_address].B3 += (beta_gamma_store.difflets[difflets_address].G1 * Xc);
							beta_gamma_store.dynamic_diffs[diffs_address].B4 += beta_gamma_store.difflets[difflets_address].G2 * Xc2;
							beta_gamma_store.dynamic_diffs[diffs_address].B5 += (beta_gamma_store.difflets[difflets_address].G2 * doubleXc);
							//}
						}
					}
				}
			}

			int blocktimesApoints = A_points * numblocks;
			long A0_B0_adj_inc = numblocks;
			long MaxA0_times_numblocks = A0_max * numblocks;
			long B0_adj_min = -B0_max * blocktimesApoints;
			long B0_adj_max = B0_max * blocktimesApoints;
			long B0_adj_inc = blocktimesApoints;

			tbb::parallel_for((long int)Bindex_i, Bindex_f, (long int)1, [&](long int Bindex)
			{
				for (long int Aindex = Aindex_i; Aindex < Aindex_f; ++Aindex) {
					double B1 = combos.B_combos[Bindex].P1;
					double B2 = combos.B_combos[Bindex].P2;
					double B3 = combos.B_combos[Bindex].P3;
					double A1 = combos.A_combos[Aindex].P1;
					double A2 = combos.A_combos[Aindex].P2;
					double A3 = combos.A_combos[Aindex].P3;

					for (long B0_adj = B0_adj_min; B0_adj <= B0_adj_max; B0_adj += B0_adj_inc) {
						long A0_B0_adj_max = MaxA0_times_numblocks + B0_adj;

						for (long A0_B0_adj = -MaxA0_times_numblocks + B0_adj; A0_B0_adj <= A0_B0_adj_max; A0_B0_adj += A0_B0_adj_inc) {

							Eigen::Matrix4d A;
							Eigen::Vector4d B;
							Eigen::Vector4d C;
							Eigen::Matrix4d D;
							Eigen::Vector4d E;

							// A and B are the only ones necessary to be zero at the start due to summations. The rest get overwritten.
							A << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
							B << 0, 0, 0, 0;

							for (int i = 0; i<numblocks; i++) {
								double beta1, double_beta2, beta3, beta4, beta5;

								long xpos = (long)(A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + A_total_drift);
								long ypos = (long)(B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + B_total_drift);

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
							

							
							unsigned long long ignored = 0;
							double function = 0;
							for (int i = 0; i<numblocks; i++){
								long xpos = (long)(A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + A_total_drift);
								long ypos = (long)(B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + B_total_drift);
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

								if (sum <= 0)	{
									printf("Negative or zero sum!\n");
								}

								function += sum;
								// Ignored is for debugging purposes. Remove in the future.
								ignored += function > results.bestdiff ? numblocks - ( i + 1 ) : 0;
							}

							// Critical section lock to prevent data overwriting.
							cs.lock();
							// Note: this "if" block doesn't seem to slow down the program at all, compared to simpler statements.
							// If we have generated a new minimum, overwrite values in results struct.
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
							results.iterations_ignored += ignored;
							// Unlock here.
							cs.unlock();
						}
					}
				}
			} );
		}
	}
	delete[] yb_arr;
	delete[] yb2_arr;
	delete[] yb3_arr;
	delete[] yb4_arr;
	delete[] yb5_arr;
	delete[] yb6_arr;
}

void argo::finalizeGridSearch() {
	delete[] beta_gamma_store.difflets;
	delete[] beta_gamma_store.dynamic_diffs;
	delete[] combos.A_combos;
	delete[] combos.B_combos;

	beta_gamma_store.difflets = 0;
	beta_gamma_store.dynamic_diffs = 0;
	combos.A_combos = 0;
	combos.B_combos = 0;
}

/**
 *	Makes the call to simplex, taking in all values calculated by grid-search previously and passing them over along with
 *	corresponding precision values to help simplex function properly. Performs fast-z or slow-z correction depending
 *	upon user inputted value for simplex_mode. Defaults to slow-z correction.
 */
void argo::performSimplexRoutine () {
	// Normalize z-correction parameters with the precision that was applied to the image.
	if ( !simplex_mode ) {
		results.bestC1 /= precision;
		results.bestC2 /= precision * precision;
		results.bestC3 /= precision * precision * precision;

		grid_best[ 8 ] = results.bestC0;
		grid_best[ 9 ] = results.bestC1;
		grid_best[ 10 ] = results.bestC2;
		grid_best[ 11 ] = results.bestC3;

		precisionArr[ 8 ] = 0.1 * data_range;
		precisionArr[ 9 ] = 0.1 * data_range / ( double ) images_store.orig_base->height;
		precisionArr[ 10 ] = 0.1 * data_range / ( double ) ( images_store.orig_base->height * images_store.orig_base->height );
		precisionArr[ 11 ] = 0.1 * data_range / ( double ) ( images_store.orig_base->height * images_store.orig_base->height * images_store.orig_base->height );

		grid_best[ 12 ] = results.bestA1;
		grid_best[ 13 ] = results.bestB1 + 1;
		grid_best[ 14 ] = results.bestC1;

		precisionArr[ 12 ] = 0.1 / ( double ) images_store.orig_base->height;
		precisionArr[ 13 ] = 0.1 / ( double ) images_store.orig_base->height;
		precisionArr[ 14 ] = 0.1 * data_range;


	}
	else {
		grid_best[ 8 ] = results.bestA1;
		grid_best[ 9 ] = results.bestB1 + 1;

		precisionArr[ 8 ] = 0.1 / ( double ) images_store.orig_base->height;
		precisionArr[ 9 ] = 0.1 / ( double ) images_store.orig_base->height;
	}



	// These are how you should fill the x array. Use this initialization 
	// if you haven't skipped the grid search.
	// Normalize parameters with regards to precision.
	grid_best[ 0 ] = results.bestA0 * precision;
	grid_best[ 1 ] = results.bestB0 * precision;
	grid_best[ 2 ] = results.bestA1;
	grid_best[ 3 ] = results.bestB1 + 1;
	grid_best[ 4 ] = results.bestA2 / precision;
	grid_best[ 5 ] = results.bestB2 / precision;
	grid_best[ 6 ] = results.bestA3 / ( precision * precision );
	grid_best[ 7 ] = results.bestB3 / ( precision * precision );


	for ( int i = 0; i < NUM_PARAMS; ++i ) {
		simplex_best[ i ] = grid_best[ i ];
	}
	

	/*
	*	Precision array determines by how much the Nelder-Mead simplex routine adjusts individual parameters.
	*	These values are all based on the changes in the parameters that will cause at most a shift of 0.1 pixel.
	*/
	precisionArr[ 0 ] = 0.1;
	precisionArr[ 1 ] = 0.1;
	precisionArr[ 2 ] = 0.1 / ( double ) images_store.orig_base->height;
	precisionArr[ 3 ] = 0.1 / ( double ) images_store.orig_base->height;
	precisionArr[ 4 ] = 0.1 / ( double ) ( images_store.orig_base->height * images_store.orig_base->height );
	precisionArr[ 5 ] = 0.1 / ( double ) ( images_store.orig_base->height * images_store.orig_base->height );
	precisionArr[ 6 ] = 0.1 / ( double ) ( images_store.orig_base->height * images_store.orig_base->height * images_store.orig_base->height );
	precisionArr[ 7 ] = 0.1 / ( double ) ( images_store.orig_base->height * images_store.orig_base->height * images_store.orig_base->height );

	simplex( images_store.orig_base, images_store.orig_sliver, grid_best, simplex_best, simplex_mode, precisionArr, simplex_reflect, simplex_contract, simplex_growth,
			simplex_iterations );
}

/**
 *	Use only for slow-z correction. This method is currently deprecated for fast-z correction. Image writing is done in 
 *	simplex routine after r and c parameters are generated.
 *
 *	Performs final warping and outputting of images as files. The operation of this method depends on the global array simplex_best
 *	being initialized and filled with the optimized parameters necessary for image correction.
 */
void argo::performImageCorrection () {
	// Perform final warp, and write the output tiff
	double aterms[ 4 ] = { simplex_best[ 0 ], simplex_best[ 2 ], simplex_best[ 4 ], simplex_best[ 6 ] };
	double bterms[ 4 ] = { simplex_best[ 1 ], simplex_best[ 3 ], simplex_best[ 5 ], simplex_best[ 7 ] };
	double cterms[ 4 ] = { simplex_best[ 8 ], simplex_best[ 9 ], simplex_best[ 10 ], simplex_best[ 11 ] };

	FImage* warp_base = new FImage( images_store.orig_base->width, images_store.orig_base->height, images_store.orig_base->metadata );
	images_store.orig_base->warpBase( warp_base, aterms, bterms, cterms, 1, interp_type );

	aterms[ 1 ] = simplex_best[ 12 ];
	bterms[ 1 ] = simplex_best[ 13 ];
	cterms[ 1 ] = simplex_best[ 14 ];

	FImage* warp_sliver = new FImage( images_store.orig_sliver->width, images_store.orig_sliver->height, images_store.orig_sliver->metadata );
	images_store.orig_sliver->warpSliver( warp_sliver, aterms, bterms, cterms, interp_type );
	
	// Write corrected images.
	warp_base->writeImage( "w_base.tif" );
	warp_sliver->writeImage( "w_sliver.tif" );

	double base_min = warp_base->getMin();
	double base_max = warp_base->getMax();
	double sliver_min = warp_sliver->getMin();
	double sliver_max = warp_sliver->getMax();
	double min_val = base_min < sliver_min ? base_min : sliver_min;
	double max_val = base_max < sliver_max ? base_max : sliver_max;
	double slope = 1.0 / ( max_val - min_val );

	// Write viewable images.
	warp_base->writeDisplayableImage( "dw_base.tif", slope, min_val );
	warp_sliver->writeDisplayableImage( "dw_sliver.tif", slope, min_val );





	// Write viewable original images.
	base_min = images_store.orig_base->getMin();
	base_max = images_store.orig_base->getMax();
	sliver_min = images_store.orig_sliver->getMin();
	sliver_max = images_store.orig_sliver->getMax();
	min_val = base_min < sliver_min ? base_min : sliver_min;
	max_val = base_max < sliver_max ? base_max : sliver_max;
	slope = 1.0 / ( max_val - min_val );

	images_store.orig_base->writeDisplayableImage( "orig_base.tif", slope, min_val );
	images_store.orig_sliver->writeDisplayableImage( "orig_sliver.tif", slope, min_val );

	delete warp_base;
	delete warp_sliver;
}

/**
 *	Non-debugging versions of argo execution.
 *
 *	@param	argc		The amount of command line arguments.
 *	@param	argv		The command line arguments.
 */
void argo::correctImages ( int argc, char* argv[] ) {
	readInputParams( argc, argv );
	readImages( false );
	initCalculatedParams();
	initCombos();
	initBetaGamma();
	performGridSearch( false );
	finalizeGridSearch();
	performSimplexRoutine();
	performImageCorrection();
}

/**
 *	Debugging version of argo execution.
 *
 *	@param	argc		The amount of command line arguments.
 *	@param	argv		The command line arguments.
 *	@param	verbose	Boolean flag indicating whether to display debugging information.
 */
void argo::correctImages ( int argc, char* argv[], bool verbose ) {
	clock_t begintime = clock();
	clock_t time0, time1;
	bool debugmode = false;

	FreeImage_Initialise();
	time0 = clock();
	readInputParams( argc, argv );
	if ( debugmode ) {
		precision = 1;
		blocksize = 8;
	}
	if ( verbose ) {
		logInputParams();
	}

	readImages( verbose );
	initCalculatedParams();
	time1 = clock();
	times.image_read_time = ( ( double ) time1 - ( double ) time0 ) / CLOCKS_PER_SEC;

	time0 = clock();
	initCombos();
	time1 = clock();
	times.combos_time = ( ( double ) time1 - ( double ) time0 ) / CLOCKS_PER_SEC;

	time0 = clock();
	initBetaGamma();
	time1 = clock();
	times.diffs_time = ( ( double ) time1 - ( double ) time0 ) / CLOCKS_PER_SEC;

	if ( verbose ) {
		logCalculatedParams();
		logBetaGammaInfo();
	}

	time0 = clock();
	performGridSearch( verbose );
	time1 = clock();
	times.grid_time = ( ( double ) time1 - ( double ) time0 ) / CLOCKS_PER_SEC;


	if ( verbose ) {
		logGridSearchInfo();
	}

	time0 = clock();
	performSimplexRoutine();
	finalizeGridSearch();
	time1 = clock();
	times.simplex_time = ( ( double ) time1 - ( double ) time0 ) / CLOCKS_PER_SEC;

	if ( verbose ) {
		logSimplexRoutineInfo();
	}

	if ( !simplex_mode ) {
		time0 = clock();
		performImageCorrection();
		time1 = clock();
		times.image_write_time = ( ( double ) time1 - ( double ) time0 ) / CLOCKS_PER_SEC;
	}

	clock_t endtime = clock();
	times.total_time = ( ( double ) endtime - ( double ) begintime ) / CLOCKS_PER_SEC;

	FreeImage_DeInitialise();
	if ( verbose ) {
		logProgramInformation();
	}
}