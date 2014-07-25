
#define __MAIN_RUN__
#define __OPENCL_TEST__


/**
 *	Calls Argo to perform its functionality.
 */
#ifdef __MAIN_RUN__
	#undef __OPENCL_TEST__
	#include <iostream>
	#include "argo.h"

	/**
	*	Main method for our program.
	*/
	int main( int argc, char *argv[] ) {
		printf( "************ Running Argo Version " VERSION " ************\n" );
		argo* program = new argo();
		program->correctImages( argc, argv, true );
		delete program;
		getchar();
	}
#endif


/**
 *	For OpenCL testing purposes.
 */
#ifdef __OPENCL_TEST__

#define _CRT_SECURE_NO_WARNINGS

#include <CL\cl.hpp>
#include <iostream>
#include <time.h>
#include <tbb\parallel_for.h>
#include <tbb\critical_section.h>
#include <fstream>
#include <vector>
#include "FImage.h"
#include "Eigen\Eigen"
#include "Kernels.h"

#define __CODEFOLD__
#define ITERATIONS 1
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


unsigned int nextPow2( unsigned int x ) {
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

int main( int argc, char *argv[] ) {

	// PROGRAM INITIALIZATION
#ifdef __CODEFOLD__
	int interpolation_type = F_CUBIC;
	int write_interp = F_BSPLINE;
	bool flipped;
	char command_arg;
	std::string base_name, sliver_name;
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
					flipped = atoi( argv[ i + 1 ] ) != 0 ? true : false;
					break;
				default:
					break;
			}
		}
	}

	// base_name = "base.tif";
	// sliver_name = "sliver.tif";

	FImage* base = new FImage( base_name, flipped );
	FImage* sliver = new FImage( sliver_name, flipped );
	double A0, A1, A2, A3, B0, B1, B2, B3;


	

	A0 = -2.000000e+000;
	A1 = -0.000000e+000;
	A2 = -8.676941e-006;
	A3 = 4.918900e-009;
	B0 = -2.000000e+000;
	B1 = -1.700680e-003;
	B2 = 5.784627e-006;
	B3 = 1.475670e-008;
	
	
	/*
	A0 = 3.000000e+000;
	A1 = -3.125000e-002;
	A2 = -1.356337e-005;
	A3 = -3.532127e-008;
	B0 = -2.000000e+000;
	B1 = -5.208333e-003;
	B2 = -3.390842e-005;
	B3 = -3.532127e-008;
	*/

	/*
	A0 = -3.012342314235345311111111111111e+000;
	A1 = 3.1250345345345001111111111111111111e-002;
	A2 = 1.3563453453453337111111111111111111111e-005;
	A3 = 3.5323453453451271111111111111111e-008;
	B0 = 2.000345435345000111111111111111111e+000;
	B1 = 5.208333345345345111111111111111111e-003;
	B2 = 3.390842345345345341111111111111111111e-005;
	B3 = 3.532127654654111111111111111111e-008;
	*/
	
	
	double* As = new double[4] { A0, A1, A2, A3 };
	double* Bs = new double[4] { B0, 1 + B1, B2, B3 };
	double* Cs = new double[4] { 0, 0, 0, 0 };
	double *base_array = new double[ sliver->width * sliver->height ];
	double *sliver_array = new double[ sliver->width * sliver->height ];

	double* double_sliv = new double[ sliver->width * sliver->height ];
	double* double_base = new double[ base->width * base->height ];
	for ( int i = 0; i < sliver->width * sliver->height; ++i ) {
		double_sliv[ i ] = sliver->data[ i ];
	}
	for ( int i = 0; i < base->width * base->height; ++i ) {
		double_base[ i ] = base->data[ i ];
	}
	

	int xbp = base->width / 2 - sliver->width / 2;
	bool exists = false;
	int maxX = sliver->width - 1;
	int sliver_minY = 0;
	int sliver_maxY = sliver->height - 1;
	int base_minY = 0;
	int base_maxY = base->height - 1;
	while ( !exists ) {
		if ( maxX - As[ 1 ] * maxX > sliver->width - 1 ) {
			maxX--;
		}
		else {
			exists = true;
		}
	}
	exists = false;
	while ( !exists ) {
		bool exist2 = false;
		if ( sliver_maxY + -( Bs[ 1 ] - 1 ) * maxX > sliver->height - 1 ) {
			sliver_maxY--;
		}
		else {
			exist2 = true;
		}
		if ( sliver_minY + -( Bs[ 1 ] - 1 ) * maxX < 0 ) {
			sliver_minY++;
		}
		else {
			if ( exist2 ) {
				exists = true;
			}
		}
	}
	exists = false;
	while ( !exists ) {
		int y2 = base_maxY * base_maxY;
		int y3 = y2 * base_maxY;
		bool exist2 = false;
		if ( Bs[ 0 ] + Bs[ 1 ] * base_maxY + Bs[ 2 ] * y2 + Bs[ 3 ] * y3 > sliver->height - 1 ) {
			base_maxY--;
		}
		else {
			exist2 = true;
		}

		y2 = base_minY * base_minY;
		y3 = y2 * base_minY;
		if ( Bs[ 0 ] + Bs[ 1 ] * base_minY + Bs[ 2 ] * y2 + Bs[ 3 ] * y3 < 0 ) {
			base_minY++;
		}
		else {
			if ( exist2 ) {
				exists = true;
			}
		}
	}

	int maxY, minY;
	double weightTop, weightBottom, weightRight;
	weightRight = ceil( ( maxX - As[ 1 ] * maxX ) ) - ( maxX - As[ 1 ] * maxX );
	if ( base_maxY < sliver_maxY ) {
		maxY = base_maxY;	// Ceil for weight
		weightTop = ceil( Bs[ 0 ] + Bs[ 1 ] * maxY + Bs[ 2 ] * maxY + Bs[ 3 ] * maxY ) - ( Bs[ 0 ] + Bs[ 1 ] * maxY + Bs[ 2 ] * maxY + Bs[ 3 ] * maxY );
	}
	else if ( base_maxY > sliver_maxY ) {
		maxY = sliver_maxY;	// Ceil for weight
		weightTop = ceil( maxY - ( Bs[ 1 ] - 1 ) * maxX ) - ( maxY - ( Bs[ 1 ] - 1 ) * maxX );
	}
	else {
		maxY = sliver_maxY;
		if ( ceil( Bs[ 0 ] + Bs[ 1 ] * maxY + Bs[ 2 ] * maxY + Bs[ 3 ] * maxY ) - ( Bs[ 0 ] + Bs[ 1 ] * maxY + Bs[ 2 ] * maxY + Bs[ 3 ] * maxY ) <
			 ceil( maxY - ( Bs[ 1 ] - 1 ) * maxX ) - ( maxY - ( Bs[ 1 ] - 1 ) * maxX ) ) {
			weightTop = ceil( Bs[ 0 ] + Bs[ 1 ] * maxY + Bs[ 2 ] * maxY + Bs[ 3 ] * maxY ) - ( Bs[ 0 ] + Bs[ 1 ] * maxY + Bs[ 2 ] * maxY + Bs[ 3 ] * maxY );
		}
		else {
			weightTop = ceil( maxY - ( Bs[ 1 ] - 1 ) * maxX ) - ( maxY - ( Bs[ 1 ] - 1 ) * maxX );
		}
	}

	if ( base_minY > sliver_minY ) {
		minY = base_minY;	// Floor for weight
		weightBottom = ( Bs[ 0 ] + Bs[ 1 ] * minY + Bs[ 2 ] * minY + Bs[ 3 ] * minY ) - floor( ( Bs[ 0 ] + Bs[ 1 ] * minY + Bs[ 2 ] * minY + Bs[ 3 ] * minY ) );
	}
	else if ( base_minY < sliver_minY ) {
		minY = sliver_minY;	// Floor for weight
		weightBottom = ( minY - ( Bs[ 1 ] - 1 ) * maxX ) - floor( ( minY - ( Bs[ 1 ] - 1 ) * maxX ) );
	}
	else {
		minY = sliver_minY;
		if ( ( minY - ( Bs[ 1 ] - 1 ) * maxX ) - floor( ( minY - ( Bs[ 1 ] - 1 ) * maxX ) ) <
			 ( Bs[ 0 ] + Bs[ 1 ] * minY + Bs[ 2 ] * minY + Bs[ 3 ] * minY ) - floor( ( Bs[ 0 ] + Bs[ 1 ] * minY + Bs[ 2 ] * minY + Bs[ 3 ] * minY ) ) ) {
			weightBottom = ( minY - ( Bs[ 1 ] - 1 ) * maxX ) - floor( ( minY - ( Bs[ 1 ] - 1 ) * maxX ) );
		}
		else {
			weightBottom = ( Bs[ 0 ] + Bs[ 1 ] * minY + Bs[ 2 ] * minY + Bs[ 3 ] * minY ) - floor( ( Bs[ 0 ] + Bs[ 1 ] * minY + Bs[ 2 ] * minY + Bs[ 3 ] * minY ) );
		}
	}

	int m = maxX + 1;
	int n = maxY - minY + 1;
	int size = m - 1 + n;
	double inv_n = 1 / ( double ) n;
	double inv_m = 1 / ( double ) m;
	double inv_mn = 1 / ( double ) ( m * n );

	double* KRCs = new double[ sliver->width + sliver->height ];
	double* RCs = new double[ sliver->width + sliver->height ];
	tbb::critical_section cs;
#endif


	// TBB (CPU THREADED) CALLS
#ifdef __CODEFOLD__


	printf("Starting TBB calls!\n");

	// Warp Sliver Bicubic
	tbb::parallel_for( 0, maxX + 1, 1, [&]( int x )	{
	// for ( int x = 0; x < maxX + 1; x++ ) {
		double sx[ 4 ]; // these are the distances for each row
		double sy[ 4 ]; // these are the distances for each column
		double ux[ 4 ]; // these are the weight for the rows
		double uy[ 4 ]; // these are the weights for the columns

		double yoffset = -( Bs[ 1 ] - 1 ) * x;
		double origX = x - As[ 1 ] * x;

		/////////////////////////////////////////////////////////////////////////////////////

		sx[ 0 ] = abs( origX - floor( origX - 1 ) ); // these get the distances for each row and column from the initial point, positive
		sy[ 0 ] = abs( yoffset - floor( yoffset - 1 ) );
		sx[ 1 ] = abs( 1 - sx[ 0 ] );
		sy[ 1 ] = abs( 1 - sy[ 0 ] );
		sx[ 2 ] = abs( 1 - sx[ 1 ] );
		sy[ 2 ] = abs( 1 - sy[ 1 ] );
		sx[ 3 ] = 1 + sx[ 2 ];
		sy[ 3 ] = 1 + sy[ 2 ];

		for ( int j = 0; j < 4; ++j ) { // these, dependent upon the distance, implement the polynomial weighting
			double x2 = sx[ j ] * sx[ j ];
			double x3 = x2 * sx[ j ];
			double y2 = sy[ j ] * sy[ j ];
			double y3 = y2 * sy[ j ];
			if ( sx[ j ] <= 1 && sx[ j ] >= 0 ) {
				ux[ j ] = ( S_INTERP_CONST_A + 2 ) * x3 - ( S_INTERP_CONST_A + 3 ) * x2 + 1;
			}
			else if ( sx[ j ] > 1 && sx[ j ] <= 2 ) {
				ux[ j ] = S_INTERP_CONST_A * x3 - 5 * S_INTERP_CONST_A * x2 + 8 * S_INTERP_CONST_A * sx[ j ] - 4 * S_INTERP_CONST_A;
			}
			else {
				ux[ j ] = 0;
			}
			if ( sy[ j ] <= 1 && sy[ j ] >= 0 ) {
				uy[ j ] = ( S_INTERP_CONST_A + 2 ) * y3 - ( S_INTERP_CONST_A + 3 ) * y2 + 1;
			}
			else if ( sy[ j ] > 1 && sy[ j ] <= 2 ) {
				uy[ j ] = S_INTERP_CONST_A * y3 - 5 * S_INTERP_CONST_A * y2 + 8 * S_INTERP_CONST_A * sy[ j ] - 4 * S_INTERP_CONST_A;
			}
			else {
				uy[ j ] = 0;
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////
		for ( int y = minY; y < maxY + 1; ++y ) {
		// tbb::parallel_for ( minY, maxY + 1, 1, [ &] ( int y ) {
			double pixelvalue = 0;
			double origY = y + yoffset;
			if ( origX < 0 || origX > sliver->width - 1 || origY < 0 || origY > sliver->height - 1 ) {
				printf( "We hit an invalid value during warping of sliver. (%d, %d) -> (%f, %f)!\n", x, y, origX, origY );
				sliver_array[ x + sliver->width * y ] = -std::numeric_limits<double>::infinity();
			}
			else {
				if ( floor ( origX ) == 0 || floor ( origY ) == 0 || floor ( origX ) >= sliver->width - 2 || floor ( origY ) >= sliver->height - 2 ) {
					pixelvalue = sliver->interpPixel( origX, origY, F_BILINEAR );
				}
				else {
					for ( int j = 0; j < 4; ++j ) {
						for ( int i = 0; i < 4; ++i ) {
							// pixelvalue += double_sliv[ ( int ) ( floor( origX ) - 1 + i + sliver->width * ( floor( origY ) - 1 + j ) ) ] * ux[ i ] * uy[ j ];
							pixelvalue += sliver->get( floor( origX ) - 1 + i, floor( origY ) - 1 + j ) * ux[ i ] * uy[ j ];
						}
					}
				}
				sliver_array[ x + sliver->width * y ] = pixelvalue;
			}
		}
	} );
	// Warp Base Bicubic
	tbb::parallel_for( minY, maxY + 1, 1, [&]( int y )	{
	// for ( int y = minY; y < maxY + 1; ++y ) {
		double sx[ 4 ]; // these are the distances for each row
		double sy[ 4 ]; // these are the distances for each column
		double ux[ 4 ]; // these are the weight for the rows
		double uy[ 4 ]; // these are the weights for the columns

		int ym2 = y * y;
		int ym3 = ym2 * y;
		double xoffset = As[ 0 ] + As[ 1 ] * y + As[ 2 ] * ym2 + As[ 3 ] * ym3;
		double origY = Bs[ 0 ] + Bs[ 1 ] * y + Bs[ 2 ] * ym2 + Bs[ 3 ] * ym3;
		/////////////////////////////////////////////////////////////////////////////////////

		sx[ 0 ] = abs( xoffset - floor( xoffset - 1 ) ); // these get the distances for each row and column from the initial point, positive
		sy[ 0 ] = abs( origY - floor( origY - 1 ) );
		sx[ 1 ] = abs( 1 - sx[ 0 ] );
		sy[ 1 ] = abs( 1 - sy[ 0 ] );
		sx[ 2 ] = abs( 1 - sx[ 1 ] );
		sy[ 2 ] = abs( 1 - sy[ 1 ] );
		sx[ 3 ] = 1 + sx[ 2 ];
		sy[ 3 ] = 1 + sy[ 2 ];
		for ( int j = 0; j < 4; ++j ) { // these, dependent upon the distance, implement the polynomial weighting
			double x2 = sx[ j ] * sx[ j ];
			double x3 = x2 * sx[ j ];
			double y2 = sy[ j ] * sy[ j ];
			double y3 = y2 * sy[ j ];
			if ( sx[ j ] <= 1 && sx[ j ] >= 0 ) {
				ux[ j ] = ( S_INTERP_CONST_A + 2 ) * x3 - ( S_INTERP_CONST_A + 3 ) * x2 + 1;
			}
			else if ( sx[ j ] > 1 && sx[ j ] <= 2 ) {
				ux[ j ] = S_INTERP_CONST_A * x3 - 5 * S_INTERP_CONST_A * x2 + 8 * S_INTERP_CONST_A * sx[ j ] - 4 * S_INTERP_CONST_A;
			}
			else {
				ux[ j ] = 0;
			}
			if ( sy[ j ] <= 1 && sy[ j ] >= 0 ) {
				uy[ j ] = ( S_INTERP_CONST_A + 2 ) * y3 - ( S_INTERP_CONST_A + 3 ) * y2 + 1;
			}
			else if ( sy[ j ] > 1 && sy[ j ] <= 2 ) {
				uy[ j ] = S_INTERP_CONST_A * y3 - 5 * S_INTERP_CONST_A * y2 + 8 * S_INTERP_CONST_A * sy[ j ] - 4 * S_INTERP_CONST_A;
			}
			else {
				uy[ j ] = 0;
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////
		for ( int x = 0; x < maxX + 1; ++x ) {
			double pixelvalue = 0;
			double origX = xbp + x + xoffset;
			if ( origX < 0 || origX > base->width - 1 || origY < 0 || origY > base->height - 1 ) {
				printf( "We hit an invalid value during warping of base!\n" );
				base_array[ x + sliver->width * y ] = -std::numeric_limits<double>::infinity();
			}
			else {
				if ( floor ( origX ) == 0 || floor ( origY ) == 0 || floor ( origX ) >= base->width - 2 || floor ( origY ) >= base->height - 2 ) {
					pixelvalue = base->interpPixel( origX, origY, F_BILINEAR );
				}
				else {
					for ( int j = 0; j < 4; ++j ) {
						for ( int i = 0; i < 4; ++i ) {
							pixelvalue += base->get( floor( origX ) - 1 + i, floor( origY ) - 1 + j ) * ux[ i ] * uy[ j ];
						}
					}
				}
				base_array[ x + sliver->width * y ] = pixelvalue;
			}
		}
	} );
	
	double skr = 0;
	double skr_error = 0;
	// Column Sums
	tbb::parallel_for( 0, m, 1, [&]( int i )	{
	// for ( int i = 1; i < m; i++ ) {
		double sum = 0;
		double sum_error = 0;
		for ( int j = 0; j < n; ++j ) {
			if ( base_array[ i + sliver->width * ( j + minY ) ] == -std::numeric_limits<double>::infinity () ) {
				printf ( "We hit an invalid on base during calculating column sums!\n" );
			}
			if ( sliver_array[ i + sliver->width * ( j + minY ) ] == -std::numeric_limits<double>::infinity () ) {
				printf ( "We hit an invalid on sliver during calculating column sums!\n" );
			}
			double diff = base_array[ i + sliver->width * ( j + minY ) ] - sliver_array[ i + sliver->width * ( j + minY ) ];
			double y = diff - sum_error;
			double t = sum + y;
			sum_error = ( t - sum ) - y;
			sum = t;
		}
		KRCs[ i ] = sum;
	} );

	// Row Sums
	tbb::parallel_for( 0, n, 1, [&]( int j )	{
	// for ( int j = 0; j < n; ++j ) {
		double sum = 0;
		double sum_error = 0;
		for ( int i = 0; i < m; ++i ) {
			if ( base_array[ i + sliver->width * ( j + minY ) ] == -std::numeric_limits<double>::infinity () ) {
				printf ( "We hit an invalid on base during calculating row sums!\n" );
			}
			if ( sliver_array[ i + sliver->width * ( j + minY ) ] == -std::numeric_limits<double>::infinity () ) {
				printf ( "We hit an invalid on sliver during calculating row sums!\n" );
			}
			double diff = sliver_array[ i + sliver->width * ( j + minY ) ] - base_array[ i + sliver->width * ( j + minY ) ];
			double y = diff - sum_error;
			double t = sum + y;
			sum_error = ( t - sum ) - y;
			sum = t;
		}
		KRCs[ m + j ] = sum;
	} );

	// KRC Sum
	tbb::parallel_for( 0, n, 1, [&]( int j )	{
	// for ( int j = 0; j < n; ++j ) {
		cs.lock();
		double y = KRCs[ m + j ] - skr_error;
		double t = skr + y;
		skr_error = ( t - skr ) - y;
		skr = t;
		cs.unlock();
	} );

	// Matrix Computation
	tbb::parallel_for( 1, m, 1, [&]( int i )	{
		//for ( int i = 1; i < m; ++i ) {
		RCs[ i - 1 ] = inv_n * ( KRCs[ i ] - KRCs[ 0 ] );
	} );

	tbb::parallel_for( m, size + 1, 1, [&]( int i ) {
		// for ( int i = m; i <= size; ++i ) {
		RCs[ i - 1 ] = -inv_n * KRCs[ 0 ] - inv_mn * skr + inv_m * KRCs[ i ];
	} );

	// Warp Sliver B-Spline
	tbb::parallel_for( 0, maxX + 1, 1, [&]( int x )	{
	// for ( int x = 0; x < maxX + 1; x++ ) {
		double sx[ 4 ]; // these are the distances for each row
		double sy[ 4 ]; // these are the distances for each column
		double ux[ 4 ]; // these are the weight for the rows
		double uy[ 4 ]; // these are the weights for the columns

		double yoffset = -( Bs[ 1 ] - 1 ) * x;
		double origX = x - As[ 1 ] * x;

		/////////////////////////////////////////////////////////////////////////////////////

		sx[ 0 ] = abs( origX - floor( origX - 1 ) ); // these get the distances for each row and column from the initial point, positive
		sy[ 0 ] = abs( yoffset - floor( yoffset - 1 ) );
		sx[ 1 ] = abs( 1 - sx[ 0 ] );
		sy[ 1 ] = abs( 1 - sy[ 0 ] );
		sx[ 2 ] = abs( 1 - sx[ 1 ] );
		sy[ 2 ] = abs( 1 - sy[ 1 ] );
		sx[ 3 ] = 1 + sx[ 2 ];
		sy[ 3 ] = 1 + sy[ 2 ];



		for ( int j = 0; j < 4; ++j ) { // these, dependent upon the distance, implement the polynomial weighting
			double x2 = sx[ j ] * sx[ j ];
			double x3 = x2 * sx[ j ];
			double y2 = sy[ j ] * sy[ j ];
			double y3 = y2 * sy[ j ];
			if ( sx[ j ] <= 1 && sx[ j ] >= 0 ) {
				ux[ j ] = ( ( 12 - ( 9 * S_INTERP_CONST_B ) - ( 6 * S_INTERP_CONST_C ) ) * x3 + ( -18 + ( 12 * S_INTERP_CONST_B ) + ( 6 * S_INTERP_CONST_C ) )*x2 + 6 - ( 2 * S_INTERP_CONST_B ) ) / 6;
			}
			else if ( sx[ j ] > 1 && sx[ j ] <= 2 ) {
				ux[ j ] = ( ( -S_INTERP_CONST_B - ( 6 * S_INTERP_CONST_C ) )*x3 + ( ( 6 * S_INTERP_CONST_B ) + ( 30 * S_INTERP_CONST_C ) )*x2 + ( -( 12 * S_INTERP_CONST_B ) - ( 48 * S_INTERP_CONST_C ) )*( sx[ j ] ) + ( 8 * S_INTERP_CONST_B ) + ( 24 * S_INTERP_CONST_C ) ) / 6;
			}
			else	{
				ux[ j ] = 0;
			}
			if ( sy[ j ] <= 1 && sy[ j ] >= 0 ) {
				uy[ j ] = ( ( 12 - ( 9 * S_INTERP_CONST_B ) - ( 6 * S_INTERP_CONST_C ) )*y3 + ( -18 + ( 12 * S_INTERP_CONST_B ) + ( 6 * S_INTERP_CONST_C ) )*y2 + 6 - ( 2 * S_INTERP_CONST_B ) ) / 6;
			}
			else if ( sy[ j ] > 1 && sy[ j ] <= 2 ) {
				uy[ j ] = ( ( -S_INTERP_CONST_B - ( 6 * S_INTERP_CONST_C ) )*y3 + ( ( 6 * S_INTERP_CONST_B ) + ( 30 * S_INTERP_CONST_C ) )*y2 + ( -( 12 * S_INTERP_CONST_B ) - ( 48 * S_INTERP_CONST_C ) )*( sy[ j ] ) + ( 8 * S_INTERP_CONST_B ) + ( 24 * S_INTERP_CONST_C ) ) / 6;
			}
			else	{
				uy[ j ] = 0;
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////
		for ( int y = minY; y < maxY + 1; y++ ) {
			double pixelvalue = 0;
			double origY = y + yoffset;
			if ( origX < 0 || origX > sliver->width - 1 || origY < 0 || origY > sliver->height - 1 ) {
				printf( "We hit an invalid value during secondary warping of sliver!\n" );
				sliver_array[ x + sliver->width * y ] = -std::numeric_limits<double>::infinity();
			}
			else {
				if ( floor ( origX ) == 0 || floor ( origY ) == 0 || floor ( origX ) >= sliver->width - 2 || floor ( origY ) >= sliver->height - 2 ) {
					pixelvalue = sliver->interpPixel( origX, origY, F_BILINEAR );
				}
				else {
					for ( int j = 0; j < 4; ++j ) {
						for ( int i = 0; i < 4; ++i ) {
							pixelvalue += sliver->get( floor( origX ) - 1 + i, floor( origY ) - 1 + j ) * ux[ i ] * uy[ j ];
						}
					}
				}
				sliver_array[ x + sliver->width * y ] = pixelvalue;
			}
		}
	} );

	// Warp Base B-Spline
	tbb::parallel_for( minY, maxY + 1, 1, [&]( int y )	{
	// for ( int y = minY; y < maxY + 1; ++y ) {
		double sx[ 4 ]; // these are the distances for each row
		double sy[ 4 ]; // these are the distances for each column
		double ux[ 4 ]; // these are the weight for the rows
		double uy[ 4 ]; // these are the weights for the columns

		int ym2 = y * y;
		int ym3 = ym2 * y;
		double xoffset = As[ 0 ] + As[ 1 ] * y + As[ 2 ] * ym2 + As[ 3 ] * ym3;
		double origY = Bs[ 0 ] + Bs[ 1 ] * y + Bs[ 2 ] * ym2 + Bs[ 3 ] * ym3;
		/////////////////////////////////////////////////////////////////////////////////////

		sx[ 0 ] = abs( xoffset - floor( xoffset - 1 ) ); // these get the distances for each row and column from the initial point, positive
		sy[ 0 ] = abs( origY - floor( origY - 1 ) );
		sx[ 1 ] = abs( 1 - sx[ 0 ] );
		sy[ 1 ] = abs( 1 - sy[ 0 ] );
		sx[ 2 ] = abs( 1 - sx[ 1 ] );
		sy[ 2 ] = abs( 1 - sy[ 1 ] );
		sx[ 3 ] = 1 + sx[ 2 ];
		sy[ 3 ] = 1 + sy[ 2 ];
		for ( int j = 0; j < 4; ++j ) { // these, dependent upon the distance, implement the polynomial weighting
			double x2 = sx[ j ] * sx[ j ];
			double x3 = x2 * sx[ j ];
			double y2 = sy[ j ] * sy[ j ];
			double y3 = y2 * sy[ j ];
			if ( sx[ j ] <= 1 && sx[ j ] >= 0 ) {
				ux[ j ] = ( ( 12 - ( 9 * S_INTERP_CONST_B ) - ( 6 * S_INTERP_CONST_C ) ) * x3 + ( -18 + ( 12 * S_INTERP_CONST_B ) + ( 6 * S_INTERP_CONST_C ) )*x2 + 6 - ( 2 * S_INTERP_CONST_B ) ) / 6;
			}
			else if ( sx[ j ] > 1 && sx[ j ] <= 2 ) {
				ux[ j ] = ( ( -S_INTERP_CONST_B - ( 6 * S_INTERP_CONST_C ) )*x3 + ( ( 6 * S_INTERP_CONST_B ) + ( 30 * S_INTERP_CONST_C ) )*x2 + ( -( 12 * S_INTERP_CONST_B ) - ( 48 * S_INTERP_CONST_C ) )*( sx[ j ] ) + ( 8 * S_INTERP_CONST_B ) + ( 24 * S_INTERP_CONST_C ) ) / 6;
			}
			else	{
				ux[ j ] = 0;
			}
			if ( sy[ j ] <= 1 && sy[ j ] >= 0 ) {
				uy[ j ] = ( ( 12 - ( 9 * S_INTERP_CONST_B ) - ( 6 * S_INTERP_CONST_C ) )*y3 + ( -18 + ( 12 * S_INTERP_CONST_B ) + ( 6 * S_INTERP_CONST_C ) )*y2 + 6 - ( 2 * S_INTERP_CONST_B ) ) / 6;
			}
			else if ( sy[ j ] > 1 && sy[ j ] <= 2 ) {
				uy[ j ] = ( ( -S_INTERP_CONST_B - ( 6 * S_INTERP_CONST_C ) )*y3 + ( ( 6 * S_INTERP_CONST_B ) + ( 30 * S_INTERP_CONST_C ) )*y2 + ( -( 12 * S_INTERP_CONST_B ) - ( 48 * S_INTERP_CONST_C ) )*( sy[ j ] ) + ( 8 * S_INTERP_CONST_B ) + ( 24 * S_INTERP_CONST_C ) ) / 6;
			}
			else	{
				uy[ j ] = 0;
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////
		for ( int x = 0; x < maxX + 1; ++x ) {
			double pixelvalue = 0;
			double origX = xbp + x + xoffset;
			if ( origX < 0 || origX > base->width - 1 || origY < 0 || origY > base->height - 1 ) {
				printf( "We hit an invalid value during secondary warping of base!\n" );
				base_array[ x + sliver->width * y ] = -std::numeric_limits<double>::infinity();
			}
			else {
				if ( floor ( origX ) == 0 || floor ( origY ) == 0 || floor ( origX ) >= base->width - 2 || floor ( origY ) >= base->height - 2 ) {
					pixelvalue = base->interpPixel( origX, origY, F_BILINEAR );
				}
				else {
					for ( int j = 0; j < 4; ++j ) {
						for ( int i = 0; i < 4; ++i ) {
							pixelvalue += base->get( floor( origX ) - 1 + i, floor( origY ) - 1 + j ) * ux[ i ] * uy[ j ];
						}
					}
				}
				base_array[ x + sliver->width * y ] = pixelvalue;
			}
		}
	} );
		
	
	// printf( "Total running time for TBB base and sliver warping at %d iterations : %f\n", ITERATIONS, ( ( double ) end - ( double ) start ) / CLOCKS_PER_SEC );
	printf( "\n" );

	double sum = 0;
	double area = 0;
	double sum_error = 0;
	double area_error = 0;


	double* old_sums = new double[ sliver->height ];
	double* old_weights = new double[ sliver->height ];
	tbb::parallel_for( 0, maxX + 1, 1, [&]( int i )	{
	// for ( int i = 0; i < maxX + 1; ++i ) {
		double c_change = i == 0 ? 0 : RCs[ i - 1 ];

		double solo_sum = 0;
		double solo_sum_error = 0;
		double solo_weight = 0;
		double solo_weight_error = 0;

		for ( int j = minY; j <= maxY; ++j ) {
			double r_change = RCs[ maxX + j - minY ];
			
			if ( base_array[ i + sliver->width * j ] == -std::numeric_limits<double>::infinity()
				 || sliver_array[ i + sliver->width * j ] == -std::numeric_limits<double>::infinity() ) {
				printf( "We hit an invalid value during difference calculation!\n" );
				continue;
			}
			double weight = 1.0;
			if ( i == maxX ) {
				weight *= weightRight;
			}
			if ( j == maxY ) {
				weight *= weightTop;
			}
			else if ( j == minY ) {
				weight *= weightBottom;
			}

			double diff = ( base_array[ i + sliver->width * j ] + r_change ) - ( sliver_array[ i + sliver->width * j ] + c_change );
			double y;
			double t;
			cs.lock();
			y = diff * diff * weight - sum_error;
			t = sum + y;
			sum_error = ( t - sum ) - y;
			sum = t;

			y = weight - area_error;
			t = area + y;
			area_error = ( t - area ) - y;
			area = t;
			cs.unlock();

			y = diff * diff * weight - solo_sum_error;
			t = solo_sum + y;
			solo_sum_error = ( t - solo_sum ) - y;
			solo_sum = t;

			y = weight - solo_weight_error;
			t = solo_weight + y;
			solo_weight_error = ( t - solo_weight ) - y;
			solo_weight = t;
		}

		old_sums[ i ] = solo_sum;
		old_weights[ i ] = solo_weight;

	} );
#endif


	// WRITE TO FILE FOR TBB (CPU THREADED)
#ifdef __CODEFOLD__

	FILE *old_file = fopen( "old_sliver.txt", "w" );
	if ( old_file == NULL ) {
		printf( "File Not opened!\n" );
	}
	else {
		for (int x = 0; x < maxX + 1; x++) {
			for (int y = minY; y < maxY + 1; y++) {
				fprintf( old_file, "(%d, %d) = %.30e\n", x, y, sliver_array[ x + sliver->width * y ] );
			}
		}
	}
	fclose( old_file );


	old_file = fopen( "old_sw.txt", "w" );
	if ( old_file == NULL ) {
		printf( "File Not opened!\n" );
	}
	else {
		fprintf( old_file, "Calculated Sums Values : \n" );
		for ( int i = 0; i < maxX + 1; i++ ) {
			fprintf( old_file, "%.30e\n", old_sums[ i ] );
		}
		fprintf( old_file, "\n\n\n" );
		fprintf( old_file, "Calculated Weights Values : \n" );
		for ( int i = 0; i < maxX + 1; i++ ) {
			fprintf( old_file, "%.30e\n", old_weights[ i ] );
		}
	}
	fclose( old_file );


	old_file = fopen ( "old_base.txt", "w" );
	if ( old_file == NULL ) {
		printf ( "File Not opened!\n" );
	}
	else {
		for ( int x = 0; x < maxX + 1; x++ ) {
			for ( int y = minY; y < maxY + 1; y++ ) {
				fprintf ( old_file, "(%d, %d) = %.30e\n", x, y, base_array[ x + sliver->width * y ] );
			}
		}
	}
	fclose ( old_file );

	
	old_file = fopen( "old_info.txt", "w" );
	if ( old_file == NULL ) {
		printf( "File Not opened!\n" );
	}
	else {
		fprintf( old_file, "MinY : %d, MaxY : %d, MaxX : %d\n", minY, maxY, maxX );
		fprintf( old_file, "Fitness : %.30e\n", sum / area );
	}
	fclose( old_file );
	

	old_file = fopen( "old_matrix.txt", "w" );
	if ( old_file != NULL ) {
		fprintf( old_file, "This is the generated C vector : \n\n" );
		for ( int i = 0; i < m - 1; ++i ) {
			fprintf( old_file, "%.30e\n", KRCs[ i ] );
		}
		fprintf( old_file, "This is the generated R vector : \n\n" );
		for ( int j = m - 1; j < size; ++j ) {
			fprintf( old_file, "%.30e\n", KRCs[ j ] );
		}
		fprintf( old_file, "\n\n\n" );
		fprintf( old_file, "This is the computed C vector : \n\n" );
		for ( int i = 0; i < m - 1; ++i ) {
			fprintf( old_file, "%.30e\n", RCs[ i ] );
		}
		fprintf( old_file, "This is the computed R vector : \n\n" );
		for ( int j = m - 1; j < size; ++j ) {
			fprintf( old_file, "%.30e\n", RCs[ j ] );
		}
	}
	fclose( old_file );
	
#endif

	// OPENCL INITIALIZATION
#ifdef __CODEFOLD__
	// Get all the platforms available
	cl::Program::Sources sources;
	std::vector<cl::Platform> all_platforms;
	cl::Platform::get(&all_platforms);
	if (all_platforms.size() == 0) {
		printf("No platforms found. Check OpenCL installation!\n");
		getchar();
		exit(1);
	}
	else {
		printf("Found %d platforms.\n", all_platforms.size());
	}
	printf("\n");
	cl::Platform default_platform = all_platforms[0];
	printf("Using platform %s\n", default_platform.getInfo<CL_PLATFORM_NAME>().c_str());
	printf("\n");

	// Get default device of the default platform.
	std::vector<cl::Device> all_devices;
	default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
	if (all_devices.size() == 0) {
		printf("No devices found. Check OpenCL installation!\n");
		getchar();
		exit(1);
	}
	else {
		printf("Found %d devices.\n", all_devices.size());
	}

	// Print device list.
	for (int i = 0; i < all_devices.size(); ++i) {
		printf("Device %d name : %s\n", i, all_devices[i].getInfo<CL_DEVICE_NAME>().c_str());
	}

	printf("\n");
	cl::Device default_device = all_devices[ 0 ];
	printf("Using Primary Device : %s\n", default_device.getInfo<CL_DEVICE_NAME>().c_str());
	printf("\n");
	cl::Context context({ default_device });


	std::vector<size_t> max_work_group_size = default_device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
	printf( "Maximum work item sizes : ( %d, %d, %d ).\n", max_work_group_size[ 0 ], max_work_group_size[ 1 ], max_work_group_size[ 2 ] );
	// max_work_group_size = 1;

	int local_mem_size = default_device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
	printf("Maximum local memory size : %d kBs.\n", local_mem_size / 1024);

	unsigned long global_mem_size = default_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
	printf("Maximum global memory size : %I64u MBs.\n", global_mem_size / 1024 / 1024);

	int number_compute_units = default_device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS> ();
	printf ( "Number of compute units : %d.\n", number_compute_units );
	
	sources.push_back( { bicubic_kernels_string.c_str(), bicubic_kernels_string.length() } );
	sources.push_back( { summation_kernel_string.c_str(), summation_kernel_string.length() } );
	sources.push_back( { bspline_kernels_string.c_str(), bspline_kernels_string.length() } );

	cl::Program program( context, sources );
	if ( program.build( { default_device } ) != CL_SUCCESS ) {
		printf( "Error building: %s\n", program.getBuildInfo<CL_PROGRAM_BUILD_LOG>( default_device ).c_str() );
		getchar();
		exit( 1 );
	}
	else {
		printf( "Kernel building success!\n" );
	}

#endif


	// BUFFER CREATION
#ifdef __CODEFOLD__
	// create buffers on the device
	cl_int* error_code;
	cl::Buffer cl_sliver( context, CL_MEM_READ_ONLY, sizeof( double ) * sliver->width * sliver->height );
	cl::Buffer cl_sliver_interp( context, CL_MEM_READ_WRITE, sizeof( double ) * sliver->width * sliver->height );
	cl::Buffer cl_base ( context, CL_MEM_READ_ONLY, sizeof ( double ) * base->width * base->height );
	cl::Buffer cl_base_interp ( context, CL_MEM_READ_WRITE, sizeof ( double ) * sliver->width * sliver->height );
	cl::Buffer cl_krcs( context, CL_MEM_READ_WRITE, sizeof( double ) * ( sliver->width + sliver->height + 1 ) );
	cl::Buffer cl_rcs( context, CL_MEM_READ_WRITE, sizeof( double ) * ( sliver->width + sliver->height ) );

	cl::Buffer cl_partial_sums( context, CL_MEM_READ_WRITE, sizeof( double ) * sliver->height );
	cl::Buffer cl_partial_weights( context, CL_MEM_READ_WRITE, sizeof( double ) * sliver->height );
	cl::Buffer cl_partial_sums_errors( context, CL_MEM_READ_WRITE, sizeof( double ) * sliver->height );
	cl::Buffer cl_partial_weights_errors( context, CL_MEM_READ_WRITE, sizeof( double ) * sliver->height );

	cl::Buffer cl_final_result( context, CL_MEM_READ_WRITE, sizeof( double ) );

	//create queue to which we will push commands for the device.
	cl::CommandQueue queue( context, default_device );

	queue.enqueueWriteBuffer( cl_sliver, CL_BLOCKING, 0, sizeof( double ) * sliver->width * sliver->height, double_sliv );
	queue.enqueueWriteBuffer( cl_base, CL_BLOCKING, 0, sizeof( double ) * base->width * base->height, double_base );
#endif


	// KERNEL CREATION
#ifdef __CODEFOLD__

	cl::Kernel sliverKernel = cl::Kernel( program, "clWarpSliverCubic" );
	cl::Kernel baseKernel = cl::Kernel( program, "clWarpBaseCubic" );

	cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double,
		cl_int, cl_int, cl_int, cl_int, cl_int >
		clWarpSliverCubic( cl::Kernel( program, "clWarpSliverCubic" ) );

	cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double,
		cl_int, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int >
		clWarpBaseCubic( cl::Kernel( program, "clWarpBaseCubic" ) );

	cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double,
		cl_int, cl_int, cl_int, cl_int, cl_int, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg& >
		clWarpSliverBspline( cl::Kernel( program, "clWarpSliverBspline" ) );

	cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double,
		cl_int, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg& >
		clWarpBaseBspline( cl::Kernel( program, "clWarpBaseBspline" ) );

	cl::make_kernel< cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int >
		clColumnSums( cl::Kernel( program, "clColumnSums" ) );
	cl::make_kernel< cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int >
		clRowSums( cl::Kernel( program, "clRowSums" ) );
	cl::make_kernel< cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl_int, cl_int, cl_int > clRCResults( cl::Kernel( program, "clRCResults" ) );


	cl::make_kernel < cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
		cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&,
		cl_int, cl_int, cl_int, cl_int,
		cl_double, cl_double, cl_double > clPartialDifference( cl::Kernel( program, "clPartialDifference" ) );

	cl::make_kernel < cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
		cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&,
		cl_int, cl_int, cl_int > clFinalDifference( cl::Kernel( program, "clFinalDifference" ) );


	
	cl::make_kernel < cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
		cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&,
		cl_int, cl_int, cl_int, cl_int,
		cl_double, cl_double, cl_double >
		clDifferenceCalculation( cl::Kernel( program, "clDifferenceCalculation" ) );
	

#endif

	// KERNEL INFO
#ifdef __CODEFOLD__ 
	size_t kernel_size;
	sliverKernel.getWorkGroupInfo<size_t>( default_device, CL_KERNEL_WORK_GROUP_SIZE, &kernel_size );
	printf( "Kernel size : %d.\n", kernel_size );

	std::vector<size_t> kernel_compile_size;
	sliverKernel.getWorkGroupInfo<std::vector<size_t>>( default_device, CL_KERNEL_COMPILE_WORK_GROUP_SIZE, &kernel_compile_size );
	printf( "Kernel compile size : ( %d, %d, %d ).\n", kernel_compile_size[ 0 ], kernel_compile_size[ 1 ], kernel_compile_size[ 2 ] );

	size_t kernel_preferred_size;
	sliverKernel.getWorkGroupInfo<size_t>( default_device, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, &kernel_preferred_size );
	printf( "Kernel preferred size : %d.\n", kernel_preferred_size );

	cl_ulong kernel_local_memory;
	sliverKernel.getWorkGroupInfo<cl_ulong>( default_device, CL_KERNEL_LOCAL_MEM_SIZE, &kernel_local_memory );
	printf( "Kernel local size : %I64u.\n", kernel_local_memory );
#endif

	// KERNEL INVOCATION
#ifdef __CODEFOLD__
	// Find largest power of 2 the dimensions fit into
	int dimX = maxX + 1;
	int dimensionX = nextPow2( dimX );
	int dimY = maxY - minY + 1;
	int dimensionY = nextPow2( dimY );
	int ratio1 = dimensionY > 64 ? ceil( dimensionY / ( double ) 64 ) : 1;
	int ratio2 = dimensionX > 64 ? ceil( dimensionX / ( double ) 128 ) : 1;
	printf( "Image Dimensions Raised to Power of 2 : ( %d, %d )\n\n", dimensionX, dimensionY );
	printf( "Kernel ratios : %d, %d \n \n", ratio1, ratio2 );
	printf( "Sliver Kernel Dimensions : ( %d, %d ) \n \n", dimX, dimensionY / ratio1 );
	printf( "Base Kernel Dimensions : ( %d, %d ) \n \n", dimensionX / ratio2, dimY );

	cl::LocalSpaceArg locals = cl::Local( sizeof( double ) * 4 );
	cl::EnqueueArgs range_sliver_args( queue, cl::NullRange, cl::NDRange( ( int ) ceil( dimX / 1.0 ) * dimensionY / ratio1 ), cl::NDRange( dimensionY / ratio1 ) );
	cl::EnqueueArgs range_base_args( queue, cl::NullRange, cl::NDRange( ( int ) ceil( dimY / 1.0 ) * dimensionX / ratio2 ), cl::NDRange( dimensionX / ratio2 ) );

	int num_sub_blocks1 = ceil( sqrt( n ) * 2 );
	// num_sub_blocks1 = nextPow2( num_sub_blocks1 );
	int num_blocks1 = 128 / num_sub_blocks1;
	int num_sub_blocks2 = ceil( sqrt( m ) * 2 );
	// num_sub_blocks2 = nextPow2( num_sub_blocks2 );
	int num_blocks2 = 256 / num_sub_blocks2;
	if ( default_device.getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU ) {
		num_sub_blocks1 = n;
		num_blocks1 = 1;
		num_sub_blocks2 = m;
		num_blocks2 = 1;
	}
	cl::EnqueueArgs range_column_sum_args( queue, cl::NullRange, cl::NDRange( m * num_sub_blocks1 ), cl::NDRange( num_sub_blocks1 * num_blocks1 ) );
	cl::EnqueueArgs range_row_sum_args( queue, cl::NullRange, cl::NDRange( n * num_sub_blocks2 ), cl::NDRange( num_sub_blocks2 * num_blocks2 ) );

	

	int size_global = nextPow2( n );
	int row_local = nextPow2( sqrt( m ) * 2 );
	cl::EnqueueArgs range_results_args( queue, cl::NullRange, cl::NDRange( size_global ), cl::NDRange( row_local ) );
	
	int partial_range = sqrt( n );
	cl::EnqueueArgs range_partial_sum_args( queue, cl::NullRange, cl::NDRange( m * partial_range ), cl::NDRange( partial_range ) );

	int final_range = sqrt( m );
	cl::EnqueueArgs range_final_results_args( queue, cl::NullRange, cl::NDRange( final_range ), cl::NDRange( final_range ) );

	int final_size = 256;
	// cl::EnqueueArgs range_final_args( queue, cl::NullRange, cl::NDRange( final_size ), cl::NDRange( final_size ) );

	clWarpSliverCubic( range_sliver_args, cl_sliver, cl_sliver_interp,
					   ( double ) As[ 1 ], ( double ) Bs[ 1 ],
					   ( int ) sliver->width, ( int ) sliver->height,
					   maxX, minY, maxY ).wait();
	clWarpBaseCubic( range_base_args, cl_base, cl_base_interp,
					 ( double ) As[ 0 ], ( double ) As[ 1 ], ( double ) As[ 2 ], ( double ) As[ 3 ],
					 ( double ) Bs[ 0 ], ( double ) Bs[ 1 ], ( double ) Bs[ 2 ], ( double ) Bs[ 3 ],
					 ( int ) sliver->width, ( int ) sliver->height, ( int ) base->width,
					 xbp, maxX, minY, maxY ).wait();
	clColumnSums( range_column_sum_args, cl_sliver_interp, cl_base_interp, cl_krcs,
				  cl::Local( sizeof( double ) * num_blocks1 * num_sub_blocks1 ), cl::Local( sizeof( double ) * num_blocks1 * num_sub_blocks1 ),
				  ( int ) sliver->width, minY, m, n, num_sub_blocks1, num_blocks1 ).wait();
	clRowSums( range_row_sum_args, cl_sliver_interp, cl_base_interp, cl_krcs,
			   cl::Local( sizeof( double ) * num_blocks2 * num_sub_blocks2 ), cl::Local( sizeof( double ) * num_blocks2 * num_sub_blocks2 ),
			   ( int ) sliver->width, minY, m, n, num_sub_blocks2, num_blocks2 ).wait();

	clRCResults( range_results_args, cl_krcs, cl_rcs, cl::Local( sizeof( double ) * row_local ), cl::Local( sizeof( double ) * row_local ), ( int ) sliver->width, m, n ).wait();

	printf("Starting OpenCL calls!\n");
	// start = clock();
	for (int i = 0; i < ITERATIONS; ++i) {
		clWarpSliverBspline( range_sliver_args, cl_sliver, cl_sliver_interp,
						   ( double ) As[ 1 ], ( double ) Bs[ 1 ],
						   ( int ) sliver->width, ( int ) sliver->height,
						   maxX, minY, maxY, locals, locals, locals, locals ).wait();
		clWarpBaseBspline( range_base_args, cl_base, cl_base_interp,
						 ( double ) As[ 0 ], ( double ) As[ 1 ], ( double ) As[ 2 ], ( double ) As[ 3 ],
						 ( double ) Bs[ 0 ], ( double ) Bs[ 1 ], ( double ) Bs[ 2 ], ( double ) Bs[ 3 ],
						 ( int ) sliver->width, ( int ) sliver->height, ( int ) base->width,
						 xbp, maxX, minY, maxY, locals, locals, locals, locals ).wait();
	}
	// end = clock();

	
	
	clPartialDifference( range_partial_sum_args, cl_partial_sums, cl_partial_weights, cl_partial_sums_errors, cl_partial_weights_errors, cl_rcs, cl_base_interp, cl_sliver_interp,
						 cl::Local( sizeof( double ) * partial_range ), cl::Local( sizeof( double ) * partial_range ),
						 cl::Local( sizeof( double ) * partial_range ), cl::Local( sizeof( double ) * partial_range ),
						 maxX, minY, maxY, ( int ) sliver->width, weightRight, weightTop, weightBottom ).wait();

	
	

	clFinalDifference( range_final_results_args, cl_partial_sums, cl_partial_weights, cl_partial_sums_errors, cl_partial_weights_errors, cl_final_result,
					   cl::Local( sizeof( double ) * row_local ), cl::Local( sizeof( double ) * row_local ),
					   cl::Local( sizeof( double ) * row_local ), cl::Local( sizeof( double ) * row_local ),
					   maxX, minY, maxY ).wait();
	

	/*

	clDifferenceCalculation( range_final_args, cl_final_result, cl_rcs,
							 cl_base_interp, cl_sliver_interp,
							 cl::Local( sizeof( double ) * final_size ), cl::Local( sizeof( double ) * final_size ),
							 cl::Local( sizeof( double ) * final_size ), cl::Local( sizeof( double ) * final_size ),
							 maxX, minY, maxY, ( int ) sliver->width, weightRight, weightTop, weightBottom ).wait();
	*/
							 

	// printf ( "Total running time for OpenCL base and sliver warping at %d iterations : %f\n", ITERATIONS, ( ( double ) end - ( double ) start ) / CLOCKS_PER_SEC );

#endif


	// WRITE TO FILE FOR OPENCL
#ifdef __CODEFOLD__

	double* sliv_array = new double[ sliver->width * sliver->height ];
	double* b_array = new double[ sliver->width * sliver->height ];
	double* sums_array = new double[ sliver->width + sliver->height + 1 ];
	double* results_array = new double[ sliver->width + sliver->height ];
	double* fitness = new double[ 1 ];

	double* sums = new double[ sliver->height ];
	double* weights = new double[ sliver->height ];

	queue.enqueueReadBuffer( cl_sliver_interp, CL_TRUE, 0, sizeof( double ) * sliver->width * sliver->height, sliv_array );
	queue.enqueueReadBuffer ( cl_base_interp, CL_TRUE, 0, sizeof ( double ) * sliver->width * sliver->height, b_array );
	queue.enqueueReadBuffer( cl_krcs, CL_TRUE, 0, sizeof( double ) * ( sliver->width + sliver->height + 1 ), sums_array );
	queue.enqueueReadBuffer( cl_rcs, CL_TRUE, 0, sizeof( double ) * ( sliver->width + sliver->height ), results_array );
	queue.enqueueReadBuffer( cl_final_result, CL_TRUE, 0, sizeof( double ), fitness );

	queue.enqueueReadBuffer( cl_partial_sums, CL_TRUE, 0, sizeof( double ) * ( sliver->width ), sums );
	queue.enqueueReadBuffer( cl_partial_weights, CL_TRUE, 0, sizeof( double ) * ( sliver->width ), weights );


	FILE *new_file = fopen( "new_sliver.txt", "w" );
	if ( new_file == NULL ) {
		printf( "File Not opened!\n" );
	}
	else {
		for ( int x = 0; x < maxX + 1; x++ ) {
			for ( int y = minY; y < maxY + 1; y++ ) {
				fprintf( new_file, "(%d, %d) = %.30e\n", x, y, sliv_array[ x + sliver->width * y ] );
			}
		}
	}
	fclose( new_file );

	new_file = fopen( "new_sw.txt", "w" );
	if ( new_file == NULL ) {
		printf( "File Not opened!\n" );
	}
	else {
		fprintf( new_file, "Calculated Sums Values : \n" );
		for ( int i = 0; i < maxX + 1; i++ ) {
			fprintf( new_file, "%.30e\n", sums[ i ] );
		}
		fprintf( new_file, "\n\n\n" );
		fprintf( new_file, "Calculated Weights Values : \n" );
		for ( int i = 0; i < maxX + 1; i++ ) {
			fprintf( new_file, "%.30e\n", weights[ i ] );
		}
	}
	fclose( new_file );

	
	new_file = fopen ( "new_base.txt", "w" );
	if ( new_file == NULL ) {
		printf ( "File Not opened!\n" );
	}
	else {
		for ( int x = 0; x < maxX + 1; x++ ) {
			for ( int y = minY; y < maxY + 1; y++ ) {
				fprintf ( new_file, "(%d, %d) = %.30e\n", x, y, b_array[ x + sliver->width * y ] );
			}
		}
	}
	fclose ( new_file );

	
	new_file = fopen( "new_info.txt", "w" );
	if ( old_file == NULL ) {
		printf( "File Not opened!\n" );
	}
	else {
		fprintf( new_file, "MinY : %d, MaxY : %d, MaxX : %d\n", minY, maxY, maxX );
		fprintf( new_file, "Fitness : %.30e\n", fitness[ 0 ] );
		// fprintf( new_file, "Sums : %f\n", fitness[ 1 ] );
		// fprintf( new_file, "Weights : %f\n", fitness[ 2 ] );
	}
	fclose( new_file );
	

	new_file = fopen( "new_matrix.txt", "w" );
	if ( old_file != NULL ) {
		fprintf( new_file, "This is the generated C vector : \n\n" );
		for ( int i = 0; i < m - 1; ++i ) {
			fprintf( new_file, "%.30e\n", sums_array[ i ] );
		}
		fprintf( new_file, "This is the generated R vector : \n\n" );
		for ( int j = 0; j < n; ++j ) {
			fprintf( new_file, "%.30e\n", sums_array[ sliver->width + 1 + j ] );
		}
		fprintf( new_file, "\n\n\n" );
		fprintf( new_file, "This is the computed C vector : \n\n" );
		for ( int i = 0; i < m - 1; ++i ) {
			fprintf( new_file, "%.30e\n", results_array[ i ] );
		}
		fprintf( new_file, "This is the computed R vector : \n\n" );
		for ( int j = 0; j < n; ++j ) {
			fprintf( new_file, "%.30e\n", results_array[ sliver->width + j ] );
		}
	}
	fclose( new_file );

	printf( "Finished!\n" );

	delete[] As;
	delete[] Bs;
	delete[] Cs;
	delete base;
	delete sliver;
	delete[] sliver_array;
	delete[] base_array;
	delete[] double_sliv;
	delete[] double_base;
	delete[] sliv_array;
	delete[] b_array;
	getchar();

#endif
}



#endif

