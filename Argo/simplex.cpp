/*Here is an implementation in C.  Using it on Rosenbrock's famous function

determineDifference(x,y)   =   (1-x)^2  +  100[ (x^2 - y)^2 ]

from the traditional starting point   (x,y) = (-1.2, +1.0)
this code calculated a minimum at     (x,y) = ( 1.000026  ,  1.000051 )
while using 177 evaluations of the function determineDifference.
(It's easy to see that Rosenbrock's function determineDifference has a minimum at (1,1)).

No warranty is expressed or implied, use at your own risk, remember that
this software is worth what you paid for it, and you paid zero.

=============================================================================
=============================================================================
*/

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <vector>
#include "Eigen\Eigen"
#include "simplex.h"


#define         PRINTEM         (100)

//#define         ALP             (1.5)   /* reflection parameter   */
//#define         BET             (0.7)   /* contraction parameter  */
//#define         GAM             (2.1)   /* expansion parameter    */
//#define         TINY            (1e-5)  /* first Simplex displace */

bool fastZ = true;

long difference_evals = 0;
bool debug = false;
int operations = 0;

FImage *base;
FImage *sliver;
int *basesize;
int *sliversize;


// Declare all necessary OpenCL devices, contexts, queues, buffers, kernels, etc.
#ifdef S_MODE_OPENCL
unsigned int nextPow2( unsigned int x ) {
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

cl::Program::Sources sources;
cl::Platform default_platform;
cl::Device default_device;
cl::Context context;
cl::Program program;
cl::Buffer cl_sliver;
cl::Buffer cl_sliver_interp;
cl::Buffer cl_base;
cl::Buffer cl_base_interp;
cl::Buffer cl_krcs;
cl::Buffer cl_rcs;

cl::Buffer cl_partial_sums;
cl::Buffer cl_partial_weights;
cl::Buffer cl_partial_sums_errors;
cl::Buffer cl_partial_weights_errors;

cl::Buffer cl_final_result;
cl::CommandQueue queue;

cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double,
	cl_int, cl_int, cl_int, cl_int, cl_int > *
	clWarpSliverCubic;

cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double,
	cl_int, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int > *
	clWarpBaseCubic;

cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double,
	cl_int, cl_int, cl_int, cl_int, cl_int, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg& > *
	clWarpSliverBspline;

cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double,
	cl_int, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg& > *
	clWarpBaseBspline;

cl::make_kernel< cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int > * clColumnSums;
cl::make_kernel< cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int > * clRowSums;
cl::make_kernel< cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl_int, cl_int, cl_int > * clRCResults;

cl::make_kernel < cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
	cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&,
	cl_int, cl_int, cl_int, cl_int,
	cl_double, cl_double, cl_double > * clPartialDifference;

cl::make_kernel < cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
	cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&,
	cl_int, cl_int, cl_int > * clFinalDifference;

/*
cl::make_kernel < cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
	cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&,
	cl_int, cl_int, cl_int, cl_int,
	cl_double, cl_double, cl_double > *
	clDifferenceCalculation;
*/

cl::LocalSpaceArg locals;
double* difference_value;


/**
 *
 *
 */
void initializeOpenCLContext() {
	double* double_sliv = new double[ sliver->width * sliver->height ];
	double* double_base = new double[ base->width * base->height ];
	for ( int i = 0; i < sliver->width * sliver->height; ++i ) {
		double_sliv[ i ] = sliver->data[ i ];
	}
	for ( int i = 0; i < base->width * base->height; ++i ) {

		double_base[ i ] = base->data[ i ];
		// fprintf( image_file, "%f\n", base->data[ i ] );
		// fprintf( image_array_file, "%f\n", double_base[ i ] );
	}

	// Get all the platforms available
	
	std::vector<cl::Platform> all_platforms;
	cl::Platform::get( &all_platforms );
	if ( all_platforms.size() == 0 ) {
		printf( "No platforms found. Check OpenCL installation!\n" );
		getchar();
		exit( 1 );
	}
	else {
		printf( "Found %d platforms.\n", all_platforms.size() );
	}
	printf( "\n" );
	default_platform = all_platforms[ 0 ];
	printf( "Using platform %s\n", default_platform.getInfo<CL_PLATFORM_NAME>().c_str() );
	printf( "\n" );

	// Get default device of the default platform.
	std::vector<cl::Device> all_devices;
	default_platform.getDevices( CL_DEVICE_TYPE_ALL, &all_devices );
	if ( all_devices.size() == 0 ) {
		printf( "No devices found. Check OpenCL installation!\n" );
		getchar();
		exit( 1 );
	}
	else {
		printf( "Found %d devices.\n", all_devices.size() );
	}

	// Print device list.
	for ( int i = 0; i < all_devices.size(); ++i ) {
		printf( "Device %d name : %s\n", i, all_devices[ i ].getInfo<CL_DEVICE_NAME>().c_str() );
	}

	printf( "\n" );
	default_device = all_devices[ 1 ];
	printf( "Using Primary Device : %s\n", default_device.getInfo<CL_DEVICE_NAME>().c_str() );
	printf( "\n" );
	context = cl::Context( { default_device } );


	std::vector<size_t> max_work_group_size = default_device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
	printf( "Maximum work item sizes : ( %d, %d, %d ).\n", max_work_group_size[ 0 ], max_work_group_size[ 1 ], max_work_group_size[ 2 ] );
	// max_work_group_size = 1;

	int local_mem_size = default_device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
	printf( "Maximum local memory size : %d kBs.\n", local_mem_size / 1024 );

	unsigned long global_mem_size = default_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
	printf( "Maximum global memory size : %I64u MBs.\n", global_mem_size / 1024 / 1024 );

	int number_compute_units = default_device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
	printf( "Number of compute units : %d.\n", number_compute_units );

	sources.push_back( { bicubic_kernels_string.c_str(), bicubic_kernels_string.length() } );
	sources.push_back( { summation_kernel_string.c_str(), summation_kernel_string.length() } );
	sources.push_back( { bspline_kernels_string.c_str(), bspline_kernels_string.length() } );

	program = cl::Program( context, sources );
	if ( program.build( { default_device } ) != CL_SUCCESS ) {
		printf( "Error building: %s\n", program.getBuildInfo<CL_PROGRAM_BUILD_LOG>( default_device ).c_str() );
		getchar();
		exit( 1 );
	}
	else {
		printf( "Kernel building success!\n" );
	}

	// create buffers on the device
	cl_int* error_code;
	cl_sliver = cl::Buffer( context, CL_MEM_READ_ONLY, sizeof( double ) * sliver->width * sliver->height );
	cl_sliver_interp = cl::Buffer( context, CL_MEM_READ_WRITE, sizeof( double ) * sliver->width * sliver->height );
	cl_base = cl::Buffer( context, CL_MEM_READ_ONLY, sizeof( double ) * base->width * base->height );
	cl_base_interp = cl::Buffer( context, CL_MEM_READ_WRITE, sizeof( double ) * sliver->width * sliver->height );
	cl_krcs = cl::Buffer( context, CL_MEM_READ_WRITE, sizeof( double ) * ( sliver->width + sliver->height + 1 ) );
	cl_rcs = cl::Buffer( context, CL_MEM_READ_WRITE, sizeof( double ) * ( sliver->width + sliver->height ) );

	cl_partial_sums = cl::Buffer( context, CL_MEM_READ_WRITE, sizeof( double ) * sliver->height );
	cl_partial_weights = cl::Buffer( context, CL_MEM_READ_WRITE, sizeof( double ) * sliver->height );
	cl_partial_sums_errors = cl::Buffer( context, CL_MEM_READ_WRITE, sizeof( double ) * sliver->height );
	cl_partial_weights_errors = cl::Buffer( context, CL_MEM_READ_WRITE, sizeof( double ) * sliver->height );

	cl_final_result = cl::Buffer( context, CL_MEM_READ_WRITE, sizeof( double ) );

	//create queue to which we will push commands for the device.
	queue = cl::CommandQueue( context, default_device );

	queue.enqueueWriteBuffer( cl_sliver, CL_BLOCKING, 0, sizeof( double ) * sliver->width * sliver->height, double_sliv );
	queue.enqueueWriteBuffer( cl_base, CL_BLOCKING, 0, sizeof( double ) * base->width * base->height, double_base );

	clWarpSliverCubic = new cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double,
		cl_int, cl_int, cl_int, cl_int, cl_int >( cl::Kernel( program, "clWarpSliverCubic" ) );

	clWarpBaseCubic = new cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double,
		cl_int, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int >
		( cl::Kernel( program, "clWarpBaseCubic" ) );

	clWarpSliverBspline = new cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double,
		cl_int, cl_int, cl_int, cl_int, cl_int, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg& >
		( cl::Kernel( program, "clWarpSliverBspline" ) );

	clWarpBaseBspline = new cl::make_kernel< cl::Buffer&, cl::Buffer&, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double, cl_double,
		cl_int, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg& >
		( cl::Kernel( program, "clWarpBaseBspline" ) );

	clColumnSums = new cl::make_kernel< cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int >
		( cl::Kernel( program, "clColumnSums" ) );
	clRowSums = new cl::make_kernel< cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl_int, cl_int, cl_int, cl_int, cl_int, cl_int >
		( cl::Kernel( program, "clRowSums" ) );
	clRCResults = new cl::make_kernel< cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl_int, cl_int, cl_int >( cl::Kernel( program, "clRCResults" ) );

	clPartialDifference = new cl::make_kernel < cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
		cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&,
		cl_int, cl_int, cl_int, cl_int,
		cl_double, cl_double, cl_double >( cl::Kernel( program, "clPartialDifference" ) );

	clFinalDifference = new cl::make_kernel < cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
		cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&,
		cl_int, cl_int, cl_int >( cl::Kernel( program, "clFinalDifference" ) );

	/*

	clDifferenceCalculation = new cl::make_kernel < cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
		cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&, cl::LocalSpaceArg&,
		cl_int, cl_int, cl_int, cl_int,
		cl_double, cl_double, cl_double >
		( cl::Kernel( program, "clDifferenceCalculation" ) );
	*/

	locals = cl::Local( sizeof( double ) * 4 );
	difference_value = new double[ 1 ];
}
#else
	double* base_array;
	double* sliver_array;
	double* KRCs;
	double* RCs;
	#ifdef S_MODE_TBB
		tbb::task_scheduler_init init;
		tbb::critical_section cs;
	#endif
#endif

/**
 *	Performs fastZ correction.
 *
 *	@param	vertex	An array containing parameter values to apply to the image.
 *	@param	writeFile	A flag indicating whether to output data to an image. Useful for the last iteration.
 *	
 *	@return	A double precision difference value of the image.
 */
double fastZDifference( double vertex[], bool writeFile ) {
	// Precomputation of bounds and weights necessary.
	double As[ 4 ] = { vertex[ 0 ], vertex[ 2 ], vertex[ 4 ], vertex[ 6 ] };
	double Bs[ 4 ] = { vertex[ 1 ], vertex[ 3 ], vertex[ 5 ], vertex[ 7 ] };

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

	// OpenCL mode of FastZ correction.
#ifdef S_MODE_OPENCL

	// Initialize NDRange arguments.
	// Find largest power of 2 the dimensions fit into
	int dimX = maxX + 1;
	int dimY = maxY - minY + 1;
	int size1 = 64;
	int size2 = 64;

	cl::EnqueueArgs range_sliver_args( queue, cl::NullRange, cl::NDRange( ( int ) ceil( dimX / 1.0 ) * size1 ), cl::NDRange( size1 ) );
	cl::EnqueueArgs range_base_args( queue, cl::NullRange, cl::NDRange( ( int ) ceil( dimY / 1.0 ) * size2 ), cl::NDRange( size2 ) );

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

	// int final_size = 256;
	// cl::EnqueueArgs range_final_args( queue, cl::NullRange, cl::NDRange( final_size ), cl::NDRange( final_size ) );


	// Kernel calls.
	clWarpSliverCubic->operator()( range_sliver_args, cl_sliver, cl_sliver_interp,
					   ( double ) As[ 1 ], ( double ) Bs[ 1 ],
					   ( int ) sliver->width, ( int ) sliver->height,
					   maxX, minY, maxY ).wait();
	clWarpBaseCubic->operator()( range_base_args, cl_base, cl_base_interp,
					 ( double ) As[ 0 ], ( double ) As[ 1 ], ( double ) As[ 2 ], ( double ) As[ 3 ],
					 ( double ) Bs[ 0 ], ( double ) Bs[ 1 ], ( double ) Bs[ 2 ], ( double ) Bs[ 3 ],
					 ( int ) sliver->width, ( int ) sliver->height, ( int ) base->width,
					 xbp, maxX, minY, maxY ).wait();
	clColumnSums->operator()( range_column_sum_args, cl_sliver_interp, cl_base_interp, cl_krcs,
				  cl::Local( sizeof( double ) * num_blocks1 * num_sub_blocks1 ),
				  cl::Local( sizeof( double ) * num_blocks1 * num_sub_blocks1 ),
				  ( int ) sliver->width, minY, m, n, num_sub_blocks1, num_blocks1 ).wait();
	clRowSums->operator()( range_row_sum_args, cl_sliver_interp, cl_base_interp, cl_krcs,
			   cl::Local( sizeof( double ) * num_blocks2 * num_sub_blocks2 ),
			   cl::Local( sizeof( double ) * num_blocks2 * num_sub_blocks2 ),
			   ( int ) sliver->width, minY, m, n, num_sub_blocks2, num_blocks2 ).wait();

	clRCResults->operator()( range_results_args, cl_krcs, cl_rcs,
				cl::Local( sizeof( double ) * row_local ), cl::Local( sizeof( double ) * row_local ),
				( int ) sliver->width, m, n ).wait();

	clWarpSliverBspline->operator()( range_sliver_args, cl_sliver, cl_sliver_interp,
							( double ) As[ 1 ], ( double ) Bs[ 1 ],
							( int ) sliver->width, ( int ) sliver->height,
							maxX, minY, maxY, locals, locals, locals, locals ).wait();
	clWarpBaseBspline->operator()( range_base_args, cl_base, cl_base_interp,
						( double ) As[ 0 ], ( double ) As[ 1 ], ( double ) As[ 2 ], ( double ) As[ 3 ],
						( double ) Bs[ 0 ], ( double ) Bs[ 1 ], ( double ) Bs[ 2 ], ( double ) Bs[ 3 ],
						( int ) sliver->width, ( int ) sliver->height, ( int ) base->width,
						xbp, maxX, minY, maxY, locals, locals, locals, locals ).wait();

	clPartialDifference->operator()( range_partial_sum_args, cl_partial_sums, cl_partial_weights, cl_partial_sums_errors, cl_partial_weights_errors, cl_rcs, cl_base_interp, cl_sliver_interp,
						 cl::Local( sizeof( double ) * partial_range ), cl::Local( sizeof( double ) * partial_range ),
						 cl::Local( sizeof( double ) * partial_range ), cl::Local( sizeof( double ) * partial_range ),
						 maxX, minY, maxY, ( int ) sliver->width, weightRight, weightTop, weightBottom ).wait();

	clFinalDifference->operator()( range_final_results_args, cl_partial_sums, cl_partial_weights, cl_partial_sums_errors, cl_partial_weights_errors, cl_final_result,
						cl::Local( sizeof( double ) * final_range ), cl::Local( sizeof( double ) * final_range ),
						cl::Local( sizeof( double ) * final_range ), cl::Local( sizeof( double ) * final_range ),
						maxX, minY, maxY ).wait();

	/*
	clDifferenceCalculation->operator()( range_final_args, cl_final_result, cl_rcs,
							 cl_base_interp, cl_sliver_interp,
							 cl::Local( sizeof( double ) * final_size ), cl::Local( sizeof( double ) * final_size ),
							 cl::Local( sizeof( double ) * final_size ), cl::Local( sizeof( double ) * final_size ),
							 maxX, minY, maxY, ( int ) sliver->width, weightRight, weightTop, weightBottom ).wait();
	*/

	queue.enqueueReadBuffer( cl_final_result, CL_TRUE, 0, sizeof( double ), difference_value );

	// If we are at the end of the program, ouput data to files.
	if ( writeFile ) {
		double Cs[ 4 ] = { 0, 0, 0, 0 };

		// Warp images.
		FImage *warp_base = new FImage( base->width, base->height, base->metadata );
		FImage *warp_sliver = new FImage( sliver->width, sliver->height, sliver->metadata );
		base->warpBase ( warp_base, As, Bs, Cs, 1, S_WRITE_INTERP );
		sliver->warpSliver ( warp_sliver, As, Bs, Cs, S_WRITE_INTERP );

		double* RCs = new double[ sliver->width + sliver->height ];
		queue.enqueueReadBuffer ( cl_rcs, CL_TRUE, 0, sizeof ( double ) * ( sliver->width + sliver->height ), RCs );

		// Apply RC Values.
		for ( int i = 1; i <= maxX; i++ ) {
			double z_change = RCs[ i - 1 ];
			for ( int j = 0; j < warp_sliver->height; ++j ) {
				if ( warp_sliver->fastGet( i, j ) == -std::numeric_limits<double>::infinity() ) continue;
				warp_sliver->fastSet( i, j, warp_sliver->fastGet( i, j ) + z_change );
			}
		}
		for ( int j = minY; j <= maxY; ++j ) {
			double z_change = RCs[ sliver->width + j - minY ];
			for ( int i = 0; i < warp_base->width; i++ ) {
				if ( warp_base->fastGet( i, j ) == -std::numeric_limits<double>::infinity() )	continue;
				warp_base->fastSet( i, j, warp_base->fastGet( i, j ) + z_change );
			}
		}

		double base_min = warp_base->getMin();
		double base_max = warp_base->getMax();
		double sliver_min = warp_sliver->getMin();
		double sliver_max = warp_sliver->getMax();
		double min_val = base_min < sliver_min ? base_min : sliver_min;
		double max_val = base_max < sliver_max ? base_max : sliver_max;
		double slope = 1 / ( max_val - min_val );

		// Write warped and corrected images.
		warp_base->writeImage( "w_base.tif" );
		warp_sliver->writeImage( "w_sliver.tif" );

		// Write viewable images.
		warp_base->writeDisplayableImage( "dw_base.tif", slope, min_val );
		warp_sliver->writeDisplayableImage( "dw_sliver.tif", slope, min_val );

		delete[] RCs;
		delete warp_base;
		delete warp_sliver;
	}

	return difference_value[ 0 ];

#else	// Non-OpenCL (Intel's Thread Building Blocks and sequential) versions of fastZ correction


	// Bicubic sliver warping
#ifdef S_MODE_TBB
	tbb::parallel_for( 0, maxX + 1, 1, [&]( int x )	{
#elif defined( S_MODE_SEQUENTIAL )
	for ( int x = 0; x < maxX + 1; x++ ) {
#endif
		double sx[ 4 ]; // these are the distances for each row
		double sy[ 4 ]; // these are the distances for each column
		double ux[ 4 ]; // these are the weight for the rows
		double uy[ 4 ]; // these are the weights for the columns

		double yoffset = -( Bs[ 1 ] - 1 ) * x;
		double origX = x - As[ 1 ] * x;

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
			double pixelvalue = 0;
			double origY = y + yoffset;
			if ( floor ( origX ) == 0 || floor ( origY ) == 0 || floor ( origX ) >= sliver->width - 2 || floor ( origY ) >= sliver->height - 2 ) {
				pixelvalue = sliver->interpPixel( origX, origY, F_BILINEAR );
			}
			else {
				for ( int j = 0; j < 4; ++j ) {
					for ( int i = 0; i < 4; ++i ) {
						pixelvalue += sliver->fastGet( floor( origX ) - 1 + i, floor( origY ) - 1 + j ) * ux[ i ] * uy[ j ];
					}
				}
			}
			sliver_array[ x + sliver->width * y ] = pixelvalue;
		}
	}
#ifdef S_MODE_TBB
	);
#endif

	// Bicubic base warping.
#ifdef S_MODE_TBB
	tbb::parallel_for( minY, maxY + 1, 1, [&]( int y )	{
#elif defined( S_MODE_SEQUENTIAL )
	for ( int y = minY; y < maxY + 1; ++y ) {
#endif
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
			if ( floor ( origX ) == 0 || floor ( origY ) == 0 || floor ( origX ) >= base->width - 2 || floor ( origY ) >= base->height - 2 ) {
				pixelvalue = base->interpPixel( origX, origY, F_BILINEAR );
			}
			else {
				for ( int j = 0; j < 4; ++j ) {
					for ( int i = 0; i < 4; ++i ) {
						pixelvalue += base->fastGet( floor( origX ) - 1 + i, floor( origY ) - 1 + j ) * ux[ i ] * uy[ j ];
					}
				}
			}
			base_array[ x + sliver->width * y ] = pixelvalue;
		}
	}
#ifdef S_MODE_TBB
	);
#endif


	double skr = 0;
	double skr_error = 0;
	// Column Sums
#ifdef S_MODE_TBB
	tbb::parallel_for( 0, m, 1, [&]( int i )	{
#elif defined( S_MODE_SEQUENTIAL )
	for ( int i = 0; i < m; i++ ) {
#endif
		double sum = 0;
		double sum_error = 0;
		for ( int j = 0; j < n; ++j ) {
			double diff = base_array[ i + sliver->width * ( j + minY ) ] - sliver_array[ i + sliver->width * ( j + minY ) ];
			double y = diff - sum_error;
			double t = sum + y;
			sum_error = ( t - sum ) - y;
			sum = t;
		}
		KRCs[ i ] = sum;
	}
#ifdef S_MODE_TBB
	);
#endif


	// Row Sums
#ifdef S_MODE_TBB
	tbb::parallel_for( 0, n, 1, [&]( int j )	{
#elif defined( S_MODE_SEQUENTIAL )
	for ( int j = 0; j < n; ++j ) {
#endif
		double sum = 0;
		double sum_error = 0;
		for ( int i = 0; i < m; ++i ) {
			double diff = sliver_array[ i + sliver->width * ( j + minY ) ] - base_array[ i + sliver->width * ( j + minY ) ];
			double y = diff - sum_error;
			double t = sum + y;
			sum_error = ( t - sum ) - y;
			sum = t;
		}
		KRCs[ m + j ] = sum;
	}
#ifdef S_MODE_TBB
	);
#endif

	// KRC Sum
// #ifdef S_MODE_TBB
//	tbb::parallel_for( 0, n, 1, [&]( int j )	{
// #elif defined( S_MODE_SEQUENTIAL )
	for ( int j = 0; j < n; ++j ) {
// #endif
//#ifdef S_MODE_TBB
//		cs.lock();
//#endif
		double y = KRCs[ m + j ] - skr_error;
		double t = skr + y;
		skr_error = ( t - skr ) - y;
		skr = t;
//#ifdef S_MODE_TBB
//		cs.unlock();
//#endif
	}
//#ifdef S_MODE_TBB
//	);
//#endif

	// Matrix Computation
//#ifdef S_MODE_TBB
//	tbb::parallel_for( 1, m, 1, [&]( int i )	{
//#elif defined( S_MODE_SEQUENTIAL )
	for ( int i = 1; i < m; ++i ) {
//#endif
		RCs[ i - 1 ] = inv_n * ( KRCs[ i ] - KRCs[ 0 ] );
	}
//#ifdef S_MODE_TBB
//	);
//#endif

	// RC Sums.
//#ifdef S_MODE_TBB
//	tbb::parallel_for( m, size + 1, 1, [&]( int i ) {
//#elif defined( S_MODE_SEQUENTIAL )
	for ( int i = m; i <= size; ++i ) {
//#endif
		RCs[ i - 1 ] = -inv_n * KRCs[ 0 ] - inv_mn * skr + inv_m * KRCs[ i ];
	}
//#ifdef S_MODE_TBB
//	);
//#endif

	// B-Spline sliver warping.
#ifdef S_MODE_TBB
	tbb::parallel_for( 0, maxX + 1, 1, [&]( int x )	{
#elif defined( S_MODE_SEQUENTIAL )
	for ( int x = 0; x < maxX + 1; x++ ) {
#endif
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
			else {
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
			if ( floor ( origX ) == 0 || floor ( origY ) == 0 || floor ( origX ) >= sliver->width - 2 || floor ( origY ) >= sliver->height - 2 ) {
				pixelvalue = sliver->interpPixel( origX, origY, F_BILINEAR );
			}
			else {
				for ( int j = 0; j < 4; ++j ) {
					for ( int i = 0; i < 4; ++i ) {
						pixelvalue += sliver->fastGet( floor( origX ) - 1 + i, floor( origY ) - 1 + j ) * ux[ i ] * uy[ j ];
					}
				}
			}
			sliver_array[ x + sliver->width * y ] = pixelvalue;
		}
	}
#ifdef S_MODE_TBB
	);
#endif

	// B-Spline base warping.
#ifdef S_MODE_TBB
	tbb::parallel_for( minY, maxY + 1, 1, [&]( int y )	{
#elif defined( S_MODE_SEQUENTIAL )
	for ( int y = minY; y < maxY + 1; ++y ) {
#endif
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
			if ( floor ( origX ) == 0 || floor ( origY ) == 0 || floor ( origX ) >= base->width - 2 || floor ( origY ) >= base->height - 2 ) {
				pixelvalue = base->interpPixel( origX, origY, F_BILINEAR );
			}
			else {
				for ( int j = 0; j < 4; ++j ) {
					for ( int i = 0; i < 4; ++i ) {
						pixelvalue += base->fastGet( floor( origX ) - 1 + i, floor( origY ) - 1 + j ) * ux[ i ] * uy[ j ];
					}
				}
			}
			base_array[ x + sliver->width * y ] = pixelvalue;
		}
	}
#ifdef S_MODE_TBB
	);
#endif

	// Sum of Differences and weights computation.
	double sum = 0;
	double area = 0;
	double sum_error = 0;
	double area_error = 0;

#ifdef S_MODE_TBB
	tbb::parallel_for( 0, maxX + 1, 1, [&]( int i )	{
#elif defined( S_MODE_SEQUENTIAL )
	for ( int i = 0; i <= maxX; ++i ) {
#endif
		double c_change = i == 0 ? 0 : RCs[ i - 1 ];
		for ( int j = minY; j <= maxY; ++j ) {
			double r_change = RCs[ maxX + j - minY ];
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
			double val = diff * diff * weight;
#ifdef S_MODE_TBB
			cs.lock();
#endif
			double y = val - sum_error;
			double t = sum + y;
			sum_error = ( t - sum ) - y;
			sum = t;

			y = weight - area_error;
			t = area + y;
			area_error = ( t - area ) - y;
			area = t;
#ifdef S_MODE_TBB
			cs.unlock();
#endif
		}
	}
#ifdef S_MODE_TBB
	);
#endif

	// If we are at the end of the program, ouput data to files.
	if ( writeFile ) {
		std::ofstream output;
		// Debugging info.
		output.open( "matrix.txt", std::ios::out );
		if ( output.is_open() ) {
			output << "MinY : " << minY << "\t MaxY : " << maxY << "\t MaxX : " << maxX << "\n";
			output << "Difference : " << sum / area << "\n";
			output << "\n\n";
			output << "This is the generated vector : \n\n";
			for ( int i = 0; i < size; ++i ) {
				output << KRCs[ i ] << "\n";
			}
			output << "\n\n";
			output << "These are the computed C's : \n\n";
			for ( int i = 0; i < m - 1; ++i ) {
				output << RCs[ i ] << "\n";
			}
			output << "\n\n";
			output << "These are the computed R's : \n\n";
			for ( int i = m - 1; i < size; ++i ) {
				output << RCs[ i ] << "\n";
			}
		}
		output.close();

		double Cs[ 4 ] = { 0, 0, 0, 0 };

		// Warp images.
		FImage *warp_base = new FImage( base->width, base->height, base->metadata );
		FImage *warp_sliver = new FImage( sliver->width, sliver->height, sliver->metadata );
		base->warpBase ( warp_base, As, Bs, Cs, 1, S_WRITE_INTERP );
		sliver->warpSliver ( warp_sliver, As, Bs, Cs, S_WRITE_INTERP );

		// Apply RC Values.
		for ( int i = 1; i <= maxX; i++ ) {
			double z_change = RCs[ i - 1 ];
			for ( int j = 0; j < warp_sliver->height; ++j ) {
				if ( warp_sliver->fastGet( i, j ) == -std::numeric_limits<double>::infinity() ) continue;
				warp_sliver->fastSet( i, j, warp_sliver->fastGet( i, j ) + z_change );
			}
		}
		for ( int j = minY; j <= maxY; ++j ) {
			double z_change = RCs[ maxX + j - minY ];
			for ( int i = 0; i < warp_base->width; i++ ) {
				if ( warp_base->fastGet( i, j ) == -std::numeric_limits<double>::infinity() )	continue;
				warp_base->fastSet( i, j, warp_base->fastGet( i, j ) + z_change );
			}
		}

		double base_min = warp_base->getMin();
		double base_max = warp_base->getMax();
		double sliver_min = warp_sliver->getMin();
		double sliver_max = warp_sliver->getMax();
		double min_val = base_min < sliver_min ? base_min : sliver_min;
		double max_val = base_max < sliver_max ? base_max : sliver_max;
		double slope = 1 / ( max_val - min_val );

		// Write warped and corrected images.
		warp_base->writeImage( "w_base.tif" );
		warp_sliver->writeImage( "w_sliver.tif" );

		// Write viewable images.
		warp_base->writeDisplayableImage( "dw_base.tif", slope, min_val );
		warp_sliver->writeDisplayableImage( "dw_sliver.tif", slope, min_val );

		delete warp_base;
		delete warp_sliver;
	}

	return ( sum / area );

#endif

}

/**
 *	This is the function that is minimized by simplex.  What it does, conceptually, is similar to performing a
 *	polynomial warp on BOTH the base image and the sliver, taking the sum of differences squared between them (normalized
 *	by area), and returning that.  In fact what it does is a little faster and more efficient.
 *
 *	First, it "warps" the sliver by calculating, for each point in the sliver, where that point would end up,
 *	IF the transformation were done in reverse.  Here the sliver points (xs, ys) are transformed to (xsp, ysp),
 *	for "prime".  It takes only the liniar portion for this transformation.
 *
 *	Then, it "warps" the base image by doing the full 3rd order polynomial warp on (xsp,ysp) to yield (xspp, yspp).
 *	The actual value is gotten from the main image by interpolation.
 *
 *	I tried writing this function using floats (instead of doubles) for some of the x and y coordinates, and found
 *	that the answer never converged beyond about 8 decimal places.  Plus, using floats instead of doubles didn't ACTUALLY
 *	make it ANY faster.  So doubles it is.
 */
double slowZDifference( double vertex[] ) {
	difference_evals++;  //increment global variable used to count function evaluations

	double sum = 0;
	double area = 0;
	int basexstart = basesize[ 0 ] / 2 - sliversize[ 0 ] / 2;
	//	float ymin=0.0,ymax=(float) sliversize[1]-1;
	double ymin = 0.0, ymax = ( double ) sliversize[ 1 ] - 1;

	double C0 = vertex[ 8 ];
	double C1 = vertex[ 9 ];
	double C2 = vertex[ 10 ];
	double C3 = vertex[ 11 ];

#ifdef S_MODE_TBB
	tbb::parallel_for( 1, sliversize[ 1 ], 1, [&]( int ys ) { //think about whether to reverse these two loops.
#else
	for ( int ys = 1; ys < sliversize[ 1 ] - 1; ++ys )	{
#endif
		//sum+= C0*C0 + (C1-1)*(C1-1) + (C2-2)*(C2-2) +(C3-3)*(C3-3);
		//for (int ys=0; ys<sliversize[1];++ys) { //think about whether to reverse these two loops.
		for ( int xs = 0; xs < sliversize[ 0 ]; ++xs ) {   //fix so not starting at zero every time, only goes as far as needed
			//reverse the linear transform to fix sliver drift
			double weight;
			double ysp = ( double ) ( ys + ( vertex[ 3 ] - 1.0 )*xs );//be sure to get these signs correct!
			double xsp = ( double ) ( xs + vertex[ 2 ] * xs );//7/2/2008; I think this should be +.
			//double ysp = (double) (ys); //- (vertex[3]-1.0)*xs);  //be sure to get these signs correct!
			//double xsp = (double) (xs); //- vertex[2]*xs);
			double ysp2 = ysp*ysp;
			double ysp3 = ysp2*ysp;
			//apply polynomial transform to account for image shift
			double xspp = ( double ) ( xsp + vertex[ 0 ] + vertex[ 2 ] * ysp + vertex[ 4 ] * ysp2 + vertex[ 6 ] * ysp3 );
			double yb = ( double ) ( vertex[ 1 ] + vertex[ 3 ] * ysp + vertex[ 5 ] * ysp2 + vertex[ 7 ] * ysp3 );
			//transform to base image coordinates
			double xb = basexstart + xspp;
			//see if new point falls on the existing base image;
			if ( ( yb >= ymin ) && ( yb <= ymax ) ) {
				weight = 1L; //
				if ( ( yb - ymin ) < 1.0 ) weight = ( yb - ymin );
				if ( ( ymax - yb ) < 1.0 ) weight = ( ymax - yb );
				//int Zs = sliver->get(rounds(xsp),rounds(ysp));
				double Zs = sliver->get( xs, ys );
				double Zm = base->interpPixel( xb, yb, F_BILINEAR );

				//double diff = sliver->get(xs,ys) - interp_pixel(base,xb,yb);
				double diff = ( C0 + Zm + C1*yb + C2*yb*yb + C3*yb*yb*yb ) - ( Zs - C1*xsp );

				//	printf("MYFUNC VALS: %d %f %f\n", Zs, Zm, diff);
				//	getchar();
#ifdef S_MODE_TBB
				cs.lock();
#endif
				sum += diff * diff * weight;
				area += weight;
#ifdef S_MODE_TBB
				cs.unlock();
#endif
			}
		}
	} 
#ifdef S_MODE_TBB
	);
#endif

	return ( sum / area );
}

/**
 *	Apply fast or slow Z correction depending upon global flag set.
 *	
 *	@param	vertex	An array storing parameter values.
 *	@return		The calculated difference.
 */
double determineDifference( double vertex[] ) {
	difference_evals++;  //increment global variable used to count function evaluations
	if ( fastZ )
		return fastZDifference( vertex, false );
	else
		return slowZDifference( vertex );
}

/**
 *	Make a copy of the base.
 *
 *	@param	_base	The FImage object to copy.
 */
void copyBase( FImage *_base ) {
	base = new FImage( _base->width, _base->height, _base->metadata );
	for ( unsigned int i = 0; i < _base->width; i++ ) {
		for ( unsigned int j = 0; j < _base->height; j++ ) {
			base->set( i, j, _base->get( i, j ) );
		}
	}
}

/**
 *	Make a copy of the sliver.
 *
 *	@param	_sliver	The FImage object to copy.
 */
void copySliver( FImage *_sliver ) {
	sliver = new FImage( _sliver->width, _sliver->height, _sliver->metadata );
	for ( unsigned int i = 0; i < _sliver->width; i++ ) {
		for ( unsigned int j = 0; j < _sliver->height; j++ ) {
			sliver->set( i, j, _sliver->get( i, j ) );
		}
	}
}

/**
 *	Performs the simplex routine. FImage files are inputted as well as a starting vector which the simplex is trying to refine. We determine
 *	what type of Z correction is being made with the _fastZ boolean flag. Depending upon the type of correction, either 8 parameters will be
*	optimized or else 12 parameters will be optimized. Various simplex parameters are provided as well for reflection, contraction
 *	etc.
 *
 *	@param	_base	Input base image.
 *	@param	_sliver	Input sliver image.
 *	@param	input_vector	Original parameter values usually given by grid-search.
 *	@param	return_vector	An array the size of input_vector to store the final calculated values.
 *	@param	_fastZ	Boolean indicator that determines if we perform fast Z or slow Z correction.
 *	@param	precision
 *	@param	alpha
 *	@param	beta
 *	@param	gamma
 *	@param	errhalt
 *	@param	maxi
 *
 *	@return		0 if success
 */
int simplex( FImage *_base, FImage *_sliver, double input_vector[], double return_vector[], bool _fastZ, double precision[], double alpha, double beta, double gamma,
			 double errhalt, int maxi ) {

	fastZ = _fastZ;
	int n;
	if ( fastZ ) {
		n = 8;
	}
	else {
		n = 12;
	}

	int i, j;
	


	double** simplex_matrix = new double*[ n + 1 ];
	for ( i = 0; i < n + 1; ++i )	simplex_matrix[ i ] = new double[ n + 1 ];
	double* differences = new double[ n + 1 ];
	double* vertex_reflection = new double[ n + 1 ];
	double* vertex_expansion = new double[ n + 1 ];
	double* centroid = new double[ n + 1 ];
	double reflection_difference, expansion_difference, contraction_difference;

	printf( "Begin Simplex Routine!\n" );
	copyBase( _base );
	copySliver( _sliver );
	sliversize = new int[ 2 ];
	sliversize[ 0 ] = sliver->width;
	sliversize[ 1 ] = sliver->height;
	basesize = new int[ 2 ];
	basesize[ 0 ] = base->width;
	basesize[ 1 ] = base->height;

#ifdef S_MODE_OPENCL
	initializeOpenCLContext();
#else
	base_array = new double[ sliver->width * sliver->height ];
	sliver_array = new double[ sliver->width * sliver->height ];
	KRCs = new double[ sliver->width + sliver->height ];
	RCs = new double[ sliver->width + sliver->height ];
	#ifdef S_MODE_TBB
	if ( !init.is_active() ) {
		init.initialize();
		printf( "TBB Thread Scheduler Initialized : %d Threads.\n", init.default_num_threads() );
	}
	else {
		printf( "TBB Thread Scheduler Already Exists : %d Threads.\n", init.default_num_threads() );
	}
	#endif
#endif

	//int ok;	Never used
	int iterations;
	int max_index, min_index;
	double c;
	double max_difference;
	double min_difference;
	double deviation;
	double mean_difference;
	double next_max_difference;

	/***************************************************/
	//fill the simplex
	for ( i = 0; i < n; i++ ) {
		simplex_matrix[ 0 ][ i ] = input_vector[ i ];
	}

	for ( i = n; i < n + 1; i++ ) {
		vertex_reflection[ i ] = 0.0;
		vertex_expansion[ i ] = 0.0;
		centroid[ i ] = 0.0;
		differences[ i ] = 0.0;
	}
	/* fill in the remaining points of the initial Simplex            */
	/* point no. k is the initial point displaced by 2x along coord k */

	for ( i = 1; i <= n; i++ ) {
		for ( j = 0; j < n; j++ ) {
			if ( ( j + 1 ) != i )
				simplex_matrix[ i ][ j ] = simplex_matrix[ 0 ][ j ];
			else {
				double temp = simplex_matrix[ 0 ][ j ];
				//if(abs(t0) < TINY) t0 = 0.5 ;
				simplex_matrix[ i ][ j ] = precision[ j ] + temp;
			}
		}
	}

	/* go through all of the points of the Simplex & find max, min */
	max_index = 0;
	min_index = 0;
	differences[ 0 ] = determineDifference( input_vector );
	max_difference = differences[ 0 ];
	min_difference = max_difference;
	iterations = 0;

	//min_difference = F_L, max_difference = F_A
	for ( i = 1; i <= n; i++ ) {
		for ( j = 0; j < n; j++ )
			vertex_reflection[ j ] = simplex_matrix[ i ][ j ];
		differences[ i ] = determineDifference( vertex_reflection );
		if ( differences[ i ] >= max_difference ) {
			max_difference = differences[ i ];
			next_max_difference = differences[ max_index ];
			//x_h[i]=max_difference;
			//		x_a[i] = differences[max_index];
			max_index = i;
		}
		if ( differences[ i ] < min_difference ) {
			min_difference = differences[ i ];
			//					x_l[i]=max_difference;
			min_index = i;
		}
	}
	for ( j = 0; j < n; j++ ) {
		centroid[ j ] = 0.00;
		for ( i = 0; i <= n; i++ ) {
			if ( i != max_index )
				centroid[ j ] += simplex_matrix[ i ][ j ];
		}
		centroid[ j ] = centroid[ j ] / ( ( double ) n );
	}
	
	mean_difference = 0.00;
	deviation = 0.00;
	for ( i = 0; i <= n; i++ ) {
		mean_difference += differences[ i ];
	}
	mean_difference = mean_difference / ( ( double ) ( n + 1 ) );
	for ( i = 0; i <= n; i++ ) {
		deviation += ( mean_difference - differences[ i ] ) * ( mean_difference - differences[ i ] );
	}
	deviation = sqrt( deviation / ( ( double ) ( n + 1 ) ) );
	// double error = ( max_difference - min_difference ) / max( max_difference, 1.0 );

	double* param_differences = new double[ n ];
	double** past_param_extrema = new double*[ n ];
	double** current_param_extrema = new double*[ n ];
	for ( i = 0; i < n; ++i ) {
		past_param_extrema[ i ] = new double[ 2 ];
		past_param_extrema[ i ][ 0 ] = INFINITY;
		past_param_extrema[ i ][ 1 ] = INFINITY;
		current_param_extrema[ i ] = new double[ 2 ];
	}
	double max_param_difference = INFINITY;
	double max_successive_difference = INFINITY;

	double* normalize_factors = new double[ n ];
	normalize_factors[ 0 ] = 1.0;
	normalize_factors[ 1 ] = 1.0;
	normalize_factors[ 2 ] = sliver->height;
	normalize_factors[ 3 ] = sliver->height;
	normalize_factors[ 4 ] = sliver->height * sliver->height;
	normalize_factors[ 5 ] = sliver->height * sliver->height;
	normalize_factors[ 6 ] = sliver->height * sliver->height * sliver->height;
	normalize_factors[ 7 ] = sliver->height * sliver->height * sliver->height;

	//while ( ( deviation >= errhalt ) && ( iterations < maxi + 1 ) ) {
	while ( ( ( max_param_difference > 1e-13 ) || ( max_successive_difference >= 1e-10 ) ) && ( iterations < maxi + 1 ) ) {
		debug = false;
		/*new code*/
		///*******************************************************///
		///************************reflection*********************///
		///*******************************************************///
		//calculate the points of the reflected simplex
		for ( i = 0; i < n; i++ ) {
			c = centroid[ i ];
			vertex_reflection[ i ] = c + alpha * ( c - simplex_matrix[ max_index ][ i ] );
		}

		reflection_difference = determineDifference( vertex_reflection );

		//if the reflected simplex is decent, use it...but not if its very good.
		//if the reflected simplex is very good, try an expansion
		if ( min_difference <= reflection_difference && reflection_difference < next_max_difference ) {
			for ( i = 0; i < n; i++ ) {
				simplex_matrix[ max_index ][ i ] = vertex_reflection[ i ];
			}
			differences[ max_index ] = reflection_difference;
			operations++;
			//printf("R");
			debug = true;
		}
		///*******************************************************///
		///************************greedy expansion***************///
		///*********************non smooth function***************///
		/*if(reflection_difference < min_difference){
		for(i=0; i<n; i++){
		c = centroid[i];
		vertex_expansion[i] = c + gamma*(vertex_reflection[i]-c);
		}
		expansion_difference = determineDifference(vertex_expansion);
		if(expansion_difference < min_difference){
		for(i =0; i<n; i++){Simplex[max_index][i]=vertex_expansion[i];}
		max_difference = expansion_difference ;
		printf("e");
		operations++;
		debug = true;
		goto do_more;
		}
		else{
		for(i =0; i<n; i++){Simplex[max_index][i]=vertex_reflection[i];}
		max_difference = reflection_difference ;
		printf("r");
		operations++;
		debug = true;
		goto do_more;
		}
		}*/

		///*******************************************************///
		///*********************greedy minimization***************///
		///************************smooth function****************///
		//if the reflected point produced a simplex with a smaller value than the current minimum,
		//check to see if an expansion produces an even better simplex
		else if ( reflection_difference < min_difference ) {		//never make this <= ~ Nathan
			for ( j = 0; j < n; j++ ) {
				c = centroid[ j ];
				vertex_expansion[ j ] = c + gamma * ( vertex_reflection[ j ] - c );
			}
			expansion_difference = determineDifference( vertex_expansion );
			//if the expanded simplex is better, use it.
			if ( expansion_difference < reflection_difference ) {
				for ( i = 0; i < n; i++ ) {
					simplex_matrix[ max_index ][ i ] = vertex_expansion[ i ];
				}
				max_difference = expansion_difference;
				//printf("e");
				operations++;
				debug = true;
			}
			//if the expanded simplex isn't better, use the reflected simplex, this will keep the
			//simplex smaller
			else if ( expansion_difference >= reflection_difference ) {
				for ( i = 0; i < n; i++ ) {
					simplex_matrix[ max_index ][ i ] = vertex_reflection[ i ];
				}
				max_difference = reflection_difference;
				//printf("r");
				operations++;
				debug = true;
			}
		}
		/*Nathan's Contraction and Shrink code to include both contraction cases*/
		///*******************************************************///
		///************************Contraction********************///
		///*******************************************************///
		//try A contraction if the reflected value is greater than the second greatest value
		//(which from the other situations it actually has to be...if statement might not be
		//needed)
		else if ( reflection_difference >= next_max_difference ) {
			//cout << "\nerror impending\nFXR: " << reflection_difference << "\nFXA: " << next_max_difference << "\nymax: " << max_difference;
			bool le = true;
			if ( reflection_difference < max_difference ) {
				le = true;
				for ( i = 0; i < n; i++ ) {
					double c = centroid[ i ];
					vertex_expansion[ i ] = c + ( beta * ( vertex_reflection[ i ] - c ) );
				}
			}
			if ( reflection_difference >= max_difference ) {
				le = false;
				for ( i = 0; i < n; i++ ) {
					double c = centroid[ i ];
					vertex_expansion[ i ] = c + ( beta * ( simplex_matrix[ max_index ][ i ] - c ) );
				}
			}

			expansion_difference = determineDifference( vertex_expansion );
			if ( le == true ) {
				if ( expansion_difference <= reflection_difference ) {
					//printf("c");
					operations++;
					debug = true;
					//cout << "case1\n";
					//getchar();
					for ( i = 0; i < n; i++ ) {
						simplex_matrix[ max_index ][ i ] = vertex_expansion[ i ];
						differences[ max_index ] = reflection_difference;
						max_difference = reflection_difference;
					}
				}
				else {
					/* contraction failed; collapse simplex points */
					/* towards the present minimum point */
					//printf("S-c");
					operations++;
					debug = true;
					for ( i = 0; i <= n; i++ ) {
						for ( j = 0; j < n; j++ ) {
							simplex_matrix[ i ][ j ] = 0.5 * ( simplex_matrix[ i ][ j ] + simplex_matrix[ min_index ][ j ] );
							vertex_reflection[ j ] = simplex_matrix[ i ][ j ];
						}
						reflection_difference = determineDifference( vertex_reflection );
						differences[ i ] = reflection_difference;
					}
				}
			}
			if ( le == false ) {
				if ( expansion_difference < max_difference ) {
					//printf("C");
					operations++;
					debug = true;
					//cout << "case2\n";
					//getchar();
					for ( i = 0; i < n; i++ ) {
						simplex_matrix[ max_index ][ i ] = vertex_expansion[ i ];
						differences[ max_index ] = expansion_difference;
						max_difference = expansion_difference;
					}
				}
				else {
					/* contraction failed; collapse simplex points */
					/* towards the present minimum point */
					//printf("S-C");
					operations++;
					debug = true;
					for ( i = 0; i <= n; i++ ) {
						for ( j = 0; j < n; j++ ) {
							simplex_matrix[ i ][ j ] = 0.5 * ( simplex_matrix[ i ][ j ] + simplex_matrix[ min_index ][ j ] );
							vertex_reflection[ j ] = simplex_matrix[ i ][ j ];
						}
						reflection_difference = determineDifference( vertex_reflection );
						differences[ i ] = reflection_difference;
					}
				}
			}
		}
		//end contraction

		/* compute the "standard error" of the Simplex */
		/* which is the termination criterion for the  */
		/* outermost (while) loop.                     */
		mean_difference = 0.00;
		deviation = 0.00;
		for ( i = 0; i <= n; i++ )
			mean_difference += differences[ i ];
		mean_difference /= ( ( double ) ( n + 1 ) );

		for ( i = 0; i <= n; i++ )
			deviation += ( mean_difference - differences[ i ] ) * ( mean_difference - differences[ i ] );
		deviation = sqrt( deviation / ( ( double ) ( n + 1 ) ) );
		// error = ( max_difference - min_difference ) / max( max_difference, 1 );

		max_index = min_index = 0;
		max_difference = min_difference = differences[ 0 ];
		for ( i = 1; i <= n; i++ ) {
			if ( differences[ i ] >= max_difference ) {
				max_difference = differences[ i ];
				//get the second highest value
				next_max_difference = differences[ max_index ];
				max_index = i;
			}
			/*if(differences[i] < min_difference)*/
			else if ( differences[ i ] < min_difference ) {
				min_difference = differences[ i ];
				min_index = i;
			}
		}

		/* compute the centroid */
		for ( j = 0; j < n; j++ ) {
			centroid[ j ] = 0.00;
			for ( i = 0; i <= n; i++ ) {
				if ( i != max_index )
					centroid[ j ] += simplex_matrix[ i ][ j ];
			}
			centroid[ j ] /= ( ( double ) n );
		}

		for ( i = 0; i < n; ++i ) {	// For each individual parameter
			param_differences[ i ] = 0.0;	// Zero out parameter differences
			current_param_extrema[ i ][ 0 ] = simplex_matrix[ 1 ][ i ];
			current_param_extrema[ i ][ 1 ] = simplex_matrix[ 1 ][ i ];
			for ( j = 1; j <= n; ++j ) {	// Go through each vertex
				if ( simplex_matrix[ j ][ i ] < current_param_extrema[ i ][ 0 ] ) {
					current_param_extrema[ i ][ 0 ] = simplex_matrix[ j ][ i ];	// Get minimum value of parameter
				}
				else if ( simplex_matrix[ j ][ i ] > current_param_extrema[ i ][ 1 ] ) {
					current_param_extrema[ i ][ 1 ] = simplex_matrix[ j ][ i ];	// Get maximum value of parameter
				}
			}
			param_differences[ i ] = ( current_param_extrema[ i ][ 1 ] - current_param_extrema[ i ][ 0 ] ) * normalize_factors[ i ];	// Param_Diff = Max - Min
		}

		max_param_difference = 0.0;
		for ( i = 0; i < n; ++i ) {
			if ( param_differences[ i ] > max_param_difference ) {
				max_param_difference = param_differences[ i ];
			}
		}

		int limiting_param = -1;

		if ( ( iterations % 100 ) == 0 ) {
			max_successive_difference = 0.0;
			for ( i = 0; i < n; ++i ) {
				double diff_mins = abs( current_param_extrema[ i ][ 0 ] - past_param_extrema[ i ][ 0 ] ) * normalize_factors[ i ];
				double diff_maxs = abs( current_param_extrema[ i ][ 1 ] - past_param_extrema[ i ][ 1 ] ) * normalize_factors[ i ];
				double diff_min_max = abs( current_param_extrema[ i ][ 0 ] - past_param_extrema[ i ][ 1 ] ) * normalize_factors[ i ];
				double diff_max_min = abs( current_param_extrema[ i ][ 1 ] - past_param_extrema[ i ][ 0 ] ) * normalize_factors[ i ];
				if ( diff_mins > max_successive_difference ) {
					max_successive_difference = diff_mins;
					limiting_param = i;
				}
				if ( diff_maxs > max_successive_difference ) {
					max_successive_difference = diff_maxs;
					limiting_param = i;
				}
				if ( diff_min_max > max_successive_difference ) {
					max_successive_difference = diff_min_max;
					limiting_param = i;
				}
				if ( diff_max_min > max_successive_difference ) {
					max_successive_difference = diff_max_min;
					limiting_param = i;
				}
			}

			for ( i = 0; i < n; ++i ) {
				past_param_extrema[ i ][ 0 ] = current_param_extrema[ i ][ 0 ];
				past_param_extrema[ i ][ 1 ] = current_param_extrema[ i ][ 1 ];
			}
		}

		if ( ( iterations % ( PRINTEM ) ) == 0 ) {
			printf( "ITERATION %4d    Best Difference = %e\n", iterations, min_difference );
			printf( "%d Operations completed.\n", operations );
			operations = 0;
			printf( "Current Best Parameters : \n" );
			printf( "%e %e %e %e\n", simplex_matrix[ min_index ][ 0 ], simplex_matrix[ min_index ][ 2 ], simplex_matrix[ min_index ][ 4 ], simplex_matrix[ min_index ][ 6 ] );
			printf( "%e %e %e %e\n", simplex_matrix[ min_index ][ 1 ], simplex_matrix[ min_index ][ 3 ], simplex_matrix[ min_index ][ 5 ], simplex_matrix[ min_index ][ 7 ] );
			// printf( "%e %e %e %e\n", simplex_matrix[ min_index ][ 8 ], simplex_matrix[ min_index ][ 9 ], simplex_matrix[ min_index ][ 10 ], simplex_matrix[ min_index ][ 11 ] );
			for ( int i = 0; i < n + 1; ++i ) {
				printf( "Difference %d : %.20e\n", i, differences[ i ] );
			}
			printf( "Average Difference : %.20e\n", mean_difference );
			printf( "Deviation : %.20e\n", deviation );
			printf( "Normalized A0 Difference : %.20e\n", param_differences[ 0 ] );
			printf( "Normalized B0 Difference : %.20e\n", param_differences[ 1 ] );
			printf( "Normalized A1 Difference : %.20e\n", param_differences[ 2 ] );
			printf( "Normalized B1 Difference : %.20e\n", param_differences[ 3 ] );
			printf( "Normalized A2 Difference : %.20e\n", param_differences[ 4 ] );
			printf( "Normalized B2 Difference : %.20e\n", param_differences[ 5 ] );
			printf( "Normalized A3 Difference : %.20e\n", param_differences[ 6 ] );
			printf( "Normalized B3 Difference : %.20e\n", param_differences[ 7 ] );
			printf( "Highest Parameter Difference : %.20e\n", max_param_difference );
			printf( "Highest Successive Difference Resulting from Param %d : %.20e\n", limiting_param, max_successive_difference );
			printf( "\n" );
		}

		if ( debug == false ) {
			printf( "\n" );
			printf( "Reflection Difference : %e\n", reflection_difference );
			printf( "Expansion Difference : %e\n", expansion_difference );
			printf( "Best Difference : %e\n", max_difference );
			printf( "Worst Difference : %e\n", min_difference );
			printf( "Second Worst Difference : %e\n", next_max_difference );
			getchar();
		}
		iterations++;
	} /* end of while loop */

	/* find biggest and smallest values of determineDifference */
	max_index = 0;
	min_index = 0;
	max_difference = differences[ 0 ];
	min_difference = max_difference;
	for ( i = 1; i <= n; i++ ) {
		if ( differences[ i ] >= max_difference ) {
			max_difference = differences[ i ];
			max_index = i;
		}
		if ( differences[ i ] < min_difference ) {
			min_difference = differences[ i ];
			min_index = i;
		}
	}

	/* put the minimum point in the return vector z */
	for ( i = 0; i < n; i++ ) {
		return_vector[ i ] = simplex_matrix[ min_index ][ i ];
	}

	if ( fastZ ) {
		fastZDifference( return_vector, true );
	}

	printf( "\n\nSimplex took %ld evaluations.\n", difference_evals );
	printf( "%e %e %e %e\n", simplex_matrix[ min_index ][ 0 ], simplex_matrix[ min_index ][ 2 ], simplex_matrix[ min_index ][ 4 ], simplex_matrix[ min_index ][ 6 ] );
	printf( "%e %e %e %e\n", simplex_matrix[ min_index ][ 1 ], simplex_matrix[ min_index ][ 3 ], simplex_matrix[ min_index ][ 5 ], simplex_matrix[ min_index ][ 7 ] );


#ifdef S_MODE_OPENCL
	delete[] difference_value;
	delete clWarpSliverCubic;
	delete clWarpBaseCubic;
	delete clWarpSliverBspline;
	delete clWarpBaseBspline;
	delete clRCResults;
	delete clPartialDifference;
	delete clFinalDifference;
	// delete clDifferenceCalculation;
#else
	delete[] base_array;
	delete[] sliver_array;
	delete[] KRCs;
	delete[] RCs;
#endif
	delete[] param_differences;
	for ( i = 0; i < n; ++i ) {
		delete[] past_param_extrema[ i ];
		delete[] current_param_extrema[ i ];
	}
	delete[] past_param_extrema;
	delete[] current_param_extrema;
	delete[] normalize_factors;

	for ( i = 0; i < n + 1; ++i )	delete[] simplex_matrix[ i ];
	delete[] simplex_matrix;
	delete[] differences;
	delete[] vertex_reflection;
	delete[] vertex_expansion;
	delete[] centroid;
	delete[] basesize;
	delete[] sliversize;
	delete base;
	delete sliver;
	return ( 0 );
}
