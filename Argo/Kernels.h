#ifndef __KERNELS__
	#define __KERNELS__

	#include <string>


	#ifndef __CODEFOLD__
		#define __CODEFOLD__
	#endif

	// BICUBIC WARPING KERNEL STRING
	#ifdef __CODEFOLD__

	const std::string bicubic_kernels_string = R"(

	#ifdef cl_amd_fp64
	#pragma OPENCL EXTENSION cl_amd_fp64 : enable
	#elif defined( cl_khr_fp64 )
	#pragma OPENCL EXTENSION cl_khr_fp64 : enable
	#endif
	#ifdef cl_amd_printf
	#pragma OPENCL EXTENSION cl_amd_printf : enable
	#endif

	#define A -0.5

	#define DIVISION_SLIVER 1
	#define DIVISION_BASE 1

	/*
	 *	Performs bicubic warping of the sliver image given Nelder-Mead optimized parameters. Makes use of local
	 *	data.
	 *
	 *	@param sliver			Global original unwarped sliver.
	 *	@param sliver_array		Global warped sliver.
	 *	@param width			Sliver width.
	 *	@param height			Sliver height.
	 *	@param maxX				The maximum possible x-value (pre-computed) given parameters.
	 *	@param minY				The minimum possible y-value (pre-computed) given parameters.
	 *	@param maxY				The maximum possible y-value (pre-computed) given parameters.
	 *
	 */
	__kernel void clWarpSliverCubic( __global __const double* sliver, __global double* sliver_array, double A1, double B1,
							int width, int height, int maxX, int minY, int maxY )	{

		__local int local_size;
		local_size = get_local_size( 0 );

		__private int local_id = get_local_id( 0 );
	
		__private int x = get_group_id( 0 ) + ( get_global_size( 0 ) / local_size ) * ( local_id * DIVISION_SLIVER / local_size );
		x = ( x > maxX + 1 ) ? maxX : x;

		local_id = local_id - ( local_id / ( local_size / DIVISION_SLIVER ) ) * ( local_size / DIVISION_SLIVER );
		__private double yoffset = -( B1 - 1 ) * x;
		__private double origX = x - A1 * x;
	

		__private double sx0 = fabs( origX - floor( origX - 1 ) ); // these get the distances for each row and column from the initial point, positive
		__private double sy0 = fabs( yoffset - floor( yoffset - 1 ) );
		__private double sx1 = fabs( 1 - sx0 );
		__private double sy1 = fabs( 1 - sy0 );
		__private double sx2 = fabs( 1 - sx1 );
		__private double sy2 = fabs( 1 - sy1 );
		__private double sx3 = 1 + sx2;
		__private double sy3 = 1 + sy2;
		__private double ux0;
		__private double ux1;
		__private double ux2;
		__private double ux3;
		__private double uy0;
		__private double uy1;
		__private double uy2;
		__private double uy3;

		__private double x2 = sx0 * sx0;
		__private double x3 = x2 * sx0;
		__private double y2 = sy0 * sy0;
		__private double y3 = y2 * sy0;
		if ( sx0 <= 1 && sx0 >= 0 ) {
			ux0 = ( A + 2 ) * x3 - ( A + 3 ) * x2 + 1;
		}
		else if ( sx0 > 1 && sx0 <= 2 ) {
			ux0 = A * x3 - 5 * A * x2 + 8 * A * sx0 - 4 * A;
		}
		else {
			ux0 = 0;
		}
		if ( sy0 <= 1 && sy0 >= 0 ) {
			uy0 = ( A + 2 ) * y3 - ( A + 3 ) * y2 + 1;
		}
		else if ( sy0 > 1 && sy0 <= 2 ) {
			uy0 = A * y3 - 5 * A * y2 + 8 * A * sy0 - 4 * A;
		}
		else {
			uy0 = 0;
		}

		x2 = sx1 * sx1;
		x3 = x2 * sx1;
		y2 = sy1 * sy1;
		y3 = y2 * sy1;
		if ( sx1 <= 1 && sx1 >= 0 ) {
			ux1 = ( A + 2 ) * x3 - ( A + 3 ) * x2 + 1;
		}
		else if ( sx1 > 1 && sx1 <= 2 ) {
			ux1 = A * x3 - 5 * A * x2 + 8 * A * sx1 - 4 * A;
		}
		else {
			ux1 = 0;
		}
		if ( sy1 <= 1 && sy1 >= 0 ) {
			uy1 = ( A + 2 ) * y3 - ( A + 3 ) * y2 + 1;
		}
		else if ( sy1 > 1 && sy1 <= 2 ) {
			uy1 = A * y3 - 5 * A * y2 + 8 * A * sy1 - 4 * A;
		}
		else {
			uy1 = 0;
		}

		x2 = sx2 * sx2;
		x3 = x2 * sx2;
		y2 = sy2 * sy2;
		y3 = y2 * sy2;
		if ( sx2 <= 1 && sx2 >= 0 ) {
			ux2 = ( A + 2 ) * x3 - ( A + 3 ) * x2 + 1;
		}
		else if ( sx2 > 1 && sx2 <= 2 ) {
			ux2 = A * x3 - 5 * A * x2 + 8 * A * sx2 - 4 * A;
		}
		else {
			ux2 = 0;
		}
		if ( sy2 <= 1 && sy2 >= 0 ) {
			uy2 = ( A + 2 ) * y3 - ( A + 3 ) * y2 + 1;
		}
		else if ( sy2 > 1 && sy2 <= 2 ) {
			uy2 = A * y3 - 5 * A * y2 + 8 * A * sy2 - 4 * A;
		}
		else {
			uy2 = 0;
		}

		x2 = sx3 * sx3;
		x3 = x2 * sx3;
		y2 = sy3 * sy3;
		y3 = y2 * sy3;
		if ( sx3 <= 1 && sx3 >= 0 ) {
			ux3 = ( A + 2 ) * x3 - ( A + 3 ) * x2 + 1;
		}
		else if ( sx3 > 1 && sx3 <= 2 ) {
			ux3 = A * x3 - 5 * A * x2 + 8 * A * sx3 - 4 * A;
		}
		else {
			ux3 = 0;
		}
		if ( sy3 <= 1 && sy3 >= 0 ) {
			uy3 = ( A + 2 ) * y3 - ( A + 3 ) * y2 + 1;
		}
		else if ( sy3 > 1 && sy3 <= 2 ) {
			uy3 = A * y3 - 5 * A * y2 + 8 * A * sy3 - 4 * A;
		}
		else {
			uy3 = 0;
		}

		for( int y = local_id + minY; y < maxY + 1; y += local_size / DIVISION_SLIVER )	{
			double origY = y + yoffset;
			int floor_origX = ( int ) origX;
			int floor_origY = ( int ) origY;
			int ceil_origX = ceil( origX );
			int ceil_origY = ceil( origY );
			int bool_flag = floor_origX == 0 || floor_origY == 0 || floor_origX >= width - 2 || floor_origY >= height - 2;

			double pixelvalue = 0;
		
			double s = origX - floor_origX;
			double t = origY - floor_origY;

			// Bilinear Interpolation
			double left_val, right_val;
			double top_val, bottom_val;
			// Interpolate across the top edge
			left_val = sliver[ floor_origX + width * floor_origY ];
			right_val = sliver[ ceil_origX + width * floor_origY ];
			// Linear interpolation recoded to only one multiply
			top_val = s * right_val + ( 1 - s ) * ( left_val );
			// Interpolate across the bottom edge
			left_val = sliver[ floor_origX + width * ceil_origY ];
			right_val = sliver[ ceil_origX + width * ceil_origY ];
			bottom_val = s * right_val + ( 1 - s ) * ( left_val );
			// Interpolate between top and bottom
			pixelvalue = ( t * bottom_val + ( 1 - t ) * top_val );
			
			// Bicubic Sum
			int temp_index = width * floor_origY;

			if( !bool_flag )	{
			
				pixelvalue = 0;

				pixelvalue += sliver[ floor_origX - 1 + temp_index - width ] * ux0 * uy0;
				pixelvalue += sliver[ floor_origX + temp_index - width ] * ux1 * uy0;
				pixelvalue += sliver[ floor_origX + 1 + temp_index - width ] * ux2 * uy0;
				pixelvalue += sliver[ floor_origX + 2 + temp_index - width ] * ux3 * uy0;

				pixelvalue += sliver[ floor_origX - 1 + temp_index ] * ux0 * uy1;
				pixelvalue += sliver[ floor_origX + temp_index ] * ux1 * uy1;
				pixelvalue += sliver[ floor_origX + 1 + temp_index ] * ux2 * uy1;
				pixelvalue += sliver[ floor_origX + 2 + temp_index ] * ux3 * uy1;

				pixelvalue += sliver[ floor_origX - 1 + temp_index + width ] * ux0 * uy2;
				pixelvalue += sliver[ floor_origX + temp_index + width ] * ux1 * uy2;
				pixelvalue += sliver[ floor_origX + 1 + temp_index + width ] * ux2 * uy2;
				pixelvalue += sliver[ floor_origX + 2 + temp_index + width ] * ux3 * uy2;

				pixelvalue += sliver[ floor_origX - 1 + temp_index + width + width ] * ux0 * uy3;
				pixelvalue += sliver[ floor_origX + temp_index + width + width ] * ux1 * uy3;
				pixelvalue += sliver[ floor_origX + 1 + temp_index + width + width ] * ux2 * uy3;
				pixelvalue += sliver[ floor_origX + 2 + temp_index + width + width ] * ux3 * uy3;

			}

			// Choose Between Bilinear and Bicubic based upon flag.
			sliver_array[ x + width * y ] = pixelvalue;
		}
	}

	/*
	 *	Performs bicubic warping of the sliver image given Nelder-Mead optimized parameters. Makes use of local
	 *	data.
	 *
	 *	@param base				Global original unwarped base.
	 *	@param base_array		Global warped base.
	 *	@param width			Sliver width.
	 *	@param height			Sliver height.
	 *	@param bwidth			Base width.
	 *	@param maxX				The maximum possible x-value (pre-computed) given parameters.
	 *	@param minY				The minimum possible y-value (pre-computed) given parameters.
	 *	@param maxY				The maximum possible y-value (pre-computed) given parameters.
	 *
	 */
	__kernel void clWarpBaseCubic( __global __const double* base, __global double* base_array, double A0, double A1, double A2, double A3, double B0, double B1, double B2, double B3,
							int width, int height, int bwidth, int xbp, int maxX, int minY, int maxY )	{
		__local int local_size;
		local_size = get_local_size( 0 );

		__private int local_id = get_local_id( 0 );
		__private int y = get_group_id( 0 ) + ( get_global_size( 0 ) / local_size ) * ( local_id * DIVISION_BASE / local_size ) + minY;
		y = ( y > maxY + 1 ) ? maxY : y;

		local_id = local_id - ( local_id / ( local_size / DIVISION_BASE ) ) * ( local_size / DIVISION_BASE );

		__private int ym2 = y * y;
		__private int ym3 = ym2 * y;
		__private double xoffset = A0 + A1 * y + A2 * ym2 + A3 * ym3;
		__private double origY = B0 + B1 * y + B2 * ym2 + B3 * ym3;

		__private double sx0 = fabs( xoffset - floor( xoffset - 1 ) ); // these get the distances for each row and column from the initial point, positive
		__private double sy0 = fabs( origY - floor( origY - 1 ) );
		__private double sx1 = fabs( 1 - sx0 );
		__private double sy1 = fabs( 1 - sy0 );
		__private double sx2 = fabs( 1 - sx1 );
		__private double sy2 = fabs( 1 - sy1 );
		__private double sx3 = 1 + sx2;
		__private double sy3 = 1 + sy2;
		__private double ux0;
		__private double ux1;
		__private double ux2;
		__private double ux3;
		__private double uy0;
		__private double uy1;
		__private double uy2;
		__private double uy3;

		__private double x2 = sx0 * sx0;
		__private double x3 = x2 * sx0;
		__private double y2 = sy0 * sy0;
		__private double y3 = y2 * sy0;
		if ( sx0 <= 1 && sx0 >= 0 ) {
			ux0 = ( A + 2 ) * x3 - ( A + 3 ) * x2 + 1;
		}
		else if ( sx0 > 1 && sx0 <= 2 ) {
			ux0 = A * x3 - 5 * A * x2 + 8 * A * sx0 - 4 * A;
		}
		else {
			ux0 = 0;
		}
		if ( sy0 <= 1 && sy0 >= 0 ) {
			uy0 = ( A + 2 ) * y3 - ( A + 3 ) * y2 + 1;
		}
		else if ( sy0 > 1 && sy0 <= 2 ) {
			uy0 = A * y3 - 5 * A * y2 + 8 * A * sy0 - 4 * A;
		}
		else {
			uy0 = 0;
		}

		x2 = sx1 * sx1;
		x3 = x2 * sx1;
		y2 = sy1 * sy1;
		y3 = y2 * sy1;
		if ( sx1 <= 1 && sx1 >= 0 ) {
			ux1 = ( A + 2 ) * x3 - ( A + 3 ) * x2 + 1;
		}
		else if ( sx1 > 1 && sx1 <= 2 ) {
			ux1 = A * x3 - 5 * A * x2 + 8 * A * sx1 - 4 * A;
		}
		else {
			ux1 = 0;
		}
		if ( sy1 <= 1 && sy1 >= 0 ) {
			uy1 = ( A + 2 ) * y3 - ( A + 3 ) * y2 + 1;
		}
		else if ( sy1 > 1 && sy1 <= 2 ) {
			uy1 = A * y3 - 5 * A * y2 + 8 * A * sy1 - 4 * A;
		}
		else {
			uy1 = 0;
		}

		x2 = sx2 * sx2;
		x3 = x2 * sx2;
		y2 = sy2 * sy2;
		y3 = y2 * sy2;
		if ( sx2 <= 1 && sx2 >= 0 ) {
			ux2 = ( A + 2 ) * x3 - ( A + 3 ) * x2 + 1;
		}
		else if ( sx2 > 1 && sx2 <= 2 ) {
			ux2 = A * x3 - 5 * A * x2 + 8 * A * sx2 - 4 * A;
		}
		else {
			ux2 = 0;
		}
		if ( sy2 <= 1 && sy2 >= 0 ) {
			uy2 = ( A + 2 ) * y3 - ( A + 3 ) * y2 + 1;
		}
		else if ( sy2 > 1 && sy2 <= 2 ) {
			uy2 = A * y3 - 5 * A * y2 + 8 * A * sy2 - 4 * A;
		}
		else {
			uy2 = 0;
		}

		x2 = sx3 * sx3;
		x3 = x2 * sx3;
		y2 = sy3 * sy3;
		y3 = y2 * sy3;
		if ( sx3 <= 1 && sx3 >= 0 ) {
			ux3 = ( A + 2 ) * x3 - ( A + 3 ) * x2 + 1;
		}
		else if ( sx3 > 1 && sx3 <= 2 ) {
			ux3 = A * x3 - 5 * A * x2 + 8 * A * sx3 - 4 * A;
		}
		else {
			ux3 = 0;
		}
		if ( sy3 <= 1 && sy3 >= 0 ) {
			uy3 = ( A + 2 ) * y3 - ( A + 3 ) * y2 + 1;
		}
		else if ( sy3 > 1 && sy3 <= 2 ) {
			uy3 = A * y3 - 5 * A * y2 + 8 * A * sy3 - 4 * A;
		}
		else {
			uy3 = 0;
		}

		for( int x = local_id; x < maxX + 1; x += local_size / DIVISION_BASE )	{
			double origX = x + xbp + xoffset;

			int floor_origX = ( int ) origX;
			int floor_origY = ( int ) origY;
			int ceil_origX = ceil( origX );
			int ceil_origY = ceil( origY );
			int bool_flag = floor_origX == 0 || floor_origY == 0 || floor_origX >= bwidth - 2 || floor_origY >= height - 2;

			double pixelvalue = 0;

			double s = origX - floor_origX;
			double t = origY - floor_origY;

			// Bilinear Interpolation
			double left_val, right_val;
			double top_val, bottom_val;
			// Interpolate across the top edge
			left_val = base[ floor_origX + bwidth * floor_origY ];
			right_val = base[ ceil_origX + bwidth * floor_origY ];
			// Linear interpolation recoded to only one multiply
			top_val = s * right_val + ( 1 - s ) * ( left_val );
			// Interpolate across the bottom edge
			left_val = base[ floor_origX + bwidth * ceil_origY ];
			right_val = base[ ceil_origX + bwidth * ceil_origY ];
			bottom_val = s * right_val + ( 1 - s ) * ( left_val );
			// Interpolate between top and bottom
			pixelvalue = ( t * bottom_val + ( 1 - t ) * top_val );

			int temp_index = bwidth * floor_origY;

			if( !bool_flag )	{

				pixelvalue = 0;

				pixelvalue += base[ floor_origX - 1 + temp_index - bwidth ] * ux0 * uy0;
				pixelvalue += base[ floor_origX + temp_index - bwidth ] * ux1 * uy0;
				pixelvalue += base[ floor_origX + 1 + temp_index - bwidth ] * ux2 * uy0;
				pixelvalue += base[ floor_origX + 2 + temp_index - bwidth ] * ux3 * uy0;

				pixelvalue += base[ floor_origX - 1 + temp_index ] * ux0 * uy1;
				pixelvalue += base[ floor_origX + temp_index ] * ux1 * uy1;
				pixelvalue += base[ floor_origX + 1 + temp_index ] * ux2 * uy1;
				pixelvalue += base[ floor_origX + 2 + temp_index ] * ux3 * uy1;

				pixelvalue += base[ floor_origX - 1 + temp_index + bwidth ] * ux0 * uy2;
				pixelvalue += base[ floor_origX + temp_index + bwidth ] * ux1 * uy2;
				pixelvalue += base[ floor_origX + 1 + temp_index + bwidth ] * ux2 * uy2;
				pixelvalue += base[ floor_origX + 2 + temp_index + bwidth ] * ux3 * uy2;

				pixelvalue += base[ floor_origX - 1 + temp_index + bwidth + bwidth ] * ux0 * uy3;
				pixelvalue += base[ floor_origX + temp_index + bwidth + bwidth ] * ux1 * uy3;
				pixelvalue += base[ floor_origX + 1 + temp_index + bwidth + bwidth ] * ux2 * uy3;
				pixelvalue += base[ floor_origX + 2 + temp_index + bwidth + bwidth ] * ux3 * uy3;

			}

			base_array[ x + width * y ] = pixelvalue;
		}
	}


	)";

	#endif

	// BSPLINE WARPING KERNEL STRING!!!
	#ifdef __CODEFOLD__
	const std::string bspline_kernels_string = R"(



	#define B 1
	#define C 0

	/*
	 *	Performs cubic b-spline warping of the base image given Nelder-Mead optimized parameters. Makes use of local
	 *	data.
	 *
	 *	@param sliver			Global original unwarped sliver.
	 *	@param sliver_array		Global warped sliver.
	 *	@param A1								
	 *	@param B1				
	 *	@param width			Sliver width.
	 *	@param height			Sliver height.
	 *	@param maxX				The maximum possible x-value (pre-computed) given parameters.
	 *	@param minY				The minimum possible y-value (pre-computed) given parameters.
	 *	@param maxY				The maximum possible y-value (pre-computed) given parameters.
	 *
	 */
	__kernel void clWarpSliverBspline( __global __const double* sliver, __global double* sliver_array, double A1, double B1,
							int width, int height, int maxX, int minY, int maxY, __local double* sx, __local double* sy, __local double* ux, __local double* uy )	{
		__local int local_size;
		__local int iterations;
		__local double yoffset;
		__local double origX;
		__local int x;
		__private int y = get_local_id( 0 );

		x = get_group_id( 0 );
		local_size = get_local_size( 0 );
		iterations = ceil( ( double ) ( maxY - minY + 1 ) / ( double ) get_local_size( 0 ) );
		yoffset = -( B1 - 1 ) * x;
		origX = x - A1 * x;

		__private double sx0 = fabs( origX - floor( origX - 1 ) ); // these get the distances for each row and column from the initial point, positive
		__private double sy0 = fabs( yoffset - floor( yoffset - 1 ) );
		__private double sx1 = fabs( 1 - sx0 );
		__private double sy1 = fabs( 1 - sy0 );
		__private double sx2 = fabs( 1 - sx1 );
		__private double sy2 = fabs( 1 - sy1 );
		__private double sx3 = 1 + sx2;
		__private double sy3 = 1 + sy2;
		__private double ux0;
		__private double ux1;
		__private double ux2;
		__private double ux3;
		__private double uy0;
		__private double uy1;
		__private double uy2;
		__private double uy3;

		__private double x2 = sx0 * sx0;
		__private double x3 = x2 * sx0;
		__private double y2 = sy0 * sy0;
		__private double y3 = y2 * sy0;
		if ( sx0 <= 1 && sx0 >= 0 ) {
			ux0 = ( ( 12 - ( 9 * B ) - ( 6 * C ) ) * x3 + ( -18 + ( 12 * B ) + ( 6 * C ) ) * x2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sx0 > 1 && sx0 <= 2 ) {
			ux0 =  ( ( -B - ( 6 * C ) ) * x3 + ( ( 6 * B ) + ( 30 * C ) ) * x2 + ( -( 12 * B ) - ( 48 * C ) )*( sx0 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			ux0 = 0;
		}
		if ( sy0 <= 1 && sy0 >= 0 ) {
			uy0 = ( ( 12 - ( 9 * B ) - ( 6 * C ) )*y3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*y2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sy0 > 1 && sy0 <= 2 ) {
			uy0 = ( ( -B - ( 6 * C ) )*y3 + ( ( 6 * B ) + ( 30 * C ) )*y2 + ( -( 12 * B ) - ( 48 * C ) )*( sy0 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			uy0 = 0;
		}

		x2 = sx1 * sx1;
		x3 = x2 * sx1;
		y2 = sy1 * sy1;
		y3 = y2 * sy1;
		if ( sx1 <= 1 && sx1 >= 0 ) {
			ux1 = ( ( 12 - ( 9 * B ) - ( 6 * C ) ) * x3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*x2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sx1 > 1 && sx1 <= 2 ) {
			ux1 =  ( ( -B - ( 6 * C ) )*x3 + ( ( 6 * B ) + ( 30 * C ) )*x2 + ( -( 12 * B ) - ( 48 * C ) )*( sx1 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			ux1 = 0;
		}
		if ( sy1 <= 1 && sy1 >= 0 ) {
			uy1 = ( ( 12 - ( 9 * B ) - ( 6 * C ) )*y3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*y2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sy1 > 1 && sy1 <= 2 ) {
			uy1 = ( ( -B - ( 6 * C ) )*y3 + ( ( 6 * B ) + ( 30 * C ) )*y2 + ( -( 12 * B ) - ( 48 * C ) )*( sy1 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			uy1 = 0;
		}

		x2 = sx2 * sx2;
		x3 = x2 * sx2;
		y2 = sy2 * sy2;
		y3 = y2 * sy2;
		if ( sx2 <= 1 && sx2 >= 0 ) {
			ux2 = ( ( 12 - ( 9 * B ) - ( 6 * C ) ) * x3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*x2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sx2 > 1 && sx2 <= 2 ) {
			ux2 =  ( ( -B - ( 6 * C ) )*x3 + ( ( 6 * B ) + ( 30 * C ) )*x2 + ( -( 12 * B ) - ( 48 * C ) )*( sx2 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			ux2 = 0;
		}
		if ( sy2 <= 1 && sy2 >= 0 ) {
			uy2 = ( ( 12 - ( 9 * B ) - ( 6 * C ) )*y3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*y2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sy2 > 1 && sy2 <= 2 ) {
			uy2 = ( ( -B - ( 6 * C ) )*y3 + ( ( 6 * B ) + ( 30 * C ) )*y2 + ( -( 12 * B ) - ( 48 * C ) )*( sy2 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			uy2 = 0;
		}

		x2 = sx3 * sx3;
		x3 = x2 * sx3;
		y2 = sy3 * sy3;
		y3 = y2 * sy3;
		if ( sx3 <= 1 && sx3 >= 0 ) {
			ux3 = ( ( 12 - ( 9 * B ) - ( 6 * C ) ) * x3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*x2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sx3 > 1 && sx3 <= 2 ) {
			ux3 = ( ( -B - ( 6 * C ) )*x3 + ( ( 6 * B ) + ( 30 * C ) )*x2 + ( -( 12 * B ) - ( 48 * C ) )*( sx3 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			ux3 = 0;
		}
		if ( sy3 <= 1 && sy3 >= 0 ) {
			uy3 = ( ( 12 - ( 9 * B ) - ( 6 * C ) )*y3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*y2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sy3 > 1 && sy3 <= 2 ) {
			uy3 = ( ( -B - ( 6 * C ) )*y3 + ( ( 6 * B ) + ( 30 * C ) )*y2 + ( -( 12 * B ) - ( 48 * C ) )*( sy3 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			uy3 = 0;
		}

		for( int it = 0; it < iterations; ++it )	{
			double origY = y * iterations + it + minY + yoffset;
			if( y * iterations + it + minY < maxY + 1 )	{
				int floor_origX = ( int ) origX;
				int floor_origY = ( int ) origY;
				int ceil_origX = ceil( origX );
				int ceil_origY = ceil( origY );
				int bool_flag = floor_origX == 0 || floor_origY == 0 || floor_origX >= width - 2 || floor_origY >= height - 2;

				double pixelvalue = 0;
		
				double s = origX - floor_origX;
				double t = origY - floor_origY;

				// Bilinear Interpolation
				double left_val, right_val;
				double top_val, bottom_val;
				// Interpolate across the top edge
				left_val = sliver[ floor_origX + width * floor_origY ];
				right_val = sliver[ ceil_origX + width * floor_origY ];
				// Linear interpolation recoded to only one multiply
				top_val = s * right_val + ( 1 - s ) * ( left_val );
				// Interpolate across the bottom edge
				left_val = sliver[ floor_origX + width * ceil_origY ];
				right_val = sliver[ ceil_origX + width * ceil_origY ];
				bottom_val = s * right_val + ( 1 - s ) * ( left_val );
				// Interpolate between top and bottom
				pixelvalue = ( t * bottom_val + ( 1 - t ) * top_val );
			
				// Bicubic Sum
				int temp_index = width * floor_origY;

				if( !bool_flag )	{

					pixelvalue = 0;

					pixelvalue += sliver[ floor_origX - 1 + temp_index - width ] * ux0 * uy0;
					pixelvalue += sliver[ floor_origX + temp_index - width ] * ux1 * uy0;
					pixelvalue += sliver[ floor_origX + 1 + temp_index - width ] * ux2 * uy0;
					pixelvalue += sliver[ floor_origX + 2 + temp_index - width ] * ux3 * uy0;

					pixelvalue += sliver[ floor_origX - 1 + temp_index ] * ux0 * uy1;
					pixelvalue += sliver[ floor_origX + temp_index ] * ux1 * uy1;
					pixelvalue += sliver[ floor_origX + 1 + temp_index ] * ux2 * uy1;
					pixelvalue += sliver[ floor_origX + 2 + temp_index ] * ux3 * uy1;

					pixelvalue += sliver[ floor_origX - 1 + temp_index + width ] * ux0 * uy2;
					pixelvalue += sliver[ floor_origX + temp_index + width ] * ux1 * uy2;
					pixelvalue += sliver[ floor_origX + 1 + temp_index + width ] * ux2 * uy2;
					pixelvalue += sliver[ floor_origX + 2 + temp_index + width ] * ux3 * uy2;

					pixelvalue += sliver[ floor_origX - 1 + temp_index + width + width ] * ux0 * uy3;
					pixelvalue += sliver[ floor_origX + temp_index + width + width ] * ux1 * uy3;
					pixelvalue += sliver[ floor_origX + 1 + temp_index + width + width ] * ux2 * uy3;
					pixelvalue += sliver[ floor_origX + 2 + temp_index + width + width ] * ux3 * uy3;

				}

				// Choose Between Bilinear and Bicubic based upon flag.
				sliver_array[ x + width * ( y * iterations + it + minY ) ] = pixelvalue;
			}
		}
	}

	/*
	 *	Performs cubic b-spline warping of the base image given Nelder-Mead optimized parameters. Makes use of local
	 *	data.
	 *
	 *	@param base				Global original unwarped base.
	 *	@param base_array		Global warped base.
	 *	@param width			Sliver width.
	 *	@param height			Sliver height.
	 *	@param bwidth			Base width.
	 *	@param maxX				The maximum possible x-value (pre-computed) given parameters.
	 *	@param minY				The minimum possible y-value (pre-computed) given parameters.
	 *	@param maxY				The maximum possible y-value (pre-computed) given parameters.
	 *
	 */
	__kernel void clWarpBaseBspline( __global __const double* base, __global double* base_array, double A0, double A1, double A2, double A3, double B0, double B1, double B2, double B3,
							int width, int height, int bwidth, int xbp, int maxX, int minY, int maxY, __local double* sx, __local double* sy, __local double* ux, __local double* uy )	{
		__local int local_size;
		__local int iterations;
		__local double xoffset;
		__local double origY;
		__local int y;
		__private int x = get_local_id( 0 );

		y = get_group_id( 0 ) + minY;
		local_size = get_local_size(0);
		iterations = ceil( ( double ) ( maxX + 1 ) / ( double ) get_local_size( 0 ) );

		__private int ym2 = y * y;
		__private int ym3 = ym2 * y;
		xoffset = A0 + A1 * y + A2 * ym2 + A3 * ym3;
		origY = B0 + B1 * y + B2 * ym2 + B3 * ym3;

		__private double sx0 = fabs( xoffset - floor( xoffset - 1 ) ); // these get the distances for each row and column from the initial point, positive
		__private double sy0 = fabs( origY - floor( origY - 1 ) );
		__private double sx1 = fabs( 1 - sx0 );
		__private double sy1 = fabs( 1 - sy0 );
		__private double sx2 = fabs( 1 - sx1 );
		__private double sy2 = fabs( 1 - sy1 );
		__private double sx3 = 1 + sx2;
		__private double sy3 = 1 + sy2;
		__private double ux0;
		__private double ux1;
		__private double ux2;
		__private double ux3;
		__private double uy0;
		__private double uy1;
		__private double uy2;
		__private double uy3;

		__private double x2 = sx0 * sx0;
		__private double x3 = x2 * sx0;
		__private double y2 = sy0 * sy0;
		__private double y3 = y2 * sy0;
		if ( sx0 <= 1 && sx0 >= 0 ) {
			ux0 = ( ( 12 - ( 9 * B ) - ( 6 * C ) ) * x3 + ( -18 + ( 12 * B ) + ( 6 * C ) ) * x2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sx0 > 1 && sx0 <= 2 ) {
			ux0 =  ( ( -B - ( 6 * C ) ) * x3 + ( ( 6 * B ) + ( 30 * C ) ) * x2 + ( -( 12 * B ) - ( 48 * C ) )*( sx0 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			ux0 = 0;
		}
		if ( sy0 <= 1 && sy0 >= 0 ) {
			uy0 = ( ( 12 - ( 9 * B ) - ( 6 * C ) )*y3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*y2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sy0 > 1 && sy0 <= 2 ) {
			uy0 = ( ( -B - ( 6 * C ) )*y3 + ( ( 6 * B ) + ( 30 * C ) )*y2 + ( -( 12 * B ) - ( 48 * C ) )*( sy0 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			uy0 = 0;
		}

		x2 = sx1 * sx1;
		x3 = x2 * sx1;
		y2 = sy1 * sy1;
		y3 = y2 * sy1;
		if ( sx1 <= 1 && sx1 >= 0 ) {
			ux1 = ( ( 12 - ( 9 * B ) - ( 6 * C ) ) * x3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*x2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sx1 > 1 && sx1 <= 2 ) {
			ux1 =  ( ( -B - ( 6 * C ) )*x3 + ( ( 6 * B ) + ( 30 * C ) )*x2 + ( -( 12 * B ) - ( 48 * C ) )*( sx1 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			ux1 = 0;
		}
		if ( sy1 <= 1 && sy1 >= 0 ) {
			uy1 = ( ( 12 - ( 9 * B ) - ( 6 * C ) )*y3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*y2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sy1 > 1 && sy1 <= 2 ) {
			uy1 = ( ( -B - ( 6 * C ) )*y3 + ( ( 6 * B ) + ( 30 * C ) )*y2 + ( -( 12 * B ) - ( 48 * C ) )*( sy1 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			uy1 = 0;
		}

		x2 = sx2 * sx2;
		x3 = x2 * sx2;
		y2 = sy2 * sy2;
		y3 = y2 * sy2;
		if ( sx2 <= 1 && sx2 >= 0 ) {
			ux2 = ( ( 12 - ( 9 * B ) - ( 6 * C ) ) * x3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*x2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sx2 > 1 && sx2 <= 2 ) {
			ux2 =  ( ( -B - ( 6 * C ) )*x3 + ( ( 6 * B ) + ( 30 * C ) )*x2 + ( -( 12 * B ) - ( 48 * C ) )*( sx2 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			ux2 = 0;
		}
		if ( sy2 <= 1 && sy2 >= 0 ) {
			uy2 = ( ( 12 - ( 9 * B ) - ( 6 * C ) )*y3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*y2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sy2 > 1 && sy2 <= 2 ) {
			uy2 = ( ( -B - ( 6 * C ) )*y3 + ( ( 6 * B ) + ( 30 * C ) )*y2 + ( -( 12 * B ) - ( 48 * C ) )*( sy2 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			uy2 = 0;
		}

		x2 = sx3 * sx3;
		x3 = x2 * sx3;
		y2 = sy3 * sy3;
		y3 = y2 * sy3;
		if ( sx3 <= 1 && sx3 >= 0 ) {
			ux3 = ( ( 12 - ( 9 * B ) - ( 6 * C ) ) * x3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*x2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sx3 > 1 && sx3 <= 2 ) {
			ux3 = ( ( -B - ( 6 * C ) )*x3 + ( ( 6 * B ) + ( 30 * C ) )*x2 + ( -( 12 * B ) - ( 48 * C ) )*( sx3 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			ux3 = 0;
		}
		if ( sy3 <= 1 && sy3 >= 0 ) {
			uy3 = ( ( 12 - ( 9 * B ) - ( 6 * C ) )*y3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*y2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sy3 > 1 && sy3 <= 2 ) {
			uy3 = ( ( -B - ( 6 * C ) )*y3 + ( ( 6 * B ) + ( 30 * C ) )*y2 + ( -( 12 * B ) - ( 48 * C ) )*( sy3 ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			uy3 = 0;
		}

		for( int it = 0; it < iterations; ++it )	{
			double origX = x * iterations + it + xbp + xoffset;

			if( x * iterations + it < maxX + 1 )	{

				int floor_origX = ( int ) origX;
				int floor_origY = ( int ) origY;
				int ceil_origX = ceil( origX );
				int ceil_origY = ceil( origY );
				int bool_flag = floor_origX == 0 || floor_origY == 0 || floor_origX >= bwidth - 2 || floor_origY >= height - 2;

				double pixelvalue = 0;

				double s = origX - floor_origX;
				double t = origY - floor_origY;

				// Bilinear Interpolation
				double left_val, right_val;
				double top_val, bottom_val;
				// Interpolate across the top edge
				left_val = base[ floor_origX + bwidth * floor_origY ];
				right_val = base[ ceil_origX + bwidth * floor_origY ];
				// Linear interpolation recoded to only one multiply
				top_val = s * right_val + ( 1 - s ) * ( left_val );
				// Interpolate across the bottom edge
				left_val = base[ floor_origX + bwidth * ceil_origY ];
				right_val = base[ ceil_origX + bwidth * ceil_origY ];
				bottom_val = s * right_val + ( 1 - s ) * ( left_val );
				// Interpolate between top and bottom
				pixelvalue = ( t * bottom_val + ( 1 - t ) * top_val );

				int temp_index = bwidth * floor_origY;

				if( !bool_flag )	{

					pixelvalue = 0;

					pixelvalue += base[ floor_origX - 1 + temp_index - bwidth ] * ux0 * uy0;
					pixelvalue += base[ floor_origX + temp_index - bwidth ] * ux1 * uy0;
					pixelvalue += base[ floor_origX + 1 + temp_index - bwidth ] * ux2 * uy0;
					pixelvalue += base[ floor_origX + 2 + temp_index - bwidth ] * ux3 * uy0;

					pixelvalue += base[ floor_origX - 1 + temp_index ] * ux0 * uy1;
					pixelvalue += base[ floor_origX + temp_index ] * ux1 * uy1;
					pixelvalue += base[ floor_origX + 1 + temp_index ] * ux2 * uy1;
					pixelvalue += base[ floor_origX + 2 + temp_index ] * ux3 * uy1;

					pixelvalue += base[ floor_origX - 1 + temp_index + bwidth ] * ux0 * uy2;
					pixelvalue += base[ floor_origX + temp_index + bwidth ] * ux1 * uy2;
					pixelvalue += base[ floor_origX + 1 + temp_index + bwidth ] * ux2 * uy2;
					pixelvalue += base[ floor_origX + 2 + temp_index + bwidth ] * ux3 * uy2;

					pixelvalue += base[ floor_origX - 1 + temp_index + bwidth + bwidth ] * ux0 * uy3;
					pixelvalue += base[ floor_origX + temp_index + bwidth + bwidth ] * ux1 * uy3;
					pixelvalue += base[ floor_origX + 1 + temp_index + bwidth + bwidth ] * ux2 * uy3;
					pixelvalue += base[ floor_origX + 2 + temp_index + bwidth + bwidth ] * ux3 * uy3;

				}

				base_array[ x * iterations + it + width * y ] = pixelvalue;
			}
		}
	}


	)";

	#endif

	// SUMMATION KERNEL STRING!!!
	#ifdef __CODEFOLD__

	const std::string summation_kernel_string = R"(

	/*
	 *	Performs the sum of differences of values within each individual column. We fix an x-value
	 *	and compute kC_x = sum_{i=0}^{n} (base(x,i) - sliver(x,i)). Each individual work group computes
	 *	num_blocks kC values. We use Kahan Summation implementation to achieve reproducibility in results.
	 *
	 *	*Note: This kernel generally achieves 80%+ kernel occupancy but it is not the most efficient implementation.*
	 *
	 *	@param sliver_array		Warped sliver image.
	 *	@param base_array		Warped base image.
	 *	@param sums_array		Storage area for sums.
	 *	@param partial_sum		Local area to store partial sums.
	 *	@param partial_error	Local area to store partial sum errors.
	 *	@param width			Sliver width.
	 *	@param minY				The minimum possible y-value (pre-computed) given parameters.
	 *	@param m				Precomputed value for matrix. Represents maximum possible x-values.
	 *	@param n				Precomputed value for matrix. Represents maximum possible y-values.
	 *	@param num_sub_blocks	The stride that each individual work unit makes.
	 *	@param num_blocks		The amount of column sums each work group computes.
	 */
	__kernel void clColumnSums( __global double* sliver_array, __global double* base_array, __global double* sums_array, __local double* partial_sum, __local double* partial_errors,
								int width, int minY, int m, int n, int num_sub_blocks, int num_blocks )	{
		__private int local_id = get_local_id( 0 );
		__private int temp = local_id / num_sub_blocks;
		__private int sub_id = local_id - temp * num_sub_blocks;		// Inline modulo since OpenCL has modulo cost.
		__private int i = ( get_group_id( 0 ) * num_blocks ) + temp > m - 1 ? m - 1 : ( get_group_id( 0 ) * num_blocks ) + temp;	// Compute kC0 as well.
		
		double sum = 0;
		double sum_error = 0;
		temp = width * minY + i;
		for ( uint j = sub_id; j < n; j += num_sub_blocks )	{
			double diff = base_array[ width * j + temp ] - sliver_array[ width * j + temp ];	// [ i + width * ( j + minY ) ]
			double y = diff - sum_error;
			double t = sum + y;
			sum_error = ( t - sum ) - y;
			sum = t;
		}

		partial_sum[ local_id ] = sum;
		partial_errors[ local_id ] = sum_error;
		barrier( CLK_LOCAL_MEM_FENCE );

		if ( sub_id == 0 ) {
			sum = 0;
			sum_error = 0;

			for ( uint j = 0; j < num_sub_blocks; ++j )	{
				sum_error += partial_errors[ local_id + j ];
			}
			for ( uint j = 0; j < num_sub_blocks; ++j )	{
				double diff = partial_sum[ local_id + j ];
				double y = diff - sum_error;
				double t = sum + y;
				sum_error = ( t - sum ) - y;
				sum = t;
			}
			sums_array[ i ] = sum;
		}
	}

	/*
	 *	Performs the sum of differences of values within each individual row. We fix a y-value
	 *	and compute kR_y = sum_{i=0}^{m} (sliver(i, y) - base(i, y)). Each individual work group computes
	 *	num_blocks kR values. We use Kahan Summation implementation to achieve reproducibility in results.
	 *
	 *	*Note: This kernel generally achieves 80%+ kernel occupancy but it is not the most efficient implementation.*
	 *
	 *	@param sliver_array		Warped sliver image.
	 *	@param base_array		Warped base image.
	 *	@param sums_array		Storage area for sums.
	 *	@param partial_sum		Local area to store partial sums.
	 *	@param partial_error	Local area to store partial sum errors.
	 *	@param width			Sliver width.
	 *	@param minY				The minimum possible y-value (pre-computed) given parameters.
	 *	@param m				Precomputed value for matrix. Represents maximum possible x-values.
	 *	@param n				Precomputed value for matrix. Represents maximum possible y-values.
	 *	@param num_sub_blocks	The stride that each individual work unit makes.
	 *	@param num_blocks		The amount of column sums each work group computes.
	 */
	__kernel void clRowSums( __global double* sliver_array, __global double* base_array, __global double* sums_array, __local double* partial_sum, __local double* partial_errors,
							int width, int minY, int m, int n, int num_sub_blocks, int num_blocks )	{
		__private int local_id = get_local_id( 0 );
		__private int temp = local_id / num_sub_blocks;
		__private int sub_id = local_id - temp * num_sub_blocks;		// Inline modulo since OpenCL has modulo cost.
		__private int j = ( get_group_id( 0 ) * num_blocks ) + temp > n - 1 ? n - 1 : ( get_group_id( 0 ) * num_blocks ) + temp;

		double sum = 0;
		double sum_error = 0;
		temp = width * j + width * minY;
		for ( uint i = sub_id; i < m; i += num_sub_blocks )	{
			double diff = sliver_array[ i + temp ] - base_array[ i + temp ];	// [ i + width * ( j + minY ) ]
			double y = diff - sum_error;
			double t = sum + y;
			sum_error = ( t - sum ) - y;
			sum = t;
		}
	
		partial_sum[ local_id ] = sum;
		partial_errors[ local_id ] = sum_error;
		barrier( CLK_LOCAL_MEM_FENCE );
	
		if ( sub_id == 0 ) {
			sum = 0;
			sum_error = 0;

			for ( uint i = 0; i < num_sub_blocks; ++i )	{
				sum_error += partial_errors[ local_id + i ];
			}

			for ( uint i = 0; i < num_sub_blocks; ++i )	{
				double diff = partial_sum[ local_id + i ];
				double y = diff - sum_error;
				double t = sum + y;
				sum_error = ( t - sum ) - y;
				sum = t;
			}
			sums_array[ width + 1 + j ] = sum;
		}
	}

	/*
	 *	Performs the closed form expression calculation for r and c values.
	 *	Where we have that r_j = 1/n * sum_{i=1}^{m-1}(kC_i) + (m-1)/(mn) * sum_{i=0}^{n-1}(kR_i) + n/(mn) * kR_j
	 *	and that c_i = 1/n * (-kC_0+kC_i)
	 *	We use Kahan Summation implementation to achieve reproducibility in results.
	 *
	 *	*Note: This kernel generally achieves 80%+ kernel occupancy but it is not the most efficient implementation.*
	 *
	 *	@param sums_array		Storage area for sums.
	 *	@param matrix_result	Storage area for r and c computation.
	 *	@param partial_sum		Local area to store partial sums.
	 *	@param partial_error	Local area to store partial sum errors.
	 *	@param width			Sliver width.
	 *	@param m				Precomputed value for matrix. Represents maximum possible x-values.
	 *	@param n				Precomputed value for matrix. Represents maximum possible y-values.
	 */
	__kernel void clRCResults( __global double* sums_array, __global double* matrix_result,
							__local double* partial_sum, __local double* partial_errors, int width, int m, int n )	{
		__private int local_id = get_local_id( 0 );
		__private int index = get_global_id( 0 );
		__private double sum = 0;
		__private double skr = 0;
		__private double sum_error = 0;
		__private double skr_error = 0;
		__private double inv_n = 1 / ( double ) n;
		__private double inv_mn = 1 / ( double ) ( m * n );
		__private double inv_m = 1 / ( double ) m;

		for( uint i = index; i < m - 1; i += get_global_size( 0 ) )	{
			// c_i = 1/n * ( -kc_0 + kc_i )
			matrix_result[ i ] = ( inv_n ) * ( sums_array[ i + 1 ] - sums_array[ 0 ] );
		}
		for( uint i = get_local_id( 0 ); i < n; i += get_local_size( 0 ) )	{
			double val = sums_array[ width + 1 + i ];
			double y = val - sum_error;
			double t = sum + y;
			sum_error = ( t - sum ) - y;
			sum = t;
		}
		partial_sum[ local_id ] = sum;
		partial_errors[ local_id ] = sum_error;
		barrier( CLK_LOCAL_MEM_FENCE );

		skr_error = 0;
		for ( uint i = 0; i < get_local_size( 0 ); ++i )	{
			skr_error += partial_errors[ i ];
		}
		for( uint i = 0; i < get_local_size( 0 ); ++i )	{
			double val = partial_sum[ i ];
			double y = val - skr_error;
			double t = skr + y;
			skr_error = ( t - skr ) - y;
			skr = t;
		}
		for( uint j = index; j < n; j += get_global_size( 0 ) )	{
			// r_i = -1/n * kc_0 - 1/mn * skr + 1/m * kr_i
			matrix_result[ width + j ] = ( -inv_n ) * sums_array[ 0 ] - ( inv_mn ) * skr + inv_m * sums_array[ width + 1 + j ];
		}

	}

	__kernel void clPartialDifference( __global double* partial_sums, __global double* partial_weights,
											__global double* partial_sums_errors, __global double* partial_weights_errors,
											__global double* matrix_results,
											__global double* base_array, __global double* sliver_array,
											__local double* local_sums, __local double* local_weights,
											__local double* local_sums_errors, __local double* local_weights_errors,
											int maxX, int minY, int maxY, int width, double weightRight, double weightTop, double weightBottom )	{
		
		int local_id = get_local_id( 0 );
		int local_size = get_local_size( 0 );
		int i = get_group_id( 0 );
		double c_change = ( i == 0 ) ? 0 : matrix_results[ i - 1 ];

		__private double sum = 0;
		__private double area = 0;
		__private double sum_error = 0;
		__private double area_error = 0;

		for( int j = local_id + minY; j < maxY + 1; j += local_size )	{
			int index = i + width * j;
			double r_change = matrix_results[ width + j - minY ];

			double weight = ( i == maxX ) ? weightRight : 1;
			weight *= ( j == maxY ) ? weightTop : 1;
			weight *= ( j == minY ) ? weightBottom : 1;


			double diff = ( base_array[ index ] + r_change ) - ( sliver_array[ index ] + c_change );	// [ i + width * ( j + minY ) ]
			double y = diff * diff * weight - sum_error;
			double t = sum + y;
			sum_error = ( t - sum ) - y;
			sum = t;

			y = weight - area_error;
			t = area + y;
			area_error = ( t - area ) - y;
			area = t;
		}

		local_sums[ local_id ] = sum;
		local_sums_errors[ local_id ] = sum_error;
		local_weights[ local_id ] = area;
		local_weights_errors[ local_id ] = area_error;
		barrier( CLK_LOCAL_MEM_FENCE );
	
		if ( local_id == 0 ) {
			sum = 0;
			area = 0;
			sum_error = 0;
			area_error = 0;

			for ( uint k = 0; k < local_size; ++k )	{
				sum_error += local_sums_errors[ k ];
				area_error += local_weights_errors[ k ];
			}

			for ( uint k = 0; k < local_size; ++k )	{
				double diff = local_sums[ k ];
				double y = diff - sum_error;
				double t = sum + y;
				sum_error = ( t - sum ) - y;
				sum = t;

				double weight = local_weights[ k ];
				y = weight - area_error;
				t = area + y;
				area_error = ( t - area ) - y;
				area = t;
			}
			partial_sums[ i ] = sum;
			partial_sums_errors[ i ] = sum_error;
			partial_weights[ i ] = area;
			partial_weights_errors[ i ] = area_error;
		}
	}


	__kernel void clFinalDifference( __global double* partial_sums, __global double* partial_weights,
									__global double* partial_sums_errors, __global double* partial_weights_errors,
									__global double* final_result,
									__local double* local_sums, __local double* local_weights,
									__local double* local_sums_errors, __local double* local_weights_errors,
									int maxX, int minY, int maxY )	{

		__private int local_id = get_local_id( 0 );
		__private double sum = 0;
		__private double weight = 0;
		__private double sum_error = 0;
		__private double weight_error = 0;

		for( uint i = local_id; i < maxX + 1; i += get_local_size( 0 ) )	{
			double val = partial_sums[ i ];
			double y = val - sum_error;
			double t = sum + y;
			sum_error = ( t - sum ) - y;
			sum = t;

			val = partial_weights[ i ];
			y = val - weight_error;
			t = weight + y;
			weight_error = ( t - weight ) - y;
			weight = t;
		}
		local_sums[ local_id ] = sum;
		local_weights[ local_id ] = weight;
		local_sums_errors[ local_id ] = sum_error;
		local_weights_errors[ local_id ] = weight_error;
		barrier( CLK_LOCAL_MEM_FENCE );

		if( local_id == 0 )	{
			sum = 0;
			sum_error = 0;
			weight = 0;
			weight_error = 0;

			for ( uint i = 0; i < maxX + 1; ++i )	{
				sum_error += partial_sums_errors[ i ];
				weight_error += partial_weights_errors[ i ];
			}
			for ( uint i = 0; i < get_local_size( 0 ); ++i )	{
				sum_error += local_sums_errors[ i ];
				weight_error += local_weights_errors[ i ];
			}
			for( uint i = 0; i < get_local_size( 0 ); ++i )	{
				double y = local_sums[ i ] - sum_error;
				double t = sum + y;
				sum_error = ( t - sum ) - y;
				sum = t;

				y = local_weights[ i ] - weight_error;
				t = weight + y;
				weight_error = ( t - weight ) - y;
				weight = t;
			}
			final_result[ 0 ] = sum / weight;
		}
	}


	__kernel void clDifferenceCalculation( __global double* final_result, __global double* matrix_results,
											__global double* base_array, __global double* sliver_array,
											__local double* partial_sums, __local double* partial_weights,
											__local double* partial_sums_errors, __local double* partial_weights_errors,
											int maxX, int minY, int maxY, int width, double weightRight, double weightTop, double weightBottom )	{
		int local_id = get_local_id( 0 );
		int work_size = get_local_size( 0 );

		double sum = 0;
		double area = 0;
		double sum_error = 0;
		double area_error = 0;

		for( int j = local_id + minY; j < maxY + 1; j += work_size )	{
			double r_change = matrix_results[ width + j - minY ];
			for( int i = 0; i < maxX + 1; i++ )	{
				int index = i + width * j;
				double c_change = ( i == 0 ) ? 0 : matrix_results[ i - 1 ];
				double weight = ( i == maxX ) ? weightRight : 1;
				weight *= ( j == maxY ) ? weightTop : 1;
				weight *= ( j == minY ) ? weightBottom : 1;

				double diff = ( base_array[ index ] + r_change ) - ( sliver_array[ index ] + c_change );

				double y = diff * diff * weight - sum_error;
				double t = sum + y;
				sum_error = ( t - sum ) - y;
				sum = t;

				y = weight - area_error;
				t = area + y;
				area_error = ( t - area ) - y;
				area = t;
			}
		}

		partial_sums[ local_id ] = sum;
		partial_sums_errors[ local_id ] = sum_error;
		partial_weights[ local_id ] = area;
		partial_weights_errors[ local_id ] = area_error;
		barrier(CLK_LOCAL_MEM_FENCE);

		if( local_id == 0 )	{
			double differences = 0;
			double areas = 0;

			double differences_error = 0;
			double areas_error = 0;

			for( int i = 0; i < work_size; ++i )	{
				differences_error += partial_sums_errors[ i ];
				areas_error += partial_weights_errors[ i ];
			}

			for( int i = 0; i < work_size; ++i )	{
				double y = partial_sums[ i ] - differences_error;
				double t = differences + y;
				differences_error = ( t - differences ) - y;
				differences = t;

				y = partial_weights[ i ] - areas_error;
				t = areas + y;
				areas_error = ( t - areas ) - y;
				areas = t;
			}
			final_result[ 0 ] = differences / areas;
		}
	}

	)";

	#endif

#endif __KERNELS__