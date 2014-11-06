#pragma OPENCL EXTENSION cl_amd_printf : enable

__kernel void clWarpSliverCubic( __global __const double* sliver, __global double* sliver_array, double A1, double B1,
						int width, int height, int minY, int maxY, __local double* sx, __local double* sy, __local double* ux, __local double* uy )	{
	__local double A;
	__local int local_size;
	__local int iterations;
	__local double yoffset;
	__local double origX;
	__local int x;
	__private int y = get_local_id( 0 );

	if( y == 0 )	{
		x = get_group_id( 0 );
		A = -0.5;
		local_size = get_local_size(0);
		iterations = ceil((double)( maxY - minY + 1 ) / (double) get_local_size(0));

		yoffset = -( B1 - 1 ) * x;
		origX = x - A1 * x;

		sx[ 0 ] = fabs( origX - floor( origX - 1 ) ); // these get the distances for each row and column from the initial point, positive
		sy[ 0 ] = fabs( yoffset - floor( yoffset - 1 ) );
		sx[ 1 ] = fabs( 1 - sx[ 0 ] );
		sy[ 1 ] = fabs( 1 - sy[ 0 ] );
		sx[ 2 ] = fabs( 1 - sx[ 1 ] );
		sy[ 2 ] = fabs( 1 - sy[ 1 ] );
		sx[ 3 ] = 1 + sx[ 2 ];
		sy[ 3 ] = 1 + sy[ 2 ];

		for ( int j = 0; j < 4; ++j ) { // these, dependent upon the distance, implement the polynomial weighting
			double x2 = sx[ j ] * sx[ j ];
			double x3 = x2 * sx[ j ];
			double y2 = sy[ j ] * sy[ j ];
			double y3 = y2 * sy[ j ];
			if ( sx[ j ] <= 1 && sx[ j ] >= 0 ) {
				ux[ j ] = ( A + 2 ) * x3 - ( A + 3 ) * x2 + 1;
			}
			else if ( sx[ j ] > 1 && sx[ j ] <= 2 ) {
				ux[ j ] = A * x3 - 5 * A * x2 + 8 * A * sx[ j ] - 4 * A;
			}
			else {
				ux[ j ] = 0;
			}
			if ( sy[ j ] <= 1 && sy[ j ] >= 0 ) {
				uy[ j ] = ( A + 2 ) * y3 - ( A + 3 ) * y2 + 1;
			}
			else if ( sy[ j ] > 1 && sy[ j ] <= 2 ) {
				uy[ j ] = A * y3 - 5 * A * y2 + 8 * A * sy[ j ] - 4 * A;
			}
			else {
				uy[ j ] = 0;
			}
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	for( __private int it = 0; it < iterations; ++it )	{
		__private double pixelvalue = 0.0;
		__private double origY = y * iterations + it + minY + yoffset;
		if( y * iterations + it + minY < maxY + 1 )	{
			if ( (int)( origX ) == 0 || (int)( origY ) == 0 || ceil( origX ) == width - 1 || ceil( origY ) == height - 1 ) {
				__private double s, t;
				__private double left_val, right_val;
				__private double top_val, bottom_val;

				// Get integer coordinates for top left corner
				__private int left_index = ( int ) origX;
				__private int top_index = ( int ) origY;
				__private int right_index = ceil( origX );
				__private int bottom_index = ceil( origY );

				// Interpolate across the top edge
				s = origX - left_index;
				t = origY - top_index;
				left_val = sliver[ left_index + width * top_index ];
				right_val = sliver[ right_index + width * top_index ];
				// Linear interpolation recoded to only one multiply
				top_val = s * right_val + ( 1 - s ) * ( left_val );

				// Interpolate across the bottom edge
				left_val = sliver[ left_index + width * bottom_index ];
				right_val = sliver[ right_index + width * bottom_index ];
				bottom_val = s * right_val + ( 1 - s ) * ( left_val );

				// Interpolate between top and bottom
				pixelvalue = ( t * bottom_val + ( 1 - t ) * top_val );
			}
			else {
				for ( int i = 0; i < 4; ++i ) {
					for ( int j = 0; j < 4; ++j ) {
						pixelvalue += sliver[ (int)( origX ) - 1 + i + width * ( (int)( origY ) - 1 + j ) ] * ux[ i ] * uy[ j ];
					}
				}
			}
			sliver_array[ x + width * ( y * iterations + it + minY ) ] = pixelvalue;
		}
	}
}

__kernel void clWarpBaseCubic( __global __const double* base, __global double* base_array, double A0, double A1, double A2, double A3, double B0, double B1, double B2, double B3,
						int width, int height, int bwidth, int xbp, int maxX, int minY, __local double* sx, __local double* sy, __local double* ux, __local double* uy )	{
	__local double A;
	__local int local_size;
	__local int iterations;
	__local double xoffset;
	__local double origY;
	__local int y;
	__private int x = get_local_id( 0 );

	if( x == 0 )	{
		y = get_group_id( 0 ) + minY;
		A = -0.5;
		local_size = get_local_size(0);
		iterations = ceil((double)( maxX + 1 ) / (double) get_local_size(0));

		__private int ym2 = y * y;
		__private int ym3 = ym2 * y;
		xoffset = A0 + A1 * y + A2 * ym2 + A3 * ym3;
		origY = B0 + B1 * y + B2 * ym2 + B3 * ym3;

		sx[ 0 ] = fabs( xoffset - floor( xoffset - 1 ) ); // these get the distances for each row and column from the initial point, positive
		sy[ 0 ] = fabs( origY - floor( origY - 1 ) );
		sx[ 1 ] = fabs( 1 - sx[ 0 ] );
		sy[ 1 ] = fabs( 1 - sy[ 0 ] );
		sx[ 2 ] = fabs( 1 - sx[ 1 ] );
		sy[ 2 ] = fabs( 1 - sy[ 1 ] );
		sx[ 3 ] = 1 + sx[ 2 ];
		sy[ 3 ] = 1 + sy[ 2 ];

		for ( int j = 0; j < 4; ++j ) { // these, dependent upon the distance, implement the polynomial weighting
			__private double x2 = sx[ j ] * sx[ j ];
			__private double x3 = x2 * sx[ j ];
			__private double y2 = sy[ j ] * sy[ j ];
			__private double y3 = y2 * sy[ j ];
			if ( sx[ j ] <= 1 && sx[ j ] >= 0 ) {
				ux[ j ] = ( A + 2 ) * x3 - ( A + 3 ) * x2 + 1;
			}
			else if ( sx[ j ] > 1 && sx[ j ] <= 2 ) {
				ux[ j ] = A * x3 - 5 * A * x2 + 8 * A * sx[ j ] - 4 * A;
			}
			else {
				ux[ j ] = 0;
			}
			if ( sy[ j ] <= 1 && sy[ j ] >= 0 ) {
				uy[ j ] = ( A + 2 ) * y3 - ( A + 3 ) * y2 + 1;
			}
			else if ( sy[ j ] > 1 && sy[ j ] <= 2 ) {
				uy[ j ] = A * y3 - 5 * A * y2 + 8 * A * sy[ j ] - 4 * A;
			}
			else {
				uy[ j ] = 0;
			}
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	for( __private int it = 0; it < iterations; ++it )	{
		__private double pixelvalue = 0.0;
		__private double origX = x * iterations + it + xbp + xoffset;
		if( x * iterations + it < maxX + 1 )	{
			if ( (int)( origX ) == 0 || (int)( origY ) == 0 || ceil( origX ) == bwidth - 1 || ceil( origY ) == height - 1 ) {
				__private double s, t;
				__private double left_val, right_val;
				__private double top_val, bottom_val;

				// Get integer coordinates for top left corner
				__private int left_index = ( int ) origX;
				__private int top_index = ( int ) origY;
				__private int right_index = ceil( origX );
				__private int bottom_index = ceil( origY );

				// Interpolate across the top edge
				s = origX - left_index;
				t = origY - top_index;
				left_val = base[ left_index + bwidth * top_index ];
				right_val = base[ right_index + bwidth * top_index ];
				// Linear interpolation recoded to only one multiply
				top_val = s * right_val + ( 1 - s ) * ( left_val );

				// Interpolate across the bottom edge
				left_val = base[ left_index + bwidth * bottom_index ];
				right_val = base[ right_index + bwidth * bottom_index ];
				bottom_val = s * right_val + ( 1 - s ) * ( left_val );

				// Interpolate between top and bottom
				pixelvalue = ( t * bottom_val + ( 1 - t ) * top_val );
			}
			else {
				for ( int i = 0; i < 4; ++i ) {
					for ( int j = 0; j < 4; ++j ) {
						pixelvalue += base[ (int)( origX ) - 1 + i + bwidth * ( (int)( origY ) - 1 + j ) ] * ux[ i ] * uy[ j ];
					}
				}
			}
			base_array[ x * iterations + it + width * y ] = pixelvalue;
		}
	}
}



__kernel void clWarpSliverBspline( __global __const double* sliver, __global double* sliver_array, double A1, double B1,
						int width, int height, int minY, int maxY, __local double* sx, __local double* sy, __local double* ux, __local double* uy )	{
	__local int B;
	__local int C;
	__local int local_size;
	__local int iterations;
	__local double yoffset;
	__local double origX;
	__local int x;
	__private int y = get_local_id( 0 );

	if( y == 0 )	{
		x = get_group_id( 0 );
		B = 1;
		C = 0;
		local_size = get_local_size(0);
		iterations = ceil((double)( maxY - minY + 1 ) / (double) get_local_size(0));

		yoffset = -( B1 - 1 ) * x;
		origX = x - A1 * x;

		sx[ 0 ] = fabs( origX - floor( origX - 1 ) ); // these get the distances for each row and column from the initial point, positive
		sy[ 0 ] = fabs( yoffset - floor( yoffset - 1 ) );
		sx[ 1 ] = fabs( 1 - sx[ 0 ] );
		sy[ 1 ] = fabs( 1 - sy[ 0 ] );
		sx[ 2 ] = fabs( 1 - sx[ 1 ] );
		sy[ 2 ] = fabs( 1 - sy[ 1 ] );
		sx[ 3 ] = 1 + sx[ 2 ];
		sy[ 3 ] = 1 + sy[ 2 ];

		for ( int j = 0; j < 4; ++j ) { // these, dependent upon the distance, implement the polynomial weighting
			double x2 = sx[ j ] * sx[ j ];
			double x3 = x2 * sx[ j ];
			double y2 = sy[ j ] * sy[ j ];
			double y3 = y2 * sy[ j ];
			if ( sx[ j ] <= 1 && sx[ j ] >= 0 ) {
				ux[ j ] = ( ( 12 - ( 9 * B ) - ( 6 * C ) ) * x3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*x2 + 6 - ( 2 * B ) ) / 6;
			}
			else if ( sx[ j ] > 1 && sx[ j ] <= 2 ) {
				ux[ j ] = ( ( -B - ( 6 * C ) )*x3 + ( ( 6 * B ) + ( 30 * C ) )*x2 + ( -( 12 * B ) - ( 48 * C ) )*( sx[ j ] ) + ( 8 * B ) + ( 24 * C ) ) / 6;
			}
			else	{
				ux[ j ] = 0;
			}
			if ( sy[ j ] <= 1 && sy[ j ] >= 0 ) {
				uy[ j ] = ( ( 12 - ( 9 * B ) - ( 6 * C ) )*y3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*y2 + 6 - ( 2 * B ) ) / 6;
			}
			else if ( sy[ j ] > 1 && sy[ j ] <= 2 ) {
				uy[ j ] = ( ( -B - ( 6 * C ) )*y3 + ( ( 6 * B ) + ( 30 * C ) )*y2 + ( -( 12 * B ) - ( 48 * C ) )*( sy[ j ] ) + ( 8 * B ) + ( 24 * C ) ) / 6;
			}
			else	{
				uy[ j ] = 0;
			}
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	for( __private int it = 0; it < iterations; ++it )	{
		__private double pixelvalue = 0;
		__private double origY = y * iterations + it + minY + yoffset;
		if( y * iterations + it + minY < maxY + 1 )	{
			if ( (int)( origX ) == 0 || (int)( origY ) == 0 || ceil( origX ) == width - 1 || ceil( origY ) == height - 1 ) {
				__private double s, t;
				__private double left_val, right_val;
				__private double top_val, bottom_val;

				// Get integer coordinates for top left corner
				__private int left_index = ( int ) origX;
				__private int top_index = ( int ) origY;
				__private int right_index = ceil( origX );
				__private int bottom_index = ceil( origY );

				// Interpolate across the top edge
				s = origX - left_index;
				t = origY - top_index;
				left_val = sliver[ left_index + width * top_index ];
				right_val = sliver[ right_index + width * top_index ];
				// Linear interpolation recoded to only one multiply
				top_val = s * right_val + ( 1 - s ) * ( left_val );

				// Interpolate across the bottom edge
				left_val = sliver[ left_index + width * bottom_index ];
				right_val = sliver[ right_index + width * bottom_index ];
				bottom_val = s * right_val + ( 1 - s ) * ( left_val );

				// Interpolate between top and bottom
				pixelvalue = ( t * bottom_val + ( 1 - t ) * top_val );
			}
			else {
				for ( int i = 0; i < 4; ++i ) {
					for ( int j = 0; j < 4; ++j ) {
						pixelvalue += sliver[ (int)( origX ) - 1 + i + width * ( (int)( origY ) - 1 + j ) ] * ux[ i ] * uy[ j ];
					}
				}
			}
			sliver_array[ x + width * ( y * iterations + it + minY ) ] = pixelvalue;
		}
	}
}

__kernel void clWarpBaseBspline( __global __const double* base, __global double* base_array, double A0, double A1, double A2, double A3, double B0, double B1, double B2, double B3,
						int width, int height, int bwidth, int xbp, int maxX, int minY, __local double* sx, __local double* sy, __local double* ux, __local double* uy )	{
	__local int B;
	__local int C;
	__local int local_size;
	__local int iterations;
	__local double xoffset;
	__local double origY;
	__local int y;
	__private int x = get_local_id( 0 );

	if( x == 0 )	{
		y = get_group_id( 0 ) + minY;
		B = 1;
		C = 0;
		local_size = get_local_size(0);
		iterations = ceil((double)( maxX + 1 ) / (double) get_local_size(0));

		__private int ym2 = y * y;
		__private int ym3 = ym2 * y;
		xoffset = A0 + A1 * y + A2 * ym2 + A3 * ym3;
		origY = B0 + B1 * y + B2 * ym2 + B3 * ym3;

		sx[ 0 ] = fabs( xoffset - floor( xoffset - 1 ) ); // these get the distances for each row and column from the initial point, positive
		sy[ 0 ] = fabs( origY - floor( origY - 1 ) );
		sx[ 1 ] = fabs( 1 - sx[ 0 ] );
		sy[ 1 ] = fabs( 1 - sy[ 0 ] );
		sx[ 2 ] = fabs( 1 - sx[ 1 ] );
		sy[ 2 ] = fabs( 1 - sy[ 1 ] );
		sx[ 3 ] = 1 + sx[ 2 ];
		sy[ 3 ] = 1 + sy[ 2 ];

		for ( int j = 0; j < 4; ++j ) { // these, dependent upon the distance, implement the polynomial weighting
			__private double x2 = sx[ j ] * sx[ j ];
			__private double x3 = x2 * sx[ j ];
			__private double y2 = sy[ j ] * sy[ j ];
			__private double y3 = y2 * sy[ j ];
			if ( sx[ j ] <= 1 && sx[ j ] >= 0 ) {
				ux[ j ] = ( ( 12 - ( 9 * B ) - ( 6 * C ) ) * x3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*x2 + 6 - ( 2 * B ) ) / 6;
			}
			else if ( sx[ j ] > 1 && sx[ j ] <= 2 ) {
				ux[ j ] = ( ( -B - ( 6 * C ) )*x3 + ( ( 6 * B ) + ( 30 * C ) )*x2 + ( -( 12 * B ) - ( 48 * C ) )*( sx[ j ] ) + ( 8 * B ) + ( 24 * C ) ) / 6;
			}
			else	{
				ux[ j ] = 0;
			}
			if ( sy[ j ] <= 1 && sy[ j ] >= 0 ) {
				uy[ j ] = ( ( 12 - ( 9 * B ) - ( 6 * C ) )*y3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*y2 + 6 - ( 2 * B ) ) / 6;
			}
			else if ( sy[ j ] > 1 && sy[ j ] <= 2 ) {
				uy[ j ] = ( ( -B - ( 6 * C ) )*y3 + ( ( 6 * B ) + ( 30 * C ) )*y2 + ( -( 12 * B ) - ( 48 * C ) )*( sy[ j ] ) + ( 8 * B ) + ( 24 * C ) ) / 6;
			}
			else	{
				uy[ j ] = 0;
			}
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	for( __private int it = 0; it < iterations; ++it )	{
		__private double pixelvalue = 0;
		__private double origX = x * iterations + it + xbp + xoffset;
		if( x * iterations + it < maxX + 1 )	{
			if ( (int)( origX ) == 0 || (int)( origY ) == 0 || ceil( origX ) == bwidth - 1 || ceil( origY ) == height - 1 ) {
				__private double s, t;
				__private double left_val, right_val;
				__private double top_val, bottom_val;

				// Get integer coordinates for top left corner
				__private int left_index = ( int ) origX;
				__private int top_index = ( int ) origY;
				__private int right_index = ceil( origX );
				__private int bottom_index = ceil( origY );

				// Interpolate across the top edge
				s = origX - left_index;
				t = origY - top_index;
				left_val = base[ left_index + bwidth * top_index ];
				right_val = base[ right_index + bwidth * top_index ];
				// Linear interpolation recoded to only one multiply
				top_val = s * right_val + ( 1 - s ) * ( left_val );

				// Interpolate across the bottom edge
				left_val = base[ left_index + bwidth * bottom_index ];
				right_val = base[ right_index + bwidth * bottom_index ];
				bottom_val = s * right_val + ( 1 - s ) * ( left_val );

				// Interpolate between top and bottom
				pixelvalue = ( t * bottom_val + ( 1 - t ) * top_val );
			}
			else {
				for ( int i = 0; i < 4; ++i ) {
					for ( int j = 0; j < 4; ++j ) {
						pixelvalue += base[ (int)( origX ) - 1 + i + bwidth * ( (int)( origY ) - 1 + j ) ] * ux[ i ] * uy[ j ];
					}
				}
			}
			base_array[ x * iterations + it + width * y ] = pixelvalue;
		}
	}
}



#ifdef cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif
#pragma OPENCL EXTENSION cl_amd_printf : enable
__kernel void clColumnSums( __global double* sliver_array, __global double* base_array, __global double* sums_array, __local double* partial_sum,
							int width, int minY, int m, int n, int num_sub_blocks, int num_blocks )	{
	__private int local_id = get_local_id( 0 );
	__private int temp = local_id / num_sub_blocks;
	__private int sub_id = local_id - temp * num_sub_blocks;		// Inline modulo since OpenCL has modulo cost.
	__private int i = ( get_group_id( 0 ) * num_blocks ) + temp > m - 1 ? m - 1 : ( get_group_id( 0 ) * num_blocks ) + temp;	// Compute kC0 as well.
	double sum = 0;
	temp = width * minY + i;
	for ( uint j = sub_id; j < n; j += num_sub_blocks )	{
		sum += base_array[ width * j + temp ] - sliver_array[ width * j + temp ];	// [ i + width * ( j + minY ) ]
	}

	partial_sum[ local_id ] = sum;
	barrier( CLK_LOCAL_MEM_FENCE );

	if ( sub_id == 0 ) {
		sum = 0;
		for ( uint t = 0; t < num_sub_blocks; ++t )
			sum += partial_sum[ local_id + t ];
		sums_array[ i ] = sum;
	}
}

__kernel void clRowSums( __global double* sliver_array, __global double* base_array, __global double* sums_array, __local double* partial_sum,
						int width, int minY, int m, int n, int num_sub_blocks, int num_blocks )	{
	__private int local_id = get_local_id( 0 );
	__private int temp = local_id / num_sub_blocks;
	__private int sub_id = local_id - temp * num_sub_blocks;		// Inline modulo since OpenCL has modulo cost.
	__private int j = ( get_group_id( 0 ) * num_blocks ) + temp > n - 1 ? n - 1 : ( get_group_id( 0 ) * num_blocks ) + temp;

	double sum = 0;
	temp = width * j + width * minY;
	for ( uint i = sub_id; i < m; i += num_sub_blocks )	{
		sum += sliver_array[ i + temp ] - base_array[ i + temp ];	// [ i + width * ( j + minY ) ]
	}
	
	partial_sum[ local_id ] = sum;
	barrier( CLK_LOCAL_MEM_FENCE );
	
	if ( sub_id == 0 ) {
		sum = 0;
		for ( uint t = 0; t < num_sub_blocks; ++t )
			sum += partial_sum[ local_id + t ];
		sums_array[ width + 1 + j ] = sum;
	}
}

__kernel void clRCResults( __global double* sums_array, __global double* matrix_result, __local double* partial_sum, int width, int m, int n )	{
	__private int local_id = get_local_id( 0 );
	__private int index = get_global_id( 0 );
	__private double sum = 0;

	for( uint i = index; i < m - 1; i += get_global_size( 0 ) )	{
		// c_i = 1/n * ( -kc_0 + kc_i )
		matrix_result[ i ] = ( 1 / ( double ) n ) * ( sums_array[ i + 1 ] - sums_array[ 0 ] );
	}
	for( uint i = get_local_id( 0 ); i < n; i += get_local_size( 0 ) )	{
		sum += sums_array[ width + 1 + i ];
	}
	partial_sum[ local_id ] = sum;
	barrier( CLK_LOCAL_MEM_FENCE );
	
	double skr;
	skr = 0;
	for( uint i = 0; i < get_local_size( 0 ); ++i )	{
		skr += partial_sum[ i ];
	}
	for( uint j = index; j < n; j += get_global_size( 0 ) )	{
		// r_i = -1/n * kc_0 - 1/mn * skr + 1/m * kr_i
		matrix_result[ width + j ] = ( -1 / ( double ) n ) * sums_array[ 0 ] - ( 1 / ( double ) ( m * n ) ) * skr + 1 / ( double ) m * sums_array[ width + 1 + j ];
	}

}


__kernel void clDifferenceCalculation( __global double* final_result, __global double* matrix_results,
										__global double* base_array, __global double* sliver_array,
										__local double* partial_sums, __local double* partial_weights,
										int problem_width, int problem_height,
										int maxX, int minY, int maxY, int width, double weightRight, double weightTop, double weightBottom )	{
	int local_id = get_local_id( 0 );
	int work_size = get_local_size( 0 );

	double sum = 0;
	double area = 0;

	for( int j = local_id + minY; j < maxY + 1; j += work_size )	{

		double r_change = matrix_results[ width + j - minY ];

		for( int i = 0; i < maxX + 1; i++ )	{

			double c_change = ( i == 0 ) ? 0 : matrix_results[ i - 1 ];

			double weight = ( i == maxX ) ? weightRight : 1;
			weight *= ( j == maxY ) ? weightTop : 1;
			weight *= ( j == minY ) ? weightBottom : 1;

			int index = i + width * j;

			double diff = ( base_array[ index ] + r_change ) - ( sliver_array[ index ] + c_change );
			sum += diff * diff * weight;
			area += weight;

		}
	}

	partial_sums[ local_id ] = sum;
	partial_weights[ local_id ] = area;
	barrier(CLK_LOCAL_MEM_FENCE);

	if( local_id == 0 )	{
		double differences = 0;
		double areas = 0;
		for( int i = 0; i < work_size; ++i )	{
			differences += partial_sums[ i ];
			areas += partial_weights[ i ];
		}
		final_result[ 0 ] = differences / areas;
	}
}
