#pragma OPENCL EXTENSION cl_amd_printf : enable

__kernel void clWarpSliverCubic( __global __const float* sliver, __global float* sliver_array, float A1, float B1,
						int width, int height, int minY, int maxY, __local float* sx, __local float* sy, __local float* ux, __local float* uy )	{
	__local float A;
	__local int local_size;
	__local int iterations;
	__local float yoffset;
	__local float origX;
	__local int x;
	__private int y = get_local_id( 0 );

	if( y == 0 )	{
		x = get_group_id( 0 );
		A = -0.5f;
		local_size = get_local_size( 0 );
		iterations = ceil( ( float ) ( maxY - minY + 1 ) / ( float ) get_local_size( 0 ) );

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
			float x2 = sx[ j ] * sx[ j ];
			float x3 = x2 * sx[ j ];
			float y2 = sy[ j ] * sy[ j ];
			float y3 = y2 * sy[ j ];
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
		__private float pixelvalue = 0.0f;
		__private float origY = y * iterations + it + minY + yoffset;
		if( y * iterations + it + minY < maxY + 1 )	{
			if ( ( int ) ( origX ) == 0 || ( int ) ( origY ) == 0 || ceil( origX ) == width - 1 || ceil( origY ) == height - 1 ) {
				__private float s, t;
				__private float left_val, right_val;
				__private float top_val, bottom_val;

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
						pixelvalue += sliver[ ( int ) ( origX ) - 1 + i + width * ( ( int ) ( origY ) - 1 + j ) ] * ux[ i ] * uy[ j ];
					}
				}
			}
			sliver_array[ x + width * ( y * iterations + it + minY ) ] = pixelvalue;
		}
	}
}

__kernel void clWarpBaseCubic( __global __const float* base, __global float* base_array, float A0, float A1, float A2, float A3, float B0, float B1, float B2, float B3,
						int width, int height, int bwidth, int xbp, int maxX, int minY, __local float* sx, __local float* sy, __local float* ux, __local float* uy )	{
	__local float A;
	__local int local_size;
	__local int iterations;
	__local float xoffset;
	__local float origY;
	__local int y;
	__private int x = get_local_id( 0 );

	if( x == 0 )	{
		y = get_group_id( 0 ) + minY;
		A = -0.5f;
		local_size = get_local_size( 0 );
		iterations = ceil( ( float ) ( maxX + 1 ) / ( float ) get_local_size( 0 ) );

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
			__private float x2 = sx[ j ] * sx[ j ];
			__private float x3 = x2 * sx[ j ];
			__private float y2 = sy[ j ] * sy[ j ];
			__private float y3 = y2 * sy[ j ];
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
		__private float pixelvalue = 0.0f;
		__private float origX = x * iterations + it + xbp + xoffset;
		if( x * iterations + it < maxX + 1 )	{
			if ( ( int ) ( origX ) == 0 || ( int ) ( origY ) == 0 || ceil( origX ) == bwidth - 1 || ceil( origY ) == height - 1 ) {
				__private float s, t;
				__private float left_val, right_val;
				__private float top_val, bottom_val;

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
						pixelvalue += base[ ( int ) ( origX ) - 1 + i + bwidth * ( ( int ) ( origY ) - 1 + j ) ] * ux[ i ] * uy[ j ];
					}
				}
			}
			base_array[ x * iterations + it + width * y ] = pixelvalue;
		}
	}
}



__kernel void clWarpSliverBspline( __global __const float* sliver, __global float* sliver_array, float A1, float B1,
						int width, int height, int minY, int maxY, __local float* sx, __local float* sy, __local float* ux, __local float* uy )	{
	__local int B;
	__local int C;
	__local int local_size;
	__local int iterations;
	__local float yoffset;
	__local float origX;
	__local int x;
	__private int y = get_local_id( 0 );

	if( y == 0 )	{
		x = get_group_id( 0 );
		B = 1;
		C = 0;
		local_size = get_local_size(0);
		iterations = ceil( ( float ) ( maxY - minY + 1 ) / ( float ) get_local_size( 0 ) );

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
			float x2 = sx[ j ] * sx[ j ];
			float x3 = x2 * sx[ j ];
			float y2 = sy[ j ] * sy[ j ];
			float y3 = y2 * sy[ j ];
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
	barrier( CLK_LOCAL_MEM_FENCE );

	for( __private int it = 0; it < iterations; ++it )	{
		__private float pixelvalue = 0.0f;
		__private float origY = y * iterations + it + minY + yoffset;
		if( y * iterations + it + minY < maxY + 1 )	{
			if ( ( int ) ( origX ) == 0 || ( int ) ( origY ) == 0 || ceil( origX ) == width - 1 || ceil( origY ) == height - 1 ) {
				__private float s, t;
				__private float left_val, right_val;
				__private float top_val, bottom_val;

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
						pixelvalue += sliver[ ( int ) ( origX ) - 1 + i + width * ( ( int ) ( origY ) - 1 + j ) ] * ux[ i ] * uy[ j ];
					}
				}
			}
			sliver_array[ x + width * ( y * iterations + it + minY ) ] = pixelvalue;
		}
	}
}

__kernel void clWarpBaseBspline( __global __const float* base, __global float* base_array, float A0, float A1, float A2, float A3, float B0, float B1, float B2, float B3,
						int width, int height, int bwidth, int xbp, int maxX, int minY, __local float* sx, __local float* sy, __local float* ux, __local float* uy )	{
	__local int B;
	__local int C;
	__local int local_size;
	__local int iterations;
	__local float xoffset;
	__local float origY;
	__local int y;
	__private int x = get_local_id( 0 );

	if( x == 0 )	{
		y = get_group_id( 0 ) + minY;
		B = 1;
		C = 0;
		local_size = get_local_size( 0 );
		iterations = ceil( ( float ) ( maxX + 1 ) / ( float ) get_local_size( 0 ) );

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
			__private float x2 = sx[ j ] * sx[ j ];
			__private float x3 = x2 * sx[ j ];
			__private float y2 = sy[ j ] * sy[ j ];
			__private float y3 = y2 * sy[ j ];
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
	barrier( CLK_LOCAL_MEM_FENCE );

	for( __private int it = 0; it < iterations; ++it )	{
		__private float pixelvalue = 0.0f;
		__private float origX = x * iterations + it + xbp + xoffset;
		if( x * iterations + it < maxX + 1 )	{
			if ( ( int ) ( origX ) == 0 || ( int ) ( origY ) == 0 || ceil( origX ) == bwidth - 1 || ceil( origY ) == height - 1 ) {
				__private float s, t;
				__private float left_val, right_val;
				__private float top_val, bottom_val;

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
						pixelvalue += base[ ( int ) ( origX ) - 1 + i + bwidth * ( ( int ) ( origY ) - 1 + j ) ] * ux[ i ] * uy[ j ];
					}
				}
			}
			base_array[ x * iterations + it + width * y ] = pixelvalue;
		}
	}
}

__kernel void clColumnSums( __global float* sliver_array, __global float* base_array, __global float* sums_array, __local float* partial_sum,
							int width, int minY, int m, int n, int num_sub_blocks, int num_blocks )	{
	__private int local_id = get_local_id( 0 );
	__private int temp = local_id / num_sub_blocks;
	__private int sub_id = local_id - temp * num_sub_blocks;		// Inline modulo since OpenCL has modulo cost.
	__private int i = ( get_group_id( 0 ) * num_blocks ) + 1 + temp;
	if( i > m - 1 )	{
		i = m - 1;
	}
	float sum = 0;
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
		sums_array[ i - 1 ] = sum;
	}
}

__kernel void clRowSums( __global float* sliver_array, __global float* base_array, __global float* sums_array, __local float* partial_sum,
						int width, int minY, int m, int n, int num_sub_blocks, int num_blocks )	{
	__private int local_id = get_local_id( 0 );
	__private int temp = local_id / num_sub_blocks;
	__private int sub_id = local_id - temp * num_sub_blocks;		// Inline modulo since OpenCL has modulo cost.
	__private int j = ( get_group_id( 0 ) * num_blocks ) + temp;
	if ( j > n - 1 )	{
		j = n - 1;
	}
	float sum = 0;
	
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
		sums_array[ width + j ] = sum;
	}
}