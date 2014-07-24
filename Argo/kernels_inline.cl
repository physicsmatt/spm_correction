#ifdef cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif
#pragma OPENCL EXTENSION cl_amd_printf : enable

__kernel void clWarpSliverCubic( __global __const double* sliver, __global double* sliver_array, double A1, double B1,
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
		ux0 = 1.5 * x3 - 2.5 * x2 + 1;
	}
	else if ( sx0 > 1 && sx0 <= 2 ) {
		ux0 = -0.5 * x3 + 2.5 * x2 + -4 * sx0 + 2;
	}
	else {
		ux0 = 0;
	}
	if ( sy0 <= 1 && sy0 >= 0 ) {
		uy0 = 1.5 * y3 - 2.5 * y2 + 1;
	}
	else if ( sy0 > 1 && sy0 <= 2 ) {
		uy0 = -0.5 * y3 + 2.5 * y2 + -4 * sy0 + 2;
	}
	else {
		uy0 = 0;
	}

	x2 = sx1 * sx1;
	x3 = x2 * sx1;
	y2 = sy1 * sy1;
	y3 = y2 * sy1;
	if ( sx1 <= 1 && sx1 >= 0 ) {
		ux1 = 1.5 * x3 - 2.5 * x2 + 1;
	}
	else if ( sx1 > 1 && sx1 <= 2 ) {
		ux1 = -0.5 * x3 + 2.5 * x2 + -4 * sx1 + 2;
	}
	else {
		ux1 = 0;
	}
	if ( sy1 <= 1 && sy1 >= 0 ) {
		uy1 = 1.5 * y3 - 2.5 * y2 + 1;
	}
	else if ( sy1 > 1 && sy1 <= 2 ) {
		uy1 = -0.5 * y3 + 2.5 * y2 + -4 * sy1 + 2;
	}
	else {
		uy1 = 0;
	}

	x2 = sx2 * sx2;
	x3 = x2 * sx2;
	y2 = sy2 * sy2;
	y3 = y2 * sy2;
	if ( sx2 <= 1 && sx2 >= 0 ) {
		ux2 = 1.5 * x3 - 2.5 * x2 + 1;
	}
	else if ( sx2 > 1 && sx2 <= 2 ) {
		ux2 = -0.5 * x3 + 2.5 * x2 + -4 * sx2 + 2;
	}
	else {
		ux2 = 0;
	}
	if ( sy2 <= 1 && sy2 >= 0 ) {
		uy2 = 1.5 * y3 - 2.5 * y2 + 1;
	}
	else if ( sy2 > 1 && sy2 <= 2 ) {
		uy2 = -0.5 * y3 + 2.5 * y2 + -4 * sy2 + 2;
	}
	else {
		uy2 = 0;
	}

	x2 = sx3 * sx3;
	x3 = x2 * sx3;
	y2 = sy3 * sy3;
	y3 = y2 * sy3;
	if ( sx3 <= 1 && sx3 >= 0 ) {
		ux3 = 1.5 * x3 - 2.5 * x2 + 1;
	}
	else if ( sx3 > 1 && sx3 <= 2 ) {
		ux3 = -0.5 * x3 + 2.5 * x2 + -4 * sx3 + 2;
	}
	else {
		ux3 = 0;
	}
	if ( sy3 <= 1 && sy3 >= 0 ) {
		uy3 = 1.5 * y3 - 2.5 * y2 + 1;
	}
	else if ( sy3 > 1 && sy3 <= 2 ) {
		uy3 = -0.5 * y3 + 2.5 * y2 + -4 * sy3 + 2;
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
			int bool_flag = floor_origX == 0 || floor_origY == 0 || ceil_origX == width - 1 || ceil_origY == height - 1;

			double pixelvalue = 0;
			double sum = 0;
		
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

			sum += sliver[ floor_origX - 1 + temp_index - width ] * ux0 * uy0;
			sum += sliver[ floor_origX + temp_index - width ] * ux1 * uy0;
			sum += sliver[ floor_origX + 1 + temp_index - width ] * ux2 * uy0;
			sum += sliver[ floor_origX + 2 + temp_index - width ] * ux3 * uy0;

			sum += sliver[ floor_origX - 1 + temp_index ] * ux0 * uy1;
			sum += sliver[ floor_origX + temp_index ] * ux1 * uy1;
			sum += sliver[ floor_origX + 1 + temp_index ] * ux2 * uy1;
			sum += sliver[ floor_origX + 2 + temp_index ] * ux3 * uy1;

			sum += sliver[ floor_origX - 1 + temp_index + width ] * ux0 * uy2;
			sum += sliver[ floor_origX + temp_index + width ] * ux1 * uy2;
			sum += sliver[ floor_origX + 1 + temp_index + width ] * ux2 * uy2;
			sum += sliver[ floor_origX + 2 + temp_index + width ] * ux3 * uy2;

			sum += sliver[ floor_origX - 1 + temp_index + width + width ] * ux0 * uy3;
			sum += sliver[ floor_origX + temp_index + width + width ] * ux1 * uy3;
			sum += sliver[ floor_origX + 1 + temp_index + width + width ] * ux2 * uy3;
			sum += sliver[ floor_origX + 2 + temp_index + width + width ] * ux3 * uy3;

			// Choose Between Bilinear and Bicubic based upon flag.
			sliver_array[ x + width * ( y * iterations + it + minY ) ] = bool_flag ? pixelvalue : sum;
		}
	}
}

__kernel void clWarpBaseCubic( __global __const double* base, __global double* base_array, double A0, double A1, double A2, double A3, double B0, double B1, double B2, double B3,
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
		ux0 = 1.5 * x3 - 2.5 * x2 + 1;
	}
	else if ( sx0 > 1 && sx0 <= 2 ) {
		ux0 = -0.5 * x3 + 2.5 * x2 + -4 * sx0 + 2;
	}
	else {
		ux0 = 0;
	}
	if ( sy0 <= 1 && sy0 >= 0 ) {
		uy0 = 1.5 * y3 - 2.5 * y2 + 1;
	}
	else if ( sy0 > 1 && sy0 <= 2 ) {
		uy0 = -0.5 * y3 + 2.5 * y2 + -4 * sy0 + 2;
	}
	else {
		uy0 = 0;
	}

	x2 = sx1 * sx1;
	x3 = x2 * sx1;
	y2 = sy1 * sy1;
	y3 = y2 * sy1;
	if ( sx1 <= 1 && sx1 >= 0 ) {
		ux1 = 1.5 * x3 - 2.5 * x2 + 1;
	}
	else if ( sx1 > 1 && sx1 <= 2 ) {
		ux1 = -0.5 * x3 + 2.5 * x2 + -4 * sx1 + 2;
	}
	else {
		ux1 = 0;
	}
	if ( sy1 <= 1 && sy1 >= 0 ) {
		uy1 = 1.5 * y3 - 2.5 * y2 + 1;
	}
	else if ( sy1 > 1 && sy1 <= 2 ) {
		uy1 = -0.5 * y3 + 2.5 * y2 + -4 * sy1 + 2;
	}
	else {
		uy1 = 0;
	}

	x2 = sx2 * sx2;
	x3 = x2 * sx2;
	y2 = sy2 * sy2;
	y3 = y2 * sy2;
	if ( sx2 <= 1 && sx2 >= 0 ) {
		ux2 = 1.5 * x3 - 2.5 * x2 + 1;
	}
	else if ( sx2 > 1 && sx2 <= 2 ) {
		ux2 = -0.5 * x3 + 2.5 * x2 + -4 * sx2 + 2;
	}
	else {
		ux2 = 0;
	}
	if ( sy2 <= 1 && sy2 >= 0 ) {
		uy2 = 1.5 * y3 - 2.5 * y2 + 1;
	}
	else if ( sy2 > 1 && sy2 <= 2 ) {
		uy2 = -0.5 * y3 + 2.5 * y2 + -4 * sy2 + 2;
	}
	else {
		uy2 = 0;
	}

	x2 = sx3 * sx3;
	x3 = x2 * sx3;
	y2 = sy3 * sy3;
	y3 = y2 * sy3;
	if ( sx3 <= 1 && sx3 >= 0 ) {
		ux3 = 1.5 * x3 - 2.5 * x2 + 1;
	}
	else if ( sx3 > 1 && sx3 <= 2 ) {
		ux3 = -0.5 * x3 + 2.5 * x2 + -4 * sx3 + 2;
	}
	else {
		ux3 = 0;
	}
	if ( sy3 <= 1 && sy3 >= 0 ) {
		uy3 = 1.5 * y3 - 2.5 * y2 + 1;
	}
	else if ( sy3 > 1 && sy3 <= 2 ) {
		uy3 = -0.5 * y3 + 2.5 * y2 + -4 * sy3 + 2;
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
			int bool_flag = ( floor_origX == 0 ) || ( floor_origY == 0 ) || ( ceil_origX == bwidth - 1 ) || ( ceil_origY == height - 1 );

			double pixelvalue = 0;
			double sum = 0;

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

			sum += base[ floor_origX - 1 + temp_index - bwidth ] * ux0 * uy0;
			sum += base[ floor_origX + temp_index - bwidth ] * ux1 * uy0;
			sum += base[ floor_origX + 1 + temp_index - bwidth ] * ux2 * uy0;
			sum += base[ floor_origX + 2 + temp_index - bwidth ] * ux3 * uy0;

			sum += base[ floor_origX - 1 + temp_index ] * ux0 * uy1;
			sum += base[ floor_origX + temp_index ] * ux1 * uy1;
			sum += base[ floor_origX + 1 + temp_index ] * ux2 * uy1;
			sum += base[ floor_origX + 2 + temp_index ] * ux3 * uy1;

			sum += base[ floor_origX - 1 + temp_index + bwidth ] * ux0 * uy2;
			sum += base[ floor_origX + temp_index + bwidth ] * ux1 * uy2;
			sum += base[ floor_origX + 1 + temp_index + bwidth ] * ux2 * uy2;
			sum += base[ floor_origX + 2 + temp_index + bwidth ] * ux3 * uy2;

			sum += base[ floor_origX - 1 + temp_index + bwidth + bwidth ] * ux0 * uy3;
			sum += base[ floor_origX + temp_index + bwidth + bwidth ] * ux1 * uy3;
			sum += base[ floor_origX + 1 + temp_index + bwidth + bwidth ] * ux2 * uy3;
			sum += base[ floor_origX + 2 + temp_index + bwidth + bwidth ] * ux3 * uy3;

			base_array[ x * iterations + it + width * y ] = bool_flag ? pixelvalue : sum;
		}
	}
}

__kernel void clWarpSliverBspline( __global __const double* sliver, __global double* sliver_array, double A1, double B1,
						int width, int height, int maxX, int minY, int maxY, __local double* sx, __local double* sy, __local double* ux, __local double* uy )	{
	__local int local_size;
	__local int iterations;
	__local double yoffset;
	__local double origX;
	__local int x;
	__local int B;
	__local int C;
	__private int y = get_local_id( 0 );

	x = get_group_id( 0 );
	local_size = get_local_size( 0 );
	iterations = ceil( ( double ) ( maxY - minY + 1 ) / ( double ) get_local_size( 0 ) );
	yoffset = -( B1 - 1 ) * x;
	origX = x - A1 * x;
	B = 1;
	C = 0;

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
			int bool_flag = floor_origX == 0 || floor_origY == 0 || ceil_origX == width - 1 || ceil_origY == height - 1;

			double pixelvalue = 0;
			double sum = 0;
		
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

			sum += sliver[ floor_origX - 1 + temp_index - width ] * ux0 * uy0;
			sum += sliver[ floor_origX + temp_index - width ] * ux1 * uy0;
			sum += sliver[ floor_origX + 1 + temp_index - width ] * ux2 * uy0;
			sum += sliver[ floor_origX + 2 + temp_index - width ] * ux3 * uy0;

			sum += sliver[ floor_origX - 1 + temp_index ] * ux0 * uy1;
			sum += sliver[ floor_origX + temp_index ] * ux1 * uy1;
			sum += sliver[ floor_origX + 1 + temp_index ] * ux2 * uy1;
			sum += sliver[ floor_origX + 2 + temp_index ] * ux3 * uy1;

			sum += sliver[ floor_origX - 1 + temp_index + width ] * ux0 * uy2;
			sum += sliver[ floor_origX + temp_index + width ] * ux1 * uy2;
			sum += sliver[ floor_origX + 1 + temp_index + width ] * ux2 * uy2;
			sum += sliver[ floor_origX + 2 + temp_index + width ] * ux3 * uy2;

			sum += sliver[ floor_origX - 1 + temp_index + width + width ] * ux0 * uy3;
			sum += sliver[ floor_origX + temp_index + width + width ] * ux1 * uy3;
			sum += sliver[ floor_origX + 1 + temp_index + width + width ] * ux2 * uy3;
			sum += sliver[ floor_origX + 2 + temp_index + width + width ] * ux3 * uy3;

			// Choose Between Bilinear and Bicubic based upon flag.
			sliver_array[ x + width * ( y * iterations + it + minY ) ] = bool_flag ? pixelvalue : sum;
		}
	}
}

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
			int bool_flag = ( floor_origX == 0 ) || ( floor_origY == 0 ) || ( ceil_origX == bwidth - 1 ) || ( ceil_origY == height - 1 );

			double pixelvalue = 0;
			double sum = 0;

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

			sum += base[ floor_origX - 1 + temp_index - bwidth ] * ux0 * uy0;
			sum += base[ floor_origX + temp_index - bwidth ] * ux1 * uy0;
			sum += base[ floor_origX + 1 + temp_index - bwidth ] * ux2 * uy0;
			sum += base[ floor_origX + 2 + temp_index - bwidth ] * ux3 * uy0;

			sum += base[ floor_origX - 1 + temp_index ] * ux0 * uy1;
			sum += base[ floor_origX + temp_index ] * ux1 * uy1;
			sum += base[ floor_origX + 1 + temp_index ] * ux2 * uy1;
			sum += base[ floor_origX + 2 + temp_index ] * ux3 * uy1;

			sum += base[ floor_origX - 1 + temp_index + bwidth ] * ux0 * uy2;
			sum += base[ floor_origX + temp_index + bwidth ] * ux1 * uy2;
			sum += base[ floor_origX + 1 + temp_index + bwidth ] * ux2 * uy2;
			sum += base[ floor_origX + 2 + temp_index + bwidth ] * ux3 * uy2;

			sum += base[ floor_origX - 1 + temp_index + bwidth + bwidth ] * ux0 * uy3;
			sum += base[ floor_origX + temp_index + bwidth + bwidth ] * ux1 * uy3;
			sum += base[ floor_origX + 1 + temp_index + bwidth + bwidth ] * ux2 * uy3;
			sum += base[ floor_origX + 2 + temp_index + bwidth + bwidth ] * ux3 * uy3;

			base_array[ x * iterations + it + width * y ] = bool_flag ? pixelvalue : sum;
		}
	}
}
