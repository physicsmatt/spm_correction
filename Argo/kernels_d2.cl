#pragma OPENCL EXTENSION cl_amd_printf : enable

__kernel void clWarpSliverCubic( __global __const double* sliver, __global double* sliver_array, double A1, double B1,
						int width, int height, int maxX, int minY, int maxY )	{
	__private int x = get_global_id( 0 );
	__private int y = get_global_id( 1 ) + minY;

	__private int bool_flag = x > maxX || y > maxY;
	x = ( maxX & bool_flag ) | ( x & ~bool_flag );
	y = ( maxY & bool_flag ) | ( y & ~bool_flag );


	__private double origX = x - A1 * x;
	__private double origY = y + -( B1 - 1 ) * x;
	__private double sx0 = fabs( origX - floor( origX - 1 ) ); // these get the distances for each row and column from the initial point, positive
	__private double sy0 = fabs( origY - floor( origY - 1 ) );
	__private double sx1 = fabs( 1 - sx0 );
	__private double sy1 = fabs( 1 - sy0 );
	__private double sx2 = fabs( 1 - sx1 );
	__private double sy2 = fabs( 1 - sy1 );
	__private double ux0;
	__private double ux1;
	__private double ux2;
	__private double ux3;
	__private double uy0;
	__private double uy1;
	__private double uy2;
	__private double uy3;

	if ( sx0 <= 1 && sx0 >= 0 ) {
		ux0 = 1.5 * sx0 * sx0 * sx0 - 2.5 * sx0 * sx0 + 1;
	}
	else if ( sx0 > 1 && sx0 <= 2 ) {
		ux0 = -0.5 * sx0 * sx0 * sx0 + 2.5 * sx0 * sx0 + -4 * sx0 + 2;
	}
	else {
		ux0 = 0;
	}
	if ( sy0 <= 1 && sy0 >= 0 ) {
		uy0 = 1.5 * sy0 * sy0 * sy0 - 2.5 * sy0 * sy0 + 1;
	}
	else if ( sy0 > 1 && sy0 <= 2 ) {
		uy0 = -0.5 * sy0 * sy0 * sy0 + 2.5 * sy0 * sy0 + -4 * sy0 + 2;
	}
	else {
		uy0 = 0;
	}

	if ( sx1 <= 1 && sx1 >= 0 ) {
		ux1 = 1.5 * sx1 * sx1 * sx1 - 2.5 * sx1 * sx1 + 1;
	}
	else if ( sx1 > 1 && sx1 <= 2 ) {
		ux1 = -0.5 * sx1 * sx1 * sx1 + 2.5 * sx1 * sx1 + -4 * sx1 + 2;
	}
	else {
		ux1 = 0;
	}
	if ( sy1 <= 1 && sy1 >= 0 ) {
		uy1 = 1.5 * sy1 * sy1 * sy1 - 2.5 * sy1 * sy1 + 1;
	}
	else if ( sy1 > 1 && sy1 <= 2 ) {
		uy1 = -0.5 * sy1 * sy1 * sy1 + 2.5 * sy1 * sy1 + -4 * sy1 + 2;
	}
	else {
		uy1 = 0;
	}

	if ( sx2 <= 1 && sx2 >= 0 ) {
		ux2 = 1.5 * sx2 * sx2 * sx2 - 2.5 * sx2 * sx2 + 1;
	}
	else if ( sx2 > 1 && sx2 <= 2 ) {
		ux2 = -0.5 * sx2 * sx2 * sx2 + 2.5 * sx2 * sx2 + -4 * sx2 + 2;
	}
	else {
		ux2 = 0;
	}
	if ( sy2 <= 1 && sy2 >= 0 ) {
		uy2 = 1.5 * sy2 * sy2 * sy2 - 2.5 * sy2 * sy2 + 1;
	}
	else if ( sy2 > 1 && sy2 <= 2 ) {
		uy2 = -0.5 * sy2 * sy2 * sy2 + 2.5 * sy2 * sy2 + -4 * sy2 + 2;
	}
	else {
		uy2 = 0;
	}

	if ( ( 1 + sx2 ) <= 1 && ( 1 + sx2 ) >= 0 ) {
		ux3 = 1.5 * ( 1 + sx2 ) * ( 1 + sx2 ) * ( 1 + sx2 ) - 2.5 * ( 1 + sx2 ) * ( 1 + sx2 ) + 1;
	}
	else if ( ( 1 + sx2 ) > 1 && ( 1 + sx2 ) <= 2 ) {
		ux3 = -0.5 * ( 1 + sx2 ) * ( 1 + sx2 ) * ( 1 + sx2 ) + 2.5 * ( 1 + sx2 ) * ( 1 + sx2 ) + -4 * ( 1 + sx2 ) + 2;
	}
	else {
		ux3 = 0;
	}
	if ( ( 1 + sy2 ) <= 1 && ( 1 + sy2 ) >= 0 ) {
		uy3 = 1.5 * ( 1 + sy2 ) * ( 1 + sy2 ) * ( 1 + sy2 ) - 2.5 * ( 1 + sy2 ) * ( 1 + sy2 ) + 1;
	}
	else if ( ( 1 + sy2 ) > 1 && ( 1 + sy2 ) <= 2 ) {
		uy3 = -0.5 * ( 1 + sy2 ) * ( 1 + sy2 ) * ( 1 + sy2 ) + 2.5 * ( 1 + sy2 ) * ( 1 + sy2 ) + -4 * ( 1 + sy2 ) + 2;
	}
	else {
		uy3 = 0;
	}

	
	__private int floor_origX = ( int ) origX;
	__private int floor_origY = ( int ) origY;
	__private int ceil_origX = ceil( origX );
	__private int ceil_origY = ceil( origY );
	bool_flag = floor_origX == 0 || floor_origY == 0 || ceil_origX == width - 1 || ceil_origY == height - 1;

	__private double pixelvalue = 0;
	__private double sum = 0;

	// Bilinear Interpolation
	__private double s = origX - floor_origX;
	__private double t = origY - floor_origY;
	__private double left_val, right_val;
	__private double top_val, bottom_val;
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
	sum += sliver[ floor_origX - 1 + width * ( floor_origY - 1 ) ] * ux0 * uy0;
	sum += sliver[ floor_origX + width * ( floor_origY - 1 ) ] * ux1 * uy0;
	sum += sliver[ floor_origX + 1 + width * ( floor_origY - 1 ) ] * ux2 * uy0;
	sum += sliver[ floor_origX + 2 + width * ( floor_origY - 1 ) ] * ux3 * uy0;
	sum += sliver[ floor_origX - 1 + width * ( floor_origY ) ] * ux0 * uy1;
	sum += sliver[ floor_origX + width * ( floor_origY ) ] * ux1 * uy1;
	sum += sliver[ floor_origX + 1 + width * ( floor_origY ) ] * ux2 * uy1;
	sum += sliver[ floor_origX + 2 + width * ( floor_origY ) ] * ux3 * uy1;
	sum += sliver[ floor_origX - 1 + width * ( floor_origY + 1 ) ] * ux0 * uy2;
	sum += sliver[ floor_origX + width * ( floor_origY + 1 ) ] * ux1 * uy2;
	sum += sliver[ floor_origX + 1 + width * ( floor_origY + 1 ) ] * ux2 * uy2;
	sum += sliver[ floor_origX + 2 + width * ( floor_origY + 1 ) ] * ux3 * uy2;
	sum += sliver[ floor_origX - 1 + width * ( floor_origY + 2 ) ] * ux0 * uy3;
	sum += sliver[ floor_origX + width * ( floor_origY + 2 ) ] * ux1 * uy3;
	sum += sliver[ floor_origX + 1 + width * ( floor_origY + 2 ) ] * ux2 * uy3;
	sum += sliver[ floor_origX + 2 + width * ( floor_origY + 2 ) ] * ux3 * uy3;

	sliver_array[ x + width * y ] = bool_flag ? pixelvalue : sum;
}

__kernel void clWarpBaseCubic( __global __const double* base, __global double* base_array, double A0, double A1, double A2, double A3, double B0, double B1, double B2, double B3,
						int width, int height, int bwidth, int xbp, int maxX, int minY, int maxY )	{
	__private int x = get_global_id( 0 );
	__private int y = get_global_id( 1 ) + minY;

	__private int bool_flag = x > maxX || y > maxY;
	x = ( maxX & bool_flag ) | ( x & ~bool_flag );
	y = ( maxY & bool_flag ) | ( y & ~bool_flag );

	__private int ym2 = y * y;
	__private int ym3 = ym2 * y;
	__private double origX = x + xbp + A0 + A1 * y + A2 * ym2 + A3 * ym3;
	__private double origY = B0 + B1 * y + B2 * ym2 + B3 * ym3;

	__private double sx0 = fabs( origX - floor( origX - 1 ) ); // these get the distances for each row and column from the initial point, positive
	__private double sy0 = fabs( origY - floor( origY - 1 ) );
	__private double sx1 = fabs( 1 - sx0 );
	__private double sy1 = fabs( 1 - sy0 );
	__private double sx2 = fabs( 1 - sx1 );
	__private double sy2 = fabs( 1 - sy1 );
	__private double ux0;
	__private double ux1;
	__private double ux2;
	__private double ux3;
	__private double uy0;
	__private double uy1;
	__private double uy2;
	__private double uy3;

	if ( sx0 <= 1 && sx0 >= 0 ) {
		ux0 = 1.5 * sx0 * sx0 * sx0 - 2.5 * sx0 * sx0 + 1;
	}
	else if ( sx0 > 1 && sx0 <= 2 ) {
		ux0 = -0.5 * sx0 * sx0 * sx0 + 2.5 * sx0 * sx0 + -4 * sx0 + 2;
	}
	else {
		ux0 = 0;
	}
	if ( sy0 <= 1 && sy0 >= 0 ) {
		uy0 = 1.5 * sy0 * sy0 * sy0 - 2.5 * sy0 * sy0 + 1;
	}
	else if ( sy0 > 1 && sy0 <= 2 ) {
		uy0 = -0.5 * sy0 * sy0 * sy0 + 2.5 * sy0 * sy0 + -4 * sy0 + 2;
	}
	else {
		uy0 = 0;
	}

	if ( sx1 <= 1 && sx1 >= 0 ) {
		ux1 = 1.5 * sx1 * sx1 * sx1 - 2.5 * sx1 * sx1 + 1;
	}
	else if ( sx1 > 1 && sx1 <= 2 ) {
		ux1 = -0.5 * sx1 * sx1 * sx1 + 2.5 * sx1 * sx1 + -4 * sx1 + 2;
	}
	else {
		ux1 = 0;
	}
	if ( sy1 <= 1 && sy1 >= 0 ) {
		uy1 = 1.5 * sy1 * sy1 * sy1 - 2.5 * sy1 * sy1 + 1;
	}
	else if ( sy1 > 1 && sy1 <= 2 ) {
		uy1 = -0.5 * sy1 * sy1 * sy1 + 2.5 * sy1 * sy1 + -4 * sy1 + 2;
	}
	else {
		uy1 = 0;
	}

	if ( sx2 <= 1 && sx2 >= 0 ) {
		ux2 = 1.5 * sx2 * sx2 * sx2 - 2.5 * sx2 * sx2 + 1;
	}
	else if ( sx2 > 1 && sx2 <= 2 ) {
		ux2 = -0.5 * sx2 * sx2 * sx2 + 2.5 * sx2 * sx2 + -4 * sx2 + 2;
	}
	else {
		ux2 = 0;
	}
	if ( sy2 <= 1 && sy2 >= 0 ) {
		uy2 = 1.5 * sy2 * sy2 * sy2 - 2.5 * sy2 * sy2 + 1;
	}
	else if ( sy2 > 1 && sy2 <= 2 ) {
		uy2 = -0.5 * sy2 * sy2 * sy2 + 2.5 * sy2 * sy2 + -4 * sy2 + 2;
	}
	else {
		uy2 = 0;
	}

	if ( ( 1 + sx2 ) <= 1 && ( 1 + sx2 ) >= 0 ) {
		ux3 = 1.5 * ( 1 + sx2 ) * ( 1 + sx2 ) * ( 1 + sx2 ) - 2.5 * ( 1 + sx2 ) * ( 1 + sx2 ) + 1;
	}
	else if ( ( 1 + sx2 ) > 1 && ( 1 + sx2 ) <= 2 ) {
		ux3 = -0.5 * ( 1 + sx2 ) * ( 1 + sx2 ) * ( 1 + sx2 ) + 2.5 * ( 1 + sx2 ) * ( 1 + sx2 ) + -4 * ( 1 + sx2 ) + 2;
	}
	else {
		ux3 = 0;
	}
	if ( ( 1 + sy2 ) <= 1 && ( 1 + sy2 ) >= 0 ) {
		uy3 = 1.5 * ( 1 + sy2 ) * ( 1 + sy2 ) * ( 1 + sy2 ) - 2.5 * ( 1 + sy2 ) * ( 1 + sy2 ) + 1;
	}
	else if ( ( 1 + sy2 ) > 1 && ( 1 + sy2 ) <= 2 ) {
		uy3 = -0.5 * ( 1 + sy2 ) * ( 1 + sy2 ) * ( 1 + sy2 ) + 2.5 * ( 1 + sy2 ) * ( 1 + sy2 ) + -4 * ( 1 + sy2 ) + 2;
	}
	else {
		uy3 = 0;
	}


	__private int floor_origX = ( int ) origX;
	__private int floor_origY = ( int ) origY;
	__private int ceil_origX = ceil( origX );
	__private int ceil_origY = ceil( origY );
	bool_flag = floor_origX == 0 || floor_origY == 0 || ceil_origX == bwidth - 1 || ceil_origY == height - 1;

	__private double pixelvalue = 0;
	__private double sum = 0;

	__private double s = origX - floor_origX;
	__private double t = origY - floor_origY;
	__private double left_val, right_val;
	__private double top_val, bottom_val;
	// Interpolate across the top edge
	s = origX - floor_origX;
	t = origY - floor_origY;
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

	sum += base[ floor_origX - 1 + bwidth * ( floor_origY - 1 ) ] * ux0 * uy0;
	sum += base[ floor_origX + bwidth * ( floor_origY - 1 ) ] * ux1 * uy0;
	sum += base[ floor_origX + 1 + bwidth * ( floor_origY - 1 ) ] * ux2 * uy0;
	sum += base[ floor_origX + 2 + bwidth * ( floor_origY - 1 ) ] * ux3 * uy0;

	sum += base[ floor_origX - 1 + bwidth * ( floor_origY ) ] * ux0 * uy1;
	sum += base[ floor_origX + bwidth * ( floor_origY ) ] * ux1 * uy1;
	sum += base[ floor_origX + 1 + bwidth * ( floor_origY ) ] * ux2 * uy1;
	sum += base[ floor_origX + 2 + bwidth * ( floor_origY ) ] * ux3 * uy1;

	sum += base[ floor_origX - 1 + bwidth * ( floor_origY + 1 ) ] * ux0 * uy2;
	sum += base[ floor_origX + bwidth * ( floor_origY + 1 ) ] * ux1 * uy2;
	sum += base[ floor_origX + 1 + bwidth * ( floor_origY + 1 ) ] * ux2 * uy2;
	sum += base[ floor_origX + 2 + bwidth * ( floor_origY + 1 ) ] * ux3 * uy2;

	sum += base[ floor_origX - 1 + bwidth * ( floor_origY + 2 ) ] * ux0 * uy3;
	sum += base[ floor_origX + bwidth * ( floor_origY + 2 ) ] * ux1 * uy3;
	sum += base[ floor_origX + 1 + bwidth * ( floor_origY + 2 ) ] * ux2 * uy3;
	sum += base[ floor_origX + 2 + bwidth * ( floor_origY + 2 ) ] * ux3 * uy3;

	// Final write.
	base_array[ x + width * y ] = bool_flag ? pixelvalue : sum;
}