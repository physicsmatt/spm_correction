#define _SCL_SECURE_NO_WARNINGS
#include "FImage.h"
#include "FreeImage.h"
#include <iostream>
#include <algorithm>

FImage::FImage( Metadata _meta ) {
	loaded = false;
	data = 0;
	width = 0;
	height = 0;
	metadata.type = _meta.type;
	metadata.format = _meta.format;
	metadata.flipped = _meta.flipped;
}

FImage::FImage ( std::string file ) {
	metadata.flipped = false;
	load( file );
}

FImage::FImage( std::string file, bool _flipped ) {
	metadata.flipped = _flipped;
	load( file );
}

FImage::FImage ( unsigned int size_x, unsigned int size_y, Metadata _meta ) {
	loaded = false;
	width = size_x;
	height = size_y;
	double* _data = new double[ width * height ];
	data = _data;
	metadata.type = _meta.type;
	metadata.format = _meta.format;
	metadata.flipped = _meta.flipped;
}

FImage::~FImage () {
	if ( data ) {
		delete[] data;
	}
}

void FImage::load ( std::string file ) {
	FreeImage_Initialise();
	if ( loaded ) {
		unload();
	}
	filename = file;
	metadata.format = FreeImage_GetFileType( filename.c_str(), 0 );
	FIBITMAP* img = FreeImage_Load( metadata.format, filename.c_str() );
	if ( !img ) {
		printf( "Error loading file: %s.\n", filename.c_str() );
		return;
	}
	if ( !FreeImage_HasPixels( img ) ) {
		printf( "Image %s does not have any pixel data.\n", filename.c_str() );
		FreeImage_Unload( img );
		return;
	}

	// color types:
	// 'uint8   : one channel 8 bit int
	// 'uint16  : one channel 16 bit int
	// 'float32 : one channel 32 bit float
	// 'rgb     : 24 bit rgb color
	// 'rgba    : 32 bit rgba color
	// 'rgb32f  : 3 x 32 bit float color  
	// 'rgba32f : 4 x 32 bit float color  
	metadata.type = FreeImage_GetImageType( img );
	unsigned int bpp = FreeImage_GetBPP( img );

	width = FreeImage_GetWidth( img );
	height = FreeImage_GetHeight( img );
	unsigned int pixels = width * height;

	loaded = true;
	if ( metadata.flipped ) {
		FreeImage_FlipVertical( img );
	}

	switch ( metadata.type ) {
		double* _data;
	case FIT_BITMAP:
		switch ( bpp ) {
			case 8:
				_data = new double[ pixels ];
				data = _data;
				for ( unsigned int y = 0; y < height; ++y ) {
					BYTE *bits = ( BYTE * ) FreeImage_GetScanLine( img, y );
					for ( unsigned int x = 0; x < width; ++x ) {
						data[ x + width * y ] = ( double ) ( bits[ x ] );
					}
				}
				break;
			default:
				loaded = false;
				printf( "File %s contains multiple channels of data, prepare a single channel file for input.\n", filename.c_str() );
		}
		break;
	case FIT_INT16:
		_data = new double[ pixels ];
		data = _data;
		for ( unsigned int y = 0; y < height; ++y ) {
			short *bits = ( short * ) FreeImage_GetScanLine( img, y );
			for ( unsigned int x = 0; x < width; ++x ) {
				data[ x + width * y ] = ( double ) ( bits[ x ] );
			}
		}
		break;
	case FIT_UINT16:
		_data = new double[ pixels ];
		data = _data;
		for ( unsigned int y = 0; y < height; ++y ) {
			unsigned short *bits = ( unsigned short * ) FreeImage_GetScanLine( img, y );
			for ( unsigned int x = 0; x < width; ++x ) {
				data[ x + width * y ] = ( double ) ( bits[ x ] );
			}
		}
		break;
	case FIT_FLOAT:
		_data = new double[ pixels ];
		data = _data;
		for ( unsigned int y = 0; y < height; ++y ) {
			float *bits = ( float * ) FreeImage_GetScanLine( img, y );
			for ( unsigned int x = 0; x < width; ++x ) {
				data[ x + width * y ] = ( double ) ( bits[ x ] );
			}
		}
		break;
	case FIT_DOUBLE:
		_data = new double[ pixels ];
		data = _data;
		for ( unsigned int y = 0; y < height; ++y ) {
			double *bits = ( double * ) FreeImage_GetScanLine( img, y );
			for ( unsigned int x = 0; x < width; ++x ) {
				data[ x + width * y ] = bits[ x ];
			}
		}
		break;
	default:
		printf( "Failed to load file %s!\n", filename.c_str() );
		loaded = false;
		break;
	}
	FreeImage_Unload( img );
	FreeImage_DeInitialise();
}

void FImage::unload () {
	if ( data ) {
		delete[] data;
	}
	data = NULL;
	loaded = false;
}

void FImage::initialize ( unsigned int size_x, unsigned int size_y, Metadata _meta ) {
	width = size_x;
	height = size_y;
	double* _data = new double[ width * height ];
	data = _data;
	metadata.type = _meta.type;
	metadata.format = _meta.format;
	metadata.flipped = _meta.flipped;
}

double FImage::get ( int x, int y ) const {
	if ( x >= 0 && x < width && y >= 0 && y < height ) {
		return ( fastGet( x, y ) );
	}
	printf( "Retreiving an invalid value!\n" );
	return ( 0 );
}

void FImage::set ( int x, int y, double value ) {
	if ( x >= 0 && x < width && y >= 0 && y < height )
		fastSet( x, y, value );
}

double FImage::getRange() {
	double* temp = new double[ width * height ];
	std::copy(data, data + width * height, temp);
	std::sort( temp, temp + width * height );
	int i = 0;
	while ( temp[ i++ ] == -std::numeric_limits<double>::infinity() );
	double min = temp[ i ];
	double max = temp[ width * height - 1 ];
	delete[] temp;
	return ( max - min );
}

double FImage::getMin() {
	double* temp = new double[ width * height ];
	std::copy( data, data + width * height, temp );
	std::sort( temp, temp + width * height );
	int i = 0;
	while ( temp[ i++ ] == -std::numeric_limits<double>::infinity() );
	double min = temp[ i ];
	delete[] temp;
	return ( min );
}

double FImage::getMax() {
	double* temp = new double[ width * height ];
	std::copy( data, data + width * height, temp );
	std::sort( temp, temp + width * height );
	double max = temp[ width * height - 1 ];
	delete[] temp;
	return ( max );
}

bool FImage::writeImage ( std::string file ) {
	FreeImage_Initialise();
	FIBITMAP *image = FreeImage_AllocateT( metadata.type, width, height );
	if ( !image ) {
		printf( "Failed to allocate memory to write image!\n" );
		FreeImage_DeInitialise();
		return ( false );
	}

	switch ( metadata.type ) {
		case FIT_BITMAP:
			for ( unsigned int y = 0; y < height; ++y ) {
				BYTE *bits = ( BYTE * ) FreeImage_GetScanLine( image, y );
				for ( unsigned int x = 0; x < width; ++x ) {
					bits[ x ] = ( BYTE ) round( data[ x + width * y ] );
				}
			}
			break;
		case FIT_INT16:
			for ( unsigned int y = 0; y < height; ++y ) {
				short *bits = ( short * ) FreeImage_GetScanLine( image, y );
				for ( unsigned int x = 0; x < width; ++x ) {
					bits[ x ] = ( short ) round( data[ x + width * y ] );
				}
			}
			break;
		case FIT_UINT16:
			for ( unsigned int y = 0; y < height; ++y ) {
				unsigned short *bits = ( unsigned short * ) FreeImage_GetScanLine( image, y );
				for ( unsigned int x = 0; x < width; ++x ) {
					bits[ x ] = ( unsigned short ) round( data[ x + width * y ] );
				}
			}
			break;
		case FIT_FLOAT:
			for ( unsigned int y = 0; y < height; ++y ) {
				float *bits = ( float * ) FreeImage_GetScanLine( image, y );
				for ( unsigned int x = 0; x < width; ++x ) {
					bits[ x ] = ( float ) data[ x + width * y ];
				}
			}
			break;
		case FIT_DOUBLE:
			for ( unsigned int y = 0; y < height; ++y ) {
				double *bits = ( double * ) FreeImage_GetScanLine( image, y );
				for ( unsigned int x = 0; x < width; ++x ) {
					bits[ x ] = data[ x + width * y ];
				}
			}
			break;
		default:
			loaded = false;
			break;
	}

	if ( metadata.flipped ) {
		FreeImage_FlipVertical( image );
	}

	if ( !FreeImage_Save( metadata.format, image, file.c_str() ) ) {
		printf( "Failed to save file!\n" );
		FreeImage_Unload( image );
		FreeImage_DeInitialise();
		return ( false );
	}

	FreeImage_Unload( image );
	FreeImage_DeInitialise();
	return ( true );
}

bool FImage::writeDisplayableImage( std::string file, double slope, double min ) {
	FreeImage_Initialise();
	FIBITMAP *image = FreeImage_AllocateT( FIT_FLOAT, width, height );
	if ( !image ) {
		printf( "Failed to allocate memory to write image!\n" );
		FreeImage_DeInitialise();
		return ( false );
	}

	for ( unsigned int y = 0; y < height; ++y ) {
		float *bits = ( float * ) FreeImage_GetScanLine( image, y );
		for ( unsigned int x = 0; x < width; ++x ) {
			if ( data[ x + width * y ] == -std::numeric_limits<double>::infinity() ) continue;
			bits[ x ] = ( float ) ( slope * ( data[ x + width * y ] - min ) );
		}
	}

	if ( metadata.flipped ) {
		FreeImage_FlipVertical( image );
	}

	if ( !FreeImage_Save( FIF_TIFF, image, file.c_str() ) ) {
		printf( "Failed to save file!\n" );
		FreeImage_Unload( image );
		FreeImage_DeInitialise();
		return ( false );
	}
	FreeImage_Unload( image );
	FreeImage_DeInitialise();
	return ( true );
}

double FImage::interpPixel( double x, double y, unsigned int type ) {
	switch ( type ) {
		case F_BILINEAR:
			return bilinear( x, y );
		case F_NEAREST:
			return nearestNeighbor( x, y );
		case F_CUBIC:
			if ( floor ( x ) == 0 || floor ( y ) == 0 || floor ( x ) >= width - 2 || floor ( y ) >= height - 2 ) {
				return bilinear( x, y );
			}
			return cubic( x, y );
		default:
			if ( floor ( x ) == 0 || floor ( y ) == 0 || floor ( x ) >= width - 2 || floor ( y ) >= height - 2 ) {
				return bilinear( x, y );
			}
			return bspline( x, y, type );
	}
}

/**
* Returns pixel value from double x,y by bilerping. Has more precision and similar runtime as interp_pixel_float.
* @param baseImage pointer to image
* @param x coordinate
* @param y coordinate
* @return pixel value as double
*/
double FImage::bilinear( double x, double y ) {
	double s, t;
	double left_val, right_val;
	double top_val, bottom_val;

	// Get integer coordinates for top left corner
	int left_index = ( int ) x, top_index = ( int ) y;
	int right_index;
	int bottom_index;

	// Get integer coordinates of bottom and right edges
	if ( left_index > width - 1 )
		right_index = left_index;             // Don't fall off the edge!
	else
		right_index = ( int ) ceil( x );             //left_index + 1;
	if ( top_index > height - 1 )
		bottom_index = top_index;             // Don't fall off the edge!
	else
		bottom_index = ( int ) ceil( y );             //top_index + 1;

	// Interpolate across the top edge
	s = x - ( double ) left_index;
	t = y - ( double ) top_index;
	left_val = fastGet( left_index, top_index );
	right_val = fastGet( right_index, top_index );
	// Linear interpolation recoded to only one multiply
	top_val = s * right_val + ( 1 - s ) * ( left_val );

	// Interpolate across the bottom edge
	left_val = fastGet( left_index, bottom_index );
	right_val = fastGet( right_index, bottom_index );
	bottom_val = s * right_val + ( 1 - s ) * ( left_val );

	// Interpolate between top and bottom
	return ( t * bottom_val + ( 1 - t ) * top_val );
}

double FImage::nearestNeighbor( double x, double y ) {
	return fastGet( round( x ), round( y ) );
}

double FImage::cubic( double x, double y ) {
	double A = -0.5;
	double sx[ 4 ]; // these are the distances for each row
	double sy[ 4 ]; // these are the distances for each column
	double ux[ 4 ]; // these are the weight for the rows
	double uy[ 4 ]; // these are the weights for the columns
	double sum = 0;

	sx[ 0 ] = abs( x - floor ( x - 1 ) ); // these get the distances for each row and column from the initial point, positive
	sy[ 0 ] = abs( y - floor( y - 1 ) );
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

	for ( int i = 0; i < 4; ++i ) {
		for ( int j = 0; j < 4; ++j ) {
			sum += fastGet( floor( x ) - 1 + i, floor( y ) - 1 + j ) * ux[ i ] * uy[ j ];
		}
	}
	return sum;
}

double FImage::bspline( double x, double y, int type ) {
	double B = 0;
	double C = 0;
	if ( type == F_BSPLINE ) {	// Cubic B-spline
		B = 1;
		C = 0;
	}
	else if ( type == F_CATMULL_BSPLINE ) {	// Catmull-Rom Spline
		B = 0;
		C = 0.5;
	}
	else {	// Mitchell-Netravali Cubic Filter
		B = 1 / 3;
		C = 1 / 3;
	}

	double sx[4]; // these are the distances for each row
	double sy[4]; // these are the distances for each column
	double ux[4]; // these are the weight for the rows
	double uy[ 4 ]; // these are the weights for the columns

	double sum = 0;

	sx[ 0 ] = abs( x - floor( x - 1 ) ); // these get the distances for each row and column from the initial point, positive
	sy[ 0 ] = abs( y - floor( y - 1 ) );
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
			ux[ j ] = ( ( 12 - ( 9 * B ) - ( 6 * C ) ) * x3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*x2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sx[ j ] > 1 && sx[ j ] <= 2 ) {
			ux[ j ] = ( ( -B - ( 6 * C ) )*x3 + ( ( 6 * B ) + ( 30 * C ) )*x2 + ( -( 12 * B ) - ( 48 * C ) )*( sx[ j ] ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			ux[ j ] = 0;
		}
		if ( sy[ j ] <= 1 && sy[ j ] >= 0 ) {
			uy[ j ] = ( ( 12 - ( 9 * B ) - ( 6 * C ) )*y3 + ( -18 + ( 12 * B ) + ( 6 * C ) )*y2 + 6 - ( 2 * B ) ) / 6;
		}
		else if ( sy[ j ] > 1 && sy[ j ] <= 2 ) {
			uy[ j ] = ( ( -B - ( 6 * C ) )*y3 + ( ( 6 * B ) + ( 30 * C ) )*y2 + ( -( 12 * B ) - ( 48 * C ) )*( sy[ j ] ) + ( 8 * B ) + ( 24 * C ) ) / 6;
		}
		else {
			uy[ j ] = 0;
		}
	}


	for ( int i = 0; i < 4; ++i )  {
		for ( int j = 0; j < 4; ++j ) {
			sum += fastGet( floor( x ) - 1 + i, floor( y ) - 1 + j ) * ux[ i ] * uy[ j ];
		}
	}

	return sum;
}

void FImage::warpSliver( FImage* _warped, double aterms[], double bterms[], double cterms[], int interp_type ) {
	int x = 0;
	int y = 0;

	double origX = 0;
	double origY = 0;
	double newz = 0;
	double pixelvalue = 0;

	for ( x = 0; x < width; x++ ) {
		double yoffset = -( bterms[ 1 ] - 1 ) * x;
		origX = x - aterms[ 1 ] * x;
		for ( y = 0; y < height; y++ ) {
			origY = y + yoffset;
			if ( origX < 0 || origX > width - 1 || origY < 0 || origY > height - 1 ) {
				_warped->set( x, y, -std::numeric_limits<double>::infinity() );
			}
			else {
				newz = -cterms[ 1 ] * origX;
				pixelvalue = interpPixel( origX, origY, interp_type ) + newz;
				_warped->set( x, y, pixelvalue );
			}
		}
	}
}

/**
* Based upon given parameter values, the base image is warped.
*
* @param aterms A parameter values
* @param bterms B parameter values
* @param cterms C parameter values
* @param dominantAxis Determines which axis (x or y) to apply the warping to.
*/
void FImage::warpBase( FImage* _warped, double aterms[], double bterms[], double cterms[], int dominant_axis, int interp_type ) {
	int x;
	int y;

	double origX = 0;
	double origY = 0;
	double newz = 0;
	double pixelvalue = 0;

	if ( dominant_axis == 1 ) {
		for ( y = 0; y < height; y++ ) {
			int y2 = y * y;
			int y3 = y2 * y;
			double xoffset = aterms[ 0 ] + aterms[ 1 ] * y + aterms[ 2 ] * y2 + aterms[ 3 ] * y3;
			origY = bterms[ 0 ] + bterms[ 1 ] * y + bterms[ 2 ] * y2 + bterms[ 3 ] * y3;
			newz = cterms[ 0 ] + cterms[ 1 ] * y + cterms[ 2 ] * y2 + cterms[ 3 ] * y3;
			for ( x = 0; x < width; x++ ) {
				origX = x + xoffset;
				if ( origX < 0 || origX > width - 1 || origY < 0 || origY > height - 1 ) {
					_warped->set( x, y, -std::numeric_limits<double>::infinity() );
				}
				else {
					pixelvalue = interpPixel( origX, origY, interp_type ) + newz;
					_warped->set( x, y, pixelvalue );
				}
			}
		}
	}
	else {
		for ( x = 0; x < width; x++ ) {
			int x2 = x * x;
			int x3 = x2 * x;
			for ( y = 0; y < height; y++ ) {
				origX = aterms[ 0 ] + ( aterms[ 1 ] + 1 ) * x + aterms[ 2 ] * x2 + aterms[ 3 ] * x3;
				origY = y + bterms[ 0 ] + ( bterms[ 1 ] - 1 ) * x + bterms[ 2 ] * x2 + bterms[ 3 ] * x3;
				if ( origX < 0 || origX > width - 1 || origY < 0 || origY > height - 1 ) {
					_warped->set( x, y, -std::numeric_limits<double>::infinity() );
				}
				else {
					pixelvalue = interpPixel( origX, origY, interp_type );
					_warped->set( x, y, pixelvalue );
				}
			}
		}
	}
}

/**
* This function down-samples an image by an integer factor.
* Assumes values of images are integers
* This is NOT a particularly fast function.
* @param factor
*/
void FImage::resample( FImage* resamp, int factor ) {
	
	unsigned int rx_size = width / factor;
	unsigned int ry_size = height / factor;
	resamp->initialize( rx_size, ry_size, metadata );

	double inv_area = 1.0 / ( factor * factor );
	for ( unsigned int ry = 0; ry < ry_size; ++ry ) {
		for ( unsigned int rx = 0; rx < rx_size; ++rx ) {
			double pixel_value = 0;
			for ( unsigned int by = ry * factor; by < ( ry + 1 ) * factor; ++by )
				for ( unsigned int bx = rx * factor; bx < ( rx + 1 ) * factor; ++bx )
					pixel_value += get( bx, by );
			resamp->set( rx, ry, pixel_value * inv_area );
		}
	}
}