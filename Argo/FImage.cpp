#include "FImage.h"

FImage::FImage ( Metadata _meta ) {
	loaded = false;
	data = 0;
	width = 0;
	height = 0;
	metadata.type = _meta.type;
	metadata.format = _meta.format;
}

FImage::FImage ( std::string file ) {
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
}

FImage::~FImage () {
	if ( data ) {
		delete[] data;
	}
}

void FImage::load ( std::string file ) {
	if ( loaded ) {
		unload();
	}
	filename = file;
	metadata.format = FreeImage_GetFileType( filename.c_str(), 0 );
	FIBITMAP* img = FreeImage_Load( metadata.format, filename.c_str() );
	if ( !img ) {
		printf( "Error loading file: %s.\n", filename );
		return;
	}
	if ( !FreeImage_HasPixels( img ) ) {
		printf( "Image %s does not have any pixel data.\n", filename );
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
				printf( "File %s contains multiple channels of data, prepare a single channel file for input.\n", filename );
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
		printf( "Failed to load file %s!\n", filename );
		loaded = false;
		break;
	}
	FreeImage_Unload( img );
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
}

double FImage::get ( int x, int y ) const {
	if ( x >= 0 && x < width && y >= 0 && y < height ) {
		return ( fastGet( x, y ) );
	}
	return ( 0 );
}

void FImage::set ( int x, int y, double value ) {
	if ( x >= 0 && x < width && y >= 0 && y < height )
		fastSet( x, y, value );
}

bool FImage::writeImage ( std::string file ) {
	FIBITMAP *image = FreeImage_AllocateT( metadata.type, width, height );
	if ( !image ) {
		printf( "Failed to allocate memory to write image!\n" );
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

	if ( !FreeImage_Save( metadata.format, image, file.c_str() ) ) {
		printf( "Failed to save file!\n" );
		return ( false );
	}
	FreeImage_Unload( image );
	return ( true );
}

/**
 * Returns pixel value from float x,y by bilerping (could use trilinear algorithm...)
 * @param baseImage pointer to image
 * @param x coordinate
 * @param y coordinate
 * @return pixel value as float
 */
float FImage::interpPixelFloat ( float x, float y ) {

	float s, t;
	int image_width = width;
	int image_height = height;
	float left_val, right_val;
	float top_val, bottom_val;

	// Get integer coordinates for top left corner
	int left_index = ( int ) x, top_index = ( int ) y;
	int right_index;
	int bottom_index;

	// Get integer coordinates of bottom and right edges
	if ( left_index > image_width )
		right_index = left_index;             // Don't fall off the edge!
	else
		right_index = ( int ) ceil( x );             //left_index + 1;
	if ( top_index > image_height )
		bottom_index = top_index;             // Don't fall off the edge!
	else
		bottom_index = ( int ) ceil( y );             //top_index + 1;

	// Interpolate across the top edge
	s = x - ( float ) left_index;
	t = y - ( float ) top_index;
	left_val = get( left_index, top_index );
	right_val = get( right_index, top_index );
	// Linear interpolation recoded to only one multiply
	top_val = s * ( float ) right_val + ( 1 - s ) * ( float ) left_val; //right_val + s * (left_val-right_val);

	// Interpolate across the bottom edge
	left_val = get( left_index, bottom_index );
	right_val = get( right_index, bottom_index );
	//	bottom_val = right_val + s * (left_val-right_val);
	//Matt changed above to below
	bottom_val = s * ( float ) right_val + ( 1 - s ) * ( float ) left_val;

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
double FImage::interpPixel ( double x, double y ) {
	double s, t;
	int image_width = width;
	int image_height = height;
	double left_val, right_val;
	double top_val, bottom_val;

	// Get integer coordinates for top left corner
	int left_index = ( int ) x, top_index = ( int ) y;
	int right_index;
	int bottom_index;

	// Get integer coordinates of bottom and right edges
	if ( left_index > image_width )
		right_index = left_index;             // Don't fall off the edge!
	else
		right_index = ( int ) ceil( x );             //left_index + 1;
	if ( top_index > image_height )
		bottom_index = top_index;             // Don't fall off the edge!
	else
		bottom_index = ( int ) ceil( y );             //top_index + 1;

	// Interpolate across the top edge
	s = x - ( double ) left_index;
	t = y - ( double ) top_index;
	left_val = get( left_index, top_index );
	right_val = get( right_index, top_index );
	// Linear interpolation recoded to only one multiply
	top_val = s * right_val + ( 1 - s ) * ( left_val ); //right_val + s * (left_val-right_val);

	// Interpolate across the bottom edge
	left_val = get( left_index, bottom_index );
	;
	right_val = get( right_index, bottom_index );
	//	bottom_val = right_val + s * (left_val-right_val);
	//Matt changed above to below
	bottom_val = s * right_val + ( 1 - s ) * ( left_val );

	// Interpolate between top and bottom
	return ( t * bottom_val + ( 1 - t ) * top_val );	//(bottom_val + t * (top_val-bottom_val));
}
