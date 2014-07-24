#ifndef _FIMAGE_H_
#define _FIMAGE_H_

#include <string>
#include "FreeImage.h"



#define F_BILINEAR 0x1000
#define F_NEAREST 0x1001
#define F_CUBIC 0x1002
#define F_BSPLINE 0x1003
#define F_CATMULL_BSPLINE 0x1004
#define F_MITCHELL_NETRAVALI_BSPLINE 0x1005

struct Metadata {
		FREE_IMAGE_FORMAT format;
		FREE_IMAGE_TYPE type;
		bool flipped;
};

class FImage {
	protected:
	private:
		bool loaded;
		std::string filename;

		double bilinear( double x, double y );
		double cubic( double x, double y );
		double bspline( double x, double y, int type );
		double nearestNeighbor( double x, double y );
	public:
		Metadata metadata;
		unsigned int width, height;
		double* data;

		FImage( Metadata _meta );

		/**
		 * Constructor
		 *
		 * @param file name of file to load.
		 */
		FImage ( std::string file );
		FImage( std::string file, bool _flipped );

		FImage ( unsigned int size_x, unsigned int size_y, Metadata _meta );

		~FImage ();

		void load ( std::string file );

		void unload ();

		void initialize ( unsigned int size_x, unsigned int size_y, Metadata _meta );

		double fastGet ( int x, int y ) const {
			return ( data[ x + width * y ] );
		}
		double get ( int x, int y ) const;
		void fastSet ( int x, int y, double value ) {
			data[ x + width * y ] = value;
		}
		void set ( int x, int y, double value );

		double interpPixel ( double x, double y, unsigned int type );

		void resample( FImage* img, int factor );

		void warpBase( FImage* img, double aterms[], double bterms[], double cterms[], int dominant_axis, int interp_type );

		void warpSliver( FImage* img, double aterms[], double bterms[], double cterms[], int interp_type );

		bool writeImage ( std::string file );

		bool writeDisplayableImage( std::string file, double slope, double min );

		double getRange();

		double getMin();

		double getMax();
};

#endif _FIMAGE_H_
