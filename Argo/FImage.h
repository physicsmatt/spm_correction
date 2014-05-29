#ifndef _FIMAGE_H_
#define _FIMAGE_H_

#include <string>
#include <iostream>
#include <fstream>
#include "FreeImage.h"

struct Metadata {
		FREE_IMAGE_FORMAT format;
		FREE_IMAGE_TYPE type;
};

class FImage {
	protected:
	private:
		double* data;
		bool loaded;
		std::string filename;
	public:
		Metadata metadata;
		unsigned int width, height;
		FImage ( Metadata _meta );

		/**
		 * Constructor
		 *
		 * @param file name of file to load.
		 */
		FImage ( std::string file );

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

		double interpPixel ( double x, double y );
		float interpPixelFloat ( float x, float y );

		bool writeImage ( std::string file );
};

#endif _FIMAGE_H_
