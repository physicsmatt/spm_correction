/*********************************************************************
*Image Master procedural image registration and stabilization routine*
*Coded by: Brian Salmons & Dr. Matthew Trawick                       *
**********************************************************************
*Dr. Matthew Trawick                                                 *
*University of Richmond                                              *
**********************************************************************
*/

/*Version 1.0
*
*Known Issues:
*-Sliverdrift is SIGNIFICANTLY slower.  Need to try to speed this up.
*-Images still have to flipped.  It seems... random.  I need to find the cause of this.
*
*-Significant speed improvements for non-sliverdrift version.
*-Fixed annoying division bug resulting from changes to the way areas are calculated
*-Program now outputs the version built on (up to the coder to change it)
*-Initial simplex values are now no longer output... they aren't needed.
*-Significant changes to the way that sliver drift computations occur.
*-Printouts from simplex supressed.
*/




//#define WINVER 0x0501 

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <math.h>
#include <time.h>
#include "Eigen\Eigen"
//#include <ctime>	No need to use this as a result of time.h being included.
#include <stdio.h>
//#include <fstream>	iostream included so why include this?
//#include <IStream>	again why?
#include <valarray>
#include <algorithm>
#include "image_basic.h"
#include "simplex.h"
#include <ppl.h>

using namespace Concurrency;
//#define SANITYCHECKS
#define VERSION "1.0"
#define codefold
//Define struct for holding the 6 Beta values used to evaluate the C coefficients.
struct beta_values{
	double B0, B1, B2, B3, B4, B5;
};
//Define struct for holding the 3 gamma values represinting the 3 points on our parabola used in the creation of diffs.
struct gamma_values{
	double G0, G1, G2;
};
//Definition of parameter Combo Struct.
struct param_combo {
	double P1, P2, P3;
};


#ifdef codefold
// Basically just return floor() for negatives and ceil() for positives. Gosh why such a crappy implementation?
int rnd(double a) {
	if (a>=0)
		return int(a + 0.5);
	else
		return int(a - 0.5);
}

bool operator<(const param_combo& a, const param_combo& b) {
	if (a.P1 != b.P1)
		return a.P1 < b.P1;
	else if (a.P2 != b.P2)
		return a.P2 < b.P2;
	else
		return a.P3 < b.P3;
}

bool operator>(const param_combo& a, const param_combo& b) {
	if (a.P1 != b.P1)
		return a.P1 > b.P1;
	else if (a.P2 != b.P2)
		return a.P2 > b.P2;
	else
		return a.P3 > b.P3;
}

bool operator<=(const param_combo& a, const param_combo& b) {
	if (a.P1 != b.P1)
		return a.P1 <= b.P1;
	else if (a.P2 != b.P2)
		return a.P2 <= b.P2;
	else
		return a.P3 <= b.P3;
}

bool operator>=(const param_combo& a, const param_combo& b) {
	if (a.P1 != b.P1)
		return a.P1 >= b.P1;
	else if (a.P2 != b.P2)
		return a.P2 >= b.P2;
	else
		return a.P3 >= b.P3;
}


// returns pixel value from float x,y by bilerping (could use trilinear algorithm...)
float interp_pixel_float(image_basic *baseImage, float x, float y) {

	float s, t;
	int image_width = baseImage->width;
	int image_height = baseImage->height;
	int left_val, right_val;
	float top_val, bottom_val;

	// Get integer coordinates for top left corner
	int left_index = (int)x, top_index = (int)y;
	int right_index;
	int bottom_index;

	// Get integer coordinates of bottom and right edges
	if (left_index > image_width)
		right_index = left_index;             // Don't fall off the edge!
	else
		right_index = (int)ceil(x);//left_index + 1;
	if (top_index > image_height)
		bottom_index = top_index;             // Don't fall off the edge!
	else
		bottom_index = (int)ceil(y);//top_index + 1;

	// Interpolate across the top edge
	s = x - (float)left_index;
	t = y - (float)top_index;
	left_val = baseImage->get(left_index, top_index);
	right_val = baseImage->get(right_index, top_index);
	// Linear interpolation recoded to only one multiply
	top_val = s*right_val + (1 - s)*(left_val);//right_val + s * (left_val-right_val);

	// Interpolate across the bottom edge
	left_val = baseImage->get(left_index, bottom_index);
	right_val = baseImage->get(right_index, bottom_index);
	//	bottom_val = right_val + s * (left_val-right_val);
	//Matt changed above to below
	bottom_val = s* right_val + (1 - s) * (left_val);

	// Interpolate between top and bottom
	return (t*bottom_val + (1 - t)*top_val);//(bottom_val + t * (top_val-bottom_val));

}

// returns pixel value from float x,y by bilerping
//  The "double" version ran just as fast as the "float" version in one of my tests, so
// I'm keeping it as a double, even if it might not be strictly necessary all the time.
double interp_pixel(image_basic *baseImage, double x, double y) {

	double s, t;
	int image_width = baseImage->width;
	int image_height = baseImage->height;
	int left_val, right_val;
	double top_val, bottom_val;

	// Get integer coordinates for top left corner
	int left_index = (int)x, top_index = (int)y;
	int right_index;
	int bottom_index;

	// Get integer coordinates of bottom and right edges
	if (left_index > image_width)
		right_index = left_index;             // Don't fall off the edge!
	else
		right_index = (int)ceil(x);//left_index + 1;
	if (top_index > image_height)
		bottom_index = top_index;             // Don't fall off the edge!
	else
		bottom_index = (int)ceil(y);//top_index + 1;

	// Interpolate across the top edge
	s = x - (double)left_index;
	t = y - (double)top_index;
	left_val = baseImage->get(left_index, top_index);
	right_val = baseImage->get(right_index, top_index);
	// Linear interpolation recoded to only one multiply
	top_val = s*right_val + (1 - s)*(left_val);//right_val + s * (left_val-right_val);

	// Interpolate across the bottom edge
	left_val = baseImage->get(left_index, bottom_index);
	right_val = baseImage->get(right_index, bottom_index);
	//	bottom_val = right_val + s * (left_val-right_val);
	//Matt changed above to below
	bottom_val = s* right_val + (1 - s) * (left_val);

	// Interpolate between top and bottom
	return (t*bottom_val + (1 - t)*top_val);//(bottom_val + t * (top_val-bottom_val));

}

void warp_image(image_basic *baseImage, int sizex, int sizey, double aterms[], double bterms[], double cterms[], int dominantAxis, image_basic *warpedImage)
{
	//we have to retain the integrity of the old image, so we have to make a new, unwarped version.
	//image_basic warpedImage(sizex,sizey,1);
	int pixelx = 0;
	int pixely = 0;

	double newx = 0;
	double newy = 0;
	double newz = 0;
	int pixelvalue = 0;

	if (dominantAxis == 1)
	{
		for (pixely = 0; pixely<sizey; pixely++)
		{
			int pixely2 = pixely*pixely;
			int pixely3 = pixely2 *pixely;
			for (pixelx = 0; pixelx<sizex; pixelx++)
			{
				//math.pow() is slightly slower than explicit multiplication.  Since this will be happening thousands of times, 
				//this IS slightly faster.  Please note that since
				//Nathan's edit (saving the value is even faster...moved calculation up to reduce redundancy by a factor of 1000.)

				newx = pixelx + aterms[0] + aterms[1] * pixely + aterms[2] * pixely2 + aterms[3] * pixely3;
				newy = bterms[0] + bterms[1] * pixely + bterms[2] * pixely2 + bterms[3] * pixely3;
				newz = cterms[0] + cterms[1] * pixely + cterms[2] * pixely2 + cterms[3] * pixely3;
				if (newx < 0 || newx >= baseImage->width || newy <0 || newy >= baseImage->height){
					warpedImage->set(pixelx, pixely, -1);
				}
				else
				{
					pixelvalue = rnd(interp_pixel(baseImage, (float)newx, (float)newy) + newz);     //where the new value is actually computed.
					//note,above, the cheap-ass implementation of rounding
					//ceil(newz);
					//pixelvalue+= newz;
					warpedImage->set(pixelx, pixely, pixelvalue);
				}
			}
		}
	}
	else
	{
		for (pixelx = 0; pixelx<sizex; pixelx++)
		{
			int pixelx2 = pixelx*pixelx;
			int pixelx3 = pixelx2*pixelx;
			for (pixely = 0; pixely<sizey; pixely++)
			{
				//math.pow() is slightly slower than explicit multiplication.  Since this will be happening thousands of times, 
				//this IS slightly faster.  Please note that since 
				//Nathan's edit. Again. see above. Reversed order of for loops 5/23/2011
				newx = aterms[0] + (aterms[1] + 1)*pixelx + aterms[2] * pixelx2 + aterms[3] * pixelx3;
				newy = pixely + bterms[0] + (bterms[1] - 1)*pixelx + bterms[2] * pixelx2 + bterms[3] * pixelx3;
				if (newx < 0 || newx >= baseImage->width || newy <0 || newy >= baseImage->height){
					warpedImage->set(pixelx, pixely, -1);
				}
				else
				{
					pixelvalue = rnd(interp_pixel(baseImage, (float)newx, (float)newy));     //where the new value is actually computed.
					warpedImage->set(pixelx, pixely, pixelvalue);

				}
			}
		}

	}
	//return warpedImage;
}

//This method makes an array of pixelvalues, rather than the image_basic inhouse valarray usage, which, while fast, is also highly confusing.
void val2Array(image_basic *baseImage, int **image)
{
	int sizex = baseImage->width;
	int sizey = baseImage->height;
	int pixelx = 0;
	int pixely = 0;

	for (pixely = 0; pixely<sizey; pixely++)
	{
		for (pixelx = 0; pixelx<sizex; pixelx++)
		{
			image[pixelx][pixely] = (int)interp_pixel(baseImage, (float)pixelx, (float)pixely);
			//note: since pixelx and pixely are integers, I have no idea why Brian is using interp_pixel for this.
		}
	}
}

//Compute the array index in an x-contiguous array.
int compute1D(int sizex, int sizey, int x, int y){
	return (x + sizex * y);
}

//This method makes an array of pixelvalues, rather than the image_basic inhouse valarray usage, which, while fast, is also highly confusing.
void val21DArray(image_basic *baseImage, int image[])
{
	int sizex = baseImage->width;
	int sizey = baseImage->height;
	int pixelx = 0;
	int pixely = 0;

	for (pixely = 0; pixely<sizey; pixely++)
	{
		for (pixelx = 0; pixelx<sizex; pixelx++)
		{
			image[compute1D(sizex, sizey, pixelx, pixely)] = baseImage->get(pixelx, pixely);
		}
	}
}

//This method makes an array of pixelvalues, rather than the image_basic inhouse valarray usage, which, while fast, is also highly confusing.
void val22DArray(image_basic *baseImage, int image[])
{
	int sizex = baseImage->width;
	int sizey = baseImage->height;
	int pixelx = 0;
	int pixely = 0;

	for (pixely = 0; pixely<sizey; pixely++)
	{
		for (pixelx = 0; pixelx<sizex; pixelx++)
		{
			image[pixelx, pixely] = baseImage->get(pixelx, pixely);
		}
	}
}

void nearestneighbor(double x, double y, double* newx, double* newy)
{
	*newx = rnd(x);
	*newy = rnd(y);
}

/* Seriously, why create min(), max() templated functions, especially when you never use them
int min(int x, int y)
{
	if (x < y)
		return x;
	else
		return y;
}


int max(int x, int y)
{
	if (x < y)
		return y;
	else
		return x;
}
*/

//This function down-samples an image by an integer factor.  
//Assumes values of images are integers
//This is NOT a particularly fast function.
void resample(image_basic *resamp, image_basic base, int factor)
{
	unsigned int rx_size = base.width / factor;
	unsigned int ry_size = base.height / factor;
	//resamp->set_color_mode(base.get_color_mode());
	resamp->initialize(rx_size, ry_size, base.get_color_mode());
	//int cm = base.get_color_mode();
	//image_basic *resamp(rx_size,ry_size,base.get_color_mode()); 

	double inv_area = 1.0 / (factor*factor);
	for (unsigned int ry = 0; ry < ry_size; ++ry) {
		for (unsigned int rx = 0; rx < rx_size; ++rx) {
			int pixel_value = 0;
			for (unsigned int by = ry*factor; by<(ry + 1)*factor; ++by)
			for (unsigned int bx = rx*factor; bx<(rx + 1)*factor; ++bx)
				pixel_value += base.fast_get(bx, by);
			resamp->set(rx, ry, rnd(pixel_value*inv_area));
		}
	}
}
/*************************************************************************************************************************/
/*main begins here*/
/*************************************************************************************************************************/

int main(int argc, char *argv[]){

#ifdef codefold
	int begintime = clock();
	bool debugmode = true;
	//Starting declarations
	int blocksize;
	char c;
	bool output;
	bool displaydiff;
	bool testind;
	clock_t time0, time1, time2, time3, time4;
	//time_t time0,time1,time2,time3;  //Brian used time_t, instead of clock_t

	int MaxA0 = 5;  //maximum constant shift, in pixels
	int MaxB0 = 5;
	double a1_Max = 0.1;  //maximum linear shift, as fraction of image size (if a1_Max = 0.02, and image height is 500 pixels, this term causes max shift at top of image of 10 pixels)
	double b1_Max = 0.1;
	double a2_Multiplier = (1.0 / 9.0);  //fractional shift, relative to a1_max.  If a1_max implies max shift of 10 pixels at top, a2_Multiplier of 0.5 causes deviation from linearity of 5 pixels at midpoint.
	double b2_Multiplier = (1.0 / 9.0);
	double a3_Multiplier = .1;//(0.0/9.0);  //fractional shift for cubic term a3, also relative to shift caused by a1_max. 
	double b3_Multiplier = .1;//(0.0/9.0);  

	//For initial grid search, image is rescaled (downsampled) by this factor.  
	//This parameter should be set to the size, in pixels, of the smallest real feature on the image.
	int precision = 1;

	//These parameters control the simplex routine.
	//These parameters were found, using trial and error on not a whole lot of different test cases.
	//To avoid being trapped in some rut, it seemse to be VERY important to have the reflection
	//parameter slightly less than 1.0, and the growth parameter less than 2.0.  
	double growthParam = 1.5;
	double contractParam = 0.5;
	double reflectParam = 0.9;
	double haltParam = 1e-11;
	int maxRefineIterationsParam = 5000;

	//blocksize=8;

	std::string rundata1, rundata2;


	printf("Running on version %s\n", VERSION);

	if ((argc < 33) || (argc > 33)) {
		//perror("Error: Usage is ImageMaster -i: ImagePath -s: SliverPath -A0: value -A1: value -A2: value -B0: value -B1: value -B2: value -p: precision -b: blocksize -g: growth -c: contract -r: reflection\n");	
		//ACTUALLY, there should be at least 31 arguments now, since I added cubic terms  
		//Also, usage is -0 A0value, -1 B0value, -2 a1value -3 b1value -4 a2multiplier -5 b2multiplier, etc.
		printf("%d", argc);
		getchar();
		//return 0;
	}

	for (int i = 1; i < argc; ++i) {
		if (argv[i][0] == '-') {
			c = argv[i][1];

			/* found an option, so look at next
			* argument to get the value of
			* the option */
			switch (c) {
			case 'i':
				rundata1 = argv[i + 1];
				for (int x = 0; x<rundata1.length(); x++){
					if (rundata1[x] == '|')
						rundata1[x] = ' ';
				}
				break;
			case 's':
				rundata2 = argv[i + 1];
				for (int x = 0; x<rundata2.length(); x++){
					if (rundata2[x] == '|')
						rundata2[x] = ' ';
				}
				break;
			case '0':
				//A0percent = atof(argv[i+1]);
				MaxA0 = atoi(argv[i + 1]);
				break;
			case '1':
				//A1percent = atof(argv[i+1]);
				MaxB0 = atoi(argv[i + 1]);
				break;
			case '2':
				//A2percent = atof(argv[i+1]);
				a1_Max = atof(argv[i + 1]);
				break;
			case '3':
				//B0percent = atof(argv[i+1]);
				b1_Max = atof(argv[i + 1]);
				break;
			case '4':
				//B1percent = atof(argv[i+1]);
				a2_Multiplier = atof(argv[i + 1]);
				break;
			case '5':
				//B2percent = atof(argv[i+1]);
				b2_Multiplier = atof(argv[i + 1]);
				break;
			case '6':
				a3_Multiplier = atof(argv[i + 1]);
				break;
			case '7':
				b3_Multiplier = atof(argv[i + 1]);
				break;
			case 'p':
				precision = atoi(argv[i + 1]);
				break;
			case 'b':
				blocksize = atoi(argv[i + 1]);
				break;
			case 'g':
				growthParam = atof(argv[i + 1]);
				break;
			case 'c':
				contractParam = atof(argv[i + 1]);
				break;
			case 'r':
				reflectParam = atof(argv[i + 1]);
				break;
			case 'h':
				haltParam = atof(argv[i + 1]);
				break;
			case 't':
				maxRefineIterationsParam = atoi(argv[i + 1]);
				break;

			default:
				break;
			}
		}
	}


	//printing out paramters:
	printf("MaxA0 %d\n", MaxA0);
	printf("MaxB0 %d\n", MaxB0);
	printf("Maxa1 %f\n", a1_Max);
	printf("Maxb1 %f\n", b1_Max);
	printf("a2 %f\n", a2_Multiplier);
	printf("b2 %f\n", b2_Multiplier);
	printf("a3 %f\n", a3_Multiplier);
	printf("b3 %f\n", b3_Multiplier);
	printf("precision %d\n", precision);
	printf("growth %f\n", growthParam);
	printf("contraction %f\n", contractParam);
	printf("reflection %f\n", reflectParam);
	printf("errhalt %f\n", haltParam);
	printf("maxi %d\n", maxRefineIterationsParam);


	//Note, this program is currently set up to request a file to open.  This may be changed to windows style
	//or .NET usage, should the need arise.
	//string rundata1="Y:\\Brian Salmons\\ImageMaster\\C++\\ImageMaster\\ImageMaster\\256abstractb.tif";
	//string rundata2="Y:\\Brian Salmons\\ImageMaster\\C++\\ImageMaster\\ImageMaster\\256abstractsliverb.tif";
	FILE *stream;
	char line[100];

	if ((stream = fopen("Y:\\Brian Salmons\\ImageMaster\\C++\\ImageMaster\\ImageMaster\\testing\\run_data.txt", "r")) != NULL)
	{
		if (fgets(line, 100, stream) == NULL)
			printf("fgets error\n");
		else
			printf("%s", line);

		if (fgets(line, 100, stream) == NULL)
			printf("fgets error\n");
		else
			printf("%s", line);
		fclose(stream);
	}



	if (argc != 1){
		for (int i = 1; i<argc; i++)
		{
			if (argv[i] == "/testdiff" || argv[i] == "testdiff"){
				//this is an obselete flag, should be removed in final v8
			}
			if (argv[i] == "/precision" || argv[i] == "precision"){

			}
			if (argv[i] == "/blocksize" || argv[i] == "blocksize"){

			}
			if (argv[i] == "/output" || argv[i] == "output"){

			}
			if (argv[i] == "/displaydiff" || argv[i] == "displaydiff"){

			}
			if (argv[i] == "/testind" || argv[i] == "testind"){

			}
		}
	}



	//This is a software level flag.  This is simply for convenience to the coder so that simple runs within
	//the IDE spawn full detail and precision.  Furthermore, all the flags can be set
	//this statement.

	if (debugmode){
		//precision=1;
		blocksize = 8;
		output = true;
		displaydiff = true;
		testind = true;
	}

	/*If precision has not been specified, this tells the algorithm to spawn the code at maximum
	*precision (see above for variable descriptions)
	*/
	if (precision == 0)
		precision = 1;

	//HANDLE FILE INPUT HERE!

	//Remember, this "image_basic" is actually quite useful, and as such, we will keep this around.
	image_basic base_Image, sliver_Image, resamp_base_Image, resamp_sliver_Image;

	//base_Image.file_read_tiff("stars_n0_s0_warpedWithZ.tif");
	//sliver_Image.file_read_tiff("stars_n0_s0_sliverWithZ.tif");
	base_Image.file_read_tiff(rundata1);
	sliver_Image.file_read_tiff(rundata2);

	//resample the image and sliver, making them smaller, according to the "precision" parameter
	resample(&resamp_base_Image, base_Image, precision);
	resample(&resamp_sliver_Image, sliver_Image, precision);

	//THIS IS AN IMPORTANT DEFINITION.  HEIGHT will only be used within diffs creation and the main grid search.
	int HEIGHT = resamp_base_Image.height;
	//string base = "base.tif";
	//string sliver = "sliver.tif";
	//resamp_base_Image.file_write_tiff(base);
	//resamp_sliver_Image.file_write_tiff(sliver);
	MaxA0 = (int)ceil((float)MaxA0 / precision);
	MaxB0 = (int)ceil((float)MaxB0 / precision);

	//	int imagesize[2];
	//	imagesize[0]=base_Image.width,
	//	imagesize[1]=base_Image.height;

	//	int sliversize[2];
	//	sliversize[0]=sliver_Image.width;
	//	sliversize[1]=sliver_Image.height;



	double bestA0, bestA1, bestA2, bestA3;
	double bestB0, bestB1, bestB2, bestB3;
	/*	Make the best diff smaller. Also, while performing diff summations, if diff goes past bestdiff, break and continue
	Honestly, bestdiff should not be so huge.	*/
	double bestdiff = 1e37;		//essentially infinity

#endif
	// START COMMENTING TO SKIP GRID SEARCH HERE


	//This is the start of code block for dealing with new, smarter way of parametrizing the variables.
	//double a1_Max = (.05 * HEIGHT) / HEIGHT ; //hard coded for now.  This is as a percentage of the image height.
	double a1_step = 1.0 / HEIGHT; //that's 1 pixel
	double a2_Max = a2_Multiplier * (a1_Max * 4 / HEIGHT);  //one third of max deviation from linear part.  This is totally arbitrary.
	double a2_step = 4.0 / (HEIGHT * HEIGHT);
	double a3_Max = a3_Multiplier * 20.78 * a1_Max / (HEIGHT * HEIGHT);
	double a3_step = 20.78 / (HEIGHT * HEIGHT * HEIGHT); //changed constant from 10.4 on 1/15/2007 <--- Nathan would like to know where magic 20.78 comes from.
	//kludge below:
	//a3_step *= 2.0;
	//a2_step *= 1.3;
	//a1_step *= 1.0;
	//force "0.0" to be one of the actual values tested for all parameters
	a1_Max = ceil(a1_Max / a1_step) * a1_step;
	a2_Max = ceil(a2_Max / a2_step) * a2_step;
	a3_Max = ceil(a3_Max / a3_step) * a3_step;

	//double b1_Max = (.05 * HEIGHT) / HEIGHT ; //hard coded for now.  This is as a percentage of the image height.
	double b1_step = 1.0 / HEIGHT; //that's 1 pixel
	double b2_Max = b2_Multiplier * (b1_Max * 4 / HEIGHT);
	double b2_step = 4.0 / (HEIGHT * HEIGHT);
	double b3_Max = b3_Multiplier * 20.78 * b1_Max / (HEIGHT * HEIGHT);
	double b3_step = 20.78 / (HEIGHT * HEIGHT * HEIGHT); //changed constant from 10.4 on 1/15/2007
	//kludge below:
	//b3_step *= 2.0;
	//b2_step *= 1.3;
	//b1_step *= 1.0;
	//force "0.0" to be one of the actual values tested for all parameters
	b1_Max = ceil(b1_Max / b1_step) * b1_step;
	b2_Max = ceil(b2_Max / b2_step) * b2_step;
	b3_Max = ceil(b3_Max / b3_step) * b3_step;

	printf("%f %f %f \n", a1_Max / a1_step, a2_Max / a2_step, a3_Max / a3_step);
	printf("%f %f %f \n", b1_Max / b1_step, b2_Max / b2_step, b3_Max / b3_step);

	//;This set of varialbes record how many pixels of drift in either axis that will be encountered by the program.
	long totaldriftA = MaxA0 + (int)ceil(1.38 * a1_Max * HEIGHT);
	long totaldriftB = MaxB0 + (int)ceil(1.38 * b1_Max * HEIGHT);
	//note that the factor of 1.38 is hard-coded for a2_multiplier = 1/3, a3_multiplier = 1/9

	printf("drifts: %d, %d\n", MaxA0, MaxB0);
	printf("total drifts: %ld, %ld\n", totaldriftA, totaldriftB);
	Eigen::Matrix4d A;
	Eigen::MatrixXd B(4, 1);
	Eigen::MatrixXd C(4, 1);
	Eigen::Matrix4d D;
	Eigen::MatrixXd E(4, 1);


	// There used to be a #endif here... Moved to right before where commenting starts to skip grid search



	//******************
	//* Block Calculations
	//*******************
	// Block size is determined by the maximum drift that may occur at a maximum.
	// This may be changed to become more dynamic in later versions, however, as of
	//	now it is separate from the loop, and so worst case must be assumed.
	double maxslope, maxslopeA, maxslopeB;
	maxslopeA = a1_Max + a2_Max*HEIGHT + .5*HEIGHT*HEIGHT*a3_Max;
	maxslopeB = b1_Max + b2_Max*HEIGHT + .5*HEIGHT*HEIGHT*b3_Max;
	if (maxslopeB > maxslopeA)
		maxslope = maxslopeB;
	else
		maxslope = maxslopeA;
	if (blocksize <= 0)
	{
		//		double maxslope=max((MaxA1+2*MaxA2*(subsize[1]-1)),MaxB1+2*MaxB2*(subsize[1]-1)); //An approximation
		//maxslope = a1_Max + a2_Max*HEIGHT + 0.5*HEIGHT*HEIGHT*a3_Max;
		//This is the correct calculation, though it leads to an inconveniently small blocksize.
		//Better to hardcode at 8 or something?
		//Note: need to adjust this in case b1's are bigger.
		blocksize = (int)(1 / maxslope);
	}
	//int remainderexists=((HEIGHT % blocksize) > 0); //we do this as an int so that it may be included in calculations
	int numblocks;//= HEIGHT/blocksize + (remainderexists);  //this adds one if there is a remainder;
	//int remainderarea = HEIGHT % blocksize;


	int Asliverdrift = (int)(maxslopeA*resamp_sliver_Image.width);
	int Bsliverdrift = (int)(maxslopeB*resamp_sliver_Image.width);
	int Adiffletpoints = (totaldriftA + Asliverdrift) * 2 + 1;  //get factor of two correct here.
	int Bdiffletpoints = (totaldriftB + Bsliverdrift) * 2 + 1;
	printf("height removed %d\n", (totaldriftB + 2 * Bsliverdrift + MaxB0));
	numblocks = (HEIGHT - totaldriftB - 2 * Bsliverdrift - MaxB0) / blocksize;

	//Hard code the block information, for now.
	//numblocks = subsize[1] / blocksize * 0.9;
	//	remainderexists = 0;
	//remainderarea = 0;

	int Apoints = totaldriftA * 2 + 1;
	int Bpoints = totaldriftB * 2 + 1;
	int num_sliver_blocks = 8;//resamp_sliver_Image.width/blocksize; //hardcode for now
	std::cout << resamp_sliver_Image.width << "\n";
	int sliver_block_width = resamp_sliver_Image.width / num_sliver_blocks;
	std::cout << sliver_block_width << "\n";
	time0 = clock();

	//*********************
	//*Diffs Array Creation
	//*********************

	//int smin=imagesize[0]/2-sliver_Image.width/2;
	int smin = resamp_base_Image.width / 2 - resamp_sliver_Image.width / 2;
	unsigned long long partial_sum1;
	unsigned long long partial_sum2;
	unsigned long long partial_sum3;
	unsigned long long partial_sum4;
	unsigned long long partial_sum5;
	unsigned long long partial_sum6;
	unsigned long long * diffs = new unsigned long long[numblocks*Apoints*Bpoints];
	long long diff;

	//Difflets Array Creation
	//*******************************************************
	//*******************************************************
	//*******************************************************
	//************************NEW THINGS*********************
	//*******************************************************
	//*******************************************************
	//*******************************************************
	//*******************************************************

	double * infoC = new double[4];
	double bestC0, bestC1, bestC2, bestC3;
	double changeInZForBlock;
	//int smin=resamp_base_Image.width /2- resamp_sliver_Image.width/2;
	//unsigned long partial_sum;
	beta_values * dynamic_diffs = new beta_values[numblocks*Apoints*Bpoints];

	//figure out maximum sliver drift, make difflets larger by that much.  (Apoints+Asliverpoints)
	gamma_values * difflets = new gamma_values[numblocks*Adiffletpoints*Bdiffletpoints*num_sliver_blocks];
	gamma_values difflet;
	printf("blocklet width: %d blockHeight: %d number of blocklets: %d number of blocks %d\n",
		resamp_sliver_Image.width / num_sliver_blocks, blocksize, num_sliver_blocks, numblocks);
	std::cout << "\nBsliverdrift: " << Bsliverdrift << "\n";

	for (long int sb = 0; sb<num_sliver_blocks; ++sb) {
		long int sliver_block_x = sb*sliver_block_width;
		for (long int j = 0; j <= 2 * (totaldriftB + Bsliverdrift); j++) {
			for (long int i = 0; i <= 2 * (totaldriftA + Asliverdrift); i++) {
				long int x_initial = smin - totaldriftA - Asliverdrift + i;
				x_initial += sliver_block_x;
				for (long int b = 0; b < numblocks; b++) {
					partial_sum1 = 0;
					partial_sum2 = 0;
					partial_sum3 = 0;
					//partial_sum4=0;
					//partial_sum5=0;
					//partial_sum6=0;
					//for(long int y=MaxB0 - totaldriftB - Bsliverdrift + j + b*blocksize, ys=MaxB0+b*blocksize; ys<MaxB0+(b+1)*blocksize; y++, ys++) {
					for (long int y = MaxB0 - totaldriftB + j + b*blocksize, ys = MaxB0 + Bsliverdrift + b*blocksize; ys<MaxB0 + (b + 1)*blocksize + Bsliverdrift; y++, ys++) {
						if (y < 0) {
							continue;
						}// just skip all negative elements.  
						//for(unsigned int x= x_initial, xs=0; xs < resamp_sliver_Image.width; x++, xs++) {
						for (long int x = x_initial, xs = sliver_block_x; xs < sliver_block_x + sliver_block_width; x++, xs++) {
							//capture 3 points with the differents in z;
							//this gives us three points on our perfect parabola
							//***************************
							///**Changed from main - sliver
							///********june 22 2010********
							//diff = resamp_sliver_Image.fast_get(xs,ys)-resamp_base_Image.fast_get(x,y)-1; //start at -1 delZ
							diff = -resamp_sliver_Image.fast_get(xs, ys) + resamp_base_Image.fast_get(x, y) - 1;
							partial_sum1 += diff*diff;
							//printf("X %d Y %d XS %d YS %d ps %lld diff %lld\n",x,y,xs,ys,partial_sum1,diff);
							//getchar();
							diff++;
							partial_sum2 += diff*diff;
							diff++;
							partial_sum3 += diff*diff;
						}
					}
					double y1 = partial_sum1;
					double y2 = partial_sum2;
					double y3 = partial_sum3;

					double g2 = (y1 + y3) / 2.0 - y2;
					double g1 = (y3 - y1) / 2.0;
					double g0 = y2;

					difflet.G0 = g0;
					difflet.G1 = g1;
					difflet.G2 = g2;

					//diffs[b+i*numblocks+j*numblocks*Apoints] = partial_sum;
					difflets[(b + i*numblocks + j*numblocks*Adiffletpoints + sb*numblocks*Adiffletpoints*Bdiffletpoints)] = difflet;

					//printf("diffs: %d, %d, %d, %lu\n",b,i,j,diffs[b+i*numblocks+j*numblocks*Apoints]);
				} //end for b
				//getchar();
			} //end for i
		} //end for j
	} //end for sb


	//unsigned long long sum0, sum1,term0,term1;
	//double function;

	//We use Matt's trick here of getting a list of indexes here, even though we have to use a for loop.  Pre-calculating the positions IS faster,
	//and Matt's trick precomputes as much as possible.

	unsigned long long count = 0;
	//float delz;

	//long * diffs_indices1=(long *)&(diffs_indices[1]);  //interestingly, using this never bought me much!
	double * yb_arr = new double[numblocks];  //for a 2nd order polynomial, this precalculation of yb and yb^2 made it slower.
	double * yb2_arr = new double[numblocks];
	double * yb3_arr = new double[numblocks];
	double * yb4_arr = new double[numblocks];
	double * yb5_arr = new double[numblocks];
	double * yb6_arr = new double[numblocks];
	int mult1, mult2, mult3;
	mult3 = numblocks*Adiffletpoints*Bdiffletpoints;
	mult2 = numblocks*Adiffletpoints;
	mult1 = numblocks*Apoints;
	//Grab the time
	time1 = clock();
	for (int i = 0; i<numblocks; i++){
		yb_arr[i] = i*blocksize + MaxB0; //+ Bsliverdrift;
		yb2_arr[i] = yb_arr[i] * yb_arr[i];
		yb3_arr[i] = yb2_arr[i] * yb_arr[i];
		yb4_arr[i] = yb3_arr[i] * yb_arr[i];
		yb5_arr[i] = yb4_arr[i] * yb_arr[i];
		yb6_arr[i] = yb5_arr[i] * yb_arr[i];
	}

	//*********************
	//*Make lists of all vaues of A1,A2,A3, B1,B2,B3 to try
	//*********************
	long int number_of_B_combos = ((int)(b1_Max / b1_step) * 2 + 1)*((int)(b2_Max / b2_step) * 2 + 1)*((int)(b3_Max / b3_step) * 2 + 1);
	//double (*B_combos)[3] = new double[number_of_B_combos][3];
	param_combo * B_combos = new param_combo[number_of_B_combos];

	long int B_combos_counter = 0;
	for (double b1 = -b1_Max; b1 <= b1_Max + b1_step / 2; b1 += b1_step){
		for (double b2 = -b2_Max; b2 <= b2_Max + b2_step / 2; b2 += b2_step){
			for (double b3 = -b3_Max; b3 <= b3_Max + b3_step / 2; b3 += b3_step){
				B_combos[B_combos_counter].P1 = b1 - b2*HEIGHT + 0.5*HEIGHT*HEIGHT*b3;
				B_combos[B_combos_counter].P2 = b2 - 1.5 * HEIGHT * b3;
				B_combos[B_combos_counter].P3 = b3;
				++B_combos_counter;
			}
		}
	}
	//	cout << "hi";
	//getchar();
	std::sort(B_combos, (B_combos + (number_of_B_combos)));

	long int number_of_A_combos = ((int)(a1_Max / a1_step) * 2 + 1)*((int)(a2_Max / a2_step) * 2 + 1)*((int)(a3_Max / a3_step) * 2 + 1);
	//double (*A_combos)[3] = new double[number_of_A_combos][3];
	param_combo * A_combos = new param_combo[number_of_A_combos];
	long int A_combos_counter = 0;
	for (double a1 = -a1_Max; a1 <= a1_Max + a1_step / 2; a1 += a1_step){
		for (double a2 = -a2_Max; a2 <= a2_Max + a2_step / 2; a2 += a2_step){
			for (double a3 = -a3_Max; a3 <= a3_Max + a3_step / 2; a3 += a3_step){
				A_combos[A_combos_counter].P1 = a1 - a2*HEIGHT + 0.5*HEIGHT*HEIGHT*a3;
				A_combos[A_combos_counter].P2 = a2 - 1.5 * HEIGHT * a3;
				A_combos[A_combos_counter].P3 = a3;
				++A_combos_counter;
			}
		}
	}
	std::sort(A_combos, (A_combos + (number_of_A_combos)));

	time2 = clock();

	//TO DO for tomorrow: think about how to do initial and final B1 sensibly and symmetrically
	//Incorporate actual sliver_width into calculation, no hard-wiring to 32.
	//Refill DIFFS

	long int Aindex_f = 0, Bindex_f = 0;

	for (long int Bindex_i = 0; Bindex_i < number_of_B_combos - 1; Bindex_i = Bindex_f) {
		double B1_i = B_combos[Bindex_i].P1;
		double B1_f = B1_i + .01;//(1.0/32);
		//search ahead for index Bi_f
		for (Bindex_f = Bindex_i + 1; (B_combos[Bindex_f].P1 < B1_f) && (Bindex_f < number_of_B_combos); Bindex_f++);

		for (long int Aindex_i = 0; Aindex_i < number_of_A_combos - 1; Aindex_i = Aindex_f) {
			double A1_i = A_combos[Aindex_i].P1;

			double A1_f = A1_i + .01;//(1.0/32);
			//search ahead for index Bi_f
			for (Aindex_f = Aindex_i + 1; (A_combos[Aindex_f].P1 < A1_f) && (Aindex_f < number_of_A_combos); Aindex_f++);

			//CALCULATE NEW DIFFS HERE
			//printf("calculating new diffs\n");
			//getchar();

			for (long int diffs_ind = 0; diffs_ind<Apoints*Bpoints*numblocks; diffs_ind++){
				dynamic_diffs[diffs_ind].B0 = 0;
				dynamic_diffs[diffs_ind].B1 = 0;
				dynamic_diffs[diffs_ind].B2 = 0;
				dynamic_diffs[diffs_ind].B3 = 0;
				dynamic_diffs[diffs_ind].B4 = 0;
				dynamic_diffs[diffs_ind].B5 = 0;
			}

			long int diffs_address;
			long int difflets_address;
			int Xc = -sliver_block_width;
			int totX;
			int totY;
			int blocktimesApoints = Apoints*numblocks;
			long A0_B0_adj_inc = numblocks;  //A0increase is now assumed to be 1.
			long MaxA0_times_numblocks = MaxA0*numblocks;
			long B0_adj_min = -MaxB0*blocktimesApoints;
			long B0_adj_max = MaxB0*blocktimesApoints;
			//long BO_adj_inc = B0increase*numblocks*Apoints;
			long BO_adj_inc = blocktimesApoints;  //B0increase is now assumed to be 1.

			int finalAddress = -mult3;
			for (int sb = 0; sb<num_sliver_blocks; ++sb) {
				//Now calculate offsets sx and symeh
				finalAddress += mult3;
				Xc += sliver_block_width;
				int Xc2 = Xc * Xc;
				int doubleXc = 2 * Xc;
				int sx = Asliverdrift - (int)(Xc * A1_i);  //which A1 do I use?  first?  last?  average?
				int sy = Bsliverdrift - (int)(Xc * B1_i);
				//printf("%f, %d, %d, %d\n",A1_i,sb,sx,sy);
				int YPart = -mult1;
				for (int dy = 0; dy<Bpoints; ++dy) {
					totY = (dy + sy)*mult2;
					YPart += mult1;
					int XPart = -numblocks;
					for (int dx = 0; dx<Apoints; ++dx) {
						totX = (dx + sx)*numblocks;
						XPart += (int)numblocks;
						for (int b = 0; b<numblocks; ++b){  //partially unroll this loop?
							//diffs[dy,dx,b] += difflets[sb,dy+sy,dx+sx,b]

							//Xc=0;
							//printf("X position: %f \n",Xc);

							//difflets_address =(b + (dx+sx)*numblocks + (dy+sy)*numblocks*Adiffletpoints + sb*numblocks*Adiffletpoints*Bdiffletpoints);
							difflets_address = (b + totX + totY + finalAddress);
							//diffs_address = b+dx*numblocks+dy*numblocks*Apoints;
							diffs_address = b + XPart + YPart;
							dynamic_diffs[diffs_address].B0 += difflets[difflets_address].G0;
							dynamic_diffs[diffs_address].B1 += difflets[difflets_address].G1;
							dynamic_diffs[diffs_address].B2 += difflets[difflets_address].G2;
							dynamic_diffs[diffs_address].B3 += (difflets[difflets_address].G1 * Xc);
							dynamic_diffs[diffs_address].B4 += difflets[difflets_address].G2 * Xc2;
							dynamic_diffs[diffs_address].B5 += (difflets[difflets_address].G2 * doubleXc);
							//}
						}
					}
				}
			}

			//printf("Betas for block 5\n B0: %g B1: %f B2 %f B3 %f B4 %f B5 %f",
			//dynamic_diffs[5].B0, dynamic_diffs[5].B1,dynamic_diffs[5].B2
			//				,dynamic_diffs[5].B3,dynamic_diffs[5].B4,dynamic_diffs[5].B5);
			//getchar();

			//#define numBytes 8 //change to 4 if you make 
			//difflets, base_diffs_adress, diffs indices, and dynamic diffs 32 bit instead of 64
			critical_section cs;
			combinable<double> sums;
			//for (long int Bindex=Bindex_i; Bindex<Bindex_f; ++Bindex) {
			parallel_for((long int)Bindex_i, Bindex_f, (long int)1, [&](long int Bindex){
				long * diffs_indices = new  long[numblocks];
				for (long int Aindex = Aindex_i; Aindex<Aindex_f; ++Aindex) {
					//					printf("combo indices: %d %d \n",Bindex,Aindex);
					double B1 = B_combos[Bindex].P1;
					double B2 = B_combos[Bindex].P2;
					double B3 = B_combos[Bindex].P3;
					double A1 = A_combos[Aindex].P1;
					double A2 = A_combos[Aindex].P2;
					double A3 = A_combos[Aindex].P3;
					//for (long int i=0;i<numblocks;i++){
					//long yb=i*blocksize + MaxB0drift;// + (blocksize/2);
					//long xpos = (long) (     A1*yb + A2*yb*yb + A3*yb*yb*yb + totaldriftA);
					//long ypos = (long) (     B1*yb + B2*yb*yb + B3*yb*yb*yb + totaldriftB);
					//int xpos = (int) (      totaldriftA);
					//int ypos = (int) (   totaldriftB);

					//long xpos = (long) (     A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + totaldriftA);
					//long ypos = (long) (     B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + totaldriftB);


					//long yb=i*blocksize + MaxB0drift;// + (blocksize/2);
					//int xpos = (int) (     A1*yb_arr[i] + A2*yb2_arr[i] + totaldriftA);
					//int ypos = (int) ( (B1-1)*yb_arr[i] + B2*yb2_arr[i] + totaldriftB);
					//diffs_indices[i]=xpos*numblocks+ypos*numblocks*Apoints+i;



					//diffs_indices[i]=(int) ( (      A1*yb + A2*yb*yb + totaldriftA  )
					//						+( (  (B1-1)*yb + B2*yb*yb + totaldriftB)  * Apoints) * numblocks);

					//}						
					//long A0_B0_adj_inc = A0increase*numblocks;  

					//for (long int i=0;i<numblocks;i++){
					//
					//	long xpos = (long) (     A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + totaldriftA);
					//	long ypos = (long) (     B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + totaldriftB);
					//
					//			
					//	diffs_indices[i]=(i + xpos*numblocks + ypos*blocktimesApoints);
					//
					// }

					for (long B0_adj = B0_adj_min; B0_adj <= B0_adj_max; B0_adj += BO_adj_inc){
						long A0_B0_adj_max = MaxA0_times_numblocks + B0_adj;
						//for (double A0=-MaxA0;A0<=MaxA0;A0+=A0increase){ 
						//Eigen::Matrix4d A;
						//Eigen::MatrixXd B(4,1);

						//A *=0;
						//B *=0;
						A << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
						B << 0, 0, 0, 0;
						//Eigen::Matrix4d D;
						//Eigen::MatrixXd E(4,1);
						for (long A0_B0_adj = -MaxA0_times_numblocks + B0_adj; A0_B0_adj <= A0_B0_adj_max; A0_B0_adj += A0_B0_adj_inc){



							for (int i = 0; i<numblocks; i++){
								double Beta0, Beta1, Beta2, Beta3, Beta4, Beta5;

								//term0 = * (unsigned long long *)((char*)base_diffs_address+(* ((long long *) (((char*)diffs_indices)+i))));
								//term1 = * (unsigned long long *)((char*)base_diffs_address+(* ((long long *) (((char*)diffs_indices)+i+numBytes))));
								//sum0 += term0;
								//sum1 += term1;
								//j + A0_B0_adj;
								//if(i<30 && i %10 ==0){	
								long xpos = (long)(A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + totaldriftA);
								long ypos = (long)(B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + totaldriftB);
								//int index = diffs_indices[i] + A0_B0_adj;//*numBytes
								int index = i + xpos*numblocks + ypos*blocktimesApoints + A0_B0_adj;
								double yb = yb_arr[i];
								double yb2 = yb2_arr[i];
								double yb3 = yb3_arr[i];
								double yb4 = yb4_arr[i];
								double yb5 = yb5_arr[i];
								double yb6 = yb6_arr[i];

								Beta0 = dynamic_diffs[index].B0;
								Beta1 = dynamic_diffs[index].B1;
								Beta2 = dynamic_diffs[index].B2;
								Beta3 = dynamic_diffs[index].B3;
								Beta4 = dynamic_diffs[index].B4;
								Beta5 = dynamic_diffs[index].B5;
								Beta2 *= 2;
								//printf("\n%f: %f: %f: %f: %f: %f\n",Beta0,Beta1,Beta2,Beta3,Beta4,Beta5);
								//printf("\n%f: %f: %f: %f: %f: %f\n",yb,yb2,yb3,yb4,yb5,yb6);
								//getchar();

								/// D << 2*Beta2		 , 2*Beta2*yb+(Beta5)				   , 2*Beta2*yb2		   , 2*Beta2*yb3, //1st row
								//	2*Beta2*yb+(Beta5) , 2*Beta2*yb2 + (2*Beta4 + 2*Beta5*yb) , 2*Beta2*yb3 +(Beta5*yb2),2*Beta2*yb4+(Beta5*yb3), //2nd Row
								//	2*Beta2*yb2	 , 2*Beta2*yb3+(Beta5*yb2)		   , 2*Beta2*yb4		   ,2*Beta2*yb5, //3rd Row
								//	2*Beta2*yb3	 , 2*Beta2*yb4+(Beta5*yb3)		   , 2*Beta2*yb5		   ,2*Beta2*yb6; //4th row

								D << Beta2, Beta2*yb + (Beta5), Beta2*yb2, Beta2*yb3, //1st row
									Beta2*yb + (Beta5), Beta2*yb2 + (2 * Beta4 + 2 * Beta5*yb), Beta2*yb3 + (Beta5*yb2), Beta2*yb4 + (Beta5*yb3), //2nd Row
									Beta2*yb2, Beta2*yb3 + (Beta5*yb2), Beta2*yb4, Beta2*yb5, //3rd Row
									Beta2*yb3, Beta2*yb4 + (Beta5*yb3), Beta2*yb5, Beta2*yb6; //4th row


								//		D << 1 , yb, yb2 , yb3, //1st row
								//			yb , yb2 , yb3 ,yb4, //2nd Row
								//			yb2	 , yb3, yb4 ,yb5, //3rd Row
								//			yb3	 , yb4, yb5 ,yb6; //4th row
								//		D = D*Beta2;
								//		Eigen::Matrix4d G;
								//		F << 0, 1,0,0,
								//			1, 2*yb,(yb2),(yb3)
								//			,0,yb2,0,0,
								//			0,yb3,0,0;
								//		F = F*Beta5;
								//		G << 0,0,0,0,
								//			0,2*Beta4,0,0,
								//			0,0,0,0
								//			,0,0,0,0;
								//		F = F+G;
								//		D = D+F;

								E << Beta1,
									Beta1*yb + Beta3,
									Beta1*yb2,
									Beta1*yb3;

								A += D;
								B -= E;


								//	}

							}
							//cout << "\n----------------------------------\n";
							//cout << A <<"\n\n";

							//	cout << A <<"\n\n";
							//cout << B <<"\n\n";


							//Eigen::Matrix4d C = A.inverse()*B;
							C = A.inverse()*B;
							//cout << C <<"\n\n";
							//getchar();
							//cout << C;
							//					cout << "I'm here";
							//getchar();
							double C0 = C.data()[0];
							double C1 = C.data()[1];
							double C2 = C.data()[2];
							double C3 = C.data()[3];
							//printf("%f %e %e %e", C0,C1,C2,C3);
							//getchar();

							cs.lock();
							double function = 0;
							for (int i = 0; i<numblocks; i++){
								long xpos = (long)(A1*yb_arr[i] + A2*yb2_arr[i] + A3*yb3_arr[i] + totaldriftA);
								long ypos = (long)(B1*yb_arr[i] + B2*yb2_arr[i] + B3*yb3_arr[i] + totaldriftB);
								//int index = diffs_indices[i] + A0_B0_adj;//*numBytes
								int index = i + xpos*numblocks + ypos*blocktimesApoints + A0_B0_adj;
								//int index = diffs_indices[i] + A0_B0_adj;
								double yb = yb_arr[i];
								double yb2 = yb2_arr[i];
								double yb3 = yb3_arr[i];

								double Beta0 = dynamic_diffs[index].B0;
								double Beta1 = dynamic_diffs[index].B1;
								double Beta2 = dynamic_diffs[index].B2;
								double Beta3 = dynamic_diffs[index].B3;
								double Beta4 = dynamic_diffs[index].B4;
								double Beta5 = dynamic_diffs[index].B5;
								double delz = C0 + C1*yb + C2*yb2 + C3*yb3;
								function += Beta0 + Beta1*delz + Beta2*delz*delz
									+ Beta3*C1 + Beta4*C1*C1 + Beta5*C1*delz;
								//sums.local() += Beta0+ Beta1*delz + Beta2*delz*delz
								//+ Beta3*C1 + Beta4*C1*C1 + Beta5*C1*delz;
							}
							//function = sums.combine(plus<double>());
							//cs.lock();
							if (function <= bestdiff){  //note: this "if" block doesn't seem to slow down the program at all, compared to simpler statements.
								bestdiff = function;
								bestA0 = (A0_B0_adj - B0_adj) / numblocks;
								bestB0 = B0_adj / blocktimesApoints;
								bestA1 = A1;
								bestA2 = A2;
								bestA3 = A3;
								bestB1 = B1;
								bestB2 = B2;
								bestB3 = B3;
								bestC0 = C0;
								bestC1 = C1;
								bestC2 = C2;
								bestC3 = C3;


								//printf("BEST: %u, %f, %f, %e, %e, %f, %f, %f\n",bestdiff,bestdiffA0,A1,A2,A3,bestdiffB0,B1,B2);
								printf("BEST: %f   COUNT: %llu\n", bestdiff, count);
								printf("           %f, %f, %e, %e\n", bestA0, A1, A2, A3);
								//printf("           a1, a2, a3: %e, %e, %e\n",a1,a2,a3);
								printf("           %f, %f, %e, %e\n", bestB0, B1, B2, B3);
								printf("           %f, %f, %e, %e\n", bestC0, bestC1, bestC2, bestC3);
								//printf("           b1, b2, b3: %e, %e, %e\n",b1,b2,b3);
								//getchar();
								//if(bestA0==-1&&bestB0==-2&&bestA1==-.028&&bestB1==-.012&&bestA2==-.00028&&bestB2==-.000064){
							}

							count++;
							cs.unlock();

						}

					}
				}
				delete diffs_indices;
			});
		}
	}
	time3 = clock();
	printf("Time 0: \t %d \n", time0);
	printf("Time 1: \t %d \n", time1);
	printf("Time 2: \t %d \n", time2);
	printf("Time 3: \t %d \n", time3);
	printf("Time for creating diffs array: \t %f \n", (time1 - time0)*0.001);
	printf("Time for creating and sorting A and B param arrays: \t %f \n", (time2 - time1)*0.001);
	double grid_time = ((double)time3 - (double)time2) / 1000;
	printf("Time for grid search: \t %lf \n", grid_time);
	printf("count: %u\n", count);
	printf("Grid search calculations took %lf clockcycles per set of parameters.\n", (grid_time*3460000000.0) / count); //!~! clock line
	// HEIGHT*=precision;
	// bestC1=-10.0/HEIGHT;
	// bestC2=-15.0/(HEIGHT*HEIGHT);
	// bestC3=5.0/(long)(HEIGHT*HEIGHT*HEIGHT);

	bestC1 /= precision;
	bestC2 /= precision*precision;
	bestC3 /= precision*precision*precision;
	//bestC0=bestC1=bestC2=bestC3=0;


	// END COMMENTING TO SKIP GRID SEARCH HERE

	printf("Done grid search\n");
	getchar();

	int n = 12;         // number of parameters to be optimized 

	//now calculate the starting parameters for the simplex routine
	double x[12];

	/*
	// Values that grid search would find (If we didn't skip it)
	x[0]  = -4.000000;
	x[1]  = -3.000000;
	x[2]  = -6.376210e-003;
	x[3]  = 9.894744e-001;
	x[4]  = -7.493367e-005;
	x[5]  = -7.493367e-005;
	x[6]  = 5.498328e-008;
	x[7]  = 5.498328e-008;
	x[8]  = -1.037122e+004;
	x[9]  = -1.696923e+001;
	x[10] = -4.265604e-002;
	x[11] = 1.345831e-005;
	*/


	// These are how you should fill the x array. Use this initialization 
	// if you haven't skipped the grid search
	x[0] = bestA0 * precision;
	x[1] = bestB0 * precision;
	x[2] = bestA1;
	x[3] = bestB1 + 1;
	x[4] = bestA2 / precision;
	x[5] = bestB2 / precision;
	x[6] = bestA3 / (precision*precision);
	x[7] = bestB3 / (precision*precision);
	x[8] = bestC0;
	x[9] = bestC1;
	x[10] = bestC2;
	x[11] = bestC3;


	// bestA0 = -3.0;
	// bestA1=-15.0/HEIGHT;
	// bestA2=-10.0/(HEIGHT*HEIGHT);
	// bestA3=0/(HEIGHT*HEIGHT*HEIGHT);
	// bestB0 = -3;
	// bestB1=(-15.0/HEIGHT);
	// bestB2=-10.0/(HEIGHT*HEIGHT);
	// bestB3=0/(HEIGHT*HEIGHT*HEIGHT);
	// x[0] = bestA0 ;          
	// x[1] = bestB0;
	// x[2] = bestA1;          
	// x[3] = bestB1+1;
	// x[4] = bestA2;
	// x[5] = bestB2;
	// x[6] = bestA3;
	// x[7] = bestB3;
	// x[8] = bestC0;
	// x[9] = bestC1;
	// x[10] = bestC2;
	// x[11] = bestC3;


	//intentionally screw up the best guess, as a torture test for the simplex refinement routine
	//x[2] *=1.2
	//x[4] *=1.2
	//x[1] +=1


	// Setup precision for simplex.  
	// These values are all based on the changes in the parameters that will cause at most a shift of 
	// 0.1 pixel.
	double z[12];
	for (int i = 0; i<n; ++i) z[i] = x[i];
	/*
	printf("Coefficients of drift were\n");
	printf("Totaldiff: %d\n",bestdiff);
	printf("A0   %f\n",z[0]);
	printf("A1   %f   %f\n",z[2],z[2]*base_Image.height);
	printf("A2   %e   %f\n",z[4],z[4]*base_Image.height*base_Image.height);
	printf("A3   %e   %f\n",z[6],z[6]*base_Image.height*base_Image.height*base_Image.height);
	printf("B0   %f\n",z[1]);
	printf("B1   %f   %f\n",z[3]-1,(z[3]-1)*base_Image.height);
	printf("B2   %e   %f\n",z[5],z[5]*base_Image.height*base_Image.height);
	printf("B3   %e   %f\n",z[7],z[7]*base_Image.height*base_Image.height*base_Image.height);
	printf("C0: %e  %f\nC1: %e  %f\nC2: %e %f\nC3: %e % f\n", z[8],z[8],z[9],z[9]*base_Image.height,z[10],z[10]*base_Image.height*base_Image.height,z[11],z[11]*base_Image.height*base_Image.height*base_Image.height);
	*/
	//printf("time: %f\n",time2-time3);

	printf("Beginning simplex\n");
	double ycorrect;
	double precisionArr[12];
	precisionArr[0] = 0.1;
	precisionArr[1] = 0.1;
	precisionArr[2] = 0.1 / base_Image.height;
	precisionArr[3] = 0.1 / base_Image.height;
	precisionArr[4] = 0.1 / pow(base_Image.height, 2.0);
	precisionArr[5] = 0.1 / pow(base_Image.height, 2.0);
	precisionArr[6] = 0.1 / pow(base_Image.height, 3.0);
	precisionArr[7] = 0.1 / pow(base_Image.height, 3.0);
	precisionArr[8] = 0.1;///base_Image.height;
	precisionArr[9] = 0.1 / base_Image.height;
	precisionArr[10] = 0.1 / pow(base_Image.height, 2.0);
	precisionArr[11] = 0.1 / pow(base_Image.height, 3.0);


	//do a final warp of the image, according to the parameters from the grid search
	double aterms_n[4] = { z[0], z[2], z[4], z[6] };
	double bterms_n[4] = { z[1], z[3], z[5], z[7] };
	double cterms_n[4] = { z[8], z[9], z[10], z[11] };
	double ctermCopy[4] = { z[8], z[9], z[10], z[11] };
	image_basic final_norefine(base_Image.width, base_Image.height, base_Image.get_color_mode());
	//warp_image(&base_Image,base_Image.width, base_Image.height,aterms_n,bterms_n,1,&final_norefine);
	std::string finalstring = "C:\\Research Data\\nathan follin\\warping\\testing with idl\\images\\final_norefine.tif";
	//final_norefine.file_write_tiff(finalstring);

	simplex(&base_Image, &sliver_Image, x, z, n, precisionArr, reflectParam, contractParam, growthParam, haltParam, maxRefineIterationsParam);
	printf("Coefficients of drift were\n");
	printf("Totaldiff: %d\n", bestdiff);
	printf("A0   %f\n", z[0]);
	printf("A1   %f   %f\n", z[2], z[2] * base_Image.height);
	printf("A2   %e   %f\n", z[4], z[4] * base_Image.height*base_Image.height);
	printf("A3   %e   %f\n", z[6], z[6] * base_Image.height*base_Image.height*base_Image.height);
	printf("B0   %f\n", z[1]);
	printf("B1   %f   %f\n", z[3] - 1, (z[3] - 1)*base_Image.height);
	printf("B2   %e   %f\n", z[5], z[5] * base_Image.height*base_Image.height);
	printf("B3   %e   %f\n", z[7], z[7] * base_Image.height*base_Image.height*base_Image.height);
	printf("C0: %e  %f\nC1: %e  %f\nC2: %e %f\nC3: %e % f\n", z[8], z[8], z[9], z[9] * base_Image.height, z[10], z[10] * base_Image.height*base_Image.height, z[11], z[11] * base_Image.height*base_Image.height*base_Image.height);

	time4 = clock();
	printf("Time for simplex refinement: \t %f seconds\n", (time4 - time3)*0.001);

	//now do a final warp, and write the output tiff
	double aterms[4] = { z[0], z[2], z[4], z[6] };  //-1 is temporary, needed to fix an offset error.
	double bterms[4] = { z[1], z[3], z[5], z[7] };
	double cterms[4] = { z[8], z[9], z[10], z[11] };
	image_basic final(base_Image.width, base_Image.height, base_Image.get_color_mode());
	warp_image(&base_Image, base_Image.width, base_Image.height, aterms, bterms, cterms, 1, &final);
	int endtime = clock();
	//finalstring="C:\\Research Data\\trawick\\testing sliver warping\\set3\\redo\\final.tif";
	finalstring = rundata1;
	std::basic_string <char> finalmarkerstring("_corrected");
	finalstring.insert(finalstring.rfind("."), finalmarkerstring);

	//printf("output string:\n");
	//cout << finalstring << endl;

	final.file_write_tiff(finalstring);
	double totalelapsed = (double)(endtime - begintime) *.001;

	//getchar();

	//string namestring="C:\\Research Data\\trawick\\testing sliver warping\\idl_tests\\resampled.tif";
	//resamp_sliver_Image.file_write_tiff(namestring);

	FILE *logfile;
	HEIGHT = final.height;
	double a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3;
	a0 = 3;
	a1 = 15.0 / HEIGHT;
	a2 = 10.0 / (HEIGHT*HEIGHT);
	a3 = 0.0 / (HEIGHT*HEIGHT*HEIGHT);

	b0 = 3;
	b1 = 15.0 / HEIGHT;
	b2 = 10.0 / (HEIGHT*HEIGHT);
	b3 = 0 / (HEIGHT*HEIGHT*HEIGHT);

	c0 = 40.0;
	c1 = 45.0 / HEIGHT;
	c2 = 90.0 / (HEIGHT*HEIGHT);
	c3 = -20.0 / (HEIGHT*HEIGHT*HEIGHT);
	c0 *= 256;
	c1 *= 256;
	c2 *= 256;
	c3 *= 256;
	//c0=c1=c2=c3=0;
	printf("Trying to output log data\n\n");
	//size_t firstpart = rundata1.find_last_of("\\");
	//size_t endpart = rundata1.find("_s0");
	//string filename = rundata1.substr(firstpart,endpart) + ".dat";
	//printf("filename  %s\n\n",filename);
	std::string outputFileName = rundata1.substr(0, rundata1.length() - 4) + ".dat";//"C:\\Research Data\\nathan follin\\warping\\testing with idl\\images\\jan2011\\" + filename;
	logfile = fopen(outputFileName.c_str(), "w");
	FILE *timefile = fopen("runtimes.txt", "a");
	if (logfile == NULL){
		printf("file Not opened");
	}
	else{


		fprintf(logfile, "%d %e %e %e %e %e %e %e %e %e %e %e %e\n", final.height, a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3);
		fprintf(logfile, "%e %e %e %e %e %e %e %e %e %e %e %e\n", z[0], z[2], z[4], z[6], z[1], z[3] - 1.0, z[5], z[7], z[8], z[9], z[10], z[11]);
		fprintf(timefile, "Total Program Time =\t%f\n", totalelapsed, 'a');
		//fprintf(logfile,"%e %e %e %e %e %e %e %e %e %e %e %e\n" ,z[0],z[2], z[4], z[6], z[1],z[3]-1.0, z[5], ctermCopy[0],ctermCopy[0],ctermCopy[1],ctermCopy[2],ctermCopy[3]);
		printf("%e %e %e %e %e %e %e %e\n", z[0], z[2], z[4], z[6], z[1], z[3] - 1.0, z[5], z[7]);
		fclose(logfile);

		printf("------End of Program------\n");
	}
	if (argc < 33) getchar();  //puts a getchar in only if not run from commandline.
}


#else


#endif
