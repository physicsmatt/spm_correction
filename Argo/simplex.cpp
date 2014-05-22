/*Here is an implementation in C.  Using it on Rosenbrock's famous function

myfunc(x,y)   =   (1-x)^2  +  100[ (x^2 - y)^2 ]

from the traditional starting point   (x,y) = (-1.2, +1.0)
this code calculated a minimum at     (x,y) = ( 1.000026  ,  1.000051 )
while using 177 evaluations of the function myfunc.
(It's easy to see that Rosenbrock's function myfunc has a minimum at (1,1)).

No warranty is expressed or implied, use at your own risk, remember that
this software is worth what you paid for it, and you paid zero.

=============================================================================
=============================================================================
*/

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>
#include "image_basic.h"
#include <iostream>
#include <time.h>
#include <vector>
#include "Eigen\Eigen"
#include <tbb\parallel_for.h>
#include <tbb\combinable.h>

//extern float interp_pixel_float(image_basic *baseImage, float x, float y);
extern double interp_pixel(image_basic *baseImage, double x, double y);
//extern double interp_pixelArray(/*vector<vector<double>> baseImage,*/int width,int sheight, double x, double y);

/* Nelder and Mead    Simplex for minimization of functions */

#define         VARS            (13)//(20)    /* # of parms to optimize */
//seems to me that VARS= 8 should be okay.  But it caused a buffer overrun error.  --Matt
#define         MAXI            (5000)  /* iteration limit        */
#define         ERRHALT         (1e-9)//(1e-9)  /* std error termination  */
//Setting this to 1e-10 stops it maybe 25 iterations sooner than 1e-11.  Big deal!  Might as well go for higher precision.
#define         PRINTEM         (500)//(10)    /* progress report period */
#define         ALP             (1.5)   /* reflection parameter   */
#define         BET             (0.7)   /* contraction parameter  */
#define         GAM             (2.1)   /* expansion parameter    */
#define         TINY            (1e-5)  /* first Simplex displace */

/* global variables */
long    funevals = 0;
double  p[VARS][VARS];
bool debug = false;
double lastfunc;
int operations = 0;
//global declaration
double Simplex[VARS][VARS];
double fvals[VARS];
double x_r[VARS];
double x_e[VARS];
double x_c[VARS];
double x_j[VARS];
//double x_l[VARS];
//double x_a[VARS];
//double x_h[VARS];
double centroid[VARS];
long double f_x_r, f_x_e, f_x_c, f_x_j, f_x_a;


int min2(int x, int y)
{
	if (x < y)
		return (x);
	else
		return (y);
}


int max2(int x, int y)
{
	if (x < y)
		return (y);
	else
		return (x);
}


/* Here is Rosenbrock's test function */
/*      the infamous parabolic valley or "banana function" */

image_basic *base;
image_basic *sliver;
int *basesize;
int *sliversize;
std::vector<std::vector<int> > base2;
std::vector<std::vector<int> > sliver2;

double interp_pixelArray(int width, int height, double x, double y, int ys) {

	double s, t;
	int image_width = width;
	int image_height = height;
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
	left_val = base2[left_index][top_index];
	right_val = base2[right_index][top_index];

	// Linear interpolation recoded to only one multiply
	top_val = s*right_val + (1 - s)*(left_val);//right_val + s * (left_val-right_val);

	// Interpolate across the bottom edge
	left_val = base2[left_index][bottom_index];
	right_val = base2[right_index][bottom_index];
	bottom_val = right_val + s * (left_val - right_val);
	//Matt changed above to below
	bottom_val = s* right_val + (1 - s) * (left_val);

	// Interpolate between top and bottom
	return (t*bottom_val + (1 - t)*top_val);//(bottom_val + t * (top_val-bottom_val));

}
double  m1yfunc(double x[])
{
	double  a, b;

	funevals++;             /* for informational purposes only */
	a = (1.0 - x[0]);
	b = 10.0 * ((x[0] * x[0]) - x[1]);
	return(((a*a) + (b*b)));

}

/* // Never gets used so commenting out for now
int rounds(double a) {
	if (a >= 0)
		return int(a + 0.5);
	else
		return int(a - 0.5);
}
*/

double myfunc(double z[])
{
	/*This is the function that is minimized by simplex.  What it does, conceptually, is similar to performing a
	polynomial warp on BOTH the base image and the sliver, taking the sum of differences squared between them (normalized
	by area), and returning that.  In fact what it does is a little faster and more efficient.

	First, it "warps" the sliver by calculating, for each point in the sliver, where that point would end up,
	IF the transformation were done in reverse.  Here the sliver points (xs, ys) are transformed to (xsp, ysp),
	for "prime".  It takes only the liniar portion for this transformation.

	Then, it "warps" the base image by doing the full 3rd order polynomial warp on (xsp,ysp) to yield (xspp, yspp).
	The actual value is gotten from the main image by interpolation.

	I tried writing this function using floats (instead of doubles) for some of the x and y coordinates, and found
	that the answer never converged beyond about 8 decimal places.  Plus, using floats instead of doubles didn't ACTUALLY
	make it ANY faster.  So doubles it is.
	*/
	funevals++;  //increment global variable used to count function evaluations

	double sum = 0.0;
	double area = 0.0;
	int basexstart = basesize[0] / 2 - sliversize[0] / 2;
	//	float ymin=0.0,ymax=(float) sliversize[1]-1;
	double ymin = 0.0, ymax = (double)sliversize[1] - 1;

	double C0 = z[8];
	double C1 = z[9];
	double C2 = z[10];
	double C3 = z[11];
	tbb::combinable<double> sums;
	tbb::combinable<double> areas;
	tbb::parallel_for(1, sliversize[1], 1, [&](int ys) { //think about whether to reverse these two loops.
		//sum+= C0*C0 + (C1-1)*(C1-1) + (C2-2)*(C2-2) +(C3-3)*(C3-3);
		//for (int ys=0; ys<sliversize[1];++ys) { //think about whether to reverse these two loops.
		for (int xs = 0; xs<sliversize[0]; ++xs) {   //fix so not starting at zero every time, only goes as far as needed
			//reverse the linear transform to fix sliver drift
			double weight;
			double ysp = (double)(ys + (z[3] - 1.0)*xs);  //be sure to get these signs correct!
			double xsp = (double)(xs + z[2] * xs);  //7/2/2008; I think this should be +.
			//double ysp = (double) (ys); //- (z[3]-1.0)*xs);  //be sure to get these signs correct!
			//double xsp = (double) (xs); //- z[2]*xs);  
			double ysp2 = ysp*ysp;
			double ysp3 = ysp2*ysp;
			//apply polynomial transform to account for image shift
			double xspp = (double)(xsp + z[0] + z[2] * ysp + z[4] * ysp2 + z[6] * ysp3);
			double yb = (double)(z[1] + z[3] * ysp + z[5] * ysp2 + z[7] * ysp3);
			//transform to base image coordinates
			double xb = basexstart + xspp;
			//see if new point falls on the existing base image;
			if ((yb >= ymin) && (yb <= ymax)) {
				weight = 1L; //
				if ((yb - ymin) < 1.0) weight = (yb - ymin);
				if ((ymax - yb) < 1.0) weight = (ymax - yb);
				//int Zs = sliver->get(rounds(xsp),rounds(ysp));
				int Zs = sliver->get(xs, ys);
				double Zm = interp_pixel(base, xb, yb);



				//double diff = sliver->get(xs,ys) - interp_pixel(base,xb,yb);
				double diff = (C0 + Zm + C1*yb + C2*yb*yb + C3*yb*yb*yb) - (Zs - C1*xsp);


				//	printf("MYFUNC VALS: %d %f %f\n", Zs, Zm, diff);
				//	getchar();

				//sum += diff*diff* weight;
				//area += weight;  //fix this so area is continuous function of x[...]
				sums.local() += diff*diff* weight;
				areas.local() += weight;
			}
		}
	});
	sum = sums.combine(std::plus<double>());
	area = areas.combine(std::plus<double>());
	/*	printf("at my function:\n");
	for(int j=0; j<8; j++) printf("    X(%3d) = %e\n",
	j, z[j]);

	*/
	return (sum / area);
}
void copyBase(image_basic *baseImage){
	base = new image_basic(baseImage->width, baseImage->height, 1);
	for (unsigned int i = 0; i<baseImage->width; i++)
	{
		for (unsigned int j = 0; j<baseImage->height; j++)
		{
			base->set(i, j, baseImage->get(i, j));
		}
	}
}

void copySliver(image_basic *sliverImage){
	sliver = new image_basic(sliverImage->width, sliverImage->height, 1);
	for (unsigned int i = 0; i<sliverImage->width; i++)
	{
		for (unsigned int j = 0; j<sliverImage->height; j++)
		{
			sliver->set(i, j, sliverImage->get(i, j));
		}
	}
}
int simplex(image_basic *baseImage, image_basic *sliverImage, double x[], double z[], int n,
	double precision[], double alpha, double beta, double gamma, double errhalt, int maxi)
{
	/*
	// Output values of x and z arrays (used to debug simplex)
	printf("Trying to output log data\n\n");
	std::string outputFileName = "output_x_array.txt";
	FILE *logfile = fopen(outputFileName.c_str(), "w");
	if (logfile == NULL)
	{
		printf("file Not opened");
	}
	else
	{
		for (int i = 0; i < 12; i++)
		{
			fprintf(logfile, "%e %f, %e %f\n", x[i], x[i], z[i], z[i]);
		}
		fclose(logfile);

		printf("------End of Program------\n");
	}
	

	// Output values of params (used to debug simplex)
	printf("Trying to output log data\n\n");
	std::string outputFileName3 = "output_params.txt";
	FILE *logfile3 = fopen(outputFileName3.c_str(), "w");
	if (logfile3 == NULL)
	{
		printf("file Not opened");
	}
	else
	{
		fprintf(logfile3, "alpha: %f\n", alpha);
		fprintf(logfile3, "beta:  %f\n", beta);
		fprintf(logfile3, "gamma: %f\n", gamma);
		fprintf(logfile3, "halt:  %f\n", errhalt);
		fprintf(logfile3, "maxi:  %i\n", maxi);
		fclose(logfile3);

		printf("------End of Program------\n");
	}
	*/

	printf("Begin Simplex Routine!\n");
	copyBase(baseImage);
	copySliver(sliverImage);
	sliversize = new int[2];
	sliversize[0] = sliver->width;
	sliversize[1] = sliver->height;
	basesize = new int[2];
	basesize[0] = base->width;

	/*

	// Output values of some random tests (used to debug simplex)
	printf("Trying to output log data\n\n");
	std::string outputFileName4 = "testoutput.txt";
	FILE *logfile4 = fopen(outputFileName4.c_str(), "w");
	if (logfile4 == NULL)
	{
		printf("file Not opened");
	}
	else
	{
		fprintf(logfile4, "BASE SIZE: %d x %d\n", base->width, base->height);
		fprintf(logfile4, "SLIVER SIZE: %d x %d\n", sliver->width, sliver->height);

		fprintf(logfile4, "EXTRA VALS: %d, %d, %d, %d\n", base->fast_get(0, 0), base->fast_get(0, 1), base->fast_get(1, 0), base->fast_get(1, 1));

		fprintf(logfile4, "VALUES FROM BASE:    %d, %d, %d, %d\n", base->fast_get(0, 0), base->fast_get(0, base->height - 1), base->fast_get(base->width - 1, 0), base->fast_get(base->width - 1, base->height - 1));
		fprintf(logfile4, "VALUES FROM SLIVER:    %d, %d, %d, %d\n", sliver->fast_get(0, 0), sliver->fast_get(0, sliver->height - 1), sliver->fast_get(sliver->width - 1, 0), sliver->fast_get(sliver->width - 1, sliver->height - 1));
		fprintf(logfile4, "VALUES FROM MYFUNC: %f\n", myfunc(z));
		fprintf(logfile4, "INTERP: %f", interp_pixel(base, 1.2, 3.3));
		fprintf(logfile4, "INTERP: %f", interp_pixel(base, 3.5, 4.3));
		fprintf(logfile4, "INTERP: %f", interp_pixel(base, 5.9, 7.1));
		fclose(logfile4);
	}
	
	*/


	basesize[1] = base->height;
	double t0;
	//int ok;	Never used
	int i, j, pass;
	int max, min;
	double c;
	double ymax;
	double ymin;
	double err;
	double meany;
	double y_a;
	n = 12;/*used to debug with smaller sample size*/

	/***************************************************/
	//fill the simplex
	for (i = 0; i<n; i++){
		Simplex[0][i] = x[i];
		printf("x[%d] = %e\n", i, x[i]);
	}
	for (i = n; i<VARS; i++){
		//				x_l[i]=0.0;
		x_r[i] = 0.0;
		//				x_a[i]=0.0;
		x_j[i] = 0.0;
		x_e[i] = 0.0;
		centroid[i] = 0.0;
		fvals[i] = 0.0;
	}
	/* fill in the remaining points of the initial Simplex            */
	/* point no. k is the initial point displaced by 2x along coord k */

	for (i = 1; i <= n; i++)
	{
		for (j = 0; j<n; j++)
		{
			if ((j + 1) != i)  Simplex[i][j] = Simplex[0][j];
			else    {
				t0 = Simplex[0][j];
				//if(abs(t0) < TINY) t0 = 0.5 ;
				Simplex[i][j] = precision[j] + t0;
			}
		}
	}

	/* go through all of the points of the Simplex & find max, min */
	max = 0; min = 0;
	fvals[0] = myfunc(x);
	ymax = fvals[0];  ymin = ymax;

	//ymin = F_L, ymax = F_A
	for (i = 1; i <= n; i++)
	{
		for (j = 0; j<n; j++) x_r[j] = Simplex[i][j];
		Simplex[i][8] = 0;
		fvals[i] = myfunc(x_r);
		if (fvals[i] >= ymax) {
			ymax = fvals[i];
			//x_h[i]=ymax;
			//		x_a[i] = fvals[max];
			max = i;
		}
		if (fvals[i] < ymin)  {
			ymin = fvals[i];
			//					x_l[i]=ymax;
			min = i;
		}
	}
	for (j = 0; j<n; j++)
	{
		centroid[j] = 0.00;
		for (i = 0; i <= n; i++)
		{
			if (i != max) centroid[j] += Simplex[i][j];
		}
		centroid[j] = centroid[j] / ((double)n);
	}
	meany = 0.00;
	for (i = 0; i <= n; i++){
		meany += fvals[i];
	}
	meany = meany / ((double)(n + 1));
	err = 0.00;
	for (i = 0; i <= n; i++){
		err += (meany - fvals[i]) * (meany - fvals[i]);

	}
	err = sqrt(err / ((double)(n + 1)));
	pass = 0;




	double myFuncResult = myfunc(z);


	while ((err >= errhalt) && (pass < maxi)){
		max = min = 0;
		debug = false;
		ymax = ymin = fvals[0];
		for (i = 1; i <= n; i++)
		{
			if (fvals[i] >= ymax) {
				ymax = fvals[i];
				//get the second highest value
				y_a = fvals[max];
				max = i;
			}
			/*if(fvals[i] < ymin)*/
			else  {
				ymin = fvals[i];
				min = i;
			}
		}

		/* compute the centroid */
		for (j = 0; j<n; j++)
		{
			centroid[j] = 0.00;
			for (i = 0; i <= n; i++)
			{
				if (i != max) centroid[j] += Simplex[i][j];
			}
			centroid[j] /= /*centroid[j]/ */ ((double)n);
		}

		if ( (pass % (PRINTEM)) == 0 )
		{
			printf("%d\n", operations);
			operations = 0;
			printf("\n ITERATION %4d     func = %20.12lf\n",
				pass, ymin);
			printf("%f %f %f %f \n", Simplex[min][0], Simplex[min][2], Simplex[min][8], Simplex[min][11]);
			//for(j=0; j<n; j++) printf("    X(%3d) = %-10.3lf\n",
			// for(j=0; j<n; j++) printf("    X(%3d) = %e\n",
			//j, Simplex[min][j]);
			printf("\n");
			//getchar();
		}
		/*new code*/
		///*******************************************************///
		///************************reflection*********************///
		///*******************************************************///
		//calculate the points of the reflected simplex
		for (i = 0; i<n; i++){
			c = centroid[i];
			x_r[i] = c + alpha*(c - Simplex[max][i]);
		}
		f_x_r = myfunc(x_r);
		//if the reflected simplex is decent, use it...but not if its very good.
		//if the reflected simplex is very good, try an expansion
		if (ymin <= f_x_r && f_x_r < y_a){
			for (i = 0; i<n; i++){ Simplex[max][i] = x_r[i]; }
			for (i = 0; i<n; i++) Simplex[max][i] = x_r[i];
			fvals[max] = f_x_r;
			operations++;
			//printf("R");
			debug = true;
			goto do_more;
		}
		///*******************************************************///
		///************************greedy expansion***************///
		///*********************non smooth function***************///	
		/*if(f_x_r < ymin){
		for(i=0; i<n; i++){
		c = centroid[i];
		x_e[i] = c + gamma*(x_r[i]-c);
		}
		f_x_e = myfunc(x_e);
		if(f_x_e < ymin){
		for(i =0; i<n; i++){Simplex[max][i]=x_e[i];}
		ymax = f_x_e ;
		printf("e");
		operations++;
		debug = true;
		goto do_more;
		}
		else{
		for(i =0; i<n; i++){Simplex[max][i]=x_r[i];}
		ymax = f_x_r ;
		printf("r");
		operations++;
		debug = true;
		goto do_more;
		}
		}*/

		///*******************************************************///
		///*********************greedy minimization***************///
		///************************smooth function****************///
		//if the reflected point produced a simplex with a smaller value than the current minimum,
		//check to see if an expansion produces an even better simplex
		if (f_x_r < ymin){		//never make this <= ~ Nathan
			for (int j = 0; j<n; j++){
				c = centroid[j];
				x_e[j] = c + gamma*(x_r[j] - c);
			}
			f_x_e = myfunc(x_e);
			//if the expanded simplex is better, use it.
			if (f_x_e < f_x_r){
				for (i = 0; i<n; i++){ Simplex[max][i] = x_e[i]; }
				ymax = f_x_e;
				//printf("e");
				operations++;
				debug = true;
				goto do_more;
			}
			//if the expanded simplex isn't better, use the reflected simplex, this will keep the
			//simplex smaller
			if (f_x_e >= f_x_r){
				for (i = 0; i<n; i++){ Simplex[max][i] = x_r[i]; }
				ymax = f_x_r;
				//printf("r");
				operations++;
				debug = true;
				goto do_more;
			}

		}
		/*Nathan's Contraction and Shrink code to include both contraction cases*/
		///*******************************************************///
		///************************Contraction********************///
		///*******************************************************///
		//try A contraction if the reflected value is greater than the second greatest value
		//(which from the other situations it actually has to be...if statement might not be 
		//needed)
		if (f_x_r >= y_a){
			//cout << "\nerror impending\nFXR: " << f_x_r << "\nFXA: " << y_a << "\nymax: " << ymax;
			bool le = true;
			if (f_x_r < ymax){
				le = true;
				for (i = 0; i < n; i++){
					double c = centroid[i];
					x_e[i] = c + (beta * (x_r[i] - c));
				}
			}
			if (f_x_r >= ymax){
				le = false;
				for (i = 0; i<n; i++){
					double c = centroid[i];
					x_e[i] = c + (beta * (Simplex[max][i] - c));
				}
			}

			f_x_e = myfunc(x_e);
			if (le == true){
				if (f_x_e <= f_x_r){
					//printf("c");
					operations++;
					debug = true;
					//cout << "case1\n";
					//getchar();
					for (i = 0; i<n; i++){
						Simplex[max][i] = x_e[i];
						fvals[max] = f_x_r;
						ymax = f_x_r;
					}
				}
				else
				{

					/* contraction failed; collapse simplex points */
					/* towards the present minimum point */
					//printf("S-c");
					operations++;
					debug = true;
					for (i = 0; i <= n; i++)
					{
						for (j = 0; j<n; j++)
						{
							Simplex[i][j] = 0.5*(Simplex[i][j] + Simplex[min][j]);
							x_r[j] = Simplex[i][j];
						}
						f_x_r = myfunc(x_r);
						fvals[i] = f_x_r;
					}
				}
			}
			if (le == false){
				if (f_x_e < ymax){
					//printf("C");
					operations++;
					debug = true;
					//cout << "case2\n";
					//getchar();
					for (i = 0; i<n; i++){
						Simplex[max][i] = x_e[i];
						fvals[max] = f_x_e;
						ymax = f_x_e;
					}
				}
				else
				{

					/* contraction failed; collapse simplex points */
					/* towards the present minimum point */
					//printf("S-C");
					operations++;
					debug = true;
					for (i = 0; i <= n; i++)
					{
						for (j = 0; j<n; j++)
						{
							Simplex[i][j] = 0.5*(Simplex[i][j] + Simplex[min][j]);
							x_r[j] = Simplex[i][j];
						}
						f_x_r = myfunc(x_r);
						fvals[i] = f_x_r;
					}
				}
			}
		}
		//end contraction


		/* compute the "standard error" of the Simplex */
		/* which is the termination criterion for the  */
		/* outermost (while) loop.                     */

	do_more:        meany = 0.00;
		for (i = 0; i <= n; i++) meany += fvals[i];
		meany /= /*meany /*/ ((double)(n + 1));
		err = 0.00;
		for (i = 0; i <= n; i++)
			err += (meany - fvals[i]) * (meany - fvals[i]);
		err = sqrt(err / ((double)(n + 1)));
		if (debug == false){
			std::cout << "\nfxc: " << f_x_c << "\nfxr: " << f_x_r << "\nymax: " <<
				ymax << "\nymin: " << ymin << "\nfxe: " << f_x_e << "\ny_a: "
				<< y_a << "\n";
			getchar();
		}
		pass++;
	}       /* end of while loop */

	/* find biggest and smallest values of myfunc */
	max = 0;  min = 0;
	ymax = fvals[0];  ymin = ymax;
	for (i = 1; i <= n; i++)
	{
		if (fvals[i] >= ymax) { ymax = fvals[i]; max = i; }
		if (fvals[i] < ymin)  { ymin = fvals[i]; min = i; }
	}

	/* put the minimum point in the return vector z */
	for (i = 0; i<n; i++) z[i] = Simplex[min][i];

	delete[] basesize;
	delete[] sliversize;
	delete base;
	delete sliver;
	printf("\n\nSIMPLEX took %ld evaluations; returned\n", funevals);
	printf("%f %f %f %f %f \n", Simplex[min][0], Simplex[min][2], Simplex[min][18], Simplex[min][19], Simplex[min][16]);

	/*
	// Output modified z array values (used to debug simplex)
	printf("Trying to output log data\n\n");
	std::string outputFileName2 = "output_z_array.txt";
	FILE *logfile2 = fopen(outputFileName2.c_str(), "w");
	if (logfile2 == NULL)
	{
		printf("file Not opened");
	}
	else
	{
		for (int i = 0; i < 12; i++)
		{
			fprintf(logfile2, "%e %f\n", z[i], z[i]);
		}
		fclose(logfile2);

		printf("------End of Program------\n");
	}

	*/

	return 1;
}
/* end of Simplex() */
