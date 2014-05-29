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
#include <iostream>
#include <time.h>
#include <vector>
#include "Eigen\Eigen"
#include <tbb\parallel_for.h>
#include <tbb\combinable.h>
#include "simplex.h"
#include "FImage.h"

/* Nelder and Mead    Simplex for minimization of functions */

#define         VARS            (13)//(20)    /* # of parms to optimize */
//seems to me that VARS= 8 should be okay.  But it caused a buffer overrun error.  --Matt
//#define         MAXI            (5000)  /* iteration limit        */
//#define         ERRHALT         (1e-9)//(1e-9)  /* std error termination  */
//Setting this to 1e-10 stops it maybe 25 iterations sooner than 1e-11.  Big deal!  Might as well go for higher precision.
#define         PRINTEM         (500)//(10)    /* progress report period */
//#define         ALP             (1.5)   /* reflection parameter   */
//#define         BET             (0.7)   /* contraction parameter  */
//#define         GAM             (2.1)   /* expansion parameter    */
//#define         TINY            (1e-5)  /* first Simplex displace */

/* global variables */
long    function_evals = 0;
//double  p[VARS][VARS];
bool debug = false;
//double lastfunc; Never used.
int operations = 0;
//global declaration
double simplex_matrix[VARS][VARS];
double fvals[VARS];
double x_r[VARS];
double x_e[VARS];
//double x_c[VARS];
//double x_j[VARS];
//double x_l[VARS];
//double x_a[VARS];
//double x_h[VARS];
double centroid[VARS];
long double f_x_r, f_x_e, f_x_c;	// f_x_j, f_x_a;

/* Here is Rosenbrock's test function */
/*      the infamous parabolic valley or "banana function" */

FImage *base;
FImage *sliver;
int *basesize;
int *sliversize;

double  m1yfunc(double x[])
{
	double  a, b;

	function_evals++;             /* for informational purposes only */
	a = (1.0 - x[0]);
	b = 10.0 * ((x[0] * x[0]) - x[1]);
	return(((a*a) + (b*b)));
}

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
double myfunc(double z[])
{
	function_evals++;  //increment global variable used to count function evaluations

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
		for (int xs = 0; xs < sliversize[0]; ++xs) {   //fix so not starting at zero every time, only goes as far as needed
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
				double Zs = sliver->get(xs, ys);
				double Zm = base->interpPixel(xb, yb);

				//double diff = sliver->get(xs,ys) - interp_pixel(base,xb,yb);
				double diff = (C0 + Zm + C1*yb + C2*yb*yb + C3*yb*yb*yb) - (Zs - C1*xsp);

				//	printf("MYFUNC VALS: %d %f %f\n", Zs, Zm, diff);
				//	getchar();

				//sum += diff*diff* weight;
				//area += weight;  //fix this so area is continuous function of x[...]
				sums.local() += diff * diff * weight;
				areas.local() += weight;
			}
		}
	});
	sum = sums.combine(std::plus<double>());
	area = areas.combine(std::plus<double>());
	return (sum / area);
}
void copyBase(FImage *_base){
	base = new FImage(_base->width, _base->height, _base->metadata);
	for (unsigned int i = 0; i < _base->width; i++)
	{
		for (unsigned int j = 0; j < _base->height; j++)
		{
			base->set(i, j, _base->get(i, j));
		}
	}
}

void copySliver(FImage *_sliver){
	sliver = new FImage(_sliver->width, _sliver->height, _sliver->metadata);
	for (unsigned int i = 0; i < _sliver->width; i++)
	{
		for (unsigned int j = 0; j < _sliver->height; j++)
		{
			sliver->set(i, j, _sliver->get(i, j));
		}
	}
}
int simplex(FImage *_base, FImage *_sliver, double input_vector[], double return_vector[], int n,
	double precision[], double alpha, double beta, double gamma, double errhalt, int maxi)
{
	printf("Begin Simplex Routine!\n");
	copyBase(_base);
	copySliver(_sliver);
	sliversize = new int[2];
	sliversize[0] = sliver->width;
	sliversize[1] = sliver->height;
	basesize = new int[2];
	basesize[0] = base->width;
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

	/***************************************************/
	//fill the simplex
	for (i = 0; i < n; i++){
		simplex_matrix[0][i] = input_vector[i];
	}
	for (i = n; i < VARS; i++){
		//				x_l[i]=0.0;
		x_r[i] = 0.0;
		//				x_a[i]=0.0;
		//				x_j[i] = 0.0;
		x_e[i] = 0.0;
		centroid[i] = 0.0;
		fvals[i] = 0.0;
	}
	/* fill in the remaining points of the initial Simplex            */
	/* point no. k is the initial point displaced by 2x along coord k */

	for (i = 1; i <= n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if ((j + 1) != i)  simplex_matrix[i][j] = simplex_matrix[0][j];
			else    {
				t0 = simplex_matrix[0][j];
				//if(abs(t0) < TINY) t0 = 0.5 ;
				simplex_matrix[i][j] = precision[j] + t0;
			}
		}
	}

	/* go through all of the points of the Simplex & find max, min */
	max = 0; min = 0;
	fvals[0] = myfunc(input_vector);
	ymax = fvals[0];  ymin = ymax;

	//ymin = F_L, ymax = F_A
	for (i = 1; i <= n; i++)
	{
		for (j = 0; j < n; j++) x_r[j] = simplex_matrix[i][j];
		simplex_matrix[i][8] = 0;
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
	for (j = 0; j < n; j++)
	{
		centroid[j] = 0.00;
		for (i = 0; i <= n; i++)
		{
			if (i != max) centroid[j] += simplex_matrix[i][j];
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

	double myFuncResult = myfunc(return_vector);

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
		for (j = 0; j < n; j++)
		{
			centroid[j] = 0.00;
			for (i = 0; i <= n; i++)
			{
				if (i != max) centroid[j] += simplex_matrix[i][j];
			}
			centroid[j] /= /*centroid[j]/ */ ((double)n);
		}

		if ((pass % (PRINTEM)) == 0)
		{
			printf("%d\n", operations);
			operations = 0;
			printf("\n ITERATION %4d     func = %20.12lf\n",pass, ymin);
			printf("%f %f %f %f \n", simplex_matrix[min][0], simplex_matrix[min][2], simplex_matrix[min][8], simplex_matrix[min][11]);
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
		for (i = 0; i < n; i++){
			c = centroid[i];
			x_r[i] = c + alpha*(c - simplex_matrix[max][i]);
		}
		f_x_r = myfunc(x_r);
		//if the reflected simplex is decent, use it...but not if its very good.
		//if the reflected simplex is very good, try an expansion
		if (ymin <= f_x_r && f_x_r < y_a){
			for (i = 0; i < n; i++){ simplex_matrix[max][i] = x_r[i]; }
			for (i = 0; i < n; i++) simplex_matrix[max][i] = x_r[i];
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
			for (int j = 0; j < n; j++){
				c = centroid[j];
				x_e[j] = c + gamma*(x_r[j] - c);
			}
			f_x_e = myfunc(x_e);
			//if the expanded simplex is better, use it.
			if (f_x_e < f_x_r){
				for (i = 0; i < n; i++){ simplex_matrix[max][i] = x_e[i]; }
				ymax = f_x_e;
				//printf("e");
				operations++;
				debug = true;
				goto do_more;
			}
			//if the expanded simplex isn't better, use the reflected simplex, this will keep the
			//simplex smaller
			if (f_x_e >= f_x_r){
				for (i = 0; i < n; i++){ simplex_matrix[max][i] = x_r[i]; }
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
				for (i = 0; i < n; i++){
					double c = centroid[i];
					x_e[i] = c + (beta * (simplex_matrix[max][i] - c));
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
					for (i = 0; i < n; i++){
						simplex_matrix[max][i] = x_e[i];
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
						for (j = 0; j < n; j++)
						{
							simplex_matrix[i][j] = 0.5*(simplex_matrix[i][j] + simplex_matrix[min][j]);
							x_r[j] = simplex_matrix[i][j];
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
					for (i = 0; i < n; i++){
						simplex_matrix[max][i] = x_e[i];
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
						for (j = 0; j < n; j++)
						{
							simplex_matrix[i][j] = 0.5*(simplex_matrix[i][j] + simplex_matrix[min][j]);
							x_r[j] = simplex_matrix[i][j];
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
	for (i = 0; i < n; i++) return_vector[i] = simplex_matrix[min][i];

	delete[] basesize;
	delete[] sliversize;
	delete base;
	delete sliver;
	printf("\n\nSIMPLEX took %ld evaluations; returned\n", function_evals);
	printf("%f %f %f %f %f \n", simplex_matrix[min][0], simplex_matrix[min][2], simplex_matrix[min][18], simplex_matrix[min][19], simplex_matrix[min][16]);

	return (1);
}
/* end of Simplex() */
