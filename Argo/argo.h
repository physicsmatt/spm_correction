#ifndef _ARGO_H_
#define _ARGO_H_
#include <string>
#include "FImage.h"

#define VERSION "1.3"
#define NUM_PARAMS 15

/**
 *	Struct for holding the 6 Beta values used to evaluate the C coefficients.
 */
struct beta_values {
		double B0, B1, B2, B3, B4, B5;
};

/**
 *	Struct for holding the 3 gamma values representing the 3 points on our parabola used in the creation of diffs.
 */
struct gamma_values {
		double G0, G1, G2;
};

/**
 *	Parameter Combo Struct.
 */
struct param_combo {
		double P1, P2, P3;
};

class argo {
	public:
		// Basic constructor.
		argo ();
		// Basic destructor.
		~argo ();

		// This struct hold various run-time values for debugging and profiling purposes.
		struct {
				double combos_time = 0, diffs_time = 0, grid_time = 0, simplex_time = 0, image_read_time = 0, image_write_time = 0, total_time = 0;
		} times;

		// File names for sliver and base images.
		std::string base_name, sliver_name;

		// Size of a block in pixels.
		int blocksize;

		// Maximum possible constant shift in pixels.
		int A0_max, B0_max;

		/**
		 *	Maximum linear shift, as fraction of image size
		 *	(if a1_Max = 0.02, and image height is 500 pixels,
		 *	this term causes max shift at top of image of 10 pixels)
		 */
		double a1_max, b1_max;

		/**
		 *	Fractional shift, relative to a1_max.
		 *	If a1_max implies max shift of 10 pixels at top, a2_Multiplier of 0.5 causes deviation
		 *	from linearity of 5 pixels at midpoint.
		 */
		double a2_mult, b2_mult;

		// Fractional shift for cubic term a3, also relative to shift caused by a1_max.
		double a3_mult, b3_mult;

		/**
		 *	For initial grid search, image is rescaled (down-sampled) by this factor.
		 *	This parameter should be set to the size, in pixels, of the smallest real feature on the image.
		 */
		int precision;

		/**
		 *	These parameters control the simplex routine.
		 */
		double simplex_growth, simplex_contract, simplex_reflect;
		int simplex_iterations;
		bool simplex_mode;		/* True for fast-z correction, false for slow-z correction. */

		/**
		 *	The dimensions of resampled images.
		 */
		int rsliver_width, rsliver_height, rbase_width, rbase_height, numblocks;

		/**
		 *	The maximum values for higher order parameters.
		 */
		double a2_max, a3_max, b2_max, b3_max;

		/**
		 *	The step values for higher order parameters.
		 */
		double a1_step, a2_step, a3_step, b1_step, b2_step, b3_step;

		long A_total_drift, B_total_drift;

		int A_sliver_drift, B_sliver_drift;
		int A_difflet_points, B_difflet_points;
		int A_points, B_points;
		int num_sliver_blocks;
		int sliver_block_width;

		unsigned long int dyn_diffs_size;
		long int A_combos_size, B_combos_size;

		// Range of the data that we have. abs(Minimum data point - Maximum data point)
		double data_range;

		/**
		 *	Struct to store original and resampled image files.
		 */
		struct {
				FImage * orig_sliver;
				FImage * orig_base;
				FImage * resamp_sliver;
				FImage * resamp_base;
				bool flipped;	// This boolean specifies if the image being read in should be read in y-down or y-up
		} images_store;

		/**
		 *	Struct to store the best computed A, B, and C parameters along with corresponding differences and grid-search debugging information.
		 */
		struct {
				double bestA0, bestA1, bestA2, bestA3;
				double bestB0, bestB1, bestB2, bestB3;
				double bestC0, bestC1, bestC2, bestC3;
				double bestdiff;

				// Grid search debugging values.
				unsigned long long count;
				unsigned long long iterations_ignored;
		} results;

		/**
		 *	Stores the beta and gamma diff arrays.
		 */
		struct {
				gamma_values * difflets;
				beta_values * dynamic_diffs;
		} beta_gamma_store;

		/**
		 *	Stores the range of precomputed A and B parameters grid search will search through.
		 */
		struct {
				param_combo * A_combos;
				param_combo * B_combos;
		} combos;

		// Number of parameters to be optimized
		double grid_best[ NUM_PARAMS ];
		double simplex_best[ NUM_PARAMS ];
		double precisionArr[ NUM_PARAMS ];

		/**
		 *	Debugging and non-debugging versions of argo execution.
		 */
		void correctImages ( int argc, char *argv[] );
		void correctImages ( int argc, char *argv[], bool verbose );

	private:
		// Deprecated variable for interpolation type. Will be refactored and removed eventually. Routine uses a bicubic/b-spline routine now.
		int interp_type = F_BSPLINE;

		/**
		 *	Various logging methods.
		 */
		void logInputParams ();
		void logCalculatedParams ();
		void logBetaGammaInfo ();
		void logCurrentBest ();
		void logGridSearchInfo ();
		void logSimplexRoutineInfo ();
		void logProgramInformation ();

		// Reads necessary parameters to run grid-search.
		void readInputParams ( int argc, char *argv[] );
		// Calculates values necessary for grid search run.
		void initCalculatedParams ();
		// Reads in image files.
		void readImages ( bool verbose );
		// Initializes A and B combos.
		void initCombos ();
		// Initializes Beta and Gamma diffs.
		void initBetaGamma ();
		// Performs the Grid Search routine.
		void performGridSearch ( bool verbose );
		// Deletes globals for Grid Search routine
		void finalizeGridSearch();
		// Performs the Simplex routine.
		void performSimplexRoutine ();
		// Performs the final image correction wtih the data calculated.
		void performImageCorrection ();
};

#endif _ARGO_H_
