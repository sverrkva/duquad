/*
 * alm.h
 *
 *  Created on: Nov 3, 2014
 *      Author: sverre
 */

/** \file
 *  ### Augmented Lagrangian Method ###
 */

#ifndef ALM_H_
#define ALM_H_

#include "head.h"
#include "math_functions.h"
#include "fgm.h"
#include "print.h"


struct Struct_ALM {

	struct Problem * prob;
	struct Options * opt;
	struct Info * info;
	struct Result * res;

	// Vectors used for calculations
	real_t * z;
	real_t * lambda;
	real_t * temp1_dim_N;
	real_t * temp2_dim_M;
	real_t * temp3_dim_M;
	real_t * z_avg;
	real_t * summ;
	real_t * pf_vec;
	real_t * A_z;

	// New for alm
	real_t * H_hat;
	real_t * A2;
	real_t * rho_At_b;

}; /**< Struct containing all necessary vectors and parameters for running ALM */

// Public functions
int32_t ALM(struct Struct_ALM *s);

#endif /* ALM_H_ */
