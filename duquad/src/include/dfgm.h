/*
 * DFGM.h
 *
 *  Created on: Sep 19, 2014
 *      Author: sverre
 */

#ifndef DFGM_H_
#define DFGM_H_

#include "head.h"
#include "math_functions.h"
#include "fgm.h"
#include "print.h"

struct Struct_DFGM {

	struct Problem * prob;
	struct Options * opt;
	struct Info * info;
	struct Result * res;

	// Vectors used for calculations
	real_t * z;
	real_t * lambda1;
	real_t * lambda2;
	real_t * temp1_dim_N;
	real_t * temp2_dim_M;
	real_t * temp3_dim_M;
	real_t * b_ub_hat;	// b + ub_hat
	real_t * b_lb_hat;	// b + lb_hat
	real_t * z_avg;
	real_t * summ;
	real_t * pf_vec;
	real_t * A_z;

	// Different from DGM
	real_t * lambda1_old;
	real_t * lambda2_old;
	real_t * y1;
	real_t * z_ds;
	real_t * y2;
	real_t * A_z_ds;

	real_t time_inner_y;
	uint32_t iterations_inner_y;

};

int32_t DFGM();

#endif /* DFGM_H_ */
