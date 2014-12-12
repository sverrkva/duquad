/*
 * DGM.h
 *
 *  Created on: 27. aug. 2014
 *      Author: Sverre
 */

#ifndef DGM_H_
#define DGM_H_

#include "head.h"
#include "math_functions.h"
#include "fgm.h"
#include "print.h"
#include "init_problem.h"


struct Struct_DGM {

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
};

int32_t DGM();


#endif /* DGM_H_ */
