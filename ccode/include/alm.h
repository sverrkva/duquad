/*
 * alm.h
 *
 *  Created on: Nov 3, 2014
 *      Author: sverre
 */

#ifndef ALM_H_
#define ALM_H_

#include "head.h"
#include "math_functions.h"
#include "fgm.h"
#include "print.h"
#include "init_problem.h"


struct Struct_ALM {

	struct Problem * prob;
	struct Options * opt;
	struct Info * info;
	struct Result * res;

	// Vectors used for calculations
	real_t * z;
	real_t * lambda;
	//real_t * lambda2;
	real_t * temp1_dim_N;
	real_t * temp2_dim_M;
	real_t * temp3_dim_M;
	//real_t * b_ub_hat;	// b + ub_hat
	//real_t * b_lb_hat;	// b + lb_hat
	real_t * z_avg;
	real_t * summ;
	real_t * pf_vec;
	real_t * A_z;

	// New for alm
	real_t * H_hat;
	real_t * A2;
	real_t * rho_At_b;

};

int32_t ALM();


#endif /* ALM_H_ */
