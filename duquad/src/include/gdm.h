/*
 * GDM.h
 *
 *  Created on: Aug 19, 2014
 *      Author: Sverre
 */

#ifndef GRADIENT_DESCENT_H_
#define GRADIENT_DESCENT_H_

#include "head.h"
#include "math_functions.h"
#include "general_functions.h"

struct Struct_GDM
{
	real_t * H;
	real_t * c;
	real_t * lb;
	real_t * ub;
	real_t * z0;

	// Optimal outputs
	real_t * zopt;
	real_t fopt;
	uint32_t exitflag;

	// INFO
	boolean lb_is_inf;
	boolean ub_is_inf;
	real_t eigH_max;
	real_t eigH_min;

	// Vectors used for calculations
	real_t * z;

	real_t * znew;

	real_t * temp1_dim_N;

	// Options
	uint32_t maxiter;
	real_t eps;
};

// Functions
uint32_t GDM();
void clean_up_GDM_C();

#endif /* GRADIENT_DESCENT_H_ */
