/*
 * fast_gradient.h
 *
 *  Created on: Aug 21, 2014
 *      Author: Sverre
 */

#ifndef FAST_GRADIENT_H_
#define FAST_GRADIENT_H_

#include "head.h"
#include <math.h>
#include "math_functions.h"
#include "print.h"

struct Struct_FGM
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
	real_t * y;
	real_t * znew;
	real_t * ynew;
	real_t * temp1_dim_N;

	// Options
	uint32_t maxiter;
	real_t eps;

};

uint32_t FGM();
void clean_up_FGM_C();

#endif /* FAST_GRADIENT_H_ */
