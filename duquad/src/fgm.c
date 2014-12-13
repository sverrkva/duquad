/*
 * FGM.c
 *
 *  Created on: Aug 21, 2014
 *      Author: Sverre
 */

#include "fgm.h"

// NOTE: z0 must be feasible

uint32_t FGM(struct Struct_FGM *s)
{

	s->exitflag = 1;

	// z = z0, y = z
	vector_copy(s->z0,s->z,N);
	vector_copy(s->z,s->y,N);

	// Initialize the function value, and difference between the new and previous function value
	real_t f = obj(s->z, s->H, s->c, s->temp1_dim_N);
	real_t fnew = 0.0;
	real_t f_diff = s->eps + 1;

	// Declaration and initialization of other necessary variables
	real_t L = s->eigH_max;
	real_t mu = s->eigH_min;
	real_t one_over_L = 1.0/L;
	real_t q = mu/L;
	real_t alpha = 0.5;
	real_t alpha_new = 0.0;
	real_t theta = 0.0;
	real_t beta = 0.0;
	real_t alpha_pow2 = 0.0;

	// iteration counter
	uint32_t niter = 0;


	/* ##################################################################################
			START WHILE-LOOP
	################################################################################## */
	while(f_diff > s->eps)
	{

		if (niter > s->maxiter){
			//fprintf(stderr, "reached maximum number of iterations in FGM (inner problem)\n");
			s->exitflag = 2;
			break;
		}

		/* *********************************************************************
	 			Finding znew
		********************************************************************* */

		// Calculating gradient f'(y) = H*y + c
		mtx_vec_mul(s->H,s->y,s->temp1_dim_N,N,N); // H * y
		vector_add(s->c,s->temp1_dim_N,s->temp1_dim_N,N); // + c

		// z_new = y - (1/L * f'(y))
		vector_scalar_mul(s->temp1_dim_N, one_over_L, s->temp1_dim_N, N); // 1/L * f'(y)
		vector_sub(s->y, s->temp1_dim_N, s->znew, N);	// y - ans

		// Projection znew on the feasible set made of lb and ub
		if(!s->ub_is_inf)
			vector_min(s->ub,s->znew,s->znew,N);
		if(!s->lb_is_inf)
			vector_max(s->lb,s->znew,s->znew,N);

		/* *********************************************************************
	 			Finding ynew
		********************************************************************* */

		// Calculating beta: beta = alpha(1-alpha) / alpha²+alpha_new,
		// where: alpha² = (1-alpha_new)*alpha² + q*alpha_new
		alpha_pow2 = alpha*alpha;
		theta = q - alpha_pow2;
		alpha_new = 0.5 * (theta + sqrt((theta*theta) + 4.0*alpha_pow2));
		beta = (alpha*(1.0-alpha)) / (alpha_pow2 + alpha_new);

		// ynew = znew + beta*(znew-z)
		vector_sub(s->znew,s->z,s->temp1_dim_N,N); // znew-z
		vector_scalar_mul(s->temp1_dim_N, beta, s->temp1_dim_N, N); // beta * (znew-z)
		vector_add(s->znew,s->temp1_dim_N,s->ynew,N); // + znew

		/* *********************************************************************
	 			Finding f_diff
		********************************************************************* */

		// Calculate the difference between the new and previous function value
		fnew = obj(s->znew, s->H, s->c, s->temp1_dim_N);
		f_diff = abs_2(fnew - f);

		/* *********************************************************************
	 			Update the variables
		********************************************************************* */

		vector_copy(s->znew,s->z,N);
		vector_copy(s->ynew,s->y,N);
		f = fnew;
		alpha = alpha_new;
		niter++;

	}	/* END WHILE-LOOP */

	vector_copy(s->znew,s->zopt,N);
	s->fopt = fnew;

	return niter;
}

// This function is not used by fgm or dfgm
void clean_up_FGM_C(struct Struct_FGM *s)
{
	free_pointer(s->H);
	free_pointer(s->c);
	free_pointer(s->lb);
	free_pointer(s->ub);
	free_pointer(s->z0);
	free_pointer(s->zopt);
	free_pointer(s->z);
	free_pointer(s->y);
	free_pointer(s->znew);
	free_pointer(s->ynew);
	free_pointer(s->temp1_dim_N);
}
