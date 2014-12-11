/*
 * falm.c
 *
 *  Created on: Nov 3, 2014
 *      Author: sverre
 */


#include "falm.h"

/* static functions declaration */

static int32_t solve_FALM();
static real_t dual_obj();
static void clean_up_inner_problem();
static void init_inner_problem();

int32_t FALM(struct Struct_FALM *s)
{
	_SDEBUG("I'm in FALM()\n");

	// Make z0 feasible
	if(!s->info->ub_is_inf)
		vector_min(s->prob->ub,s->prob->z0,s->prob->z0,N);
	if(!s->info->lb_is_inf)
		vector_max(s->prob->lb,s->prob->z0,s->prob->z0,N);


	// Initialize the "inner" problem
	struct Struct_FGM p_in;
	init_inner_problem(s, &p_in);
	if (solve_FALM(s,&p_in) == -1){
		ERROR("Error in solving solveFALM()\n");
		clean_up_inner_problem(&p_in);
		return -1;
	};
	clean_up_inner_problem(&p_in);
	_SDEBUG("End of FALM()\n");
	return 0;
}

static int32_t solve_FALM(struct Struct_FALM *s, struct Struct_FGM *p_in)
{
	// calculating alpha
	real_t alpha = s->opt->rho;

	// Must get right b from mex!!!!!!!!!!!!!!

	// Making constant vectors and matrices
	//A2 = A'*A; ???????
	//mtx_vec_mul(s->prob->A_t,s->prob->b,s->temp1_dim_N,N,M);
	//vector_scalar_mul(s->temp1_dim_N,rho,s->rho_At_b,N); // rho_At_b = rho * A' * b:

	// New H (H_hat). NOTE: when rho in not constant H_hat has to be calculated every time
	// Get H_hat from mex !!!!!


	// Dynamc Arrays for ds and pf
	struct Array ds_array, pf_array;
	initArray(&ds_array, 100);
	initArray(&pf_array, 100);

	// Declaration and initialization of other necessary variables
	real_t pf = s->opt->eps_pf + 1;
	uint32_t niter_feasible_ds = 0;
	uint32_t niter_feasible_pf = 0;
	uint32_t last_eps_ds = 0;
	uint32_t niter = 1;
	uint32_t niter_inner = 0;
	real_t dual_value_new = 0.0;
	real_t dual_value_diff = s->opt->eps_ds + 1;
	s->res->exitflag = 1;
	s->res->out->exitflag_inner = 1;
	s->res->out->num_exceeded_max_niter_inner = 0;
	// new:
	real_t theta = 1.0;
	real_t theta_new = 0.0;
	real_t beta = 0.0;
	real_t S = theta;

	// Clock
	clock_t tic, toc, tic_in, toc_in;
	tic = clock();

	// solve inner problem first time
	// with: c_in = c + A'*y1 - A'*y2 = c + A'(y1-y2) = c + A'(lambda-lambda2)
	// NOTE: As long as lambda = lambda2 = y1 = y2 0, c_in = c
	vector_sub(s->prob->c,s->rho_At_b,p_in->c,N);
	vector_copy(s->prob->z0,p_in->z0,N); // Initialize z0

	tic_in = clock();
	niter_inner += FGM(p_in);
	toc_in = clock();

	s->res->out->time_tot_inner = (real_t)(toc_in-tic_in);
	s->time_inner_y = (real_t)(toc_in-tic_in);
	if(p_in->exitflag == 2){
		s->res->out->exitflag_inner = 2;
		s->res->out->num_exceeded_max_niter_inner++;
	}
	vector_copy(p_in->zopt,s->z,N);
	vector_copy(p_in->zopt,s->z_ds,N);
	mtx_vec_mul(s->prob->A,s->z_ds,s->A_z_ds,M,N); // used in dual_obj
	vector_copy(s->z,s->summ,N);
	vector_scalar_mul(s->summ,(theta/S),s->summ,N);

	// find the value of the dual function (Lagrangian) first time
	real_t dual_value = dual_obj(s);

	/* ##################################################################################
			START WHILE-LOOP
	################################################################################## */

	while (dual_value_diff > s->opt->eps_ds || pf > s->opt->eps_pf)
	{

		if (niter > s->opt->maxiter_outer){
			//printf("reached maximum number of iterations in FALM\n");
			niter_feasible_ds--;
			niter_feasible_pf--;
			s->res->exitflag = 2;
			break;
		}
		if (dual_value_diff < s->opt->eps_ds){
			niter_feasible_ds++;
			last_eps_ds = 1;
		}
		else{
			last_eps_ds = 0;
		}
		if (pf < s->opt->eps_pf)
			niter_feasible_pf++;

		/* *********************************************************************
				Finding the next lambda
		********************************************************************* */

		// lambda += alpha * (A*z - b)
		// NOTE: re-use the 'lambda' vector instead of creating 'lambda_new'

		mtx_vec_mul(s->prob->A,s->z,s->A_z,M,N); 	// A*z

		// lambda = y1 + alpha*(A*z-b_ub_hat)
		vector_sub(s->A_z,s->prob->b,s->temp3_dim_M,M);	// ANS - b
		vector_scalar_mul(s->temp3_dim_M,alpha,s->temp3_dim_M,M);	// ANS * alpha
		vector_add(s->temp3_dim_M,s->y1,s->lambda,M);	// lambda = y1 + ANS
		// NOTE: No projection


	    /* *********************************************************************
				Finding the next y
		********************************************************************* */

		// y = lambda + beta * (lambda-lambda_old)

		// calculating beta
		theta_new = 0.5 * (1.0 + sqrt(1.0 + 4.0*(theta*theta)));
		beta = (theta-1.0)/theta_new;

		// y1
		vector_sub(s->lambda,s->lambda_old,s->temp2_dim_M,M); // lambda - lambda_old
		vector_scalar_mul(s->temp2_dim_M,beta,s->temp2_dim_M,M); // ANS * beta
		vector_add(s->lambda,s->temp2_dim_M,s->y1,M); 	// y1 = ANS + lambda


		/* *********************************************************************
				Solving the inner problem
		********************************************************************* */

		// First calculate the new c (c_hat) for the inner problem,
		// c_hat = p_in->c = c + A' * (y1_new - y2_new)

		vector_copy(s->prob->c,p_in->c,N);

		// p_in->c = c + A' * y1
		mtx_vec_mul(s->prob->A_t,s->y1,s->temp1_dim_N,N,M); // A' * y1
		vector_add(p_in->c,s->temp1_dim_N,p_in->c,N);	// p_in->c = ANS + c
		vector_sub(p_in->c,s->rho_At_b,p_in->c,N);	// - rho *A'*b

		// Compute the new optimal z (s->z) from the inner problem
		vector_copy(s->z, p_in->z0, N); // Warm start

		tic_in = clock();
		niter_inner += FGM(p_in);
		toc_in = clock();
		s->res->out->time_tot_inner += (real_t)(toc_in-tic_in);
		s->time_inner_y += (real_t)(toc_in-tic_in);
		if(p_in->exitflag == 2){
			s->res->out->exitflag_inner = 2;
			s->res->out->num_exceeded_max_niter_inner++;
		}
		vector_copy(p_in->zopt,s->z,N);

		/* *********************************************************************
				Finding dual_value_diff
		********************************************************************* */

		// calculate the difference between the new and previousdual function value (ds)
		// NOTE: must calculate the 'inner problem' one more time, with lambda instead of y

		vector_copy(s->prob->c,p_in->c,N);


		// p_in->c = c + A' * lambda
		mtx_vec_mul(s->prob->A_t,s->lambda,s->temp1_dim_N,N,M); // A' * lambda
		vector_add(p_in->c,s->temp1_dim_N,p_in->c,N);	// p_in->c = ANS + c
		vector_sub(p_in->c,s->rho_At_b,p_in->c,N);	// - rho *A'*b


		 // Finding new z_ds
		vector_copy(s->z_ds, p_in->z0, N); // Warm start
		tic_in = clock();
		niter_inner += FGM(p_in);
		toc_in = clock();
		s->res->out->time_tot_inner += (real_t)(toc_in-tic_in);

		if(p_in->exitflag == 2){
			s->res->out->exitflag_inner = 2;
			s->res->out->num_exceeded_max_niter_inner++;
		}
		vector_copy(p_in->zopt,s->z_ds,N);
		mtx_vec_mul(s->prob->A,s->z_ds,s->A_z_ds,M,N); // used in dual_obj

		// Calculate dual_value_new
		dual_value_new = dual_obj(s);	// z_ds is used in this function
		dual_value_diff = abs_2(dual_value_new - dual_value);
		dual_value = dual_value_new;

		/* *********************************************************************
				Finding pf
		********************************************************************* */

		// pf = norm2(  max([A*z_pf - b_ub_hat ; -A*z_pf + b_lb_hat]) , 0)  )
		// NOTE: Algorithm 3 (last) => we use z_pf = z_ds
		// 		 Algorithm 4 (average) => we use z_pf = z_avg = summ_k / (niter+1)

		if (s->opt->algorithm == 7) {
		// *** LAST z in stopping criteria ***

			vector_sub(s->A_z_ds,s->prob->b,s->pf_vec,M); // (A*z) - b_ub_hat
			// NOTE: No projection

		}
		else if (s->opt->algorithm == 8) {
		// *** AVERAGE z in pf stopping criteria ***

			// Calculating z_avg
			S = S + theta_new;
			vector_scalar_mul(s->z,theta_new,s->temp1_dim_N,N); // S * theta_new
			vector_add(s->summ,s->temp1_dim_N,s->summ,N); // summ += ANS
			vector_scalar_mul(s->summ,1.0/S,s->z_avg,N); //z_avg = summ / S
			mtx_vec_mul(s->prob->A,s->z_avg,s->temp2_dim_M,M,N); // A * z_pf = A * z_avg

			vector_sub(s->temp2_dim_M,s->prob->b,s->pf_vec,M); // (A*z) - b
			// NOTE: No projection

		}
		else {
			ERROR("Choose algorithm 7 or 8 when running FALM()\n");
    		return -1;
		}

	    // Take the norm
	    pf = vector_norm_2(s->pf_vec, s->info->pf_vec_length); // norm2

		// store the result of the stopping criteria
		insertArray(&ds_array, dual_value_diff);
		insertArray(&pf_array, pf);

		/* *********************************************************************
			 			Update the variables
		********************************************************************* */

		vector_copy(s->lambda,s->lambda_old,M);
		//vector_copy(s->lambda2,s->lambda2_old,M);
		theta = theta_new;
		niter++;

	} 	/* #### END WHILE-LOOP #### */


	/* *********************************************************************
				 		Setting the result
	********************************************************************* */

	// Clock
	toc = clock();
	s->res->out->time = (real_t)(toc - tic) / CLOCKS_PER_SEC;
	s->res->out->time_tot_inner /=  CLOCKS_PER_SEC;
	s->time_inner_y /=  CLOCKS_PER_SEC;
	//printf("time inner y: %f\n", s->time_inner_y);

	if (s->opt->algorithm == 7) {
		s->res->zopt = s->z_ds;
		s->res->fopt = obj(s->z_ds, s->prob->H, s->prob->c, s->temp1_dim_N);
	}
	else if (s->opt->algorithm == 8) {
		s->res->zopt = s->z_avg;
		s->res->fopt = obj(s->z_avg, s->prob->H, s->prob->c, s->temp1_dim_N);
	}
	else {
		ERROR("Choose algorithm 7 or 8 when running FALM()\n");
		return -1;
	}

	s->res->lambda1 = s->lambda;
	//s->res->lambda2 = s->lambda2;
	s->res->out->iterations = --niter;
	s->res->out->iterations_inner_tot = niter_inner;
	s->res->out->niter_feasible_ds = ++niter_feasible_ds;
	s->res->out->niter_feasible_pf = ++niter_feasible_pf;

	if (last_eps_ds == 1)
		s->res->out->flag_last_satisfied = 2;
	else
		s->res->out->flag_last_satisfied = 1;

	// Copy ds_vector and pf_vector
	s->res->out->ds_vector = (real_t *)realloc(s->res->out->ds_vector, (niter) * sizeof(real_t));
	s->res->out->pf_vector = (real_t *)realloc(s->res->out->pf_vector, (niter) * sizeof(real_t));
	vector_copy(ds_array.array,s->res->out->ds_vector,niter);
	vector_copy(pf_array.array,s->res->out->pf_vector,niter);
	freeArray(&ds_array);
	freeArray(&pf_array);

	return 0;
}


static real_t dual_obj(struct Struct_FALM *s)
{
	/*return: 0.5z'Hz + c'z + lambda'*(A*z - b) + 0.5*rho*norm(A*z-b) */

	// f_value = 0.5*z'Hz + c'z;
	real_t f_value = obj(s->z_ds, s->prob->H, s->prob->c, s->temp1_dim_N);

	// f_value += lambda'*(A*z-b)

	vector_sub(s->A_z_ds,s->prob->b,s->temp3_dim_M,M); // (A*z) - b
	f_value += vector_mul(s->lambda,s->temp3_dim_M,M);

	real_t nor = vector_norm_2(s->temp3_dim_M,M);
	f_value += 0.5 * s->opt->rho * (nor*nor);

	return f_value;
}

static void init_inner_problem(const struct Struct_FALM *s, struct Struct_FGM *p_in)
{
	// Setting the inner problem to point to the outer problem except for c and z0
	// Setting the inner problem to point to the same as the outer saves some memory
	p_in->H = s->H_hat;
	p_in->lb = s->prob->lb;
	p_in->ub = s->prob->ub;
    p_in->lb_is_inf = s->info->lb_is_inf;
	p_in->ub_is_inf = s->info->ub_is_inf;
	p_in->eigH_max = s->info->eigH_max;
	p_in->eigH_min = s->info->eigH_min;

	// Make own allocation for c (This is changing for each outer iteration)
	p_in->c = vector_alloc(N);
	p_in->z0 = vector_alloc(N);

	// Allocate other necessary vectors
	p_in->zopt = vector_alloc(N);
	p_in->z = vector_alloc(N);
	p_in->y = vector_alloc(N);
	p_in->znew = vector_alloc(N);
	p_in->ynew = vector_alloc(N);
	p_in->temp1_dim_N = vector_alloc(N);

	// Initialize options
	p_in->maxiter = s->opt->maxiter_inner;
	p_in->eps = s->opt->eps_inner;
}

static void clean_up_inner_problem(struct Struct_FGM *s)
{
	free_pointer(s->c);
	free_pointer(s->z0);
	free_pointer(s->zopt);
	free_pointer(s->z);
	free_pointer(s->y);
	free_pointer(s->znew);
	free_pointer(s->ynew);
	free_pointer(s->temp1_dim_N);
}
