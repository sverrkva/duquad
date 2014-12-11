#include "gdm.h"

// This function is not maintained

//static void gradientGD();

uint32_t GDM(struct Struct_GDM *s)
{

	s->exitflag = 1;

	// z = z0
	vector_copy(s->z0,s->z,N);


	// Initialize the function value, and difference between the new and previous function value
	real_t f = obj(s->z, s->H, s->c, s->temp1_dim_N);
	real_t fnew = 0.0;
	real_t f_diff = s->eps + 1;

	// Declaration and initialization of other necessary variables
	real_t L = s->eigH_max;
	// real_t mu = s->eigH_min;
	real_t one_over_L = 1.0/L;

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

		// calculating gradient
		//gradientGD(s->z,s->H,s->c,s->temp1_dim_N);
		mtx_vec_mul(s->H,s->z,s->temp1_dim_N,N,N); // H * y
		vector_add(s->c,s->temp1_dim_N,s->temp1_dim_N,N); // + c



		/* *********************************************************************
	 			Finding znew
		********************************************************************* */
//		for(i=0;i<N;i++){
//			s->znew[i] = s->z[i] - (opt->alpha * s->temp1_dim_N[i]);
//		}

		vector_scalar_mul(s->temp1_dim_N, one_over_L, s->temp1_dim_N, N); // 1/L * f'(y)
		vector_sub(s->z, s->temp1_dim_N, s->znew, N);	// y - ans



		// projection znew on the feasible set made of lb and ub
		if(!s->ub_is_inf)
			vector_min(s->ub,s->znew,s->znew,N);
		if(!s->lb_is_inf)
			vector_max(s->lb,s->znew,s->znew,N);

		/* *********************************************************************
	 			Finding f_diff
		********************************************************************* */

		// calculate the difference between the new and previous function value
		fnew = obj(s->znew, s->H, s->c, s->temp1_dim_N);
		f_diff = abs_2(fnew - f);

		/* *********************************************************************
	 			Update the variables
		********************************************************************* */
//		for(i=0;i<N;i++)
//			s->z[i] = s->znew[i];
		vector_copy(s->znew,s->z,N);
		f = fnew;
		niter = niter+1;
	}

	vector_copy(s->znew,s->zopt,N);
	s->fopt = fnew;

	return niter;
}


// This function is not used by fgm or dfgm
void clean_up_GDM_C(struct Struct_GDM *s)
{
	free_pointer(s->H);
	free_pointer(s->c);
	free_pointer(s->lb);
	free_pointer(s->ub);
	free_pointer(s->z0);
	free_pointer(s->zopt);
	free_pointer(s->z);
	free_pointer(s->znew);
	free_pointer(s->temp1_dim_N);
}
