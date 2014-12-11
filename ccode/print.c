/*
 * print.c
 *
 *  Created on: Aug 21, 2014
 *      Author: Sverre
 */
#include "print.h"

void print_vector(double *v, int length)	// Maybe make private function
{
	int i;
	printf("[%.8f, ", v[0]);
	for(i=1;i<length-1;i++){
			printf("%.8f, ", v[i]);
	}
	printf("%.8f", v[length-1]);
	printf("]'\n");
}

void print_matrix(const real_t *mtx, const uint32_t rows, const uint32_t cols)	// Maybe make private function
{
	printf("\n");
	uint32_t i;	/* row number */
	uint32_t j; /* column number */
	uint32_t k = 0; /* matrix index (row * column) */
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			printf("%.4f, ", mtx[k]);
			k++;
		}
		printf("\n");
	}
}

void print_Problem(struct Problem * s)
{
	printf("\n");
	printf("\n-- Problem --\n");
	printf("H: "); print_matrix(s->H,N,N);
	printf("\n");
	printf("c: "); print_vector(s->c,N);
	printf("A: "); print_matrix(s->A,M,N);
	printf("\n");
	printf("b: ");	print_vector(s->b,M);
	printf("lb_hat: "); print_vector(s->lb_hat,M);
	printf("ub_hat: ");	print_vector(s->ub_hat,M);
	printf("lb: "); print_vector(s->lb,N);
	printf("ub: ");	print_vector(s->ub,N);
	printf("z0: ");	print_vector(s->z0,N);
}

void print_Options(struct Options * opt)
{
	printf("\n");
	printf("\n-- Options --\n");
	printf("maxiter_outer: %d\n",opt->maxiter_outer);
	printf("maxiter_inner: %d\n",opt->maxiter_inner);
	printf("eps_ds: %f\n",opt->eps_ds);
	printf("eps_pf: %f\n",opt->eps_pf);
	printf("eps_inner: %f\n",opt->eps_inner);
	printf("algorithm: %d\n",opt->algorithm);
}

void print_Info(struct Info * info)
{
	printf("\n");
	printf("\n-- Info --\n");
	printf("lb_is_inf: %d\n",info->lb_is_inf);
	printf("ub_is_inf: %d\n",info->ub_is_inf);
	printf("lb_hat_is_inf: %d\n",info->lb_hat_is_inf);
	printf("ub_hat_is_inf: %d\n",info->ub_hat_is_inf);
	printf("eigH_max: %f\n",info->eigH_max);
	printf("eigH_min: %f\n",info->eigH_min);
	printf("Ld: %f\n",info->Ld);
	printf("problem_case: %d\n",info->problem_case);
	printf("pf_vec_length: %d\n",info->pf_vec_length);
}

void print_Result(struct Result *s)
{
	printf("\n");
	printf("\n-- RESULT --\n");
	//printf("zopt: "); print_vector(s->zopt,N);
	//printf("lambda1: "); print_vector(s->lambda1,N);
	//printf("lambda2: "); print_vector(s->lambda2,N);
	printf("fopt: %f\n", s->fopt);
	printf("iterations: %d\n", s->out->iterations);
	printf("totale iterations inner: %d\n", s->out->iterations_inner_tot);
	printf("niter_feasible_ds: %d\n", s->out->niter_feasible_ds);
	printf("niter_feasible_pf: %d\n", s->out->niter_feasible_pf);
	printf("flag_last_satisfied: %d\n", s->out->flag_last_satisfied);
	//printf("ds vector: "); print_vector(s->out->ds_vector,s->out->iterations);
	//printf("pf vector: "); print_vector(s->out->pf_vector,s->out->iterations);
	printf("time: %f\n", s->out->time);
	printf("time_tot_inner: %f\n", s->out->time_tot_inner);
	printf("exitflag: %d\n", s->exitflag);
	printf("exitflag_inner: %d\n", s->out->exitflag_inner);
	printf("num reached max iter inner: %d\n", s->out->num_exceeded_max_niter_inner);
	printf("---\n");
}

void print_problem_FGM(struct Struct_FGM *s)
{
	printf("\n");
	printf("\n-- Problem --\n");
	printf("H: "); print_matrix(s->H,N,N);
	printf("c: "); print_vector(s->c,N);
	printf("lb: "); print_vector(s->lb,N);
	printf("ub: ");	print_vector(s->ub,N);
	printf("z0: ");	print_vector(s->z0,N);
	printf("\n");
	printf("lb_is_inf: %d\n", s->lb_is_inf );
	printf("ub_is_inf: %d\n", s->ub_is_inf );
	printf("eigH_min: %f\n", s->eigH_min );
	printf("eigH_max: %f\n", s->eigH_max );
	printf("\n");
}

void print_result_FGM(struct Struct_FGM *s, int niter)
{
	printf("\n");
	printf("\n-- RESULT --\n");
	printf("zopt: "); print_vector(s->zopt,N);
	printf("fopt: %f\n", s->fopt);
	printf("niter: %d\n", niter);
	printf("---\n");
}

void print_result_GDM(struct Struct_GDM *s, int niter)
{
	printf("\n");
	printf("\n-- RESULT --\n");
	printf("zopt: "); print_vector(s->zopt,N);
	printf("fopt: %f\n", s->fopt);
	printf("niter: %d\n", niter);
	printf("---\n");
}
