/*
 * qp_structs.h
 *
 *  Created on: Oct 1, 2014
 *      Author: sverre
 */

/** \file
 * Contains the declaration of the structs that are common for the algorithms
 */
#ifndef QP_STRUCTS_H_
#define QP_STRUCTS_H_
#include "typedefs.h"

struct Problem {
	real_t * H;	/**< The Hessian matrix. Dimensions n x n, and has to be positive definite */
	real_t * c; /**< The gradient vector */
	real_t * A; /**< Linear constraints matrix Dimensions m x n */
	real_t * A_t; /**< A transposed */
	real_t * b; /**< Linear constraints vector */
	real_t * lb_hat; /**< The lower bound for the linear constraints */
	real_t * ub_hat; /**< The upper bound for the linear constraints */
	real_t * lb; /**< The lower bound for optimization variable z */
	real_t * ub; /**< The upper bound for optimization variable z */
	real_t * z0; /**< The initial point */
}; /**< Contains all the matrices and vectors used to describe the general QP */

struct Options {
	uint32_t maxiter_outer; /**< Maximum number of iterations in the outer loop */
	uint32_t maxiter_inner; /**< Maximum number of iterations in the inner loop */
	real_t eps_ds; /**< Tolerance for dual suboptimality */
	real_t eps_pf; /**< Tolerance for primal feasibility */
	real_t eps_inner; /**< Tolerance for primal feasibility in the inner problem */
	uint32_t algorithm; /**< Specifies the algorithm used to solve the problem. Values: 1: DGM last, 2: DGM avg, 3: DFGM last, 4: DFGM avg, 5: ALM last, 6: ALM avg, 7: FALM last, 8: FALM avg */
	real_t rho; /**< Penalty parameter used in ALM and FALM */
}; /**< Option specified by the user. Default values can be found in the user manual */

struct Info {
	boolean lb_is_inf; /**< true: no lower bound on optimization variable */
	boolean ub_is_inf; /**< true: no upper bound on optimization variable */
	boolean lb_hat_is_inf; /**< true: no lower bound on linear constraints */
	boolean ub_hat_is_inf; /**< true: no upper bound on linear constraint */
	real_t eigH_max; /**< Largest eigenvalue of the Hessian H */
	real_t eigH_min; /**< Smallest eigenvalue of the Hessian H */
	real_t Ld; /**< Lipschitz constant */
	uint32_t problem_case; /**< case 1: lb_hat != ub_hat, case 2: lb_hat == ub_hat, case 3: lb_hat = -inf, case 4: ub_hat = inf */
	uint32_t pf_vec_length;
}; /**< parameters that are calculated automatically off-line in duquad.m  */

struct Output {
	uint32_t iterations; /**< Number of outer iterations */
	uint32_t iterations_inner_tot; /**< Total number of iterations for the inner problem  */
	real_t time; /**< Runtime of the algorithm after all initialization is done */
	real_t time_tot_inner; /**< Total time spent on solving the inner problem */
	uint32_t flag_last_satisfied; /**< Flag spesifies which stopping criteria was resolved last. Value: 0 = dual suboptimality, 1 = primal feasibility  */
	uint32_t niter_feasible_ds; /**< Number of iterations the criterion for dual suboptimality was satisfied */
	uint32_t niter_feasible_pf; /**< Number of iterations the criterion for primal feasibility was satisfied */
	uint32_t exitflag_inner; /**< Exitflag for the inner problem. Values: 1 = feasible point found, 2 = Maximum number of iterations exceeded */
	uint32_t num_exceeded_max_niter_inner; /**< Total number of times the inner problem exceeded the number of iterations */
	real_t * ds_vector; /**< Vector storing all the value of the dual suboptimality every iteration */
	real_t * pf_vector; /**< Vector storing all the value of the primal feasibility every iteration */
}; /**< Important results are collected in the Output struct */

struct Result {
	real_t * zopt; /**< Optimal point */
	real_t fopt; /**< Optimal value */
	uint32_t exitflag; /**< Values: 1 = optimal point found, 2 = maximum number of iterations exceeded, -1 = error */
	real_t * lambda1; /**< Set of lagrangian multipliers */
	real_t * lambda2; /**< Set of lagrangian multipliers */
	struct Output * out; /**< Sruct containing other results */
}; /**< Outputs of the algorithms */

struct Array{
  real_t * array;
  uint32_t used;
  uint32_t size;
}; /**< for internal use */


#endif /* QP_STRUCTS_H_ */
