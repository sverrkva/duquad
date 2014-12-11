/*
 * math_functions.h
 *
 *  Created on: Aug 21, 2014
 *      Author: Sverre
 */

#ifndef MATH_FUNCTIONS_H_
#define MATH_FUNCTIONS_H_

#include "head.h"

void mtx_vec_mul(); // res = mtx * v
void vector_min();
void vector_max();
void vector_sub();		// res = v1-v2
void vector_add();   		// res = v1+v2
void vector_scalar_mul(); // res = v1 * scalar
void vector_elements_to_zero();
uint32_t vector_is_equal();
void vector_copy();	// copy v1 to v2
void vector_max_with_zero();
real_t vector_norm_2();

real_t abs_2();


void mtx_transpose();

real_t obj();
real_t vector_mul();

#endif /* MATH_FUNCTIONS_H_ */
