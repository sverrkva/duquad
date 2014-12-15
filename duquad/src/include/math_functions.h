/*
 * math_functions.h
 *
 *  Created on: Aug 21, 2014
 *      Author: Sverre
 */

/** \file
 * Header for the internal math library. Contains basic vector and matrix operations. The library is tried to be as optimized as possible.
 */

#ifndef MATH_FUNCTIONS_H_
#define MATH_FUNCTIONS_H_

#include "head.h"


/** matrix-vector multiplication: res = mtx * v */
void mtx_vec_mul(const real_t *mtx, const real_t *v, real_t *res, const uint32_t rows, const uint32_t cols);

/** Computes the transpose of the matrix: mtx_t = transpose(mtx) */
void mtx_transpose(const real_t * mtx, real_t * mtx_t, const uint32_t rows, const uint32_t cols);

/** compare two vectors element by element, and returns the smallest element of the two: res = min(v1,v2) */
void vector_min(const real_t *v1, const real_t *v2, real_t *res, const uint32_t length);

/** compare two vectors element by element, and returns the largest element of the two: res = max(v1,v2) */
void vector_max(const real_t *v1, const real_t *v2, real_t *res, const uint32_t length);

/** vector subtract: res = v1 - v2 */
void vector_sub(const real_t *v1, const real_t *v2, real_t *res, const uint32_t length);

/** vector addition res = v1 + v2 */
void vector_add(const real_t *v1, const real_t *v2, real_t *res, const uint32_t length);

/** vector multiplication: res = v1 * v2 */
real_t vector_mul(const real_t *v1,const real_t *v2, const uint32_t length);

/** Vector times a scalar: res = v1 * scalar */
void vector_scalar_mul(const real_t *v1, const real_t scalar, real_t *res, const uint32_t length);

/** set all elements in vector v to zero */
void vector_elements_to_zero(real_t *v, const uint32_t length);// may not be optimal

/** returns true if all elements in the two vectors are equal */
uint32_t vector_is_equal(const real_t *v1, const real_t *v2, const uint32_t length);

/** copy a vector element by element: v2 = v1 */
void vector_copy(const real_t *v1, real_t *v2, const uint32_t length);

/** project a vector on all positive numbers: v = max(v,0.0) */
void vector_max_with_zero(real_t *v, const uint32_t length);

/** compute the euclidean norm of the vector v */
real_t vector_norm_2(real_t *v, const uint32_t length);

/** return the absolute value of the input */
real_t abs_2(const real_t a);

/** compute the primal objective function */
real_t obj(const real_t *z, const real_t *H, const real_t *c, real_t *temp);

#endif /* MATH_FUNCTIONS_H_ */
