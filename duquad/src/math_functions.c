/*
 * functions.c
 *
 *  Created on: Aug 21, 2014
 *      Author: Sverre
 */

#include "math_functions.h"

static void mtx_vec_mul_small(const real_t *mtx, const real_t *v, real_t *res, const uint32_t rows, const uint32_t cols)
{
	uint32_t i;	/* row number */
	uint32_t j; /* column number */
	uint32_t k = 0; /* matrix index (row * col) */
	for (i=0;i<rows;i++) {
		res[i] = 0.0;
		for (j=0;j<cols;j++) {
			res[i] += mtx[k++] * v[j];
		}
	}
}

void mtx_vec_mul(const real_t *mtx, const real_t *v, real_t *res, const uint32_t rows, const uint32_t cols)
{
	if(rows < 10 || cols < 10){
		mtx_vec_mul_small(mtx,v,res,rows,cols);
	}
	else{
		uint32_t i;	/* row number */
		uint32_t j; /* column number */
		uint32_t k;
		real_t temp;
		k = 0;
		for (i=0;i<rows;i++) {
			temp = 0.0;
			for (j=0;j<=cols-4;j+=4) {
				temp += (mtx[k] * v[j] +
						 mtx[k+1] * v[j+1] +
						 mtx[k+2] * v[j+2] +
						 mtx[k+3] * v[j+3] );
				k+=4;
			}
			for (; j<cols; j++){
				temp += mtx[k++] * v[j];
			}
			res[i] = temp;
		}
	}
}

void mtx_transpose(const real_t * mtx, real_t * mtx_t, const uint32_t rows, const uint32_t cols)
{
	uint32_t i;	/* row number */
	uint32_t j; /* column number */
	uint32_t k = 0; /* matrix index (row * column) */
	uint32_t kT; /* matrix transpose index */
	for (i=0;i<rows;i++) {
		for (j=0;j<cols;j++) {
			kT = j * rows + i;
			mtx_t[kT] = mtx[k];
			k++;
		}
	}
}

void vector_min(const real_t *v1, const real_t *v2, real_t *res, const uint32_t length)
{
	uint32_t i;
	for(i=0;i<length;i++){
		if(v1[i] < v2[i])
			res[i] = v1[i];
		else
			res[i] = v2[i];
	}
}

void vector_max(const real_t *v1, const real_t *v2, real_t *res, const uint32_t length)
{
	uint32_t i;
	for(i=0;i<length;i++){
		if(v1[i] > v2[i])
			res[i] = v1[i];
		else
			res[i] = v2[i];
	}
}

void vector_sub(const real_t *v1, const real_t *v2, real_t *res, const uint32_t length)
{
	uint32_t i;
	for(i=0;i<length;i++){
		res[i] = v1[i] - v2[i];
	}
}

void vector_add(const real_t *v1, const real_t *v2, real_t *res, const uint32_t length)
{
	uint32_t i;
	for(i=0;i<length;i++){
		res[i] = v1[i] + v2[i];
	}
}

real_t vector_mul(const real_t *v1,const real_t *v2, const uint32_t length)
{
	uint32_t i;
	real_t val = 0.0;
	for(i=0;i<length;i++){
		val += v1[i] * v2[i];
	}
	return val;
}

void vector_scalar_mul(const real_t *v1, const real_t scalar, real_t *res, const uint32_t length)
{
	uint32_t i;
	for(i=0;i<length;i++){
		res[i] = v1[i] * scalar;
	}
}

void vector_elements_to_zero(real_t *v, const uint32_t length)// may not be optimal
{
	uint32_t i;
	for(i=0;i<length;i++){
		v[i] = 0.0;
	}
}

uint32_t vector_is_equal(const real_t *v1, const real_t *v2, const uint32_t length)
{
	uint32_t i;
	for(i=0;i<length;i++){
		if (v1[i] != v2[i])
			return FALSE;
	}
	return TRUE;
}

void vector_copy(const real_t *v1, real_t *v2, const uint32_t length)
{
	uint32_t i;
	for(i=0;i<length;i++){
		v2[i] = v1[i];
	}
}

void vector_max_with_zero(real_t *v, const uint32_t length)
{
	uint32_t i;
	for(i=0;i<length;i++){
		if(v[i] < 0.0)
			v[i] = 0.0;
	}
}

real_t vector_norm_2(real_t *v, const uint32_t length)
{

	uint32_t i;
	real_t m_sum = 0.0;
	for (i=0;i<length;i++)
	{
	    m_sum  += v[i] * v[i];
	}
	return sqrt(m_sum);	// might overflow
//	A safer algorithm would be:
//	 y=max(abs(x))
//	for (uint32_t i = 1; i<=n;i++){
//	      m_sum+=(x(i)/y)^2;
//	}
//	 return y*sqrt(m_sum);

}

real_t abs_2(const real_t a)
{
	if(a < 0.0)
		return a * -1.0;
	return a;
}

real_t obj(const real_t *z, const real_t *H, const real_t *c, real_t *temp)
{
	// Calculate: f_value = z'*H*z + c'z
	real_t f_value = 0.0;
	mtx_vec_mul(H,z,temp,N,N);	// z' * H
	f_value += vector_mul(temp,z,N);// * z
	f_value *= 0.5; 				// * 0.5
	f_value += vector_mul(c,z,N); 	// + c'z
	return f_value;
}



