/*
 * general_functions.c
 *
 *  Created on: Oct 4, 2014
 *      Author: sverre
 */

#include "general_functions.h"

real_t * vector_alloc(uint32_t size)
{
	return (real_t*)calloc(size, sizeof(real_t));
}


void free_pointer(real_t * pointer)
{
	free((void*)pointer);
	pointer = NULL;

}

void initArray(struct Array *a, uint32_t initialSize) {
  a->array = (real_t *)malloc(initialSize * sizeof(real_t));
  a->used = 0;
  a->size = initialSize;
}

void insertArray(struct Array *a, real_t element) {
  if (a->used == a->size) {
    a->size *= 2;
    a->array = (real_t *)realloc(a->array, a->size * sizeof(real_t));
  }
  a->array[a->used++] = element;
}

void freeArray(struct Array *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}
