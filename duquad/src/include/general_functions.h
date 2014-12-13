/*
 * general_functions.h
 *
 *  Created on: Oct 4, 2014
 *      Author: sverre
 */

/** \file
 * Contains some general functions that are common for most other files.
 */

#ifndef GENERAL_FUNCTIONS_H_
#define GENERAL_FUNCTIONS_H_

#include "head.h"

real_t * vector_alloc();
void free_pointer(real_t * pointer);

void initArray();
void insertArray();
void freeArray();

#endif /* GENERAL_FUNCTIONS_H_ */
