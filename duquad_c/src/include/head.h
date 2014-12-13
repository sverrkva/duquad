/*
 * head.h
 *
 *  Created on: Aug 19, 2014
 *      Author: Sverre
 */

/** \file
 * Contains system libraries, mex libraries, global constants, global variables and some debugging macros.
 */
#ifndef HEAD_H_
#define HEAD_H_

// *** System header files ***
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>

#include "typedefs.h"
#include "qp_structs.h"

#include <matrix.h>
#include <mex.h>
#include <string.h>

// *** Global constants ***
#define TRUE 1
#define FALSE 0

// Global variables
uint32_t N;	/**< Dimension of the Hessian matrix */
uint32_t M;	/**< Dimension of the linear constraint matrix */
uint32_t ALGORITHM;

// *** Macros ***
#ifdef DEBUG
	#define _DEBUG(fmt, args...) printf("%s:%s:%d: "fmt, __FILE__, __FUNCTION__, __LINE__, args)
#else
	#define _DEBUG(fmt, args...)
#endif
#ifdef DEBUG2
	#define _DEBUG2(fmt, args...) printf("%s:%s:%d: "fmt, __FILE__, __FUNCTION__, __LINE__, args)
#else
	#define _DEBUG2(fmt, args...)
#endif
/* Use: _DEBUG("hello world", arg)
add <-D -DEBUG> in gcc comand to compiler
e.g.: mex -O CFLAGS="\$CFLAGS -D DEBUG -std=c99" main.c calculations.c
Note: will not work if no 'arg' is used */

#define YO printf("YOYO\n")
#define ERROR(fmt) fprintf(stderr,"%s:%s:%d: "fmt, __FILE__, __FUNCTION__, __LINE__);

// A Simpler debug (without arguments)
#ifdef SDEBUG
	#define _SDEBUG(fmt) printf("%s:%s:%d: "fmt, __FILE__, __FUNCTION__, __LINE__)
#else
	#define _SDEBUG(fmt)
#endif



#endif /* HEAD_H_ */
