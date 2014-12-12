/*
 * typedefs.h
 *
 *  Created on: Sep 26, 2014
 *      Author: sverre
 */

#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

/* MISRA C 2004 compliant numeric typedef */

typedef char char_t;

///* Probably already defined by inttypes.h or types.h*/
//# ifndef __int8_t_defined
//#  define __int8_t_defined
//typedef signed char int8_t;
//typedef signed short int16_t;
//typedef signed int int32_t;
//#endif
//
//#ifndef _SYS_TYPES_H
//#if !defined(_INT64_T) && !defined(INT64_MAX)
//#define _INT64_T
//typedef signed long int64_t;
//#endif
//#endif


typedef signed char int8_t;
typedef signed short int16_t;
typedef signed int int32_t;
//typedef signed long int64_t;

typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;
//typedef unsigned long uint64_t;

typedef float float32_t;
typedef double float64_t;
//typedef long double float128_t;

typedef uint32_t boolean;
typedef float64_t real_t;


#endif /* TYPEDEFS_H_ */
