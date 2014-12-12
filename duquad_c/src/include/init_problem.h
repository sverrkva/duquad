/*
 * init_problem.h
 *
 *  Created on: Aug 19, 2014
 *      Author: Sverre
 */

#ifndef INIT_PROBLEM_H_
#define INIT_PROBLEM_H_

#include "head.h"
#include "fgm.h"
#include "dgm.h"
#include "dfgm.h"
#include "general_functions.h"


void init_problem_from_file_DGM();
void init_problem_from_file_DFGM();

void clean_up_DGM_C();
void clean_up_DFGM_C();



//void init_optDGM();
//void init_opt_DFGM();




void init_problem_FGM();
void create_problem_FGM_1();
void init_problem_GDM();
void create_problem_GDM_1();


#endif /* INIT_PROBLEM_H_ */
