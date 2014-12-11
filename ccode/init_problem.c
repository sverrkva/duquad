/*
 * init_problem.c
 *
 *  Created on: Aug 19, 2014
 *      Author: Sverre
 */


#include "init_problem.h"

static void allocate_Problem(struct Problem *s);
static void clean_up_Problem(struct Problem *s);
static void clean_up_Output(struct Output *s);
static void allocate_DGM(struct Struct_DGM *s);
static void allocate_DFGM(struct Struct_DFGM *s);
static void vector_read_from_file(real_t *v, FILE * file, const uint32_t length);
static void matrix_read_from_file(real_t *mtx, FILE * file, const uint32_t rows, const uint32_t cols);

void init_problem_from_file_DGM(struct Struct_DGM *s)
{

	FILE *file;
	if((file = fopen("../../matlab/problem.txt", "r")) == NULL){
		printf("could not open file\n");
		return;
	}

	uint32_t itemp;
	real_t dtemp;
	fscanf(file,"%d", &itemp); N = itemp;
	fscanf(file,"%d", &itemp); M = itemp;

	static struct Info info_global;
	fscanf(file,"%d", &itemp); info_global.lb_is_inf = itemp;
	fscanf(file,"%d", &itemp); info_global.ub_is_inf = itemp;
	fscanf(file,"%d", &itemp); info_global.problem_case = itemp;
	fscanf(file,"%lf", &dtemp); info_global.eigH_min = dtemp;
	fscanf(file,"%lf", &dtemp); info_global.eigH_max = dtemp;
	fscanf(file,"%lf", &dtemp); info_global.Ld = dtemp;

    switch (info_global.problem_case) {
    case 1:
    	info_global.lb_hat_is_inf = FALSE;
    	info_global.ub_hat_is_inf = FALSE;
    	info_global.pf_vec_length = 2*M;
        break;
    case 2:
    	info_global.lb_hat_is_inf = FALSE;
    	info_global.ub_hat_is_inf = FALSE;
    	info_global.pf_vec_length = M;
        break;
    case 3:
    	info_global.lb_hat_is_inf = TRUE;
    	info_global.ub_hat_is_inf = FALSE;
    	info_global.pf_vec_length = M;
        break;
    case 4:
    	info_global.lb_hat_is_inf = FALSE;
    	info_global.ub_hat_is_inf = TRUE;
    	info_global.pf_vec_length = M;
        break;
    default:
    	printf("Could not resolve problem case\n");
        break;
    }

    s->info = &info_global;

    static struct Options opt_global;
	fscanf(file,"%d", &itemp); opt_global.maxiter_outer = itemp;
	fscanf(file,"%d", &itemp); opt_global.maxiter_inner = itemp;
	fscanf(file,"%lf", &dtemp); opt_global.eps_ds = dtemp;
	fscanf(file,"%lf", &dtemp); opt_global.eps_pf = dtemp;
	fscanf(file,"%lf", &dtemp); opt_global.eps_inner = dtemp;
	fscanf(file,"%d", &itemp); opt_global.algorithm = itemp;

	s->opt = &opt_global;

	static struct Problem prob_global;
	allocate_Problem(&prob_global);

	matrix_read_from_file(prob_global.H,file,N,N);
	vector_read_from_file(prob_global.c,file, N);
	matrix_read_from_file(prob_global.A,file,M,N);
	vector_read_from_file(prob_global.b,file, M);
	if (s->info->problem_case == 1 || s->info->problem_case == 4)
		vector_read_from_file(prob_global.lb_hat,file, M);
	if (s->info->problem_case !=4)
		vector_read_from_file(prob_global.ub_hat,file, M);
	vector_read_from_file(prob_global.lb,file, N);
	vector_read_from_file(prob_global.ub,file, N);
	vector_read_from_file(prob_global.z0,file, N);

	mtx_transpose(prob_global.A,prob_global.A_t,M,N);

	s->prob = & prob_global;

	s->res = malloc(sizeof * s->res);
	s->res->fopt = 0;
	s->res->out = malloc(sizeof * s->res->out);

	allocate_DGM(s);

	fclose(file);
}

// Almost identical to init_problem_from_file_DGM
void init_problem_from_file_DFGM(struct Struct_DFGM *s)
{

	FILE *file;
	if((file = fopen("../../matlab/problem.txt", "r")) == NULL){
		printf("could not open file\n");
		return;
	}

	uint32_t itemp;
	real_t dtemp;
	fscanf(file,"%d", &itemp); N = itemp;
	fscanf(file,"%d", &itemp); M = itemp;

	static struct Info info_global;
	fscanf(file,"%d", &itemp); info_global.lb_is_inf = itemp;
	fscanf(file,"%d", &itemp); info_global.ub_is_inf = itemp;
	fscanf(file,"%d", &itemp); info_global.problem_case = itemp;
	fscanf(file,"%lf", &dtemp); info_global.eigH_min = dtemp;
	fscanf(file,"%lf", &dtemp); info_global.eigH_max = dtemp;
	fscanf(file,"%lf", &dtemp); info_global.Ld = dtemp;

	switch (info_global.problem_case) {
	case 1:
		info_global.lb_hat_is_inf = FALSE;
		info_global.ub_hat_is_inf = FALSE;
		info_global.pf_vec_length = 2*M;
		break;
	case 2:
		info_global.lb_hat_is_inf = FALSE;
		info_global.ub_hat_is_inf = FALSE;
		info_global.pf_vec_length = M;
		break;
	case 3:
		info_global.lb_hat_is_inf = TRUE;
		info_global.ub_hat_is_inf = FALSE;
		info_global.pf_vec_length = M;
		break;
	case 4:
		info_global.lb_hat_is_inf = FALSE;
		info_global.ub_hat_is_inf = TRUE;
		info_global.pf_vec_length = M;
		break;
	default:
		printf("Could not resolve problem case\n");
		break;
	}

	s->info = &info_global;

	static struct Options opt_global;
	fscanf(file,"%d", &itemp); opt_global.maxiter_outer = itemp;
	fscanf(file,"%d", &itemp); opt_global.maxiter_inner = itemp;
	fscanf(file,"%lf", &dtemp); opt_global.eps_ds = dtemp;
	fscanf(file,"%lf", &dtemp); opt_global.eps_pf = dtemp;
	fscanf(file,"%lf", &dtemp); opt_global.eps_inner = dtemp;
	fscanf(file,"%d", &itemp); opt_global.algorithm = itemp;

	s->opt = &opt_global;

	static struct Problem prob_global;
	allocate_Problem(&prob_global);

	matrix_read_from_file(prob_global.H,file,N,N);
	vector_read_from_file(prob_global.c,file, N);
	matrix_read_from_file(prob_global.A,file,M,N);
	vector_read_from_file(prob_global.b,file, M);
	if (s->info->problem_case == 1 || s->info->problem_case == 4)
		vector_read_from_file(prob_global.lb_hat,file, M);
	if (s->info->problem_case !=4)
		vector_read_from_file(prob_global.ub_hat,file, M);
	vector_read_from_file(prob_global.lb,file, N);
	vector_read_from_file(prob_global.ub,file, N);
	vector_read_from_file(prob_global.z0,file, N);

	mtx_transpose(prob_global.A,prob_global.A_t,M,N);

	s->prob = &prob_global;

	s->res = malloc(sizeof * s->res);
	s->res->fopt = 0;
	s->res->out = malloc(sizeof * s->res->out);

	allocate_DFGM(s);

	fclose(file);
}

void clean_up_DGM_C(struct Struct_DGM *s)
{
	clean_up_Problem(s->prob);
//	free(s->opt); s->opt = NULL;
//	free(s->info); s->info = NULL;
	free_pointer(s->z);
	free_pointer(s->lambda1);
	free_pointer(s->lambda2);
	free_pointer(s->temp1_dim_N);
	free_pointer(s->temp2_dim_M);
	free_pointer(s->temp3_dim_M);
	free_pointer(s->b_ub_hat);
	free_pointer(s->b_lb_hat);
	free_pointer(s->z_avg);
	free_pointer(s->summ);
	free_pointer(s->pf_vec);
	free_pointer(s->A_z);
	clean_up_Output(s->res->out);

	free(s->res->out);
	s->res->out = NULL;
	free(s->res);
	s->res = NULL;
}

void clean_up_DFGM_C(struct Struct_DFGM *s)
{
	clean_up_Problem(s->prob);
	free_pointer(s->z);
	free_pointer(s->lambda1);
	free_pointer(s->lambda2);
	free_pointer(s->temp1_dim_N);
	free_pointer(s->temp2_dim_M);
	free_pointer(s->temp3_dim_M);
	free_pointer(s->b_ub_hat);
	free_pointer(s->b_lb_hat);
	free_pointer(s->z_avg);
	free_pointer(s->summ);
	free_pointer(s->pf_vec);
	free_pointer(s->A_z);
	clean_up_Output(s->res->out);
	free(s->res->out);
	s->res->out = NULL;
	free(s->res);
	s->res = NULL;

	// Different from DGM()
	free_pointer(s->lambda1_old);
	free_pointer(s->lambda2_old);
	free_pointer(s->y1);
	free_pointer(s->y2);
	free_pointer(s->z_ds);
	free_pointer(s->A_z_ds);
}



static void vector_read_from_file(real_t *v, FILE * file, const uint32_t length)
{
	uint32_t i;
	for (i=0;i<length;i++){
//		fscanf(file,"%lf,",&v[i]);
		fscanf(file,"%lf,",&v[i]);
	}
}

static void matrix_read_from_file(real_t *mtx, FILE * file, const uint32_t rows, const uint32_t cols)
{
	uint32_t i;	/* row number */
	uint32_t j; /* column number */
	uint32_t k = 0; /* matrix index (row * col) */
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			fscanf(file,"%lf,",&mtx[k]);
			k++;
		}
	}
}

static void allocate_DGM(struct Struct_DGM *s)
{
	s->z = vector_alloc(N);
	s->lambda1 = vector_alloc(M);
	s->lambda2 = vector_alloc(M);
	s->temp1_dim_N = vector_alloc(N);
	s->temp2_dim_M = vector_alloc(M);
	s->temp3_dim_M = vector_alloc(M);
	s->b_ub_hat = vector_alloc(M);
	s->b_lb_hat = vector_alloc(M);
	s->z_avg = vector_alloc(N);
	s->summ = vector_alloc(N);
	s->pf_vec = vector_alloc(s->info->pf_vec_length);
	s->A_z = vector_alloc(M);
	s->res->out->ds_vector  = vector_alloc(2);
	s->res->out->pf_vector  = vector_alloc(2);
}

static void allocate_DFGM(struct Struct_DFGM *s)
{
	s->z = vector_alloc(N);
	s->lambda1 = vector_alloc(M);
	s->lambda2 = vector_alloc(M);
	s->temp1_dim_N = vector_alloc(N);
	s->temp2_dim_M = vector_alloc(M);
	s->temp3_dim_M = vector_alloc(M);
	s->b_ub_hat = vector_alloc(M);
	s->b_lb_hat = vector_alloc(M);
	s->z_avg = vector_alloc(N);
	s->summ = vector_alloc(N);
	s->pf_vec = vector_alloc(s->info->pf_vec_length);
	s->A_z = vector_alloc(M);
	s->res->out->ds_vector  = vector_alloc(2);
	s->res->out->pf_vector  = vector_alloc(2);

	// Different form DGM()
	s->lambda1_old = vector_alloc(M);
	s->lambda2_old = vector_alloc(M);
	s->y1 = vector_alloc(M);
	s->y2 = vector_alloc(M);
	s->z_ds = vector_alloc(N);
	s->A_z_ds = vector_alloc(M);
}

static void allocate_Problem(struct Problem *s)
{
	s->H = vector_alloc(N*N);
	s->c = vector_alloc(N);
	s->A = vector_alloc(M*N);
	s->A_t = vector_alloc(N*M);
	s->b = vector_alloc(M);
	s->lb_hat = vector_alloc(M);
	s->ub_hat = vector_alloc(M);
	s->lb = vector_alloc(N);
	s->ub = vector_alloc(N);
	s->z0 = vector_alloc(N);
}

static void clean_up_Problem(struct Problem *s)
{
	free_pointer(s->H);
	free_pointer(s->c);
	free_pointer(s->A);
	free_pointer(s->A_t);
	free_pointer(s->b);
	free_pointer(s->lb_hat);
	free_pointer(s->ub_hat);
	free_pointer(s->lb);
	free_pointer(s->ub);
	free_pointer(s->z0);
	s = NULL;
}

static void clean_up_Output(struct Output *s)
{
	free_pointer(s->ds_vector);
	free_pointer(s->pf_vector);
	s = NULL;
}



// Fast Gradient

void init_problem_FGM(struct Struct_FGM *s)
{
	N = 3;

	s->H = 				vector_alloc(N*N);
	s->c = 				vector_alloc(N);
	s->lb = 			vector_alloc(N);
	s->ub = 			vector_alloc(N);
	s->z0 = 			vector_alloc(N);

	s->zopt = 			vector_alloc(N);
	s->fopt = 			0;

	s->lb_is_inf = 		TRUE;
	s->ub_is_inf = 		TRUE;
	s->eigH_min = 		0.;
	s->eigH_max = 		0.;

	s->z = 				vector_alloc(N);
	s->y = 				vector_alloc(N);
	s->znew = 			vector_alloc(N);
	s->ynew = 			vector_alloc(N);
	s->temp1_dim_N = 	vector_alloc(N);

	s->maxiter = 10000;
	s->eps = 0.00001;
}

void init_problem_GDM(struct Struct_GDM *s)
{
	N = 3;

	s->H = 				vector_alloc(N*N);
	s->c = 				vector_alloc(N);
	s->lb = 			vector_alloc(N);
	s->ub = 			vector_alloc(N);
	s->z0 = 			vector_alloc(N);

	s->zopt = 			vector_alloc(N);
	s->fopt = 			0;

	s->lb_is_inf = 		TRUE;
	s->ub_is_inf = 		TRUE;
	s->eigH_min = 		0.;
	s->eigH_max = 		0.;

	s->z = 				vector_alloc(N);
	s->znew = 			vector_alloc(N);
	s->temp1_dim_N = 	vector_alloc(N);

	s->maxiter = 10000;
	s->eps = 0.00001;
}

void create_problem_FGM_1(struct Struct_FGM *s)
{

	s->H[0] = 1.2204;
	s->H[1] = 1.1123;
	s->H[2] = -3.8935;
	s->H[3] = 1.1123;
	s->H[4] = 3.5821;
	s->H[5] = -3.3333;
	s->H[6] = -3.8935;
	s->H[7] = -3.3333;
	s->H[8] = 19.6174;

	s->c[0] = 4;
	s->c[1] = 3;
	s->c[2] = -2;

	s->lb[0] = -1;
	s->lb[1] = -1;
	s->lb[2] = -1;
	s->lb_is_inf = FALSE;

	s->ub[0] = 1;
	s->ub[1] = 1;
	s->ub[2] = 1;
	s->ub_is_inf = FALSE;

	s->z0[0] = 0;
	s->z0[1] = 0;
	s->z0[2] = 0;

	s->eigH_min = 0.3605;
	s->eigH_max = 21.1022;

}

void create_problem_GDM_1(struct Struct_GDM *s)
{

	s->H[0] = 1.2204;
	s->H[1] = 1.1123;
	s->H[2] = -3.8935;
	s->H[3] = 1.1123;
	s->H[4] = 3.5821;
	s->H[5] = -3.3333;
	s->H[6] = -3.8935;
	s->H[7] = -3.3333;
	s->H[8] = 19.6174;

	s->c[0] = 4;
	s->c[1] = 3;
	s->c[2] = -2;

	s->lb[0] = -1;
	s->lb[1] = -1;
	s->lb[2] = -1;
	s->lb_is_inf = FALSE;

	s->ub[0] = 1;
	s->ub[1] = 1;
	s->ub[2] = 1;
	s->ub_is_inf = FALSE;

	s->z0[0] = 0;
	s->z0[1] = 0;
	s->z0[2] = 0;

	s->eigH_min = 0.3605;
	s->eigH_max = 21.1022;

}
