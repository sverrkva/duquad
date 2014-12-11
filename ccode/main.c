

#include "head.h"
#include "gdm.h"
#include "fgm.h"
#include "dgm.h"
#include "print.h"
#include "dfgm.h"
#include "alm.h"
#include "falm.h"



// Inputs
#define H_IN prhs[0]
#define C_IN prhs[1]
#define A_IN prhs[2]
#define B_IN prhs[3]
#define LB_HAT_IN prhs[4]
#define UB_HAT_IN prhs[5]
#define LB_IN prhs[6]
#define UB_IN prhs[7]
#define Z0_IN prhs[8]
#define OPT_IN prhs[9]
#define INFO_IN prhs[10]

// Used in augmented cases
#define H_HAT_IN prhs[11]
#define A2_IN prhs[12]
#define RHO_AT_B_IN prhs[13]

// Outputs
#define ZOPT_OUT plhs[0]
#define FOPT_OUT plhs[1]
#define EXITFLAG_OUT plhs[2]
#define OUTPUT_STRUCT_OUT plhs[3]
#define LAMBDA1_OUT plhs[4]
#define LAMBDA2_OUT plhs[5]

// Private functions
static void init_Info();
static void init_Options();
static void init_Problem();

static void make_output();
static void make_output_fail_safe();
static void make_Output_struct();

static void clean_up_Output_M();    // Same function in init_problem.h
static void clean_up_Result();

static void init_DGM();
static void allocate_DGM();
static void clean_up_DGM_matlab();

static void init_DFGM();
static void allocate_DFGM();
static void clean_up_DFGM_matlab();

static void init_ALM();
static void allocate_ALM();
static void clean_up_ALM_matlab();

static void init_FALM();
static void allocate_FALM();
static void clean_up_FALM_matlab();


// z0 might be a problem
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    _SDEBUG("*** Main start**\n");
    

    int algorithm = mxGetScalar(mxGetField(OPT_IN, 0, "algorithm"));
    ALGORITHM = algorithm;
    
    if (algorithm == 1 || algorithm == 2)
    {
        struct Struct_DGM s;
        init_DGM(nrhs, prhs, &s);
       
        if (DGM(&s) == -1){
            ERROR("Error in DGM()\n");
            make_output_fail_safe(nlhs, plhs);
            clean_up_DGM_matlab(&s);
            return;
        }
 //       print_Options(s.opt);
//         print_Info(s.info);
//         print_Result(s.res);
//         print_Problem(s.prob);
        make_output(nlhs, plhs, s.res);
        clean_up_DGM_matlab(&s);
    }
    else if (algorithm == 3 || algorithm == 4)
    {
        struct Struct_DFGM s;
        init_DFGM(nrhs, prhs, &s);
        if (DFGM(&s) == -1){
            ERROR("Error in DFGM()\n");
            make_output_fail_safe(nlhs, plhs);
            clean_up_DFGM_matlab(&s);
            return;
        }
//         print_Options(s.opt);
//         print_Info(s.info);
//         print_Result(s.res);
//         print_Problem(s.prob);
        make_output(nlhs, plhs, s.res);
        clean_up_DFGM_matlab(&s);
    }
    else if (algorithm == 5 || algorithm == 6)
    {
        struct Struct_ALM s;
        init_ALM(nrhs, prhs, &s);
        if (ALM(&s) == -1){
            ERROR("Error in ALM()\n");
            make_output_fail_safe(nlhs, plhs);
            clean_up_ALM_matlab(&s);
            return;
        }
        make_output(nlhs, plhs, s.res);
        clean_up_ALM_matlab(&s);
    }
    
    else if (algorithm == 7 || algorithm == 8)
    {
        struct Struct_FALM s;
        init_FALM(nrhs, prhs, &s);
        if (FALM(&s) == -1){
            ERROR("Error in FALM()\n");
            make_output_fail_safe(nlhs, plhs);
            clean_up_FALM_matlab(&s);
            return;
        }
        make_output(nlhs, plhs, s.res);
        clean_up_FALM_matlab(&s);
    }
    

    _SDEBUG("*** Main finish**\n");
    return;
}

static void init_Info(int nrhs, const mxArray *prhs[], struct Info *p)
{
    p->lb_is_inf = mxGetScalar(mxGetField(INFO_IN, 0, "lb_is_inf"));
    p->ub_is_inf = mxGetScalar(mxGetField(INFO_IN, 0, "ub_is_inf"));
    p->problem_case = mxGetScalar(mxGetField(INFO_IN, 0, "problem_case"));
    p->eigH_min = mxGetScalar(mxGetField(INFO_IN, 0, "eigH_min"));
    p->eigH_max = mxGetScalar(mxGetField(INFO_IN, 0, "eigH_max"));
    p->Ld = mxGetScalar(mxGetField(INFO_IN, 0, "Ld"));

    switch (p->problem_case) {
    case 1:
        p->lb_hat_is_inf = FALSE;
        p->ub_hat_is_inf = FALSE;
        p->pf_vec_length = 2*M;
        break;
    case 2:
        p->lb_hat_is_inf = FALSE;
        p->ub_hat_is_inf = FALSE;
        p->pf_vec_length = M;
        break;
    case 3:
        p->lb_hat_is_inf = TRUE;
        p->ub_hat_is_inf = FALSE;
        p->pf_vec_length = M;
        break;
    case 4:
        p->lb_hat_is_inf = FALSE;
        p->ub_hat_is_inf = TRUE;
        p->pf_vec_length = M;
        break;
    default:
    	printf("Could not resolve problem case\n");
        break;
    }   
}

static void init_Options(int nrhs, const mxArray *prhs[], struct Options *opt)
{
    // Get parameters from the struct opt
    opt->maxiter_outer = mxGetScalar(mxGetField(OPT_IN, 0, "maxiter_outer"));
    opt->maxiter_inner = mxGetScalar(mxGetField(OPT_IN, 0, "maxiter_inner"));
    opt->eps_ds = mxGetScalar(mxGetField(OPT_IN, 0, "eps_ds"));
    opt->eps_pf = mxGetScalar(mxGetField(OPT_IN, 0, "eps_pf"));
    opt->eps_inner = mxGetScalar(mxGetField(OPT_IN, 0, "eps_inner"));
    opt->algorithm = mxGetScalar(mxGetField(OPT_IN, 0, "algorithm")); 
}

static void init_Problem(int nrhs, const mxArray *prhs[], struct Problem *p)
{
   
    mxArray *b1,*b2;
    b1 = mxDuplicateArray(H_IN);
    b2 = mxDuplicateArray(A_IN);
    p->H = mxGetPr(b1);
    p->A = mxGetPr(b2);
    
    p->A_t = vector_alloc(N*M);
    mtx_transpose(p->A,p->A_t,M,N);
    
    // Other vectors
    mxArray *a1,*a2,*a3,*a4,*a5,*a6,*a7;
    a1 = mxDuplicateArray(C_IN);
    a2 = mxDuplicateArray(B_IN);
    a3 = mxDuplicateArray(LB_HAT_IN);
    a4 = mxDuplicateArray(UB_HAT_IN);
    a5 = mxDuplicateArray(LB_IN);
    a6 = mxDuplicateArray(UB_IN);
    a7 = mxDuplicateArray(Z0_IN);
    p->c = mxGetPr(a1);
    p->b = mxGetPr(a2);
    p->lb_hat = mxGetPr(a3);
    p->ub_hat = mxGetPr(a4);
    p->lb = mxGetPr(a5);
    p->ub = mxGetPr(a6);
    p->z0 = mxGetPr(a7);
    
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

static void allocate_ALM(struct Struct_ALM *s)
{
    s->z = vector_alloc(N);
	s->lambda = vector_alloc(M);
	s->temp1_dim_N = vector_alloc(N);
	s->temp2_dim_M = vector_alloc(M);
	s->temp3_dim_M = vector_alloc(M);
	s->z_avg = vector_alloc(N);
	s->summ = vector_alloc(N);
    s->pf_vec = vector_alloc(s->info->pf_vec_length);
    s->A_z = vector_alloc(M);
    s->res->out->ds_vector  = vector_alloc(2);
	s->res->out->pf_vector  = vector_alloc(2);
}

static void allocate_FALM(struct Struct_FALM *s)
{
    s->z = vector_alloc(N);
	s->lambda = vector_alloc(M);
	s->temp1_dim_N = vector_alloc(N);
	s->temp2_dim_M = vector_alloc(M);
	s->temp3_dim_M = vector_alloc(M);
	s->z_avg = vector_alloc(N);
	s->summ = vector_alloc(N);
    s->pf_vec = vector_alloc(s->info->pf_vec_length);
    s->A_z = vector_alloc(M);
    s->res->out->ds_vector  = vector_alloc(2);
	s->res->out->pf_vector  = vector_alloc(2);
    
   	// Different form DGM()
	s->lambda_old = vector_alloc(M);
	s->y1 = vector_alloc(M);
	s->z_ds = vector_alloc(N);
    s->A_z_ds = vector_alloc(M);
}

static void init_DGM(int nrhs, const mxArray *prhs[], struct Struct_DGM *s) 
{
    // Get parameters from the struct INFO
    N = mxGetScalar(mxGetField(INFO_IN, 0, "N"));
    M = mxGetScalar(mxGetField(INFO_IN, 0, "M"));
    
    s->info = malloc(sizeof * s->info);
    init_Info(nrhs, prhs, s->info);
    
    s->opt = malloc(sizeof * s->opt);
    init_Options(nrhs, prhs, s->opt);
    
    s->prob = malloc(sizeof * s->prob);
    init_Problem(nrhs, prhs, s->prob);
     
    s->res = malloc(sizeof * s->res);
	s->res->fopt = 0.0;
	s->res->out = malloc(sizeof * s->res->out); 
    
    allocate_DGM(s);
}

static void init_DFGM(int nrhs, const mxArray *prhs[], struct Struct_DFGM *s) 
{
   
    // Get parameters from the struct INFO
    N = mxGetScalar(mxGetField(INFO_IN, 0, "N"));
    M = mxGetScalar(mxGetField(INFO_IN, 0, "M"));
    
    s->info = malloc(sizeof * s->info);
    init_Info(nrhs, prhs, s->info);
    
    s->opt = malloc(sizeof * s->opt);
    init_Options(nrhs, prhs, s->opt);
    
    s->prob = malloc(sizeof * s->prob);
    init_Problem(nrhs, prhs, s->prob);
     
    s->res = malloc(sizeof * s->res);
	s->res->fopt = 0.0;
	s->res->out = malloc(sizeof * s->res->out); 
    
    allocate_DFGM(s);
 }

static void init_ALM(int nrhs, const mxArray *prhs[], struct Struct_ALM *s) 
{
    // Get parameters from the struct INFO
    N = mxGetScalar(mxGetField(INFO_IN, 0, "N"));
    M = mxGetScalar(mxGetField(INFO_IN, 0, "M"));
    
    s->info = malloc(sizeof * s->info);
    init_Info(nrhs, prhs, s->info);
    
    s->opt = malloc(sizeof * s->opt);
    init_Options(nrhs, prhs, s->opt);
    s->opt->rho = mxGetScalar(mxGetField(OPT_IN, 0, "rho")); 
    
    s->prob = malloc(sizeof * s->prob);
    init_Problem(nrhs, prhs, s->prob);
     
    s->res = malloc(sizeof * s->res);
	s->res->fopt = 0.0;
	s->res->out = malloc(sizeof * s->res->out); 
    
    allocate_ALM(s);
    
    // Take in some other vectors spesial for augmented case
    mxArray *a1,*a2,*a3;
    a1 = mxDuplicateArray(H_HAT_IN);
    a2 = mxDuplicateArray(A2_IN);
    a3 = mxDuplicateArray(RHO_AT_B_IN);
    
    s->H_hat = mxGetPr(a1);
    s->A2 = mxGetPr(a2);
    s->rho_At_b = mxGetPr(a3);
}

static void init_FALM(int nrhs, const mxArray *prhs[], struct Struct_FALM *s) 
{
    // Get parameters from the struct INFO
    N = mxGetScalar(mxGetField(INFO_IN, 0, "N"));
    M = mxGetScalar(mxGetField(INFO_IN, 0, "M"));
    
    s->info = malloc(sizeof * s->info);
    init_Info(nrhs, prhs, s->info);
    
    s->opt = malloc(sizeof * s->opt);
    init_Options(nrhs, prhs, s->opt);
    s->opt->rho = mxGetScalar(mxGetField(OPT_IN, 0, "rho")); 
    
    s->prob = malloc(sizeof * s->prob);
    init_Problem(nrhs, prhs, s->prob);
     
    s->res = malloc(sizeof * s->res);
	s->res->fopt = 0.0;
	s->res->out = malloc(sizeof * s->res->out); 
    
    allocate_FALM(s);
    
    // Take in some other vectors spesial for augmented case
    mxArray *a1,*a2,*a3;
    a1 = mxDuplicateArray(H_HAT_IN);
    a2 = mxDuplicateArray(A2_IN);
    a3 = mxDuplicateArray(RHO_AT_B_IN);
    
    s->H_hat = mxGetPr(a1);
    s->A2 = mxGetPr(a2);
    s->rho_At_b = mxGetPr(a3);
}

static void clean_up_Output_M(struct Output *s)
{
	free_pointer(s->ds_vector);
	free_pointer(s->pf_vector);
    free(s);
 	s = NULL;
}

static void clean_up_Result(struct Result *res)
{   
    clean_up_Output_M(res->out);
    res->zopt = NULL; // Just a pointer (newer allocated space)
    res->lambda1 = NULL; // Just a pointer
    res->lambda2 = NULL; // Just a pointer
    free(res);
 	res = NULL;
}

static void clean_up_DGM_matlab (struct Struct_DGM *s)
{
	free_pointer(s->prob->A_t);
    free(s->opt); s->opt = NULL;
 	free(s->info); s->info = NULL;
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
    clean_up_Result(s->res); 
}

static void clean_up_DFGM_matlab (struct Struct_DFGM *s)
{
	free_pointer(s->prob->A_t);
    free(s->opt); s->opt = NULL;
 	free(s->info); s->info = NULL;
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
    clean_up_Result(s->res); 
    
    // Different from DGM()
	free_pointer(s->lambda1_old);
	free_pointer(s->lambda2_old);
	free_pointer(s->y1);
	free_pointer(s->y2);
	free_pointer(s->z_ds);
    free_pointer(s->A_z_ds);
}

static void clean_up_ALM_matlab (struct Struct_ALM *s)
{
	free_pointer(s->prob->A_t);
    free(s->opt); s->opt = NULL;
 	free(s->info); s->info = NULL;
 	free_pointer(s->z);
	free_pointer(s->lambda);
	free_pointer(s->temp1_dim_N);
	free_pointer(s->temp2_dim_M);
	free_pointer(s->temp3_dim_M);
	free_pointer(s->z_avg);
	free_pointer(s->summ);
	free_pointer(s->pf_vec);
    free_pointer(s->A_z);
    clean_up_Result(s->res); 
}

static void clean_up_FALM_matlab (struct Struct_FALM *s)
{
	free_pointer(s->prob->A_t);
    free(s->opt); s->opt = NULL;
 	free(s->info); s->info = NULL;
 	free_pointer(s->z);
	free_pointer(s->lambda);
	free_pointer(s->temp1_dim_N);
	free_pointer(s->temp2_dim_M);
	free_pointer(s->temp3_dim_M);
	free_pointer(s->z_avg);
	free_pointer(s->summ);
	free_pointer(s->pf_vec);
    free_pointer(s->A_z);
    clean_up_Result(s->res); 
    
    // Different from DGM()
	free_pointer(s->lambda_old);
	free_pointer(s->y1);
	free_pointer(s->z_ds);
    free_pointer(s->A_z_ds);
}

static void make_output(int nlhs, mxArray *plhs[], const struct Result * res)
{
    uint32_t i;

    // zopt
    mxArray *MEX_zopt;
    MEX_zopt = ZOPT_OUT = mxCreateDoubleMatrix(N,1,mxREAL);
    real_t *ZOPT = mxGetPr(MEX_zopt);
    for (i=0;i<N;i++)
        ZOPT[i] = res->zopt[i];

    // fopt
    mxArray *MEX_fopt;
    MEX_fopt=FOPT_OUT = mxCreateDoubleMatrix(1,1,mxREAL);
    real_t *FOPT = mxGetPr(MEX_fopt);
    *FOPT = res->fopt;

    // exitflag
    mxArray *MEX_exitflag;
    MEX_exitflag = EXITFLAG_OUT = mxCreateDoubleMatrix(1,1,mxREAL);
    real_t *EXITFLAG = mxGetPr(MEX_exitflag);
    *EXITFLAG = res->exitflag;

    // output (struct)
    make_Output_struct(nlhs, plhs, res->out);

    // Lambda 1
    mxArray *MEX_lambda1;
    MEX_lambda1 = LAMBDA1_OUT = mxCreateDoubleMatrix(N,1,mxREAL);
    real_t *LAMBDA1 = mxGetPr(MEX_lambda1);
    for (i=0;i<N;i++)
        LAMBDA1[i] = res->lambda1[i];

    if (ALGORITHM <= 4){
        // Lambda 2
        mxArray *MEX_lambda2;
        MEX_lambda2 = LAMBDA2_OUT = mxCreateDoubleMatrix(N,1,mxREAL);
        real_t *LAMBDA2 = mxGetPr(MEX_lambda2);
        for (i=0;i<N;i++)
            LAMBDA2[i] = res->lambda2[i];
    }
}

static void make_output_fail_safe(int nlhs, mxArray *plhs[])
{
    uint32_t i;
    
    // zopt
    mxArray *MEX_zopt;
    MEX_zopt = ZOPT_OUT = mxCreateDoubleMatrix(1,1,mxREAL);
    real_t *ZOPT = mxGetPr(MEX_zopt);
    *ZOPT = -1;

    // fopt
    mxArray *MEX_fopt;
    MEX_fopt=FOPT_OUT = mxCreateDoubleMatrix(1,1,mxREAL);
    real_t *FOPT = mxGetPr(MEX_fopt);
    *FOPT = -1;

    // exitflag
    mxArray *MEX_exitflag;
    MEX_exitflag = EXITFLAG_OUT = mxCreateDoubleMatrix(1,1,mxREAL);
    real_t *EXITFLAG = mxGetPr(MEX_exitflag);
    *EXITFLAG = -1;

    // output (struct)
    mxArray *MEX_output;
    MEX_output = OUTPUT_STRUCT_OUT = mxCreateDoubleMatrix(1,1,mxREAL);
    real_t *OUTPUT = mxGetPr(MEX_output);
    *OUTPUT = -1;

    // Lambda 1
    mxArray *MEX_lambda1;
    MEX_lambda1 = LAMBDA1_OUT = mxCreateDoubleMatrix(N,1,mxREAL);
    real_t *LAMBDA1 = mxGetPr(MEX_lambda1);
    *LAMBDA1 = -1;

    // Lambda 2
    mxArray *MEX_lambda2;
    MEX_lambda2 = LAMBDA2_OUT = mxCreateDoubleMatrix(N,1,mxREAL);
    real_t *LAMBDA2 = mxGetPr(MEX_lambda2);
    *LAMBDA2 = -1;
    
}


static void make_Output_struct(int nlhs, mxArray *plhs[], const struct Output * out)
{
    #define NUM_FIELDS 11   // Number of fields in struct
    #define NAME_LENGT 30   

    #define ITERATION 0
    #define ITERATION_INNER_TOT 1
    #define TIME 2
    #define TIME_TOT_INNER 3
    #define FLAG_LAST_SATISFIED 4
    #define NITER_FEASIBLE_DS 5
    #define NITER_FEASIBLE_PF 6
    #define EXITFLAG_INNER 7
    #define NUM_EXCEEDED_MAX_NITER_INNER 8 
    #define DS_VECTOR 9
    #define PF_VECTOR 10
    uint32_t i; 

    //Allocate memory for the fieldnames
    const char *fieldnames[NUM_FIELDS]; //This will hold field names.
    for (i=0;i<NUM_FIELDS;i++){
        fieldnames[i] = (char*)mxMalloc(NAME_LENGT);
    }       

    // Get all the information
    memcpy((void*)fieldnames[ITERATION],"iterations",sizeof(char)*NAME_LENGT);
    mxArray *MEX_iterations = mxCreateDoubleMatrix (1,1,mxREAL);
    real_t *iterations = mxGetPr(MEX_iterations);
    iterations[0] = out->iterations;

    memcpy((void*)fieldnames[ITERATION_INNER_TOT],"iterations_inner_tot",sizeof(char)*NAME_LENGT);
    mxArray *MEX_iterations_inner_tot = mxCreateDoubleMatrix (1,1,mxREAL);
    real_t *iterations_inner_tot = mxGetPr(MEX_iterations_inner_tot);
    iterations_inner_tot[0] = out->iterations_inner_tot;

    memcpy((void*)fieldnames[TIME],"time",sizeof(char)*NAME_LENGT);
    mxArray *MEX_time = mxCreateDoubleMatrix (1,1,mxREAL);
    real_t *time = mxGetPr(MEX_time);
    time[0] = out->time;

    memcpy((void*)fieldnames[TIME_TOT_INNER],"time_tot_inner",sizeof(char)*NAME_LENGT);
    mxArray *MEX_time_tot_inner = mxCreateDoubleMatrix (1,1,mxREAL);
    real_t *time_tot_inner = mxGetPr(MEX_time_tot_inner);
    time_tot_inner[0] = out->time_tot_inner;

    memcpy((void*)fieldnames[FLAG_LAST_SATISFIED],"flag_last_satisfied",sizeof(char)*NAME_LENGT);
    mxArray *MEX_flag_last_satisfied = mxCreateDoubleMatrix (1,1,mxREAL);
    real_t *flag_last_satisfied = mxGetPr(MEX_flag_last_satisfied);
    flag_last_satisfied[0] = out->flag_last_satisfied;

    memcpy((void*)fieldnames[NITER_FEASIBLE_DS],"niter_feasible_ds",sizeof(char)*NAME_LENGT);
    mxArray *MEX_niter_feasible_ds = mxCreateDoubleMatrix (1,1,mxREAL);
    real_t *niter_feasible_ds = mxGetPr(MEX_niter_feasible_ds);
    niter_feasible_ds[0] = out->niter_feasible_ds;

    memcpy((void*)fieldnames[NITER_FEASIBLE_PF],"niter_feasible_pf",sizeof(char)*NAME_LENGT);
    mxArray *MEX_niter_feasible_pf = mxCreateDoubleMatrix (1,1,mxREAL);
    real_t *niter_feasible_pf = mxGetPr(MEX_niter_feasible_pf);
    niter_feasible_pf[0] = out->niter_feasible_pf;

    memcpy((void*)fieldnames[EXITFLAG_INNER],"exitflag_inner",sizeof(char)*NAME_LENGT);
    mxArray *MEX_exitflag_inner = mxCreateDoubleMatrix (1,1,mxREAL);
    real_t *exitflag_inner = mxGetPr(MEX_exitflag_inner);
    exitflag_inner[0] = out->exitflag_inner;

    memcpy((void*)fieldnames[NUM_EXCEEDED_MAX_NITER_INNER],"num_exceeded_max_niter_inner",sizeof(char)*NAME_LENGT);
    mxArray *MEX_num_exceeded_max_niter_inner = mxCreateDoubleMatrix (1,1,mxREAL);
    real_t *num_exceeded_max_niter_inner = mxGetPr(MEX_num_exceeded_max_niter_inner);
    num_exceeded_max_niter_inner[0] = out->num_exceeded_max_niter_inner;

    memcpy((void*)fieldnames[PF_VECTOR],"pf_vector",sizeof(char)*NAME_LENGT);
    mxArray *MEX_pf_vector = mxCreateDoubleMatrix (out->iterations,1,mxREAL);
    real_t *pf_vector = mxGetPr(MEX_pf_vector);
    for (i=0;i<out->iterations;i++)
        pf_vector[i] = out->pf_vector[i];

    memcpy((void*)fieldnames[DS_VECTOR],"ds_vector",sizeof(char)*NAME_LENGT);
    mxArray *MEX_ds_vector = mxCreateDoubleMatrix (out->iterations,1,mxREAL);
    real_t *ds_vector = mxGetPr(MEX_ds_vector);
    for (i=0;i<out->iterations;i++)
        ds_vector[i] = out->ds_vector[i];

    // Create the struct       
    OUTPUT_STRUCT_OUT = mxCreateStructMatrix(1,1,NUM_FIELDS,fieldnames);

    //Deallocate memory for the fieldnames
    for (i=0;i<NUM_FIELDS;i++){
        mxFree((void*) fieldnames[i] );
    }

    // Set the info into the struct
    mxSetFieldByNumber(OUTPUT_STRUCT_OUT,0,ITERATION, MEX_iterations);
    mxSetFieldByNumber(OUTPUT_STRUCT_OUT,0,ITERATION_INNER_TOT, MEX_iterations_inner_tot);
    mxSetFieldByNumber(OUTPUT_STRUCT_OUT,0,TIME, MEX_time);
    mxSetFieldByNumber(OUTPUT_STRUCT_OUT,0,TIME_TOT_INNER, MEX_time_tot_inner);
    mxSetFieldByNumber(OUTPUT_STRUCT_OUT,0,FLAG_LAST_SATISFIED, MEX_flag_last_satisfied);
    mxSetFieldByNumber(OUTPUT_STRUCT_OUT,0,NITER_FEASIBLE_DS, MEX_niter_feasible_ds);
    mxSetFieldByNumber(OUTPUT_STRUCT_OUT,0,NITER_FEASIBLE_PF, MEX_niter_feasible_pf);
    mxSetFieldByNumber(OUTPUT_STRUCT_OUT,0,EXITFLAG_INNER, MEX_exitflag_inner);
    mxSetFieldByNumber(OUTPUT_STRUCT_OUT,0,NUM_EXCEEDED_MAX_NITER_INNER, MEX_num_exceeded_max_niter_inner);
    mxSetFieldByNumber(OUTPUT_STRUCT_OUT,0,DS_VECTOR, MEX_ds_vector);
    mxSetFieldByNumber(OUTPUT_STRUCT_OUT,0,PF_VECTOR, MEX_pf_vector);
}