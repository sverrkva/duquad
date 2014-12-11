% KVAMME

% algorith:
% 1 = DGM last
% 2 = DGM avg
% 3 = DFGM last
% 4 = DFGM avg

addpath('/home/sverre/Dropbox/NTNU/Master/code/thesis/MEX')
%load('problem_case_2_DFGM_not_equal_to_matlab_2','H','c','A','b','lb_hat','ub_hat','lb','ub','z0');
%print_problem_to_file(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt)


% Parameters
maxiter_outer =5000;
maxiter_inner=100;
eps_ds=0.0001;
eps_pf = 0.05;


% Set options for DGM()
opt_DGM.maxiter_outer = maxiter_outer;
opt_DGM.maxiter_inner = maxiter_inner;
opt_DGM.eps_ds = eps_ds;
opt_DGM.eps_pf = eps_pf;
opt_DGM.eps_inner = 0.001;

% Set options for DFGM()
opt_DFGM.maxiter_outer = maxiter_outer;
opt_DFGM.maxiter_inner = maxiter_inner;
opt_DFGM.eps_ds = eps_ds;
opt_DFGM.eps_pf = eps_pf;
opt_DFGM.eps_inner = 0.000001;
opt_DFGM.rho = 1;



%% *** DGM last ***

fprintf('\n****** DGM last ******');
opt_DGM.algorithm = 1;

% C
fprintf('\n-- C --\n');
tic;
[zopt,fopt,exitflag,output]...
    = QP_toolbox_C(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DGM);
time = toc;
print_result(zopt,fopt,exitflag,output,time);

dgm_last_C = output;

% Matlab
fprintf('\n-- Matlab --\n')
tic;
[zopt,fopt,exitflag,output]...
    = DGM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DGM);
time = toc;
print_result(zopt,fopt,exitflag,output,time);

dgm_last_M = output;

% close all;
% figure(); hold on;grid on;
% title('pf')
% plot(dgm_last_M.pf_vector,'r')
% plot(dgm_last_C.pf_vector,'g')
% figure(); hold on;grid on;
% title('ds')
% plot(dgm_last_M.ds_vector,'r')
% plot(dgm_last_C.ds_vector,'g')



%% *** DGM avg ***

fprintf('\n\n****** DGM avg ******');
opt_DGM.algorithm = 2;

% C
fprintf('\n-- C --\n');
tic;
[zopt,fopt,exitflag,output]...
    = QP_toolbox_C(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DGM);
time = toc;

print_result(zopt,fopt,exitflag,output,time);

% Matlab
fprintf('\n-- Matlab --\n')
tic;
[zopt,fopt,exitflag,output]...
    = DGM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DGM);
time = toc;
print_result(zopt,fopt,exitflag,output,time);



%% *** DFGM last ***
fprintf('\n\n****** DFGM last ******');
opt_DFGM.algorithm = 3;

% C
fprintf('\n-- C --\n');
tic;
[zopt,fopt,exitflag,output]...
    = QP_toolbox_C(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;

print_result(zopt,fopt,exitflag,output,time);

% Matlab
fprintf('\n-- Matlab --\n')
tic;
[zopt,fopt,exitflag,output]...
    = DFGM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;
print_result(zopt,fopt,exitflag,output,time);


%% *** DFGM avg ***
fprintf('\n\n****** DFGM avg ******');
opt_DFGM.algorithm = 4;

% C
fprintf('\n-- C --\n');
tic;
[zopt,fopt,exitflag,output]...
    = QP_toolbox_C(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;

print_result(zopt,fopt,exitflag,output,time);

% Matlab
fprintf('\n-- Matlab --\n')
tic;
[zopt,fopt,exitflag,output]...
    = DFGM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;
print_result(zopt,fopt,exitflag,output,time);


%% *** ALM last ***




fprintf('\n\n****** ALM last ******');
opt_DFGM.algorithm = 5;

% C
fprintf('\n-- C --\n');
tic;
[zopt,fopt,exitflag,output]...
    = QP_toolbox_C(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;

print_result(zopt,fopt,exitflag,output,time);

% Matlab
fprintf('\n-- Matlab --\n')
tic;
[zopt,fopt,exitflag,output]...
    = ALM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;
print_result(zopt,fopt,exitflag,output,time);

%% *** ALM avg ***
fprintf('\n\n****** ALM avg ******');
opt_DFGM.algorithm = 6;

% C
fprintf('\n-- C --\n');
tic;
[zopt,fopt,exitflag,output]...
    = QP_toolbox_C(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;

print_result(zopt,fopt,exitflag,output,time);

% Matlab
fprintf('\n-- Matlab --\n')
tic;
[zopt,fopt,exitflag,output]...
    = ALM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;
print_result(zopt,fopt,exitflag,output,time);

%% *** Gurobi ***
fprintf('\n\n****** GUROBI ******');

% Options Gurobi
clear params;
params.outputflag = 0;
params.Threads = 1;

clear model;
model.Q = sparse(H*0.5);
model.obj = c;
model.A = sparse(A);
model.rhs = b+ub_hat;
model.sense = '=';
model.lb = lb;
model.ub = ub;
model.start = z0;

tic;
res_gurobi = gurobi(model,params)
time = toc;
fprintf('\n*** Gurobi finished ***\n');
fprintf('fopt: %f\n',res_gurobi.objval);
fprintf('exitflag: %s\n',res_gurobi.status);
fprintf('runtime: %f\n',time); 
res_gurobi



% *** FALM last ***

% Parameters
maxiter_outer =300;
maxiter_inner=100;
eps_ds=0.0001;
eps_pf = 0.01;

% Set options for DFGM()
opt_DFGM.maxiter_outer = maxiter_outer;
opt_DFGM.maxiter_inner = maxiter_inner;
opt_DFGM.eps_ds = eps_ds;
opt_DFGM.eps_pf = eps_pf;
opt_DFGM.eps_inner = 0.00001;
opt_DFGM.rho = 1;



fprintf('\n\n****** FALM last ******');
opt_DFGM.algorithm = 7;

% C
fprintf('\n-- C --\n');
tic;
[zopt,fopt,exitflag,output]...
    = QP_toolbox_C(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;

print_result(zopt,fopt,exitflag,output,time);

% Matlab
fprintf('\n-- Matlab --\n')
tic;
[zopt,fopt,exitflag,output]...
    = FALM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;
print_result(zopt,fopt,exitflag,output,time);

%% *** FALM avg ***
fprintf('\n\n****** FALM avg ******');
opt_DFGM.algorithm = 8;

% C
fprintf('\n-- C --\n');
tic;
[zopt,fopt,exitflag,output]...
    = QP_toolbox_C(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;

print_result(zopt,fopt,exitflag,output,time);

% Matlab
fprintf('\n-- Matlab --\n')
tic;
[zopt,fopt,exitflag,output]...
    = FALM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_DFGM);
time = toc;
print_result(zopt,fopt,exitflag,output,time);



