% problems from repository

addpath('/home/sverre/Dropbox/NTNU/Master/QP-Test-Problems-master/MAT_Files')
load('CVXQP1_S'); % f = 0.115907181D+05


% OBJECTIVE
H   = full(Q);
c = c;
% CONSTRAINTS
% linear Inequalities
A       = full(A); 
b       = zeros(size(A,1),1);
lb_hat  = rl;
ub_hat  = ru;
% decision Variable Bounds
lb  = lb;
ub  = ub;

z0 = zeros(size(H,1),1);


opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
    = quadprog(H,c,[A;-A],[b+ub_hat;-b-lb_hat],[],[],lb,ub,z0,opts);
time_quad = toc;
fprintf('quadprog finished\n');
fprintf('f: %f\n',fopt_quadprog);
fprintf('iter: %d\n\n',OUTPUT.iterations);

clear model;
model.Q = sparse(H*0.5);
model.obj = c;
model.A = sparse([A;-A]);
model.rhs = [b+ub_hat;-b-lb_hat];
model.sense = '<';
model.lb = lb;
model.ub = ub;
model.start = z0;


res_gurobi = gurobi(model,params)
fprintf('f: %0.8f\n',res_gurobi.objval);


% Parameters
maxiter_outer =2000;
maxiter_inner=10000;
eps_DS_outer=0.0001;
eps_PF_outer = 0.8;
eps_inner=0.0001;

% *** DGM last ***
opt = set_options_dual(...
'maxiter_outer',maxiter_outer,...
'maxiter_inner',maxiter_inner,...
'eps_DS_outer',eps_DS_outer,...
'eps_PF_outer',eps_PF_outer,...
'eps_inner',eps_inner,...
'algorithm',(5)...
);

opt.rho = 1.5;
tic;
[zopt_last_DGM,fopt_last_DGM,niter_last_DGM]...
    = ALM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt);
time_last = toc;
fprintf('DGM last finished\n');
fprintf('f: %f\n',fopt_last_DGM);
fprintf('iter: %d\n\n',niter_last_DGM);








