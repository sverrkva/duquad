% Testing the memory

% Parameters
maxiter_outer =5000;
maxiter_inner=100;
eps_ds=0.0001;
eps_pf = 0.05;


% Set options
opt.maxiter_outer = maxiter_outer;
opt.maxiter_inner = maxiter_inner;
opt.eps_ds = eps_ds;
opt.eps_pf = eps_pf;
opt.eps_inner = 0.000001;
opt.rho = 1;


%% QUADPROG
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
fprintf('\n-- Quadprog --\n');
tic;
[zopt,fopt,exitflag,output,LAMBDA]...
    = quadprog(H,c,[],[],A,b+ub_hat,lb,ub,z0,opts);
time = toc;
print_result(zopt,fopt,exitflag,output,time);




%% *** FALM last ***

fprintf('\n\n****** FALM last ******');
opt.algorithm = 7;

% C
fprintf('\n-- C --\n');
tic;
[zopt,fopt,exitflag,output]...
    = QP_toolbox_C(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt);
time = toc;

print_result(zopt,fopt,exitflag,output,time);

% Matlab
fprintf('\n-- Matlab --\n')
tic;
[zopt,fopt,exitflag,output]...
    = FALM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt);
time = toc;
print_result(zopt,fopt,exitflag,output,time);