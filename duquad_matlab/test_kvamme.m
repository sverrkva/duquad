


%%
% KVAMME


% Parameters
maxiter_outer =1000;
maxiter_inner=1000;
eps_ds=0.001;
eps_pf = 0.01;
eps_inner=0.0001;

% *** DGM last ***
opt_last_DGM = set_options_dual(...
'maxiter_outer',maxiter_outer,...
'maxiter_inner',maxiter_inner,...
'eps_ds',eps_ds,...
'eps_pf',eps_pf,...
'eps_inner',eps_inner,...
'algorithm',(1)...
);

tic;
[zopt,fopt,exitflag,output]...
    = DGM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_last_DGM);
time_last = toc;
fprintf('DGM last finished\n');
fprintf('f: %f\n',fopt);
fprintf('exitflag: %d\n',exitflag);
fprintf('iterations: %d\n',output.iterations);
fprintf('iterations inner tot: %d\n',output.iterations_inner_tot);
fprintf('niter_feasible_ds: %d\n',output.niter_feasible_ds);
fprintf('niter_feasible_pf: %d\n',output.niter_feasible_pf);

%zopt

% clear model;
% model.Q = sparse(H*0.5);
% model.obj = c;
% model.A = sparse(A);
% model.rhs = [b];
% model.sense = '=';
% model.lb = lb;
% model.ub = ub;
% model.start = z0;
% clear params;
% params.outputflag = 0;
% params.BarIterLimit =10000;
% res = gurobi(model,params)


print_problem_to_file(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_last_DGM)



%%
% *** DGM avg ***
opt_avg_DGM = set_options_dual(...
'maxiter_outer',maxiter_outer,...
'maxiter_inner',maxiter_inner,...
'eps_ds',eps_ds,...
'eps_pf',eps_pf,...
'eps_inner',eps_inner,...
'algorithm',(2)...
);

tic;
[zopt_avg_DGM,fopt_avg_DGM,niter_avg_DGM]...
     = DGM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_avg_DGM);
time_avg = toc;
fprintf('DGM avg finished\n');
fprintf('f: %f\n',fopt_avg_DGM);
fprintf('iter: %d\n\n',niter_avg_DGM);

%print_problem_to_file(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_avg_DGM)

%%

% Parameters
maxiter_outer =10000;
maxiter_inner=1000;
eps_ds=0.0001;
eps_pf = 0.0001;
eps_inner=0.000001;

% *** DFGM last ***
opt_last_DFGM = set_options_dual(...
'maxiter_outer',maxiter_outer,'maxiter_inner',maxiter_inner,...
'eps_ds',eps_ds,'eps_pf',eps_pf,'eps_inner',eps_inner,...
'algorithm',(3)...
);

print_problem_to_file(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_last_DFGM)
tic;
[zopt_last_DFGM,fopt_last_DFGM,niter_last_DFGM]...
     = DFGM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_last_DFGM);
time_last_DFGM = toc;
fprintf('DFGM last finished\n');
fprintf('f: %f\n',fopt_last_DFGM);
fprintf('iter: %d\n\n',niter_last_DFGM);


%%
% *** DFGM avg ***
opt_avg_DFGM = set_options_dual(...
'maxiter_outer',maxiter_outer,'maxiter_inner',maxiter_inner,...
'eps_ds',eps_ds,'eps_pf',eps_pf,'eps_inner',eps_inner,...
'algorithm',(4)...
);
tic;
[zopt_avg_DFGM,fopt_avg_DFGM,niter_avg_DFGM]...
     = DFGM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_avg_DFGM);
time_avg_DFGM = toc;
fprintf('DFGM avg finished\n');
fprintf('f: %f\n',fopt_avg_DFGM);
fprintf('iter: %d\n\n',niter_avg_DFGM);


%% *** ALM last ***
% opt_last_ALM = set_options_dual(...
% 'maxiter_outer',maxiter_outer,...
% 'maxiter_inner',maxiter_inner,...
% 'eps_ds',eps_ds,...
% 'eps_pf',eps_pf,...
% 'eps_inner',eps_inner,...
% 'algorithm',(5)...
% );
% 
% opt_last_ALM.rho = 1;
% 
% tic;
% [zopt_last_ALM,fopt_last_ALM,niter_last_ALM]...
%     = ALM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_last_ALM);
% time_last = toc;
% fprintf('ALM last finished\n');
% fprintf('f: %f\n',fopt_last_ALM);
% fprintf('iter: %d\n\n',niter_last_ALM);
% %zopt_last_DGM
% 
% % *** FALM last ***
% opt_last_FALM = set_options_dual(...
% 'maxiter_outer',maxiter_outer,...
% 'maxiter_inner',maxiter_inner,...
% 'eps_ds',eps_ds,...
% 'eps_pf',eps_pf,...
% 'eps_inner',eps_inner,...
% 'algorithm',(7)...
% );
% 
% opt_last_FALM.rho = 15;
% 
% tic;
% [zopt_last_FALM,fopt_last_FALM,niter_last_FALM]...
%     = FALM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt_last_FALM);
% time_last = toc;
% fprintf('FALM last finished\n');
% fprintf('f: %f\n',fopt_last_FALM);
% fprintf('iter: %d\n\n',niter_last_FALM);
% %zopt_last_DGM
% 
% 



