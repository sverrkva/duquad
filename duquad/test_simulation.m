% testing the simulation function
clear all


par.run_quadprog = 0;   % Matlab
par.run_gurobi = 1;     % Gurobi
par.a1 = 1;             % DGM last
par.a2 = 0;             % DGM avg
par.a3 = 1;             % DFGM last
par.a4 = 0;             % DFGM avg
% 
% par.run_quadprog = 1;   % Matlab
% par.run_gurobi = 1;     % Gurobi
% par.a1 = 1;             % DGM last
% par.a2 = 0;             % DGM avg
% par.a3 = 1;             % DFGM last
% par.a4 = 1;             % DFGM avg

% Options
maxiter_outer = 4000;
maxiter_inner=100;
eps_ds=0.001;
eps_pf = 0.5;

% Set options for DGM()
opt_DGM.maxiter_outer = maxiter_outer;
opt_DGM.maxiter_inner = maxiter_inner;
opt_DGM.eps_ds = eps_ds;
opt_DGM.eps_pf = eps_pf;
opt_DGM.eps_inner = 0.0001;
par.opt_DGM = opt_DGM;

% Set options for DFGM()
opt_DFGM.maxiter_outer = maxiter_outer;
opt_DFGM.maxiter_inner = maxiter_inner;
opt_DFGM.eps_ds = eps_ds;
opt_DFGM.eps_pf = eps_pf;
opt_DFGM.eps_inner = 0.000001;
par.opt_DFGM = opt_DFGM;

% Set problem parameters
par.n = 1200;
par.m = round(2*par.n/3);
par.N_simulations = 1;
par.gamma = 10;


% ---- Simulate ----
fprintf('\n------------ START NEW SIMULATION --------------\n');
% Matlab code

test = simulate_matlab(par,[],false);

% problem = test{1}.problem;
% % C code
% test = simulate_c(par,problem,false);



fprintf('\n\n ------------ Simulation finished ----------------\n\n');
play_sound(0.3)


% for i=1:1
%    
%    test{i} = simulate(par);
%    fprintf('\n ------- finished with n = %d --------\n', par.n);
%    save('results/test_n');
% end
% 
% for i=1:length(test)
%    info{i} = extract_info(test{i});
% end
% 
% for i=1:length(test)-1
%     info{i} = extract_info(test{i});
% end

