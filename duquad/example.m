% *** DuQuad Examplefile ***

% HELP: type 'help duquad' in console.

% Make a QP problem with constraints

H = [11 4 ; 4 22];  % Hessian matrix
c = [3 ; 4];        % gradient vector
A = [1 1;2 1];      % linear constraints matrix
b = [2 ; 3];        % linear constraints vector
lb_hat = [-2 ; -2]; % lower bound for the linear constraints
ub_hat = [2 ; 2];   % upper bound for the linear constraints
lb = [-1 ; -2];     % lower bound for optimization variable z
ub = [0.5 ; 2];     % upper bound for optimization variable z
z0 = [0.5 ; -0.5];  % initial point


% Set options for the toolbox. (default values in documentaion)

opt.maxiter_outer = 1000;   % Maximum number of iterations in the outer loop
opt.maxiter_inner = 100;    % Maximum number of iterations in the inner loop
opt.eps_ds = 0.000001;       %  olerance for dual suboptimality
opt.eps_pf = 0.001;         % Tolerance for primal feasibility
opt.eps_inner = 0.000001;   % Tolerance for primal feasibility in the inner problem
opt.rho = 1;                % Penalty parameter used in ALM and FALM
opt.algorithm = 3;          % Spesifies the algoritm used to solve the problem. 

%   Algorithm values:
%   1: DGM last
%   2: DGM avg
%   3: DFGM last
%   4: DFGM avg
%   5: ALM last
%   6: ALM avg
%   7: FALM last
%   8: FALM avg


% Run DuQuad

[zopt,fopt,exitflag,output,lambda1,lambda2] = duquad(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt);

fprintf('\nf: %f\n',fopt);
fprintf('iterations: %d\n',output.iterations);
output

%% Check result with quadprog

opts = optimoptions('quadprog','Algorithm','active-set','Display','off');

[zopt,fopt,EXITFLAG,OUTPUT]...
    = quadprog(H,c,[A;-A],[b+ub_hat;-b-lb_hat],[],[],lb,ub,z0,opts);

fprintf('\nQuadprog result:\n');
fprintf('f: %f\n',fopt);

EXITFLAG
