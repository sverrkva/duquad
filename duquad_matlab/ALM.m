function [zopt,fopt,exitflag,output] = ALM(varargin)
% H has to be symetric and positiv definite
% lb_hat must be equal to ub_hat
% if Az = b, lab_hat = ub_hat = [0;0;0...]

if nargin < 10
    opt = [];
    if nargin < 9
        z0 = [];
        if nargin < 8
            ub = [];
            if nargin < 7
                lb = [];
                if nargin < 6
                    ub_hat = [];
                    if nargin < 5
                        lb_hat = [];
                        if nargin < 4
                            b = [];
                            if nargin < 3
                                A = [];
                                if nargin < 2
                                    error('Incorrect number of input arguments in ALM().')                 
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

[H,c] = deal(varargin{1:2});

if nargin > 2
    A = varargin{3};
end    
if nargin > 3
    b = varargin{4};
end    
if nargin > 4
    lb_hat = varargin{5};
end  
if nargin > 5
    ub_hat = varargin{6};
end
if nargin > 6
    lb = varargin{7};
end    
if nargin > 7
    ub = varargin{8};
end  
if nargin > 8
    z0 = varargin{9};
end
if nargin > 9
    opt = varargin{10};
end
if nargin > 10
    error('Incorrect number of input arguments in ALM.')
end    

n = size(H,2);
m = size(A,1);

if isempty(A)
    error('No boint in using DGM without defining any linear_inequalities')
end
if isempty(b)
    b = zeros(m,1);     % not optimal to calculate with this value when empty
end
if isempty(lb)
    lb_inf = true;
else
    lb_inf = false;
end
if isempty(ub)
    ub_inf = true;
else
    ub_inf = false;
end
if isempty(z0)
    z0 = zeros(n,1);
end
if isempty(opt)
    % Set default values
    opt.maxiter_outer = 10000;
    opt.maxiter_inner = 100;
    opt.eps_ds = 0.001;
    opt.eps_pf = 0.05;
    opt.eps_inner = 0.0001;
    opt.algorithm = 5;
    opt.rho = 1;
end

% initial point must be feasible
if ~ub_inf
    z0 = min(ub,z0);
end
if ~lb_inf
    z0 = max(lb,z0);
end

% Calculating alpha
rho = opt.rho;
%Ld = rho;
alpha = rho;

% Setting options for the inner problem
optFG.maxiter = opt.maxiter_inner;
optFG.eps = opt.eps_inner;

% initialize lambda
lambda = zeros(m,1);

% Making constant vectors and matrices
% b = b + lb_hat 
if isempty(ub_hat) || isempty(lb_hat)
    warning('lb_hat and/or ub_hat is unconstrained. Running the case Az = b')
elseif (ub_hat ~= lb_hat)
    error('AML will only run with linear equalities');
else
    b = b + lb_hat;       % lb_hat sholud be zero
end

A_t = A';
A2 = A_t*A;
rho_At_b = rho * A_t * b;


% New H (H_hat). NOTE: when rho in not constant H_hat has to be calculated every time 
H_hat = H + rho*A2; % Condtant rho
H_eig = eig(H_hat);
info.L = max(H_eig);
info.mu = min(H_eig);

% solve inner problem first time - with H_hat and c_hat
c_hat = c - rho_At_b; %c + A'*lambda - rho*Atb;
[z,~, iterations_inner] = FGM(H_hat,c_hat,lb,ub,z0,optFG,info);
iterations_inner_tot = iterations_inner;
A_z = A*z; % used later

% find the value of the Lagrangian first time
dual_value = dual_obj_c2(z,H,c,lambda,A_z,b,rho);
dual_value_diff = opt.eps_ds + 1;

% initialize other necessary variables
pf = opt.eps_pf + 1;
summ = z;
niter_feasible_ds = 0;
niter_feasible_pf = 0;
exitflag = 1;
niter = 1;


% #########################################################################
%			START WHILE-LOOP
% #########################################################################
while or(dual_value_diff > opt.eps_ds , pf > opt.eps_pf) 

    % check if the maximum number of iterations is reached
    if niter >= opt.maxiter_outer
        warning('reached maximum number of iterations in ALM\n')
        niter_feasible_ds = niter_feasible_ds - 1;
        niter_feasible_pf = niter_feasible_pf - 1;
        exitflag = 2;
        break;
    end
    
    if dual_value_diff < opt.eps_ds
        niter_feasible_ds = niter_feasible_ds + 1;
        last_eps_ds = 1;
    else
        last_eps_ds = 0;
    end

    if pf < opt.eps_pf
        niter_feasible_pf = niter_feasible_pf + 1;
    end

    % *********************************************************************
	%			Finding the next lambda
	% *********************************************************************
    
    lambda = lambda + alpha * (A_z - b);
    % NOTE: no projectiond of lambda

    % *********************************************************************
	%			Solving the inner problem
	% *********************************************************************
    % obtaining z from the inner problem
    c_hat = c + A_t*lambda - rho_At_b; % rho *A' * b is constant
    [z,~,iterations_inner] = FGM(H_hat,c_hat,lb,ub,z,optFG,info);
    iterations_inner_tot = iterations_inner_tot + iterations_inner;
    A_z = A*z; % only do this operation here
    
    % *********************************************************************
	%			Finding dual_value_diff
	% *********************************************************************
    dual_value_new = dual_obj_c2(z,H,c,lambda,A_z,b,rho);
    dual_value_diff = abs(dual_value_new - dual_value);
    dual_value = dual_value_new;

    % *********************************************************************
	%			Finding pf
	% *********************************************************************
    % *** LAST z in stopping criteria ***
    if opt.algorithm == 5
        pf_vec = A_z - b;
        % NOTE: no projectiond

    % *** AVERAGE z in Pf stopping criteria ***
    elseif opt.algorithm == 6
        summ = summ + z;
        z_avg = summ/(niter+1);
        pf_vec = A*z_avg - b;
        % NOTE: no projectiond
    else
        error('Choose opt.algorithm = 5 / 6 when running ALM')
    end

    % take the norm
    pf = norm2(pf_vec);
    
    % set ds and pf vectors for output
    output.ds_vector(niter) = dual_value_diff;
    output.pf_vector(niter) = pf;
    
    % *********************************************************************
	%			Update the variables
	% *********************************************************************
    niter = niter+1;
    
end % #### END WHILE-LOOP ####

if opt.algorithm == 5
    zopt = z;
elseif opt.algorithm == 6
    zopt = z_avg;
else
    error('Choose opt.algorithm = 5 / 6 when running ALM')
end

fopt = obj(zopt,H,c);   
output.iterations = niter-1;
output.iterations_inner_tot = iterations_inner_tot;
output.niter_feasible_ds = niter_feasible_ds+1;
output.niter_feasible_pf = niter_feasible_pf+1;  
     
end
    
% Local functions
function value = obj(z, H, c) 
    value = 0.5* (z'*H*z) + c'*z;
end


function value = dual_obj_c2(z,H,c,lambda,A_z,b,rho)
    Az_b = A_z - b;
    value = 0.5* (z'*H*z) + c'*z...
        + lambda'*(Az_b)...
        + 0.5*rho*(norm2(Az_b)^2);
end


function printv(v,length)
    fprintf('[ ');
    for i=1:length
        fprintf('%f  ', v(i));
    end
    fprintf(']\n');
end
