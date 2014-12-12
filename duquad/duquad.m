function [zopt,fopt,exitflag,output,lambda1,lambda2] = duquad(varargin)
%
% *** DuQuad - Quadratic Programming Optimization ***
% 
%  Attempts to solve the quadratic programming problem:
%
%  --------------------------------------
%  |   min f(z) = 0.5*z'*H*z + c'z,     |
%  |    z                               |
%  |       s.t.                         |
%  |       lb_hat <= Az - b <= ub       |
%  |       lb <= z <= ub                |
%  --------------------------------------         
%
% where z = [z1 z2 ... zN]^T
%
% INPUTS:
% H:        Hessian matrix (must be positive definite and symmetric)
% c:        gradient vector
% A:        linear constraints matrix
% b:        linear constraints vector
% lb_hat:   lower bound for the linear constraints
% ub_hat:   upper bound for the linear constraints
% lb:       lower bound for optimization variable z
% ub:       upper bound for optimization variable z
% z0:       initial point
% opt:      struct containing options, see OPTIONS
% 
% OUTPUTS:
% zopt:     optimal solution
% fopt:     optimal value (f(xopt))
% niter:    number of iteration for finding the solution
% output:   struct containing info from the problem solving
% 
% USAGE:
% [zopt,fopt,niter] = duquad(H,c);
% [zopt,fopt,niter] = duquad(H,c,A,b);
% [zopt,fopt,niter] = duquad(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt);
% [zopt,fopt,niter] = duquad(H,c,A,b,[],[],[],[],[],opt);
% 
% OPTIONS:
% opt.maxiter_outer:  Maximum number of iterations in the outer loop
% opt.maxiter_inner:  Maximum number of iterations in the inner loop
% opt.eps_ds:         Tolerance for dual suboptimality
% opt.eps_pf:         Tolerance for primal feasibility
% opt.eps_inner:      Tolerance for primal feasibility in the inner problem
% opt.rho:            Penalty parameter used in ALM and FALM
% opt.algorithm:      Spesifies the algoritm used to solve the problem. Values: 
%                     1: DGM last
%                     2: DGM avg
%                     3: DFGM last
%                     4: DFGM avg
%                     5: ALM last
%                     6: ALM avg
%                     7: FALM last
%                     8: FALM avg
%
% AUTHOR:
% Sverre Kvamme
% 
% Project webpage and more information:
% http://sverrkva.github.io/duquad/


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
                                    error('Incorrect number of input arguments.')                 
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
    error('Incorrect number of input arguments.')
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
    lb_is_inf = true;
else
    lb_is_inf = false;
end
if isempty(ub)
    ub_is_inf = true;
else
    ub_is_inf = false;
end
if isempty(z0)
    z0 = zeros(n,1);
end
if isempty(opt)
    % Set default values
    opt.maxiter_outer = 1000;
    opt.maxiter_inner = 100;
    opt.eps_ds = 0.001;
    opt.eps_pf = 0.05;
    opt.eps_inner = 0.000001;
    opt.algorithm = 3;
    opt.rho = 1;
end


% Desciding the category/case of the problem 
if isempty(ub_hat) && isempty(lb_hat)
    error('both lb_hat and ub_hat is unconstrained')
elseif isempty(ub_hat)
    problem_case = 4;
elseif isempty(lb_hat)
    problem_case = 3;
elseif lb_hat == ub_hat
    problem_case = 2;
else
    problem_case = 1;
end

INFO = struct(...
'N',[]...
,'M',[]...
,'lb_is_inf',[]...
,'ub_is_inf',[]...
,'problem_case',[]...
,'eigH_min',[]...
,'eigH_max',[]...
,'Ld',[]...
);

INFO.N = n;
INFO.M = m;
INFO.lb_is_inf = lb_is_inf;
INFO.ub_is_inf = ub_is_inf;
INFO.problem_case = problem_case;

A_t = A'; % only do this operation once

if opt.algorithm <= 4
    eigH = eig(H);
    eigH_min = min(eigH);
    eigH_max = max(eigH);
    INFO.eigH_min = eigH_min;
    INFO.eigH_max = eigH_max;

    % check for errors in H and A
    % check if all elements in lb <= ub

    switch problem_case

        %*************************************************
        case 1  % ub_hat and lb_hat exists
        %*************************************************   

            Ld = norm([A;-A],2)^2 / eigH_min;

        % *************************************************    
        case 2  % ub_hat = lb_hat
        % *************************************************

            Ld = norm(A,2)^2 / eigH_min;

        %*************************************************
        case 3  % lb_hat = -inf
        %*************************************************

             Ld = norm(A,2)^2 / eigH_min;

        %*************************************************
        case 4  % ub_hat = inf
        %*************************************************

             Ld = norm(A,2)^2 / eigH_min;

    end

        INFO.Ld = Ld;
        [zopt,fopt,exitflag,output,lambda1,lambda2] = main(H,c,A_t,b,lb_hat,ub_hat,lb,ub,z0,opt,INFO);
        % NOTE: Sending in transposed of A, A', because when reading input in
        % mex this is more easy. Maybe do something about this
        if (exitflag == -1)
            error('Error trying to run the program. See stderr for error messages.');
        elseif (exitflag == 2)
            fprintf('Algorithm %d reached maximum number of iterations.\n', opt.algorithm);
        end
else
    % Making constant vectors and matrices
    % b = b + lb_hat 
    if isempty(ub_hat) || isempty(lb_hat)
        warning('lb_hat and/or ub_hat is unconstrained. Running the case Az = b')
    elseif (ub_hat ~= lb_hat)
        error('ALM and FALM will only run with linear equalities');
    else
        b = b + lb_hat;       % lb_hat sholud be zero
    end

    A2 = A_t * A;
    rho_At_b = opt.rho * A_t * b; % Does inverting A more than one time
    H_hat = H + opt.rho*A2; % Condtant rho
    H_eig = eig(H_hat);
    INFO.eigH_min = min(H_eig);
    INFO.eigH_max = max(H_eig);
    [zopt,fopt,exitflag,output,lambda1]...
        = main(H,c,A_t,b,lb_hat,ub_hat,lb,ub,z0,opt,INFO,H_hat,A2,rho_At_b);
    % NOTE: Sending in transposed of A, A', because when reading input in
    % mex this is more easy. Maybe do something about this
    if (exitflag == -1)
        error('Error trying to run the program. See stderr for error messages.');
    elseif (exitflag == 2)
        fprintf('Algorithm %d reached maximum number of iterations.\n', opt.algorithm);
    end
    
end

end



