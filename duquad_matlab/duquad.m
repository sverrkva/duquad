function [zopt,fopt,exitflag,output] = duquad(varargin)
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

addpath('src');

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
n = size(H,2);

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

if isempty(lb)
    lb = [];
end
if isempty(ub)
    ub = [];
end
if isempty(z0)
    z0 = zeros(n,1);
end
if isempty(A)
    warning('No linear constraints. Running FGM')
    [zopt,fopt,output.iterations] = FGM(H,c,lb,ub);
    exitflag = 1;
    return;
end

m = size(A,1);

if isempty(b)
    b = zeros(m,1);
end
if isempty(lb_hat)
    lb_hat  = [];
end
if isempty(ub_hat)
    ub_hat = [];
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


% Start the algorithms
if opt.algorithm <= 2
    [zopt,fopt,exitflag,output]...
        = DGM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt);
    
elseif opt.algorithm <= 4
    [zopt,fopt,exitflag,output]...
        = DFGM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt);
    
elseif opt.algorithm <= 6
    [zopt,fopt,exitflag,output]...
        = ALM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt);
    
elseif opt.algorithm <= 8
    [zopt,fopt,exitflag,output]...
        = FALM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt);
    
else
    error('Can not descide the problem case');
end

if (exitflag == -1)
    error('Error trying to run the program');
end

end



