function [zopt,fopt,niter] = GDM(varargin)
% *** Descend Gradient method ***
%
% Solves the problem: f(z) = 0.5*z'*H*z + c'z,
% where z = [z1 z2 ... zN]^T
%
% INPUTS:
% H: input matrix (must be positive definite and symmetric)
% c: input vector
% lb: lower bound vector
% ub: upper bound vector
% z0: initial guess
% 
% OUTPUTS:
% zopt: optimal solution
% fopt: optimal value (f(zopt))
% niter: number of iteration for finding zopt
% 
% USAGE:
% [zopt,fopt,niter] = fast_gradient(H,c,lb,ub);
% [zopt,fopt,niter] = fast_gradient(H,c,lb,ub,z0);
% [zopt,fopt,niter] = fast_gradient(H,c,[],[],z0); - unconstrained case
% 
% AUTHOR:
% Sverre Kvamme

% TODO: 
% - Kommentere og gj�re finere
% - ta innput med maxiter osv
% - make default for lb and ub
% - sholud express if one of the stoppingcriterias are met
% - Do so that you sont have to push lb and ub
% - check for typos
% - make help text perfect

if nargin < 7
    info = [];
    if nargin < 6
        opt = [];
        if nargin < 5
            z0 = [];
            if nargin < 4
                ub = [];
                if nargin < 3
                    lb = [];
                    if nargin < 2
                        error('Incorrect number of input arguments.')                 
                    end
                end
            end
        end
    end
end

[H,c] = deal(varargin{1:2});

if nargin > 2
    lb = varargin{3};
end    
if nargin > 3
    ub = varargin{4};
end    
if nargin > 4
    z0 = varargin{5};
end  
if nargin > 5
    opt = varargin{6};
end
if nargin > 6
    info = varargin{7};
end
if nargin > 7
    error('Incorrect number of input arguments.')
end   

n = size(H,2);

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
    opt = set_options('default');
end
if isempty(info)
    H_eig = eig(H);
    info.L =  max(H_eig);
    info.mu = min(H_eig);
end

alpha = 1/(info.L);

% initialize other necessary variables and vectors
diff_function_value = opt.eps+1; 
f_value = obj(z0,H,c);
z = z0;

niter = 0;

while diff_function_value > opt.eps
    
    if niter > opt.maxiter
        %warning('reached maximum number of iterations in fast_gradient'
        break;
    end
    
    % calculating next z  (z_{k+1})
    znew = z - alpha * grad(z,H,c);
   
    % projection znew on the feasible set made of lb and ub,
    % but do not project if problem is unconstrained 
    if ~ub_inf
        znew = min(ub,znew);
    end
    if ~lb_inf
        znew = max(lb,znew);
    end
    
    % calculate the difference between the new and previous function value
    f_value_new = obj(znew,H,c);
    diff_function_value = abs(f_value_new - f_value);
    
    % update the variables
    z = znew;
    f_value = f_value_new;
    niter = niter+1;
end

zopt = z;
fopt = obj(zopt,H,c);

end

function f = obj(z, H, c) 
f = 0.5* z'*H*z + c'*z;
end

function g = grad(z, H, c)    
g = H*z + c;   
end
