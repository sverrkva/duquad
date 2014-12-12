function [zopt,fopt,niter] = FGM(varargin)
%
% *** Fast Gradient method ***
%
% Solves the problem: 
%
%      min 0.5*z'*H*z + c'*z,
%       z
%          s.t.
%          lb <= z <= ub,
%           
% INPUTS:
% H: Hessian matrix (must be positive definite and symmetric)
% c: gradient vector
% lb: lower bound for optimization variable z
% ub: upper bound for optimization variable z
% z0: initial point
% opt: struct containing algorithm options
% info: contains smallest and largest eigenvalue of H
% 
% OUTPUTS:
% zopt: optimal solution
% fopt: optimal value (f(zopt))
% niter: number of iteration for finding the solution
% 
% USAGE:
% [zopt,fopt,niter] = FGM(H,c);
% [zopt,fopt,niter] = FGM(H,c,lb,ub);
% [zopt,fopt,niter] = FGM(H,c,lb,ub,z0);
% [zopt,fopt,niter] = FGM(H,c,lb,ub,z0,opt);
% [zopt,fopt,niter] = FGM(H,c,[],[],[],opt); 
% 
% AUTHOR:
% Sverre Kvamme

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
    opt.maxiter = 10000;
    opt.eps = 0.001;
end
if isempty(info)
    H_eig = eig(H);
    info.L =  max(H_eig);
    info.mu = min(H_eig);
end
    
q = info.mu / info.L;
one_over_L = 1/info.L;
alpha = 0.5;

% initialize other necessary variables and vectors
diff_function_value = opt.eps+1; 
f_value = obj(z0,H,c);
z = z0;
y = z;

niter = 0;


while diff_function_value > opt.eps
   
    if niter > opt.maxiter
        %warning('reached maximum number of iterations in fast_gradient'
        break;
    end
    
    % calculating next z = y - (1/L)(H*y+c)
    znew = y - one_over_L * (H*y + c);%grad(y,H,c);
     
    % projection xnew on the feasible set made of lb and ub,
    % but do not project if problem is unconstrained 
    if ~ub_inf
        znew = min(ub,znew);
    end
    if ~lb_inf
        znew = max(lb,znew);
    end

    % calculating next y 
    alpha_pow2 = alpha*alpha;
    theta = q - alpha_pow2;
    alpha_new = 0.5 * (theta + sqrt(theta^2 + 4*alpha_pow2));
    beta = (alpha*(1-alpha)) / (alpha_pow2 + alpha_new);
    y = znew + beta*(znew - z);
        
    % calculate the difference between the new and previous function value
    f_value_new = obj(znew,H,c);
    diff_function_value = abs(f_value_new - f_value);

    % update the variables
    z = znew;
    f_value = f_value_new;
    alpha = alpha_new;
    niter = niter+1;
   
end

zopt = z;
fopt = f_value;

end

function f = obj(z, H, c)
f = 0.5 * (z'*H*z) + c'*z;
end