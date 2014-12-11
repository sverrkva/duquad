function [zopt,fopt,niter] = FGM(varargin)
% *** Fast Gradient method ***
%
% Solves the problem: f(x) = 0.5*x'*H*x + c'x,
% where x = [x1 x2 ... xN]^T
%
% INPUTS:
% H: input matrix (must be positive definite and symmetric)
% c: input vector
% lb: lower bound vector
% ub: upper bound vector
% z0: initial guess
% opt: struct containing algorithm options
% 
% OUTPUTS:
% xopt: optimal solution
% fopt: optimal value (f(xopt))
% niter: number of iteration for finding xopt
% 
% USAGE:
% [xopt,fopt,niter] = FG(H,c);
% [xopt,fopt,niter] = FG(H,c,lb,ub);
% [xopt,fopt,niter] = FG(H,c,lb,ub,z0);
% [xopt,fopt,niter] = FG(H,c,lb,ub,z0,opt);
% [xopt,fopt,niter] = FG(H,c,[],[],[],opt); 
% 
% AUTHOR:
% Sverre Kvamme

% TODO:
% - check if matrix dim er korrekte
% - sholud express if one of the stoppingcriterias are met
% - check for typos
% - make help text perfect
% - take the time
% - update the gradient_descent with the same changes

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
%fprintf('niter inner: %d\n', niter)
zopt = z;
fopt = f_value;

end

function f = obj(z, H, c)
f = 0.5 * (z'*H*z) + c'*z;
end