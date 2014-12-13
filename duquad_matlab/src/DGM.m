function [zopt,fopt,exitflag,output,lambda1,lambda2] = DGM(varargin)

[H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt] = deal(varargin{1:10});

n = size(H,2);
m = size(A,1);

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

% initial point must be feasible
if ~ub_inf
    z0 = min(ub,z0);
end
if ~lb_inf
    z0 = max(lb,z0);
end

% Desciding the category/case of the problem 
if isempty(ub_hat) && isempty(lb_hat)
    error('both lb_hat and ub_hat is unconstrained in DGM()')
elseif isempty(ub_hat)
    problem_case = 4;
elseif isempty(lb_hat)
    problem_case = 3;
elseif lb_hat == ub_hat
    problem_case = 2;
else
    problem_case = 1;
end

% calculating alpha
H_eig = eig(H);
sigma = min(H_eig);
switch problem_case
    case 1   
        Ld = norm([A;-A],2)^2 / sigma;
    case 2 
        Ld = norm(A,2)^2 / sigma;
    case 3  
        Ld = norm(A,2)^2 / sigma;
    case 4  
        Ld = norm(A,2)^2 / sigma;
end
alpha = 1/(2*Ld);

% Setting options and info for the inner problem
optFG.maxiter = opt.maxiter_inner;
optFG.eps = opt.eps_inner;
info.L = max(H_eig);
info.mu = sigma;

% initialize lambda
lambda1 = zeros(m,1);
lambda2 = zeros(m,1);

% making A', b+lb_hat and b+ub_hat, because this is an constant operations
A_t = A';
switch problem_case
    case 1   
        b_ub_hat = b + ub_hat;
        b_lb_hat = b + lb_hat;
    case 2 
        b_ub_hat = b + ub_hat;
        b_lb_hat = [];
    case 3  
        b_ub_hat = b + ub_hat;
        b_lb_hat = [];
    case 4  
        b_ub_hat = [];
        b_lb_hat = b + lb_hat;
end


% solve inner problem first time
c_hat = c;  % lambda = 0;
[z,~, iterations_inner] = FGM(H,c_hat,lb,ub,z0,optFG,info);
iterations_inner_tot = iterations_inner;
A_z = A*z; % used later

% find the value of the dual function (Lagrangian) first time
dual_value = dual_obj(z,H,c,lambda1,lambda2,A_z,b_ub_hat,b_lb_hat,problem_case);
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
while or(dual_value_diff >= opt.eps_ds , pf > opt.eps_pf) 

    % check if the maximum number of iterations is reached
    if niter >= opt.maxiter_outer
        warning('reached maximum number of iterations in DGM\n')
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
    switch problem_case
        case 1   
            % calculating next lambda  (lambda_{lambda+1})
            lambda1 = lambda1 + alpha * (A_z - b_ub_hat);
            lambda2 = lambda2 + alpha * (b_lb_hat - A_z);

            % project lambda_new on R+
            lambda1 = max(0, lambda1);
            lambda2 = max(0, lambda2);
        case 2 
            lambda1 = lambda1 + alpha * (A_z - b_ub_hat);
            % NOTE: No projection
        case 3  
            lambda1 = lambda1 + alpha * (A_z - b_ub_hat);
            lambda1 = max(0, lambda1);
        case 4  
            lambda2 = lambda2 + alpha * (b_lb_hat - A_z);
            lambda2 = max(0, lambda2);
    end

    % *********************************************************************
	%			Solving the inner problem
	% *********************************************************************
    switch problem_case
        case 1   
            c_hat = c + A_t*(lambda1 - lambda2);
        case 2 
            c_hat = c + A_t*lambda1;
        case 3  
            c_hat = c + A_t*lambda1;
        case 4  
            c_hat = c - A_t*lambda2;
    end
    
    [z,~,iterations_inner] = FGM(H,c_hat,lb,ub,z,optFG,info);
    iterations_inner_tot = iterations_inner_tot + iterations_inner;
    A_z = A*z; % only do this operation here
    
    % *********************************************************************
	%			Finding dual_value_diff
	% *********************************************************************
    dual_value_new = dual_obj(z,H,c,lambda1,lambda2,A_z,b_ub_hat,b_lb_hat,problem_case);
    dual_value_diff = abs(dual_value_new - dual_value);
    dual_value = dual_value_new;

    % *********************************************************************
	%			Finding pf
	% *********************************************************************
    % *** LAST z in stopping criteria ***
    if opt.algorithm == 1
        switch problem_case
            case 1   
                pf_vec = [A_z - b_ub_hat ; b_lb_hat - A_z];
                pf_vec = max(0, pf_vec);
            case 2 
                pf_vec = A_z - b_ub_hat;
                % NOTE: No projection
            case 3  
                pf_vec = A_z - b_ub_hat;
                pf_vec = max(0, pf_vec);
            case 4  
                pf_vec = b_lb_hat - A_z;
                pf_vec = max(0, pf_vec);
        end

    % *** AVERAGE z in pf stopping criteria ***
    elseif opt.algorithm == 2
        summ = summ + z;
        z_avg = summ /(niter+1);
        A_z_avg = A*z_avg;
        switch problem_case
            case 1   
                pf_vec = [A_z_avg - b_ub_hat ; b_lb_hat - A_z_avg];
                pf_vec = max(0, pf_vec);
            case 2 
                pf_vec = A_z_avg - b_ub_hat;
                % NOTE: No projection
            case 3  
                pf_vec = A_z_avg - b_ub_hat;
                pf_vec = max(0, pf_vec);
            case 4  
                pf_vec = b_lb_hat - A_z_avg;
                pf_vec = max(0, pf_vec);
        end
        
    else
        error('Choose opt.algorithm = 1 / 2 when running DGM')
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


if opt.algorithm == 1
    zopt = z;
elseif opt.algorithm == 2
    zopt = z_avg;
else
    error('Choose opt.algorithm = 1 / 2 when running DGM')
end

fopt = obj(zopt,H,c);
output.iterations = niter-1;
output.iterations_inner_tot = iterations_inner_tot;
output.niter_feasible_ds = niter_feasible_ds+1;
output.niter_feasible_pf = niter_feasible_pf+1;


if opt.algorithm == 1
    output.algorithm = 'DGM last';
else 
    output.algorithm = 'DGM avg';
end
   
end



% *************************************************************************
%			Local functions
% *************************************************************************
function value = obj(z, H, c) 
    value = 0.5* (z'*H*z) + c'*z;
end

function value = dual_obj(z,H,c,lambda1,lambda2,A_z,b_ub_hat,b_lb_hat,problem_case)
  
    switch problem_case
        case 1  
            value = 0.5* (z'*H*z) + c'*z... 
               + lambda1'*(A_z-b_ub_hat)...
               + lambda2'*(b_lb_hat-A_z);
       case 2 
            value = 0.5* (z'*H*z) + c'*z... 
               + lambda1'*(A_z-b_ub_hat);
        case 3  

            value = 0.5* (z'*H*z) + c'*z... 
               + lambda1'*(A_z-b_ub_hat);
        case 4  
            value = 0.5* (z'*H*z) + c'*z... 
               + lambda2'*(b_lb_hat-A_z);
    end
end