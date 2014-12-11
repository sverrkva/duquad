function print_problem_to_file(varargin)

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

if any(eig(H)<=0)
    error('Non-convex. (Quadratic form is not positive definite.)')
end
if any(H~=H')
    error('Quadratic form must be symmetric.')
end
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
    opt = set_options_DGM('default');
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

eigH = eig(H);
eigH_min = min(eigH);
eigH_max = max(eigH);

switch problem_case
    case 1  % ub_hat and lb_hat exists
        Ld = norm([A;-A],2)^2 / eigH_min;
    case 2  % ub_hat = lb_hat       
        Ld = norm(A,2)^2 / eigH_min;
    case 3  % lb_hat = -inf
        Ld = norm(A,2)^2 / eigH_min;
    case 4  % ub_hat = inf
        Ld = norm(A,2)^2 / eigH_min;       
end

% print problem to file on this format
% n
% m
% lb_is_inf
% ub_is_inf
% problem_case
% eigH_min
% eigH_max
% Ld
% maxiter_outer
% maxiter_inner
% eps_ds
% eps_pf
% eps_inner
% algorithm
% H
% c
% A
% b
% lb_hat
% ub_hat
% lb
% ub
% z0
% opt

filname = '/home/sverre/Dropbox/NTNU/Master/code/thesis/matlab/problem.txt';
fid = fopen(filname,'w');

fprintf(fid,'%d\n',n);
fprintf(fid,'%d\n',m);
fprintf(fid,'%d\n',lb_is_inf);
fprintf(fid,'%d\n',ub_is_inf);
fprintf(fid,'%d\n',problem_case);
fprintf(fid,'%f\n',eigH_min);
fprintf(fid,'%f\n',eigH_max);
fprintf(fid,'%f\n',Ld);
fprintf(fid,'%d\n',opt.maxiter_outer);
fprintf(fid,'%d\n',opt.maxiter_inner);
fprintf(fid,'%f\n',opt.eps_ds);
fprintf(fid,'%f\n',opt.eps_pf);
fprintf(fid,'%f\n',opt.eps_inner);
fprintf(fid,'%d\n',opt.algorithm);

for i=1:n
    for j=1:n
        fprintf(fid,'%f,',H(i,j));
    end
    fprintf(fid,'\n');
end
for i=1:n
    fprintf(fid,'%f,',c(i));
end
fprintf(fid,'\n');
for i=1:m
    for j=1:n
        fprintf(fid,'%f,',A(i,j));
    end
    fprintf(fid,'\n');
end
for i=1:m
    fprintf(fid,'%f,',b(i));
end
fprintf(fid,'\n');
if (problem_case == 1 || problem_case == 4)
    for i=1:m
        fprintf(fid,'%f,',lb_hat(i));
    end
    fprintf(fid,'\n');
end
if (problem_case ~= 4)
    for i=1:m
        fprintf(fid,'%f,',ub_hat(i));
    end
    fprintf(fid,'\n');
end
for i=1:n
    fprintf(fid,'%f,',lb(i));
end
fprintf(fid,'\n');
for i=1:n
    fprintf(fid,'%f,',ub(i));
end
fprintf(fid,'\n');
for i=1:n
    fprintf(fid,'%f,',z0(i));
end
fprintf(fid,'\n');

fclose(fid);

end