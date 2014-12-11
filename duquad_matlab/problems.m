%% ex 1 - SIMPLE 1
clear all; clc

% OBJECTIVE
H = [11,4 ; 4,22];
c = [3;4];

% CONSTRAINTS
% linear Inequalities
A = [1 1;2 1];
b = [2;3];
lb_hat = []; % can not be inf
ub_hat = [0;0];
% decision Variable Bounds
lb = [-99999;-999999];
ub = [99999;999999];
% initial point
z0 = [30;40];

% QUADPROG
%fprintf('quadprog start\n');
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
    = quadprog(H,c,A,b,[],[],lb,ub,z0,opts);
time_quad = toc;
fprintf('quadprog finish\n');
fprintf('quadprog finished\n');
fprintf('f: %f\n',fopt_quadprog);
fprintf('iter: %d\n\n',OUTPUT.iterations);



%% ex 2 - SIMPLE 2  (constrained)
clear all
clc

% OBJECTIVE
H = [30,0,0 ; 0,6,0 ; 0,0,20];
c = [4;3;20];
% CONSTRAINTS
% linear Inequalities
A = [1,4,5;2,1,3;7,8,9];
b = [2;3;4];
%lb_hat = [-5;-6;-7;-7]; % can not be inf
lb_hat = [-10;-10;-10];%[-1;-2;-3];
ub_hat = [51;20;30];
% decision Variable Bounds
%lb = [-20;-20;-20];
%ub = [20;20;20];
lb = [-99;-99;-99]; %[0;0;0];
ub = [99;99;99];
% initial point
z0 = [30;40;20];


%load('problem','H','c','A','b','lb_hat','ub_hat','lb','ub','z0');

% QUADPROG
%fprintf('quadprog start\n');
% opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
% tic;
% [zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
%     = quadprog(H,c,A,b+ub_hat,[],[],lb,ub,z0,opts);
%     %= quadprog(H,c,[A;-A],[b+ub_hat;-b-lb_hat],[],[],lb,ub,z0,opts);
% time_quad = toc;
% fprintf('quadprog finish\n');
% fprintf('quadprog finished\n');
% fprintf('f: %f\n',fopt_quadprog);
% fprintf('iter: %d\n\n',OUTPUT.iterations);
% 
% addpath('andrei');
% 
% tol_outer_Pf = 0.01;
% tol_inner=0.00001;
% bb = -b;
% 
% [ps_ist pf_ist k]...
%    =dual_fastgradient_lastiter(H,c,0,[0;0;0],[0;0;0],A,bb,0,tol_outer_Pf,tol_inner,10000,z0);
% 
% test_kvamme

fprintf('quadprog start\n');
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
    = quadprog(H,c,[A;-A],[b+ub_hat;-b-lb_hat],[],[],lb,ub,z0,opts);
time_quad = toc;
fprintf('quadprog finished\n');
fprintf('f: %f\n',fopt_quadprog);
fprintf('iter: %d\n\n',OUTPUT.iterations);


%% ex 3 - RANDOM PROBLEM 
clear all
clc

n=10;
gamma = 1; % much better performance if this number is high
m=round(n/2); 

% OBJECTIVE
H = randn(n,n-5); H= H*H' + gamma * eye(n);
c = 5 * randn(n,1);
% CONSTRAINTS
% linear Inequalities
A = randn(m,n); 
b = randn(m,1);
%lb_hat = [];
lb_hat = [];
ub_hat = zeros(m,1);
% decision Variable Bounds
lb = [];
ub = [];
% initial point
z0 = [];




%% ex 4 - RANDOM PROBLEM (case 1: constrained) 
clear all;  clc


n=20;
gamma = 5;
m=round(3*n / 2);

% OBJECTIVE
H   = 2 * randn(n,n-5); H= H*H' + gamma * eye(n);
c   = randn(n,1)*5;
% CONSTRAINTS
% linear Inequalities
A       = randn(m,n); 
b       = randn(m,1);    % This parameter has a lot to say for finding the answer
lb_hat  = -b - 10*rand(m,1);  %rand(m,1)*-10;
ub_hat  = - b + 10*rand(m,1);  %rand(m,1)*10;
% decision Variable Bounds
lb  = ones(n,1)*-2;
ub  = ones(n,1)*2;
% initial point
z0  = randn(n,1)*-30000;

%load('problem','H','c','A','b','lb_hat','ub_hat','lb','ub','z0');

fprintf('Problem generated\n');

fprintf('quadprog start\n');
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
    = quadprog(H,c,[A;-A],[b+ub_hat;-b-lb_hat],[],[],lb,ub,z0,opts);
time_quad = toc;
fprintf('quadprog finished\n');
fprintf('f: %f\n',fopt_quadprog);
fprintf('iter: %d\n\n',OUTPUT.iterations);

% GUROBI

% clear model;
% model.Q = sparse(H*0.5);
% model.obj = c;
% model.A = sparse([A;-A]);
% model.rhs = [b+ub_hat;-b-lb_hat];
% model.sense = '<';
% model.lb = lb;
% model.ub = ub;
% model.start = z0;
% clear params;
% params.outputflag = 1;
% params.BarIterLimit =10000;
% res = gurobi(model,params)



%% ex 5 - RANDOM PROBLEM (case 2: ub_hat = lb_hat)  
clear all;clc
%make

n=1000;
gamma = 2;
m=round(n/2);
%m=round(3*n / 2);

% OBJECTIVE
H   = 2 * randn(n,n-5); H= H*H' + gamma * eye(n);
c   = randn(n,1)*5;
% CONSTRAINTS
% linear Inequalities
A       = randn(m,n); 
b       = randn(m,1);    % This parameter has a lot to say for finding the answer
lb_hat  = [];% -rand(m,1); %-b + 10*rand(m,1);  %rand(m,1)*-10;
ub_hat  = [];%lb_hat;  %rand(m,1)*10;
% decision Variable Bounds
lb  = ones(n,1)*-2;
ub  = ones(n,1)*2;
% initial point
z0  =  randn(n,1)*-30000;

%load('problem_equal_2','H','c','A','b','lb_hat','ub_hat','lb','ub','z0');

% QUADPROG
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
    = quadprog(H,c,[],[],A,b,lb,ub,z0,opts);
time_quad = toc;
fprintf('quadprog finished\n');
fprintf('f: %f\n',fopt_quadprog);
fprintf('iter: %d\n\n',OUTPUT.iterations);




%% ex 6 - RANDOM PROBLEM (case 3: lb_hat = -inf) 
clear all
clc

n=30;
gamma = 1;
%m=round(n/2);
m=round(3*n / 2);

% OBJECTIVE
H   = 2 * randn(n,n-5); H= H*H' + gamma * eye(n);
c   = randn(n,1)*5;
% CONSTRAINTS
% linear Inequalities
A       = randn(m,n); 
b       = randn(m,1);    % This parameter has a lot to say for finding the answer
lb_hat  = [];%-b - 10*rand(m,1);%[];
ub_hat  = - b + 10*rand(m,1);  %rand(m,1)*10;
% decision Variable Bounds
lb  = ones(n,1)*-2;
ub  = ones(n,1)*2;
% initial point
z0  = randn(n,1)*-30000;

%load('problem_lb_inf_2','H','c','A','b','lb_hat','ub_hat','lb','ub','z0');

% QUADPROG
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
    = quadprog(H,c,A,b+ub_hat,[],[],lb,ub,z0,opts);
time_quad = toc;
fprintf('quadprog finished\n');
fprintf('f: %f\n',fopt_quadprog);
fprintf('iter: %d\n\n',OUTPUT.iterations);




%% ex 7 - RANDOM PROBLEM (case 4: ub_hat = inf) 
clear all
clc


n=10;
gamma = 1;
%m=round(n/2);
m=round(3*n / 2);

% OBJECTIVE
H   = 2 * randn(n,n-5); H= H*H' + gamma * eye(n);
c   = randn(n,1)*5;
% CONSTRAINTS
% linear Inequalities
A       = randn(m,n); 
b       = randn(m,1);    % This parameter has a lot to say for finding the answer
lb_hat  = -b - 10*rand(m,1);
ub_hat  = [];
% decision Variable Bounds
lb  = ones(n,1)*-2;
ub  = ones(n,1)*2;
% initial point
z0  = randn(n,1)*-30000;

%load('problem_ub_inf_1','H','c','A','b','lb_hat','ub_hat','lb','ub','z0');

% QUADPROG
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
    = quadprog(H,c,-A,-b-lb_hat,[],[],lb,ub,z0,opts);
time_quad = toc;
fprintf('quadprog finished\n');
fprintf('f: %f\n',fopt_quadprog);
fprintf('iter: %d\n\n',OUTPUT.iterations);

%% ex 10 - compare to OPTI

clear all
clc

% OBJECTIVE
H = [1 -1; -1  2];          %Objective Function (min 0.5x'Hx + f'x)
c = -[2 6]';                
% CONSTRAINTS
% linear Inequalities
A = [1,1; -1,2; 2, 1];      %Linear Inequality Constraints (Ax <= b)
b = [2;2;3];    
lb_hat = []; % can not be inf
ub_hat = zeros(3,1);
% decision Variable Bounds
lb = [0;0];  
ub = [];
% initial point
z0 = [30;40];


% Create OPTI Object
Opt = opti('qp',H,c,'ineq',A,b,'lb',lb)

% Solve the QP problem
[x,fval,exitflag,info] = solve(Opt)

% KVAMME
[xopt,fopt,niter] = DGM(H,c,A,b,lb_hat,ub_hat,lb,ub,z0);
xopt
fopt
niter

%% Problem to compare to Andrei's code


clear all
clc

n=10;
gamma = 2;
%m=round(n/2);
m=n;%round(3*n / 2);

% OBJECTIVE
H   = 2 * randn(n,n-5); H= H*H' + gamma * eye(n);
c   = randn(n,1)*5;
% CONSTRAINTS
% linear Inequalities
A       = randn(m,n); 
b       = randn(m,1);    % This parameter has a lot to say for finding the answer
lb_hat  = [];
ub_hat  = zeros(m,1); 
% decision Variable Bounds
lb  = zeros(n,1);
ub  = [];
% initial point
z0  = randn(n,1)*-30000;

%load('problem_lb_inf_2','H','c','A','b','lb_hat','ub_hat','lb','ub','z0');

% QUADPROG
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
    = quadprog(H,c,A,b+ub_hat,[],[],lb,ub,z0,opts);
time_quad = toc;
fprintf('quadprog finished\n');

maxiter_outer =100000;
maxiter_inner=10000;
tol_outer_ds=0.00001;
tol_inner=0.000001;
tol_outer_Pf = 0.01;

addpath('andrei');

[ps_ist pf_ist k]...
   =dual_fastgradient_lastiter(H,c,0,zeros(n,1),zeros(n,1),A,-b,0,tol_outer_Pf,tol_inner,100000,z0,tol_outer_ds);



%% Problem for FGM

clear all
clc

% OBJECTIVE
H = [1.2204,1.1123,-3.8935 ; 1.1123,3.5821,-3.3333 ; -3.8935,-3.3333,19.6174];
c = [4;3;-2];
% CONSTRAINTS
% decision Variable Bounds
lb = [-1;-1;-1];
ub = [1;1;1];
% initial point
z0 = [10;30;-10];


% QUADPROG
%fprintf('quadprog start\n');
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
    = quadprog(H,c,[],[],[],[],lb,ub,z0,opts);
time_quad = toc;
fprintf('quadprog finish\n');
fprintf('quadprog finished\n');
fprintf('f: %f\n',fopt_quadprog);
fprintf('iter: %d\n\n',OUTPUT.iterations);




% *** FGM ***
opt_FG = set_options_FG(...
'maxiter',10000,...
'eps',0.00001);

tic;
[zopt_FGM,fopt_FGM,niter_FGM]...
    = FGM(H,c,lb,ub,z0,opt_FG);
time_last = toc;
fprintf('DGM last finished\n');
fprintf('f: %f\n',fopt_FGM);
fprintf('iter: %d\n\n',niter_FGM);
zopt_FGM

%% GUROBI

% OBJECTIVE
H = [30,0,0 ; 0,6,0 ; 0,0,20];
c = [4;3;20];
% CONSTRAINTS
% linear Inequalities
A = [1,4,5;2,1,3;7,8,9];
b = [2;3;4];
%lb_hat = [-5;-6;-7;-7]; % can not be inf
lb_hat = [];%[-1;-2;-3];
ub_hat = [0; 0; 0];
% decision Variable Bounds
%lb = [-20;-20;-20];
%ub = [20;20;20];
lb = [0;-99;-99]; %[0;0;0];
ub = [99;99;99];
% initial point
z0 = [30;40;20];

clear model;
model.Q = sparse(H*0.5);
model.obj = c;
model.A = sparse(A);
model.rhs = b;
model.sense = '<';
model.lb = lb;
model.ub = ub;

res = gurobi(model)
res.x 
res.objval
fprintf('quadprog start\n');
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
    = quadprog(H,c,A,b,[],[],lb,ub,z0,opts);
time_quad = toc;
fprintf('quadprog finished\n');
fprintf('f: %f\n',fopt_quadprog);
fprintf('iter: %d\n\n',OUTPUT.iterations);

zopt_quadprog


%% Test DALM: Az = b

clear all;clc
%make

n=40;
gamma = 2;
m=round(n/2);
%m=round(3*n / 2);

% OBJECTIVE
H   = 2 * randn(n,n-5); H= H*H' + gamma * eye(n);
c   = randn(n,1)*5;
% CONSTRAINTS
% linear Inequalities
A       = randn(m,n); 
b       = randn(m,1); 
lb_hat  = -rand(m,1); %-b + 10*rand(m,1);  %rand(m,1)*-10;
ub_hat  = lb_hat;  %rand(m,1)*10;% decision Variable Bounds
lb  = ones(n,1)*-2;
ub  = ones(n,1)*2;
% initial point
z0  = randn(n,1)*-30000;

%load('problem_equal_2','H','c','A','b','lb_hat','ub_hat','lb','ub','z0');

% QUADPROG
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog,EXITFLAG,OUTPUT,LAMBDA]...
    = quadprog(H,c,[],[],A,b,lb,ub,z0,opts);
time_quad = toc;
fprintf('quadprog finished\n');
fprintf('f: %f\n',fopt_quadprog);
fprintf('iter: %d\n\n',OUTPUT.iterations);





