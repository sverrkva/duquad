% INFO: 
% s008: gdm vs fgd
% increasing n
% plot: iterations, convergence


clear all
close all
clc

n=20;
gamma = 0.5;
m=round(3*n / 2);

% OBJECTIVE
H   = 2 * randn(n,n-5); H= H*H' + gamma * eye(n);
c   = randn(n,1)*50;

lb  = ones(n,1)*-2;
ub  = ones(n,1)*2;

z0  = lb;


opt.eps = 0.001;
opt.maxiter = 5000;
opt.alpha = 1;
        
[zopt,fopt,iterations_inner] = FGM(H,c,lb,ub,z0,opt);
fopt
zopt

opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic;
[zopt_quadprog,fopt_quadprog]...
    = quadprog(H,c,[],[],[],[],lb,ub,z0,opts);

fopt_quadprog
% 
% figure()
% hold on;
% legend_array = [];
% cstring = [];
% p=plot(par.n_vec, fgm.iter,'r');
% legend_array = [legend_array p];
% cstring{end+1} = 'FGM';
% 
% p=plot(par.n_vec, gdm.iter,'b');
% legend_array = [legend_array p];
% cstring{end+1} = 'GDM';
% 
% set_legend(legend_array,cstring);
% 
% x='$n$';
% y='iterations';
% tit='Number of iterations';
% set_labels(x,y,tit);
% %%
% 
% opt.eps = 0.00;
% opt.maxiter = 2000;
% opt.alpha = 1;
% 
% par.n = 200;
% par.m = round(par.n/2.2);% round(2*par.n/3);
% 
% [H,c,A,b,lb_hat,ub_hat,lb,ub,z0] = generate_problem_case1(par);
% problem.H = H;
% problem.c = c;
% problem.A = A;
% problem.b = b;
% problem.lb_hat = lb_hat;
% problem.ub_hat = ub_hat;
% problem.lb = lb;
% problem.ub = ub;
% problem.z0 = z0;
% 
% [zopt,fopt,iterations_inner,fgm_fdiff] = FGM(H,c,lb,ub,z0,opt);    
% [zopt,fopt,iterations_inner,gfm_fdiff] = GDM(H,c,lb,ub,z0,opt);
% 
% figure()
% hold on;
% grid on;
% legend_array = [];
% cstring = [];
% p=plot(fgm_fdiff,'r');
% legend_array = [legend_array p];
% cstring{end+1} = 'FGM';
% 
% p=plot(gfm_fdiff,'b');
% legend_array = [legend_array p];
% cstring{end+1} = 'GDM';
% 
% set_legend(legend_array,cstring);
% 
% x='$n$';
% y='iterations';
% tit='Number of iterations';
% set_labels(x,y,tit);
% ylim([-0.5 2])
