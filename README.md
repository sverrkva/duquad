DuQuad
======


*** DuQuad - Quadratic Programming Optimization ***

 Attempts to solve the quadratic programming problem:


   min f(z) = 0.5*z'*H*z + c'z,
    z
       s.t.
       lb_hat <= Az - b <= ub
       lb <= z <= ub

where z = [z1 z2 ... zN]^T

INPUTS:
H:        Hessian matrix (must be positive definite and symmetric)
c:        gradient vector
A:        linear constraints matrix
b:        linear constraints vector
lb_hat:   lower bound for the linear constraints
ub_hat:   upper bound for the linear constraints
lb:       lower bound for optimization variable z
ub:       upper bound for optimization variable z
z0:       initial point
opt:      struct containing options, see OPTIONS

OUTPUTS:
zopt:     optimal solution
fopt:     optimal value (f(xopt))
niter:    number of iteration for finding the solution
output:   struct containing info from the problem solving

USAGE:
[zopt,fopt,niter] = duquad(H,c);
[zopt,fopt,niter] = duquad(H,c,A,b);
[zopt,fopt,niter] = duquad(H,c,A,b,lb_hat,ub_hat,lb,ub,z0,opt);
[zopt,fopt,niter] = duquad(H,c,A,b,[],[],[],[],[],opt);

OPTIONS:
opt.maxiter_outer:  Maximum number of iterations in the outer loop
opt.maxiter_inner:  Maximum number of iterations in the inner loop
opt.eps_ds:         Tolerance for dual suboptimality
opt.eps_pf:         Tolerance for primal feasibility
opt.eps_inner:      Tolerance for primal feasibility in the inner problem
opt.rho:            Penalty parameter used in ALM and FALM
opt.algorithm:      Spesifies the algoritm used to solve the problem. Values: 
                    1: DGM last
                    2: DGM avg
                    3: DFGM last
                    4: DFGM avg
                    5: ALM last
                    6: ALM avg
                    7: FALM last
                    8: FALM avg

AUTHOR:
Sverre Kvamme

Project webpage and more information:
http://sverrkva.github.io/duquad/

