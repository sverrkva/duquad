DuQuad
======

### Quadratic Programming Optimization 

DuQuad attempts to solve the quadratic programming problem:

min 0.5\*z'\*H\*z + c'\*z,  
&nbsp;z  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;s.t.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lb_hat <= Az - b <= ub_hat   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lb <= z <= ub  

### duquad/
This folder contains the program and an example file. The program is written in C but utilizes a Matlab interface. The mex file is compiled of linux.

### duquad_matlab/
This folder contains...

AUTHOR:
Sverre Kvamme

Project webpage and more information:
http://sverrkva.github.io/duquad/

