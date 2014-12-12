DuQuad
======

The DuQuad optimization toolbox attempts to solve the quadratic programming problem on the form:

min 0.5\*z'\*H\*z + c'\*z,  
&nbsp;&nbsp;z  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;s.t.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lb_hat <= Az - b <= ub_hat   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;lb <= z <= ub  

The program is written in C but utilizes a Matlab interface. 

### duquad_c/
This folder contains  the c-code for the algorithms, a makefile (make.m) for compiling a mex file, and a matlab-file (duquad.m) which run the program. The program has only been compiled in linux. An example file is also included to get started. 

### duquad_matlab/
This folder contains the duquad program written in matlab. An example file is included to get started

### Info
Project webpage and more information:
http://sverrkva.github.io/duquad/

