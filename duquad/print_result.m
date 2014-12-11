function print_result(zopt,fopt,exitflag,output,time) 

fprintf('f: %f\n',fopt);
fprintf('exitflag: %d\n',exitflag);
fprintf('iterations: %d\n',output.iterations);
fprintf('iterations inner tot: %d\n',output.iterations_inner_tot);
fprintf('niter_feasible_ds: %d\n',output.niter_feasible_ds);
fprintf('niter_feasible_pf: %d\n',output.niter_feasible_pf);
fprintf('runtime: %f\n',time);

end

