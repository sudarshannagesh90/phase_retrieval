function [x_PM] = PhaseMax(y,A,niters)
y = sqrt(y);

x0 = init_wirt( A, y, true, true ); 
%  Solve phasemax using a gradient based method
[x_PM,outs] = solve_phasemax_grad(x0, A, y,niters);
