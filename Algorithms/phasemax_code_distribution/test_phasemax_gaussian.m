%  Build a random phase retrieval problem and solve it using phasemax.
%  Then, test whether exact reconstruction occured.
%  This script uses a random Guassian measurement matrix, which is
%  explicitly formed and handed to the PhaseMax solver as an argument.  For
%  an example using function handles and fast Fourier transforms, see the
%  script test_phasemax_fourier.
%
%  Copyright Goldstein & Studer, 2016.  For more details, visit 
%  https://www.cs.umd.edu/~tomg/projects/phasemax/


rng(1)             % seed the random number generator
addpath('solvers')

n=100;             % number of unknowns
m = 8*n;           % number of measurements
isComplex = true;  % use complex matrices? or just stick to real?
fprintf('Testing PhaseMax with %d unknowns and %d Gaussian measurements\n', n, m);

%  Build a random test problem
xt = randn(n,1)+isComplex*randn(n,1)*1i; % the true solution vector
A = randn(m,n)+isComplex*randn(m,n)*1i;  % measurment matrix
b = abs(A*xt); % the unsigned phase measurments

%% Try to recover x0
%  generate an initial guess using a spectral method
x0 = init_wirt( A, b, true, true ); 
%  Solve phasemax using a gradient based method
[x,outs] = solve_phasemax_grad(x0, A, b);

fprintf('totalIters = %d\n',outs.totalIterations);

% Determine the optimal phase rotation so that the recovered solution
% matches the true solution as well as possible.  
alpha = (x'*xt)/(x'*x);
x = alpha*x;

% Determine the relative reconstruction error.  If the true signal was 
% recovered, the error should be very small - on the order of the numerical
% accuracy of the solver.
recon_error = norm(xt-x)/norm(xt);
fprintf('relative recon error = %d\n',recon_error);

