%  Build a phase retrieval problem using Fourier measurements, and solve it
%  using phasemax.  Then, test whether exact reconstruction occured.
%  Note:  this script uses function handles to specify the measurement
%  operator, as opposed to dense matrices.
%
%  Copyright Goldstein & Studer, 2016.  For more details, visit 
%  https://www.cs.umd.edu/~tomg/projects/phasemax/

rng(1)     % seed the random number generator
n=10000;   % number of unknowns
m = 8*n;   % number of measurements
fprintf('Testing PhaseMax with %d unknowns and %d Fourier measurements\n', n, m);

%  The true solution
xt = randn(n,1)+randn(n,1)*1i; % the true solution vector

%% Define the Fourier measurements operators
%  For mapping the entries of a n-vector into random spots in an m-vector
takers = randperm(m,n);
expand = sparse(takers,1:n,1,m,n);
A = @(x) fft(expand*x);
At = @(x) expand'*ifft(x)*m;

%% Compute the phaseless measurements
b = abs(A(xt)); 

%% Try to recover x0
%  generate an initial guess using a spectral method
x0 = init_wirt( A, b, true, true, n, At ); 
%  Solve phasemax using a gradient based method
[x,outs] = solve_phasemax_grad(x0, A, b, At);

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

