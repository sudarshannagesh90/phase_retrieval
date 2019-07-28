function [x_PL] = PhaseLift(y,A,niter)
% function implementing PhaseLift algorithm from Candes et al. 2011. 
% y: Linear phase-less measurements
% A: Gaussian measurement matrix
addpath(genpath('sparsepr/'))
[M,N]  = size(A);
x_init = randn(N,N);
%% Solve the SDP 
lambda = 5e-2;
opts.maxits = niter;
opts.tol    = 1e-10;
opts.restart= 200;
tfocsLinOp = initializeLinopPR({[N,N],M},A); 
recoveredMat = solver_TraceLS( tfocsLinOp, y, lambda, x_init, opts );
%% Compute the leading eigenvector of the matrix X
[recoveredSig,eVal] = eig(recoveredMat);
x_PL = sqrt(eVal(end))*recoveredSig(:,end);
