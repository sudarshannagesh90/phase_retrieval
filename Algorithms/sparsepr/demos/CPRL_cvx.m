% Copyright (c) 2014 Michigan State University and the CHARMS research
% group
% 
% This file is part of the SparsePR software package
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the "Software"), to deal in
% the Software without restriction, including without limitation the rights to
% use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
% of the Software, and to permit persons to whom the Software is furnished to do
% so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%% CPRL.m
% 
% Simple implementation of Compressive Phase Retrieval via Lifting
% Uses the CVX optimization software package
%
% For see details, see
%   Compressive Phase Retrieval From Squared Output Measurements Via 
%   Semidefinite Programming
%   Ohlsson, Henrik; Yang, Allen Y.; Dong, Roy; Shankar Sastry, S.
%   eprint arXiv:1111.6323
%   2011
%

clear all; close all; clc
addpath ../src

% For repeatable results
rng(102014);


%% Problem parameters
n = 2^05;                           % Dimension of the problem
s = 002;                            % Sparsity
m = 30;                             % Total no. of phaseless measurements

% Noisy data?
addNoise = false;
snr = 40;                       % Desired noise level (SNR, in dB)

% Test function
testFncType = 'sparse';         % Type of test function

% Type of measurements
measurementType = 'randGauss';  % random Gaussian

% Print out parameters to standard output
fprintf('\nProblem dimension - %d\n', n );
fprintf('Sparsity - %d\n', s );
fprintf('No. of measurements - %d\n', m );


%% Generate test signal and measurements
% Test signal
testFnc = generateTestFunction( testFncType, n, s );

% Measurements
[measurements, measurementMat, noise] = generateCPRLMeasurements( ...
                        measurementType, testFnc, m, addNoise, snr );

                    
%% CPRL
%   Solve the Sparse Phase Retrieval problem using
%       min      || b - A(X) ||_1 + mu*|| X ||_1
%   subject to          X >= 0
% Note: We are using a variant of PhaseLift for noisy data
% Standard PhaseLift problem formulation is
%       min             Tr(X)
%   subject to         A(X) = b
%                       X >= 0
% 
% For noisy data, [1] recommends solving
%       min         || b - A(X) ||_1
%   subject to          X >= 0
% The CPRL formulation we solve is based on this implementation of
% PhaseLift
%
% [1]   Solving Quadratic Equations via PhaseLift when There Are About As 
%       Many Equations As Unknowns
%       Emmanuel J. Candes and Xiaodong Li
%       http://arxiv.org/abs/1208.6247v2
%


% Regularization parameter
% Change this with added noise level
mu = 5e-2;

% Measurements matrix
A = measurementMat;

% Use CVX to solve the problem
cvx_begin sdp
    % If you do not have access to the MOSEK solver (MOSEK ApS or CVX 
    % Professional license), comment this line
    cvx_solver mosek
    variable X(n,n) hermitian
    
    minimize norm( diag(A*X*A') - measurements, 1 ) + mu*norm( X(:), 1 )
    subject to
        X >= 0;
cvx_end

% The above SDP recovers the matrix X = xx*; we need to extract x
% Since X is low-rank (ideally rank-1), choose solution to be (scaled) 
% leading eigenvector
[recoveredSig, eVal] = eig(X);
recoveredSig = sqrt(eVal(end))*recoveredSig(:,end);
                        
% Correct for global phase factor
phaseFac = exp( 1i* angle( (recoveredSig'*testFnc)/(testFnc'*testFnc) ) );
recoveredSig = recoveredSig*phaseFac;


%% Plot results
% Compute error 
err_2 = norm( testFnc - recoveredSig )/ norm(testFnc);
err_dB = 10*log10( norm(testFnc-recoveredSig,2)^2 / norm(testFnc,2)^2 );

% Print errors to standard output
fprintf( '\nError in recovery (2-norm) is %3.3e \n', err_2 );
fprintf( 'Error in recovery (dB) is %3.3e \n', err_dB );

% Plot reconstructions
% ... real part
figure(1); subplot(1,2,1); 
stem( real(testFnc), 'r-+', 'linewidth', 2 ); hold on
stem( real(recoveredSig), 'b-o', 'linewidth', 2 );
xlabel 'k'; ylabel 'Re(x[k])'; grid; title 'Real Part'
legend( 'True', 'Recovered' )

% ... imaginary part
subplot(1,2,2); 
stem( imag(testFnc), 'r-+', 'linewidth', 2 ); hold on
stem( imag(recoveredSig), 'b-o', 'linewidth', 2 );
xlabel 'k'; ylabel 'Im(x[k])'; grid; title 'Imaginary Part'
legend( 'True', 'Recovered' )
