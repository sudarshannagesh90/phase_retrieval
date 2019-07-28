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

%% sparsePR.m
% 
% Sparse Phase Retrieval made easy - recovering sparse vectors from
% phaseless measurements using a two step process:
%   1. solve a reduced dimension non-sparse phase retrieval problem (we use
%   PhaseLift)
%   2. solve a standard CS recovery problem using the result from step 1.
%   as measurements.
%
% For see details, see arXiv preprint
%   Robust Sparse Phase Retrieval Made Easy
%   Mark Iwen, Aditya Viswanathan, Yang Wang
%   2014
%
% Note that this method requires measurement matrices of the form 
% M = PC, where P \in C^{m x mCS} and C \in C^{mCS x n} are admissible 
% phase retrieval and CS matrices respectively.
%

clear all; close all; clc
addpath ../src
addpath ../third/TFOCS

% For repeatable results
rng(102014);


%% Problem parameters
n = 2^10;                       % Dimension of the problem
s = 005;                        % Sparsity
mCS = round(1.75*s*log(n/s));   % No. of CS measurements
m = round(1.75*mCS*log(mCS));   % Total no. of phaseless measurements

% Noisy data?
addNoise = false;
snr = 40;                       % Desired noise level (SNR, in dB)

% Test function
testFncType = 'sparse';         % Type of test function

% Type of measurements
% We will assume that both the phase retrieval and CS matrices P and C
% respectively are random Gaussian
measurementType = 'randGauss';

% Other options
checkIntSoln = true;           % Check intermediate PhaseLift soln?

% Print out parameters to standard output
fprintf('\nProblem dimension - %d\n', n );
fprintf('Sparsity - %d\n', s );
fprintf('No. of measurements - %d\n', m );
fprintf('CS problem dimension - %d\n', mCS );


%% Generate test signal and measurements
% Test signal
testFnc = generateTestFunction( testFncType, n, s );

% Measurements
[measurements, measurementMat, noise] = generateMeasurements( ...
                        measurementType, testFnc, m, mCS, addNoise, snr );

                    
% Sparse Phase Retrieval Made Easy
%% Step I
%   PhaseLift
%   Solve the Phase Retrieval problem using m measurements, obtaining an 
%   mCS length intermediate vector
%
%   In this script, we will use PhaseLift
%   We solve
%       min      0.5 || b - A(X) ||_2^2 + lambda*trace(X)
%   subject to          X >= 0
fprintf( '\n\n Step 1: Now solving the PhaseLift problem ... \n\n' );

% Use TFOCS to solve the problem
% Set TFOCS parameters
opts.maxIts = 1e3;
opts.tol = 1e-10;
opts.restart = 200;
% opts.largescale = true;

% Regularization parameter
lambda = 5e-2;

% Initial guess
initGuess = zeros(mCS);

% TFOCS requires special initialization of the measurement model
% We use the Phase Retrieval matrix, P
tfocsLinOp = initializeLinopPR( {[mCS, mCS], m}, measurementMat.P ); 

recoveredMat = solver_TraceLS( tfocsLinOp, measurements, ...
                                            lambda, initGuess, opts );

% Intermediate solution
% Choose intermediate solution to be (scaled) leading eigenvector
% Use eigs for large problems ...
[intSoln, eVal] = eig(recoveredMat);
intSoln = sqrt(eVal(end))*intSoln(:,end);


if( checkIntSoln )
    % Check whether intermediate solution is correct
    trueIntSoln = measurementMat.C*testFnc;
    % Correct for global phase factor
    checkSoln = intSoln;        % work on copy of int. soln.
    phaseFac = exp(1i*angle((checkSoln'*trueIntSoln)/(trueIntSoln'*trueIntSoln)));
    checkSoln = checkSoln*phaseFac;
    % Compute error in intermediate solution
    err_2 = norm( trueIntSoln - checkSoln )/ norm(trueIntSoln);
    err_dB = 10*log10( norm(trueIntSoln-checkSoln,2)^2 / norm(trueIntSoln,2)^2 );

    % Print errors to standard output
    fprintf( '\nIntermediate Soln: Error in recovery is %3.3e \n', err_2 );
    fprintf( 'Intermediate Soln: Error in recovery (dB) is %3.3e \n', err_dB );
end

%% Step 2
% CS recovery 
% We will use Basis Pursuit De-Noising (BPDN)
%       min         || x ||_1 
%   subject to     || Ax - b ||_2 <= eps
fprintf( '\n\n Step 2: Now solving the CS problem ... \n\n' );

% De-Noising parameter, eps
if( addNoise )
    epsilon = min(1e-3, 1e-3*noise.noise_power); % Noisy reconstruction
else
    epsilon = 1e-10;                        % No Noise
end

% Use TFOCS to solve the problem
% TFOCS solves the "smoothed" version of the problem parameterized by mu. 
% Typically, mu is small; when mu=0, we solve the exact BPDN problem.
% Set TFOCS parameters
mu = 5e-2;

% TFOCS requires special initialization for complex linear measurements
% This time, use the CS matrix C
tfocsLinOp = linop_matrix( measurementMat.C, 'C2C' );

recoveredSig = solver_sBPDN( tfocsLinOp, intSoln, ...
                                    epsilon, mu, [], [], opts );

                        
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
