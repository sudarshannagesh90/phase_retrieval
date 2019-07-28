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

%% PhaseLift.m
% 
% Recover signal from magnitude measurements using PhaseLift
% Uses the TFOCS optimization package
%

clear all; close all; clc
addpath ../src
addpath ../third/TFOCS

% For repeatable results
rng(102014);


%% Problem parameters
n = 2^05;                       % Dimension of the problem
m = round(2.00*n*log(n));       % Total no. of phaseless measurements

% Noisy data?
addNoise = false;
snr = 40;                       % Desired noise level (SNR, in dB)

% Test function
% testFncType = 'sparse';         % Sparse functions
testFncType = 'randgauss';      % Random Gaussian
% testFncType = 'unirand';        % Uniform random
% testFncType = 'sinusoid';       % Sinusoidal

% % No. of components
s = 5;                         % For sparse or sinusoidal signals

% Type of measurements
measurementType = 'randGauss';  % Random Gaussian

% Print out parameters to standard output
fprintf('\nProblem dimension - %d\n', n );
fprintf('No. of measurements - %d\n\n', m );


%% Generate test signal and measurements
% Test signal
if( strcmpi(testFncType, 'sparse') || strcmpi(testFncType, 'sinusoid') )
    testFnc = generateTestFunction( testFncType, n, s );
else
    testFnc = generateTestFunction( testFncType, n );
end

% Measurements
[measurements, measurementMat, noise] = generateCPRLMeasurements( ...
                        measurementType, testFnc, m, addNoise, snr );

                    
%%   PhaseLift
%   Solve the Phase Retrieval problem using m measurements, using PhaseLift
%
%   We solve
%       min      0.5 || b - A(X) ||_2^2 + lambda*trace(X)
%   subject to          X >= 0

% Use TFOCS to solve the problem
% Set TFOCS parameters
opts.maxIts = 1e3;
opts.tol = 1e-10;
opts.restart = 200;
% opts.largescale = true;

% Regularization parameter
lambda = 5e-2;

% Initial guess
initGuess = zeros(n);

% TFOCS requires special initialization of the measurement model
tfocsLinOp = initializeLinopPR( {[n, n], m}, measurementMat ); 

recoveredMat = solver_TraceLS( tfocsLinOp, measurements, ...
                                            lambda, initGuess, opts );

% The above SDP recovers the matrix X = xx*; we need to extract x
% Since X is low-rank (ideally rank-1), choose solution to be (scaled) 
% leading eigenvector                                        
[recoveredSig, eVal] = eig(recoveredMat);
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
