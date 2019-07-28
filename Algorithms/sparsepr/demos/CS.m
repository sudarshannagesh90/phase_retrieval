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

%% CS.m
% 
% Simple implementation of Compressive Sensing reconstruction for the 
% recovery of sparse signals using TFOCS
%

clear all; close all; clc
addpath ../src
addpath ../third/TFOCS

% For repeatable results
rng(102014);


%% Problem parameters
n = 2^10;                   % Dimension of the problem
s = 050;                    % Sparsity
m = round(1.75*s*log(n/s)); % No. of measurements

% Noisy data?
addNoise = false;
snr = 40;                   % Desired noise level (SNR, in dB)

% Test function
testFncType = 'sparse';     % Sparse test function

% Type of measurements
measurementType = 'randGauss';      % Random Gaussian measurements


%% Generate test signal and measurements
% Test signal
testFnc = generateTestFunction( testFncType, n, s );

% Measurements
[measurements, measurementMat, noise] = generateCSMeasurements( ...
                        measurementType, testFnc, m, addNoise, snr );

                    
%% CS recovery 
% We will use Basis Pursuit De-Noising (BPDN)
%       min         || x ||_1 
%   subject to     || Ax - b ||_2 <= eps

% De-Noising parameter, eps
if( addNoise )
    epsilon = 1e-2*noise.noise_power;        % Noisy reconstruction
else
    epsilon = 1e-10;                        % No Noise
end

% Use TFOCS to solve the problem
% TFOCS solves the "smoothed" version of the problem parameterized by mu. 
% Typically, mu is small; when mu=0, we solve the exact BPDN problem.
% Set TFOCS parameters
mu = 2e-1;
opts.maxIts = 1000;
opts.tol = 1e-10;
opts.restart = 100;

% TFOCS requires special initialization for complex linear measurements
tfocsLinOp = linop_matrix( measurementMat, 'C2C' );
[recoveredSig, out, opts] = solver_sBPDN( tfocsLinOp, measurements, ...
                            epsilon, mu, [], [], opts );
                        
                                            
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
