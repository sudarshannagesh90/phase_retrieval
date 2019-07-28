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

%% sparsePR_robustness.m
% 
% Sparse Phase Retrieval - Robustness study
%   We generate a plot of the noise performance of the 2-step sparse phase
%   retrieval algorithm
%
% For see details, see arXiv preprint
%   Robust SParse Phase Retrieval Made Easy
%   Mark Iwen, Aditya Viswanathan, Yang Wang
%   2014
%
% See sparsePR.m in the demos/ folder for a simple implementation of the
% method.
%

clear all; close all; clc
addpath ../src
addpath ../third/TFOCS

% Use Matlab's parallel computing parfor to speed up the script
% We will do several trials to generate the robustness plot which 
% can be split across process
matlabpool('local', 5);

% For repeatable results
rngSeed = 102014;

%% Problem parameters
n = 2^10;                       % Dimension of the problem
s = 005;                        % Sparsity

% No. of measurements
mCS = round(2.25*s*log(n/s));   % No. of CS measurements
m = round(7.00*mCS);            % Total no. of phaseless measurements

% Robutness plot parameters 
addNoise = true;
snrVals = (20:10:60).';         % Added noise level (SNR, in dB)

% Regularization parameters
% ... contains the reg. parameter lambda for the PhaseLift problem
% each row corresponds to parameters for a different SNR value
regPars = [5e-0; 5e-2; 2.5e-1; 1e-2; 1e-3].';

% ... reg. parameter mu for the CS problem
% TFOCS solves the "smoothed" version of the problem parameterized by mu. 
% Typically, mu is small; when mu=0, we solve the exact BPDN problem.
mu = 1e-2;

% No. of trials
nTrials = 100;                  % No. of trials to generate each plot 
                                % point

% Store 2-norm error (in dB) for plotting
errVals = zeros( nTrials, length(snrVals) );

% Test function
testFncType = 'sparse';         % Type of test function

% Type of measurements
% We will assume that both the phase retrieval and CS matrices are 
% random Gaussian
measurementType = 'randGauss';

% Set TFOCS parameters
opts.maxIts = 1e3;
opts.tol = 1e-10;
opts.restart = 200;
opts.printEvery = 0;
% opts.largescale = true;

% Print out parameters to standard output
fprintf('\nProblem dimension - %d\n', n );
fprintf('Sparsity - %d\n', s );
fprintf('CS Problem dimension - %d\n', mCS );
fprintf('Total no. of measurements - %d\n', m );

fprintf( '\nTrial SNR=20dB SNR=30dB SNR=40dB SNR=50dB SNR=60dB' );
fprintf( '\n----- -------- -------- -------- -------- --------' );

parfor idxTrials = 1:nTrials

% Record errors for this trial in this array
errTrial = zeros(1, length(snrVals) );

%% Generate test signal and measurements
% Test signal
testFnc = generateTestFunction( testFncType, n, s );

% For each value of SNR
for idxSNR = 1:length(snrVals)

% Noise level
snr = snrVals(idxSNR);

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

% Use TFOCS to solve the problem

% Regularization parameter
lambda = regPars(idxSNR);

% Initial guess
initGuess = zeros(mCS);

% TFOCS requires special initialization of the measurement model
tfocsLinOp = initializeLinopPR( {[mCS, mCS], m}, measurementMat.P ); 

recoveredMat = solver_TraceLS( tfocsLinOp, measurements, ...
                                            lambda, initGuess, opts );

% Choose intermediate solution to be (scaled) leading eigenvector
[intSoln, eVal] = eig(recoveredMat);
intSoln = sqrt(eVal(end))*intSoln(:,end);


%% Step 2
% CS recovery 
% We will use Basis Pursuit De-Noising (BPDN)
%       min         || x ||_1 
%   subject to     || Ax - b ||_2 <= eps

% De-Noising parameter, eps
epsilon = min(1e-3, 1e-3*noise.noise_power);

% Use TFOCS to solve the problem

% TFOCS requires special initialization for complex linear measurements
tfocsLinOp = linop_matrix( measurementMat.C, 'C2C' );

recoveredSig = solver_sBPDN( tfocsLinOp, intSoln, ...
                                    epsilon, mu, [], [], opts );

                        
% Correct for global phase factor
phaseFac = exp( 1i* angle( (recoveredSig'*testFnc)/(testFnc'*testFnc) ) );
recoveredSig = recoveredSig*phaseFac;


%% Compute error 
err_dB = 10*log10( norm(testFnc-recoveredSig,2)^2 / norm(testFnc,2)^2 );
errTrial(idxSNR) = err_dB;

end

% End of a trial; write out and print all errors
errVals(idxTrials,:) = errTrial;

fprintf('\n %03d   %03.3f  %03.3f  %03.3f  %03.3f  %03.3f', idxTrials, ...
    errTrial(1), errTrial(2), errTrial(3), ...
                    errTrial(4), errTrial(5) );
                
end

%% Plot SNR vs error curve

% Average SNR from all trials
avgSNR = mean(errVals);

% Prin out averaged SNRs
fprintf('\n\nAverage SNR across all trials is');
fprintf('\n       %3.3f  %3.3f  %3.3f  %3.3f  %3.3f\n', ...
        avgSNR(1), avgSNR(2), avgSNR(3), avgSNR(4), avgSNR(5) );

% Save variables to file
save( 'sparsePR_robustness.mat', 'avgSNR', 'snrVals', 'errVals' );

% Plot
figure('visible', 'off'); plot( snrVals, avgSNR, 'r-+', 'linewidth', 2 ); grid
xlabel( 'SNR (dB)', 'interpreter', 'latex', 'fontsize', 15 );
ylabel( 'Recontruction Error (dB)', 'interpreter', 'latex', 'fontsize', 15 );
title( 'Robustness to Additive Noise: $N=1024$, sparsity $s=5$', ...
                            'interpreter', 'latex', 'fontsize', 15 );
axis([15 65 -65 -15])
print( gcf, '-dpng', 'sparsePR_robustness.png' );
print( gcf, '-depsc', 'sparsePR_robustness.eps' );
print( gcf, '-dpdf', 'sparsePR_robustness.pdf' );
