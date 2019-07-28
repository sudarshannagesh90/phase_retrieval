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

%% CPRL_runtime_measurements.m
% 
% Generate plots illustrating number of measurements required and runtime
% for CPRL as a function of sparsity
%
% For details, see
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
n = 2^06;                           % Dimension of the problem

% Noisy data?
addNoise = false;
snr = 40;                       % Desired noise level (SNR, in dB)

% Test function
testFncType = 'sparse';         % Type of test function

% Type of measurements
measurementType = 'randGauss';  % Random Gaussian

% No. of trials for each data point in plot
nTrials = 100;

% Regularization parameter
mu = 5e-2;

% Sparsity values
sVals = (1:5).';

% We will generate plots as a function of sparsity for fixed n
runtime = zeros( nTrials, length(sVals) );
nMeasurements = zeros( length(sVals), 1 );
error = zeros( nTrials, length(sVals) );

% We will find no. of measurements for 95% successful reconstruction. Use
% this array to keep track of unsuccessful reconstructions
failure = zeros( nTrials, length(sVals) );

% Print out problem parameters
fprintf( '\n Problem dimension: n=%d\n', n );

% Sparsity loop ...
for idxSVals = 1:length(sVals)

% Set sparsity
s = sVals(idxSVals);

% Start with no. of measurements equal to ...
% (using plot from reference paper as guide)
m = 10*s;

% Print status message
fprintf( '\nEvaluating measurements reqd. for sparsity %d; trying m=%d\n', ...
                        s, m );
fprintf( '\nTrial    Error (dB)    Success?    Runtime (in secs.)' );
fprintf( '\n-----    ----------    --------    ------------------' );

% Initialize trial counter and failure tracker
idxTrials = 1;
failure(:, idxSVals) = 0;
    
% Trials loop ...
while(idxTrials <= nTrials)

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
% Note: We are using a variant of PhaseLift for noisy data with added l1
% norm constraints

% Timer
timer = tic;

% Use CVX to solve the problem
cvx_begin sdp quiet
    cvx_solver mosek
    variable X(n,n) hermitian
    
    minimize norm(diag(measurementMat*X*measurementMat')-measurements,1) + mu*norm( X(:), 1 )
    subject to
        X >= 0;
cvx_end

% Choose solution to be (scaled) leading eigenvector
[recoveredSig, eVal] = eig(X);
recoveredSig = sqrt(eVal(end))*recoveredSig(:,end);
                        
% Correct for global phase factor
phaseFac = exp( 1i* angle( (recoveredSig'*testFnc)/(testFnc'*testFnc) ) );
recoveredSig = recoveredSig*phaseFac;

%% Record error and execution times
% Record execution time
runtime( idxTrials, idxSVals ) = toc(timer);

% Compute error 
err_dB = 10*log10( norm(testFnc-recoveredSig,2)^2 / norm(testFnc,2)^2 );
error( idxTrials, idxSVals ) = err_dB;

% Check if we have successful reconstruction
if( err_dB > -100 )
    % Unsuccessful reconstruction
    failure( idxTrials, idxSVals ) = 1;
    % Print out status
    fprintf('\n %03d      %+08.3f         %d             %+08.3f', ...
        idxTrials, err_dB, false, runtime(idxTrials,idxSVals) );
    % Update trial counter
    idxTrials = idxTrials+1;    
    
    if( mean( failure(:,idxSVals) ) > 0.05 )
        % If more than 5% of trials are unsuccessful, we need more 
        % measurements
        m = m+1;
        % Reset counters
        idxTrials = 1;
        failure(:,idxSVals) = 0;
        % Print message
        fprintf('\n\n Too many unsuccessful reconstructions.' );
        fprintf(' Need ... more ... measurements!' );
        fprintf('\n Now trying m=%d\n\n', m );
        fprintf( '\nEvaluating measurements reqd. for sparsity %d; trying m=%d\n', ...
                                s, m );
        fprintf( '\nTrial    Error (dB)    Success?    Runtime (in secs.)' );
        fprintf( '\n-----    ----------    --------    ------------------' );       
    end
else
    % Print out status
    fprintf('\n %03d      %+08.3f         %d             %+08.3f', ...
        idxTrials, err_dB, true, runtime(idxTrials,idxSVals) );
    % Update trial counter
    idxTrials = idxTrials+1;    
end


end

% Successfull reconstruction in at least 95% trials
% Store no. of measurements required
nMeasurements( idxSVals) = m;

% Print no. of measurements required
fprintf( '\n Success!' );
fprintf( '\n No. of measurements required at sparsity %d is %d', s, m ); 
fprintf( ' Avg. runtime at sparsity %d is % 3.3f', s, ...
                                            mean(runtime(:,idxSVals)) );
                                
end


% Save results to file
save( 'CPRL_runtime_measurements.mat', 'sVals', 'n', 'runtime', ...
                        'nMeasurements', 'error' );


% Plot runtime and no. of measurements
figure('visible', 'on'); 
plot( sVals, mean(runtime), 'r-+', 'linewidth', 2 );
xlabel( 'Sparsity, $s$', 'interpreter', 'latex', 'fontsize', 15 );
ylabel( 'Avg. runtime (in secs.)', 'interpreter', 'latex', 'fontsize', 15 );
grid; title( 'CPRL Runtime: N=64, Noiseless Measurements', ...
     'interpreter', 'latex', 'fontsize', 15 );
 
% Save plots
print( gcf, '-dpng', 'CPRL_runtime.png' );
print( gcf, '-depsc', 'CPRL_runtime.eps' );
print( gcf, '-dpdf', 'CPRL_runtime.pdf' );
 
figure('visible', 'on'); 
plot( sVals, nMeasurements, 'r-+', 'linewidth', 2 );
xlabel( 'Sparsity, $s$', 'interpreter', 'latex', 'fontsize', 15 );
ylabel( 'No. of measurements', 'interpreter', 'latex', 'fontsize', 15 );
grid; title( 'CPRL Measurements: N=64, Noiseless Measurements', ...
     'interpreter', 'latex', 'fontsize', 15 );

% Save plots
print( gcf, '-dpng', 'CPRL_measurements.png' );
print( gcf, '-depsc', 'CPRL_measurements.eps' );
print( gcf, '-dpdf', 'CPRL_measurements.pdf' ); 
