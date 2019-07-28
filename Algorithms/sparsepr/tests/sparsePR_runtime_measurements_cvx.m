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

%% sparsePR_runtime_measurements.m
% 
% Generate plots showing number of measurements required and runtime
% (as a funtion of sparsity) for the two-step sparse phase retrieval 
% algorithm
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

% We fix the CS problem dimension
% (obtained using experiments on general CS problems)
mCS = max( 10, round(1.75*s.*log(n./s)) );

% Start with 3.75 times the CS dimension for total no. of phaseless 
% measurements
m = round(3.75*mCS);

% Print status message
fprintf( '\nEvaluating measurements reqd. for sparsity %d; trying m=%d\n', ...
                        s, m );
fprintf( '\nTrial    Error (dB)    Int. Soln. Error (dB)    Success?    Runtime (in secs.)' );
fprintf( '\n-----    ----------    ---------------------    --------    ------------------' );

% Initialize trial counter and failure tracker
idxTrials = 1;
failure(:, idxSVals) = 0;
    
% Trials loop ...
while(idxTrials <= nTrials)

%% Generate test signal and measurements
% Test signal
testFnc = generateTestFunction( testFncType, n, s );

% Measurements
[measurements, measurementMat, noise] = generateMeasurements( ...
                        measurementType, testFnc, m, mCS, addNoise, snr );
P = measurementMat.P;       % Pahse retrieval matrix
C = measurementMat.C;       % Compressive sensing matrix
                    
                    
% Sparse Phase Retrieval Made Easy
%% Step I
%   PhaseLift
%   Solve the Phase Retrieval problem using m measurements, obtaining an 
%   mCS length intermediate vector
%
%   In this script, we will use PhaseLift
%   We solve
%       min         || b - A(X) ||_1
%   subject to          X >= 0

% Initialize timer
timer = tic;

cvx_begin sdp quiet
    cvx_solver mosek;
    variable X(mCS,mCS) hermitian
    
    minimize norm( diag(P*X*P') - measurements, 1 )
    subject to
        X >= 0;
cvx_end

% Choose solution to be (scaled) leading eigenvector
[intSoln, eVal] = eig(X);
intSoln = sqrt(eVal(end))*intSoln(:,end);

% Don't count time while printing out the intermediate solution error
step1_time = toc(timer);

% Check intermediate solution error
trueIntSoln = C*testFnc;
% Correct for global phase factor
checkSoln = intSoln;        % work on copy of int. soln.
phaseFac = exp(1i*angle((checkSoln'*trueIntSoln)/(trueIntSoln'*trueIntSoln)));
checkSoln = checkSoln*phaseFac;
% Compute error in intermediate solution
err_int_dB = 10*log10( norm(trueIntSoln-checkSoln,2)^2 / ...
                                            norm(trueIntSoln,2)^2 );


%% Step 2
% CS recovery 
% We will use Basis Pursuit De-Noising (BPDN)
%       min         || x ||_1 
%   subject to     || Ax - b ||_2 <= eps

% De-Noising parameter, eps
epsilon = 1e-10;                        % Error tolerance

% Time to compute step 2
timer = tic;

cvx_begin quiet
    cvx_solver mosek;
    variable recoveredSig(n,1) complex;
    
    minimize norm(recoveredSig,1)
    subject to
        norm(C*recoveredSig-intSoln) <= epsilon;
cvx_end

% Correct for global phase factor
phaseFac = exp( 1i* angle( (recoveredSig'*testFnc)/(testFnc'*testFnc) ) );
recoveredSig = recoveredSig*phaseFac;

%% Record error and execution times
% Record execution time
runtime( idxTrials, idxSVals ) = step1_time + toc(timer);

% Compute error 
err_dB = 10*log10( norm(testFnc-recoveredSig,2)^2 / norm(testFnc,2)^2 );
error( idxTrials, idxSVals ) = err_dB;

% Check if we have successful reconstruction
if( err_dB > -100 )
    % Unsuccessful reconstruction
    failure( idxTrials, idxSVals ) = 1;
    % Print out status
    fprintf('\n %03d      %+08.3f           %+08.3f              %d             %+08.3f', ...
        idxTrials, err_dB, err_int_dB, false, runtime(idxTrials,idxSVals) );
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
        fprintf( '\nTrial    Error (dB)    Int. Soln. Error (dB)    Success?    Runtime (in secs.)' );
        fprintf( '\n-----    ----------    ---------------------    --------    ------------------' );
    end
else
    % Print out status
    fprintf('\n %03d      %+08.3f           %+08.3f              %d             %+08.3f', ...
        idxTrials, err_dB, err_int_dB, true, runtime(idxTrials,idxSVals) );
    % Update trial counter
    idxTrials = idxTrials+1;    
end


end

% Successfull reconstruction in at least 95% trials
% Store no. of measurements required
nMeasurements( idxSVals) = m;

% Print no. of measurements required
fprintf( '\n\n Success!' );
fprintf( '\n No. of measurements required at sparsity %d is %d', s, m ); 
fprintf( ' Avg. runtime (in secs.) at sparsity %d is % 3.3f\n', s, ...
                                            mean(runtime(:,idxSVals)) );
                                
end


% Save results to file
save( 'sparsePR_runtime_measurements.mat', 'sVals', 'n', 'runtime', ...
                        'nMeasurements', 'error' );


% Plot runtime and no. of measurements
figure('visible', 'on'); 
plot( sVals, mean(runtime), 'r-+', 'linewidth', 2 );
xlabel( 'Sparsity, $s$', 'interpreter', 'latex', 'fontsize', 15 );
ylabel( 'Avg. runtime (in secs.)', 'interpreter', 'latex', 'fontsize', 15 );
grid; title( 'SparsePR Runtime: N=64, Noiseless Measurements', ...
     'interpreter', 'latex', 'fontsize', 15 );
 
% Save plots
print( gcf, '-dpng', 'sparsePR_runtime.png' );
print( gcf, '-depsc', 'sparsePR_runtime.eps' );
print( gcf, '-dpdf', 'sparsePR_runtime.pdf' );
 
figure('visible', 'on'); 
plot( sVals, nMeasurements, 'r-+', 'linewidth', 2 );
xlabel( 'Sparsity, $s$', 'interpreter', 'latex', 'fontsize', 15 );
ylabel( 'No. of measurements', 'interpreter', 'latex', 'fontsize', 15 );
grid; title( 'SparsePR Measurements: N=64, Noiseless Measurements', ...
     'interpreter', 'latex', 'fontsize', 15 );

% Save plots
print( gcf, '-dpng', 'sparsePR_measurements.png' );
print( gcf, '-depsc', 'sparsePR_measurements.eps' );
print( gcf, '-dpdf', 'sparsePR_measurements.pdf' ); 
