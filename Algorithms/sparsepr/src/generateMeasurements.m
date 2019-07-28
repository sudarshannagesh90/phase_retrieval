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

function [measurements, measurementMat, noise] = generateMeasurements( ...
                    measurementType, trueSignal, m, mCS, addNoise, snr )

% generateMeasurements Generate magnitude measurements for Phase Retrieval
%   Generate magnitude measurements 
%           b_i = |<p_i, x>|^2 + n_i,       i = 1,...,m
%               where   p_i \in C^n is a measurement vector
%                       n_i \in R is measurement noise
%                       x \in C^n is the unknown signal
%                       b_i \in R is the generated measurement
%   to simulate and test the Phase Retrieval problem. 
%
%   Inputs:
%       measurementType - type of measurements to generate; choices are
%                           randgauss - random Gaussian
%       trueSignal      - the underlying true signal
%       m               - No. of measurements
%       mCS             - CS problem dimension
%       addNoise        - Do we add noise to the measurements?
%       snr             - If we add noise, at what snr (dB)?
%
%   Outputs:
%       measurements    - the generated measurements
%       measurementMat  - Matrix used to generate measurements
%                           this is a structure containing the measurement
%                           matrix as well as its components P and C
%       noise           - Structure containing noise parameters
%


% Length of true signal
n = length(trueSignal);

% Generate measurement matrix
switch lower(measurementType)
    case 'randgauss'
        % Random (complex) Gaussian measurements
        C = randn(mCS, n) + 0i*randn(mCS, n);   % CS matrix
        P = randn(m, mCS) + 1i*randn(m, mCS);   % Phase Retrieval matrix
        M = P*C;
end

% Store components of measurement matrix in a structure
measurementMat.M = M;
measurementMat.P = P;
measurementMat.C = C;

% Generate magnitude measurements
measurements = abs( M*trueSignal ).^2;

% Add Noise
if( addNoise )        
    % Calculate signal power
    % This is just the norm^2 of the measurements divided by dimension
    signal_power = (1/m)*norm( measurements )^2;

    % Add noise
    % Noise model: Additive uniform random noise
    % b_i^ = b_i + n_i, where
    %   b_i is ith entry of measurement vector
    %   n_i is the added noise
    %   b_i^ is the ith entry of noise corrupted measurement
    % n_i ~ U[0,a] is iid uniform random noise
    % Note: mean(n) = a/2, var(n) = (a^2)/12
    % We choose a as follows:
    %   mean^2 + var^2 = noise_power
    %   (1/4 + 1/12)*a^2 = a^2/3 = noise_power
    % Noise power
    noise_power = signal_power / ( 10^(snr/10) );

    % Choose uniform distribution parameter a
    alpha = sqrt( 3*noise_power );
    noiseVec = alpha*rand( size(measurements) );
    
    % Additive noise
    measurements = measurements + noiseVec;    
end

% Store noise parameters in structure
noise.addNoise = addNoise;
noise.snr = snr;
if( addNoise )
    noise.signal_power = signal_power;
    noise.noise_power = noise_power;
    noise.noiseVec = noiseVec;
end

end
