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

function [measurements, measurementMat, noise] = generateCSMeasurements( ...
                        measurementType, trueSignal, m, addNoise, snr )

% generateCSMeasurements Generate linear measurements for testing CS
%   Generate linear measurements for simulating compressive sensing
%
%   Inputs:
%       measurementType - type of measurements to generate; choices are
%                           randgauss - random Gaussian
%       trueSignal      - the underlying true signal
%       m               - No. of measurements
%       addNoise        - Do we add noise to the measurements?
%       snr             - If we add noise, at what snr (dB)?
%
%   Outputs:
%       measurements    - the generated measurements
%       measurementMat  - Matrix used to generate measurements
%       noise           - Structure containing noise parameters


% Length of true signal
n = length(trueSignal);

% Generate measurement matrix
switch lower(measurementType)
    case 'randgauss'
        % Random Gaussian measurements
        measurementMat = randn(m,n);
        
        % Scale to have unit-normed columns
        measurementMat = bsxfun( @times, measurementMat, ...
                                    1./sqrt(sum(measurementMat.^2,1)) );

end

% Generate measurements
measurements = measurementMat*trueSignal;

% Add Noise
if( addNoise )        
    % Calculate signal power
    % This is just the norm^2 of the measurements divided by dimension
    signal_power = (1/m)*norm( measurements )^2;

    % Add iid additive zero mean Gaussian noise
    % Noise power
    noise_power = signal_power / ( 10^(snr/10) );

    % Additive noise
    noiseVec = sqrt(noise_power) * randn(m,1);
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
