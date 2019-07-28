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

function testFnc = generateTestFunction( testFncType, n, varargin )

% generateTestFunction Generate test functions for Phase Retrieval problem
%   Generate test functions to simulate and test the Phase Retrieval
%   problem. 
%
%   Inputs:
%       testFncType - type of test function to generate; possible choices
%                       randgauss - random Gaussian
%                       unirand   - uniform random
%                       sinusoid  - sinusoidal
%                       sparse    - sparse signal
%       n           - length of test signal
%   Optional:
%       nComps      - No. of components in test signal (for sinusoids and
%                     sparse signals)
%
%   Output(s):
%       testFnc     - the generated test function


% Process optional input arguments
if( nargin == 3 )
    nComps = varargin{1};
end

% Generate test function
switch lower(testFncType)
    case 'randgauss'
        % Random (complex) Gaussian test function
        testFnc = randn(n,1) + 1i*randn(n,1);
        
    case 'unirand'
        % Uniform (complex) Random test function in [-0.5,0.5]
        testFnc = (rand(n,1) - 0.5) + 1i*(rand(n,1) - 0.5);
        
    case 'sparse'
        % Sparse (complex) test function
        % There are 'nComp' non-zero entries
        % Each is a (complex) uniform random entry in [-0.5,0.5]
        testFnc = zeros(n,1);
        testFnc( randperm(n,nComps) ) = (rand(nComps,1) - 0.5) + ...
                1i*(rand(nComps,1) - 0.5);
                
    case 'sinusoid'
        % Trigonometric test function
        % There are 'nComp' frequencies
        % Amplitude of each component is uniform random (and independent)
        % in [-0.5,0.5]
        % Frequencies randomly selected in [1,floor(n/2)]
        testFnc = zeros(n,1);
        for idx =1:nComps
            testFnc = testFnc + ...
                (rand-0.5)*exp(2i*pi*randi(floor(n/2))*(0:n-1)/n).';
        end
end

% Normalize so that test functions have unit norm
testFnc = testFnc/norm(testFnc);

end
