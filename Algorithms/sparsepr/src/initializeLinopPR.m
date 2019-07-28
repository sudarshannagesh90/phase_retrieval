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

function tfocsLinOp = initializeLinopPR( size , PhMat )

% initializeLinopPR Initialize Phase Retrieval measurement operator
%   Initialize a linear operator for use with tfocs given a phase retrieval
%   measurement matrix
%
%   Inputs:
%       size        - dimensions of the Phase Retrieval matrix
%       PhMat       - Phase retrieval matrix
%
%   Output(s):
%       tfocsLinOp  - linear operator for use with tfocs
%
%   Adapted from 
%   Code for 'Phase Retrieval from Coded Diffraction Patterns'
%   http://web.stanford.edu/~mahdisol/PRcode.html
%   Mahdi Soltanolkotabi
%   Oct 2014
%

tfocsLinOp = @(X, mode) PRFullMatrixMeasurements(size, PhMat, X, mode);
    
end

% As required by tfocs, use anonymous functions to return one of three
% modes
function y = PRFullMatrixMeasurements(size, PhMat, X, mode)

switch mode,
    % Return the size
    case 0, 
        y = size;
        
    % Return the operator
    case 1, 
        % Note that the squared magnitude measurements may be expressed in
        % the following form (where X = xx*)
        y = diag(PhMat*X*PhMat');

    % Return 
    case 2,
        y = PhMat' * diag(X) * PhMat;
end

end
