% Solve the PhaseMax signal reconstruction problem
%         maximize <c,x>
%         subject to |Ax|<=b
% Where 'c' is an nx1 initial guess chosen by the user, 'A' is an mxn 
% measurement matrix, and 'b' is an mx1 vector of non-negative real 
% measurements.
% 
% The problem is solved by approximately enforcing the constraints using a
% quadratic barrier function.  A continuation method is used to increase
% the strength of the barrier function until a high level of numerical
% accuracy is reached.
%   The objective plus the quadratic barrier has the form
%     <-c,x> + 0.5*max{|Ax|-b,0}^2.
%
%  Inputs:
%  A: The linear operator 'A' may be either a matrix or a function 
%    handle.  If 'A' is a function handle, then the user must supply the
%    adjoint of 'A' using the argument 'At'.  If 'A' is a matrix, then the 
%    value of 'At' is discarded. 
%  c: An estimate of the true solution. 'c' can be a column vector, or it 
%    can be a 2D array in the case of recovering an image from phaseless 
%    measurements.
%  b: A column vector of real, non-negative measurements
%  At:  The adjoint of 'A'.  This is only used if 'A' is a function handle.
%
%  Copyright Goldstein & Studer, 2016.  For more details, visit 
%  https://www.cs.umd.edu/~tomg/projects/phasemax/

function [x, outs] = solve_phasemax_grad( c, A, b, niters )

%% Check preconditions
assert(size(c,2)==1,'Input must be column vector');
assert(norm(imag(b))==0,'b must be real');
assert(sum(b<0)==0,'b must be non-negative');
assert(size(b,2)==1,'b must be a column vector');
assert(norm(c)>0==1,'guess vector cannot be zero');
m = size(b,1);
n = size(c,1);

%% Set parameters
maxIters = niters;  % maximum iterations per continuation step
tol = 1e-6;        % The relative tolerance of the numerical method. 

if ~exist('isNonNeg','var')
    isNonNeg = false;
end

%  Normalize the initial guess relative to the number of measurements
c = (c/norm(c))*mean(b)*(m/n)*100;
%  If 'A' is numeric (i.e., a matrix), then normalize the rows
%  Also, convert 'A' to a function handle so that we only need to address
%  the function handle case.
if isnumeric(A)
    % normalize A
    normals = sqrt(sum(A.*conj(A),2));
    A = A./(normals*ones(1,n));
    b = b./normals;
    At = @(x) A'*x;
    A = @(x) A*x;
end

%  re-scale the initial guess so that it approximately satisfies |Ax|=b
x = c.* min(b./abs(A(c)));

%%  Define objective function components for the gradient descent method FASTA
f = @(z) 0.5*norm(max(abs(z)-b,0))^2;     %  f(z) = 0.5*max{|x|-b,0}^2 : This is the quadratic penalty function
gradf = @(z)  (sign(z).*max(abs(z)-b,0)); % The gradient of the quadratic penalty
opts = [];          % Options to hand to FASTA
opts.verbose=1;
opts.recordObjective=true;
opts.maxIters=maxIters;
opts.stopRule = 'residual';
iterCount = 0;     %  A counter to maintain the total iterations across continuation steps
for homotopy = 1:200    % Iterate over continuation steps
    g = @(x) -real(c'*x);   %  The linear part of the objective
    proxg = @(x,t) x+t*c;   %  The proximal operator of the linear objective
    opts.tol=norm(c)/100;   % use a tighter tolerance when the solution is more exact
    [x,outs]= fasta(A, At, f, gradf, g, proxg, x, opts );  % Call a solver to minimize the quadratic barrier problem
    opts.tau = outs.stepsizes(end);     % Record the most recent stepsize for recycling.
    iterCount = iterCount+outs.iterationCount;  % Record total iterations used so far
    c = c/10;  % do continuation - this makes the quadratic penalty stronger
    if max(abs(A(x))-b) <tol;   % test how accurately the constraints are satisfied
        break;
    end
end
outs.totalIterations = iterCount;

end

