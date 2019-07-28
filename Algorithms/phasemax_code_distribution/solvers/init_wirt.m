%  Compute the intializer proposed for truncated Wirtinger flow.
%  The user may choose whether to use truncation by setting the boolean
%  argument 'isTruncated'.  The user may also turn on/off a least-squares
%  based method for determining the optimal scale of the initializer.
%     Finally, the measurement matrix 'A' may be either numeric/matrix
%  valued, or a function handle.  When a function handle is used, the
%  value of 'N' (the length of the unknown signal) and 'At' (a function 
%  handle for the adjoint of 'A') must be supplied.  When 'A' is numeric,
%  the values of 'At' and 'N' are ignored and inferred from the arguments.
%
%  Copyright Goldstein & Studer, 2016.  For more details, visit 
%  https://www.cs.umd.edu/~tomg/projects/phasemax/

function [ x ] = init_wirt( A, b0, isTruncated, isScaled, N, At )

M = size(b0,1);
% Truncated Wirtinger flow initialization
alphay = 3; %3 (4 also works fine)
lambda0 = sqrt(1/M*sum(b0.^2));
idx = b0<=alphay*lambda0;
if ~isTruncated
    idx = ones(size(b0));
end

if isnumeric(A)
    N = size(A,2);
    [V,sigma] = eigs(1/M*A'*diag(idx.*b0.^2)*A,1);
    [~,maxidx] = max(sigma);
    x = V(:,maxidx)*lambda0*sqrt(M*N/norm(A,'fro')^2); 
    A = @(x) A*x;  % create function handle for later use
else
    x = randn(N,1);
    x0 = x+1;
    op = @(x) 1/M*At((idx.*b0.^2).*A(x));
    count = 0;
    while count<100 && norm(x-x0)/norm(x0)>1e-2
        x0 = x;
        x = op(x);
        x = x/norm(x);
        count = count + 1;
    end
end 


% rescale the solution to have approximately the correct magnitude 
if isScaled
 b = b0(idx);
 Ax = abs(A(x));
 Ax = Ax(idx,:);
 % solve min_s || s|Az| - b ||
 s = Ax'*b/(Ax'*Ax);
 x = x*s;
end

end

