function [z] = WF(y, A, iters, tau0)

[m, n] = size(A);
y = y.^2;

%% Initialization
npower_iter = 50;       % Number of power iterations
z0 = randn(n, 1); 
z0 = z0 / norm(z0, 'fro');
			% Initial guess
for tt = 1:npower_iter
	z0 = A' * (y.* (A * z0));
	z0 = z0 / norm(z0, 'fro');
end

normest = sqrt(sum(y) / numel(y)); 
			% Estimate norm to scale eigenvector
z = normest * z0;	% Apply scaling

%% Loop
mu = @(t) min(1 - exp(-t / tau0), 0.2);
			% Schedule for step size

for t = 1:iters
	yz = A*z;
	grad = 1/m * A' * ((abs(yz).^2 - y) .* yz);
			% Wirtinger gradient
	z = z - mu(t) / normest^2 * grad;
			% Gradient update
end

