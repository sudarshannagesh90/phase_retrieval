%% Script for using phase-retrieval algorithms.
clear; clc; close all;
addpath('Algorithms/');

%% Problem parameters
n = 16^2;
p = 12 * n;
sigma_w = 1e-8;                       % noise standard-deviation
meas_type = 1;                        % 1: 0, 1 measurements
                                      % 2: -1, 1 measurements
                                      % 3: Gaussian measurements
                                      
%% Algorithm parameters
alg_type = 1;                         % 1: Gerchberg-Saxton algorithm
                                      % 2: prVAMP
n_iters = 500;                        % Number of iterations

%% Make measurements
x_o = randn(n, 1) + 1j * randn(n, 1);
switch meas_type
    case 1
        A = round(rand(p, n));
    case 2
        A = 2 * round(rand(p, n)) - 1;
    case 3
        A = randn(p, n) + 1j * randn(p, n);
end
noise_vec = sigma_w * (1/sqrt(2) * randn(p, 1) + 1j * 1/sqrt(2) * randn(p, 1));
y = abs(A * x_o + noise_vec);

%% Retrive phase from magnitude measurements
switch alg_type
    case 1
        A_dag = pinv(A);
        x_recovered = GS(y, A, A_dag, n_iters);
end

%% Unwrap phase of the solution
x_unwrapped = x_recovered * exp(-1j * angle(trace(x_o' * x_recovered)));
x_nrmse = norm(x_o - x_unwrapped, 'fro') / norm(x_o, 'fro');
fprintf('NRMSE: %f', x_nrmse)

%%