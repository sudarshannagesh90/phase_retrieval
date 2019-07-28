function X = prsamp_demo_original(H, Y, x_init)
%% This function solves Y=|HX| for binary H and complex X

[m, n] = size(H);

%% Parameters
delta = 1e-4; % noise variation parameter  
damp = 0.01; % how fast the algorithm converges to a local minima     
maxIter = n/2;
vnf = 0.5; % variance normalization factor   
prec = 1e-5; % convergence criterion. minimum difference between two consequent x
%%
init_c = ones(n,1)*0.1;

opt = struct;
opt.delta = delta;
opt.priorPrms = [1, 0, 1/sqrt(n)];
opt.maxIter = maxIter;
opt.prec = prec;
opt.display = 0;
opt.damp = damp;
opt.vnf = vnf;

% init_a = complex(normrnd(0, 1/n, n, 1), normrnd(0, 1/n, n, 1));
init_a = x_init;
opt.initState = [init_a; init_c];
H = sparse(H);
%X = SwAMP_Calib(Y,H,opt);%Has only been compiled for Windows
X = myprSwAMP(Y,H,opt);%Has only been compiled for Linux
