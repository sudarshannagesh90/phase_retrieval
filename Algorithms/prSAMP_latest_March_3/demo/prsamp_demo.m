function X = prsamp_demo(H, Y,x_init,noise_var,damp_seq,damp_par,par_frac)
%% This function solves Y=|HX| for binary H and complex X

[m, n] = size(H);

%% Parameters
delta=noise_var;
%delta = 1e-4; % noise variation parameter  
% damp = 0.01; % how fast the algorithm converges to a local minima    
% damp = .9;
% damp = 0.1; % how fast the algorithm converges to a local minima    
maxIter = n/2;
% maxIter = 25;
vnf = 0.5; % variance normalization factor   
prec = 1e-5; % convergence criterion. minimum difference between two consequent x
%%
% init_c = ones(n,1)*0.1;
init_c = var(x_init)*ones(n,1);

opt = struct;
opt.delta = delta;
opt.priorPrms = [1, 0, 1/sqrt(n)];
opt.maxIter = maxIter;
opt.prec = prec;
opt.display = 1;
opt.damp_seq = damp_seq;
opt.damp_par = damp_par;
opt.par_frac = par_frac;
opt.vnf = vnf;

% init_a = complex(normrnd(0, 1/n, n, 1), normrnd(0, 1/n, n, 1));
init_a=x_init;
%init_a=[1:length(x_init)]';
opt.initState = [init_a; init_c];
H = sparse(H);
% X = PartialprSwAMP(Y,H,opt);
X = ArbMatprSwAMP(Y,H,opt);
