function X = prsamp_demo_newchannel(H, Y,x_init,noise_var,damp_seq,damp_par,par_frac,maxIter)
%% This function solves Y=|HX| for binary H and complex X

[m, n] = size(H);

%% Parameters
if length(noise_var(:))==1
    delta=ones(m,1)*noise_var;
else
    delta=noise_var;%noise_var should be vector of length m
end
%delta = 1e-4; % noise variation parameter  
% damp = 0.01; % how fast the algorithm converges to a local minima    
% damp = .9;
% damp = 0.1; % how fast the algorithm converges to a local minima    
% maxIter = n/2;
% maxIter=100;
% maxIter = 25;
vnf = 0.5; % variance normalization factor   
prec = 1e-5; % convergence criterion. minimum difference between two consequent x
%%
% init_c = ones(n,1)*0.1;
init_c = var(x_init)*ones(n,1);

opt = struct;
opt.delta = delta;
opt.priorPrms = [1, 0, var(x_init)];
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
%X = PartialprSwAMP(Y,H,opt);
X = SwAMP_Calib(Y,H,opt);
