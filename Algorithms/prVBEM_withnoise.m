function [x_hat,KLdiv]=prVBEM_withnoise(y,D,sigma,niters,opt)
      

% Mean-Field approximation for phase retrieval
%*********************************************
% This matlab function implements the phase recovery algorithm described in
% "Phase recovery from a Bayesian point of view: the variational approach"
% by A. Drémeau and F. Krzakala
%
%
% AUTHORS: A. Drémeau
%
%
% INPUTS:
% y: vector of observations
% D: dictionary
% opt: options
% .var_a: variance of the non-zero coefficients in x
% .var_n: variance of Gaussian noise n (default: 10^-8))
% .pas_est: used in the stopping criterion to force the decrease of the
% estimated noise variance (default: 0.1)
% .flag_est_n: if 'on', the noise variance is estimated at each
% iteration (default:'on')
%
%
% OUTPUTS:
% x_hat: recovered signal
% KLdiv: value of the Kullback-Leibler divergence at stopping point
%




% Initialization
%*****************                       
[dim_y,dim_s]=size(D);

% default parameters
if nargin < 5,
    tmp_x  = pinv(D)*y;
	var_a  = max(abs(tmp_x))^2;
    phase = exp(1i*2*pi*rand(dim_y,1));
	x0 = pinv(D)*(y.*phase);
    
    opt.var_a      = var_a;
    opt.var_n      = sigma;
    opt.x0         = x0;
    opt.flag_est_n = 'on';
    opt.niter      = niters;
    opt.pas_est    = 0.1;
    opt.flag_cv    = 'KL';
end

var_a      = opt.var_a;
var_n      = opt.var_n;
x0         = opt.x0;
flag_est_n = opt.flag_est_n;
flag_cv    = opt.flag_cv;
pas_est    = opt.pas_est;
niter      = opt.niter;
if var_n==0
    var_n = 10^(-6);
    warning('Noise variance set to a small value for implementation reasons.')
end

moy_x   = zeros(dim_s,1);
var_x   = zeros(dim_s,1);

if ~isempty(x0)
    moy_x(:)=x0;
end
ybar=y;
z=D'*ybar;
H=D'*D;

var_x(:)=var_a*var_n./(var_n+var_a*diag(H));
Dm=D*moy_x;

var_n_true = var_n;


% Iterative process
%*******************
% Convergence criteria
OK_outer=1;
compt=1;
KLdiv_old=100000;

while OK_outer && compt<niter
        
    % Estimation of var_n  
    %*********************
    if strcmp(flag_est_n,'on') || strcmp(flag_est_n,'on_off') 
        var_n=inv(dim_y)*(Dm'*Dm+var_x'*diag(H)+y'*y-2*real(ybar'*Dm));
        if ~isreal(var_n)
            error('var_n is not real')
        end
    end
    
	% Update of the q(x_i)=G(moy_x,var_x)
    %************************************
    for k=1:1:dim_s
        
        val_tmp=z(k)-H(k,:)*moy_x+H(k,k)*moy_x(k);
        moy_x(k)=(var_a/(var_n+var_a*H(k,k)))*val_tmp;
        var_x(k)=var_n*var_a/(var_n+var_a*H(k,k));
        
    end 
    
    t=conj(y).*(D*moy_x);
    phi = t./abs(t);
    I0=besseli(zeros(dim_y,1),(2/var_n)*abs(t));
    I1=besseli(ones(dim_y,1),(2/var_n)*abs(t));
    fac_bessel=I1./I0;
    if ~isempty(find(isnan(fac_bessel)==1,1))
        fac_bessel(isnan(I1./I0)==1)=1;
    end
    ybar=y.*phi.*fac_bessel;
    z=D'*ybar;
    Dm=D*moy_x;    
    
    % Convergence checks
    %*******************            
    if strcmp(flag_cv,'KL')

        % Computation of the KL divergence
        % int log p(y,x) = int log p(y|x) + int log p(x)
        frst = (dim_y)*log(var_n)+(1/var_n)*(Dm'*Dm+var_x'*diag(H)+y'*y-2*real(ybar'*Dm));
        scnd = (dim_s)*log(var_a)+(1/var_a)*sum((var_x(:,1)+abs(moy_x(:,1)).^2));
        
        % int log q(x)
        sxth = sum(log(var_x));
        vec_tmp = abs(t).*fac_bessel;
        I0_tmp=I0;
        svth = (2/var_n)*sum(vec_tmp)-sum(log(nonzeros(I0_tmp)));

        KLdiv = frst + scnd - sxth + svth;
        diff  = KLdiv_old - KLdiv;
        
        if isnan(diff) || isinf(abs(diff))
            warning('Stopping criterion changed to fixed number of iterations.')
            flag_cv = 'iter';
        end

        if diff<1e-8
            if (var_n/var_n_true)>1.1
                var_n = var_n - pas_est*(var_n-var_n_true);
                flag_est_n = 'off';

                % Computation of the KL divergence
                % int log p(y,x) = int log p(y|x) + int log p(x)
                frst = (dim_y)*log(var_n)+(1/var_n)*(Dm'*Dm+var_x'*diag(H)+y'*y-2*real(ybar'*Dm));
                scnd = (dim_s)*log(var_a)+(1/var_a)*sum((var_x(:,1)+abs(moy_x(:,1)).^2));

                % int log q(x)
                sxth = sum(log(var_x));
                vec_tmp = abs(t).*fac_bessel;
                I0_tmp=I0;
                svth = (2/var_n)*sum(vec_tmp)-sum(log(nonzeros(I0_tmp)));

                KLdiv = frst + scnd - sxth + svth;

            else
                OK_outer=0;
            end              
        end
        
        KLdiv_old = KLdiv;
        
    end
    compt=compt+1;
end
x_hat = moy_x;
KLdiv = KLdiv_old;


