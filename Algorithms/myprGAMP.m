function [x_hat,residual] = myprGAMP(Beta,wvar,x_init,x_o,y,iters,width,height,denoiser1,stepsize,gout,prmts,A_func,At_func)
% function x_hat = DAMP(y,iters,width,height,denoiser,M_func,Mt_func)
% this function implements D-AMP based on any denoiser present in the
% denoise function
% Input:
%       y       : the measurements 
%       iters   : the number of iterations
%       width   : width of the sampled signal
%       height  : height of the sampeled signal. height=1 for 1D signals
%       denoiser: string that determines which denosier to use. e.g.
%       denoiser='BM3D'
%       M_func  : function handle that projects onto M. Or a matrix M.
%       Mt_func : function handle that projects onto M'. Or no entry
%Output:
%       x_hat   : the recovered signal.

if nargin==14%function
    error('Function handles not supported yet');
    A=@(x) A_func(x);
    At=@(z) At_func(z);
else%Matrix
    A_norm2=norm(A_func,'fro')^2;
    A=@(x)A_func*x;
    At=@(z)A_func'*z;
end

if isequal(gout,'phaseless')
    eta_out=@(phat,pvar,y,sigma2_w) eta_out_phaseless(phat,pvar,y,sigma2_w);
else
    eta_out=@(phat,pvar,y,sigma2_w) eta_out_gaussian(phat,pvar,y,sigma2_w);
end
% denoi_in=@(noisy,sigma_hat) mydenoise(noisy,sigma_hat,width,height,denoiser1);

n=width*height;
m=length(y);

% x=mean(x_o(:))*ones(n,1);%From GAMP paper: This should be the expected value of x based on it's priors
x=x_init;
sigma2_x=var(x);%This may not be this value is actually set. Should be around 2.5e3
s=eps*ones(m,1);
sigma2_p=eps;
x_bar=x_init;
sigma2_s=eps;

for i=1:iters
    sigma2_p=Beta*A_norm2/m*mean(abs(sigma2_x))+(1-Beta)*sigma2_p;
    sigma2_p=max(sigma2_p,eps);
    
    p=A(x)-sigma2_p*s;
    [g,dg] = eta_out(p,sigma2_p,y,wvar);
    s=Beta*g+(1-Beta)*s;
    sigma2_s=-Beta*dg+(1-Beta)*sigma2_s;
    stepx=stepsize;%This shouldn't be neceessary.  Leave it for now then figure out how to remove it.
    sigma2_r=stepx/(A_norm2/n*sigma2_s);
    x_bar=Beta*x+(1-Beta)*x_bar;
    r=x_bar+sigma2_r*At(s);
    [x,sigma2_x]=ComplexBournoulliGaussian( r, sigma2_r, prmts );
%     eta=randn(1,n);
%     epsilon_in=max(abs(r))/1000+eps;
%     div_in=eta*((denoi_in(real(r)+epsilon_in*eta',sqrt(var_used))-x)/epsilon_in);
%     sigma2_x=max(sigma2_r*div_in/n,0);
    residual(i) = norm(y(:)-abs(A_func*x(:)))/norm(y);
end
x_hat=reshape(x,[height width]);
x_hat=conj(x_hat');
end
function [g,dg] = eta_out_phaseless(phat,pvar,y,sigma2_w)
        var0=sigma2_w;
        y_abs=y;
        m=length(y);
        phat_abs=abs(phat);
        B = 2*y_abs.*phat_abs./(var0+pvar);
        I1overI0 = min( B./(sqrt(B.^2+4)), ...
            B./(0.5+sqrt(B.^2+0.25)) );%upper bounds (11)&(16) from Amos
        y_sca = y_abs./(1+var0./pvar);
        phat_sca = phat_abs./(1+pvar./var0);
        zhat = (phat_sca + y_sca.*I1overI0).*sign(phat);
        sigma2_z = y_sca.^2 + phat_sca.^2 ...
            + (1+B.*I1overI0)./(1./var0+1./pvar) ...
            - abs(zhat).^2;
        g=1/pvar*(zhat-phat);
        dg=1/pvar*(mean(sigma2_z)/pvar-1);%Using uniform variance.
end
function [g,dg] = eta_out_gaussian(phat,pvar,y,sigma2_w)
    g=(y-phat).*(pvar+sigma2_w).^-1;
    dg=-(pvar+sigma2_w).^-1;
end
% function [a, c] = ComplexBournoulliGaussian( r, sig, prmts )
%     rho = prmts(1);%Sparsity
%     pr_mean = prmts(2);
%     pr_var = prmts(3);
% 
%     isv = 1 ./ (pr_var + sig);
%     rsc = .5 .* (pr_mean - r) .* (pr_mean - r) .* isv;
%     eff = (pr_mean .* sig + r .* pr_var) .* isv;
%     vrp = pr_var .* sig .* isv;
% 
%     gamma = ((1. - rho) / rho) .* sqrt(pr_var ./ vrp) .* ...
%         exp(-.5 * r .* r ./ sig + rsc);
% 
%     a = eff ./ (1 + gamma);
%     c = bsxfun( @max, gamma .* a .^ 2 + vrp ./ (1 + gamma), 1e-19 );
% end

function [xhat,xvar]=ComplexGaussian(r, sigma2_r, prmts )
            uhat0 = prmts(2);
            uvar0 = prmts(3); 
            % Compute posterior mean and variance
            gain = uvar0./(uvar0+sigma2_r);
            xhat = gain.*(r-uhat0)+uhat0;
            xvar = gain.*sigma2_r;
end
function [xhat,xvar]=ComplexBournoulliGaussian(r, sigma2_r, prmts )
            py1=prmts(1);
            py0=1-prmts(1);
            x0_mean=0;%Harcoded, 0 with prob py0
            [xhat1,xvar1]=ComplexGaussian(r,sigma2_r,prmts);
            xhat = py1.*xhat1 + py0.*x0_mean;
            xvar = py1.*(abs(xhat1).^2 + xvar1) + py0.*(abs(x0_mean).^2)...
                - abs(xhat).^2;
end
