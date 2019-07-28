addpath(genpath('.'));
n=8^2;
m=12*n;
x_o=complex(normrnd(0, 1/n, n, 1), normrnd(0, 1/n, n, 1));

A=round(rand(m,n));
sigma=1e-4;
y=abs(A*x_o+sigma*(sqrt(2)/2*randn(m,1)+1i*sqrt(2)/2*randn(m,1)));

noise_var=sigma^2;

%prVBEM
x_hat_prVBEM=prVBEM(y,A);

%prSAMP
iters=100;
x_hat_prSAMP_latest=prsamp_demo_newchannel(A, y,x_hat_prVBEM,noise_var,.9,.9,1e5,iters);


prVBEM_error=norm(disambiguate(x_hat_prVBEM,x_o)-x_o)
prSAMP_error=norm(disambiguate(x_hat_prSAMP_latest,x_o)-x_o)