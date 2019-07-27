function [x_hat]=GS(y,A,Adag,iters)
x_hat=Adag*y;
for i=1:iters
    z=y.*exp(1i*angle(A*x_hat));
    x_hat=Adag*z;
end
end
