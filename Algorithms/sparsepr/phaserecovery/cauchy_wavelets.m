function [phi,psi] = cauchy_wavelets (N,nb_filters,p)

% Creates a family of Cauchy wavelets (in Fourier domain)
% N : size of each wavelet
% nb_filters : number of wavelets
% p : exponent

lambda = 0.01;
epsilon = lambda.^(1/p) / exp(1);
m = -log(epsilon) + log(-log(epsilon));
omega = (0:N-1)*m*1.5/N ;
S = zeros(N,1);

for j=0:nb_filters-2
    psi(:,j+1) = (2^j * omega).^(p).* exp(-p*(2^j * omega-1));
    S = S + abs(psi(:,j+1).^2);
end
psi = psi ./ sqrt(max(S));

a = 0.45 ;

phi = zeros(1,N);
phi(N/2+1) = 0;
for n=N/2:-1:1
    phi(n) = phi(n+1) +(2^(nb_filters-1-a) * omega(n)).^(p) ...
        .* exp(-p*(2^(nb_filters-1-a) * omega(n)-1));
end
phi = phi ./phi(1);
phi(N/2+2:N) = fliplr(phi(2:N/2)) ;
phi = phi.' ;