function [M,G] = relax (A,b,a_priori)

% Calculates the matrix M associated to the relaxed problem.
% The maxcut algorithm will return the minimizor of Tr(MX).

% The matrix G expands a vector of phases whose coefficients at
% which phases are known to be equal have been merged.
% Reason : when some phases are known, b is given with these known
% phases. One thus has to reconstruct a vector of phases, whose
% coefficients at the known places are all equal, to 1. So these
% coefficients are merged and only one is reconstructed for all of
% them. After reconstruction, the matrix G will be used to recover
% all coefficients from the calculated ones.

n = size(A,1) ; % Number of measurements
N = size(A,2) ; % Size of signal

G = eye(n) ;
G(a_priori.known_indexes,a_priori.known_indexes) = 1 ;
for k=1:n
    k_equals = find(G(k,:)) ;
    G(:,k_equals(2:end)) = [] ;
end

if a_priori.sym
    S = [eye(N/2);fliplr(eye(N/2))] ;
    A = A*S ;
end

if a_priori.real
    Ar = [real(A);imag(A)] ;
    Ar_dag = pinv(Ar,1e-10) ;
    Br = [diag(real(b)),-diag(imag(b)) ;
          diag(imag(b)),diag(real(b))] ;
    Gr = [G,zeros(size(G));zeros(size(G)),G] ;
    M = (Ar*Ar_dag-eye(2*n))*Br*Gr ;
    M = M'*M ;
else
    A_dag = pinv(A,1e-10) ;
    M = (A*A_dag-eye(n))*diag(b)*G ;
    M = M'*M ;
end