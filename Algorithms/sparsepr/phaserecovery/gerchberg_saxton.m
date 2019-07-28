function f_rec = gerchberg_saxton (f_init,A,b,proj,a_priori,nb_its)

% Gerchberg-Saxton algorithm

% f_init is the initial guess
% A is the matrix of the linear transformation applied to the entry signal
% b is a vector containing the moduli
% proj is the projection onto Range(A)
% a_priori contains the informations about the signal (only the
% indexes of the known phases will be used)
% nb_its is the number of projections which will be performed (optional)

% At each step, the linear transformation (named hj) of the current
% guess is calculated. It is then projected onto the space of
% vectors whose moduli are correct, then onto the image of the
% linear application A, which gives a new guess for the signal.

% Remark : the symmetry constraint is not handled (but if it proved
% useful, it could easily be done)

if (nargin<6)
    nb_its = 3000 ;
end
n = size(A,1) ; % Number of measurements
N = size(A,2) ; % Size of signal

hj = A*f_init ;

for k=1:nb_its
   
    % Adjustment of moduli and known phases
    hj = hj +(abs(hj)<10^(-6)*max(abs(hj))) ;
    hj = hj.*b./abs(hj) ;
    hj(a_priori.known_indexes) = b(a_priori.known_indexes) ;
    
    % Projection onto Range(A)
    hj = proj(hj) ;
    
end

[hj,f_rec] = proj(hj) ;