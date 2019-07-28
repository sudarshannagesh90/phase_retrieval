function [X,phases] = rec_maxcut (M,a_priori,opts)

% X is the maximizer of the convex problem
% phases is the projection of X's main eigenvector onto the set of
%      vectors whose coefficients are all unitary complexes.
% opts is facultative ; it may contain :
%      digits : precision
%      disp : 0 to prevent the computer from displaying information

n = size(M,1) ;
if a_priori.real ; n = n/2 ; end

warning('off') ;

if (nargin<3) ; opts = [] ; end
X = rec_maxcut_impl (M,a_priori,opts) ;           % Solve
    
[phases,d] = eigs(X,1,'lm',struct('disp',0)) ;    % Get main eigenvector
phases = phases(1:n) + i*phases(n+1:2*n) ;
phases = phases + (abs(phases)<10^(-4)) ;
phases = phases./abs(phases) ;             % Renormalize

if not(a_priori.real)
    X = (X(1:n,1:n)+X(n+1:2*n,n+1:2*n) + ...
         i*(-X(1:n,n+1:2*n)+X(n+1:2*n,1:n))) ...
        /2 ;
end