function X = rec_maxcut_impl (M,a_priori,opts)

% Solves the problem :
%   Find max(Trace(MX)) s.t. X >= 0, diag(X)=1
% The dual problem is :
%   Find min(scal(y,ones)) s.t. M <= Diag(y)
% The method is as described in "An interior-point method for
% semidefinite programming", by Helmberg, Rendl, Vanderbei and
% Wolkowicz and follows the code reproduced on page 21.

% Display
if (nargin>2) && (isfield(opts,'disp'))
    disp_b = opts.disp ;
else
    disp_b = 1 ;
end

% Precision
if (nargin>2) && (isfield(opts,'digits'))
    digits = opts.digits ;
else
    digits = 6 ;
end

if a_priori.real
    n = size(M,1)/2 ;
else
    % Come down to the case where M is real
    n = size(M,1) ;
    M = [real(M),-imag(M) ; imag(M),real(M)] ;
end
    
Mscale = norm(M) ;
b = ones(n,1);            % Any b>0 works just as well

% Initial admissible values for primal and dual problems
X = eye(2*n) ;
if a_priori.real ; X = X/2 ; end
ybuf = sum(abs(M))' * 1.1 ;
mbuf = max(ybuf) ; ybuf = ybuf+0.01*mbuf ; % To avoid problems if
                                           % some coefficients of M
                                           % are zero

if a_priori.real                           % initial y
    y = ybuf(1:n)+ybuf(n+1:2*n) ;
else
    y = ybuf(1:n) ;                        % only half of y,
                                           % because of symmetries
end
Z = diag([y;y]) - M ;

if a_priori.real                           % initial dual
    phi = b'*y ;
else
    phi = 2*b'*y ;
end
psi = M(:)' * X(:) ;                       % and primal costs
mu = Z(:)'*X(:)/(4*n) ;                    % initial complementarity

iter = 0 ;                                 % iteration count

if disp_b
    disp(['iter    alphap  alphad    gap       lower      ' ...
          'upper     digit']); drawnow ;
end

while phi-psi > max([1,Mscale]) * 10^(-digits)
        
    % Start a new iteration
    iter = iter + 1 ;
    
    [U,S,V] = svd(Z) ;
    Zi = V*diag(diag(S).^(-1))*U' ;
    Zi = (Zi + Zi')/2 ;
    Zi_diag = diag(Zi) ;
    Zi_diag = Zi_diag(1:n)+Zi_diag(n+1:2*n) ;
    Hm = X.*Zi ;
    Hm = Hm(1:n,:) + Hm(n+1:2*n,:) ; Hm = Hm(:,1:n) + Hm(:,n+1:2*n) ;
    dy =  Hm \ (mu * Zi_diag - b) ;       % Solve for dy
    dZ =  diag ([dy;dy]) ;
    dX = - Zi * dZ * X + mu * Zi - X ;    % Back substitute for dX
    dX = (dX + dX') / 2 ;
    
    % Line search on primal
    alphap = 1;                                % initial steplength
    [dummy,posdef] = chol( X + alphap * dX );  % test if pos.def 
    while posdef > 0,
        alphap = alphap * .7;
        [dummy,posdef] = chol( X + alphap * dX );
    end ;
    alphap = alphap * .98 ;                    % stay away from boundary
    % Line search on dual
    alphad = 1;                                % initial steplength
    [dummy,posdef] = chol( Z + alphad * dZ );
    while posdef > 0;
        alphad = alphad * .7;
        [dummy,posdef] = chol( Z + alphad * dZ );
    end;
    alphad = alphad * .98 ;
    
    % Update
    X = X + alphap * dX;
    y = y + alphad * dy;
    Z = Z + alphad * dZ;
    mu = X(:)' * Z(:) / (4*n) ;
    if alphap + alphad > 1.8, mu = mu/2 ; end ; % Speed up for long steps
    if a_priori.real
        phi = b'*y ;
    else
        phi = 2*b'*y ;
    end
    psi = M(:)' * X(:) ;

    % Display iteration count
    if disp_b && (mod(iter,10)==0)
        fprintf('%2.0f\t%.3f\t%.3f\t  %.1e   %.1e   %.1e   %.1e\n',iter,alphap, ...
                alphad,(phi-psi),psi,phi, ...
                -log((phi-psi)/max(1,Mscale))/log(10)) ;
    end

end      % End of main loop
