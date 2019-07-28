function f = reconstruct (A,b,phi,a_priori)

% Reconstructs the signal.
% A is the considered linear application.
% b is the vector of moduli.
% phi is the vector of reconstructed phases.
% a_priori contains the informations about the signal

n = size(A,1) ; % Number of measurements
N = size(A,2) ; % Size of signal

if a_priori.sym
    S = [eye(N/2);fliplr(eye(N/2))] ;
    A = A*S ;
end

% Determination of f such that Af is as close to b.*phi as possible.
if a_priori.real
    Ar = [real(A);imag(A)] ;
    rec_meas = b.*phi ;
    rec_meas = [real(rec_meas);imag(rec_meas)] ;
    Ar_dag = pinv(Ar,1e-10) ;
    f = Ar_dag*rec_meas ;
else
    A_dag = pinv(A,1e-10) ;
    f = A_dag*(b.*phi) ;
end

if a_priori.sym
    f = [f;flipud(f)] ;
end