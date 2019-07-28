function [x_PC] = PhaseCut(b,A,niter)
b = sqrt(b);
[nb_meas,~] = size(A);
nb_its_gs = niter;
 
a_priori.real = false ; % The signal is real
a_priori.sym = false ; % The signal is symmetric
a_priori.known_indexes = [] ; % Indexes where the phases are known
                              % For wavelet transform / random
                              % filters, if the phase corresponding
                              % to the last wavelet / filter is
                              % known, replace [] by :
                              % [nb_meas-N+1:nb_meas]
string_filters = 'Gaussian_measurements' ; % Name of the filters
string_function = 'Gaussian_noise' ; % Name of the test function
                                     % See files
                                     % 'create_function.m' and
                                     % 'create_A_matrix' for a list
                                     % of admissible names.
% Make the signal compatible with the a priori informations
if a_priori.real, f0 = real(f0) ; end ;
if a_priori.sym, f0 = 0.5*(f0+flipud(f0)) ; end ;

mod_indexes = [1:nb_meas] ;
mod_indexes(a_priori.known_indexes) = [] ;

[M,G] = relax(A,b,a_priori) ;

% Apply the maxcut algorithm
tic;[X,phi]=rec_maxcut(-M,a_priori);toc;

% Lengthen the reconstructed phases to a full-sized vector (see
% relax.m for details)
phi = G*phi ;
% Calculate the reconstructed function
f_rec = reconstruct(A,b,phi,a_priori) ;

% Refine the result with the Gerchberg-Saxton algorithm
proj = create_projector(A,a_priori,string_filters) ;
f_rec = gerchberg_saxton(f_rec,A,b,proj,a_priori,nb_its_gs) ;

x_PC = f_rec;

