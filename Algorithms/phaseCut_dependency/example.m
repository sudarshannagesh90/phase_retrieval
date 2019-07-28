clear all ;

N = 2^6 ; % Size of the signal
nb_meas = 4*N ; % Number of measurements
nb_its_gs = 3000 ;

string_filters = 'Gaussian_measurements' ; % Name of the filters
string_function = 'Gaussian_noise' ; % Name of the test function
                                     % See files
                                     % 'create_function.m' and
                                     % 'create_A_matrix' for a list
                                     % of admissible names.

% Define the function to reconstruct
f0 = create_function (N,string_function) ;
% Define the matrix associated with the filters
A = create_A_matrix (N,nb_meas,string_filters) ;

% Additional information about the signal to reconstruct
a_priori.real = true ; % The signal is real
a_priori.sym = false ; % The signal is symmetric
a_priori.known_indexes = [] ; % Indexes where the phases are known
                              % For wavelet transform / random
                              % filters, if the phase corresponding
                              % to the last wavelet / filter is
                              % known, replace [] by :
                              % [nb_meas-N+1:nb_meas]

% Make the signal compatible with the a priori informations
if a_priori.real, f0 = real(f0) ; end ;
if a_priori.sym, f0 = 0.5*(f0+flipud(f0)) ; end ;

% Indexes where the phase is unknown
mod_indexes = [1:nb_meas] ;
mod_indexes(a_priori.known_indexes) = [] ;
% Calculate the moduli
b = A*f0 ;
b(mod_indexes) = abs(b(mod_indexes)) ;
    
% Construct the relaxation matrix (see file relax.m for details)
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
% Multiply by a unitary complex
if not(a_priori.real)
    theta = angle(scal(f_rec,f0)/scal(f0,f0)) ;
    f_rec = f_rec * exp(i*theta) ;
else
    if (scal(f_rec,f0)<scal(-f_rec,f0))
        f_rec = -f_rec ;
    end
end
    
% Display the results
diff_f = abs(f0-f_rec) ;
diff_f_norm = sqrt(scal(diff_f,diff_f)./scal(f0,f0)) ;
diff_b = abs(abs(b)-abs(A*f_rec)) ;
diff_b_norm = sqrt(scal(diff_b,diff_b)./scal(b,b)) ;
eigenvalues = sort(eigs(X,2,'lm',struct('disp',0))) ;
fprintf(['Second over first eigenvalue : %f\n' ...
         'Differences in L2 norm :\n' ...
         '  Between initial and reconstructed functions : %f\n' ...
         '  Between initial and reconstructed moduli : %f\n'], ...
        eigenvalues(1)/eigenvalues(2), diff_f_norm, diff_b_norm) ;

% Initial and reconstructed functions
figure (1) ; clf ; hold on ;
plot (real(f0),'b-o') ;
plot (real(f_rec),'r') ;
legend ('Initial function (real part)', ...
        'Reconstructed function (real part)') ;

% If of interest, initial and reconstructed moduli
if (strcmp(string_filters,'Cauchy_wavelets')) ...
    || (strcmp(string_filters,'Cauchy_wavelets_non_circ')) ...
    || (strcmp(string_filters,'Cauchy_wavelets_oversampled')) ...
    || (strcmp(string_filters,'Gaussian_random_filters')) ...
    || (strcmp(string_filters,'Binary_random_filters')) ...
    
    figure (2) ; clf ;
    b_init = A*f0 ;
    b_rec = A*f_rec ;
    b_diff = abs(b_init)-abs(b_rec) ;
    if (strcmp(string_filters,'Cauchy_wavelets_oversampled'))
        nb_filters = nb_meas/(2*N) ;
    else
        nb_filters = nb_meas/N ;
    end
    filters_size = nb_meas/nb_filters ;
    for k=1:nb_filters
        subplot (nb_filters,1,k) ; hold on ;
        to_disp = [abs(b_init((k-1)*filters_size+1:k*filters_size)), ...
                   abs(b_rec((k-1)*filters_size+1:k*filters_size)), ...
                   abs(b_diff((k-1)*filters_size+1:k*filters_size))] ;
        m = max(max(to_disp)) ;
        axis([0,filters_size,0,m*1.1]) ;
        plot (to_disp(:,1),'b-o') ;
        plot (to_disp(:,2),'r') ;
        plot (to_disp(:,3),'g--') ;
        if (k==1)
            legend ('Initial moduli', ...
                    'Reconstructed moduli', 'Difference') ;
        end
    end
        
end