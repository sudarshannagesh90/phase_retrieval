function proj = create_projector (A,a_priori,string)

% Calculates the orthogonal projection onto range(A).

% A is the linear application
% a_priori contains the informations about the signal
% string is the name of the chosen linear application

[nb_meas,N] = size(A) ;

if (strcmp(string,'Cauchy_wavelets'))
    % Circular convolution with Cauchy wavelets.
    % Complexity : O ( N * log(N) * filters' number )

    [phi,psi] = cauchy_wavelets (N,nb_meas/N,5) ;
    filters = [psi,phi] ;

    somme = sum(abs(filters).^2,2) ;
    if a_priori.real
        somme = somme + [somme(1);flipud(somme(2:N))] ;
        proj = @(hj) (proj_cauchy_wavelets_real(N,nb_meas,filters,somme,hj)) ;
    else
        proj = @(hj) (proj_cauchy_wavelets(N,nb_meas,filters,somme,hj)) ;
    end
            
elseif (strcmp(string,'Cauchy_wavelets_oversampled'))
    
    % Convolution with Cauchy wavelets, oversampled by a factor 2.
    % Complexity : O ( N * log(N) * filters' number )
    
    [phi,psi] = cauchy_wavelets (N,nb_meas/(2*N),5) ;
    filters = [[psi;zeros(N,nb_meas/(2*N)-1)], ...
               [phi(1:N/2);zeros(N,1);phi(N/2+1:N)]] ;

    somme = sum(abs(filters).^2,2) ;
    somme = (somme(1:N)+somme(N+1:2*N))/2 ;
    if a_priori.real
        somme = somme + [somme(1);flipud(somme(2:N))] ;
        proj = @(hj) (proj_cauchy_wavelets_oversampled_real ...
                      (N,nb_meas,filters,somme,hj)) ;
    else
        proj = @(hj) (proj_cauchy_wavelets_oversampled ...
                      (N,nb_meas,filters,somme,hj)) ;
    end
    
elseif ((strcmp(string,'Gaussian_random_filters')) ...
    || (strcmp(string,'Binary_random_filters')))
    
    % Fourier transform after multiplication by a filters
    % Complexity : O ( N * log(N) * filters' number )

    % Get filters back from A
    filters = A(1+N*[0:nb_meas/N-1],:).' ;

    somme = sum(abs(filters).^2,2) ;
    proj = @(hj) ...
           (proj_filters(N,nb_meas,filters,somme,a_priori.real,hj)) ;
    
elseif ((strcmp(string,'Gaussian_measurements')) ...
    || (strcmp(string,'Cauchy_wavelets_non_circ')))

    % It is not possible to use the structure of A
    % Complexity : O ( nb_meas^2 )

    if a_priori.real
        A_pinv = pinv([real(A);imag(A)],1e-10) ;
    else
        A_pinv = pinv(A,1e-10) ;
    end
    proj = @(hj) (general_projector(N,nb_meas,A,A_pinv,a_priori.real,hj)) ...
           ;
    
else
   
    error(['create_A_matrix : ', string , ' is not a valid option.']) ;
 
end



function [hj2,f_rec] = proj_cauchy_wavelets_real(N,nb_meas,filters,somme,hj)
    hj = reshape(hj,N,nb_meas/N) ;
    f_rec_fourier = sum(fft(hj).*conj(filters),2) ;
    f_rec_fourier = f_rec_fourier + ...
        conj([f_rec_fourier(1);flipud(f_rec_fourier(2:N))]) ;
    f_rec_fourier = f_rec_fourier./somme ;
    hj2 = [] ;
    for j=1:size(filters,2)
        hj2 = [hj2;ifft(f_rec_fourier.*filters(:,j))] ;
    end
    f_rec = real(ifft(f_rec_fourier)) ;



function [hj2,f_rec] = proj_cauchy_wavelets(N,nb_meas,filters,somme,hj)
    hj = reshape(hj,N,nb_meas/N) ;
    f_rec_fourier = sum(fft(hj).*conj(filters),2)./somme ;
    hj2 = [] ;
    for j=1:size(filters,2)
        hj2 = [hj2;ifft(f_rec_fourier.*filters(:,j))] ;
    end
    f_rec = ifft(f_rec_fourier) ;



function [hj2,f_rec] = proj_cauchy_wavelets_oversampled(N,nb_meas,filters,somme,hj)
    hj = reshape(hj,2*N,nb_meas/(2*N)) ;
    A_star_hj = sum(fft(hj).*conj(filters),2) ;
    A_star_hj = (A_star_hj(1:N)+A_star_hj(N+1:2*N))/2 ;
    f_rec_fourier = A_star_hj./somme ;
    hj2 = [] ;
    for j=1:size(filters,2)
        hj2 = [hj2;ifft([f_rec_fourier;f_rec_fourier].*filters(:,j))] ;
    end
    f_rec = ifft(f_rec_fourier) ;
    
    
    
function [hj2,f_rec] = proj_cauchy_wavelets_oversampled_real ...
        (N,nb_meas,filters,somme,hj)
    hj = reshape(hj,2*N,nb_meas/(2*N)) ;
    A_star_hj = sum(fft(hj).*conj(filters),2) ;
    A_star_hj = (A_star_hj(1:N)+A_star_hj(N+1:2*N))/2 ;
    A_star_hj = A_star_hj + ...
        conj([A_star_hj(1);flipud(A_star_hj(2:N))]) ;
    f_rec_fourier = A_star_hj./somme ;
    hj2 = [] ;
    for j=1:size(filters,2)
        hj2 = [hj2;ifft([f_rec_fourier;f_rec_fourier].*filters(:,j))] ;
    end
    f_rec = real(ifft(f_rec_fourier)) ;

function [hj2,f_rec] = proj_filters (N,nb_meas,filters,somme, ...
                                     real_bool,hj)
    hj = reshape(hj,N,nb_meas/N) ;
    f_rec = sum(ifft(hj).*conj(filters),2)./somme ;
    if (real_bool)
        f_rec = real(f_rec) ;
    end
    hj2 = [] ;
    for j=1:size(filters,2)
    hj2 = [hj2;fft(f_rec.*filters(:,j))] ; 
    end
    
function [hj2,f_rec] = general_projector(N,nb_meas,A,A_pinv, ...
                                         real_bool,hj)
    if real_bool
        f_rec = A_pinv*[real(hj);imag(hj)] ;
    else
        f_rec = A_pinv * hj ;
    end
    hj2 = A*f_rec ;    