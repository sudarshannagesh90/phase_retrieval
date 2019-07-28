function A = create_A_matrix (N,nb_meas,string)

% Generates the matrix associated to the considered linear application.

if (strcmp(string,'Cauchy_wavelets'))
    % Circular convolution with Cauchy wavelets.

    if (mod(nb_meas,N)~=0)
        error (['create_A_matrix : the number of measurements must ' ...
                'be a multiple of the signals size.']) ;
    end

    [phi,psi] = cauchy_wavelets (N,nb_meas/N,5) ;
    filters = [psi,phi] ;

    [x_ind,y_ind] = meshgrid([0:N-1]) ;
    indexes = 1+mod(y_ind-x_ind,N) ;
    A = [] ;
    for j=1:nb_meas/N
        filter = ifft(filters(:,j)) ;
        A = [A;filter(indexes)] ;
    end
    
elseif (strcmp(string,'Cauchy_wavelets_non_circ'))
    
    % Convolution with Cauchy wavelets.
    % It is non circular : the signal is implicitly lengthened by
    % zero to the left and to the right.
    
    if (mod(nb_meas,N)~=0)
        error (['create_A_matrix : the number of measurements must ' ...
                'be a multiple of the signals size.']) ;
    end
    
    [phi,psi] = cauchy_wavelets (N,nb_meas/N,5) ;
    filters = [psi,phi] ;
    
    [x_ind,y_ind] = meshgrid([0:N-1]) ;
    indexes = N+1+y_ind-x_ind ;
    A = [] ;
    for j=1:nb_meas/N
        filter = ifft(filters(:,j)) ;
        double_filter = [zeros(N/2,1) ; filter(N/2+1:N) ; ...
                         filter(1:N/2) ; zeros(N/2,1)] ;
        A = [A ; double_filter(indexes)] ;        
    end
    
elseif (strcmp(string,'Cauchy_wavelets_oversampled'))
    
    % Convolution with Cauchy wavelets, oversampled by a factor 2.
    
    if (mod(nb_meas,2*N)~=0)
        error (['create_A_matrix : the number of measurements must ' ...
                'be a multiple of the double of the signals size.']) ;
    end

    [phi,psi] = cauchy_wavelets (N,nb_meas/(2*N),5) ;
    filters = [[psi;zeros(N,nb_meas/(2*N)-1)], ...
               [phi(1:N/2);zeros(N,1);phi(N/2+1:N)]] ;
    
    [x_ind,y_ind] = meshgrid([0:N-1],[0:2*N-1]) ;
    indexes = mod(y_ind-2*x_ind,2*N)+1 ;
    A = [] ;
    for j=1:nb_meas/(2*N)
        filter = ifft(filters(:,j)) ;
        A = [A ; filter(indexes)] ;
    end
    
elseif (strcmp(string,'Gaussian_random_filters'))
    
    % Fourier transform after multiplication by a gaussian random
    % filter.
    
    if (mod(nb_meas,N)~=0)
        error (['create_A_matrix : the number of measurements must ' ...
                'be a multiple of the signals size.']) ;
    end

    filters = randn(N,nb_meas/N) + i*randn(N,nb_meas/N) ;
    
    A = [] ;
    TF = exp(-2*pi*i/N*[0:(N-1)]'*[0:(N-1)]) ;
    for j=1:nb_meas/N
        F = ones(N,1)*filters(:,j).' ;
        A = [A ; F.*TF] ;
    end
    
elseif (strcmp(string,'Binary_random_filters'))
    
    % Fourier transform after multiplication by a binary random
    % filter.
    % The filters are such that there is no point of the entry
    % signal which every filter multiplies by zero.
        
    if (mod(nb_meas,N)~=0)
        error (['create_A_matrix : the number of measurements must ' ...
                'be a multiple of the signals size.']) ;
    end

    filters = +(rand(N,nb_meas/N)<0.5) ;
    always_zero = find(prod(1-filters,2)) ;
    indexes_to_change = always_zero + ...
        N * floor(nb_meas/N*rand(size(always_zero))) ;
    filters(indexes_to_change) = 1 ;
    
    A = [] ;
    TF = exp(-2*pi*i/N*[0:(N-1)]'*[0:(N-1)]) ;
    for j=1:nb_meas/N
        F = ones(N,1)*filters(:,j).' ;
        A = [A ; F.*TF] ;
    end

elseif (strcmp(string,'Gaussian_measurements'))

    % Gaussian random measurements
    
    A = randn(nb_meas,N) + i*randn(nb_meas,N) ;
    
else
   
    error(['create_A_matrix : ', string , ' is not a valid option.']) ;
 
end