function f0 = create_function(N,string)

% Generates a test signal of size N

if strcmp(string,'Gaussian_noise')

    f0 = randn(N,1) + i*randn(N,1) ;

elseif strcmp(string,'Uniform_noise')
    
    f0 = 2 * (rand(N,1) + i*rand(N,1)) - (1 + i) ;
    
elseif strcmp(string,'Sinusoids')
    
    f0_fourier = zeros(N,1) ;
    while (sum(abs(f0_fourier))<1)
        f0_fourier = rand(N,1)<0.05 ;    
        f0 = ifft(f0_fourier.*exp(2*pi*i*rand(N,1))) ;
    end

elseif strcmp(string,'Scan-line')
    
    f0 = scanline(N) + i*scanline(N) ;
    
else
    
    error(['create_function : ', string , ' is not a valid option.']) ;

end