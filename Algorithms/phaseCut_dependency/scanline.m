function f = scanline (N)

folder = dir ('images') ;
images_number = size(folder,1)-2 ;
title = folder(3+floor(images_number*rand(1))).name ;

im = imread(['images/',title]) ;
im = sum(im,3)/765 ;
[m1,m2] = size(im) ;

f1 = im(floor(m1*rand(1))+1,:) ;

if (m2>=N)
    %f = f1(round(m2*[1:N]/N)) ;
    f1_f = fft(f1) ;
    f1_f = [f1_f(1:N-round(N/2)),f1_f(m2-round(N/2)+1:m2)] ;
    f = real(ifft(f1_f))' ;
else
    f1_f = fft(f1) ;
    f1_f = [f1_f(1:round(m2/2)),zeros(N-m2),f1_f(round(m2/2+1):m2)] ;
    f = real(ifft(f1_f))' ;
end