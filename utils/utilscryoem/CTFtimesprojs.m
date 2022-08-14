function z = CTFtimesprojs(projs, ctfs, defocusID)
tmpa = projs;
tmpb = fft2(ifftshift(tmpa)); 
tmpb = fftshift(tmpb) .* ctfs(:,:, defocusID); 
z    = fftshift(ifft2(ifftshift(tmpb)));
%z    = real(z); 