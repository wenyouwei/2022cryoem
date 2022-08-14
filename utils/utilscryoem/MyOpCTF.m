% The code is based on the following source codes
%% precomp.m by Lanhui Wang(2012).
%% cryo_projection by Gene Katsevich and Yoel Shkolnisky(2014).
%% Written by Prof.You-Wen Wen (wenyouwei@gmail.com) on April 2018. 

classdef MyOpCTF
    properties
        %rotmatrix
        %invrotmatrix
        ndim
        NOproj
        ctfs
        defocusID
    end
    
    methods
        %% construct the class MyOpCTF
        function obj = MyOpCTF(n, K, defocus)
            if nargin == 2, defocus = [1.4 1.75 2]; end
            alpha  = 0.07;    res = 3.36;  
            lambda = 0.0251;  Cs  = 2;      B = 100;
            %if mod(n,2)==0, siz=n+1; else siz = n; end
            siz = n + 1; 
            n_d  = length(defocus);  
            ctfs = zeros(n,n,n_d);

            for k = 1:n_d
                ctf = CTF(siz, res, lambda, defocus(k), Cs, B, alpha);
                ctfs(:,:,k) = ctf(1:end-1,1:end-1);
            end
            
            obj.ctfs    = ctfs; 
            obj.NOproj  = K;
            obj.ndim    = n; 
            obj.defocusID = mod(0:K-1,n_d) + 1;
        end
        
        %% 
        function y = mtimes(A,x)
            tmpa = x;
            tmpb = fft2(ifftshift(tmpa)); 
            tmpb = fftshift(tmpb) .* A.ctfs(:,:, A.defocusID); 
            z    = fftshift(ifft2(ifftshift(tmpb)));
            y    = z; 
            %y    = real(z); 
        end
        
    %%   methods end  
    end
end
