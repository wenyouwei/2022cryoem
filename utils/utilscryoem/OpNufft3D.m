% The code is based on the following source codes
%% precomp.m by Lanhui Wang(2012).
%% cryo_projection by Gene Katsevich and Yoel Shkolnisky(2014).
%% Written by Prof.You-Wen Wen (wenyouwei@gmail.com) on April 2018. 

classdef OpNufft3D
    properties
        ndim
        projnumber
        omega
        idxcub
        transpose
        range
        st
    end
    
    methods
        %% construct the class MyOpNuff3D
        function obj = OpNufft3D(rotmatrix,n)
            invrotmatrix   = permute(rotmatrix,[2,1,3]);
            obj.transpose  = 0; 
            obj.projnumber = size(rotmatrix,3);
            %obj.rotmatrix  = rotmatrix; % rotation matrix
            obj.ndim       = n; 
            
            %range    = -(n-1)/2:(n-1)/2;    % grid of polar coordinates 
            if  mod(n,2)==0
                range = -fix(n/2):fix(n/2)-1;
            else
                range = -(n-1)/2:(n-1)/2;    % grid of polar coordinates
            end
            
            n_x(:,:) = invrotmatrix(:,1,:); 
            n_y(:,:) = invrotmatrix(:,2,:);
            %[I,J]    = ndgrid(range,range); % grid of cartesian coordinates
            [I,J] = meshgrid(range,range);% grid of cartesian coordinates
            I = I(:);            J = J(:);
            omega = zeros(n^2*obj.projnumber,3);
            for k=1:obj.projnumber
                P = I * n_x(:,k)' + J * n_y(:,k)'; % polar --> cartesian
                omega((k-1)*n^2+1:k*n^2,:) = P; % projection coordinates
            end
            % points in the Fourier cube
            idxcub = max(omega,[],2)<fix(n/2) & min(omega,[],2)>-fix(n/2);
            
            obj.idxcub = idxcub; 
            obj.omega  = 2 * pi /n * omega; 
            obj.range  = range; 
                        
            Jd = [4 4 4];
            Nd = [n n n];
            Ld = [2^11 2^11 2^11];  % use table for large matrix
            Kd = 2 * Nd;
            n_shift = fix(Nd/2); % stress it

            ktype = 'minmax:kb';
            st = nufft_init(obj.omega(idxcub,:), Nd, Jd, Kd, n_shift, 'table', Ld, ktype);
            obj.st = st; 
        end
        
        %% 
        function [y,pfs] = mtimes(A,x)
            n = size(x,1);
            P = A.omega; 
            %% compute A' x ==>back projection
            if A.transpose == 1  % x is the projected image                  
                xfft = fftshift(fft2(ifftshift(x))); % spatial --> frequency  
                pfs  = xfft(A.idxcub); 
                % pfs: images' fourier coefficient  A.st: direction of
                % projections. 
                y    = nufft_adj(pfs, A.st)/n^2; 
                %y    = anufft3(pfs, P', [n n n]);
                y    = real(y); 

                return;
            else % x is the volume, A.transpose == 0
                %% compute A x ==> forward projection
                proj_fourier = zeros(size(A.idxcub,1),1); 
                proj_fourier(A.idxcub) = nufft(x,A.st);
                
                %proj_fourier = reshape(proj_fourier, n, n, A.projnumber);   
                %proj_fourier = ifftshift(proj_fourier);
                %y   = fftshift(ifft2(proj_fourier));
                
                
                if mod(n,2)==0
                    [I,J]     = meshgrid(A.range,A.range); %
                    %[I,J]     = ndgrid(A.range,A.range); 
                    IandJ  = I(:) + J(:);
                    tmp = exp(2*pi/(2*n)*1i*(IandJ-1));
                    tmp = repmat(tmp,[A.projnumber 1]);   
                    proj_fourier = proj_fourier.*exp(1i.*sum(-P,2)/2);
                    proj_fourier = proj_fourier .* tmp; 
                end
            
                proj_fourier = reshape(proj_fourier, n, n, A.projnumber);   
                proj_fourier = ifftshift(proj_fourier);
                y   = fftshift(ifft2(proj_fourier));
 
                if mod(n,2)==0
                    tmp = reshape(exp(2*pi/(2*n)*1i*(IandJ)),n,n); 
                    tmp = repmat(tmp,[1 1 A.projnumber]); 
                    y = y .* tmp;
                end
                y = real(y); 
%                 if isa(x,'single')
%                     imagtol=5.0e-7;
%                 elseif isa(x,'double')
%                     imagtol=5.0e-13;
%                 end
%                         
%                 if norm(imag(projection(:)))/norm(projection(:))>imagtol
%                     error('GCAR:imaginaryComponents','projection has imag components');
%                 end
%                 y = real(projection);
            end
        end
        
        %% 
        function y = ctranspose(A)
            y = A; 
            if A.transpose == 0 
                y.transpose = 1; 
            else
                y.transpose = 0; 
            end
        end
        
        %% 
        function z = compkernel(A,ctf,defocusID)
            switch nargin 
                case 1
                    weights = ones(A.st.M,1); 
                case 2 
                    weights = ctf.^2; 
                    weights = weights(A.idxcub);
                case 3
                    weights = ctf(:,:,defocusID).^2; 
                    weights = weights(A.idxcub);
            end
 
            n     = A.ndim; 
            P     = A.omega; 
            P     = P(A.idxcub,:); 
            idx   = 2*n:-1:2; 
            toepker = zeros(2*n,2*n,2*n); 
            %myshift = [-n -n 0;-n 0 0; 0 -n 0; 0 0 0]' + fix((n)/2); 
            myshift = [n n 0;n 0 0; 0 n 0; 0 0 0]' ;
            
            st1     = A.st;   st2 = A.st; st3 = A.st; st4 = A.st; 
            st1.phase_shift   = exp(1i*P*myshift(:,1));
            st2.phase_shift   = exp(1i*P*myshift(:,2));
            st3.phase_shift   = exp(1i*P*myshift(:,3));
            st4.phase_shift   = exp(1i*P*myshift(:,4));
                      
            toepker(1:n,1:n,n+1:2*n)         = nufft_adj(weights, st1);
            toepker(1:n,n+1:2*n,n+1:2*n)     = nufft_adj(weights, st2);
            toepker(n+1:2*n,1:n,n+1:2*n)     = nufft_adj(weights, st3);
            toepker(n+1:2*n,n+1:2*n,n+1:2*n) = nufft_adj(weights, st4);
            toepker(2:end,2:end,2:n)         = conj(toepker(idx,idx,2*n:-1:2+n)); 
            z = toepker(2:end,2:end,2:end)/n^2;
        end
        
        function B = AtCtCA(A,C)
            B = compkernel(A,C.ctfs,C.defocusID);
        end
        
    %%   methods end  
    end
end
