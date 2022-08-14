function [z, sigma] = ProjAddNoise(projs, SNR)

randn('state', 10000); 
K = size(projs,3); 
n = size(projs,1); 
z = projs; 
sigma = zeros(K,1); 
%tmp  = 0;
for k = 1:K
    proj  = projs(:,:,k);
    sigma(k) = sqrt(var(proj(:))/SNR);
    noise    = randn(n);
    
    noise = noise/std(noise(:));
    noise = noise * sigma(k);
    z(:,:,k) = proj + noise;
    %tmpa = [norm(noise(:))^2 n*n*sigma(k)^2]
    %tmp = tmp + norm(noise(:))^2;
end
%tmp

