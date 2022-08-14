clear; close all; clc;
load('gaussian_ph.mat');
load p500;

K = 1; 
%q = qrand(K);
%proja  = cryo_project(ph,q);
%R      = reshape(quat2rotm(q),9,K);
R = inv_rot_matrices(:,:,1:K); proja = projections(:,:,1:K)'; 
% R = reshape(inv(R'),9,K); 
% projs  = CryoForwardProj(ph,R,64,'single'); 
% norm(projs(:)-proja(:))/norm(proja(:))
[v, v_b, kernel ] = recon3d_firm(proja,R,[], 1e-4, 1, zeros(64,64,64));

return; 
q2 = rotm2quat(reshape(R,3,3,K))';
tic; proja = cryo_project(ph,q2); toc; 
norm(projs(:)-proja(:))/norm(projs(:))

% R2 = reshape(R,3,3,K); 
% [v, v_b, kernel ] = recon3d_firm(projs,R2,[], 1e-4, 500, zeros(64,64,64));
% projs_fourier = FFT_projections( projs);
% v_b2          = nufft_adj(projs_fourier(idxcub), sta); %toc;
% norm(v_b(:)-v_b2(:))/norm(v_b(:))
% z = CryoEM_Atimesx(kernel,ph);
% norm(v_b(:)-z(:))/norm(z(:))


% q2 = rotm2quat(reshape(R,3,3,K))';
% norm(q2-q)

for i = 1:K
    tmp = (quat2rotm(q(:,i)))';    
    R(:,i) = tmp(:);
end
RotMatrixinv    = reshape(R, 9, K);

maskR = size(ph,1)/2-1; 
tic; projaa = OpForwardProj(ph, R, maskR); toc
norm(projaa(:)-proja(:))/norm(proja(:))



return;
%% %%
load p500;
K = 2; 
R = inv_rot_matrices(:,:,1:K); 
projaa = projections(:,:,1:K); 
q2 = rotm2quat(reshape(R,3,3,K))';
tic; proja = cryo_project(ph,q2); toc; 
norm(projaa(:)-proja(:))/norm(projaa(:))


for i = 1:K
    tmp2  = inv_rot_matrices(:,:,i);  %  
    q(:,i) = rotm2quat(tmp2)';
end
proj2 = cryo_project(ph,q); 
norm(projaa(:)-proj2(:))/norm(projaa(:))
