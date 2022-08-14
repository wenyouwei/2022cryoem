clear; close all; clc;
load('gaussian_ph.mat');
load p500;

n = size(ph,1); K = 500; 


% q = qrand(K);
% proja = cryo_project(ph,q);  
% R     = quat2rotm(q);
% inv_rot_matrices = R; 
R = inv_rot_matrices; 

projs = CryoForwardProj(ph, reshape(R,9,K),n,'single'); 
inv_rot_matrices = R; 
projections = projs; 
projections=permute(projections,[2 1 3]);

% fprintf('Reconstruction from clean centered projections \n')
% [ v, v_b, kernel ,err, iter, flag] = recon3d_firm( projections,inv_rot_matrices,[], 1e-6, 1000, zeros(64,64,64));
% fprintf('The relative error of reconstruction is %f.\n',...
%     norm(v(:)-ph(:))/norm(ph(:)));
% 

% for k = 1:K, 
%     tmp = inv_rot_matrices(:, :, k); 
%     tmp = inv(tmp); 
%     RotMatrix(:,k) =  tmp(:); 
% end   
% RotMatrix    = reshape(inv_rot_matrices, 9, K);
% projs  = CryoForwardProj(ph,RotMatrix,n,'single'); 
% projs2 = permute(projs,[2 1 3]);
% norm(projs(:)-projections(:))
% 
% norm(projs2(:)-projections(:))

projectionsold = projections; 
%projections = projections + max(projections(:))/100 * randn(size(projections)); 

fprintf('Reconstruction from clean centered projections \n')
[ v, v_b, kernel ,err, iter, flag,v2] = recon3d_firm(projections,inv_rot_matrices,[], 1e-4, 500, zeros(64,64,64));
fprintf('The relative error of reconstruction is %f.\n',norm(v(:)-ph(:))/norm(ph(:)));



return;