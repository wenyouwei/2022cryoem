clear; close all; clc;
load('gaussian_ph.mat');
load p500;

n = size(ph,1); K = 1; 

% test the paper's data
R = inv_rot_matrices(:,:,1:K); 
projs  = projections(:,:,1:K); 
projs  = permute(projs,[2 1 3]);

maskR  = size(ph,1)/2-1; 
R2 = R(:); ph = permute(ph, [ 3 1 2]); 
projaa = CryoForwardProj(ph, reshape(R2,9,K),n,'single'); 
norm(projaa(:)-projs(:))/norm(projs(:))
R3 = R'; R3 = R3(:); ph = permute(ph, [ 3 1 2]); 
projaa = CryoForwardProj(ph, reshape(R3,9,K),n,'single'); 
norm(projaa(:)-projs(:))/norm(projs(:))
R3 = inv(R); R3 = R3(:); ph = permute(ph, [ 3 1 2]); 
projaa = CryoForwardProj(ph, reshape(R3,9,K),n,'single'); 
norm(projaa(:)-projs(:))/norm(projs(:))
R3 = inv(R'); R3 = R3(:); ph = permute(ph, [ 3 1 2]); 
projaa = CryoForwardProj(ph, reshape(R3,9,K),n,'single'); 
norm(projaa(:)-projs(:))/norm(projs(:))


return; 

fprintf('Reconstruction from clean centered projections \n')
[va, v_b, kernel ,err, iter, flag] = recon3d_firm(projections,inv_rot_matrices,[], 1e-4, 500, zeros(64,64,64));
fprintf('The relative error of reconstruction is %f.\n',norm(va(:)-ph(:))/norm(ph(:)));

[At_b, st] = OpBackProjNUfft(projections, inv_rot_matrices);
err1 =  norm(At_b(:)-v_b(:)); 
disp('*************************************************')
fprintf('check subprogram nufft:, e err=%1.2e\n', err1);

z = CryoEM_Atimesx(kernel,ph);
err2  = norm(v_b(:)-z(:))/norm(z(:)); 
fprintf('check kernel:, err=%1.2e\n', err2);


disp('*************************************************')


disp(''); disp('');
disp('******************************************')
disp('         Nufft vs. interpolation ');
q = qrand(K);
% proja = cryo_project(ph,q);  
R     = quat2rotm(q);
projs = CryoForwardProj(ph, reshape(R,9,K),n,'single'); 

Rintp = zeros(9,K); 
for i = 1:K
    tmp  = inv(R(:,:,i));
    Rintp(:,i) = tmp(:);
end
maskR = size(ph,1)/2-1; 
tic; projaa = OpForwardProj(ph, Rintp, maskR); toc  % interpolation 

disp('***** Forward Projection *******'); 
norm(projaa(:)-projs(:))/norm(projs(:))

[ va, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,R,[], 1e-4, 500, zeros(64,64,64));
[At_b, st] = OpBackProjNUfft(projs, R);
err1 =  norm(At_b(:)-v_b(:)); 
fprintf('check subprogram nufft:, e err=%1.2e\n', err1);


norm(va(:)-ph(:))/norm(va(:))

return; 
vb2 = OpBackProj(projaa, Rintp, maskR);
disp('***** Back Projection *******'); 
norm(vb2(:)-v_b(:))/norm(vb2(:))