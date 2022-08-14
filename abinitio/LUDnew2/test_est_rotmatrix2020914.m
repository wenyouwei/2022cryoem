
clear all; close all; clc;

%% Generate simulated projections
% Generate 200 simulated projections of size 65x65.
% For simplicity, the projections are centered.
load cleanrib;
figure;isosurface(real(volref),max(real(volref(:)))/5);title('OR');
k = size(volref,1);
n = 129;  %65
V = zeros(n,n,n);
V(1:k,1:k,1:k)= volref;
volref = V;
K   = 500;
SNR = 1/16;
%SNR=1000; % No noise

a              = qrand(K);    %
rotmatrices    = quat2rotm(a);
rotmatricesinv = permute(rotmatrices, [2 1 3]);

A     = OpNufft3D(rotmatrices,n); % projection operator
projs = A * volref;      % projected images

[noisy_projs, sigma] = ProjAddNoise(projs, SNR); 

ref_rot = rotmatrices; 


%[projs,noisy_projs,~,ref_rot]=cryo_gen_projections(n,K,SNR);
figure;viewstack(noisy_projs,5,5); % Show some noisy projections
masked_r = 45;
masked_projs=mask_fuzzy(noisy_projs,masked_r); % Applly circular mask
figure;viewstack(masked_projs,5,5); % Show some noisy projections

% Compute polar Fourier transform, using radial resolution n_r and angular
% resolution n_theta. n_theta is the same as above.
n_theta = 360; %360;%72
n_r = 100;     %100;    %33
[npf,sampling_freqs]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections   

% Find common lines from projections
max_shift=0;
shift_step=1;
common_lines_matrix = commonlines_gaussian(npf,max_shift,shift_step);
C = clstack2C( common_lines_matrix,n_theta );
% Find reference common lines and compare
% [ref_clstack,~]=clmatrix_cheat_qq(ref_rot,n_theta);
[ref_clstack,~]=clmatrix_cheat(ref_rot,n_theta);
p = comparecl( common_lines_matrix, ref_clstack, n_theta, 10 );
fprintf('Percentage of correct common lines: %f%%\n\n',p*100);

%% Methods for est_rotmatrix 
Err_v = zeros(100,1); 
% 1-Model p =2, q =2 
k = 1;
tic;
est_rots = R2_PEG_p2q2(C,n_theta, ref_rot);
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
% fprintf('Reconstruction projections  from clean centered projections \n');
% [ Vol_1, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(Vol_1),max(real(Vol_1(:)))/5);title('R2PEGp2q2');
% Err_v(k) = norm(volref(:) - real(Vol_1(:)))/norm(volref(:));

k = 2;
tic;
est_rots = R2_BLSPG_p2q2(C,n_theta, ref_rot);
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
% fprintf('Reconstruction projections  from clean centered projections \n');
% [ Vol_2, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(Vol_2),max(real(Vol_2(:)))/4);title('R2BLSPGp2q2');
% Err_v(k) = norm(volref(:) - real(Vol_2(:)))/norm(volref(:));


k = 3;
W = ones(2*K);
tic;
for j = 1:11
    est_rots = R2_PGM_p2q2w(W, C,n_theta, ref_rot);
    [W, res] = W_weights(est_rots, C);
     MSEss(j) = check_MSE(est_rots, ref_rot);
end
Time(k) = toc;
figure; plot(1:j,MSEss);title('MSEs');
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
% fprintf('Reconstruction projections  from clean centered projections \n');
% [ Vol_3, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(Vol_3),max(real(Vol_3(:)))/5);title('R2PGMp2q2w');
% Err_v(k) = norm(volref(:) - real(Vol_3(:)))/norm(volref(:));

k = 4;
tic;
est_rots = R3_PGM_p2q2w(C, n_theta, ref_rot);
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
% fprintf('Reconstruction projections  from clean centered projections \n');
% [ Vol_4, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(Vol_4),max(real(Vol_4(:)))/5);title('R3PGMp2q2w');
% Err_v(k) = norm(volref(:) - real(Vol_4(:)))/norm(volref(:));

k = 5;
tic;
est_rots = R2_PEG_p1q1(C,n_theta,ref_rot);
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
fprintf('Reconstruction projections  from clean centered projections \n');
% [ Vol_5, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(Vol_5),max(real(Vol_5(:)))/5);title('R2PEGp1q1');
% Err_v(k) = norm(volref(:) - real(Vol_5(:)))/norm(volref(:)); 
 
fprintf('SNR = %f, pp = %f\n',  SNR, p);
fprintf('K = %d, L = %d, pp = %f\n', K, n_theta, p);
fprintf('-----------------------------------------\n')
fprintf(' Exp  Method         MSE            Time  Err_v \n')
fprintf('  1   R2_PEG_p2q2    %1.5f    %6.2f    %6.2f\n',   MSEs(1),  Time(1), Err_v(1));
fprintf('  2   R2_BLSPG_p2q2  %1.5f    %6.2f    %6.2f\n',   MSEs(2),  Time(2), Err_v(2));
fprintf('  3   R2_PGM_p2q2w   %1.5f    %6.2f    %6.2f\n',   MSEs(3),  Time(3), Err_v(3));
fprintf('  4   R3_PGM_p2q2w   %1.5f    %6.2f    %6.2f\n',   MSEs(4),  Time(4), Err_v(4));
fprintf('  5   R2_PEG_p1q1    %1.5f    %6.2f    %6.2f\n',   MSEs(5),  Time(5), Err_v(5));
fprintf('-------------------------------------------\n');


%% 2-Model p = 2, q = 2;
k = 1;
tic;
est_rots = R_PGM_p2q2(C,n_theta, ref_rot);
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot); 
fprintf('Reconstruction projections  from clean centered projections \n');
[ Vol_1, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(Vol_1),max(real(Vol_1(:)))/5);title('RPGMp2q2');
Err_v(k) = norm(volref(:) - real(Vol_1(:)))/norm(volref(:));

k = 2;
tic;
est_rots = R_APGM_p2q2(C,n_theta, ref_rot);
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
fprintf('Reconstruction projections  from clean centered projections \n');
[ Vol_2, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(Vol_2),max(real(Vol_2(:)))/5);title('RAPGMp2q2');
Err_v(k) = norm(volref(:) - real(Vol_2(:)))/norm(volref(:));

k=3;
tic;
est_rots = R_BLSPG_p2q2(C,n_theta, ref_rot);
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
fprintf('Reconstruction projections  from clean centered projections \n');
[ Vol_3, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(Vol_3),max(real(Vol_3(:)))/5);title('RBLSPGp2q2');
Err_v(k) = norm(volref(:) - real(Vol_3(:)))/norm(volref(:));

k = 4;
tic;
est_rots = R_PEGM_p2q2(C,n_theta, ref_rot);
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
fprintf('Reconstruction projections  from clean centered projections \n');
[ Vol_4, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(Vol_4),max(real(Vol_4(:)))/5);title('RPEGMp2q2');
Err_v(k) = norm(volref(:) - real(Vol_4(:)))/norm(volref(:));

%fprintf('SNR = %f, pp = %f\n',  SNR, p);
fprintf('SNR = %f,K = %d, L = %d, pp = %f\n',SNR, K, n_theta, p);
fprintf('-----------------------------------------\n')
fprintf(' Exp  Method         MSE      Time   Err_v\n')
fprintf('  1   R_PGM_p2q2     %1.5f   %6.2f   %6.2f\n',   MSEs(1), Time(1),Err_v(1));
fprintf('  2   R_APGM_p2q2    %1.5f   %6.2f   %6.2f\n',   MSEs(2), Time(2),Err_v(1));
fprintf('  3   R_BLSPG_p2q2   %1.5f   %6.2f   %6.2f\n',   MSEs(3), Time(3),Err_v(1));
fprintf('  4   R_PEGM_p2q2    %1.5f   %6.2f   %6.2f\n',   MSEs(4), Time(4),Err_v(1));
fprintf('-------------------------------------------\n')


%%  Model p = 2, q = 1/2

k = 1;
W = ones(2*K);
tic;
for j = 1:5
    est_rots = R_PGM_p2q2w(W, C,n_theta, ref_rot);
    [W, res] = W_weights(est_rots, C);
end
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
fprintf('Reconstruction projections  from clean centered projections \n');
[ Vol_1, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(Vol_1),max(real(Vol_1(:)))/5);title('RPGMp2q2w');
Err_v(k) = norm(volref(:) - real(Vol_1(:)))/norm(volref(:));


k = 2;
W = ones(2*K);
tic;
for j = 1:5
    est_rots = R_APGM_p2q2w(W,C,n_theta, ref_rot);
    [W, res] = W_weights(est_rots, C);
end
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
fprintf('Reconstruction projections  from clean centered projections \n');
[ Vol_2, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(Vol_2),max(real(Vol_2(:)))/5);title('RAPGMp2q2w');
Err_v(k) = norm(volref(:) - real(Vol_2(:)))/norm(volref(:));

k=3;
W = ones(2*K);
tic;
for j = 1:5
    est_rots = R_BLSPG_p2q2w(W, C,n_theta, ref_rot);
    [W, res] = W_weights(est_rots, C);
end
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
fprintf('Reconstruction projections  from clean centered projections \n');
[ Vol_3, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(Vol_3),max(real(Vol_3(:)))/5);title('RBLSPGp2q2w');
Err_v(k) = norm(volref(:) - real(Vol_2(:)))/norm(volref(:));


k = 4;
W = ones(2*K);
tic;
W = ones(2*K);
tic;
for j = 1:5
    est_rots = R_PEGM_p2q2w(W,C,n_theta, ref_rot);
    [W, res] = W_weights(est_rots, C);
end
Time(k) = toc;
[MSEs(k),est_inv_rots, err]= check_MSE(est_rots, ref_rot);
fprintf('Reconstruction projections  from clean centered projections \n');
[ Vol_4, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(Vol_4),max(real(Vol_4(:)))/5);title('RPEGMp2q2w');
Err_v(k) = norm(volref(:) - real(Vol_4(:)))/norm(volref(:));


%%  Model p = 2, q = 1/2
fprintf('SNR = %f, pp = %f\n',  SNR, p);
fprintf('K = %d, L = %d, pp = %f\n', K, n_theta, p);
fprintf('-----------------------------------------\n')
fprintf(' Exp  Method         MSE         Time     Err_v\n')
fprintf('  1   R_PGM_p2q2w    %1.5f   %6.2f   %6.2f\n',   MSEs(1), Time(1),Err_v(1));
fprintf('  2   R_APGM_p2q2w   %1.5f   %6.2f   %6.2f\n',   MSEs(2), Time(2),Err_v(1));
fprintf('  3   R_BLSPG_p2q2w  %1.5f   %6.2f   %6.2f\n',   MSEs(3), Time(3),Err_v(1));
fprintf('  4   R_PEGM_p2q2w   %1.5f   %6.2f   %6.2f\n',   MSEs(4), Time(4),Err_v(1));
fprintf('-------------------------------------------\n')


%% Model p = 2, q = 1;
k = 1;
tic;
est_inv_rots = R_PDM_p2q1(C,n_theta, ref_rot);
Time(k) = toc;
MSE(k) = check_MSE(est_inv_rots, ref_rot);


k = 2;
tic;
est_inv_rots = R_PEGM_p2q1(C,n_theta, ref_rot);
Time(k) = toc;
MSE(k) = check_MSE(est_inv_rots, ref_rot);

k = 3;
tic;
est_inv_rots = R_PGBLS_p2q1(C,n_theta, ref_rot);
Time(k) = toc;
MSE(k) = check_MSE(est_inv_rots, ref_rot);


fprintf('SNR = %f, pp = %f\n',  SNR, p);
fprintf('K = %d, L = %d, pp= %f \n', K, n_theta,p);
fprintf('-----------------------------------------\n')
fprintf(' Exp  Method          MSE         Time\n')
fprintf('  1   R_PDM_p2q1     %1.5f      %6.2f\n',   MSE(1),  Time(1));
fprintf('  2   R_PEGM_p2q1    %1.5f      %6.2f\n',   MSE(2),  Time(2));
fprintf('  3   R_PGBLS_p2q1   %1.5f      %6.2f\n',   MSE(3),  Time(3));
fprintf('-------------------------------------------\n')









%% Test est_orientations_LS + SDPLR
% tic;
% est_inv_rots = est_orientations_LS(common_lines_matrix, n_theta);
% Time(1) = toc;
% MSEs(1) = check_MSE(est_inv_rots, ref_rot);
% fprintf('Reconstruction projections  from clean centered projections \n')
% 
% [ v_LS, v_b, kernel ,err, iter, flag] = recon3d_firm(masked_projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(v_LS),max(real(v_LS(:)))/5);


%% With spectral norm constraint + ADMM
k = 1;
pars.alpha = 2/3;
tic;
est_inv_rots = est_orientations_LS(common_lines_matrix, n_theta, pars,ref_rot);
Time(k) = toc;
[MSEs(k), est_inv_rots]= check_MSE(est_inv_rots, ref_rot);
[ v_LSc, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(v_LSc),max(real(v_LSc(:)))/5);title('LSc');
Err_v(k) = norm(volref(:) - real(v_LSc(:)))/norm(volref(:));

%% Test est_orientations_LUD
k=2;
tic;
est_inv_rots = est_orientations_LUD(common_lines_matrix,n_theta);
Time(k) = toc;
[MSEs(k), est_inv_rots]= check_MSE(est_inv_rots, ref_rot);
[ v_LUDADMM, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(v_LUDADMM),max(real(v_LUDADMM(:)))/5);title('ADMM');
Err_v(k) = norm(volref(:) - real(v_LUDADMM(:)))/norm(volref(:));

% ADMM with spectral norm constraint
k = 3;
pars.alpha = 2/3;
pars.solver = 'ADMM';
tic;
est_inv_rots = est_orientations_LUD(common_lines_matrix,n_theta, pars);
Time(k) = toc;
[MSEs(k), est_inv_rots] = check_MSE(est_inv_rots, ref_rot);
[ v_LUDADMMc, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(v_LUDADMMc),max(real(v_LUDADMMc(:)))/5);title('ADMMc');
Err_v(k)= norm(volref(:) - real(v_LUDADMMc(:)))/norm(volref(:));

%% IRLS
%  pars.solver = 'IRLS';
% pars.alpha = 0;
% tic;
% est_inv_rots = est_orientations_LUD(common_lines_matrix,n_theta, pars);
% Time(5) = toc;
% [MSEs(5), est_inv_rots] = check_MSE(est_inv_rots, ref_rot);
% [ v_IRLS, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(v_IRLS),max(real(v_IRLS(:)))/5);
% Err_v(5) = norm(volref(:) - real(v_IRLS(:)))/norm(volref(:))
% MSEs(5)

% IRLS with spectral norm constraint
k= 4;
pars.solver = 'IRLS';
pars.alpha = 2/3;
tic;
est_inv_rots = est_orientations_LUD(common_lines_matrix,n_theta, pars);
Time(k) = toc;
[MSEs(k), est_inv_rots] = check_MSE(est_inv_rots, ref_rot);
[ v_IRLSc, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(v_IRLSc),max(real(v_IRLSc(:)))/4);title('IRLSc');
Err_v(k) = norm(volref(:) - real(v_IRLSc(:)))/norm(volref(:));




%% Print the MSEs and cost time of the results

fprintf('SNR = %f, pp = %f\n',  SNR, p);
fprintf('K = %d, n_r = %d, L = %d,masked_r = %d \n', K,n_r, n_theta,masked_r);
fprintf('-----------------------------------------\n')
fprintf('  1   LSc     %1.5f    %1.5f  %6.2f\n',   MSEs(1), Err_v(1), Time(1));
fprintf('  2   LUD     %1.5f    %1.5f  %6.2f\n',   MSEs(2), Err_v(2), Time(2));
fprintf('  3   LUDc    %1.5f    %1.5f  %6.2f\n',   MSEs(3), Err_v(3), Time(3));
fprintf('  4  IRLSc    %1.5f    %1.5f  %6.2f\n',    MSEs(4),Err_v(4), Time(4));
fprintf('-------------------------------------------\n')




% [FSC_PD, spatialFrequency, meanIntensity] = FourierShellCorrelate(v_PD,volref,n,n);
% [FSC_PDc, spatialFrequency, meanIntensity] = FourierShellCorrelate(v_PDc,volref,n,n);
% [FSC_LSc, spatialFrequency, meanIntensity] = FourierShellCorrelate(v_LSc,volref,n,n);
% [FSC_ADMM, spatialFrequency, meanIntensity] = FourierShellCorrelate(v_LUDADMM,volref,n,n);
% [FSC_ADMMc, spatialFrequency, meanIntensity] = FourierShellCorrelate(v_LUDADMMc,volref,n,n);
% [FSC_IRLSc, spatialFrequency, meanIntensity] = FourierShellCorrelate(v_IRLSc,volref,n,n);
% [FSC_eig, spatialFrequency, meanIntensity] = FourierShellCorrelate(v_eig,volref,n,n);
% 
% 
% figure;    
% tick_spacing = 1;     
% plot(FSC_PD,'b','LineWidth',2);hold on     
% plot(FSC_LSc,'c--','LineWidth',2);hold on ;   
% plot(FSC_ADMM,'g','LineWidth',2);hold on ;  
% plot(FSC_ADMMc,'b--','LineWidth',2);  hold on ;
% plot(FSC_IRLSc,'r--','LineWidth',2);  hold on ;
% plot(FSC_eig,'y','LineWidth',2);  
% ylim([0 0.1]); %xlim([1 n/2]) 
% set(gca,'YTick',[0,0.143,0.5,1])
% f_ind_s = num2str([0:0.01:0.05]');    
% set(gca,'XTickLabel',f_ind_s);   
% ylabel('FSC');   
% xlabel('Spatial frequency (1/pixel size)')   
% grid on   
% legend('firm','Tikfix','Tik','TV'); set(gca, 'FontSize',18);  
% % print('-dpng',strcat('emd8523FSC',num2str(1/SNR),num2str(-K)));
% 
% 
% 






