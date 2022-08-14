
clear all; close all; clc;

%% Objection:
%   min || RiCij - RjCji ||_{p}^{q}

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
K = 100;
SNR = 1;
%SNR=1000; % No noise
[projs,noisy_projs,~,ref_rot]=cryo_gen_projections(n,K,SNR);
figure;viewstack(noisy_projs,5,5); % Show some noisy projections
masked_r = 35;
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

% Find reference common lines and compare
% [ref_clstack,~]=clmatrix_cheat_qq(ref_rot,n_theta);
[ref_clstack,~]=clmatrix_cheat(ref_rot,n_theta);
p = comparecl( common_lines_matrix, ref_clstack, n_theta, 10 );
fprintf('Percentage of correct common lines: %f%%\n\n',p*100);


MSEs = zeros(9,1);
Err_v = zeros(9,1);
Time = zeros(9,1);

opts.tol = 1e-3;
% proximal Method 
est_rots = L2_PG(common_lines_matrix,n_theta, ref_rot);
% Accelerated proximal Method 
est_rots = L2_APG(common_lines_matrix,n_theta, ref_rot);

est_rots = L2_BLSPG(common_lines_matrix,n_theta, ref_rot);
[MSEs(1) est_inv_rots] = check_MSE(est_rots, ref_rot);
inv_est_rot = permute(est_rots, [2 1 3]);
[MSEs(2) est_inv_rots]= check_MSE(inv_est_rot, ref_rot);
[ v_PG, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(v_PG),max(real(v_PG(:)))/5);title('APG');
Err_v(1) = norm(volref(:) - real(v_PG(:)))/norm(volref(:));


pars.solver = 'PDM1';
tic;
% est_inv_rots = cryoEMPDM(common_lines_matrix,n_theta,opts);
[est_inv_rots ]= cryoEMPDMM(common_lines_matrix,n_theta,pars,opts );
Time(1) = toc;
[MSEs(1) est_inv_rots] = check_MSE(est_inv_rots, ref_rot);
MSEs(1)
% fprintf('Reconstruction projections  from clean centered projections \n')
% [ v_PD, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(v_PD),max(real(v_PD(:)))/5);title('PDM');
% Err_v(1) = norm(volref(:) - real(v_PD(:)))/norm(volref(:))

pars.solver = 'PDMc';
tic;
[est_inv_rots ]= cryoEMPDMM(common_lines_matrix,n_theta,pars,opts );
% est_inv_rots = cryoEMPDMc(common_lines_matrix,n_theta, opts);
Time(2) = toc;
[MSEs(2), est_inv_rots] = check_MSE(est_inv_rots, ref_rot);
MSEs(2)
fprintf('Reconstruction projections  from clean centered projections \n')
% [ v_PDc, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 200, zeros(n,n,n));
% figure;isosurface(real(v_PDc),max(real(v_PDc(:)))/5);title('PDMc');
% Err_v(7) = norm(volref(:) - real(v_PDc(:)))/norm(volref(:))

pars.solver = 'PDMc2';
tic;
% est_inv_rots = cryoEMPDM(common_lines_matrix,n_theta,opts);
[est_inv_rots ]= cryoEMPDMM(common_lines_matrix,n_theta,pars,opts );
Time(3) = toc;
[MSEs(3) est_inv_rots] = check_MSE(est_inv_rots, ref_rot);
MSEs(3)
% fprintf('Reconstruction projections  from clean centered projections \n')
% [ v_PD, v_b, kernel, err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(v_PD),max(real(v_PD(:)))/5);title('PDMc2');
% Err_v(9) = norm(volref(:) - real(v_PD(:)))/norm(volref(:))

%% Test est_orientations_LS + SDPLR
% tic;
% est_inv_rots = est_orientations_LS(common_lines_matrix, n_theta);
% Time(1) = toc;
% MSEs(1) = check_MSE(est_inv_rots, ref_rot);
% fprintf('Reconstruction projections  from clean centered projections \n')
% [ v_LS, v_b, kernel ,err, iter, flag] = recon3d_firm(masked_projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(v_LS),max(real(v_LS(:)))/5);


%% With spectral norm constraint + ADMM
pars.alpha = 2/3;
tic;
est_inv_rots = est_orientations_LS(common_lines_matrix, n_theta, pars);
Time(4) = toc;
[MSEs(4), est_inv_rots]= check_MSE(est_inv_rots, ref_rot);
MSEs(4)
% [ v_LSc, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(v_LSc),max(real(v_LSc(:)))/5);title('LSc');
% Err_v(2) = norm(volref(:) - real(v_LSc(:)))/norm(volref(:))


%% Test est_orientations_LUD
% ADMM
tic;
est_inv_rots = est_orientations_LUD(common_lines_matrix,n_theta);
Time(5) = toc;
[MSEs(5), est_inv_rots]= check_MSE(est_inv_rots, ref_rot);
MSEs(5)
% [ v_LUDADMM, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(v_LUDADMM),max(real(v_LUDADMM(:)))/5);title('ADMM');
% Err_v(3) = norm(volref(:) - real(v_LUDADMM(:)))/norm(volref(:))

% ADMM with spectral norm constraint
pars.alpha = 2/3;
pars.solver = 'ADMM';
tic;
est_inv_rots = est_orientations_LUD(common_lines_matrix,n_theta, pars);
Time(6) = toc;
[MSEs(6), est_inv_rots] = check_MSE(est_inv_rots, ref_rot);
MSEs(6)
% [ v_LUDADMMc, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(v_LUDADMMc),max(real(v_LUDADMMc(:)))/5);title('ADMMc');
% Err_v(4)= norm(volref(:) - real(v_LUDADMMc(:)))/norm(volref(:))


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
pars.solver = 'IRLS';
pars.alpha = 2/3;
tic;
est_inv_rots = est_orientations_LUD(common_lines_matrix,n_theta, pars);
Time(7) = toc;
[MSEs(7), est_inv_rots] = check_MSE(est_inv_rots, ref_rot);
MSEs(7)
% [ v_IRLSc, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(v_IRLSc),max(real(v_IRLSc(:)))/5);title('IRLSc');
% Err_v(5) = norm(volref(:) - real(v_IRLSc(:)))/norm(volref(:))

%% Estimate orientations using sychronization./ eig
% tic;
% S = cryo_syncmatrix_vote(common_lines_matrix,n_theta);
% inv_rot_eig = cryo_syncrotations(S,ref_rot);
% Time(6) = toc;
% [MSEs(6), est_inv_rots] = check_MSE(inv_rot_eig, ref_rot);
% [ v_eig, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
% figure;isosurface(real(v_eig),max(real(v_eig(:)))/5);title('Eig');
% Err_v(6) = norm(volref(:) - real(v_eig(:)))/norm(volref(:))
% MSEs(6)
%fprintf('MSE of the estimated rotations: %f\n\n',mse);
%fprintf('MSE of the estimated rotations: %f\n\n',check_MSE(rotations,refq));

%l1 norm  by PH
%------------------------------------------------------------
tic;
est_inv_rot_l1 = l1_norm_rotatmatix(common_lines_matrix,n_theta, ref_rot);
Time(8) = toc;
MSEs(8) = check_MSE(est_inv_rot_l1, ref_rot);


tic;
est_inv_rots = l1_norm_rotatmatixC2(common_lines_matrix,n_theta, ref_rot);
Time(8) = toc;
est_inv_rots = permute(est_inv_rots, [2 1 3]);
[MSEs,est_inv_rots,~] = check_MSE(est_inv_rots, ref_rot);
[ v_eig, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 200, zeros(n,n,n));
figure;isosurface(real(v_eig),max(real(v_eig(:)))/8);title('l1');

tic;
est_inv_rots = l1_norm_rotatmatixC2W(common_lines_matrix,n_theta, ref_rot);
Time(8) = toc;
est_inv_rots = permute(est_inv_rots, [2 1 3]);
[MSEs,est_inv_rots,~] = check_MSE(est_inv_rots, ref_rot);
[ v_eig, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 200, zeros(n,n,n));
figure;isosurface(real(v_eig),max(real(v_eig(:)))/6);title('l1');



tic;
est_inv_rot_l1 = l1_norm_rotatmatixC(common_lines_matrix,n_theta, ref_rot);
Time(8) = toc;
MSEs(8) = check_MSE(est_inv_rot_l1, ref_rot);


tic;
est_inv_rot_l1 = l1_norm_rotatmatixC1(common_lines_matrix,n_theta, ref_rot);
Time(8) = toc;
[MSEs(8), est_inv_rots] = check_MSE(est_inv_rot_l1, ref_rot);
[ v_L1TV, v_b, kernel ,err, iter, flag] = recon3d_firm(projs,est_inv_rots,[], 1e-4, 500, zeros(n,n,n));
figure;isosurface(real(v_L1TV),max(real(v_L1TV(:)))/5);
err_v = norm(volref(:) - real(v_L1TV(:)))/norm(volref(:))
% MSEs(8)




%% Print the MSEs and cost time of the results

fprintf('SNR = %f, p = %f\n',  SNR, p);
fprintf('K = %d, n_r = %d, L = %d,masked_r = %d \n', K,n_r, n_theta,masked_r);
fprintf('-----------------------------------------\n')
fprintf(' Exp  Method   MSE        Err_v     Time\n')
fprintf('  1   PDM     %1.5f    %1.5f  %6.2f\n',   MSEs(1), Err_v(1), Time(1));
fprintf('  2   PDMc    %1.5f    %1.5f  %6.2f\n',   MSEs(2), Err_v(2), Time(2));
fprintf('  3   PDMc2   %1.5f    %1.5f  %6.2f\n',   MSEs(3), Err_v(3), Time(3));
fprintf('  4   LSc     %1.5f    %1.5f  %6.2f\n',   MSEs(4), Err_v(4), Time(4));
fprintf('  5   LUD     %1.5f    %1.5f  %6.2f\n',   MSEs(5), Err_v(5), Time(5));
fprintf('  6   LUDc    %1.5f    %1.5f  %6.2f\n',   MSEs(6), Err_v(6), Time(6));
fprintf('  7   IRLSc   %1.5f    %1.5f  %6.2f\n',   MSEs(7), Err_v(7), Time(7));
% fprintf('  7   eig     %1.5f    %1.5f  %6.2f\n',   MSEs(6), Err_v(6), Time(6));
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

% fprintf('K = %d, n_r = %d, L = %d,masked_r = %d, SNR = %f, p = %f\n', K,n_r, n_theta,masked_r, SNR, p);
% fprintf('------------------------------------------------------------------\n')
% fprintf(' Exp PDM  LS  LUD  alpha  ADMM  IRLS  eig  MSE       Err_v  Time\n')
% fprintf('  1   Y                                    %1.5f  %1.5f  %6.2f\n',   MSEs(1), Err_v(1), Time(1));
% fprintf('  2       Y         2/3    Y               %1.5f  %1.5f  %6.2f\n',   MSEs(2), Err_v(2), Time(2));
% fprintf('  3            Y           Y               %1.5f  %1.5f  %6.2f\n',   MSEs(3), Err_v(3), Time(3));
% fprintf('  4            Y    2/3    Y               %1.5f  %1.5f  %6.2f\n',   MSEs(4), Err_v(4), Time(4));
% fprintf('  5            Y    2/3           Y        %1.5f  %1.5f  %6.2f\n',   MSEs(5), Err_v(5), Time(5));
% fprintf('  6                                    Y   %1.5f  %1.5f  %6.2f\n',   MSEs(6), Err_v(6), Time(6));
% fprintf('------------------------------------------------------------------\n')





% fprintf('K = %d,  n_r = %d,     L = %d,   SNR = %f, masked_r = %d   p = %f\n', K,n_r, n_theta, SNR,masked_r, p);
% fprintf('=======================================================================================\n')
% fprintf(' Exp PDM   LS    LUD    alpha   SDPLR  ADMM  IRLS  Eig      MSE          Time\n')
% fprintf('  1   Y                                                    %1.5f       %6.2f\n',   MSEs(1), Time(1));
% fprintf('  2        Y             2/3             Y                       %1.5f     %6.2f\n',   MSEs(2), Time(2));
% fprintf('  3              Y                     Y                       %1.5f     %6.2f\n',   MSEs(3), Time(3));
% fprintf('  4              Y       2/3             Y                     %1.5f       %6.2f\n',   MSEs(4), Time(4));
% fprintf('  5              Y                            Y                %1.5f     %6.2f\n',   MSEs(5), Time(5));
% fprintf('  6              Y       2/3                    Y                %1.5f     %6.2f\n',   MSEs(6), Time(6));
% fprintf('  7                                           Y           %1.5f     %6.2f\n',   MSEs(7), Time(7));
% fprintf('  8                                                 Y    %1.5f      %6.2f\n',   MSEs(8), Time(8));
% fprintf('=======================================================================================\n')
% 
% %% 3D reconstruction
% 

% fprintf('Reconstrprojectionsuction from clean centered projections \n')
% [ v, v_b, kernel ,err, iter, flag] = recon3d_firm(masked_projs,inv_rot_matrices,[], 1e-4, 500, zeros(n,n,n));
% %fprintf('The relative error of reconstruction is %f.\n',norm(v(:)-ph(:))/norm(ph(:)));


