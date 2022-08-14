% Simulation experment
% 2020-9-11 By Panhuan
clear all; close all; clc;
%addpath ./LUD_new2022/fig
%addpath ./LUD_new2022/fig/PGDMSE

% 1.1-download 3D data set
load cleanrib;

k = size(volref,1);  % data set size
n = 129; %129;  %65 ;129
V = NewSizeVol(volref,n);
volref = V;


Result_LS = [];
Result_LUD = [];
K = 500; % Generate K simulated projections

for SNR =  1/64 %[1 1/4  1/8  1/16  1/32]
    result_MSE_p2q2 = [];
    result_MSE_p2q1 = [];
    
    for i = 1:1
        
        %% 1.2-Generate projections
        
        a          = qrand(K);    %
        ref_rot    = q_to_rot(a); %quat2rotm(a);
        A     = OpNufft3D(ref_rot,n); % projection operatorz
        projs = A * volref;
        % figure;viewstack(projs,5,5); % Show some noisy projections
        [noisy_projs, sigma] = ProjAddNoise(projs, SNR);
        
        % figure;viewstack(noisy_projs,5,5); % Show some noisy projections
        masked_r = 45; % mask
        masked_projs=mask_fuzzy(noisy_projs,masked_r); % Applly circular mask
        % 1.3-Compute polar Fourier transform of projections
        n_theta = 360; %360;%72
        n_r = 100;     %100;    %33
        [npf,sampling_freqs]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections
        
        % 1.4-Find common lines from Fourier projections in 3D
        max_shift=0;
        shift_step=1;
        common_lines_matrix = commonlines_gaussian(npf,max_shift,shift_step);
        C = clstack2C( common_lines_matrix,n_theta ); % common lines matrix
        % Find reference common lines and compare
        % [ref_clstack,~]=clmatrix_cheat_qq(ref_rot,n_theta);
        [ref_clstack,~]=clmatrix_cheat(ref_rot,n_theta);
        p = comparecl( common_lines_matrix, ref_clstack, n_theta, 10 ); %  detection correct rate of common lines
        
        %%   %% 2.2- PGM for LS model
        %% 2.1- PGM for LS model        
        Param.OrigRot = ref_rot;
        [est_rots,  MSEiter, REiter]= ProjGradRotLS(C, Param);
        
        
        MSE_p2q2 = check_MSE(est_rots, ref_rot);
        result_MSE_p2q2 = [result_MSE_p2q2; MSE_p2q2];
        
        
        %%  part II    %% 2.2- PGM for LUD model
        
        
        Param.OrigRot = ref_rot; Param.InitFlag = 2;
        %[est_rots,  MSEiter, REiter,ObjV]= ProjGradRotLUD(C, Param);
        [est_rots,  MSEiter, REiter]= ProjGradRotIterWLS(C, Param);
        MSE_p2q1 = check_MSE(est_rots, ref_rot);
        result_MSE_p2q1 = [result_MSE_p2q1; MSE_p2q1 ];
    end
    
    Result_LS = [Result_LS result_MSE_p2q2];
    Result_LUD = [Result_LUD result_MSE_p2q1];
    
end
figure; semilogy(ObjV);
% writematrix(Result_LS,'MES_LSK500SNR.csv');
% writematrix(Result_LUD,'MES_LUDK500SNR.csv');
% 
