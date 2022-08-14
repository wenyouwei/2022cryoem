
clear all; close all; clc;

%% Generate simulated projections
% Generate 200 simulated projections of size 65x65.
% For simplicity, the projections are centered.
load cleanrib;
k = size(volref,1);
n = 129;  %65
V = zeros(n,n,n);
V(1:k,1:k,1:k)= volref;
volref = V;
K   = 6; 500;
SNR = 1/16;

a              = qrand(K);    %
rotmatrices    = quat2rotm(a);
rotmatricesinv = permute(rotmatrices, [2 1 3]);

A     = OpNufft3D(rotmatrices,n); % projection operator
projs = A * volref;      % projected images

[noisy_projs, sigma] = ProjAddNoise(projs, SNR); 

ref_rot = rotmatrices; 


masked_r = 45;
masked_projs=mask_fuzzy(noisy_projs,masked_r); % Applly circular mask

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


a1 = squeeze(C(:,1,2:end));
a2 = [-C(:,2,1) zeros(2,K-2)]; 
a3 = [zeros(2,1) -C(:,3,1) zeros(2,K-3)];
a4 = [zeros(2,2) -C(:,4,1) zeros(2,K-4)];
a5 = [zeros(2,3) -C(:,5,1) zeros(2,K-5)];
a6 = [zeros(2,4) -C(:,6,1) zeros(2,K-6)];
%%
A1 = [a1; a2; a3; a4; a5; a6]; 

a1 = [zeros(2,K-2); squeeze(C(:,2,3:end))];
a2 = [-C(:,3,2) zeros(2,K-3)]; 
a3 = [zeros(2,1) -C(:,4,2) zeros(2,K-4)];
a4 = [zeros(2,2) -C(:,5,2) zeros(2,K-5)];
a5 = [zeros(2,3) -C(:,6,2) zeros(2,K-6)];
%%
A2 = [a1; a2; a3; a4; a5]; 

a1 = [zeros(4,K-3); squeeze(C(:,3,4:end))];
a2 = [-C(:,4,3) zeros(2,K-4)]; 
a3 = [zeros(2,1) -C(:,5,3) zeros(2,K-5)];
a4 = [zeros(2,2) -C(:,6,3) zeros(2,K-6)];
A3 = [a1; a2; a3; a4]; 


a1 = [zeros(6,K-4); squeeze(C(:,4,5:end))];
a2 = [-C(:,5,4) zeros(2,K-5)]; 
a3 = [zeros(2,1) -C(:,6,4) zeros(2,K-6)];
A4 = [a1; a2; a3]; 

a1 = [zeros(8,K-5); squeeze(C(:,5,6:end))];
a2 = [-C(:,6,5) zeros(2,K-6)]; 
A5 = [a1; a2]; 

A  = [A1 A2 A3 A4 A5];
[rank(A)  rank(A(:,1:12))]
