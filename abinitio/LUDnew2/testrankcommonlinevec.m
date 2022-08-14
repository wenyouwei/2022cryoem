
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
K   = 400; %500;
SNR = 1/16;
for iii = 1:1
a              = qrand(K);    %
%rotmatrices    = quat2rotm(a);
rotmatrices    = q_to_rot(a); %quat2rotm(a); 
rotmatricesinv = permute(rotmatrices, [2 1 3]);
A     = OpNufft3D(rotmatrices,n); % projection operator
projs = A * volref;      % projected images
[noisy_projs, sigma] = ProjAddNoise(projs, SNR); 
ref_rot = rotmatrices; 
masked_r = 45;
masked_projs=mask_fuzzy(noisy_projs,masked_r); % Applly circular mask

%% Compute polar Fourier transform, using radial resolution n_r and angular
% resolution n_theta. n_theta is the same as above.
n_theta = 360; %360;%72
n_r = 100;     %100;    %33
[npf,sampling_freqs]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections   

%% Find common lines from projections
max_shift=0;
shift_step=1;
common_lines_matrix = commonlines_gaussian(npf,max_shift,shift_step);
C = clstack2C( common_lines_matrix,n_theta );


%% compute  column rank (ZZ) and prove column rank (ZZ) is full
% common_lines is polluted by noisy
% 2022-7-17 

ZZ = zeros( 2*K, K*(K-1)/2 );
for i = 1
   for  j = i+1:K
       ZZ(2*i-1:2*i,j-1) =  C(:,i,j);
       ZZ(2*j-1:2*j,j-1)= - C(:,j,i);
   end 
end

%  when i>=2, There is a law in the column elements of the ZZ \

for i = 2:K-1
   for  j = i+1:K 
       ZZ(2*i-1:2*i,j-1+(i-1)*K-i*(i+1)/2+1) =  C(:,i,j);
       ZZ(2*j-1:2*j,j-1+(i-1)*K-i*(i+1)/2+1)= - C(:,j,i);
   end 
end
rank(ZZ)
% A = ZZ(:,1:2*K); %whos A
% if rank(A)<12
%     disp('wrong');rank(A)
% end
end
return;
% 
% for i = 2
%    for  j = i+1:K 
%        ZZ(2*i-1:2*i,j-1+K-i) =  C(:,i,j)
%        ZZ(2*j-1:2*j,j-1+K-i)= - C(:,j,i)
%    end 
% end
% 
% for i = 3
%    for  j = i+1:K 
%        ZZ(2*i-1:2*i,j-1+K-2+K-i) =  C(:,i,j)
%        ZZ(2*j-1:2*j,j-1+K-2+K-i)= - C(:,j,i)
%    end 
% end
% 
% 
% for i = 4
%    for  j = i+1:K 
%        ZZ(2*i-1:2*i,j-1+K-2+K-3+K-i) =  C(:,i,j)
%        ZZ(2*j-1:2*j,j-1+K-2+K-3+K-i)= - C(:,j,i)
%    end 
% end
% 
% for i = 5
%    for  j = i+1:K 
%        ZZ(2*i-1:2*i,j-1+K-2+K-3+K-4+K-i) =  C(:,i,j)
%        ZZ(2*j-1:2*j,j-1+K-2+K-3+K-4+K-i)= - C(:,j,i)
%    end 
% end
% 
% 
% for i = 6
%    for  j = i+1:K 
%        ZZ(2*i-1:2*i,j-1+K-2+K-3+K-4+K-5+K-i) =  C(:,i,j)
%        ZZ(2*j-1:2*j,j-1+K-2+K-3+K-4+K-5+K-i)= - C(:,j,i)
%    end 
% end
% 
% rank(ZZ)


%% method 1 by Wen 2022-07-17
%% i = 1; j = 2:K
a1 = squeeze(C(:,1,2:end));
a2 = [-C(:,2,1) zeros(2,K-2)]; 
a3 = [zeros(2,1) -C(:,3,1) zeros(2,K-3)];
a4 = [zeros(2,2) -C(:,4,1) zeros(2,K-4)];
a5 = [zeros(2,3) -C(:,5,1) zeros(2,K-5)];
a6 = [zeros(2,4) -C(:,6,1) zeros(2,K-6)];
a7 = [zeros(2,5) -C(:,7,1) zeros(2,K-7)];
a8 = [zeros(2,6) -C(:,8,1) zeros(2,K-8)];
A1 = [a1; a2; a3; a4; a5; a6;a7;a8]; 

%% i = 2; j = 3:K
a1 = [zeros(2,K-2); squeeze(C(:,2,3:end))];
a2 = [-C(:,3,2) zeros(2,K-3)]; 
a3 = [zeros(2,1) -C(:,4,2) zeros(2,K-4)];
a4 = [zeros(2,2) -C(:,5,2) zeros(2,K-5)];
a5 = [zeros(2,3) -C(:,6,2) zeros(2,K-6)];
a6 = [zeros(2,4) -C(:,7,2) zeros(2,K-7)];
a7 = [zeros(2,5) -C(:,8,2) zeros(2,K-8)];
%%
A2 = [a1; a2; a3; a4; a5;a6;a7]; 

%% i = 3; j = 4:K
a1 = [zeros(4,K-3); squeeze(C(:,3,4:end))];
a2 = [-C(:,4,3) zeros(2,K-4)]; 
a3 = [zeros(2,1) -C(:,5,3) zeros(2,K-5)];
a4 = [zeros(2,2) -C(:,6,3) zeros(2,K-6)];
a5 = [zeros(2,3) -C(:,7,3) zeros(2,K-7)];
a6 = [zeros(2,4) -C(:,8,3) zeros(2,K-8)];
A3 = [a1; a2; a3; a4;a5;a6]; 

%% i = 4; j = 5:K
a1 = [zeros(6,K-4); squeeze(C(:,4,5:end))];
a2 = [-C(:,5,4) zeros(2,K-5)]; 
a3 = [zeros(2,1) -C(:,6,4) zeros(2,K-6)];
a4 = [zeros(2,2) -C(:,7,4) zeros(2,K-7)];
a5 = [zeros(2,3) -C(:,8,4) zeros(2,K-8)];
A4 = [a1; a2; a3;a4;a5]; 

%% i = 5; j = 6:K
a1 = [zeros(8,K-5); squeeze(C(:,5,6:end))];
a2 = [-C(:,6,5) zeros(2,K-6)]; 
a3 = [zeros(2,1) -C(:,7,5) zeros(2,K-7)]; 
a4 = [zeros(2,2) -C(:,8,5) zeros(2,K-8)]; 
A5 = [a1; a2;a3;a4]; 

%% i = 6; j = 7:K
a1 = [zeros(10,K-6); squeeze(C(:,6,7:end))];
a2 = [-C(:,7,6) zeros(2,K-7)]; 
a3 = [zeros(2,1) -C(:,8,6) zeros(2,K-8)]; 
A6 = [a1; a2;a3]; 

%% i = 7; j = 8:K 
a1 = [zeros(12,K-7); squeeze(C(:,7,8:end))];
a2 = [-C(:,8,7) zeros(2,K-8)]; 

A7 = [a1; a2]; 

A  = [A1 A2 A3 A4 A5 A6 A7]; 
[rank(A)  rank(A(:,1:28))] 

%%





