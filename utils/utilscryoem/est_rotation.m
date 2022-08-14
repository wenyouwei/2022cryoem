function inv_rots_aligned = est_rotation(ims,n_r,n_theta,rots_true,K); 



%%% Estimate viewing angles %%%
inv_rots_true = invert_rots(rots_true);
% Calculate polar Fourier transforms of images.
pf = cryo_pft(ims, n_r, n_theta);

% Estimate common-lines matrix.
clstack_est = cryo_clmatrix(pf);

% Calculate the "true" common-lines matrix from true rotations.
clstack_true = clmatrix_cheat(rots_true, n_theta);

% Construct syncronization matrix from common lines.
S_est = cryo_syncmatrix_vote(clstack_est, n_theta);

% Estimate rotations using synchronization.
inv_rots_est = cryo_syncrotations(S_est);

% Align our estimates rotations to the true ones. This lets us compute the MSE,
% but also allows us to more easily compare the reconstructed volume.
inv_rots_aligned = align_rots(inv_rots_est, inv_rots_true);


% Calculate proportion of correctly identified common lines.
cl_prop = comparecl(clstack_est, clstack_true, n_theta, 10);

% Calculate MSE of rotations.
mse_rots = tnorm(inv_rots_aligned-inv_rots_true)^2/K;
fprintf('%-40s%20g\n', 'Proportion of correct common lines:', cl_prop);
fprintf('%-40s%20g\n', 'MSE of rotations:', mse_rots);

