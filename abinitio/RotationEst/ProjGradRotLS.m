function [est_rots,MSEiter, REiter,ObjValue]= ProjGradRotLS(C, Param)
% 2019-10-14
%% Object:
%  \sum_{ij} \min 1/2*|| RiCij - RjCji||_2^2
%  s.t  Ri'Ri = I
%  Ri^{k+1} = Ri - t*(RiCij - RjCji)*Cij^{T}
%% Proximal Gradient with Backtracking Line-Search

MaxIter = 1000;   SolRE = 1e-6;
K  = size(C,3);

if nargin == 2
    if isfield(Param,'OrigRot'),    truerot   = Param.OrigRot;    end
    if isfield(Param,'InitRot'),    R_old     = Param.InitRot;    end
    if isfield(Param,'MaxIter'),    MaxIter   = Param.MaxIter;    end
    if isfield(Param,'SolRE'),      SolRE     = Param.SolRE;      end
end

if ~exist('R_old','var')    % initial the rotation matrix value
    R_old   = rand(3,3,K);
    for i = 1:K
        [U,~,V]   = svd(R_old(:,:,i));
        R_old(:,:,i) = U * V';
    end
end


TOL   = 1.0e-14;
MSEiter = 0;

tk    = 1; %[0.99,0.7,0.618,0.5,0.2,0.08]
alpha = 0.618; %0.618;%0.99; % SNR = 1/16 alpha =0.618
R_new = R_old;

for iter = 1:MaxIter
    flag = 1;
    kk   = 0;
    [GradRi,JRK] = FunCostAndGradp2q2(R_old(:,1:2,:),C);
    Grad_Ri      = zeros(3,3,K);
    Grad_Ri(:,1:2,:) = GradRi;
    
    tk = min(0.618,tk*2.5);
    while flag
        kk = kk+1;
        GradObj = 0;
        for i =1:K
            temp    = R_old(:,:,i) - tk * Grad_Ri(:,:,i);
            [U,~,V] = svd(temp,0);
            R_new(:,:,i) = U*V';
            
            tmp     = Grad_Ri(:,:,i).* (R_new(:,:,i)-R_old(:,:,i));
            GradObj = GradObj + sum(tmp(:)); % Grad_R*(Rnew-Rold)
        end
        
        [~,JR] = FunCostAndGradp2q2(R_new(:,1:2,:),C);
        JRKnew = JRK + GradObj + (norm(R_new(:)-R_old(:))^2)/(2*tk);
        
        if JR > JRKnew && kk<20
            tk = alpha*tk;
        else
            flag = 0;
        end
    end
    %fprintf('step=%2e\n',tk);
    
    REiter(iter) = norm(R_new(:) - R_old(:))/norm(R_new(:));
    
    R_old = R_new;
        
    ObjValue(iter) = JR; 

    %% Make sure that we got true rotations.
    if exist('truerot','var')        
        est_rots = zeros(3,3,K);
        for k=1:K
            est_rots(:,:,k) = [R_new(:,1:2,k),cross(R_new(:,1,k),R_new(:,2,k))];
            R = est_rots(:,:,k);
            erro = norm(R*R.'-eye(3));
            if erro > TOL || abs(det(R)-1)> TOL
                [U,~,V] = svd(R);
                est_rots(:,:,k) = U*V.';
            end
        end
        [MSEiter(iter),~,~] = check_MSE(est_rots, truerot);
    end
    if REiter(iter) <SolRE, fprintf('iter = %d\n',iter); break; end
    
end
est_rots = zeros(3,3,K);
for k=1:K
    est_rots(:,:,k) = [R_new(:,1:2,k),cross(R_new(:,1,k),R_new(:,2,k))];
    R = est_rots(:,:,k);
    erro = norm(R*R.'-eye(3));
    if erro > TOL || abs(det(R)-1)> TOL
        [U,~,V] = svd(R);
        est_rots(:,:,k) = U*V.';
    end
end
% figure(9); subplot(2,1,1); semilogy(REiter); subplot(2,1,2); semilogy(MSEiter);