function [est_rots,MSEiter, REiter,B]= ProjGradRotIterWLS(C, Param)
% 2019-10-14
%% Object:
%  \sum_{ij} \min 1/2*wij*|| RiCij - RjCji||_2^2
%  s.t  Ri'Ri = I
%  Ri^{k+1} = Ri - t*wij*(RiCij - RjCji)*Cij^{T}
%% Proximal Gradient with Backtracking Line-Search

MaxIter = 100;   SolRE  = 1e-6; %Rflag = 2;
K  = size(C,3);
MSEiter = [];     REiter = [];
W     = ones(K,K);

if nargin == 2
    if isfield(Param,'OrigRot'),    truerot   = Param.OrigRot;    end
    if isfield(Param,'InitRot'),    R_old     = Param.InitRot;    end
    %if isfield(Param,'InitFlag'),   Rflag     = Param.InitFlag;   end
    if isfield(Param,'MaxIter'),    MaxIter   = Param.MaxIter;    end
    if isfield(Param,'SolRE'),      SolRE     = Param.SolRE;      end
end

if ~exist('R_old','var')    % initial the rotation matrix value
    %     if Rflag == 1
    R_old   = rand(3,3,K);
    for i = 1:K
        [U,~,V]   = svd(R_old(:,:,i));
        R_old(:,:,i) = U * V';
    end
    MSEiter = []; REiter = [];
    %     else
    %         Param.SolRE = 1e-3;
    %         [R_old,  MSEiter, REiter]= ProjGradRotLS(C, Param);
    %         [GradRi,JRK,W] = FunCostAndGradp2q2w(R_old(:,1:2,:),C,W);
    %     end
    
end

TOL   = 1.0e-14;


tk    = 1; %[0.99,0.7,0.618,0.5,0.2,0.08]
alpha = 0.618; %0.618;%0.99; % SNR = 1/16 alpha =0.618
R_newa = R_old;

for iter = 1:MaxIter
    %     flag = 1;
    %     kk   = 0;
    %     Param1  = Param;
    %     Param1.InitRot = R_old;  Param1.MaxIter = 10; Param1.SolRE = 1e-4;
    %
    %
    %
    %     [R_new,  MSEiter2, REiter2, W]= ProjGradRotWLS(C,W, Param1);
    %     MSEiter = [MSEiter MSEiter2];
    %     REiter  = [REiter REiter2];
    
    
    R_old2 = R_old;
    clear MSEiter2 REiter2 ObjValue
    
    for ii = 1:100
        flag = 1;
        kk   = 0;
        [GradRi,JRK,B] = FunCostAndGradp2q2w(R_old(:,1:2,:),C,W);
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
            
            [~,JR,B] = FunCostAndGradp2q2w(R_new(:,1:2,:),C,W);
            JRKnew = JRK + GradObj + (norm(R_new(:)-R_old(:))^2)/(2*tk);
            
            if JR > JRKnew && kk<20
                tk = alpha*tk;
            else
                flag = 0;
            end
        end
        %fprintf('step=%2e\n',tk);
        %W = B;
        REiter2(ii) = norm(R_new(:) - R_old(:))/norm(R_new(:));
        
        R_old = R_new;
        
        ObjValue(ii) = JR;
        
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
            [MSEiter2(ii),~,~] = check_MSE(est_rots, truerot);
        end
        if REiter2(ii) <1e-4,  break; end
        %if iter > 2, break; end
    end
    W       = B;
    MSEiter = [MSEiter MSEiter2];
    REiter  = [REiter REiter2];
    tmp   = R_new(:,1:2,:) - R_old2(:,1:2,:);
    
    REiterb  = norm(tmp(:))/norm(R_new(:));
    
    R_old = R_new;
    
    
    
    if REiterb  <SolRE, fprintf('iter = %d\n',iter); break; end
    
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
%
% MSEiter = [MSEiter1 MSEiter];
% REiter  = [REiter1 REiter];