
function inv_rotations_l1 = l1_norm_rotatmatixC(clstack,n_theta, rot_ture)

Kproj  = size(clstack,1);
%% commom matrix C is different from c
K = Kproj;
L = n_theta;
C=zeros(3,K,K);
for k1=1:K
    k2=(k1+1):K;
    l1 = clstack(k1,k2)-1;
    l2 = clstack(k2,k1)-1;
%     l1=l1(:);
%     l2=l2(:);
    
    x12 = cos(2*pi*l1/L);
    y12 = sin(2*pi*l1/L);
    x21 = cos(2*pi*l2/L);
    y21 = sin(2*pi*l2/L);
    
    C(1,k1,k2)=x12;
    C(2,k1,k2)=y12;
    C(1,k2,k1)=x21;
    C(2,k2,k1)=y21;
    
end


%% test optimal s,t
%kk =0;
% for t = 0.0001:0.0005:0.01
%     kk =kk+1; iter =0; i =0; 
    
 s = 0.05;%1.618;
 t = 0.008;%1.618;
% initial value of R
%Rold  = zeros(3,3,Kproj);
Rold  = zeros(3,3,Kproj);
Rold(1,1,:) = 1; Rold(2,2,:) = 1; Rold(3,3,:) = 1;
zold  = zeros(3, Kproj, Kproj); 
zold(3,:,:)  = 1;
tic;
for iter = 1:100
    for i = 1:Kproj
        % -- zij
        for j = i+1:Kproj
            tmpa = zold(:,i,j) + s*(Rold(:,:,i)*C(:,i,j)-Rold(:,:,j)*C(:,j,i));
            tmpb = zold(:,j,i) + s*(Rold(:,:,j)*C(:,j,i)-Rold(:,:,i)*C(:,i,j));
            tmpa(tmpa>1) = 1; tmpa(tmpa<-1) = -1;
            tmpb(tmpb>1) = 1; tmpb(tmpb<-1) = -1;
            zhat(:,i,j) = tmpa;
            zhat(:,j,i) = tmpb;            
        end
    end
    
    %% -- update R_ij
    tmp =  ones(3,3);
    for i =1:Kproj
        for j = 1:Kproj
            %tmp = tmp + c(:,i,j) * (zhat(:,i,j) - zhat(:,j,i))';
            tmp = tmp +   (zhat(:,i,j) - zhat(:,j,i))*C(:,i,j)'; % 2019-6-26
        end
      % Rnew(:,:,i) =  tmp;
         
        Rnew(:,:,i) = Rold(:,:,i) - t * tmp;
        [U,~,V] = svd(Rnew(:,1:2,i),0);
        R = U * V';
        est_inv_rots(:,:,i) = [R(:,1) R(:,2) cross(R(:,1),R(:,2))];
        
    end
    
%% -- update z_ij
    for i = 1:Kproj
        for j = i+1:Kproj
            tmpa = zold(:,i,j) + s*(Rnew(:,:,i)*C(:,i,j)-Rnew(:,:,j)*C(:,j,i));
            tmpb = zold(:,j,i) + s*(Rnew(:,:,j)*C(:,j,i)-Rnew(:,:,i)*C(:,i,j));
         
            tmpa(tmpa>1) = 1; tmpa(tmpa<-1) = -1;
            tmpb(tmpb>1) = 1; tmpb(tmpb<-1) = -1;      
            zhat(:,i,j) = tmpa;            
            zhat(:,j,i) = tmpb;            
        end
    end
%% 
% V1(:,:) = Rnew(:,1,:); V2(:,:) =Rnew(:,2,:);
% A = ATA_solver(V1,V2,K);% V1*A'=R1 and V2*A'=R2
% R1 = A*V1;
% R2 = A*V2;
% R1(:,:) = Rnew(:,1,:);
% R2(:,:) = Rnew(:,2,:);
% est_inv_rots = zeros(3,3,K);
% 
% for k=1:K
%     R = [R1(:,k) R2(:,k)];
%     [U,~,V] = svd(R,0);
%     R = U*V';
%     est_inv_rots(:,:,k) = [R(:,1) R(:,2) cross(R(:,1),R(:,2))];
% end

    ER_l1a1(iter) = check_MSE(Rnew,rot_ture);
    ER_l1c1(iter) = check_MSE(est_inv_rots,rot_ture);
    inv_rotations_l1 = permute(est_inv_rots,[2,1,3]);
    ER_l1b1(iter) = check_MSE(inv_rotations_l1,rot_ture);
     

    zold = zhat;
    Rnew = est_inv_rots;
    Rold = Rnew;   
end
figure;plot(1:iter,ER_l1a1);
figure;plot(1:iter,ER_l1c1);
figure;plot(1:iter,ER_l1b1);
fprintf('t =%1.5f', t);
time = toc



% test optimal s,t
% ER_l1aa(kk)= ER_l1a(iter); ER_l1cc(kk)= ER_l1c(iter); ER_l1bb(kk)=ER_l1b(iter);
% end
% figure;plot(1:kk,ER_l1aa);
% figure;plot(1:kk,ER_l1cc);
% figure;plot(1:kk,ER_l1bb);
 

