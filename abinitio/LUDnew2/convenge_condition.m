function [JR,JRK] = convenge_condition (R_new,R_old,C,tk)

% JR & JRk =\sum_{i,j}||RiCij - RjCji||

K = size(R_new,3);

%% JR, JRk
JR = 0;
JRk = 0;
for  i = 1:K
    for j = 1:K
        JR  = JR  + (norm(R_new(:,:,i)*C(:,i,j) - R_new(:,:,j)*C(:,j,i),2))^2;
        JRk = JRk + (norm(R_old(:,:,i)*C(:,i,j) - R_old(:,:,j)*C(:,j,i),2))^2;
    end
end


%% Grad_R * (R_{k+1}-R_k)
Grad_R = 0;
for i = 1:K
    Grad_Ri = zeros(3,3);
    for j = 1:K
        if i ~=j
%             Grad_Ri = Grad_Ri + R_old(:,:,j) *C(:,j,i)*C(:,i,j)';
            Grad_Ri = Grad_Ri + (R_old(:,:,i) *C(:,i,j)-R_old(:,:,j) *C(:,j,i))*C(:,i,j)';

        end
    end
    Grad_R = Grad_R + sum(dot(Grad_Ri,R_new(:,:,i)-R_old(:,:,i)));
end

% || R-R_old ||_F^2
R_Rk = 0;
for i = 1:K
    R_Rk = R_Rk + (norm(R_new(:,:,i)-R_old(:,:,i), 'fro'))^2;
end
R_Rk  =  R_Rk /(2*tk);
if R_Rk ==Inf
    R_Rk =0;
end
JRK = JRk + Grad_R + R_Rk;