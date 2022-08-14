function OutPut = CryoEMPrimalDual3DTV(H, b, alpha, voltrue,LSDF)


MaxIter = 100;              SolRE  = 1e-4; tol = 1e-4; 
n       = size(b,1); 
re     = 1; k = 0; 
s      = min(1/alpha,.618); %s = alpha
t      = 5e-2; %1/32/s;  % t = beta


Heig = CryoEMToeplitz2CirEig(H);
%A       = fftnCryoEMkernel(H); 

start_time = cputime;
vol    = zeros(n,n,n); 
vxi    = GradVol3D(vol);
beta   = t * alpha; 
T = @(x)CryoEMToeplitzXvol(Heig,x); 
while (k<MaxIter)&&(re>SolRE)
    k = k+1;  
    
    %% update the primal variable 
    v      = vol - t * Divz3D(vxi);
    rhs    = beta * b + v; 
    %volnew = pcg(@(x)AtAReg(x),rhs(:),tol, 100,[],[],vol(:)); 
    T2     = @(x)(beta * CryoEMToeplitzXvol(Heig,x)+ x); %(A'A + alpha*I)*x
    Minvmat = @(x)(x); 
    volnew = pcgCryoEM(T2,Minvmat, rhs, tol, 100, vol); % A'A x = Atb
    %Costdatafitting(Heig, volnew, rhs);
    %volnew(volnew<0)=0;
    
    %if mod(k,5)==0
    %    tau = LSDF(volnew); 
    %    fprintf('k = %5d  vol error:%f,    tau =%f\n',k,norm(volnew(:)-voltrue(:))/norm(voltrue(:)),tau);
    %end
    volhat = volnew + volnew - vol;
        
    dvol   = GradVol3D(volhat); 
    vxi    = VolGradProjLinfty(vxi - s * dvol);  
    
    %% relative error
    re   = norm(volnew(:)-vol(:))/norm(volnew(:));
    vol  = volnew;
    %t  = t * 1.05; 
    
    %if isvar('xtrue'),      ISNR(k) = fun_ISNR(xtrue,g,f);     end
    %iter_time(k) = cputime - start_time;
    %regpar(k)    = mu/t;         
end

OutPut.Sol      = vol;
OutPut.vxi      = vxi;


function q = VolGradProjLinfty(p)
% email: wenyouwei@gmail.com

tmp = sqrt(p(1,:,:,:).^2+p(2,:,:,:).^2+p(3,:,:,:).^2);
tmp(tmp<1) = 1; 

q(1,:,:,:) = p(1,:,:,:)./tmp; 
q(2,:,:,:) = p(2,:,:,:)./tmp; 
q(3,:,:,:) = p(3,:,:,:)./tmp; 


function z = GradVol3D(u)

n  = size(u, 1); 
z  = zeros(3,n, n,n); 
z(1,1:end-1,:,:) = u(2:end,:,:) - u(1:end-1,:,:);  %kernel [0;-1;1]  
z(2,:,1:end-1,:) = u(:,2:end,:) - u(:,1:end-1,:);  %kernel [0 -1 1]
z(3,:,:,1:end-1) = u(:,:,2:end) - u(:,:,1:end-1);    

% z.dx(n,:,:) = zeros(n,n); 
% z.dy(:,n,:) = zeros(n,n);
% z.dz(:,:,n) = zeros(n,n);

function w = Divz3D(z)

zx(:,:,:) = z(1,:,:,:); % x-axis
zy(:,:,:) = z(2,:,:,:); % y-axis
zz(:,:,:) = z(3,:,:,:); % z-axis
n  = size(zx, 1); 

tmp1 = zx(1,:,:); tmp2 = zy(:,1,:); tmp3 = zz(:,:,1); 

tmp1(2:n-1,:,:) = zx(2:n-1,:, :) - zx(1:n-2, :, :);  %kernel [0;-1;1]  
tmp2(:,2:n-1,:) = zy(:,2:n-1, :) - zy(:,1:n-2,  :);  %kernel [0;-1;1]  
tmp3(:,:,2:n-1) = zz(:, :,2:n-1) - zz(:,  :,1:n-2);  

tmp1(n,:,:) = -zx(n-1,:,:); 
tmp2(:,n,:) = -zy(:,n-1,:); 
tmp3(:,:,n) = -zz(:,:,n-1); 

w = tmp1 + tmp2 + tmp3; 
