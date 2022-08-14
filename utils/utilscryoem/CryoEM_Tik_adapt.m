function [unew,alpha] = CryoEM_Tik_adapt(T,b,eta,eta0,u,alpha)

if nargin == 4, alpha = 1e-5;  u = 0 * b;  end
if nargin == 5, alpha = 1e-5; end

   
Minvmat = @(x)(x);   %without preconditioner
unew = u; 
k    = 0; 
while 1
    k     = k + 1; 
    Treg  = @(x)(T(x) + x/alpha); 
    unew  = pcgCryoEM(Treg,Minvmat, b + u/alpha,  1e-6, 100,unew); 
    tmp   = conj(unew) .* (T(unew) - 2*b);
    Jdata = real(sum(tmp(:)) + eta0); 
    rtmp  = Jdata/eta; p = 2; 
    if rtmp<1.5 && rtmp>0.9, p = 20; end
    alpha = real(rtmp^p * alpha);
    %fprintf('%d-th iter, ratio: %1.2f , alpha=%1.2e\n',k, rtmp, alpha);
    if abs(rtmp-1)<0.0015 || k>5,break; end
    %if k>20, break; end
end

beta  = alpha; 

unew = NewtonTikDPRegParam(T,b,eta,eta0,u,beta);
