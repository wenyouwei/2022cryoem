function [unew, beta] = NewtonTikDPRegParam(T,b,eta,eta0,u,beta)
%%
%%

if nargin == 4, beta = 1e-5;  u = 0 * b;  end
if nargin == 5, beta = 1e-5; end

Minvmat = @(x)(x);   %without preconditioner
unew  = u; 
k = 0;  cont = 1; 
while cont
    k    = k + 1; 
    Treg  = @(x)(beta * T(x) + x); 
    unew  = pcgCryoEM(Treg,Minvmat, beta * b + u,  1e-5, 100,unew); 
    unew  = real(unew); 
    tmp   = conj(unew) .* (T(unew) - 2*b);
    Jdata = real(sum(tmp(:)) + eta0);    
    %fprintf('eta: %1.2e, cost: %1.2e,  reg:%1.2e\n',eta,Jdata,beta); 
    a    = b - real(T(unew)); 
    z    = pcgCryoEM(@(x)(beta * T(x) + x),Minvmat, a,  1e-6, 100,unew); 
    z    = real(z); 
    dPhi = - 2 * (a.*z); dPhi = real(sum(dPhi(:))); 
    betanew = beta - min((Jdata-eta)/dPhi, beta*0.618); 
    beta = betanew;
    cont = abs(Jdata-eta)/eta>1e-2 && k<=20;
end
