    function [x, error, iter] = pcgCryoEM(A,Minvmat, b, tol, max_it, x0, xtrue)
% Fastest version.
%
% Input:
%   A       convolution kernel
%   b       backprojection of the imagesguess
%   tol     error toleration
%   max_it  maximun number of iterations
%   x0      initial 
% Output:
%   x       solution
%   err     estimated error
%   iter    number of iterations
%   flag    INTEGER: 0 = solution found to tolerance
%                    1 = no convergence given max_it
%
n = size(b,1);

if isempty(tol), tol=1e-3; end
if isempty(max_it), max_it=n^3; end
if nargin<6, x0 = zeros(n,n,n); end

%%
% initialization
x    = x0;
iter = 0;
bnrm2 = norm( b(:) ); 
if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end 


%% CG method for reconstruction
r = b- A(x0); %SystemRegMatrixtimesx(A,lambda, x0);% what's
error=norm(r(:));
if ( norm(b(:)) == 0 ),    return, end
if ( error < tol ),    return, end
z = Minvmat(r);  % what's    

p = r; 
rzold = r(:)' * z(:); 
err=zeros(1,max_it);
for iter = 1:max_it                       % begin iteration
    Ap    = A(p); %SystemRegMatrixtimesx(A,lambda, p); 
    alpha = rzold /(p(:)' * Ap(:)); 
    x     = x + alpha * p;  
    x = real(x); 
    r     = r - alpha * Ap; 
    error = norm(r(:));               % check convergence   
    if error < tol, break; end
    z     = Minvmat(r); 
    rznew = r(:)' * z(:);
    beta  = rznew / rzold; 
    p     = r + beta * p; 
    rzold = rznew; 
    err(iter)=error/bnrm2;
    
    if ( error/bnrm2 <= tol ),
        %os = sprintf(['CPCG converged at iteration %d to a solution with '...
        %    'relative residual %0.2g'], iter,error/bnrm2);
        %disp(os);
        break, 
    end
    if isvar('xtrue'), reerr(iter) = norm(x(:)-xtrue(:))/norm(xtrue(:)); end
    
    %if mod(iter,50)==1
    %    fprintf('iter %d  vol error:%f\n',iter,err(iter));
    %end
    
end
plot(err); 
%fprintf('iter %d  vol error:%f\n',iter,err(iter));

% if ( error/bnrm2 > tol ),
%     os = sprintf(['CPCG stopped after %d iteration'...
%         'and has relative residual %0.2g'],max_it,error/bnrm2);
%     flag = 1; disp(os);                   % no convergence
% end