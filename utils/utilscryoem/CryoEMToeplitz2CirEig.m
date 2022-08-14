function a = CryoEMToeplitz2CirEig(A)

n = size(A,1); 
n = (n+1)/2; 

a = zeros(2*n,2*n,2*n);
a(1:n,1:n,1:n)         = A(n:2*n-1,n:2*n-1,n:2*n-1);
a(n+1:2*n,1:n,1:n)     = A([n 1:n-1],n:2*n-1,n:2*n-1);
a(1:n,n+1:2*n,1:n)     = A(n:2*n-1,[n 1:n-1],n:2*n-1);
a(n+1:2*n,n+1:2*n,1:n) = A([n 1:n-1],[n 1:n-1],n:2*n-1);
a(1:n,1:n,n+1:2*n)     = A(n:2*n-1,n:2*n-1,[n 1:n-1]);
a(n+1:2*n,1:n,n+1:2*n) = A([n 1:n-1],n:2*n-1,[n 1:n-1]);
a(1:n,n+1:2*n,n+1:2*n) = A(n:2*n-1,[n 1:n-1],[n 1:n-1]);
a(n+1:2*n,n+1:2*n,n+1:2*n)=A([n 1:n-1],[n 1:n-1],[n 1:n-1]);
a = fftn(a);  %fast Fourier transform