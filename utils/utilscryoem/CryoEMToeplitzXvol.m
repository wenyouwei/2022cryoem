function y = CryoEMToeplitzXvol(a,x)

n = size(a,1)/2; 
x = reshape(x, n, n, n); 
z = zeros(2*n,2*n,2*n);
z(1:n,1:n,1:n) = x;
Ax = ( ifftn( a .* fftn(z) ) ); 

y = Ax(1:n,1:n,1:n);
%y = real(y); 
