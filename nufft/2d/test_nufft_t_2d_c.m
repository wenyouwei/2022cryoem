function test_nufft_t_2d_c
%
% Test the function nufft_t_2d_optimized.
% See test_nufft_t_2d_a for more information.
%
% Yoel Shkolnisky, December 2009.

m=128;
n_tests=50;

for k=1:n_tests
    n=min(2^k,512);
    fprintf('Test for size n=%4d...',n)
    x=rand(m,2); % sampling points
    beta=rand(n);
    
    tic
    v1=nufft_t_2d(beta,x);
    t1=toc;
        
    tic;
    v2=nufft_t_2d_optimized(beta,x,'double');
    t2=toc;
    
    d=norm(v1-v2)/norm(v1);

    if t1<1.e-12
        t1=0;
        t2=1;
    end

    if d<1.0e-5
        fprintf('OK    d=%7.4e    r=%7.4f\n',d,t1/t2);
    else
        fprintf('ERROR\n');
    end
end