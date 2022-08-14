function  theta = CC2theta(phi,C,mu)

for i=1:K
    for j=i+1:K
      t=0;
      for k=1:2
      theta[j*2*K+i*2+k] = C[j*2*K+i*2+k]-*mu*(phi[2*j*2*K+2*i+k]*C[i*2*K+j*2]+phi[(2*j+1)*2*K+2*i+k]*C[i*2*K+j*2+1]);
      t = t + theta[j*2*K+i*2+k]* theta[j*2*K+i*2+k];
      end
      t=sqrt(t);
      for k=1:2
      theta[j*2*K+i*2+k] = theta[j*2*K+i*2+k]/t;
      S[2*j*2*K+2*i+k] = theta[j*2*K+i*2+k]*C[i*2*K+j*2];
      S[(2*j+1)*2*K+2*i+k] = theta[j*2*K+i*2+k]*C[i*2*K+j*2+1];
      end
    end
 end
