clear
K = 6; 
theta = 1/K:1/K:1;     

noise = randn(size(theta));
theta2 = theta * pi.*(1+0.001*noise);

a1 = cos(theta2); b1 = sin(theta2); aa = [a1; b1];

c1 =[aa(:,1);-aa(:,1);zeros(8,1)];
c2 =[aa(:,2);0;0;-aa(:,2);zeros(6,1)];
c3 =[aa(:,3);zeros(4,1);-aa(:,3);zeros(4,1)];
c4 =[aa(:,4);zeros(6,1);-aa(:,4);zeros(2,1)];
c5 =[aa(:,5);zeros(8,1);-aa(:,5)];

A0  =[c1 c2 c3 c4 c5]; 

% noise = randn(size(theta));
% theta2 = theta * pi.*(1+0.01*noise);
% 
% a1 = cos(theta2); b1 = sin(theta2); aa = [a1; b1];
% 
% c1 =[aa(:,1);-aa(:,1);zeros(8,1)];
% c2 =[aa(:,2);0;0;-aa(:,2);zeros(6,1)];
% c3 =[aa(:,3);zeros(4,1);-aa(:,3);zeros(4,1)];
% c4 =[aa(:,4);zeros(6,1);-aa(:,4);zeros(2,1)];
% c5 =[aa(:,5);zeros(8,1);-aa(:,5)];

A  =[c1 c2 c3 c4 c5]; 


B  = zeros(12,4); B(3:end,:) = A(1:end-2,1:4); 


A  =[c1 c2 c3 c4 c5]; 


C  = zeros(12,3); C(5:end,:) = A(1:end-4,1:3); 

  
D  = zeros(12,2); D(7:end,:) = A(1:end-6,1:2); 
 
A  =[c1 c2 c3 c4 c5]; 

E  = zeros(12,1); E(9:end,:) = A(1:end-8,1); 

%AA = [A0(:,1:4) B(:,1:3) C D E];

AA = [A0(:,1:4) B(:,1:3) C(:,1:2) D(:,1)];
%AA = [A0(:,1:3) B(:,1:2) C(:,1)];
rank(AA)
