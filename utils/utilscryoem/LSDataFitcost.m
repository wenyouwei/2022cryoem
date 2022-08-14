function tau = LSDataFitcost(ToepCirEig,u,back_projs,obsimg,upbound)
% this function is to calculate the value of the data fitting term.

tmpb = CryoEMToeplitzXvol(ToepCirEig,u) - 2 * back_projs; 
tmpb = real(u .* tmpb); 
Phiu = sum(tmpb(:)) + norm(obsimg(:))^2;
tau  = Phiu /upbound; 
 
