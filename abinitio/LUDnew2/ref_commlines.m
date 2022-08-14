function common_lines_matrix= ref_commlines(rot_matrices, n_theta,p)

if nargin < 3
    p = 0;
end
L = n_theta;
K = size(rot_matrices ,3);

inv_rot_matrices = permute(rot_matrices,[2 1 3]);
n_x(:,:) = inv_rot_matrices(:,1,:);
n_y(:,:) = inv_rot_matrices(:,2,:);
n_z(:,:) = inv_rot_matrices(:,3,:);

fprintf('Computing common lines... ');
% Find common lines -- the intersection points of the circles
common_lines_matrix = zeros(K);
eqs_matrix = zeros(3,4);


   for i=1:K
     for j=(i+1):K
        % 3 equations in 4 unknowns varibles 
        % Cij = [Ci1,Ci2,0];
        % Cji = [Cj1, Cj2,0];
        % Ri*Cij = Rj*Cji ==> Ri*Cij - Rj*Cji =0 ==>:
        eqs_matrix = [n_x(1,i), n_y(1,i), -n_x(1,j), -n_y(1,j) ; ...
            n_x(2,i), n_y(2,i), -n_x(2,j), -n_y(2,j) ; ...
            n_x(3,i), n_y(3,i), -n_x(3,j), -n_y(3,j)];
        Z = 2*null(eqs_matrix);   
        
        phi_1 = atan2(Z(2),Z(1));
        if (phi_1 < 0)
            phi_1 = phi_1 + 2*pi;
        end
        
        phi_2 = atan2(Z(4),Z(3));
        if (phi_2 < 0)
            phi_2 = phi_2 + 2*pi;
        end
        
        dtheta = L/(2*pi);
        l1 = mod(round(phi_1*dtheta)+L-1,L)+1;
        l2 = mod(round(phi_2*dtheta)+L-1,L)+1;
        common_lines_matrix(i,j) = l1;
        common_lines_matrix(j,i) = l2;
     end
  end

%% Perturb the common lines matrix
if p~= 0
    is_perturbed = zeros(K); % to verify the success of the algorithm we store which common lines were perturbed
    for i=1:K
        for j=(i+1):K
            r = rand;
            if (r > p) % with probability 1-p the common line needs to be perturbed
            is_perturbed(i,j) = 1;
            common_lines_matrix(i,j) = floor(rand*L)+1;
            common_lines_matrix(j,i) = floor(rand*L)+1;
            end
        end
    end
end
fprintf('Finished!\n');