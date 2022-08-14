%% Experiments with simulated rotation
% 2022-7-22
clc; clear all; close all
inneriter = 10;

fp = fopen('./abinitio/RotationEst/effectinitial.txt','a+')
%% parameter set
for K = 1000  %[ 500 1000 2000, 3000]
    fprintf(fp,'\n\nK=%4d\n',K);
    n_theta = 360; %[ 72, 90, 120, 144, 180, 360] %360;%72 120 180
    ref_rot = rand_rots(K);
    inv_rot_matrices = permute(ref_rot,[2 1 3]);
 
    fprintf('Computing common lines... ');
    is_perturbed = zeros(K); % to verify the success of the algorithm we store which common lines were perturbed
    
    for pp = [0.9 0.7 0.5   0.3] % the proportion of correctly detected common lines
        fprintf(fp,'\np=%4f\n',pp);
        common_lines_matrix= ref_commlines(ref_rot, n_theta,pp);
        
        %% Test different orientation determination algorithms
        tic;
        C = clstack2C( common_lines_matrix,n_theta );
        tt = toc;
        
        %% random initial value for LS
        for j = 1:inneriter
            Param.OrigRot = ref_rot;
            [est_rots,  MSEiter, REiter]= ProjGradRotLS(C, Param);
            figure(100); semilogy(MSEiter);hold on; pause(0.1)
            ppMSELS(j) = MSEiter(end);
        end
        fprintf(fp,'LS: %1.4e %1.4e %1.4e %1.4e',mean(ppMSELS),min(ppMSELS),max(ppMSELS),std(ppMSELS));
           
        ylabel('MSE');xlabel('Iteration');
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'PaperSize',[11 8]);
        fig = gcf;
        fig.PaperUnits = 'centimeters';
        fig.PaperPosition = [0 0 11 8];
        fig.Units = 'centimeters';
        fig.PaperSize=[11 8];
        fig.Units = 'centimeters';
        print(fig,'-dpdf',strcat('./abinitio/RotationEst/K',num2str(K),'p',num2str(10*pp),'MESLS'));
        
        
        %% LS initial value for LUD
        for j = 1:inneriter
            Param.OrigRot = ref_rot; Param.InitFlag = 2;
            
            %W = ones(K,K);
            [est_rots,  MSEiter, REiter]= ProjGradRotIterWLS(C, Param);
            figure(102);semilogy(MSEiter);hold on; pause(0.1)
            ppMSELUD(j) = MSEiter(end);
        end
        
        fprintf(fp,'\nLUD: %1.4e %1.4e %1.4e %1.4e',mean(ppMSELUD),min(ppMSELUD),max(ppMSELUD),std(ppMSELUD));
        
        ylabel('MSE');xlabel('Iteration');
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'PaperSize',[11 8]);
        fig = gcf;
        fig.PaperUnits = 'centimeters';
        fig.PaperPosition = [0 0 11 8];
        fig.Units = 'centimeters';
        fig.PaperSize=[11 8];
        fig.Units = 'centimeters';
        print(fig,'-dpdf',strcat('./abinitio/RotationEst/K',num2str(K),'p',num2str(10*pp),'MESLUD'));
        
        close all;
        
    end
end


