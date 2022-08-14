close all;
N = 10;
x = randn(100,N);
x_mu = mean(x,1);
x_std = std(x,0,1);
% Now make the plot:
figure();
% mean +/- std:
x_max = x_mu+x_std;
x_min = x_mu-x_std;
% XData and YData of the error bar lines:
xd = (1:N)+([-0.1; 0.1; 0; 0; -0.1; 0.1; NaN]);
yd = [x_max([1 1 1],:); x_min([1 1 1],:); NaN(1,N)];
% create the error bar lines:
line(xd(:),yd(:),'Color','k');
% create a line for marking the mean of each column of x:
line(1:N,x_mu, ...
    'Marker','s', ...
    'MarkerSize',12, ...
    'LineStyle','none', ...
    'Color','k', ...
    'MarkerFaceColor','w');
% set the grid and xlim:
grid on
xlim([0 N+1]);