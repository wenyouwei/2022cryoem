function [ W ]=script_find_graph_weights_v3(X, c, n, k)
%Function for using linear programming to construct graph that has same
%weight in each node. This function requires cvx tool box.
%Need modification
%   Input: 
%       X: distance between nearest neighbors. size: number of edges x 1.
%       c: nearest neighbor list. Matrix size: number of edges x 2.
%       n: size of the graph.
%       k: weight for each node.
%   Output: 
%       W: weight on each edge ij.

% cd /scratch/summer2011/post_classaveraging/
% setup;
% cd /u/zhizhenz/matlab_cryo

N=length(X);

cvx_begin
    variable W(N)
    minimize( W'*X );
    subject to
        W<=1;
        W>=0;
        for i=1:n
            id1=find(c(:, 1)==i);
            id2=find(c(id1, 2)>i);
            id3=find(c(:, 2)==i);
            id4=find(c(id3, 1)<i);
            id=[id1(id2); id3(id4)];
            sum(W(id))==k;
        end;
cvx_end
