%% This code generates a cycle graph with N nodes 
%% Inputs: 
% N: the number of nodes in the cycle
% Outputs: 
% G: Adjacency matrix representing the graph  
% 
% Note: requires config_model.m to generate graph
%%
function G = cycle(N)

v = ones(1, N-1); 
adj = diag(v, 1);
adj(N, 1) = 1;
adj = adj + adj.';

G = adj;
