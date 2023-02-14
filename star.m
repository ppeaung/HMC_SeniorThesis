%% This code generates a star graph with N nodes 
%% Inputs: 
% N: the number of nodes in the cycle
% Outputs: 
% G: Adjacency matrix representing the graph  
% 
% Note: requires config_model.m to generate graph
%%
function G = star(N)

adj = zeros(N);
adj(:, ceil(N/2)) = ones(N,1); 
adj = adj - diag(diag(adj)); % Remove self loops
adj = adj + adj.'; 

G = adj; 


