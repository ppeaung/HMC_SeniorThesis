%% This code generates a star graph with N nodes 
%% Inputs: 
% N: the number of nodes in the cycle
% Outputs:
% P: the probability for establishing an edge 
% G: Adjacency matrix representing the graph  
% 
% Note: requires config_model.m to generate graph
%%
function G = erdos_renyi(N,p)

adj = rand(N);
probMatrix = ones(N)*p;
adj = adj < probMatrix;
adj = adj - diag(diag(adj)); % Remove self loops
adj = triu(adj);
adj = adj + adj'; 

G = adj; 