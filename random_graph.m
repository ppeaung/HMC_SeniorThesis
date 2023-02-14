%% This code generates a random graph with N nodes 
%% Inputs: 
% N: the number of nodes in the random graph
% Outputs: 
% G: Adjacency matrix representing the graph  
% 
% Note: requires config_model.m to generate graph
%%
function G = random_graph(N)

graph = config_model(N, 'Poisson', 4); % note: graph is not adj matrix obj
adj = adjacency(graph);
adj = adj - diag(diag(adj)); % Remove self loops
G = adj;
