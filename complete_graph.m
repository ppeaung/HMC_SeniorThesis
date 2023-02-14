%% This code generates a complete graph with N nodes 
%% Inputs: 
% N: the number of nodes in the complete graph
% Outputs: 
% G: Adjacency matrix representing the complete graph  
%%
function G = complete_graph(N)

adj = ones(N);
adj = adj - diag(diag(adj)); % Remove self loops

G = adj;



