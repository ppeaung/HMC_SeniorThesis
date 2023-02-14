%% This code generates a k-partite graph with N nodes 
%% Inputs: 
% N: the number of nodes in the cycle
% Outputs: 
% G: Adjacency matrix representing the graph  
% 
% Note: requires config_model.m to generate graph
%%
function G = k_partite(N, k, v)

% alg outline 
% make adj nan matrix size k x k 
% make zero matrices of size k x k 
% loop through v pairwise k-1 times, make ones matrices of size v_i x v_j

% adj will have entry 
% adj(k,k) = zeros(k) 
% adj(k, k+1) = ones(v_i, v_j)

adj = nan(k);
for i = 1:k
    adj(i,i) = zeros(i);
    adj(i, i+1) = ones(v(i), v(i+1));
end

%remove self loops 

adj = adj + adj.';
G = adj;


% make all the 1's and 0's 



% loop through k and make adj matrix 


% revist this again 


% Expand this to generalize to k partite problem 
% Complete Bipartite 
% Assume N even. Have the two partite sets equal 

                %{
                A = ones(N/2);
                B = zeros(N/2);
                AdjacencyMatrix = [B, A; A, B.'];
                %}
                
                % Ask about complete k-partite / explore more 