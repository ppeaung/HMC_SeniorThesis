%% This code runs the Adaptive HK model  
%% Inputs: 
%
% A: Adjacency matrix 
% C: Confidence Bound
% B: Beta (Homophily) Bound
% BAILOUT_TIME: Bailout time
% T: Number of Trials to Run
%%
function [X,Y, CONVERGED_TIME] = AdaptiveHK_Bailout(A, C, B, BAILOUT_TIME, T, alpha)

N = size(A, 1); % num nodes 

X = nan(N, T);
%fixed initial condition for N = 50 for testing
%X(:, 1) = [0.266284734031702;0.443434385780392;0.604934179178457;0.176540936299639;0.393020579253233;0.297511756985254;0.215711194961260;0.947739971868184;0.364007267158769;0.233742172913786;0.649220477663479;0.624389742642062;0.745741503756702;0.952381560419347;0.582549454810714;0.813708628242106;0.490194149193398;0.656263963789855;0.487039036624862;0.969342809898414;0.812469515083649;0.218062052652031;0.234913288600567;0.649817559819567;0.999300341564186;0.509650781455572;0.299597493343715;0.611423368268675;0.264517311617740;0.374170374185495;0.427728725311291;0.675666222060413;0.597440012095185;0.907773311105797;0.698254195156154;0.00863993180236089;0.0973585763834047;0.0811965673153248;0.562633787061655;0.266662760772722;0.253147160110808;0.533572261662614;0.998001894818611;0.285084524559171;0.200842568136013;0.379003574720872;0.700258078409084;0.100454234851962;0.492963643715687;0.564676448639263];
X(:, 1) = rand(N, 1);
eps_vec = zeros(N, 1) + 10^-4;
CURR_TRIAL = 1;

% Do first run manually 
CURR_OPINS = X(:, CURR_TRIAL);
 
% Add new edges based off distance 
    
DIST_MAT = 1 - abs(CURR_OPINS.'-CURR_OPINS);
PROB_MAT = rand(N); 
PROB_MAT = tril(PROB_MAT) + triu(PROB_MAT', 1); % Make symmetric 
    
PROB_MAT = PROB_MAT * alpha; 

ADD_EDGE = PROB_MAT < DIST_MAT;
   
ADD_EDGE = ADD_EDGE - diag(diag(ADD_EDGE)); % Make sure not adding self loops
  
% Keep a copy of previous
A_PREV = A;

A = A | ADD_EDGE; 
   

% Update opinions 
DIST_MAT = abs(CURR_OPINS.'-CURR_OPINS) <= C;
    
NUM = ((A+eye(N)).*DIST_MAT)*CURR_OPINS;
DEN = sum((A+eye(N)).*DIST_MAT)';
X(:, CURR_TRIAL + 1) = NUM./DEN; % Put new opinion into next 
CURR_TRIAL = CURR_TRIAL + 1;    

% Dissolve connections & update adj matrix % as a note this was first in
% the trials you did beforehand 
DISS_MAT = abs(CURR_OPINS.'-CURR_OPINS) >= B;
A = A & ~DISS_MAT;





    
% keep going if we didn't reach bailout time, if opinions aren't equal, or
% if adjacency matrix changed 

% opinion convergence is a stronger condition it seems
while ((CURR_TRIAL < BAILOUT_TIME)) && ((isequal(abs(X(:,CURR_TRIAL) - X(:,CURR_TRIAL - 1)) <= eps_vec, zeros(N,1)))) %% ~isequal(A_PREV, A))
    
    if CURR_TRIAL >= size(X, 2) % Double the size of X if we reach that point
        Z = nan(N, size(X,2)*2); 
        Z(:,1:size(X,2)) = X;
        X = Z;
    end 
        
    CURR_OPINS = X(:, CURR_TRIAL);
    
    % Add new edges based off distance 
    DIST_MAT = 1 - abs(CURR_OPINS.'-CURR_OPINS);
    PROB_MAT = rand(N); 
    PROB_MAT = tril(PROB_MAT) + triu(PROB_MAT', 1); % Make symmetric 
    
    PROB_MAT = PROB_MAT * alpha; 

    ADD_EDGE = PROB_MAT < DIST_MAT;
   
    ADD_EDGE = ADD_EDGE - diag(diag(ADD_EDGE)); % Make sure not adding self loops
    
    A_PREV = A;
    
    A = A | ADD_EDGE; 
    
   
    % Update opinions 
    DIST_MAT = abs(CURR_OPINS.'-CURR_OPINS) <= C;
    
    NUM = ((A+eye(N)).*DIST_MAT)*CURR_OPINS;
    DEN = sum((A+eye(N)).*DIST_MAT)';
    X(:, CURR_TRIAL + 1) = NUM./DEN; % Put new opinion into next 
    CURR_TRIAL = CURR_TRIAL + 1;


    % keep track of previous adj matrx to compare 
    % Dissolve connections & update adj matrix 
    
    DISS_MAT = abs(CURR_OPINS.'-CURR_OPINS) >= B;
    A = A & ~DISS_MAT;


    
    
    % checking graphs after each timestep
    %%{
    figure(20)
    plot(graph(A), 'MarkerSize', 5, 'LineWidth', 1, 'EdgeColor', 'k', 'NodeLabel', [], 'NodeCData', CURR_OPINS);
    c = colorbar;
    c.FontSize = 20; c.Label.String = 'Opinion state at convergence';
    caxis([0 1])        
    %}
    
end
CONVERGED_TIME = CURR_TRIAL;
Y = A;
end