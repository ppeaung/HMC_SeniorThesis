%% This code runs the Sigmoidal Opinion Dynamic Model via the 4th order RK 
%% Inputs: 
%
% A: Adjacency matrix 
% X: array of opinions 
% delta: Confidence bound
% gamma: controls "sharpness" of dependence of nodes i and j, as gamma 
% approaches infinity it approaches a stepwise function (classical HK
% model) 
% Z: the set of Zealots or nodes that are not persaudable

function Y = sigmoidal(A, X, delta, gamma, Z) 

% fixed number of iterations
numIterations = 5000; 
Y = nan(size(X, 1), numIterations); % size size(X,1) x numIterations 

% first index is the size of 

Xnew = X; % initialize X 

h = 0.01; % What do I choose as h here? was it possible that h was too small? 

% simulation via RK4 
for i = 1:numIterations
    Xold = Xnew; % reset old value
    a = h*sigmoidalDE(A, Xold, gamma, delta); 
    b = h*sigmoidalDE(A, Xold + (1/2)*a, gamma, delta); 
    c = h*sigmoidalDE(A, Xold + (1/2)*b, gamma, delta);
    d = h*sigmoidalDE(A, Xold + c, gamma, delta);
    
    % modify so that there are no changes in zealots 
    if size(Z,1) > 0
        for j = 1:size(Z,1)
        a(Z(j)) = 0;
        b(Z(j)) = 0; 
        c(Z(j)) = 0;
        d(Z(j)) = 0; 
        end
    end
    Xnew = Xold + (1/6)*(a + 2*b + 2*c + d);
    Y(:,i) = Xnew; 
end

% Helper functions 
% Input: adj matrix A, vector of opinions X, and parameters gamma and delta
% Returns: function representing DE for sigmoidal DE that's coupled 
    function f = sigmoidalDE(A, X, gamma, delta)
        % calculate weight matrix 
        W = 1./(1 + exp(gamma*(X-X').^2 - gamma*delta)); 
        W = W - diag(diag(W)); % remove diagonal weights (should we do this)?  
        if sum(A.*W, 1) == 0 % don't divide by 0, but I don't think this does what you want it to do 
            f = zeros(size(X,2));
        else
            num = sum(A.*W.*(X-X'), 1)'; 
            den = sum(A.*W, 1)'; % does a node relate to itself? 
            f = num./den;
            f(isnan(f)) = 0; % set nan values to 0 
        end
    end
end
