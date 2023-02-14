%% Inputs: 
%
% A: Adjacency matrix 
% X: array of opinions 
% delta: Confidence bound
% gamma: controls "sharpness" of dependence of nodes i and j, as gamma 
% approaches infinity it approaches a stepwise function (classical HK
% model) 
% BAILOUT_TIME: Bailout time
% T: Number of Trials to Run
%%

% using inf makes inf and nan values appear, seems using "big" values may
% be solution 
% ok so maybe implement Euler's method and then implement a fourth order
% Runge-Kutta method... 
% or maybe just use the damn ode solver in matlab 

function Y = sigmoidalHK(A, X, delta, gamma) 

% fixed number of iterations
numIterations = 50; 
Y = nan(size(X, 2), numIterations); 

Xnew = X; % initialize X 

h = 0.001; % What do I choose as h here? 


for i = 1:numIterations
    Xold = Xnew; % reset old value
    a = h*sigmoidalDE(A, weight(Xold, gamma, delta), Xold); 
    b = h*sigmoidalDE(A, weight(Xold + (1/2)*a, gamma, delta), Xold + (1/2)*a); 
    c = h*sigmoidalDE(A, weight(Xold + (1/2)*b, gamma, delta), Xold + (1/2)*b);
    d = h*sigmoidalDE(A, weight(Xold + (1/2)*c, gamma, delta), Xold + (1/2)*c);
    Xnew = Xold + (1/6)*(a + 2*b + 2*c + d);
    Y(:, i) = Xnew; 
end

% Input: vector of opinions X and parameters gamma and delta
% Returns: matrix of the corresponding weights
function W = weight(X, gamma, delta)
    W = (1 + exp(gamma*(X-X').^2 - gamma*delta)).^(-1); 
    W = W - diag(diag(W)); % remove diagonal weights 
end

% Input: adj matrix A, weight matrix W, and vector of opinions X 
% Returns: function representing DE for sigmoidal DE 
function f = sigmoidalDE(A, W, X)
    f = (sum(A.*(W.*(X-X')), 2))./sum(A.*W,2); 
end
    %end
end
