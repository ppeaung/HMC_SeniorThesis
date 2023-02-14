%% Trial runs of Sigmoidal 
%
% Sigmoidal Bounded-Confidence Model
%
% PP 2022
%
% This script will simulate a sigmoidal Hegselmann model on a
% network generated with a configuration model. This simulation uses 
% Runge-Kutta's method in solving DE's here, using a fixed number of 
% iterations. Produces a time series of opinion states
% and a visualization of the graph + opinion states at convergence and the
% initial graph.
% 
% Note: this requires sigmoidal.m
% 
% Parameters:   
%
%% Global variables 
% TimeSteps = 50;
%

clear 
clc

rng shuffle 

%% Test cases on complete graphs with gamma --> inf (HK model): 
% 1) Singleton graph 
% 2) K_2 
% 3) K_3 
% 4) K_100

% To do's:
% 1) Check to see if your implementation truly works
    % it looks functional or aligns with intuitions with gamma \to \infty
    % and delta = 1 but what about delta = 0? 
    % Check with Karate Club implementation and compare your graph with
    % what pc, hzb, and mp did on theirs 
% 2) reread forecast election paper that prof. volkening did 
% 3) reask prof. deford about how to get the parameters again and polling
% data... 
% 4) could check to see how the model performs under various graph
% topologies again... 
% 5) mention of stochasticity here! 


% doing weird things, such as exceeding opinion value of 1 or being
% negative 

gamma = 2.5;
delta = 1; % works for $delta = 1$ but doesn't match intuition with $delta = 0$ here? Or will it not exactly be the same... 
numTrials = 2;
numIterations = 5000; % from sigmoidal 
Zempty = int16.empty; 

% to do: check with karate club data set and compare with pc, hzb, mp

% Singleton 
%{
adjSingle = ones(1); 
opinSingle = rand(1, numTrials); 
dataSingle = nan(size(opinSingle,1), numIterations, numTrials); 
for n = 1:numTrials
    Y = sigmoidal(adjSingle, opinSingle(n), delta, gamma, Zempty); % matrix of size size(opinSingle) x numIterations 
    dataSingle(:,:,n) = Y; 
end

meanSingle = mean(dataSingle, 'all'); 

    figure(1)
    x = linspace(1,numIterations, numIterations);
    y = dataSingle(1,:, 1);  %arbitrarily, the first trial  
    plot(x, y)
%}

%{
% K2
adjK2 = ones(2);
adjK2 = adjK2 - eye(2); 
% opinK2 = rand(2, numTrials); since random values aren't working rn, do
% fixed ones 

opinK2 = [0, 0; 1, 1]; % if opinions aren't the same, things go weird. If they're same, invariant 
dataK2 = nan(size(opinK2,1), numIterations, numTrials); 
for n = 1:numTrials
    Y = sigmoidal(adjK2, opinK2(:, n), delta, gamma, Zempty); % matrix of size size(opinSingle) x numIterations 
    dataK2(:,:,n) = Y; % need sanity checks to make sure values stay within bounds 
end


meanK2 = mean(dataK2, [1 2]); 

    figure(2)
    x = linspace(1,numIterations, numIterations);
    y = dataK2(1,:, 1);  %arbitrarily, the first trial  
    plot(x, y)

    hold on 
    plot(x, dataK2(2,:,1));
    hold off 
%}

%% K_3 
%{
adjK3 = ones(3);
adjK3 = adjK3 - eye(3); 
opinK3 = rand(3, numTrials); 
dataK3 = nan(size(opinK3,1), numIterations, numTrials); 
for n = 1:numTrials
    Y = sigmoidal(adjK3, opinK3(:, n), delta, gamma, Zempty); % matrix of size size(opinSingle) x numIterations 
    dataK3(:,:,n) = Y; 
end

    figure(3)
    x = linspace(1,numIterations, numIterations);
    y = dataK3(1,:, 1);  %arbitrarily, the first trial  
    plot(x, y)

    hold on 
    for i=2:3
        plot(x, dataK3(i,:,1));
    end
    hold off 
meanK3 = mean(dataK3, [1 2]); 
%}


%% K_100 
%{
adjK100 = ones(100);
adjK100 = adjK100 - eye(100); 
opinK100 = rand(100, numTrials); 
dataK100 = nan(size(opinK100,1), numIterations, numTrials); 
for n = 1:numTrials
    Y = sigmoidal(adjK100, opinK100(:, n), delta, gamma, Zempty); %matrix of size size(opinSingle) x numIterations 
    dataK100(:,:,n) = Y; 
end

meanK100 = mean(dataK100, [1 2]); 

    figure(100)
    x = linspace(1,numIterations, numIterations);
    y = dataK100(1,:, 1);  %arbitrarily, the first trial  
    plot(x, y)

    hold on 
    for i=2:100
        plot(x, dataK100(i,:,1));
    end
    hold off 
%}

%% Erdos Renyi 
%{
j = 100;
p = 0.3; 

adjER = rand(j);
probMatrix = ones(j)*p;
adjER = adjER < probMatrix;
adjER = adjER - diag(diag(adjER)); % Remove self loops
adjER = triu(adjER);
adjER = adjER + adjER'; 

opinER = rand(j, numTrials); 
dataER = nan(size(opinER,1), numIterations, numTrials);

for n = 1:numTrials
    Y = sigmoidal(adjER, opinER(:, n), delta, gamma, Zempty); %matrix of size size(opinSingle) x numIterations 
    dataER(:,:,n) = Y; 
end

figure(1000)
    plot(graph(adjER), 'MarkerSize', 5, 'LineWidth', 1, 'EdgeColor', 'k', 'NodeLabel', [], 'NodeCData', dataER(:, 1, 1));
    c = colorbar;
    c.FontSize = 20; c.Label.String = 'Opinion state at convergence';
    caxis([0 1])    

figure(1001)
    plot(graph(adjER), 'MarkerSize', 5, 'LineWidth', 1, 'EdgeColor', 'k', 'NodeLabel', [], 'NodeCData', dataER(:, end, 1));
    c = colorbar;
    c.FontSize = 20; c.Label.String = 'Opinion state at convergence';
    caxis([0 1])    

figure(101)
    x = linspace(1,numIterations, numIterations);
    y = dataER(1,:, 1);  %arbitrarily, the first trial  
    plot(x, y)

    hold on 
    for i=2:j
        plot(x, dataER(i,:,1));
    end
    hold off 
%}

%% Karate Club 
karate = load('karate.mat');
adjK = karate.Problem.A; 
opinKarate = zeros(size(full(karate.Problem.A), 1), 1); 
opinKarate(1) = 1; 
opinKarate(34) = -1; 
ZKarate = [1; 34]; 
dataKarate = nan(size(opinKarate,1), numIterations, numTrials); 

for n = 1:numTrials
    Y = sigmoidal(full(adjK), opinKarate, delta, gamma, ZKarate); %matrix of size size(opinSingle) x numIterations 
    dataKarate(:,:,n) = Y; 
end

% node 1 and 34 are zealots here 


figure(2000)
    plot(graph(full(adjK)), 'MarkerSize', 5, 'LineWidth', 1, 'EdgeColor', 'k', 'NodeLabel', [], 'NodeCData', dataKarate(:,end,1));
    c = colorbar;
    c.FontSize = 20; c.Label.String = 'Opinion state at convergence';
    caxis([-1 1])   

figure(2001)
    x = linspace(1,numIterations, numIterations);
    y = dataKarate(1,:, 1);  %arbitrarily, the first trial  
    plot(x, y)

    hold on 
    for i=2:size(opinKarate, 1)
        plot(x, dataKarate(i,:,1));
    end
    hold off 
