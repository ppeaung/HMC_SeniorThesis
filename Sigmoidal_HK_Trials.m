%% Trial Sigmoidal 
%
% Sigmoidal Bounded-Confidence Model
%
% PP 2022
%
% This script will simulate a sigmoidal model on a network. 
% This simulation uses the 4th order Runge-Kutta's method in solving DE's here, using a fixed number of 
% iterations. Produces a time series of opinion states and a visualization of the graph + opinion states at convergence and the
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

gamma = 2.5;
delta = 1; % works for $delta = 1$ but doesn't match intuition with $delta = 0$ here? Or will it not exactly be the same... 
numTrials = 2;
numIterations = 5000; % from sigmoidal 
Zempty = int16.empty; 

% Singleton 

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


%% K_100 
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

%% Erdos Renyi 
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
