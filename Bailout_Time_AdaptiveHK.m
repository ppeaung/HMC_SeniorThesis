%% Trial runs of Adaptive 
%
% Adaptive Hegselmann--Krause Bounded-Confidence Model
%
% PP 2021
%
% This script will simulate an adaptive Hegselmann--Krause (HK) model on a
% network generated with a configuration model. The script terminates after
% a fixed number of time steps. Produces a time series of opinion states
% and a visualization of the graph + opinion states at convergence and the
% initial graph.
% 
% Parameters:   timeSteps, number of time steps for when the model is run
%               numTrials, number of trials we iterate for some C and Beta 
%               CONFIDENCE_BOUND_RANGE, max range of C paramter
%               BETA_BOUND_RANGE, max value of Beta parameter
%               N, number of nodes that will appear in the network
%
%

%% Global variables 
% TimeSteps = 50;
clear
clc

rng shuffle 

NUM_TRIALS = 10; % 50 as default 
C_SPACING = 0.1;
B_SPACING = 0.1;
A_SPACING = 0.1; 
C_START_PARAMETER = 0; 
B_START_PARAMETER = 0; 
A_START_PARAMETER = 0;
CONFIDENCE_BOUND_RANGE = 1;
BETA_BOUND_RANGE = 1;
A_BOUND_RANGE = 1;
BAILOUT_TIME = 500;
CURR_TRIAL = 0;
p = 0.6; % probability for erdos renyi model 
T = 100; % guess of num trials; parameterizes guess of opinion matrix
BEG_NUMNODES = 100;
MAX_NUMNODES = 100;
NUM_NODES_SPACING = 100;
GRAPH_TYPE = 'erdos_renyi'; % 'random_graph', 'complete_graph', 'cycle'
TEST_INIT_GRAPH = false;
TEST_CONV_GRAPH = false;

% advance by rows then columns 


meanNumComponentMatrix = nan([fix((CONFIDENCE_BOUND_RANGE - C_START_PARAMETER) ...
    / C_SPACING) + 1, fix((BETA_BOUND_RANGE - B_START_PARAMETER) / B_SPACING) ...
    + 1, fix((A_BOUND_RANGE - A_START_PARAMETER) / A_SPACING) ...
    + 1]); 
meanConvergenceMatrix = nan([fix((CONFIDENCE_BOUND_RANGE - C_START_PARAMETER) ...
    / C_SPACING) + 1, fix((BETA_BOUND_RANGE - B_START_PARAMETER) / B_SPACING) ...
    + 1, fix((A_BOUND_RANGE - A_START_PARAMETER) / A_SPACING) ...
    + 1]);
meanDifferenceMetricMatrix = nan([fix((CONFIDENCE_BOUND_RANGE - C_START_PARAMETER) ...
    / C_SPACING) + 1, fix((BETA_BOUND_RANGE - B_START_PARAMETER) / B_SPACING) ...
    + 1, fix((A_BOUND_RANGE - A_START_PARAMETER) / A_SPACING) ...
    + 1]);
meanAssortativityCoefficientMatrix = nan([fix((CONFIDENCE_BOUND_RANGE - C_START_PARAMETER) ...
    / C_SPACING) + 1, fix((BETA_BOUND_RANGE - B_START_PARAMETER) / B_SPACING) ...
    + 1, fix((A_BOUND_RANGE - A_START_PARAMETER) / A_SPACING) ...
    + 1]);

for i = C_START_PARAMETER:C_SPACING:CONFIDENCE_BOUND_RANGE 
    for j = B_START_PARAMETER:B_SPACING: BETA_BOUND_RANGE
        for k = BEG_NUMNODES:NUM_NODES_SPACING:MAX_NUMNODES
            for l = A_START_PARAMETER:A_SPACING:A_BOUND_RANGE
            N = k;
            tollMeanComponents = 0; % reset counters
            tollMeanConvergence = 0;
            tollMeanDifferenceMetric = 0;
            tollMeanAssort = 0;
            for t = 1: NUM_TRIALS
                % Initalize the model
                if strcmp(GRAPH_TYPE, 'random_graph')
                    AdjacencyMatrix = random_graph(N);
                elseif strcmp(GRAPH_TYPE, 'complete_graph')
                    AdjacencyMatrix = complete_graph(N);
                elseif strcmp(GRAPH_TYPE, 'cycle')
                    AdjacencyMatrix = cycle(N);
                elseif strcmp(GRAPH_TYPE, 'square_lattice')
                    l = 10;
                    AdjacencyMatrix = square_lattice(l); 
                    % AdjacencyMatrix = square_lattice(...);
                elseif strcmp(GRAPH_TYPE, 'k_partite')
                    % AdjacencyMatrix = k_partite(...);
                elseif strcmp(GRAPH_TYPE, 'erdos_renyi')
                    AdjacencyMatrix = erdos_renyi(N,p);
                elseif strcmp(GRAPH_TYPE, 'barbell')
                    AdjacencyMatrix = barbell(N);
                elseif strcmp(GRAPH_TYPE, 'shortcut')
                    % AdjacencyMatrix = shortcut(...);
                elseif strcmp(GRAPH_TYPE, 'star')
                    AdjacencyMatrix = star(N);
                end
                
                if TEST_INIT_GRAPH == true
                    X(:, 1) = rand(N, 1); 
                    figure (10) 
                    plot(graph(AdjacencyMatrix), 'MarkerSize', 5, 'LineWidth', 1, 'EdgeColor', 'k', 'NodeLabel', [], 'NodeCData', X(:, 1))
                end 
            
                [X_1,Y_1, CONVERGED_TIME] = AdaptiveHK_Bailout(AdjacencyMatrix, i, j, BAILOUT_TIME, T, l);
           
                X_1 = X_1(:, 1:CONVERGED_TIME); % Remove nan values
            
                if TEST_CONV_GRAPH == true
                    figure(11)
                    plot(graph(Y_1), 'MarkerSize', 5, 'LineWidth', 1, 'EdgeColor', 'k', 'NodeLabel', [], 'NodeCData', X_1(:,end))
                    c = colorbar;
                    c.FontSize = 20; c.Label.String = 'Opinion state at convergence';
                    caxis([0 1])
                end 

                numComponents = spectral_gap(Y_1);

                tollMeanComponents = numComponents + tollMeanComponents;
                tollMeanConvergence = CONVERGED_TIME + tollMeanConvergence; 

                % Calculate how different the opinions are at the end 
                x = X_1(:, CONVERGED_TIME);
                % Round to nearest dec place 
                x = round(x, 2);
                difference_metric = numComponents / length(unique(x));

                tollMeanDifferenceMetric = difference_metric + tollMeanDifferenceMetric;


                % Calculuate the assort coefficient
                % A is vector of opinions 
                % Y_1 is the adjacency matrix that was converged  
            
                % Testing: It should be in [-1, 1] 
                % Perfectly assortive = 1 perfectly dissassortive = -1

                %%{
                degree_vec = sum(Y_1); %entry i represents degree for node
                degree_vec = degree_vec'; 
                edge_count = sum(Y_1, 'all');
                 
                if edge_count == 0
                    assort_coef = 0;
                else
                    num_mat1 = Y_1.*(x*x'); 
                    num_mat2 = ((degree_vec*degree_vec').*(x*x'))/edge_count;
                    den_mat1 = degree_vec.*(x.^2);
                    a = num_mat1 - num_mat2;
                    num = sum(a, 'all');
                    den = sum(den_mat1, 'all') - sum(num_mat2, 'all');
                    if num==0 && den==0
                        assort_coef = 0;
                    else 
                        assort_coef = num/den;
                    end
                 end
                tollMeanAssort = assort_coef + tollMeanAssort;

            end 
            meanNumComponentMatrix(round((i-C_START_PARAMETER) / C_SPACING) + 1, round((j-B_START_PARAMETER) / B_SPACING) + 1, round((l-A_START_PARAMETER) / A_SPACING) + 1) = tollMeanComponents / NUM_TRIALS;
            meanConvergenceMatrix(round((i-C_START_PARAMETER) / C_SPACING) + 1, round((j-B_START_PARAMETER) / B_SPACING) + 1, round((l-A_START_PARAMETER) / A_SPACING) + 1) = tollMeanConvergence / NUM_TRIALS;
            meanDifferenceMetricMatrix(round((i-C_START_PARAMETER) / C_SPACING) + 1, round((j-B_START_PARAMETER) / B_SPACING) + 1, round((l-A_START_PARAMETER) / A_SPACING) + 1) = tollMeanDifferenceMetric / NUM_TRIALS;
            meanAssortativityCoefficientMatrix(round((i-C_START_PARAMETER) / C_SPACING) + 1, round((j-B_START_PARAMETER) / B_SPACING) + 1, round((l-A_START_PARAMETER) / A_SPACING) + 1) = tollMeanAssort / NUM_TRIALS;
            end
        end
    end
end

%% Plots

% set up for plotx = linspace(0,1);
c = zeros(1, fix((CONFIDENCE_BOUND_RANGE - C_START_PARAMETER)/C_SPACING)+1);

for i=1:(CONFIDENCE_BOUND_RANGE/C_SPACING)+1
    c(i) = (i-1)* fix((CONFIDENCE_BOUND_RANGE - C_START_PARAMETER)/C_SPACING);
end

% plot average converged time
%{
figure(2)
set(0,'defaultaxeslinestyleorder',{'-*',':o','-*'})
hold on
for Beta = 2:size(meanNumComponentMatrix, 2)
    txt = ['Beta = ',num2str(SPACING*(Beta-1))];
    plot(c, meanNumComponentMatrix(:,Beta), 'DisplayName',txt);
end
lgd = legend; 
lgd.NumColumns = 3;
xlabel('Confidence Bound C'); ylabel('Average Num Clusters')
hold off
%}

% plot convergence time 
%{
figure(3)
set(0,'defaultaxeslinestyleorder',{'-*',':o','-*'})
hold on
for Beta = 2:size(meanNumComponentMatrix, 2)
    txt = ['Beta = ',num2str(SPACING*(Beta-1))];
    plot(c, meanConvergenceMatrix(:,Beta), 'DisplayName',txt);
end
lgd = legend; 
lgd.NumColumns = 3;
xlabel('Confidence Bound C'); ylabel('Average Convergence Time')
hold off
%}
modifiedMeanNumComponentMatrix = meanNumComponentMatrix(:,(2:11), :);

prob_var = num2str(p);
prob_var = strrep(prob_var, '.', '');
c_space_var = num2str(CONFIDENCE_BOUND_RANGE);
c_space_var = strrep(c_space_var, '.', '');
b_space_var = num2str(B_SPACING);
b_space_var = strrep(b_space_var, '.', '');




savedir_data = '/Users/loaner/Desktop/Fall 2021/MATLAB/MATH 197/_Data';
savedir_figures = '/Users/loaner/Desktop/Fall 2021/MATLAB/MATH 197/_Figures';

%if strcmp(GRAPH_TYPE, 'erdos_renyi')
%    save(fullfile(savedir_data, append(GRAPH_TYPE, ' p=', prob_var, ' N=', ...
%        int2str(MAX_NUMNODES), ...
%        ' C_spacing =', c_space_var, ' B_spacing =', ...
%        b_space_var)),"modifiedMeanNumComponentMatrix", ...
%        "meanNumComponentMatrix", "meanConvergenceMatrix", "BEG_NUMNODES", ...
%        "MAX_NUMNODES", "NUM_NODES_SPACING", "NUM_TRIALS", ...
%        "C_START_PARAMETER", "C_SPACING", "B_START_PARAMETER", ...
%        "B_SPACING", "meanDifferenceMetricMatrix");
%else 
%    save(fullfile(savedir_data, append(GRAPH_TYPE, ' N=', int2str(MAX_NUMNODES), ...
%        ' C_spacing =', c_space_var, ' B_spacing =', ...
%        b_space_var)),"modifiedMeanNumComponentMatrix", ...
%       "meanNumComponentMatrix", "meanConvergenceMatrix", ...
%        "BEG_NUMNODES", "MAX_NUMNODES", "NUM_NODES_SPACING", ...
%        "NUM_TRIALS", "C_START_PARAMETER", "C_SPACING", ...
%        "B_START_PARAMETER", "B_SPACING", "meanDifferenceMetricMatrix");
%end

%need to check what im saving to make sure it's fine 


% plot components wrt heatmap
%{
figure(4)
heatmap(modifiedMeanNumComponentMatrix(:,:,1));
colormap winter;
title('Average Number of Clusters');
ylabel('Confidence Bound C'); xlabel('Homophilly Parameter \beta');
%}

% plot convergence wrt imagesc

figure(5)
colormap default;
f = imagesc(meanConvergenceMatrix(:,:,1));
title('Convergence Time for varying parameters');
ylabel('Confidence Bound C'); xlabel('Homophilly Parameter \beta');
%xticklabels({'0', '0.03', '0.06', '0.09', '0.12', '0.15', '0.18', '0.21', '0.24', '0.27', '.3'});
%yticklabels({'0', '0.03', '0.06', '0.09', '0.12', '0.15', '0.18', '0.21', '0.24', '0.27', '.3'});
xticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'});
yticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'})
h = colorbar;
ylabel(h, 'Timestep for convergence');
 saveas(h, fullfile(savedir_figures, append(GRAPH_TYPE, ' p=', prob_var, ' N=', ...
        int2str(MAX_NUMNODES), ...
        ' C_spacing =', int2str(C_SPACING), ' B_spacing =', ...
        int2str(B_SPACING), ' Av_Conv_Time')), 'epsc');



%{
figure(6)
heatmap(meanConvergenceMatrix(:,:,1));
%}
% Using this one in thesis 

figure(7)
colormap default;
g = imagesc(modifiedMeanNumComponentMatrix(:,:,1));
title('Average Number of Clusters');
ylabel('Confidence Bound C'); xlabel('Homophilly Parameter \beta');
%xticklabels({'0.03', '0.06', '0.09', '0.12', '0.15', '0.18', '0.21', '0.24', '0.27', '.3'});
%yticklabels({'0', '0.03', '0.06', '0.09', '0.12', '0.15', '0.18', '0.21', '0.24', '0.27', '.3'});
xticklabels({'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'});
yticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'})
h = colorbar;
ylabel(h, 'Number of Clusters');
saveas(g, fullfile(savedir_figures, append(GRAPH_TYPE, ' p=', prob_var, ' N=', ...
        int2str(MAX_NUMNODES), ...
        ' C_spacing =', int2str(C_SPACING), ' B_spacing =', ...
        int2str(B_SPACING), ' Av_Num_Cl')), 'epsc');


figure(8)
colormap default;
j = imagesc(meanDifferenceMetricMatrix(:,2:end, 1));
title('Difference Metric for varying parameters');
ylabel('Confidence Bound C'); xlabel('Homophilly Parameter \beta');
xticklabels({'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'});
yticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'});
%xticklabels({'0.03', '0.06', '0.09', '0.12', '0.15', '0.18', '0.21', '0.24', '0.27', '.3'});
%yticklabels({'0', '0.03', '0.06', '0.09', '0.12', '0.15', '0.18', '0.21', '0.24', '0.27', '.3'});
h = colorbar;
ylabel(h, 'Difference Metric Value');
saveas(j, fullfile(savedir_figures, append(GRAPH_TYPE, ' p=', prob_var, ' N=', ...
        int2str(MAX_NUMNODES), ...
        ' C_spacing =', int2str(C_SPACING), ' B_spacing =', ...
        int2str(B_SPACING), 'mean_Difference_Coef')), 'epsc');

figure(9)
colormap default;
k = imagesc(meanAssortativityCoefficientMatrix(:,2:end,1));
title('Average Scalar Assortativity Coefficient for varying parameters');
ylabel('Confidence Bound C'); xlabel('Homophilly Parameter \beta');
xticklabels({'0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'});
yticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1'});
%xticklabels({'0.03', '0.06', '0.09', '0.12', '0.15', '0.18', '0.21', '0.24', '0.27', '.3'});
%yticklabels({'0', '0.03', '0.06', '0.09', '0.12', '0.15', '0.18', '0.21', '0.24', '0.27', '.3'});
h = colorbar;
ylabel(h, 'Scalar Assortitively Coefficient');


saveas(k, fullfile(savedir_figures, append(GRAPH_TYPE, ' p=', prob_var, ' N=', ...
        int2str(MAX_NUMNODES), ...
        ' C_spacing =', int2str(C_SPACING), ' B_spacing =', ...
        int2str(B_SPACING), 'Av_Assort_Coef')), 'epsc');
