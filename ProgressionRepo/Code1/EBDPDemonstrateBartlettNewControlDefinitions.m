% Script to explain the event-based disease progression model when applied
% to Jonathan 

clear all
close all
clc

% First setting up paths
% All functions are in repos/Code1 and repos/Code1/netlab

addpath /cs/research/vision/camino1/camino/fonteijn/Progression/Code1
addpath /cs/research/vision/camino1/camino/fonteijn/Progression/Code1/netlab

% Reading in Jonathan's data
file_data = '/cs/research/vision/camino1/camino/fonteijn/Progression/BartlettData/2011_06_09ADNINewClassificationsXLS98.xls';
[data_num, data_txt, data_raw] = xlsread(file_data);
name_variables = data_txt(1, 2:(end-1));
% Splitting off real data from control/patient identifiers which are in the
% last column
id_controls = data_num(:, end);
data = data_num(:, 2:(end-1));

% Throwing out all patients which contain missing data
I_missing = find(sum(isnan(data), 2));
data(I_missing, :) = [];
id_controls(I_missing) = [];

% Not sure what to do with the CDR scores and avtot, so I'll throw them out for now
I_keep = [1:7];
nr_events = length(I_keep);
data = data(:, I_keep);
name_variables = name_variables(I_keep);

[nr_subj, nr_var] = size(data);
data_controls = data(id_controls == 1, :);
data_patients = data(id_controls == 0, :);

% Checking the distributions of all variables that remain in
% controls(green) and patients (red)

figure(1), clf
for i_var = 1:nr_var,
    
    subplot(3, 4, i_var), hold on
    [hist_c, x_c] = ksdensity(data_controls(:, i_var));
    plot(x_c, hist_c, 'g');
    [hist_p, x_p] = ksdensity(data_patients(:, i_var));
    plot(x_p, hist_p, 'r');
    title(name_variables{i_var})
    
end

% Visual inspection of the histograms shows that some variables are higher
% in patients than in controls. I'm going to keep the convention in the
% modelig that values should be lower in patients than in controls

I_swap = [1 3:7];
data_controls(:, I_swap) = -data_controls(:, I_swap);
data_patients(:, I_swap) = -data_patients(:, I_swap);
I_log_pos = [];
data_controls(:, I_log_pos) = log(data_controls(:, I_log_pos));
data_patients(:, I_log_pos) = log(data_patients(:, I_log_pos));
I_log_neg = [1 3];
data_controls(:, I_log_neg) = -log(-data_controls(:, I_log_neg));
data_patients(:, I_log_neg) = -log(-data_patients(:, I_log_neg));

figure(2), clf
for i_var = 1:nr_var,
    
    subplot(3, 4, i_var), hold on
    [hist_c, x_c] = ksdensity(data_controls(:, i_var));
    plot(x_c, hist_c, 'g');
    [hist_p, x_p] = ksdensity(data_patients(:, i_var));
    plot(x_p, hist_p, 'r');
    title(name_variables{i_var})
    
end

% Now fitting the EBDP model to the data
% We only need two functions for this: EBDPComputeLikelihood.m and EBDPMCMC.m

% First computing the likelihood of the values given events have happened
% or not
% The histograms suggest that a mixture of Gaussians would do quite nicely
% First trying out option 1, which fits a single Gaussian to the control
% distribution, then fits a mixture to the whole data distribution while
% keeping the control distribution fixed

version_likelihood = 4;
[likelihood, gmix_struct] = ...
    EBDPComputeLikelihood(data_patients, data_controls, version_likelihood);

% Checking how the mixture distributions fit the real distributions
% Waling through one variable at a time
% The left plot is the true data distribution, the right plot is the
% mixture likelihood
for i_var = 1:nr_var,
    
    figure(1), clf
    p_controls = gmmprob(gmix_struct.gmix_patients{i_var}, data_controls(:, i_var));
    p_patients = gmmprob(gmix_struct.gmix_patients{i_var}, data_patients(:, i_var));
    [hist_c, x_c] = ksdensity(data_controls(:, i_var));
    [hist_p, x_p] = ksdensity(data_patients(:, i_var));
    
    subplot(121), hold on
    plot(x_c, hist_c, 'g');
    plot(x_p, hist_p, 'r');
    
    subplot(122), hold on
    plot(data_controls(:, i_var), p_controls, 'g.')
    plot(data_patients(:, i_var), p_patients, 'r.');
    
    fprintf('%s\n', name_variables{i_var})
    pause
    
end

%% Performing MCMC

% Performing mcmc
parm_mcmc.nr_gradient_ascent = 5;
parm_mcmc.nr_it_gradient_ascent = 2e3;
parm_mcmc.nr_it_burnin = 1e5;
parm_mcmc.nr_it_mcmc = 1e5;
parm_mcmc.interval_display = 1e2;
parm_mcmc.flag_sum = 2; % This should always be set to 2.
parm_mcmc.nr_it_check_p_false = 1e2;
parm_mcmc.range_p_false = [0.05 0.1
    0.05 0.1];
parm_mcmc.std_p_false = [1e-4 1e-4]';
parm_mcmc.idx_clinevent = []; 

parm_struct = EBDPMCMCTest(likelihood, ...
    parm_mcmc);

% Computing positional variance
hist2_mat = zeros([nr_events nr_events]);
[ord inds_av] = sort(mean(parm_struct.event_order_mcmc, 2));
for roi1 = 1:nr_events,
    
    for roi2 = 1:nr_events,
        
        hist2_mat(roi1, roi2) = ...
            sum(parm_struct.event_order_mcmc(inds_av(roi1), :) == roi2);
        
    end
    
end

% Visualizing positional variance...
position_roi_mean = mean(parm_struct.event_order_mcmc, 2);
position_roi_std = std(parm_struct.event_order_mcmc, [], 2);

[d, order_roi] = sort(position_roi_mean, 'ascend');
% Making a annotated version of the histogram
max_length_str = 0;
for roi = 1:nr_events,
    
    max_length_str = max([max_length_str length(name_variables{roi})]);
    
end
mat_labels = zeros(nr_events, max_length_str);
for roi = 1:nr_events,
    
    length_str = length(name_variables{order_roi(roi)});
    mat_labels(roi, 1:length_str) = name_variables{order_roi(roi)};
    
end
figure(3), clf, 
imagesc(hist2_mat)
map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
colormap(map_inverted)
axis square
set(gca, ...
    'YTick', [1:nr_events], 'YTickLabel', char(mat_labels), ...
    'YGrid', 'on')
hXLabel = xlabel('model stage');
hYLabel = ylabel('region');
set([hXLabel hYLabel], ...
    'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperPositionMode', 'auto');






