function robustnessTest2Demo(div_size,version_likelihood,n_repeat)

%This function tests the robustness of the algoirthmn to the outliers in
%the control cohort. i.e. how consistent is the ordering of the events
%given that there may be some outliners in control i.e. patients classified
%as control

%To do this I will 
% - randomely select a set of controls (equal to 1/3 of the control size)
% - remove this from the control group
% - run the experiment
% - compare the ordering
% - repeat the expermint several times

%add path
addpath('/Users/jaghajan/Documents/Experiments/ADNI_BIOMARKERS_Svn/netlab');

%% ====================== Load data ======================
% Reading in Jonathan's data
%file_data = '/cs/research/vision/camino1/camino/fonteijn/Progression/BartlettData/data_biomarkers_ADNI.mat';
file_data = 'BartlettData/data_biomarkers_ADNI.mat';
load(file_data);

name_variables = data_struct.name_var;
data = data_struct.data;
% This is somewhat cumbersome: the first variable is the subject id, which
% is actually a string and therefore not included in the data matrix. The
% data belonging to a variable with variable name name_variables{x} can 
% therefore be found in data(:, x-1)
% I'm correcting for that here
I_keep = [2:length(name_variables)];
name_variables = name_variables(I_keep);

% There are quite a lot of patients in which one of the measurements is missing.
% Throwing out all of these patients
I_missing = find(sum(isnan(data), 2));
data(I_missing, :) = [];

% There are control definitions: the last column are the clinically defined
% controls (controls == 0, patients == 1), the column before that the
% CSF-defined controls
id_controls_clin = data(:, end);
id_controls_csf = data(:, end-1);
% Splitting off these variables to avoid confusion
data(:, (end-1):end) = [];

% Not sure what to do with the CDR scores and avtot, so I'll throw them out for now
% The CDR scores are either 0, 0.5 and 1, so fitting continuous distributions to those
% is problematic ***jania:***  I also removed the mmscore as it is discrete
% and a MoG does not do a good job of fitting it
I_keep = [1:4 7:12];
data = data(:, I_keep);
name_variables = name_variables(I_keep);

% Enforcing zero mean and unit standard deviation
% data = zscore(data');
[nr_subj, nr_var] = size(data);
% Continuing with the clinically defined controls
data_controls = data(id_controls_csf == 1, :);
data_patients = data(id_controls_csf == 0, :);

% Visual inspection of the histograms shows that some variables are higher
% in patients than in controls. I'm going to keep the convention in the
% modelig that values should be lower in patients than in controls

%I_swap = [1 3 4 7 9 10 11];
I_swap = [1 3 4 6 8 9 10];

data_controls(:, I_swap) = -data_controls(:, I_swap);
data_patients(:, I_swap) = -data_patients(:, I_swap);

n_control=size(data_controls,1);
n_patient=size(data_patients,1);

% %% randomely select some patinets to be part of control
%div_size=0.3333; 
rand_size=round(n_control*div_size);
stored_rand_perm=zeros(n_repeat,n_control);

if exist(['Mat_Data/test2_stored_rand_perm' num2str(n_repeat) 'Repeat' num2str(div_size) 'Portion.mat'])
    load(['Mat_Data/test2_stored_rand_perm' num2str(n_repeat) 'Repeat' num2str(div_size) 'Portion.mat'],'stored_rand_perm');
else
    
    for i=1:n_repeat
        stored_rand_perm(i,:)=randperm(n_control);
    end
    save(['Mat_Data/test2_stored_rand_perm' num2str(n_repeat) 'Repeat' num2str(div_size) 'Portion.mat'],'stored_rand_perm');
end

% 
% %% ====================== fit the model ======================
% 
% % Now fitting the EBDP model to the data repeatedly
% % We only need two functions for this: EBDPComputeLikelihood.m and EBDPMCMC.m
% 
% for cRep=1:n_repeat
%     cRep
% this_rand_select_indx= stored_rand_perm(cRep,1:rand_size);
% this_control_indx=setdiff([1:n_control],this_rand_select_indx);
% 
% this_data_controls=data_controls(this_control_indx,:);
% this_data_patients=data_patients;
% 
% %% fit the model
% [likelihood, gmix_struct] = ...
%     EBDPComputeLikelihood(this_data_patients, this_data_controls, version_likelihood);
% 
% %Performing mcmc
% parm_mcmc.nr_gradient_ascent = 5;
% parm_mcmc.nr_it_gradient_ascent = 2e3;
% parm_mcmc.nr_it_burnin = 1e5;
% parm_mcmc.nr_it_mcmc = 1e5;
% parm_mcmc.interval_display = 1e2;
% parm_mcmc.flag_sum = 2; % This should always be set to 2.
% parm_mcmc.nr_it_check_p_false = 1e2;
% parm_mcmc.range_p_false = [0.05 0.1
%     0.05 0.1];
% parm_mcmc.std_p_false = [1e-4 1e-4]';
% parm_mcmc.idx_clinevent = [];
% 
% parm_struct = EBDPMCMCTest(likelihood, ...
%     parm_mcmc);
% 
% 
% save(['Mat_Data/test2_results_versionLike' num2str(version_likelihood) 'Rep' num2str(cRep) 'Repeat' num2str(div_size) 'Portion.mat'],'this_data_patients','this_data_controls','likelihood', 'gmix_struct','parm_struct','this_rand_select_indx','this_control_indx');
% 
% end

% %% ====================== Compare the results ======================
mcmcSize=100000;
event_order_mcmc_all=zeros(nr_var,mcmcSize*n_repeat);
position_roi_mean_all=zeros(nr_var,n_repeat);
for cRep=1:n_repeat
    cRep
    close all
    
    load(['Mat_Data/test2_results_versionLike' num2str(version_likelihood) 'Rep' num2str(cRep) 'Repeat' num2str(div_size) 'Portion.mat'],'this_data_patients','this_data_controls','likelihood', 'gmix_struct','parm_struct','this_rand_select_indx','this_control_indx');
    
    
    %% visualizing the event ordering
    % Computing positional variance
    nr_events = size(parm_struct.event_order_mcmc, 1);
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
    figure(1);
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
    set(gcf,'Color',[1 1 1])
%     %%  ====================== visualizing the patient control distribution for this random selection ======================
%     figure(2);
%     for i_var = 1:nr_var,
%         
%         subplot(3, 4, i_var), hold on
%         [hist_c, x_c] = ksdensity(this_data_controls(:, i_var));
%         plot(x_c, hist_c, 'g');
%         [hist_p, x_p] = ksdensity(this_data_patients(:, i_var));
%         plot(x_p, hist_p, 'r');
%         title(name_variables{i_var})     
%         
%     end
%     set(gcf,'Color',[1 1 1])
%     
     event_order_mcmc_all(:,(cRep-1)*mcmcSize+1:cRep*mcmcSize)=parm_struct.event_order_mcmc;
     position_roi_mean_all(:,cRep)=position_roi_mean;
%     set(gcf,'Color',[1 1 1])
%     pause
    
end
%  
% % %% ====================== Draw all of the results in one plot ======================
% % figure;
% % for cRep=1:n_repeat    
% % load(['Mat Data/test2_results_versionLike' num2str(version_likelihood) 'Rep' num2str(cRep) 'Repeat' num2str(div_size) 'Portion.mat'],'this_data_patients','this_data_controls','likelihood', 'gmix_struct','parm_struct','this_rand_select_indx','this_control_indx');
% %     
% %     nr_events = size(parm_struct.event_order_mcmc, 1);
% %     hist2_mat = zeros([nr_events nr_events]);
% %     [ord inds_av] = sort(mean(parm_struct.event_order_mcmc, 2));
% %     for roi1 = 1:nr_events,
% %         
% %         for roi2 = 1:nr_events,
% %             
% %             hist2_mat(roi1, roi2) = ...
% %                 sum(parm_struct.event_order_mcmc(inds_av(roi1), :) == roi2);
% %             
% %         end
% %         
% %     end
% %     
% %     % Visualizing positional variance...
% %     position_roi_mean = mean(parm_struct.event_order_mcmc, 2);
% %     position_roi_std = std(parm_struct.event_order_mcmc, [], 2);
% %     
% %     [d, order_roi] = sort(position_roi_mean, 'ascend');
% %     % Making a annotated version of the histogram
% %     max_length_str = 0;
% %     for roi = 1:nr_events,
% %         
% %         max_length_str = max([max_length_str length(name_variables{roi})]);
% %         
% %     end
% %     mat_labels = zeros(nr_events, max_length_str);
% %     for roi = 1:nr_events,
% %         
% %         length_str = length(name_variables{order_roi(roi)});
% %         mat_labels(roi, 1:length_str) = name_variables{order_roi(roi)};
% %         
% %     end
% %     subplot(2,ceil(n_repeat/2),cRep)
% %     imagesc(hist2_mat)
% %     hold on;
% %     map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
% %     colormap(map_inverted)
% %     axis square
% %     set(gca, ...
% %         'YTick', [1:nr_events], 'YTickLabel', char(mat_labels), ...
% %         'YGrid', 'on')
% %     hXLabel = xlabel('model stage');
% %     hYLabel = ylabel('region');
% %     set([hXLabel hYLabel], ...
% %         'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
% %     set(gcf, 'PaperPositionMode', 'auto');
% %     set(gcf,'Color',[1 1 1])
% % end
% 
% 
%% ======================Draw the resulting event ordering as an average position of the results of n repetiotions ======================

% visualizing the event ordering
%Computing positional variance
nr_events = size(event_order_mcmc_all, 1);
hist2_mat = zeros([nr_events nr_events]);
[ord inds_av] = sort(mean(event_order_mcmc_all, 2));
for roi1 = 1:nr_events,
    
    for roi2 = 1:nr_events,
        
        hist2_mat(roi1, roi2) = ...
            sum(event_order_mcmc_all(inds_av(roi1), :) == roi2);
        
    end
    
end

% Visualizing positional variance...
position_roi_mean = mean(event_order_mcmc_all, 2);
position_roi_std = std(event_order_mcmc_all, [], 2);

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
figure(4);
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
set(gcf,'Color',[1 1 1])


