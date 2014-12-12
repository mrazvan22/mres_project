function robustnessTest2PCADemo

%This function tests the robustness of the algoirthmn to the outliers in
%the control cohort. i.e. how consistent is the ordering of the events
%given that there may be some outliners in control i.e. patients classified
%as control

% The idea is to repeat the model fitting with various subsets of controls gradually eliminating those closest to the patient
%distribution

%To do this I will 
% - Do PCA on control and patient data seperately
% - Sort the controls by their didtance to the patient distribution (closest first)
% - remove the first n% of cntrols from the sprted control group
% - run the experiment
% - compare the ordering
% - repeat the expermint several times

%% ============== Initialize ==============
close all
clear

%add path
addpath('/Users/jaghajan/Documents/Experiments/ADNI_BIOMARKERS_Svn/netlab')

version_likelihood=8;
div_size=[0.1:0.1:0.9];

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


%% ============= DO PCA on control data ==============

% subtract the mean
[coefs,scores,variances,t2] = princomp(data_controls);
X_c_transformed_low_dim=scores(:,1:2);
X_c_transformed_low_dim=X_c_transformed_low_dim';

plot(X_c_transformed_low_dim(1,:),X_c_transformed_low_dim(2,:),'gx')
hold on

%% ============= DO PCA on patient data ==============

% subtract the mean
[coefs,scores,variances,t2] = princomp(data_patients);
X_p_transformed_low_dim=scores(:,1:2);
X_p_transformed_low_dim=X_p_transformed_low_dim';

plot(X_p_transformed_low_dim(1,:),X_p_transformed_low_dim(2,:),'ro')

%% ======== find Mahalanobis distacne of each control point from mean of patient distributon in PCA space

mean_p=mean(X_p_transformed_low_dim,2);
cov_p=cov(X_p_transformed_low_dim');
for c_control=1:n_control
    this_control=X_c_transformed_low_dim(:,c_control);
    mahal_dist(c_control)=(this_control-mean_p)'*inv(cov_p)*(this_control-mean_p);
end
i

[s_control s_indx]=sort(mahal_dist);

% Visulize the closest control points
for c_control=1:n_control
    plot(X_c_transformed_low_dim(1,s_indx(c_control)),X_c_transformed_low_dim(2,s_indx(c_control)),'bx')
    pause
end

%% ======= fit the model to subset of controls
for c_size=1:length(div_size)
    c_size
    this_size=div_size(c_size);
    rmv_size=round(n_control*this_size);
    this_control_indx=s_indx(rmv_size+1:end);
    this_data_controls=data_controls(this_control_indx,:);
    this_data_patients=data_patients;
    
    %fit the model
    [likelihood, gmix_struct] = ...
        EBDPComputeLikelihood(this_data_patients, this_data_controls, version_likelihood);
    
    %Performing mcmc
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
    
    
    save(['Mat_Data/test2_PCA_results_versionLike' num2str(version_likelihood) 'Portion' num2str(this_size) '.mat'],'this_data_patients','this_data_controls','likelihood', 'gmix_struct','parm_struct','this_control_indx');
    
    
end



