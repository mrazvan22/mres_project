clear
close all
% First setting up paths
addpath('/Users/jaghajan/Documents/Experiments/ADNI_BIOMARKERS_Svn/netlab')
addpath('/Users/jaghajan/Documents/Experiments/ADNI_BIOMARKERS_Svn/Archive')

% All functions are in repos/Code1 and repos/Code1/netlab

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

[nr_subj, nr_var] = size(data);

% Continuing with the csf defined controls
data_controls = data(id_controls_csf == 1, :);
data_patients = data(id_controls_csf == 0, :);

% Checking the distributions of all variables that remain in
% controls(green) and patients (red)

figure(1), clf
for i_var = 1:nr_var
    
    subplot(3, 4, i_var), hold on
    [hist_c, x_c] = ksdensity(data_controls(:, i_var));
    plot(x_c, hist_c, 'g');
    [hist_p, x_p] = ksdensity(data_patients(:, i_var));
    plot(x_p, hist_p, 'r');
    title(name_variables{i_var})
    set(gcf,'Color',[1 1 1]);
    
end

% Visual inspection of the histograms shows that some variables are higher
% in patients than in controls. I'm going to keep the convention in the
% modelig that values should be lower in patients than in controls

%I_swap = [1 3 4 7 9 10 11];
I_swap = [1 3 4 6 8 9 10];

data_controls(:, I_swap) = -data_controls(:, I_swap);
data_patients(:, I_swap) = -data_patients(:, I_swap);
data(:,I_swap) = -data(:,I_swap);
figure(2), clf
for i_var = 1:nr_var,
    
    subplot(3, 4, i_var), hold on
    [hist_c, x_c] = ksdensity(data_controls(:, i_var));
    plot(x_c, hist_c, 'g');
    [hist_p, x_p] = ksdensity(data_patients(:, i_var));
    plot(x_p, hist_p, 'r');
    title(name_variables{i_var})
    set(gcf,'Color',[1 1 1]);
    
end

%cov type parameter
%covType = 0;  %unconstrained - individual covariance matrices
%covType = 1;  %unconstrained - shared covariance matrices
%covType = 2;  %diagonal - individual covariance matrices
%covType = 3;  %diagonal - shared covariance matrices
covType = 4;  %spherical - individual covariance matrices
%covType = 5;  %spherical - shared covariance matrices


%% Perform unsupervised clustering of biomarker data

cluster_version=3;
%1 = EM MoG on single biomarkrs
%12 = version 1.2. same as 1 but with random initialization
%2 = EM MoG on entire data
%22 = same as two but only on selected biomarkers
%3 = EM on PCA Dimensioanlity reduecd data
%4 = K-means clustering

switch cluster_version
  
         
    case 3
        % ============= PCA dimensionality reduced data =============
        [coefs,scores,variances,t2] = princomp(data);
        X_transformed_low_dim=scores(:,1:2);
        X_transformed_low_dim=X_transformed_low_dim';
        
%         
%                  figure
%                  plot3(X_transformed_low_dim(1,id_controls_csf == 1),X_transformed_low_dim(2,id_controls_csf == 1),X_transformed_low_dim(3,id_controls_csf == 1),'bx');
%                  hold on
%                  plot3(X_transformed_low_dim(1,id_controls_csf == 0),X_transformed_low_dim(2,id_controls_csf == 0),X_transformed_low_dim(2,id_controls_csf == 0),'rx');
%         
        
        
        figure
        plot(X_transformed_low_dim(2,id_controls_csf == 0),X_transformed_low_dim(1,id_controls_csf == 0),'rx');
        hold on
        plot(X_transformed_low_dim(2,id_controls_csf == 1),X_transformed_low_dim(1,id_controls_csf == 1),'gx');
        
        [label, model] = EMMixT(X_transformed_low_dim, 2,covType);
        save('Clustering_Results/MixT_PCA_Clustering_Results.mat','label','model');
        
end

%% Visualize
% csf_contorl_indx=find(id_controls_clin == 1);
% csf_patient_index=find(id_controls_clin == 0);
%
% est_class1_indx=find(label==1); %patient cluster
% est_class2_indx=find(label==2); %control cluster
%
% %% Now fit the EBDP model with new control and patient groups
% new_data_controls = data(est_class2_indx, :);
% new_data_patients = data(est_class1_indx, :);
%
%
% EBDPDemonstrateBartlettGivenData(new_data_patients,new_data_controls,name_variables)
%


