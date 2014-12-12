function common_control_indx=findCommonControlIndices(three_flag,rand_flag)

close all

%cluster_version=1;
%1 = EM MoG on single biomarkrs
%2 = EM MoG on entire data
%3 = EM on PCA Dimensioanlity reduecd data
%31 = For visualization only --> compare the PCA clustering with ND All data clustering

% First setting up paths
addpath('/Users/jaghajan/Documents/Experiments/ADNI_BIOMARKERS_Svn/netlab')
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

for roi = 1:10
    if rand_flag
        load(['Clustering_Results/Clustering_Results_with_' name_variables{roi} '.mat'],'this_data','init_labels','label','model','llh','init_centres');
        
    else
        load(['Clustering_Results/Clustering_Results_with_' name_variables{roi} 'Random_init.mat'],'this_data','label','model','llh');
    end
    roi1_model=model;
    roi1_label=label;
    roi1_llh=llh;
    roi1_thisdata=this_data;
    roi1_init_centres=init_centres;
    
    roi1_patient_indx=find(roi1_label==1);
    roi1_control_indx=find(roi1_label==2);
    
    all_roi_controls{roi}=roi1_control_indx;
end

if three_flag % select the first three bio9markers
    common_control_indx=mintersect(all_roi_controls{1},all_roi_controls{2},all_roi_controls{3});
else
common_control_indx=mintersect(all_roi_controls{1},all_roi_controls{2},all_roi_controls{3}, ...
    all_roi_controls{4},all_roi_controls{5},all_roi_controls{6}, ...
    all_roi_controls{8},all_roi_controls{9});
end
