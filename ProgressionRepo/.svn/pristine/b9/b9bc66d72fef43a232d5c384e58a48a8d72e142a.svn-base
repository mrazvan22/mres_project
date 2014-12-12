% Script to convert Marc's Jabocbian maps into regionally averaged
% Jacobians using Freesurfer segmentation
% You need to add spm5 into your path (can be found at:
% /cs/research/vision/camino1/camino/fonteijn/forDannyAndHubert/spm5

clear all
close all

dir_work = spm_select(1, 'dir', 'Give data dir');
% Controls are cxxtxx, patients are pxxtxx, you can find both in the
% freesurfer directory
dir_control = spm_select([1 100], 'dir', 'Give controls');
dir_patient = spm_select([1 100], 'dir', 'Give patients');
file_diff_t = 'diff_t2baseline.mat';
file_cat_at_scan = 'cat_at_scan.mat';

for tp = 1:size(dir_patient, 1),
    
    patient_id_str(tp, :) = dir_patient(tp, (end-6):(end-1));
    patient_id(tp) = str2num(patient_id_str(tp, 2:3));
    patient_tp(tp) = str2num(patient_id_str(tp, 5:6));
    
end
for tp = 1:size(dir_control, 1),
    
    control_id_str(tp, :) = dir_control(tp, (end-6):(end-1));
    control_id(tp) = str2num(control_id_str(tp, 2:3));
    control_tp(tp) = str2num(control_id_str(tp, 5:6));
    
end

load(file_diff_t);
load(file_cat_at_scan);
diff_t_control = diff_t(1:size(dir_control, 1));
diff_t_patient = diff_t((size(dir_control, 1)+1):end);
cat_control = cat_at_scan(1:size(dir_control, 1));
cat_patient = cat_at_scan((size(dir_control, 1)+1):end);

% Removing all first time points from the list.
I_baseline_patient = find(diff_t_patient == 0);
I_baseline_control = find(diff_t_control == 0);
patient_id(I_baseline_patient) = [];
patient_tp(I_baseline_patient) = [];
control_id(I_baseline_control) = [];
control_tp(I_baseline_control) = [];
diff_t_control(I_baseline_control) = [];
diff_t_patient(I_baseline_patient) = [];
cat_control(I_baseline_control) = [];
cat_patient(I_baseline_patient) = [];

% FreesurferColorLUT contains the Freesurfer segment names and their values 
file_lut = [dir_work '/FreeSurferColorLUT.txt'];
[label_id, label_name, d1,d2, d3, d4] = ...
    textread(file_lut, '%d%s%d%d%d%d');

id_roi_select = [17 53 1001:1034 2001:2034]; % Freesurfer indicators that are hippocampus (17 = left, 53 = right) and the cortical areas
nr_roi_select = length(id_roi_select);
clear roi_id*
for roi = 1:nr_roi_select,
    
    id_roi = find(label_id == id_roi_select(roi));
    name_roi{roi} = label_name{id_roi};
    
end

for flag_reg = [1 2],
    
    [data_control{flag_reg}, data_patient{flag_reg}, I_label_select] = ...
        ad_network(dir_control, dir_patient, dir_work, id_roi_select, ...
        flag_reg);
    
end

data_struct.data_control = data_control;
data_struct.data_patient = data_patient;
data_struct.cat_patient = cat_patient;
data_struct.cat_control = cat_control;
data_struct.flagreg{1} = 'Modat';
data_struct.flagreg{2} = 'Freeborough';
data_struct.diff_t_control = diff_t_control;
data_struct.diff_t_patient = diff_t_patient;
data_struct.id_roi_select = id_roi_select;
data_struct.name_roi = name_roi;
data_struct.I_label_select = I_label_select;
data_struct.patient_id = patient_id;
data_struct.patient_tp = patient_tp;
data_struct.control_id = control_id;
data_struct.control_tp = control_tp;

save data_AD data_struct
