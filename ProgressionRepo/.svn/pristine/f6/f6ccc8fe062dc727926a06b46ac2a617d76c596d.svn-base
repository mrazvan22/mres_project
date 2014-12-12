clear all
close all

file_data = spm_select(1, 'mat');
load(file_data);
%% Create the precedence probability matrix.
atrophy_patient = data_struct.data_jacc_patient_select';
atrophy_control = data_struct.data_jacc_control_select';
I_label_select = data_struct.I_label_select;

[nr_roi nr_pat] = size(atrophy_patient);

atrophy_patient = data_struct.data_patient{1}.data_nomixt;
atrophy_control = data_struct.data_control{1}.data_nomixt;
roi_id = data_struct.roi_id;
I_label_select = data_struct.I_label_select;
name_roi = data_struct.name_roi;

atrophy_patient = atrophy_patient(roi_id, :);
atrophy_control = atrophy_control(roi_id, :);

p_atrophy = zeros(nr_roi, nr_pat);
for pat = 1:nr_pat,
    
    for roi = 1:nr_roi,
        
        [h, p_atrophy(roi, pat)] = ...
            ttest(atrophy_control(roi, :), ...
            atrophy_patient(roi, pat), [], 'left');
        
    end
    
end