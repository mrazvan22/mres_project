% Script to convert Marc's Jabocbian maps into regionally averaged
% Jacobians using Freesurfer segmentation
% You need to add spm5 into your path (can be found at:
% /cs/research/vision/camino1/camino/fonteijn/forDannyAndHubert/spm5

clear all
close all

dir_work = spm_select(1, 'dir', 'Give data dir');
% Controls are cxxtxx, patients are pxxtxx, you can find both in the
% freesurfer directory
file_data = spm_select([1 100], 'image', 'Give baseline data');
nr_subj = size(file_data, 1);
 
% FreesurferColorLUT contains the Freesurfer segment names and their values 
file_lut = [dir_work '/FreeSurferColorLUT.txt'];
[label_id, label_name, d1,d2, d3, d4] = ...
    textread(file_lut, '%d%s%d%d%d%d');
nr_labels = length(label_id);
   
[p, f, d1, d2] = fileparts(deblank(file_data(1, :)));
I1 = strfind(f, 'patient');
I2 = strfind(f, 'grp');
id_patient = str2num(f((I1+7):(I2-2)));
file_segm_nii = sprintf('%sfreesurfer/%s/mri/aparc+aseg.nii', ...
    dir_work, f);
if ~exist(file_segm_nii, 'file'),

    file_segm_mgz = sprintf('%sfreesurfer/%s/mri/aparc+aseg.mgz', ...
        dir_work, f);
    unix(sprintf('mri_convert %s %s', file_segm_mgz, file_segm_nii));

end
V_segm = spm_vol(file_segm_nii);
img_segm = spm_read_vols(V_segm);
for label = 1:nr_labels,

    I_label{label} = find(img_segm == label_id(label));

end

roi_select = [11 12 13 1001:1034 50 51 52 2001:2034]; % Freesurfer indicators that are hippocampus (17 = left, 53 = right) and the cortical areas
nr_roi_select = length(roi_select);
clear roi_id*
for id = 1:nr_roi_select,
    
    roi_id(id) = find(label_id == roi_select(id));
    name_roi{id} = label_name{roi_id(id)};
    I_label_select{id} = I_label{roi_id(id)};
    
end

[data_jacc_fluid, data_jacc_ffd] = ...
    hd_network(file_data, dir_work, roi_id);
id_group = zeros(nr_subj, 1);
for subj = 1:nr_subj,
    
    [p, f, d1, d2] = fileparts(deblank(file_data(subj, :)));
    I = strfind(f, 'grp');
    id_group(subj) = str2num(f(I+3));
    
end

data_struct.data_jacc_ffd = data_jacc_ffd;
data_struct.data_jacc_fluid = data_jacc_fluid;
data_struct.roi_id = roi_id;
data_struct.name_roi = name_roi;
data_struct.roi_select = roi_select;
data_struct.I_label_select = I_label_select;
data_struct.V_segm = V_segm;
data_struct.id_group = id_group;

save data_hd data_struct