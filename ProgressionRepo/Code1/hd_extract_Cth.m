% Script to convert Marc's Jabocbian maps into regionally averaged
% Jacobians using Freesurfer segmentation
% You need to add spm5 into your path (can be found at:
% /cs/research/vision/camino1/camino/fonteijn/forDannyAndHubert/spm5

clear all
close all

dir_work = spm_select(1, 'dir', 'Give data dir');
% Controls are cxxtxx, patients are pxxtxx, you can find both in the
% freesurfer directory
dir_data = spm_select([1 1000], 'dir', 'Give all data');
nr_data = size(dir_data, 1);
 
% FreesurferColorLUT contains the Freesurfer segment names and their values 
file_lut = [dir_work '/FreeSurferColorLUT.txt'];
[label_id, label_name, d1,d2, d3, d4] = ...
    textread(file_lut, '%d%s%d%d%d%d');
nr_labels = length(label_id);
   
file_segm_nii = sprintf('%smri/aparc+aseg.nii', ...
    deblank(dir_data(1, :)));
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

[data_struct] = ...
    hd_network_Cth(dir_data);
id_fu = zeros(nr_data, 1);
id_patient = zeros(nr_data, 1);
id_grp = zeros(nr_data, 1);
for d = 1:nr_data,
    
    dir_local = deblank(dir_data(d, :));
    I = strfind(dir_local, 'freesurfer');;
    str_file = dir_local((I+11):(end-1));
    if findstr(str_file, 'bl'),
        
        id_fu(d) = 0;
        
    elseif findstr(str_file, 'fu12'),
        
        id_fu(d) = 1;
        
    elseif findstr(str_file, 'fu24'),
        
        id_fu(d) = 2;
        
    end
    
    I1 = findstr(str_file, 'patient');
    I2 = findstr(str_file((I1+7):end), '_');
    id_patient(d) = str2num(str_file((I1+7):(I1 + 7 + I2 - 2)));
    
    I1 = findstr(str_file, 'grp');
    id_grp(d) = str2num(str_file((I1 + 3):end));
    
end
        
data_struct.I_label = I_label;
data_struct.label_name = label_name;
data_struct.V_segm = V_segm;
data_struct.id_fu = id_fu;
data_struct.id_patient = id_patient;
data_struct.id_grp = id_grp;

save data_HD_CthVol data_struct


save data_hd data_struct