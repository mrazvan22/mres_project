% Script that copies the bare essentials from the processed ADNI data
% It only couples the cortical/subcortical segmentation of the Freesurfer
% analysis and the jacobian images of the ffd registration
%% Copying all data
clear all
close all

% starting from a file and a directory
file_data_location = 'subjects_for_interim_testing.txt';
dir_data = '/cluster/project1/drc-freesurfer/matt/forDannyAndHubert/adni-freesurfer/';
dir_work = pwd;

% Now reading in the directory names of the baseline images, the follow-up
% images and the patient status
[dir_baseline, dir_fu, status] = textread(file_data_location, '%s%s%s');
nr_subj = size(dir_baseline, 1);
str_tar = 'tar -cf';

for subj = 1:nr_subj,
    
    dir_local_bl = ...
        sprintf('%s/%s', dir_work, deblank(dir_baseline{subj}));
    unix(sprintf('mkdir %s', dir_local_bl));
    dir_local_fu = sprintf('%s/%s', dir_work, deblank(dir_fu{subj}));
    unix(sprintf('mkdir %s', dir_local_fu));
    
    status_local = status{subj};
    eval(sprintf('save %s/status status_local', dir_local))
    
    % Copying the segmented image from freesurfer...
    unix(sprintf('scp -pr fonteijn@comic100:/%s%s/mri/aparc+aseg.mgz %s', ...
        dir_data, dir_baseline{subj}, dir_local_bl))
    unix(sprintf('scp -pr fonteijn@comic100:/%s%s/stats %s', ...
        dir_data, dir_baseline{subj}, dir_local_bl));
    unix(sprintf('scp -pr fonteijn@comic100:/%s%s/stats %s', ...
        dir_data, dir_fu{subj}, dir_local_fu));

    unix(sprintf('scp -pr fonteijn@comic100:/%snifty_reg/jacobian/%s_r_to_b_nrr_jac.nii.gz %s', ...
        dir_data, dir_baseline{subj}(1:10), dir_local_bl));
    
end
   
%% Reading in all data
work_dir = pwd;
file_lut = [work_dir '/FreeSurferColorLUT.txt'];
[label_id, label_name, d1,d2, d3, d4] = ...
    textread(file_lut, '%d%s%d%d%d%d');

id_roi_select = [10:13 17:18 49:54 1001:1034 2001:2034];
nr_roi_select = length(id_roi_select);
for roi = 1:nr_roi_select,
    
    id_local = find(label_id == id_roi_select(roi));
    name_roi_select{roi} = label_name{id_local};
    
end

data_mean = zeros(nr_roi_select, nr_subj);
data_median = zeros(nr_roi_select, nr_subj);
for subj = 1:nr_subj,
    
    dir_local = sprintf('%s/%s', dir_work, deblank(dir_baseline{subj}));
    file_segm_nii = [dir_local '/aparc+aseg.nii'];
    if ~exist(file_segm_nii, 'file'),
        
        file_segm_mgz = [dir_local '/aparc+aseg.mgz'];
        unix(sprintf('mri_convert %s %s', file_segm_mgz, file_segm_nii));
        
    end    
    V_segm = spm_vol(file_segm_nii);
    img_segm = spm_read_vols(V_segm);
        
    for roi = 1:nr_roi_select,
        
        I_roi{roi} = find(img_segm == id_roi_select(roi));
        
    end
    file_jacc_t = sprintf('%s/nifti_reg_new/%s_r_to_b_nrr_jac.nii', ...
        dir_work, dir_baseline{subj});
    if ~exist(file_jacc_t),
        
        unix(sprintf('gunzip %s/nifti_reg_new/%s_r_to_b_nrr_jac.nii.gz', ...
            dir_work, dir_baseline{subj}));
        
    end
    
    if exist(file_jacc_t, 'file'),
        
        img_jacc = spm_read_vols(spm_vol(file_jacc_t));
        
        for roi = 1:nr_roi_select
            
            data_local = img_jacc(I_roi{roi});
            data_local(find(isnan(data_local))) = [];
            data_mean(roi, subj) = mean(data_local);
            data_median(roi, subj) = median(data_local);
            
        end
        unix(sprintf('rm %s', file_segm_nii));
        fprintf('subj: %d in %d subjects\n', subj, nr_subj);
    end
        
           
end

clear dir_subj_bl dir_subj_fu
for subj = 1:nr_subj,
    
    dir_subj_bl(subj, :) = sprintf('%s/%s/', dir_work, deblank(dir_baseline{subj}));
    dir_subj_fu(subj, :) = sprintf('%s/%s/', dir_work, deblank(dir_fu{subj}));
    
end
data_Cthvol_bl = hd_network_Cth(dir_subj_bl);
data_Cthvol_fu = hd_network_Cth(dir_subj_fu);

status_mat = zeros(nr_subj, 1);
for subj = 1:nr_subj,
    
    if status{subj}(1) == 'N',
        
        status_mat(subj) = 0;
        
    elseif status{subj}(1) == 'M',
        
        status_mat(subj) = 1;
        
    elseif status{subj}(1) == 'A',
        
        status_mat(subj) = 2;
        
    end
    
end

data_controls.data_mean = data_mean(:, status_mat == 0);
data_controls.data_median = data_median(:, status_mat == 0);
data_patients.data_mean = data_mean(:, status_mat > 0);
data_patients.data_median = data_median(:, status_mat > 0);
data_patients.status = status_mat(status_mat > 0);

data_struct.data_ffd_mean = data_mean;
data_struct.data_ffd_median = data_median;
data_struct.data_Cthvol_bl = data_Cthvol_bl;
data_struct.data_Cthvol_fu = data_Cthvol_fu;
data_struct.nr_roi_select = nr_roi_select;
data_struct.name_roi_select = name_roi_select;
data_struct.V_segm = V_segm;
data_struct.I_roi = I_roi;
data_struct.status_mat = status_mat;
save data_adni.mat data_struct