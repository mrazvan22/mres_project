function [data_control, data_patient, I_label] = ...
    ad_network(dir_control, dir_patient, dir_work, ...
    id_roi_select, flag_reg)

if flag_reg == 1,
    
    fprintf('Fluid registration\n')
    
elseif flag_reg == 2,
    
    fprintf('drc fluids\n')
    
else
    
    fprintf('No recognized registration\n')
    
end
data_control = ...
    extract_data(dir_control, 'c', id_roi_select, ...
    dir_work, flag_reg);
data_patient = ...
    extract_data(dir_patient, 'p', id_roi_select, ...
    dir_work, flag_reg);

file_segm_nii = [dir_control(1, :) 'mri/aparc+aseg.nii'];
V_segm = spm_vol(file_segm_nii);
img_segm = spm_read_vols(V_segm);
file_hippocampus_lh_nii = sprintf('%s/data/c01t01_leftHippo.nii', ...
    dir_work);
file_hippocampus_rh_nii = sprintf('%s/data/c01t01_rightHippo.nii', ...
    dir_work);
V_hippocampus_lh = spm_vol(file_hippocampus_lh_nii);
V_hippocampus_rh = spm_vol(file_hippocampus_rh_nii);
img_hippocampus_lh = spm_read_vols(V_hippocampus_lh);
img_hippocampus_rh = spm_read_vols(V_hippocampus_rh);
I = find(img_hippocampus_lh > 0);
I_hippocampus{1} = convert_indimg(I, V_hippocampus_lh, V_segm);
I = find(img_hippocampus_rh > 0);
I_hippocampus{2} = convert_indimg(I, V_hippocampus_rh, V_segm);

for roi = 1:length(id_roi_select),

    I_label{roi} = find(img_segm == id_roi_select(roi));

end
I_label{id_roi_select == 17} = I_hippocampus{2};
I_label{id_roi_select == 53} = I_hippocampus{1};

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function data = extract_data(dir_subj, str_corp, ...
    id_roi_select, dir_work, flag_reg)

opts_gmm = zeros(1, 14);
opts_gmm(1) = -1;
opts_gmm(3) = 1e-4;
opts_gmm(14) = 1e3;

nr_roi_select = length(id_roi_select);
for tp = 1:size(dir_subj, 1),
    
    subj_id_str(tp, :) = dir_subj(tp, (end-6):(end-1));
    subj_id(tp) = str2num(subj_id_str(tp, 2:3));
    subj_tp(tp) = str2num(subj_id_str(tp, 5:6));
    
end

subj_id_unique = unique(subj_id);
data_mean = zeros(nr_roi_select, size(dir_subj, 1) - length(subj_id_unique));
data_median = zeros(nr_roi_select, size(dir_subj, 1) - length(subj_id_unique));
data_std = zeros(nr_roi_select, size(dir_subj, 1) - length(subj_id_unique));
data_nr_vox = zeros(nr_roi_select, size(dir_subj, 1) - length(subj_id_unique));
data_meanneg = zeros(nr_roi_select, size(dir_subj, 1) - length(subj_id_unique));
data_medianneg = zeros(nr_roi_select, size(dir_subj, 1) - length(subj_id_unique));
data_medianneg_nw = zeros(nr_roi_select, size(dir_subj, 1) - length(subj_id_unique));
cnt_jacc = 1;
for pid = subj_id_unique,
    
    I_subj = find(subj_id == pid);
    file_segm_nii = [dir_subj(I_subj(1), :) 'mri/aparc+aseg.nii'];
    if ~exist(file_segm_nii, 'file'),
        
        file_segm_mgz = [dir_subj(I_subj(1), :) 'mri/aparc+aseg.mgz'];
        unix(sprintf('mri_convert %s %s', file_segm_mgz, file_segm_nii));
        
    end    
    V_segm = spm_vol(file_segm_nii);
    img_segm = spm_read_vols(V_segm);
    
    file_hippocampus_lh_nii = sprintf('%s/data/%s%.2dt%.2d_leftHippo.nii', ...
        dir_work, str_corp, pid, 1);
    file_hippocampus_rh_nii = sprintf('%s/data/%s%.2dt%.2d_rightHippo.nii', ...
        dir_work, str_corp, pid, 1);
    if ~exist(file_hippocampus_lh_nii, 'file'),
        
        file_hippocampus_lh_mgz = sprintf('%s/data/%s%.2dt%.2d_leftHippo.mgz', ...
            dir_work, str_corp, pid, 1);
        unix(sprintf('mri_convert %s %s', file_hippocampus_lh_mgz, ...
            file_hippocampus_lh_nii));

    end
    if ~exist(file_hippocampus_rh_nii, 'file'),
        
        file_hippocampus_rh_mgz = sprintf('%s/data/%s%.2dt%.2d_rightHippo.mgz', ...
            dir_work, str_corp, pid, 1);
        unix(sprintf('mri_convert %s %s', file_hippocampus_rh_mgz, ...
            file_hippocampus_rh_nii));
        
    end
    V_hippocampus_lh = spm_vol(file_hippocampus_lh_nii);
    V_hippocampus_rh = spm_vol(file_hippocampus_rh_nii);    
    img_hippocampus_lh = spm_read_vols(V_hippocampus_lh);
    img_hippocampus_rh = spm_read_vols(V_hippocampus_rh);
    I = find(img_hippocampus_lh > 0);
    I_hippocampus{1} = convert_indimg(I, V_hippocampus_lh, V_segm);
    I = find(img_hippocampus_rh > 0);
    I_hippocampus{2} = convert_indimg(I, V_hippocampus_rh, V_segm);
    
    for roi = 1:nr_roi_select,,
        
        I_roi{roi} = find(img_segm == id_roi_select(roi));
        
    end
    subj_tp_local = subj_tp(I_subj);
    for t = 2:max(subj_tp_local),
        
        if flag_reg == 1,
            
            file_jacc_t = sprintf('%s/fluidreg/nonRigid/%s%.2d_jac_%.2d-to-01.nii', ...
                dir_work, str_corp, pid, t);
            
        elseif flag_reg == 2,
            
            file_jacc_t = sprintf('%s/drcFluid/%s%.2dt01/%s%.2dt01-%s%.2dt%.2d_unsep.nii', ...
                dir_work, str_corp, pid, str_corp, pid, str_corp, pid, t);
            
        end

        img_jacc = spm_read_vols(spm_vol(file_jacc_t));
        if flag_reg == 2,
            
            img_jacc = exp(img_jacc);
            
        end
       
        for roi = 1:nr_roi_select
            
            data_local = img_jacc(I_roi{roi});
            data_mean(roi, cnt_jacc) = mean(data_local);
            data_median(roi, cnt_jacc) = median(data_local);
            data_std(roi, cnt_jacc) = std(data_local);
            data_nr_vox(roi, cnt_jacc) = length(data_local);
            I_neg = find(data_local < 1);
            data_meanneg(roi, cnt_jacc) = mean(data_local(I_neg))*...
                (length(I_neg)/length(data_local));
            data_medianneg(roi, cnt_jacc) = median(data_local(I_neg))*...
                (length(I_neg)/length(data_local));
            data_medianneg_nw(roi, cnt_jacc) = median(data_local(I_neg));
            
            
        end
        
        for i_hc = 1:2,
            
            data_local = img_jacc(I_hippocampus{i_hc});
            data_hc_mean(i_hc, cnt_jacc) = ...
                mean(data_local);
            data_hc_median(i_hc, cnt_jacc) = ...
                median(data_local);
            data_hc_std(i_hc, cnt_jacc) = ...
                std(data_local);
            data_hc_nr_vox(i_hc, cnt_jacc) = length(data_local);
            I_neg = find(data_local < 1);
            data_hc_meanneg(i_hc, cnt_jacc) = mean(data_local(I_neg))*...
                (length(I_neg)/length(data_local));
            data_hc_medianneg(i_hc, cnt_jacc) = median(data_local(I_neg))*...
                (length(I_neg)/length(data_local));
            data_hc_medianneg_nw(i_hc, cnt_jacc) = median(data_local(I_neg));

            
        end
        fprintf('Data set: %d\n', ...
            cnt_jacc)
        cnt_jacc = cnt_jacc + 1;
        
    end
    
end
    
data_mean(id_roi_select == 17, :) = data_hc_mean(2, :);
data_mean(id_roi_select == 53, :) = data_hc_mean(1, :);
data_meanneg(id_roi_select == 17, :) = data_hc_meanneg(2, :);
data_meanneg(id_roi_select == 53, :) = data_hc_meanneg(1, :);
data_median(id_roi_select == 17, :) = data_hc_median(2, :);
data_median(id_roi_select == 53, :) = data_hc_median(1, :);
data_medianneg(id_roi_select == 17, :) = data_hc_medianneg(2, :);
data_medianneg(id_roi_select == 53, :) = data_hc_medianneg(1, :);
data_medianneg_nw(id_roi_select == 17, :) = data_hc_medianneg_nw(2, :);
data_medianneg_nw(id_roi_select == 53, :) = data_hc_medianneg_nw(1, :);
data_std(id_roi_select == 17, :) = data_hc_std(2, :);
data_std(id_roi_select == 53, :) = data_hc_std(1, :);
data_nr_vox(id_roi_select == 17, :) = data_hc_nr_vox(2, :);
data_nr_vox(id_roi_select == 53, :) = data_hc_nr_vox(1, :);

data.data_mean = data_mean;
data.data_meanneg = data_meanneg;
data.data_median = data_median;
data.data_medianneg = data_medianneg;
data.data_medianneg_nw = data_medianneg_nw;
data.data_std = data_std;
data.data_nr_vox = data_nr_vox;
