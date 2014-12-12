function [data_control, data_patient, I_label] = ...
    ad_network2(dir_control, dir_patient, dir_work, ...
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
data_nomixt = zeros(nr_roi_select, size(dir_subj, 1) - length(subj_id_unique));
data_mixt = zeros(nr_roi_select, size(dir_subj, 1) - length(subj_id_unique));
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
        if flag_reg == 1,
            
            img_jacc = log(img_jacc);
            
        end
        %         subplot(121), hist(img_jacc(I_hippocampus{1}))
        %         fprintf('mean: %d\n', mean(img_jacc(I_hippocampus{1})))
        %         fprintf('median: %d\n', median(img_jacc(I_hippocampus{1})))
        %         data_local = img_jacc(I_hippocampus{1});
        %         fprintf('measure hubfon %d\n', mean(data_local(data_local < 0))*length(find(data_local < 0))/length(data_local));
        %
        %         subplot(122), hist(img_jacc(I_roi{8}))
        %         fprintf('mean: %d\n', mean(img_jacc(I_roi{8})))
        %         fprintf('median: %d\n', median(img_jacc(I_roi{8})))
        %         data_local = img_jacc(I_roi{8});
        %         fprintf('measure hubfon %d\n', mean(data_local(data_local < 0))*length(find(data_local < 0))/length(data_local));
       
        for roi = 1:nr_roi_select
            
            %             if sum(img_jacc(I_roi{roi})),
            %
            %                 data_local = img_jacc(I_roi{roi});
            %                 BIC_local = zeros(3, 10);
            %                 for c = 1:3,
            %
            %                     for it = 1:10,
            %
            %                         gmix{c, it} = gmm(1, c, 'full');
            %                         gmix{c, it} = ...
            %                             gmmem_fsem(gmix{c, it}, data_local, opts_gmm);
            %                         BIC_local(c, it) = gmmem_BIC(gmix{c, it}, data_local);
            %
            %                     end
            %
            %                 end
            %                 I_opt = find(BIC_local == min(BIC_local(:)));
            %                 [c_opt, it_opt] = ind2sub(size(BIC_local), I_opt);
            %                 gmix_opt = gmix{c_opt(1), it_opt(1)};
            %                 I_min = find(gmix_opt.centres == min(gmix_opt.centres));
            %                 if gmix_opt.priors(I_min(1)) > 0.2,
            %
            %                     data_mixt(roi, cnt_jacc) = ...
            %                         gmix_opt.centres(I_min(1));
            %
            %                 else
            %
            %                     data_mixt(roi, cnt_jacc) = ...
            %                         mean(img_jacc(I_label{label}));
            %
            %                 end
            %
            %
            %             else
            %
            %                 data_mixt(roi, cnt_jacc) = ...
            %                     mean(img_jacc(I_roi{roi}));
            %
            %             end
            data_local = img_jacc(I_roi{roi});
            data_nomixt(roi, cnt_jacc) = ...
                mean(data_local(data_local < 0))*(length(find(data_local < 0))/length(data_local));
            
        end
        
        for i_hc = 1:2,
            
            data_local = img_jacc(I_hippocampus{i_hc});
            %             BIC_local = zeros(3, 10);
            %             for c = 1:3,
            %
            %                 for it = 1:10,
            %
            %                     gmix{c, it} = gmm(1, c, 'full');
            %                     gmix{c, it} = ...
            %                         gmmem_fsem(gmix{c, it}, data_local, opts_gmm);
            %                     BIC_local(c, it) = gmmem_BIC(gmix{c, it}, data_local);
            %
            %                 end
            %
            %             end
            %             I_opt = find(BIC_local == min(BIC_local(:)));
            %             [c_opt, it_opt] = ind2sub(size(BIC_local), I_opt);
            %             gmix_opt = gmix{c_opt(1), it_opt(1)};
            %             I_min = find(gmix_opt.centres == min(gmix_opt.centres));
            %             if gmix_opt.priors(I_min(1)) > 0.2,
            %
            %                 data_hc_mixt(i_hc, cnt_jacc) = ...
            %                     gmix_opt.centres(I_min(1));
            %
            %             else
            %
            %                 data_hc_mixt(i_hc, cnt_jacc) = ...
            %                     mean(data_local);
            %
            %             end
            data_hc_nomixt(i_hc, cnt_jacc) = ...
                mean(data_local(data_local < 0))*(length(find(data_local < 0))/length(data_local));
            
        end
        fprintf('Data set: %d\n', ...
            cnt_jacc)
        cnt_jacc = cnt_jacc + 1;
        
    end
    
end
    
% data_mixt(id_roi_select == 17, :) = data_hc_mixt(2, :);
% data_mixt(id_roi_select == 53, :) = data_hc_mixt(1, :);
data_nomixt(id_roi_select == 17, :) = data_hc_nomixt(2, :);
data_nomixt(id_roi_select == 53, :) = data_hc_nomixt(2, :);

% data.data_mixt = data_mixt;
data.data_nomixt = data_nomixt;

