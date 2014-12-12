%% Fluidreg
clear all
close all

work_dir = '/cs/research/vision/camino1/camino/fonteijn/forDannyAndHubert';
dir_jac = '/cs/research/vision/camino1/camino/fonteijn/forDannyAndHubert/fluidreg/nonRigid/';
dir_control = dir([work_dir '/freesurfer/c*']);
nr_control = length(dir_control);
dir_patient = dir([work_dir '/freesurfer/p*']);
nr_patient = length(dir_patient);
kern_smooth = [8 8 8];

id_c = zeros(nr_control, 2);
for c = 1:nr_control,
    
    id_c(c, 1) = str2num(dir_control(c).name(2:3));
    id_c(c, 2) = str2num(dir_control(c).name(5:6));
    
end
id_p = zeros(nr_patient, 2);
for c = 1:nr_patient,
    
    id_p(c, 1) = str2num(dir_patient(c).name(2:3));
    id_p(c, 2) = str2num(dir_patient(c).name(5:6));
    
end

id_c_unique = unique(id_c(:, 1));
id_p_unique = unique(id_p(:, 1));
file_template = '/cs/research/vision/camino1/camino/fonteijn/spm8/templates/T1.nii';
V_template = spm_vol(file_template);
flags_write.vox = [4 4 4];
cnt = 1;
for c = 1:length(id_c_unique),
    
    I_c = find(id_c(:, 1) == id_c_unique(c));
    t_c = id_c(I_c, 2);
    file_T1 = [work_dir '/freesurfer/' dir_control(I_c(1)).name '/mri/T1.nii'];
    if ~exist(file_T1, 'file'),
        
        file_T1_mgz = [deblank(dir_subj(subj, :)) 'mri/T1.mgz'];
        unix(sprintf('mri_convert %s %s', file_T1_mgz, file_T1));
        
    end
    V_T1 = spm_vol(file_T1);
    spm_normalise(V_template, V_T1, 'norm_parm.mat');
    clear file_rewrite file_norm file_smooth
    for t = 2:length(t_c),
        
        I_t = find(t_c == t);
        file_rewrite(t-1, :) = sprintf('%sc%.2d_jac_%.2d-to-01.nii', ...
            dir_jac, id_c_unique(c), t);
        file_norm(t-1, :) = sprintf('%swc%.2d_jac_%.2d-to-01.nii', ...
            dir_jac, id_c_unique(c), t);
        file_smooth(t-1, :) = sprintf('%sswc%.2d_jac_%.2d-to-01.nii', ...
            dir_jac, id_c_unique(c), t);
        file_control(cnt, :) = file_smooth(t-1, :);
        cnt = cnt + 1;
        
    end
    spm_write_sn(file_rewrite, 'norm_parm.mat', flags_write);
    for c2 = 1:size(file_norm, 1),
        
        spm_smooth(file_norm(c2, :), file_smooth(c2, :), kern_smooth);
        
    end
    fprintf('control: %d\n', c)
    
end
        
cnt = 1;
for p = 1:length(id_p_unique),
    
    I_p = find(id_p(:, 1) == id_p_unique(p));
    t_p = id_p(I_p, 2);
    file_T1 = [work_dir '/freesurfer/' dir_patient(I_p(1)).name '/mri/T1.nii'];
    if ~exist(file_T1, 'file'),
        
        file_T1_mgz = [deblank(dir_subj(subj, :)) 'mri/T1.mgz'];
        unix(sprintf('mri_convert %s %s', file_T1_mgz, file_T1));
        
    end
    V_T1 = spm_vol(file_T1);
    spm_normalise(V_template, V_T1, 'norm_parm.mat');
    clear file_rewrite file_norm file_smooth
    for t = 2:length(t_p),
        
        I_t = find(t_p == t);
        file_rewrite(t-1, :) = sprintf('%sp%.2d_jac_%.2d-to-01.nii', ...
            dir_jac, id_p_unique(p), t);
        file_norm(t-1, :) = sprintf('%swp%.2d_jac_%.2d-to-01.nii', ...
            dir_jac, id_p_unique(p), t);
        file_smooth(t-1, :) = sprintf('%sswp%.2d_jac_%.2d-to-01.nii', ...
            dir_jac, id_p_unique(p), t);
        file_patient(cnt, :) = file_smooth(t-1, :);
        cnt = cnt + 1;
        
    end
    spm_write_sn(file_rewrite, 'norm_parm.mat', flags_write);
    for p2 = 1:size(file_norm, 1),
        
        spm_smooth(file_norm(p2, :), file_smooth(p2, :), kern_smooth);
        
    end
    fprintf('patient: %d\n', p)
    
end

file_aal = '/cs/research/vision/camino1/camino/fonteijn/aal_toolbox/ROI_MNI_V4.nii';
V_aal = spm_vol(file_aal);
img_aal = spm_read_vols(V_aal);
I_aal_hr = find(img_aal ~= 0);
V = spm_vol(file_control(1, :));
I_aal_lr = convert_indimg(I_aal_hr, V_aal, V);
nr_vox = length(I_aal_lr);
data_control = zeros(nr_vox, size(file_control, 1));
data_patient = zeros(nr_vox, size(file_patient, 1));
for c = 1:size(file_control, 1),
    
    V = spm_vol(file_control(c, :));
    img = spm_read_vols(V);
    data_control(:, c) = img(I_aal_lr);
    
end
for p = 1:size(file_patient, 1),
    
    V = spm_vol(file_patient(p, :));
    img = spm_read_vols(V);
    data_patient(:, p) = img(I_aal_lr);
    
end

data_clust = [data_control data_patient];
nr_clust = 20;
idx = kmeans(data_clust, nr_clust, 'Replicates', 10, 'Display', 'iter');
atrophy_control = zeros(nr_clust, size(data_control, 2));
atrophy_patient = zeros(nr_clust, size(data_patient, 2));
for clust = 1:20,
    
    I_clust = find(idx == clust);
    atrophy_control(clust, :) = mean(data_control(I_clust, :));
    atrophy_patient(clust, :) = mean(data_patient(I_clust, :));
    
end

p_A_D = ...
    compute_p_A_D(atrophy_control, ...
    atrophy_patient, 1);

nr_it_mcmc = 1e4;
nr_it_burnin = 1e3;
nr_it_hillclimb = 1e3;
thinning = 1e0;
nr_hillclimb = 5;
[parm_struct, diag_struct] = ...
    AtrophyModelMCMC2(p_A_D, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb);

img_order = zeros(V.dim);
for clust = 1:nr_clust,
    
    I_clust = find(idx == clust);
    img_order(I_aal_lr(I_clust)) = parm_struct.order_events_max(clust);
    
end
V_order= V;
V_order.fname = 'img_order_fluidreg_clust.nii';
spm_create_vol(V_order);
spm_write_vol(V_order, img_order);

%% drcreg
clear all
close all

work_dir = '/cs/research/vision/camino1/camino/fonteijn/forDannyAndHubert';
dir_jac = '/cs/research/vision/camino1/camino/fonteijn/forDannyAndHubert/drcFluid/';
dir_control = dir([work_dir '/freesurfer/c*']);
nr_control = length(dir_control);
dir_patient = dir([work_dir '/freesurfer/p*']);
nr_patient = length(dir_patient);
kern_smooth = [8 8 8];

id_c = zeros(nr_control, 2);
for c = 1:nr_control,
    
    id_c(c, 1) = str2num(dir_control(c).name(2:3));
    id_c(c, 2) = str2num(dir_control(c).name(5:6));
    
end
id_p = zeros(nr_patient, 2);
for c = 1:nr_patient,
    
    id_p(c, 1) = str2num(dir_patient(c).name(2:3));
    id_p(c, 2) = str2num(dir_patient(c).name(5:6));
    
end

id_c_unique = unique(id_c(:, 1));
id_p_unique = unique(id_p(:, 1));
file_template = '/cs/research/vision/camino1/camino/fonteijn/spm8/templates/T1.nii';
V_template = spm_vol(file_template);
flags_write.vox = [4 4 4];
cnt = 1;
for c = 1:length(id_c_unique),
    
    I_c = find(id_c(:, 1) == id_c_unique(c));
    t_c = id_c(I_c, 2);
    file_T1 = [work_dir '/freesurfer/' dir_control(I_c(1)).name '/mri/T1.nii'];
    if ~exist(file_T1, 'file'),
        
        file_T1_mgz = [deblank(dir_subj(subj, :)) 'mri/T1.mgz'];
        unix(sprintf('mri_convert %s %s', file_T1_mgz, file_T1));
        
    end
    V_T1 = spm_vol(file_T1);
    spm_normalise(V_template, V_T1, 'norm_parm.mat');
    clear file_rewrite file_norm file_smooth
    for t = 2:length(t_c),
        
        I_t = find(t_c == t);
        file_rewrite(t-1, :) = sprintf('%sc%.2dt01/c%.2dt01-c%.2dt%.2d.img', ...
            dir_jac, id_c_unique(c), id_c_unique(c), id_c_unique(c),t);
        file_norm(t-1, :) = sprintf('%sc%.2dt01/wc%.2dt01-c%.2dt%.2d.img', ...
            dir_jac, id_c_unique(c), id_c_unique(c), id_c_unique(c),t);
        file_smooth(t-1, :) = sprintf('%sc%.2dt01/swc%.2dt01-c%.2dt%.2d.img', ...
            dir_jac, id_c_unique(c), id_c_unique(c), id_c_unique(c),t);
        file_control(cnt, :) = file_smooth(t-1, :);
        cnt = cnt + 1;
        
    end
    spm_write_sn(file_rewrite, 'norm_parm.mat', flags_write);
    for c2 = 1:size(file_norm, 1),
        
        spm_smooth(file_norm(c2, :), file_smooth(c2, :), kern_smooth);
        
    end
    fprintf('control: %d\n', c)
    
end
        
cnt = 1;
for p = 1:length(id_p_unique),
    
    I_p = find(id_p(:, 1) == id_p_unique(p));
    t_p = id_p(I_p, 2);
    file_T1 = [work_dir '/freesurfer/' dir_patient(I_p(1)).name '/mri/T1.nii'];
    if ~exist(file_T1, 'file'),
        
        file_T1_mgz = [deblank(dir_subj(subj, :)) 'mri/T1.mgz'];
        unix(sprintf('mri_convert %s %s', file_T1_mgz, file_T1));
        
    end
    V_T1 = spm_vol(file_T1);
    spm_normalise(V_template, V_T1, 'norm_parm.mat');
    clear file_rewrite file_norm file_smooth
    for t = 2:length(t_c),
        
        I_t = find(t_c == t);
        file_rewrite(t-1, :) = sprintf('%sp%.2dt01/p%.2dt01-p%.2dt%.2d.img', ...
            dir_jac, id_c_unique(c), id_c_unique(c), id_c_unique(c),t);
        file_norm(t-1, :) = sprintf('%sp%.2dt01/wp%.2dt01-p%.2dt%.2d.img', ...
            dir_jac, id_c_unique(c), id_c_unique(c), id_c_unique(c),t);
        file_smooth(t-1, :) = sprintf('%sp%.2dt01/swp%.2dt01-p%.2dt%.2d.img', ...
            dir_jac, id_c_unique(c), id_c_unique(c), id_c_unique(c),t);
        file_control(cnt, :) = file_smooth(t-1, :);
        cnt = cnt + 1;
        
    end
 
    spm_write_sn(file_rewrite, 'norm_parm.mat', flags_write);
    for p2 = 1:size(file_norm, 1),
        
        spm_smooth(file_norm(p2, :), file_smooth(p2, :), kern_smooth);
        
    end
    fprintf('patient: %d\n', p)
    
end

file_aal = '/cs/research/vision/camino1/camino/fonteijn/aal_toolbox/ROI_MNI_V4.nii';
V_aal = spm_vol(file_aal);
img_aal = spm_read_vols(V_aal);
I_aal_hr = find(img_aal ~= 0);
V = spm_vol(file_control(1, :));
I_aal_lr = convert_indimg(I_aal_hr, V_aal, V);
nr_vox = length(I_aal_lr);
data_control = zeros(nr_vox, size(file_control, 1));
data_patient = zeros(nr_vox, size(file_patient, 1));
for c = 1:size(file_control, 1),
    
    V = spm_vol(file_control(c, :));
    img = spm_read_vols(V);
    data_control(:, c) = img(I_aal_lr);
    
end
for p = 1:size(file_patient, 1),
    
    V = spm_vol(file_patient(p, :));
    img = spm_read_vols(V);
    data_patient(:, p) = img(I_aal_lr);
    
end

data_clust = [data_control data_patient];
nr_clust = 20;
idx = kmeans(data_clust, nr_clust, 'Replicates', 10, 'Display', 'iter');
atrophy_control = zeros(nr_clust, size(data_control, 2));
atrophy_patient = zeros(nr_clust, size(data_patient, 2));
for clust = 1:20,
    
    I_clust = find(idx == clust);
    atrophy_control(clust, :) = mean(data_control(I_clust, :));
    atrophy_patient(clust, :) = mean(data_patient(I_clust, :));
    
end

p_A_D = ...
    compute_p_A_D(atrophy_control, ...
    atrophy_patient, 1);

nr_it_mcmc = 1e4;
nr_it_burnin = 1e3;
nr_it_hillclimb = 1e3;
thinning = 1e0;
nr_hillclimb = 5;
[parm_struct, diag_struct] = ...
    AtrophyModelMCMC2(p_A_D, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb);

img_order = zeros(V.dim);
for clust = 1:nr_clust,
    
    I_clust = find(idx == clust);
    img_order(I_aal_lr(I_clust)) = parm_struct.order_events_max(clust);
    
end
V_order= V;
V_order.fname = 'img_order_drcreg_clust.nii';
spm_create_vol(V_order);
spm_write_vol(V_order, img_order);

%% Segmentation-based

load job_segment
job_segment = matlabbatch;
dir_subj = spm_select([1 200], 'dir');
nr_subj = size(dir_subj, 1);
for subj = 1:nr_subj,
    
    file_T1 = [deblank(dir_subj(subj, :)) 'mri/T1.nii'];
    if ~exist(file_T1, 'file'),
        
        file_T1_mgz = [deblank(dir_subj(subj, :)) 'mri/T1.mgz'];
        unix(sprintf('mri_convert %s %s', file_T1_mgz, file_T1));
        
    end
    job_segment{1}.spm.spatial.preproc.data{1} = ...
        file_T1;
    matlabbatch = job_segment;
    save job_segment matlabbatch
    spm_jobman('run', 'job_segment.mat')
    fprintf('subject: %d in %d subjects\n', ...
        subj, nr_subj)
    
end

kern_smooth = [12 12 12];
for subj = 1:nr_subj,
    
    file_segm = [deblank(dir_subj(subj, :)) 'mri/mwc1T1.nii'];
    file_smooth = [deblank(dir_subj(subj, :)) 'mri/smwc1T1.nii'];
    spm_smooth(file_segm, file_smooth, kern_smooth)
    fprintf('subject: %d in %d subjects\n', ...
        subj, nr_subj)
    
end
    
    
