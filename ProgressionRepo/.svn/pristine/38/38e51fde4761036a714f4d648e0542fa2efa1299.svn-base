%% Loading data and computing atrophy
clear all
close all

file_data = spm_select(1, 'mat');
load(file_data);
id_roi_select = data_struct.roi_id;
I_label_select = data_struct.I_roi_select;
name_roi = data_struct.name_roi;
I_lh = 1:38;
I_rh = 39:76;
I_hem{1} = I_lh;
I_hem{2} = I_rh;
id_group = data_struct.id_group;
nr_event = 38;

atrophy_fu12.data_mean = mean(cat(3, log(data_struct.data_jacc_ffd_new.data_median(I_lh, :, 1)), ...
    log(data_struct.data_jacc_ffd_new.data_median(I_rh, :, 1))), 3);
atrophy_fu24.data_mean = mean(cat(3, log(data_struct.data_jacc_ffd_new.data_median(I_lh, :, 2)), ...
    log(data_struct.data_jacc_ffd_new.data_median(I_rh, :, 2))), 3);

atrophy_control_fu12.data_mean = atrophy_fu12.data_mean(:, id_group == 0);
I_purge = find(atrophy_control_fu12.data_mean(1, :) == -Inf);
atrophy_control_fu12.data_mean(:, I_purge) = [];
atrophy_control_fu24.data_mean = atrophy_fu24.data_mean(:, id_group == 0);
I_purge = find(atrophy_control_fu24.data_mean(1, :) == -Inf);
atrophy_control_fu24.data_mean(:, I_purge) = [];
atrophy_patient_fu12.data_mean = atrophy_fu12.data_mean(:, id_group > 0);
I_purge = find(atrophy_patient_fu12.data_mean(1, :) == -Inf);
atrophy_patient_fu12.data_mean(:, I_purge) = [];
atrophy_patient_fu24.data_mean = atrophy_fu24.data_mean(:, id_group > 0);
I_purge = find(atrophy_patient_fu24.data_mean(1, :) == -Inf);
atrophy_patient_fu24.data_mean(:, I_purge) = [];

flag_mixt = 0;
flag_vis = 0;
flag_filt = 0;
X_fu12 = ...
    compute_p_A_D(atrophy_control_fu12, ...
    atrophy_patient_fu12, flag_mixt, flag_filt, flag_vis);
X_fu24 = ...
    compute_p_A_D(atrophy_control_fu24, ...
    atrophy_patient_fu24, flag_mixt, flag_filt, flag_vis);
X = cat(2, X_fu12, X_fu24);
thr_sig = 0.95;
X(X >= thr_sig) = 1;
X(X < thr_sig) = 0;
%% Performing MCMC
nr_it_mcmc = 1e4;
nr_it_burnin = 2e4;
nr_it_hillclimb = 2e3;
thinning = 1e0;
nr_hillclimb = 5;
version_like = 2;

roi_select = [1:nr_event];
X_select = X(roi_select, :);

[parm_struct, diag_struct] = ...
    EventOnsetModel(X_select, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb);

subplot(4, 2, 1), imagesc(X_select)
subplot(4, 2, 2), plot(diag_struct.logp_hillclimb)
subplot(4, 2, 3), imagesc(parm_struct.model_progress_mean)
subplot(4, 2, 4), plot(diag_struct.logp_mcmc)
subplot(4, 2, 5), imagesc(parm_struct.model_progress_max)
subplot(4, 2, 6), imagesc(parm_struct.mu)
subplot(4, 2, 8), imagesc(parm_struct.sigma)
