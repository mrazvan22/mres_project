%% Setting up simulated data
clear all
close all

cnt_fig = 1;
version_mcmc = 2;
flag_mixt = 3;
if flag_mixt == 1,
    
    flag_filt = input('Filtered controls? \n');
    
else
    
    flag_filt = 0;
    
end

file_data = spm_select(1, 'mat');
load(file_data);
for flag_reg = 1:2,
    
    atrophy_patient{flag_reg}.data_mean = data_struct.data_patient{flag_reg}.data_median;
    atrophy_control{flag_reg}.data_mean = data_struct.data_control{flag_reg}.data_median;
    
end
id_roi_select = data_struct.id_roi_select;
I_label_select = data_struct.I_label_select;
name_roi = data_struct.name_roi;
[nr_roi nr_pat] = size(atrophy_patient{1}.data_mean);
nr_roi_lhrh = nr_roi/2;

nr_events = nr_roi;
p_A_D = zeros(nr_events, nr_pat, 2);  
prob_pat = zeros(nr_events, nr_pat, 2);  
prob_con = zeros(nr_events, nr_pat, 2);  
span = zeros(nr_events, 2, 2);
for flag_reg = 1:2,
    
    [p_A_D(1:nr_roi, :, flag_reg), prob_pat(1:nr_roi, :, flag_reg), ...
        prob_con(1:nr_roi, :, flag_reg), span(1:nr_roi, :, flag_reg)] = ...
        compute_p_A_D(atrophy_control{flag_reg}, ...
        atrophy_patient{flag_reg}, flag_mixt, flag_filt);
    
end
for pat = 1:nr_pat,
    
    p_A_D((nr_roi+1):(nr_roi + data_struct.cat_patient(pat)), pat, :) = 1;
    prob_pat((nr_roi+1):(nr_roi + data_struct.cat_patient(pat)), pat, :) = 1;
    
end
prob_con((nr_roi+1):(nr_roi+3), :, :) = 1-prob_pat((nr_roi+1):(nr_roi+3), :, :);
span((nr_roi_lhrh+1):(nr_roi_lhrh+3), 1, :) = 0;
span((nr_roi+1):(nr_roi+3), 2, :) = 1;
idx_clinevent = zeros(nr_roi + 3, 1);
idx_clinevent((nr_roi+1):(nr_roi+3)) = 1;


%% Performing mcmc
nr_it_mcmc = 1e5;
nr_it_burnin = 1e4;
nr_it_hillclimb = 1e4;
thinning = 1e0;
nr_hillclimb = 5;
version_like = 2;

flag_reg = 1;

roi_select = [1 2 9];
p_A_D_select = p_A_D(roi_select, :, flag_reg);
prob_pat_select = prob_pat(roi_select, :, flag_reg);
prob_con_select = prob_con(roi_select, :, flag_reg);
prob_pat_select(prob_pat_select == 0) = 1e-3;
prob_con_select(prob_con_select == 0) = 1e-3;
span_select = span(roi_select, :, flag_reg);
idx_clinevent_select = idx_clinevent(roi_select);

[parm_struct, diag_struct] = ...
    AtrophyModelMCMC2Corrected(prob_pat_select, prob_con_select, span_select, ...
    idx_clinevent_select, nr_it_hillclimb, nr_it_burnin, nr_it_mcmc, ...
    thinning, nr_hillclimb, version_like);

subplot(3, 2, 1), imagesc(p_A_D_select)
subplot(3, 2, 3), imagesc(parm_struct.order_events)
subplot(3, 2, 2), plot(diag_struct.logp_hillclimb)
subplot(3, 2, 4), plot(diag_struct.logp_order_events)
subplot(3, 2, 6), plot(parm_struct.f_uniform')


nr_it_burnin = 1e5;
nr_it_hillclimb = 1e4;
nr_hillclimb = 5;
nr_it_mcmc = 1e3;

thr_vec = [1e-4 1e-2 0.2 0.4 0.6];
for thr1 = 1:length(thr_vec),
    
    for thr2 = 1:length(thr_vec),
        
        p_A_D_thr = p_A_D(:, :, 1);
        p_A_D_thr(p_A_D_thr < thr_vec(thr1)) = thr_vec(thr1);
        p_A_D_thr(p_A_D_thr > 1-thr_vec(thr2)) = 1-thr_vec(thr2);
        [parm_struct{thr1, thr2}, diag_struct{thr1, thr2}] = ...
            AtrophyModelMCMC2a(p_A_D_thr, nr_it_hillclimb, ...
            nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_mcmc);
        fprintf('thr1: %d thr2: %d\n', thr1, thr2)
        
    end
    
end

for thr1 = 1:length(thr_vec),
    
    for thr2 = 1:length(thr_vec),
        
        logp(thr1, thr2) = diag_struct{thr1, thr2}.logp_order_events_max;
        
    end
    
end






