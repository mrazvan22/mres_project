clear all
close all

file_data = spm_select(1, 'mat');
load(file_data);
flag_reg = input('Which registration method (Modat = 1/Freeborough = 2)? ');
flag_mixt = input('Mixtures (1)/no mixtures (0)? ');
if flag_mixt == 1,
    
    atrophy_patient = data_struct.data_patient{flag_reg}.data_median;
    atrophy_control = data_struct.data_control{flag_reg}.data_median;
    
elseif flag_mixt == 0,
    
    atrophy_patient = data_struct.data_patient{flag_reg}.data_median;
    atrophy_control = data_struct.data_control{flag_reg}.data_median;
    
end
    
I_label_select = data_struct.I_label_select;
name_roi = data_struct.name_roi;

I_lh = [1 3:36];
I_rh = [2 37:70];
atrophy_patient_lhrh = cat(4, atrophy_patient(I_lh, :, :), ...
    atrophy_patient(I_rh, :, :));
atrophy_patient_lhrh = mean(atrophy_patient_lhrh, 4);
atrophy_control_lhrh = cat(4, atrophy_control(I_lh, :, :), ...
    atrophy_control(I_rh, :, :));
atrophy_control_lhrh = mean(atrophy_control_lhrh, 4);
atrophy_patient = atrophy_patient_lhrh;
atrophy_control = atrophy_control_lhrh;
[nr_roi nr_pat] = size(atrophy_patient);

p_atrophy = zeros(nr_roi, nr_pat);
for pat = 1:nr_pat,
    
    for roi = 1:nr_roi,
        
        [h, p] = ...
            ttest(atrophy_control(roi, :), ...
            atrophy_patient(roi, pat), [], 'left');
        p_atrophy(roi, pat) = p;
        
    end
    
end
for pat = 1:nr_pat,
    
    p_atrophy((nr_roi+1):(nr_roi + data_struct.cat_patient(pat)), pat, :) = 1;

end
nr_roi = nr_roi + 3;

%% Performing MCMC
nr_it_mcmc = 1e5;
nr_it_burnin = 1e5;
thinning = 100;

order_events_current = randperm(nr_roi);
atrophy_model_current = zeros(nr_roi, nr_roi);
for roi = 1:nr_roi,
    
    atrophy_model_current(order_events_current==roi, roi:end) = 1;
    
end

p_atrophy_model_current = zeros(nr_roi, nr_pat);
for pat = 1:nr_pat,
    
    p_atrophy_local = p_atrophy(:, pat);
    p_atrophy_model_local = zeros(nr_roi, nr_roi);
    for roi = 1:nr_roi,
        
        I_atrophy = find(atrophy_model_current(:, roi) == 1);
        I_noatrophy = find(atrophy_model_current(:, roi) == 0);
        p_atrophy_model_local(I_atrophy, roi) = ...
            p_atrophy_local(I_atrophy);
        p_atrophy_model_local(I_noatrophy, roi) = ...
            1 - p_atrophy_local(I_noatrophy);
        
    end
    p_atrophy_model_local(p_atrophy_model_local == 0) = eps;
    I_maxp = find(sum(log(p_atrophy_model_local)) == max(sum(log(p_atrophy_model_local))));
    p_atrophy_model_current(:, pat) = p_atrophy_model_local(:, I_maxp(1));
    
end
logp_current = sum(log(p_atrophy_model_current(:)));

order_events = zeros(nr_roi, nr_it_mcmc/thinning);
logp_order_events = zeros(nr_it_mcmc/thinning, 1);
cnt_it_mcmc = 1:500:(nr_it_mcmc + nr_it_burnin);
logp_order_events_burnin = zeros(nr_it_burnin, 1);
cnt = 1;
logp_max = logp_current;
order_events_max = order_events_current;
for it_mcmc = 1:(nr_it_mcmc + nr_it_burnin),
    
    % Swapping regions around
    roi_swap = randperm(nr_roi);
    roi_swap = roi_swap(1:2);
    order_events_new = order_events_current;
    order_events_new(roi_swap(1)) = order_events_current(roi_swap(2));
    order_events_new(roi_swap(2)) = order_events_current(roi_swap(1));
    
    atrophy_model_new = zeros(nr_roi, nr_roi);
    for roi = 1:nr_roi,
        
        atrophy_model_new(order_events_new==roi, roi:end) = 1;
        
    end
    
    p_atrophy_model_new = zeros(nr_roi, nr_pat);
    for pat = 1:nr_pat,

        p_atrophy_local = p_atrophy(:, pat);
        p_atrophy_model_local = zeros(nr_roi, nr_roi);
        for roi = 1:nr_roi,

            I_atrophy = find(atrophy_model_new(:, roi) == 1);
            I_noatrophy = find(atrophy_model_new(:, roi) == 0);
            p_atrophy_model_local(I_atrophy, roi) = ...
                p_atrophy_local(I_atrophy);
            p_atrophy_model_local(I_noatrophy, roi) = ...
                1 - p_atrophy_local(I_noatrophy);

        end
        p_atrophy_model_local(p_atrophy_model_local == 0) = eps;
        I_maxp = find(sum(log(p_atrophy_model_local)) == max(sum(log(p_atrophy_model_local))));
        p_atrophy_model_new(:, pat) = p_atrophy_model_local(:, I_maxp(1));

    end
    logp_new = sum(log(p_atrophy_model_new(:)));

    % Accept new ordering?
    if it_mcmc <= nr_it_burnin,
        
        alpha = (logp_new > logp_current);
        
    else
        
        alpha = exp(logp_new - logp_current);
        
    end
    %     alpha = exp(logp_new - logp_current);
    if alpha > rand,
        
        order_events_current = order_events_new;
        logp_current = logp_new;
        if logp_current > logp_max,
            
            logp_max = logp_current;
            order_events_max = order_events_current;
            
        end
        
    end
    if (it_mcmc > nr_it_burnin) & ~mod(it_mcmc, thinning),
        
        order_events(:, cnt) = order_events_current;
        logp_order_events(cnt) = logp_current;
        fprintf('alpha: %d\n', alpha);
        cnt = cnt + 1;
        
    end
    if it_mcmc <= nr_it_burnin,
        
        logp_order_events_burnin(it_mcmc) = logp_current;
        if find(cnt_it_mcmc == it_mcmc),
            
            fprintf('Burnin it: %d in %d burnin its\n', ...
                it_mcmc, nr_it_burnin);
            figure(1)
            plot(logp_order_events_burnin(1:it_mcmc))
            
        end
        
    elseif it_mcmc > nr_it_burnin,
        
        if find(cnt_it_mcmc == (it_mcmc - nr_it_burnin)),
            
            fprintf('MCMC it: %d in %d mcmc its\n', ...
                it_mcmc-nr_it_burnin, nr_it_mcmc)
            
        end
        
    end
    
end
hist2_mat = zeros(nr_roi);
[ord inds_av] = sort(mean(order_events, 2));
for roi1 = 1:nr_roi,
    
    for roi2 = 1:nr_roi,
        
        hist2_mat(roi1, roi2) = sum(order_events(inds_av(roi1), :) == roi2);
        
    end
    
end
subplot(131), plot(logp_order_events_burnin)
subplot(132), plot(logp_order_events)
subplot(133), imagesc(hist2_mat)

%% Post-processing
I_maxp_mcmc = find(logp_order_events == max(logp_order_events));
order_events_max = order_events(:, I_maxp_mcmc(1));
atrophy_model_max = zeros(nr_roi, nr_roi);
for roi = 1:nr_roi,
    
    atrophy_model_max(order_events_max==roi, roi:end) = 1;
    
end

p_atrophy_model_max = zeros(nr_roi, nr_pat);
class_pat = zeros(nr_pat, 1);
for pat = 1:nr_pat,
    
    p_atrophy_local = p_atrophy(:, pat);
    p_atrophy_model_local = zeros(nr_roi, nr_roi);
    for roi = 1:nr_roi,
        
        I_atrophy = find(atrophy_model_max(:, roi) == 1);
        I_noatrophy = find(atrophy_model_max(:, roi) == 0);
        p_atrophy_model_local(I_atrophy, roi) = ...
            p_atrophy_local(I_atrophy);
        p_atrophy_model_local(I_noatrophy, roi) = ...
            1 - p_atrophy_local(I_noatrophy);
        
    end
    p_atrophy_model_local(p_atrophy_model_local == 0) = eps;
    class_pat(pat) = find(sum(log(p_atrophy_model_local)) == max(sum(log(p_atrophy_model_local))));
    p_atrophy_model_current(:, pat) = p_atrophy_model_local(:, class_pat(pat));
    
end

order_events_mean = mean(order_events, 2);
file_segm = 'img_segm.nii';
V_segm = spm_vol(file_segm);
img_segm = spm_read_vols(V_segm);
img_order = zeros(V_segm.dim);
for roi = 1:nr_roi,
    
    img_order(I_label_select{roi}) = ...
        round(100*order_events_mean(roi));
    
end

[p, file_data_f, d1, d2] = fileparts(file_data);
V_order = V_segm;
V_order.fname = sprintf('img_order_%s_reg%d_mixt%d.nii', ...
    file_data_f, flag_reg, flag_mixt);
spm_create_vol(V_order);
spm_write_vol(V_order, img_order);

results_struct.img_order = img_order;
results_struct.order_events = order_events;
results_struct.logp_order_events = logp_order_events;
results_struct.order_events_max = order_events_max;
results_struct.class_pat = class_pat;
results_struct.flag_reg = flag_reg;
results_struct.data_struct = data_struct;

eval(sprintf('save results_%s_reg%d_mixt%d results_struct', ...
    file_data_f, flag_reg, flag_mixt))

