clear all
close all

file_data = spm_select(1, 'mat');
load(file_data);

atrophy_patient(:, :, 1) = ...
    squeeze(mean(data_struct.data_Cth_patient, 2));
atrophy_control(:, :, 1) = ...
    squeeze(mean(data_struct.data_Cth_control, 2));
atrophy_patient(:, :, 2) = ...
    squeeze(mean(data_struct.data_vol_patient, 2));
atrophy_control(:, :, 2) = ...
    squeeze(mean(data_struct.data_vol_control, 2));
  
patient_id = data_struct.patient_id;
patient_tp = data_struct.patient_tp;
for flg_Cthvol = 1:2,
    
    cnt = 1;
    for id = unique(patient_id),
        
        I_pat = find(patient_id == id);
        I_baseline = find(patient_tp(I_pat) == 1);
        atrophy_baseline = atrophy_patient(:, I_pat(I_baseline), flg_Cthvol);
        for t = 1:length(I_pat),
            
            if t ~= I_baseline,
                
                atrophy_patient_rel2baseline(:, cnt, flg_Cthvol) = ...
                    (atrophy_patient(:, I_pat(t), flg_Cthvol) - ...
                    atrophy_patient(:, I_pat(I_baseline), flg_Cthvol))./ ...
                    atrophy_patient(:, I_pat(I_baseline), flg_Cthvol);
                dummy(cnt) = id;
                cnt = cnt + 1;

            end
            
        end
        
    end
    control_id = data_struct.control_id;
    control_tp = data_struct.control_tp;
    cnt = 1;
    for id = unique(control_id),
        
        I_pat = find(control_id == id);
        I_baseline = find(control_tp(I_pat) == 1);
        atrophy_baseline = atrophy_control(:, I_pat(I_baseline));
        for t = 1:length(I_pat),
            
            if t ~= I_baseline,
                
                atrophy_control_rel2baseline(:, cnt, flg_Cthvol) = ...
                    (atrophy_control(:, I_pat(t), flg_Cthvol) - ...
                    atrophy_control(:, I_pat(I_baseline), flg_Cthvol))./ ...
                    atrophy_control(:, I_pat(I_baseline), flg_Cthvol);
                cnt = cnt + 1;
                
            end
            
        end
        
    end
    
end
      
[nr_roi nr_pat, d] = size(atrophy_patient_rel2baseline);

p_atrophy = zeros(nr_roi, nr_pat, 2);
for flg_Cthvol = 1:2,
    
    for pat = 1:nr_pat,
        
        for roi = 1:nr_roi,
            
            [h, p_atrophy(roi, pat, flg_Cthvol)] = ...
                ttest(squeeze(atrophy_control_rel2baseline(roi, :, flg_Cthvol)), ...
                atrophy_patient_rel2baseline(roi, pat, flg_Cthvol), [], 'left');
            
        end
        
    end
    
end

%% Performing MCMC

nr_it_mcmc = 1e5;
nr_it_burnin = 1e5;
nr_it_hillclimb = 1e5;
thinning = 1e2;

hist2_mat = zeros([nr_roi nr_roi 4 2]);
order_events = zeros(nr_roi, nr_it_mcmc/thinning, 4, 2);
order_events_max = zeros(nr_roi, 4, 2);
logp_order_events = zeros(nr_it_mcmc/thinning, 4, 2);
logp_order_events_max = zeros(4, 2);
for flg_Cthvol = 1:2, 
    
    for l = [2 4],
        
        [order_events(:, :, l, flg_Cthvol), logp_order_events(:, l, flg_Cthvol), ...
            logp_order_events_max(l, flg_Cthvol), order_events_max(:, l, flg_Cthvol), ...
            diagnostics{l, flg_Cthvol}] = AtrophyModelMCMC(p_atrophy(:, :, flg_Cthvol), ...
            l, nr_it_hillclimb, nr_it_burnin, nr_it_mcmc, thinning);
        [ord inds_av] = sort(mean(order_events(:, :, l, flg_Cthvol), 2));
        for roi1 = 1:nr_roi,
            
            for roi2 = 1:nr_roi,
                
                hist2_mat(roi1, roi2, l, flg_Cthvol) = sum(order_events(inds_av(roi1), :, l, flg_Cthvol) == roi2);
                
            end
            
        end
        
    end
    
end

save dataAD_Cthvol_Results
%% Investigating spread of positions
cnt_fig = 1;
for flg_Cthvol = 1:2,

    for l = [2 4],
        
        subplot(2, 2, cnt_fig)
        imagesc(hist2_mat(:, :, l, flg_Cthvol))
        axis square, colorbar
        cnt_fig = cnt_fig + 1;
        
    end
    
end
        

%% Making images...
file_segm = 'img_segm.nii';
V_segm = spm_vol(file_segm);
img_segm = spm_read_vols(V_segm);
l = 2;
[p, f, d1, d2] = fileparts(file_data);
for flg_Cthvol = 1:2,
    
    img_order = zeros(size(img_segm));
    for roi1 = 1:nr_roi,
        
        cnt = 1;
        for roi2 = 1:length(data_struct.label_name),
            
            if strfind(data_struct.label_name{roi2}, ...
                    data_struct.name_regions{roi1}),
                
                idx_roi1(cnt) = roi2;
                cnt = cnt + 1;
                
            end
            
        end
        
        img_order(data_struct.I_label{idx_roi1(1)}) = ...
            round(mean(order_events(roi1, :, l, flg_Cthvol), 2));
        
    end
    V_order = V_segm;
    V_order.fname = sprintf('img_order_%s_flagreg%d.nii', ...
        f, flg_Cthvol);
    spm_create_vol(V_order);
    spm_write_vol(V_order, img_order);
    
end

%% Looking at hemispheric symmetry in results

I_lh = [1 3:36];
I_rh = [2 37:70];

figure(1), clf
cnt_fig = 1;
for flg_Cthvol = 1:2,
    
    for l = [2 4],
    
        mean_order_lh = mean(order_events(I_lh, :, l, flg_Cthvol), 2);
        std_order_lh = std(order_events(I_lh, :, l, flg_Cthvol), [], 2);
        
        mean_order_rh = mean(order_events(I_rh, :, l, flg_Cthvol), 2);
        std_order_rh = std(order_events(I_rh, :, l, flg_Cthvol), [], 2);
    
        subplot(2, 2, cnt_fig), hold on
        errorbar(mean_order_lh, mean_order_rh, std_order_rh, '.')
        herrorbar(mean_order_lh, mean_order_rh, std_order_lh, '.')
        cnt_fig = cnt_fig + 1;
        
    end
    
end

%% Looking at correspondence between registration methods

figure(2), clf, hold on
for a = 1:nr_analysis,
    
    order_events = results{a}.order_events;
    mean_order(:, a) = mean(order_events, 2);
    std_order(:, a) = std(order_events, [], 2);
    
end
errorbar(mean_order(:, 1), mean_order(:, 2), std_order(:, 2), '.k')
herrorbar(mean_order(:, 1), mean_order(:, 2), std_order(:, 1), '.k')

%% Making histograms...

figure(3), clf
inds_av = zeros(nr_roi, nr_analysis);
hist2_mat = zeros([nr_roi nr_roi nr_analysis]);
for a = 1:nr_analysis,
    
    [ord inds_av(:, a)] = sort(mean_order(:, a));
    for roi1 = 1:nr_roi,

        for roi2 = 1:nr_roi,

            hist2_mat(roi1, roi2, a) = sum(results{a}.order_events(inds_av(roi1, a), :) == roi2);

        end

    end
    subplot(1, nr_analysis, a), imagesc(log(hist2_mat(:, :, a)))
    
end

%% Checking ordering of patients within the model
patient_id = results{1}.data_struct.patient_id;
patient_tp = results{1}.data_struct.patient_tp;
figure(4), clf
for a = 1:nr_analysis,
    
    class_pat = results{a}.class_pat;
    for pid = unique(patient_id),
        
        I_pat = find(patient_id == pid);
        plot(patient_tp(I_pat), class_pat(I_pat))
        axis([min(patient_tp) max(patient_tp) 1 70])
        fprintf('Reg: %d\tPatient: %d\n', a, pid)
        pause
        
    end
    
end
    




