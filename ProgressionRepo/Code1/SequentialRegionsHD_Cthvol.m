clear all
close all

file_data = spm_select(1, 'mat');
load(file_data);

atrophy = cat(1, data_struct.data_subcort(8:10, :), ...
    squeeze(data_struct.data_vol(2:end, 1, :)), ...
    data_struct.data_subcort(28:29, :), ...
    squeeze(data_struct.data_vol(2:end, 2, :)));

group_id = data_struct.id_grp;
fu_id = data_struct.id_fu;
patient_id = data_struct.id_patient;

cnt_control = 1;
cnt_patient = 1;
patient_id_unique = unique(patient_id);
for it_id = 1:length(patient_id_unique),
    
    id = patient_id_unique(it_id);
    I_pat = find(patient_id == id);
    I_baseline = find(fu_id(I_pat) == 0);
    atrophy_baseline = atrophy(:, I_pat(I_baseline));
    
    for t = 1:length(I_pat),
        
           
            atrophy_rel2baseline_local = ...
                (atrophy(:, I_pat(t)) - ...
                atrophy(:, I_pat(I_baseline)))./ ...
                atrophy(:, I_pat(I_baseline));
            
            if group_id(I_pat(t)) == 0,
                
                atrophy_control_rel2baseline(:, cnt_control) = ...
                    atrophy_rel2baseline_local;
                fu_id_control(cnt_control) = fu_id(I_pat(t));
                group_id_control(cnt_control) = group_id(I_pat(t));
                cnt_control = cnt_control + 1;
                
            elseif group_id(I_pat(t)) > 0,
                
                atrophy_patient_rel2baseline(:, cnt_patient) = ...
                    atrophy_rel2baseline_local;
                fu_id_patient(cnt_patient) = fu_id(I_pat(t));
                group_id_patient(cnt_patient) = group_id(I_pat(t));
                cnt_patient = cnt_patient + 1;
                
            end
            
        end
        
    end
    
end

atrophy_control.data_mean = atrophy_control_rel2baseline;
atrophy_patient.data_mean = atrophy_patient_rel2baseline;
[nr_roi nr_pat, d] = size(atrophy_patient_rel2baseline);
flag_mixt = 3;
flag_filt = 1;
[p_A_D, prob_pat, prob_con] = ...
    compute_p_A_D(atrophy_control, ...
    atrophy_patient, flag_mixt, flag_filt);
p_A_D(p_A_D == 0) = eps;
p_A_D(p_A_D == 1) = 1-eps;
idx_clinevent = zeros(nr_roi, 1);

%% Performing MCMC
nr_it_mcmc = 1e4;
nr_it_burnin = 1e4;
nr_it_hillclimb = 2e3;
thinning = 1e1;
nr_hillclimb = 5;
version_mcmc = 2;

hist2_mat = zeros([nr_roi nr_roi]);
[parm_struct, diag_struct] = ...
    AtrophyModelMCMC2d(prob_pat, prob_con, idx_clinevent, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_mcmc);

save dataHD_Cthvol_Results
%% Investigating spread of positions
cnt_fig = 1;
for flag_reg = 1:2,

    for l = [2 4],
        
        subplot(2, 2, cnt_fig)
        imagesc(hist2_mat(:, :, l, flag_reg))
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
for flag_reg = 1:2,
    
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
            round(mean(order_events(roi1, :, l, flag_reg), 2));
        
    end
    V_order = V_segm;
    V_order.fname = sprintf('img_order_%s_flagreg%d.nii', ...
        f, flag_reg);
    spm_create_vol(V_order);
    spm_write_vol(V_order, img_order);
    
end

%% Looking at hemispheric symmetry in results

I_lh = [1 3:36];
I_rh = [2 37:70];

figure(1), clf
cnt_fig = 1;
for flag_reg = 1:2,
    
    for l = [2 4],
    
        mean_order_lh = mean(order_events(I_lh, :, l, flag_reg), 2);
        std_order_lh = std(order_events(I_lh, :, l, flag_reg), [], 2);
        
        mean_order_rh = mean(order_events(I_rh, :, l, flag_reg), 2);
        std_order_rh = std(order_events(I_rh, :, l, flag_reg), [], 2);
    
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
    

