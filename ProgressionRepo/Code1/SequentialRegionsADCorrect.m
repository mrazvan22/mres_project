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

I_lh = [1 3:36];
I_rh = [2 37:70];
I_hem{1} = I_lh;
I_hem{2} = I_rh;

for flag_reg = 1:2,
        
    for roi = 1:nr_roi_lhrh,
        
        for pat = 1:nr_pat,
            
            mean_local = [atrophy_patient{flag_reg}.data_mean(I_lh(roi), pat) ...
                atrophy_patient{flag_reg}.data_mean(I_rh(roi), pat)];
            %             std_local = [atrophy_patient{flag_reg}.data_std(I_lh(roi), pat) ...
            %                 atrophy_patient{flag_reg}.data_std(I_rh(roi), pat)];
            %             N_local = [atrophy_patient{flag_reg}.data_nr_vox(I_lh(roi), pat) ...
            %                 atrophy_patient{flag_reg}.data_nr_vox(I_rh(roi), pat)];
            atrophy_patient_lhrh{flag_reg}.data_mean(roi, pat) = ...
                mean(mean_local);
            %             atrophy_patient_lhrh{flag_reg}.data_std(roi, pat) = ...
            %                 sqrt((sum(N_local.*((std_local.^2) + (mean_local.^2)))/sum(N_local)) - ...
            %                 (atrophy_patient_lhrh{flag_reg}.data_mean(roi, pat).^2));
            %             atrophy_patient_lhrh{flag_reg}.data_nr_vox(roi, pat) = sum(N_local);
 
            
        end
        
    end
    
end
for flag_reg = 1:2,
    
    for roi = 1:nr_roi_lhrh,
        
        for con = 1:size(atrophy_control{1}.data_mean, 2),
            
            mean_local = [atrophy_control{flag_reg}.data_mean(I_lh(roi), con) ...
                atrophy_control{flag_reg}.data_mean(I_rh(roi), con)];
            %             std_local = [atrophy_control{flag_reg}.data_std(I_lh(roi), con) ...
            %                 atrophy_control{flag_reg}.data_std(I_rh(roi), con)];
            %             N_local = [atrophy_control{flag_reg}.data_nr_vox(I_lh(roi), con) ...
            %                 atrophy_control{flag_reg}.data_nr_vox(I_rh(roi), con)];
            atrophy_control_lhrh{flag_reg}.data_mean(roi, con) = ...
                mean(mean_local);
            %             atrophy_control_lhrh{flag_reg}.data_std(roi, con) = ...
            %                 sqrt((sum(N_local.*((std_local.^2) + (mean_local.^2)))/sum(N_local)) - ...
            %                 (atrophy_control_lhrh{flag_reg}.data_mean(roi, con).^2));
            %             atrophy_control_lhrh{flag_reg}.data_nr_vox(roi, con) = sum(N_local);
            
        end
        
    end
    
end
for flag_reg = 1:2,
    
    for hem = 1:2,
        
        atrophy_control_hem{hem}{flag_reg}.data_mean = ...
            atrophy_control{flag_reg}.data_mean(I_hem{hem}, :);
        atrophy_patient_hem{hem}{flag_reg}.data_mean = ...
            atrophy_patient{flag_reg}.data_mean(I_hem{hem}, :);
        atrophy_control_hem{hem}{flag_reg}.data_median = ...
            atrophy_control{flag_reg}.data_mean(I_hem{hem}, :);
        atrophy_patient_hem{hem}{flag_reg}.data_median = ...
            atrophy_patient{flag_reg}.data_mean(I_hem{hem}, :);
        
    end
    
end

nr_events = nr_roi;
nr_events_lhrh = nr_roi_lhrh + 3;
nr_events_hem = nr_events_lhrh;

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
span((nr_roi+1):(nr_roi+3), 1, :) = 0;
span((nr_roi+1):(nr_roi+3), 2, :) = 1;
idx_clinevent = zeros(nr_events, 1);
idx_clinevent((nr_roi+1):(nr_roi+3)) = 1;
% p_A_D(p_A_D < 1e-4) = 1e-4;
% p_A_D(p_A_D > 1-1e-4) = 1-1e-4;

p_A_D_lhrh = zeros(nr_events_lhrh, nr_pat, 2);
prob_pat_lhrh = zeros(nr_events_lhrh, nr_pat, 2);  
prob_con_lhrh = zeros(nr_events_lhrh, nr_pat, 2);  
span_lhrh = zeros(nr_events_lhrh, 2, 2);
for flag_reg = 1:2,
    
    [p_A_D_lhrh(1:nr_roi_lhrh, :, flag_reg), prob_pat_lhrh(1:nr_roi_lhrh, :, flag_reg), ...
        prob_con_lhrh(1:nr_roi_lhrh, :, flag_reg), span_lhrh(1:nr_roi_lhrh, :, flag_reg)] = ...
        compute_p_A_D(atrophy_control_lhrh{flag_reg}, ...
        atrophy_patient_lhrh{flag_reg}, flag_mixt, flag_filt);
    
end
for pat = 1:nr_pat,
    
    p_A_D_lhrh((nr_roi_lhrh+1):(nr_roi_lhrh + data_struct.cat_patient(pat)), pat, :) = 1;
    prob_pat_lhrh((nr_roi_lhrh+1):(nr_roi_lhrh + data_struct.cat_patient(pat)), pat, :) = 1;
    
end
prob_con_lhrh((nr_roi_lhrh+1):(nr_roi_lhrh+3), :, :) = 1-prob_pat_lhrh((nr_roi_lhrh+1):(nr_roi_lhrh+3), :, :);
span_lhrh((nr_roi_lhrh+1):(nr_roi_lhrh+3), 1, :) = 0;
span_lhrh((nr_roi_lhrh+1):(nr_roi_lhrh+3), 2, :) = 1;
idx_clinevent_lhrh = zeros(nr_events_lhrh, 1);
idx_clinevent_lhrh((nr_roi_lhrh+1):(nr_roi_lhrh+3)) = 1;
% p_A_D_lhrh(p_A_D_lhrh < 1e-4) = 1e-4;
% p_A_D_lhrh(p_A_D_lhrh > 1-1e-4) = 1-1e-4;

p_A_D_hem = zeros(nr_events_hem, nr_pat, 2, 2);
prob_pat_hem = zeros(nr_events_hem, nr_pat, 2, 2);  
prob_con_hem = zeros(nr_events_hem, nr_pat, 2, 2);  
span_hem = zeros(nr_events_hem, 2, 2, 2);
for flag_reg = 1:2,
    
    for hem = 1:2,
        
        [p_A_D_hem(1:nr_roi_lhrh, :, flag_reg, hem), prob_pat_hem(1:nr_roi_lhrh, :, flag_reg, hem), ...
            prob_con_hem(1:nr_roi_lhrh, :, flag_reg, hem), ...
            span_hem(1:nr_roi_lhrh, :, flag_reg, hem)] = ...
            compute_p_A_D(atrophy_control_hem{hem}{flag_reg}, ...
            atrophy_patient_hem{hem}{flag_reg}, flag_mixt, flag_filt);
        
    end
    
end
for pat = 1:nr_pat,
    
    p_A_D_hem((nr_roi_lhrh+1):(nr_roi_lhrh + data_struct.cat_patient(pat)), pat, :, :) = 1;
    prob_pat_hem((nr_roi_lhrh+1):(nr_roi_lhrh + data_struct.cat_patient(pat)), pat, :, :) = 1;
    
end
prob_con_hem((nr_roi_lhrh+1):(nr_roi_lhrh+3), :, :, :) = ...
    1-prob_pat_hem((nr_roi_lhrh+1):(nr_roi_lhrh+3), :, :, :);
span_hem((nr_roi_lhrh+1):(nr_roi_lhrh+3), 1, :, :) = 0;
span_hem((nr_roi_lhrh+1):(nr_roi_lhrh+3), 2, :, :) = 1;
idx_clinevent_hem = zeros(nr_events_lhrh, 1);
idx_clinevent_hem((nr_roi_lhrh+1):(nr_roi_lhrh+3)) = 1;
%% Performing mcmc
nr_it_mcmc = 1e5;
nr_it_burnin = 1e5;
nr_it_hillclimb = 1e4;
thinning = 1e0;
nr_hillclimb = 5;
version_mcmc = 2;

% for flag_reg = 1:2,
%     
%     [parm_struct{flag_reg}, diag_struct{flag_reg}] = ...
%         AtrophyModelMCMC2(p_A_D(:, :, flag_reg), nr_it_hillclimb, ...
%         nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_mcmc);
%     [parm_struct_lhrh{flag_reg}, diag_struct_lhrh{flag_reg}] = ...
%         AtrophyModelMCMC2(p_A_D_lhrh(:, :, flag_reg), nr_it_hillclimb, ...
%         nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_mcmc);
%     
%     %     for hem = 1:2,
%     %
%     %         [parm_struct_hem{flag_reg}{hem}, diag_struct_hem{flag_reg}{hem}] = ...
%     %             AtrophyModelMCMC2c(prob_pat_hem(:, :, flag_reg, hem), ...
%     %             prob_con_hem(:, :, flag_reg, hem), span_hem(:, :, flag_reg, hem), nr_it_hillclimb, ...
%     %             nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_mcmc);
%     %
%     %     end
%     
% end
for flag_reg = 1:2,
    
    [parm_struct{flag_reg}, diag_struct{flag_reg}] = ...
        AtrophyModelMCMC2d(prob_pat(:, :, flag_reg), prob_con(:, :, flag_reg), ...
        span(:, :, flag_reg), idx_clinevent, nr_it_hillclimb, ...
        nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_mcmc);
    [parm_struct_lhrh{flag_reg}, diag_struct_lhrh{flag_reg}] = ...
        AtrophyModelMCMC2d(prob_pat_lhrh(:, :, flag_reg), prob_con_lhrh(:, :, flag_reg), ...
        span_lhrh(:, :, flag_reg), idx_clinevent_lhrh, nr_it_hillclimb, ...
        nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_mcmc);  
    for hem = 1:2,
        
        [parm_struct_hem{flag_reg}{hem}, diag_struct_hem{flag_reg}{hem}] = ...
            AtrophyModelMCMC2d(prob_pat_hem(:, :, flag_reg, hem), ...
            prob_con_hem(:, :, flag_reg, hem), span_hem(:, :, flag_reg, hem), ... 
            idx_clinevent_hem, nr_it_hillclimb, nr_it_burnin, nr_it_mcmc, ...
            thinning, nr_hillclimb, version_mcmc);
        
    end       
    
end

for flag_reg = 1:2,
    
    atrophy_model = zeros(nr_events);
    for roi = 1:nr_events,
        
        atrophy_model(parm_struct{flag_reg}.order_events_max==roi, roi:end) = 1;
        
    end
    [logLik(flag_reg), logLikmat(:, :, flag_reg)] = ...
        logLikelihoodAtrophyModel(p_A_D(:, :, flag_reg), atrophy_model, version_mcmc);
    
end


hist2_mat = zeros([nr_events nr_events 2]);
hist2_mat_lhrh = zeros([nr_events_lhrh nr_events_lhrh 2]);
for hem = 1:2,
    
    hist2_mat_hem{hem} = zeros([nr_events_hem nr_events_hem 2]);
    
end
for flag_reg = 1:2,
    
    [ord inds_av] = sort(mean(parm_struct{flag_reg}.order_events, 2));
    for roi1 = 1:nr_events,
        
        for roi2 = 1:nr_events,
            
            hist2_mat(roi1,roi2, flag_reg) = ...
                sum(parm_struct{flag_reg}.order_events(inds_av(roi1), :) == roi2);
            
        end
        
    end
    max_order = max(parm_struct{flag_reg}.order_events(:));
    hist2_mat(max_order:end, max_order:end, flag_reg) = 0.1*nr_it_mcmc;
    [ord inds_av] = sort(mean(parm_struct_lhrh{flag_reg}.order_events, 2));
    for roi1 = 1:nr_events_lhrh,
        
        for roi2 = 1:nr_events_lhrh,
            
            hist2_mat_lhrh(roi1, roi2, flag_reg) = ...
                sum(parm_struct_lhrh{flag_reg}.order_events(inds_av(roi1), :) == roi2);
            
        end
        
    end
    max_order = max(parm_struct_lhrh{flag_reg}.order_events(:));
    hist2_mat_lhrh(max_order:end, max_order:end, flag_reg) = 0.1*nr_it_mcmc;

    for hem = 1:2,
        
        [ord inds_av] = sort(mean(parm_struct_hem{flag_reg}{hem}.order_events, 2));
        for roi1 = 1:nr_events_hem,
            
            for roi2 = 1:nr_events_hem,
                
                hist2_mat_hem{hem}(roi1, roi2, flag_reg) = ...
                    sum(parm_struct_hem{flag_reg}{hem}.order_events(inds_av(roi1), :) == roi2);
                
            end
            
        end
        max_order = max(parm_struct_hem{flag_reg}{hem}.order_events(:));
        hist2_mat_hem{hem}(max_order:end, max_order:end, flag_reg) = 0.1*nr_it_mcmc;
        
    end
    
end 
eval(sprintf('save dataADResults_Squeezed', flag_filt, version_mcmc));

%% Investigating spread of positions

subplot(2, 2, 1),
imagesc(hist2_mat(:, :, 1))
title('70 regions Modat-registration')
subplot(2, 2, 2),
imagesc(hist2_mat(:, :, 2))
title('70 regions Freeborough-registration')
subplot(2, 2, 3),
imagesc(hist2_mat_lhrh(:, :, 1))
title('35 regions Modat-registration')
subplot(2, 2, 4),
imagesc(hist2_mat_lhrh(:, :, 2))
title('35 regions Freeborough-registration')

        
%% Classifying patients...
flag_reg = 1;
atrophy_model = zeros(nr_events);
[i1, events_order_mean] = sort(mean(parm_struct{flag_reg}.order_events, 2), 'ascend');
[i2, order_events_mean] = sort(events_order_mean);
for roi = 1:nr_events,
    
    atrophy_model(order_events_mean==roi, roi:end) = 1;
    
end
[logLik, logLikmat] = ...
    logLikelihoodAtrophyModel(p_A_D(:, :, flag_reg), atrophy_model, version_mcmc);

class_pat = zeros(nr_pat, 1);
for pat = 1:nr_pat,
    
    class_pat(pat) = find(logLikmat(pat, :) == max(logLikmat(pat, :)));
    
end

%% Making images...
file_segm = 'img_segm.nii';
V_segm = spm_vol(file_segm);
img_segm = spm_read_vols(V_segm);
[p, f, d1, d2] = fileparts(file_data);
flag_reg = 1;

img_order = zeros(size(img_segm));
for roi = 1:nr_roi_lhrh,
    
    img_order(I_label_select{I_lh(roi)}) = ...
        round(mean(parm_struct_lhrh{flag_reg}.order_events(roi, :), 2));
    
end
V_order = V_segm;
V_order.fname = sprintf('img_order_lhrh_%s_flagreg%d.nii', ...
    f, flag_reg);
spm_create_vol(V_order);
spm_write_vol(V_order, img_order);
    
 
%% Making progression images for mean hemispheres.
flag_reg = 1;
progression_vec = [4 11 18 25 32];

file_segm = 'img_segm.nii';
V_segm = spm_vol(file_segm);
img_segm = spm_read_vols(V_segm);

img_stage = zeros(V_segm.dim);
for roi = 1:nr_roi_lhrh,
    
    img_stage(I_label_select{I_lh(roi)}) = 1;
    
end

pos_roi_mean = round(mean(parm_struct_lhrh{flag_reg}.order_events, 2));
pos_roi_mean = pos_roi_mean(1:nr_roi_lhrh);
[d, order_roi_lhrh] = ...
    sort(pos_roi_mean, 'ascend');

for stage = 1:nr_roi_lhrh,
    
    img_stage(I_label_select{I_lh(order_roi_lhrh(stage))}) = 2;
    if find(progression_vec == stage),
        
        V_stage = V_segm;
        V_stage.fname = sprintf('img_order_lhrh__flagreg%d_stage%d.nii', ...
            flag_reg, stage);
        spm_create_vol(V_stage);
        spm_write_vol(V_stage, img_stage);
        
    end
    
end
%% Looking at hemispheric symmetry in results
figure(1), clf, 
cnt_fig = 1;
for flag_reg = 1:2,
    
    mean_order = zeros(nr_events_hem, 2);
    for hem = 1:2,
        
        mean_order(:, hem) = ...
            mean(parm_struct_hem{flag_reg}{hem}.order_events, 2);
        
    end    
    subplot(1, 2, cnt_fig), hold on
    plot(mean_order(:, 1), mean_order(:, 2), '.k')
    xlabel('Right hemisphere')
    ylabel('Left hemisphere')
    cnt_fig = cnt_fig + 1;
    
end
subplot(1, 2, 1), title('Modat-registration'), axis square
subplot(1, 2, 2), title('Freeborough-registration'), axis square

%% Looking at correspondence between registration methods

figure(2), clf, hold on, cnt_fig = cnt_fig + 1;
mean_order = zeros(nr_events, 2);
std_order = zeros(nr_events, 2);
for flag_reg = 1:2,
    
    mean_order(:, flag_reg) = ...
        mean(parm_struct{flag_reg}.order_events, 2);
    std_order(:, flag_reg) = ...
        std(parm_struct{flag_reg}.order_events, [], 2);
    
end
subplot(1, 2, 1),
scatter(mean_order(:, 1), mean_order(:, 2), 'k.')
title('Both hemispheres'), axis square
xlabel('Freeborough registration')
ylabel('Modat registration')

mean_order = zeros(nr_events_lhrh, 2);
std_order = zeros(nr_events_lhrh, 2);
for flag_reg = 1:2,
    
    mean_order(:, flag_reg) = ...
        mean(parm_struct_lhrh{flag_reg}.order_events, 2);
    std_order(:, flag_reg) = ...
        std(parm_struct_lhrh{flag_reg}.order_events, [], 2);
    
end
subplot(1, 2, 2)
scatter(mean_order(:, 1), mean_order(:, 2), '.k')
title('mean hemispheres'), axis square
xlabel('Freeborough registration')
ylabel('Modat registration')

%% Checking ordering of patients within the model
patient_id = data_struct.patient_id;
patient_tp = data_struct.patient_tp;
figure(1), clf, hold on
for pid = unique(patient_id),
    
    I_pat = find(patient_id == pid);
    plot(patient_tp(I_pat), class_pat(I_pat), 'k', 'LineWidth', 2)
    axis([min(patient_tp) max(patient_tp) 1 nr_events]) 
        
end
hxlabel = xlabel('follow-up scan');
hylabel = ylabel('disease stage')'
set([hxlabel hylabel], ...
    'FontName', 'Helvetica', ...
    'FontSize', 12)
grid on

%% Making goose plots 
flag_reg = 1;

position_roi_mean = mean(parm_struct{flag_reg}.order_events, 2);
position_roi_std = std(parm_struct{flag_reg}.order_events, [], 2);
[d, order_roi] = sort(position_roi_mean, 'ascend');
position_roi_lhrh_mean = mean(parm_struct_lhrh{flag_reg}.order_events, 2);
position_roi_lhrh_std = std(parm_struct_lhrh{flag_reg}.order_events, [], 2);
[d, order_roi_lhrh] = sort(position_roi_lhrh_mean, 'ascend');
for hem = 1:2,
    
    position_roi_hem_mean(:, hem) = ...
        mean(parm_struct_hem{flag_reg}{hem}.order_events, 2);
    position_roi_hem_std(:, hem) = ...
        std(parm_struct_hem{flag_reg}{hem}.order_events, [], 2);
    [d, order_roi_hem(:, hem)] = ...
        sort(position_roi_hem_mean(:, hem), 'ascend');
    
end

% First making standard goose plots (without standard deviation)...
labels{1} = sprintf('%02iL', 1);
labels{2} = sprintf('%02iR', 1);
for roi = 2:35,
    
    labels{roi+1} = sprintf('%02iL', roi);
    labels{roi+35} = sprintf('%02iR', roi);
    
end
labels{nr_roi+1} = 'A';
labels{nr_roi+2} = 'B';
labels{nr_roi+3} = 'C';

for roi = 1:(nr_roi_lhrh),
    
    labels_lhrh{roi} = sprintf('%02i', roi);
    
end
labels_lhrh{nr_roi_lhrh+1} = 'A';
labels_lhrh{nr_roi_lhrh+2} = 'B';
labels_lhrh{nr_roi_lhrh+3} = 'C';

for roi = 1:70,
    
    I_str = strfind(name_roi{roi}, 'ctx-');
    if ~isempty(I_str),
        
        name_roi{roi} = name_roi{roi}(5:end);
        
    end
    I_str = strfind(name_roi{roi}, 'Left');
    if ~isempty(I_str),
        
        name_roi{roi} = ['lh-' name_roi{roi}(6:end)];
        
    end
    I_str = strfind(name_roi{roi}, 'Right');
    if ~isempty(I_str),
        
        name_roi{roi} = ['rh-' name_roi{roi}(7:end)];
        
    end    
    
end
name_roi{roi+1} = 'pre-MCI';
name_roi{roi+2} = 'MCI';
name_roi{roi+3} = 'AD';

name_roi_lhrh{1} = name_roi{1}(4:end);
for roi = 2:(nr_roi_lhrh),
    
    name_roi_lhrh{roi} = name_roi{roi+1}(4:end);
    
end
name_roi_lhrh{nr_roi_lhrh+1} = 'pre-MCI';
name_roi_lhrh{nr_roi_lhrh+2} = 'MCI';
name_roi_lhrh{nr_roi_lhrh+3} = 'AD';
           
t = 0:pi/120:2*pi;
colourRGB_regions = [30 144 255]/255;
colourRGB_diagnosis = [233 150 122]/255;
colourRGB = cat(1, repmat(colourRGB_regions, nr_roi, 1), ...
    repmat(colourRGB_diagnosis, nr_events-nr_roi, 1));
edgecolour = [repmat([0 0 1], nr_roi, 1);
    repmat([1 0 0], nr_events-nr_roi, 1)];

%% For whole brain

switch_baseline = [1:5:nr_events];
r_x = 1;
r_y = 1;
x_pos = 0;
y_pos0 = 0;
x_multiply = [0 1 1.5 2 2.5];
y_multiply = [0 1 -1 2 -2];  
dx_pos = 2;
dy_pos = 2;
figure(cnt_fig), clf, cnt_fig = cnt_fig + 1;
for roi = 1:nr_events,
    
    if find(switch_baseline == roi),
        
        cnt_y = 1;
        x_pos0 = x_pos + dx_pos;
        
    end
    x_pos = x_pos0 + x_multiply(cnt_y)*dx_pos;
    y_pos = y_pos0 + y_multiply(cnt_y)*dy_pos;
    cnt_y = cnt_y + 1;
    px = r_x*cos(t) + x_pos;
    py = r_y*sin(t) + y_pos;
    pp = patch(px, py, colourRGB(order_roi(roi), :), ...
        'EdgeColor', edgecolour(order_roi(roi), :), 'LineWidth', 2);
    text(x_pos, y_pos, labels{order_roi(roi)}, ...
        'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, 'FontSize', 8, 'FontWeight', 'demi');

end
hold on
axis equal, axis off

% Making a annotated version of the histogram
max_length_str = 0;
for roi = 1:nr_events,
    
    max_length_str = max([max_length_str length(name_roi{roi})]);
    
end
mat_labels = zeros(nr_roi, max_length_str);
for roi = 1:nr_events,
    
    length_str = length(name_roi{order_roi(roi)});
    mat_labels(roi, 1:length_str) = name_roi{order_roi(roi)};
    
end
figure(cnt_fig), clf, cnt_fig = cnt_fig + 1;
imagesc(hist2_mat(:, :, flag_reg))
map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
colormap(map_inverted)
axis square
set(gca, ...
    'YTick', [1:nr_events], 'YTickLabel', char(mat_labels), ...
    'YGrid', 'on')
hXLabel = xlabel('Model place');
hYLabel = ylabel('Region');
set([hXLabel hYLabel], ...
    'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');

    
%% For mean hemispheres
colourRGB_regions = [30 144 255]/255;
colourRGB_diagnosis = [233 150 122]/255;
colourRGB = cat(1, repmat(colourRGB_regions, nr_roi_lhrh, 1), ...
    repmat(colourRGB_diagnosis, nr_events_lhrh-nr_roi_lhrh, 1));
edgecolour = [repmat([0 0 1], nr_roi_lhrh, 1);
    repmat([1 0 0], nr_events_lhrh-nr_roi_lhrh, 1)];

switch_baseline = [1:5:nr_events_lhrh];
r_x = 1;
r_y = 1;
x_pos = 0;
y_pos0 = 0;
x_multiply = [0 1 1.5 2 2.5];
y_multiply = [0 1 -1 2 -2];  
dx_pos = 2;
dy_pos = 2;
figure(cnt_fig), clf, cnt_fig = cnt_fig +1;
for roi = 1:nr_events_lhrh,
    
    if find(switch_baseline == roi),
        
        cnt_y = 1;
        x_pos0 = x_pos + dx_pos;
        
    end
    x_pos = x_pos0 + x_multiply(cnt_y)*dx_pos;
    y_pos = y_pos0 + y_multiply(cnt_y)*dy_pos;
    cnt_y = cnt_y + 1;
    px = r_x*cos(t) + x_pos;
    py = r_y*sin(t) + y_pos;
    pp = patch(px, py, colourRGB(order_roi_lhrh(roi), :), ...
        'EdgeColor', edgecolour(order_roi_lhrh(roi), :), 'LineWidth', 2);
    text(x_pos, y_pos, labels_lhrh{order_roi_lhrh(roi)}, ...
        'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, 'FontSize', 12, 'FontWeight', 'demi');

end
hold on
axis equal, axis off

% Making a annotated version of the histogram
max_length_str = 0;
for roi = 1:nr_events_lhrh,
    
    max_length_str = max([max_length_str length(name_roi_lhrh{roi})]);
    
end
mat_labels = zeros(nr_events_lhrh, max_length_str);
for roi = 1:nr_events_lhrh,
    
    length_str = length(name_roi_lhrh{order_roi_lhrh(roi)});
    mat_labels(roi, 1:length_str) = name_roi_lhrh{order_roi_lhrh(roi)};
    
end
figure(cnt_fig), clf, cnt_fig = cnt_fig + 1;
imagesc(log(hist2_mat_lhrh(:, :, flag_reg)))
map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
colormap(map_inverted)
axis square
set(gca, ...
    'YTick', [1:nr_events_lhrh], 'YTickLabel', char(mat_labels), ...
    'YGrid', 'on')
hXLabel = xlabel('Model place');
hYLabel = ylabel('Region');
set([hXLabel hYLabel], ...
    'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
%% For seperately fitted left and right hemispheres
colourRGB_regions = [30 144 255]/255;
colourRGB_diagnosis = [233 150 122]/255;
colourRGB = cat(1, repmat(colourRGB_regions, nr_roi_hem, 1), ...
    repmat(colourRGB_diagnosis, nr_events_hem-nr_roi_hem, 1));
edgecolour = [repmat([0 0 1], nr_roi_hem, 1);
    repmat([1 0 0], nr_events_hem-nr_roi_hem, 1)];

switch_baseline = [1:5:nr_events_hem];
r_x = 1;
r_y = 1;
x_multiply = [0 1 1.5 2 2.5];
y_multiply = [0 1 -1 2 -2];  
dx_pos = 2;
dy_pos = 2;
figure(cnt_fig), clf, cnt_fig = cnt_fig + 1;
% Making a annotated version of the histogram
name_roi_hem = name_roi_lhrh;
for hem = 1:2,
    
    max_length_str = 0;
    for roi = 1:nr_events_hem,
        
        max_length_str = max([max_length_str length(name_roi_hem{roi})]);
        
    end
    mat_labels = zeros(nr_events_hem, max_length_str);
    for roi = 1:nr_events_lhrh,
        
        length_str = length(name_roi_hem{order_roi_hem(roi, hem)});
        mat_labels(roi, 1:length_str) = ...
            name_roi_hem{order_roi_hem(roi, hem)};
        
    end
    subplot(1, 2, hem)
    imagesc(hist2_mat_hem{hem}(:, :, flag_reg))
    map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
    colormap(map_inverted)
    axis square
    set(gca, ...
        'YTick', [1:nr_events_hem], 'YTickLabel', char(mat_labels), ...
        'YGrid', 'on')
    hXLabel = xlabel('Model place');
    hYLabel = ylabel('Region');
    set([hXLabel hYLabel], ...
        'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
    
end

%% Checking dependence of ordering on std atrophy of controls
std_atrophy_control = squeeze(std(atrophy_control, [], 2));
[d, orderingRegionsStd] = ...
    sort(std_atrophy_control(:, 1), 'ascend');
[d, placeRegionsStd] = sort(orderingRegionsStd, 'ascend');
placeRegionsAlgorithm = order_events_max;
flag_reg = 1;
figure(cnt_fig), clf, cnt_fig = cnt_fig +1;
subplot(121)
scatter(placeRegionsAlgorithm(1:70, flag_reg), placeRegionsStd, '.k')
axis square, title('Effect of std control atrophy on ordering')

% To check let's look at the influence of the average atrophy over all
% patients and its influence on the regions' place in the ordering
sum_p_atrophy = squeeze(sum(p_atrophy, 2));
[d, orderingRegionsSumAtrophy] = sort(sum_p_atrophy(:, flag_reg), 'descend');
[d, placeRegionsSumAtrophy] = sort(orderingRegionsSumAtrophy, 'ascend');
subplot(122)
scatter(placeRegionsAlgorithm(:, flag_reg), placeRegionsSumAtrophy, 'k.')
axis square, title('Effect of mean patient atrophy on ordering')








