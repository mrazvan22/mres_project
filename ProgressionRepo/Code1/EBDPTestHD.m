% Loading data

clear all
close all

file_data = spm_select(1, 'mat');
load(file_data);
id_roi_select = data_struct.roi_id;
I_label_select = data_struct.I_roi_select;
name_roi = data_struct.name_roi;
I_lh = 1:38;
I_rh = 39:76;
id_group = data_struct.id_group;
nr_roi = length(I_lh);

atrophy_fu12.data_mean = data_struct.data_jacc_ffd_new.data_median(:, :, 1);
I_purge = find(atrophy_fu12.data_mean(1, :) == 0);
atrophy_fu12.data_mean(:, I_purge) = [];
id_group_fu12 = id_group;
id_group_fu12(I_purge) = [];
atrophy_fu24.data_mean = data_struct.data_jacc_ffd_new.data_median(:, :, 2);
I_purge = find(atrophy_fu24.data_mean(1, :) == 0);
atrophy_fu24.data_mean(:, I_purge) = [];
id_group_fu24 = id_group;
id_group_fu24(I_purge) = [];

atrophy_fu12_lhrh.data_mean = mean(cat(3, atrophy_fu12.data_mean(I_lh, :), ...
    atrophy_fu12.data_mean(I_rh, :)), 3);
atrophy_fu24_lhrh.data_mean = mean(cat(3, atrophy_fu24.data_mean(I_lh, :), ...
    atrophy_fu24.data_mean(I_rh, :)), 3);

atrophy_control_fu12.data_mean = atrophy_fu12.data_mean(:, id_group_fu12 == 0);
atrophy_control_fu24.data_mean = atrophy_fu24.data_mean(:, id_group_fu24 == 0);
atrophy_patient_fu12.data_mean = atrophy_fu12.data_mean(:, id_group_fu12 > 0);
atrophy_patient_fu24.data_mean = atrophy_fu24.data_mean(:, id_group_fu24 > 0);
atrophy_control_fu12_lhrh.data_mean = ...
    atrophy_fu12_lhrh.data_mean(:, id_group_fu12 == 0);
atrophy_control_fu24_lhrh.data_mean = ...
    atrophy_fu24_lhrh.data_mean(:, id_group_fu24 == 0);
atrophy_patient_fu12_lhrh.data_mean = ...
    atrophy_fu12_lhrh.data_mean(:, id_group_fu12 > 0);
atrophy_patient_fu24_lhrh.data_mean = ...
    atrophy_fu24_lhrh.data_mean(:, id_group_fu24 > 0);

atrophy_control.data_mean = cat(2, atrophy_control_fu12.data_mean, atrophy_control_fu24.data_mean);
atrophy_patient.data_mean = cat(2, atrophy_patient_fu12.data_mean, atrophy_patient_fu24.data_mean);
atrophy_control_lhrh.data_mean = cat(2, atrophy_control_fu12_lhrh.data_mean, ...
    atrophy_control_fu24_lhrh.data_mean);
atrophy_patient_lhrh.data_mean = cat(2, atrophy_patient_fu12_lhrh.data_mean, ...
    atrophy_patient_fu24_lhrh.data_mean);

%% Computing likelihood
version_mixt = 4;
% [likelihood_events_fu12, gmix_fu12] = ...
%     EBDPComputeLikelihood(atrophy_patient_fu12.data_mean', ...
%     atrophy_control_fu12.data_mean', version_mixt);
% [likelihood_events_fu24, gmix_fu24] = ...
%     EBDPComputeLikelihood(atrophy_patient_fu24.data_mean', ...
%     atrophy_control_fu24.data_mean', version_mixt);
% likelihood_events = cat(2, likelihood_events_fu12, likelihood_events_fu24);
% % likelihood_events = likelihood_events_fu24;
[likelihood_events, gmix] = ...
    EBDPComputeLikelihood(atrophy_patient.data_mean', ...
    atrophy_control.data_mean', version_mixt);
[likelihood_events_lhrh, gmix_lhrh] = ...
    EBDPComputeLikelihood(atrophy_patient_lhrh.data_mean', ...
    atrophy_control_lhrh.data_mean', version_mixt);
idx_clinevent = [];

%% Performing mcmc
parm_mcmc.nr_gradient_ascent = 5;
parm_mcmc.nr_it_gradient_ascent = 2e3;
parm_mcmc.nr_it_burnin = 1e5;
parm_mcmc.nr_it_mcmc = 1e6;
parm_mcmc.interval_display = 1e2;
parm_mcmc.nr_it_check_p_false = 1e2;
parm_mcmc.range_p_false = [0.099 0.1
    0.099 0.1];
parm_mcmc.std_p_false = [1e-4 1e-4]';
parm_mcmc.idx_clinevent = idx_clinevent;
parm_mcmc.flag_sum = 2;

[parm_struct] = EBDPMCMCTest(likelihood_events, parm_mcmc);
[parm_struct_lhrh] = EBDPMCMCTest(likelihood_events_lhrh, parm_mcmc);
%% Visualization
close all
figure(1), clf
subplot(221), plot(parm_struct_lhrh.log_likelihood_gradient_ascent)
subplot(222), plot(parm_struct_lhrh.log_likelihood_mcmc)
subplot(223), imagesc(parm_struct_lhrh.event_order_mcmc);
subplot(224), plot(parm_struct_lhrh.p_false_mcmc')

nr_events = nr_roi;
hist2_mat = zeros([nr_roi nr_roi]);
[ord inds_av] = sort(mean(parm_struct_lhrh.event_order_mcmc, 2));
for roi1 = 1:nr_events,
    
    for roi2 = 1:nr_events,
        
        hist2_mat(roi1, roi2) = ...
            sum(parm_struct_lhrh.event_order_mcmc(inds_av(roi1), :) == roi2);
        
    end
    
end

position_roi_mean = mean(parm_struct_lhrh.event_order_mcmc, 2);
position_roi_std = std(parm_struct_lhrh.event_order_mcmc, [], 2);
[d, order_roi] = sort(position_roi_mean, 'ascend');
% First making standard goose plots (without standard deviation)...
for roi = 1:nr_roi,
    
    labels{roi} = sprintf('%02i', roi);
    
end
for roi = 1:nr_roi,
    
    I_str = strfind(data_struct.name_roi{roi}, 'Left');
    if isempty(I_str),
        
        I_str = strfind(data_struct.name_roi{roi}, 'Right');
        
        if isempty(I_str),
            
            I_str = strfind(data_struct.name_roi{roi}, 'ctx');
            name_roi_lhrh{roi} = data_struct.name_roi{roi}((I_str+7):end);
            
        else
            
            name_roi_lhrh{roi} = data_struct.name_roi{roi}((I_str+6):end);
            
            
        end
        
    else
        
        name_roi_lhrh{roi} = data_struct.name_roi{roi}((I_str+5):end);
        
    end
    
end

figure(2), clf
t = 0:pi/120:2*pi;
colourRGB = [30 144 255]/255;
switch_baseline = [1:5:nr_roi];
r_x = 1;
r_y = 1;
x_pos = 0;
y_pos0 = 0;
x_multiply = [0 1 1.5 2 2.5];
y_multiply = [0 1 -1 2 -2];  
dx_pos = 2;
dy_pos = 2;
for roi = 1:nr_roi,
    
    if find(switch_baseline == roi),
        
        cnt_y = 1;
        x_pos0 = x_pos + dx_pos;
        
    end
    x_pos = x_pos0 + x_multiply(cnt_y)*dx_pos;
    y_pos = y_pos0 + y_multiply(cnt_y)*dy_pos;
    cnt_y = cnt_y + 1;
    px = r_x*cos(t) + x_pos;
    py = r_y*sin(t) + y_pos;
    pp = patch(px, py, colourRGB, 'EdgeColor', [0 0 1], 'LineWidth', 2);
    text(x_pos, y_pos, labels{order_roi(roi)}, ...
        'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, 'FontSize', 12, 'FontWeight', 'demi');

end
hold on
axis equal, axis off
set(gcf, 'PaperPositionMode', 'auto');

% Making a annotated version of the histogram
max_length_str = 0;
for roi = 1:nr_roi,
    
    max_length_str = max([max_length_str length(name_roi_lhrh{roi})]);
    
end
mat_labels = zeros(nr_roi, max_length_str);
for roi = 1:nr_roi,
    
    length_str = length(name_roi_lhrh{order_roi(roi)});
    mat_labels(roi, 1:length_str) = name_roi_lhrh{order_roi(roi)};
    
end
figure(3), clf
imagesc(log(hist2_mat(:, :)))
map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
colormap(map_inverted)
axis square
set(gca, ...
    'YTick', [1:nr_roi], 'YTickLabel', char(mat_labels), ...
    'YGrid', 'on')
hXLabel = xlabel('model stage');
hYLabel = ylabel('region');
set([hXLabel hYLabel], ...
    'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperPositionMode', 'auto');

%% Making progression images for mean hemispheres.
flag_reg = 1;
progression_vec = [1:5:38];

file_segm = 'img_segm.nii';
V_segm = spm_vol(file_segm);
img_segm = spm_read_vols(V_segm);

img_stage = zeros(V_segm.dim);
for roi = 1:nr_roi,
    
    img_stage(I_label_select{I_lh(roi)}) = 1;
    
end

pos_roi_mean = round(mean(parm_struct_lhrh.event_order_mcmc, 2));
[d, order_roi_lhrh] = ...
    sort(pos_roi_mean, 'ascend');

for stage = 1:nr_roi,
    
    img_stage(I_label_select{I_lh(order_roi_lhrh(stage))}) = 2;
    if find(progression_vec == stage),
        
        V_stage = V_segm;
        V_stage.fname = sprintf('img_order_lhrh_HD_%s_stage%d.nii', ...
            date, stage);
        spm_create_vol(V_stage);
        spm_write_vol(V_stage, img_stage);
        
    end
    
end

%% Make comparison between hemispheres...
figure(4), clf, hold on
mean_order = zeros(nr_roi, 2);
std_order = zeros(nr_roi, 2);
mean_order(:, 1) = ...
    mean(parm_struct.event_order_mcmc(I_lh, :), 2);
mean_order(:, 2) = ...
    mean(parm_struct.event_order_mcmc(I_rh, :), 2);
std_order(:, 1) = ...
    std(parm_struct.event_order_mcmc(I_lh, :), [], 2);
std_order(:, 2) = ...
    std(parm_struct.event_order_mcmc(I_rh, :), [], 2);

hData = line(mean_order(:, 1), mean_order(:, 2));
herrorbar1 = herrorbar(mean_order(:, 1), mean_order(:, 2), ...
    std_order(:, 1), '.');
herrorbar2 = errorbar(mean_order(:, 1), mean_order(:, 2), ...
    std_order(:, 2), '.');

axis square
hXlabel = xlabel('Event position right hemisphere');
hYlabel = ylabel('Event position left hemisphere');
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:10:70   , ...
    'XTick'       , 0:10:70   , ...
    'LineWidth'   , 1         );
set(hData, ...
    'LineStyle'       , 'none'      , ...
    'Marker'          , 'o', ...
    'MarkerSize'      , 4           , ...
    'MarkerEdgeColor' , 'none'      , ...
    'MarkerFaceColor' , [0 0 0] );
set([herrorbar1; herrorbar2], ...
    'Color', [0 0 0]);
set([hXlabel hYlabel], ...
    'FontName', 'Helvetica', ...
    'FontSize', 12)
set(gcf, 'PaperPositionMode', 'auto');
eval(sprintf('print -dpng comparehemispheres_HD_%s.png', date));
