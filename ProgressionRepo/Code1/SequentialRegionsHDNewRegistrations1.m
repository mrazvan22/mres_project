clear all
close all

flag_mixt = 3;
if flag_mixt == 1,
    
    flag_filt = input('Filtered controls? ');
    
else
    
    flag_filt = 1;
    
end

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

atrophy_fu12.data_mean = mean(cat(3, data_struct.data_jacc_ffd_new.data_median(I_lh, :, 1), ...
    data_struct.data_jacc_ffd_new.data_median(I_rh, :, 1)), 3);
atrophy_fu24.data_mean = mean(cat(3, data_struct.data_jacc_ffd_new.data_median(I_lh, :, 2), ...
    data_struct.data_jacc_ffd_new.data_median(I_rh, :, 2)), 3);

atrophy_control_fu12.data_mean = atrophy_fu12.data_mean(:, id_group == 0);
I_purge = find(atrophy_control_fu12.data_mean(1, :) == 0);
atrophy_control_fu12.data_mean(:, I_purge) = [];
atrophy_control_fu24.data_mean = atrophy_fu24.data_mean(:, id_group == 0);
I_purge = find(atrophy_control_fu24.data_mean(1, :) == 0);
atrophy_control_fu24.data_mean(:, I_purge) = [];
atrophy_patient_fu12.data_mean = atrophy_fu12.data_mean(:, id_group > 0);
I_purge = find(atrophy_patient_fu12.data_mean(1, :) == 0);
atrophy_patient_fu12.data_mean(:, I_purge) = [];
atrophy_patient_fu24.data_mean = atrophy_fu24.data_mean(:, id_group > 0);
I_purge = find(atrophy_patient_fu24.data_mean(1, :) == 0);
atrophy_patient_fu24.data_mean(:, I_purge) = [];

[p_A_D_fu12, prob_pat_fu12, ...
    prob_con_fu12, span_fu12] = ...
    compute_p_A_D(atrophy_control_fu12, ...
    atrophy_patient_fu12, flag_mixt, flag_filt);
[p_A_D_fu24, prob_pat_fu24, ...
    prob_con_fu24, span_fu24] = ...
    compute_p_A_D(atrophy_control_fu24, ...
    atrophy_patient_fu12, flag_mixt, flag_filt);
p_A_D = cat(2, p_A_D_fu12, p_A_D_fu24);
prob_pat = cat(2, prob_pat_fu12, prob_pat_fu24);
prob_con = cat(2, prob_con_fu12, prob_con_fu24);
span = cat(2, span_fu12, span_fu24);
idx_clinev = zeros(size(span_fu12, 1), 1);

%% Performing MCMC
nr_it_mcmc = 1e5;
nr_it_burnin = 1e5;
nr_it_hillclimb = 1e4;
thinning = 1e1;
nr_hillclimb = 10;
version_mcmc = 2;

[parm_struct, diag_struct] = ...
    AtrophyModelMCMC2d(prob_pat, prob_con, span, idx_clinev, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_mcmc);

nr_roi = size(p_A_D, 1);
nr_events = nr_roi;
hist2_mat = zeros([nr_roi nr_roi]);
[ord inds_av] = sort(mean(parm_struct.order_events, 2));
for roi1 = 1:nr_events,
    
    for roi2 = 1:nr_events,
        
        hist2_mat(roi1, roi2, flag_reg) = ...
            sum(parm_struct.order_events(inds_av(roi1), :) == roi2);
        
    end
    
end
eval(sprintf('save dataHD_Results'));

%% Investigating spread of positions
cnt_fig = 1;
for flag_reg = 1:2,

    subplot(2, 2, cnt_fig)
    imagesc(hist2_mat(:, :, flag_reg))
    axis square, colorbar
    cnt_fig = cnt_fig + 1;
    
end
        

%% Making images...
file_segm = 'img_segm.nii';
V_segm = spm_vol(file_segm);
img_segm = spm_read_vols(V_segm);
[p, f, d1, d2] = fileparts(file_data);
for flag_reg = 1:2,
    
    img_order = zeros(size(img_segm));
    for roi = 1:nr_roi,
    
        img_order(I_label_select{roi}) = ...
            round(mean(parm_struct{flag_reg}.order_events(roi, :), 2));
        
    end
    V_order = V_segm;
    V_order.fname = sprintf('img_order_%s_filt%d_flagreg%d.nii', ...
        f, flag_filt, flag_reg);
    spm_create_vol(V_order);
    spm_write_vol(V_order, img_order);
    
end

%% Looking at correspondence between registration methods

figure(2), clf, hold on
for flag_reg = 1:2,
    
    order_events = parm_struct{flag_reg}.order_events;
    mean_order(:, flag_reg) = mean(order_events, 2);
    std_order(:, flag_reg) = std(order_events, [], 2);
    
end
errorbar(mean_order(:, 1), mean_order(:, 2), std_order(:, 2), '.k')
herrorbar(mean_order(:, 1), mean_order(:, 2), std_order(:, 1), '.k')

    
%% Making goose plots 
position_roi_mean = mean(parm_struct.order_events, 2);
position_roi_std = std(parm_struct.order_events, [], 2);
[d, order_roi] = sort(position_roi_mean, 'ascend');
% First making standard goose plots (without standard deviation)...
nr_events = length(position_roi_mean);
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
        'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, 'FontSize', 8, 'FontWeight', 'demi');

end
hold on
axis equal, axis off

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
imagesc(hist2_mat(:, :, flag_reg))
map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
colormap(map_inverted)
axis square
set(gca, ...
    'YTick', [1:nr_roi], 'YTickLabel', char(mat_labels), ...
    'YGrid', 'on')
hXLabel = xlabel('Model place');
hYLabel = ylabel('Region');
set([hXLabel hYLabel], ...
    'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');



