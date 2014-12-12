%% Setting up simulated data
clear all
close all

cnt_fig = 1;
file_data = spm_select(1, 'mat');
load(file_data);

roi_id = data_struct.roi_id;
I_label_select = data_struct.I_label_select;
name_roi = data_struct.name_roi;
id_group = data_struct.id_group;

atrophy_control(:, :, 1) = data_struct.data_jacc_ffd(:, id_group == 0, 2);
atrophy_control(:, :, 2) = data_struct.data_jacc_fluid(:, id_group == 0, 2);
atrophy_patient(:, :, 1) = data_struct.data_jacc_ffd(:, id_group > 0, 2);
atrophy_patient(:, :, 2) = data_struct.data_jacc_fluid(:, id_group > 0, 2);

I_purge = find(atrophy_control(1, :, 1) == 0);
atrophy_control(:, I_purge, :) = [];
I_purge = find(atrophy_patient(1, :, 1) == 0);
atrophy_patient(:, I_purge, :) = [];

I_lh = [1:37];
I_rh = [38:74];
I_hem{1} = I_lh;
I_hem{2} = I_rh;
atrophy_patient_lhrh = cat(4, atrophy_patient(I_lh, :, :), ...
    atrophy_patient(I_rh, :, :));
atrophy_patient_lhrh = mean(atrophy_patient_lhrh, 4);
atrophy_control_lhrh = cat(4, atrophy_control(I_lh, :, :), ...
    atrophy_control(I_rh, :, :));
atrophy_control_lhrh = mean(atrophy_control_lhrh, 4);

[nr_roi nr_pat, d] = size(atrophy_patient);
nr_roi_lhrh = size(atrophy_patient_lhrh, 1);
nr_roi_hem = nr_roi_lhrh;
nr_events = nr_roi;
nr_events_lhrh = nr_roi_lhrh;
nr_events_hem = nr_events_lhrh;

p1 = zeros([nr_events nr_events 2]);
p2 = zeros([nr_events nr_events 2]);
p1_lhrh = zeros([nr_events_lhrh nr_events_lhrh 2]);
p2_lhrh = zeros([nr_events_lhrh nr_events_lhrh 2]);
p1_hem = zeros([nr_events_hem nr_events_hem 2 2]);
p2_hem = zeros([nr_events_hem nr_events_hem 2 2]);

for flag_reg = 1:2,
    
    p1(:, :, flag_reg) = PrecedenceMatrix1(atrophy_patient(:, :, flag_reg)');
    p2(:, :, flag_reg) = ...
        PrecedenceMatrix2(atrophy_patient(:, :, flag_reg)', ...
        atrophy_control(:, :, flag_reg)');
    p2(p2 == 0) = eps;
    p1_lhrh(:, :, flag_reg) = ...
        PrecedenceMatrix1(atrophy_patient_lhrh(:, :, flag_reg)');
    p2_lhrh(:, :, flag_reg) = ...
        PrecedenceMatrix2(atrophy_patient_lhrh(:, :, flag_reg)', ...
        atrophy_control_lhrh(:, :, flag_reg)');
    p2_lhrh(p2_lhrh == 0) = eps;
    for hem = 1:2,
        
        p1_hem(:, :, flag_reg, hem) = ...
            PrecedenceMatrix1(atrophy_patient(I_hem{hem}, :, flag_reg)');
        p2_hem(:, :, flag_reg, hem) = ...
            PrecedenceMatrix2(atrophy_patient(I_hem{hem}, :, flag_reg)', ...
            atrophy_control(I_hem{hem}, :, flag_reg)');       
        
    end
    p2_hem(p2_hem == 0) = eps;
    
end

%% Performing mcmc
nr_it_mcmc = 1e5;
nr_it_burnin = 1e4;
nr_it_hillclimb = 1e4;
thinning = 1e2;
nr_hillclimb = 5;

for flag_reg = 1:2,
    
    [parm_struct_lhrh{flag_reg}, diag_struct_lhrh{flag_reg}] = ...
        AtrophyModelMCMCOldStyle(p1_lhrh(:, :, flag_reg), ...
        nr_it_hillclimb, nr_hillclimb, nr_it_burnin, nr_it_mcmc, thinning);
    for hem = 1:2,
        
        [parm_struct_hem{hem}{flag_reg}, ...
            diag_struct_hem{hem}{flag_reg}] = ...
            AtrophyModelMCMCOldStyle(p1_hem(:, :, flag_reg, hem), ...
            nr_it_hillclimb, nr_hillclimb, nr_it_burnin, nr_it_mcmc, thinning);

    end
    [parm_struct{flag_reg}, diag_struct{flag_reg}] = ...
        AtrophyModelMCMCOldStyle(p1(:, :, flag_reg), ...
        nr_it_hillclimb, nr_hillclimb, nr_it_burnin, nr_it_mcmc, thinning);

end

hist2_mat = zeros([nr_events nr_events 2]);
hist2_mat_lhrh = zeros([nr_events_lhrh nr_events_lhrh 2]);
hist2_mat_hem= zeros([nr_events_hem nr_events_hem 2 2]);
for flag_reg = 1:2,
    
    [ord inds_av] = sort(mean(parm_struct{flag_reg}.order_events, 2));
    for roi1 = 1:nr_events,
        
        for roi2 = 1:nr_events,
            
            hist2_mat(roi1, roi2, flag_reg) = ...
                sum(parm_struct{flag_reg}.order_events(inds_av(roi1), :) == roi2);
            
        end
        
    end
    [ord inds_av] = sort(mean(parm_struct_lhrh{flag_reg}.order_events, 2));
    for roi1 = 1:nr_events_lhrh,
        
        for roi2 = 1:nr_events_lhrh,
            
            hist2_mat_lhrh(roi1, roi2, flag_reg) = ...
                sum(parm_struct_lhrh{flag_reg}.order_events(inds_av(roi1), :) == roi2);
            
        end
        
    end
    for hem = 1:2,
        
        [ord inds_av] = sort(mean(parm_struct_hem{hem}{flag_reg}.order_events, 2));
        for roi1 = 1:nr_events_hem,
            
            for roi2 = 1:nr_events_hem,
                
                hist2_mat_hem(roi1, roi2, flag_reg, hem) = ...
                    sum(parm_struct_hem{hem}{flag_reg}.order_events(inds_av(roi1), :) == roi2);
                
            end
            
        end
        
    end
    
end 
% eval(sprintf('save dataADResults_puolamaki_mixt%d', flag_mixt));

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
atrophy_model_max = zeros(nr_events_lhrh, nr_events_lhrh);
for roi = 1:nr_events_lhrh,
    
    atrophy_model_max(order_events_max_lhrh(:, 1)==roi, roi:end) = 1;
    
end

p_atrophy_model_max = zeros(nr_events_lhrh, nr_pat);
class_pat = zeros(nr_pat, 1);
for pat = 1:nr_pat,
    
    p_atrophy_local = p_atrophy_lhrh(:, pat);
    p_atrophy_model_local = zeros(nr_events_lhrh, nr_events_lhrh);
    for roi = 1:nr_events_lhrh,
        
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


%% Making images...
file_segm = 'img_segm.nii';
V_segm = spm_vol(file_segm);
img_segm = spm_read_vols(V_segm);
fname = 'data_HD_OldStyle';
for flag_reg = 1:2,
    
    img_order = zeros(size(img_segm));
    for roi = 1:nr_roi,
    
        img_order(I_label_select{roi}) = ...
            round(mean(parm_struct{flag_reg}.order_events(roi, :), 2));
        
    end
    V_order = V_segm;
    V_order.fname = sprintf('img_order_%s_flagreg%d.nii', ...
        fname, flag_reg);
    spm_create_vol(V_order);
    spm_write_vol(V_order, img_order);
    img_order = zeros(size(img_segm));
    for roi = 1:nr_roi_lhrh,
    
        img_order(I_label_select{I_lh(roi)}) = ...
            round(mean(parm_struct_lhrh{flag_reg}.order_events(roi, :), 2));
        
    end
    V_order = V_segm;
    V_order.fname = sprintf('img_order_lhrh_%s_flagreg%d.nii', ...
        fname, flag_reg);
    spm_create_vol(V_order);
    spm_write_vol(V_order, img_order);
    
    for hem = 1:2,
        
        img_order = zeros(size(img_segm));
        for roi = 1:nr_roi_hem,
            
            img_order(I_label_select{I_hem{hem}(roi)}) = ...
                round(mean(parm_struct_hem{flag_reg}{hem}.order_events(roi, :), 2));
            
        end
        V_order = V_segm;
        V_order.fname = sprintf('img_order_%s_hem%d_flagreg%d.nii', ...
            fname, hem, flag_reg);
        spm_create_vol(V_order);
        spm_write_vol(V_order, img_order);
        
    end

    
end

%% Looking at hemispheric symmetry in results
figure(cnt_fig), clf, cnt_fig = cnt_fig +1;
cnt_fig = 1;
for flag_reg = 1:2,
    
    mean_order_lh = mean(parm_struct_hem{1}{flag_reg}.order_events, 2);
    std_order_lh = std(parm_struct_hem{1}{flag_reg}.order_events, [], 2);   
    mean_order_rh = mean(parm_struct_hem{2}{flag_reg}.order_events, 2);
    std_order_rh = std(parm_struct_hem{2}{flag_reg}.order_events, [], 2);
   
    subplot(1, 2, cnt_fig), hold on
    errorbar(mean_order_lh, mean_order_rh, std_order_rh, '.')
    herrorbar(mean_order_lh, mean_order_rh, std_order_lh, '.')
    xlabel('Right hemisphere')
    ylabel('Left hemisphere')
    cnt_fig = cnt_fig + 1;
    
end
subplot(1, 2, 1), title('Modat-registration'), axis square
subplot(1, 2, 2), title('Freeborough-registration'), axis square

%% Looking at correspondence between registration methods

figure(cnt_fig), clf, hold on, cnt_fig = cnt_fig + 1;
mean_order = zeros(nr_events, 2);
std_order = zeros(nr_events, 2);
for flag_reg = 1:2,
    
    mean_order(:, flag_reg) = ...
        mean(parm_struct{flag_reg}.order_events, 2);
    std_order(:, flag_reg) = ...
        std(parm_struct{flag_reg}.order_events, [], 2);
    
end
subplot(1, 2, 1),
errorbar(mean_order(:, 1), mean_order(:, 2), std_order(:, 2), '.k')
herrorbar(mean_order(:, 1), mean_order(:, 2), std_order(:, 1), '.k')
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
errorbar(mean_order(:, 1), mean_order(:, 2), std_order(:, 2), '.k')
herrorbar(mean_order(:, 1), mean_order(:, 2), std_order(:, 1), '.k')
title('mean hemispheres'), axis square
xlabel('Freeborough registration')
ylabel('Modat registration')

%% Checking ordering of patients within the model
patient_id = results{1}.data_struct.patient_id;
patient_tp = results{1}.data_struct.patient_tp;
figure(cnt_fig), clf, cnt_fig = cnt_fig + 1;
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
        mean(parm_struct_hem{hem}{flag_reg}.order_events, 2);
    position_roi_hem_std(:, hem) = ...
        std(parm_struct_hem{hem}{flag_reg}.order_events, [], 2);
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

labels_lhrh{1} = labels{1};
for roi = 2:(nr_roi_lhrh+1),
    
    labels_lhrh{roi} = labels{roi+1};
    
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
name_roi{roi+1} = 'Phase 1';
name_roi{roi+2} = 'Phase 2';
name_roi{roi+3} = 'Phase 3';

name_roi_lhrh{1} = name_roi{1}(4:end);
for roi = 2:(nr_roi_lhrh+1),
    
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
imagesc(hist2_mat(:, :, l, flag_reg))
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
        'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, 'FontSize', 8, 'FontWeight', 'demi');

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
imagesc(hist2_mat_lhrh(:, :, flag_reg))
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








