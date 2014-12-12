%% Loading data
% This stuff that isn't immediately important, I hope the variable names
% are self-explanatory anyway
clear all
close all

cnt_fig = 1;
file_data = 'd:\temp\Progression\data_fAD\data_AD.mat'; % This needs to changed...
load(file_data);

id_roi_select = data_struct.id_roi_select;
I_label_select = data_struct.I_label_select;
name_roi = data_struct.name_roi;
[nr_roi nr_pat] = size(data_struct.data_patient{1}.data_median);
nr_roi_lhrh = nr_roi/2;
nr_events = nr_roi + 3;
nr_events_lhrh = nr_roi_lhrh +3;
I_lh = [1 3:36];
I_rh = [2 37:70];

% There are two sorts of registration that were used on this data: 1 = the
% Modat method, 2 = the Freeborough method. That's what flag_reg
% indicates...
for flag_reg = 1:2,
    
    atrophy_patient{flag_reg}.data_mean = ...
        data_struct.data_patient{flag_reg}.data_median;
    atrophy_control{flag_reg}.data_mean = ...
        data_struct.data_control{flag_reg}.data_median;
    for roi = 1:nr_roi_lhrh,
        
        for pat = 1:nr_pat,
            
            mean_local = [atrophy_patient{flag_reg}.data_mean(I_lh(roi), pat) ...
                atrophy_patient{flag_reg}.data_mean(I_rh(roi), pat)];
            atrophy_patient_lhrh{flag_reg}.data_mean(roi, pat) = ...
                mean(mean_local);
            
        end
        
    end
    for roi = 1:nr_roi_lhrh,
        
        for con = 1:size(atrophy_control{flag_reg}.data_mean, 2),
            
            mean_local = [atrophy_control{flag_reg}.data_mean(I_lh(roi), con) ...
                atrophy_control{flag_reg}.data_mean(I_rh(roi), con)];
            atrophy_control_lhrh{flag_reg}.data_mean(roi, con) = ...
                mean(mean_local);
            
        end
        
    end

end

%% Compute data Likelihood given event/no event
version_mixt = 4; % Choosing the mixture of the Gaussian and uniform distributions
opts = foptions;
opts(14) = 1e3;
for flag_reg = 1:2,
    
    likelihood_events{flag_reg} = zeros(nr_events, nr_pat, 2);
    % likelihood_events{flag_reg}(:, :, 1) is the likelihood given that no
    % event has occurred, likelihood_events{flag_reg}(:, :, 2) is the
    % likelihood given that an event has occurred
    % imagesc(likelihood_events{flag_reg}(:, :,
    % 2)./sum(likelihood_events{flag_reg}, 3) gives something equivalent
    % to a posterior probability that an event has occurred
    
    [likelihood_events{flag_reg}(1:nr_roi, :, :), gmix_roi] = ...
        EBDPComputeLikelihood(atrophy_patient{flag_reg}.data_mean', ...
        atrophy_control{flag_reg}.data_mean', version_mixt);
    
    % The following is a bit complicated. We're dealing with clinical
    % events and atrophy events here. The likelihood of atrophy events,
    % given an event has/has not occurred is determined with EBDPComputeLikelihood
    % For the clinical events, no such likelihood exists and we have to
    % make up some values to get the MCMC fitting going. To this end, I'm
    % determining some 'typical' values from the distribution of likelihood
    % values, given that events have or have not occurred. I'm doing this
    % by fitting a mixture of Gaussians to these values.
    % If you don't like this, I'm also fitting the EBDP model to the
    % atrophy values only, with very similar results.
    
    gmix_noevents = gmm(1, 2, 'full');
    l_local = likelihood_events{flag_reg}(:, :, 1);
    gmix_noevents = gmmem(gmix_noevents, l_local(:), opts);
    max_noevents = max(gmix_noevents.centres);
    min_noevents = min(gmix_noevents.centres);
    gmix_events = gmm(1, 2, 'full');
    l_local = likelihood_events{flag_reg}(:, :, 2);
    max_events = mean(l_local(l_local ~= 0));
    min_events = min_noevents;
    for pat = 1:nr_pat,
        
        likelihood_events{flag_reg}((nr_roi + data_struct.cat_patient(pat) + 1):end, ...
            pat, 1) = max_noevents;
        likelihood_events{flag_reg}((nr_roi + data_struct.cat_patient(pat) + 1):end, ...
            pat, 2) = min_events;
        likelihood_events{flag_reg}((nr_roi + 1):(nr_roi + data_struct.cat_patient(pat)), ...
            pat, 1) = min_noevents;
        likelihood_events{flag_reg}((nr_roi + 1):(nr_roi + data_struct.cat_patient(pat)), ...
            pat, 2) = max_events;
        likelihood_events{flag_reg}(nr_roi + 1, pat, 1) = 0;
        likelihood_events{flag_reg}(nr_roi + 1, pat, 2) = 1;
        
    end

    % The variables with _lhrh are likelihood values when left and right
    % hemisphere atrophy values are averaged. This is mainly done to reduce
    % the number of stages in the model and make the results somewhat more
    % presentable
    likelihood_events_lhrh{flag_reg} = zeros(nr_events_lhrh, nr_pat, 2);
    [likelihood_events_lhrh{flag_reg}(1:nr_roi_lhrh, :, :), gmix_roi] = ...
        EBDPComputeLikelihood(atrophy_patient_lhrh{flag_reg}.data_mean', ...
        atrophy_control_lhrh{flag_reg}.data_mean', version_mixt);
    gmix_noevents = gmm(1, 2, 'full');
    l_local = likelihood_events_lhrh{flag_reg}(:, :, 1);
    gmix_noevents = gmmem(gmix_noevents, l_local(:), opts);
    max_noevents = max(gmix_noevents.centres);
    min_noevents = min(gmix_noevents.centres);
    gmix_events = gmm(1, 2, 'full');
    l_local = likelihood_events_lhrh{flag_reg}(:, :, 2);
    max_events = mean(l_local(l_local ~= 0));
    min_events = min_noevents;
    for pat = 1:nr_pat,
         
        likelihood_events_lhrh{flag_reg}((nr_roi_lhrh + ...
            data_struct.cat_patient(pat) + 1):end, pat, 1) = max_noevents;
        likelihood_events_lhrh{flag_reg}((nr_roi_lhrh + ...
            data_struct.cat_patient(pat) + 1):end, pat, 2) = min_events;
        likelihood_events_lhrh{flag_reg}((nr_roi_lhrh+1):...
            (nr_roi_lhrh + data_struct.cat_patient(pat)), pat, 1) = min_noevents;
        likelihood_events_lhrh{flag_reg}((nr_roi_lhrh+1):...
            (nr_roi_lhrh + data_struct.cat_patient(pat)), pat, 2) = max_events;
        likelihood_events_lhrh{flag_reg}(nr_roi_lhrh + 1, pat, 1) = 0;
        likelihood_events_lhrh{flag_reg}(nr_roi_lhrh + 1, pat, 2) = 1;
        
    end
    
end
idx_clinevent = []; % This variable is outdated

% This was apparently needed to make things work properly. Some of the
% likelihood values are zero, which the MCMC algorithm doesn't respond to
% particularly well...
likelihood_events{2}(likelihood_events{2} < 1e-3) = 1e-3;
likelihood_events_lhrh{2}(likelihood_events_lhrh{2} < 1e-3) = 1e-3;

% This is fitting the model to the atrophy likelihoods only
likelihood_events_noclin = likelihood_events_lhrh{1}(1:nr_roi_lhrh, :, :);
%% Performing mcmc
parm_mcmc.nr_gradient_ascent = 5;
parm_mcmc.nr_it_gradient_ascent = 2e3;
parm_mcmc.nr_it_burnin = 1e5;
parm_mcmc.nr_it_mcmc = 1e5;
parm_mcmc.interval_display = 1e2;
parm_mcmc.flag_sum = 2; % This should always be set to 2.

% The following parameters are used in a deprecated part of the algorithm.
% This part is effectively switched off when idx_clinevent = []
parm_mcmc.nr_it_check_p_false = 1e2;
parm_mcmc.range_p_false = [0.05 0.1
    0.05 0.1];
parm_mcmc.std_p_false = [1e-4 1e-4]';
parm_mcmc.idx_clinevent = idx_clinevent; 

for flag_reg = 1:2,
    
    parm_struct{flag_reg} = EBDPMCMCTest(likelihood_events{flag_reg}, ...
        parm_mcmc);
    parm_struct_lhrh{flag_reg} = EBDPMCMCTest(likelihood_events_lhrh{flag_reg}, ...
        parm_mcmc);
    
end
parm_struct_noclin = EBDPMCMCTest(likelihood_events_noclin, ...
    parm_mcmc);

flag_reg = 1;
parm_struct_hem{1} = EBDPMCMCTest(likelihood_events{flag_reg}(I_lh, :, :), ...
    parm_mcmc);
parm_struct_hem{2} = EBDPMCMCTest(likelihood_events{flag_reg}(I_rh, :, :), ...
    parm_mcmc);

  
% it_vec = 1:1e2:1e5;
% for flag_reg = 1:2,
%     
%     event_order = parm_struct{flag_reg}.event_order_mcmc;
%     event_order_unique = unique(event_order', 'rows')';
%     nr_event_order_unique = zeros(size(event_order_unique, 2), 1);
%     for it = 1:size(event_order_unique, 2),
%         
%         mat_comp = repmat(event_order_unique(:, it), 1, 1e5);
%         mat_comp = event_order == mat_comp;
%         nr_event_order_unique(it) = sum(min(mat_comp, [], 1));
%         if find(it_vec == it),
%             
%             fprintf('it: %d in %d it\n', it, size(event_order_unique, 2))
%             
%         end
%         
%     end
%     event_order_mode{flag_reg} = event_order_unique(:, ...
%         find(nr_event_order_unique == max(nr_event_order_unique)));
%     
%     event_order = parm_struct_lhrh{flag_reg}.event_order_mcmc;
%     event_order_unique = unique(event_order', 'rows')';
%     nr_event_order_unique = zeros(size(event_order_unique, 2), 1);
%     for it = 1:size(event_order_unique, 2),
%         
%         mat_comp = repmat(event_order_unique(:, it), 1, 1e5);
%         mat_comp = event_order == mat_comp;
%         nr_event_order_unique(it) = sum(min(mat_comp, [], 1));
%         if find(it_vec == it),
%             
%             fprintf('it: %d in %d it\n', it, size(event_order_unique, 2))
%             
%         end
%         
%     end
%     event_order_lhrh_mode{flag_reg} = event_order_unique(:, ...
%         find(nr_event_order_unique == max(nr_event_order_unique)));
%     
%     
% end
%% Visualizing results
% Some diagnostics first
flag_reg = 1;

close all
cnt_fig = 1;
figure(1), clf
subplot(221), plot(parm_struct_lhrh{flag_reg}.log_likelihood_gradient_ascent)
subplot(222), plot(parm_struct_lhrh{flag_reg}.log_likelihood_mcmc)
subplot(223), imagesc(parm_struct_lhrh{flag_reg}.event_order_mcmc);
subplot(224), plot(parm_struct_lhrh{flag_reg}.p_false_mcmc')
cnt_fig = cnt_fig + 1;   
    
% Making histogram
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

% Making goose plot
position_roi_lhrh_mean = mean(parm_struct_lhrh{flag_reg}.event_order_mcmc, 2);
position_roi_lhrh_std = std(parm_struct_lhrh{flag_reg}.event_order_mcmc, [], 2);

t = 0:pi/120:2*pi;
colourRGB_regions = [30 144 255]/255;
colourRGB_diagnosis = [233 150 122]/255;
colourRGB = cat(1, repmat(colourRGB_regions, nr_roi_lhrh, 1), ...
    repmat(colourRGB_diagnosis, nr_events_lhrh-nr_roi_lhrh, 1));
edgecolour = [repmat([0 0 1], nr_roi_lhrh, 1);
    repmat([1 0 0], nr_events_lhrh-nr_roi_lhrh, 1)];
[d, order_roi_lhrh] = ...
    sort(position_roi_lhrh_mean, 'ascend');

switch_baseline = [1:5:nr_events_lhrh];
r_x = 1;
r_y = 1;
x_pos = 0;
y_pos0 = 0;
x_multiply = [0 1 1.5 2 2.5];
y_multiply = [0 1 -1 2 -2];  
dx_pos = 2;
dy_pos = 2;
figure(2), clf, 
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
set(gcf, 'PaperPositionMode', 'auto');

hist2_mat_lhrh = zeros([nr_events_lhrh nr_events_lhrh]);
[ord inds_av] = sort(mean(parm_struct_lhrh{flag_reg}.event_order_mcmc, 2));
for roi1 = 1:nr_events_lhrh,
    
    for roi2 = 1:nr_events_lhrh,
        
        hist2_mat_lhrh(roi1, roi2) = ...
            sum(parm_struct_lhrh{flag_reg}.event_order_mcmc(inds_av(roi1), :) == roi2);
        
    end
    
end

[d, order_roi_lhrh] = sort(position_roi_lhrh_mean, 'ascend');
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
figure(3), clf, 
imagesc(log(hist2_mat_lhrh(:, :, flag_reg)))
map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
colormap(map_inverted)
axis square
set(gca, ...
    'YTick', [1:nr_events_lhrh], 'YTickLabel', char(mat_labels), ...
    'YGrid', 'on')
hXLabel = xlabel('model stage');
hYLabel = ylabel('region');
set([hXLabel hYLabel], ...
    'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'PaperPositionMode', 'auto');

%% Making progression images for mean hemispheres.
flag_reg = 1;
progression_vec = 5:5:38;

file_segm = 'img_segm.nii';
V_segm = spm_vol(file_segm);
img_segm = spm_read_vols(V_segm);

img_stage = zeros(V_segm.dim);
for roi = 1:nr_roi_lhrh,
    
    img_stage(I_label_select{I_lh(roi)}) = 1;
    
end

[d, order_roi_lhrh] = ...
    sort(position_roi_lhrh_mean(1:nr_roi_lhrh), 'ascend');
for stage = 1:nr_roi_lhrh,
    
    img_stage(I_label_select{I_lh(order_roi_lhrh(stage))}) = 2;
    if find(progression_vec == stage),
        
        V_stage = V_segm;
        V_stage.fname = sprintf('img_order_fAD_%s_lhrh_stage%d.nii', date, ...
            stage);
        spm_create_vol(V_stage);
        spm_write_vol(V_stage, img_stage);
        
    end
    
end

%% Looking at correspondence between registration methods
hist2_mat_comparereg = zeros([nr_events_lhrh nr_events_lhrh 2]);
flag_reg_follow = 1;
[ord inds_av] = sort(mean(parm_struct_lhrh{flag_reg_follow}.event_order_mcmc, 2));
for flag_reg = 1:2,
    
    for roi1 = 1:nr_events_lhrh,
        
        for roi2 = 1:nr_events_lhrh,
            
            hist2_mat_comparereg(roi1, roi2, flag_reg) = ...
                sum(parm_struct_lhrh{flag_reg}.event_order_mcmc(inds_av(roi1), :) == roi2);
            
        end
        
    end
    
end

position_roi_comparereg_mean = ...
    mean(parm_struct_lhrh{flag_reg_follow}.event_order_mcmc, 2);
[d, order_roi_comparereg] = sort(position_roi_comparereg_mean, 'ascend');
% Making a annotated version of the histogram
max_length_str = 0;
for roi = 1:nr_events_lhrh,
    
    max_length_str = max([max_length_str length(name_roi_lhrh{roi})]);
    
end
mat_labels_comparereg = zeros(nr_events_lhrh, max_length_str);
for roi = 1:nr_events_lhrh,
    
    length_str = length(name_roi_lhrh{order_roi_comparereg(roi)});
    mat_labels_comparereg(roi, 1:length_str) = name_roi_lhrh{order_roi_comparereg(roi)};
    
end

figure(4), clf
map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
for flag_reg = 1:2,
    
    subplot(1, 2, flag_reg),
    imagesc(log(hist2_mat_comparereg(:, :, flag_reg)))
    colormap(map_inverted)
    axis square
    set(gca, ...
        'YTick', [1:nr_events_lhrh], 'YTickLabel', char(mat_labels_comparereg), ...
        'YGrid', 'on')
    hXLabel = xlabel('model stage');
    hYLabel = ylabel('region');
    set([hXLabel hYLabel], ...
        'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
    set(gcf, 'PaperPositionMode', 'auto');
    
end

figure(5), clf, hold on, 

mean_order = zeros(nr_events_lhrh, 2);
std_order = zeros(nr_events_lhrh, 2);
for flag_reg = 1:2,
    
    mean_order(:, flag_reg) = ...
        mean(parm_struct_lhrh{flag_reg}.event_order_mcmc, 2);
    std_order(:, flag_reg) = ...
        std(parm_struct_lhrh{flag_reg}.event_order_mcmc, [], 2);
    
end
hData = line(mean_order(:, 1), mean_order(:, 2));
herrorbar1 = herrorbar(mean_order(:, 1), mean_order(:, 2), std_order(:, 1), '.');
herrorbar2 = errorbar(mean_order(:, 1), mean_order(:, 2), std_order(:, 2), '.');

axis square
hXlabel = xlabel('Event position Freeborough registration');
hYlabel = ylabel('Event position Modat registration');
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
    'YTick'       , 0:10:35   , ...
    'XTick'       , 0:10:35   , ...
    'LineWidth'   , 1         );
set(hData, ...
    'LineStyle'       , 'none'      , ...
    'Marker'          , 'o', ...
    'MarkerSize'      , 4           , ...
    'MarkerEdgeColor' , 'none'      , ...
    'MarkerFaceColor' , [0 0 0] );
set([herrorbar1; herrorbar2], ...
    'Color' , [0 0 0] );

set(hTitle, ...
    'FontName', 'Helvetica', ...
    'FontSize', 14, ...
    'FontWeight', 'Bold')
set([hXlabel hYlabel], ...
    'FontName', 'Helvetica', ...
    'FontSize', 12)
set(gcf, 'PaperPositionMode', 'auto');
eval(sprintf('print -dpng compareregistrations_fAD_%s.png', date));

%% Looking at hemispheric correspondence...
figure(6), clf, hold on
flag_reg = 1;
nr_roi_hem = nr_roi_lhrh;

mean_order = zeros(nr_roi_lhrh, 2);
std_order = zeros(nr_roi_lhrh, 2);
mean_order(:, 1) = ...
    mean(parm_struct{flag_reg}.event_order_mcmc(I_lh, :), 2);
mean_order(:, 2) = ...
    mean(parm_struct{flag_reg}.event_order_mcmc(I_rh, :), 2);
std_order(:, 1) = ...
    std(parm_struct{flag_reg}.event_order_mcmc(I_lh, :), [], 2);
std_order(:, 2) = ...
    std(parm_struct{flag_reg}.event_order_mcmc(I_rh, :), [], 2);

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
    'Color', [0 0 0])
set([hXlabel hYlabel], ...
    'FontName', 'Helvetica', ...
    'FontSize', 12)
set(gcf, 'PaperPositionMode', 'auto');
eval(sprintf('print -dpng comparehemispheres_fAD_%s.png', date));


%% Checking ordering of patients within the model
patient_id = data_struct.patient_id;
patient_tp = data_struct.patient_tp;
class_pat = parm_struct_lhrh{1}.class_pat;
figure(6), clf, hold on
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
set(gcf, 'PaperPositionMode', 'auto');
eval(sprintf('print -dpng patientclassifications_fAD_%s.png', date));

%% Making figure for results with no clinical events
position_roi_noclin_mean = mean(parm_struct_noclin.event_order_mcmc, 2);
position_roi_noclin_std = std(parm_struct_noclin.event_order_mcmc, [], 2);

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
    'Color', [0 0 0])
set([hXlabel hYlabel], ...
    'FontName', 'Helvetica', ...
    'FontSize', 12)
set(gcf, 'PaperPositionMode', 'auto');

t = 0:pi/120:2*pi;
colourRGB_regions = [30 144 255]/255;
colourRGB_diagnosis = [233 150 122]/255;
colourRGB = cat(1, repmat(colourRGB_regions, nr_roi_lhrh, 1), ...
    repmat(colourRGB_diagnosis, nr_events_lhrh-nr_roi_lhrh, 1));
edgecolour = [repmat([0 0 1], nr_roi_lhrh, 1);
    repmat([1 0 0], nr_events_lhrh-nr_roi_lhrh, 1)];
[d, order_roi_noclin] = ...
    sort(position_roi_noclin_mean, 'ascend');

switch_baseline = [1:5:nr_events_lhrh];
r_x = 1;
r_y = 1;
x_pos = 0;
y_pos0 = 0;
x_multiply = [0 1 1.5 2 2.5];
y_multiply = [0 1 -1 2 -2];  
dx_pos = 2;
dy_pos = 2;
figure(7), clf, 
for roi = 1:nr_roi_lhrh,
    
    if find(switch_baseline == roi),
        
        cnt_y = 1;
        x_pos0 = x_pos + dx_pos;
        
    end
    x_pos = x_pos0 + x_multiply(cnt_y)*dx_pos;
    y_pos = y_pos0 + y_multiply(cnt_y)*dy_pos;
    cnt_y = cnt_y + 1;
    px = r_x*cos(t) + x_pos;
    py = r_y*sin(t) + y_pos;
    pp = patch(px, py, colourRGB(order_roi_noclin(roi), :), ...
        'EdgeColor', edgecolour(order_roi_noclin(roi), :), 'LineWidth', 2);
    text(x_pos, y_pos, labels_lhrh{order_roi_noclin(roi)}, ...
        'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, 'FontSize', 12, 'FontWeight', 'demi');

end
hold on
axis equal, axis off
set(gcf, 'PaperPositionMode', 'auto');

%% Leave on out experiment
version_mixt = 4;
opts = foptions;
opts(14) = 1e3;
idx_clinevent = [];

parm_mcmc.nr_gradient_ascent = 5;
parm_mcmc.nr_it_gradient_ascent = 2e3;
parm_mcmc.nr_it_burnin = 1e5;
parm_mcmc.nr_it_mcmc = 1e5;
parm_mcmc.interval_display = 1e2;
parm_mcmc.nr_it_check_p_false = 1e2;
parm_mcmc.range_p_false = [0.05 0.1
    0.05 0.1];
parm_mcmc.std_p_false = [1e-4 1e-4]';
parm_mcmc.idx_clinevent = idx_clinevent;
parm_mcmc.flag_sum = 2;

class_pat_loo = zeros(nr_pat, 1);

flag_reg = 1;
nr_pat_loo = nr_pat -1;
likelihood_events_tot = zeros(nr_events_lhrh, nr_pat, 2);
[likelihood_events_tot(1:nr_roi_lhrh, :, :), gmix_roi] = ...
    EBDPComputeLikelihood(atrophy_patient_lhrh{flag_reg}.data_mean', ...
    atrophy_control_lhrh{flag_reg}.data_mean', version_mixt);
gmix_noevents = gmm(1, 2, 'full');
l_local = likelihood_events_tot(:, :, 1);
gmix_noevents = gmmem(gmix_noevents, l_local(:), opts);
max_noevents = max(gmix_noevents.centres);
min_noevents = min(gmix_noevents.centres);
gmix_events = gmm(1, 2, 'full');
l_local = likelihood_events_tot(:, :, 2);
max_events = mean(l_local(l_local ~= 0));
min_events = min_noevents;
for pat = 1:nr_pat,
    
    likelihood_events_tot((nr_roi_lhrh + ...
        data_struct.cat_patient(pat) + 1):end, pat, 1) = max_noevents;
    likelihood_events_tot((nr_roi_lhrh + ...
        data_struct.cat_patient(pat) + 1):end, pat, 2) = min_events;
    likelihood_events_tot((nr_roi_lhrh+1):...
        (nr_roi_lhrh + data_struct.cat_patient(pat)), pat, 1) = min_noevents;
    likelihood_events_tot((nr_roi_lhrh+1):...
        (nr_roi_lhrh + data_struct.cat_patient(pat)), pat, 2) = max_events;
    likelihood_events_tot(nr_roi_lhrh + 1, pat, 1) = 0;
    likelihood_events_tot(nr_roi_lhrh + 1, pat, 2) = 1;
    
end

for pat_loo = 1:nr_pat,
        
    likelihood_events_loo = likelihood_events_tot;
    likelihood_events_loo(:, pat_loo, :) = [];
    likelihood_events_test = likelihood_events_tot(:, pat_loo, :);
    
    parm_struct_loo{pat_loo} = EBDPMCMCTest(likelihood_events_loo, ...
        parm_mcmc, likelihood_events_test);
    class_pat_loo(pat_loo) = parm_struct_loo{pat_loo}.class_pat_test;
    
    save progress_loo.txt pat_loo -ASCII
    save class_pat_loo.txt class_pat_loo -ASCII
    
end
    
parm_struct_tot = EBDPMCMCTest(likelihood_events_tot, parm_mcmc);
class_pat_tot = parm_struct_tot.class_pat;

figure(1), 
clf, hold on
for pid = unique(patient_id),
    
    I_pat = find(patient_id == pid);
    plot(patient_tp(I_pat), class_pat_tot(I_pat), 'k', 'LineWidth', 2)
    axis([min(patient_tp) max(patient_tp) 1 nr_events_lhrh]) 
        
end
hxlabel = xlabel('follow-up scan');
hylabel = ylabel('disease stage')'
set([hxlabel hylabel], ...
    'FontName', 'Helvetica', ...
    'FontSize', 12)
grid on
set(gcf, 'PaperPositionMode', 'auto');

%% Looking at hemispheric correspondence...
flag_reg = 1;
nr_roi_hem = nr_roi_lhrh;

mean_order = zeros(nr_roi_lhrh, 2);
std_order = zeros(nr_roi_lhrh, 2);
mean_order(:, 1) = ...
    mean(parm_struct_hem{1}.event_order_mcmc, 2);
mean_order(:, 2) = ...
    mean(parm_struct_hem{2}.event_order_mcmc, 2);
std_order(:, 1) = ...
    std(parm_struct_hem{1}.event_order_mcmc, [], 2);
std_order(:, 2) = ...
    std(parm_struct_hem{2}.event_order_mcmc, [], 2);

figure(6), clf, hold on
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
    'Color', [0 0 0])
set([hXlabel hYlabel], ...
    'FontName', 'Helvetica', ...
    'FontSize', 16)
set(gcf, 'PaperPositionMode', 'auto');
