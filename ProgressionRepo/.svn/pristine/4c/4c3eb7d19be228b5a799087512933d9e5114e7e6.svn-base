%% This function runs the model on a subset of regins for fast training and testing

%% Loading data
close all

cnt_fig = 1;
file_data = '/Users/jaghajan/Documents/Experiments/Mat Data/MAPER_DATA.mat'; % This needs to changed...
%file_data = '/Users/jaghajan/Documents/Experiments/Mat Data/MAPER_ICV_NORM_DATA.mat'; % This needs to changed...

load(file_data);

roi_select=[1,2,3,4,5,6,9,10,34,35,50,51];

roi_id=data_struct.roi_id(roi_select)/256;
roi_label=data_struct.roi_id(roi_select);
roi_name=data_struct.roi_name(roi_select);
nr_roi=length(roi_select);
nr_pat = size(data_struct.data_patient.avg_vol,2);
nr_con = size(data_struct.data_control.avg_vol,2);

nr_roi_lhrh = floor(nr_roi/2);
nr_events = nr_roi;
nr_events_lhrh = nr_roi_lhrh;
I_lh = [2:2:nr_roi];
I_rh = [1:2:nr_roi];

%% Inilialize volumes for patient and control and for the combined left and right side
avg_vol_patient=data_struct.data_patient.avg_vol(roi_select,:);
avg_vol_control=data_struct.data_control.avg_vol(roi_select,:);

avg_vol_patient_lhrh=zeros(nr_roi_lhrh,nr_pat);
avg_vol_control_lhrh=zeros(nr_roi_lhrh,nr_con);

%for each region average over the left and right side
for roi = 1:nr_roi_lhrh,
    
    for pat = 1:nr_pat,
        
        mean_local = [data_struct.data_patient.avg_vol(roi_select(I_lh(roi)), pat) ...
            data_struct.data_patient.avg_vol(roi_select(I_rh(roi)), pat)];
        avg_vol_patient_lhrh(roi, pat) = mean(mean_local);
        
    end
    
    
    for con = 1:nr_con,
        
        mean_local = [data_struct.data_control.avg_vol(roi_select(I_lh(roi)), con) ...
            data_struct.data_control.avg_vol(roi_select(I_rh(roi)), con)];
        avg_vol_control_lhrh(roi, con) = mean(mean_local);
        
    end
    
end

%% Compute data Likelihood given event/no event
version_mixt = 7; % Choosing the mixture of the Gaussian and uniform distributions
opts = foptions;
opts(14) = 1e3;

likelihood_events = zeros(nr_events, nr_pat, 2);
% likelihood_events(:, :, 1) is the likelihood given that no
% event has occurred, likelihood_events(:, :, 2) is the
% likelihood given that an event has occurred
% imagesc(likelihood_events(:, :,
% 2)./sum(likelihood_events, 3) gives something equivalent
% to a posterior probability that an event has occurred

[likelihood_events, gmix_roi] = EBDPComputeLikelihood(avg_vol_patient', ...
    avg_vol_control', version_mixt);

% The variables with _lhrh are likelihood values when left and right
% hemisphere atrophy values are averaged. This is mainly done to reduce
% the number of stages in the model and make the results somewhat more
% presentable
likelihood_events_lhrh = zeros(nr_events_lhrh, nr_pat, 2);
[likelihood_events_lhrh, gmix_roi] = EBDPComputeLikelihood(avg_vol_patient_lhrh', ...
    avg_vol_control_lhrh', version_mixt);


idx_clinevent = []; % This variable is outdated

% This was apparently needed to make things work properly. Some of the
% likelihood values are zero, which the MCMC algorithm doesn't respond to
% particularly well...
likelihood_events(likelihood_events < 1e-3) = 1e-3;
likelihood_events_lhrh(likelihood_events_lhrh < 1e-3) = 1e-3;

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

    
parm_struct = EBDPMCMCTest(likelihood_events, ...
    parm_mcmc);
parm_struct_lhrh = EBDPMCMCTest(likelihood_events_lhrh, ...
    parm_mcmc);

parm_struct_hem{1} = EBDPMCMCTest(likelihood_events(I_lh, :, :), ...
    parm_mcmc);
parm_struct_hem{2} = EBDPMCMCTest(likelihood_events(I_rh, :, :), ...
    parm_mcmc);


%% Visualizing results
% Some diagnostics first
flag_reg = 1;

close all
cnt_fig = 1;
figure(1), clf
subplot(221), plot(parm_struct_lhrh.log_likelihood_gradient_ascent)
subplot(222), plot(parm_struct_lhrh.log_likelihood_mcmc)
subplot(223), imagesc(parm_struct_lhrh.event_order_mcmc);
subplot(224), plot(parm_struct_lhrh.p_false_mcmc')
cnt_fig = cnt_fig + 1;

% Making histogram
% Making histogram
for roi = 1:(nr_roi_lhrh),
    
    labels_lhrh{roi} = sprintf('%02i', roi);
    
end

for roi = 1:(nr_roi_lhrh),
    
    roi_name_lhrh{roi} = roi_name{2*roi}(1:end-5);
    
end

% Making goose plot
position_roi_lhrh_mean = mean(parm_struct_lhrh.event_order_mcmc, 2);
position_roi_lhrh_std = std(parm_struct_lhrh.event_order_mcmc, [], 2);

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
[ord inds_av] = sort(mean(parm_struct_lhrh.event_order_mcmc, 2));
for roi1 = 1:nr_events_lhrh,
    
    for roi2 = 1:nr_events_lhrh,
        
        hist2_mat_lhrh(roi1, roi2) = ...
            sum(parm_struct_lhrh.event_order_mcmc(inds_av(roi1), :) == roi2);
        
    end
    
end

[d, order_roi_lhrh] = sort(position_roi_lhrh_mean, 'ascend');
% Making a annotated version of the histogram
max_length_str = 0;
for roi = 1:nr_events_lhrh,
    
    max_length_str = max([max_length_str length(roi_name_lhrh{roi})]);
    
end
mat_labels = zeros(nr_events_lhrh, max_length_str);
for roi = 1:nr_events_lhrh,
    
    length_str = length(roi_name_lhrh{order_roi_lhrh(roi)});
    mat_labels(roi, 1:length_str) = roi_name_lhrh{order_roi_lhrh(roi)};
    
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
