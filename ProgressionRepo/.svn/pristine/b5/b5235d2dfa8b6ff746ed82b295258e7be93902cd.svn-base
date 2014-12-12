%% Setting up simulated data
clear all
close all

cnt_fig = 1;
version_mcmc = 2;
flag_mixt = 0;
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
X = zeros(nr_events, nr_pat, 2);  
for flag_reg = 1:2,
    
    [X(1:nr_roi, :, flag_reg)] = ...
        compute_p_A_D(atrophy_control{flag_reg}, ...
        atrophy_patient{flag_reg}, flag_mixt, flag_filt);
    
end
for pat = 1:nr_pat,
    
    X((nr_roi+1):(nr_roi + data_struct.cat_patient(pat)), pat, :) = 1;
    
end

%% Performing mcmc
nr_it_mcmc = 1e4;
nr_it_burnin = 1e4;
nr_it_hillclimb = 2e3;
thinning = 1e0;
nr_hillclimb = 5;

flag_reg = 1;

roi_select = [1 9 71 72 73];
X_select = X(roi_select, :, flag_reg);
thr_sig = 0.95;
X_select(X_select >= thr_sig) = 1;
X_select(X_select < thr_sig) = 0;

[parm_struct, diag_struct] = ...
    EventOnsetModel(X_select, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb);

subplot(4, 2, 1), imagesc(X_select)
subplot(4, 2, 2), plot(diag_struct.logp_hillclimb)
subplot(4, 2, 3), imagesc(parm_struct.model_progress_mean)
subplot(4, 2, 4), plot(diag_struct.logp_mcmc)
subplot(4, 2, 5), imagesc(parm_struct.model_progress_max)
subplot(4, 2, 6), imagesc(parm_struct.mu)
subplot(4, 2, 8), imagesc(parm_struct.sigma)

position_roi_mean = mean(parm_struct.mu, 2);
[d, order_roi] = sort(position_roi_mean, 'ascend');

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
figure(2), clf
imagesc(parm_struct.model_progress_mean(order_roi, :))
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


