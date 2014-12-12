%% Preparing data and calculating probability of atrophy in each region in
%% each patient
clear all
close all

load data_AD
atrophy_control = data_struct.data_control{1}.data_medianneg;
atrophy_patient = data_struct.data_patient{1}.data_medianneg;
diff_t_control = data_struct.diff_t_control;
diff_t_patient = data_struct.diff_t_patient;
control_id = data_struct.control_id;
control_tp = data_struct.control_tp;
patient_id = data_struct.patient_id;
patient_tp = data_struct.patient_tp;

[nr_roi, nr_con] = size(atrophy_control);
nr_pat = size(atrophy_patient, 2);

atrophy_sub_control = zeros(nr_roi, nr_con);
for id = unique(control_id),
    
    I_s = find(control_id == id);
    atrophy_sub_control(:, I_s(1)) = atrophy_control(:, I_s(1))/...
        diff_t_control(I_s(1));
    if max(control_tp(I_s)) > 2,
        
        for t = 3:max(control_tp(I_s)),
            
            I_t1 = find(control_tp(I_s) == t);
            I_t0 = find(control_tp(I_s) == t-1);
            atrophy_sub_control(:, I_s(I_t1)) = ...
                (atrophy_control(:, I_s(I_t1))-atrophy_control(:, I_s(I_t0)))/...
                (diff_t_control(I_s(I_t1)) - diff_t_control(I_s(I_t0)));
            
        end
        
    end
    
end
atrophy_sub_patient = zeros(nr_roi, nr_pat);
for id = unique(patient_id),
    
    I_s = find(patient_id == id);
    atrophy_sub_patient(:, I_s(1)) = atrophy_patient(:, I_s(1))/...
        diff_t_patient(I_s(1));
    if max(patient_tp(I_s)) > 2,
        
        for t = 3:max(patient_tp(I_s)),
            
            I_t1 = find(patient_tp(I_s) == t);
            I_t0 = find(patient_tp(I_s) == t-1);
            atrophy_sub_patient(:, I_s(I_t1)) = ...
                (atrophy_patient(:, I_s(I_t1))-atrophy_patient(:, I_s(I_t0)))/...
                (diff_t_patient(I_s(I_t1)) - diff_t_patient(I_s(I_t0)));
            
        end
        
    end
    
end

p_A_D = compute_p_A_D_aar(atrophy_patient, atrophy_control, atrophy_sub_patient, ...
    atrophy_sub_control, diff_t_patient, diff_t_control);

for id = unique(patient_id),
    
    I_s = find(patient_id == id);
    if max(patient_tp(I_s)) > 2,
        
        for t = 3:max(patient_tp(I_s)),
            
            I_t1 = find(patient_tp(I_s) == t);
            I_t0 = find(patient_tp(I_s) == t-1);
            p_A_D(:, I_s(I_t1)) = ...
                max([p_A_D(:, I_s(I_t1)) p_A_D(:, I_s(I_t0))], [], 2);
            
        end
        
    end        
    
end
p_A_D(p_A_D == 0) = eps;
p_A_D(p_A_D == 1) = 1-eps;

%% Performing MCMC

nr_it_mcmc = 1e5;
nr_it_burnin = 1e5;
nr_it_hillclimb = 2e3;
thinning = 1e0;
nr_hillclimb = 20;
version_mcmc = 2;
[parm_struct, diag_struct] = ...
    AtrophyModelMCMC2(p_A_D, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_mcmc);
hist2_mat = zeros(nr_roi, nr_roi);
    
[ord inds_av] = sort(mean(parm_struct.order_events, 2));
for roi1 = 1:nr_roi,
    
    for roi2 = 1:nr_roi,
        
        hist2_mat(roi1,roi2) = ...
            sum(parm_struct.order_events(inds_av(roi1), :) == roi2);
        
    end
    
end
%% Visualize results

figure(1), clf
subplot(121), plot(diag_struct.logp_hillclimb), 
xlabel('it'), ylabel('logp'), title('hill-climb')
subplot(122), plot(diag_struct.logp_order_events), 
xlabel('it'), ylabel('logp'), title('MCMC chain')

position_roi_mean = mean(parm_struct.order_events, 2);
position_roi_std = std(parm_struct.order_events, [], 2);
[d, order_roi] = sort(position_roi_mean, 'ascend');
name_roi = data_struct.name_roi;

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

% Making a annotated version of the histogram
max_length_str = 0;
for roi = 1:nr_roi,
    
    max_length_str = max([max_length_str length(name_roi{roi})]);
    
end
mat_labels = zeros(nr_roi, max_length_str);
for roi = 1:nr_roi,
    
    length_str = length(name_roi{order_roi(roi)});
    mat_labels(roi, 1:length_str) = name_roi{order_roi(roi)};
    
end
figure(2), clf
imagesc(hist2_mat)
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

    

