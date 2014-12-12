%% Loading data
clear all
close all

files = spm_select([1 10], 'mat');
nr_analysis = size(files, 1);
for a = 1:nr_analysis,
    
    load(deblank(files(a, :)));
    results{a} = results_struct;
    
end
nr_roi = length(results{1}.data_struct.roi_id);
nr_pat = size(results{1}.data_struct.data_patient{1}.data_mixt, 2);
%% Looking at hemispheric symmetry in results

I_lh = [1 3:36];
I_rh = [2 37:70];

figure(1), clf
for a = 1:nr_analysis,
    
    order_events = results{a}.order_events;
    mean_order_lh = mean(order_events(I_lh, :), 2);
    std_order_lh = std(order_events(I_lh, :), [], 2);
    
    mean_order_rh = mean(order_events(I_rh, :), 2);
    std_order_rh = std(order_events(I_rh, :), [], 2);
    
    subplot(1, nr_analysis, a), hold on
    errorbar(mean_order_lh, mean_order_rh, std_order_rh, '.')
    herrorbar(mean_order_lh, mean_order_rh, std_order_lh, '.')
    
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
    



