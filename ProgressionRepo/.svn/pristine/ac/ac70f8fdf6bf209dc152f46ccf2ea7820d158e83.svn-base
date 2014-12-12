close all
clear all

% load the learned parameters
load('JaniaAllRegionParamsMoGNormICV.mat')
 %data_patients=data_struct.data_patient.avg_vol;
 %data_controls=data_struct.data_control.avg_vol;
%for lhrh regions combined
data_patients=avg_vol_patient_lhrh;
data_controls=avg_vol_control_lhrh;


% Checking how the mixture distributions fit the real distributions
% Waling through one variable at a time
% The left plot is the true data distribution, the right plot is the
% mixture likelihood

for i_var = 1:nr_events_lhrh,
    
    data_tot = cat(2, data_controls(i_var,:), data_patients(i_var,:))';
    figure;
    %p_controls = gmmprob(gmix_roi_all.gmix_controls{i_var}, data_controls(i_var,:)');
    %p_patients = gmmprob(gmix_roi_all.gmix_patients{i_var}, data_patients(i_var,:)');
    %for lhrh regions combined
    p_controls = gmmprob(gmix_roi.gmix_controls{i_var}, data_controls(i_var,:)');
    p_patients = gmmprob(gmix_roi.gmix_patients{i_var}, data_patients(i_var,:)');
    
    %         p_controls = gmix_struct.gmix_roi{i_var}.prob_N;
    %         p_patients = gmix_struct.gmix_roi{i_var}.prob_U;
    [hist_c, x_c] = ksdensity(data_controls(i_var,:));
    [hist_p, x_p] = ksdensity(data_patients(i_var,:));
    
    subplot(121), hold on
    plot(x_c, hist_c, 'g');
    plot(x_p, hist_p, 'r');
        axis([min([min(x_p),min(x_c),min(data_controls(i_var,:)), min(data_patients(i_var,:))]), ...
        max([max(x_p),max(x_c), max(data_controls(i_var,:)), max(data_patients(i_var,:))]), ...
        min([min(hist_p),min(hist_c), min(p_controls), min(p_patients)]), ...
        max([max(hist_p),max(hist_c), max(p_controls), max(p_patients)])])
    
    subplot(122), hold on
    plot(data_controls(i_var,:), p_controls, 'g.')
    plot(data_patients(i_var,:), p_patients, 'r.');
    %plot(data_controls(:, i_var), p_controls(id_controls_csf == 1), 'g.')
    %plot(data_patients(:, i_var), p_patients(id_controls_csf == 0), 'r.');
    axis([min([min(x_p),min(x_c),min(data_controls(i_var,:)), min(data_patients(i_var,:))]), ...
        max([max(x_p),max(x_c), max(data_controls(i_var,:)), max(data_patients(i_var,:))]), ...
        min([min(hist_p),min(hist_c), min(p_controls), min(p_patients)]), ...
        max([max(hist_p),max(hist_c), max(p_controls), max(p_patients)])])
    set(gcf,'Color',[1 1 1])
    
    fprintf('%s\n', roi_name_lhrh{i_var})
    pause
    
end
