function robustnessTest2ImagingSubGroups(div_size,version_mixt,n_repeat,data_struct,cGroup)

%This function tests the robustness of the algoirthmn to the outliers in
%the control cohort. i.e. how consistent is the ordering of the events
%given that there may be some outliners in control i.e. patients classified
%as control

%To do this I will 
% - randomely select a set of controls (equal to 1/3 of the control size)
% - remove this from the control group
% - run the experiment
% - compare the ordering
% - repeat the expermint several times

cnt_fig = 1;
roi_id=data_struct.roi_id/256;
roi_label=data_struct.roi_id;
roi_name=data_struct.roi_name;
[nr_roi nr_pat] = size(data_struct.data_patient.avg_vol);
nr_con = size(data_struct.data_control.avg_vol,2);

nr_roi_lhrh = floor(nr_roi/2);
nr_events = nr_roi;
nr_events_lhrh = nr_roi_lhrh;
I_lh = [2:2:83];
I_rh = [1:2:17 21:2:83];

% Inilialize volumes for patient and control and for the combined left and right side
avg_vol_patient=data_struct.data_patient.avg_vol;
avg_vol_control=data_struct.data_control.avg_vol;

avg_vol_patient_lhrh=zeros(nr_roi_lhrh,nr_pat);
avg_vol_control_lhrh=zeros(nr_roi_lhrh,nr_con);

%for each region average over the left and right side
for roi = 1:nr_roi_lhrh,
    
    for pat = 1:nr_pat,
        
        mean_local = [data_struct.data_patient.avg_vol(I_lh(roi), pat) ...
            data_struct.data_patient.avg_vol(I_rh(roi), pat)];
        avg_vol_patient_lhrh(roi, pat) = mean(mean_local);
        
    end
    
    
    for con = 1:nr_con,
        
        mean_local = [data_struct.data_control.avg_vol(I_lh(roi), con) ...
            data_struct.data_control.avg_vol(I_rh(roi), con)];
        avg_vol_control_lhrh(roi, con) = mean(mean_local);
        
    end
    
end
%% randomely select some patinets to be part of control
%div_size=0.3333; 
rand_size=round(nr_con*div_size);
stored_rand_perm=zeros(n_repeat,nr_con);

if exist(['Mat Data/test2_stored_cGroup_' num2str(cGroup) 'rand_perm' num2str(n_repeat) 'Repeat' num2str(div_size) 'Portion.mat'])
    load(['Mat Data/test2_stored_cGroup_' num2str(cGroup) 'rand_perm' num2str(n_repeat) 'Repeat' num2str(div_size) 'Portion.mat'],'stored_rand_perm');
else
    
    for i=1:n_repeat
        stored_rand_perm(i,:)=randperm(nr_con);
    end
    save(['Mat Data/test2_stored_cGroup_' num2str(cGroup) 'rand_perm' num2str(n_repeat) 'Repeat' num2str(div_size) 'Portion.mat'],'stored_rand_perm');
end


%% ====================== fit the model ======================

%Now fitting the EBDP model to the data repeatedly
%We only need two functions for this: EBDPComputeLikelihood.m and EBDPMCMC.m

for cRep=1:n_repeat
    cRep
this_rand_select_indx= stored_rand_perm(cRep,1:rand_size);
this_control_indx=setdiff([1:nr_con],this_rand_select_indx);

this_data_controls=avg_vol_control_lhrh(:,this_control_indx);
this_data_patients=avg_vol_patient_lhrh;

%% fit the model
 opts = foptions;
 opts(14) = 1e3;
 
% likelihood_events = zeros(nr_events, nr_pat, 2);
% % likelihood_events(:, :, 1) is the likelihood given that no
% % event has occurred, likelihood_events(:, :, 2) is the
% % likelihood given that an event has occurred
% % imagesc(likelihood_events(:, :,
% % 2)./sum(likelihood_events, 3) gives something equivalent
% % to a posterior probability that an event has occurred
% 
% [likelihood_events, gmix_roi] = EBDPComputeLikelihood(avg_vol_patient', ...
%     avg_vol_control', version_mixt);
% 
 % The variables with _lhrh are likelihood values when left and right
 % hemisphere atrophy values are averaged. This is mainly done to reduce
 % the number of stages in the model and make the results somewhat more
 % presentable
 likelihood_events_lhrh = zeros(nr_events_lhrh, nr_pat, 2);
 [likelihood_events_lhrh, gmix_roi] = EBDPComputeLikelihood(this_data_patients', ...
     this_data_controls', version_mixt);
 
 
 idx_clinevent = []; % This variable is outdated
 
 % This was apparently needed to make things work properly. Some of the
 % likelihood values are zero, which the MCMC algorithm doesn't respond to
 % particularly well...
 %likelihood_events(likelihood_events < 1e-3) = 1e-3;
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

    
% parm_struct = EBDPMCMCTest(likelihood_events, ...
%     parm_mcmc);
 parm_struct_lhrh = EBDPMCMCTest(likelihood_events_lhrh, ...
     parm_mcmc);
% 
% parm_struct_hem{1} = EBDPMCMCTest(likelihood_events(I_lh, :, :), ...
%     parm_mcmc);
% parm_struct_hem{2} = EBDPMCMCTest(likelihood_events(I_rh, :, :), ...
%     parm_mcmc);

%save the parameters for this model
save(['Mat Data/test2_results_cGroup_' num2str(cGroup) 'versionLike' num2str(version_likelihood) 'Rep' num2str(cRep) 'Repeat' num2str(div_size) 'Portion.mat'],'this_data_patients','this_data_controls','likelihood', 'gmix_struct','parm_struct','this_rand_select_indx','this_control_indx');

end

% % %% ====================== Compare the results ======================
% mcmcSize=100000;
% event_order_mcmc_all=zeros(nr_var,mcmcSize*n_repeat);
% position_roi_mean_all=zeros(nr_var,n_repeat);
% for cRep=1:n_repeat
%     cRep
%     close all
%     
%     load(['Mat Data/test2_results_versionLike' num2str(version_likelihood) 'Rep' num2str(cRep) 'Repeat' num2str(div_size) 'Portion.mat'],'this_data_patients','this_data_controls','likelihood', 'gmix_struct','parm_struct','this_rand_select_indx','this_control_indx');
%     
%     
%     %% visualizing the event ordering
%     % Computing positional variance
%     nr_events = size(parm_struct.event_order_mcmc, 1);
%     hist2_mat = zeros([nr_events nr_events]);
%     [ord inds_av] = sort(mean(parm_struct.event_order_mcmc, 2));
%     for roi1 = 1:nr_events,
%         
%         for roi2 = 1:nr_events,
%             
%             hist2_mat(roi1, roi2) = ...
%                 sum(parm_struct.event_order_mcmc(inds_av(roi1), :) == roi2);
%             
%         end
%         
%     end
%     
%     % Visualizing positional variance...
%     position_roi_mean = mean(parm_struct.event_order_mcmc, 2);
%     position_roi_std = std(parm_struct.event_order_mcmc, [], 2);
%     
%     [d, order_roi] = sort(position_roi_mean, 'ascend');
%     % Making a annotated version of the histogram
%     max_length_str = 0;
%     for roi = 1:nr_events,
%         
%         max_length_str = max([max_length_str length(name_variables{roi})]);
%         
%     end
%     mat_labels = zeros(nr_events, max_length_str);
%     for roi = 1:nr_events,
%         
%         length_str = length(name_variables{order_roi(roi)});
%         mat_labels(roi, 1:length_str) = name_variables{order_roi(roi)};
%         
%     end
%     figure(1);
%     imagesc(hist2_mat)
%     map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
%     colormap(map_inverted)
%     axis square
%     set(gca, ...
%         'YTick', [1:nr_events], 'YTickLabel', char(mat_labels), ...
%         'YGrid', 'on')
%     hXLabel = xlabel('model stage');
%     hYLabel = ylabel('region');
%     set([hXLabel hYLabel], ...
%         'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
%     set(gcf, 'PaperPositionMode', 'auto');
%     set(gcf,'Color',[1 1 1])
% %     %%  ====================== visualizing the patient control distribution for this random selection ======================
% %     figure(2);
% %     for i_var = 1:nr_var,
% %         
% %         subplot(3, 4, i_var), hold on
% %         [hist_c, x_c] = ksdensity(this_data_controls(:, i_var));
% %         plot(x_c, hist_c, 'g');
% %         [hist_p, x_p] = ksdensity(this_data_patients(:, i_var));
% %         plot(x_p, hist_p, 'r');
% %         title(name_variables{i_var})     
% %         
% %     end
% %     set(gcf,'Color',[1 1 1])
% %     
%      event_order_mcmc_all(:,(cRep-1)*mcmcSize+1:cRep*mcmcSize)=parm_struct.event_order_mcmc;
%      position_roi_mean_all(:,cRep)=position_roi_mean;
% %     set(gcf,'Color',[1 1 1])
% %     pause
%     
% end
%  
% % %% ====================== Draw all of the results in one plot ======================
% % figure;
% % for cRep=1:n_repeat    
% % load(['Mat Data/test2_results_versionLike' num2str(version_likelihood) 'Rep' num2str(cRep) 'Repeat' num2str(div_size) 'Portion.mat'],'this_data_patients','this_data_controls','likelihood', 'gmix_struct','parm_struct','this_rand_select_indx','this_control_indx');
% %     
% %     nr_events = size(parm_struct.event_order_mcmc, 1);
% %     hist2_mat = zeros([nr_events nr_events]);
% %     [ord inds_av] = sort(mean(parm_struct.event_order_mcmc, 2));
% %     for roi1 = 1:nr_events,
% %         
% %         for roi2 = 1:nr_events,
% %             
% %             hist2_mat(roi1, roi2) = ...
% %                 sum(parm_struct.event_order_mcmc(inds_av(roi1), :) == roi2);
% %             
% %         end
% %         
% %     end
% %     
% %     % Visualizing positional variance...
% %     position_roi_mean = mean(parm_struct.event_order_mcmc, 2);
% %     position_roi_std = std(parm_struct.event_order_mcmc, [], 2);
% %     
% %     [d, order_roi] = sort(position_roi_mean, 'ascend');
% %     % Making a annotated version of the histogram
% %     max_length_str = 0;
% %     for roi = 1:nr_events,
% %         
% %         max_length_str = max([max_length_str length(name_variables{roi})]);
% %         
% %     end
% %     mat_labels = zeros(nr_events, max_length_str);
% %     for roi = 1:nr_events,
% %         
% %         length_str = length(name_variables{order_roi(roi)});
% %         mat_labels(roi, 1:length_str) = name_variables{order_roi(roi)};
% %         
% %     end
% %     subplot(2,ceil(n_repeat/2),cRep)
% %     imagesc(hist2_mat)
% %     hold on;
% %     map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
% %     colormap(map_inverted)
% %     axis square
% %     set(gca, ...
% %         'YTick', [1:nr_events], 'YTickLabel', char(mat_labels), ...
% %         'YGrid', 'on')
% %     hXLabel = xlabel('model stage');
% %     hYLabel = ylabel('region');
% %     set([hXLabel hYLabel], ...
% %         'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
% %     set(gcf, 'PaperPositionMode', 'auto');
% %     set(gcf,'Color',[1 1 1])
% % end
% 
% 
% %% ======================Draw the resulting event ordering as an average position of the results of n repetiotions ======================
% 
% % visualizing the event ordering
% %Computing positional variance
% nr_events = size(event_order_mcmc_all, 1);
% hist2_mat = zeros([nr_events nr_events]);
% [ord inds_av] = sort(mean(event_order_mcmc_all, 2));
% for roi1 = 1:nr_events,
%     
%     for roi2 = 1:nr_events,
%         
%         hist2_mat(roi1, roi2) = ...
%             sum(event_order_mcmc_all(inds_av(roi1), :) == roi2);
%         
%     end
%     
% end
% 
% % Visualizing positional variance...
% position_roi_mean = mean(event_order_mcmc_all, 2);
% position_roi_std = std(event_order_mcmc_all, [], 2);
% 
% [d, order_roi] = sort(position_roi_mean, 'ascend');
% % Making a annotated version of the histogram
% max_length_str = 0;
% for roi = 1:nr_events,
%     
%     max_length_str = max([max_length_str length(name_variables{roi})]);
%     
% end
% mat_labels = zeros(nr_events, max_length_str);
% for roi = 1:nr_events,
%     
%     length_str = length(name_variables{order_roi(roi)});
%     mat_labels(roi, 1:length_str) = name_variables{order_roi(roi)};
%     
% end
% figure(4);
% imagesc(hist2_mat)
% map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
% colormap(map_inverted)
% axis square
% set(gca, ...
%     'YTick', [1:nr_events], 'YTickLabel', char(mat_labels), ...
%     'YGrid', 'on')
% hXLabel = xlabel('model stage');
% hYLabel = ylabel('region');
% set([hXLabel hYLabel], ...
%     'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
% set(gcf, 'PaperPositionMode', 'auto');
% set(gcf,'Color',[1 1 1])
% 
% 
