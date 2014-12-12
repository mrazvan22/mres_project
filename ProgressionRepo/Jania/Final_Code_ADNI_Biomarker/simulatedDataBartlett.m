function simulatedDataBartlett

%% Initialization
%close all

%read the data
file_name='BartlettData/simulated_data_2011_09_27testData1JWB.xls';
[xls_data name_variables]=xlsread(file_name);

name_variables={name_variables{1:end-1}};

data=xls_data(:,1:end-1);
id_control=xls_data(:,end);
if strcmp(file_name,'BartlettData/Simulation_data_for_Danny.xls')
    id_control(id_control==0)=100;
    id_control(id_control==1)=0;
    id_control(id_control==100)=1;
end

data_controls=data(id_control==1,:);
data_patients=data(id_control==0,:);

[nr_patients,nr_events] = size(data_patients)
%% visualize the data

% Checking the distributions of all variables that remain in
% controls(green) and patients (red)

figure(1), clf
for i_var = 1:nr_events,
    
    subplot(3, 4, i_var), hold on
    [hist_c, x_c] = ksdensity(data_controls(:, i_var));
    plot(x_c, hist_c, 'g');
    [hist_p, x_p] = ksdensity(data_patients(:, i_var));
    plot(x_p, hist_p, 'r');
    title(name_variables{i_var})
    
end

% %noticed some rois need swapping to have patient values lower than control
% I_swap=[1,5:9];
% data_controls(:,I_swap)=-data_controls(:,I_swap);
% data_patients(:,I_swap)=-data_patients(:,I_swap);
% data(:,I_swap)=-data(:,I_swap);
% 
% figure(2), clf
% for i_var = 1:nr_events,
%     
%     subplot(3, 4, i_var), hold on
%     [hist_c, x_c] = ksdensity(data_controls(:, i_var));
%     plot(x_c, hist_c, 'g');
%     [hist_p, x_p] = ksdensity(data_patients(:, i_var));
%     plot(x_p, hist_p, 'r');
%     title(name_variables{i_var})
%     
% end

% 
% % draw L-plots
% figure
% pCount=0;
% for roi1=1:nr_events
%     for roi2=roi1+1:nr_events
%         pCount=pCount+1;
%         if mod(pCount,12)==1
%             figure;
%             pCount=1;
%         end
%         h1=subplot(3,4,pCount);plot(data_patients(:,roi1),data_patients(:,roi2),'rx');
%         hold on; plot(data_controls(:,roi1),data_controls(:,roi2),'gx');
%         axis([h1],'square')
%         set(gcf,'Color',[ 1 1 1]); 
%         proceding_region=findWhichLPlotQuarter(data_patients(:,roi1),mean(data_controls(:,roi1)),data_patients(:,roi2),mean(data_controls(:,roi2)));
%         title(proceding_region)
%         xlabel(name_variables{roi1});
%         ylabel(name_variables{roi2});
%     end
% end

%% fit EBP model
version_likelihood = 13;
threshold_flag=0;
nAttempts=5;
[likelihood, gmix_struct] = ...
    EBDPComputeLikelihood(data_patients, data_controls, version_likelihood,threshold_flag,nAttempts);


% for i_var = 1:nr_events,
%     
%     data_tot = cat(1, data_controls(:, i_var), data_patients(:, i_var));
%     figure;
%     p_controls = gmmprob(gmix_struct.gmix_controls{i_var}, data_controls(:, i_var));
%     p_patients = gmmprob(gmix_struct.gmix_patients{i_var}, data_patients(:, i_var));
%     %         p_controls = gmix_struct.gmix_roi{i_var}.prob_N;
%     %         p_patients = gmix_struct.gmix_roi{i_var}.prob_U;
%     [hist_c, x_c] = ksdensity(data_controls(:, i_var));
%     [hist_p, x_p] = ksdensity(data_patients(:, i_var));
%     
%     subplot(121), hold on
%     plot(x_c, hist_c, 'g');
%     plot(x_p, hist_p, 'r');
%     
%     subplot(122), hold on
%     plot(data_controls(:, i_var), p_controls, 'g.')
%     plot(data_patients(:, i_var), p_patients, 'r.');
%     %plot(data_controls(:, i_var), p_controls(id_controls_csf == 1), 'g.')
%     %plot(data_patients(:, i_var), p_patients(id_controls_csf == 0), 'r.');
%     set(gcf,'Color',[1 1 1])
%     
%     fprintf('%s\n', name_variables{i_var})
%     pause
%     
% end

%findNEventOrdering(likelihood)

%% draw the likelihood ratios
%close all
for roi=1:5
subplot(1,5,roi);plot(data_patients(:,roi),log(likelihood(roi,:,2)./likelihood(roi,:,1)),'r.');
hold on
plot(data_patients(:,roi),log(likelihood(roi,:,1)./likelihood(roi,:,2)),'g.');
end

%% Performing MCMC

% Performing mcmc
parm_mcmc.nr_gradient_ascent = 5;
parm_mcmc.nr_it_gradient_ascent = 2e3;
parm_mcmc.nr_it_burnin = 1e5;
parm_mcmc.nr_it_mcmc = 1e5;
parm_mcmc.interval_display = 1e2;
parm_mcmc.flag_sum = 2; % This should always be set to 2.
parm_mcmc.nr_it_check_p_false = 1e2;
parm_mcmc.range_p_false = [0.05 0.1
    0.05 0.1];
parm_mcmc.std_p_false = [1e-4 1e-4]';
parm_mcmc.idx_clinevent = [];

parm_struct = EBDPMCMCTest(likelihood, ...
    parm_mcmc);

% Computing positional variance
nr_events = size(parm_struct.event_order_mcmc, 1);
hist2_mat = zeros([nr_events nr_events]);
[ord inds_av] = sort(mean(parm_struct.event_order_mcmc, 2));
for roi1 = 1:nr_events,

    for roi2 = 1:nr_events,

        hist2_mat(roi1, roi2) = ...
            sum(parm_struct.event_order_mcmc(inds_av(roi1), :) == roi2);

    end

end

% Visualizing positional variance...
position_roi_mean = mean(parm_struct.event_order_mcmc, 2);
position_roi_std = std(parm_struct.event_order_mcmc, [], 2);

[d, order_roi] = sort(position_roi_mean, 'ascend');
% Making a annotated version of the histogram
max_length_str = 0;
for roi = 1:nr_events,

   max_length_str = max([max_length_str length(name_variables{roi})]);

end
mat_labels = zeros(nr_events, max_length_str);
for roi = 1:nr_events,

    length_str = length(name_variables{order_roi(roi)});
    mat_labels(roi, 1:length_str) = name_variables{order_roi(roi)};
end
figure(3), clf,
imagesc(hist2_mat)
map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
colormap(map_inverted)
axis square
set(gca, ...
    'YTick', [1:nr_events], 'YTickLabel', char(mat_labels), ...
    'YGrid', 'on')
hXLabel = xlabel('model stage');
hYLabel = ylabel('region');
set([hXLabel hYLabel], ...
    'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'Color',[1 1 1]);
