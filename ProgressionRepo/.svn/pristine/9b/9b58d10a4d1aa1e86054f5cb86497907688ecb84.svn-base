function dataSimulationAllPermutationsDemo
clear
close all

%% Initialization
nr_events=3;        % select the nuber of events
nr_patients=300;    % select the nuber of patients
nr_controls=300;     % select the nuber of controls
sim_version=1;

% rand('seed',14);
% randn('seed',14);

%% Define the idstribututions
switch sim_version
    
    case 1 %% case 1: Simple case where control and patient distributions are the same for all biomarkers
        %Set the mean and sd of all biomarker to be the same
        for roi=1:nr_events
            gmix_struct_est(roi).mean_control=0;
            gmix_struct_est(roi).sd_control=1;
            gmix_struct_est(roi).mean_patient=-1;
            gmix_struct_est(roi).sd_patient=1;
        end
    case 2 %% case 2: Biomarkers distributions have different means but same variance
        gmix_struct_est(1).mean_control=5;
        gmix_struct_est(1).sd_control=2.5;
        gmix_struct_est(1).mean_patient=-5;
        gmix_struct_est(1).sd_patient=2.5;
        
        for roi=2:nr_events
            gmix_struct_est(roi).mean_control=gmix_struct_est(1).mean_control+randn*2;
            gmix_struct_est(roi).sd_control=2.5;
            gmix_struct_est(roi).mean_patient=gmix_struct_est(1).mean_patient+randn*2;
            gmix_struct_est(roi).sd_patient=2.5;
        end
        
    case 3 %% case 2: Biomarkers distributions the same means but higher variance for patient dist
        for roi=1:nr_events
            gmix_struct_est(roi).mean_control=5;
            gmix_struct_est(roi).sd_control=2.5;
            gmix_struct_est(roi).mean_patient=-5;
            gmix_struct_est(roi).sd_patient=5.5;
        end
        
    case 4 %% case 4: Biomarkers distributions have different means and higher variance for patient dist (same for all roi)
        gmix_struct_est(1).mean_control=5;
        gmix_struct_est(1).sd_control=2.5;
        gmix_struct_est(1).mean_patient=-5;
        gmix_struct_est(1).sd_patient=5.5;
        
        for roi=2:nr_events
            rand_dist=abs(randn*5);
            gmix_struct_est(roi).mean_control=gmix_struct_est(1).mean_control+rand_dist;
            gmix_struct_est(roi).sd_control=gmix_struct_est(1).sd_control;
            gmix_struct_est(roi).mean_patient=gmix_struct_est(1).mean_patient+rand_dist;
            gmix_struct_est(roi).sd_patient=gmix_struct_est(1).sd_patient;
        end
        
    case 5 %% case 4: Biomarkers distributions have different means and higher variance for patient dist diff for each roi
        
        gmix_struct_est(1).mean_control=-67.5;
        gmix_struct_est(1).sd_control=27.68;
        gmix_struct_est(1).mean_patient=-136.33;
        gmix_struct_est(1).sd_patient=54.14;
        
        gmix_struct_est(2).mean_control=213.22;
        gmix_struct_est(2).sd_control=52.42;
        gmix_struct_est(2).mean_patient=130.45;
        gmix_struct_est(2).sd_patient=24.24;
        
        gmix_struct_est(3).mean_control=-0.0230;
        gmix_struct_est(3).sd_control=0.0116;
        gmix_struct_est(3).mean_patient=-0.0418;
        gmix_struct_est(3).sd_patient=0.0151;
        
        gmix_struct_est(4).mean_control= 0.0035;
        gmix_struct_est(4).sd_control=4.1151e-04;
        gmix_struct_est(4).mean_patient=0.0028;
        gmix_struct_est(4).sd_patient=5.3181e-04;
        
        gmix_struct_est(5).mean_control= -6.1745;
        gmix_struct_est(5).sd_control=6.1680;
        gmix_struct_est(5).mean_patient= -13.4318;
        gmix_struct_est(5).sd_patient=9.0717;
        
        
end





% assume a disease stage for each patient and allow equal number of
% patients in each stage
patient_stage=[];
for i=1:nr_events
    patient_stage=[patient_stage ones(1,nr_patients/nr_events)*i];
end
name_variables={'Biomarker1','Biomarker2','Biomarker3','Biomarker4','Biomarker5'}

%% Choose a random ordering for the biomarkers
true_order=perms([1:nr_events]);
n_perms=size(true_order,1);

for c_order=3:n_perms
    tic
    this_order=true_order(c_order,:);
    
    %% Data Generation
    if sim_version==0  % The very basic case with step like fixed values for has/hasnot occured
        has_occr_val=-5;
        has_not_occ_val=5;
        control_val=5;
        for k=1:nr_events
            indx=find(patient_stage==k);
            for roi=1:k
                data_patients(indx,this_order(roi))=has_occr_val;
            end
            
            for roi=k+1:nr_events
                data_patients(indx,this_order(roi))=has_not_occ_val;
            end
        end
        
        %generate controls
        for roi=1:nr_events
            data_controls(:,this_order(roi))=ones(nr_controls,1)*control_val;
        end
    else
        %generate patients
        % for each stage
        for k=1:nr_events
            indx=find(patient_stage==k);
            for roi=1:k
                data_patients(indx,this_order(roi))=gmix_struct_est(this_order(roi)).mean_patient+randn(1,length(indx))*gmix_struct_est(this_order(roi)).sd_patient;
            end
            
            for roi=k+1:nr_events
                data_patients(indx,this_order(roi))=gmix_struct_est(this_order(roi)).mean_control+randn(1,length(indx))*gmix_struct_est(this_order(roi)).sd_control;
            end
        end
        
        %generate controls
        for roi=1:nr_events
            data_controls(:,this_order(roi))=gmix_struct_est(this_order(roi)).mean_control+randn(1,nr_controls)*gmix_struct_est(this_order(roi)).sd_control;
        end
    end
    %% Data Visualization
%     
%     figure;
%     for i_var = 1:nr_events
%         
%         subplot(3, 4, i_var), hold on
%         [hist_c, x_c] = ksdensity(data_controls(:, i_var));
%         plot(x_c, hist_c, 'g');
%         [hist_p, x_p] = ksdensity(data_patients(:, i_var));
%         plot(x_p, hist_p, 'r');
%         title(['Biomarker ' num2str(i_var)])
%         set(gcf,'Color',[1 1 1])
%     end
    
    %% fit the model
    
%     I_select=[1:5];
%     data_controls=data_controls(:,I_select);
%     data_patients=data_patients(:,I_select);
%     
    version_likelihood = 13;  %version_likelihood==10 is the same as 4 but with threshold on likelihood ratio
    nrAttempts=5;
    threshold_flag=0;

    [likelihood, gmix_struct] = ...
        EBDPComputeLikelihood(data_patients, data_controls, version_likelihood,threshold_flag,nrAttempts);
    
    %% Since the number of events is small we can get the closed form solution
    est_order(c_order,:)=findNEventOrdering(likelihood);
    
    
    %% Set the likelihood of event occured for all vallues>control.mean to zero
    % for roi=1:5
    % indx=find(data_patients(:,roi)>gmix_struct.gmix_controls{roi}.centres);
    % like2=likelihood(roi,:,2);
    % like2(indx)=0;
    % likelihood(roi,:,2)=like2;
    % end
    
    % % Checking how the mixture distributions fit the real distributions
    % % Waling through one variable at a time
    % % The left plot is the true data distribution, the right plot is the
    % % mixture likelihood
    % for i_var = 1:nr_var,
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
end

%% calculate percentage of ecorrect estimation of ordering
for i=1:n_perms
diff(i)=sum(abs(true_order(i,:)-est_order(i,:)));
end
pct_correct=length(find(diff==0))/n_perms*100

i