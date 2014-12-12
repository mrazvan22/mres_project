function DataSimulationWithParamsDemo1k
clear
%close all

%% Initialization
nr_events=5;        % select the nuber of events
nr_patients=500;    % select the nuber of patients
nr_controls=500;     % select the nuber of controls
sim_version=1;
stp=1;              %stp is a number to multiply by the sd to make the samples more spread out
threshold_flag=0;
mixT_flag=0;
noise_flag=0;
nois_pct=[0:0.05:0.95];
   rand('seed',1);
   randn('seed',1);

for i=1:nr_events
    name_variables{i}=['Biomarker' num2str(i)];
end
% assume a disease stage for each patient and allow equal number of
%patients in each stage
patient_stage=[];
for i=1:nr_patients/nr_events
    patient_stage=[patient_stage ones(1,nr_patients/nr_events)*i];
end
%patient_stage=ceil(rand(1,nr_patients)*nr_events);

%% Define the distribututions
switch sim_version
    
    case 1 %% case 1: Simple case where control and patient distributions are the same for all biomarkers
        %Set the mean and sd of all biomarker to be the same
        for roi=1:nr_events
            gmix_struct(roi).mean_control=0;
            gmix_struct(roi).sd_control=1;
            gmix_struct(roi).mean_patient=-1;
            gmix_struct(roi).sd_patient=1;
        end
    case 2 %% case 2: Biomarkers distributions have different means but same variance
        gmix_struct(1).mean_control=5;
        gmix_struct(1).sd_control=2.5;
        gmix_struct(1).mean_patient=-5;
        gmix_struct(1).sd_patient=2.5;
        
        for roi=2:nr_events
            gmix_struct(roi).mean_control=gmix_struct(1).mean_control+randn*5;
            gmix_struct(roi).sd_control=2.5;
            gmix_struct(roi).mean_patient=gmix_struct(1).mean_patient+randn*5;
            gmix_struct(roi).sd_patient=2.5;
        end
        
    case 3 %% case 2: Biomarkers distributions the same means but higher variance for patient dist
        for roi=1:nr_events
            gmix_struct(roi).mean_control=5;
            gmix_struct(roi).sd_control=2.5;
            gmix_struct(roi).mean_patient=-5;
            gmix_struct(roi).sd_patient=5.5;
        end
        
    case 4 %% case 4: Biomarkers distributions have different means and higher variance for patient dist (same for all roi)
        gmix_struct(1).mean_control=6;
        gmix_struct(1).sd_control=2.5;
        gmix_struct(1).mean_patient=-6;
        gmix_struct(1).sd_patient=5.5;
        
        for roi=2:nr_events
            rand_num=randn*2;
            gmix_struct(roi).mean_control=gmix_struct(1).mean_control+rand_num;
            gmix_struct(roi).sd_control=gmix_struct(1).sd_control;
            gmix_struct(roi).mean_patient=gmix_struct(1).mean_patient+rand_num;
            gmix_struct(roi).sd_patient=gmix_struct(1).sd_patient;
        end
        
    case 5 %% case 4: Biomarkers distributions have different means and higher variance for patient dist diff for each roi
        
        gmix_struct(1).mean_control=-67.5;
        gmix_struct(1).sd_control=27.68;
        gmix_struct(1).mean_patient=-136.33;
        gmix_struct(1).sd_patient=54.14;
        
        gmix_struct(2).mean_control=213.22;
        gmix_struct(2).sd_control=52.42;
        gmix_struct(2).mean_patient=130.45;
        gmix_struct(2).sd_patient=24.24;
        
        gmix_struct(3).mean_control=-0.0230;
        gmix_struct(3).sd_control=0.0116;
        gmix_struct(3).mean_patient=-0.0418;
        gmix_struct(3).sd_patient=0.0151;
        
        gmix_struct(4).mean_control= 0.0035;
        gmix_struct(4).sd_control=4.1151e-04;
        gmix_struct(4).mean_patient=0.0028;
        gmix_struct(4).sd_patient=5.3181e-04;
        
        gmix_struct(5).mean_control= -6.1745;
        gmix_struct(5).sd_control=6.1680;
        gmix_struct(5).mean_patient= -13.4318;
        gmix_struct(5).sd_patient=9.0717;
        
        
        
        
        
    case 6 %% case generate from T-distributions, same mean and same dof
        %create data
        
        N_DIM = 1;
        pr_dist=1;
        
        for roi=1:nr_events
            %fill in actual mean, prior and covariance for controls
            tDist_controls(roi).mean = 5;
            tDist_controls(roi).prior = pr_dist;
            tDist_controls(roi).dof = 3;
            nDimNoise = randperm(N_DIM);
            F = randn(N_DIM,nDimNoise(1));
            %tDist_controls(roi).cov = F*F'+0.1*eye(N_DIM)+diag(0.2*rand(N_DIM,1));
            tDist_controls(roi).cov=2.5*2.5;
            
            %fill in actual mean, prior and covariance for patinets
            tDist_patients(roi).mean = -5;
            tDist_patients(roi).prior = pr_dist;
            tDist_patients(roi).dof = 3;
            nDimNoise = randperm(N_DIM);
            F = randn(N_DIM,nDimNoise(1));
            %tDist_patients(roi).cov = F*F'+0.1*eye(N_DIM)+diag(0.2*rand(N_DIM,1));
            tDist_patients(roi).cov=5*5;
        end
        
    case 7 %% case generate from T-distributions, same mean and diff dof
        %create data
        
        N_DIM = 1;
        prDist=1;
        
        for (roi = 1:nr_events)
            tDist_controls(roi).mean = 5+rand(N_DIM,1);
            tDist_controls(roi).prior = prDist;
            tDist_controls(roi).dof = rand(1)*2+1.0;
            nDimNoise = randperm(N_DIM);
            F = randn(N_DIM,nDimNoise(1));
            tDist_controls(roi).cov = F*F'+0.1*eye(N_DIM)+diag(0.2*rand(N_DIM,1));
            
            tDist_patients(roi).mean = -5+rand(N_DIM,1);
            tDist_patients(roi).prior = prDist;
            tDist_patients(roi).dof = rand(1)*5+1.0;
            nDimNoise = randperm(N_DIM);
            F = randn(N_DIM,nDimNoise(1));
            tDist_patients(roi).cov = F*F'+0.1*eye(N_DIM)+diag(0.2*rand(N_DIM,1));
        end
end


%% Data Generation
if noise_flag>0
    itr=length(nois_pct);
else
    itr=1;
end

for cItr=1:itr
    %close all
    if sim_version==0  % The very basic case with step like fixed values for has/hasnot occured
        has_occr_val=-5;
        has_not_occ_val=5;
        control_val=10;
        for k=1:nr_events
            indx=find(patient_stage==k);
            for roi=1:k
                data_patients(indx,roi)=has_occr_val;
            end
            
            for roi=k+1:nr_events
                data_patients(indx,roi)=has_not_occ_val;
            end
        end
        %generate controls
        for roi=1:nr_events
            data_controls(:,roi)=ones(nr_controls,1)*control_val;
        end
        
        
    elseif sim_version>=6 %Generate samples form a T-distribution
        % for each stage
        for k=1:nr_events
            indx=find(patient_stage==k);
            for roi=1:k
                data_patients(indx,roi)=mixTGenerate(tDist_patients(roi),length(indx));
            end
            
            for roi=k+1:nr_events
                data_patients(indx,roi)=mixTGenerate(tDist_controls(roi),length(indx));
            end
        end
        %generate controls
        for roi=1:nr_events
            data_controls(:,roi)=mixTGenerate(tDist_controls(roi),length(indx));
        end
        
        
    else
        %generate patients
        % for each stage
        for k=1:nr_events
            indx=find(patient_stage==k);
            for roi=1:k
                data_patients(indx,roi)=gmix_struct(roi).mean_patient+randn(1,length(indx))*gmix_struct(roi).sd_patient;
            end
            
            for roi=k+1:nr_events
                data_patients(indx,roi)=gmix_struct(roi).mean_control+randn(1,length(indx))*gmix_struct(roi).sd_control;
            end
        end
        %generate controls
        for roi=1:nr_events
            data_controls(:,roi)=gmix_struct(roi).mean_control+randn(1,nr_controls)*gmix_struct(roi).sd_control;
        end
        
        switch noise_flag
            case 0 % Do not add noise
            case 1 % Put some patients into the control distribution
                perm=randperm(nr_controls);
                noise_indx=perm(1:round(nois_pct(cItr)*nr_controls));
                for roi=1:nr_events
                    %                     perm=randperm(nr_controls);
                    %                     noise_indx=perm(1:round(nois_pct(cItr)*nr_controls));
                    data_controls(noise_indx,roi)=gmix_struct(roi).mean_patient+randn(length(noise_indx),1)*gmix_struct(roi).sd_patient;
                end
                
            case 2 % Add some patients with random uniform distribution
                for roi=1:nr_events
                    perm=randperm(nr_patients);
                    noise_indx=perm(1:round(nois_pct(cItr)*nr_patients));
                    a=gmix_struct(roi).mean_patient-3*gmix_struct(roi).sd_patient;
                    b=gmix_struct(roi).mean_patient+3*gmix_struct(roi).sd_patient;
                    %data_patients(noise_indx,roi)=a+(b-a).*rand(length(noise_indx),1);
                    data_patients(noise_indx,roi)=unifrnd(a,b,length(noise_indx),1);
                end
        end
        
    end
end
%% Data Visualization

figure;
for i_var = 1:nr_events
    
    subplot(3, 4, i_var), hold on
    [hist_c, x_c] = ksdensity(data_controls(:, i_var));
    plot(x_c, hist_c, 'g');
    [hist_p, x_p] = ksdensity(data_patients(:, i_var));
    plot(x_p, hist_p, 'r');
    title(['Biomarker ' num2str(i_var)])
    set(gcf,'Color',[1 1 1])
end

%% fit the model

I_select=[1:nr_events];
data_controls=data_controls(:,I_select);
data_patients=data_patients(:,I_select);

%% Performing MCMC

% Performing mcmc
parm_mcmc.nr_gradient_ascent = 5;
parm_mcmc.nr_it_gradient_ascent = 2e3;
parm_mcmc.nr_it_burnin = 1e4;
parm_mcmc.nr_it_mcmc = 1e4;
parm_mcmc.interval_display = 1e3;
parm_mcmc.flag_sum = 2; % This should always be set to 2.
parm_mcmc.nr_it_check_p_false = 1e2;
parm_mcmc.range_p_false = [0.05 0.1
    0.05 0.1];
parm_mcmc.std_p_false = [1e-4 1e-4]';
parm_mcmc.idx_clinevent = [];


parm_struct = EBDPMCMCTestWithParams1k(data_controls, ... 
    data_patients, parm_mcmc);
 
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

   max_length_str = max([max_length_str length(name_variables{I_select(roi)})]);

end
mat_labels = zeros(nr_events, max_length_str);
for roi = 1:nr_events,

    length_str = length(name_variables{I_select(order_roi(roi))});
    mat_labels(roi, 1:length_str) = name_variables{I_select(order_roi(roi))};
end
figure, clf,
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

figure
for i=1:5
subplot(2,3,i);hist(parm_struct.event_order_mcmc(i,:));title(['Position ' num2str(i)]); set(gcf,'Color',[ 1 1 1]);
end
i
