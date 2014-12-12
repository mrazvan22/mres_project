function dataSimulationAllPermutationsWithParamsClosedFormDemo
clear
close all

%% Initialization
nr_events=3;        % select the nuber of events
nr_patients=500;    % select the nuber of patients
nr_controls=500;     % select the nuber of controls
sim_version=1;
threshold_flag=0;

rand('seed',1);
randn('seed',1);

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

for c_order=1:n_perms
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
    
    
    version_likelihood = 13;  %version_likelihood==13 fit a Mixture of Gaussians several
    %times with different random initializations
    nrAttempts=5;              %number of times to fit MoG
    [likelihood, gmix_struct] = ...
        EBDPComputeLikelihood(data_patients, data_controls, version_likelihood,threshold_flag,nrAttempts);
    
    %% Closed form solution without prior over params
    est_order(c_order,:)=findNEventOrdering(likelihood);
    
    %% closed form solution with prior over params
    [est_order_params(c_order,:),  ~, ~]= EBDPMCMCTestWithParams1kCorrectClosedForm(data_controls, data_patients)
    
end

%% calculate percentage of ecorrect estimation of ordering
for i=1:n_perms
    diff(i)=sum(abs(true_order(i,:)-est_order(i,:)));
end
pct_correct=length(find(diff==0))/n_perms*100

i