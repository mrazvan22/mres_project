% Script to explain the event-based disease progression model when applied
% to Jonathan

clear
close all
clc

% First setting up paths
addpath('/Users/jaghajan/Documents/Experiments/ADNI_BIOMARKERS_Svn/netlab')
% All functions are in repos/Code1 and repos/Code1/netlab

% Reading in Jonathan's data
%file_data = '/cs/research/vision/camino1/camino/fonteijn/Progression/BartlettData/data_biomarkers_ADNI.mat';
file_data = 'BartlettData/data_biomarkers_ADNI.mat';
load(file_data);

name_variables = data_struct.name_var;
data = data_struct.data;
% This is somewhat cumbersome: the first variable is the subject id, which
% is actually a string and therefore not included in the data matrix. The
% data belonging to a variable with variable name name_variables{x} can
% therefore be found in data(:, x-1)
% I'm correcting for that here
I_keep = [2:length(name_variables)];
name_variables = name_variables(I_keep);

% There are quite a lot of patients in which one of the measurements is missing.
% Throwing out all of these patients
I_missing = find(sum(isnan(data), 2));
data(I_missing, :) = [];

% There are control definitions: the last column are the clinically defined
% controls (controls == 0, patients == 1), the column before that the
% CSF-defined controls
id_controls_clin = data(:, end);
id_controls_csf = data(:, end-1);
% Splitting off these variables to avoid confusion
data(:, (end-1):end) = [];

% Not sure what to do with the CDR scores and avtot, so I'll throw them out for now
% The CDR scores are either 0, 0.5 and 1, so fitting continuous distributions to those
% is problematic ***jania:***  I also removed the mmscore as it is discrete
% and a MoG does not do a good job of fitting it
I_keep = [1:4 7:12];
data = data(:, I_keep);
name_variables = name_variables(I_keep);

[nr_subj, nr_var] = size(data);

cluster_version=3;
%1 = EM MoG on single biomarkrs
%12 = version 1.2 same as 1 but with random initialization
%2 = EM MoG on entire data
%22 = same as two but on selected biomarkers
%3 = EM on PCA Dimensioanlity reduecd data
%4 = k-means clustering using all of the biomarkers

switch cluster_version
    
    case 1 %Clustered using a single biomarker at a time
        
        three_flag=1; %if three_flag=1 select controls common to the first three biomarkers only
        rand_flag=0;
        if three_flag
            outName='First_three_Single_Biomarker_Cluster_SinglePlusTwoGaussRsults.mat';
        else
            outName='Single_Biomarker_Cluster_SinglePlusTwoGaussRsults.mat';
        end
        
        % ************************************************************************
        % *** Visual inspection of that two of the biomarkers 'hippovol'
        % *** and 'hippoatrophy' do not discriminate between controls and
        % *** patients hence their clustering results are removed
        % ************************************************************************
        
        
        % ************************************************************************
        % *** Choose only those individuals as controls that were clustered as
        % *** controls in all remaining 8 biomarkers
        % ************************************************************************
        
        common_control_indx=findCommonControlIndices(three_flag,rand_flag);
        patient_indx=setdiff([1:nr_subj],common_control_indx);
        
        data_controls=data(common_control_indx,:);
        data_patients=data(patient_indx,:);
        
        
        figure(1), clf
        for i_var = 1:nr_var,
            
            subplot(3, 4, i_var), hold on
            [hist_c, x_c] = ksdensity(data_controls(:, i_var));
            plot(x_c, hist_c, 'g');
            [hist_p, x_p] = ksdensity(data_patients(:, i_var));
            plot(x_p, hist_p, 'r');
            title(name_variables{i_var})
            
        end
        set(gcf,'Color',[1 1 1])
        
        % Visual inspection of the histograms shows that some variables are higher
        % in patients than in controls. I'm going to keep the convention in the
        % modelig that values should be lower in patients than in controls
        
        %I_swap = [1 3 4 7 9 10 11];
        I_swap = [1 3 4 6 8 9 10];
        
        data_controls(:, I_swap) = -data_controls(:, I_swap);
        data_patients(:, I_swap) = -data_patients(:, I_swap);
        figure(2), clf
        for i_var = 1:nr_var,
            
            subplot(3, 4, i_var), hold on
            [hist_c, x_c] = ksdensity(data_controls(:, i_var));
            plot(x_c, hist_c, 'g');
            [hist_p, x_p] = ksdensity(data_patients(:, i_var));
            plot(x_p, hist_p, 'r');
            title(name_variables{i_var})
            
        end
        set(gcf,'Color',[1 1 1])
        
    case 12 %Clustered using a single biomarker at a time with random initialization
        
        three_flag=1; %if three_flag=1 select controls common to the first three biomarkers only
        rand_flag=1; % with random initialization
        if three_flag
            outName='First_three_Single_Bio_RandInit_Cluster_SinglePlusTwoGaussRsults.mat';
        else
            outName='Single_Bio_RandInit_Cluster_SinglePlusTwoGaussRsults.mat';
        end
        
        % ************************************************************************
        % *** Visual inspection of that two of the biomarkers 'hippovol'
        % *** and 'hippoatrophy' do not discriminate between controls and
        % *** patients hence their clustering results are removed
        % ************************************************************************
        
        
        % ************************************************************************
        % *** Choose only those individuals as controls that were clustered as
        % *** controls in all remaining 8 biomarkers
        % ************************************************************************
        
        common_control_indx=findCommonControlIndices(three_flag,rand_flag);
        patient_indx=setdiff([1:nr_subj],common_control_indx);
        
        data_controls=data(common_control_indx,:);
        data_patients=data(patient_indx,:);
        
        
        figure(1), clf
        for i_var = 1:nr_var,
            
            subplot(3, 4, i_var), hold on
            [hist_c, x_c] = ksdensity(data_controls(:, i_var));
            plot(x_c, hist_c, 'g');
            [hist_p, x_p] = ksdensity(data_patients(:, i_var));
            plot(x_p, hist_p, 'r');
            title(name_variables{i_var})
            
        end
        set(gcf,'Color',[1 1 1])
        
        % Visual inspection of the histograms shows that some variables are higher
        % in patients than in controls. I'm going to keep the convention in the
        % modelig that values should be lower in patients than in controls
        
        %I_swap = [1 3 4 7 9 10 11];
        I_swap = [1 3 4 6 8 9 10];
        
        data_controls(:, I_swap) = -data_controls(:, I_swap);
        data_patients(:, I_swap) = -data_patients(:, I_swap);
        figure(2), clf
        for i_var = 1:nr_var,
            
            subplot(3, 4, i_var), hold on
            [hist_c, x_c] = ksdensity(data_controls(:, i_var));
            plot(x_c, hist_c, 'g');
            [hist_p, x_p] = ksdensity(data_patients(:, i_var));
            plot(x_p, hist_p, 'r');
            title(name_variables{i_var})
            
        end
        set(gcf,'Color',[1 1 1])
        
    case 2 % Clustered using all of the biomarkers
        outName='ND_AllData_Cluster_SinglePlusTwoGaussRsults.mat';
        %outName='ND_AllData_Cluster_SinglePlusTwoGaussRsults_ClinControl.mat';
        
        load('Clustering_Results/Clustering_Results_AllData.mat');
        %load('Clustering_Results/Clustering_Results_AllData_ClinControl.mat');
        %noticed in PCA clustering 1st component is control and the second
        %component is patient
        patient_indx=find(label==1);
        control_indx=find(label==2);
        data_controls=data(control_indx,:);
        data_patients=data(patient_indx,:);
        
        % Checking the distributions of all variables that remain in
        % controls(green) and patients (red)
        
        figure(1), clf
        for i_var = 1:nr_var,
            
            subplot(3, 4, i_var), hold on
            [hist_c, x_c] = ksdensity(data_controls(:, i_var));
            plot(x_c, hist_c, 'g');
            [hist_p, x_p] = ksdensity(data_patients(:, i_var));
            plot(x_p, hist_p, 'r');
            title(name_variables{i_var})
            
        end
        % ************************************************************************
        % *** Visual inspection of the histograms shows that we do n ot need
        % *** to flip any variables
        % ************************************************************************
        
        
    case 22 % Clustered using all of the biomarkers
        outName='ND_selected_biomarkers_Data_Cluster_SinglePlusTwoGaussRsults.mat';
        %outName='ND_AllData_Cluster_SinglePlusTwoGaussRsults_ClinControl.mat';
        
        load('Clustering_Results/Clustering_Results_selected_biomarkers_Data.mat');
        %load('Clustering_Results/Clustering_Results_AllData_ClinControl.mat');
        %noticed in PCA clustering 1st component is control and the second
        %component is patient
        patient_indx=find(label==1);
        control_indx=find(label==2);
        data_controls=data(control_indx,:);
        data_patients=data(patient_indx,:);
        
        % Checking the distributions of all variables that remain in
        % controls(green) and patients (red)
        
        figure(1), clf
        for i_var = 1:nr_var,
            
            subplot(3, 4, i_var), hold on
            [hist_c, x_c] = ksdensity(data_controls(:, i_var));
            plot(x_c, hist_c, 'g');
            [hist_p, x_p] = ksdensity(data_patients(:, i_var));
            plot(x_p, hist_p, 'r');
            title(name_variables{i_var})
            
        end
        % ************************************************************************
        % *** Visual inspection of the histograms shows that we do n ot need
        % *** to flip any variables
        % ************************************************************************
        
        
    case 3  % Clustered using PCA dimensionality reduced data
        outName='PCA_Cluster_MixTResults.mat';
        load('Clustering_Results/MixT_PCA_Clustering_Results.mat');
        %noticed in PCA clustering 1st component is patient and the second
        %component is control
        patient_indx=find(label==1);
        control_indx=find(label==2);
        data_controls=data(control_indx,:);
        data_patients=data(patient_indx,:);
        
        % Checking the distributions of all variables that remain in
        % controls(green) and patients (red)
        
        figure(1), clf
        for i_var = 1:nr_var,
            
            subplot(3, 4, i_var), hold on
            [hist_c, x_c] = ksdensity(data_controls(:, i_var));
            plot(x_c, hist_c, 'g');
            [hist_p, x_p] = ksdensity(data_patients(:, i_var));
            plot(x_p, hist_p, 'r');
            title(name_variables{i_var})
            
        end
        % ************************************************************************
        % *** Visual inspection of the histograms shows that we do n ot need
        % *** to flip any variables
        % ************************************************************************
       
        %I_swap = [1 3 4 7 9 10 11];
        I_swap = [1 3 4 6 8 9 10];
        
        data_controls(:, I_swap) = -data_controls(:, I_swap);
        data_patients(:, I_swap) = -data_patients(:, I_swap);
        data(:,I_swap)=-data(:,I_swap);
        
        figure(2), clf
        for i_var = 1:nr_var,
            
            subplot(3, 4, i_var), hold on
            [hist_c, x_c] = ksdensity(data_controls(:, i_var));
            plot(x_c, hist_c, 'g');
            [hist_p, x_p] = ksdensity(data_patients(:, i_var));
            plot(x_p, hist_p, 'r');
            title(name_variables{i_var})
            
        end
        set(gcf,'Color',[1 1 1])
    
    case 4
        outName='K-means_Cluster_SinglePlusTwoGaussResults.mat'
        load('Clustering_Results/K-means_Clustering_Results.mat');
        patient_indx=find(label==1);
        control_indx=find(label==2);
        data_controls=data(control_indx,:);
        data_patients=data(patient_indx,:);
        
        % Checking the distributions of all variables that remain in
        % controls(green) and patients (red)
        
        figure(1), clf
        for i_var = 1:nr_var,
            
            subplot(3, 4, i_var), hold on
            [hist_c, x_c] = ksdensity(data_controls(:, i_var));
            plot(x_c, hist_c, 'g');
            [hist_p, x_p] = ksdensity(data_patients(:, i_var));
            plot(x_p, hist_p, 'r');
            title(name_variables{i_var})
            
        end
        
end






%%  Now fitting the EBDP model to the data only using the selected controls

%cov type parameter
%covType = 0;  %unconstrained - individual covariance matrices
%covType = 1;  %unconstrained - shared covariance matrices
%covType = 2;  %diagonal - individual covariance matrices
%covType = 3;  %diagonal - shared covariance matrices
covType = 4;  %spherical - individual covariance matrices
%covType = 5;  %spherical - shared covariance matrices


like_version=1;
%like_version=1 fit a mixture of T-distributions
%like_version=2 fit a mixture of Uniform and T-distributions

%% transpose the data to fit the mix-T 
data=data';
data_controls=data_controls';
data_patients=data_patients';
prior_control=size(data_controls,2)/nr_subj;
nr_events=size(data_controls,1);

threshold_flag=1;
min_prior_flag=1;
for roi=1:nr_events
    
    switch like_version
        case 1 % fit a mixture of T-distributions 
            figure
            tEst_control = mixTFit(data_controls(roi,:),1,covType);
            figure
            tEst_tot = mixTFitFixedComp(data(roi,:),2,covType,tEst_control(1).mean,prior_control,min_prior_flag);
            tEst_patient=tEst_tot(2);
            tEst_patient.prior=1;
            
            % find likelihood event_occurred and event_not_occured
            likelihood(roi, :, 1) = exp(getLogTProb(data_patients(roi,:),tEst_control.mean,tEst_control.cov,tEst_control.dof));
            likelihood(roi, :, 2) = exp(getLogTProb(data_patients(roi,:),tEst_patient.mean,tEst_patient.cov,tEst_patient.dof));
            
            
            
            % Visualize how well the T-distributions fit the data
            
            p_controls=exp(getLogTProb(data_controls(roi,:),tEst_control.mean,tEst_control.cov,tEst_control.dof));
            p_patients=exp(getLogTProb(data_patients(roi,:),tEst_patient.mean,tEst_patient.cov,tEst_patient.dof));
            p_all_MoT=getExpMixTProb(data(roi,:),tEst_tot);
            
            [hist_c, x_c] = ksdensity(data_controls(roi,:));
            [hist_p, x_p] = ksdensity(data_patients(roi,:));
            
            
            figure
            subplot(121), hold on
            plot(x_c, hist_c, 'g');
            plot(x_p, hist_p, 'r');
            
            subplot(122), hold on
            plot(data_controls(roi,:), p_controls, 'g.')
            plot(data_patients(roi,:), p_patients, 'r.');
            plot(data(roi,:), p_all_MoT, 'b.');
            
            set(gcf,'Color',[1 1 1])
            pause
            
        case 2 % fit a mixture of Uniform and T-distributions 
            figure
            tEst_control = mixTFit(data_controls(roi,:),1,covType);
            
            parm_mixtUnif.mu_init = tEst_control.mean;
            parm_mixtUnif.sig2_init = tEst_control.cov;
            parm_mixtUnif.dof_init=tEst_control.dof;
            
            parm_mixtUnif.a = min(data(roi,:));
            parm_mixtUnif.b = tEst_control.mean;
            parm_mixtUnif.nr_it = 1e3;
            
            parm = fitFixedTDistUnif(data(roi,:), parm_mixtUnif);
            parm.data = data(roi,:);
            gmix_struct.gmix_roi{roi} = parm;
            
            likelihood(roi, :, 1) = parm.prob_T(I_patients);
            likelihood(roi, :, 2) = parm.prob_U(I_patients);
            
    end
    
    %% set a maximumn likelilhood ratio
    if threshold_flag
        T=[3.2];
        ratio=abs(log(likelihood(roi, :, 1))-log(likelihood(roi, :, 2)));
        indx=find(ratio>log(T));
        
        for cIndx=1:length(indx)
            thisLike=likelihood(roi,indx(cIndx),:);
            %if you are in the far left of the distributrion
            if data_patients(roi,indx(cIndx))<tEst_control.mean
                minIndx=find(thisLike==min(thisLike));
                mxIndx=find(thisLike==max(thisLike));
                thisLike(minIndx)=thisLike(mxIndx)/T;
            else
                % if you are in the far right of the distribution
                minIndx=find(thisLike==min(thisLike));
                mxIndx=find(thisLike==max(thisLike));
                thisLike(mxIndx)=thisLike(minIndx)*T;
            end
            likelihood(roi,indx(cIndx),:)=thisLike;
        end
    end
end


%% 
% %visuaize the likelihoods
% for roi=1:length(I_select)
% figure
% plot(data_patients(roi,:),log(likelihood(roi,:,2))-log(likelihood(roi,:,1)),'r.');
% hold on
% plot(data_patients(roi,:),log(likelihood(roi,:,1))-log(likelihood(roi,:,2)),'g.');
% set(gcf,'Color',[1 1 1])
% pause
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

save(['Mat_Data/' outName])




