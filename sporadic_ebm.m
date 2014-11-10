load('alex_data/ADNIdata_Baseline.mat')

% ADNIdataBL raw data as received from the ADNI website
%EBMdataBL - pre-processed data for use in the EBM model
% EMBLdxBL - the level of impairment for each patient (1 - Cognitively normal, 2 - MCI 3 - AD)
% EBMevents - labels of the EBM events

[nr_patients, nr_biomarkers] = size(EBMdataBL)

control_indices = find(EBMdxBL == 1)
patient_indices = find(EBMdxBL > 1)

% pXgEnE(i,j,1) = p(x_ij | Ei) (patients)
% pXgEnE(i,j,2) = p(x_ij | not Ei) (controls)
pXgEnE = zeros(nr_patients, nr_biomarkers,2);

for biomk=1:nr_biomarkers
   control_levels = EBMdataBL(control_indices, biomk);
   patient_levels = EBMdataBL(patient_indices, biomk);
    
   mu_control = mean(control_levels);
   sigma_control = std(control_levels);
   
   mu_patient = mean(patient_levels);
   sigma_patient = std(patient_levels);
   
   %pXgEnE(:, biomk, 2) = normpdf(control_levels, mu_control, sigma_control);
   %pXgEnE(:, biomk, 1) = normpdf(patient_levels, mu_patient, sigma_patient);
   
   
   all_levels = EBMdataBL(:, biomk);
   % fit two gaussians on all the data
   minA = min(all_levels)
   maxA = max(all_levels)
   
   OPTIONS = statset('MaxIter',500,'Display','final')
   GMModel = gmdistribution.fit(all_levels, 2, 'Options', OPTIONS);
   
   %plot(pdf(GMModel, (minA:0.01:maxA)'))
   %hold on
   hist(all_levels,50)
   
end