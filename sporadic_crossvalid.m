function [] = sporadic_crossvalid()


% ADNIdataBL raw data as received from the ADNI website
% EBMdataBL - pre-processed data for use in the EBM model
% EMBLdxBL - the level of impairment for each patient (1 - Cognitively normal, 2 - MCI 3 - AD)
% EBMevents - labels of the EBM events
load('alex_data/ADNIdata_Baseline.mat')

RAND_SHUFFLE= 1;
BOOTSTRAPS = 100;
REPLACEMENT = true;
[NR_SUBJECTS, NR_BIOMK] = size(EBMdataBL);
global DEBUG;
DEBUG = false;

[muBoot, sigmaBoot, piBoot] = deal(zeros(BOOTSTRAPS, NR_BIOMK, 2), zeros(BOOTSTRAPS, NR_BIOMK, 2), zeros(BOOTSTRAPS, NR_BIOMK, 2));
[bootData, bootDiag] = deal(zeros(BOOTSTRAPS, NR_SUBJECTS, NR_BIOMK), zeros(BOOTSTRAPS, NR_SUBJECTS));

%% find the Gaussian distribution parameters
for b=1:BOOTSTRAPS
  b
  sampleIndices = randsample(NR_SUBJECTS,NR_SUBJECTS,REPLACEMENT);
  bootData(b,:,:) = EBMdataBL(sampleIndices, :);
  bootDiag(b,:) = EBMdxBL(sampleIndices, :);
  [muBoot(b,:,:), sigmaBoot(b,:,:), piBoot(b,:,:)] = calc_gaussian_parameters(...
    squeeze(bootData(b,:,:)), squeeze(bootDiag(b,:)), @fmincon_fit);
end

save('cross_params.mat', 'muBoot', 'sigmaBoot', 'piBoot', 'bootData', 'bootDiag'); 

%% Find the maximum likelihood sequence
load('cross_params.mat');
bootSeq = zeros(BOOTSTRAPS, NR_BIOMK);
for b=1:BOOTSTRAPS
  b
  [bootSeq(b,:),max_lik_grad,final_sequences_grad,final_lik_grad] = ...
    char_seq_grad_asc(squeeze(bootData(b,:,:)), squeeze(muBoot(b,:,:)), ...
    squeeze(sigmaBoot(b,:,:)), squeeze(piBoot(b,:,:)), RAND_SHUFFLE);
end
save('boot_samples.mat', 'bootSeq');

%% Plot the positional variance matrix
load('boot_samples.mat');
noboot = load('grad_asc.mat');

[boot_all_matrix] = calc_variance_diagrams(bootSeq,noboot.max_seq_grad); % bootstrap matrix for all subjects

h = plotVarMatrix(boot_all_matrix, noboot.max_seq_grad);
saveas(h, 'figures/bootPosVarAll', 'png');
hgexport(h, 'figures/bootPosVarAll.eps');

end



