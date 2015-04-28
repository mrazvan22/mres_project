function [] = sporadic_ebm()

% ADNIdataBL raw data as received from the ADNI website
% EBMdataBL - pre-processed data for use in the EBM model
% EMBLdxBL - the level of impairment for each patient (1 - Cognitively normal, 2 - MCI 3 - AD)
% EBMevents - labels of the EBM events
load('alex_data/ADNIdata_Baseline.mat')

%[mu_mix, sigma_mix, pi_mix] = calc_gaussian_parameters(EBMdataBL, EBMdxBL, @em_mix)
%[mu_mix, sigma_mix, pi_mix] = calc_gaussian_parameters(EBMdataBL, EBMdxBL, @fmincon_fit)
%save('sporadic_params.mat', 'mu_mix', 'sigma_mix', 'pi_mix'); 

load('sporadic_params.mat');
a = load('alex_params.mat');
sum(sum(abs(a.mu_mix - mu_mix)))
sum(sum(abs(a.sigma_mix - sigma_mix)))
sum(abs(a.sigma_mix - sigma_mix),2)

randShuffle = 1;
%[max_seq_grad,max_lik_grad,final_sequences_grad,final_lik_grad] = char_seq_grad_asc(EBMdataBL, mu_mix, sigma_mix, pi_mix, randShuffle);
%save('grad_asc.mat', 'max_seq_grad','max_lik_grad','final_sequences_grad','final_lik_grad');

load('grad_asc.mat');
a_grad = load('alex_grad_asc.mat')
% max likelihood labels
EBMevents(max_seq_grad)
sporadic_likelihood = calcLikelihoodFast(EBMdataBL, max_seq_grad, mu_mix, sigma_mix, pi_mix)
alex_likelihood = calcLikelihoodFast(EBMdataBL, a_grad.max_seq_grad, mu_mix, sigma_mix, pi_mix)

%[samples, acceptance_rate] = MCMC_sampling(EBMdataBL, mu_mix, sigma_mix,pi_mix, max_seq_grad, randShuffle);
%save('MCMC100k.mat','samples','acceptance_rate');

load('MCMC.mat');
[all_pop_matrix] = calc_variance_diagrams(samples,max_seq_grad);

norm_mat = normaliseVarMatrix(all_pop_matrix);
h = plotVarMatrix(all_pop_matrix, max_seq_grad);
saveas(h, 'figures/posVarianceMatrix', 'png')

patientStages = calcPatientStages(EBMdataBL, max_seq_grad, mu_mix, sigma_mix);
h = plotStagesHist(patientStages, EBMdxBL);
hgexport(h, 'figures/patientStages.eps')

end


