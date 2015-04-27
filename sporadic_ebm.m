function [] = sporadic_ebm()

% ADNIdataBL raw data as received from the ADNI website
% EBMdataBL - pre-processed data for use in the EBM model
% EMBLdxBL - the level of impairment for each patient (1 - Cognitively normal, 2 - MCI 3 - AD)
% EBMevents - labels of the EBM events
load('alex_data/ADNIdata_Baseline.mat')

%[mu_mix, sigma_mix, pi_mix] = calc_gaussian_parameters(EBMdataBL, EBMdxBL, @em_mix)
[mu_mix, sigma_mix, pi_mix] = calc_gaussian_parameters(EBMdataBL, EBMdxBL, @fmincon_fit)
save('sporadic_params.mat', 'mu_mix', 'sigma_mix', 'pi_mix'); 

load('sporadic_params.mat');
a = load('alex_params.mat');
sum(sum(abs(a.mu_mix - mu_mix)))
sum(sum(abs(a.sigma_mix - sigma_mix)))

load('sporadic_params.mat')
%[max_seq_grad,max_lik_grad,final_sequences_grad,final_lik_grad] = char_seq_grad_asc(EBMdataBL, a.mu_mix, sigma_mix, pi_mix);
[max_seq_grad,max_lik_grad,final_sequences_grad,final_lik_grad] = char_seq_grad_asc(EBMdataBL, mu_mix, sigma_mix, pi_mix);
save('grad_asc.mat', 'max_seq_grad','max_lik_grad','final_sequences_grad','final_lik_grad');

load('grad_asc.mat');
% max likelihood labels
EBMevents(max_seq_grad)

%[samples, acceptance_rate] = MCMC_sampling(EBMdataBL, mu_mix, sigma_mix,pi_mix, max_seq_grad);
%save('MCMC2.mat','samples','acceptance_rate');
%[~, max_lik_event_positions] = sort(max_seq_grad) % compare these to Alex's results, not the max_seq_grad which is the permutation


load('MCMC.mat');
[all_pop_matrix] = calc_variance_diagrams(samples,max_seq_grad);

%norm_mat = normaliseVarMatrix(all_pop_matrix);
h = plotVarMatrix(all_pop_matrix, max_seq_grad);


end


