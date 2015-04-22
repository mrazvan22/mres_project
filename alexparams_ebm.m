function [] = alexparams_ebm()

% ADNIdataBL raw data as received from the ADNI website
% EBMdataBL - pre-processed data for use in the EBM model
% EMBLdxBL - the level of impairment for each patient (1 - Cognitively normal, 2 - MCI 3 - AD)
% EBMevents - labels of the EBM events
load('alex_data/ADNIdata_Baseline.mat')

%[mu_mix, sigma_mix, pi_mix] = calc_gaussian_parameters(EBMdataBL, EBMdxBL)

load('alex_params.mat')
%[max_seq_grad,max_lik_grad,final_sequences_grad,final_lik_grad] = char_seq_grad_asc(EBMdataBL, mu_mix, sigma_mix, pi_mix);
%save('alex_grad_asc.mat', 'max_seq_grad','max_lik_grad','final_sequences_grad','final_lik_grad');

load('alex_grad_asc.mat');
% max likelihood labels
EBMevents(max_seq_grad)

%[samples, acceptance_rate] = MCMC_sampling(EBMdataBL, mu_mix, sigma_mix, pi_mix, max_seq_grad);
%save('alex_MCMC2_adj.mat','samples','acceptance_rate');

load('alex_MCMC.mat');
[all_pop_matrix] = calc_variance_diagrams(samples,max_seq_grad);

%norm_mat = normaliseVarMatrix(all_pop_matrix);
h = plotVarMatrix(all_pop_matrix, max_seq_grad);
%saveas(h, 'figures/alex_mat_max_like', 'png')

%           [1,2,3,4,5,6,7 ,8,9,10,11,12,13,14]
%matrixSeqPerm = [1,2,3,5,7,4,11,6,8,9, 12,13,14,10]
%matrixSeq = max_seq_grad(matrixSeqPerm)

char_seq = calc_char_sequence(samples);
[matrix_char_seq] = calc_variance_diagrams(samples,char_seq);
%norm_mat_char_seq = normaliseVarMatrix(matrix_char_seq);

h = plotVarMatrix(matrix_char_seq, char_seq);
%saveas(h, 'figures/alex_mat_char_seq', 'png')

calc_likelihood(EBMdataBL, max_seq_grad, mu_mix, sigma_mix, pi_mix)
calc_likelihood(EBMdataBL, char_seq, mu_mix, sigma_mix, pi_mix)
end
