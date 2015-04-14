function [] = sporadic_ebm()

% ADNIdataBL raw data as received from the ADNI website
% EBMdataBL - pre-processed data for use in the EBM model
% EMBLdxBL - the level of impairment for each patient (1 - Cognitively normal, 2 - MCI 3 - AD)
% EBMevents - labels of the EBM events
load('alex_data/ADNIdata_Baseline.mat')

%[mu_mix, sigma_mix, pi_mix] = calc_gaussian_parameters(EBMdataBL, EBMdxBL)

load('sporadic_params.mat')
%[max_seq_grad,max_lik_grad,final_sequences_grad,final_lik_grad] = char_seq_grad_asc(EBMdataBL, mu_mix, sigma_mix);
%save('grad_asc.mat', 'max_seq_grad','max_lik_grad','final_sequences_grad','final_lik_grad');

load('grad_asc.mat');
%[samples, acceptance_rate] = MCMC_sampling(EBMdataBL, mu_mix, sigma_mix, max_seq_grad);
%save('MCMC2.mat','samples','acceptance_rate');
%[~, max_lik_event_positions] = sort(max_seq_grad) % compare these to Alex's results, not the max_seq_grad which is the permutation

load('MCMC.mat');
[char_sequence,char_perm,all_pop_matrix] = calc_variance_diagrams(samples);

calc_likelihood(EBMdataBL, 1:14, mu_mix, sigma_mix)

end

function [char_sequence,char_perm,all_pop_matrix] = calc_variance_diagrams(samples)
    
biomk_avg_ordering = mean(samples);
[~, char_sequence] = sort(biomk_avg_ordering); % for comparing with Alex's seq of events in the matrix
char_sequence
[~, char_perm] = sort(char_sequence); % for permuting the matrix

NR_BIOMK = length(biomk_avg_ordering);
all_pop_matrix = zeros(NR_BIOMK, NR_BIOMK); % all_pop_matrix(event, position);

for s=1:size(samples,1)
  for biomk = 1:NR_BIOMK
    all_pop_matrix(biomk, samples(s,biomk)) = all_pop_matrix(biomk, samples(s,biomk)) + 1;
  end
end

%permute the matrix using the characteristic/max_likelihood ordering
 
all_pop_matrix = all_pop_matrix(char_perm,:);

end
