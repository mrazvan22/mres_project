function [mu_old, sigma_old, pi]  = ...
    gradient_ascent_fit(data, mu_init, sigma_init, pi_init)

% data is a 1d array of biomarker levels which is assumed to be generated
% from a mixture of 2 gaussian distributions

mu_old = mu_init;
sigma_old = max_sigma;
min_sigma = max_sigma / 5;

% pi_control is the weight of the control gaussian
% pi(1) - control pi(2) - patient
pi = [0.5 0.5];

N = length(data); 
K = 2;

assert(length(mu_old) == K && length(sigma_old) == K && length(pi) == K);

eval_matrix = zeros(N,K);

log_likely = inf;
thresh = 0.0001;

%iterations = 1000;
old_likelihood = singleBiomkLikelihood(data, mu_old, sigma_old, pi);
new_likelihood = old_likelihood + 1;

while old_likelihood < new_likelihood
  
  
end

end