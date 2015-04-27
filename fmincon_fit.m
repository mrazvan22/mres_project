function [mu, sigma, pi]  = ...
    fmincon_fit(data, mu_init, sigma_init, pi_init)

% data is a 1d array of biomarker levels which is assumed to be generated
% from a mixture of 2 gaussian distributions

%mu_init(2) = 27

max_sigma = sigma_init;
%min_sigma = max_sigma / 5;
min_sigma = max_sigma/5;

if (mu_init(1) < mu_init(2))
  min_mu = [0 mu_init(2)];
  max_mu = [mu_init(1) inf];
else
  min_mu = [mu_init(1) 0];
  max_mu = [inf mu_init(2)];
end

% pi_control is the weight of the control gaussian
% pi(1) - control pi(2) - patient
pi_init = [0.5 0.5];

%N = length(data); 
K = 2;

assert(length(mu_init) == K && length(sigma_init) == K && length(pi_init) == K);

lb = [min_mu(1), min_mu(2), min_sigma(1), min_sigma(2), 0];
ub = [max_mu(1), max_mu(2), max_sigma(1), max_sigma(2), 1]; 

fminconOptions = optimset('MaxFunEvals', 20000, 'Algorithm', 'interior-point',...
    'TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'iter');

startX = [mu_init, sigma_init, pi_init(1)];
x = fmincon(@singleBiomkLikelihood, startX, [],[],[],[],lb, ub, [], fminconOptions, data);

[mu, sigma, pi_control] = deal(x(1:2), x(3:4), x(5));

pi = [pi_control, 1-pi_control];

end