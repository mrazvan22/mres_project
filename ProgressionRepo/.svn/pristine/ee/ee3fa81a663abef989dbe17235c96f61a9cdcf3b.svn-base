function mix_struct = fitFixedTDistUnif(x, parm)
% mix_struct = fitFixedGaussUnif(x, parm)
% This function fits a mixture of a Gaussian and a Uniform to the data
% The parameters of the Gaussian are kept fixed, the only free parameter is
% therefore the volume fraction of the uniform distribution (the range of
% the Uniform distribution is fixed as well).

a = parm.a;
b = parm.b;
nr_it = parm.nr_it;
mu = parm.mu_init; % I'm keeping the Gaussian parameters fixed here
sig2 = parm.sig2_init;
dof=parm.dof_init;
f_T = rand;
f_U = 1 - f_T;
nr_dp = length(x);

for it = 1:nr_it,
    
    prob_U = zeros(nr_dp, 1);
    prob_U((x >= a) & (x <= b)) = 1/(b-a);
    tEst_T.mean = mu;
    tEst_T.cov = sig2;
    tEst_T.dof=dof;
    prob_T = exp(getLogTProb(x,tEst_T.mean,tEst_T.cov,tEst_T.dof));
    prob_T=prob_T';
    
    norm_fact = (f_U*prob_U) + (f_T*prob_T);
    post_U = (f_U*prob_U)./norm_fact;
    post_T = (f_T*prob_T)./norm_fact;
    
    f_U = sum(post_U)/nr_dp;
    f_T = sum(post_T)/nr_dp;
    %     mu = mean(post_N.*x)/mean(post_N);
    %     sig2 = mean(post_N.*((x - mu).^2))/mean(post_N);
    
end

mix_struct.mu = mu;
mix_struct.sig2 = sig2;
mix_struct.f_U = f_U;
mix_struct.f_T = f_T;
mix_struct.prob_U = prob_U;
mix_struct.prob_T = prob_T;
mix_struct.post_U = post_U;
mix_struct.post_T = post_T;