function [samples,acceptance_rate] = MCMC_sampling(X, mu_mix, sigma_mix, init_seq)

[nr_subjects,nr_biomk] = size(X);


[nr_biomk2,~] = size(mu_mix);
[nr_biomk3,~] = size(sigma_mix);

assert(nr_biomk == nr_biomk2 && nr_biomk == nr_biomk3);

% initialise curr_seq to initial sequence
curr_seq = init_seq;

debug = 0;

if(debug == 0)
  burnout_iterations = 10^5;
  actual_iterations = 10^6;
else
  burnout_iterations = 1000;
  actual_iterations = 100;
end

old_likelihood = calc_likelihood(X, curr_seq, mu_mix, sigma_mix);

s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);

samples = zeros(actual_iterations, nr_biomk);

acceptance_rate = 0;

for i=1:(burnout_iterations + actual_iterations)
    i
    p1 = ceil(rand * 14);
    p2 = ceil(rand * 14);
    while (p1 == p2)
        p2 = ceil(rand * 14);
    end
    new_seq = curr_seq;
    
    % swap the 2 events 
    tmp = new_seq(p1);
    new_seq(p1) = new_seq(p2);
    new_seq(p2) = tmp;
    
    new_likelihood = calc_likelihood(X, new_seq, mu_mix, sigma_mix);
    likely_ratio = new_likelihood - old_likelihood; % subtract because it's log likelihood
    
    if(log(rand) < likely_ratio) % take log of rand because the likelihood was in logspace
        % swapt them with probability min(likely_ratio, 1)
        curr_seq = new_seq;
        old_likelihood = new_likelihood;
        acceptance_rate = acceptance_rate + 1;
    end
    
    % if burnout is over record the samples
    if (i > burnout_iterations)
        samples(i-burnout_iterations,:) = curr_seq;
    end
    
    if(mod(i,1000) == 0)
       display('1000 iterations') 
    end
    
end

acceptance_rate = acceptance_rate/(burnout_iterations + actual_iterations);


end