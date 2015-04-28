function [mu, sigma, pi]  = ...
    em_mix(data, mu_init, sigma_init, pi_init)

% data is a 1d array of biomarker levels which is assumed to be generated
% from a mixture of 2 gaussian distributions

%mu_init(2) = 28.359

mu = mu_init;
sigma = sigma_init;
%pi = pi_init;
max_sigma = sigma_init;
min_sigma = max_sigma / 5;

% pi_control is the weight of the control gaussian
% pi(1) - control pi(2) - patient
pi = [0.5 0.5];


if (mu_init(1) < mu_init(2))
  min_mu = [0 mu_init(2)];
  max_mu = [mu_init(1) inf];
else
  min_mu = [mu_init(1) 0];
  max_mu = [inf mu_init(2)];
end

N = length(data); 
K = 2;

assert(length(mu) == K && length(sigma) == K && length(pi) == K);

znk = zeros(N,K);
eval_matrix = zeros(N,K);

log_likely = inf;
thresh = 0.0001;

iterations = 1000;

Nk = [-1 -1];
%% keep updating until convergence
for iter=1:iterations
    
    % E-step
    for k=1:K
        znk(:,k) = pi(k) * normpdf(data, mu(k), sigma(k));
    end
% 
    % normalise the znk rows
    for n=1:N
        if(sum(znk(n,:)) ~= 0)
            znk(n,:) = znk(n,:) ./ sum(znk(n,:));
        else
            display('Warning - znk is zero')
            znk(n,:) = [0.5 0.5];
        end
    end

    % M-step

    Nk = sum(znk);
    
    if(~all(Nk) || isnan(Nk(1)) || isnan(Nk(2)))
        break
    end
    
    for k=1:K 
       mu(k) =  sum(znk(:,k) .* data) /Nk(k);
       
       if (mu(k) > max_mu(k))
         %display('max mu constraint')
         mu(k) = max_mu(k);
       end
       
       if (mu(k) < min_mu(k))
         %display('min mu constraint')
         mu(k) = min_mu(k);
       end
       
       % take sqrt because we need std dev, not variance
       sigma(k) = sqrt(sum(znk(:,k) .* (data - mu(k)) .* (data - mu(k)) ) / Nk(k)); 
       
       % set the constraint that alex implemented in the paper: i.e. the sigma
       % of each distribution should not be greater than the sigma of the CN
       % and AD groups taken separately
       if (sigma(k) > max_sigma(k))
           %display('max sigma constraint')
           sigma(k) = max_sigma(k);
       end
      
       % also make sure the sigma doesn't fall below a certain threshold
       if (sigma(k) < min_sigma(k))
           %display('min sigma constraint')
           sigma(k) = min_sigma(k);
       end
       
       pi = Nk / N;
       
       % normalise the pi, although it shouldn't normally need to be
       % normalised
       pi = pi ./ sum(pi);

%        X = min(data):0.005:max(data);
%        hist(data,25)
%        hold on
%        plot(X, eval_mix_gaussians(X, mu, sigma, pi))

    end
    
    %log_likely = singleBiomkLikelihood(data, mu, sigma, pi);  

    %[mu, sigma, pi(1)]
    if(~all(sigma) || isnan(sigma(1)) || isnan(sigma(2)))
        break
    end
end
    

end