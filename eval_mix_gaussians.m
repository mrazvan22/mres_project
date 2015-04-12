function Y = eval_mix_gaussians(X, mu, sigma, pi)
% X is a vector of data points
% mu is a vector of means
% sigma is a vector of std deviations

K = length(mu);
Y = zeros(size(X));

for k=1:K
   if (sigma(k) ~= 0) 
       Y = Y + pi(k) * normpdf(X, mu(k), sigma(k)); 
   end
end

end