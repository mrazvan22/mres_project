function log_likely = singleBiomkLikelihood(x, data)

mu = x(1:2);
sigma = x(3:4);
pi_control = x(5);
pi = [pi_control, 1-pi_control];

% Evaluate the log_likelihood
K = 2;
N = length(data); 
eval_matrix = zeros(N,K);
for k=1:K
   eval_matrix(:,k) = pi(k) * normpdf(data, mu(k), sigma(k));
end

sumK = sum(eval_matrix, 2);
assert(length(sumK) == N);
log_likely = -sum(log(sumK));

end