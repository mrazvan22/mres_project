function log_likelihood=getMixGaussLogLikeWithParams(data,est_gauss,prior_mean,prior_sigma,a,b)

log_likelihood=log((getGaussProb(data,est_gauss(1).mean,est_gauss(1).cov).*est_gauss(1).prior.*prior_mean(1).*prior_sigma(1)) ...
    +(getGaussProb(data,est_gauss(2).mean,est_gauss(2).cov).*est_gauss(2).prior.*prior_mean(2).*prior_sigma(2)));

log_likelihood=sum(log_likelihood);