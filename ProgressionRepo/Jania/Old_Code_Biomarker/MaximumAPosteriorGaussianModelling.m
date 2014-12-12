function MaximumAPosteriorGaussianModelling
close all
% define the likelihood distribution
like_mean=0;
like_sd=1;

%gerenrate data from the defined distribution
x=randn(1,100);
y=getGaussProb(x,like_mean,like_sd);


%plot the likelihood distribution
subplot(1,5,1);plot(x,y,'b.');
title('Likelihood Gaussian Distribution');

%define priors over the mean and the variance of the gaussian dist
prior_mean=like_mean;

%set the sigma of the prior normal distribution such that
%min and max values of current biomarker are within + / - 3sigma
minX=min(x);
maxX=max(x);
sigmaMin=(prior_mean-minX)/3;
sigmaMax=(prior_mean)/3;
prior_sigma=(sigmaMin+sigmaMax)/2;

% set the alpha and beta hyperparameters of tge gamma distribution
%such that alpha*beta=sigma and 90% of populationof the gamma
%distribution is between 0 and sigma
a=(1/25);
b=like_sd*25;

prior_thisMean=normpdf(like_mean,prior_mean,prior_sigma);
prior_thisCov=gampdf(like_sd,a,b);

%plot the prior distriutions
hold on
subplot(1,5,2);plot([-5:0.1:5],normpdf([-5:0.1:5],prior_mean,prior_sigma),'r.');
title('Prior Gaussian Distribution');
subplot(1,5,3);plot([0:0.1:10],gampdf([0:0.1:10],a,b),'g.');
title('Prior Gamma Distribtion')

%plot the posterior
posterior=x.*prior_thisMean.*prior_thisCov;
%subplot(1,5,4);plot(x,posterior,'m.');
%title('Posterior Gaussian Distribution');
subplot(1,5,4);plot(x,posterior,'m.');
plot([-5:0.1:5],getGaussProb([-5:0.1:5],like_mean,like_sd).*normpdf([-5:0.1:5],prior_mean,prior_sigma).*gampdf([0:0.1:10],a,b));
i