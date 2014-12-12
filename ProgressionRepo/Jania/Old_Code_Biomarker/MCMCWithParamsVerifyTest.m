function MCMCWithParamsVerifyTest

% This function to the do a sanity check on the code that performs
% MCMCwithParams to see if the code is working properly

% The idea is to generate data form a mixture of Gaussian model and see
% whether the MCMC code will recover the true means and variances

%% Initialize the Gaussian parameters
gauss(1).mean=0;
gauss(1).cov=1;
gauss(1).prior=0.25;

gauss(2).mean=-3;
gauss(2).cov=2.5;
gauss(2).prior=0.75;

nData=500;
prGaussCum = cumsum([gauss(:).prior]);

nr_it_burnin =1e4;
nr_it_mcmc = 1e4;
it_vec=1:100:nr_it_mcmc + nr_it_burnin;
gauss_mcmc_mean=zeros(2,nr_it_mcmc);
gauss_mcmc_cov=zeros(2,nr_it_mcmc);


%% Gerenate data
for (cData = 1:nData)
    %choose which Gaussian
    randVal = rand(1);
    decision = prGaussCum-randVal;
    %retain gaussians where cum prob is greater
    decision(find(decision<0))=1.1;
    %choose closest Gaussian
    gaussChosen = find(decision == min(decision));
    %generate this point and store
    data(cData) = gauss(gaussChosen).mean+gauss(gaussChosen).cov*randn;
    membership(cData)=gaussChosen;
end;

%visualize the generated dataset
figure
[f,xi] = ksdensity(data(membership==1)); 
plot(xi,f,'g');

hold on
[f2,xi2] = ksdensity(data(membership==2)); 
plot(xi2,f2,'r');


%% Model fitting
%initialize the means and variances keeping the component prior fixed as
%the true prior
covType=0;
nAttempts=1;
[init_gauss,pr_data]=mixGaussFit(data,2,covType,nAttempts);

init_gauss(1).mean=1.5;
init_gauss(1).cov=0.3;
init_gauss(1).prior=0.25;

init_gauss(2).mean=-5;
init_gauss(2).cov=0.9;
init_gauss(2).prior=0.75;


means=[init_gauss(:).mean]
class1Indx=find(means==max(means));
class2Indx=find(means==min(means));

est_gauss(1).mean=init_gauss(class1Indx).mean;
est_gauss(1).cov=init_gauss(class1Indx).cov;
est_gauss(1).prior=gauss(1).prior;

est_gauss(2).mean=init_gauss(class2Indx).mean;
est_gauss(2).cov=init_gauss(class2Indx).cov;
est_gauss(2).prior=gauss(2).prior;

%save the initial estimate of the mean
init_gauss=est_gauss;


%% Set the priors over the parameters
% initialize the priors over the parameters
%set the sigma of the prior normal distribution such that
%min and max values of current biomarker are within + / - 3sigma
minX=min(data);
maxX=max(data);
sigmaMin=(init_gauss(1).mean-minX)/3;
sigmaMax=(maxX-init_gauss(1).mean)/3;
prior_sigma(1)=max([sigmaMin,sigmaMax]);

sigmaMin=(init_gauss(2).mean-minX)/3;
sigmaMax=(maxX-init_gauss(2).mean)/3;
prior_sigma(2)=max([sigmaMin,sigmaMax]);

% set the alpha and beta hyperparameters of tge gamma distribution
%such that alpha*beta=sigma and 90% of populationof the gamma
%distribution is between 0 and sigma
a(1)=50;
b(1)=init_gauss(1).cov/50;

a(2)=50;
b(2)=init_gauss(2).cov/50;


prior_thisMean(1)=normpdf(est_gauss(1).mean,init_gauss(1).mean,prior_sigma(1));
prior_thisMean(2)=normpdf(est_gauss(2).mean,init_gauss(2).mean,prior_sigma(2));
prior_thisCov(1)=gampdf(est_gauss(1).cov,a(1),b(1));
prior_thisCov(2)=gampdf(est_gauss(2).cov,a(2),b(2));


%get likelihood of current Model
log_likelihood_current=getMixGaussLogLikeWithParams(data,est_gauss,prior_thisMean,prior_thisCov,a,b)
%log_likelihood_current=sum(mixGaussLogProb(data,est_gauss));

%% Visualize thel ikelihood and the prior distributions
close all
%plot for no-event
like_class1_dist=normpdf(data(membership==1),est_gauss(1).mean,sqrt(est_gauss(1).cov));
xmu0=linspace(min(data(membership==1)),max(data(membership==1)));
prior_dist_mean_class1=normpdf(xmu0,gauss(1).mean,prior_sigma(1));
xsig0=[0.1:0.1:10];
prior_dist_sigma_class1=gampdf(xsig0,a(1),b(1));
cum_dist_prior_sigma_class1=gamcdf(xsig0,a(1),b(1));

subplot(1,4,1);plot(data(membership==1),like_class1_dist,'g.');title('Likelihood Class1');
hold on
subplot(1,4,2);plot(xmu0,prior_dist_mean_class1,'b.');title('Prior Dist over Mean')
subplot(1,4,3);plot(xsig0,prior_dist_sigma_class1,'m.');title('Prior Dist over Sigma');
subplot(1,4,4);plot(xsig0,cum_dist_prior_sigma_class1,'m.-');title('Cumulative Dist Prior Sigma');
set(gcf,'Color',[ 1 1 1]);

%plot for event
figure
like_class2_dist=normpdf(data(membership==2),est_gauss(2).mean,sqrt(est_gauss(2).cov));
xmu0=linspace(min(data(membership==2)),max(data(membership==2)));
prior_dist_mean_class2=normpdf(xmu0,gauss(2).mean,prior_sigma(2));
xsig0=[0.1:0.1:10];
prior_dist_sigma_class2=gampdf(xsig0,a(2),b(2));
cum_dist_prior_sigma_class2=gamcdf(xsig0,a(2),b(2));

subplot(1,4,1);plot(data(membership==2),like_class2_dist,'r.');title('Likelihood Class2');
hold on
subplot(1,4,2);plot(xmu0,prior_dist_mean_class2,'b.');title('Prior Dist over Mean')
subplot(1,4,3);plot(xsig0,prior_dist_sigma_class2,'m.');title('Prior Dist over Sigma');
subplot(1,4,4);plot(xsig0,cum_dist_prior_sigma_class2,'m.-');title('Cumulative Dist Prior Sigma');
set(gcf,'Color',[ 1 1 1]);


% compar the estimated dis vs true dist
figure
p_class1 = like_class1_dist;
p_class2 = like_class2_dist;
[hist_class1, x_class1] = ksdensity(data(membership==1));
[hist_class2, x_class2] = ksdensity(data(membership==2));

subplot(121), hold on
plot(x_class1, hist_class1, 'g');
plot(x_class2, hist_class2, 'r');

subplot(122), hold on
plot(data(membership==1), p_class1, 'g.')
plot(data(membership==2), p_class2, 'r.');
set(gcf,'Color',[ 1 1 1]);

%% Start MCMC
log_likelihood_max = log_likelihood_current;
max_gauss=est_gauss;

change_step=0.01;
min_sigma=0.09;
change_dir=[1,-1];

% calculate a proposal distribution in this case based on the variance of thetotal data
sigma=std(data);
if sigma<min_sigma
    sigma=min_sigma;
end
proposal_dist=sigma*change_step;


new_gauss=est_gauss;

for it_mcmc = 1:(nr_it_burnin + nr_it_mcmc)
    
    %Vary the means and the variances
    rand_dir=randperm(2);
    new_gauss(1).mean=new_gauss(1).mean+proposal_dist*randn;
    new_gauss(2).mean=new_gauss(2).mean+ proposal_dist*randn;
    new_gauss(1).cov=new_gauss(1).cov+ proposal_dist*randn;
    new_gauss(2).cov=new_gauss(2).cov+ proposal_dist*randn;
    
    prior_thisMean(1)=normpdf(new_gauss(1).mean,init_gauss(1).mean,prior_sigma(1));
    prior_thisMean(2)=normpdf(new_gauss(2).mean,init_gauss(2).mean,prior_sigma(2));
    prior_thisCov(1)=gampdf(new_gauss(1).cov,a(1),b(1));
    prior_thisCov(2)=gampdf(new_gauss(2).cov,a(2),b(2));
    
    % constrain the Mean Class1 to always be greater than Mean Class2
    if new_gauss(2).mean>new_gauss(1).mean
        prior_thisMean(1)=0;
        prior_thisMean(2)=0;
    end
    
    %Compute the new likelihood after changing the mu and sigma
    log_likelihood_new=getMixGaussLogLikeWithParams(data,new_gauss,prior_thisMean,prior_thisCov,a,b);
    %log_likelihood_new=sum(mixGaussLogProb(data,new_gauss));
    
    
    alpha = exp(log_likelihood_new - log_likelihood_current);
    if alpha > rand,
        
        est_gauss = new_gauss;
        log_likelihood_current = log_likelihood_new;
        
    end
    if log_likelihood_current > log_likelihood_max,
        
        max_gauss = est_gauss;
        log_likelihood_max = log_likelihood_current;
        
    end
    

    if it_mcmc > nr_it_burnin,
        gauss_mcmc_mean(:,it_mcmc-nr_it_burnin)=[est_gauss.mean]';
        gauss_mcmc_cov(:,it_mcmc-nr_it_burnin)=[est_gauss.cov]';
        log_likelihood_mcmc(it_mcmc-nr_it_burnin) = log_likelihood_current;
        
    end
    
    if find(it_vec == it_mcmc),
        
        if it_mcmc < nr_it_burnin,
            
            fprintf('burnin it: %d\n', it_mcmc)
            
        else
            
            fprintf('mcmc it: %d\n', it_mcmc-nr_it_burnin)
            
        end
        
    end
    
end % end MCMC
       

i

