function parm_struct=findBestOrderingAndParams(data_controls,data_patients,parm_mcmc)

%% Initialization
data_tot=[data_controls;data_patients];
[nr_subj,nr_events]=size(data_tot);
nr_it_burnin = parm_mcmc.nr_it_burnin;
nr_it_mcmc = parm_mcmc.nr_it_mcmc;

%% Get an initial estimate of the means and variances from the data
 version_likelihood = 13;  %version_likelihood==13 is the same as 8 but with a different GMM package
 nAttempts=5;
 threshold_flag=0;
    [likelihood, gmix_struct_est] = ...
        EBDPComputeLikelihood(data_patients, data_controls, version_likelihood,threshold_flag,nAttempts);

%% Get an initial ordering with this means and variances
[event_order_current]=findNEventOrdering(likelihood);

%% Evaluate the priors for current Mu and sigma
for roi=1:nr_events
    thisMean_control(roi)=gmix_struct_est.gmix_controls{roi}.mean;
    thisMean_patient(roi)=gmix_struct_est.gmix_patients{roi}.mean;
    thisCov_control(roi)=gmix_struct_est.gmix_controls{roi}.cov;
    thisCov_patient(roi)=gmix_struct_est.gmix_patients{roi}.cov;
    
    prior_thisMean_control(roi)=normpdf(thisMean_control_init,0,1);
    prior_thisMean_patient(roi)=normpdf(thisMean_patient_init,0,1);
    prior_thisCov_control(roi)=gampdf(thisCov_control_init,1,2);
    prior_thisCov_patient(roi)=gampdf(thisCov_patient_init,1,2);
    
end
%find the ikelihood of current ordering
log_likelihood_current=findLikelihoodCurrentOrderingAndParams(likelihood,event_order_current,prior_thisMean_control,prior_thisMean_patient,prior_thisCov_control,prior_thisCov_patient);
    

event_order_max=event_order_current;
log_likelihood_max = log_likelihood_current;

event_order_mcmc = zeros(nr_events, nr_it_mcmc);
log_likelihood_mcmc = zeros(nr_it_mcmc, 1);
change_step=0.01;
change_dir=[1,-1];

for it_mcmc=1:nr_it_mcmc
    it_mcmc
    %swap the order of two events
    swap=randperm(nr_events);
    event_order_new=event_order_current; 
    event_order_new(swap(1))=event_order_current(swap(2));
    event_order_new(swap(2))=event_order_current(swap(1));  
    
    %vary the means and the variances 
    for roi=1:nr_events
        rand_dir=randperm(2);
        thisMean_control(roi)=thisMean_control(roi)+ change_dir(rand_dir(1))*thisMean_control(roi)*change_step;
        thisMean_patient(roi)=thisMean_patient(roi)+ thisMean_patient(roi)*change_step;
        thisCov_control(roi)=thisCov_control(roi)+ change_dir(rand_dir(1))*thisCov_control(roi)*change_step;
        thisCov_patient(roi)=thisCov_patient(roi)+ thisCov_patient(roi)*change_step;
        prior_thisMean_control(roi)=normpdf(thisMean_control(roi),0,1);
        prior_thisMean_patient(roi)=normpdf(thisMean_patient(roi),0,1);
        prior_thisCov_control(roi)=gampdf(thisCov_control(roi),1,2);
        prior_thisCov_patient(roi)=gampdf(thisCov_patient(roi),1,2);
        
        
        %find the current likelihood event no event
        new_likelihood(roi, :, 1) = getGaussProb(data_tot(:, roi)',thisMean_control(roi),sqrt(thisCov_control(roi)));
        new_likelihood(roi, :, 2) = getGaussProb(data_tot(:, roi)',thisMean_patient(roi),sqrt(thisCov_patient(roi))); 
        
    end
    
 
    %find likelihood of the new ordering ordering
    log_likelihood_new=findLikelihoodCurrentOrderingAndParams(new_likelihood,event_order_new,prior_thisMean_control,prior_thisMean_patient,prior_thisCov_control,prior_thisCov_patient);
    
    alpha = exp(log_likelihood_new - log_likelihood_current);
    if alpha > rand,
        
        event_order_current = event_order_new;
        log_likelihood_current = log_likelihood_new;
        
    end
    if log_likelihood_current > log_likelihood_max,
        
        event_order_max = event_order_current;
        log_likelihood_max = log_likelihood_current;
        
    end
    %*********** NOTE I HAVE NOT ADDED BURNING YET *************
    event_order_mcmc(:, it_mcmc) = event_order_current;
    log_likelihood_mcmc(it_mcmc) = log_likelihood_current;
    
    
end


parm_struct.event_order_mcmc = event_order_mcmc;
parm_struct.log_likelihood_mcmc = log_likelihood_mcmc;


