% p_false(1) = c = probability false positives
% p_false(2) = d = probability false negatives (see Puolamaki)

function [parm_struct] = EBDPMCMCTestWithParamsFast(data_controls, data_patients, parm_mcmc)
        
[nr_pat, nr_events] = size(data_patients);
nr_gradient_ascent = parm_mcmc.nr_gradient_ascent;
nr_it_gradient_ascent = parm_mcmc.nr_it_gradient_ascent;
nr_it_burnin = parm_mcmc.nr_it_burnin;
nr_it_mcmc = parm_mcmc.nr_it_mcmc;
nr_it_check_p_false = parm_mcmc.nr_it_check_p_false;
range_p_false = parm_mcmc.range_p_false;
std_p_false = parm_mcmc.std_p_false;
flag_sum = parm_mcmc.flag_sum;

interval_display = parm_mcmc.interval_display;
it_vec = 1:interval_display:(nr_it_mcmc + nr_it_burnin);
data_tot=[data_controls;data_patients];
[nr_subj, nr_events]=size(data_tot);

% Get an initial estimate of the means and variances from the data
 version_likelihood = 13;  %version_likelihood==13 is the same as 8 but with a different GMM package
 nAttempts=1;
 threshold_flag=0;
    [likelihood_events, gmix_struct_est] = ...
        EBDPComputeLikelihood(data_patients, data_controls, version_likelihood,threshold_flag,nAttempts);

% Evaluate the priors for current Mu and sigma
for roi=1:nr_events
    thisMean_noevent(roi)=gmix_struct_est.gmix_controls{roi}.mean;
    thisMean_event(roi)=gmix_struct_est.gmix_patients{roi}.mean;
    thisCov_noevent(roi)=gmix_struct_est.gmix_controls{roi}.cov;
    thisCov_event(roi)=gmix_struct_est.gmix_patients{roi}.cov;
    
    %set the hyperparameters of the priors over Mu and Sigma
    %set the prior normal distribution centered over current mean
    init_mean_noevent(roi)=thisMean_noevent(roi);
    init_mean_event(roi)=thisMean_event(roi);
    
    %set the sigma of the prior normal distribution such that
    %min and max values of current biomarker are within + / - 3sigma
    minX=min(data_tot(:,roi));
    maxX=max(data_tot(:,roi));
    sigmaMin=(init_mean_noevent(roi)-minX)/3;
    sigmaMax=(maxX-init_mean_noevent(roi))/3; 
    prior_sigma_noevent(roi)=max([sigmaMin,sigmaMax]);
    
    sigmaMin=(init_mean_event(roi)-minX)/3;
    sigmaMax=(maxX-init_mean_event(roi))/3; 
    prior_sigma_event(roi)=max([sigmaMin,sigmaMax]);
    
    % set the alpha and beta hyperparameters of tge gamma distribution
    %such that alpha*beta=sigma and 90% of populationof the gamma
    %distribution is between 0 and sigma
    a_noevent(roi)=50;
    b_noevent(roi)=thisCov_noevent(roi)/50;
    
    a_event(roi)=50;
    b_event(roi)=thisCov_event(roi)/50;
        
    
    prior_thisMean_noevent(roi)=normpdf(thisMean_noevent(roi),init_mean_noevent(roi),prior_sigma_noevent(roi));
    prior_thisMean_event(roi)=normpdf(thisMean_event(roi),init_mean_event(roi),prior_sigma_event(roi));
    prior_thisCov_noevent(roi)=gampdf(thisCov_noevent(roi),a_noevent(roi),b_noevent(roi));
    prior_thisCov_event(roi)=gampdf(thisCov_event(roi),a_event(roi),b_event(roi));
    
    
    %% Visualize thel ikelihood and the prior distributions
    close all
    %plot for no-event
    like_noevent_dist=normpdf(data_controls(:,roi),thisMean_noevent(roi),sqrt(thisCov_noevent(roi)));
    xmu0=linspace(min(data_controls(:,roi)),max(data_controls(:,roi)));
    prior_dist_mean_noevent=normpdf(xmu0,init_mean_noevent(roi),prior_sigma_noevent(roi));
    xsig0=[0.1:0.1:10];
    prior_dist_sigma_noevent=gampdf(xsig0,a_noevent(roi),b_noevent(roi));
    cum_dist_prior_sigma_noevent=gamcdf(xsig0,a_noevent(roi),b_noevent(roi));
    
    subplot(1,4,1);plot(data_controls(:,roi),like_noevent_dist,'g.');title('Likelihood Dist No-event');
    hold on
    subplot(1,4,2);plot(xmu0,prior_dist_mean_noevent,'b.');title('Prior Dist over Mean')
    subplot(1,4,3);plot(xsig0,prior_dist_sigma_noevent,'m.');title('Prior Dist over Sigma');
    subplot(1,4,4);plot(xsig0,cum_dist_prior_sigma_noevent,'m.-');title('Cumulative Dist Prior Sigma');
    set(gcf,'Color',[ 1 1 1]);
    
    %plot for event
    figure
    like_event_dist=normpdf(data_patients(:,roi),thisMean_event(roi),sqrt(thisCov_event(roi)));
    xmu0=linspace(min(data_patients(:,roi)),max(data_patients(:,roi)));
    prior_dist_mean_event=normpdf(xmu0,init_mean_event(roi),prior_sigma_event(roi));
    xsig0=[0.1:0.1:10];
    prior_dist_sigma_event=gampdf(xsig0,a_event(roi),b_event(roi));
    cum_dist_prior_sigma_event=gamcdf(xsig0,a_event(roi),b_event(roi));
    
    subplot(1,4,1);plot(data_patients(:,roi),like_event_dist,'r.');title('Likelihood Dist Event');
    hold on
    subplot(1,4,2);plot(xmu0,prior_dist_mean_event,'b.');title('Prior Dist over Mean')
    subplot(1,4,3);plot(xsig0,prior_dist_sigma_event,'m.');title('Prior Dist over Sigma');
    subplot(1,4,4);plot(xsig0,cum_dist_prior_sigma_event,'m.-');title('Cumulative Dist Prior Sigma');
    set(gcf,'Color',[ 1 1 1]);
    
    
    % compar the estimated dis vs true dist
    figure
    p_controls = like_noevent_dist;
    p_patients = like_event_dist;
    [hist_c, x_c] = ksdensity(data_controls(:, roi));
    [hist_p, x_p] = ksdensity(data_patients(:, roi));
    
    subplot(121), hold on
    plot(x_c, hist_c, 'g');
    plot(x_p, hist_p, 'r');
    
    subplot(122), hold on
    plot(data_controls(:, roi), p_controls, 'g.')
    plot(data_patients(:, roi), p_patients, 'r.');
    set(gcf,'Color',[ 1 1 1]);
    
end
    
    
flag_sum=2;
% First gradient ascent
log_likelihood_gradient_ascent = ...
    zeros(nr_it_gradient_ascent, nr_gradient_ascent);
event_order_gradient_ascent = zeros(nr_events, nr_gradient_ascent);
p_false_gradient_ascent = zeros(2, nr_gradient_ascent);
for it_gradient_ascent = 1:nr_gradient_ascent,
    
    event_order_current = randperm(nr_events);
    
    log_likelihood_current=compute_loglikelihood( ...
        likelihood_events,event_order_current,flag_sum);

    for it = 1:nr_it_gradient_ascent,
        
        event_swap = randperm(nr_events);
        event_swap = event_swap(1:2);

        event_order_new = event_order_current;
        event_order_new(event_swap(1)) = event_order_current(event_swap(2));
        event_order_new(event_swap(2)) = event_order_current(event_swap(1));
        
        
%         log_likelihood_new =findLikelihoodCurrentOrdering( ...
%             likelihood_events,event_order_new);
        
        log_likelihood_new=compute_loglikelihood(likelihood_events, ...
            event_order_new, flag_sum);

        if log_likelihood_new > log_likelihood_current,
            
            event_order_current = event_order_new;
            log_likelihood_current = log_likelihood_new;
            
        end

        log_likelihood_gradient_ascent(it, it_gradient_ascent) = ...
            log_likelihood_current;
        
    end
    event_order_gradient_ascent(:, it_gradient_ascent) = ...
        event_order_current;
    
end
it_gradient_ascent_max = find(log_likelihood_gradient_ascent(end, :) == ...
    max(log_likelihood_gradient_ascent(end, :)));
it_gradient_ascent_max = it_gradient_ascent_max(1);

% Now MCMC
event_order_mcmc = zeros(nr_events, nr_it_mcmc);
p_false_mcmc = zeros(2, nr_it_mcmc);
log_likelihood_mcmc = zeros(nr_it_mcmc, 1);

event_order_current = event_order_gradient_ascent(:, it_gradient_ascent_max);
event_order_max = event_order_gradient_ascent(:, it_gradient_ascent_max);
p_false_current  = p_false_gradient_ascent(:, it_gradient_ascent_max);
p_false_max = p_false_current;
  
%Note I replace the log_likelihood_current of the gradient ascent with that
%of when also including the prior over parameters
% Also for this I need to find the likelihood over all of the data not just
% patients hence:
for roi=1:nr_events
    %find the current likelihood event no event
    likelihood_events(roi, 1:nr_subj, 1) = getGaussProb(data_tot(:, roi)',thisMean_noevent(roi),sqrt(thisCov_noevent(roi)));
    likelihood_events(roi, 1:nr_subj, 2) = getGaussProb(data_tot(:, roi)',thisMean_event(roi),sqrt(thisCov_event(roi)));    
end

log_likelihood_current=compute_loglikelihood_and_paramsZeroK(likelihood_events, ...
    event_order_current, flag_sum, prior_thisMean_noevent, ...
        prior_thisMean_event,prior_thisCov_noevent,prior_thisCov_event);

% log_likelihood_current = ...
%     log_likelihood_gradient_ascent(end, it_gradient_ascent_max);

log_likelihood_max = log_likelihood_current;
change_step=0.01;
min_sigma=0.09;
change_dir=[1,-1];

% calculate a proposal distribution in this case based on the variance of thetotal data
for roi=1:nr_events
    sigma=std(data_tot(:,roi));
    if sigma<min_sigma 
        sigma=min_sigma;
    end
    proposal_dist(roi)=sigma*change_step;
end



%start MCMC
for it_mcmc = 1:(nr_it_burnin + nr_it_mcmc),
    
    % Swapping regions
    event_swap = randperm(nr_events);
    event_swap = event_swap(1:2);    
    event_order_new = event_order_current;
    event_order_new(event_swap(1)) = event_order_current(event_swap(2));
    event_order_new(event_swap(2)) = event_order_current(event_swap(1));    
    
    % varying the means and the variances
    for roi=1:nr_events
        rand_dir=randperm(2);
        thisMean_noevent(roi)=thisMean_noevent(roi)+ change_dir(rand_dir(1))*proposal_dist(roi);
        thisMean_event(roi)=thisMean_event(roi)+ change_dir(rand_dir(1))*proposal_dist(roi);
        thisCov_noevent(roi)=thisCov_noevent(roi)+ change_dir(rand_dir(1))*proposal_dist(roi);
        thisCov_event(roi)=thisCov_event(roi)+ change_dir(rand_dir(1))*proposal_dist(roi);    
        
        prior_thisMean_noevent(roi)=normpdf(thisMean_noevent(roi),init_mean_noevent(roi),prior_sigma_noevent(roi));
        prior_thisMean_event(roi)=normpdf(thisMean_event(roi),init_mean_event(roi),prior_sigma_event(roi));
        prior_thisCov_noevent(roi)=gampdf(thisCov_noevent(roi),a_noevent(roi),b_noevent(roi));
        prior_thisCov_event(roi)=gampdf(thisCov_event(roi),a_event(roi),b_event(roi));
        
        
        % constrain the Mean Event to always be smaller than Mean No-event
        if thisMean_event(roi)>thisMean_noevent(roi)
            prior_thisMean_noevent(roi)=0;
            prior_thisMean_event(roi)=0;
        end
        
        
        %find the current likelihood event no event
        likelihood_events(roi, 1:nr_subj, 1) = getGaussProb(data_tot(:, roi)',thisMean_noevent(roi),sqrt(thisCov_noevent(roi)));
        likelihood_events(roi, 1:nr_subj, 2) = getGaussProb(data_tot(:, roi)',thisMean_event(roi),sqrt(thisCov_event(roi)));
        
    end
    
    %Compute the new likelihood after swapping some events and changing the
    %mu and sigma

    log_likelihood_new=compute_loglikelihood_and_paramsZeroK(likelihood_events, ...
    event_order_new, flag_sum, prior_thisMean_noevent, ...
        prior_thisMean_event,prior_thisCov_noevent,prior_thisCov_event);
    
    alpha = exp(log_likelihood_new - log_likelihood_current);
    if alpha > rand,
        
        event_order_current = event_order_new;
        log_likelihood_current = log_likelihood_new;
        
    end
    if log_likelihood_current > log_likelihood_max,
        
        event_order_max = event_order_current;
        log_likelihood_max = log_likelihood_current;
        
    end
    

    if it_mcmc > nr_it_burnin,
        
        event_order_mcmc(:, it_mcmc-nr_it_burnin) = event_order_current;
        log_likelihood_mcmc(it_mcmc-nr_it_burnin) = log_likelihood_current;
        
    end
    
    if find(it_vec == it_mcmc),
        
        if it_mcmc < nr_it_burnin,
            
            fprintf('burnin it: %d\n', it_mcmc)
            
        else
            
            fprintf('mcmc it: %d\n', it_mcmc-nr_it_burnin)
            
        end
        
    end
    
end % end MCMC iter

% [log_likelihood_dummy, class_pat] = compute_loglikelihood(likelihood_events, ...
%     event_order_max, p_false_max, flag_sum);
% if nargin == 3,
%     
%     [log_likelihood_dummy, class_pat_test] = compute_loglikelihood(likelihood_events_test, ...
%         event_order_max, p_false_max, flag_sum);
%     
% end

parm_struct.log_likelihood_gradient_ascent = log_likelihood_gradient_ascent;
parm_struct.event_order_gradient_ascent = event_order_gradient_ascent;
parm_struct.event_order_mcmc = event_order_mcmc;
parm_struct.log_likelihood_mcmc = log_likelihood_mcmc;
parm_struct.event_order_max = event_order_max;
parm_struct.log_likelihood_max = log_likelihood_max;
%parm_struct.class_pat = class_pat;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [log_likelihood_tot, class_pat] = compute_loglikelihood(likelihood_events, ...
    event_order, flag_sum)

[nr_events, nr_pat] = size(likelihood_events(:, :, 1));

likelihood_nonclinevents = likelihood_events;
likelihood_nonclinevents(likelihood_nonclinevents == 0) = ...
    min(likelihood_nonclinevents(likelihood_nonclinevents ~= 0));
log_likelihood_nonclinevents = log(likelihood_nonclinevents);

nr_phase = nr_events;
event_pattern = zeros(nr_events, nr_events);
for event = 1:nr_events,
    
    event_pattern(event_order==event, event:end) = 1;
    
end

events = event_pattern == 1;
events_nonclin = events;
noevents = event_pattern == 0;
noevents_nonclin = noevents;
log_likelihood = zeros(nr_phase, nr_pat);
for phase = 1:nr_phase,
   
    log_likelihood(phase, :) = ...
        sum(log_likelihood_nonclinevents(events_nonclin(:, phase), :, 2), 1) + ...
        sum(log_likelihood_nonclinevents(noevents_nonclin(:, phase), :, 1), 1);
        
end

if flag_sum == 1,
    
    sum_log_likelihood = 0;
    for pat = 1:nr_pat,
        
        I_pat = find(log_likelihood(:, pat) == max(log_likelihood(:, pat)));
        sum_log_likelihood = sum_log_likelihood + log_likelihood(I_pat(1), pat);
        
    end
    log_likelihood_tot = sum_log_likelihood;
    
elseif flag_sum == 2,
    
    log_likelihood_tot = sum(log(sum(exp(log_likelihood), 1)));
    
end

if nargout == 2,
    
    class_pat = zeros(nr_pat, 1);
    for pat = 1:nr_pat,
        
        I_pat = find(log_likelihood(:, pat) == max(log_likelihood(:, pat)));
        class_pat(pat) = I_pat(1);
        
    end
    
end



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [log_likelihood_tot, class_pat] = compute_loglikelihood_and_params(likelihood_events, ...
    event_order, flag_sum, prior_thisMean_noevent, ...
        prior_thisMean_event,prior_thisCov_noevent,prior_thisCov_event)

[nr_events, nr_pat] = size(likelihood_events(:, :, 1));

likelihood_nonclinevents = likelihood_events;
likelihood_nonclinevents(likelihood_nonclinevents == 0) = ...
    min(likelihood_nonclinevents(likelihood_nonclinevents ~= 0));
log_likelihood_nonclinevents = log(likelihood_nonclinevents);

nr_phase = nr_events;
event_pattern = zeros(nr_events, nr_events);
for event = 1:nr_events,
    
    event_pattern(event_order==event, event:end) = 1;
    
end

events = event_pattern == 1;
events_nonclin = events;
noevents = event_pattern == 0;
noevents_nonclin = noevents;
log_likelihood = zeros(nr_phase, nr_pat);
for phase = 1:nr_phase,
   

    logPrEvent=log_likelihood_nonclinevents(events_nonclin(:, phase), :, 2);
    logPrMuEvent=log(prior_thisMean_event(events_nonclin(:, phase)'));
    logPrCovEvent=log(prior_thisCov_event(events_nonclin(:, phase)'));
    
    logPrNoEvent=log_likelihood_nonclinevents(noevents_nonclin(:, phase), :, 1);
    logPrMuNoEvent=log(prior_thisMean_noevent(noevents_nonclin(:, phase)'));
    logPrCovNoEvent=log(prior_thisCov_noevent(noevents_nonclin(:, phase)'));
    
    log_likelihood(phase, :)=sum(logPrEvent+...
        repmat(logPrMuEvent',1,size(logPrEvent,2))+ ...
        repmat(logPrCovEvent',1,size(logPrEvent,2)),1)+...
        sum(logPrNoEvent+repmat(logPrMuNoEvent', 1, size(logPrNoEvent,2))+ ...
        repmat(logPrCovNoEvent',1,size(logPrNoEvent,2)),1);
        
end

if flag_sum == 1,
    
    sum_log_likelihood = 0;
    for pat = 1:nr_pat,
        
        I_pat = find(log_likelihood(:, pat) == max(log_likelihood(:, pat)));
        sum_log_likelihood = sum_log_likelihood + log_likelihood(I_pat(1), pat);
        
    end
    log_likelihood_tot = sum_log_likelihood;
    
elseif flag_sum == 2,
    
    log_likelihood_tot = sum(log(sum(exp(log_likelihood), 1)));
    
end

if nargout == 2,
    
    class_pat = zeros(nr_pat, 1);
    for pat = 1:nr_pat,
        
        I_pat = find(log_likelihood(:, pat) == max(log_likelihood(:, pat)));
        class_pat(pat) = I_pat(1);
        
    end
    
end