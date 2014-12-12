function [parm_struct] = VerifyTestWithParamsFast(data_controls, data_patients, parm_mcmc)

[nr_events, nr_pat] = size(data_patients);
nr_gradient_ascent = parm_mcmc.nr_gradient_ascent;
nr_it_gradient_ascent = parm_mcmc.nr_it_gradient_ascent;
nr_it_burnin = parm_mcmc.nr_it_burnin;
nr_it_mcmc = parm_mcmc.nr_it_mcmc;
nr_it_check_p_false = parm_mcmc.nr_it_check_p_false;
range_p_false = parm_mcmc.range_p_false;
std_p_false = parm_mcmc.std_p_false;
flag_sum = parm_mcmc.flag_sum;

thisMean_noevent_mcmc=zeros(nr_it_mcmc,1);
thisMean_event_mcmc=zeros(nr_it_mcmc,1);
thisCov_noevent_mcmc=zeros(nr_it_mcmc,1);
thisCov_event_mcmc=zeros(nr_it_mcmc,1);

thisMean_noevent_burn=zeros(nr_it_burnin,1);
thisMean_event_burn=zeros(nr_it_burnin,1);
thisCov_noevent_burn=zeros(nr_it_burnin,1);
thisCov_event_burn=zeros(nr_it_burnin,1);


pMNE=zeros(nr_it_mcmc,1);
pME=zeros(nr_it_mcmc,1);
pCNE=zeros(nr_it_mcmc,1);
pCE=zeros(nr_it_mcmc,1);



interval_display = parm_mcmc.interval_display;
it_vec = 1:interval_display:(nr_it_mcmc + nr_it_burnin);
data_tot=[data_controls data_patients];
[nr_events, nr_subj]=size(data_tot);

%Get an initial estimate of the means and variances from the data
opts = foptions;
opts(5) = 1;
opts(14) = 1e3;
opts(19) = 10; % nr_components_add_max
opts(20) = 0.6; % reduction_covars
opts(21) = 1e2; % #repetitions for adding components
gmix_struct_est = gmm(1, 2, 'full');
gmix_struct_est = gmmem(gmix_struct_est, [data_controls'; data_patients'], opts);


% Evaluate the priors for current Mu and sigma
for roi=1:nr_events
    %% set means using gmem
    thisMean_noevent=gmix_struct_est.centres(1);
    thisMean_event=gmix_struct_est.centres(2);
    thisCov_noevent=gmix_struct_est.covars(1);
    thisCov_event=gmix_struct_est.covars(2);
    
    
    %% set the hyperparameters of the priors over Mu and Sigma
    
    %set the prior normal distribution centered over current mean
    init_mean_noevent=thisMean_noevent;
    init_mean_event=thisMean_event;
    init_cov_noevent=thisCov_noevent;
    init_cov_event=thisCov_event;
    
    %     %% **** Option1: set the sima for both event and noevent based on minimum and maximum
    %     %  ****  value of the entire data:  THIS IS NOT GOOD BECAUSE IT CAUSED THE
    %     %  **** NO-EVENT DISTRIBUTION TO HAVE A LARGER VARIANCE THAN THAT OF EVENT WHICH
    %     %  **** IS NOT WHAT WE ASSUME IN THE MODEL
    %
    %         %set the sigma of the prior normal distribution such that
    %         %min and max values of current biomarker are within + / - 3sigma
    %         minX=min(data_tot);
    %         maxX=max(data_tot);
    %         sigmaMin=(init_mean_noevent-minX)/3;
    %         sigmaMax=(maxX-init_mean_noevent)/3;
    %         prior_sigma_noevent=max([sigmaMin,sigmaMax]);
    %
    %         sigmaMin=(init_mean_event-minX)/3;
    %         sigmaMax=(maxX-init_mean_event)/3;
    %         prior_sigma_event=max([sigmaMin,sigmaMax]);
    %
    
    
    %% **** Option2: set the sigma for event and noevent based on minimum and maximum
    %  ****  value of the data_atient and data_control distributions separately
    
    %set the sigma of the prior normal distribution such that
    %min and max values of current biomarker are within + / - 3sigma
    minC=min(data_controls);
    maxC=max(data_controls);
    sigmaMin=(init_mean_noevent-minC)/3;
    sigmaMax=(maxC-init_mean_noevent)/3;
    prior_sigma_noevent=max([sigmaMin,sigmaMax]);
    
    minP=min(data_patients);
    maxP=max(data_patients);
    sigmaMin=(init_mean_event-minP)/3;
    sigmaMax=(maxP-init_mean_event)/3;
    prior_sigma_event=max([sigmaMin,sigmaMax]);
    
    
    %% SET THE ALPHA AND BETA HYPERPARAMETERS OF THE PRIOR DISTRIBUTION
    % set the alpha and beta hyperparameters of tge gamma distribution
    %such that the mean i.e. alpha*beta=sigma and the standard deviation is
    %10 times the mean
    
    b_noevent=100;
    a_noevent=thisCov_noevent/100;
    
    b_event=100;
    a_event=thisCov_event/100;
    
    
    prior_thisMean_noevent=normpdf(thisMean_noevent,init_mean_noevent,prior_sigma_noevent);
    prior_thisMean_event=normpdf(thisMean_event,init_mean_event,prior_sigma_event);
    prior_thisCov_noevent=gampdf(thisCov_noevent,a_noevent,b_noevent);
    prior_thisCov_event=gampdf(thisCov_event,a_event,b_event);
    
    
end

% Now MCMC
log_likelihood_mcmc = zeros(nr_it_mcmc, 1);
log_likelihood_burn=zeros(nr_it_burnin,1);


weight_noevent=0.25;
weight_event=0.75;

%find the current likelihood event no event
likelihood_events(:, 1) = getGaussProb(data_tot,thisMean_noevent,sqrt(thisCov_noevent));
likelihood_events(:, 2) = getGaussProb(data_tot,thisMean_event,sqrt(thisCov_event));


likelihood_data= likelihood_events(:, 1).*weight_noevent*prior_thisMean_noevent*prior_thisCov_noevent+ ...
    likelihood_events(:,2).*weight_event*prior_thisMean_event*prior_thisCov_event;
log_likelihood_current=sum(log(likelihood_data));


log_likelihood_max = log_likelihood_current;
change_step_muNE=0.0004;
change_step_muE=0.0009;
change_step_sigNE=0.0003;
change_step_sigE=0.0001;
min_sigma=0.09;

% calculate a proposal distribution in this case based on the variance of thetotal data
sigma=std(data_tot).^2;
if sigma<min_sigma
    sigma=min_sigma;
end

thisMean_noevent_current=thisMean_noevent;
thisMean_event_current=thisMean_event;
thisCov_noevent_current= thisCov_noevent;
thisCov_event_current= thisCov_event;

thisMean_noevent_new=thisMean_noevent_current;
thisMean_event_new=thisMean_event_current;
thisCov_noevent_new= thisCov_noevent_current;
thisCov_event_new= thisCov_event_current;


%start MCMC
for it_mcmc = 1:(nr_it_burnin + nr_it_mcmc),
    
    % varying the means and the variances
    for roi=1:nr_events
        noise=randn;
        thisMean_noevent_new=thisMean_noevent_new+sigma*change_step_muNE*noise;
        noise=randn;
        thisMean_event_new=thisMean_event_new+ sigma*change_step_muE*noise;
        noise=randn;
        thisCov_noevent_new=thisCov_noevent_new+ sigma*change_step_sigNE*noise;
        noise=randn;
        thisCov_event_new=thisCov_event_new+ sigma*change_step_sigE*noise;
        
        prior_thisMean_noevent_new=normpdf(thisMean_noevent_new,init_mean_noevent,prior_sigma_noevent);
        prior_thisMean_event_new=normpdf(thisMean_event_new,init_mean_event,prior_sigma_event);
        prior_thisCov_noevent_new=gampdf(thisCov_noevent_new,a_noevent,b_noevent);
        prior_thisCov_event_new=gampdf(thisCov_event_new,a_event,b_event);
        
        
        % constrain the Mean Event to always be smaller than Mean No-event
        if thisMean_event_new>thisMean_noevent_new
            prior_thisMean_noevent_new=0;
            prior_thisMean_event_new=0;
        end
        
        
        %find the current likelihood event no event
        likelihood_events_new(:,1) = getGaussProb(data_tot,thisMean_noevent_new,sqrt(thisCov_noevent_new));
        likelihood_events_new(:,2) = getGaussProb(data_tot,thisMean_event_new,sqrt(thisCov_event_new));
        
        
    end
    
    %Compute the new likelihood after swapping some events and changing the
    %mu and sigma  
    
    likelihood_data_new= likelihood_events_new(:, 1).*weight_noevent*prior_thisMean_noevent_new*prior_thisCov_noevent_new+ ...
        likelihood_events_new(:,2).*weight_event*prior_thisMean_event_new*prior_thisCov_event_new;
    log_likelihood_new=sum(log(likelihood_data_new));
    
    
    
    alpha= exp(log_likelihood_new - log_likelihood_current);
    u=rand;
    
    if alpha> u
        
        thisMean_noevent_current=thisMean_noevent_new;
        thisMean_event_current=thisMean_event_new;
        thisCov_noevent_current= thisCov_noevent_new;
        thisCov_event_current= thisCov_event_new;
        log_likelihood_current = log_likelihood_new;
        
    end
    if log_likelihood_current > log_likelihood_max,
        
        thisMean_noevent_max=thisMean_noevent_current;
        thisMean_event_max=thisMean_event_current;
        thisCov_noevent_max= thisCov_noevent_current;
        thisCov_event_max= thisCov_event_current;
        log_likelihood_max = log_likelihood_current;
        
    end
    
    
    if it_mcmc > nr_it_burnin,
        
        
        thisMean_noevent_mcmc(it_mcmc-nr_it_burnin)=thisMean_noevent_current;
        thisMean_event_mcmc(it_mcmc-nr_it_burnin)=thisMean_event_current;
        thisCov_noevent_mcmc(it_mcmc-nr_it_burnin)= thisCov_noevent_current;
        thisCov_event_mcmc(it_mcmc-nr_it_burnin)= thisCov_event_current;
        log_likelihood_mcmc(it_mcmc-nr_it_burnin) = log_likelihood_current;
    else
        
        thisMean_noevent_burn(it_mcmc)=thisMean_noevent_current;
        thisMean_event_burn(it_mcmc)=thisMean_event_current;
        thisCov_noevent_burn(it_mcmc)= thisCov_noevent_current;
        thisCov_event_burn(it_mcmc)= thisCov_event_current;
        log_likelihood_burn(it_mcmc) = log_likelihood_current;
    end
    
    if find(it_vec == it_mcmc),
        
        if it_mcmc < nr_it_burnin,
            
            fprintf('burnin it: %d\n', it_mcmc)
            
        else
            
            fprintf('mcmc it: %d\n', it_mcmc-nr_it_burnin)
            
        end
        
    end
    
    
end % end MCMC iter


parm_struct.log_likelihood_mcmc = log_likelihood_mcmc;
%parm_struct.p_false_max = p_false_max;
parm_struct.log_likelihood_max = log_likelihood_max;
%parm_struct.class_pat = class_pat;
parm_struct.thisMean_noevent_mcmc=thisMean_noevent_mcmc;
parm_struct.thisMean_event_mcmc=thisMean_event_mcmc;
parm_struct.thisCov_noevent_mcmc=thisCov_noevent_mcmc;
parm_struct.thisCov_event_mcmc=thisCov_event_mcmc;

parm_struct.log_likelihood_burn = log_likelihood_burn;
parm_struct.thisMean_noevent_burn=thisMean_noevent_burn;
parm_struct.thisMean_event_burn=thisMean_event_burn;
parm_struct.thisCov_noevent_burn=thisCov_noevent_burn;
parm_struct.thisCov_event_burn=thisCov_event_burn;

save(['verifyDataSimulationMCMC/VerifyMCMCResultWith' num2str(10) '-' num2str(log10(nr_it_mcmc)) 'Itr' num2str(change_step) 'Steps.mat'],'parm_struct');
i
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