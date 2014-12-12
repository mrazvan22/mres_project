% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [log_likelihood_tot, class_pat] = compute_loglikelihood_and_paramsZeroK(likelihood_events, ...
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

%% add a stage zero to account for controls
event_pattern=[zeros(nr_events,1), event_pattern];

events = event_pattern == 1;
events_nonclin = events;
noevents = event_pattern == 0;
noevents_nonclin = noevents;
log_likelihood = zeros(nr_phase, nr_pat);
for phase = 1:nr_phase+1,
   

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