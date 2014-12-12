function [log_likelihood_tot, class_pat] = compute_loglikelihood_Jania(likelihood_events, ...
    event_order, flag_sum,prior_thisMean_noevent, ...
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
   
    logPrMuEvent=log(prior_thisMean_event(events_nonclin(:, phase)'));
    logPrCovEvent=log(prior_thisCov_event(events_nonclin(:, phase)'));
   
    logPrMuNoEvent=log(prior_thisMean_noevent(noevents_nonclin(:, phase)'));
    logPrCovNoEvent=log(prior_thisCov_noevent(noevents_nonclin(:, phase)'));

    logPriorMu(phase,:)=[logPrMuEvent logPrMuNoEvent];
    logPriorCov(phase,:)=[logPrCovEvent logPrCovNoEvent];

    logPrEvent=log_likelihood_nonclinevents(events_nonclin(:, phase), :, 2);
    
    logPrNoEvent=log_likelihood_nonclinevents(noevents_nonclin(:, phase), :, 1);
    
    log_likelihood(phase, :)=sum(logPrEvent,1)+sum(logPrNoEvent,1);
        
end


    
if flag_sum == 2,
    
    log_likelihood_tot = sum(log(sum(exp(log_likelihood), 1))) + sum(logPriorMu(:))+sum(logPriorCov(:));
    
end

if nargout == 2,
    
    class_pat = zeros(nr_pat, 1);
    for pat = 1:nr_pat,
        
        I_pat = find(log_likelihood(:, pat) == max(log_likelihood(:, pat)));
        class_pat(pat) = I_pat(1);
        
    end
    
end