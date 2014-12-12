% p_false(1) = c = probability false positives
% p_false(2) = d = probability false negatives (see Puolamaki)

function [parm_struct] = EBDPMCMCTest(likelihood_events, parm_mcmc, ...
    likelihood_events_test)
        
[nr_events, nr_pat] = size(likelihood_events);
nr_gradient_ascent = parm_mcmc.nr_gradient_ascent;
nr_it_gradient_ascent = parm_mcmc.nr_it_gradient_ascent;
nr_it_burnin = parm_mcmc.nr_it_burnin;
nr_it_mcmc = parm_mcmc.nr_it_mcmc;
nr_it_check_p_false = parm_mcmc.nr_it_check_p_false;
range_p_false = parm_mcmc.range_p_false;
std_p_false = parm_mcmc.std_p_false;
idx_clinevent = parm_mcmc.idx_clinevent;
flag_sum = parm_mcmc.flag_sum;

interval_display = parm_mcmc.interval_display;
it_vec = 1:interval_display:(nr_it_mcmc + nr_it_burnin);

% First gradient ascent
log_likelihood_gradient_ascent = ...
    zeros(nr_it_gradient_ascent, nr_gradient_ascent);
event_order_gradient_ascent = zeros(nr_events, nr_gradient_ascent);
p_false_gradient_ascent = zeros(2, nr_gradient_ascent);
for it_gradient_ascent = 1:nr_gradient_ascent,
    
    event_order_current = randperm(nr_events);
    p_false_current = zeros(2, 1);
    for it_p = 1:2,
        
        p_false_current(it_p) = ...
            (range_p_false(it_p, 2)-range_p_false(it_p, 1))*rand + ...
            range_p_false(it_p, 1);
        
    end
    log_likelihood_current = compute_loglikelihood(likelihood_events, ...
        event_order_current, p_false_current, idx_clinevent, flag_sum);
    for it = 1:nr_it_gradient_ascent,
        
        event_swap = randperm(nr_events);
        event_swap = event_swap(1:2);

        event_order_new = event_order_current;
        event_order_new(event_swap(1)) = event_order_current(event_swap(2));
        event_order_new(event_swap(2)) = event_order_current(event_swap(1));
        
        log_likelihood_new = compute_loglikelihood(likelihood_events, ...
            event_order_new, p_false_current, idx_clinevent, flag_sum);
        if log_likelihood_new > log_likelihood_current,
            
            event_order_current = event_order_new;
            log_likelihood_current = log_likelihood_new;
            
        end

        for it_p = 1:2,
            
            p_false_new = p_false_current;
            p_false_new(it_p) = max(range_p_false(it_p, 1), ...
                min(range_p_false(it_p, 2), ...
                normrnd(p_false_current(it_p), std_p_false(it_p))));      
            
            log_likelihood_new = compute_loglikelihood(likelihood_events, ...
                event_order_current, p_false_new, idx_clinevent, flag_sum);
            if log_likelihood_new > log_likelihood_current,
                
                p_false_current = p_false_new;
                log_likelihood_current = log_likelihood_new;
                
            end
            
        end
        log_likelihood_gradient_ascent(it, it_gradient_ascent) = ...
            log_likelihood_current;
        
    end
    event_order_gradient_ascent(:, it_gradient_ascent) = ...
        event_order_current;
    p_false_gradient_ascent(:, it_gradient_ascent) = p_false_current;
    
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
log_likelihood_current = ...
    log_likelihood_gradient_ascent(end, it_gradient_ascent_max);
log_likelihood_max = log_likelihood_current;
alpha_store = zeros(2, nr_it_check_p_false);
it_alpha_store = 1;
for it_mcmc = 1:(nr_it_burnin + nr_it_mcmc),
    
    % Swapping regions
    event_swap = randperm(nr_events);
    event_swap = event_swap(1:2);    
    event_order_new = event_order_current;
    event_order_new(event_swap(1)) = event_order_current(event_swap(2));
    event_order_new(event_swap(2)) = event_order_current(event_swap(1));    
    log_likelihood_new = compute_loglikelihood(likelihood_events, ...
        event_order_new, p_false_current, idx_clinevent, flag_sum);
    alpha = exp(log_likelihood_new - log_likelihood_current);
    if alpha > rand,
        
        event_order_current = event_order_new;
        log_likelihood_current = log_likelihood_new;
        
    end
    if log_likelihood_current > log_likelihood_max,
        
        event_order_max = event_order_current;
        p_false_max = p_false_current;
        log_likelihood_max = log_likelihood_current;
        
    end
    
    % Changing probabilities false positives/negatives
    for it_p = 1:2,
        
        p_false_new = p_false_current;
        while p_false_new(it_p) == p_false_current(it_p),
            
            p_false_new(it_p) = max(range_p_false(it_p, 1), ...
                min(range_p_false(it_p, 2), ...
                normrnd(p_false_current(it_p), std_p_false(it_p))));
            
        end
        log_likelihood_new = compute_loglikelihood(likelihood_events, ...
            event_order_current, p_false_new, idx_clinevent, flag_sum);

        alpha = min(1, exp(log_likelihood_new - log_likelihood_current));
        if alpha > rand,
            
            p_false_current = p_false_new;
            log_likelihood_current = log_likelihood_new;
            
        end
    
        if log_likelihood_current > log_likelihood_max,
            
            event_order_max = event_order_current;
            p_false_max = p_false_current;
            log_likelihood_max = log_likelihood_current;
            
        end

        alpha_store(it_p, it_alpha_store) = alpha;
        if it_alpha_store == nr_it_check_p_false,
            
            mean_alpha = mean(alpha_store(it_p, :));
            if mean_alpha < 0.15,
                
                std_p_false(it_p) = std_p_false(it_p) - 0.1*std_p_false(it_p);
                
            elseif mean_alpha > 0.5,
                
                std_p_false(it_p) = std_p_false(it_p) + 0.1*std_p_false(it_p);
                
            end
            alpha_store(it_p, :) = zeros(nr_it_check_p_false, 1);
            
        else
            
            it_alpha_store = it_alpha_store + 1;
            
        end
        
    end
    if it_alpha_store == nr_it_check_p_false,
        
        it_alpha_store = 1;
        
    end

    if it_mcmc > nr_it_burnin,
        
        event_order_mcmc(:, it_mcmc-nr_it_burnin) = event_order_current;
        p_false_mcmc(:, it_mcmc-nr_it_burnin) = p_false_current;
        log_likelihood_mcmc(it_mcmc-nr_it_burnin) = log_likelihood_current;
        
    end
    
    if find(it_vec == it_mcmc),
        
        if it_mcmc < nr_it_burnin,
            
            fprintf('burnin it: %d\n', it_mcmc)
            
        else
            
            fprintf('mcmc it: %d\n', it_mcmc-nr_it_burnin)
            
        end
        
    end
    
end

[log_likelihood_dummy, class_pat] = compute_loglikelihood(likelihood_events, ...
    event_order_max, p_false_max, idx_clinevent, flag_sum);
if nargin == 3,
    
    [log_likelihood_dummy, class_pat_test] = compute_loglikelihood(likelihood_events_test, ...
        event_order_max, p_false_max, idx_clinevent, flag_sum);
    
end

parm_struct.log_likelihood_gradient_ascent = log_likelihood_gradient_ascent;
parm_struct.event_order_gradient_ascent = event_order_gradient_ascent;
parm_struct.event_order_mcmc = event_order_mcmc;
parm_struct.p_false_mcmc = p_false_mcmc;
parm_struct.log_likelihood_mcmc = log_likelihood_mcmc;
parm_struct.event_order_max = event_order_max;
parm_struct.p_false_max = p_false_max;
parm_struct.log_likelihood_max = log_likelihood_max;
parm_struct.class_pat = class_pat;
if nargin == 3,
    
    parm_struct.class_pat_test = class_pat_test;
    
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [log_likelihood_tot, class_pat] = compute_loglikelihood(likelihood_events, ...
    event_order, p_false, idx_clinevent, flag_sum)

[nr_events, nr_pat] = size(likelihood_events(:, :, 1));
nr_events_nonclin = nr_events - length(idx_clinevent);
nr_events_clin = length(idx_clinevent);
idx_nonclinevent = 1:nr_events;
idx_nonclinevent(idx_clinevent) = [];

likelihood_clinevents = likelihood_events(idx_clinevent, :, :);
likelihood_nonclinevents = likelihood_events;
likelihood_nonclinevents(idx_clinevent, :, :) = [];
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
events_nonclin(idx_clinevent, :) = [];
events_clin = events(idx_clinevent, :);
noevents = event_pattern == 0;
noevents_nonclin = noevents;
noevents_nonclin(idx_clinevent, :) = [];
noevents_clin = noevents(idx_clinevent, :);
log_likelihood = zeros(nr_phase, nr_pat);
for phase = 1:nr_phase,
    
    log_likelihood_nonclinevents_phase = ...
        sum(log_likelihood_nonclinevents(events_nonclin(:, phase), :, 2), 1) + ...
        sum(log_likelihood_nonclinevents(noevents_nonclin(:, phase), :, 1), 1);
    
    if sum(noevents_clin(:, phase)),
        
        alpha = sum(repmat(noevents_clin(:, phase), 1, nr_pat).*...
            likelihood_clinevents(:, :, 2), 1)./...
            repmat(sum(noevents_clin(:, phase)), 1, nr_pat);
        beta = sum(repmat(noevents_clin(:, phase), 1, nr_pat).*...
            likelihood_clinevents(:, :, 1), 1)./...
            repmat(sum(noevents_clin(:, phase)), 1, nr_pat);
        
    else
        
        alpha = 0;
        beta = 0;
        
    end
    if sum(events_clin(:, phase)),
        
        gamma = sum(repmat(events_clin(:, phase), 1, nr_pat).*...
            likelihood_clinevents(:, :, 1), 1)./...
            repmat(sum(events_clin(:, phase)), 1, nr_pat);
        delta = sum(repmat(events_clin(:, phase), 1, nr_pat).*...
            likelihood_clinevents(:, :, 2), 1)./...
            repmat(sum(events_clin(:, phase)), 1, nr_pat);
        
    else
        
        gamma = 0;
        delta = 0;
        
    end    
    log_likelihood_clinevents_phase = alpha*log(p_false(1)) + ...
        beta*log(1-p_false(1)) + ...
        gamma*log(p_false(2)) + ...
        delta*log(1-p_false(2));
    log_likelihood(phase, :) = log_likelihood_nonclinevents_phase + ...
        log_likelihood_clinevents_phase;
        
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


% events_clin_match = zeros(nr_events_clin, nr_pat);
% for pat = 1:nr_pat,
%     
%     match_mat = events_clin.*repmat(likelihood_clinevents(:, pat, 2), 1, nr_phase);
%     phase_pat = find(sum(match_mat, 1) == max(sum(match_mat, 1)));
%     events_clin_match(:, pat) = events_clin(:, phase_pat(1));
%     
% end
% 
% alpha = sum(sum((1 - events_clin_match).*likelihood_clinevents(:, :, 2)));
% beta = sum(sum((1 - events_clin_match).*likelihood_clinevents(:, :, 1)));
% gamma = sum(sum(events_clin_match.*likelihood_clinevents(:, :, 1)));
% delta = sum(sum(events_clin_match.*likelihood_clinevents(:, :, 2)));
% log_likelihood_clinevents = alpha*log(p_false(1)) + ...
%     beta*log(1-p_false(1)) + ...
%     gamma*log(p_false(2)) + ...
%     delta*log(1-p_false(2));
% log_likelihood_nonclinevents_tot = ...
%     sum(log(max(likelihood_nonclinevents_phase, [], 1)));
% log_likelihood = log_likelihood_nonclinevents_tot + log_likelihood_clinevents;
      
        