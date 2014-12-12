function [parm_struct] = EBDPMCMC(likelihood_events, parm_mcmc)
% parm_struct = EBDPMCMC(likelihood_events, parm_mcmc)
% The MCMC procedure consists of three phases:
% 1. A gradient ascent phase. In this phase the normal MCMC procedure is
% followed, but only steps which increase the likelihood are accepted.
% This procedure is repeated several times from different starting
% positions to get close to the maximum likelihood sequence
% 2. A burnin phase
% 3. The MCMC phase


[nr_events, nr_pat] = size(likelihood_events);
nr_gradient_ascent = parm_mcmc.nr_gradient_ascent;
nr_it_gradient_ascent = parm_mcmc.nr_it_gradient_ascent; 
nr_it_burnin = parm_mcmc.nr_it_burnin;
nr_it_mcmc = parm_mcmc.nr_it_mcmc;
nr_it_check_f0 = parm_mcmc.nr_it_check_f0;
f0_range = parm_mcmc.f0_range;
std_f0 = parm_mcmc.std_f0;
data_range = parm_mcmc.data_range;

interval_display = parm_mcmc.interval_display;
it_vec = 1:interval_display:(nr_it_mcmc + nr_it_burnin);

% First gradient ascent
log_likelihood_gradient_ascent = ...
    zeros(nr_it_gradient_ascent, nr_gradient_ascent);
event_order_gradient_ascent = zeros(nr_events, nr_gradient_ascent);
f0_gradient_ascent = zeros(nr_gradient_ascent, 1);
for it_gradient_ascent = 1:nr_gradient_ascent,
    
    event_order_current = randperm(nr_events);
    f0_current = (f0_range(2)-f0_range(1))*rand + f0_range(1);
    log_likelihood_current = compute_loglikelihood(likelihood_events, ...
        event_order_current, f0_current, data_range);
    for it = 1:nr_it_gradient_ascent,
        
        event_swap = randperm(nr_events);
        event_swap = event_swap(1:2);

        event_order_new = event_order_current;
        event_order_new(event_swap(1)) = event_order_current(event_swap(2));
        event_order_new(event_swap(2)) = event_order_current(event_swap(1));
        
        log_likelihood_new = compute_loglikelihood(likelihood_events, ...
            event_order_new, f0_current, data_range);
        if log_likelihood_new > log_likelihood_current,
            
            event_order_current = event_order_new;
            log_likelihood_current = log_likelihood_new;
            
        end

        f0_new = max(f0_range(1), min(f0_range(2), normrnd(f0_current, std_f0)));
        log_likelihood_new = compute_loglikelihood(likelihood_events, ...
            event_order_current, f0_new, data_range);        
        if log_likelihood_new > log_likelihood_current,

            f0_current = f0_new;
            log_likelihood_current = log_likelihood_new;

        end

        log_likelihood_gradient_ascent(it, it_gradient_ascent) = ...
            log_likelihood_current;
        
    end
    event_order_gradient_ascent(:, it_gradient_ascent) = ...
        event_order_current;
    f0_gradient_ascent(it_gradient_ascent) = f0_current;
    
end
it_gradient_ascent_max = find(log_likelihood_gradient_ascent(end, :) == ...
    max(log_likelihood_gradient_ascent(end, :)));
it_gradient_ascent_max = it_gradient_ascent_max(1);

% Now MCMC
event_order_mcmc = zeros(nr_events, nr_it_mcmc);
f0_mcmc = zeros(nr_it_mcmc, 1);
log_likelihood_mcmc = zeros(nr_it_mcmc, 1);

event_order_current = event_order_gradient_ascent(:, it_gradient_ascent_max);
event_order_max = event_order_gradient_ascent(:, it_gradient_ascent_max);
f0_current  = f0_gradient_ascent(it_gradient_ascent_max);
f0_max = f0_current;
log_likelihood_current = log_likelihood_gradient_ascent(end, it_gradient_ascent_max);
log_likelihood_max = log_likelihood_current;
alpha_store = zeros(nr_it_check_f0, 1);
it_alpha_store = 1;
for it_mcmc = 1:(nr_it_burnin + nr_it_mcmc),
    
    % Swapping regions
    event_swap = randperm(nr_events);
    event_swap = event_swap(1:2);    
    event_order_new = event_order_current;
    event_order_new(event_swap(1)) = event_order_current(event_swap(2));
    event_order_new(event_swap(2)) = event_order_current(event_swap(1));    
    log_likelihood_new = compute_loglikelihood(likelihood_events, ...
        event_order_new, f0_current, data_range);
    alpha = exp(log_likelihood_new - log_likelihood_current);
    if alpha > rand,
        
        event_order_current = event_order_new;
        log_likelihood_current = log_likelihood_new;
        
    end
    if log_likelihood_current > log_likelihood_max,
        
        event_order_max = event_order_current;
        f0_max = f0_current;
        log_likelihood_max = log_likelihood_current;
        
    end
    
    % Changing volume fractions
    f0_new = max(f0_range(1), min(f0_range(2), normrnd(f0_current, std_f0)));
    log_likelihood_new = compute_loglikelihood(likelihood_events, ...
        event_order_current, f0_new, data_range);
    alpha = exp(log_likelihood_new - log_likelihood_current);
    if alpha > rand,
        
        f0_current = f0_new;
        log_likelihood_current = log_likelihood_new;
        
    end
    
    if log_likelihood_current > log_likelihood_max,
        
        event_order_max = event_order_current;
        f0_max = f0_current;
        log_likelihood_max = log_likelihood_current;
        
    end

    alpha_store(it_alpha_store) = alpha;
    if it_alpha_store == nr_it_check_f0,
        
        it_alpha_store = 1;
        mean_alpha = mean(alpha_store);
        if mean_alpha < 0.15,

            std_f0 = std_f0 - 0.1*std_f0;

        elseif mean_alpha > 0.5,

            std_f0 = std_f0 + 0.1*std_f0;

        end
        alpha_store = zeros(nr_it_check_f0, 1);

    else
        
        it_alpha_store = it_alpha_store + 1;
        
    end

    if it_mcmc > nr_it_burnin,
        
        event_order_mcmc(:, it_mcmc-nr_it_burnin) = event_order_current;
        f0_mcmc(it_mcmc-nr_it_burnin) = f0_current;
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

parm_struct.log_likelihood_gradient_ascent = log_likelihood_gradient_ascent;
parm_struct.event_order_gradient_ascent = event_order_gradient_ascent;
parm_struct.event_order_mcmc = event_order_mcmc;
parm_struct.f0_mcmc = f0_mcmc;
parm_struct.log_likelihood_mcmc = log_likelihood_mcmc;
parm_struct.event_order_max = event_order_max;
parm_struct.f0_max = f0_max;
parm_struct.log_likelihood_max = log_likelihood_max;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function log_likelihood = compute_loglikelihood(likelihood_events, ...
    event_order, f0, data_range)

[nr_events, nr_pat] = size(likelihood_events(:, :, 1));
for event = 1:nr_events,
    
    likelihood_events(event, :, 1) = (f0/data_range(event)) + ...
        (1-f0)*likelihood_events(event, :, 1);
    likelihood_events(event, :, 2) = (f0/data_range(event)) + ...
        (1-f0)*likelihood_events(event, :, 2);
    
end

nr_phase = nr_events;
event_pattern = zeros(nr_events, nr_events + 1);
for event = 1:nr_events,
    
    event_pattern(event_order==event, (event+1):end) = 1;
    
end

events = event_pattern == 1;
noevents = event_pattern == 0;
likelihood = zeros(nr_phase+1, nr_pat);
for phase = 1:(nr_phase+1),
    
    likelihood(phase, :) = prod(likelihood_events(events(:, phase), :, 2), 1).*...
        prod(likelihood_events(noevents(:, phase), :, 1), 1);
    
end
log_likelihood = sum(log(sum(likelihood, 1)));
      
        