function [parm_struct, diag_struct] = ...
    AtrophyModelMCMC2c(p_A_D, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_like)

sig_swap = 3;

nr_roi_org = size(p_A_D, 1);
I_zero = find(sum(p_A_D, 2) == 0);
I_nonzero = find(sum(p_A_D, 2) ~= 0);
p_A_D(I_zero, :) = [];
[nr_roi, nr_pat] = size(p_A_D);

thr_p_lim = [1e-6 0.8];

% sum_p_atrophy = sum(p_atrophy, 2);
% [d, I_sort] = sort(sum_p_atrophy, 'descend');
% [d, order_events_current] = sort(I_sort);

logp_hillclimb = zeros(nr_it_hillclimb, nr_hillclimb);
order_events_hillclimb = zeros(nr_roi, nr_hillclimb);
thr_p_hillclimb = zeros(nr_hillclimb, 1);
logp_order_events_burnin = zeros(nr_it_burnin, 1);
order_events = zeros(nr_roi, nr_it_mcmc/thinning);
thr_p = zeros(nr_it_mcmc/thinning, 1);
logp_order_events = zeros(nr_it_mcmc/thinning, 1);
cnt_it_mcmc = 1:100:(nr_it_mcmc + nr_it_burnin + nr_it_hillclimb);
cnt = 1;
accept_randomswap = 0.5;
for it_hillclimb = 1:nr_hillclimb,
    
    order_events_current = randperm(nr_roi);    
    thr_p_current = (thr_p_lim(2) - thr_p_lim(1))*rand + thr_p_lim(1);
    logp_current = logLikelihoodAtrophyModel_thrp(p_A_D, ...
        order_events_current, thr_p_current, version_like);
    
    for it = 1:nr_it_hillclimb,
        
        choose_swap = rand;
        if choose_swap < accept_randomswap,
                        
            roi_swap = randperm(nr_roi);
            roi_swap = roi_swap(1:2);

        else
            
            roi_swap(1) = ceil(nr_roi*rand);
            pos_swap(1) = order_events_current(roi_swap(1));
            pos_swap(2) = ceil(normrnd(pos_swap(1), sig_swap));
            pos_swap(2) = min(nr_roi, max(1, pos_swap(2)));
            roi_swap(2) = find(order_events_current == pos_swap(2));
            
        end
        order_events_new = order_events_current;
        order_events_new(roi_swap(1)) = order_events_current(roi_swap(2));
        order_events_new(roi_swap(2)) = order_events_current(roi_swap(1));
        
        logp_new = logLikelihoodAtrophyModel(p_A_D, ...
            order_events_new, version_like);
    
        if logp_new > logp_current,
            
            order_events_current = order_events_new;
            logp_current = logp_new;
            
        end
        logp_hillclimb(it, it_hillclimb) = logp_current;
        
    end
    order_events_hillclimb(:, it_hillclimb) = order_events_current;
    fprintf('hill-climb it: %d in %d iterations\n', it_hillclimb, nr_hillclimb);
    
end

it_hillclimb_max = find(logp_hillclimb(end, :) == max(logp_hillclimb(end, :)));
order_events_current = order_events_hillclimb(:, it_hillclimb_max(1));
logp_current = logp_hillclimb(end, it_hillclimb_max(1));
logp_order_events_max = logp_current;
order_events_max = order_events_current;
cnt = 1;
cnt_alpha = 1;
accept_randomswap = 0.3;
alpha_store = zeros(100, 1);
sig_swap_store = sig_swap;
for it_mcmc = 1:(nr_it_mcmc + nr_it_burnin),
    
    % Swapping regions around
    choose_swap = rand;
    if choose_swap < accept_randomswap,
        
        roi_swap = randperm(nr_roi);
        roi_swap = roi_swap(1:2);
        flag_swap = 1;
        
    else
        
        flag_swap = 2;
        roi_swap(1) = ceil(nr_roi*rand);
        pos_swap(1) = order_events_current(roi_swap(1));
        pos_swap(2) = ceil(normrnd(pos_swap(1), sig_swap));
        pos_swap(2) = min(nr_roi, max(1, pos_swap(2)));
        roi_swap(2) = find(order_events_current == pos_swap(2));
        
    end
    order_events_new = order_events_current;
    order_events_new(roi_swap(1)) = order_events_current(roi_swap(2));
    order_events_new(roi_swap(2)) = order_events_current(roi_swap(1));
         
    logp_new = logLikelihoodAtrophyModel(p_A_D, ...
        order_events_new, version_like);
    
    % Accept new ordering?
    alpha = exp(logp_new - logp_current);
    if flag_swap == 2,
        
        alpha_store(cnt_alpha) = alpha;
        cnt_alpha = cnt_alpha + 1;
        if cnt_alpha > 100,
            
            mean_alpha = mean(alpha_store);
            if mean_alpha < 0.15,
                
                sig_swap = sig_swap - 0.1*sig_swap;
                
            elseif mean_alpha > 0.5,
                
                sig_swap = sig_swap + 0.1*sig_swap;
                
            end
            sig_swap_store = [sig_swap_store sig_swap];
            cnt_alpha = 1;
            alpha_store = zeros(100, 1);

        end
        
    end
    if alpha > rand,
        
        order_events_current = order_events_new;
        logp_current = logp_new;
        if logp_current > logp_order_events_max,
            
            logp_order_events_max = logp_current;
            order_events_max = order_events_current;
            
        end
        
    end
    if (it_mcmc > nr_it_burnin) && ~mod(it_mcmc, thinning),
        
        order_events(:, cnt) = order_events_current;
        logp_order_events(cnt) = logp_current;
        cnt = cnt + 1;
        
    end
    
      
    if it_mcmc <= nr_it_burnin,
        
        logp_order_events_burnin(it_mcmc) = logp_current;
        if find(cnt_it_mcmc == it_mcmc),
            
            fprintf('Burnin it: %d in %d burnin its\n', ...
                it_mcmc, nr_it_burnin);
            
        end
        
    elseif it_mcmc > nr_it_burnin,
        
        if find(cnt_it_mcmc == it_mcmc),
            
            fprintf('MCMC it: %d in %d mcmc its\n', ...
                it_mcmc-nr_it_burnin, nr_it_mcmc)
            
        end
        
    end
    
end

[logp_max, logLik_max, model_bestfit] = logLikelihoodAtrophyModel(p_A_D, ...
    order_events_max, version_like);

order_events_org = zeros(nr_roi_org, size(order_events, 2));
order_events_org(I_nonzero, :) = order_events;
order_events_org(I_zero, :) = nr_roi + 1;
order_events_max_org = zeros(nr_roi_org, 1);
order_events_max_org(I_nonzero) = order_events_max;
order_events_max_org(I_zero) = nr_roi+1;
order_events_hillclimb_org = zeros(nr_roi_org, size(order_events_hillclimb, 2));
order_events_hillclimb_org(I_nonzero, :) = order_events_hillclimb;
order_events_hillclimb_org(I_zero) = nr_roi+1;

parm_struct.order_events = order_events_org;
parm_struct.order_events_max = order_events_max_org;
parm_struct.logLik_max = logLik_max;
parm_struct.model_bestfit = model_bestfit;

diag_struct.logp_order_events = logp_order_events;
diag_struct.logp_order_events = logp_order_events;
diag_struct.logp_order_events_max = logp_order_events_max;
diag_struct.logp_hillclimb = logp_hillclimb;
diag_struct.order_events_hillclimb = order_events_hillclimb_org;
diag_struct.logp_order_events_burnin = logp_order_events_burnin;
    
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------
