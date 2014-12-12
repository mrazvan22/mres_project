function [parm_struct, diag_struct] = ...
    AtrophyModelMCMC3(p_A_D, thr, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb)

[nr_roi, nr_pat] = size(p_A_D);
p_A_D = p_A_D > thr;
sig_sig_error = 0.03;

logp_hillclimb = zeros(nr_it_hillclimb, nr_hillclimb);
order_events_hillclimb = zeros(nr_roi, nr_hillclimb);
sig_error_hillclimb = zeros(nr_hillclimb);
logp_order_events_burnin = zeros(nr_it_burnin, 1);
order_events = zeros(nr_roi, nr_it_mcmc/thinning);
sig_error = zeros(nr_it_mcmc/thinning, 1);
logp_order_events = zeros(nr_it_mcmc/thinning, 1);
cnt_it_mcmc = 1:100:(nr_it_mcmc + nr_it_burnin + nr_it_hillclimb);
cnt = 1;
for it_hillclimb = 1:nr_hillclimb,
    
    order_events_current = randperm(nr_roi);
    sig_error_current = 3*rand;
    
    logp_current = mcmc_it(p_A_D, order_events_current, sig_error_current);
    for it = 1:nr_it_hillclimb,
        
        roi_swap = randperm(nr_roi);
        roi_swap = roi_swap(1:2);
        order_events_new = order_events_current;
        order_events_new(roi_swap(1)) = order_events_current(roi_swap(2));
        order_events_new(roi_swap(2)) = order_events_current(roi_swap(1));
        
        logp_new = mcmc_it(p_A_D, order_events_new, sig_error_current);    
        if logp_new > logp_current,
            
            order_events_current = order_events_new;
            logp_current = logp_new;
            
        end
        flg_continue = 1;
        while flg_continue
            
            sig_error_new = normrnd(sig_error_current, sig_sig_error);
            if (sig_error_new > 0),
                
                flg_continue = 0;
                
            end
            
        end
        logp_new = mcmc_it(p_A_D, order_events_current, sig_error_new);
        if logp_new > logp_current,
            
            sig_error_current = sig_error_new;
            logp_current = logp_new;
            
        end
        logp_hillclimb(it, it_hillclimb) = logp_current;
        
    end
    order_events_hillclimb(:, it_hillclimb) = order_events_current;
    sig_error_hillclimb(it_hillclimb) = sig_error_current;
    fprintf('hill-climb it: %d in %d iterations\n', it_hillclimb, nr_hillclimb);
    
end

it_hillclimb_max = find(logp_hillclimb(end, :) == max(logp_hillclimb(end, :)));
order_events_current = order_events_hillclimb(:, it_hillclimb_max(1));
sig_error_current = sig_error_hillclimb(it_hillclimb_max(1));
logp_current = logp_hillclimb(end, it_hillclimb_max(1));
logp_order_events_max = logp_current;
order_events_max = order_events_current;
sig_error_max = sig_error_current;
for it_mcmc = 1:(nr_it_mcmc + nr_it_burnin),
    
    % Swapping regions around
    roi_swap(1) = ceil(nr_roi*rand);
    order_events_new = order_events_current;
    pos_swap(1) = order_events_new(roi_swap(1));
    flg_continue = 1;
    while flg_continue
        
        pos_swap(2) = ceil(normrnd(pos_swap(1), 3));
        if (pos_swap(2) > 0) & (pos_swap(2) <= nr_roi),
            
            flg_continue = 0;
            
        end
        
    end
    roi_swap(2) = find(order_events_new == pos_swap(2));
    order_events_new(roi_swap(1)) = order_events_current(roi_swap(2));
    order_events_new(roi_swap(2)) = order_events_current(roi_swap(1));
             
    logp_new = mcmc_it(p_A_D, order_events_new, sig_error_current);
    % Accept new ordering?
    alpha = exp(logp_new - logp_current);
    if alpha > rand,
        
        order_events_current = order_events_new;
        logp_current = logp_new;
        if logp_current > logp_order_events_max,
            
            logp_order_events_max = logp_current;            
            order_events_max = order_events_current;
            sig_error_max = sig_error_current;
            
        end
        
    end
    
    flg_continue = 1;
    while flg_continue
        
        sig_error_new = normrnd(sig_error_current, sig_sig_error);
        if (sig_error_new > 0),
            
            flg_continue = 0;
            
        end
        
    end
    [logp_new, pat_sim] = mcmc_it(p_A_D, order_events_current, sig_error_new);
    alpha = exp(logp_new - logp_current);
    if alpha > rand,
        
        sig_error_current = sig_error_new;
        logp_current = logp_new;
        if logp_current > logp_order_events_max,
            
            logp_order_events_max = logp_current;
            order_events_max = order_events_current;
            sig_error_max = sig_error_current;            
            
        end
        
    end
    
    if (it_mcmc > nr_it_burnin) && ~mod(it_mcmc, thinning),
        
        order_events(:, cnt) = order_events_current;
        sig_error(cnt) = sig_error_current;
        logp_order_events(cnt) = logp_current;

        cnt = cnt + 1;
        
    end
    
      
    if it_mcmc <= nr_it_burnin,
        
        logp_order_events_burnin(it_mcmc) = logp_current;
        if find(cnt_it_mcmc == it_mcmc),
            
            fprintf('Burnin it: %d in %d burnin its logLik: %3.3f\n', ...
                it_mcmc, nr_it_burnin, logp_current);            
            
        end
        
    elseif it_mcmc > nr_it_burnin,
        
        if find(cnt_it_mcmc == it_mcmc),
            
            fprintf('MCMC it: %d in %d mcmc its logLik: %3.3f\n', ...
                it_mcmc-nr_it_burnin, nr_it_mcmc, logp_current)
            
        end
        
    end
    
end

[logp_dummy, model_bestfit] = mcmc_it(p_A_D, order_events_max, sig_error_max);

parm_struct.order_events = order_events;
parm_struct.order_events_max = order_events_max;
parm_struct.sig_error = sig_error;
parm_struct.sig_error_max = sig_error_max;
parm_struct.model_bestfit = model_bestfit;

diag_struct.logp_order_events = logp_order_events;
diag_struct.logp_order_events = logp_order_events;
diag_struct.logp_order_events_max = logp_order_events_max;
diag_struct.logp_hillclimb = logp_hillclimb;
diag_struct.order_events_hillclimb = order_events_hillclimb;
diag_struct.sig_error_hillblimc = sig_error_hillclimb;
diag_struct.logp_order_events_burnin = logp_order_events_burnin;
    
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [logp_new, pat_sim] = mcmc_it(p_A_D, order_events, sig_error)

[nr_roi, nr_pat] = size(p_A_D);
atrophy_model = zeros(nr_roi, nr_roi);
for roi = 1:nr_roi,
    
    atrophy_model(order_events==roi, roi:end) = 1;
    
end

overlap_mat = zeros(nr_pat, nr_roi);
for pat = 1:nr_pat,
    
    for roi = 1:nr_roi,
        
        overlap_mat(pat, roi) = sum(p_A_D(:, pat) == atrophy_model(:, order_events(roi)))/nr_roi;
        
    end
    
end
I_pat = zeros(nr_pat, 1);
for pat = 1:nr_pat,
    
    I_pat_local = find(overlap_mat(pat, :) == max(overlap_mat(pat, :)));
    I_pat(pat) = I_pat_local(1);
    
end
pat_sim = atrophy_model(:, order_events(I_pat));
error = p_A_D - pat_sim;
logp_new = logp_Npdf(error(:), sig_error);
    
    
