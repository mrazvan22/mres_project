function [parm_struct, diag_struct] = ...
    AtrophyModelMCMCPuolamaki(atrophy_data, c_lim, ...
    d_lim, nr_it_hillclimb, nr_it_burnin, nr_it_mcmc, ...
    thinning, nr_hillclimb)

nr_roi_org = size(atrophy_data, 1);
I_zero = find(sum(atrophy_data, 2) == 0);
I_nonzero = find(sum(atrophy_data, 2) ~= 0);
atrophy_data(I_zero, :) = [];
[nr_roi, nr_pat] = size(atrophy_data);
order_events_current = randperm(nr_roi);
% sum_p_atrophy = sum(p_atrophy, 2);
% [d, I_sort] = sort(sum_p_atrophy, 'descend');
% [d, order_events_current] = sort(I_sort);

logp_hillclimb = zeros(nr_it_hillclimb, nr_hillclimb);
order_events_hillclimb = zeros(nr_roi, nr_hillclimb);
c_hillclimb = zeros(nr_hillclimb, 1);
d_hillclimb = zeros(nr_hillclimb, 1);
logp_order_events_burnin = zeros(nr_it_burnin, 1);
c = zeros(nr_it_mcmc/thinning, 1);
d = zeros(nr_it_mcmc/thinning, 1);
order_events = zeros(nr_roi, nr_it_mcmc/thinning);
logp_order_events = zeros(nr_it_mcmc/thinning, 1);
cnt_it_mcmc = 1:100:(nr_it_mcmc + nr_it_burnin + nr_it_hillclimb);
cnt = 1;
for it_hillclimb = 1:nr_hillclimb,
        
    c_current = 0.1;
    d_current = 0.1;
    order_events_current = randperm(nr_roi);
    logp_current = logLikelihoodAtrophyModelPuolamaki(atrophy_data, ...
        order_events_current, c_current, d_current);
    
    for it = 1:nr_it_hillclimb,
        
        roi_swap = randperm(nr_roi);
        roi_swap = roi_swap(1:2);
        order_events_new = order_events_current;
        order_events_new(roi_swap(1)) = order_events_current(roi_swap(2));
        order_events_new(roi_swap(2)) = order_events_current(roi_swap(1));
        %         c_new = normrnd(c_current, 0.03);
        %         d_new = normrnd(d_current, 0.03);
        %         if (c_new < c_lim(1)) || (c_new > c_lim(2)),
        %
        %             c_new = c_current;
        %
        %         end
        %         if (d_new < d_lim(1)) || (d_new > d_lim(2)),
        %
        %             d_new = d_current;
        %
        %         end
        c_new = c_current;
        d_new = d_current;
        logp_new = logLikelihoodAtrophyModelPuolamaki(atrophy_data, ...
            order_events_new, c_new, d_new);
    
        if logp_new > logp_current,
            
            order_events_current = order_events_new;
            c_current = c_new;
            d_current = d_new;
            logp_current = logp_new;
            
        end
        logp_hillclimb(it, it_hillclimb) = logp_current;
        
    end
    order_events_hillclimb(:, it_hillclimb) = order_events_current;
    c_hillclimb(it_hillclimb) = c_current;
    d_hillclimb(it_hillclimb) = d_current;
    fprintf('hill-climb it: %d in %d iterations\n', it_hillclimb, nr_hillclimb);
    
end

it_hillclimb_max = find(logp_hillclimb(end, :) == max(logp_hillclimb(end, :)));
order_events_current = order_events_hillclimb(:, it_hillclimb_max(1));
logp_current = logp_hillclimb(end, it_hillclimb_max(1));
logp_order_events_max = logp_current;
order_events_max = order_events_current;
c_max = c_current;
d_max = d_current;
for it_mcmc = 1:(nr_it_mcmc + nr_it_burnin),
    
    % Swapping regions around
    roi_swap = randperm(nr_roi);
    roi_swap = roi_swap(1:2);
    order_events_new = order_events_current;
    order_events_new(roi_swap(1)) = order_events_current(roi_swap(2));
    order_events_new(roi_swap(2)) = order_events_current(roi_swap(1));
    logp_new = logLikelihoodAtrophyModelPuolamaki(atrophy_data, ...
        order_events_new, c_current, d_current);
    
    alpha = exp(logp_new - logp_current);
    if alpha > rand,
        
        order_events_current = order_events_new;
        if logp_current > logp_order_events_max,
            
            logp_order_events_max = logp_current;
            order_events_max = order_events_current;
            c_max = c_current;
            d_max = d_current;
            
        end
        
    end
        
    c_new = normrnd(c_current, 0.03);
    if (c_new < c_lim(1)) || (c_new > c_lim(2)),
        
        c_new = c_current;
        
    end
    d_new = normrnd(d_current, 0.03);
    if (d_new < d_lim(1)) || (d_new > d_lim(2)),
        
        d_new = d_current;
        
    end
    logp_new = logLikelihoodAtrophyModelPuolamaki(atrophy_data, ...
        order_events_current, c_new, d_new);

    % Accept new ordering?
    alpha = exp(logp_new - logp_current);
    if alpha > rand,
        
        c_current = c_new;
        d_current = d_new;
        logp_current = logp_new;
        if logp_current > logp_order_events_max,
            
            logp_order_events_max = logp_current;
            order_events_max = order_events_current;
            c_max = c_current;
            d_max = d_current;
            
        end
        
    end
    if (it_mcmc > nr_it_burnin) && ~mod(it_mcmc, thinning),
        
        order_events(:, cnt) = order_events_current;
        logp_order_events(cnt) = logp_current;
        c(cnt) = c_current;
        d(cnt) = d_current;
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

[logp_new, model_bestfit] = logLikelihoodAtrophyModelPuolamaki(atrophy_data, ...
    order_events_max, c_max, d_max);

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
parm_struct.c = c;
parm_struct.d = d;
parm_struct.c_max = c_max;
parm_struct.d_max = d_max;
parm_struct.model_bestfit = model_bestfit;

diag_struct.logp_order_events = logp_order_events;
diag_struct.logp_order_events = logp_order_events;
diag_struct.logp_order_events_max = logp_order_events_max;
diag_struct.logp_hillclimb = logp_hillclimb;
diag_struct.order_events_hillclimb = order_events_hillclimb_org;
diag_struct.c_hillclimb = c_hillclimb;
diag_struct.d_hillclimb = d_hillclimb;
diag_struct.logp_order_events_burnin = logp_order_events_burnin;
    
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------
function [logLik, model_bestfit] = logLikelihoodAtrophyModelPuolamaki(data_atrophy, ...
    order_events, c, d)

[nr_roi, nr_pat] = size(data_atrophy);
atrophy_model = zeros(nr_roi, nr_roi);
for roi = 1:nr_roi,
    
    atrophy_model(order_events==roi, roi:end) = 1;
    
end

% logLikmat = zeros(nr_roi, nr_pat);
logLikmat = zeros(nr_pat, 1);

I_pat = zeros(nr_pat, 1);
overlap_mat = zeros(nr_pat, nr_roi);
for pat = 1:nr_pat,
    
    for roi = 1:nr_roi,
        
        overlap_mat(pat, roi) = sum(data_atrophy(:, pat) == atrophy_model(:, order_events(roi)))/(nr_roi);
        
    end
    
end

for pat = 1:nr_pat,
    
    I_roi = find(overlap_mat(pat, :) == max(overlap_mat(pat, :)));
    I_roi = I_roi(1);
    I_pat(pat) = I_roi;
    I_alpha = find((data_atrophy(:, pat) == 1) & ...
        (atrophy_model(:, order_events(I_roi)) == 0));
    p_alpha = length(I_alpha)*log(c);
    I_beta = find((data_atrophy(:, pat) == 1) & ...
        (atrophy_model(:, order_events(I_roi)) == 1));
    p_beta = length(I_beta)*log(1-c);
    I_gamma = find((data_atrophy(:, pat) == 0) & ...
        (atrophy_model(:, order_events(I_roi)) == 1));
    p_gamma = length(I_gamma)*log(d);
    I_delta = find((data_atrophy(:, pat) == 0) & ...
        (atrophy_model(:, order_events(I_roi)) == 0));
    p_delta = length(I_delta)*log(1-d);

    logLikmat(pat) = p_alpha + p_beta + p_gamma + p_delta;
    
end
% logLik = sum(log(sum(exp(logLikmat), 1)));
logLik = sum(logLikmat);
model_bestfit = atrophy_model(:, order_events(I_pat));

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
