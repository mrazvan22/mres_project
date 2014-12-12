function [parm_struct, diag_struct] = ...
    AtrophyModelMCMC2b(p_A_D, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_like)

[nr_roi, nr_pat] = size(p_A_D);

% sum_p_atrophy = sum(p_atrophy, 2);
% [d, I_sort] = sort(sum_p_atrophy, 'descend');
% [d, order_events_current] = sort(I_sort);

logp_hillclimb = zeros(nr_it_hillclimb, nr_hillclimb);
order_events_hillclimb = zeros(nr_roi, nr_hillclimb);
logp_order_events_burnin = zeros(nr_it_burnin, 1);
order_events = zeros(nr_roi, nr_it_mcmc/thinning);
logp_order_events = zeros(nr_it_mcmc/thinning, 1);
cnt_it_mcmc = 1:100:(nr_it_mcmc + nr_it_burnin + nr_it_hillclimb);
cnt = 1;
for it_hillclimb = 1:nr_hillclimb,
    
    order_events_current = randperm(nr_roi);    
    logp_current = logLikelihoodAtrophyModel(p_A_D, ...
        order_events_current, version_like);
    
    for it = 1:nr_it_hillclimb,
        
        roi_swap = randperm(nr_roi);
        roi_swap = roi_swap(1:2);
        order_events_new = order_events_current;
        order_events_new(roi_swap(1)) = order_events_current(roi_swap(2));
        order_events_new(roi_swap(2)) = order_events_current(roi_swap(1));
        
        logp_new = logLik_localfunc(p_A_D, ...
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
for it_mcmc = 1:(nr_it_mcmc + nr_it_burnin),
    
    % Swapping regions around
    roi_swap = randperm(nr_roi);
    roi_swap = roi_swap(1:2);
    order_events_new = order_events_current;
    order_events_new(roi_swap(1)) = order_events_current(roi_swap(2));
    order_events_new(roi_swap(2)) = order_events_current(roi_swap(1));
             
    logp_new = logLik_localfunc(p_A_D, ...
        order_events_new, version_like);
    
    % Accept new ordering?
    alpha = exp(logp_new - logp_current);
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

[logp_max, logLik_max, model_bestfit] = logLik_localfunc(p_A_D, ...
    order_events_max, version_like);

parm_struct.order_events = order_events;
parm_struct.order_events_max = order_events_max;
parm_struct.logLik_max = logLik_max;
parm_struct.model_bestfit = model_bestfit;

diag_struct.logp_order_events = logp_order_events;
diag_struct.logp_order_events = logp_order_events;
diag_struct.logp_order_events_max = logp_order_events_max;
diag_struct.logp_hillclimb = logp_hillclimb;
diag_struct.order_events_hillclimb = order_events_hillclimb;
diag_struct.logp_order_events_burnin = logp_order_events_burnin;
    
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------


function [logLik, logLikmat, model_bestfit] = logLik_localfunc(p_A_D, ...
    order_events, version)

[nr_roi, nr_pat] = size(p_A_D);
atrophy_model = zeros(nr_roi, nr_roi);
for roi = 1:nr_roi,
    
    atrophy_model(order_events==roi, roi:end) = 1;
    
end

logp_A_D = log(p_A_D);
p_noA_D = 1- p_A_D;
p_noA_D(p_noA_D == 0) = eps;
logp_noA_D = log(p_noA_D);
atrophy_model = atrophy_model == 1;
noatrophy_model = atrophy_model == 0;

if ~exist('version', 'var'),
    
    version = 2;
    
end
if version ~= 3,
    
    logLikmat = zeros(nr_pat, nr_roi);
    for roi = 1:nr_roi,
        
        for pat = 1:nr_pat,
            
            logLikmat(pat, roi) = logLikmat(pat, roi) + ...
                sum(logp_A_D(atrophy_model(:, order_events(roi)), pat));
            logLikmat(pat, roi) = logLikmat(pat, roi) + ...
                sum(logp_noA_D(noatrophy_model(:, order_events(roi)), pat));
            
        end
        
    end
    I_pat = zeros(nr_pat, 1);
    for pat = 1:nr_pat,
        
        I_pat_local = find(logLikmat(pat, :) == max(logLikmat(pat, :)));
        I_pat(pat) = I_pat_local(1);
        
    end

end
if version == 1,
    
    logLikmat_filt = zeros(nr_pat, 1);
    for pat = 1:nr_pat,
        
        logLikmat_filt(pat) = logLikmat(pat, I_pat(pat));;
        
    end
    logLik = sum(sum(logLikmat_filt));
    model_bestfit = atrophy_model(:, order_events(I_pat));
    
elseif version == 2,

    logLik = sum(log(sum(exp(logLikmat), 2)));
    model_bestfit = atrophy_model(:, order_events(I_pat));
    
elseif version == 3,
    
    logLik_stage = zeros(nr_roi, 1);
    for stage = 1:nr_roi,
        
        
        logLik_stage(stage) = sum(prod(p_A_D(atrophy_model(:, stage), :), 1).*...
            prod(p_noA_D(noatrophy_model(:, stage), :), 1));
        
    end
    logLik = log(sum(logLik_stage));
    if nargout == 3,
        
        logLikmat = zeros(nr_pat, nr_roi);
        for roi = 1:nr_roi,
            
            for pat = 1:nr_pat,
                
                logLikmat(pat, roi) = logLikmat(pat, roi) + ...
                    sum(logp_A_D(atrophy_model(:, order_events(roi)), pat));
                logLikmat(pat, roi) = logLikmat(pat, roi) + ...
                    sum(logp_noA_D(noatrophy_model(:, order_events(roi)), pat));
                
            end
            
        end
        I_pat = zeros(nr_pat, 1);
        for pat = 1:nr_pat,
            
            I_pat_local = find(logLikmat(pat, :) == max(logLikmat(pat, :)));
            I_pat(pat) = I_pat_local(1);
            
        end
        model_bestfit = atrophy_model(:, order_events(I_pat));
        
    end
    
end      

