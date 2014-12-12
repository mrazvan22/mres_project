function [parm_struct, diag_struct] = ...
    EventOnsetModel(X, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb)

nr_event_org = size(X, 1);
I_zero = find(sum(X, 2) == 0);
I_nonzero = find(sum(X, 2) ~= 0);
X(I_zero, :) = [];
[nr_event, nr_pat] = size(X);

% sum_p_atrophy = sum(p_atrophy, 2);
% [d, I_sort] = sort(sum_p_atrophy, 'descend');
% [d, order_events_current] = sort(I_sort);

logp_hillclimb = zeros(nr_it_hillclimb, nr_hillclimb);
mu_hillclimb = zeros(nr_event, nr_hillclimb);
sigma_hillclimb = zeros(nr_event, nr_hillclimb);
p_fpfn_hillclimb = zeros(2, nr_hillclimb); %fpfn = false positive false negative
logp_burnin = zeros(nr_it_burnin, 1);
mu = zeros(nr_event, nr_it_mcmc/thinning);
sigma = zeros(nr_event, nr_it_mcmc/thinning);
p_fpfn = zeros(2, nr_it_mcmc/thinning);
logp_mcmc = zeros(nr_it_mcmc/thinning, 1);
cnt_it_mcmc = 1:100:(nr_it_mcmc + nr_it_burnin + nr_it_hillclimb);
cnt = 1;
accept_randomswap = 0.5;
sig_swap = 1;
std_jump_sigma = 0.5*ones(nr_event, 1);
std_jump_p_fpfn = 0.03;

sigma_range = [1e-6 50];
p_fpfn_range = [1e-6 0.3];
for it_hillclimb = 1:nr_hillclimb,
    
    mu_current = randperm(nr_event)';
    sigma_current = 0.5*rand(nr_event, 1);
    p_fpfn_current = 0.1*rand(1, 2);
    logp_current = logLikFun(X, mu_current, sigma_current, ...
        p_fpfn_current);
    

    for it = 1:nr_it_hillclimb,
        
        choose_swap = rand;
        if choose_swap < accept_randomswap,
            
            roi_swap = randperm(nr_event);
            roi_swap = roi_swap(1:2);
            flag_swap = 1;
            
        else
            
            flag_swap = 2;
            roi_swap(1) = ceil(nr_event*rand);
            pos_swap(1) = mu_current(roi_swap(1));
            pos_swap(2) = ceil(normrnd(pos_swap(1), sig_swap));
            pos_swap(2) = min(nr_event, max(1, pos_swap(2)));
            roi_swap(2) = find(mu_current == pos_swap(2));
            
        end
        mu_new = mu_current;
        mu_new(roi_swap(1)) = mu_current(roi_swap(2));
        mu_new(roi_swap(2)) = mu_current(roi_swap(1));
        
        logp_new = logLikFun(X, mu_new, sigma_current, ...
            p_fpfn_current);
        if logp_new > logp_current,
            
            mu_current = mu_new;
            logp_current = logp_new;
            
        end
        
        for e = 1:nr_event,
                       
            flg_continue = 1;
            while flg_continue,
                
                sigma_new = sigma_current;
                sigma_new(e) = normrnd(sigma_new(e), std_jump_sigma(e));
                if (sigma_new(e) > sigma_range(1)) && ...
                        (sigma_new(e) < sigma_range(2)),
                    
                    flg_continue = 0;
                    
                end
                
            end
            logp_new = logLikFun(X, mu_current, sigma_new, ...
                p_fpfn_current);
            if logp_new > logp_current,
                
                sigma_current = sigma_new;
                logp_current = logp_new;
                
            end
            
        end
            
        flg_continue = 1;
        while flg_continue,
            
            p_fpfn_new = normrnd(p_fpfn_current, std_jump_p_fpfn);
            if (min(p_fpfn_new) > p_fpfn_range(1)) && ...
                    (max(p_fpfn_new) < p_fpfn_range(2)),
                
                flg_continue = 0;
                
            end
            
        end
        logp_new = logLikFun(X, mu_current, sigma_current, ...
            p_fpfn_new);
        if logp_new > logp_current,
            
            p_fpfn_current = p_fpfn_new;
            logp_current = logp_new;
            
        end
        logp_hillclimb(it, it_hillclimb) = logp_current;
        %         mu_current
        %         sigma_current
        %         p_fpfn_current
        %         logp_current
        %         pause
        
    end
    mu_hillclimb(:, it_hillclimb) = mu_current;
    sigma_hillclimb(:, it_hillclimb) = sigma_current;
    p_fpfn_hillclimb(:, it_hillclimb) = p_fpfn_current;
    fprintf('hill-climb it: %d in %d iterations\n', it_hillclimb, nr_hillclimb);
    
end

it_hillclimb_max = find(logp_hillclimb(end, :) == max(logp_hillclimb(end, :)));
mu_current = mu_hillclimb(:, it_hillclimb_max(1));
sigma_current = sigma_hillclimb(:, it_hillclimb_max(1));
p_fpfn_current = p_fpfn_hillclimb(:, it_hillclimb_max(1));
logp_current = logp_hillclimb(end, it_hillclimb_max(1));
logp_max = logp_current;
mu_max = mu_current;
sigma_max = sigma_current;
p_fpfn_max = p_fpfn_current;
cnt = 1;
cnt_alpha = 1;
accept_randomswap = 0.3;
alpha_store = zeros(100, 1);
sig_swap_store = sig_swap;
for it_mcmc = 1:(nr_it_mcmc + nr_it_burnin),
    
    choose_swap = rand;
    if choose_swap < accept_randomswap,
        
        roi_swap = randperm(nr_event);
        roi_swap = roi_swap(1:2);
        flag_swap = 1;
        
    else
        
        flag_swap = 2;
        roi_swap(1) = ceil(nr_event*rand);
        pos_swap(1) = mu_current(roi_swap(1));
        pos_swap(2) = ceil(normrnd(pos_swap(1), sig_swap));
        pos_swap(2) = min(nr_event, max(1, pos_swap(2)));
        roi_swap(2) = find(mu_current == pos_swap(2));
        
    end
    mu_new = mu_current;
    mu_new(roi_swap(1)) = mu_current(roi_swap(2));
    mu_new(roi_swap(2)) = mu_current(roi_swap(1));
    
    logp_new = logLikFun(X, mu_new, sigma_current, ...
        p_fpfn_current);
    
    % Accept new ordering?
    alpha = exp(logp_new - logp_current);
    if alpha > rand,
        
        mu_current = mu_new;
        logp_current = logp_new;
        if logp_current > logp_max,
            
            logp_max = logp_current;
            mu_max = mu_current;
            
        end
        
    end
    
    for e = 1:nr_event,
        
        flg_continue = 1;
        while flg_continue,
            
            sigma_new = sigma_current;
            sigma_new(e) = normrnd(sigma_new(e), std_jump_sigma(e));
            if (sigma_new(e) > sigma_range(1)) && ...
                    (sigma_new(e) < sigma_range(2)),
                
                flg_continue = 0;
                
            end
            
        end
        logp_new = logLikFun(X, mu_current, sigma_new, ...
            p_fpfn_current);
        
        alpha = exp(logp_new - logp_current);
        if alpha > rand,
            
            sigma_current = sigma_new;
            logp_current = logp_new;
            if logp_current > logp_max,
                
                logp_max = logp_current;
                sigma_max = sigma_current;
                
            end
            
        end

        if logp_new > logp_current,
            
            sigma_current = sigma_new;
            logp_current = logp_new;
            
        end
        
    end
    for e = 1:2,
        
        flg_continue = 1;
        while flg_continue,
            
            p_fpfn_new = p_fpfn_current;
            p_fpfn_new(e) = normrnd(p_fpfn_new(e), std_jump_p_fpfn);
            if (p_fpfn_new(e) > p_fpfn_range(1)) && ...
                    (p_fpfn_new(e) < p_fpfn_range(2)),
                
                flg_continue = 0;
                
            end
            
        end
        logp_new = logLikFun(X, mu_current, sigma_current, ...
            p_fpfn_new);
        alpha = exp(logp_new - logp_current);
        if alpha > rand,
            
            p_fpfn_current = p_fpfn_new;
            logp_current = logp_new;
            if logp_current > logp_max,
                
                logp_max = logp_current;
                p_fpfn_max = p_fpfn_current;
                
            end
            
        end
        
    end


    if (it_mcmc > nr_it_burnin) && ~mod(it_mcmc, thinning),
        
        mu(:, cnt) = mu_current;
        sigma(:, cnt) = sigma_current;
        p_fpfn(:, cnt) = p_fpfn_current;
        logp_mcmc(cnt) = logp_current;
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

mu_org = zeros(nr_event_org, size(mu, 2));
mu_org(I_nonzero, :) = mu;
mu_org(I_zero, :) = nr_event + 1;
mu_max_org = zeros(nr_event_org, 1);
mu_max_org(I_nonzero) = mu_max;
mu_max_org(I_zero) = nr_event+1;
mu_hillclimb_org = zeros(nr_event_org, size(mu_hillclimb, 2));
mu_hillclimb_org(I_nonzero, :) = mu_hillclimb;
mu_hillclimb_org(I_zero) = nr_event+1;
sigma_org = zeros(nr_event_org, size(sigma, 2));
sigma_org(I_nonzero, :) = sigma;
sigma_org(I_zero, :) = mean(sigma(:));
sigma_max_org = zeros(nr_event_org, 1);
sigma_max_org(I_nonzero) = sigma_max;
sigma_max_org(I_zero) = mean(sigma_max);
sigma_hillclimb_org = zeros(nr_event_org, size(sigma_hillclimb, 2));
sigma_hillclimb_org(I_nonzero, :) = sigma_hillclimb;
sigma_hillclimb_org(I_zero) = mean(sigma_hillclimb(:));

nr_stage = nr_event_org + 1;
stage_mat = repmat([1:nr_stage], nr_event_org, 1);
mu_mat = repmat(mean(mu_org, 2), 1, nr_stage);
sigma_mat = repmat(mean(sigma_org, 2), 1, nr_stage);
model_progress_mean = normpdf(stage_mat, mu_mat, sigma_mat);
mu_mat = repmat(mu_max_org, 1, nr_stage);
sigma_mat = repmat(sigma_max_org, 1, nr_stage);
model_progress_max = normpdf(stage_mat, mu_mat, sigma_mat);

parm_struct.mu = mu_org;
parm_struct.mu_max = mu_max_org;
parm_struct.sigma = sigma_org;
parm_struct.sigma_max = sigma_max_org;
parm_struct.p_fpfn = p_fpfn;
parm_struct.model_progress_mean = model_progress_mean;
parm_struct.model_progress_max = model_progress_max;

diag_struct.logp_mcmc = logp_mcmc;
diag_struct.logp_max = logp_max;
diag_struct.logp_hillclimb = logp_hillclimb;
diag_struct.mu_hillclimb = mu_hillclimb_org;
diag_struct.sigma_hillclimb = sigma_hillclimb_org;
diag_struct.p_fpfn_hillclimb = p_fpfn_hillclimb;
diag_struct.logp_order_events_burnin = logp_order_events_burnin;
    
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [logLik, model_progress] = logLikFun(X, mu, sigma, p_fpfn)

c = p_fpfn(1);
d = p_fpfn(2);
[nr_event, nr_pat] = size(X);
nr_stage = nr_event + 1;
likmat = zeros(nr_event, nr_stage, nr_pat);
p_AgivenX = zeros(nr_event, nr_pat);
p_AgivenX(X==1) = 1-c;
p_AgivenX(X==0) = d;
p_noAgivenX = zeros(nr_event, nr_pat);
p_noAgivenX(X==1) = c;
p_noAgivenX(X==0) = 1-d;

stage_mat = repmat([1:nr_stage], nr_event, 1);
mu_mat = repmat(mu, 1, nr_stage);
sigma_mat = repmat(sigma, 1, nr_stage);
model_progress = normcdf(stage_mat, mu_mat, sigma_mat);
model_progress(model_progress < 0) = 0;
model_progress(model_progress > 1) = 1;

for pat = 1:nr_pat,
    
    likmat(:, :, pat) = model_progress.*repmat(p_AgivenX(:, pat), 1, nr_stage) + ...
        (1-model_progress).*repmat(p_noAgivenX(:, pat), 1, nr_stage);
    
end
logLik = sum(log(sum(prod(likmat, 1), 2)));
    