function [parm_struct, diag_struct] = ...
    AtrophyModelMCMC2Corrected(prob_pat, prob_con, span, idx_clinevent, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_like)

nr_roi_org = size(prob_pat, 1);
I_zero = find(sum(prob_pat, 2) == 0);
I_nonzero = find(sum(prob_pat, 2) ~= 0);
prob_pat(I_zero, :) = [];
prob_con(I_zero, :) = [];
span(I_zero, :) = [];
idx_clinevent(I_zero) = [];
[nr_roi, nr_pat] = size(prob_pat);
% logf_uniform_span = [-10 log(0.2)];

I_clinevent = find(idx_clinevent == 1);
nr_clinevent = length(I_clinevent);
I_nonclinevent = find(idx_clinevent == 0);
nr_nonclinevent = length(I_nonclinevent);

sig_swap = 3;
% sig_f_uniform = 0.3*ones(2, 1);
diff_span = span(:, 2) - span(:, 1);

% sum_p_atrophy = sum(p_atrophy, 2);
% [d, I_sort] = sort(sum_p_atrophy, 'descend');
% [d, order_events_current] = sort(I_sort);

logp_hillclimb = zeros(nr_it_hillclimb, nr_hillclimb);
order_events_hillclimb = zeros(nr_roi, nr_hillclimb);
% logf_uniform_hillclimb = zeros(2, nr_hillclimb, 1);
logp_order_events_burnin = zeros(nr_it_burnin, 1);
order_events = zeros(nr_roi, nr_it_mcmc/thinning);
% logf_uniform = zeros(2, nr_it_mcmc/thinning, 1);
logp_order_events = zeros(nr_it_mcmc/thinning, 1);
cnt_it_mcmc = 1:100:(nr_it_mcmc + nr_it_burnin + nr_it_hillclimb);
accept_randomswap = 0.5;
% p_uniform = 1./repmat(diff_span, 1, nr_pat);
for it_hillclimb = 1:nr_hillclimb,
    
    order_events_current = randperm(nr_roi);    
    %     logf_uniform_current = (logf_uniform_span(2) - logf_uniform_span(1))*...
    %         rand(2, 1) + logf_uniform_span(1);
    %     logf_uniform_current = zeros(2, 1);

    %     f_local = zeros(nr_roi, nr_pat);
    %     f_local(I_nonclinevent, :) = repmat(exp(logf_uniform_current(1)), nr_nonclinevent, nr_pat);
    %     f_local(I_clinevent, :) = repmat(exp(logf_uniform_current(2)), nr_clinevent, nr_pat);
    %     prob_pat_current = ((1-f_local).*prob_pat + f_local.*p_uniform);
    %     prob_con_current = ((1-f_local).*prob_con + f_local.*p_uniform);

    logp_current = logLikelihoodAtrophyModelCorrect(prob_pat, ...
        prob_con, order_events_current);
    
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
            
        %         f_local = zeros(nr_roi, nr_pat);
        %         f_local(I_nonclinevent, :) = repmat(exp(logf_uniform_current(1)), nr_nonclinevent, nr_pat);
        %         f_local(I_clinevent, :) = repmat(exp(logf_uniform_current(2)), nr_clinevent, nr_pat);
        %         prob_pat_current = ((1-f_local).*prob_pat + f_local.*p_uniform);
        %         prob_con_current = ((1-f_local).*prob_con + f_local.*p_uniform);
        
        logp_new = logLikelihoodAtrophyModelCorrect(prob_pat, ...
            prob_con, order_events_new);

        if logp_new > logp_current,
            
            order_events_current = order_events_new;
            logp_current = logp_new;
            
        end
        %         for roi= 1:2,
        %
        %             flg_continue = 1;
        %             while flg_continue,
        %
        %                 logf_uniform_new = logf_uniform_current;
        %                 logf_uniform_new(roi) = normrnd(logf_uniform_new(roi), ...
        %                     sig_f_uniform(roi));
        %                 if (logf_uniform_new(roi) > logf_uniform_span(1)) && ...
        %                         (logf_uniform_new(roi) < logf_uniform_span(2)),
        %
        %                     flg_continue = 0;
        %
        %                 end
        %
        %             end
        %
        %             f_local = zeros(nr_roi, nr_pat);
        %             f_local(I_nonclinevent, :) = repmat(exp(logf_uniform_new(1)), nr_nonclinevent, nr_pat);
        %             f_local(I_clinevent, :) = repmat(exp(logf_uniform_new(2)), nr_clinevent, nr_pat);
        %             prob_pat_current = ((1-f_local).*prob_pat + f_local.*p_uniform);
        %             prob_con_current = ((1-f_local).*prob_con + f_local.*p_uniform);
        %
        %             logp_new = logLikelihoodAtrophyModelCorrect(prob_pat_current, ...
        %                 prob_con_current, order_events_current);
        %
        %             if logp_new > logp_current,
        %
        %                 logf_uniform_current = logf_uniform_new;
        %                 logp_current = logp_new;
        %
        %             end
        %
        %         end
        logp_hillclimb(it, it_hillclimb) = logp_current;
        
    end
    order_events_hillclimb(:, it_hillclimb) = order_events_current;
    %     logf_uniform_hillclimb(:, it_hillclimb) = logf_uniform_current;
    fprintf('hill-climb it: %d in %d iterations\n', it_hillclimb, nr_hillclimb);
    
end

it_hillclimb_max = find(logp_hillclimb(end, :) == max(logp_hillclimb(end, :)));
order_events_current = order_events_hillclimb(:, it_hillclimb_max(1));
% logf_uniform_current = logf_uniform_hillclimb(:, it_hillclimb_max(1));
logp_current = logp_hillclimb(end, it_hillclimb_max(1));
logp_order_events_max = logp_current;
order_events_max = order_events_current;
% logf_uniform_max = logf_uniform_current;
cnt = 1;
cnt_alpha = 1;
cnt_alpha_roi_swap = 1;
accept_randomswap = 0.3;
alpha_store = zeros(nr_roi+1, 100);
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
            
    %     f_local = zeros(nr_roi, nr_pat);
    %     f_local(I_nonclinevent, :) = repmat(exp(logf_uniform_current(1)), nr_nonclinevent, nr_pat);
    %     f_local(I_clinevent, :) = repmat(exp(logf_uniform_current(2)), nr_clinevent, nr_pat);
    %     p_A_D = ((1-f_local).*prob_pat + f_local.*p_uniform)./...
    %         ((1-f_local).*prob_pat + (1-f_local).*prob_con + 2*f_local.*p_uniform);

    logp_new = logLikelihoodAtrophyModelCorrect(prob_pat, prob_con, ...
        order_events_new);
    
    % Accept new ordering?
    alpha = exp(logp_new - logp_current);
    if flag_swap == 2,
        
        alpha_store(nr_roi+1, cnt_alpha_roi_swap) = alpha;
        cnt_alpha_roi_swap = cnt_alpha_roi_swap + 1;
        if cnt_alpha_roi_swap == 100,
            
            mean_alpha = mean(alpha_store(nr_roi + 1, :));
            if mean_alpha < 0.15,
                
                sig_swap = sig_swap - 0.1*sig_swap;
                
            elseif mean_alpha > 0.5,
                
                sig_swap = sig_swap + 0.1*sig_swap;
                
            end
            cnt_alpha_roi_swap = 1;
            alpha_store(nr_roi +1, :) = zeros(100, 1);

        end
        
    end
    if alpha > rand,
        
        order_events_current = order_events_new;
        logp_current = logp_new;
        if logp_current > logp_order_events_max,
            
            logp_order_events_max = logp_current;
            order_events_max = order_events_current;
            %             logf_uniform_max = logf_uniform_current;
            
        end
        
    end
    
    % Changing thresholds...
    %     for roi= 1:2,
    %
    %         flg_continue = 1;
    %         while flg_continue,
    %
    %             logf_uniform_new = logf_uniform_current;
    %             logf_uniform_new(roi) = normrnd(logf_uniform_new(roi), ...
    %                 sig_f_uniform(roi));
    %             if (logf_uniform_new(roi) > logf_uniform_span(1)) && ...
    %                     (logf_uniform_new(roi) < logf_uniform_span(2)),
    %
    %                 flg_continue = 0;
    %
    %             end
    %
    %         end
    %
    %         f_local = zeros(nr_roi, nr_pat);
    %         f_local(I_nonclinevent, :) = repmat(exp(logf_uniform_new(1)), nr_nonclinevent, nr_pat);
    %         f_local(I_clinevent, :) = repmat(exp(logf_uniform_new(2)), nr_clinevent, nr_pat);
    %         p_A_D = ((1-f_local).*prob_pat + f_local.*p_uniform)./...
    %             ((1-f_local).*prob_pat + (1-f_local).*prob_con + 2*f_local.*p_uniform);
    %
    %         logp_new = logLikelihoodAtrophyModel(p_A_D, ...
    %             order_events_current, version_like);
    %         alpha = exp(logp_new - logp_current);
    %         alpha_store(roi, cnt_alpha) = alpha;
    %         if alpha > rand,
    %
    %             logf_uniform_current = logf_uniform_new;
    %             logp_current = logp_new;
    %             if logp_current > logp_order_events_max,
    %
    %                 logp_order_events_max = logp_current;
    %                 order_events_max = order_events_current;
    %                 logf_uniform_max = logf_uniform_current;
    %
    %             end
    %
    %         end
    %
    %     end
    %     if cnt_alpha == 100,
    %
    %         for roi = 1:nr_roi,
    %
    %             mean_alpha = mean(alpha_store(roi, :));
    %             if mean_alpha < 0.15,
    %
    %                 sig_f_uniform(roi) = sig_f_uniform(roi) - 0.1*sig_f_uniform(roi);
    %
    %             elseif mean_alpha > 0.5,
    %
    %                 sig_f_uniform(roi) = sig_f_uniform(roi) + 0.1*sig_f_uniform(roi);
    %
    %             end
    %
    %         end
    %         cnt_alpha = 1;
    %         alpha_store(1:nr_roi, :) = zeros(nr_roi, 100);
    %
    %     else
    %
    %         cnt_alpha = cnt_alpha + 1;
    %
    %     end

    if (it_mcmc > nr_it_burnin) && ~mod(it_mcmc, thinning),
        
        order_events(:, cnt) = order_events_current;
        %         logf_uniform(:, cnt) = logf_uniform_current;
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

% f_local = zeros(nr_roi, nr_pat);
% f_local(I_nonclinevent, :) = repmat(exp(logf_uniform_max(1)), nr_nonclinevent, nr_pat);
% f_local(I_clinevent, :) = repmat(exp(logf_uniform_max(2)), nr_clinevent, nr_pat);
% p_A_D = ((1-f_local).*prob_pat + f_local.*p_uniform)./...
%     ((1-f_local).*prob_pat + (1-f_local).*prob_con +
%     2*f_local.*p_uniform);
[logp_max, logLik_max, model_bestfit] = logLikelihoodAtrophyModelCorrect(prob_pat, ...
    prob_con, order_events_max);

order_events_org = zeros(nr_roi_org, size(order_events, 2));
order_events_org(I_nonzero, :) = order_events;
order_events_org(I_zero, :) = nr_roi + 1;
order_events_max_org = zeros(nr_roi_org, 1);
order_events_max_org(I_nonzero) = order_events_max;
order_events_max_org(I_zero) = nr_roi+1;
order_events_hillclimb_org = zeros(nr_roi_org, size(order_events_hillclimb, 2));
order_events_hillclimb_org(I_nonzero, :) = order_events_hillclimb;
order_events_hillclimb_org(I_zero) = nr_roi+1;
% f_uniform_org = exp(logf_uniform);
% f_uniform_max_org = exp(logf_uniform_max);

parm_struct.order_events = order_events_org;
parm_struct.order_events_max = order_events_max_org;
% parm_struct.f_uniform = f_uniform_org;
% parm_struct.f_uniform_max = f_uniform_max_org;
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
