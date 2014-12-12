function [order_events, logp_order_events, logp_order_events_max, ...
    order_events_max, diagnostics] = AtrophyModelMCMC(p_atrophy, ...
    likelihoodVersion, nr_it_hillclimb, nr_it_burnin, nr_it_mcmc, thinning)

[nr_roi, nr_pat] = size(p_atrophy);
sum_p_atrophy = sum(p_atrophy, 2);
[d, I_sort] = sort(sum_p_atrophy, 'descend');
[d, order_events_current] = sort(I_sort);

% order_events_current = randperm(nr_roi);
atrophy_model_current = zeros(nr_roi, nr_roi);

for roi = 1:nr_roi,
    
    atrophy_model_current(order_events_current==roi, roi:end) = 1;
    
end

if likelihoodVersion ~= 5,
    
    eval(sprintf('logp_current = logLikelihoodAtrophyModel%d(p_atrophy, atrophy_model_current);', ...
        likelihoodVersion));
    
else
    
    c = zeros(nr_it_mcmc/thinning, 1);
    d = zeros(nr_it_mcmc/thinning, 1);
    logp_current = logLikelihoodAtrophyModel5(p_atrophy, atrophy_model_current, c);
    
end

order_events = zeros(nr_roi, nr_it_mcmc/thinning);
logp_order_events = zeros(nr_it_mcmc/thinning, 1);
if nargout == 5,
    
    logp_order_events_hillclimb = zeros(nr_it_hillclimb, 1);
    logp_order_events_burnin = zeros(nr_it_burnin, 1);
    
end

cnt_it_mcmc = 1:100:(nr_it_mcmc + nr_it_burnin + nr_it_hillclimb);
cnt = 1;

logp_order_events_max = logp_current;
order_events_max = order_events_current;
for it_mcmc = 1:(nr_it_mcmc + nr_it_burnin + nr_it_hillclimb),
    
    % Swapping regions around
    roi_swap = randperm(nr_roi);
    roi_swap = roi_swap(1:2);
    order_events_new = order_events_current;
    order_events_new(roi_swap(1)) = order_events_current(roi_swap(2));
    order_events_new(roi_swap(2)) = order_events_current(roi_swap(1));
    
    atrophy_model_new = zeros(nr_roi, nr_roi);
    for roi = 1:nr_roi,
        
        atrophy_model_new(order_events_new==roi, roi:end) = 1;
        
    end
    
    if likelihoodVersion ~= 5,
        eval(sprintf('logp_new = logLikelihoodAtrophyModel%d(p_atrophy, atrophy_model_new);', ...
            likelihoodVersion));
    end
    
    % Accept new ordering?
    if it_mcmc <= nr_it_hillclimb,
        
        % hill climbing ...
        alpha = (logp_new > logp_current);
        
    else
        
        alpha = exp(logp_new - logp_current);
        
    end
    if alpha > rand,
        
        order_events_current = order_events_new;
        logp_current = logp_new;
        if logp_current > logp_order_events_max,
            
            logp_order_events_max = logp_current;
            order_events_max = order_events_current;
            
        end
        
    end
    if (it_mcmc > (nr_it_burnin + nr_it_hillclimb)) && ~mod(it_mcmc, thinning),
        
        order_events(:, cnt) = order_events_current;
        logp_order_events(cnt) = logp_current;
        cnt = cnt + 1;
        
    end
    
    if it_mcmc <= nr_it_hillclimb,
        
        if nargout == 5,
            
            logp_order_events_hillclimb(it_mcmc) = logp_current;
            
        end
        if find(cnt_it_mcmc == it_mcmc),
            
            fprintf('Hillclimb it: %d in %d hillclimb its\n', ...
                it_mcmc, nr_it_hillclimb);
            
        end
        
    elseif it_mcmc <= (nr_it_hillclimb + nr_it_burnin),
        
        if nargout == 5,
            
            logp_order_events_burnin(it_mcmc) = logp_current;
            
        end
        if find(cnt_it_mcmc == it_mcmc),
            
            fprintf('Burnin it: %d in %d burnin its\n', ...
                it_mcmc - nr_it_hillclimb, nr_it_burnin);
            
        end
        
    elseif it_mcmc > (nr_it_burnin + nr_it_hillclimb),
        
        if find(cnt_it_mcmc == it_mcmc),
            
            fprintf('MCMC it: %d in %d mcmc its\n', ...
                it_mcmc-nr_it_burnin-nr_it_hillclimb, nr_it_mcmc)
            
        end
        
    end
    
end
if nargout == 5,
    
    diagnostics.logp_order_events_hillclimb = logp_order_events_hillclimb;
    diagnostics.logp_order_events_burnin = logp_order_events_burnin;
    
end
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------
    function logLik = logLikelihoodAtrophyModel1(p_atrophy, atrophy_model)
    
    [nr_roi, nr_pat] = size(p_atrophy);
    p_atrophy_model = zeros(nr_roi, nr_roi, nr_pat);
    for pat = 1:nr_pat,
        
        p_atrophy_local = p_atrophy(:, pat);
        p_atrophy_model_local = zeros(nr_roi, nr_roi);
        for roi = 1:nr_roi,
            
            I_atrophy = find(atrophy_model(:, roi) == 1);
            I_noatrophy = find(atrophy_model(:, roi) == 0);
            p_atrophy_model_local(I_atrophy, roi) = ...
                p_atrophy_local(I_atrophy);
            p_atrophy_model_local(I_noatrophy, roi) = ...
                1 - p_atrophy_local(I_noatrophy);
            
        end
        p_atrophy_model_local(p_atrophy_model_local == 0) = eps;
        p_atrophy_model(:, :, pat) = p_atrophy_model_local;
        
    end
    logLik = sum(log(sum(sum(p_atrophy_model, 1))));
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    function logLik = logLikelihoodAtrophyModel2(p_atrophy, atrophy_model)
        
        [nr_roi, nr_pat] = size(p_atrophy);
        p_atrophy_model = zeros(nr_roi, nr_roi, nr_pat);
        for pat = 1:nr_pat,
            
            p_atrophy_local = p_atrophy(:, pat);
            p_atrophy_model_local = zeros(nr_roi, nr_roi);
            for roi = 1:nr_roi,
                
                I_atrophy = find(atrophy_model(:, roi) == 1);
                I_noatrophy = find(atrophy_model(:, roi) == 0);
                p_atrophy_model_local(I_atrophy, roi) = ...
                    p_atrophy_local(I_atrophy);
                p_atrophy_model_local(I_noatrophy, roi) = ...
                    1 - p_atrophy_local(I_noatrophy);
                
            end
            p_atrophy_model_local(p_atrophy_model_local == 0) = eps;
            p_atrophy_model(:, :, pat) = p_atrophy_model_local;
            
        end
        logLik = sum(log(sum(prod(p_atrophy_model, 1))));
        % -------------------------------------------------------------------------
        % -------------------------------------------------------------------------
        
        function logLik = logLikelihoodAtrophyModel3(p_atrophy, atrophy_model)
            
            [nr_roi, nr_pat] = size(p_atrophy);
            p_atrophy_model = zeros(nr_roi, nr_pat);
            for pat = 1:nr_pat,
                
                p_atrophy_local = p_atrophy(:, pat);
                p_atrophy_model_local = zeros(nr_roi, nr_roi);
                for roi = 1:nr_roi,
                    
                    I_atrophy = find(atrophy_model(:, roi) == 1);
                    I_noatrophy = find(atrophy_model(:, roi) == 0);
                    p_atrophy_model_local(I_atrophy, roi) = ...
                        p_atrophy_local(I_atrophy);
                    p_atrophy_model_local(I_noatrophy, roi) = ...
                        1 - p_atrophy_local(I_noatrophy);
                    
                end
                p_atrophy_model_local(p_atrophy_model_local == 0) = eps;
                I_maxp = find(sum(p_atrophy_model_local) == max(sum(p_atrophy_model_local)));
                p_atrophy_model(:, pat) = p_atrophy_model_local(:, I_maxp(1));
                
            end
            logLik = sum(log(sum(p_atrophy_model, 1)));
            
            % -------------------------------------------------------------------------
            % -------------------------------------------------------------------------
            
            function logLik = logLikelihoodAtrophyModel4(p_atrophy, atrophy_model)
                
                [nr_roi, nr_pat] = size(p_atrophy);
                p_atrophy_model = zeros(nr_roi, nr_pat);
                for pat = 1:nr_pat,
                    
                    p_atrophy_local = p_atrophy(:, pat);
                    p_atrophy_model_local = zeros(nr_roi, nr_roi);
                    for roi = 1:nr_roi,
                        
                        I_atrophy = find(atrophy_model(:, roi) == 1);
                        I_noatrophy = find(atrophy_model(:, roi) == 0);
                        p_atrophy_model_local(I_atrophy, roi) = ...
                            p_atrophy_local(I_atrophy);
                        p_atrophy_model_local(I_noatrophy, roi) = ...
                            1 - p_atrophy_local(I_noatrophy);
                        
                    end
                    p_atrophy_model_local(p_atrophy_model_local == 0) = eps;
                    I_maxp = find(sum(log(p_atrophy_model_local)) == max(sum(log(p_atrophy_model_local))));
                    p_atrophy_model(:, pat) = p_atrophy_model_local(:, I_maxp(1));
                    
                end
                logLik = sum(log(p_atrophy_model(:)));
                
                % -------------------------------------------------------------------------
                % -------------------------------------------------------------------------
                
                % function logLik = logLikelihoodAtrophyModel5(p_atrophy, ...
                %     atrophy_model, c, d, p_thr)
                %
                % Xnm = p_atrophy > p_thr;
                % cc = corr(Xnm, atrophy_modl);
                %
                % [nr_roi, nr_pat] = size(p_atrophy);
                % p_atrophy_model = zeros(nr_roi, nr_pat);
                %
                % for pat = 1:nr_pat,
                %
                %     p_atrophy_local = p_atrophy(:, pat);
                %     p_atrophy_model_local = zeros(nr_roi, nr_roi);
                %     for roi = 1:nr_roi,
                %
                %         I_atrophy = find(atrophy_model(:, roi) == 1);
                %         I_noatrophy = find(atrophy_model(:, roi) == 0);
                %         p_atrophy_model_local(I_atrophy, roi) = ...
                %             p_atrophy_local(I_atrophy);
                %         p_atrophy_model_local(I_noatrophy, roi) = ...
                %             1 - p_atrophy_local(I_noatrophy);
                %
                %     end
                %     p_atrophy_model_local(p_atrophy_model_local == 0) = eps;
                %     I_maxp = find(sum(log(p_atrophy_model_local)) == max(sum(log(p_atrophy_model_local))));
                %     p_atrophy_model(:, pat) = p_atrophy_model_local(:, I_maxp(1));
                %
                % end
                % logLik = sum(log(p_atrophy_model(:)));
